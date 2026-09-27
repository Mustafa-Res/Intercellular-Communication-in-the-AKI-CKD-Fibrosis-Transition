"""
Pillar-1 EV-specificity enrichment engine.

Provides the shared machinery used by Analyses A/B, the negative-control ladder,
and the covariate-matched permutation test:

  - run_gsea_prerank : GSEA on a genome-wide ranked log2FC vector over a gene set
                       (PRIMARY enrichment; governs Gate 1). Wraps gseapy.prerank.
  - run_ora          : hypergeometric ORA over a DEG set (SECONDARY/supportive).
  - matched_permutation_pvalue : empirical p that EV-set enrichment exceeds random
                       non-EV gene sets matched on expression level, gene length,
                       and detection rate (B=10,000 by default).

Design notes
------------
* GSEA statistic = normalized enrichment score (NES) and FDR from gseapy.prerank.
* The ranking metric is the per-gene log2FC (Analysis A) or pooled log2FC
  (Analysis B), so enrichment reflects coordinated directional shift, not a DEG
  cutoff. This is the pre-registered PRIMARY test.
* The permutation null answers the specific reviewer question "is EV enrichment
  greater than a covariate-matched background?" using the SAME enrichment
  statistic (mean signed rank-score of the set) for EV vs each sampled null set.
* All randomness is seeded for reproducibility.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Iterable


# ----------------------------------------------------------------------------
# GSEA (primary)
# ----------------------------------------------------------------------------
def run_gsea_prerank(ranked: pd.Series, gene_sets: dict, seed: int = 0,
                     min_size: int = 5, max_size: int = 5000,
                     permutation_num: int = 1000) -> pd.DataFrame:
    """GSEA preranked over one or more gene sets.

    Parameters
    ----------
    ranked : pd.Series indexed by gene symbol, values = ranking metric (log2FC).
             Should contain the full genome-wide axis (no DEG cutoff).
    gene_sets : {set_name: iterable_of_genes}
    Returns a tidy DataFrame: term, es, nes, pval, fdr, size, leading_edge.
    """
    import gseapy as gp

    rnk = ranked.dropna().sort_values(ascending=False)
    # gseapy expects unique gene index
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    gs = {k: sorted(set(v) & set(rnk.index)) for k, v in gene_sets.items()}
    gs = {k: v for k, v in gs.items() if len(v) >= min_size}
    if not gs:
        return pd.DataFrame(columns=["term", "es", "nes", "pval", "fdr", "size", "leading_edge"])

    res = gp.prerank(rnk=rnk.reset_index().rename(columns={"index": 0, rnk.name or 0: 1}) if False else rnk,
                     gene_sets=gs, min_size=min_size, max_size=max_size,
                     permutation_num=permutation_num, seed=seed, no_plot=True,
                     outdir=None, threads=4, verbose=False)
    r = res.res2d.copy()
    # normalize column names across gseapy versions
    colmap = {c: c.lower().replace(" ", "_") for c in r.columns}
    r = r.rename(columns=colmap)
    def pick(*names):
        for n in names:
            if n in r.columns:
                return n
        return None
    out = pd.DataFrame({
        "term": r[pick("term", "name")],
        "es": pd.to_numeric(r[pick("es")], errors="coerce"),
        "nes": pd.to_numeric(r[pick("nes")], errors="coerce"),
        "pval": pd.to_numeric(r[pick("nom_p-val", "pval", "p-val")], errors="coerce"),
        "fdr": pd.to_numeric(r[pick("fdr_q-val", "fdr", "fdr_qval")], errors="coerce"),
    })
    sz = pick("gene_%", "tag_%", "size")
    le = pick("lead_genes", "leading_edge", "ledge_genes")
    out["size"] = [len(gs.get(t, [])) for t in out["term"]]
    out["leading_edge"] = r[le].values if le else ""
    return out.reset_index(drop=True)


# ----------------------------------------------------------------------------
# ORA (secondary / supportive)
# ----------------------------------------------------------------------------
def run_ora(deg_genes: Iterable[str], gene_set: Iterable[str],
            universe: Iterable[str]) -> dict:
    """Hypergeometric over-representation of a gene set in a DEG list.

    Returns dict with k (overlap), K (set in universe), n (DEGs in universe),
    N (universe size), fold_enrichment, and p (survival function).
    """
    from scipy.stats import hypergeom

    universe = set(universe)
    deg = set(deg_genes) & universe
    gs = set(gene_set) & universe
    N = len(universe)
    K = len(gs)
    n = len(deg)
    k = len(deg & gs)
    if N == 0 or K == 0 or n == 0:
        return dict(k=k, K=K, n=n, N=N, fold_enrichment=np.nan, p=np.nan)
    expected = n * K / N
    fold = (k / expected) if expected > 0 else np.nan
    # P(X >= k)
    p = hypergeom.sf(k - 1, N, K, n)
    return dict(k=k, K=K, n=n, N=N, fold_enrichment=fold, p=float(p))


# ----------------------------------------------------------------------------
# Covariate-matched permutation
# ----------------------------------------------------------------------------
def _set_score(ranked: pd.Series, genes: set) -> float:
    """Enrichment statistic for the permutation null: mean of the signed
    standardized rank-score over set members. Uses rankdata so it is robust and
    directly comparable between the EV set and matched random sets."""
    from scipy.stats import rankdata
    g = [x for x in genes if x in ranked.index]
    if not g:
        return np.nan
    # standardized rank in [-1, 1]: high positive logFC -> +1
    r = pd.Series(rankdata(ranked.values), index=ranked.index)
    r = (r - r.mean()) / (r.max() - r.min() + 1e-12)  # ~[-0.5,0.5]
    return float(r.loc[g].mean())


def build_covariate_table(expr: pd.DataFrame, gene_lengths: pd.Series | None = None) -> pd.DataFrame:
    """Build per-gene covariates for matching from an expression matrix
    (genes x samples, linear or log scale):
      - expr_level : mean expression across samples
      - detection  : fraction of samples with expression > 0 (or > min)
      - length     : gene length if provided (else NaN, dropped from matching)
    """
    lvl = expr.mean(axis=1)
    det = (expr > expr.min().min()).mean(axis=1)
    cov = pd.DataFrame({"expr_level": lvl, "detection": det})
    if gene_lengths is not None:
        cov["length"] = gene_lengths.reindex(cov.index)
    return cov


def matched_permutation_pvalue(ranked: pd.Series, ev_genes: Iterable[str],
                               covariates: pd.DataFrame, B: int = 10000,
                               n_bins: int = 10, seed: int = 0,
                               exclude: Iterable[str] | None = None) -> dict:
    """Empirical p that EV-set enrichment exceeds covariate-matched random sets.

    For each of B iterations, sample a null gene set of the SAME size as the EV
    set, matched to the EV set's joint covariate distribution via stratified
    binning (quantile bins per covariate), drawn from NON-EV genes. Compare the
    observed EV set score to the null distribution of scores.

    Returns dict: observed, null_mean, null_sd, z, p_greater, B_effective.
    """
    rng = np.random.default_rng(seed)
    ranked = ranked.dropna()
    genes_axis = set(ranked.index) & set(covariates.index)
    cov = covariates.loc[sorted(genes_axis)].copy()
    cov = cov.dropna(axis=1, how="all")  # drop length if absent

    ev = sorted(set(ev_genes) & genes_axis)
    exclude = set(exclude or set()) | set(ev_genes)
    pool = [g for g in cov.index if g not in exclude]
    if len(ev) < 3 or len(pool) < len(ev):
        return dict(observed=np.nan, null_mean=np.nan, null_sd=np.nan, z=np.nan,
                    p_greater=np.nan, B_effective=0)

    # Assign each gene to a joint covariate bin
    bincols = []
    for c in cov.columns:
        try:
            b = pd.qcut(cov[c].rank(method="first"), q=n_bins, labels=False, duplicates="drop")
        except ValueError:
            b = pd.cut(cov[c], bins=n_bins, labels=False)
        bincols.append(b.astype("Int64").astype(str))
    joint = bincols[0]
    for b in bincols[1:]:
        joint = joint.str.cat(b, sep="_")
    joint = pd.Series(joint.values, index=cov.index)

    # Pool genes grouped by joint bin (non-EV only)
    pool_by_bin: dict[str, list] = {}
    for g in pool:
        pool_by_bin.setdefault(joint.loc[g], []).append(g)

    ev_bins = joint.loc[ev].values
    from collections import Counter
    bin_need = Counter(ev_bins)

    obs = _set_score(ranked, set(ev))
    null = np.empty(B)
    feasible = True
    for i in range(B):
        picked = []
        ok = True
        for b, cnt in bin_need.items():
            cand = pool_by_bin.get(b, [])
            if len(cand) >= cnt:
                idx = rng.choice(len(cand), size=cnt, replace=False)
                picked.extend([cand[j] for j in idx])
            else:
                # fall back: sample with replacement from whole pool for this deficit
                idx = rng.choice(len(pool), size=cnt, replace=True)
                picked.extend([pool[j] for j in idx])
                feasible = False
        null[i] = _set_score(ranked, set(picked))

    null_mean = float(np.nanmean(null))
    null_sd = float(np.nanstd(null, ddof=1))
    z = (obs - null_mean) / (null_sd + 1e-12)
    # one-sided empirical p (EV enrichment greater, i.e. more up-shifted)
    p_greater = float((np.sum(null >= obs) + 1) / (B + 1))
    return dict(observed=obs, null_mean=null_mean, null_sd=null_sd, z=z,
                p_greater=p_greater, B_effective=B, strict_matching=feasible)
