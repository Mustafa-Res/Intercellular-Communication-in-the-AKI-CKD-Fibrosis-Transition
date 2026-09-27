#!/usr/bin/env python
"""Covariate-matched permutation test (B=10,000) for the EV set + negative-control
ladder, in both contexts (native, transplant-meta). Writes a tidy summary TSV.

For each gene set, compares its signed rank-score (from the context's log2FC
ranking) to a null of 10,000 covariate-matched random non-EV sets. One-sided
empirical p for enrichment greater than matched background.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path("/mnt/results/ev_aki_ckd_project")
sys.path.insert(0, str(REPO / "01_bulk_signature" / "code"))
sys.path.insert(0, str(REPO / "00_preregistration" / "code"))
from ev_genesets import ev_genesets as evg
import enrichment as enr

OUT = REPO / "01_bulk_signature/data"
B = 10000
SEED = 0

GS = evg.build_ev_genesets()
NC = evg.build_negative_controls()
LADDER = {
    "EV_TIER1_PRIMARY": set(GS["EV_TIER1_PRIMARY"]),
    "EV_TIER2_CARGO": set(GS["EV_TIER2_CARGO"]),
    "NONEV_SECRETED_BACKGROUND": set(GS["NONEV_SECRETED_BACKGROUND"]),
    "CTRL_LYSO_ENDOSOMAL": set(NC["CTRL_LYSO_ENDOSOMAL"]),
    "CTRL_MITOCHONDRIAL": set(NC["CTRL_MITOCHONDRIAL"]),
    "CTRL_RIBOSOMAL": set(NC["CTRL_RIBOSOMAL"]),
}
# EV membership to exclude from every null draw (so background is truly non-EV)
EV_ALL = set(GS["EV_TIER1_PRIMARY"]) | set(GS["EV_TIER2_CARGO"]) | set(GS["EV_TIER1_MACHINERY"]) | set(GS["EV_MVB_PATHWAY"])

# Pre-registered expected direction of enrichment
EXPECT = {
    "EV_TIER1_PRIMARY": "up (strong)",
    "EV_TIER2_CARGO": "supportive (noisier)",
    "NONEV_SECRETED_BACKGROUND": "less specific",
    "CTRL_LYSO_ENDOSOMAL": "some overlap, weaker",
    "CTRL_MITOCHONDRIAL": "injury-assoc, not EV-specific",
    "CTRL_RIBOSOMAL": "null (not enriched)",
}


def context_covariates_from_de(de_paths):
    series = []
    for p in de_paths:
        d = pd.read_csv(p, sep="\t")
        s = d.set_index("gene")["AveExpr"]
        s = s[~s.index.duplicated()]
        series.append(s)
    idx = set(series[0].index)
    for s in series[1:]:
        idx &= set(s.index)
    idx = sorted(idx)
    avg = sum(s.loc[idx] for s in series) / len(series)
    return pd.DataFrame({"expr_level": avg})


def ranked_from(de_path=None, meta_path=None):
    if meta_path:
        m = pd.read_csv(meta_path, sep="\t")
        r = pd.Series(m.pooled_logFC.values, index=m.gene)
    else:
        d = pd.read_csv(de_path, sep="\t")
        r = pd.Series(d.logFC.values, index=d.gene)
    return r[~r.index.duplicated()].dropna()


CONTEXTS = {
    "native_GSE139061": dict(
        ranked=ranked_from(de_path=OUT / "per_cohort_de/GSE139061_de.tsv"),
        cov=context_covariates_from_de([OUT / "per_cohort_de/GSE139061_de.tsv"]),
    ),
    "transplant_meta": dict(
        ranked=ranked_from(meta_path=OUT / "analysisB_transplant_meta.tsv"),
        cov=context_covariates_from_de([
            OUT / "per_cohort_de/GSE30718_de.tsv",
            OUT / "per_cohort_de/GSE53769_de.tsv",
        ]),
    ),
}

rows = []
for cname, cfg in CONTEXTS.items():
    ranked, cov = cfg["ranked"], cfg["cov"]
    print(f"[{cname}] axis={len(ranked)} genes, covariates={cov.shape}", flush=True)
    for sname, sset in LADDER.items():
        res = enr.matched_permutation_pvalue(
            ranked=ranked, ev_genes=sset, covariates=cov,
            B=B, n_bins=10, seed=SEED, exclude=EV_ALL if sname.startswith("EV_") else set(sset),
        )
        res.update(dict(context=cname, term=sname, expected=EXPECT[sname]))
        rows.append(res)
        print(f"  {sname:28s} obs={res['observed']:+.4f} null={res['null_mean']:+.4f} "
              f"z={res['z']:+.2f} p={res['p_greater']:.4g} strict={res.get('strict_matching')}", flush=True)

df = pd.DataFrame(rows)[
    ["context", "term", "expected", "observed", "null_mean", "null_sd", "z", "p_greater", "B_effective", "strict_matching"]
]
df.to_csv(OUT / "permutation_null_summary.tsv", sep="\t", index=False)
print("\nWrote", OUT / "permutation_null_summary.tsv")
