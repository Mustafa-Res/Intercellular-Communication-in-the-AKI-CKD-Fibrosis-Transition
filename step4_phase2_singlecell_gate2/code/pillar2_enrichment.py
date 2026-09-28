"""Pillar 2 Step 4 — Sender EV-specificity test using the FROZEN Pillar-1 engine.

For each sender compartment (PT = primary, tubular_other = secondary) and each
contrast, on the limma-voom (primary engine) logFC ranking:
  - GSEA prerank over the full gene-set ladder (EV_TIER1_PRIMARY governs Gate 2A/2B).
  - ORA (hypergeometric) on FDR<0.05 DEGs (supporting).
  - Covariate-matched permutation (B=10,000): EV_TIER1_PRIMARY vs matched non-EV
    background, plus the negative-control ladder (ribosomal must be NULL).
  - Effect-size floor: |NES|>1 AND median EV_TIER1 logFC>0 with |median|>=0.20.

Gene axis: DE tables use feature_id; gene sets are SYMBOL-based -> map via gene_map.
Covariates (expr_level, detection) built from the compartment pseudobulk counts.
Gene-length matching omitted by design (documented Pillar-1 caveat; mygene/GTF absent).
"""
import sys, os, json
import numpy as np
import pandas as pd

sys.path.insert(0, "/mnt/results/ev_aki_ckd_project/01_bulk_signature/code")
from enrichment import (run_gsea_prerank, run_ora, matched_permutation_pvalue,
                        build_covariate_table)

DURABLE = "/mnt/shared-workspace/shared/pillar2"
DE_DIR = f"{DURABLE}/de_out"
IN = f"{DURABLE}/de_inputs"
OUT = "/mnt/results/ev_aki_ckd_project/02_singlecell_specificity/data"
os.makedirs(OUT, exist_ok=True)

GENESETS = json.load(open("/mnt/results/ev_aki_ckd_project/00_preregistration/feasibility/ev_genesets_frozen_v1.3.json"))
GENESETS = {k: v for k, v in GENESETS.items() if not k.startswith("_")}
LADDER = ["EV_TIER1_PRIMARY", "EV_TIER2_CARGO", "NONEV_SECRETED_BACKGROUND",
          "CTRL_LYSO_ENDOSOMAL", "CTRL_MITOCHONDRIAL", "CTRL_RIBOSOMAL"]
gene_map = pd.read_csv(f"{DURABLE}/gene_map.tsv", sep="\t").set_index("feature_id")["feature_name"]

SENDERS = ["PT", "tubular_other"]
CONTRASTS = ["AKI_vs_normal", "CKD_vs_normal", "AKI_vs_CKD"]
B_PERM = 10000
SEED = 0

def to_symbol_ranked(de):
    """feature_id-indexed DE table -> symbol-indexed logFC series (primary=voom logFC)."""
    s = de.copy()
    s["symbol"] = s["feature_id"].map(gene_map)
    s = s.dropna(subset=["symbol"])
    # collapse dup symbols by max |logFC|
    s["absfc"] = s["logFC"].abs()
    s = s.sort_values("absfc", ascending=False).drop_duplicates("symbol")
    return pd.Series(s["logFC"].values, index=s["symbol"].values, name="logFC")

def load_pseudobulk_symbols(comp):
    """compartment counts (genes x samples, feature_id) -> symbol-indexed CPM-ish matrix for covariates."""
    cnt = pd.read_csv(f"{IN}/counts_{comp}.tsv.gz", sep="\t", index_col=0)
    cnt.index = cnt.index.map(gene_map)
    cnt = cnt[~cnt.index.isna()]
    cnt = cnt.groupby(level=0).sum()
    # CPM
    cpm = cnt / (cnt.sum(axis=0) + 1e-9) * 1e6
    return cpm

gsea_rows, ora_rows, perm_rows, effect_rows = [], [], [], []

for comp in SENDERS:
    cpm = load_pseudobulk_symbols(comp)
    cov = build_covariate_table(cpm)  # expr_level, detection
    for cn in CONTRASTS:
        fpath = f"{DE_DIR}/de_{comp}__{cn}.tsv"
        if not os.path.exists(fpath):
            print(f"  missing {fpath}, skip"); continue
        de = pd.read_csv(fpath, sep="\t")
        ranked = to_symbol_ranked(de)
        print(f"\n=== {comp} | {cn} : ranked {len(ranked)} symbols ===")

        # ---- GSEA over ladder ----
        gsea = run_gsea_prerank(ranked, {k: GENESETS[k] for k in LADDER},
                                seed=SEED, permutation_num=1000)
        gsea.insert(0, "compartment", comp); gsea.insert(1, "contrast", cn)
        gsea_rows.append(gsea)
        nes_t1 = gsea.loc[gsea["term"] == "EV_TIER1_PRIMARY", "nes"]
        fdr_t1 = gsea.loc[gsea["term"] == "EV_TIER1_PRIMARY", "fdr"]
        nes_t1 = float(nes_t1.iloc[0]) if len(nes_t1) else np.nan
        fdr_t1 = float(fdr_t1.iloc[0]) if len(fdr_t1) else np.nan
        print(f"  EV_TIER1 GSEA NES={nes_t1:.3f} FDR={fdr_t1:.4g}")

        # ---- ORA on FDR<0.05 DEGs (up-regulated, since injury induction) ----
        de2 = de.copy(); de2["symbol"] = de2["feature_id"].map(gene_map)
        degs_up = de2.loc[(de2["FDR"] < 0.05) & (de2["logFC"] > 0), "symbol"].dropna().tolist()
        universe = de2["symbol"].dropna().tolist()
        for term in ["EV_TIER1_PRIMARY", "CTRL_RIBOSOMAL", "NONEV_SECRETED_BACKGROUND"]:
            o = run_ora(degs_up, GENESETS[term], universe)
            o.update(compartment=comp, contrast=cn, term=term, n_degs_up=len(degs_up))
            ora_rows.append(o)

        # ---- matched permutation: EV_TIER1 + negative controls ----
        for term in ["EV_TIER1_PRIMARY", "CTRL_RIBOSOMAL", "CTRL_MITOCHONDRIAL"]:
            pr = matched_permutation_pvalue(ranked, GENESETS[term], cov,
                                            B=B_PERM, seed=SEED)
            pr.update(compartment=comp, contrast=cn, term=term)
            perm_rows.append(pr)
            print(f"  perm {term}: obs={pr['observed']:.4f} p_greater={pr['p_greater']:.4g} "
                  f"z={pr['z']:.2f} strict={pr.get('strict_matching')}")

        # ---- effect-size floor ----
        ev_syms = [g for g in GENESETS["EV_TIER1_PRIMARY"] if g in ranked.index]
        med_lfc = float(np.median(ranked.loc[ev_syms])) if ev_syms else np.nan
        floor_pass = (abs(nes_t1) > 1) and (med_lfc > 0) and (abs(med_lfc) >= 0.20)
        effect_rows.append(dict(compartment=comp, contrast=cn,
                                nes_ev_tier1=nes_t1, fdr_ev_tier1=fdr_t1,
                                median_ev_tier1_logFC=med_lfc, n_ev_genes=len(ev_syms),
                                effect_size_floor_pass=bool(floor_pass)))
        print(f"  effect-size floor: median EV logFC={med_lfc:.3f} |NES|>1={abs(nes_t1)>1} -> PASS={floor_pass}")

# ---- save ----
pd.concat(gsea_rows, ignore_index=True).to_csv(f"{OUT}/sender_ev_enrichment_gsea.tsv", sep="\t", index=False)
pd.DataFrame(ora_rows).to_csv(f"{OUT}/sender_ev_enrichment_ora.tsv", sep="\t", index=False)
pd.DataFrame(perm_rows).to_csv(f"{OUT}/permutation_null_sender.tsv", sep="\t", index=False)
pd.DataFrame(effect_rows).to_csv(f"{OUT}/sender_effect_size_floor.tsv", sep="\t", index=False)
print("\nSaved sender enrichment outputs to", OUT)
print(pd.DataFrame(effect_rows).to_string())
