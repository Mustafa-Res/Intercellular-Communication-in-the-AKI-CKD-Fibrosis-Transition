"""Pillar 2 Step 5 (robustness GSEA) + Step 7 (receiver descriptive GSEA).
Reuses FROZEN Pillar-1 enrichment engine (run_gsea_prerank). All rankings by SYMBOL.

Step 5 outputs:
  per_dataset_pt_ev_nes.tsv        - EV_TIER1 GSEA NES/FDR in each atlas (a12ccb9b, dea717d4)
  leave_one_dataset_out.tsv        - EV_TIER1 GSEA NES/FDR dropping each dataset (+ full)
  downsampled_control_sensitivity.tsv - EV_TIER1 GSEA over the 1000 down-sampled logFC medians
Step 7 outputs (NON-GATING, descriptive):
  receiver_descriptive_gsea.tsv    - EV ladder GSEA in receiver compartments, AKI_vs_normal
"""
import os, sys, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, "/mnt/results/ev_aki_ckd_project/01_bulk_signature/code")
from enrichment import run_gsea_prerank

DURABLE = "/mnt/shared-workspace/shared/pillar2"
ROB = f"{DURABLE}/robust_out"
OUT = "/mnt/results/ev_aki_ckd_project/02_singlecell_specificity/data"
GS = json.load(open("/mnt/results/ev_aki_ckd_project/00_preregistration/feasibility/ev_genesets_frozen_v1.3.json"))
GENE_SETS = {k: v for k, v in GS.items() if k != "_provenance"}
gene_map = pd.read_csv(f"{DURABLE}/gene_map.tsv", sep="\t").set_index("feature_id")["feature_name"]
SEED = 1; PERM = 1000
LADDER = ["EV_TIER1_PRIMARY","EV_TIER2_CARGO","NONEV_SECRETED_BACKGROUND",
          "CTRL_LYSO_ENDOSOMAL","CTRL_MITOCHONDRIAL","CTRL_RIBOSOMAL"]

def rank_from_logfc(df, lfc_col="logFC"):
    """collapse feature_id->symbol by max|logFC|, return Series indexed by symbol."""
    d = df.copy()
    if "symbol" not in d.columns:
        d["symbol"] = d["feature_id"].map(gene_map)
    d = d.dropna(subset=["symbol", lfc_col])
    d["abs"] = d[lfc_col].abs()
    d = d.sort_values("abs", ascending=False).drop_duplicates("symbol")
    return d.set_index("symbol")[lfc_col].sort_values(ascending=False)

def ev_nes(ranked):
    g = run_gsea_prerank(ranked, {k: GENE_SETS[k] for k in ["EV_TIER1_PRIMARY","CTRL_RIBOSOMAL"]},
                         seed=SEED, min_size=5, permutation_num=PERM)
    g = g.set_index("term")
    out = {}
    for t in ["EV_TIER1_PRIMARY","CTRL_RIBOSOMAL"]:
        out[f"{t}_nes"] = float(g.loc[t,"nes"]) if t in g.index else np.nan
        out[f"{t}_fdr"] = float(g.loc[t,"fdr"]) if t in g.index else np.nan
    return out

# ---------- Step 5a: per-dataset ----------
rows = []
for f in sorted(glob.glob(f"{ROB}/perdataset_PT_AKIvsnormal_*.tsv")):
    ds = os.path.basename(f).replace("perdataset_PT_AKIvsnormal_","").replace(".tsv","")
    r = rank_from_logfc(pd.read_csv(f, sep="\t"))
    d = ev_nes(r); d.update(dict(dataset=ds, n_genes=len(r)))
    rows.append(d)
pd.DataFrame(rows).to_csv(f"{OUT}/per_dataset_pt_ev_nes.tsv", sep="\t", index=False)
print("=== Step 5a per-dataset EV_TIER1 NES ===")
print(pd.DataFrame(rows)[["dataset","EV_TIER1_PRIMARY_nes","EV_TIER1_PRIMARY_fdr",
                          "CTRL_RIBOSOMAL_nes"]].to_string(index=False))

# ---------- Step 5b: leave-one-dataset-out ----------
rows = []
for f in sorted(glob.glob(f"{ROB}/loo_PT_AKIvsnormal_drop_*.tsv")):
    dr = os.path.basename(f).replace("loo_PT_AKIvsnormal_drop_","").replace(".tsv","")
    r = rank_from_logfc(pd.read_csv(f, sep="\t"))
    d = ev_nes(r); d.update(dict(dropped=dr, n_genes=len(r)))
    rows.append(d)
pd.DataFrame(rows).to_csv(f"{OUT}/leave_one_dataset_out.tsv", sep="\t", index=False)
print("\n=== Step 5b leave-one-dataset-out EV_TIER1 NES ===")
print(pd.DataFrame(rows)[["dropped","EV_TIER1_PRIMARY_nes","EV_TIER1_PRIMARY_fdr",
                          "CTRL_RIBOSOMAL_nes"]].to_string(index=False))

# ---------- Step 5c: down-sampled control sensitivity ----------
dsf = f"{ROB}/downsample_PT_AKIvsnormal_logFC.tsv.gz"
if os.path.exists(dsf):
    m = pd.read_csv(dsf, sep="\t")
    itc = [c for c in m.columns if c.startswith("iter")]
    # median logFC per gene across down-sampled iterations -> single ranking
    m["logFC"] = m[itc].median(axis=1)
    r = rank_from_logfc(m[["feature_id","symbol","logFC"]])
    med = ev_nes(r)
    # also per-iteration EV_TIER1 median log2FC distribution (effect-size stability)
    ev_syms = [g for g in GENE_SETS["EV_TIER1_PRIMARY"] if g in set(m["symbol"])]
    ev_rows = m[m["symbol"].isin(ev_syms)]
    per_iter_median = ev_rows[itc].median(axis=0)
    summ = dict(ev_tier1_nes_median_ranking=med["EV_TIER1_PRIMARY_nes"],
                ev_tier1_fdr_median_ranking=med["EV_TIER1_PRIMARY_fdr"],
                ctrl_ribosomal_nes_median_ranking=med["CTRL_RIBOSOMAL_nes"],
                ev_tier1_median_logFC_mean=float(per_iter_median.mean()),
                ev_tier1_median_logFC_q05=float(per_iter_median.quantile(0.05)),
                ev_tier1_median_logFC_q95=float(per_iter_median.quantile(0.95)),
                frac_iters_median_logFC_gt_0=float((per_iter_median>0).mean()),
                frac_iters_median_logFC_ge_0p20=float((per_iter_median>=0.20).mean()),
                n_iters=len(itc))
    pd.DataFrame([summ]).to_csv(f"{OUT}/downsampled_control_sensitivity.tsv", sep="\t", index=False)
    print("\n=== Step 5c down-sampled control sensitivity ===")
    for k,v in summ.items(): print(f"  {k}: {v}")
else:
    print("\n[Step 5c] down-sample matrix not found yet:", dsf)

# ---------- Step 7: receiver descriptive GSEA (NON-GATING) ----------
rows = []
for comp in ["receiver_endothelial","receiver_myeloid","receiver_stromal"]:
    f = f"{DURABLE}/de_out/de_{comp}__AKI_vs_normal.tsv"
    if not os.path.exists(f): continue
    r = rank_from_logfc(pd.read_csv(f, sep="\t"))
    g = run_gsea_prerank(r, GENE_SETS, seed=SEED, min_size=5, permutation_num=PERM).set_index("term")
    for t in LADDER:
        if t in g.index:
            rows.append(dict(compartment=comp, term=t, nes=float(g.loc[t,"nes"]),
                             fdr=float(g.loc[t,"fdr"]), size=int(g.loc[t,"size"])))
pd.DataFrame(rows).to_csv(f"{OUT}/receiver_descriptive_gsea.tsv", sep="\t", index=False)
print("\n=== Step 7 receiver descriptive GSEA (AKI_vs_normal, NON-GATING) ===")
rec = pd.DataFrame(rows)
print(rec[rec.term.isin(["EV_TIER1_PRIMARY","EV_TIER2_CARGO","CTRL_RIBOSOMAL"])]
      .to_string(index=False))
print("\nStep 5 GSEA + Step 7 complete ->", OUT)
