"""Pillar 2 FULL pull — iterate over 6 kidney 10x 3' v3 datasets, pull each,
map ontology->compartment, build donor x compartment pseudobulk (raw count sums),
accumulate, and checkpoint to durable storage.

Census default X is RAW integer counts (verified). edgeR/limma-voom-ready.
Memory-safe: one dataset in memory at a time; only tiny pseudobulk matrices accumulate.
"""
import time, resource, gc, json
import numpy as np
import pandas as pd
import scanpy as sc
import cellxgene_census
from scipy.sparse import csr_matrix, coo_matrix

t0 = time.time()
CENSUS_VERSION = "2025-01-30"
DISEASES = ["normal", "chronic kidney disease", "acute kidney failure"]
ASSAY = "10x 3' v3"
OUT_DIR = "/workspace/pillar2_out"          # local (fast, random-access)
DURABLE = "/mnt/shared-workspace/shared/pillar2"  # checkpoint
import os; os.makedirs(OUT_DIR, exist_ok=True)

def peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6

# ---- compartment mapping (FROZEN) ----
SENDER_PT = {"epithelial cell of proximal tubule"}
SENDER_TUBULAR_OTHER = {
    "kidney loop of Henle thick ascending limb epithelial cell",
    "kidney distal convoluted tubule epithelial cell",
    "kidney connecting tubule epithelial cell",
    "kidney collecting duct principal cell",
    "kidney collecting duct intercalated cell",
    "renal intercalated cell", "renal principal cell",
    "kidney loop of Henle thin descending limb epithelial cell",
    "kidney loop of Henle thin ascending limb epithelial cell",
    "kidney loop of Henle epithelial cell",
    "renal alpha-intercalated cell", "renal beta-intercalated cell",
    "kidney cortex tubule cell", "epithelial cell of nephron",
}
RECEIVER_MYELOID = {
    "macrophage", "kidney interstitial alternatively activated macrophage",
    "monocyte", "mononuclear phagocyte", "conventional dendritic cell",
    "myeloid cell", "plasmacytoid dendritic cell", "dendritic cell",
    "kidney resident macrophage",
}
RECEIVER_STROMAL = {
    "kidney interstitial fibroblast", "myofibroblast cell",
    "renal interstitial pericyte", "fibroblast", "mesenchymal cell",
    "kidney interstitial cell",
}
RECEIVER_ENDO = {
    "endothelial cell", "endothelial cell of vascular tree",
    "capillary endothelial cell", "kidney capillary endothelial cell",
    "glomerular endothelial cell", "kidney glomerular endothelial cell",
}
INJURY_CONTEXT = {"podocyte", "parietal epithelial cell"}

CMAP = {}
for s, lab in [(SENDER_PT, "PT"), (SENDER_TUBULAR_OTHER, "tubular_other"),
               (RECEIVER_MYELOID, "receiver_myeloid"), (RECEIVER_STROMAL, "receiver_stromal"),
               (RECEIVER_ENDO, "receiver_endothelial"), (INJURY_CONTEXT, "injury_context")]:
    for ct in s:
        CMAP[ct] = lab

def map_compartment(ct):
    return CMAP.get(ct, "other")

# ---- discover datasets ----
dmap = pd.read_csv(f"{DURABLE}/v3_donor_dataset_map.csv")
datasets = sorted(dmap["dataset_id"].unique())
print(f"[{time.time()-t0:.0f}s] {len(datasets)} datasets to process")

pb_blocks = []      # list of (pseudobulk count DataFrame: samples x genes)
meta_blocks = []    # list of sample metadata
var_ref = None      # gene axis reference (feature_id -> feature_name)

for di, ds in enumerate(datasets):
    vf = (f"tissue_general == 'kidney' and is_primary_data == True "
          f"and disease in {DISEASES!r} and assay == {ASSAY!r} and dataset_id == '{ds}'")
    print(f"\n[{time.time()-t0:.0f}s] ({di+1}/{len(datasets)}) pulling {ds} ...")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        ad = cellxgene_census.get_anndata(
            census, organism="Homo sapiens", obs_value_filter=vf,
            obs_column_names=["disease","donor_id","dataset_id","assay","cell_type",
                              "tissue","suspension_type"],
            var_column_names=["feature_id","feature_name"],
        )
    if ad.n_obs == 0:
        print("   (empty; skipping)"); del ad; gc.collect(); continue
    ad.obs["compartment"] = ad.obs["cell_type"].astype(str).map(map_compartment)

    # gene axis: use feature_id (stable Census gene ID) as index
    ad.var_names = ad.var["feature_id"].astype(str).values
    if var_ref is None:
        var_ref = ad.var[["feature_id","feature_name"]].copy()

    # pseudobulk per (donor x compartment)
    ad.obs["dc"] = ad.obs["donor_id"].astype(str) + "||" + ad.obs["compartment"].astype(str)
    X = ad.X.tocsr() if hasattr(ad.X,"tocsr") else csr_matrix(ad.X)
    groups = ad.obs["dc"].values
    uniq = pd.unique(groups)
    gidx = {g:i for i,g in enumerate(uniq)}
    rows = np.fromiter((gidx[g] for g in groups), dtype=np.int64, count=len(groups))
    ind = coo_matrix((np.ones(len(rows)), (rows, np.arange(len(rows)))),
                     shape=(len(uniq), X.shape[0])).tocsr()
    pb = np.asarray((ind @ X).todense())               # summed raw counts
    ncell = np.asarray(ind.sum(axis=1)).ravel()

    pb_df = pd.DataFrame(pb, index=uniq, columns=ad.var_names)
    m = pd.DataFrame({"dc": uniq, "n_cells": ncell})
    m[["donor_id","compartment"]] = m["dc"].str.split(r"\|\|", expand=True)
    dcov = ad.obs.drop_duplicates("donor_id").set_index("donor_id")[
        ["disease","dataset_id","assay","suspension_type"]]
    m = m.join(dcov, on="donor_id")
    # unique sample id across datasets
    m["sample_id"] = ds[:8] + "|" + m["dc"]
    pb_df.index = m["sample_id"].values

    pb_blocks.append(pb_df)
    meta_blocks.append(m)
    print(f"   cells={ad.n_obs}  pseudobulk_samples={pb_df.shape[0]}  peak={peak_gb():.1f}GB")
    print(f"   compartments: {ad.obs['compartment'].value_counts().to_dict()}")
    del ad, X, pb, ind; gc.collect()

# ---- align gene axes (union, fill 0) and concat ----
print(f"\n[{time.time()-t0:.0f}s] aligning gene axes + concatenating ...")
all_genes = sorted(set().union(*[set(b.columns) for b in pb_blocks]))
pb_full = pd.concat([b.reindex(columns=all_genes, fill_value=0) for b in pb_blocks], axis=0)
meta_full = pd.concat(meta_blocks, axis=0, ignore_index=True)
meta_full = meta_full.set_index("sample_id").loc[pb_full.index].reset_index()
pb_full = pb_full.astype(np.int32)  # raw counts fit int32

print(f"   pseudobulk FULL: {pb_full.shape} (samples x genes)  peak={peak_gb():.1f}GB")
print(f"   samples per compartment (all): \n{meta_full.groupby('compartment').size().to_string()}")

# ---- donor floors + primary/exploratory flags (per compartment x disease) ----
CELL_FLOOR, DONOR_FLOOR = 25, 10
meta_full["pass_cell_floor"] = meta_full["n_cells"] >= CELL_FLOOR
kept = meta_full[meta_full["pass_cell_floor"]].copy()
donor_counts = (kept.groupby(["compartment","disease"])["donor_id"].nunique()
                .unstack(fill_value=0))
donor_counts.to_csv(f"{OUT_DIR}/donor_counts_by_compartment.tsv", sep="\t")
print(f"\n[{time.time()-t0:.0f}s] donor counts per (compartment x disease), >= {CELL_FLOOR} cells:")
print(donor_counts.to_string())
# primary if all 3 diseases (or at least AKI+normal) >= DONOR_FLOOR
def status_row(r):
    aki = r.get("acute kidney failure",0); ckd = r.get("chronic kidney disease",0); nrm = r.get("normal",0)
    if aki>=DONOR_FLOOR and nrm>=DONOR_FLOOR:
        return "primary"
    return "exploratory"
donor_counts["status"] = donor_counts.apply(status_row, axis=1)
print("\n  compartment status (AKI>=10 & normal>=10 -> primary):")
print(donor_counts["status"].to_string())

# ---- save checkpoints (write local first, then copy durable) ----
print(f"\n[{time.time()-t0:.0f}s] saving checkpoints ...")
# counts as parquet (samples x genes) + meta tsv + var map
pb_full.reset_index().rename(columns={"index":"sample_id"}).to_parquet(f"{OUT_DIR}/pseudobulk_counts.parquet")
meta_full.to_csv(f"{OUT_DIR}/pseudobulk_meta.tsv", sep="\t", index=False)
var_ref.to_csv(f"{OUT_DIR}/gene_map.tsv", sep="\t", index=False)
donor_counts.reset_index().to_csv(f"{OUT_DIR}/donor_counts_status.tsv", sep="\t", index=False)

import shutil
for f in ["pseudobulk_counts.parquet","pseudobulk_meta.tsv","gene_map.tsv",
          "donor_counts_by_compartment.tsv","donor_counts_status.tsv"]:
    shutil.copy(f"{OUT_DIR}/{f}", f"{DURABLE}/{f}")

summary = {
    "n_samples": int(pb_full.shape[0]), "n_genes": int(pb_full.shape[1]),
    "n_datasets": len(pb_blocks), "cell_floor": CELL_FLOOR, "donor_floor": DONOR_FLOOR,
    "peak_gb": round(peak_gb(),2), "runtime_s": round(time.time()-t0,1),
    "compartment_status": donor_counts["status"].to_dict(),
}
with open(f"{DURABLE}/pull_summary.json","w") as fh: json.dump(summary, fh, indent=2)
print(json.dumps(summary, indent=2))
print(f"\n[{time.time()-t0:.0f}s] FULL PULL DONE. peak RSS = {peak_gb():.1f} GB")
