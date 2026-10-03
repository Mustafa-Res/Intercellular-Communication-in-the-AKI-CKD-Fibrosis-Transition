"""Pillar 3 CELL-LEVEL pull — SAME frozen Census query as Pillar 2, but RETAIN the
cell x gene matrix (CellChat + LIANA both need cell-level expression + cell-type labels).

Scope (per approved PLAN.md sec 3): only the TWO replication atlases (a12ccb9b, dea717d4)
and only the FOUR Pillar-3 compartments (PT sender + 3 receivers). Written as one
compressed .h5ad shard per (dataset x disease). Raw counts kept in .X (Census default).

Memory-safe: one (dataset x disease) shard in memory at a time; subset to 4 compartments
BEFORE materializing to cut ~half the cells. Peak expected < 64 GB.

Frozen query: census_version=2025-01-30, Homo sapiens, tissue_general=kidney,
is_primary_data=True, assay=10x 3' v3, disease in {normal, CKD, AKF}.
"""
import time, resource, gc, json, os
import numpy as np
import pandas as pd
import scanpy as sc
import cellxgene_census

t0 = time.time()
CENSUS_VERSION = "2025-01-30"
ASSAY = "10x 3' v3"
DISEASES = ["normal", "chronic kidney disease", "acute kidney failure"]
DCODE = {"normal": "normal", "chronic kidney disease": "CKD", "acute kidney failure": "AKI"}
# Only the two atlases that support within-dataset AKI-vs-normal for all 4 cell types
ATLASES = {"a12ccb9b": None, "dea717d4": None}  # dataset_id[:8] -> full id (resolved below)

OUT_LOCAL = "/workspace/pillar3_out"       # random-access, fast
DURABLE = "/mnt/shared-workspace/shared/pillar3"
os.makedirs(OUT_LOCAL, exist_ok=True)
os.makedirs(DURABLE, exist_ok=True)

def peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6

# ---- compartment mapping (FROZEN — identical to Pillar 2) ----
SENDER_PT = {"epithelial cell of proximal tubule"}
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
CMAP = {}
for s, lab in [(SENDER_PT, "PT"), (RECEIVER_MYELOID, "receiver_myeloid"),
               (RECEIVER_STROMAL, "receiver_stromal"), (RECEIVER_ENDO, "receiver_endothelial")]:
    for ct in s:
        CMAP[ct] = lab
KEEP_COMPARTMENTS = {"PT", "receiver_myeloid", "receiver_stromal", "receiver_endothelial"}

# ---- resolve full dataset ids from the frozen donor-dataset map ----
dmap = pd.read_csv("/mnt/shared-workspace/shared/pillar2/v3_donor_dataset_map.csv")
for full_id in dmap["dataset_id"].unique():
    k = full_id[:8]
    if k in ATLASES:
        ATLASES[k] = full_id
missing = [k for k, v in ATLASES.items() if v is None]
assert not missing, f"could not resolve dataset ids for {missing}"
print(f"[{time.time()-t0:.0f}s] resolved atlases: {ATLASES}")

manifest = []
var_ref = None

for short_id, full_id in ATLASES.items():
    for disease in DISEASES:
        dcode = DCODE[disease]
        vf = (f"tissue_general == 'kidney' and is_primary_data == True "
              f"and disease == {disease!r} and assay == {ASSAY!r} "
              f"and dataset_id == '{full_id}'")
        print(f"\n[{time.time()-t0:.0f}s] pulling {short_id} / {disease} ({dcode}) ...")
        with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
            ad = cellxgene_census.get_anndata(
                census, organism="Homo sapiens", obs_value_filter=vf,
                obs_column_names=["disease", "donor_id", "dataset_id", "assay",
                                  "cell_type", "tissue", "suspension_type"],
                var_column_names=["feature_id", "feature_name"],
            )
        if ad.n_obs == 0:
            print("   (empty; skipping)")
            del ad; gc.collect(); continue

        # map compartment, subset to the 4 Pillar-3 compartments BEFORE anything else
        ad.obs["compartment"] = ad.obs["cell_type"].astype(str).map(lambda c: CMAP.get(c, "other"))
        keep = ad.obs["compartment"].isin(KEEP_COMPARTMENTS).values
        ad = ad[keep].copy()
        if ad.n_obs == 0:
            print("   (no cells in 4 compartments; skipping)")
            del ad; gc.collect(); continue

        # gene axis = feature_id (stable), keep symbol in var for LR matching
        ad.var_names = ad.var["feature_id"].astype(str).values
        ad.var["feature_name"] = ad.var["feature_name"].astype(str).values
        if var_ref is None:
            var_ref = ad.var[["feature_id", "feature_name"]].copy()

        # annotate shard identity
        ad.obs["disease_code"] = dcode
        ad.obs["dataset_short"] = short_id
        ad.obs["compartment"] = ad.obs["compartment"].astype(str)
        # X stays RAW counts (Census default). ensure int-compatible dtype preserved as float32 sparse.
        vc = ad.obs["compartment"].value_counts().to_dict()
        print(f"   cells={ad.n_obs}  genes={ad.n_vars}  peak={peak_gb():.1f}GB")
        print(f"   compartments: {vc}")
        print(f"   donors: {ad.obs['donor_id'].nunique()}")

        shard = f"{short_id}__{dcode}.h5ad"
        path_local = f"{OUT_LOCAL}/{shard}"
        ad.write_h5ad(path_local, compression="gzip")
        sz = os.path.getsize(path_local) / 1e6
        print(f"   wrote {shard} ({sz:.0f} MB)")
        manifest.append({
            "shard": shard, "dataset_short": short_id, "dataset_id": full_id,
            "disease": disease, "disease_code": dcode,
            "n_cells": int(ad.n_obs), "n_genes": int(ad.n_vars),
            "n_donors": int(ad.obs["donor_id"].nunique()),
            "compartments": json.dumps(vc), "size_mb": round(sz, 1),
        })
        del ad; gc.collect()

# ---- write gene map + manifest, copy shards to durable via shell cp (done outside) ----
var_ref.to_csv(f"{OUT_LOCAL}/gene_map_pillar3.tsv", sep="\t", index=False)
man = pd.DataFrame(manifest)
man.to_csv(f"{OUT_LOCAL}/pull_manifest.tsv", sep="\t", index=False)
print(f"\n[{time.time()-t0:.0f}s] MANIFEST:")
print(man.to_string())
summary = {
    "n_shards": len(manifest),
    "total_cells": int(man["n_cells"].sum()) if len(man) else 0,
    "peak_gb": round(peak_gb(), 2),
    "runtime_s": round(time.time() - t0, 1),
    "census_version": CENSUS_VERSION,
}
with open(f"{OUT_LOCAL}/pull_summary_pillar3.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print(json.dumps(summary, indent=2))
print(f"\n[{time.time()-t0:.0f}s] CELL-LEVEL PULL DONE. peak RSS = {peak_gb():.1f} GB")
print("NEXT: shell-cp shards + gene_map + manifest from OUT_LOCAL to DURABLE")
