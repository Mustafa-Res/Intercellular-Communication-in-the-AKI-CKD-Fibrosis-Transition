"""Pillar 2 PILOT — pull ONE AKI dataset from Census, map compartments, build
donor x compartment pseudobulk. Measures memory/runtime + validates logic before
the full pull. Reports peak RSS.
"""
import time, resource, gc
import numpy as np
import pandas as pd
import scanpy as sc
import cellxgene_census

t0 = time.time()
CENSUS_VERSION = "2025-01-30"
PILOT_DATASET = "0b75c598-0893-4216-afe8-5414cab7739d"  # smallest AKI-containing dataset

def peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # KB->GB on Linux

vf = (
    "tissue_general == 'kidney' and is_primary_data == True "
    f"and dataset_id == '{PILOT_DATASET}'"
)

print(f"[{time.time()-t0:.0f}s] opening Census, pulling dataset {PILOT_DATASET} ...")
with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="Homo sapiens",
        obs_value_filter=vf,
        obs_column_names=[
            "disease", "donor_id", "dataset_id", "assay", "cell_type",
            "tissue", "suspension_type",
        ],
        var_column_names=["feature_id", "feature_name"],
    )
print(f"[{time.time()-t0:.0f}s] pulled AnnData: {adata.shape}  peak={peak_gb():.1f}GB")
print("  X dtype:", adata.X.dtype, "| sparse:", hasattr(adata.X, 'tocsr'))
print("  diseases:", adata.obs['disease'].value_counts().to_dict())
print("  suspension:", adata.obs['suspension_type'].value_counts().to_dict())
print("  n cell types:", adata.obs['cell_type'].nunique())

# --- compartment mapping (frozen ontology -> compartment) ---
SENDER_PT = {"epithelial cell of proximal tubule"}
SENDER_TUBULAR_OTHER = {
    "kidney loop of Henle thick ascending limb epithelial cell",
    "kidney distal convoluted tubule epithelial cell",
    "kidney connecting tubule epithelial cell",
    "kidney collecting duct principal cell",
    "kidney collecting duct intercalated cell",
    "renal intercalated cell",
    "renal principal cell",
    "kidney loop of Henle thin descending limb epithelial cell",
    "kidney loop of Henle thin ascending limb epithelial cell",
    "kidney loop of Henle epithelial cell",
}
RECEIVER_MYELOID = {
    "macrophage", "kidney interstitial alternatively activated macrophage",
    "monocyte", "mononuclear phagocyte", "conventional dendritic cell",
    "myeloid cell", "plasmacytoid dendritic cell", "dendritic cell",
}
RECEIVER_STROMAL = {
    "kidney interstitial fibroblast", "myofibroblast cell",
    "renal interstitial pericyte", "fibroblast", "mesenchymal cell",
}
RECEIVER_ENDO = {
    "endothelial cell", "endothelial cell of vascular tree",
    "capillary endothelial cell", "kidney capillary endothelial cell",
}
INJURY_CONTEXT = {"podocyte", "parietal epithelial cell"}

def map_compartment(ct):
    if ct in SENDER_PT: return "PT"
    if ct in SENDER_TUBULAR_OTHER: return "tubular_other"
    if ct in RECEIVER_MYELOID: return "receiver_myeloid"
    if ct in RECEIVER_STROMAL: return "receiver_stromal"
    if ct in RECEIVER_ENDO: return "receiver_endothelial"
    if ct in INJURY_CONTEXT: return "injury_context"
    return "other"

adata.obs["compartment"] = adata.obs["cell_type"].astype(str).map(map_compartment)
print("\n  compartment cell counts:")
print(adata.obs["compartment"].value_counts().to_string())

# --- pseudobulk: sum raw counts per (donor x compartment) ---
print(f"\n[{time.time()-t0:.0f}s] building donor x compartment pseudobulk ...")
adata.obs["dc"] = adata.obs["donor_id"].astype(str) + "||" + adata.obs["compartment"].astype(str)
from scipy.sparse import csr_matrix
X = adata.X.tocsr() if hasattr(adata.X, "tocsr") else csr_matrix(adata.X)
groups = adata.obs["dc"].values
uniq = pd.unique(groups)
gidx = {g: i for i, g in enumerate(uniq)}
rows = np.array([gidx[g] for g in groups])
# indicator (n_groups x n_cells) @ X  -> (n_groups x n_genes) summed counts
from scipy.sparse import coo_matrix
ind = coo_matrix((np.ones(len(rows)), (rows, np.arange(len(rows)))), shape=(len(uniq), X.shape[0])).tocsr()
pb = ind @ X  # summed counts
pb = np.asarray(pb.todense())
n_cells_per = np.asarray(ind.sum(axis=1)).ravel()

meta = pd.DataFrame({"dc": uniq})
meta[["donor_id", "compartment"]] = meta["dc"].str.split(r"\|\|", expand=True)
meta["n_cells"] = n_cells_per
# attach donor-level covariates
donor_cov = adata.obs.drop_duplicates("donor_id").set_index("donor_id")[
    ["disease", "dataset_id", "assay", "suspension_type"]
]
meta = meta.join(donor_cov, on="donor_id")

print(f"[{time.time()-t0:.0f}s] pseudobulk matrix: {pb.shape}  peak={peak_gb():.1f}GB")
print("  pseudobulk samples per compartment (>=25 cells):")
print(meta[meta.n_cells >= 25].groupby("compartment").size().to_string())
print("\n  PILOT pseudobulk sample table (PT only):")
print(meta[meta.compartment == "PT"][["donor_id", "n_cells", "disease", "suspension_type"]].to_string(index=False))

print(f"\n[{time.time()-t0:.0f}s] PILOT DONE. peak RSS = {peak_gb():.1f} GB")
print(f"  per-1000-cells memory estimate: {peak_gb()/adata.shape[0]*1000:.3f} GB/1k cells")
print(f"  extrapolated to 654k cells (naive linear): {peak_gb()/adata.shape[0]*653845:.1f} GB")
