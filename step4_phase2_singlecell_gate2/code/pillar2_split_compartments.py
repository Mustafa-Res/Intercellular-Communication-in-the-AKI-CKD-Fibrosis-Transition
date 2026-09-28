"""Pillar 2 Step 3a — split full pseudobulk into per-compartment count matrices
for edgeR. Writes genes x samples counts (gzipped TSV) + coldata per compartment.
Handles pandas categorical -> str so R design matrices are clean.

Only the 3 kidney diseases are retained (defensive; pull already restricted).
Only samples passing the >=25 cell floor are written.
"""
import os, json
import numpy as np
import pandas as pd

DURABLE = "/mnt/shared-workspace/shared/pillar2"
OUT = "/workspace/pillar2_de_inputs"
os.makedirs(OUT, exist_ok=True)

DISEASES = ["normal", "chronic kidney disease", "acute kidney failure"]
CELL_FLOOR = 25
# compartments to run DE on (senders + receivers + injury_context; skip 'other')
COMPARTMENTS = ["PT", "tubular_other", "receiver_myeloid", "receiver_stromal",
                "receiver_endothelial", "injury_context"]

print("loading pseudobulk parquet + meta ...")
counts = pd.read_parquet(f"{DURABLE}/pseudobulk_counts.parquet")  # sample_id col + genes
counts = counts.set_index("sample_id")
meta = pd.read_csv(f"{DURABLE}/pseudobulk_meta.tsv", sep="\t")
# the pull script wrote sample_id as the (unnamed) first col -> read as 'index'
if "sample_id" not in meta.columns:
    meta = meta.rename(columns={"index": "sample_id"})
# force str (kill categorical levels)
for c in ["disease", "dataset_id", "assay", "suspension_type", "donor_id", "compartment"]:
    meta[c] = meta[c].astype(str)
meta = meta.set_index("sample_id")

# align
meta = meta.loc[counts.index]
print(f"counts {counts.shape}  meta {meta.shape}")

# short disease codes for R friendliness
dcode = {"normal": "normal", "chronic kidney disease": "CKD", "acute kidney failure": "AKI"}
meta["disease_code"] = meta["disease"].map(dcode)

summary = {}
for comp in COMPARTMENTS:
    sel = meta[(meta["compartment"] == comp)
               & (meta["disease"].isin(DISEASES))
               & (meta["n_cells"] >= CELL_FLOOR)].copy()
    if sel.shape[0] == 0:
        print(f"  {comp}: 0 samples, skip"); continue
    cmat = counts.loc[sel.index]                     # samples x genes
    # drop all-zero genes within this compartment (edgeR filterByExpr will refilter,
    # but this shrinks files a lot)
    nz = (cmat.sum(axis=0) > 0)
    cmat = cmat.loc[:, nz]
    genes_x_samples = cmat.T                          # genes x samples
    genes_x_samples.index.name = "feature_id"
    # write
    fc = f"{OUT}/counts_{comp}.tsv.gz"
    fco = f"{OUT}/coldata_{comp}.tsv"
    genes_x_samples.to_csv(fc, sep="\t")
    cols = ["donor_id", "disease", "disease_code", "dataset_id", "suspension_type",
            "n_cells", "compartment"]
    sel[cols].reset_index().to_csv(fco, sep="\t", index=False)
    # per-disease donor counts
    dcts = sel.groupby("disease_code")["donor_id"].nunique().to_dict()
    summary[comp] = {"n_samples": int(sel.shape[0]), "n_genes_nonzero": int(cmat.shape[1]),
                     "donors": dcts}
    print(f"  {comp}: {sel.shape[0]} samples, {cmat.shape[1]} nonzero genes, donors={dcts}")

with open(f"{OUT}/split_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)

# copy durable
import shutil
os.makedirs(f"{DURABLE}/de_inputs", exist_ok=True)
for f in os.listdir(OUT):
    shutil.copy(f"{OUT}/{f}", f"{DURABLE}/de_inputs/{f}")
print("\nwrote per-compartment DE inputs -> durable/de_inputs/")
print(json.dumps(summary, indent=2))
