"""Pillar 3 per-shard: normalize + subsample (cap 2000/compartment, seed=1),
run LIANA rank_aggregate (Python engine = CellPhoneDB + NATMI + Connectome +
SingleCellSignalR + logFC + RRA consensus), and export the IDENTICAL subsampled
cells to a Seurat-ingestible bundle for the CellChat engine.

Both engines therefore see the same cells, same genes (symbols), same normalization.

Usage:
    python pillar3_prep_and_liana.py <shard.h5ad> <out_dir> [--cap 2000] [--seed 1]

Writes into <out_dir>:
    liana_<shard>.tsv                      LIANA per-LR table (PT->receivers filtered downstream)
    forR_<shard>/matrix.mtx.gz             normalized data (genes x cells) for Seurat
    forR_<shard>/genes.tsv  barcodes.tsv   axes (gene = SYMBOL)
    forR_<shard>/counts.mtx.gz             raw counts (genes x cells) for Seurat
    forR_<shard>/meta.tsv                  per-cell metadata incl compartment
    prep_<shard>.json                      cell/gene counts, subsample record
"""
import sys, os, json, gzip, argparse, time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse, io as sio

RECEIVERS = ["receiver_myeloid", "receiver_stromal", "receiver_endothelial"]
CELLTYPES = ["PT"] + RECEIVERS

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("shard")
    ap.add_argument("out_dir")
    ap.add_argument("--cap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--gene_map", default="/mnt/shared-workspace/shared/pillar3/gene_map_pillar3.tsv")
    args = ap.parse_args()
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    tag = os.path.basename(args.shard).replace(".h5ad", "")
    rng = np.random.default_rng(args.seed)

    a = sc.read_h5ad(args.shard)
    # restrict to the 4 compartments (should already be), and to cells with a compartment label
    a = a[a.obs["compartment"].isin(CELLTYPES)].copy()

    # map gene axis (feature_id) -> SYMBOL for LR matching; collapse dup symbols by max total count
    gm = pd.read_csv(args.gene_map, sep="\t").set_index("feature_id")["feature_name"].to_dict()
    syms = pd.Index([gm.get(fid, fid) for fid in a.var_names])
    a.var["symbol"] = syms.values
    # drop genes with missing/empty symbol
    keep_g = a.var["symbol"].notna() & (a.var["symbol"].astype(str).str.len() > 0)
    a = a[:, keep_g.values].copy()
    # collapse duplicate symbols: keep the copy with highest total expression
    if a.var["symbol"].duplicated().any():
        tot = np.asarray(a.X.sum(axis=0)).ravel()
        order = np.argsort(-tot)
        seen, keep_idx = set(), []
        sym_arr = a.var["symbol"].values
        for i in order:
            s = sym_arr[i]
            if s not in seen:
                seen.add(s); keep_idx.append(i)
        keep_idx = np.sort(np.array(keep_idx))
        a = a[:, keep_idx].copy()
    a.var_names = a.var["symbol"].astype(str).values
    a.var_names_make_unique()

    # keep raw counts, then normalize (lib-size 1e4 + log1p) on a working layer
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)   # a.X now log-normalized

    # subsample cap per compartment (seed fixed)
    idx_keep = []
    rec = {}
    for ct in CELLTYPES:
        ci = np.where(a.obs["compartment"].values == ct)[0]
        rec[ct] = {"n_before": int(ci.size)}
        if ci.size > args.cap:
            ci = rng.choice(ci, size=args.cap, replace=False)
        rec[ct]["n_after"] = int(ci.size)
        idx_keep.append(ci)
    idx_keep = np.sort(np.concatenate(idx_keep))
    a = a[idx_keep].copy()
    a.obs["compartment"] = pd.Categorical(a.obs["compartment"].astype(str), categories=CELLTYPES)

    print(f"[{time.time()-t0:.0f}s] {tag}: {a.n_obs} cells x {a.n_vars} genes after prep")
    print("  subsample:", {k: (v["n_before"], v["n_after"]) for k, v in rec.items()})

    # ---------------- LIANA (Python engine) ----------------
    import liana as li
    # rank_aggregate uses the log-normalized .X ; label col = compartment
    li.mt.rank_aggregate(
        a, groupby="compartment", resource_name="consensus",
        expr_prop=0.1, use_raw=False, verbose=False, seed=args.seed,
        n_perms=100,
    )
    lr = a.uns["liana_res"].copy()
    lr.to_csv(f"{args.out_dir}/liana_{tag}.tsv", sep="\t", index=False)
    print(f"  LIANA: {lr.shape[0]} LR rows; cols={list(lr.columns)[:8]}...")

    # ---------------- export identical cells for CellChat (R) ----------------
    fr = f"{args.out_dir}/forR_{tag}"
    os.makedirs(fr, exist_ok=True)
    # Seurat wants genes x cells. write normalized (data) + counts.
    data_gc = a.X.T.tocsr() if sparse.issparse(a.X) else sparse.csr_matrix(a.X.T)
    counts_gc = a.layers["counts"].T.tocsr() if sparse.issparse(a.layers["counts"]) else sparse.csr_matrix(a.layers["counts"].T)
    with gzip.open(f"{fr}/matrix.mtx.gz", "wb") as fh:
        sio.mmwrite(fh, data_gc)
    with gzip.open(f"{fr}/counts.mtx.gz", "wb") as fh:
        sio.mmwrite(fh, counts_gc)
    pd.Series(a.var_names).to_csv(f"{fr}/genes.tsv", sep="\t", index=False, header=False)
    pd.Series(a.obs_names).to_csv(f"{fr}/barcodes.tsv", sep="\t", index=False, header=False)
    meta = a.obs[["compartment", "donor_id", "disease_code", "dataset_short"]].copy()
    meta.index.name = "barcode"
    meta.to_csv(f"{fr}/meta.tsv", sep="\t")

    with open(f"{args.out_dir}/prep_{tag}.json", "w") as fh:
        json.dump({"tag": tag, "n_cells": int(a.n_obs), "n_genes": int(a.n_vars),
                   "subsample": rec, "cap": args.cap, "seed": args.seed,
                   "liana_rows": int(lr.shape[0]),
                   "runtime_s": round(time.time()-t0, 1)}, fh, indent=2)
    print(f"[{time.time()-t0:.0f}s] {tag} DONE -> liana + forR bundle")

if __name__ == "__main__":
    main()
