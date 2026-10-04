"""Build gse210622__AKI / gse210622__normal h5ad shards for the A6 external replication.

v2: raw 10x matrices (6.79M barcodes each) -> donor assignment by full-barcode
majority match against metadata prefix groups (unique suffixes dominate; cross-donor
collisions ~10 cells are negligible), then column-subset to QC-passing cells only.
AKI arm = Non-COVID AKI donors only; normal arm = all 6 controls.
Shard format mirrors pillar3 Census shards (var_names=Ensembl, obs$compartment).
"""
import gzip, os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
import fast_matrix_market as fmm

RAW = "/workspace/gse210622_raw"
OUT = "/workspace/gse210622_shards"
os.makedirs(OUT, exist_ok=True)

CT_MAP = {"PT": "PT", "Leuk": "receiver_myeloid", "IntC": "receiver_stromal", "EC": "receiver_endothelial"}

# ---- metadata (mixed 7/8-field rows) ----
meta = {}
with gzip.open("/workspace/gse210622_meta.gz", "rt") as fh:
    fh.readline()
    for ln in fh:
        f = ln.split()
        if len(f) == 7:
            bc, orig, nc, nf, mt, ct, gr = f
        elif len(f) == 8:
            bc, orig, _x, nc, nf, mt, ct, gr = f
        else:
            continue
        meta[bc] = (orig, ct, gr)
print("metadata cells:", len(meta), flush=True)

# prefix -> {raw_barcode_suffix: full_barcode}
prefix_groups = {}
for bc in meta:
    p, suf = bc.split("_", 1)
    prefix_groups.setdefault(p, {})[suf] = bc
prefix_suffix_sets = {p: set(d) for p, d in prefix_groups.items()}
print("prefix groups:", {p: len(d) for p, d in prefix_groups.items()}, flush=True)

samples = sorted({fn[:-len(".matrix.mtx.gz")] for fn in os.listdir(RAW) if fn.endswith(".matrix.mtx.gz")})

keep = {"AKI": [], "normal": []}
gene_map_written = False
assignments = {}
for s in samples:
    with gzip.open(f"{RAW}/{s}.barcodes.tsv.gz", "rt") as fh:
        bcs = [l.strip() for l in fh]
    # raw files contain the full 10x whitelist -> restrict to barcodes with
    # non-zero total counts BEFORE matching (empty whitelist entries match everything)
    X = fmm.mmread(f"{RAW}/{s}.matrix.mtx.gz").tocsr()  # genes x barcodes
    tot = np.asarray(X.sum(axis=0)).ravel()
    # real cells have high totals; ambient whitelist droplets sit near 0-5.
    # >=500 cleanly separates the true donor (~10k cells) from cross-donor noise (<200).
    nonzero = set(np.flatnonzero(tot >= 500).tolist())
    # donor assignment: match metadata suffixes against NONZERO barcodes only
    scores = {}
    for p, sufs in prefix_suffix_sets.items():
        scores[p] = sum(1 for i in nonzero if bcs[i] in sufs)
    donor = max(scores, key=scores.get)
    assignments[s] = (donor, scores[donor])
    runner = sorted(scores.items(), key=lambda kv: -kv[1])[1]
    print(f"{s}: donor={donor} matched={scores[donor]} (runner-up {runner[0]}={runner[1]}); nonzero_barcodes={len(nonzero)}", flush=True)

    # cells to keep = barcodes in the assigned donor's metadata group
    keep_set = prefix_suffix_sets[donor]
    keep_idx = np.array([i for i in sorted(nonzero) if bcs[i] in keep_set], dtype=np.int64)
    if keep_idx.size == 0:
        print(f"  SKIP {s}: no cells"); continue

    X = X[:, keep_idx].T.tocsr()  # cells x genes
    feats = pd.read_csv(f"{RAW}/{s}.features.tsv.gz", sep="\t", header=None, names=["feature_id", "feature_name", "x"])
    if not gene_map_written:
        feats[["feature_id", "feature_name"]].to_csv(f"{OUT}/gene_map_gse210622.tsv", sep="\t", index=False)
        gene_map_written = True
    a = ad.AnnData(X=X)
    a.var_names = feats["feature_id"].values
    a.obs_names = [prefix_groups[donor][bcs[i]] for i in keep_idx]
    m = pd.DataFrame([meta[b] for b in a.obs_names], columns=["orig_ident", "celltype", "group"], index=a.obs_names)
    a.obs = m
    a.obs["donor"] = donor
    a.obs["sample"] = s
    a = a[a.obs["celltype"].isin(CT_MAP)].copy()
    a.obs["compartment"] = a.obs["celltype"].map(CT_MAP)
    groups = set(a.obs["group"])
    arm = "AKI" if "Non-COVID" in groups else ("normal" if "Control" in groups else None)
    if arm is None:
        print(f"  SKIP {s} ({donor}): COVID arm excluded"); continue
    a = a[a.obs["group"] == ("Non-COVID" if arm == "AKI" else "Control")].copy()
    keep[arm].append(a)
    print(f"  -> {arm}: {a.n_obs} cells in 4 compartments", flush=True)
    del X

pd.DataFrame([(s, d, n) for s, (d, n) in assignments.items()],
             columns=["sample", "donor_prefix", "matched_barcodes"]).to_csv(f"{OUT}/sample_donor_map.tsv", sep="\t", index=False)

for arm, cond in [("AKI", "AKI"), ("normal", "normal")]:
    A = ad.concat(keep[arm], axis=0, join="outer", fill_value=0)
    A.X = sparse.csr_matrix(A.X)
    out = f"{OUT}/gse210622__{cond}.h5ad"
    A.write_h5ad(out)
    print(f"WROTE {out}: {A.n_obs} cells x {A.n_vars} genes; donors={sorted(A.obs['donor'].unique())}")
    print("  compartments:", A.obs["compartment"].value_counts().to_dict(), flush=True)
print("SHARDS_DONE")
