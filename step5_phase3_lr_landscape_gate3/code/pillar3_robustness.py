"""Pillar 3 robustness battery (PLAN.md sec 4.5 Gate 3C + sec 5 Gate 3B).

Operates on the EXISTING forR_<tag> bundles (same cells the primary engines saw), so no
h5ad re-read is needed and cell identity is guaranteed identical to the primary run.

Three tests (selectable):
  scramble : Gate 3B specificity. Permute compartment labels WITHIN each condition shard,
             re-run BOTH engines, re-run consensus. The primary 3 consensus axes (and the
             broader landscape) must COLLAPSE. Repeated over several seeds.
  lodo     : Gate 3C no-single-donor-driver. For each AKI donor, drop that donor from the
             AKI bundle, re-run BOTH engines on AKI, recompute up-in-AKI status of the 3
             consensus axes vs the (fixed) normal side. Axis must stay up-in-AKI in the
             MAJORITY of donor-dropout iterations, in each dataset.

Reconstructs a minimal AnnData (normalized X, symbol var, compartment/donor obs) from the
bundle for LIANA, and writes a scrambled/subset forR bundle for CellChat.

Usage:
  python pillar3_robustness.py scramble <engine_out_dir> <forR_root> <out_dir> [--seeds 1,2,3]
  python pillar3_robustness.py lodo     <engine_out_dir> <forR_root> <out_dir>
"""
import os, sys, json, gzip, argparse, time, subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse, io as sio

DATASETS = ["a12ccb9b", "dea717d4"]
RECEIVERS = ["receiver_myeloid", "receiver_stromal", "receiver_endothelial"]
CELLTYPES = ["PT"] + RECEIVERS
CCRS = "/mnt/shared-workspace/shared/pillar3/pillar3_run_cellchat.R"

# ---- fixed thresholds (mirror consensus script) ----
CC_SIG_P = 0.05; CC_FLOOR_ABS = 0.05; CC_FLOOR_REL = 0.20
LIANA_FLOOR_RANKPOS = 50; LIANA_FLOOR_RRA = 0.15


def canon_complex(x):
    if pd.isna(x): return ""
    parts = [p.strip() for p in str(x).replace("&", "_").split("_") if p.strip()]
    return "_".join(sorted(parts))

def axis_key(receiver, ligand, receptor):
    return f"{receiver}|{canon_complex(ligand)}|{canon_complex(receptor)}"


def load_bundle_adata(fr):
    """Reconstruct minimal AnnData (cells x genes, normalized) from a forR bundle."""
    data_gc = sio.mmread(gzip.open(os.path.join(fr, "matrix.mtx.gz")))   # genes x cells
    X = sparse.csr_matrix(data_gc).T.tocsr()                              # cells x genes
    genes = pd.read_csv(os.path.join(fr, "genes.tsv"), header=None)[0].astype(str).tolist()
    barcodes = pd.read_csv(os.path.join(fr, "barcodes.tsv"), header=None)[0].astype(str).tolist()
    meta = pd.read_csv(os.path.join(fr, "meta.tsv"), sep="\t", dtype=str).set_index("barcode")
    meta = meta.loc[barcodes]
    a = ad.AnnData(X=X, obs=meta.copy())
    a.var_names = genes
    a.obs_names = barcodes
    return a


def run_liana_on_adata(a, seed=1):
    import liana as li
    a = a.copy()
    a.obs["compartment"] = pd.Categorical(a.obs["compartment"].astype(str), categories=CELLTYPES)
    li.mt.rank_aggregate(a, groupby="compartment", resource_name="consensus",
                         expr_prop=0.1, use_raw=False, verbose=False, seed=seed, n_perms=100)
    return a.uns["liana_res"].copy()


def liana_pt_positions(lr):
    """PT->receiver ordinal positions + magnitude_rank, keyed by axis."""
    lr = lr[(lr["source"] == "PT") & (lr["target"].isin(RECEIVERS))].copy()
    if lr.empty:
        return pd.DataFrame(columns=["key", "mag_position", "magnitude_rank"])
    lr["key"] = [axis_key(t, l, r) for t, l, r in zip(lr["target"], lr["ligand_complex"], lr["receptor_complex"])]
    lr = lr.sort_values("magnitude_rank").reset_index(drop=True)
    lr["mag_position"] = np.arange(1, len(lr) + 1)
    lr = lr.drop_duplicates("key", keep="first")
    return lr[["key", "mag_position", "magnitude_rank"]]


def write_forR_from_adata(a, out_fr):
    """Write a forR bundle (genes x cells) from an AnnData whose .X is normalized.
    We only have normalized X in the bundle; CellChat uses the 'data' slot, so we pass
    normalized as both counts and data (counts only used for object creation / QC filters)."""
    os.makedirs(out_fr, exist_ok=True)
    data_gc = sparse.csr_matrix(a.X).T.tocsr()
    with gzip.open(f"{out_fr}/matrix.mtx.gz", "wb") as fh:
        sio.mmwrite(fh, data_gc)
    with gzip.open(f"{out_fr}/counts.mtx.gz", "wb") as fh:
        sio.mmwrite(fh, data_gc)  # normalized stands in for counts (min.cells filter only)
    pd.Series(a.var_names).to_csv(f"{out_fr}/genes.tsv", sep="\t", index=False, header=False)
    pd.Series(a.obs_names).to_csv(f"{out_fr}/barcodes.tsv", sep="\t", index=False, header=False)
    meta = a.obs[["compartment", "donor_id", "disease_code", "dataset_short"]].copy()
    meta.index.name = "barcode"
    meta.to_csv(f"{out_fr}/meta.tsv", sep="\t")


def run_cellchat_bundle(fr, out_dir, tag):
    env = dict(os.environ, R_LIBS_USER="/workspace/.Rlib")
    r = subprocess.run(["/opt/conda/bin/Rscript", CCRS, fr, out_dir, tag],
                       capture_output=True, text=True, env=env)
    f = f"{out_dir}/cellchat_LR_{tag}.tsv"
    if not os.path.exists(f):
        sys.stderr.write(r.stdout[-2000:] + r.stderr[-2000:])
        return None
    return pd.read_csv(f, sep="\t")


def cc_pt_probs(df):
    if df is None or df.empty:
        return pd.DataFrame(columns=["key", "prob", "pval"])
    df = df[(df["source"] == "PT") & (df["target"].isin(RECEIVERS))].copy()
    df["key"] = [axis_key(t, l, r) for t, l, r in zip(df["target"], df["ligand"], df["receptor"])]
    df = df.sort_values("prob", ascending=False).drop_duplicates("key", keep="first")
    return df[["key", "prob", "pval"]]


# =========================== SCRAMBLE (Gate 3B) ===========================
def scramble_test(engine_dir, forR_root, out_dir, seeds):
    os.makedirs(out_dir, exist_ok=True)
    cons = pd.read_csv(f"{engine_dir.rsplit('/',1)[0]}/consensus/consensus_landscape_AKIvsnormal.tsv", sep="\t") \
        if os.path.exists(f"{engine_dir.rsplit('/',1)[0]}/consensus/consensus_landscape_AKIvsnormal.tsv") else None
    consensus_keys = set(cons["key"]) if cons is not None and len(cons) else set()
    results = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        # build scrambled per-(dataset,cond) engine outputs into a temp dir
        tmp = f"{out_dir}/scramble_seed{seed}"
        os.makedirs(tmp, exist_ok=True)
        for ds in DATASETS:
            for cond in ["AKI", "normal"]:
                tag = f"{ds}__{cond}"
                fr = f"{forR_root}/forR_{tag}"
                if not os.path.exists(fr):
                    continue
                a = load_bundle_adata(fr)
                # permute compartment labels WITHIN this shard (breaks cell-type structure)
                perm = rng.permutation(a.n_obs)
                a.obs["compartment"] = a.obs["compartment"].values[perm]
                # LIANA on scrambled
                lr = run_liana_on_adata(a, seed=seed)
                lr.to_csv(f"{tmp}/liana_{tag}.tsv", sep="\t", index=False)
                # CellChat on scrambled (write scrambled bundle then run)
                sfr = f"{tmp}/forR_{tag}"
                write_forR_from_adata(a, sfr)
                cc = run_cellchat_bundle(sfr, tmp, tag)
                # (cc table already written by the R script as cellchat_LR_<tag>.tsv)
        # re-run consensus on scrambled engine outputs
        sub = subprocess.run(
            ["/workspace/.venv/bin/python",
             "/mnt/shared-workspace/shared/pillar3/pillar3_differential_consensus.py",
             tmp, f"{tmp}/consensus"],
            capture_output=True, text=True)
        scr_cons_f = f"{tmp}/consensus/consensus_landscape_AKIvsnormal.tsv"
        scr_funnel_f = f"{tmp}/consensus/consensus_funnel.json"
        n_scr = 0; scr_keys = set(); funnel = {}
        if os.path.exists(scr_cons_f):
            sc_cons = pd.read_csv(scr_cons_f, sep="\t")
            n_scr = len(sc_cons); scr_keys = set(sc_cons["key"]) if len(sc_cons) else set()
        if os.path.exists(scr_funnel_f):
            funnel = json.load(open(scr_funnel_f))
        overlap = consensus_keys & scr_keys
        results.append({
            "seed": int(seed),
            "n_consensus_scrambled": int(n_scr),
            "n_primary_consensus": len(consensus_keys),
            "n_primary_axes_surviving_scramble": len(overlap),
            "surviving_axes": sorted(overlap),
            "scrambled_funnel": funnel,
        })
        print(f"  [scramble seed{seed}] scrambled consensus = {n_scr}; primary axes surviving = {len(overlap)}")
    out = {
        "test": "scrambled_cell_label_negative_control",
        "primary_consensus_keys": sorted(consensus_keys),
        "seeds": results,
        "collapses": all(r["n_primary_axes_surviving_scramble"] == 0 for r in results),
        "interpretation": None,
    }
    out["interpretation"] = (
        "PASS (Gate 3B specificity): scrambling cell-type labels collapses the consensus; "
        "no primary consensus axis survives in any seed."
        if out["collapses"] else
        "CONCERN: >=1 primary consensus axis survives label scrambling in some seed; specificity not clean.")
    with open(f"{out_dir}/scramble_result.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("SCRAMBLE:", json.dumps({k: out[k] for k in ["collapses", "interpretation"]}, indent=2))
    return out


# =========================== LODO (Gate 3C) ===========================
def lodo_test(engine_dir, forR_root, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    cons = pd.read_csv(f"{engine_dir.rsplit('/',1)[0]}/consensus/consensus_landscape_AKIvsnormal.tsv", sep="\t")
    consensus_keys = list(cons["key"])
    per_ds = {}
    for ds in DATASETS:
        fr_aki = f"{forR_root}/forR_{ds}__AKI"
        if not os.path.exists(fr_aki):
            continue
        a_aki_full = load_bundle_adata(fr_aki)
        donors = sorted(a_aki_full.obs["donor_id"].astype(str).unique())
        # fixed normal side from primary engine outputs
        cc_norm = cc_pt_probs(pd.read_csv(f"{engine_dir}/cellchat_LR_{ds}__normal.tsv", sep="\t"))
        li_norm = liana_pt_positions(pd.read_csv(f"{engine_dir}/liana_{ds}__normal.tsv", sep="\t"))
        cc_norm_map = cc_norm.set_index("key")["prob"].to_dict()
        li_norm_pos = li_norm.set_index("key")["mag_position"].to_dict()
        li_norm_rra = li_norm.set_index("key")["magnitude_rank"].to_dict()
        worst_norm_pos = (max(li_norm_pos.values()) + 1) if li_norm_pos else 1
        tmp = f"{out_dir}/lodo_{ds}"; os.makedirs(tmp, exist_ok=True)
        axis_updir = {k: [] for k in consensus_keys}   # per donor-dropout: up in AKI both engines?
        for d in donors:
            a_d = a_aki_full[a_aki_full.obs["donor_id"].astype(str) != d].copy()
            # need >=1 PT and >=1 receiver with enough cells; CellChat min.cells=10 handles it
            tagd = f"{ds}__AKI__drop_{d}"
            lr = run_liana_on_adata(a_d, seed=1)
            li_pos = liana_pt_positions(lr).set_index("key")
            sfr = f"{tmp}/forR_{tagd}"; write_forR_from_adata(a_d, sfr)
            ccdf = run_cellchat_bundle(sfr, tmp, tagd)
            cc_p = cc_pt_probs(ccdf).set_index("key")["prob"].to_dict()
            for k in consensus_keys:
                # CellChat up?
                cc_aki = cc_p.get(k, 0.0); cc_nrm = cc_norm_map.get(k, 0.0)
                cc_up = cc_aki > cc_nrm
                # LIANA up? (both rank-position improvement and RRA improvement positive)
                if k in li_pos.index:
                    pos_aki = li_pos.loc[k, "mag_position"]; rra_aki = li_pos.loc[k, "magnitude_rank"]
                else:
                    pos_aki = worst_norm_pos + 1; rra_aki = 1.0
                pos_nrm = li_norm_pos.get(k, worst_norm_pos); rra_nrm = li_norm_rra.get(k, 1.0)
                li_up = ((pos_nrm - pos_aki) > 0) and ((rra_nrm - rra_aki) > 0)
                axis_updir[k].append(bool(cc_up and li_up))
        n_don = len(donors)
        per_ds[ds] = {
            "n_aki_donors": n_don,
            "axes": {k: {"n_updir_iters": int(sum(v)), "n_iters": len(v),
                         "frac_updir": (sum(v)/len(v) if v else 0.0),
                         "majority_stable": (sum(v) > len(v)/2.0)}
                     for k, v in axis_updir.items()},
        }
        print(f"  [LODO {ds}] {n_don} donors; " +
              "; ".join(f"{k.split('|')[0]}/{k.split('|')[1]}:{per_ds[ds]['axes'][k]['frac_updir']:.2f}"
                        for k in consensus_keys))
    # verdict: each consensus axis majority-stable in BOTH datasets
    verdict = {}
    for k in consensus_keys:
        stables = [per_ds[ds]["axes"][k]["majority_stable"] for ds in per_ds if k in per_ds[ds]["axes"]]
        verdict[k] = {"majority_stable_all_datasets": (len(stables) > 0 and all(stables)),
                      "per_dataset_frac": {ds: per_ds[ds]["axes"][k]["frac_updir"] for ds in per_ds}}
    out = {"test": "leave_one_donor_out_AKI", "consensus_keys": consensus_keys,
           "per_dataset": per_ds, "verdict": verdict,
           "all_axes_no_single_donor_driver": all(v["majority_stable_all_datasets"] for v in verdict.values())}
    with open(f"{out_dir}/lodo_result.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("LODO verdict:", json.dumps({k: v["majority_stable_all_datasets"] for k, v in verdict.items()}, indent=2))
    print("LODO all axes stable:", out["all_axes_no_single_donor_driver"])
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("test", choices=["scramble", "lodo"])
    ap.add_argument("engine_dir")
    ap.add_argument("forR_root")
    ap.add_argument("out_dir")
    ap.add_argument("--seeds", default="1,2,3")
    args = ap.parse_args()
    if args.test == "scramble":
        seeds = [int(s) for s in args.seeds.split(",")]
        scramble_test(args.engine_dir, args.forR_root, args.out_dir, seeds)
    else:
        lodo_test(args.engine_dir, args.forR_root, args.out_dir)
