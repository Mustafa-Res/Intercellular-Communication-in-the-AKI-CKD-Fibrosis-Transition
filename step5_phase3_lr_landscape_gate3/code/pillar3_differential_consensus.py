"""Pillar 3 differential + cross-engine/cross-dataset consensus (PLAN.md sec 4.2-4.3, sec 5).

Reads per-shard engine outputs:
  liana_<dataset>__<COND>.tsv        (LIANA, all cell-type pairs)
  cellchat_LR_<dataset>__<COND>.tsv  (CellChat, PT->receiver already filtered)
for dataset in {a12ccb9b, dea717d4}, COND in {AKI, normal, CKD}.

Produces, per engine per dataset, an AKI-vs-normal DIFFERENTIAL with SEPARATE columns:
  direction (up_in_AKI bool), significance flag, effect size(s).
Then the cross-engine + cross-dataset CONSENSUS with the FIXED, ABSOLUTE Gate-3 floor.

Axis identity = (receiver, ligand_symbol, receptor_symbol) with PT as fixed sender.
Ligand/receptor complexes normalized to a canonical sorted-underscore string so CellChat
and LIANA complex naming can be matched (sec A3).

Gate-3 fixed thresholds (a priori, data-INDEPENDENT):
  CellChat significant: permutation pval < 0.05 (AKI side)
  CellChat floor:      dprob >= 0.05 OR >= 20% relative increase
  LIANA significant / floor: PT->receiver ordinal rank improvement >= 50 positions
                             OR magnitude_rank (RRA) improvement >= 0.15  (AKI lower than normal)
Consensus (up-in-AKI landscape) requires:
  (1) up_in_AKI in BOTH engines AND significant in >=1 engine
  (2) same in BOTH datasets
  (3) clears the absolute effect floor on BOTH engines' scales
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd
pd.set_option("future.no_silent_downcasting", True)

DATASETS = ["a12ccb9b", "dea717d4"]
RECEIVERS = ["receiver_myeloid", "receiver_stromal", "receiver_endothelial"]

# ---- fixed thresholds (a priori) ----
CC_SIG_P = 0.05
CC_FLOOR_ABS = 0.05
CC_FLOOR_REL = 0.20
LIANA_FLOOR_RANKPOS = 50
LIANA_FLOOR_RRA = 0.15


def canon_complex(x):
    """Canonicalize a complex string: split on '_' or '&', sort, rejoin with '_'."""
    if pd.isna(x):
        return ""
    parts = str(x).replace("&", "_").split("_")
    parts = [p.strip() for p in parts if p.strip()]
    return "_".join(sorted(parts))


def axis_key(receiver, ligand, receptor):
    return f"{receiver}|{canon_complex(ligand)}|{canon_complex(receptor)}"


# ---------------- LIANA ----------------
def load_liana(indir, dataset, cond):
    f = f"{indir}/liana_{dataset}__{cond}.tsv"
    if not os.path.exists(f):
        return None
    lr = pd.read_csv(f, sep="\t")
    lr = lr[(lr["source"] == "PT") & (lr["target"].isin(RECEIVERS))].copy()
    if lr.empty:
        return lr
    lr["receiver"] = lr["target"]
    lr["ligand"] = lr["ligand_complex"]
    lr["receptor"] = lr["receptor_complex"]
    lr["key"] = [axis_key(r, l, rc) for r, l, rc in zip(lr["receiver"], lr["ligand"], lr["receptor"])]
    # integer ordinal position among PT->receiver axes (rank 1 = strongest = lowest magnitude_rank)
    lr = lr.sort_values("magnitude_rank").reset_index(drop=True)
    lr["mag_position"] = np.arange(1, len(lr) + 1)
    # dedupe key keeping best (lowest magnitude_rank)
    lr = lr.drop_duplicates("key", keep="first")
    return lr[["key", "receiver", "ligand", "receptor", "magnitude_rank",
               "specificity_rank", "mag_position", "lr_means", "cellphone_pvals"]]


def liana_differential(indir, dataset):
    aki = load_liana(indir, dataset, "AKI")
    nrm = load_liana(indir, dataset, "normal")
    if aki is None or nrm is None:
        return None
    m = aki.merge(nrm, on="key", how="outer", suffixes=("_AKI", "_normal"))
    # fill positions: axes absent in a condition get worst position + 1
    worst_aki = (aki["mag_position"].max() if len(aki) else 0) + 1
    worst_nrm = (nrm["mag_position"].max() if len(nrm) else 0) + 1
    m["mag_position_AKI"] = m["mag_position_AKI"].fillna(worst_aki)
    m["mag_position_normal"] = m["mag_position_normal"].fillna(worst_nrm)
    m["magnitude_rank_AKI"] = m["magnitude_rank_AKI"].fillna(1.0)
    m["magnitude_rank_normal"] = m["magnitude_rank_normal"].fillna(1.0)
    # coalesce annotation cols
    for c in ["receiver", "ligand", "receptor"]:
        m[c] = m[f"{c}_AKI"].fillna(m[f"{c}_normal"])
    # effect sizes
    m["rank_improvement"] = m["mag_position_normal"] - m["mag_position_AKI"]      # positive = up in AKI
    m["rra_improvement"] = m["magnitude_rank_normal"] - m["magnitude_rank_AKI"]   # positive = up in AKI
    # direction + significance/floor (LIANA floor doubles as its significance flag per plan)
    m["up_in_AKI_liana"] = (m["rank_improvement"] > 0) & (m["rra_improvement"] > 0)
    m["liana_floor"] = (m["rank_improvement"] >= LIANA_FLOOR_RANKPOS) | (m["rra_improvement"] >= LIANA_FLOOR_RRA)
    m["liana_sig"] = m["liana_floor"]  # LIANA statistical support == meeting the fixed rank/RRA thresholds
    m["dataset"] = dataset
    return m[["key", "dataset", "receiver", "ligand", "receptor",
              "mag_position_AKI", "mag_position_normal", "rank_improvement",
              "magnitude_rank_AKI", "magnitude_rank_normal", "rra_improvement",
              "up_in_AKI_liana", "liana_floor", "liana_sig"]]


# ---------------- CellChat ----------------
def load_cellchat(indir, dataset, cond):
    f = f"{indir}/cellchat_LR_{dataset}__{cond}.tsv"
    if not os.path.exists(f):
        return None
    df = pd.read_csv(f, sep="\t")
    if df.empty:
        return df
    # CellChat cols: source,target,ligand,receptor,prob,pval,interaction_name,pathway_name,...
    df = df[df["source"] == "PT"].copy()
    df = df[df["target"].isin(RECEIVERS)]
    df["receiver"] = df["target"]
    df["key"] = [axis_key(r, l, rc) for r, l, rc in zip(df["receiver"], df["ligand"], df["receptor"])]
    # dedupe key keeping strongest prob
    df = df.sort_values("prob", ascending=False).drop_duplicates("key", keep="first")
    return df[["key", "receiver", "ligand", "receptor", "prob", "pval", "pathway_name"]]


def cellchat_differential(indir, dataset):
    aki = load_cellchat(indir, dataset, "AKI")
    nrm = load_cellchat(indir, dataset, "normal")
    if aki is None or nrm is None:
        return None
    m = aki.merge(nrm, on="key", how="outer", suffixes=("_AKI", "_normal"))
    m["prob_AKI"] = m["prob_AKI"].fillna(0.0)
    m["prob_normal"] = m["prob_normal"].fillna(0.0)
    m["pval_AKI"] = m["pval_AKI"].fillna(1.0)
    for c in ["receiver", "ligand", "receptor"]:
        m[c] = m[f"{c}_AKI"].fillna(m[f"{c}_normal"])
    m["pathway_name"] = m["pathway_name_AKI"].fillna(m.get("pathway_name_normal"))
    m["dprob"] = m["prob_AKI"] - m["prob_normal"]
    with np.errstate(divide="ignore", invalid="ignore"):
        m["rel_increase"] = np.where(m["prob_normal"] > 0,
                                     (m["prob_AKI"] - m["prob_normal"]) / m["prob_normal"],
                                     np.inf)  # from 0 -> present = infinite relative increase
    m["up_in_AKI_cc"] = m["dprob"] > 0
    m["cc_sig"] = (m["pval_AKI"] < CC_SIG_P) & m["up_in_AKI_cc"]
    m["cc_floor"] = m["up_in_AKI_cc"] & ((m["dprob"] >= CC_FLOOR_ABS) | (m["rel_increase"] >= CC_FLOOR_REL))
    m["dataset"] = dataset
    return m[["key", "dataset", "receiver", "ligand", "receptor", "pathway_name",
              "prob_AKI", "prob_normal", "dprob", "rel_increase",
              "pval_AKI", "up_in_AKI_cc", "cc_sig", "cc_floor"]]


# ---------------- consensus ----------------
def build_consensus(indir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    liana_diffs, cc_diffs = [], []
    for ds in DATASETS:
        ld = liana_differential(indir, ds)
        cd = cellchat_differential(indir, ds)
        if ld is not None:
            ld.to_csv(f"{out_dir}/differential_liana_{ds}.tsv", sep="\t", index=False)
            liana_diffs.append(ld)
        if cd is not None:
            cd.to_csv(f"{out_dir}/differential_cellchat_{ds}.tsv", sep="\t", index=False)
            cc_diffs.append(cd)
    if not liana_diffs or not cc_diffs:
        print("Missing engine outputs; cannot build consensus yet.")
        return None
    L = pd.concat(liana_diffs, ignore_index=True)
    C = pd.concat(cc_diffs, ignore_index=True)

    # per (key, dataset) join engines
    per = C.merge(L, on=["key", "dataset", "receiver"], how="outer", suffixes=("", "_l"))
    # coalesce ligand/receptor
    for c in ["ligand", "receptor"]:
        per[c] = per[c].fillna(per[f"{c}_l"])
    # fill engine flags where an engine had no row for that axis in that dataset
    for col, val in [("up_in_AKI_cc", False), ("cc_sig", False), ("cc_floor", False),
                     ("up_in_AKI_liana", False), ("liana_sig", False), ("liana_floor", False)]:
        per[col] = per[col].fillna(val).astype(bool)
    # per-axis-per-dataset: direction both engines + sig>=1 + floor both
    per["dir_both"] = per["up_in_AKI_cc"] & per["up_in_AKI_liana"]
    per["sig_ge1"] = per["cc_sig"] | per["liana_sig"]
    per["floor_both"] = per["cc_floor"] & per["liana_floor"]
    per["passes_in_dataset"] = per["dir_both"] & per["sig_ge1"] & per["floor_both"]
    per.to_csv(f"{out_dir}/per_axis_per_dataset.tsv", sep="\t", index=False)

    # cross-dataset: require passes_in_dataset in BOTH datasets
    g = per.groupby("key")
    n_ds_pass = g["passes_in_dataset"].sum()
    n_ds_dirboth = g["dir_both"].sum()
    n_ds_dirsig = g.apply(lambda s: int(((s["dir_both"]) & (s["sig_ge1"])).sum()), include_groups=False) \
        if hasattr(g, "apply") else None
    keys_replicated = n_ds_pass[n_ds_pass >= 2].index  # both datasets

    # ---- consensus attrition funnel (for the Gate 3 memo; transparent drop-out) ----
    n_datasets = per["dataset"].nunique()
    funnel = {
        "n_datasets_present": int(n_datasets),
        "n_unique_axes_any_engine_any_dataset": int(per["key"].nunique()),
        "n_axes_dir_both_in_ge1_dataset": int((n_ds_dirboth >= 1).sum()),
        "n_axes_dir_both_in_both_datasets": int((n_ds_dirboth >= 2).sum()),
        "n_axes_passes(dir+sig+floor)_in_ge1_dataset": int((n_ds_pass >= 1).sum()),
        "n_axes_passes(dir+sig+floor)_in_both_datasets": int((n_ds_pass >= 2).sum()),
        "thresholds": {
            "CC_SIG_P": CC_SIG_P, "CC_FLOOR_ABS": CC_FLOOR_ABS, "CC_FLOOR_REL": CC_FLOOR_REL,
            "LIANA_FLOOR_RANKPOS": LIANA_FLOOR_RANKPOS, "LIANA_FLOOR_RRA": LIANA_FLOOR_RRA,
        },
    }
    with open(f"{out_dir}/consensus_funnel.json", "w") as fh:
        json.dump(funnel, fh, indent=2)
    print("CONSENSUS FUNNEL:", json.dumps(funnel, indent=2))
    # assemble consensus table
    rows = []
    for key in keys_replicated:
        sub = per[per["key"] == key]
        rows.append({
            "key": key,
            "receiver": sub["receiver"].iloc[0],
            "ligand": sub["ligand"].iloc[0],
            "receptor": sub["receptor"].iloc[0],
            "pathway_name": sub["pathway_name"].dropna().iloc[0] if sub["pathway_name"].notna().any() else "",
            "n_datasets_pass": int(sub["passes_in_dataset"].sum()),
            "cc_dprob_mean": float(sub["dprob"].mean()),
            "cc_sig_any": bool(sub["cc_sig"].any()),
            "liana_rank_improvement_mean": float(sub["rank_improvement"].mean()),
            "liana_rra_improvement_mean": float(sub["rra_improvement"].mean()),
            "liana_sig_any": bool(sub["liana_sig"].any()),
        })
    cons = pd.DataFrame(rows)
    if len(cons):
        cons = cons.sort_values(["receiver", "cc_dprob_mean"], ascending=[True, False]).reset_index(drop=True)
    cons.to_csv(f"{out_dir}/consensus_landscape_AKIvsnormal.tsv", sep="\t", index=False)
    print(f"Consensus axes (replicated, both datasets, both engines, floor): {len(cons)}")
    if len(cons):
        print(cons.head(20).to_string())
    return cons


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("indir")
    ap.add_argument("out_dir")
    args = ap.parse_args()
    build_consensus(args.indir, args.out_dir)
