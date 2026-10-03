"""Gate 3D biological-coherence sub-check: CKD-persistence of the 3 consensus axes.

The 3 consensus axes were discovered in the AKI-vs-normal contrast. Gate 3D asks
whether they are biologically coherent; one of its criteria is "CKD persistence" --
do the same axes remain UP in the CKD-vs-normal contrast (fibrotic chronic state),
consistent with a durable AKI->CKD fibrosis signal rather than an acute-only blip?

This reuses the EXACT validated differential logic from
pillar3_differential_consensus.py but with the case condition = CKD instead of AKI.
It reports, per consensus axis per dataset:
  - CellChat: prob_CKD, prob_normal, dprob (CKD-normal), up_in_CKD, cc_sig, cc_floor
  - LIANA: mag_position_CKD/normal, rank_improvement, rra_improvement, up_in_CKD, floor
and a per-axis 'persists' summary (up in CKD by BOTH engines in a dataset; replicated
across datasets = strong persistence).

CKD here is CKD-vs-normal *within each atlas* -- this is the secondary chronic contrast,
NOT pooling CKD with AKI. Reported as supportive coherence only (association language).

Usage: python pillar3_ckd_persistence.py <engine_out_dir> <consensus_tsv> <out_dir>
"""
import os, sys, json
import numpy as np
import pandas as pd

# Reuse the frozen, unit-tested primitives from the primary consensus script.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pillar3_differential_consensus as P

pd.set_option("future.no_silent_downcasting", True)

DATASETS = P.DATASETS
RECEIVERS = P.RECEIVERS
CC_SIG_P = P.CC_SIG_P
CC_FLOOR_ABS = P.CC_FLOOR_ABS
CC_FLOOR_REL = P.CC_FLOOR_REL
LIANA_FLOOR_RANKPOS = P.LIANA_FLOOR_RANKPOS
LIANA_FLOOR_RRA = P.LIANA_FLOOR_RRA


def liana_case_vs_normal(indir, dataset, case="CKD"):
    """LIANA differential, case vs normal, mirroring P.liana_differential but generic case."""
    case_df = P.load_liana(indir, dataset, case)
    nrm = P.load_liana(indir, dataset, "normal")
    if case_df is None or nrm is None:
        return None
    m = case_df.merge(nrm, on="key", how="outer", suffixes=("_case", "_normal"))
    worst_case = (case_df["mag_position"].max() if len(case_df) else 0) + 1
    worst_nrm = (nrm["mag_position"].max() if len(nrm) else 0) + 1
    m["mag_position_case"] = m["mag_position_case"].fillna(worst_case)
    m["mag_position_normal"] = m["mag_position_normal"].fillna(worst_nrm)
    m["magnitude_rank_case"] = m["magnitude_rank_case"].fillna(1.0)
    m["magnitude_rank_normal"] = m["magnitude_rank_normal"].fillna(1.0)
    for c in ["receiver", "ligand", "receptor"]:
        m[c] = m[f"{c}_case"].fillna(m[f"{c}_normal"])
    m["rank_improvement"] = m["mag_position_normal"] - m["mag_position_case"]
    m["rra_improvement"] = m["magnitude_rank_normal"] - m["magnitude_rank_case"]
    m["up_in_case_liana"] = (m["rank_improvement"] > 0) & (m["rra_improvement"] > 0)
    m["liana_floor"] = (m["rank_improvement"] >= LIANA_FLOOR_RANKPOS) | (m["rra_improvement"] >= LIANA_FLOOR_RRA)
    m["dataset"] = dataset
    return m[["key", "dataset", "receiver", "ligand", "receptor",
              "mag_position_case", "mag_position_normal", "rank_improvement",
              "magnitude_rank_case", "magnitude_rank_normal", "rra_improvement",
              "up_in_case_liana", "liana_floor"]]


def cellchat_case_vs_normal(indir, dataset, case="CKD"):
    case_df = P.load_cellchat(indir, dataset, case)
    nrm = P.load_cellchat(indir, dataset, "normal")
    if case_df is None or nrm is None:
        return None
    m = case_df.merge(nrm, on="key", how="outer", suffixes=("_case", "_normal"))
    m["prob_case"] = m["prob_case"].fillna(0.0)
    m["prob_normal"] = m["prob_normal"].fillna(0.0)
    m["pval_case"] = m["pval_case"].fillna(1.0)
    for c in ["receiver", "ligand", "receptor"]:
        m[c] = m[f"{c}_case"].fillna(m[f"{c}_normal"])
    m["pathway_name"] = m["pathway_name_case"].fillna(m.get("pathway_name_normal"))
    m["dprob"] = m["prob_case"] - m["prob_normal"]
    with np.errstate(divide="ignore", invalid="ignore"):
        m["rel_increase"] = np.where(m["prob_normal"] > 0,
                                     (m["prob_case"] - m["prob_normal"]) / m["prob_normal"],
                                     np.inf)
    m["up_in_case_cc"] = m["dprob"] > 0
    m["cc_sig"] = (m["pval_case"] < CC_SIG_P) & m["up_in_case_cc"]
    m["cc_floor"] = m["up_in_case_cc"] & ((m["dprob"] >= CC_FLOOR_ABS) | (m["rel_increase"] >= CC_FLOOR_REL))
    m["dataset"] = dataset
    return m[["key", "dataset", "receiver", "ligand", "receptor", "pathway_name",
              "prob_case", "prob_normal", "dprob", "rel_increase",
              "pval_case", "up_in_case_cc", "cc_sig", "cc_floor"]]


def main(indir, consensus_tsv, out_dir, case="CKD"):
    os.makedirs(out_dir, exist_ok=True)
    cons = pd.read_csv(consensus_tsv, sep="\t")
    consensus_keys = list(cons["key"])
    print(f"Consensus axes to test for {case}-persistence: {len(consensus_keys)}")
    for k in consensus_keys:
        print("  ", k)

    L = pd.concat([liana_case_vs_normal(indir, ds, case) for ds in DATASETS], ignore_index=True)
    C = pd.concat([cellchat_case_vs_normal(indir, ds, case) for ds in DATASETS], ignore_index=True)

    rows = []
    for key in consensus_keys:
        for ds in DATASETS:
            cc = C[(C["key"] == key) & (C["dataset"] == ds)]
            li = L[(L["key"] == key) & (L["dataset"] == ds)]
            cc_up = bool(cc["up_in_case_cc"].iloc[0]) if len(cc) else False
            cc_dprob = float(cc["dprob"].iloc[0]) if len(cc) else np.nan
            cc_prob_case = float(cc["prob_case"].iloc[0]) if len(cc) else 0.0
            cc_prob_norm = float(cc["prob_normal"].iloc[0]) if len(cc) else 0.0
            cc_pval = float(cc["pval_case"].iloc[0]) if len(cc) else np.nan
            cc_floor = bool(cc["cc_floor"].iloc[0]) if len(cc) else False
            li_up = bool(li["up_in_case_liana"].iloc[0]) if len(li) else False
            li_rank = float(li["rank_improvement"].iloc[0]) if len(li) else np.nan
            li_rra = float(li["rra_improvement"].iloc[0]) if len(li) else np.nan
            li_floor = bool(li["liana_floor"].iloc[0]) if len(li) else False
            rows.append({
                "key": key, "dataset": ds,
                f"cc_prob_{case}": cc_prob_case, "cc_prob_normal": cc_prob_norm,
                "cc_dprob": cc_dprob, "cc_pval": cc_pval,
                f"cc_up_in_{case}": cc_up, "cc_floor": cc_floor,
                f"liana_rank_improvement": li_rank, "liana_rra_improvement": li_rra,
                f"liana_up_in_{case}": li_up, "liana_floor": li_floor,
                "up_both_engines": cc_up and li_up,
            })
    detail = pd.DataFrame(rows)
    detail.to_csv(f"{out_dir}/ckd_persistence_detail.tsv", sep="\t", index=False)

    # per-axis summary: persists in a dataset = up in both engines; replicated = both datasets
    summ = []
    for key in consensus_keys:
        sub = detail[detail["key"] == key]
        n_ds_up_both = int(sub["up_both_engines"].sum())
        n_ds_cc_up = int(sub[f"cc_up_in_{case}"].sum())
        n_ds_li_up = int(sub[f"liana_up_in_{case}"].sum())
        summ.append({
            "key": key,
            "n_datasets_up_both_engines": n_ds_up_both,
            "n_datasets_cc_up": n_ds_cc_up,
            "n_datasets_liana_up": n_ds_li_up,
            "persists_both_datasets_both_engines": n_ds_up_both == len(DATASETS),
            "persists_ge1_dataset_both_engines": n_ds_up_both >= 1,
        })
    summary = pd.DataFrame(summ)
    summary.to_csv(f"{out_dir}/ckd_persistence_summary.tsv", sep="\t", index=False)

    result = {
        "case": case,
        "n_consensus_axes": len(consensus_keys),
        "n_persist_both_datasets_both_engines": int(summary["persists_both_datasets_both_engines"].sum()),
        "n_persist_ge1_dataset_both_engines": int(summary["persists_ge1_dataset_both_engines"].sum()),
        "per_axis": summary.to_dict(orient="records"),
    }
    with open(f"{out_dir}/ckd_persistence_result.json", "w") as fh:
        json.dump(result, fh, indent=2)
    print("\n=== CKD-PERSISTENCE DETAIL ===")
    print(detail.to_string())
    print("\n=== CKD-PERSISTENCE SUMMARY ===")
    print(summary.to_string())
    print("\n=== RESULT ===")
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    indir = sys.argv[1]
    consensus_tsv = sys.argv[2]
    out_dir = sys.argv[3]
    case = sys.argv[4] if len(sys.argv) > 4 else "CKD"
    main(indir, consensus_tsv, out_dir, case)
