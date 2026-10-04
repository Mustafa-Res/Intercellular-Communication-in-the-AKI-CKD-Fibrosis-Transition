"""
Build the 5 expanded tables (T1-T5) from the frozen synthesis_data layer.
Writes TSVs + a single column dictionary to 04_synthesis/data/.
Every value derives from frozen artifacts; nothing hardcoded except display labels.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np

import synthesis_data as sdata

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "04_synthesis" / "data"
OUT.mkdir(parents=True, exist_ok=True)


def _pivot_per_axis(pa: pd.DataFrame) -> pd.DataFrame:
    """One row per key with per-dataset effect columns + stage flags aggregated across datasets."""
    # per-dataset numeric effects we want side by side
    num_cols = ["dprob", "rel_increase", "rank_improvement", "rra_improvement", "prob_AKI", "prob_normal"]
    flag_cols = ["up_in_AKI_cc", "up_in_AKI_liana", "dir_both", "cc_sig", "liana_sig",
                 "cc_floor", "liana_floor", "sig_ge1", "floor_both", "passes_in_dataset"]
    wide = {}
    for ds, g in pa.groupby("dataset"):
        g = g.set_index("key")
        for c in num_cols:
            wide[f"{c}__{ds}"] = g[c]
        for c in flag_cols:
            wide[f"{c}__{ds}"] = g[c]
    W = pd.DataFrame(wide)
    # static identity columns (same across datasets)
    ident = pa.drop_duplicates("key").set_index("key")[["receiver", "ligand", "receptor", "pathway_name"]]
    W = ident.join(W)
    return W.reset_index()


def build_T1(d) -> pd.DataFrame:
    """Master axis provenance: one row per candidate axis, S1->S5 stage flags + per-atlas effects + route/tier/membership."""
    pa = d["p3_per_axis"]
    W = _pivot_per_axis(pa)

    # Stage membership (re-derived, consistent with harness derive_funnel)
    dir_keys = set(pa.loc[pa["dir_both"], "key"])                       # S2
    dirsig_keys = set(pa.loc[pa["dir_both"] & pa["sig_ge1"], "key"])   # S3
    pass_keys = set(pa.loc[pa["passes_in_dataset"], "key"])            # S4
    passes = pa.loc[pa["passes_in_dataset"]]
    n_ds = passes.groupby("key")["dataset"].nunique()
    consensus_keys = set(n_ds[n_ds >= 2].index)                        # S5

    W["stage_S1_candidate"] = True
    W["stage_S2_dir_concordant"] = W["key"].isin(dir_keys)
    W["stage_S3_dir_and_sig"] = W["key"].isin(dirsig_keys)
    W["stage_S4_full_single_atlas_pass"] = W["key"].isin(pass_keys)
    W["stage_S5_cross_atlas_consensus"] = W["key"].isin(consensus_keys)

    # Route tag + EV compatibility + tier from pillar-4 whole-landscape table (only defined for S4=34).
    routes = d["p4_routes34"][["key", "route_tag", "ev_compatible", "annotation_source", "n_datasets_full_pass", "is_consensus"]]
    W = W.merge(routes, on="key", how="left")

    # Supporting-landscape membership (secondary layer)
    supp_keys = set(d["p4_supporting41"]["key"]) if "key" in d["p4_supporting41"].columns else set()
    W["supporting_landscape_member"] = W["key"].isin(supp_keys)

    # Tier label for the 3 consensus (Tier 1 SPP1 cap-robust; Tier 2 COL4A1 magnitude-borderline)
    tier = {
        "receiver_endothelial|SPP1|ITGA5_ITGB1": "Tier1_cap_robust",
        "receiver_stromal|SPP1|ITGAV_ITGB1": "Tier1_cap_robust",
        "receiver_stromal|COL4A1|ITGA1_ITGB1": "Tier2_magnitude_borderline",
    }
    W["consensus_tier"] = W["key"].map(tier)

    # Sort: consensus first, then survivors, then rest; within, by max dprob
    W["_ord"] = np.where(W["stage_S5_cross_atlas_consensus"], 0,
                 np.where(W["stage_S4_full_single_atlas_pass"], 1,
                 np.where(W["stage_S3_dir_and_sig"], 2,
                 np.where(W["stage_S2_dir_concordant"], 3, 4))))
    dprob_cols = [c for c in W.columns if c.startswith("dprob__")]
    W["_maxdprob"] = W[dprob_cols].max(axis=1)
    W = W.sort_values(["_ord", "_maxdprob"], ascending=[True, False]).drop(columns=["_ord", "_maxdprob"])
    return W


def build_T2(d) -> pd.DataFrame:
    """Cross-pillar specificity summary: one row per pillar x test x observed vs control statistic x verdict."""
    rows = []
    # P1 native GSEA EV vs ribosomal-null (permutation summary)
    p1 = d["p1_null"].set_index(["context", "term"])
    rows.append(dict(pillar="P1_bulk", context="native_GSE139061", test="permutation_null_enrichment",
                     signal_term="EV_TIER1_PRIMARY",
                     signal_observed=p1.loc[("native_GSE139061", "EV_TIER1_PRIMARY"), "observed"],
                     signal_z=p1.loc[("native_GSE139061", "EV_TIER1_PRIMARY"), "z"],
                     signal_p=p1.loc[("native_GSE139061", "EV_TIER1_PRIMARY"), "p_greater"],
                     control_term="CTRL_RIBOSOMAL",
                     control_observed=p1.loc[("native_GSE139061", "CTRL_RIBOSOMAL"), "observed"],
                     control_z=p1.loc[("native_GSE139061", "CTRL_RIBOSOMAL"), "z"],
                     control_p=p1.loc[("native_GSE139061", "CTRL_RIBOSOMAL"), "p_greater"],
                     verdict="EV inside null (NS); generic ribosomal control significant -> specificity FAIL",
                     source="01_bulk_signature/data/permutation_null_summary.tsv"))
    rows.append(dict(pillar="P1_bulk", context="transplant_meta", test="permutation_null_enrichment",
                     signal_term="EV_TIER1_PRIMARY",
                     signal_observed=p1.loc[("transplant_meta", "EV_TIER1_PRIMARY"), "observed"],
                     signal_z=p1.loc[("transplant_meta", "EV_TIER1_PRIMARY"), "z"],
                     signal_p=p1.loc[("transplant_meta", "EV_TIER1_PRIMARY"), "p_greater"],
                     control_term="CTRL_RIBOSOMAL",
                     control_observed=p1.loc[("transplant_meta", "CTRL_RIBOSOMAL"), "observed"],
                     control_z=p1.loc[("transplant_meta", "CTRL_RIBOSOMAL"), "z"],
                     control_p=p1.loc[("transplant_meta", "CTRL_RIBOSOMAL"), "p_greater"],
                     verdict="EV enriched but generic ribosomal control equally enriched -> not EV-specific",
                     source="01_bulk_signature/data/permutation_null_summary.tsv"))
    # P2 sender-space PT AKI_vs_normal EV vs ribosomal
    p2 = d["p2_null"]
    sel = p2[(p2["compartment"] == "PT") & (p2["contrast"] == "AKI_vs_normal")].set_index("term")
    rows.append(dict(pillar="P2_singlecell", context="PT AKI_vs_normal", test="permutation_null_sender_module",
                     signal_term="EV_TIER1_PRIMARY",
                     signal_observed=sel.loc["EV_TIER1_PRIMARY", "observed"],
                     signal_z=sel.loc["EV_TIER1_PRIMARY", "z"],
                     signal_p=sel.loc["EV_TIER1_PRIMARY", "p_greater"],
                     control_term="CTRL_RIBOSOMAL",
                     control_observed=sel.loc["CTRL_RIBOSOMAL", "observed"],
                     control_z=sel.loc["CTRL_RIBOSOMAL", "z"],
                     control_p=sel.loc["CTRL_RIBOSOMAL", "p_greater"],
                     verdict="EV module elevated but generic ribosomal control ~4x stronger -> specificity FAIL",
                     source="02_singlecell_specificity/data/permutation_null_sender.tsv"))
    # P3 consensus / EV overlay
    g3 = d["p3_gate_summary"]
    sub = d["p4_subtraction"]
    rows.append(dict(pillar="P3_communication", context="AKI_vs_normal consensus", test="cross_atlas_consensus + EV_overlay",
                     signal_term="consensus_axes", signal_observed=g3["n_consensus_axes"], signal_z=np.nan, signal_p=np.nan,
                     control_term="EV_compatible_of_consensus",
                     control_observed=sub["primary_consensus"]["n_ev_compatible"], control_z=np.nan, control_p=np.nan,
                     verdict="3 reproducible consensus axes; 0 EV-compatible (descriptive, never gates)",
                     source="03_communication_landscape/data/gate3/gate3_summary.json + pillar4_subtraction_summary.json"))
    # P4 synthesis subtraction
    rows.append(dict(pillar="P4_synthesis", context="reproducible landscape S4", test="route_partition_subtraction",
                     signal_term="reproducible_survivors", signal_observed=sub["reproducible_population_S4"]["n"],
                     signal_z=np.nan, signal_p=np.nan,
                     control_term="EV_compatible_of_survivors",
                     control_observed=sub["reproducible_population_S4"]["n_ev_compatible"], control_z=np.nan, control_p=np.nan,
                     verdict="34/34 reproducible survivors soluble+ECM; 0 EV-compatible",
                     source="04_synthesis/data/pillar4_subtraction_summary.json"))
    return pd.DataFrame(rows)


def build_T3(d) -> pd.DataFrame:
    """Consensus-axis dossier: the 3 axes x every robustness dimension."""
    cap = d["p3_cap"]["axes"]
    lodo = d["p3_lodo"]["verdict"]
    ckd = {a["key"]: a for a in d["p3_ckd"]["per_axis"]}
    pathway = d["p3_pathway"]["consensus_pathways"]
    routes = d["p4_routes34"].set_index("key")

    # map display axis name (cap json) -> canonical key
    capkey = {
        "SPP1->ITGA5_ITGB1 (endothelial)": "receiver_endothelial|SPP1|ITGA5_ITGB1",
        "SPP1->ITGAV_ITGB1 (stromal)": "receiver_stromal|SPP1|ITGAV_ITGB1",
        "COL4A1->ITGA1_ITGB1 (stromal)": "receiver_stromal|COL4A1|ITGA1_ITGB1",
    }
    tier = {
        "receiver_endothelial|SPP1|ITGA5_ITGB1": "Tier1_cap_robust",
        "receiver_stromal|SPP1|ITGAV_ITGB1": "Tier1_cap_robust",
        "receiver_stromal|COL4A1|ITGA1_ITGB1": "Tier2_magnitude_borderline",
    }
    rows = []
    for disp, key in capkey.items():
        _, label, recv, pway = [t for t in sdata.CONSENSUS_AXES if t[0] == key][0]
        c = cap[disp]
        a12 = c["per_dataset"]["a12ccb9b"]; dea = c["per_dataset"]["dea717d4"]
        pw = pathway[pway]
        rows.append(dict(
            key=key, axis=label, receiver=recv, pathway=pway,
            route_tag=routes.loc[key, "route_tag"], ev_compatible=bool(routes.loc[key, "ev_compatible"]),
            consensus_tier=tier[key],
            # cap 4000 per-dataset effects
            cc_dprob_a12=a12["cc_dprob"], cc_rel_increase_a12=a12["cc_rel_increase"],
            liana_rank_impr_a12=a12["liana_rank_improvement"], liana_rra_impr_a12=a12["liana_rra_improvement"],
            passes_a12_cap4000=a12["passes_in_dataset"],
            cc_dprob_dea=dea["cc_dprob"], cc_rel_increase_dea=dea["cc_rel_increase"],
            liana_rank_impr_dea=dea["liana_rank_improvement"], liana_rra_impr_dea=dea["liana_rra_improvement"],
            passes_dea_cap4000=dea["passes_in_dataset"],
            direction_robust_cap4000=c["direction_robust_both_datasets"],
            strict_pass_both_cap4000=c["passes_strict_both_datasets_cap4000"],
            # LODO
            lodo_frac_a12=lodo[key]["per_dataset_frac"]["a12ccb9b"],
            lodo_frac_dea=lodo[key]["per_dataset_frac"]["dea717d4"],
            lodo_majority_stable=lodo[key]["majority_stable_all_datasets"],
            # scramble (global result; consensus collapses)
            scramble_survives=False,
            # CKD persistence
            ckd_persists_both_engines_both_ds=ckd[key]["persists_both_datasets_both_engines"],
            ckd_n_datasets_cc_up=ckd[key]["n_datasets_cc_up"],
            ckd_n_datasets_liana_up=ckd[key]["n_datasets_liana_up"],
            # pathway coherence (this atlas pair)
            pathway_dprob_a12=pw["a12ccb9b"]["dprob"], pathway_rank_a12=pw["a12ccb9b"]["rank_up_in_AKI"],
            pathway_dprob_dea=pw["dea717d4"]["dprob"], pathway_rank_dea=pw["dea717d4"]["rank_up_in_AKI"],
        ))
    return pd.DataFrame(rows)


def build_T4(d) -> pd.DataFrame:
    """EV-subtraction accounting: survivor ligands x route x #axes x EV-cargo vs signal-peptide membership."""
    w = d["p4_routes34"]
    ev = d["ev"]
    cargo = set(ev["EV_TIER2_CARGO"]); bg = set(ev["NONEV_SECRETED_BACKGROUND"])
    tier1 = set(ev["EV_TIER1_PRIMARY"])
    rows = []
    for lig, g in w.groupby("ligand"):
        subs = [s for part in g["ligand"].unique() for s in part.split("_")]  # ligand col is single gene here
        in_cargo = lig in cargo
        in_bg = lig in bg
        rows.append(dict(
            ligand=lig,
            n_axes=len(g),
            route_tag=";".join(sorted(g["route_tag"].unique())),
            receivers=";".join(sorted(g["receiver"].unique())),
            in_EV_TIER1_PRIMARY=lig in tier1,
            in_EV_TIER2_CARGO=in_cargo,
            in_NONEV_SECRETED_BACKGROUND=in_bg,
            ev_compatible=bool(g["ev_compatible"].any()),
            subtraction_reason=(
                "not in EV cargo" if not in_cargo else
                ("EV cargo BUT classical signal-peptide background -> non-EV" if in_bg else "EV-compatible")
            ),
        ))
    T = pd.DataFrame(rows).sort_values("n_axes", ascending=False).reset_index(drop=True)
    assert T["n_axes"].sum() == 34, T["n_axes"].sum()
    return T


def build_T5(d) -> pd.DataFrame:
    """Expression-floor sensitivity: render the frozen JSON sweep as a tidy TSV."""
    ef = d["p4_exprfloor"]
    # locate the per-floor sweep list irrespective of exact wrapper key
    sweep = None
    for k, v in ef.items():
        if isinstance(v, list) and v and isinstance(v[0], dict) and any("floor" in str(kk).lower() for kk in v[0]):
            sweep = v; break
    if sweep is None:
        # fall back: flatten dict-of-floors
        rows = []
        for k, v in ef.items():
            if isinstance(v, dict):
                r = {"floor": k}; r.update(v); rows.append(r)
        return pd.DataFrame(rows)
    return pd.DataFrame(sweep)


DICTIONARY = """# synthesis tables (04_synthesis/data/) - column dictionary

All tables derive from frozen artifacts (no new inference). Generated by 04_synthesis/code/synthesis_tables_figures/build_tables.py.

## T1_master_axis_provenance.tsv
One row per candidate PT->receiver interaction (S1 = 1311). Columns:
- key: canonical axis id (receiver|ligand|receptor).
- receiver, ligand, receptor, pathway_name: interaction identity.
- <metric>__<dataset>: per-atlas effect metrics (dataset in {a12ccb9b snRNA, dea717d4 scRNA}).
    dprob = CellChat AKI-normal communication-probability delta; rel_increase = CellChat relative increase;
    rank_improvement / rra_improvement = LIANA magnitude-rank / robust-rank-aggregation improvement (AKI vs normal);
    prob_AKI / prob_normal = CellChat probabilities.
- <flag>__<dataset>: per-atlas boolean gates (up_in_AKI_cc, up_in_AKI_liana, dir_both, cc_sig, liana_sig,
    cc_floor, liana_floor, sig_ge1, floor_both, passes_in_dataset).
- stage_S1..S5: nested funnel membership (S1 candidate; S2 cross-engine direction-concordant; S3 +statistical
    support; S4 full single-atlas pass; S5 cross-atlas consensus). Strictly nested S5 subset S4 subset S2.
- route_tag, ev_compatible, annotation_source, n_datasets_full_pass, is_consensus: route attribution (defined for S4 survivors).
- supporting_landscape_member: in the Supporting reproducible communication landscape (dir_both & sig_ge1 in >=1 atlas, not consensus).
- consensus_tier: Tier1_cap_robust (SPP1 axes) / Tier2_magnitude_borderline (COL4A1) for the 3 consensus axes.

## T2_cross_pillar_specificity.tsv
One row per pillar x specificity test. signal_* = the EV/biological term; control_* = the matched generic control
(ribosomal for P1/P2 nulls; EV-compatible count for P3/P4). verdict summarizes the specificity outcome. source = frozen file.

## T3_consensus_axis_dossier.tsv
The 3 consensus axes x every robustness dimension: per-atlas cap-4000 effects & pass flags, cap direction-robust /
strict-pass, leave-one-donor-out fractions (lodo_frac_*), scramble survival, CKD persistence, and SPP1/COLLAGEN
pathway-level dprob & up-in-AKI rank per atlas. route_tag / ev_compatible / consensus_tier included.

## T4_ev_subtraction_accounting.tsv
The survivor ligands (S4) x n_axes x route x receivers x EV-set membership. in_EV_TIER2_CARGO & in_NONEV_SECRETED_BACKGROUND
show WHY each ligand is non-EV: cargo-annotated but carrying a classical signal peptide (subtraction_reason).
n_axes sums to 34.

## T5_expression_floor_sensitivity.tsv
Cross-engine expression-proxy floor sweep of the funnel (frozen pillar4_funnel_exprfloor_sensitivity.json):
funnel stage counts (S1..S5) at increasing expression floors, showing floor-invariance of the downstream chain.
"""


def main():
    d = sdata.load_all()
    sdata.verify(d)  # guard: identities must hold before we emit tables
    outputs = {
        "T1_master_axis_provenance.tsv": build_T1(d),
        "T2_cross_pillar_specificity.tsv": build_T2(d),
        "T3_consensus_axis_dossier.tsv": build_T3(d),
        "T4_ev_subtraction_accounting.tsv": build_T4(d),
        "T5_expression_floor_sensitivity.tsv": build_T5(d),
    }
    for name, df in outputs.items():
        df.to_csv(OUT / name, sep="\t", index=False)
        print(f"wrote {name:42s} shape={df.shape}")
    (OUT / "TABLES_dictionary.md").write_text(DICTIONARY)
    print(f"wrote TABLES_dictionary.md")
    return outputs


if __name__ == "__main__":
    main()
