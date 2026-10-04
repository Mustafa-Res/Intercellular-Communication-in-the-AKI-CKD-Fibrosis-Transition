"""
Pillar-synthesis expanded-visualization data layer.

Single source of truth for every number rendered in the expanded figures/tables.
- Loads ONLY frozen artifacts (no new inference, no new data).
- Exposes a `load_all()` that returns a dict of parsed objects.
- Exposes `verify()` which asserts the re-derivation identities the plan (§11.5)
  requires, and returns a list of (source_file -> value) provenance rows for
  printing or logging.

Run standalone:  python3 04_synthesis/code/synthesis_tables_figures/synthesis_data.py
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import pandas as pd

# Project root = parents[3] of this file: 04_synthesis/code/synthesis_tables_figures/synthesis_data.py -> project root
ROOT = Path(__file__).resolve().parents[3]
RES = ROOT  # pillar data dirs are addressed repo-relatively


def _tsv(p):
    return pd.read_csv(RES / p, sep="\t")


def _json(p):
    with open(RES / p) as f:
        return json.load(f)


def load_all() -> dict:
    d = {}
    # Pillar 1
    d["p1_null"] = _tsv("01_bulk_signature/data/permutation_null_summary.tsv")
    d["p1_gate"] = _tsv("01_bulk_signature/data/gate1_verdicts.tsv")
    d["p1_native_gsea"] = _tsv("01_bulk_signature/data/ev_enrichment_native_gsea.tsv")
    d["p1_transplant_meta"] = _tsv("01_bulk_signature/data/analysisB_transplant_meta.tsv")
    # Pillar 2
    d["p2_null"] = _tsv("02_singlecell_specificity/data/permutation_null_sender.tsv")
    d["p2_sev"] = _tsv("02_singlecell_specificity/data/severity_within_PT.tsv")
    d["p2_sev_corr"] = _tsv("02_singlecell_specificity/data/severity_within_PT_correlation.tsv")
    d["p2_nes"] = _tsv("02_singlecell_specificity/data/per_dataset_pt_ev_nes.tsv")
    d["p2_gate"] = _tsv("02_singlecell_specificity/data/gate2_verdicts.tsv")
    # Pillar 3
    d["p3_per_axis"] = _tsv("03_communication_landscape/data/consensus/per_axis_per_dataset.tsv")
    d["p3_consensus"] = _tsv("03_communication_landscape/data/consensus/consensus_landscape_AKIvsnormal.tsv")
    d["p3_lodo"] = _json("03_communication_landscape/data/gate3/lodo_result.json")
    d["p3_cap"] = _json("03_communication_landscape/data/gate3/cap_sensitivity_result.json")
    d["p3_scramble"] = _json("03_communication_landscape/data/gate3/scramble_result.json")
    d["p3_ckd"] = _json("03_communication_landscape/data/gate3/ckd_persistence_result.json")
    d["p3_pathway"] = _json("03_communication_landscape/data/gate3/pathway_coherence_result.json")
    d["p3_gate_summary"] = _json("03_communication_landscape/data/gate3/gate3_summary.json")
    # Pillar 4
    d["p4_routes34"] = _tsv("04_synthesis/data/pillar4_whole_landscape_routes.tsv")
    d["p4_supporting41"] = _tsv("04_synthesis/data/pillar4_supporting_landscape.tsv")
    d["p4_primary3"] = _tsv("04_synthesis/data/pillar4_route_partition_primary.tsv")
    d["p4_subtraction"] = _json("04_synthesis/data/pillar4_subtraction_summary.json")
    d["p4_exprfloor"] = _json("04_synthesis/data/pillar4_funnel_exprfloor_sensitivity.json")
    # Frozen EV gene sets
    d["ev"] = _json("00_preregistration/feasibility/ev_genesets_frozen_v1.3.json")
    return d


# The three frozen consensus axes, canonical order + display labels.
CONSENSUS_AXES = [
    ("receiver_endothelial|SPP1|ITGA5_ITGB1", "SPP1->ITGA5_ITGB1", "endothelial", "SPP1"),
    ("receiver_stromal|SPP1|ITGAV_ITGB1",     "SPP1->ITGAV_ITGB1", "stromal",     "SPP1"),
    ("receiver_stromal|COL4A1|ITGA1_ITGB1",   "COL4A1->ITGA1_ITGB1", "stromal",   "COLLAGEN"),
]

# Survivor ligands (from whole-landscape route table); membership drives T4.
SURVIVOR_LIGANDS = ["SPP1", "COL4A1", "COL4A2", "LAMA5", "LAMB1", "LAMC1", "THBS1", "PROS1", "SLIT2"]


def derive_funnel(d: dict) -> dict:
    """Re-derive the nested funnel counts from per_axis_per_dataset (NOT hardcoded)."""
    pa = d["p3_per_axis"]
    S1 = pa["key"].nunique()  # all candidate axes
    # S2 = cross-engine direction-concordant in >=1 dataset (dir_both true in any dataset)
    S2 = pa.loc[pa["dir_both"], "key"].nunique()
    # S3 = + statistical support (dir_both & sig_ge1) in >=1 dataset
    S3 = pa.loc[pa["dir_both"] & pa["sig_ge1"], "key"].nunique()
    # S4 = full single-atlas pass (passes_in_dataset) in >=1 dataset
    S4 = pa.loc[pa["passes_in_dataset"], "key"].nunique()
    # S5 = cross-atlas consensus (pass in BOTH datasets)
    passes = pa.loc[pa["passes_in_dataset"]]
    n_ds = passes.groupby("key")["dataset"].nunique()
    S5 = int((n_ds >= 2).sum())
    return {"S1": S1, "S2": S2, "S3": S3, "S4": S4, "S5": S5}


def verify(d: dict) -> list[tuple[str, str, object]]:
    """
    Assert plan §11.5 identities. Returns provenance rows (label, source_file, value).
    Raises AssertionError if any identity fails.
    """
    prov: list[tuple[str, str, object]] = []

    def rec(label, src, val):
        prov.append((label, src, val))
        return val

    # --- Funnel re-derivation ---
    f = derive_funnel(d)
    rec("Funnel S1 (candidate axes)", "03_communication_landscape/data/consensus/per_axis_per_dataset.tsv", f["S1"])
    rec("Funnel S2 (dir-concordant)", "03_communication_landscape/data/consensus/per_axis_per_dataset.tsv", f["S2"])
    rec("Funnel S3 (dir+sig)", "03_communication_landscape/data/consensus/per_axis_per_dataset.tsv", f["S3"])
    rec("Funnel S4 (full single-atlas pass)", "03_communication_landscape/data/consensus/per_axis_per_dataset.tsv", f["S4"])
    rec("Funnel S5 (cross-atlas consensus)", "03_communication_landscape/data/consensus/per_axis_per_dataset.tsv", f["S5"])
    assert f["S1"] == 1311, f["S1"]
    assert f["S2"] == 44, f["S2"]
    assert f["S3"] == 44, f["S3"]        # S2==S3: statistical support removes nothing
    assert f["S4"] == 34, f["S4"]
    assert f["S5"] == 3, f["S5"]

    # --- Whole-landscape route counts (S4=34) ---
    w = d["p4_routes34"]
    rc = w["route_tag"].value_counts().to_dict()
    rec("S4 route counts", "04_synthesis/data/pillar4_whole_landscape_routes.tsv", rc)
    assert len(w) == 34
    assert rc.get("ECM", 0) == 24 and rc.get("Soluble-classical", 0) == 10
    assert int(w["ev_compatible"].sum()) == 0, "expected 0 EV-compatible in S4"
    assert 24 + 10 == len(w)

    # --- Subtraction summary cross-check ---
    sub = d["p4_subtraction"]
    rec("S4 n / route (subtraction json)", "04_synthesis/data/pillar4_subtraction_summary.json",
        (sub["reproducible_population_S4"]["n"], sub["reproducible_population_S4"]["route_tag_counts"]))
    rec("Consensus n / route", "04_synthesis/data/pillar4_subtraction_summary.json",
        (sub["primary_consensus"]["n"], sub["primary_consensus"]["route_tag_counts"]))
    rec("Supporting landscape n / route", "04_synthesis/data/pillar4_subtraction_summary.json",
        (sub["supporting_reproducible_communication_landscape"]["n"],
         sub["supporting_reproducible_communication_landscape"]["route_tag_counts"]))
    assert sub["reproducible_population_S4"]["n"] == 34
    assert sub["primary_consensus"]["n"] == 3
    assert sub["supporting_reproducible_communication_landscape"]["n"] == 41
    assert sub["reproducible_population_S4"]["n_ev_compatible"] == 0
    assert sub["primary_consensus"]["n_ev_compatible"] == 0
    assert sub["supporting_reproducible_communication_landscape"]["n_ev_compatible"] == 0

    # --- Nesting / partition identities ---
    supp = d["p4_supporting41"]
    prim = d["p4_primary3"]
    rec("Supporting = 44 - 3 = 41", "04_synthesis/data/pillar4_supporting_landscape.tsv", len(supp))
    assert len(supp) == 41
    assert len(prim) == 3
    # consensus keys subset of survivors subset of dir-concordant
    pa = d["p3_per_axis"]
    keys_S2 = set(pa.loc[pa["dir_both"], "key"])
    keys_S4 = set(pa.loc[pa["passes_in_dataset"], "key"])
    cons_keys = set(w.loc[w["is_consensus"], "key"])
    rec("Consensus keys (n)", "04_synthesis/data/pillar4_whole_landscape_routes.tsv", len(cons_keys))
    assert len(cons_keys) == 3
    assert cons_keys <= keys_S4 <= keys_S2, "nesting 3 subset 34 subset 44 violated"
    # supporting and consensus disjoint
    supp_keys = set(supp["key"]) if "key" in supp.columns else set()
    if supp_keys:
        assert cons_keys.isdisjoint(supp_keys), "consensus and supporting must be disjoint"

    # --- EV membership mechanism (T4): 8/9 survivor ligands are cargo AND background ---
    ev = d["ev"]
    cargo = set(ev["EV_TIER2_CARGO"]); bg = set(ev["NONEV_SECRETED_BACKGROUND"])
    both = [g for g in SURVIVOR_LIGANDS if g in cargo and g in bg]
    rec("Survivor ligands in cargo AND background", "ev_genesets_frozen_v1.3.json", both)
    assert set(both) == set(SURVIVOR_LIGANDS) - {"COL4A1"}, both  # COL4A1 not in cargo
    assert "COL4A1" not in cargo

    # --- LODO fractions present for all 3 axes both atlases ---
    lodo = d["p3_lodo"]["verdict"]
    assert len(lodo) == 3
    for k in lodo:
        assert set(lodo[k]["per_dataset_frac"].keys()) == {"a12ccb9b", "dea717d4"}

    # --- Cap sensitivity: 3 direction-robust, 2 strict-pass at cap4000 ---
    cap = d["p3_cap"]
    rec("cap4000 direction-robust / strict-pass", "03_communication_landscape/data/gate3/cap_sensitivity_result.json",
        (cap["n_axes_direction_robust_cap4000"], cap["n_axes_passes_strict_both_datasets_cap4000"]))
    assert cap["n_axes_direction_robust_cap4000"] == 3
    assert cap["n_axes_passes_strict_both_datasets_cap4000"] == 2

    # --- Scramble collapses; CKD persistence = 1 axis both-engines-both-ds ---
    assert d["p3_scramble"]["collapses"] is True
    for s in d["p3_scramble"]["seeds"]:
        assert s["n_consensus_scrambled"] == 0
    rec("Scramble consensus survivors (per seed)", "03_communication_landscape/data/gate3/scramble_result.json",
        [s["n_consensus_scrambled"] for s in d["p3_scramble"]["seeds"]])
    ckd = d["p3_ckd"]
    rec("CKD persist both-engines-both-datasets", "03_communication_landscape/data/gate3/ckd_persistence_result.json",
        ckd["n_persist_both_datasets_both_engines"])
    assert ckd["n_persist_both_datasets_both_engines"] == 1

    # --- Gate 3 overall FULL PASS ---
    assert d["p3_gate_summary"]["OVERALL_VERDICT"] == "FULL PASS"

    # --- Severity correlation present (all_donors strong positive) ---
    corr = d["p2_sev_corr"].set_index("subset")
    rho_all = float(corr.loc["all_donors", "spearman_rho"])
    rec("Within-PT EV-module vs injury rho (all donors)", "02_singlecell_specificity/data/severity_within_PT_correlation.tsv", rho_all)
    assert rho_all > 0.5

    return prov


if __name__ == "__main__":
    data = load_all()
    provenance = verify(data)
    print(f"[synthesis_data] loaded {len(data)} frozen objects; {len(provenance)} identity checks PASSED\n")
    print(f"{'CHECK':52s} {'SOURCE':52s} VALUE")
    print("-" * 130)
    for label, src, val in provenance:
        print(f"{label:52s} {src:52s} {val}")
