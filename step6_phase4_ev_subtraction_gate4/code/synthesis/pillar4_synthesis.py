#!/usr/bin/env python3
"""
Pillar 4 - Constructive synthesis: whole-landscape EV subtraction.

Reuses ONLY frozen Pillar-3 artifacts (no new inference, no new Census pull).
Produces:
  - pillar4_funnel.json                  (strictly-nested evidence-refinement funnel, real counts)
  - pillar4_route_partition_primary.tsv  (3 consensus axes x route tag + evidence)
  - pillar4_whole_landscape_routes.tsv   (all reproducible S4 axes x route tag)
  - pillar4_supporting_landscape.tsv     ("Supporting reproducible communication landscape", frozen secondary rule)
  - pillar4_subtraction_summary.json     (route-tag counts, fraction EV-compatible, soluble-control-gate status)

Route-tag rule (FROZEN, reused, evaluated in this precedence):
  EV-compatible  iff (ligand in EV_TIER2_CARGO) AND (ligand NOT in NONEV_SECRETED_BACKGROUND)
  else CellChat annotation: Secreted Signaling -> soluble-classical
                            ECM-Receptor       -> ECM
                            Cell-Cell Contact  -> contact
  else (no frozen membership + no CellChat annotation) -> unclassified  (reported, not dropped)

SECONDARY (descriptive) rule for the Supporting reproducible communication landscape:
  cross-engine direction agreement (dir_both) in >=1 atlas
  AND statistically supported by >=1 engine (sig_ge1) in >=1 atlas
  AND reproduced in >=1 atlas
  AND NOT a primary consensus axis (i.e. does NOT pass full chain in BOTH atlases)
  == the S3-in->=1-atlas set minus the 3 consensus axes.
"""
import json, hashlib
from pathlib import Path
import pandas as pd

SEED = 1
PROJ = Path("/mnt/results/ev_aki_ckd_project")
P3 = PROJ / "03_communication_landscape" / "data"
P4 = PROJ / "04_synthesis" / "data"
EO = Path("/mnt/shared-workspace/shared/pillar3/engine_out")
P4.mkdir(parents=True, exist_ok=True)

CONSENSUS_KEYS = [
    "receiver_endothelial|SPP1|ITGA5_ITGB1",
    "receiver_stromal|SPP1|ITGAV_ITGB1",
    "receiver_stromal|COL4A1|ITGA1_ITGB1",
]

# ---------------------------------------------------------------- frozen sets
ev = json.load(open(PROJ / "00_preregistration" / "feasibility" / "ev_genesets_frozen_v1.3.json"))
EV_TIER2_CARGO = set(ev["EV_TIER2_CARGO"])
NONEV_BG = set(ev["NONEV_SECRETED_BACKGROUND"])

def ev_compatible(ligand: str) -> bool:
    """Frozen EV route rule. Multi-subunit ligands: EV-compatible iff ANY subunit qualifies."""
    subs = str(ligand).replace("_", " ").split()
    return any((s in EV_TIER2_CARGO) and (s not in NONEV_BG) for s in subs)

# ---------------------------------------------- CellChat annotation lookup
# authoritative source for Secreted/ECM/Contact tag (the differential tables dropped this col)
ann_rows = []
for f in sorted(EO.glob("cellchat_LR_*.tsv")):
    d = pd.read_csv(f, sep="\t")
    ann_rows.append(d[["ligand", "receptor", "pathway_name", "annotation"]])
ann = pd.concat(ann_rows, ignore_index=True).drop_duplicates()

# key on (ligand, receptor); receptor names in CellChat use A_B for complexes (matches our keys)
ann_lookup = {}
for _, r in ann.iterrows():
    ann_lookup.setdefault((str(r["ligand"]), str(r["receptor"])), str(r["annotation"]))
# fallback lookup on ligand alone (majority annotation) for axes not matched on the pair
lig_ann = (ann.groupby("ligand")["annotation"]
             .agg(lambda s: s.value_counts().index[0]).to_dict())

ANN_MAP = {
    "Secreted Signaling": "Soluble-classical",
    "ECM-Receptor": "ECM",
    "Cell-Cell Contact": "Contact",
}

def route_tag(ligand: str, receptor: str):
    """Return (route_tag, ev_compatible_bool, annotation_source)."""
    if ev_compatible(ligand):
        return "EV-compatible", True, "frozen_ev_rule"
    a = ann_lookup.get((str(ligand), str(receptor)))
    src = "cellchat_pair"
    if a is None:
        a = lig_ann.get(str(ligand))
        src = "cellchat_ligand"
    if a in ANN_MAP:
        return ANN_MAP[a], False, src
    return "unclassified", False, "none"

# ------------------------------------------------------- per-axis flag table
pa = pd.read_csv(P3 / "consensus" / "per_axis_per_dataset.tsv", sep="\t")
assert bool((( pa["dir_both"] & pa["sig_ge1"] & pa["floor_both"]) == pa["passes_in_dataset"]).all())

g = pa.groupby("key")
n_dir  = g["dir_both"].sum()                      # datasets with dir concordance
n_ds   = g.apply(lambda d: int((d["dir_both"] & d["sig_ge1"]).sum()))   # dir+sig
n_full = g["passes_in_dataset"].sum()             # dir+sig+floor

set_S2 = set(n_dir[n_dir >= 1].index)
set_S3 = set(n_ds[n_ds >= 1].index)
set_S4 = set(n_full[n_full >= 1].index)
set_S5 = set(n_full[n_full >= 2].index)           # consensus (both atlases)

S1 = pa["key"].nunique()
funnel = {
    "cap": 2000,
    "note": "Strictly-nested; all stages at '>=1 atlas' scope until the explicit cross-atlas step (S5).",
    "stages": [
        {"id": "S1", "label": "Candidate PT->receiver interactions", "n": int(S1),
         "criterion": "unique axis in any engine, any atlas"},
        {"id": "S2", "label": "Cross-engine direction agreement", "n": int(len(set_S2)),
         "criterion": "up-in-AKI in BOTH engines (CellChat & LIANA), >=1 atlas"},
        {"id": "S3", "label": "+ Statistical support", "n": int(len(set_S3)),
         "criterion": "additionally significant in >=1 engine, >=1 atlas"},
        {"id": "S4", "label": "+ Effect-size floor", "n": int(len(set_S4)),
         "criterion": "additionally passes a-priori absolute+relative floor, >=1 atlas (full single-atlas pass)"},
        {"id": "S5", "label": "Cross-atlas replication (consensus)", "n": int(len(set_S5)),
         "criterion": "full evidence chain in BOTH atlases"},
    ],
}
# nesting + monotonic guards
assert set_S3 <= set_S2 and set_S4 <= set_S3 and set_S5 <= set_S4, "funnel stages not nested"
seq = [s["n"] for s in funnel["stages"]]
assert all(seq[i] >= seq[i+1] for i in range(len(seq)-1)), "funnel not monotonic"
assert set(CONSENSUS_KEYS) == set_S5, "consensus set mismatch"

# ------------------------------------------- route-tag the whole S4 landscape
cons_land = pd.read_csv(P3 / "consensus" / "consensus_landscape_AKIvsnormal.tsv", sep="\t")

def key_parts(k):
    recv, lig, rec = k.split("|")
    return recv.replace("receiver_", ""), lig, rec

rows = []
for k in sorted(set_S4):
    recv, lig, rec = key_parts(k)
    tag, evc, src = route_tag(lig, rec)
    rows.append({"key": k, "receiver": recv, "ligand": lig, "receptor": rec,
                 "route_tag": tag, "ev_compatible": evc, "annotation_source": src,
                 "n_datasets_full_pass": int(n_full[k]), "is_consensus": k in set_S5})
land = pd.DataFrame(rows).sort_values(["is_consensus", "route_tag", "key"], ascending=[False, True, True])
land.to_csv(P4 / "pillar4_whole_landscape_routes.tsv", sep="\t", index=False)

# ---------------------------------------------------- primary (3 consensus)
prim = land[land["is_consensus"]].copy()
# attach the frozen consensus effect-size columns for the paper table
prim = prim.merge(
    cons_land[["key", "pathway_name", "cc_dprob_mean", "liana_rank_improvement_mean",
               "liana_rra_improvement_mean"]],
    on="key", how="left")
prim["tier"] = prim["key"].map({
    "receiver_endothelial|SPP1|ITGA5_ITGB1": "Tier 1 (cap-robust)",
    "receiver_stromal|SPP1|ITGAV_ITGB1": "Tier 1 (cap-robust)",
    "receiver_stromal|COL4A1|ITGA1_ITGB1": "Tier 2 (magnitude-borderline at cap4000)",
})
prim.to_csv(P4 / "pillar4_route_partition_primary.tsv", sep="\t", index=False)

# ------------------------------ Supporting reproducible communication landscape
# frozen secondary rule: S3 (dir+sig, >=1 atlas) minus consensus
supp_keys = sorted(set_S3 - set_S5)
srows = []
for k in supp_keys:
    recv, lig, rec = key_parts(k)
    tag, evc, src = route_tag(lig, rec)
    srows.append({"key": k, "receiver": recv, "ligand": lig, "receptor": rec,
                  "route_tag": tag, "ev_compatible": evc, "annotation_source": src,
                  "n_datasets_dir": int(n_dir[k]),
                  "n_datasets_dir_sig": int(n_ds[k]),
                  "n_datasets_full_pass": int(n_full[k]),
                  "reaches_floor_ge1_atlas": k in set_S4})
supp = pd.DataFrame(srows).sort_values(["route_tag", "key"])
supp.to_csv(P4 / "pillar4_supporting_landscape.tsv", sep="\t", index=False)

# ----------------------------------------------- subtraction summary + gate
def route_counts(df):
    return df["route_tag"].value_counts().to_dict()

neg_gate = json.load(open(P3 / "overlay" / "ev_route_negcontrol_gate.json"))

summary = {
    "seed": SEED,
    "cap": 2000,
    "ev_geneset_version": "v1.3",
    "ev_geneset_sha16": hashlib.sha256(
        json.dumps(ev, sort_keys=True).encode()).hexdigest()[:16],
    "reproducible_population_S4": {
        "n": int(len(set_S4)),
        "route_tag_counts": route_counts(land),
        "n_ev_compatible": int(land["ev_compatible"].sum()),
        "frac_ev_compatible": round(float(land["ev_compatible"].mean()), 4),
    },
    "primary_consensus": {
        "n": int(len(prim)),
        "route_tag_counts": route_counts(prim),
        "n_ev_compatible": int(prim["ev_compatible"].sum()),
    },
    "supporting_reproducible_communication_landscape": {
        "n": int(len(supp)),
        "route_tag_counts": route_counts(supp),
        "n_ev_compatible": int(supp["ev_compatible"].sum()),
        "definition": "dir_both & sig_ge1 in >=1 atlas, NOT consensus (frozen secondary rule)",
    },
    "soluble_negative_control_gate": {
        "n_controls": neg_gate["n_soluble_negative_controls"],
        "n_scored_ev_compatible": neg_gate["n_controls_scored_EV_compatible"],
        "informative": neg_gate["ev_compatibility_annotation_informative"],
    },
}
json.dump(summary, open(P4 / "pillar4_subtraction_summary.json", "w"), indent=2)
json.dump(funnel, open(P4 / "pillar4_funnel.json", "w"), indent=2)

# ----------------------------------------------------------------- report
print("=== FUNNEL (strictly nested, cap=2000) ===")
for s in funnel["stages"]:
    print(f"  {s['id']}: {s['n']:>5}  {s['label']}")
print(f"  --> route: EV-compatible among consensus = {summary['primary_consensus']['n_ev_compatible']}")
print(f"  --> soluble/ECM survivors among consensus = "
      f"{sum(v for k,v in summary['primary_consensus']['route_tag_counts'].items() if k in ('Soluble-classical','ECM'))}")
print()
print("=== WHOLE REPRODUCIBLE LANDSCAPE (S4 =", len(set_S4), ") route tags ===")
for k, v in summary["reproducible_population_S4"]["route_tag_counts"].items():
    print(f"  {k:>18}: {v}")
print(f"  EV-compatible fraction: {summary['reproducible_population_S4']['frac_ev_compatible']}")
print()
print("=== PRIMARY consensus route tags ===", summary["primary_consensus"]["route_tag_counts"])
print("=== SUPPORTING landscape (n =", len(supp), ") route tags ===",
      summary["supporting_reproducible_communication_landscape"]["route_tag_counts"])
print("=== soluble neg-control gate: ",
      f"{summary['soluble_negative_control_gate']['n_scored_ev_compatible']}/"
      f"{summary['soluble_negative_control_gate']['n_controls']} flagged ===")
print("\nWROTE:", *[p.name for p in sorted(P4.glob('pillar4_*'))])
