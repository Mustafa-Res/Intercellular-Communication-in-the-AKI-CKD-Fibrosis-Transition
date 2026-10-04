#!/usr/bin/env python3
"""
Gate 4 - synthesis-integrity gate for Pillar 4 (constructive synthesis).

NOT a hypothesis-rejection gate. PASS means the synthesis is internally sound and
faithful to the frozen artifacts and to the locked Gate-1/2/3 verdicts.

  4A  Faithful subtraction   : frozen EV rule + soluble-control gate (0/38); no post-hoc reclassification.
  4B  Primary/secondary wall : 3 consensus axes carry the primary claim; the "Supporting reproducible
                               communication landscape" is a disjoint, separately-admitted set; naming constraint.
  4C  Quantitative dominance : among reproducible survivors soluble+ECM dominate, EV-compatible == 0.
  4D  Cross-pillar arc       : P1/P2 share the ribosomal specificity-failure mode; P2 effect-floor fails;
                               P3 is FULL PASS yet 0/3 EV-compatible; P4 survivors are soluble+ECM.
                               Checked against the ACTUAL gate1/2/3 artifacts, not memory.

Emits: 04_synthesis/data/gate4/gate4_verdicts.tsv, 04_synthesis/data/gate4/gate4_summary.json
"""
import json, re
from pathlib import Path
import pandas as pd

PROJ = Path("/mnt/results/ev_aki_ckd_project")
P4 = PROJ / "04_synthesis" / "data"
G4 = P4 / "gate4"
G4.mkdir(parents=True, exist_ok=True)

CONSENSUS_KEYS = {
    "receiver_endothelial|SPP1|ITGA5_ITGB1",
    "receiver_stromal|SPP1|ITGAV_ITGB1",
    "receiver_stromal|COL4A1|ITGA1_ITGB1",
}
FORBIDDEN_NAMES = ["the surviving program", "surviving program"]
REQUIRED_NAME = "Supporting reproducible communication landscape"

# ---- load Pillar-4 outputs
summary = json.load(open(P4 / "pillar4_subtraction_summary.json"))
funnel = json.load(open(P4 / "pillar4_funnel.json"))
land = pd.read_csv(P4 / "pillar4_whole_landscape_routes.tsv", sep="\t")
prim = pd.read_csv(P4 / "pillar4_route_partition_primary.tsv", sep="\t")
supp = pd.read_csv(P4 / "pillar4_supporting_landscape.tsv", sep="\t")

# ---- load prior gate artifacts (authoritative, on disk)
g1 = pd.read_csv(PROJ / "01_bulk_signature" / "data" / "gate1_verdicts.tsv", sep="\t")
g2 = pd.read_csv(PROJ / "02_singlecell_specificity" / "data" / "gate2_verdicts.tsv", sep="\t")
g3s = json.load(open(PROJ / "03_communication_landscape" / "data" / "gate3" / "gate3_summary.json"))

verdicts = []
def add(cid, status, evidence):
    verdicts.append({"criterion": cid, "status": "PASS" if status else "FAIL",
                     "evidence": evidence})

# ---------------------------------------------------------------- 4A
soluble_gate_ok = (summary["soluble_negative_control_gate"]["n_scored_ev_compatible"] == 0
                   and summary["soluble_negative_control_gate"]["informative"])
# no post-hoc reclassification: every ev_compatible tag must be derivable from the frozen rule flag
#   (annotation_source == 'frozen_ev_rule' whenever ev_compatible True)
no_posthoc = bool(((~land["ev_compatible"]) | (land["annotation_source"] == "frozen_ev_rule")).all())
uses_frozen = summary["ev_geneset_version"] == "v1.3"
a_ok = soluble_gate_ok and no_posthoc and uses_frozen
add("4A_faithful_subtraction", a_ok,
    f"soluble-control gate {summary['soluble_negative_control_gate']['n_scored_ev_compatible']}/"
    f"{summary['soluble_negative_control_gate']['n_controls']} flagged & informative={soluble_gate_ok}; "
    f"no post-hoc reclassification={no_posthoc}; frozen EV geneset v1.3={uses_frozen} "
    f"(sha {summary['ev_geneset_sha16']}).")

# ---------------------------------------------------------------- 4B
prim_keys = set(prim["key"])
supp_keys = set(supp["key"])
primary_ok = prim_keys == CONSENSUS_KEYS
disjoint = len(prim_keys & supp_keys) == 0
# naming constraint: scan all Pillar-4 text deliverables (memo + mini-paper) + this repo's outputs.
# Mini-paper lives at 04_synthesis/data/pillar4_mini_paper.md (co-located with the pillar).
naming_violations = []
naming_scope = (list(P4.glob("*.tsv")) + list(P4.glob("*.json"))
                + [P4 / "GATE4_memo.md", P4 / "pillar4_mini_paper.md"])
for f in naming_scope:
    if not f.exists():
        continue
    t = f.read_text(errors="ignore").lower()
    for bad in FORBIDDEN_NAMES:
        if bad in t:
            naming_violations.append(f"{f.name}:{bad}")
naming_ok = len(naming_violations) == 0
b_ok = primary_ok and disjoint and naming_ok
add("4B_primary_secondary_wall", b_ok,
    f"primary==3 consensus axes={primary_ok}; primary/secondary disjoint={disjoint} "
    f"(|primary|={len(prim_keys)}, |supporting|={len(supp_keys)}); "
    f"naming-constraint-clean (tables+json+memo+mini-paper)={naming_ok} (violations={naming_violations}).")

# ---------------------------------------------------------------- 4C
s4 = summary["reproducible_population_S4"]
rc = s4["route_tag_counts"]
soluble_ecm = rc.get("Soluble-classical", 0) + rc.get("ECM", 0)
ev_zero = s4["n_ev_compatible"] == 0
dominate = soluble_ecm == s4["n"] and ev_zero  # all survivors soluble/ECM, EV 0
prim_ev_zero = summary["primary_consensus"]["n_ev_compatible"] == 0
c_ok = dominate and prim_ev_zero
add("4C_quantitative_dominance", c_ok,
    f"S4 reproducible survivors n={s4['n']}: soluble+ECM={soluble_ecm}, EV-compatible={s4['n_ev_compatible']} "
    f"(frac_ev={s4['frac_ev_compatible']}); primary consensus EV-compatible={summary['primary_consensus']['n_ev_compatible']}. "
    f"Route counts={rc}.")

# ---------------------------------------------------------------- 4D
def verdict_of(df, key_substr, col_status):
    row = df[df.iloc[:, 0].astype(str).str.contains(key_substr, case=False, regex=False)]
    return row

# P1: ribosomal null clean must be FAIL (ribosomal VIOLATED) and overall 1C FAIL
g1_ribo = g1[g1["criterion"].str.contains("neg-control ladder", case=False)]
p1_ribo_fail = (g1_ribo["verdict"].iloc[0] == "FAIL") if len(g1_ribo) else False
g1_1c = g1[g1["criterion"].str.contains("GATE 1C", case=False)]
p1_1c_fail = (g1_1c["verdict"].iloc[0] == "FAIL") if len(g1_1c) else False

# P2: ribosomal specificity FAIL (2B) and effect_size_floor FAIL and overall FAIL
def g2_row(sub):
    r = g2[g2["criterion"].str.contains(sub, case=False, regex=False)]
    return (str(r["pass_"].iloc[0]).strip().lower() == "true") if len(r) else None
p2_2b = g2_row("2B_specificity")            # True == passed; we need this False (ribosomal violated)
p2_floor = g2_row("effect_size_floor")      # need False
p2_overall = g2_row("OVERALL_VERDICT")      # need False
p2_ribo_violated = (p2_2b is False)
p2_floor_fail = (p2_floor is False)
p2_fail = (p2_overall is False)

# shared failure mode: both P1 and P2 fail specifically on the ribosomal negative control
shared_ribo_mode = p1_ribo_fail and p2_ribo_violated

# P3: FULL PASS but 0/3 EV-compatible
p3_fullpass = g3s["OVERALL_VERDICT"] == "FULL PASS"
p3_ev0 = "0 of 3" in g3s["ev_overlay_note"]

# P4: survivors soluble+ECM (from 4C)
p4_survivors_soluble_ecm = dominate

d_ok = all([p1_ribo_fail, p1_1c_fail, p2_ribo_violated, p2_floor_fail, p2_fail,
            shared_ribo_mode, p3_fullpass, p3_ev0, p4_survivors_soluble_ecm])
add("4D_cross_pillar_arc", d_ok,
    f"P1 ribosomal-null FAIL={p1_ribo_fail} & 1C FAIL={p1_1c_fail}; "
    f"P2 ribosomal VIOLATED={p2_ribo_violated} & effect-floor FAIL={p2_floor_fail} & overall FAIL={p2_fail}; "
    f"SHARED ribosomal specificity-failure mode P1&P2={shared_ribo_mode}; "
    f"P3 FULL PASS={p3_fullpass} yet 0/3 EV-compatible={p3_ev0}; "
    f"P4 survivors soluble+ECM & EV=0={p4_survivors_soluble_ecm}. "
    f"No contradiction across locked gate verdicts.")

# ---------------------------------------------------------------- overall
vt = pd.DataFrame(verdicts)
overall = bool((vt["status"] == "PASS").all())
vt = pd.concat([vt, pd.DataFrame([{
    "criterion": "GATE4_OVERALL", "status": "PASS" if overall else "FAIL",
    "evidence": "Gate 4 = 4A AND 4B AND 4C AND 4D (synthesis-integrity)."}])],
    ignore_index=True)
vt.to_csv(G4 / "gate4_verdicts.tsv", sep="\t", index=False)

gate4_summary = {
    "gate": "Gate 4 (synthesis-integrity; constructive, not hypothesis-rejection)",
    "4A_faithful_subtraction": "PASS" if a_ok else "FAIL",
    "4B_primary_secondary_wall": "PASS" if b_ok else "FAIL",
    "4C_quantitative_dominance": "PASS" if c_ok else "FAIL",
    "4D_cross_pillar_arc": "PASS" if d_ok else "FAIL",
    "OVERALL_VERDICT": "PASS" if overall else "FAIL",
    "headline": ("Once EV-specific candidates are removed under a frozen, pre-registered rule, "
                 "the reproducible AKI PT->receiver communication that survives is soluble-classical "
                 "(SPP1) + ECM (collagen-IV/laminin) signaling; 0 of the reproducible survivors are "
                 "EV-compatible above the soluble-secretome background."),
    "naming_constraint": f'secondary layer named "{REQUIRED_NAME}" (never "the surviving program").',
    "language_ceiling": ("association/candidacy only; 'dominant' is scoped to transcriptomically-resolvable "
                         "reproducible PT->receiver axes; does not exclude EV routing at the protein level "
                         "nor claim these are the only communication routes."),
}
json.dump(gate4_summary, open(G4 / "gate4_summary.json", "w"), indent=2)

print(vt.to_string(index=False))
print()
print("GATE 4 OVERALL:", gate4_summary["OVERALL_VERDICT"])
print("wrote:", G4 / "gate4_verdicts.tsv", "|", G4 / "gate4_summary.json")
