"""Pillar 3 route-attribution overlay (PLAN.md sec 4.4) — SECONDARY, descriptive, NOT gated.

For each ligand in the consensus landscape, classify plausible communication route(s):
  - Soluble-classical: classical signal peptide (proxy: member of NONEV_SECRETED_BACKGROUND)
                       AND CellChatDB "Secreted Signaling"
  - Contact-dependent: CellChatDB annotation == "Cell-Cell Contact"
  - ECM:               CellChatDB annotation == "ECM-Receptor"
  - EV-compatible (SPECIFICITY-GATED, guarded): ligand in EV_TIER2_CARGO
                       AND NOT better explained by classical secretion
                       (proxy: NOT in NONEV_SECRETED_BACKGROUND) -> lacks classical SP,
                       AND topology-consistent (curated EV cargo membership encodes this).

MANDATORY negative-control gate: the 38 frozen soluble cytokines/chemokines MUST NOT score
EV-compatible. If >=1 does, the EV-compatibility annotation is declared UNINFORMATIVE.

A ligand may carry multiple route tags; the goal is to DISCRIMINATE, not to force EV.
Uses ONLY frozen artifacts (EV gene sets v1.3, soluble neg-control freeze) + CellChatDB
annotation carried in the CellChat LR tables. No fresh fetches -> fully reproducible.

Usage: python pillar3_ev_route_overlay.py <consensus_dir> <engine_out_dir> <out_dir>
"""
import os, sys, json
import pandas as pd

EV_JSON = "/mnt/results/ev_aki_ckd_project/00_preregistration/feasibility/ev_genesets_frozen_v1.3.json"
NEGCTRL_JSON = "/mnt/results/ev_aki_ckd_project/03_communication_landscape/data/soluble_negative_controls_frozen.json"
DATASETS = ["a12ccb9b", "dea717d4"]
CONDS = ["AKI", "normal", "CKD"]


def canon_parts(x):
    if pd.isna(x):
        return []
    return [p.strip() for p in str(x).replace("&", "_").split("_") if p.strip()]


def load_frozen():
    ev = json.load(open(EV_JSON))
    cargo = set(ev["EV_TIER2_CARGO"])
    bg = set(ev["NONEV_SECRETED_BACKGROUND"])
    nc = json.load(open(NEGCTRL_JSON))["genes"]
    return cargo, bg, set(nc), ev.get("_provenance", {}).get("plan_version", "v1.3")


def build_cellchatdb_annotation(engine_dir):
    """Map ligand-symbol -> set of CellChatDB annotation categories, harvested from all
    CellChat LR tables (which carry the 'annotation' column per interaction)."""
    lig2ann = {}
    for ds in DATASETS:
        for cond in CONDS:
            f = f"{engine_dir}/cellchat_LR_{ds}__{cond}.tsv"
            if not os.path.exists(f):
                continue
            df = pd.read_csv(f, sep="\t")
            if "annotation" not in df.columns:
                continue
            for lig, ann in zip(df["ligand"], df["annotation"]):
                for part in canon_parts(lig):
                    lig2ann.setdefault(part, set()).add(str(ann))
    return lig2ann


def route_tags_for_ligand(lig, cargo, bg, lig2ann):
    """Return (tags list, ev_compatible bool) for one ligand symbol (may be a complex)."""
    parts = canon_parts(lig)
    anns = set()
    for p in parts:
        anns |= lig2ann.get(p, set())
    tags = []
    # CellChatDB-driven route tags
    if "Cell-Cell Contact" in anns:
        tags.append("Contact-dependent")
    if "ECM-Receptor" in anns:
        tags.append("ECM")
    # Soluble-classical: any subunit in classical-secreted background AND secreted-signaling annotation
    in_bg = any(p in bg for p in parts)
    if in_bg and ("Secreted Signaling" in anns):
        tags.append("Soluble-classical")
    elif in_bg and not anns:
        tags.append("Soluble-classical")  # classical SP even if annotation missing
    # EV-compatible (guarded): in cargo AND lacks classical SP (not in bg)
    in_cargo = any(p in cargo for p in parts)
    ev_compatible = in_cargo and (not in_bg)
    if ev_compatible:
        tags.append("EV-compatible")
    if not tags:
        tags.append("Unclassified")
    return tags, ev_compatible


def negative_control_gate(cargo, bg, negctrl, lig2ann):
    """The 38 soluble controls MUST NOT score EV-compatible."""
    violations = []
    for g in sorted(negctrl):
        _, evc = route_tags_for_ligand(g, cargo, bg, lig2ann)
        if evc:
            violations.append(g)
    informative = (len(violations) == 0)
    return informative, violations


def main(consensus_dir, engine_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    cargo, bg, negctrl, evver = load_frozen()
    lig2ann = build_cellchatdb_annotation(engine_dir)

    # ---- MANDATORY negative-control gate FIRST ----
    informative, violations = negative_control_gate(cargo, bg, negctrl, lig2ann)
    gate = {
        "ev_geneset_version": evver,
        "n_soluble_negative_controls": len(negctrl),
        "n_controls_scored_EV_compatible": len(violations),
        "violating_controls": violations,
        "ev_compatibility_annotation_informative": informative,
        "rule": "EV-compatible iff (ligand in EV_TIER2_CARGO) AND (ligand NOT in NONEV_SECRETED_BACKGROUND [lacks classical signal peptide]).",
        "note": ("EV-compatibility annotation is INFORMATIVE: 0/38 soluble controls scored EV-compatible."
                 if informative else
                 "EV-compatibility annotation is UNINFORMATIVE: >=1 soluble control scored EV-compatible; do NOT interpret EV-compatible tags."),
    }
    with open(f"{out_dir}/ev_route_negcontrol_gate.json", "w") as fh:
        json.dump(gate, fh, indent=2)
    print("NEG-CONTROL GATE:", json.dumps(gate, indent=2))

    # ---- annotate consensus axes (if any) ----
    cons_f = f"{consensus_dir}/consensus_landscape_AKIvsnormal.tsv"
    if not os.path.exists(cons_f):
        print("No consensus file; overlay writes gate only.")
        return
    cons = pd.read_csv(cons_f, sep="\t")
    if cons.empty:
        print("Consensus landscape empty; overlay writes gate + empty annotation.")
        cons["route_tags"] = []
        cons["ev_compatible"] = []
        cons.to_csv(f"{out_dir}/consensus_with_routes.tsv", sep="\t", index=False)
        return
    rt, evc = [], []
    for lig in cons["ligand"]:
        tags, e = route_tags_for_ligand(lig, cargo, bg, lig2ann)
        rt.append(";".join(tags))
        evc.append(e)
    cons["route_tags"] = rt
    cons["ev_compatible"] = evc
    # honesty flag: EV-compatible only meaningful if gate is informative
    cons["ev_compatible_reportable"] = [e and informative for e in evc]
    cons.to_csv(f"{out_dir}/consensus_with_routes.tsv", sep="\t", index=False)

    n_evc = int(cons["ev_compatible"].sum())
    summary = {
        "n_consensus_axes": int(len(cons)),
        "n_EV_compatible_axes": n_evc,
        "n_EV_compatible_axes_reportable": int(cons["ev_compatible_reportable"].sum()),
        "route_tag_counts": pd.Series(";".join(cons["route_tags"]).split(";")).value_counts().to_dict(),
        "ev_annotation_informative": informative,
        "interpretation": (
            f"{n_evc} of {len(cons)} consensus injury-linked PT->receiver axes are EV-compatible above the "
            f"soluble-secretome background." if informative else
            "EV-compatibility uninformative (neg-control gate failed); EV tags suppressed."),
    }
    with open(f"{out_dir}/ev_route_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print("EV ROUTE SUMMARY:", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
