"""Freeze the soluble-cytokine/chemokine NEGATIVE-CONTROL list for the Pillar-3
EV-route specificity gate (PLAN.md sec 4.4). These are canonical CLASSICAL
signal-peptide-bearing soluble mediators. The EV-compatibility annotation MUST NOT
flag any of them EV-compatible; if it does, the annotation is declared uninformative.

Defined BEFORE any inference. Written with a SHA256 hash to durable storage.
The list is curated (classical secreted cytokines/chemokines with N-terminal signal
peptides, canonical soluble action) and then validated against the frozen gene sets:
  - each SHOULD be in NONEV_SECRETED_BACKGROUND (genuine soluble secreted), and
  - we RECORD whether each also appears in EV_TIER2_CARGO (the Vesiclepedia
    false-positive surface we are guarding against).
"""
import json, hashlib, datetime

FROZEN = "/mnt/results/ev_aki_ckd_project/00_preregistration/feasibility/ev_genesets_frozen_v1.3.json"
OUT = "/mnt/results/ev_aki_ckd_project/03_communication_landscape/data/soluble_negative_controls_frozen.json"

# --- curated canonical soluble cytokines/chemokines (classical signal peptide) ---
# Chemokines (CC + CXC families): archetypal soluble, signal-peptide-secreted.
CHEMOKINES = [
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8", "CCL19", "CCL20", "CCL21",
    "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL9", "CXCL10",
    "CXCL11", "CXCL12", "CXCL13", "CXCL16",
]
# Classical pro-inflammatory / immunomodulatory interleukins + TNF/IFN (soluble forms).
CYTOKINES = [
    "IL1B", "IL6", "IL10", "IL15", "IL18", "IL34",
    "TNF", "IFNG", "CSF1", "CSF2", "CSF3",
    "TGFB1", "TGFB2",  # secreted (latent complex) — classical signal peptide
]
# Growth factors with unambiguous classical secretion (soluble action controls).
GROWTH_FACTORS = ["HGF", "IGF1", "VEGFA", "PDGFB"]

NEG_CONTROLS = sorted(set(CHEMOKINES + CYTOKINES + GROWTH_FACTORS))

with open(FROZEN) as f:
    gs = json.load(f)
nonev = set(gs["NONEV_SECRETED_BACKGROUND"])
ev_cargo = set(gs["EV_TIER2_CARGO"])

records = []
for g in NEG_CONTROLS:
    records.append({
        "gene": g,
        "in_NONEV_SECRETED_BACKGROUND": g in nonev,   # expect True (genuine soluble)
        "in_EV_TIER2_CARGO": g in ev_cargo,           # the false-positive surface
    })

n_in_nonev = sum(r["in_NONEV_SECRETED_BACKGROUND"] for r in records)
n_in_evcargo = sum(r["in_EV_TIER2_CARGO"] for r in records)

payload = {
    "description": ("Soluble cytokine/chemokine NEGATIVE CONTROLS for the Pillar-3 "
                    "EV-route specificity gate. These classical signal-peptide soluble "
                    "mediators MUST NOT be flagged EV-compatible; if any is, the "
                    "EV-compatibility annotation is declared UNINFORMATIVE."),
    "frozen_date": datetime.date.today().isoformat(),
    "n_controls": len(NEG_CONTROLS),
    "genes": NEG_CONTROLS,
    "validation": {
        "n_in_NONEV_SECRETED_BACKGROUND": n_in_nonev,
        "n_of_total": len(NEG_CONTROLS),
        "n_also_in_EV_TIER2_CARGO": n_in_evcargo,
        "note": ("Controls also appearing in EV_TIER2_CARGO illustrate the "
                 "Vesiclepedia false-positive surface: presence in an EV catalog "
                 "must NOT by itself yield an EV-compatible call — the specificity "
                 "gate (signal-peptide + topology + over-secretome) must exclude them."),
    },
    "records": records,
    "source_genesets": {"file": FROZEN,
                        "sha256_16": hashlib.sha256(open(FROZEN, "rb").read()).hexdigest()[:16]},
}

payload_str = json.dumps(payload, indent=2, sort_keys=True)
payload["self_sha256_16"] = hashlib.sha256(payload_str.encode()).hexdigest()[:16]

import os
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    json.dump(payload, f, indent=2)

print(f"Froze {len(NEG_CONTROLS)} soluble negative controls -> {OUT}")
print(f"  in NONEV_SECRETED_BACKGROUND: {n_in_nonev}/{len(NEG_CONTROLS)}")
print(f"  ALSO in EV_TIER2_CARGO (false-positive surface): {n_in_evcargo}/{len(NEG_CONTROLS)}")
print(f"  self sha256[:16] = {payload['self_sha256_16']}")
# show the ones NOT in nonev background (curation check) and the ones in ev cargo
not_nonev = [r["gene"] for r in records if not r["in_NONEV_SECRETED_BACKGROUND"]]
in_cargo = [r["gene"] for r in records if r["in_EV_TIER2_CARGO"]]
print(f"  NOT in NONEV background (check curation): {not_nonev}")
print(f"  controls in EV cargo catalog: {in_cargo}")
