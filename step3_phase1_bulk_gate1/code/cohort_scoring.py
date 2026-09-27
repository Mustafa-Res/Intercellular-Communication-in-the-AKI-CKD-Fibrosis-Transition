"""
Cohort confidence scoring (A/B/C) for Pillar-1 (user-required, PLAN v1.4).

Transparent 5-criterion score -> tier that GOVERNS which cohorts touch the
primary result. Prevents weak cohorts from silently influencing the headline.

Criteria (each +1):
  1. kidney_tissue_confirmed  : genuine kidney tissue (not blood/urine/organoid/cell line)
  2. phenotype_clear          : unambiguous disease label at sample level
  3. control_appropriate      : has an appropriate control group (healthy/donor/nephrectomy/protocol)
  4. platform_quality         : RNA-seq or a well-established array platform, adequate n
  5. clinical_metadata        : disease + at least basic covariate metadata present

Tier:
  A (primary discovery)  : score == 5  AND control_appropriate AND phenotype_clear
  B (sensitivity)        : score 3-4, or score 5 with a tagged caveat (e.g. transplant/FFPE/tumor-adjacent)
  C (exploratory)        : score <= 2, or non-standard contrast (progression-only, no healthy control)
"""
from __future__ import annotations

import pandas as pd


def score_row(r: dict) -> dict:
    s = {
        "kidney_tissue_confirmed": int(r.get("kidney_tissue", False)),
        "phenotype_clear": int(r.get("phenotype_clear", False)),
        "control_appropriate": int(r.get("control_appropriate", False)),
        "platform_quality": int(r.get("platform_quality", False)),
        "clinical_metadata": int(r.get("clinical_metadata", False)),
    }
    total = sum(s.values())

    caveat = r.get("caveat", "")
    if total == 5 and s["control_appropriate"] and s["phenotype_clear"] and not caveat:
        tier = "A"
    elif r.get("progression_only") or (not s["control_appropriate"]) or total <= 2:
        tier = "C"
    else:
        tier = "B"

    return {**s, "score": total, "tier": tier, "caveat": caveat}


# Finalized, hand-curated cohort facts (verified via sample-level GEO metadata this session)
COHORTS = [
    # AKI discovery
    dict(GSE="GSE139061", phenotype="AKI", role="discovery", platform="RNA-seq (GPL20301)",
         n_case=39, n_control=9, control_type="nephrectomy_margin",
         kidney_tissue=True, phenotype_clear=True, control_appropriate=True,
         platform_quality=True, clinical_metadata=True, caveat="",
         notes="Native-kidney AKI/ATN; canonical human AKI biopsy cohort"),
    dict(GSE="GSE30718", phenotype="AKI", role="discovery", platform="microarray (GPL570)",
         n_case=28, n_control=19, control_type="nephrectomy+pristine_protocol",
         kidney_tissue=True, phenotype_clear=True, control_appropriate=True,
         platform_quality=True, clinical_metadata=True, caveat="",
         notes="Kidney-allograft ATN vs nephrectomy + pristine protocol biopsy (Famulski/Halloran)"),
    dict(GSE="GSE53769", phenotype="AKI", role="discovery", platform="microarray (GPL16686, Affymetrix HuGene-2.0-ST)",
         n_case=8, n_control=10, control_type="protocol_biopsy",
         kidney_tissue=True, phenotype_clear=True, control_appropriate=True,
         platform_quality=True, clinical_metadata=True, caveat="transplant ATN context; small n",
         notes="Transplant ATN vs protocol biopsy w/o pathology"),
    # CKD discovery (ERCB tubulointerstitial + glomerular)
    dict(GSE="GSE104954", phenotype="CKD", role="discovery", platform="microarray (ERCB TI)",
         n_case=62, n_control=5, control_type="nephrectomy_margin",
         kidney_tissue=True, phenotype_clear=True, control_appropriate=True,
         platform_quality=True, clinical_metadata=True, caveat="few controls; tag by diagnosis+compartment",
         notes="ERCB (European Renal cDNA Bank) tubulointerstitial multi-disease (CKD/DKD/IgAN)"),
    dict(GSE="GSE99340", phenotype="CKD", role="discovery", platform="microarray (ERCB)",
         n_case=120, n_control=21, control_type="living_donor+nephrectomy",
         kidney_tissue=True, phenotype_clear=True, control_appropriate=True,
         platform_quality=True, clinical_metadata=True, caveat="tag compartment (glom+TI)",
         notes="Large ERCB multi-disease; living-donor controls"),
    dict(GSE="GSE47185", phenotype="CKD", role="discovery", platform="microarray (ERCB)",
         n_case=34, n_control=27, control_type="living_donor+nephrectomy",
         kidney_tissue=True, phenotype_clear=True, control_appropriate=True,
         platform_quality=True, clinical_metadata=True, caveat="tag compartment",
         notes="ERCB multi-disease microarray"),
    # Replication - DKD
    dict(GSE="GSE142153", phenotype="DKD", role="replication", platform="microarray (GPL6480, Agilent)",
         n_case=23, n_control=10, control_type="healthy",
         kidney_tissue=False, phenotype_clear=True, control_appropriate=True,
         platform_quality=True, clinical_metadata=True,
         caveat="PBMC (blood), not kidney tissue",
         notes="Blood-compartment DKD study; ineligible for the kidney-tissue replication role; never executed"),
    dict(GSE="GSE294519", phenotype="DKD", role="replication", platform="RNA-seq",
         n_case=23, n_control=13, control_type="tumor_adjacent",
         kidney_tissue=True, phenotype_clear=True, control_appropriate=True,
         platform_quality=True, clinical_metadata=True, caveat="tumor-adjacent controls (tag)",
         notes="DKD renal tubules vs paraneoplastic normal"),
    # Replication - IgAN (ERCB glomerular)
    dict(GSE="GSE104948", phenotype="IgAN", role="replication", platform="microarray (ERCB glom)",
         n_case=27, n_control=5, control_type="tumor_nephrectomy",
         kidney_tissue=True, phenotype_clear=True, control_appropriate=True,
         platform_quality=True, clinical_metadata=True, caveat="tag diagnosis; glomerular compartment",
         notes="ERCB glomerular; IgAN=27, also DKD/FSGS/lupus etc."),
    # Progression-specific (exploratory)
    dict(GSE="GSE339456", phenotype="CKD-progression", role="progression", platform="RNA-seq (FFPE)",
         n_case=18, n_control=17, control_type="stable_vs_progressive",
         kidney_tissue=True, phenotype_clear=True, control_appropriate=False,
         platform_quality=True, clinical_metadata=True, caveat="no healthy control; FFPE",
         progression_only=True,
         notes="Stable vs progressive CKD (progression contrast, not case/control)"),
]


def build_table() -> pd.DataFrame:
    rows = []
    for c in COHORTS:
        sc = score_row(c)
        rows.append({
            "GSE": c["GSE"], "phenotype": c["phenotype"], "role": c["role"],
            "platform": c["platform"], "n_case": c["n_case"], "n_control": c["n_control"],
            "control_type": c["control_type"],
            "score": sc["score"], "tier": sc["tier"], "caveat": sc["caveat"],
            "notes": c["notes"],
        })
    df = pd.DataFrame(rows).sort_values(["tier", "phenotype", "GSE"]).reset_index(drop=True)
    return df


if __name__ == "__main__":
    df = build_table()
    with pd.option_context("display.max_colwidth", 60, "display.width", 200):
        print(df.to_string(index=False))
    print("\nTier counts:", dict(df["tier"].value_counts()))
    print("Tier-A by phenotype:", dict(df[df.tier == "A"].groupby("phenotype").size()))
