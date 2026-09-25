  Tier-1        : curated literature-standard (van Niel 2022 NRMCB; Colombo 2014)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# Self-contained: read frozen gene-set source JSONs from the repo's data_manifests.
# 00_preregistration/ = parents[2] (00_preregistration/code/ev_genesets/ev_genesets.py)
RAW = Path(__file__).resolve().parents[2] / "genesets" / "geneset_sources"


# ---------------------------------------------------------------- loaders
def _load_go_sets() -> dict:
    return json.loads((RAW / "go_ev_sets.json").read_text())


def _load_tier1() -> dict:
    return json.loads((RAW / "tier1_machinery.json").read_text())


def _load_secreted() -> set:
    return set(json.loads((RAW / "uniprot_secreted.json").read_text()))


def _load_negative_controls() -> dict:
    return json.loads((RAW / "negative_controls.json").read_text())


def build_negative_controls() -> dict:
    """Negative-control gene sets for the EV-specificity test (PLAN v1.3).

    C1 general secreted  -> reuse NONEV_SECRETED_BACKGROUND (built in build_ev_genesets)
    C2 mitochondrial     -> expect NOT specifically enriched
    C3 ribosomal         -> expect NOT enriched
    C4 lysosomal/endosomal -> biologically adjacent (shared trafficking); key discriminator
    """
    nc = _load_negative_controls()
    return {
        "CTRL_MITOCHONDRIAL": set(nc["MITOCHONDRIAL"]),
        "CTRL_RIBOSOMAL": set(nc["RIBOSOMAL"]),
        "CTRL_LYSO_ENDOSOMAL": set(nc["LYSO_ENDOSOMAL"]),
    }


# ---------------------------------------------------------------- builder
def build_ev_genesets() -> dict:
    """Return the frozen hierarchical EV gene-set dictionary (sets of HGNC symbols)."""
    go = _load_go_sets()
    tier1_cat = _load_tier1()
    secreted = _load_secreted()

    tier1 = {g for v in tier1_cat.values() for g in v}  # strict machinery

    tier2_cargo = set(go["GO:0070062"]["genes"])  # broad GO exosome cargo
    mvb_pathway = (
        set(go["GO:0005771"]["genes"])
        | set(go["GO:0036258"]["genes"])
        | set(go["GO:0071985"]["genes"])
    )
    # Tier-1 primary set is machinery + the mechanistic MVB pathway genes
    ev_primary = tier1 | mvb_pathway

    # Background = secreted minus the STRICT primary core only (not broad cargo)
    nonev_background = secreted - ev_primary

    return {
        "EV_TIER1_PRIMARY": ev_primary,          # headline specificity test
        "EV_TIER1_MACHINERY": tier1,             # machinery-only (subset, for reporting)
        "EV_MVB_PATHWAY": mvb_pathway,
        "EV_TIER2_CARGO": tier2_cargo,           # broad GO, sensitivity only
        "NONEV_SECRETED_BACKGROUND": nonev_background,
        "_tier1_categories": tier1_cat,          # provenance
    }


def summary() -> dict:
    s = build_ev_genesets()
    return {k: len(v) for k, v in s.items() if not k.startswith("_")}


# ---------------------------------------------- covariate-matched null sampler
def match_background(
    ev_genes,
    background_genes,
    covariates: pd.DataFrame,
    n_per_ev: int = 20,
    cov_cols=("mean_expr", "gene_length", "detection_freq"),
    n_bins: int = 10,
    random_state: int = 0,
):
    """Draw an expression/length/detection-matched background set.

    For each EV gene present in `covariates`, sample `n_per_ev` background genes
    from the same multivariate covariate stratum (quantile bins per covariate).
    Returns a de-duplicated list of matched background gene symbols.

    covariates : DataFrame indexed by gene symbol with columns in `cov_cols`
                 (use whichever covariates are available; missing cols are skipped).
    """
    rng = np.random.default_rng(random_state)
    cols = [c for c in cov_cols if c in covariates.columns]
    if not cols:
        raise ValueError("None of the requested covariate columns are present.")

    cov = covariates.dropna(subset=cols).copy()
    # Build quantile bin labels per covariate
    for c in cols:
        cov[c + "_bin"] = pd.qcut(cov[c].rank(method="first"), n_bins, labels=False)
    bin_cols = [c + "_bin" for c in cols]
    cov["stratum"] = cov[bin_cols].astype(str).agg("|".join, axis=1)

    ev_in = [g for g in ev_genes if g in cov.index]
    bg_in = [g for g in background_genes if g in cov.index]
    bg_set = set(bg_in)

    strata_to_bg = (
        cov.loc[cov.index.isin(bg_set)]
        .groupby("stratum")
        .apply(lambda d: list(d.index))
        .to_dict()
    )

    matched = []
    for g in ev_in:
        stratum = cov.at[g, "stratum"]
        pool = [x for x in strata_to_bg.get(stratum, []) if x != g]
        if not pool:
            continue
        take = min(n_per_ev, len(pool))
        matched.extend(rng.choice(pool, size=take, replace=False).tolist())

    return sorted(set(matched)), {
        "n_ev_matched": len(ev_in),
        "n_background_pool": len(bg_in),
        "n_matched_drawn": len(set(matched)),
        "covariates_used": cols,
    }


if __name__ == "__main__":
    from pprint import pprint

    sets = build_ev_genesets()
    pprint(summary())

    # --- sanity checks (fail loudly) ---
    assert not (sets["NONEV_SECRETED_BACKGROUND"] & sets["EV_TIER1_PRIMARY"]), \
        "Background overlaps Tier-1 primary!"
    for m in ["CD9", "CD63", "CD81", "TSG101", "PDCD6IP", "RAB27A", "CHMP4B"]:
        assert m in sets["EV_TIER1_MACHINERY"], f"{m} missing from machinery"
    # Conventional-secretion proteins MUST stay in the null under the new design
    assert "ALB" in sets["NONEV_SECRETED_BACKGROUND"], "ALB should be in null (v1.2 design)"
    assert "INS" in sets["NONEV_SECRETED_BACKGROUND"], "INS should be in null"
    # Tier-2 broad cargo should be a superset-ish of machinery cargo overlap
    overlap = len(sets["EV_TIER1_MACHINERY"] & sets["EV_TIER2_CARGO"])
    print(f"Tier1 machinery genes also in broad GO cargo: {overlap}/{len(sets['EV_TIER1_MACHINERY'])}")
    print("ALB in null:", "ALB" in sets["NONEV_SECRETED_BACKGROUND"])
    print("All assertions passed.")
