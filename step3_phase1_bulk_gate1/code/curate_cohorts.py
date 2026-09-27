"""
Curate GEO candidate cohorts for the Pillar-1 consensus EV-injury signature.

Resolves per-sample case/control structure (from characteristics AND sample
titles), tags control type, records platform/tissue, and applies the
pre-registered inclusion/exclusion criteria (PLAN_pillar1.md).

Output: a candidate table for USER SIGN-OFF before any DE is run.
"""
from __future__ import annotations

import os
import re
import warnings
from collections import Counter

import GEOparse
import pandas as pd

warnings.filterwarnings("ignore")
CACHE = "/workspace/geo_cache"
os.makedirs(CACHE, exist_ok=True)

# phenotype keyword -> normalized condition
CASE_KW = {
    "aki": "AKI", "acute kidney": "AKI", "atn": "AKI", "ischemia": "AKI",
    "ckd": "CKD", "chronic kidney": "CKD", "fibrosis": "fibrosis",
    "hypertensive neph": "CKD", "dkd": "DKD", "diabetic": "DKD",
    "iga": "IgAN",
}
CONTROL_KW = {
    "healthy": "healthy", "living donor": "living_donor", "donor": "living_donor",
    "nephrectomy": "nephrectomy_margin", "tumor-adjacent": "tumor_adjacent",
    "paraneoplastic": "tumor_adjacent", "adjacent": "tumor_adjacent",
    "normal": "normal_unspecified", "control": "control_unspecified",
}


def _pheno_text(gsm):
    m = gsm.metadata
    parts = (
        m.get("source_name_ch1", [])
        + m.get("characteristics_ch1", [])
        + m.get("title", [])
    )
    return " ".join(parts).lower()


def _classify(text):
    """Return (is_case, condition, is_control, control_type) from one sample's text."""
    cond = None
    for kw, c in CASE_KW.items():
        if kw in text:
            cond = c
            break
    ctrl_type = None
    for kw, t in CONTROL_KW.items():
        if kw in text:
            ctrl_type = t
            break
    # A sample is a control if it has a control keyword AND no disease keyword,
    # OR its text clearly says normal/control without a disease term.
    is_control = ctrl_type is not None and cond in (None,)
    is_case = cond is not None and not is_control
    return is_case, cond, is_control, ctrl_type


def curate(gse: str) -> dict:
    try:
        g = GEOparse.get_GEO(geo=gse, destdir=CACHE, silent=True, include_data=False)
    except Exception as e:  # noqa
        return {"gse": gse, "error": str(e)[:80]}

    plat = ";".join(sorted({p for p in g.gpls})) if g.gpls else "?"
    rows = []
    for name, gsm in g.gsms.items():
        org = " ".join(gsm.metadata.get("organism_ch1", [])).lower()
        if "sapiens" not in org:
            continue
        text = _pheno_text(gsm)
        is_case, cond, is_ctrl, ctype = _classify(text)
        rows.append((cond, is_case, is_ctrl, ctype))

    n_human = len(rows)
    cases = [r for r in rows if r[1]]
    ctrls = [r for r in rows if r[2]]
    unresolved = n_human - len(cases) - len(ctrls)
    cond_counts = Counter(r[0] for r in cases if r[0])
    ctype_counts = Counter(r[3] for r in ctrls if r[3])

    return {
        "gse": gse,
        "platform": plat,
        "n_human": n_human,
        "n_case": len(cases),
        "n_control": len(ctrls),
        "n_unresolved": unresolved,
        "case_conditions": dict(cond_counts),
        "control_types": dict(ctype_counts),
    }


if __name__ == "__main__":
    import json
    import sys

    gses = sys.argv[1:] or [
        "GSE139061", "GSE271143", "GSE339456", "GSE294519",
    ]
    out = [curate(g) for g in gses]
    print(json.dumps(out, indent=2))
