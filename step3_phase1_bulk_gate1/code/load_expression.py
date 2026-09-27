"""
Load expression matrices + resolve per-sample case/control labels for Pillar-1.

Produces, per cohort, an aligned (genes x samples) expression DataFrame and a
sample-metadata table with columns: sample, group ('case'/'control'), context,
control_type. Handles GEO series-matrix files and the GSE139061 processed CSV.

Gene-ID harmonization to Ensembl is done downstream (map_to_ensembl.py); this
module keeps the platform-native probe/gene IDs.
"""
from __future__ import annotations

import gzip
import io
import os
import re

import GEOparse
import numpy as np
import pandas as pd

MATRIX_DIR = "/workspace/geo_matrices"
CACHE = "/workspace/geo_cache"


# ----------------------------- label resolution -----------------------------
def resolve_labels(gse: str) -> pd.DataFrame:
    """Return per-sample DataFrame: sample, group, context, control_type, raw_diag."""
    g = GEOparse.get_GEO(geo=gse, destdir=CACHE, silent=True, include_data=False)
    rows = []
    for name, gsm in g.gsms.items():
        m = gsm.metadata
        org = " ".join(m.get("organism_ch1", [])).lower()
        if "sapiens" not in org:
            continue
        src = " ".join(m.get("source_name_ch1", [])).lower()
        char = " ".join(m.get("characteristics_ch1", [])).lower()
        title = " ".join(m.get("title", [])).lower()
        blob = f"{src} | {char} | {title}"
        rows.append({"sample": name, "source": src, "char": char,
                     "title": title, "blob": blob})
    df = pd.DataFrame(rows)
    return df


# Cohort-specific label rules (explicit, auditable). Each returns group + control_type.
def _label_gse139061(df):
    # native kidney: 'disease: AKI' vs control (Reference/nephrectomy/normal)
    grp, ct = [], []
    for b in df["char"]:
        if "aki" in b or "acute kidney" in b:
            grp.append("case"); ct.append(None)
        else:
            grp.append("control"); ct.append("nephrectomy_margin")
    df["group"] = grp; df["control_type"] = ct; df["context"] = "native"
    return df


def _label_gse30718(df):
    grp, ct = [], []
    for b in df["char"]:
        if "acute kidney injury" in b or "aki" in b:
            grp.append("case"); ct.append(None)
        elif "nephrectomy" in b:
            grp.append("control"); ct.append("nephrectomy")
        elif "protocol" in b or "pristine" in b:
            grp.append("control"); ct.append("protocol_biopsy")
        else:
            grp.append("other"); ct.append(None)
    df["group"] = grp; df["control_type"] = ct; df["context"] = "transplant_ATN"
    return df


def _label_gse53769(df):
    # 3 biopsy types: ATN (case), post-TX protocol biopsy w/o pathology (control),
    # zero-hour/pre-implantation biopsy (EXCLUDE - distinct baseline state, would confound).
    grp, ct = [], []
    for b in df["blob"]:
        if "atn" in b or "acute tubular" in b:
            grp.append("case"); ct.append(None)
        elif "zero-hour" in b or "zero hour" in b or "pre-impl" in b:
            grp.append("exclude_zerohour"); ct.append(None)
        elif "protocol" in b and "without pathology" in b:
            grp.append("control"); ct.append("protocol_biopsy")
        elif "protocol" in b:
            grp.append("control"); ct.append("protocol_biopsy")
        else:
            grp.append("other"); ct.append(None)
    df["group"] = grp; df["control_type"] = ct; df["context"] = "transplant_ATN"
    return df


LABELERS = {
    "GSE139061": _label_gse139061,
    "GSE30718": _label_gse30718,
    "GSE53769": _label_gse53769,
}


# ----------------------------- matrix loading -----------------------------
def _read_series_matrix(path: str) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as fh:
        txt = fh.read()
    lines = txt.splitlines()
    start = end = None
    for i, l in enumerate(lines):
        if l.startswith("!series_matrix_table_begin"):
            start = i + 1
        elif l.startswith("!series_matrix_table_end"):
            end = i
            break
    if start is None:
        raise ValueError("no data table in series matrix")
    table = "\n".join(lines[start:end])
    df = pd.read_csv(io.StringIO(table), sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip('"')
    df.columns = [c.strip('"') for c in df.columns]
    return df


def load_matrix(gse: str) -> pd.DataFrame:
    """Return genes/probes x samples expression matrix (platform-native IDs)."""
    if gse == "GSE139061":
        path = os.path.join(MATRIX_DIR, "GSE139061_Eadon_processed_QN_101419.csv.gz")
        df = pd.read_csv(path, index_col=0)
        df.index = df.index.astype(str)
        return df
    path = os.path.join(MATRIX_DIR, f"{gse}_series_matrix.txt.gz")
    return _read_series_matrix(path)


def _cohort_gse139061():
    """GSE139061 processed CSV encodes group in column names: AKI* vs REF*.
    Values are linear RNA-seq abundances (quantile-normalized). Returns (expr, meta).
    """
    expr = load_matrix("GSE139061")
    cols = list(expr.columns)
    grp, ct = [], []
    keep = []
    for c in cols:
        cu = c.upper()
        if cu.startswith("AKI"):
            keep.append(c); grp.append("case"); ct.append(None)
        elif cu.startswith("REF") or cu.startswith("CONTROL") or cu.startswith("NORM"):
            keep.append(c); grp.append("control"); ct.append("nephrectomy_margin")
    expr = expr[keep]
    meta = pd.DataFrame({"sample": keep, "group": grp, "control_type": ct,
                         "context": "native"}).set_index("sample")
    return expr, meta


def load_cohort(gse: str):
    """Return (expr, meta) aligned on samples present in BOTH."""
    if gse == "GSE139061":
        return _cohort_gse139061()
    expr = load_matrix(gse)
    meta = resolve_labels(gse)
    meta = LABELERS[gse](meta)
    # align: keep only case/control samples present in expr columns
    meta = meta[meta["group"].isin(["case", "control"])].copy()
    common = [s for s in meta["sample"] if s in expr.columns]
    meta = meta[meta["sample"].isin(common)].set_index("sample").loc[common]
    expr = expr[common]
    return expr, meta


if __name__ == "__main__":
    import sys
    for gse in sys.argv[1:] or ["GSE139061", "GSE30718", "GSE53769"]:
        expr, meta = load_cohort(gse)
        vals = expr.to_numpy(dtype=float)
        print(f"=== {gse} ===")
        print(f"  expr shape: {expr.shape}  (genes/probes x samples)")
        print(f"  groups: {meta['group'].value_counts().to_dict()}")
        print(f"  context: {meta['context'].iloc[0]}")
        print(f"  value range: min={np.nanmin(vals):.2f} max={np.nanmax(vals):.2f} "
              f"median={np.nanmedian(vals):.2f}  (log-scale if max<~30)")
        print(f"  NaN fraction: {np.isnan(vals).mean():.3f}")
