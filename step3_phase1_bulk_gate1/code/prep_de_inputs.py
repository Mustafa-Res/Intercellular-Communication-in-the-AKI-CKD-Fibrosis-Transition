"""
Prepare per-cohort inputs for the limma DE R script.

For each cohort: load expr+meta (load_expression.py), log2-transform linear data,
build a feature->gene_symbol map appropriate to the platform, and write three
TSVs to a per-cohort work dir under /workspace/de_inputs/<GSE>/.

Gene mapping strategy per platform:
  - GSE139061 : index already gene symbols -> identity map
  - GSE30718  (GPL570) : 'Gene Symbol' column from platform annotation (split '///')
  - GSE53769  (GPL16686): mapped via R annotation (hugene20sttranscriptcluster.db),
                          handled in R; here we emit GB_ACC map as fallback.
Others (CKD/DKD/IgAN) added as their platforms are resolved.
"""
from __future__ import annotations

import os
import warnings

import GEOparse
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from load_expression import load_cohort  # noqa: E402

CACHE = "/workspace/geo_cache"
WORK = "/workspace/de_inputs"


def _is_linear(vals: np.ndarray) -> bool:
    """Heuristic: log-scale arrays top out < ~30; linear counts go far higher."""
    return np.nanmax(vals) > 100


def _platform_gene_map(gse: str, gpl: str) -> pd.DataFrame | None:
    g = GEOparse.get_GEO(geo=gse, destdir=CACHE, silent=True, include_data=False)
    plat = g.gpls.get(gpl)
    if plat is None:
        return None
    tbl = plat.table
    if "Gene Symbol" in tbl.columns:  # GPL570 style
        m = tbl[["ID", "Gene Symbol"]].copy()
        m.columns = ["feature", "gene_symbol"]
        # take first symbol when multi-mapped ('DDR1 /// MIR4640')
        m["gene_symbol"] = (m["gene_symbol"].astype(str)
                            .str.split("///").str[0].str.strip())
        m = m[(m["gene_symbol"] != "") & (m["gene_symbol"].str.lower() != "nan")]
        return m
    return None  # needs R annotation (e.g. GPL16686)


PLATFORMS = {
    "GSE30718": "GPL570",
    "GSE53769": "GPL16686",
}


def prep(gse: str) -> dict:
    expr, meta = load_cohort(gse)
    vals = expr.to_numpy(dtype=float)
    logged = False
    if _is_linear(vals):
        # log2(x+1) on linear abundances
        expr = np.log2(expr.astype(float) + 1.0)
        logged = True
    out_dir = os.path.join(WORK, gse)
    os.makedirs(out_dir, exist_ok=True)

    expr_path = os.path.join(out_dir, "expr.tsv")
    meta_path = os.path.join(out_dir, "meta.tsv")
    map_path = os.path.join(out_dir, "genemap.tsv")

    expr.to_csv(expr_path, sep="\t")
    meta[["group", "control_type", "context"]].to_csv(meta_path, sep="\t")

    # gene map
    gpl = PLATFORMS.get(gse)
    map_written = None
    if gse == "GSE139061":
        # index already symbols
        gm = pd.DataFrame({"feature": expr.index, "gene_symbol": expr.index})
        gm.to_csv(map_path, sep="\t", index=False)
        map_written = map_path
    elif gpl:
        gm = _platform_gene_map(gse, gpl)
        if gm is not None:
            gm.to_csv(map_path, sep="\t", index=False)
            map_written = map_path
        else:
            map_written = None  # R will annotate (GPL16686)

    return {
        "gse": gse, "logged": logged, "n_feat": expr.shape[0],
        "n_samp": expr.shape[1], "groups": meta["group"].value_counts().to_dict(),
        "context": meta["context"].iloc[0],
        "expr": expr_path, "meta": meta_path,
        "map": map_written, "gpl": gpl,
    }


if __name__ == "__main__":
    import json
    import sys
    for gse in sys.argv[1:] or ["GSE139061", "GSE30718", "GSE53769"]:
        print(json.dumps(prep(gse), indent=2, default=str))
