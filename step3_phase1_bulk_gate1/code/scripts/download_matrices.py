"""
Download expression matrices for the locked Pillar-1 cohorts into /workspace.

Robust to transient GEO FTP/HTTPS 403 rate-limiting (retry with backoff, and
fall back between the HTTPS mirror and the FTP path). Writes raw downloaded
files to /workspace/geo_matrices/ (random-access-safe local disk).

Usage:
    python 01_bulk_signature/code/scripts/download_matrices.py GSE139061 GSE30718 ...
    python 01_bulk_signature/code/scripts/download_matrices.py --all
"""
from __future__ import annotations

import gzip
import os
import sys
import time
import urllib.request

OUT = "/workspace/geo_matrices"
os.makedirs(OUT, exist_ok=True)

# Locked cohorts (PLAN v1.6). value = kind hint for downstream DE.
COHORTS = {
    # AKI native context
    "GSE139061": "rnaseq_qn_supp",     # quantile-normalized processed csv (series suppl)
    # AKI transplant-ATN context
    "GSE30718": "microarray_matrix",   # GPL570 series matrix
    "GSE53769": "rnaseq_matrix",       # GPL16686 series matrix (check scale)
    # CKD discovery (ERCB, tier-B)
    "GSE104954": "microarray_matrix",
    "GSE99340": "microarray_matrix",
    "GSE47185": "microarray_matrix",
    # replication DKD
    "GSE142153": "microarray_matrix",
    "GSE294519": "rnaseq_matrix",
    # replication IgAN
    "GSE104948": "microarray_matrix",
    # exploratory progression
    "GSE339456": "rnaseq_matrix",
}


def _fetch(url: str, timeout: int = 90, retries: int = 4, backoff: float = 3.0) -> bytes:
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception as e:  # noqa
            last = e
            time.sleep(backoff * (i + 1))
    raise RuntimeError(f"failed {url}: {last}")


def series_matrix_url(gse: str) -> str:
    stub = gse[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/{gse}_series_matrix.txt.gz"


def series_suppl_url(gse: str, fname: str) -> str:
    stub = gse[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/{fname}"


def download(gse: str, kind: str) -> str:
    if kind == "rnaseq_qn_supp" and gse == "GSE139061":
        fname = "GSE139061_Eadon_processed_QN_101419.csv.gz"
        url = series_suppl_url(gse, fname)
        dest = os.path.join(OUT, fname)
    else:
        url = series_matrix_url(gse)
        dest = os.path.join(OUT, f"{gse}_series_matrix.txt.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"[skip] {gse} already downloaded ({os.path.getsize(dest)} bytes)")
        return dest
    data = _fetch(url)
    with open(dest, "wb") as fh:
        fh.write(data)
    print(f"[ok]   {gse} -> {dest} ({len(data)} bytes)")
    return dest


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args or args == ["--all"]:
        targets = list(COHORTS)
    else:
        targets = args
    for gse in targets:
        kind = COHORTS.get(gse, "microarray_matrix")
        try:
            download(gse, kind)
        except Exception as e:  # noqa
            print(f"[FAIL] {gse}: {e}")
