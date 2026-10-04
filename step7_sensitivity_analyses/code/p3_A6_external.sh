#!/bin/bash
# A6: external replication in GSE210622 (non-COVID AKI vs control snRNA).
# Same prep (cap 2000/compartment, seed 1), same two engines, same consensus thresholds.
set -e
PY=/workspace/.venv/bin/python3
R=/opt/conda/bin/Rscript
P3=/mnt/shared-workspace/shared/pillar3
OUT=/workspace/p3_gse210622
mkdir -p "$OUT"
for tag in gse210622__AKI gse210622__normal; do
  if [ ! -f "$OUT/liana_${tag}.tsv" ]; then
    echo "=== A6 prep+LIANA $tag ==="
    env -i PATH=/opt/conda/bin:/usr/bin:/bin HOME=/root "$PY" "$P3/pillar3_prep_and_liana.py" "$P3/${tag}.h5ad" "$OUT" --cap 2000 --seed 1 --gene_map "$P3/gene_map_gse210622.tsv"
  fi
  if [ ! -f "$OUT/cellchat_LR_${tag}.tsv" ]; then
    echo "=== A6 CellChat $tag ==="
    "$R" "$P3/pillar3_run_cellchat.R" "$OUT/forR_${tag}" "$OUT" "${tag}"
  fi
done
echo "=== A6 consensus (single-dataset) ==="
env -i PATH=/opt/conda/bin:/usr/bin:/bin HOME=/root "$PY" /workspace/p3_consensus_gse.py "$OUT" "$OUT/consensus"
echo A6_DONE
