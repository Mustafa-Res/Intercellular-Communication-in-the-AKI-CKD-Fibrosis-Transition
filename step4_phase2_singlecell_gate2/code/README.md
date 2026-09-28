# Pillar 2 — code

.

## Data source

* **CZ CELLxGENE Census**, version `2025-01-30`, accessed via the
`cellxgene-census` Python API (`open\_soma` / `get\_anndata`).
* Cell filter: `tissue\_general == 'kidney'`, `is\_primary\_data == True`,
10x 3' v3 assay, disease in {AKI, CKD, normal}, iterated per `dataset\_id`.
* Donor->dataset map: `v3\_donor\_dataset\_map.csv` (in the intermediate dir).

## Run order

1. `pillar2\_census\_probe.py` — metadata-only probe: enumerate 10x 3' v3 kidney
datasets + cell counts to plan the pull. (No matrices pulled.)
2. `pillar2\_pilot\_pull.py` — pull ONE AKI dataset; validate compartment mapping
and pseudobulk logic; report peak RSS before committing to the full pull.
3. `pillar2\_full\_pull.py` — pull all 10x 3' v3 kidney datasets, map ontology ->
compartment, build donor x compartment pseudobulk (raw count sums),
checkpoint to the intermediate dir.
4. `pillar2\_split\_compartments.py` — split the full pseudobulk into per-compartment
genes x samples count matrices (+ coldata) for the R DE step.
5. **Differential expression (pseudobulk):**

   * `pillar2\_de.R` — PRIMARY limma-voom + `duplicateCorrelation` (block = donor),
with an edgeR-QLF donor-collapsed sensitivity arm. 3 contrasts
(AKI-vs-normal primary; CKD-vs-normal; AKI-vs-CKD).
   * `pillar2\_edger\_de.R` — edgeR-QLF pseudobulk DE (independent engine).
6. `pillar2\_enrichment.py` — sender EV-specificity test using the **frozen
Pillar-1 enrichment engine** (`run\_gsea\_prerank`) on the PT logFC ranking
(EV Tier-1 vs mito / ribo / secreted-background controls).
7. **Robustness (Step 5):**

   * `pillar2\_robust\_de.R` — per-dataset within-atlas and leave-one-dataset-out
(LODO) DE rankings for PT AKI-vs-normal.
   * `pillar2\_downsample\_fix.R` — down-sampled-control robustness arm.
   * `pillar2\_robust\_gsea.py` — robustness GSEA over those rankings (+ receiver
descriptive GSEA, Step 7).
8. `pillar2\_composition\_severity.py` — compartment abundance shift, severity-
within-PT association, and per-donor pseudobulk EV module score
(supporting / non-gating).
9. `pillar2\_gate2\_eval.py` — evaluate the pre-registered Gate 2 criteria (2A-2D)
and write `02\_singlecell\_specificity/data/gate2\_verdicts.tsv`.

## Paths / regenerability

* `OUT = 02\_singlecell\_specificity/data/` — **durable** frozen outputs (gate verdicts, DE
summaries, enrichment, robustness tables). Committed in this repo.
* `DURABLE = /mnt/shared-workspace/shared/pillar2` — a working volume holding the
**regenerable intermediates**: `pseudobulk\_counts.parquet` (\~82 MB),
`pseudobulk\_meta.tsv`, `gene\_map.tsv`, and `de\_inputs/` / `de\_out/` /
`robust\_out/`. These are **not** committed; they are fully regenerable by
re-running steps 3-5 against Census `2025-01-30`. Repoint `DURABLE` to any
writable scratch location if regenerating.
* The frozen EV gene sets are read from `00\_preregistration/feasibility/ev\_genesets\_frozen\_v1.3.json`.

## Reuse of frozen Pillar-1 machinery

Steps 6-7 import the frozen Pillar-1 enrichment engine
(`sys.path` -> `01\_bulk\_signature/code`, `run\_gsea\_prerank`) and `00\_preregistration/code/plotting\_style.py`,
so the EV-specificity test is identical across pillars.

