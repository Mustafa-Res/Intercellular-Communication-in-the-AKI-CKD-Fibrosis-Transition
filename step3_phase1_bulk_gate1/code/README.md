# Pillar 1 — code

Python/R pipeline for Pillar 1. this directory documents *how*
the frozen results in `01\_bulk\_signature/data/` were produced.

## Data source

* Ten public human kidney bulk-transcriptome cohorts from NCBI GEO
(accessions, platforms, case/control structure, and tiers:
`00\_preregistration/feasibility/pillar1\_cohorts\_selected.csv`).

## Run order

1. `scripts/download\_matrices.py` — download GEO series/platform matrices
(retry/backoff against GEO rate limiting).
2. `curate\_cohorts.py` — resolve per-sample case/control structure, tag control
type, apply the pre-registered inclusion/exclusion criteria.
3. `cohort\_scoring.py` — transparent 5-criterion A/B/C cohort-confidence score
governing which cohorts touch the primary result.
4. `load\_expression.py` — library: aligned genes x samples expression +
sample metadata per cohort (series-matrix and processed-CSV handling).
5. `prep\_de\_inputs.py` — log2-transform, build feature→gene maps, write
per-cohort DE input TSVs.
6. `scripts/map\_gpl16686.R` — GPL16686 probe→gene-symbol map (org.Hs.eg.db).
7. `scripts/run\_de\_limma.R` — per-cohort limma differential expression.
8. `scripts/run\_meta\_transplant.R` — Analysis B: transplant-ATN random-effects
meta-analysis (DerSimonian–Laird primary, REML sensitivity).
9. `enrichment.py` — library: EV-specificity enrichment engine
(`run\_gsea\_prerank`, negative-control ladder machinery) shared with Pillar 2.
10. `scripts/run\_permutation\_ladder.py` — covariate-matched permutation test
(B = 10,000) for the EV set + negative-control ladder, both contexts.

## Notes

* Absolute paths recorded inside scripts reflect the original analysis
environment; adapt them before re-running.

