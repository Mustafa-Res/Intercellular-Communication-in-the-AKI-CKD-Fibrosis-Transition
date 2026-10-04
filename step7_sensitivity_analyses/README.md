# Step 7 — Sensitivity analyses

Cross-cutting robustness analyses that stress-test the headline results. The two
analyses whose code and outputs live **here** are the ones that span or sit
outside the four phases:

- **Shared-donor exclusion** (`code/p3_A1_shared_donor.py`): the two kidney
  atlases share 10 AKI donors (same Census collection). Each atlas's AKI arm is
  reduced in turn (6 snRNA / 4 scRNA AKI donors remain), both engines are
  re-run, and the 3 consensus axes are re-evaluated against the fixed normal
  side under the full per-dataset pass criteria. Result:
  `data/TS77_shared_donor_exclusion.tsv` — only SPP1→ITGAV_ITGB1 survives both
  reduced arms.
- **External replication in GSE210622** (`code/p3_A6_external.sh`,
  `code/build_gse210622_shards.py`, `code/p3_consensus_gse.py`): an independent
  non-COVID AKI snRNA cohort, prepared identically (cap 2000 cells per
  compartment, seed 1), scored with the same two engines and the same consensus
  thresholds. Result: `data/TS78_external_replication_gse210622.tsv` — the SPP1
  axes replicate at direction level with CellChat significance; COL4A1 is
  discordant; no axis passes the full criteria externally.

`data/shared_donors.json` — the 10 shared donor IDs used by the exclusion
analysis.

All other sensitivity artifacts live with the phase that produced them; see
`SENSITIVITY_INDEX.md` for the complete repository-wide index.
