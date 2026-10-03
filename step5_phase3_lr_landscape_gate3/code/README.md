# Pillar 3 — code 
Python/R pipeline for Pillar 3: build the proximal-tubule (PT) sender →  
receiver ligand–receptor landscape in human kidney single-cell atlases with two  
independent engines (LIANA rank-aggregate + CellChat v2) on identical cells,  
and test reproducibility, specificity, and robustness.

 Data source

- CZ CELLxGENE Census `2025-01-30`, same frozen query as Pillar 2, restricted  
  to the two multi-disease atlases (snRNA `a12ccb9b…`, scRNA `dea717d4…`) and  
  the four Pillar-3 compartments (PT sender + 3 receivers).

 Run order

1. `pillar3_celllevel_pull.py` — cell-level Census pull (cell x gene matrices  
   retained), one compressed `.h5ad` shard per dataset x condition.
2. `freeze_soluble_negctrl.py` — freeze the soluble-cytokine negative-control  
   panel **before any inference** (SHA-256 recorded).
3. `pillar3_prep_and_liana.py` — per-shard normalize + subsample  
   (cap 2,000 cells/compartment, seed = 1), run LIANA `rank_aggregate`  
   (CellPhoneDB + NATMI + Connectome + SingleCellSignalR + logFC, RRA  
   consensus), and export the identical subsampled cells for the CellChat run.
4. `pillar3_run_cellchat.R` — CellChat v2 engine per shard on the identical  
   cells (Seurat ingest of the `forR_<tag>` bundles).
5. `pillar3_differential_consensus.py` — differential (AKI-vs-normal,  
   CKD-vs-normal) + cross-engine / cross-dataset consensus; this is the  
   evidence-refinement funnel (1,311 → 44 → 34 → 3).
6. `pillar3_robustness.py` — robustness battery on the identical cells:  
   label-scramble specificity (Gate 3B) and related arms (Gate 3C).
7. `pillar3_ev_route_overlay.py` — route-attribution overlay  
   (soluble-classical / contact-dependent / EV-compatible), secondary and  
   descriptive, not gated.
8. `pillar3_ckd_persistence.py` — Gate 3D CKD-persistence sub-check for the  
   3 consensus axes.
9. `pillar3_gate3_eval.py` — evaluate the pre-registered Gate 3 criteria  
   (3A–3D) and write `gate3_verdicts.tsv`.


## Notes

- Engine thresholds (CellChat p < 0.05, absolute/relative probability floors,  
  LIANA rank/RRA floors) are defined in `pillar3_differential_consensus.py`  
  and were frozen before evaluation.
- Absolute paths recorded inside scripts reflect the original analysis  
  environment; adapt them before re-running.
