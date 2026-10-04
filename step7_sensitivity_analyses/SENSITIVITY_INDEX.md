# Sensitivity & robustness index

Every sensitivity/robustness artifact in the repository, what it tests, and
where it lives. Headline results are unchanged by these analyses unless noted.

## Step 4 — Phase 2 (single-cell sender EV-specificity)
| Artifact | What it tests |
|---|---|
| `step4_phase2_singlecell_gate2/data/downsampled_control_sensitivity.tsv` | Gate-2 verdict under a downsampled normal-control arm |
| `step4_phase2_singlecell_gate2/data/sender_effect_size_floor.tsv` | Sender specificity under an effect-size floor |
| `step4_phase2_singlecell_gate2/data/leave_one_dataset_out.tsv` | Leave-one-dataset-out stability of the Phase-2 result |
| `step4_phase2_singlecell_gate2/code/pillar2_robust_de.R`, `pillar2_robust_gsea.py` | Alternative DE / enrichment engines for the sender test |

## Step 5 — Phase 3 (ligand–receptor landscape)
| Artifact | What it tests |
|---|---|
| `step5_phase3_lr_landscape_gate3/data/robustness/cap_sensitivity_result.json` | Sender-cap sweep (1000/2000/4000/8000 cells per compartment) — landscape invariant |
| `step5_phase3_lr_landscape_gate3/data/robustness/scramble_result.json` | 100 label-scramble seeds → 0 consensus axes (null calibration) |
| `step5_phase3_lr_landscape_gate3/data/robustness/lodo_result.json` | Leave-one-donor-out stability |
| `step5_phase3_lr_landscape_gate3/data/gate3/` | Gate-3 evaluation copies of the above, plus CKD-persistence and pathway-coherence checks |
| `step5_phase3_lr_landscape_gate3/figures/F5_scramble_collapse.svg/.png` | Scramble-collapse figure |
| `step5_phase3_lr_landscape_gate3/code/pillar3_robustness.py` | Cap-sweep / scramble / LODO driver |

## Step 6 — Phase 4 (EV subtraction + synthesis)
| Artifact | What it tests |
|---|---|
| `step6_phase4_ev_subtraction_gate4/data/T5_expression_floor_sensitivity.tsv` | Expression-floor sweep (0/10/25/50th percentile): the EV-exclusive set (S5) is invariant; S2/S4 shrink at the 50th percentile (44→36, 34→26) |
| `step6_phase4_ev_subtraction_gate4/data/gate4_and_rule_sensitivity.tsv` | Gate-4 verdict under the alternative AND-rule configuration |
| `step6_phase4_ev_subtraction_gate4/data/pillar4_funnel_exprfloor_sensitivity.json` | Evidence funnel under the floor sweep |
| `step6_phase4_ev_subtraction_gate4/figures/FD_consensus_axis_robustness.svg/.png` | Consensus-axis robustness summary |
| `step6_phase4_ev_subtraction_gate4/figures/F_FUNNEL_SUPP_exprfloor_sensitivity.svg/.png` | Funnel under the floor sweep |

## Step 7 — this folder
| Artifact | What it tests |
|---|---|
| `data/TS77_shared_donor_exclusion.tsv` | Consensus axes after excluding the 10 AKI donors shared between atlases — only SPP1→ITGAV_ITGB1 survives both reduced arms |
| `data/TS78_external_replication_gse210622.tsv` | External replication in GSE210622 — SPP1 axes direction-level with CellChat significance; COL4A1 discordant; no full-criteria pass |
