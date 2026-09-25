# Pre-Registration Protocol — Intercellular Communication in the AKI→CKD Fibrosis Transition

## 1\. Scientific rationale

Extracellular vesicles (EVs) mediate intercellular communication in kidney injury by transferring protein/RNA cargo between resident and infiltrating cell types. Generic ligand–receptor cell–cell communication (CCC) tools cannot detect EV cargo transfer. The AKI→CKD fibrosis transition (\~30–40% of AKI progresses) is driven in part by tubular-epithelial EV signaling to macrophages/fibroblasts. This project builds a multi-evidence, convergence-based case that a *specific, reproducible, genetically-anchored, patient-validated* EV-cargo communication program participates in the AKI→CKD transition — distinguished from prior art by integrating genetics × EV-cargo × urinary-EV × spatial adjacency plus an EV-vs-background specificity test.

## 2\. Design philosophy

* **PROGRESSIVE 5-pillar design:** build each layer only if the prior gate passes ("saves months"). Not a fixed 3- or 5-pillar commitment.

## 3\. Phenotype scope

* **PRIMARY / discovery:** AKI→CKD fibrosis transition.
* **Replication cohorts (SEPARATE, never pooled):** diabetic kidney disease (DKD) and IgA nephropathy (IgAN).
* Approved phrasing: **"AKI-associated and CKD-associated fibrosis signatures representing the AKI-to-CKD transition spectrum."** Temporal/transition inference is reserved for integration with single-cell (Census), spatial (KPMP), and established biology — never claimed from bulk cross-sectional data alone.

## 4\. Hypotheses

* **H1 (primary, falsifiable):** EV genes are coordinately dysregulated as a reproducible AKI→CKD program whose enrichment exceeds a matched non-EV secreted background.
* **H2:** EV-aware CCC identifies specific sender→receiver axes (hypothesized: injured proximal tubule→macrophage; PT→fibroblast; endothelial→immune).
* **H3:** Spatially, sender/receiver pairs are adjacent beyond chance.
* **H4:** EV-program genes are genetically enriched (GWAS/eQTL/pQTL) above background.
* **H5:** EV cargo candidates are detectable in urinary-EV proteomes (detection filter, not case/control).
* **H6 (exploratory):** The EV program is pharmacologically reversible (drug repurposing).

## 5\. Datasets, inclusion/exclusion, and data sources (per pillar)

### 5.1 Pillar 1 — Bulk (discovery + replication)

* Source: GEO human kidney-tissue bulk RNA-seq (microarray in a tagged sensitivity arm).
* Availability verified: AKI=67, CKD=112, fibrosis=93, DKD=59, IgAN=15 human RNA-seq series.
* Detailed inclusion/exclusion and cohort curation: see `PLAN\\\_pillar1.md`.
* **COHORT CONFIDENCE TIERING — governs which cohorts touch the primary result.** Each cohort scored on 5 criteria (+1 each): kidney tissue confirmed; phenotype clear at sample level; appropriate control group; platform quality/adequate n; clinical metadata present. Tiers: **A = primary discovery** (score 5, appropriate control, clear phenotype, no caveat); **B = sensitivity** (score 3–4, OR score 5 with a tagged caveat e.g. transplant/FFPE/tumor-adjacent/few-controls); **C = exploratory** (score ≤2 OR non-standard contrast e.g. progression-only/no healthy control). Weak cohorts must NOT silently influence the headline: the primary meta-analysis is built on Tier-A cohorts; Tier-B enters pre-specified sensitivity analyses; Tier-C is exploratory only. Implemented in `src/cohort\\\_scoring.py`; table saved to `02\\\_bulk\\\_signature/cohorts\\\_selected.csv`.
* **Locked selection:** Tier-A discovery = GSE139061 + GSE30718 (both native-kidney AKI, independent platforms) → satisfies Gate-1 ≥2-cohort reproducibility on the AKI axis. CKD discovery = ERCB cohorts GSE104954/GSE99340/GSE47185 (Tier-B, compartment/control caveats). Replication (separate, never pooled): DKD = GSE142153 (Tier-A) + GSE294519 (Tier-B); IgAN = GSE104948 (Tier-B). Progression contrast = GSE339456 (Tier-C, exploratory). NOTE: GSE246206 does NOT provide a 2nd clean AKI cohort (its mRNA arm = GSE139061 + GSE246204 minimal-change-disease); GSE271143 excluded (no sample-level phenotype).
* **Claim framing:** AKI cohorts = acute-injury state; CKD cohorts = failed-repair/fibrosis state; the AKI-to-CKD transition is INFERRED by integration (single-cell/spatial/known biology), not claimed from longitudinal samples.

### 5.2 Pillar 2 — Single-cell (PRIMARY source = CELLxGENE Census, pinned `2025-01-30`)

* `value\\\_filter = "tissue\\\_general == 'kidney' and is\\\_primary\\\_data == True"`.
* Disease labels: `'normal'`, `'chronic kidney disease'`, `'acute kidney failure'`, `'type 2 diabetes mellitus'`.
* Verified counts: normal=597,379; CKD=266,542; AKI=125,291; T2DM=11,509 (1,367,163 total across 21 datasets).
* Prefer the 2–3 datasets carrying AKI+CKD+normal *together* (within-dataset case/control → reduces batch confounding). Top: `a12ccb9b-4fbe-457d-8590-ac78053259ef` (snRNA, KPMP atlas: AKI 54,549 / CKD 127,996 / normal 122,444), `dea717d4-7bc0-4e46-950f-fd7e1cc8df7d` (scRNA companion: AKI 51,167 / CKD 130,630 / normal 43,380).
* KPMP raw controlled-access files = supplementary only (no public REST API).

### 5.3 Pillar 3 — Spatial

* KPMP Atlas (Visium/CODEX/IMC): 203 Visium + 164 CODEX + 81 IMC datasets.

### 5.4 Pillar 4 — Genetics

* GWAS: CKDGen portal (eGFR/CKD summary stats; portal reachable).
* eQTL (kidney): Susztaklab/NephQTL/GTEx kidney — confirm access at Pillar 4.
* pQTL (plasma): UKB-PPP/deCODE — confirm access at Pillar 4.

### 5.5 Pillar 5 — Urinary-EV validation + drugs

* PRIDE urinary-EV proteomes (all verified live): PXD050720 (ADTKD-HNF1β), PXD065380 (CKD/transplant), PXD042533 (ADPKD), PXD009784 (nephronophthisis), PXD025547 (medullary sponge kidney). Framed as a cargo-DETECTION filter, not matched case/control.
* Drug repurposing: LINCS/CMap + Open Targets (only if Gate 3 passes).

## 6\. Gate criteria

* **Gate 1:** Tier-1 EV enrichment FDR<0.05 AND > matched non-EV background AND reproducible ≥2 cohorts AND consistent direction AND significant pathway enrichment.
* **Gate 2:** ≥1 EV-aware CCC axis + cargo detected + reproducible in an independent dataset + spatial adjacency permutation FDR<0.05 + ≥5 cargo candidates carry forward.
* **Gate 3:** genetic enrichment above background + ≥3 of 4 orthogonal supports (GWAS, eQTL/pQTL, urinary-EV, independent replication).
* **Gate 4:** claim rises above "EVs participate" to a specific, novel, genetically-anchored, patient-validated program with a tractable intervention.
* **General:** BH-FDR<0.05 default; Bonferroni for small confirmatory families; feature selection inside CV folds; EPV≥10 for regressions.

## 7\. QC procedures

* Bulk: library-size/outlier checks, PCA, sex-check, batch diagnostics.
* Single-cell: per-dataset QC (mito%, doublets, min genes/counts), harmonized cell-type labels.
* Spatial: spot/segmentation QC.
* Genetics: allele harmonization, ancestry-matched LD, LD-mismatch diagnostics for fine-mapping.

## 8\. Statistical analysis plan

* Pillar 1: limma-voom (RNA-seq)/limma (microarray) per-cohort DE, design `\\\~condition + batch + covariates` → random-effects meta-analysis (DerSimonian–Laird primary; REML sensitivity via metafor) → ORA (hypergeometric) + GSEA over frozen EV sets vs. matched background. Control-type heterogeneity → tagged sensitivity meta-analysis.
* **Meta-analysis effect-size definition:** per gene, effect size = **log2 fold change**; variance = **(standard error from the limma model)²**. Meta-analysis outputs per gene: **pooled log2FC, 95% CI, Cochran's Q, I²**. Gene ranking: **PRIMARY = FDR-adjusted pooled p-value; SECONDARY = absolute pooled effect size** (prevents statistically-significant-but-biologically-tiny genes from dominating).
* **Negative-control gene sets — run BEFORE any EV-specificity claim.** Test the same enrichment statistic (with the same covariate-matched null) on: (C1) general secreted proteins; (C2) mitochondrial proteins; (C3) ribosomal proteins; (C4) lysosomal/endosomal proteins. Interpretation: if only EV core is enriched → strong specificity; if all trafficking systems rise together → broader stress response.
* Pillar 2: EV-aware CCC (miRTalk primary; EV-Net/NicheNet-adapted); CellChat/CellPhoneDB as generic comparator only.
* Pillar 3: neighborhood/adjacency permutation tests.
* Pillar 4: genetic enrichment vs. background; selective MR/coloc on converged candidates.
* Multiple testing, EPV, CV discipline as in §6.

## 9\. Integration strategy (convergence-tier scoring)

* Common Ensembl gene axis (proteins via UniProt→gene; miRNAs→validated targets).
* **Tier 1 = ≥3-of-4 evidence axes + spatially adjacent → advances to MR/coloc.**
* Convergence score aggregates: bulk meta-signature, single-cell CCC, spatial adjacency, genetics, urinary-EV detection.

## 10\. Validation strategy

* Multi-evidence convergence with association language.
* Formal causal inference (MR/colocalization) on top 5–10 fully-converged candidates only.
* Urinary-EV proteomes as an orthogonal cargo-detection filter.

