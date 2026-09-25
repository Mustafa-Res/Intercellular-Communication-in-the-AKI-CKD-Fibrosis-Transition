# Pillar 1 Execution Plan — Consensus Bulk EV-Injury Signature + Gate 1

## 1\. Objective

Test **H1 (primary, falsifiable):** EV genes are coordinately dysregulated as a *reproducible* AKI→CKD transcriptional program whose enrichment **exceeds a matched non-EV secreted background**. Discovery = AKI + CKD-fibrosis (kidney tissue). Replication = DKD and IgAN, analyzed **separately** (never pooled).

\---

## 2\. Cohort selection

GEO screening (58 human kidney candidates) confirmed most hits are mechanistic/organoid/transplant-blood studies unfit for a disease-signature meta-analysis. Manual curation against inclusion criteria is required.

### 2.1 Inclusion criteria (pre-registered)

* Human **kidney tissue** (biopsy, nephrectomy, or micro-dissected tubulo-interstitium/glomerulus). NOT blood/plasma/urine/organoid/cell-line.
* Bulk RNA-seq (or microarray in a tagged sensitivity arm).
* Contains a **disease-vs-control contrast** with ≥3 per group (case AND control in the same series → within-series DE, avoids cross-study batch as the primary comparison).
* Raw or processed counts/intensities retrievable (GEO supplementary or SRA).

### 2.2 Exclusion criteria

* Non-kidney tissue; pure perturbation studies with no human disease/control axis; organoid/iPSC-only; transplant rejection-only without native disease; n<3/group; species = mouse only.

\---

## 3\. Per-cohort differential expression

* RNA-seq: `limma-voom`; microarray: `limma`. Design `\~ condition + known covariates` (sex, batch, platform where available). EPV respected.
* Gene ID harmonization → Ensembl gene IDs on a common axis; collapse duplicates by max mean expression.
* Control-type tagged per cohort (healthy nephrectomy margin vs. living-donor vs. "normal-adjacent") for the heterogeneity sensitivity meta-analysis.

## 4\. Meta-analysis

* Per-gene random-effects meta-analysis of standardized effect sizes (DerSimonian–Laird primary; REML sensitivity) via `metafor`.
* Heterogeneity: I², τ². Control-type sensitivity meta-analysis (exclude non-healthy controls) as pre-registered.
* Multiple testing: BH-FDR across genes.

## 5\. EV enrichment — the specificity test (Gate 1 core)

* **Primary:** is the Tier-1 EV set (152 genes) enriched among the reproducible meta-signature **beyond** the covariate-matched non-EV secreted background? Test via:

  * ORA (hypergeometric) of Tier-1 EV vs. background universe, AND
  * a **permutation/bootstrap null**: draw expression/length/detection-matched background gene sets (the `match\_background` sampler, already built), recompute the enrichment statistic B=10,000×, report empirical p and effect size.
* **Sensitivity:** repeat with Tier-2 broad GO cargo (2,115). Report both; headline = Tier-1.
* GSEA over EV sets on the ranked meta-statistic as an orthogonal enrichment view.
* Direction consistency: EV genes' effect-size sign agreement across cohorts.

## 6\. Gate 1 criteria

PASS requires:

1. Tier-1 EV enrichment FDR < 0.05, AND
2. enrichment > matched non-EV background (empirical permutation p < 0.05), AND
3. reproducible in ≥2 independent discovery cohorts (direction-consistent), AND
4. significant pathway enrichment (GO/Reactome exosome/MVB) on the meta-signature.
Convergence over raw significance: a direction-consistent EV program reproduced across cohorts + beating background outweighs a single-cohort low p-value.

