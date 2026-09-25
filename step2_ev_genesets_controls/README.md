# Step 2 — Two-tier EV gene set and negative controls

The study's shared inputs, built and frozen before any hypothesis test:

- `genesets/geneset_sources/` — self-contained source JSONs: Tier-1 curated EV
  machinery, Tier-2 GO exosome-cargo sets, the curated biogenesis panel, the
  UniProt secreted background, and the four negative controls (mitochondrial,
  ribosomal, lysosomal/endosomal, general secreted).
- `genesets/README.md` — accessions, query parameters, and construction notes.
- `code/ev_genesets/` — the builder/sampler module: assembles the tiers and
  draws covariate-matched null gene sets (used by every downstream phase).
- `frozen/ev_genesets_frozen_v1.3.json` — the content-hashed frozen gene sets
  (v1.3) that all analyses load.
