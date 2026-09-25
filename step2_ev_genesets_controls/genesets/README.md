# `00_preregistration/genesets/`

Records that make external-data retrieval reproducible. **No controlled-access
data are stored here** — only accession identifiers, query parameters, and the
frozen gene-set *source* JSONs from which the analysis gene sets were built.

- `geneset_sources/` — the raw inputs that define the EV gene sets, kept so the
  frozen sets can be rebuilt from scratch:
  - `tier1_machinery.json` — curated Tier-1 EV biogenesis/secretion machinery.
  - `curated_biogenesis.json` — curated EV-biogenesis support set.
  - `go_ev_sets.json` — Tier-2 broad GO exosome/EV-cargo sets (QuickGO).
  - `uniprot_secreted.json` — UniProt secreted-protein background.
  - `negative_controls.json` — mitochondrial / ribosomal / lysosomal-endosomal /
    general-secreted negative-control gene sets.

The *frozen, versioned* gene-set artifact actually consumed by the pillars is
`00_preregistration/feasibility/ev_genesets_frozen_v1.3.json` (content-hashed).
