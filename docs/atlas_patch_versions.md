# Atlas Patch Versions

Patch labels use `AtlasvMAJOR.MINOR.PATCH` followed by a short description. They identify repository patches and analysis interface changes; they are not claims that datasets or scientific results are immutable releases.

## Atlasv0.0.01 - Context-preserving Phase 1 PK ingestion

Adds a one-command PK refresh for AtlasSPD Phase 1. The workflow preserves every source-specific PK context row, writes a separately selected representative record, joins external values under `pk_context_*` columns, and does not recompute `spd_exposure_label`. SPD is primary; NCATS FRDB and PK-DB availability is recorded explicitly; DailyMed/openFDA requests are serial, cached, and low-confidence until numeric extraction is reviewed; DrugBank remains optional BYOL.

## Atlasv0.0.02 - Verified PK source health and target representation plan

Adds reproducible NCATS archive validation and PK-DB API consistency checks, a record-level DailyMed/openFDA source-text review queue, and quarantine of unreviewed SPL numeric values from representative ML features. The PK-DB health check distinguishes a live server defect from true zero coverage. It also documents the Phase 1 target-representation path: deterministic pocket descriptors, frozen ESM-2 sequence embeddings, pocket-residue pooling, and sequence/pocket-cluster validation.

DrugBank remains an optional licensed input. Academic access is requested through https://go.drugbank.com/academic_research; the current release page reports that academic downloads are temporarily paused.

## Atlasv0.0.03 - Reviewed SPL PK and target pocket features

Adjudicates every cached SPL Cmax candidate against its source excerpt, rejects known SD/AUC/other-analyte extraction failures, and emits a context-specific accepted table without promoting any value to a universal per-drug Cmax. PK-DB health now distinguishes the working PostgreSQL study metadata from the broken Elasticsearch output index, and NCATS retrieval tries every current/archive URL published on the official download page.

Adds a reproducible one-row-per-PDB target pocket feature table using label-independent residue composition, simple geometry, size, and metal-context descriptors. ESM-2 remains a frozen sequence-representation ablation; PocketVec remains a separate inverse-docking baseline rather than a name for the deterministic descriptors.
