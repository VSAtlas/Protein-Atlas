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


## Atlasv0.0.04 - Stage SPD add-ons and matched controls

Adds the versioned SPD positive-add-on workflow, same-source control staging, and
fail-closed score-scale/provenance handling for small-library additions.

## Atlasv0.0.05 - Build the auditable docking-atlas foundation

Adds release-manifest ingestion, failure-complete receptor-ligand cells, immutable
SQLite provenance, explicit primary-score ranking contracts, and bidirectional
target/drug exports without inferring scientific approvals.

## Atlasv0.0.06 - Ship the conference atlas release builder

Adds `atlas publish audit|build`, a path-redacted public database and flat
downloads, receptor/drug/pair static exploration, readiness evidence, guided local
analysis, and deterministic conference image planning.
## Atlasv0.0.07 - Safe reusable SPD add-on integration

Adds a one-command SPD add-on merge with canonical FDA identity reconciliation,
duplicate/conflict gates, training quarantine, metadata refresh, and mandatory
post-ingestion audits. Adds `atlas ml score-addons` for exact staged PDB-ligand
pairs in a frozen full-library comparison-run Vina context. The scorer rejects reference
exhaustiveness below 2, caps workers at 32, refuses concurrent Atlas runs using
manifest-owned Slurm IDs, local process evidence, and a 72-hour fail-closed
manifest window, and holds an exclusive add-on lock. Stale historical lifecycle
records no longer block forever, and fast or DUD-only runtime shortcuts are not
accepted.

## Atlasv0.0.08 - Harden qualified and scalable docking-atlas releases

Adds hashed scientific annotation sources, separate receptor-quality and
native-redock gates, rerun-preserving normalized completion records, compact public
SQLite projection, deterministic bundle integrity verification, paginated static
exploration, and a bounded read-only Worker/R2 architecture prototype.

Corrupt completion inputs, duplicate canonical master results, unsafe edge output
locations, forged overwrite markers, untracked in-site verification reports, and
private-path leaks now fail closed. Pair artifacts are visible in static pages and
missing image artifacts become explicit preflight gaps; the full 760,878-cell edge
export remains blocked on a future SQLite-streaming implementation.

## Atlasv0.0.09 - Complete contextual PK and add-on discovery

Adds a source-specific NCATS FRDB PK adapter with unit-safe Cmax conversion,
structured dose/route/regimen/formulation and fraction-unbound context, reusable
source-text-adjudicated DailyMed scenarios, FDA-only openFDA query filtering, and
atomic replacement of stale derived PK columns. The SPD exposure label remains
unchanged. Timestamped SPD add-on merges can now select the newest completed
score-ready directory while rejecting and recording newer incomplete branches.

## Atlasv0.0.10 - Stream and qualify auditable docking-atlas delivery

Adds bounded-memory SQLite export with coarse pair shards, provider-neutral
download projection and deployment preflight, typed pair-level artifact roles,
fail-closed retry selection and causal-lineage readiness, exact prepared-receptor
APO/HOLO evidence with metal-audit acceptance/rejection provenance, selected-pose
validation gates, and the conference top-valid/best-invalid image plan.

## Atlasv0.0.12 - Classify and modularize the docking atlas

Classifies all frozen-run primary scores with exact row-bound legacy-decoy-Z
provenance, separates APO exploratory ranks from HOLO headline ranks, and keeps
both tracks unranked when exact chemistry or score-to-pose evidence is missing.
Public projections now omit private, unselected, unverified, and unhashed
artifacts. A symmetry-aware native-redock RMSD evaluator replaces the unsafe
ligand-fitted fallback but qualifies only with independently verified common-atom
receptor-frame evidence and exact all-model Vina score ordering; unsupported
evidence and unspecified potential stereochemistry remain pending. Fail-closed
single-stage repeat planning covers Vina, GNINA, SCORCH, MM/GBSA, and pose
validation. Full-matrix execution plus the stereo, resonance-group, aggregate-pose,
known-pair, and partial-sampling policies stay explicit publication gaps.
