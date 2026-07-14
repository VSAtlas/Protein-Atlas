# Docking Atlas Release Workflow

This workflow turns existing Atlas run artifacts into an auditable SQLite release,
a path-redacted public snapshot, flat analysis tables, and a self-contained static
explorer. It does not choose receptors, map PDB structures to proteins, select known
pairs, or decide whether a native redock passes a scientific qualification policy.
Those choices must be recorded explicitly in the release manifest and in the
hashed scientific annotation sources that it references. Annotation CSV/JSON/YAML
files may remain external, versioned inputs; Atlas hashes their content and stores
the source record provenance during the build.
Track unresolved approvals in
[`docking_atlas_decision_register.md`](docking_atlas_decision_register.md) before
freezing a release.
For the completed historical SPD90 matrix, follow the guarded
[`docking_atlas_historical_backfill.md`](docking_atlas_historical_backfill.md)
handoff; its pose-validation command is explicitly a not-yet-implemented interface.

Copy [the example manifest](examples/docking_atlas_release.example.yaml) to a
release-specific, versioned path and replace every angle-bracket placeholder in
the copy. Keep the example unchanged and do not publish a placeholder manifest.

## Audit and build

Audit inputs before building:

```bash
atlas publish audit \
  --manifest <RELEASE_MANIFEST.yaml> \
  --out-dir outputs/data/<RELEASE_ID>
```

The audit always returns nonzero for a missing run manifest or an invalid
completion JSON/shape. Add `--strict` to also return nonzero when any run lacks
master rows or valid completion records needed for a failure-complete import. It
writes `audit.json` and does not create a database. It scans completion and
master-result identities before database creation, so duplicate successful
attempts, duplicate result rows, and unmatched explicit selectors are visible in
the audit and make strict mode fail. This command checks run inputs and annotation
declarations, but it does not parse annotation files or validate their context
keys; the build performs those annotation checks during ingestion and aborts on an
invalid record.

Generate exact prepared-receptor chemistry evidence separately. This records
observations and exact-context metal-audit acceptance or rejection; it never
assigns receptor-quality or redocking qualification:

```bash
atlas publish receptor-evidence \
  --manifest <RELEASE_MANIFEST.yaml> \
  --out <RECEPTOR_ANNOTATIONS.yaml> \
  --expected-context-count <FROZEN_COUNT> \
  --expected-inventory-sha256 <FROZEN_CONTEXT_KEYSET_SHA256>
```

Build the private database, public downloads, and static site:

```bash
atlas publish build \
  --manifest <RELEASE_MANIFEST.yaml> \
  --out-dir outputs/data/<RELEASE_ID>
```

Use `--overwrite` only to replace the generated site and private database already
present in that output directory. Use `--no-parquet` when a Parquet engine is
unavailable. The build does not launch docking or repair missing scientific
inputs.

## Output layout

```text
outputs/data/<RELEASE_ID>/
├── audit.json                         # input audit, when run separately
├── build_summary.json                 # public build summary
├── build_summary.private.json         # canonical import summary
├── docking_atlas.private.sqlite       # canonical private evidence database
└── site/
    ├── _headers                       # Cloudflare Pages cache/security policy
    ├── index.html                     # portable explorer entry point
    ├── analysis.html
    ├── proteins/                      # receptor-context target pages
    ├── drugs/  pairs/                 # drug and pair detail pages
    ├── assets/                        # local CSS, JS, and analysis JSON
    ├── site_manifest.json
    ├── release_inventory.json         # sorted size/hash/media-type inventory
    ├── release_checksums.sha256       # deterministic payload checksums
    ├── release_verification.json      # offline integrity report
    └── downloads/
        ├── docking_atlas.sqlite       # compact public, path-redacted projection
        ├── database_redaction.json    # private/public hashes and redaction counts
        ├── release_browser.json       # browser payload
        ├── release_readiness.json     # descriptive coverage checklist
        ├── image_plan.json            # deterministic screenshot plan and gaps
        ├── pairs.csv                  # every imported pair cell
        ├── pairs.parquet              # optional equivalent table
        └── export_summary.json
```

The site has no server-side database dependency: publish the complete `site/`
directory without moving individual files. Keep `docking_atlas.private.sqlite`
private and canonical. The downloadable SQLite database is a derived public copy
whose machine-local paths are redacted. To keep a full matrix downloadable, raw
per-row result JSON and completion payload JSON are omitted from this projection;
their normalized score/source columns and completion hashes remain. Do not replace
the full private provenance ledger with the compact public copy.

## Offline verification and deployment handoff

`atlas publish build` creates the Pages headers, deterministic inventory, checksum
file, and a passing verification report. Recheck the exact directory immediately
before upload:

```bash
atlas publish verify --site-dir outputs/data/<RELEASE_ID>/site
```

Use `--report <PATH>` to write the report outside the deployable `site/`
directory; only the canonical `site/release_verification.json` report may live
inside the bundle. A failed verification returns nonzero. The verifier requires no
credentials or network and detects:

- missing local HTML or manifest download links;
- files missing from the inventory or checksum list;
- changed file sizes or SHA-256 hashes;
- symlinks and private build files inside the public tree;
- machine-local `/stor/`, `/home/`, or `/tmp/` paths in text and SQLite fields.

The checksum and inventory files exclude themselves and the verification report to
avoid recursive hashes; they cover the deployable payload, including `_headers`.
Regenerate them through `atlas publish build`, not by editing release files.

### Cloudflare Pages

1. Run `atlas publish build`, then `atlas publish verify` locally.
2. Inspect `release_inventory.json` before choosing Pages. As checked on
   2026-07-13, Pages permits at most 20,000 files on Free and up to 100,000 files
   on paid plans; the paid limit requires `PAGES_WRANGLER_MAJOR_VERSION=4`. Every
   individual Pages asset must also be no larger than 25 MiB.
3. Use Pages only when both the file-count and per-file checks pass. Then select the
   contents of `outputs/data/<RELEASE_ID>/site/` as the asset/output directory. Do
   not upload its parent, the private SQLite database, or private build summaries.
4. Preserve `_headers`; it applies a restrictive same-origin security policy,
   revalidation for pages, and immutable caching for release downloads.
5. Treat a deployed bundle as immutable. Because download basenames are stable, use
   a release-specific Pages project/hostname or an immutable Pages deployment URL;
   do not replace files in place behind the same long-lived URL.
6. Retain `release_inventory.json`, `release_checksums.sha256`, and
   `release_verification.json` with the deployed site as the handoff evidence.

The historical publication candidate's compact public SQLite projection is
approximately 439 MB, so it cannot be uploaded as a Cloudflare Pages asset. Its CSV
and other downloads must also be checked individually. A passing Atlas integrity
report does not imply that a bundle satisfies provider limits.

### Object storage or another download host

R2 is not required when the selected static host can serve every generated file.
For Pages, however, any download above 25 MiB requires object storage or another
host. Mirror the exact verified files under `site/downloads/` to an immutable prefix
such as `releases/<RELEASE_ID>/downloads/`. Preserve the content types from
`release_inventory.json`, set immutable cache metadata, and never overwrite an
existing release prefix. No credentials, bucket operations, CORS policy, or public
URLs are inferred by Atlas.

The generated static site continues to use its local `downloads/...` links, so an R2
copy alone is an archive/mirror and does not make an otherwise oversized Pages
bundle deployable. Do not hand-edit the verified static payload. For the bounded
edge explorer, use the supported `--download-base-url` projection below; it validates
all declared hashes and emits release-qualified HTTPS links without changing the
static fallback. A future static-site overlay would be a separate build contract.

Atlas provides two non-deploying edge projections. Preserve the existing
browser-JSON path for bounded conference/intermediate snapshots:

```bash
atlas publish edge-bundle \
  --site-dir outputs/data/<RELEASE_ID>/site \
  --out-dir outputs/data/<RELEASE_ID>/edge \
  --download-base-url https://downloads.example.org
```

For larger matrices, stream the compact public SQLite snapshot into coarse pair
shards without generating one pair object or static HTML page per cell:

```bash
atlas publish edge-bundle \
  --database outputs/data/<RELEASE_ID>/site/downloads/docking_atlas.sqlite \
  --source-site-dir outputs/data/<RELEASE_ID>/site \
  --out-dir outputs/data/<RELEASE_ID>/edge-streamed \
  --batch-rows 1000 \
  --coarse-shard-rows 2000 \
  --download-base-url https://downloads.example.org
```

Run the same provider-neutral or provider-specific integrity preflight on either
result:

```bash
atlas publish preflight \
  --edge-dir outputs/data/<RELEASE_ID>/edge-streamed \
  --site-dir outputs/data/<RELEASE_ID>/site \
  --provider generic
```

Both paths keep the static site unchanged and emit one dependency-free SPA shell,
a GET/HEAD-only Worker, and an immutable JSON object tree. See
[the optional edge-delivery handoff](docking_atlas_edge_delivery.md) for routes,
object layouts, limits, cost tiers, and manual deployment safeguards.

The browser-JSON path remains capped at a 128 MiB payload and 50,000 pair cells.
The SQLite path uses bounded fetch batches, a disposable on-disk work database,
SQLite ranking windows, and configurable coarse shards. It supports the engineering
scale required by the audited 760,878-cell candidate, but that full export has not
yet been executed or deployment-preflighted. Do not generate the per-pair static
site solely to feed the streaming path.

The conference-safe fallback is to keep the complete verified `site/` directory
together on a static HTTP host whose file-count and asset-size limits have been
checked, or retain it as a local backup. Pages is suitable only for a bounded
bundle that passes its provider limits. For external downloads, use object storage
or another download host together with a supported link projection; a dynamic
application server remains optional.

Until an authoritative PDB-to-protein mapping is supplied, entries under
`proteins/` are receptor-context target pages identified by run, PDB, variant, and
pH—not claims of one canonical biological protein per page. The static generator
creates one page per receptor context, drug, and pair. A full
FDA-scale matrix should use the coarse-sharded streaming edge projection rather
than pre-generating every pair page; the final host and any additional pagination
or caching policy remain deployment decisions.

## Scientific and ranking contract

Every expected receptor–ligand cell remains visible. Invalid, unsuccessful,
unvalidated, and missing-result cells are not dropped, and structured status and
failure fields remain available in SQLite, CSV, Parquet, and pair pages.

`final_score` is the only primary ranking field. There is no fallback to
`atlas_score`, `consensus_score`, or a raw docking score. Every non-null primary
score also carries `final_score_source`, so consensus-Z, SCORCH-Z, and other
materialization families are not silently conflated. A pair is rank-eligible only
when all of the following are true:

- `receptor_classification` is exactly the controlled token `HOLO`, based on
  retained chemistry in the exact prepared receptor rather than a run-directory
  label;
- receptor quality has one of the controlled explicit qualification tokens;
- native redocking has one of the controlled explicit qualification tokens;
- `final_score` is present;
- `pose_valid == 1`;
- `pose_validation_scope == selected_final_score_pose`;
- the pair is neither a native control nor a decoy.

The requested run variant remains visible for provenance, but it cannot override
observed receptor content. An observed `APO` receptor is excluded even when it was
produced under a run labeled HOLO. Missing or unresolved chemistry is also excluded.
All other cells remain visible without a rank and carry an eligibility reason such
as `apo_receptor`, `receptor_classification_missing`,
`receptor_quality_not_qualified`, `native_redock_not_qualified`,
`missing_final_score`, `pose_invalid`, `pose_validation_missing`,
`pose_validation_scope_unqualified`, `native_control`, or `decoy`.

Legacy `pose_valid_any` values are retained as
`atlas_posebusters_any_stage_legacy_v1` with scope
`any_stage_for_ligand`. The browser discloses the checks performed by
`tools/pose_bust.py`: molecule loading, sanitization, valence, internal clash,
bond-length, bond-angle, and the protein-clash/relative-distance rule. It shows the
0.92 CLI default but labels the exact historical run cutoff as unresolved when it
was not recorded. A legacy any-stage pass does not prove that the pose associated
with `final_score` passed, so it cannot enter headline rankings until that
alignment is supplied. `atlas_score` and its source remain available as secondary
normalized evidence, not a substitute primary rank.

## Release manifest

The manifest is a versioned scientific contract and is stored with the release.
It must contain a non-empty `runs` list and non-empty policy records for:

- `receptor_selection`: the externally chosen receptor-quality and inclusion rule;
- `native_redocking`: the externally chosen qualification method and thresholds;
- `failure_handling`: how the complete matrix and failure taxonomy are defined;
- `normalization`: the declared score fields, directions, and comparison scope.

Atlas preserves these mappings as data; their presence does not prove that a policy
was satisfied. Receptor quality and native redocking are counted as explicitly
qualified only when their respective stored status is one of the controlled
engineering tokens `qualified`, `explicitly_qualified`, or
`qualification_passed`. Assign such a token only after the manifest's scientific
policy has actually been applied. An observed valid control result is not
automatically a qualification.

Each run may be a run ID string or a mapping with `run_id` and optional path
overrides for `run_manifest`, `master_rows`, `docked`, and `archive_index`.
Overrides are integration locations, not scientific replacements.

### Explicit scientific annotation sources

The optional `annotations` mapping accepts three independent CSV, JSON, or YAML
sources. A value may be a path string or `{path: ...}`; relative paths resolve
from the release manifest directory.

For a first release, copy the three header-complete templates in `docs/examples/`:
`docking_atlas_proteins.example.csv`, `docking_atlas_receptors.example.csv`, and
`docking_atlas_known_pairs.example.csv`. Replace every angle-bracket placeholder;
the examples intentionally contain no biological mapping, classification, pass/fail
outcome, threshold result, or known-pair choice. Use an empty CSV containing only
the header when one annotation category is deliberately not yet supplied.

- `proteins` records require `protein_key` and at least one of `uniprot_id`,
  `gene_symbol`, or `display_name`.
- `receptors` records require the exact `run_id`, `pdb_id`, `variant`, and
  `ph_label` fields of an imported receptor context. `variant` and `ph_label`
  must be present even when the imported value is blank. Each row must provide
  at least one of the optional explicit fields `protein_key`,
  `receptor_classification`, `qualification_status`, `qualification_reason`,
  `native_redock_status`, `native_redock_reason`, and `native_redock_rmsd`.
- `known_pairs` records require the same exact receptor-context key plus
  `ligand_canonical_id`. Optional fields include `selection_label` and
  `evidence_reference`.

YAML/JSON sources may be a list or a mapping containing `records: [...]`.
Every stored row retains the source path, source SHA-256, one-based record index,
original record JSON, and optional `provenance` value. JSON/YAML may supply a
structured provenance object; CSV values, including JSON-looking text, are
preserved as strings rather than parsed. Protein identities,
receptor contexts, receptor annotations, and known-pair selections remain separate
relational records. Atlas rejects duplicate context selections, unknown protein
keys, incomplete context keys, and context keys that do not match exactly. It does
not fill an omitted variant or pH, infer a PDB-to-UniProt mapping, reinterpret APO
as HOLO, decide that a redock passed, or choose a known pair.

The private and public databases store these records in `protein_identities`,
`receptor_annotations`, and `known_pair_selections`. Public copies path-redact
annotation provenance just like run provenance. A known-pair selection may name a
ligand absent from that receptor's matrix; it is retained with a null `pair_cell_id`
so downstream image planning reports the gap rather than silently replacing it.

### Known-pair syntax

Frozen public releases must use `annotations.known_pairs` as described above.
`atlas publish build` loads known-pair selections only from
`known_pair_selections` rows ingested from that hashed external source; duplicate
database-backed context selections reject the build. Legacy inline
`image_plan.known_pairs` values may remain in older manifests and are covered by
the manifest hash, but the current public builder ignores them. The advanced Python
callable `write_release_image_plan(..., known_pairs=...)` accepts explicit
compatibility input for review calls; it is not used by the public build.
If an annotation-backed selection is absent, invalid, ambiguous, or unrenderable,
the image plan records a structured gap such as `known_pair_not_supplied` rather
than falling back or choosing a pair.

`atlas publish build` generates `downloads/image_plan.json`, hashes it, and adds it
to the explorer download manifest automatically. The plan records deterministic
`atlas screenshot` argument vectors for a unique, valid native-control pose with a
recorded native-redock status, up to five top valid scored ligands, the best-scoring
invalid ligand when present, and the supplied known pair when it is not already
selected. It does not render images or invent missing structures or scores.
Missing prerequisites become structured gaps. Every selected row also receives a
`pair_artifact_not_indexed` or `pair_artifact_not_verified` gap until a
pair-linked artifact with the controlled `docking_pose` role is verified. This is
an engineering identity/hash preflight, not a scientific pose-validity decision;
execute screenshot commands only after checking the resolved receptor and pose
inputs. The advanced Python callable
`analysis.atlas_database.write_release_image_plan` can regenerate a plan for review,
but a public release should be rebuilt with `atlas publish build` so the download
hashes and site manifest remain consistent.

`release_readiness.json` is generated automatically during a site build. It is a
descriptive audit, not a scientific verdict. Readiness schema 5 reports
pair/failure completeness, qualified-headline score coverage, PoseBusters coverage,
selected-final-score-pose alignment, controlled APO/HOLO coverage, chemistry
evidence state and run-label conflicts, explicitly stored receptor-quality and
native-redock statuses, artifacts/hashes/verification, rank-exclusion reasons, and
score-source classification. It also blocks publication when a selected result
lacks a structured causal link to its matching selected completion attempt. A
controlled label is evidence, not a receptor-quality or native-redock
qualification.

Completion manifests are normalized in the private database. Each source file and
its JSON payload are stored once in `completion_records`; pair-level
`docking_attempts` reference that row by `completion_record_id`, so retries with
the same engine, stage, and chunk remain separate attempts. Every master CSV row is
stored separately in `result_attempts` with its input-table path/hash, one-based CSV
line number, canonical row hash, parsed fields, and original row JSON. The selected
row is then materialized into `pair_cells`; CSV order never chooses a silent winner.

When a pair has multiple successful completion attempts or multiple result rows,
the build fails closed unless that run's `attempt_selections` entry identifies one
candidate exactly. Completion selection requires `completion_relpath` relative to
`paths.docked` plus `completion_sha256`. Result selection requires
`result_row_number` (the logical CSV record number, with the header counted as
record 1) plus `result_sha256`, computed from Atlas's sorted compact JSON
representation of the CSV row. A selector may contain either pair of fields or both.
The receptor key
(`pdb_id`, exact `variant`, exact `ph_label`) and `ligand_canonical_id` are always
required. Unique candidates are selected automatically and record that fact;
ambiguous, unmatched, unsafe-path, partial-hash, and duplicate selectors abort.
An ambiguity error prints each candidate's safe relative path or logical record
number together with the SHA-256 values needed for the manifest entry.
Both selected and unselected attempts retain `selected_for_release`,
`selection_method`, and the one-based selector index for auditability.
Explicitly choosing one completion and one result does not prove that the selected
completion produced that result. The nullable result-attempt fields
`completion_record_id`, `completion_link_method`, and
`completion_link_evidence_json` preserve this distinction; no link is inferred
from matching names, row order, or independent selectors.

Corrupt, non-object, or structurally unsupported completion JSON aborts the input
audit/build instead of being counted as failure-complete. The private database
retains raw completion and result-attempt JSON. The compact public projection omits
those raw payloads and redacts their source paths while keeping hashes, row numbers,
statuses, and selection provenance.

Artifact archive indexes are also fail-closed release inputs. Every entry needs a
controlled `artifact_role`, content SHA-256, non-negative size, file type, archive
path, and member name. Controlled roles are `source_receptor`,
`prepared_receptor`, `source_ligand`, `prepared_ligand`, `docking_pose`,
`pose_image`, `receptor_validation_report`,
`native_redock_validation_report`, `pose_validation_report`,
`search_box_definition`, `run_configuration`, `raw_stdout`, `raw_stderr`, and
`software_environment`. The role determines a `run`, `receptor`, `ligand`, or
`pair` scope. Ligand- and pair-scoped entries must explicitly name
`ligand_canonical_id`; Atlas no longer guesses association from a filename. A
`raw_stdout` or `raw_stderr` artifact is pair-scoped so the private ledger can
reconstruct one receptor-ligand calculation rather than attaching a log to a
whole receptor by inference. A pair-scoped artifact for an unscheduled pair
remains preserved with a null
`pair_cell_id` and never creates a phantom matrix cell. Role typing records artifact
identity only: it does not decide scientific validity or whether that role may be
published.

Database storage is therefore proportional to the number and size of completion
files plus the number of master result rows, rather than multiplying each completion
payload by the number of expected ligands in the file.

## Unresolved prerequisites before a public scientific release

The following data must be supplied or completed; the builder does not infer them:

1. A frozen run list and an authoritative PDB-to-protein/target mapping for the
   labels intended for public display.
2. An explicit experimentally determined receptor-selection policy and its recorded
   outcome for every included receptor context.
3. Native-ligand redocking evidence for every receptor, including the declared
   qualification calculation and an explicitly stored outcome.
4. A `pose_valid` result for every pair that may be ranked; unvalidated pairs remain
   visible but are never ranked.
5. Materialized `final_score` and `final_score_source` for every pose-valid drug pair
   intended for ranking. Failure cells may legitimately have no score.
6. A complete expected-pair universe from completion manifests, with terminal status
   and a structured reason for every unsuccessful or invalid calculation.
7. Pair-linked pose/input artifacts with SHA-256 hashes and verification state for
   every artifact promised for download or reconstruction.
8. Ligand identity/state provenance—including protonation, tautomer, stereochemistry,
   formal charge, and source/prepared hashes—where the release claims those fields.
9. One explicitly reviewed known-pair selection per receptor context where the
   conference image set requires it.
10. Restorable receptor/pose artifacts and verified variant/pH resolution before
    executing image-plan screenshot commands.
11. Causal completion linkage for every selected result attempt.
12. A production download origin/provider, immutable release prefix, and measured
    full-scale streaming rehearsal. The provider-neutral coarse-shard architecture
    is selected; account creation alone does not authorize upload or deployment.

Treat every readiness blocker as a concrete missing-coverage item. Resolve it in the
source run artifacts or manifest, rebuild the private database, and regenerate the
public projection rather than editing public outputs by hand.
