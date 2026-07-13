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
writes `audit.json` and does not create a database. This command checks run inputs
and annotation declarations, but it does
not parse annotation files or validate their context keys; the build performs
those annotation checks during ingestion and aborts on an invalid record.

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

The generated site continues to use its local `downloads/...` links, so an R2 copy
is an archive/mirror rather than an automatic link rewrite. It does not make an
otherwise oversized Pages bundle deployable. Serving downloads directly from an R2
or other domain requires a future supported URL-projection build option and an
origin-specific CORS decision; do not hand-edit the verified payload. Until that
projection exists, serve the complete verified site from a host that accepts its
files or retain it as a local conference backup.

For bounded conference/intermediate snapshots, Atlas also provides a separate,
non-deploying edge delivery projection:

```bash
atlas publish edge-bundle \
  --site-dir outputs/data/<RELEASE_ID>/site \
  --out-dir outputs/data/<RELEASE_ID>/edge
```

This keeps the static site unchanged and emits one dependency-free SPA shell, a
GET/HEAD-only Worker, and an immutable R2 JSON object tree. See
[the optional Cloudflare Worker + R2 handoff](docking_atlas_edge_delivery.md) for
routes, object layout, limits, cost tiers, and manual deployment safeguards.

The current command is capped at a 128 MiB browser payload and 50,000 pair cells;
it explicitly does not support the audited 760,878-cell publication candidate.
Publication scale requires a future SQLite-streaming exporter that bypasses both
the monolithic `release_browser.json` and per-pair static HTML generation. Do not
run the static build solely to feed the bounded edge command.

The conference-safe fallback is to keep the complete verified `site/` directory
together on a static HTTP host whose file-count and asset-size limits have been
checked, or retain it as a local backup. Pages is suitable only for a bounded
bundle that passes its provider limits. For external downloads, use object storage
or another download host together with a supported link projection; a dynamic
application server remains optional.

Until an authoritative PDB-to-protein mapping is supplied, entries under
`proteins/` are receptor-context target pages identified by run, PDB, variant, and
pH—not claims of one canonical biological protein per page. The current generator
creates one page per receptor context, drug, and pair. Whether a full FDA-scale
matrix should continue using pre-generated pages, use pagination, or move to
client-side/on-demand rendering is a future engineering and hosting decision.

## Scientific and ranking contract

Every expected receptor–ligand cell remains visible. Invalid, unsuccessful,
unvalidated, and missing-result cells are not dropped, and structured status and
failure fields remain available in SQLite, CSV, Parquet, and pair pages.

`final_score` is the only primary ranking field. There is no fallback to
`atlas_score`, `consensus_score`, or a raw docking score. A pair is rank-eligible
only when all of the following are true:

- the imported receptor variant is exactly `HOLO`;
- `receptor_classification` is blank while the controlled vocabulary and
  run/source conflict policy await user approval;
- receptor quality has one of the controlled explicit qualification tokens;
- native redocking has one of the controlled explicit qualification tokens;
- `final_score` is present;
- `pose_valid == 1`;
- the pair is neither a native control nor a decoy.

All other cells are shown without a rank and carry an eligibility reason such as
`apo_receptor`, `receptor_variant_not_holo`,
`receptor_classification_policy_pending`, `receptor_quality_not_qualified`,
`native_redock_not_qualified`, `missing_final_score`, `pose_invalid`,
`pose_validation_missing`, `native_control`, or `decoy`. `atlas_score` and its
source remain available as secondary normalized evidence, not a substitute
primary rank.

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
recorded native-redock status, up to five top valid scored ligands, and the supplied
known pair when it is not already selected. It does not render images, does not
select the best invalid ligand, and does not invent missing structures or scores.
Missing prerequisites become structured gaps. Every selected row also receives a
`pair_artifact_not_indexed` or `pair_artifact_not_verified` gap until a
pair-linked artifact is verified. That artifact flag is an engineering preflight
signal, not proof that the artifact is the required pose member; execute screenshot
commands only after checking the resolved receptor and pose inputs. The advanced
Python callable
`analysis.atlas_database.write_release_image_plan` can regenerate a plan for review,
but a public release should be rebuilt with `atlas publish build` so the download
hashes and site manifest remain consistent.

`release_readiness.json` is generated automatically during a site build. It is a
descriptive audit, not a scientific verdict. It reports pair/failure completeness,
qualified-headline score coverage, pose validation, distinct explicitly stored
receptor-quality and native-redock statuses, receptor variants,
artifacts/hashes/verification, rank-exclusion reasons, score-source classification,
and the schema-3 `receptor_classification_policy_approval` blocker for any
non-empty classification recorded before the controlled vocabulary and conflict
policy are approved.

Completion manifests are normalized in the private database. Each source file and
its JSON payload are stored once in `completion_records`; pair-level
`docking_attempts` reference that row by `completion_record_id`, so retries with
the same engine, stage, and chunk remain separate attempts. Corrupt, non-object, or
structurally unsupported completion JSON aborts the input audit/build instead of
being counted as failure-complete. Duplicate canonical receptor-ligand result keys
in a master CSV also abort the build until an explicit attempt/result selection is
provided; CSV order never chooses a silent winner.

Database storage is therefore proportional to the number and size of completion
files, rather than multiplying each payload by the number of expected ligands in
the file. The source path and SHA-256 remain available for provenance, while the
public projection redacts the private path and raw payload.

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
11. A hosting and scaling decision, including whether the complete static page set is
    acceptable for the final matrix size and browser/device targets.

Treat every readiness blocker as a concrete missing-coverage item. Resolve it in the
source run artifacts or manifest, rebuild the private database, and regenerate the
public projection rather than editing public outputs by hand.
