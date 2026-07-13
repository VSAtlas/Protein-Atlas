# Docking Atlas Release Workflow

This workflow turns existing Atlas run artifacts into an auditable SQLite release,
a path-redacted public snapshot, flat analysis tables, and a self-contained static
explorer. It does not choose receptors, map PDB structures to proteins, select known
pairs, or decide whether a native redock passes a scientific qualification policy.
Those choices must be made explicitly in the release manifest.

Start from [the example manifest](examples/docking_atlas_release.example.yaml) and
replace every angle-bracket placeholder. Do not publish a placeholder manifest.

## Audit and build

Audit inputs before building:

```bash
atlas publish audit \
  --manifest docs/examples/docking_atlas_release.example.yaml \
  --out-dir outputs/data/<RELEASE_ID>
```

Add `--strict` to return nonzero when any run lacks the manifest, master rows, or
completion records needed for a failure-complete import. The audit writes
`audit.json` and does not create a database.

Build the private database, public downloads, and static site:

```bash
atlas publish build \
  --manifest <RELEASE_MANIFEST.yaml> \
  --out-dir outputs/data/<RELEASE_ID>
```

Use `--overwrite` only to replace an existing generated site. Use `--no-parquet`
when a Parquet engine is unavailable. The build does not launch docking or repair
missing scientific inputs.

## Output layout

```text
outputs/data/<RELEASE_ID>/
├── audit.json                         # input audit, when run separately
├── build_summary.json                 # public build summary
├── build_summary.private.json         # canonical import summary
├── docking_atlas.private.sqlite       # canonical private evidence database
└── site/
    ├── index.html                     # portable explorer entry point
    ├── analysis.html
    ├── proteins/                      # receptor-context target pages
    ├── drugs/  pairs/                 # drug and pair detail pages
    ├── assets/                        # local CSS, JS, and analysis JSON
    ├── site_manifest.json
    └── downloads/
        ├── docking_atlas.sqlite       # public, path-redacted copy
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
whose machine-local paths and stored JSON paths are redacted. Do not replace the
private provenance ledger with the redacted copy.

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

- `final_score` is present;
- `pose_valid == 1`;
- the pair is neither a native control nor a decoy.

All other cells are shown without a rank and carry an eligibility reason such as
`missing_final_score`, `pose_invalid`, `pose_validation_missing`,
`native_control`, or `decoy`. `atlas_score` and its source remain available as
secondary normalized evidence, not a substitute primary rank.

## Release manifest

The manifest is a versioned scientific contract and is stored with the release.
It must contain a non-empty `runs` list and non-empty policy records for:

- `receptor_selection`: the externally chosen receptor-quality and inclusion rule;
- `native_redocking`: the externally chosen qualification method and thresholds;
- `failure_handling`: how the complete matrix and failure taxonomy are defined;
- `normalization`: the declared score fields, directions, and comparison scope.

Atlas preserves these mappings as data; their presence does not prove that a policy
was satisfied. Native redocking is counted as explicitly qualified only when the
stored receptor status is one of the controlled engineering tokens `qualified`,
`explicitly_qualified`, or `qualification_passed`. Assign such a token only after
the manifest's scientific policy has actually been applied. An observed valid
control result is not automatically a qualification.

Each run may be a run ID string or a mapping with `run_id` and optional path
overrides for `run_manifest`, `master_rows`, `docked`, and `archive_index`.
Overrides are integration locations, not scientific replacements.

### Known-pair syntax

Conference image planning accepts `image_plan.known_pairs`. Supply at most one
representative pair per receptor context. Identify the context with either
`receptor_context_id`, or with `run_id`, `pdb_id`, and optional `variant` and
`ph_label`. Identify the ligand with `canonical_id` (the aliases `ligand_id` and
`drug_id` are also accepted). Atlas reports absent, ambiguous, invalid, or duplicate
selections as gaps; it does not choose a known pair.

`atlas publish build` generates `downloads/image_plan.json`, hashes it, and adds it
to the explorer download manifest automatically. The plan records deterministic
`atlas screenshot` argument vectors for the native control, up to five top valid
scored ligands, and the supplied known pair when it is not already selected. It does
not render images, does not select the best invalid ligand, and does not invent
missing structures or scores. The advanced Python callable
`analysis.atlas_database.write_release_image_plan` can regenerate a plan for review,
but a public release should be rebuilt with `atlas publish build` so the download
hashes and site manifest remain consistent.

`release_readiness.json` is generated automatically during a site build. It is a
descriptive audit, not a scientific verdict. It reports pair/failure completeness,
valid-drug score coverage, pose validation, explicitly stored native-redock status,
receptor variants, artifacts/hashes/verification, and score-source classification.

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
