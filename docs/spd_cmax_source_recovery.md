# SPD Cmax Source Recovery

Atlas recovers provenance for the human Cmax values in Sutherland et al.'s
Secondary Pharmacology Database (SPD) without inventing study context that is
not present in the public supplement.

## Rebuild

```bash
python -m analysis.cli.recover_spd_cmax_context \
  --out-dir data/external/spd/cmax_source_recovery
```

The command verifies the SPD workbook and the pinned public Smit source files
with SHA-256 checksums. It writes the measured Cmax inventory, the SPD-to-Smit
mapping, all contributing Smit source observations, an administration-context
review table, and a manifest recording counts and policy.

`refresh_phase1_pk_context` runs this recovery automatically when the cache is
absent. A previously reviewed context can be supplied explicitly:

```bash
python -m analysis.cli.refresh_phase1_pk_context \
  --spd-cmax-administration-context \
    data/external/spd/cmax_source_recovery/spd_cmax_administration_context.csv \
  --out-dir data/AtlasSPD_phase1/pk_context_<version>
```

## Source Contract

The public SPD workbook identifies the provider for each Cmax value but does
not include row-level PMID/DOI, dose, route, regimen, or formulation fields.
The source categories therefore have different admissibility:

- `NIBR curation`: the paper describes primary-literature curation at the
  highest approved dose, but the public row-level citation manifest is absent.
- `Pharmapendium`: SPD uses a third-quartile summary over heterogeneous records;
  it is not one study arm and upstream access restrictions may apply.
- `Smit et al.`: the public source can reproduce the SPD drug-level Cmax values,
  but those values are medians assembled from one or more observations.

Parser-derived dose, route, regimen, and formulation values are stored only in
`candidate_*` columns. Even a single contributing Smit row remains
`single_source_context_candidate`; it is not automatically training eligible.
This protects against upstream metadata errors and against assuming that a
number mentioned in prose belongs to the exact Cmax scenario.

Only an explicitly reviewed row with all of the following may be merged into
the primary SPD context:

- `context_status = verified_single_source_observation`
- `citation_status = verified_original_source_record`
- `training_allowed = true`
- exact DrugCentral ID, nonempty InChIKey, and Cmax-scenario identity
- a nonempty `original_source_record_id`, study identifier, and citation
- confirmed dose/route context and source-use rights

All other administration values remain audit or sensitivity context. Missing
context remains missing and is never inferred from another PK scenario.

## Clean Model-Table Merge

Merge a reviewed PK release into the pre-PK model freeze using a chemically
corroborated key. Do not merge into an already materialized older PK table and
do not overwrite primary conflicts.

```bash
python -m analysis.cli.merge_pk_context \
  --primary <pre_pk_model_ready.csv> \
  --pk-enriched <pk_release.csv> \
  --out <clean_merge>/model_ready.csv \
  --audit-dir <clean_merge>/pk_merge_audit \
  --key ligand_base,inchikey,target_uniprot,pdb_id

python -m analysis.cli.refresh_ml_feature_metadata \
  --run-dir <clean_merge> \
  --dataset <clean_merge>/model_ready.csv \
  --out-dir <feature_complete>

python -m analysis.cli.materialize_pk_mechanistic_features \
  --input <feature_complete>/model_ready.csv \
  --out <feature_complete>/model_ready_pk_mechanistic.csv \
  --target-endpoint spd_free_cmax_um
```

Use `--target-endpoint pk_context_free_cmax_um` only for a separately named
external-context sensitivity table. Endpoint-specific `pk_mech_*_training_allowed`
columns, not the generic availability of contextual fields, determine whether
administration values are eligible for a model.
