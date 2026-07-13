System Prompt: Atlas2 Dev Agent
1. CORE CONSTRAINTS & SAFETY
Scope: Edit ONLY code/protein_automation/ and docs/.

Protected: tests/, ci/, input_pdbs/, docked/.

Root: /stor/home/mpg2352/atlas2/code/protein_automation (or atlas/...).

File Ops:

NO: cat large files (.pdb, .sdf, .csv), ls -R, find ., tree.

YES: head -n 5, csvstat, REPO_MAP.txt.

Temp: Create temp dirs ONLY under /stor/home/mpg2352 (NEVER /tmp).

Env: direnv active. Ignore lockfile errors. Run micromamba activate docking-env if needed.

External tools: install third-party tools under `/stor/work/VDS_Beckham/atlas/tools` (or a user-specified external tools root), not inside this repo's source, data, or runtime output trees. For open-source tools, make the install reproducible through the existing tool installer bundle under `tools/installers/` so a fresh user can clone/download the repo and run one documented install command. Proprietary/BYOL tools must remain unbundled and be documented as explicit local path/config inputs.

Canonical Atlas workflows:

- Fresh-clone/public smoke: `atlas smoke public` (installer path: `bash tools/installers/install_atlas_publication_stack.sh --public-smoke`).
- First real run planning: `atlas first-run --uniprot P04637 --ligands fda --small-library --fast --dry-run`; use `--yes --no-dry-run` only when launching is explicitly intended.
- Report/analysis handoff: `atlas analysis report <RUN_ID> --status-html`; this wraps master export and report generation.
- Benchmark/throughput analysis: use `atlas analysis dud-eval <RUN_ID>`, `atlas analysis throughput integrity <RUN_ID>`, and `atlas analysis interactions export <RUN_ID>` instead of raw `python -m analysis.cli.*` commands.
- Reproducibility handoff: `atlas reproduce bundle <RUN_ID>`; reuse `outputs/data/<RUN_ID>/`, do not create a separate `outputs/repro` tree.
- Slurm workflows: use `atlas slurm submit|progress|finalize`, not raw `tools/slurm/*.sh` or `tools/finalize_distributed_run.py`, unless debugging the wrapper itself.
- Run debugging: use `atlas debug <RUN_ID> --deep` or `atlas status <RUN_ID> --errors --explain`; do not inspect raw logs.
- Artifact sizing: use `atlas artifacts measure --run-id <RUN_ID>`.
- Developer verification: use `atlas dev verify --full --fix --smoke` as the canonical gate. If repo-wide debt blocks it, report the failing components and run focused checks for touched files.

2. DEV PROTOCOL
Use multiple agents for complex or multistep tasks

Versioned commit/push handoff:

- For each completed, verified Atlas patch, assign the next `AtlasvMAJOR.MINOR.PATCH` label and append a concise description to `docs/atlas_patch_versions.md`.
- Commit only files changed for the current patch. Never stage `data/`, `outputs/`, caches, downloaded external archives, protected directories, backup files, or unrelated dirty-worktree changes.
- Use the commit subject `AtlasvMAJOR.MINOR.PATCH: <concise description>`.
- Push the current branch after verification when ongoing push was requested. Report the commit hash, branch, and push result.
- Do not commit or push a patch whose required focused checks or smoke workflow failed; report the blocker instead.

Explore (Low-Token):

Read REPO_MAP.txt ? _skeletons/ (grep definitions) ? Source (only if strictly necessary).

Edit (Surgical):

No wide refactors. Import src.path_router (never hardcode paths).

Respect APO/HOLO and pH.

When a file grows beyong 100kb, do not add more lines of code, attempt to either split the file into more reasonable chunks or simply write another file and wire it in.

do NOT use more than 32 cores at a time for docking runs

Verify (MANDATORY):

4. **Verification**
  - Prefer `atlas dev verify --full --fix --smoke` for the full canonical gate.
  - The full gate includes `ruff check --fix .`, `vulture .`, `tach check`, `python tools/architecture_gate.py`, `mypy .`, and the smoke workflow.
  - If the full gate fails on known repo-wide debt, run focused equivalents for touched files and report the full-gate blockers clearly.
  - Run the most relevant tests.
  - Report: commands run, pass/fail, and where outputs are written.

3. DEBUGGING & DATA
Logs: Use `atlas debug <RUN_ID> --deep` or `atlas status <RUN_ID> --errors --explain`. DO NOT read raw logs.

Data: Use csvstat or csvgrep.

Bugs: Found a bug? 1. Create regression test. 2. Fix bug. 3. Verify pass.

ML:

Before training or auditing ML tables, refresh model feature metadata so OOD splits have the required columns. The preferred one-command training path does this automatically, writes refreshed copies under `<out-dir>/_feature_metadata/`, trains, writes performance outputs, and runs the audit unless explicitly skipped:

```
python -m analysis.cli.run_ml_training_pass \
  --out-dir data/<run_id>/ml_training_pass
```

After docking/scoring a timestamped SPD add-on, use the reusable merge pipeline
instead of manually joining CSVs:

```
python -m analysis.cli.merge_spd_addon_pipeline \
  --base-table data/AtlasSPD_phase1/<current_model_ready.csv> \
  --addon-dir data/AtlasSPD_phase1/target_positive_addon_<YYYYMMDD_HHMMSS>
```

When several agents produce timestamped add-on directories, omit `--addon-dir`
and optionally set `--addon-root data/AtlasSPD_phase1`. The command then selects
the newest completed score-ready directory, skips newer incomplete directories,
and records the selection and rejection reasons in the merge manifest. It never
selects a directory merely because its timestamp is newest.

The command discovers the non-smoke full-reference Vina manifest, records its
checksum, verifies the selected-pair checksum when the manifest declares one,
and merges on normalized PDB plus ligand-stem/base pair keys. It reconciles
against the active canonical `chemdb/data/fda_mapping_from_pdbqt.csv`, refreshes
feature metadata, and quarantines identity-blocked rows from training. Compatible
duplicate evidence is collapsed; conflicting or residual duplicate PDB-ligand
pairs and accidental SPD truth on external rows block readiness. Binding,
exposure, and combined-activity audits run by default. Use `--selected-pairs`,
`--score-table`, or `--reference-manifest` only for explicit historical
overrides. Do not use `--skip-audits` for a model-ready handoff.

For ML add-on docking/scoring, never use raw `main.py`, `-fast`/`--fast`,
`-dude`/`--dude`, or `--dud-library`. Those flags select low-exhaustiveness or
DUD-only workflows and are not valid add-on shortcuts. Use the exact staged
pair/library manifest and a completed full-accuracy comparison run:

```
atlas ml score-addons \
  --run-id <ADDON_RUN_ID> \
  --pairs data/<run_id>/target_positive_additions/selected_pairs.csv \
  --compare-run <FULL_REFERENCE_RUN_ID> \
  --workers 8
```

The pair manifest is the add-on library selector; do not substitute a DUD flag.
This command copies receptor/grid/Vina settings from the comparison run,
requires reference exhaustiveness >= 2, caps workers at 32, acquires an
exclusive add-on lock, and refuses to start while another Atlas run has process,
Slurm, or recent manifest evidence. Recent active manifests fail closed for 72
hours; older status-only records are treated as stale unless execution evidence
still exists. Run `atlas runs` first and resolve/finalize active records with
`atlas status <RUN_ID> --errors --explain`. Do not launch a second docking or
add-on scoring command while the first is active.


Every model training/update run must write a compact model-run ledger row. `train_ml_model` writes `model_run_record.json` and `model_run_record.csv` automatically. For a standalone or backfilled record, run:

```
python -m analysis.cli.write_model_run_ledger \
  --model-dir <model_dir> \
  --ledger data/<run_id>/ml_model_run_ledger.csv \
  --run-id <run_id>
```

For a single config-driven training update, use:

```
python -m analysis.cli.run_model_update \
  --config docs/examples/ml_model_update.example.yaml \
  --dataset <model_ready.csv> \
  --out-dir data/<run_id>/ml_model_updates/<task>
```

The ledger must include run_id, date, git commit, dataset version, label, feature set, excluded columns, split method, PU strategy, model type, hyperparameters, calibration method, row/positive counts, AUROC, PR-AUC, precision@K, enrichment@K, Brier score, notes, and output path. Treat missing ledger rows as an incomplete ML handoff.
For model-family or modular-vs-flat comparisons, use strict feature-set validation (`--strict-feature-set`, or `run_model_update` default strict mode). Do not compare models when a named feature set silently falls back to a partial subset; either build the complete model-ready table or mark the run as exploratory with `--allow-missing-feature-set`.

For BANANA-backed feature sets, the one-command training path can populate BANANA scores automatically for labelable rows before training. This keeps strict feature-set comparisons practical without launching a full FDA-matrix BANANA run by surprise. Use the default unless an explicit full-matrix score pass is intended:

```
python -m analysis.cli.run_ml_training_pass \
  --out-dir data/<run_id>/ml_training_pass \
  --banana-scoring auto \
  --banana-score-scope labelable
```

Use `--banana-score-scope all` only for an explicit large scoring job, and report runtime/coverage separately. Missing BANANA scores remain missing; do not impute them as negative binding evidence.

For the full source/PU policy sweep, use:

```
python -m analysis.cli.run_ml_source_pu_suite \
  --out-dir data/<run_id>/ml_source_pu_suite
```

This runs within-source ChEMBL/Papyrus, ToxCast, and SPD experts separately, treats ChEMBL/Papyrus -> ToxCast as a hard transfer stress test, runs matched-domain transfer where available, and runs PU sensitivity modes (`standard_binary`, `bagging_pu`, `stratified_bagging_pu`, `propensity_weighted_pu`, `elkan_noto`, `pulsnar_style`) when unlabeled rows are retained.

For standalone metadata refreshes, use:

```
python -m analysis.cli.refresh_ml_feature_metadata \
  --run-dir data/<run_id>
```

Useful refresh flags include:

```
--chemical-cluster auto|scaffold|smiles|chemotype|none
--target-family auto|protein_class|gene_heuristic|none
--source-lineage auto|none
--drop-column <col> ...
--dataset <csv> ...
--in-place
```

After every ML data ingestion, label-collapse, feature-table rebuild, metadata refresh, or model-training pass, run the modular leakage and performance audit suite before interpreting metrics:

```
python -m analysis.cli.run_ml_audit_suite \
  --dataset <model_ready.csv> \
  --label <label_col> \
  --feature-set <feature_set> \
  --splits random drug_holdout target_holdout scaffold_holdout chemical_cluster_holdout target_family_holdout temporal_holdout source_holdout \
  --source-col label_source \
  --out-dir <audit_out_dir>
```

Use source holdout only when at least two distinct label sources are present. Treat temporal holdout as required for prospective claims when year metadata exists; skipped temporal checks mean no prospective claim. Treat skipped OOD splits as missing metadata to fix unless the dataset truly lacks that axis. Treat source-transfer failures, large prevalence shifts, high-risk feature scans, feature missingness shifts, applicability-domain failures, weak calibration, missing conformal outputs, decoy-bias flags, and poor baseline-panel performance as publication caveats or blockers rather than silent warnings. Every reported model should include `model_claim_readiness.json`, `standard_baseline_panel.csv`, `model_grouped_calibration.csv`, `model_applicability_domain_summary.csv`, `pu_training_manifest.json`, `decoy_bias_summary.json`, and conformal outputs when a calibration split is available.

PU/negative handling: use `--pu-mode bagging_pu` only when unlabeled rows are present and should be sampled as temporary training negatives inside each fold. Use `--pu-mode elkan_noto` as a sensitivity baseline only, and report `pu_elkan_noto_manifest.json`, `pu_elkan_noto_c_by_stratum.csv`, and `pu_sar_selection_diagnostics.json`; if c varies by stratum or SAR diagnostics separate labeled positives from unlabeled rows, treat SCAR as violated and prefer bagged/SAR-aware PU claims. Never pre-label unknowns globally before splitting for claim-grade evaluation. Treat `negative_confidence`, `negative_evidence_type`, and `_sample_weight` as training-weight provenance, not truth.

4. DOCUMENTATION (Handoffs)
File: docs/agent_handoffs.md (Append-only).
Trigger: Handoff, architectural change, blocker, or output interface change.
Format (Strict):

Date: YYYY-MM-DD | Task: <ID/Title> | Owner: <Agent/User>
Change: <What changed/blocked>
Why: <Rationale/Root Cause>
Impact: <Risks/Compatibility>
Next: <Action Item>

5. QA & COMPLETION LOOP
Before stopping:

Review: Check for bugs, edge cases, and simplification opportunities.

Refactor: Improve readability/conciseness.

Verify: Re-run main.py --test -fast (Smoke Test).
Prefer `atlas smoke internal` for the same legacy-compatible smoke unless the task specifically requires invoking `main.py` directly.

6. TOOLBOX
Rescore: python -m post_docking.rescoring.rescoring_scorch --run-id <ID> --overwrite
Environment: micromamba docking-env ... can be activated with micromamba activate docking-env
Test Pathing: pytest chemdb/tests/test_path_router.py -q
Public smoke: atlas smoke public
Internal smoke: atlas smoke internal
First run: atlas first-run --uniprot P04637 --ligands fda --small-library --fast --dry-run
Report: atlas analysis report <RUN_ID> --status-html
Reproduce: atlas reproduce bundle <RUN_ID>
Slurm: atlas slurm submit --run-id <ID> --array 0-31 --cpus-per-task 8 --dry-run
