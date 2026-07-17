# Atlas CLI Reference

Run:

```bash
atlas --help
```

The help output is organized into these sections:

- Inputs
- Run control
- Docking modes
- Ligand/library modes
- pH/APO/HOLO options
- Reporting
- Debug/developer options

Common commands:

```bash
atlas init
atlas init --write-example-inputs
atlas demo --status-html
atlas smoke public
atlas setup-report
atlas new-run --pdb 1BN1 --small-library --fast --dry-run
atlas first-run --uniprot P04637 --ligands fda --small-library --fast --dry-run
atlas status atlas_demo
atlas status <RUN_ID> --errors --explain
atlas status <RUN_ID> --watch 30
atlas status <RUN_ID> --html
atlas runs
atlas report <RUN_ID>
atlas analysis report <RUN_ID> --status-html
atlas analysis dud-eval <RUN_ID>
atlas analysis throughput integrity <RUN_ID>
atlas analysis throughput acceptance <BASELINE_RUN_ID> <CANDIDATE_RUN_ID>
atlas analysis interactions export <RUN_ID>
atlas analysis dataset-eda --dataset <CSV_OR_PARQUET> --label-col <LABEL>
atlas analysis dataset-eda --dataset-name AtlasSPD_phase1 --list-matches
atlas ml prepare-run --run-id <RUN_ID> --sample-rows 5000
atlas ml train --run-id <RUN_ID>
atlas ml train --run-id <RUN_ID> --sample-rows 5000
atlas ml doctor --run-id <RUN_ID> --strict --write-report
atlas ml splits lock --run-id <RUN_ID>
atlas ml splits validate --run-id <RUN_ID>
atlas ml audit --run-id <RUN_ID> --dataset outputs/data/<RUN_ID>/ml/training_pass/<expert>/model/model_predictions.csv --label <LABEL> --splits random
atlas ml source-pu --run-id <RUN_ID>
atlas ml hpo --run-id <RUN_ID> --dataset <CSV> --label <LABEL> --trials 25
atlas ml leaderboard --run-id <RUN_ID>
atlas ml external-eval --run-id <RUN_ID> --model-dir <MODEL_DIR> --dataset <CSV> --label <LABEL> --name <BENCHMARK>
atlas ml chemprop --run-id <RUN_ID> --dataset <CSV> --label <LABEL>
atlas ml tdc --run-id <RUN_ID> --name hERG --tdc-module single_pred --tdc-class Tox
atlas ml reinvent --run-id <RUN_ID> --generated-smiles generated.csv
atlas ml gen-bench --run-id <RUN_ID> --suite guacamol --generated-smiles generated.csv
atlas ml status --run-id <RUN_ID>
atlas throughput bench --profile micro --mode local
atlas throughput bench --profile smoke --mode local
atlas throughput bench --profile smoke --mode slurm-sim --array 0-2%3 --cpus-per-task 4
atlas throughput bench --profile medium --mode slurm-dry-run --array 0-31%4 --cpus-per-task 8
atlas reproduce bundle <RUN_ID>
atlas debug <RUN_ID> --deep
atlas artifacts measure --run-id <RUN_ID>
atlas slurm submit --run-id <RUN_ID> --array 0-31%4 --cpus-per-task 8 --dry-run
atlas slurm progress <RUN_ID>
atlas slurm finalize <RUN_ID> --reconcile-only
atlas screenshot <RUN_ID> --pdb 1BN1 --top 20
atlas screenshot <RUN_ID> --pdb 1BN1 --representatives --gallery-html
atlas screenshot <RUN_ID> --pdb 1BN1 --ligand imatinib --without-ligand
atlas targets panels
atlas targets species
atlas targets guide
atlas targets genes --genes EGFR ABL1 --gene-type KINASE --out analysis/gene_list/kinase_genes.csv
atlas targets from-genes EGFR ABL1 --out analysis/gene_list/kinase_targets.csv
atlas targets search "BRCA DNA repair" --out analysis/gene_list/brca_topic_targets.csv
atlas targets install --panel kinases --max-return 3
atlas targets install BRCA1 BRCA2
atlas targets install Trp53 --species mouse
atlas targets install --uniprot P04637
atlas targets install EGFR --ligand ATP --dedupe sequence-identity 90 --quality publication
atlas targets install EGFR --no-dedupe
atlas targets install --genes EGFR ABL1 --species "Homo sapiens" --taxonomy-id 9606 --methods X-RAY --resolution-max 2.8 --gene-type KINASE
atlas run-panel kinases --ligands chembl --limit-ligands 1000 --fast
atlas ligands install chembl --limit 1000 --dry-run
atlas ligands audit fda --mapping-csv chemdb/data/fda_mapping_from_pdbqt.csv
atlas ligands audit fda --source-sdf <fda.sdf> --approval-manifest <fda_filter_manifest.json> --library-dir <prepared_fda_pdbqts>
atlas ligands audit stereo --library-id fda --csv <ligands.csv> --csv-id-column rdk_id --csv-parent-id-column canonical_parent_id --csv-structure-column smiles --out-dir <audit_dir>
atlas doctor --pdb 1BN1 --ligands fda
atlas dev verify --full --fix --smoke
atlas --pdb 1BN1 --fast
atlas --pdbs "1BN1,2OJ9" --single imatinib
atlas --run-id pilotstudy --pdb TEST --fast
atlas --verify-tools
atlas --doctor
atlas --print-effective-config
```

Key behavior:

- `atlas init` creates a local config and performs best-effort tool path detection.
- `atlas init --write-example-inputs` also writes lightweight input setup notes under `docs/examples/`.
- `atlas demo` writes a no-docking report fixture under `outputs/data/atlas_demo/`; `--status-html` also writes the static status dashboard.
- `atlas smoke public` runs the publication-safe smoke fixture used by the installer.
- `atlas setup-report` summarizes install readiness, detected tool paths, and optional/BYOL gaps.
- `atlas new-run --pdb <ID> ... --dry-run` validates first-run inputs and prints a launchable command.
- `atlas first-run --uniprot <ACCESSION> --ligands <library> ... --dry-run` selects a first target, plans ligand-library prep, and prints the launch command.
- `atlas status [RUN_ID]` summarizes manifests, stage progress, ETA, Slurm state when available, and report/output paths.
- `atlas status RUN_ID --errors --explain` adds structured root-cause summaries and plain-language next actions.
- `atlas status RUN_ID --watch 30` refreshes the dashboard every 30 seconds.
- `atlas status RUN_ID --html` writes `outputs/data/<RUN_ID>/status.html` for new runs, while still reading legacy `data/<RUN_ID>/` outputs when present.
- `atlas runs` lists recent run IDs, progress, failures, age, and report availability.
- `atlas report RUN_ID` exports `outputs/data/<RUN_ID>/master_rows.csv` if needed, then runs the maintained report generator for a completed run.
- `atlas analysis report RUN_ID --status-html` is the discoverable report workflow: master export, report HTML/YAML, and optional status HTML.
- `atlas analysis dud-eval RUN_ID` runs the maintained DUD/decoy benchmark evaluator using run-scoped output roots.
- `atlas analysis holo-integrity RUN_ID --strict` certifies the post-prep receptor-state partition and writes `outputs/data/<RUN_ID>/holo_integrity/holo_integrity.csv` plus `holo_integrity.json`. For the frozen 93-target SPD panel, the expected evidence partition is 76 strict holo targets and 17 apo fallbacks; an apo fallback remains apo and is never relabeled as experimentally holo.
- `atlas analysis throughput integrity RUN_ID` writes throughput integrity JSON/CSV under `outputs/data/<RUN_ID>/`.
- `atlas analysis throughput acceptance BASELINE_RUN_ID CANDIDATE_RUN_ID` compares throughput KPIs for two runs.
- `atlas analysis interactions export RUN_ID` converts run reporting rows to the maintained partitioned interaction dataset.
- `atlas analysis dataset-eda` runs read-only tabular EDA on CSV/TSV/Parquet/JSONL/JSON datasets and writes summaries, optional ydata-profiling HTML/JSON, Phi_K and dython association matrices, scikit-learn mutual information rankings when `--label-col` is provided, and NetworkX/PyVis association-network outputs under `outputs/data/dataset_eda/` by default. Use `--dataset-name ... --list-matches` to resolve fuzzy names before reading a dataset.
- `atlas ml prepare-run` and `atlas ml ready` run ML doctor, metadata-aligned split locking, deterministic sampled training, and leaderboard refresh for an imported run.
- `atlas ml train RUN_ID` wraps the one-command four-expert ML pass, refreshes/caches feature metadata, writes models under `outputs/data/<RUN_ID>/ml/training_pass/`, and records dataset/model/artifact manifests. Use `--sample-rows` for deterministic fast iteration and `--split-manifest` or `--expert-split-manifest expert=PATH` to reuse locked splits.
- `atlas ml doctor` preflights ML labels, source/scaffold/target-family/temporal metadata, source imbalance, tiny class counts, unlabeled-row provenance, and feature leakage; `--strict` turns warnings into blockers and `--write-report` writes `ml/readiness/ml_doctor_report.json`.
- `atlas ml splits lock|validate` writes and verifies metadata-refreshed split manifests under `outputs/data/<RUN_ID>/ml/splits/`.
- `atlas ml audit` runs the modular leakage/performance audit suite; use `--claim-mode publication` to block missing OOD/source/temporal/calibration/conformal evidence.
- `atlas ml source-pu` runs source-transfer and PU sensitivity workflows under `outputs/data/<RUN_ID>/ml/source_pu_suite/`.
- `atlas ml hpo` runs the optional Optuna sweep wrapper and writes `outputs/data/<RUN_ID>/ml/optuna_sweep/` by default; it also accepts `--split-manifest` and `--pruner median|successive_halving|hyperband|none`.
- `atlas ml leaderboard` writes `ml/leaderboard/model_leaderboard.csv` and `.json` across model, HPO, adapter, and external-eval manifests.
- `atlas ml external-eval` scores a trained Atlas model on an external labeled table and writes predictions, metrics, reliability, subgroup metrics, and `external_eval_manifest.json`.
- `atlas ml chemprop` stages an Atlas model-ready CSV for Chemprop and can launch Chemprop when `--run` and tool paths are available.
- `atlas ml tdc` stages PyTDC or local TDC-style CSV benchmarks under `outputs/data/<RUN_ID>/ml/tdc/` by default.
- `atlas ml reinvent` captures or explicitly launches REINVENT4 generation outputs under `outputs/data/<RUN_ID>/ml/reinvent_generation/`.
- `atlas ml gen-bench` prepares and captures GuacaMol/MOSES benchmark artifacts under run-scoped ML output directories.
- `atlas ml status` summarizes run-scoped ML manifests and trained/skipped/failed expert statuses.
- `python -m analysis.cli.run_reinvent_generation` captures REINVENT4 generated SMILES/CSV/JSON outputs into `generated_smiles.csv` plus `reinvent_generation_manifest.json`; external execution requires explicit `--execute --command ...`.
- `python -m analysis.cli.benchmark_generated_smiles` prepares GuacaMol/MOSES one-SMILES-per-line inputs and captures benchmark result manifests/metric summaries without importing those packages.
- `atlas throughput bench ...` runs the throughput harness without changing production protein-prep policy. It writes `bench_config.json`, local/simulated Slurm summaries, Slurm submit previews, `prep_cache_summary.json`, `slurm_readiness.json`, optional util-bench output, and optional acceptance reports under `outputs/data/<RUN_ID>/throughput_bench/`. The `micro` profile is for fast scheduler/reporting development and requires cached receptor prep unless `--allow-fresh-prep` is used; use `smoke`, `medium`, and real Slurm runs for performance evidence.
- `atlas reproduce bundle RUN_ID` writes `outputs/data/<RUN_ID>/reproducibility_bundle.json` and status JSON beside the run report without creating a separate top-level output tree.
- `atlas debug RUN_ID --deep` wraps the structured status/error summarizer and avoids raw log inspection.
- `atlas artifacts measure --run-id RUN_ID` measures run-scoped outputs and writes reports under `outputs/data/<RUN_ID>/artifact_reports/` by default.
- `atlas slurm submit ...` submits the maintained Slurm worker/finalizer scripts with CLI overrides; choose `--array` concurrency and `--cpus-per-task` for the actual cluster allocation and queue policy. Use `atlas slurm submit --bench2-canary --run-id <RUN_ID> --dry-run` to preview the low-SU bench2-fast utilization canary; it requests three exclusive whole-node array tasks for 15 minutes and uses node CPU auto-detection unless you explicitly pass `--cpus-per-task`.
- `atlas slurm progress RUN_ID` wraps the status dashboard with live Slurm probes.
- `atlas slurm finalize RUN_ID` wraps `tools/finalize_distributed_run.py`.
- `atlas screenshot RUN_ID --pdb <ID> --top 20` restores selected archived poses when needed and writes active-site PNGs under `outputs/data/<RUN_ID>/screenshots/`. The default view shows the full protein as a cartoon with the ligand in sticks; use `--view-context pocket` for a zoomed active-site view or `--protein-style surface`, `mesh`, `lines`, or `none` for alternate receptor context.
- `atlas screenshot RUN_ID --pdb <ID> --representatives --gallery-html` renders representative best, median, and worst scored ligands and writes a self-contained `gallery.html` with embedded PNGs for run-level visual triage.
- `atlas targets guide` prompts for a built-in panel or custom genes, organism/species filters, resolution, and whether to write only the gene CSV, query candidates, dry-run install, or download PDBs.
- `atlas targets genes ...` writes a reusable `gene,category` CSV for custom genes or a built-in panel.
- `atlas targets from-genes ...` wraps the RCSB gene-to-PDB selector for wetlab-first target selection.
- `atlas targets search "BRCA DNA repair"` runs a separate RCSB full-text topic search and ranks the matching structures with the same docking-readiness filters.
- `atlas targets panels` lists built-in gene panels such as `kinases`, `gpcrs`, `ion-channels`, `nuclear-receptors`, and `adme`.
- `atlas targets species` lists organism aliases such as `human`, `mouse`, `rat`, `zebrafish`, `fly`, `worm`, `yeast`, and `ecoli`.
- `atlas targets install BRCA1 BRCA2` accepts gene symbols directly; `--genes`, `--genes-file`, `--uniprot`, and `--panel` remain available for scripted workflows.
- `atlas targets install ...` selects candidate PDB IDs, downloads PDB files into the configured input root, and writes a target install manifest.
- Target queries de-duplicate results at 90% sequence identity by default. Use `--dedupe sequence-identity 70` for stricter diversity or `--no-dedupe` when all near-identical structures are needed.
- `--ligand ATP` / `--contains-ligand HEM` and `--quality publication` can be used on target query/install commands when a run needs co-crystal ligand or publication-style quality filters.
- Target gene search defaults to human proteins, all experimental methods, up to 10 PDBs per gene, ligand-bound structures, and the resolution cutoff for entries where RCSB reports resolution. Use `--species mouse`, `--organism zebrafish`, `--taxonomy-id`, `--entity-type`, `--methods`, `--resolution-max`, `--gene-type`, and `--allow-apo` for custom target sets.
- `atlas run-panel <panel> --ligands <library>` installs panel targets, installs/preps the ligand library, and runs Atlas with those selected PDBs.
- `atlas ligands install <library> --dry-run` prints the planned cache and prepared-library paths without fetching or preparing files.
- `atlas ligands audit stereo --library-id <ID> --sdf <SOURCE.sdf> --out-dir <DIR>` audits source-defined and unresolved potential stereochemistry without changing or assigning structures. SDF input uses `_Name` as its ligand ID by default; use `--sdf-id-property` and `--sdf-parent-id-property` when the source identifiers are stored in other properties. CSV input requires `--csv-id-column` and `--csv-structure-column`, accepts SMILES, InChI, or MolBlock structures through `--csv-structure-format`, and can retain a separate canonical-parent key with `--csv-parent-id-column`.
- The stereo audit writes a row ledger, JSONL evidence, and a summary with source-file/record hashes, RDKit method/version, specified and potential-but-unassigned element counts, parse status, and unresolved status. Optional `--csv-prepared-path-column` and `--csv-prepared-sha256-column` bind and verify prepared bytes, but PDBQT is never accepted as stereochemical source evidence. `--fail-on-parse-error`, `--fail-on-unspecified`, and `--fail-on-unresolved` provide automation gates. These are audit gates only: unspecified or unknown stereo is not automatically assigned, interpreted as invalid, or allowed to qualify a rank/native redock until the scientific policy is approved.
- `atlas ligands audit fda` performs an offline, non-mutating structure/name concordance and repair audit. In addition to row evidence, review/quarantine tables, an active-parent eligibility manifest, the source catalog, `fda_approval_cross_reference.csv`, and the JSON summary, it writes four repair artifacts: `fda_mapping_repaired.csv` preserves one ledger row per original mapping with original identity fields, corrected manifest-backed DrugCentral identity, prepared-file path and byte SHA256; `fda_score_reuse_manifest.csv` records explicit historical-score joins and their claim-grade provenance status; `fda_redock_delta.csv` contains one preparation/docking request per missing unique approved parent; and `fda_repair_quarantine.csv` records unresolved or unsafe-to-collapse mapping/source records. Outputs default to `outputs/data/fda_identity_audit/`; inputs and existing score tables are never rewritten.
- The cross-reference assigns every mapping row to exactly one of seven categories: confirmed FDA active ingredient, confirmed FDA salt whose parent was docked, retained salt/counterion, collapsed combination product, additive/excipient, incorrect name–structure mapping, or unverifiable. Pass the exact-filtered DrugCentral FDA SDF with `--source-sdf`, its hash-bearing `--approval-manifest` (auto-discovered as `fda_filter_manifest.json` beside the SDF), and the corresponding prepared library with `--library-dir`. FDA confirmation requires the validated manifest, a stable exact or unique full-parent InChIKey/source-record link, and a corroborated name from that linked record; 14-character connectivity, ordinal position, formula, and names alone cannot establish approval. Coordinate-reconstructed connectivity can establish whether the prepared file retained the mapped full form or parent form, but it cannot establish FDA approval or docking eligibility. Mapping PubChem/DrugCentral labels remain review evidence until independently linked. `--fail-on-blockers` returns a nonzero status for probable identity mismatches.
- Historical score reconciliation is opt-in: repeat `--score-csv <path>` for exact score tables or `--score-run-id <RUN_ID>` to select that run's `fda_master_rows.csv` (falling back deterministically to `master_rows.csv`). The command never scans the repository or guesses a run. `chemistry_redock_required` answers only whether the approved parent lacks a compatible prepared structure; missing historical hashes do not change that chemistry decision. Claim-grade score reuse is a separate, fail-closed status requiring exact historical ligand bytes plus receptor, effective configuration, engine/input, completion, and pose/score provenance. Legacy rows lacking those hashes remain `candidate_context_unverified`. Raw per-engine scores may be candidates for reuse only when this context is verified; ranks, z-scores, FDR values, consensus scores, SCORCH values, and other cohort-derived statistics must be recomputed unless the cohort identity is also proven unchanged.
- `atlas ligands repair fda --repaired-mapping-csv <CSV> --redock-delta-csv <CSV> --source-sdf <SDF> --approval-manifest <JSON> --output-dir <DIR>` builds a separate human-readable FDA canonical-parent library without modifying legacy PDBQTs or launching docking. It copies the 1,264 compatible approved parents byte-for-byte, prepares only the unique approved-parent delta with the shared Meeko backend, and writes a per-parent manifest, chemistry/source sidecars, byte checksums, quarantine table, canonical parent SDF, summary, and library index. Filenames use sanitized preferred drug names; DrugCentral and parent-key hash suffixes are added only if names collide. The manifest, never the filename, is the identity authority.
- Use `atlas ligands repair fda ... --plan-only` to validate the complete partition, source-record parent chemistry, stereochemistry, names, checksums, and suitability policy without writing files or preparing ligands. The current manifest-backed plan contains 1,863 approved canonical parents: 1,264 compatible legacy parents and 599 delta parents. Thirteen delta rows remain in the plan but are `quarantined_non_dockable` and produce no placeholder PDBQT: six degenerate parent graphs, six nonstandard-element cases that require explicit preparation-and-engine capability evidence, and nitric oxide as a tiny carbon-free/radical parent. The finalized library at `prepped_ligands/fda_named_parent_library_20260712_final/` contains 1,850 ready PDBQTs and 13 quarantine rows with zero preparation failures. Deutetrabenazine is ready after preserving explicit isotope H/D while adding its missing implicit hydrogens; nitric oxide remains quarantined. The builder itself did not launch docking; the later separately authorized delta run is bound below. A consumable directory contains only rows with `materialization_status=ready`; approval does not imply docking suitability.
- `atlas ligands repair fda-delta --named-library-dir <FINAL_NAMED_LIBRARY> --output-dir <DIR>` creates the separately consumable 586-PDBQT `prepare_parent + ready` subset without rerunning Meeko. It validates the finalized 1,863/1,264/599/586/13 partition and every source/output checksum, prefers hardlinks with a byte-copy fallback, preserves all 599 delta rows and 13 quarantines in provenance manifests, writes a 586-entry library index, and fails on drift or changed resume artifacts. Add `--copy` to disable hardlinks.
- `atlas ligands repair fda-legacy --original-mapping-csv <CSV> --repaired-mapping-csv <CSV> --redock-delta-csv <CSV> --named-library-manifest-csv <CSV> --output-dir <DIR>` writes a versioned reconciliation bundle without changing the canonical legacy map, PDBQTs, or score tables. It emits a full-schema mapping v2, exact identity-change ledger, byte-verified file inventory, score identity/checksum joins, unresolved second-pass queue, canonical-parent redock linkage, and a hash-bearing summary. Repeat `--score-csv <path>` to name score inputs explicitly. If omitted, the command uses only `SPD_addons/fda_master_rows.csv`, `SPDaddon2/fda_master_rows.csv`, and `SPDsparsefloor/fda_master_rows.csv`; it never recursively discovers or selects newer files.
- `atlas ligands repair fda-terminal-v3` resolves every v2 legacy row structure-first against the full DrugCentral SDF, validated CID and exact-InChIKey PubChem caches, Drugs@FDA ingredients, repaired-map connectivity evidence, the approved-source SDF partition, the named canonical manifest, and the explicit repair quarantine. All inputs are required and hash-bound. The command writes mapping v3, terminal dispositions, a ready FDA subset, an excluded/non-ready subset, and a completion/hash summary without mutating inputs. It fails completion when any usable row lacks a structure-validated preferred identity, any FDA row lacks a canonical or explicit unsafe-quarantine link, any PDBQT checksum changes, or an unexpected terminal value appears. Score reuse remains fail-closed for ambiguous connectivity, missing files, and unsafe parent collapse.
- Build or verify the exact offline PubChem caches with `python tools/fetch_pubchem_cid_cache.py --input-csv <mapping-v2.csv> --cache-dir <CID_CACHE>` and `python tools/fetch_pubchem_inchikey_cache.py --input-csv <mapping-v2.csv> --cache-dir <INCHIKEY_CACHE>`; add `--offline` to forbid network access and require verified cached content. Before promotion, run the independent bundle verifier explicitly:

  ```bash
  python tools/verify_fda_terminal_bundle.py \
    --mapping <fda_legacy_mapping_v3.csv> \
    --terminal-dispositions <fda_legacy_terminal_dispositions_v3.csv> \
    --named-manifest <fda_named_library_manifest.csv> \
    --pubchem-manifest <pubchem_cache_manifest.json> \
    --repair-quarantine <fda_repair_quarantine.csv> \
    --fda-source-sdf <fda_source.sdf> \
    --output <fda_terminal_bundle_verification.json>
  ```

- `atlas ligands promote fda-v3 --mapping-v3-csv <CSV> --verifier-json <JSON> --output-dir <DIR>` promotes only a terminal-resolver v3 mapping whose SHA-256 matches a passing independent `atlas.fda-terminal-bundle-verification.v1` report with all checks passed. It writes the immutable, hash-bearing `fda_mapping_v3_promotion_manifest.json` as a pointer and never copies the mapping or changes a default. Repeat `--score-csv <CSV>` to write `fda_mapping_v3_score_identity_joins.csv`; identity fields are attached only when both the exact RDK ID and a historical/input PDBQT checksum match, and every score row is stamped with the mapping SHA-256. Accepted score evidence fields are `ligand_pdbqt_sha256`, `pdbqt_sha256`, `ligand_sha256`, `input_ligand_sha256`, and `historical_pdbqt_sha256`; score-side `current_pdbqt_sha256` and `selected_pdbqt_sha256` are current context, not historical input proof, and cannot authorize a join. The join propagates `legacy_score_reuse_status` and `prepared_connectivity_state`, but `score_reuse_authorized` is true only for the explicit matching coordinate-exact or coordinate-parent status/state pairs; ambiguous connectivity, missing files, unsafe parent collapse (including the bismuth override), missing fields, and inconsistent pairs fail closed without suppressing an otherwise exact identity attachment. `--update-config <CONFIG_PATH>` is the only config mutation path: it atomically updates `FDA_MAPPING_CSV` in that explicitly named existing file and records pre/post hashes; omitting it leaves every config untouched. The verified terminal-v3 mapping currently used by defaults and `config.txt` is `chemdb/data/fda_mapping_from_pdbqt.csv`; the hash-bound source copy remains under `outputs/data/fda_repair_20260712/terminal_reconciliation_v3/`.
- The finalized bundle under `outputs/data/fda_repair_20260712/legacy_reconciliation_v2/` preserves 4,729 usable legacy files and one unusable row, applies exactly 945 manifest-backed identity corrections in the new mapping, retains 3,459 mapping rows as unresolved, and preserves all 342 input score rows. Of those score rows, 228 join by `rdk_id` with a verified current checksum attached but lack a historical ligand checksum, so they remain `candidate_context_unverified`; 114 have unparseable/non-RDK ligand identifiers and remain unresolved rather than being labeled non-FDA. The 599-parent linkage resolves to 586 prepared-ready and 13 quarantined parents. Its `docking_launched=false` values describe the immutable pre-launch plan state; they are intentionally not rewritten when a later run starts. Cohort-dependent ranks, FDR, consensus, and SCORCH statistics still require recomputation.
- `atlas ligands repair fda-legacy-bind-run --reconciliation-dir <DIR> --run-id <RUN_ID> --target-dir <DIR> --target-manifest-csv <CSV> --delta-summary-json <JSON> --run-manifest-yaml <YAML> --cpu-count <N>` writes only `fda_legacy_run_binding_v2.json` beside the reconciliation bundle. It verifies the exact ready-parent/hash set, every current delta PDBQT checksum, the frozen target selection and PDB bytes, run ID/status/library/holo mode, and point-in-time run-manifest hash. For `postFDAfixSPDrun`, the artifact binds 586 ready parents to the 93-target manifest and the active holo run with declared CPU count 16; `--observed-worker-concurrency 16` records the separately observed process concurrency. The run manifest remains mutable while active, so its SHA-256 is explicitly a capture-time snapshot, not a final-run hash.
- If a supplied PDBQT lacks a SMILES/InChIKey remark and authoritative source sidecar, the audit reconstructs element connectivity from coordinates, with Open Babel as a fallback for files RDKit cannot parse. It reports `connectivity_consistent_exact`, `connectivity_consistent_parent`, `connectivity_disagrees`, or `connectivity_unavailable` separately from high-fidelity bond-order/stereo evidence. Connectivity-only agreement is useful for map-to-file reconciliation but can never establish FDA approval or docking eligibility.
- FDA preparation also validates that an existing PDBQT library was built from the current filtered SDF and filter manifest. A stale or legacy unfiltered library fails closed; rerun with `atlas ligands prep fda --force-prep` to preserve the old directory under an `.unvalidated` name and rebuild it.
- `atlas doctor --pdb <ID> --ligands <library>` validates the prepared receptor/ligand PDBQT inputs, detects multi-model or empty ligand files, and reports the prep command to rerun when conversion failed.
- `atlas stages repeat --run-id <RUN_ID> --stage pose-validation --all-receptors --all-ligands --all-contexts --fraction 1 --database <release.sqlite>` writes an auditable planner-only ledger for one stage. Supported stages are `vina`, `gnina`, `scorch`, `mmgbsa`, and `pose-validation`; no workload launches in this release. Pose-consuming plans reject unknown score-to-pose semantics, and aggregate scores require the complete declared set of verified contributing-pose hashes.
- Partial `--fraction` or `--count` plans require `--selection-strategy hash|top-score`. Deterministic `hash` is implemented. `top-score` fails closed until a score source, direction, provenance, and tie policy are approved. SQLite sources above 25,000 pair cells also fail closed pending the streaming planner and exact-pair executors.
- `atlas dev verify --full --fix --smoke` runs the canonical repo gate: full ruff, vulture, tach, architecture gate, full mypy, and smoke.
- `--pdb` is the preferred repeatable input selector.
- `--run-id` names the run and overrides `ATLAS_RUN_ID`.
- `--single` restricts docking to one ligand selector.
- `--ph-ligand-mode` controls pH-aware ligand selection when pH ensemble modes are active.
- `--verify-tools` performs detailed external-tool checks and exits.
- `atlas doctor` or `--doctor` performs lightweight Python/environment/tool-resolution diagnostics and exits.

For a complete fresh-install path, see `docs/quickstart_new_users.md`.

### ML data-gap report

Use this after label ingestion, feature refresh, or a model-training pass to identify which target, family, source, scaffold, and chemistry strata block claim-grade ML generalization:

```bash
python -m analysis.cli.run_ml_data_gap_report \
  --dataset data/<run_id>/model_ready/<table>.csv \
  --out-dir data/<run_id>/evaluation_ready/data_gap_report
```

Primary outputs are `ml_data_gap_report.csv`, `holdout_claim_readiness.csv`, `source_calibration_readiness.csv`, and `data_collection_priorities.csv`. Treat provenance-derived rates as audit/calibration context, not clean predictive features.

## AtlasSPD Phase 1 PK context

Refresh structured exposure context and join it to the frozen Phase 1 table:

```bash
python -m analysis.cli.refresh_phase1_pk_context \
  --out-dir data/AtlasSPD_phase1/pk_context_v0_0_15
```

The command writes `pk_context_long.csv`, `pk_context_representative.csv`,
`AtlasSPD_phase1_pk_enriched.csv`, endpoint-specific clearance and absolute-
bioavailability views, source coverage, and a fail-closed
`pk_context_validation.json`. Primary dose fields are populated only when the
dose is bound to the selected Cmax scenario. Maximum labeled recommended adult
dose is a separately named sensitivity artifact; highest-studied and
maximum-tolerated doses never enter the primary context. Existing
`spd_exposure_label` values are never recomputed from external PK.

DailyMed/openFDA numeric extraction is cached. Cmax rows enter training only
through the source-text-curated context registry. Clearance, absolute
bioavailability, and maximum-dose candidates require same-clause semantic gates
and an explicit canonical-candidate decision file:

```bash
python -m analysis.cli.adjudicate_spl_pk_candidates \
  --candidates <spl_pk_candidate_review.csv> \
  --out-dir <spl_adjudication_dir> \
  --review-decisions <spl_pk_review_decisions.csv>
```

Machine acceptance alone is not training approval. Apparent oral clearance,
plasma/systemic clearance, renal clearance, and other endpoint types remain
separate columns and are never pooled. Machine-ambiguous absolute
bioavailability can enter only through a structured review that confirms the
value, qualifier, analyte, extravascular route, endpoint, and absolute-reference
basis. Relative-only evidence stays excluded; the review decision and source
hashes remain attached to every admitted row.

When PK-DB's normalized output endpoint is unavailable, recover contextual rows
from its documented study endpoint and public study source TSVs:

```bash
python -m analysis.cli.recover_pkdb_context \
  --dataset <model_ready.csv> \
  --out-dir data/external/pkdb/recovered_phase1
```

Public PK-DB access does not grant original-source ML rights. Recovered rows stay
quarantined unless `--source-rights-manifest` supplies an affirmative per-study
`training_allowed` decision and a rights reference; the same manifest can be
passed to the one-command refresh as `--pkdb-source-rights-manifest`. Use
`--drugbank-cmax` and `--drugbank-protein-binding` only with licensed local
exports. Apply for academic DrugBank access at
https://go.drugbank.com/academic_research and verify current release availability
at https://go.drugbank.com/releases/latest.

### Pair-residual binding validation

Use grouped OOF validation to separate ligand propensity from pair-specific evidence and optionally test model-ready pocket descriptors:

```bash
python -m analysis.cli.run_pair_residual_binding \
  --dataset data/<run_id>/model_ready.csv \
  --out-dir data/<run_id>/pair_residual_grouped_oof \
  --ligand-feature-set ligand_primary_functional_branch_no_qed \
  --pair-feature final_score \
  --pair-feature banana_score_normalized \
  --protein-feature pocket_atom_count \
  --protein-feature pocket_residue_count \
  --shortcut-feature rdkit_mol_wt \
  --variant pair-only \
  --variant ligand-only \
  --variant combined \
  --variant protein-only \
  --variant pair-protein \
  --variant combined-all \
  --variant shortcut-reduced-all \
  --outer-group drug_id \
  --outer-group target_id \
  --outer-group drug_id+target_id \
  --outer-group drug_id+scaffold_key
```

With the explicit `--variant` options shown, the output compares `protein_only`, `pair_protein`, `combined_all`, and `shortcut_reduced_all` with pair-only, ligand-only, and combined variants on one complete-case cohort. `pair_residual_pooled_metrics.csv` includes PR-AUC, AUROC, Brier skill against the prevalence-only predictor, adaptive ECE, and calibration intercept/slope. These sigmoid-transformed OOF logits remain diagnostic scores until an independent grouped calibration and final holdout are used.
