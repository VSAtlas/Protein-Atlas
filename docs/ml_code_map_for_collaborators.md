# Atlas ML Code Map

This note is for collaborators who want to inspect the current Atlas ML implementation without reading the whole repository.

## Current one-command run

```bash
atlas ml train --run-id pilotstudy
```

Default experts are defined in `analysis/cli/run_ml_training_pass.py`. The `atlas ml` wrapper defaults outputs to `outputs/data/<RUN_ID>/ml/`.

| Expert | Label | Default feature set | Current role |
|---|---|---|---|
| Clean Binding expert | `spd_binding_label` | `spd_binding_nonleaky` | Uses `consensus_score`, `banana_score_normalized`, structure/class fields, and RDKit descriptors. Aggregate priors such as `binding_expert_score` are excluded. |
| Frozen Binding prior / baseline | none | `frozen_binding_prior_baseline` | Rule-based Atlas/BANANA/SCORCH-style prior. Outputs `banana_atlas_blend_score` and `binding_expert_score`; use as a baseline or compact downstream prior, not as clean Binding expert input. |
| Direct SPD Exposure expert | `spd_exposure_label` | `spd_exposure_nonleaky` | Uses `binding_expert_score`, structure/class fields, and RDKit descriptors. AC50, Cmax, free Cmax, and exposure margin are excluded. |
| Estimated-PK Exposure model | `spd_ac50_uM`, `free_cmax_um`, final `spd_exposure_label` | explicit `--run-estimated-pk-exposure` | Deferred sensitivity path. It is no longer run by default; out-of-fold predictions are required before downstream use. |
| Tissue expression / site scorer | none by default | `tissue_site_score_only` | Deterministic expression/site scorer using expression features only; output composite is `site_relevance_score`. |
| Pair-level Tissue/Site relevance expert | independent `tissue_site_label` if available | `spd_tissue_site_nonleaky` | Combines binding/exposure/tissue priors for pair-level tissue relevance. Do not treat as clean supervised ML if labels are generated from the same expression rules. |
| Clean Mechanism expert | `mechanism_ml_label_clean` preferred; `mechanism_ml_label` with caveats | `spd_mechanism_nonleaky` | Uses binding/exposure/tissue priors plus structure/class/chemistry. Graph/pathway/target-safety evidence is excluded. |
| Mechanism evidence / sensitivity model | `mechanism_ml_label` | `spd_mechanism_evidence_sensitivity` | Evidence-augmented sensitivity model using graph, pathway, drug-ADR, target-ADR, and target-safety evidence; requires explicit `--allow-label-definition-features`. |
| Graph/KGE mechanism model candidate | held-out graph edges or panel-specific mechanism labels | graph nodes/edges, not a tabular feature set | Optional sparse-graph model family using PyKEEN KGE or PyG/DGL heterograph models. Keep separate from tabular experts until edge-split leakage, relation leakage, and negative-edge sampling are audited. |
| Hypergraph mechanism model candidate | held-out mechanism hyperedges or panel-specific labels | hyperedge incidence graph, not a tabular feature set | Optional PyG/DGL path for drug-target-ADR/tissue/pathway tuples. Requires hyperedge split manifests and careful negative hyperedge sampling. |
| KAN / DeepADR-style model candidate | same label as selected tabular expert | same nonleaky tabular features as selected expert | Experimental nonlinear tabular candidate for sparse ADR labels. Use nested validation/HPO and compare against logistic/RF/XGBoost before making claims. |

## Main files

- `analysis/cli/run_ml_training_pass.py`: one-command four-expert runner.
- `analysis/ml/train_classifier.py`: model fitting, splits, PU modes, metrics, calibration, feature importance, split manifests.
- `analysis/ml/feature_sets.py`: named feature sets and label-definition exclusions.
- `analysis/ml/leakage_checks.py`: hard leakage guards.
- `analysis/ml/splits.py`: random, drug, target, scaffold, chemical-cluster, target-family, and temporal splits.
- `analysis/ml/source_benchmark_tables.py`: curated bioactivity, ToxCast, source-transfer, and SPD table builder.
- `analysis/cli/run_ml_audit_suite.py`: modular leakage/performance audit suite.
- `analysis/cli/stage_ml_label_sources.py`, `analysis/cli/build_four_state_evidence.py`, and `analysis/cli/build_pk_table.py`: external source cataloging, four-state label collapse, and PK/exposure context staging for ChEMBL, BindingDB, ToxCast/Tox21, PubChem BioAssay, IUPHAR/GtoPdb, and NCATS Inxight.
- `analysis/cli/run_ml_source_pu_suite.py`: source-specific and PU sensitivity suite.
- `analysis/cli/run_ml_optuna_sweep.py` and `analysis/ml/optuna_sweep.py`: optional Optuna HPO wrapper; Optuna is imported only when the sweep is invoked.
- `analysis/ml/preflight.py`, `analysis/ml/split_locking.py`, `analysis/ml/leaderboard.py`, and `analysis/ml/external_eval.py`: imported-run readiness checks, locked splits, leaderboard aggregation, and external benchmark evaluation.
- `analysis/ml/notebook_summary.py`: notebook-safe manifest, registry, and leaderboard summary helpers for `outputs/data/<RUN_ID>/ml/`.
- `analysis/ml/dataset_eda.py` and `analysis/cli/run_dataset_eda.py`: read-only first-pass EDA for model-ready or staged tabular datasets, including ydata-profiling, Phi_K/phik, dython, mutual information, and association-network artifacts when optional dependencies are installed.
- `analysis/ml/graph_model_backends.py` and `analysis/cli/audit_graph_model_backends.py`: optional PyKEEN/PyG/DGL backend availability audit for graph/KGE mechanism-model experiments.
- `analysis/ml/mechanism_graph_splits.py` and `analysis/cli/build_mechanism_graph_training_splits.py`: KGE/hypergraph-ready mechanism graph nodes, reified triples, hyperedges, edge split manifests, and negative-edge sampling policy files.

## ID leakage policy

Raw IDs are metadata only. They are used for joins, split manifests, prediction output, and holdout grouping, but not as predictors.

Hard-blocked predictive features include:

```text
drug_id
target_id
pdb_id
literature_supported_label
```

The guard lives in `analysis/ml/feature_sets.py` and `analysis/ml/leakage_checks.py`. Prediction files reattach IDs after scoring so results remain interpretable.

## Current output files

For each trained expert under `outputs/data/pilotstudy/ml/training_pass/<expert>/model/`:

- `dataset_version_manifest.json`
- `model_metrics.csv`
- `model_predictions.csv`
- `model_metric_bootstrap_ci.csv`
- `model_reliability_table.csv`
- `model_grouped_calibration.csv`
- `model_subgroup_metrics.csv`
- `model_decision_metrics.csv`
- `standard_baseline_panel.csv`
- `hpo_trials.csv`
- `hpo_manifest.json`
- `conformal_summary.json`
- `feature_importance.csv`
- `split_manifest.csv`
- `split_manifest.json`
- `model_claim_readiness.json`
- `model_card.json`
- `model_artifacts.json`
- `experiment_runs.jsonl`
- `registry.json`
- `trained_model.pkl`

Optuna sweep runs add `hpo/optuna/optuna_sweep_manifest.json`, `hpo/optuna/optuna_trials.csv`, `hpo/optuna/best_params.json`, and a final model under `model/` when `--skip-final-train` is not used. Use `--pruner median|successive_halving|hyperband|none` to record and activate an Optuna pruner. The Optuna manifest records `sampler`, `pruner`, `trial_state_counts`, and `pruning_status` so notebooks and status wrappers can distinguish completed, failed, and pruned trials.

The run-level manifest is:

```text
outputs/data/pilotstudy/ml/training_pass/ml_training_pass_manifest.json
```

## Notebook

Use:

```text
notebooks/atlas_ml_exploration.ipynb
```

It loads the current run, checks selected features, plots basic distributions, summarizes label/source balance, and compares logistic regression, random forest, gradient boosting, SVM, and KNN on the existing split manifests. For the current expert architecture table, import `display_expert_feature_architecture` from `analysis.ml.notebook_summary`.


## ML preflight, readiness, and claim modes

Use `atlas ml prepare-run --run-id <RUN_ID>` after importing a run. It runs doctor, locks metadata-refreshed splits, performs deterministic sampled training, and refreshes `ml/leaderboard/model_leaderboard.csv`. Use `atlas ml doctor --run-id <RUN_ID> --strict --write-report` before publication-oriented reruns. Training and audit commands accept `--claim-mode exploratory|publication`; publication mode blocks missing source/OOD/temporal evidence, missing calibration/conformal outputs, and missing claim metadata instead of leaving them as caveats.

MLflow logging is enabled by default with a local file store under `outputs/data/<RUN_ID>/ml/mlruns` unless `--no-mlflow` or `ATLAS_DISABLE_MLFLOW=1` is used. W&B and ClearML remain opt-in through `ATLAS_WANDB_PROJECT` and `ATLAS_CLEARML_PROJECT`. The local registry artifacts are always written. For a real MLflow model registry, set `ATLAS_MLFLOW_REGISTRY_URI`, `ATLAS_MLFLOW_REGISTER_MODEL=1`, and optionally `ATLAS_MLFLOW_REGISTERED_MODEL_NAME`; this keeps registry-grade model versioning explicit instead of accidentally creating local file-only registry state.

## Adapter wrappers

Useful run-scoped wrappers:

```bash
atlas ml prepare-run --run-id pilotstudy --sample-rows 5000
atlas ml splits lock --run-id pilotstudy
atlas ml train --run-id pilotstudy --sample-rows 5000
atlas ml hpo --run-id pilotstudy --dataset model_ready.csv --label y --trials 25
atlas ml leaderboard --run-id pilotstudy
atlas ml external-eval --run-id pilotstudy --model-dir outputs/data/pilotstudy/ml/training_pass/binding/model --dataset external.csv --label y --name benchmark
atlas ml chemprop --run-id pilotstudy --dataset model_ready.csv --label y
atlas ml tdc --run-id pilotstudy --name hERG --tdc-module single_pred --tdc-class Tox
atlas ml reinvent --run-id pilotstudy --generated-smiles generated.csv
atlas ml gen-bench --run-id pilotstudy --suite guacamol --generated-smiles generated.csv
atlas analysis dataset-eda --dataset outputs/data/pilotstudy/ml/training_pass/binding/model/training_design_matrix.csv --label-col spd_binding_label
python -m analysis.cli.audit_graph_model_backends --out outputs/data/pilotstudy/ml/graph_model_backends.csv
python -m analysis.cli.build_mechanism_graph_training_splits --mechanism-pu-table data/pilotstudy/negative_evidence/mechanism_four_state_labels_deep_sources_v7_adrecs_positive_source_balanced.csv --out-dir outputs/data/pilotstudy/ml/graph_splits --panel all --label-col mechanism_ml_label
```

The exploration notebook defaults to `outputs/data/<RUN_ID>/ml/`, defines `REPO_ROOT`, and preserves the old `data/pilotstudy/ml_training_pass_current_defaults/` fallback for archived pilot outputs. It also includes run-scoped manifest and leaderboard helpers: `display_ml_run_summary(RUN_ID, repo_root=REPO_ROOT)` and `display_ml_leaderboard(RUN_ID, repo_root=REPO_ROOT)`.


## Graph/KGE Manifests

Graph/KGE and hypergraph experiments should start from explicit mechanism labels, not unlabeled Atlas rows. The graph split builder writes:

- `graph_nodes.csv`
- `mechanism_hyperedges.csv`
- `kge_reified_triples.csv`
- `negative_edges.csv`
- `negative_sampling_policy.json`
- `split_<mode>/edge_split_manifest.csv`
- `split_<mode>/kge_reified_triples.split.csv`

Production negative edges are only explicit label `0` rows. Unknown rows are not converted into negatives, and corruption-based negative sampling is disabled by default unless a future experiment marks it as benchmark/sensitivity-only and generates corruptions inside each training fold.


## Expanded Tabular Model Families

The shared classifier factory supports these tabular model types in addition to the original logistic/RF/XGBoost/SVM/KNN set:

- `elastic_net_logistic`
- `explainable_boosting_machine` or `ebm`
- `catboost`
- `lightgbm` or `shallow_lightgbm`

Use the same dataset, feature set, split manifest, and leakage flags when comparing model families. The helper installer `bash tools/installers/install_tabular_ml_backends.sh` installs LightGBM, CatBoost, InterpretML core, and the optional dataset EDA stack (`ydata-profiling`, `phik`, `dython`, `networkx`, `pyvis`) while restoring `numpy==1.26.4` for RDKit compatibility.
