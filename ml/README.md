# ML Workflow (CPU-only)

This folder contains a lightweight BigBind v1.5 workflow for:

- downloading and extracting BigBind
- training a logistic-regression classifier (active vs inactive)
- evaluating holdout ranking/classification using:
  - `PR_AUC` from `analysis/dud_eval.py::pr_auc`
  - `EF@{1,2,5,10}%` from `analysis/dud_eval.py::ef_at_fractions`
  - derived `precision@{1,2,5,10}%`

No GPU is required.

Optional dependency for Bayesian search mode:

```bash
pip install optuna
```

Environment note (RDKit/NumPy ABI):

- Keep `numpy` on `1.26.4` in this repo environment.
- Use `rdkit` and avoid installing `rdkit-pypi` into the same env.
- If you see `_ARRAY_API not found`, run:

```bash
python -m pip uninstall -y rdkit-pypi
python -m pip install --force-reinstall rdkit==2025.3.6 numpy==1.26.4 Pillow==11.3.0
```

## 1) Download BigBind v1.5

Run from repo root:

```bash
python ml/download_bigbind.py \
  --url "<BIGBIND_V15_TAR_URL>" \
  --out-dir ml/data/bigbind
```

After extraction, `config.bigbind_dir` should point at the extracted `BigBindV1.5` directory.

## 2) Configure

Default config is `ml/config.json`.

Important fields:

- `bigbind_dir`: path to extracted BigBind directory.
- `train_pdb`: holdout target string.
  - Match order:
    1. case-insensitive match against `ex_rec_pdb`
    2. fallback to case-insensitive match against `pocket` if no `ex_rec_pdb` match exists
- `splits`: which `activities_<split>.csv` files to include.
- `max_rows`: optional dev-only limit.
- `model.class_weight`: defaults to `"balanced"`.
- `atlas_cfg_path`: optional Atlas config override for fpocket path resolution.
- `fpocket_center_columns`: columns used for active-site center lookup.
- `fpocket_variant_column` / `fpocket_ph_column`: optional row columns passed to fpocket loader.
- `features.vina_score`: optional placeholder slot only (no docking/import in this PR).

## 3) Train + Evaluate

```bash
python ml/train_and_eval.py --config ml/config.json
```

Feature ablation experiments:

```bash
python ml/experiment_runner.py --config ml/config.json
```

AutoML search (grid or Optuna):

```bash
python ml/run_automl.py --config ml/config.json --mode grid
python ml/run_automl.py --config ml/config.json --mode optuna --max_trials 50
```

Outputs land in `ml/outputs/<run_id>/`:

- `model.joblib`
- `featurizer_metadata.json`
- `config_snapshot.json`
- `metrics_report.json`
- `metrics_report.csv`
- `holdout_predictions.csv`

AutoML additionally writes:

- `automl_trials.csv`
- `best_trial.json`
- `best_model.joblib`
- `best_holdout_report.json`
- `best_holdout_report.csv`
- `best_holdout_predictions.csv`
- `top_k_holdout_report.csv`

## Feature set

Ligand (RDKit):

- descriptors: MW, LogP, HBD, HBA, TPSA, rotatable bonds, ring count, formal charge
- Bemis-Murcko scaffold (reported in `holdout_predictions.csv`)
- Morgan fingerprint bits (`features.ligand_morgan_fp_bits`, default 2048)

Pocket:

- sourced from fpocket summaries via `druggability_orchestrator.load_fpocket_metrics_for_ml`
- features: `fpocket_druggability`, `fpocket_volume`, `fpocket_openness`,
  `fpocket_polar_fraction`, `fpocket_has_metal`, `fpocket_tier`
- fpocket output root is resolved by `get_fpocket_output_root(atlas_cfg)`
  (default repo path is `protein_automation/fpocket`)

Ensure fpocket has been run for the receptor set you train/evaluate on before
running ML experiments.
