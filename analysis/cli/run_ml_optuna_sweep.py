from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

from analysis.io import load_config
from analysis.ml.optuna_sweep import OPTUNA_METRICS, SUPPORTED_OPTUNA_PU_MODES, run_optuna_sweep
from analysis.ml.train_classifier import train_ml_model


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run an Optuna-backed Atlas ML hyperparameter sweep, then train the selected model on the outer split."
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--feature-set", default="full_nonleaky")
    parser.add_argument("--model", default="logistic_regression", choices=["logistic_regression", "elastic_net_logistic", "random_forest", "xgboost", "gradient_boosted_trees", "lightgbm", "shallow_lightgbm", "catboost", "explainable_boosting_machine", "ebm"])
    parser.add_argument("--split", default="drug_holdout")
    parser.add_argument("--split-column", default=None)
    parser.add_argument("--train-values", nargs="*", default=None)
    parser.add_argument("--test-values", nargs="*", default=None)
    parser.add_argument("--temporal-year-col", default=None)
    parser.add_argument("--temporal-cutoff-year", type=int, default=None)
    parser.add_argument("--temporal-max-missing-fraction", type=float, default=0.20)
    parser.add_argument("--validation-fold-col", default=None)
    parser.add_argument("--validation-fold-value", default=None)
    parser.add_argument("--validation-fraction", type=float, default=0.15)
    parser.add_argument("--trials", type=int, default=25)
    parser.add_argument("--timeout", type=int, default=None)
    parser.add_argument("--metric", choices=OPTUNA_METRICS, default="AUPRC")
    parser.add_argument("--storage", default=None, help="Optional Optuna storage URL, for example sqlite:///path/to/study.db")
    parser.add_argument("--pruner", choices=["median", "successive_halving", "hyperband", "none"], default="median")
    parser.add_argument("--split-manifest", type=Path, default=None)
    parser.add_argument("--no-mlflow", action="store_true")
    parser.add_argument("--study-name", default=None)
    parser.add_argument("--load-if-exists", action="store_true")
    parser.add_argument("--bootstraps", type=int, default=200)
    parser.add_argument("--class-weight", choices=["balanced", "none"], default="balanced")
    parser.add_argument("--pu-mode", choices=sorted(SUPPORTED_OPTUNA_PU_MODES), default="standard_binary")
    parser.add_argument("--exclude-features", nargs="*", default=None)
    parser.add_argument(
        "--allow-label-definition-features",
        action="store_true",
        help="Allow label-definition fields only for explicit leakage sensitivity sweeps.",
    )
    parser.add_argument(
        "--allow-partial-rescoring-features",
        action="store_true",
        help="Allow SCORCH/final rescoring fields only on fully rescored subset datasets.",
    )
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--skip-final-train", action="store_true", help="Run only the Optuna study and write hpo/optuna artifacts.")
    args = parser.parse_args(argv)

    if args.no_mlflow:
        os.environ["ATLAS_DISABLE_MLFLOW"] = "1"
    else:
        os.environ.setdefault("ATLAS_ENABLE_MLFLOW", "1")
    config = load_config(args.config)
    seed = int(args.seed if args.seed is not None else config.get("project", {}).get("random_seed", 42))
    class_weight = None if args.class_weight == "none" else args.class_weight
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    sweep = run_optuna_sweep(
        dataset_path=args.dataset,
        label_col=args.label,
        feature_set=args.feature_set,
        model_type=args.model,
        split_mode=args.split,
        out_dir=out,
        seed=seed,
        split_column=args.split_column,
        train_values=args.train_values,
        test_values=args.test_values,
        temporal_year_col=args.temporal_year_col,
        temporal_cutoff_year=args.temporal_cutoff_year,
        temporal_max_missing_fraction=args.temporal_max_missing_fraction,
        validation_fold_col=args.validation_fold_col,
        validation_fold_value=args.validation_fold_value,
        validation_fraction=args.validation_fraction,
        n_trials=args.trials,
        timeout=args.timeout,
        metric=args.metric,
        class_weight=class_weight,
        pu_mode=args.pu_mode,
        exclude_features=args.exclude_features,
        split_manifest=args.split_manifest,
        allow_label_definition_features=args.allow_label_definition_features,
        allow_partial_rescoring_features=args.allow_partial_rescoring_features,
        study_name=args.study_name,
        storage=args.storage,
        load_if_exists=args.load_if_exists,
        pruner=args.pruner,
    )

    result: dict[str, Any] = {
        "status": "sweep_complete",
        "out_dir": str(out),
        "hpo_manifest": sweep.get("manifest_path"),
        "hpo_trials": sweep.get("trials_path"),
        "best_params": sweep.get("best_params", {}),
    }
    if not args.skip_final_train:
        model_dir = out / "model"
        hpo_metadata = {
            "backend": "optuna",
            "study_name": sweep.get("study_name"),
            "selection_metric": sweep.get("selection_metric"),
            "n_trials": sweep.get("n_trials"),
            "trials_path": sweep.get("trials_path"),
            "manifest_path": sweep.get("manifest_path"),
            "best_trial_number": sweep.get("best_trial_number"),
            "best_value": sweep.get("best_value"),
            "best_params": sweep.get("best_params", {}),
            "policy": "Optuna selected hyperparameters on fixed validation data; train_ml_model owns final outer-holdout scoring.",
        }
        trained = train_ml_model(
            args.dataset,
            args.label,
            args.feature_set,
            args.model,
            args.split,
            model_dir,
            seed,
            split_column=args.split_column,
            train_values=args.train_values,
            test_values=args.test_values,
            temporal_year_col=args.temporal_year_col,
            temporal_cutoff_year=args.temporal_cutoff_year,
            n_bootstraps=args.bootstraps,
            pu_mode=args.pu_mode,
            class_weight=class_weight,
            exclude_features=args.exclude_features,
            allow_label_definition_features=args.allow_label_definition_features,
            allow_partial_rescoring_features=args.allow_partial_rescoring_features,
            split_manifest=args.split_manifest,
            validation_fold_col=args.validation_fold_col,
            validation_fold_value=args.validation_fold_value,
            validation_fraction=args.validation_fraction,
            nested_model_selection=False,
            temporal_max_missing_fraction=args.temporal_max_missing_fraction,
            claim_mode=args.claim_mode,
            repo_root=args.repo_root,
            model_params=sweep.get("best_params", {}),
            hpo_metadata=hpo_metadata,
        )
        result.update(
            {
                "status": "trained",
                "model_dir": str(model_dir),
                "metrics": trained.get("metrics", {}),
                "n_train": trained.get("n_train"),
                "n_test": trained.get("n_test"),
                "claim_readiness": trained.get("claim_readiness", {}),
            }
        )

    summary_path = out / "ml_optuna_sweep_manifest.json"
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True, default=str), encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
