from __future__ import annotations

import argparse
import os
from pathlib import Path

from analysis.io import load_config
from analysis.ml.train_classifier import train_ml_model


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Train leakage-controlled Atlas ML model.")
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--feature-set", default="full_nonleaky")
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument("--model-n-jobs", type=int, default=0, help="Threads/jobs for model families that support parallelism; 0 keeps model defaults.")
    parser.add_argument("--split", default="drug_holdout")
    parser.add_argument("--split-column", default=None)
    parser.add_argument("--train-values", nargs="*", default=None)
    parser.add_argument("--test-values", nargs="*", default=None)
    parser.add_argument("--temporal-year-col", default=None)
    parser.add_argument("--temporal-cutoff-year", type=int, default=None)
    parser.add_argument("--validation-fold-col", default=None)
    parser.add_argument("--validation-fold-value", default=None)
    parser.add_argument("--validation-fraction", type=float, default=0.15)
    parser.add_argument(
        "--nested-model-selection",
        action="store_true",
        help="Tune supported model hyperparameters on an inner validation partition before scoring the outer holdout.",
    )
    parser.add_argument("--bootstraps", type=int, default=200)
    parser.add_argument(
        "--class-weight",
        choices=["balanced", "none"],
        default="balanced",
        help="Class weighting for supervised classifiers; use none for unweighted logistic regression.",
    )
    parser.add_argument(
        "--pu-mode",
        choices=[
            "standard_binary",
            "positive_unlabeled_weighted",
            "case_control_matched",
            "bagging_pu",
            "stratified_bagging_pu",
            "propensity_weighted_pu",
            "elkan_noto",
            "pulsnar_style",
        ],
        default="standard_binary",
    )
    parser.add_argument("--pu-bags", type=int, default=50)
    parser.add_argument("--pu-unlabeled-ratio", type=float, default=1.0)
    parser.add_argument(
        "--temporal-max-missing-fraction",
        type=float,
        default=0.20,
        help="Block temporal validation if more than this fraction of labeled rows lacks year metadata.",
    )
    parser.add_argument(
        "--exclude-features",
        nargs="*",
        default=None,
        help="Remove listed feature columns from the selected feature set for sensitivity models.",
    )
    parser.add_argument(
        "--allow-label-definition-features",
        action="store_true",
        help=(
            "Allow fields used to define the requested label. Use only for explicitly named "
            "sensitivity models, not primary nonleaky models."
        ),
    )
    parser.add_argument(
        "--allow-partial-rescoring-features",
        action="store_true",
        help=(
            "Allow SCORCH/final rescoring fields. Use only on explicit rescored-subset "
            "datasets where every supervised row has rescoring, so availability is not a feature."
        ),
    )
    parser.add_argument(
        "--group-reweight-cols",
        nargs="*",
        default=None,
        help="Capped inverse-frequency reweighting over source/family strata, e.g. source_family target_family.",
    )
    parser.add_argument(
        "--group-reweight-include-label",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include the supervised label in group-reweight strata so rare positive/negative cells are balanced.",
    )
    parser.add_argument("--group-reweight-max-factor", type=float, default=5.0)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--split-manifest", type=Path, default=None)
    parser.add_argument("--calibration", choices=["none", "raw", "sigmoid", "isotonic"], default="none")
    parser.add_argument("--strict-feature-set", action="store_true", help="Fail if any non-excluded feature in the named feature set is missing from the dataset.")
    parser.add_argument("--no-mlflow", action="store_true")
    args = parser.parse_args(argv)
    if args.no_mlflow:
        os.environ["ATLAS_DISABLE_MLFLOW"] = "1"
    else:
        os.environ.setdefault("ATLAS_ENABLE_MLFLOW", "1")
    config = load_config(args.config)
    seed = int(config.get("project", {}).get("random_seed", 42))
    model_params = None
    if int(args.model_n_jobs) > 0:
        model_params = {"n_jobs": int(args.model_n_jobs), "thread_count": int(args.model_n_jobs)}
    train_ml_model(
        args.dataset,
        args.label,
        args.feature_set,
        args.model,
        args.split,
        args.out_dir,
        seed,
        split_column=args.split_column,
        train_values=args.train_values,
        test_values=args.test_values,
        temporal_year_col=args.temporal_year_col,
        temporal_cutoff_year=args.temporal_cutoff_year,
        n_bootstraps=args.bootstraps,
        pu_mode=args.pu_mode,
        class_weight=None if args.class_weight == "none" else args.class_weight,
        exclude_features=args.exclude_features,
        allow_label_definition_features=args.allow_label_definition_features,
        allow_partial_rescoring_features=args.allow_partial_rescoring_features,
        validation_fold_col=args.validation_fold_col,
        validation_fold_value=args.validation_fold_value,
        validation_fraction=args.validation_fraction,
        nested_model_selection=args.nested_model_selection,
        temporal_max_missing_fraction=args.temporal_max_missing_fraction,
        pu_bags=args.pu_bags,
        pu_unlabeled_ratio=args.pu_unlabeled_ratio,
        claim_mode=args.claim_mode,
        repo_root=args.repo_root,
        split_manifest=args.split_manifest,
        calibration_method=args.calibration,
        strict_feature_set=args.strict_feature_set,
        model_params=model_params,
        group_reweight_cols=args.group_reweight_cols,
        group_reweight_include_label=bool(args.group_reweight_include_label),
        group_reweight_max_factor=float(args.group_reweight_max_factor),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
