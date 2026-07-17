from __future__ import annotations

import argparse

from analysis.ml.grouped_cv_stability import run_grouped_cv_stability


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run repeated grouped-CV stability diagnostics for ML tables.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--feature-set", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--models",
        nargs="+",
        default=["logistic_regression", "random_forest", "lightgbm"],
        help="Model families to evaluate.",
    )
    parser.add_argument(
        "--group-cols",
        nargs="+",
        default=["target_id", "target_family"],
        help="Columns used for repeated grouped CV.",
    )
    parser.add_argument(
        "--leave-one-group-cols",
        nargs="*",
        default=["target_family"],
        help="Columns used for leave-one-group-out diagnostics.",
    )
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--min-test-positives", type=int, default=10)
    parser.add_argument("--min-test-negatives", type=int, default=10)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--bootstraps", type=int, default=200)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--class-weight", default="balanced")
    parser.add_argument("--pu-mode", default="standard_binary")
    parser.add_argument("--exclude-features", nargs="*", default=None)
    parser.add_argument("--strict-feature-set", action="store_true")
    parser.add_argument(
        "--applicability-domain",
        action="store_true",
        help="Compute fold-specific nearest-neighbor applicability distances for held-out rows.",
    )
    parser.add_argument("--stratified-target-holdout-by-family", action="store_true")
    parser.add_argument("--stratified-target-holdout-repeats", type=int, default=None)
    parser.add_argument("--stratified-target-holdout-fraction", type=float, default=0.2)
    args = parser.parse_args(argv)
    run_grouped_cv_stability(
        args.dataset,
        label_col=args.label,
        feature_set=args.feature_set,
        out_dir=args.out_dir,
        model_types=args.models,
        group_cols=args.group_cols,
        leave_one_group_cols=args.leave_one_group_cols,
        n_splits=args.n_splits,
        repeats=args.repeats,
        min_test_positives=args.min_test_positives,
        min_test_negatives=args.min_test_negatives,
        top_k=args.top_k,
        n_bootstraps=args.bootstraps,
        n_permutations=args.permutations,
        seed=args.seed,
        class_weight=None if str(args.class_weight).lower() in {"none", "false", "0"} else args.class_weight,
        pu_mode=args.pu_mode,
        exclude_features=args.exclude_features,
        strict_feature_set=args.strict_feature_set,
        compute_applicability_domain=args.applicability_domain,
        run_stratified_target_holdout_by_family=args.stratified_target_holdout_by_family,
        stratified_target_holdout_repeats=args.stratified_target_holdout_repeats,
        stratified_target_holdout_fraction=args.stratified_target_holdout_fraction,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
