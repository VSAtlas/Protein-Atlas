from __future__ import annotations

import argparse

from analysis.ml.grouped_cv_ablation_suite import run_grouped_cv_ablation_suite


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run SPD/addon-aware grouped OOF ablations, rank-based failure analysis, "
            "score provenance checks, and the explicit exposure-margin model."
        )
    )
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--run-dir", default=None)
    parser.add_argument(
        "--feature-table",
        action="append",
        default=[],
        help="Authoritative score/metadata table; repeat for multiple tables.",
    )
    parser.add_argument("--no-feature-refresh", action="store_true")
    parser.add_argument("--models", nargs="+", default=["logistic_regression", "lightgbm"])
    parser.add_argument(
        "--group-cols",
        nargs="+",
        default=["drug_id", "chemical_cluster", "target_id", "target_family"],
    )
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--min-test-positives", type=int, default=10)
    parser.add_argument("--min-test-negatives", type=int, default=10)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument(
        "--error-top-k",
        type=int,
        default=None,
        help="FP/FN operating K; default matches K to the number of positives.",
    )
    parser.add_argument("--bootstraps", type=int, default=200)
    parser.add_argument("--permutations", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--model-n-jobs", type=int, default=4)
    parser.add_argument("--no-applicability-domain", action="store_true")
    parser.add_argument("--skip-exposure", action="store_true")
    parser.add_argument("--exposure-models", nargs="+", default=["ridge", "lightgbm"])
    args = parser.parse_args(argv)
    run_grouped_cv_ablation_suite(
        args.dataset,
        args.out_dir,
        model_types=args.models,
        group_cols=args.group_cols,
        n_splits=args.n_splits,
        repeats=args.repeats,
        min_test_positives=args.min_test_positives,
        min_test_negatives=args.min_test_negatives,
        top_k=args.top_k,
        error_top_k=args.error_top_k,
        n_bootstraps=args.bootstraps,
        n_permutations=args.permutations,
        seed=args.seed,
        model_n_jobs=args.model_n_jobs,
        compute_applicability_domain=not args.no_applicability_domain,
        run_exposure=not args.skip_exposure,
        exposure_models=args.exposure_models,
        refresh_feature_metadata=not args.no_feature_refresh,
        run_dir=args.run_dir,
        extra_feature_tables=args.feature_table,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
