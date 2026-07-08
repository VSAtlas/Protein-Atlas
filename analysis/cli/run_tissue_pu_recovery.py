from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.tissue_pu_recovery import run_tissue_pu_recovery_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run tissue/site positive-unlabeled recovery. This ranks held-out "
            "tissue positives against held-out unknown/background rows and reports "
            "top-K recovery instead of treating unknowns as confirmed negatives."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", default="tissue_site_label")
    parser.add_argument("--feature-set", default="spd_tissue_site_nonleaky")
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument("--split", default="target_holdout")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--pu-mode",
        choices=["bagging_pu", "stratified_bagging_pu"],
        default="stratified_bagging_pu",
    )
    parser.add_argument("--pu-bags", type=int, default=50)
    parser.add_argument("--pu-unlabeled-ratio", type=float, default=1.0)
    parser.add_argument("--top-k", nargs="+", type=int, default=[10, 20, 50])
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--max-eval-background", type=int, default=10000)
    parser.add_argument("--class-weight", choices=["balanced", "none"], default="balanced")
    parser.add_argument(
        "--allow-tissue-expression-features",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Allow expression/site features for an externally sourced tissue label. "
            "Disable this when auditing rule-derived tissue labels."
        ),
    )
    parser.add_argument("--strict-feature-set", action="store_true")
    parser.add_argument("--case-study-rows", type=int, default=25)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args(argv)

    result = run_tissue_pu_recovery_report(
        args.dataset,
        label_col=args.label,
        feature_set=args.feature_set,
        model_type=args.model,
        split_mode=args.split,
        out_dir=args.out_dir,
        seed=args.seed,
        pu_mode=args.pu_mode,
        pu_bags=args.pu_bags,
        pu_unlabeled_ratio=args.pu_unlabeled_ratio,
        top_k=args.top_k,
        test_fraction=args.test_fraction,
        max_eval_background=args.max_eval_background,
        class_weight=None if args.class_weight == "none" else args.class_weight,
        allow_tissue_expression_features=bool(args.allow_tissue_expression_features),
        strict_feature_set=bool(args.strict_feature_set),
        case_study_rows=int(args.case_study_rows),
        repo_root=args.repo_root,
    )
    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
