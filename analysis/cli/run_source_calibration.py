from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.source_calibration import run_source_calibration


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run source-aware calibration for cross-source ML transfer.")
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--feature-set", default="pilot_nonleaky")
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument("--train-sources", nargs="+", required=True)
    parser.add_argument("--calibration-source", required=True)
    parser.add_argument("--source-col", default="label_source", help="Column used to identify train/calibration source groups.")
    parser.add_argument("--calibration-method", choices=["logistic", "isotonic"], default="logistic")
    parser.add_argument("--class-weight", choices=["balanced", "none"], default="balanced")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--group-reweight-cols",
        nargs="*",
        default=None,
        help="Optional capped inverse-frequency reweighting over source/family strata.",
    )
    parser.add_argument(
        "--group-reweight-include-label",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include the supervised label in group-reweight strata.",
    )
    parser.add_argument("--group-reweight-max-factor", type=float, default=5.0)
    args = parser.parse_args(argv)
    run_source_calibration(
        args.dataset,
        args.label,
        args.feature_set,
        args.model,
        args.train_sources,
        args.calibration_source,
        args.out_dir,
        calibration_method=args.calibration_method,
        seed=args.seed,
        class_weight=None if args.class_weight == "none" else args.class_weight,
        group_reweight_cols=args.group_reweight_cols,
        group_reweight_include_label=bool(args.group_reweight_include_label),
        group_reweight_max_factor=float(args.group_reweight_max_factor),
        source_col=args.source_col,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
