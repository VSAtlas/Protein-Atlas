from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.data_gap_report import DEFAULT_GROUP_COLS, DEFAULT_LABELS, run_data_gap_report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Write claim-oriented Atlas ML data-gap reports showing which "
            "families, sources, scaffolds, and holdout axes block robust claims."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--labels", nargs="*", default=DEFAULT_LABELS)
    parser.add_argument("--group-cols", nargs="*", default=DEFAULT_GROUP_COLS)
    parser.add_argument("--min-group-labelable", type=int, default=20)
    parser.add_argument("--min-group-positives", type=int, default=10)
    parser.add_argument("--min-group-negatives", type=int, default=10)
    parser.add_argument("--min-holdout-test-positives", type=int, default=10)
    parser.add_argument("--min-holdout-test-negatives", type=int, default=10)
    parser.add_argument("--min-holdout-train-positives", type=int, default=20)
    parser.add_argument("--min-holdout-train-negatives", type=int, default=20)
    parser.add_argument("--min-calibration-positives", type=int, default=20)
    parser.add_argument("--min-calibration-negatives", type=int, default=20)
    parser.add_argument("--severe-positive-rate", type=float, default=0.85)
    parser.add_argument("--sparse-positive-rate", type=float, default=0.05)
    parser.add_argument("--max-unknown-fraction", type=float, default=0.90)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = run_data_gap_report(
        args.dataset,
        args.out_dir,
        labels=args.labels,
        group_cols=args.group_cols,
        min_group_labelable=args.min_group_labelable,
        min_group_positives=args.min_group_positives,
        min_group_negatives=args.min_group_negatives,
        min_holdout_test_positives=args.min_holdout_test_positives,
        min_holdout_test_negatives=args.min_holdout_test_negatives,
        min_holdout_train_positives=args.min_holdout_train_positives,
        min_holdout_train_negatives=args.min_holdout_train_negatives,
        min_calibration_positives=args.min_calibration_positives,
        min_calibration_negatives=args.min_calibration_negatives,
        severe_positive_rate=args.severe_positive_rate,
        sparse_positive_rate=args.sparse_positive_rate,
        max_unknown_fraction=args.max_unknown_fraction,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
