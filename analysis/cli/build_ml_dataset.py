from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.build_ml_dataset import build_ml_dataset


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build Atlas ML training dataset.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--feature-set", default="full_nonleaky")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--deduplicate-drug-target",
        action="store_true",
        help="Keep one labeled row per drug-target pair before ML training.",
    )
    parser.add_argument(
        "--group-cols",
        nargs="*",
        default=None,
        help="Optional grouping columns for deduplication; defaults to drug_id target_id.",
    )
    parser.add_argument(
        "--eligible-col",
        default=None,
        help="Optional column that must match --eligible-value before ML rows are kept.",
    )
    parser.add_argument(
        "--eligible-value",
        default="1",
        help="Required value for --eligible-col; default: 1.",
    )
    parser.add_argument(
        "--include-unlabeled-as-background",
        action="store_true",
        help="Keep unlabeled rows as background negatives with sample weights for PU-style training.",
    )
    parser.add_argument(
        "--unlabeled-weight",
        type=float,
        default=0.25,
        help="Sample weight for unlabeled/background rows when PU mode is enabled.",
    )
    parser.add_argument(
        "--label-source-include",
        nargs="*",
        default=None,
        help="Keep rows whose semicolon-separated label_source contains any listed source.",
    )
    parser.add_argument(
        "--require-year-col",
        default=None,
        help="Keep only rows with a nonmissing numeric value in this year column.",
    )
    parser.add_argument(
        "--require-nonmissing-cols",
        nargs="*",
        default=None,
        help="Keep only rows with nonmissing numeric values in all listed columns.",
    )
    parser.add_argument(
        "--exclude-features",
        nargs="*",
        default=None,
        help="Remove listed feature columns from the selected feature set before writing the ML table.",
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
        "--keep-all-columns",
        action="store_true",
        help=(
            "Write all source columns to the model-ready table while still excluding non-selected "
            "or explicitly excluded columns from the training feature set."
        ),
    )
    parser.add_argument(
        "--exclude-flag-cols",
        nargs="*",
        default=None,
        help="Drop rows where any listed boolean flag column is true.",
    )
    args = parser.parse_args(argv)
    build_ml_dataset(
        args.pair_table,
        args.label,
        args.feature_set,
        args.out,
        deduplicate_drug_target=args.deduplicate_drug_target,
        group_cols=args.group_cols,
        eligible_col=args.eligible_col,
        eligible_value=args.eligible_value,
        include_unlabeled_as_background=args.include_unlabeled_as_background,
        unlabeled_weight=args.unlabeled_weight,
        label_source_include=args.label_source_include,
        require_year_col=args.require_year_col,
        require_nonmissing_cols=args.require_nonmissing_cols,
        exclude_flag_cols=args.exclude_flag_cols,
        exclude_features=args.exclude_features,
        allow_label_definition_features=args.allow_label_definition_features,
        keep_all_columns=args.keep_all_columns,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
