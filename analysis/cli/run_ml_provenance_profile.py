from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.provenance_profile_audit import (
    DEFAULT_EXTENDED_GROUP_COLS,
    DEFAULT_LABEL_COLS,
    run_provenance_profile_audit,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run Atlas ML provenance, shortcut-bias, profiling, association, drift, "
            "validation, and permutation-importance audits."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--labels", nargs="*", default=DEFAULT_LABEL_COLS)
    parser.add_argument("--group-cols", nargs="*", default=DEFAULT_EXTENDED_GROUP_COLS)
    parser.add_argument("--profile-label")
    parser.add_argument("--sample-size", type=int, default=10_000)
    parser.add_argument("--profile-sample-size", type=int, default=5_000)
    parser.add_argument("--max-association-columns", type=int, default=80)
    parser.add_argument("--max-category-levels", type=int, default=200)
    parser.add_argument("--max-mi-features", type=int, default=300)
    parser.add_argument("--permutation-sample-size", type=int, default=5_000)
    parser.add_argument("--permutation-repeats", type=int, default=5)
    parser.add_argument("--permutation-model", default="fast_logistic")
    parser.add_argument(
        "--drift-splits",
        nargs="*",
        default=["target_holdout", "chemical_cluster_holdout", "target_family_holdout"],
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--skip-profile", action="store_true")
    parser.add_argument("--skip-associations", action="store_true")
    parser.add_argument("--skip-permutation-importance", action="store_true")
    parser.add_argument("--skip-drift", action="store_true")
    parser.add_argument("--skip-validation", action="store_true")
    parser.add_argument("--no-augmented-dataset", action="store_true")
    parser.add_argument("--fail-on-missing", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    manifest = run_provenance_profile_audit(
        args.dataset,
        args.out_dir,
        labels=args.labels,
        group_cols=args.group_cols,
        profile_label=args.profile_label,
        run_profile=not args.skip_profile,
        run_associations=not args.skip_associations,
        run_permutation_importance=not args.skip_permutation_importance,
        run_drift=not args.skip_drift,
        run_validation=not args.skip_validation,
        write_augmented_dataset=not args.no_augmented_dataset,
        sample_size=args.sample_size,
        profile_sample_size=args.profile_sample_size,
        max_association_columns=args.max_association_columns,
        max_category_levels=args.max_category_levels,
        max_mi_features=args.max_mi_features,
        permutation_sample_size=args.permutation_sample_size,
        permutation_repeats=args.permutation_repeats,
        permutation_model=args.permutation_model,
        drift_splits=args.drift_splits,
        seed=args.seed,
        test_fraction=args.test_fraction,
        fail_on_missing=args.fail_on_missing,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
