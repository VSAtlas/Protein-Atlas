from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.leakage_overlap_audit import DEFAULT_SPLITS, audit_ml_leakage_overlap


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audit ML train/test split overlap, duplicate labels, and forbidden predictive features."
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--feature-set")
    parser.add_argument("--exclude-features", nargs="*", default=[])
    parser.add_argument("--splits", nargs="*", default=DEFAULT_SPLITS)
    parser.add_argument("--source-col", default="label_source")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument(
        "--allow-partial-rescoring-features",
        action="store_true",
        help="Allow SCORCH/final fields only for explicit rescored-subset audits.",
    )
    args = parser.parse_args(argv)
    manifest = audit_ml_leakage_overlap(
        args.dataset,
        args.label,
        args.out_dir,
        feature_set=args.feature_set,
        exclude_features=args.exclude_features,
        split_modes=args.splits,
        source_col=args.source_col,
        seed=args.seed,
        test_fraction=args.test_fraction,
        allow_partial_rescoring_features=args.allow_partial_rescoring_features,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
