from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.target_chemotype_bias import DEFAULT_GROUP_COLS, audit_target_chemotype_bias


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Audit target/chemotype shortcut bias in an Atlas ML table by combining "
            "group prevalence diagnostics with metadata-only shortcut models."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--group-cols", nargs="*", default=DEFAULT_GROUP_COLS)
    parser.add_argument(
        "--splits",
        nargs="*",
        default=["target_holdout", "target_family_holdout", "chemical_cluster_holdout"],
    )
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--min-group-n", type=int, default=10)
    parser.add_argument("--min-positive", type=int, default=5)
    parser.add_argument("--min-negative", type=int, default=5)
    parser.add_argument("--dominance-threshold", type=float, default=0.25)
    parser.add_argument("--model-n-jobs", type=int, default=4)
    args = parser.parse_args(argv)

    manifest = audit_target_chemotype_bias(
        args.dataset,
        label_col=args.label,
        out_dir=args.out_dir,
        group_cols=args.group_cols,
        splits=args.splits,
        model_type=args.model,
        seed=args.seed,
        top_k=args.top_k,
        min_group_n=args.min_group_n,
        min_positive=args.min_positive,
        min_negative=args.min_negative,
        dominance_threshold=args.dominance_threshold,
        model_n_jobs=args.model_n_jobs,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
