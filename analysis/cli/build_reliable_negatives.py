from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.reliable_negatives import build_reliable_negative_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Mine fold-local reliable negatives using external similarity/descriptors only."
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--selection-features", nargs="+", required=True)
    parser.add_argument("--split", default="drug_holdout")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--n-per-positive", type=float, default=1.0)
    parser.add_argument("--max-negatives", type=int, default=None)
    args = parser.parse_args(argv)
    build_reliable_negative_table(
        args.dataset,
        args.out,
        label_col=args.label,
        selection_features=args.selection_features,
        split_mode=args.split,
        seed=args.seed,
        test_fraction=args.test_fraction,
        n_per_positive=args.n_per_positive,
        max_negatives=args.max_negatives,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
