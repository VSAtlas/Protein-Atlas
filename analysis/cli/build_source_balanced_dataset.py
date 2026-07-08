from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.source_balance import build_source_balanced_dataset


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build a source-balanced ML sensitivity dataset.")
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--source-col", default="negative_source")
    parser.add_argument("--negative-label", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-per-source", type=int, default=None)
    args = parser.parse_args(argv)
    build_source_balanced_dataset(
        args.dataset,
        args.out,
        label_col=args.label,
        source_col=args.source_col,
        negative_label=args.negative_label,
        seed=args.seed,
        max_per_source=args.max_per_source,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
