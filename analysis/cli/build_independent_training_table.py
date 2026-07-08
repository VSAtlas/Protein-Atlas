from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.dataset_independence_audit import build_independent_training_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a training table with held-out benchmark/calibration pairs removed."
    )
    parser.add_argument("--train", required=True, type=Path)
    parser.add_argument("--heldout", nargs="+", required=True, help="Held-out tables as NAME=PATH or PATH.")
    parser.add_argument("--exclude-source-token", nargs="*", default=[])
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    build_independent_training_table(
        args.train,
        args.heldout,
        args.out,
        exclude_source_tokens=args.exclude_source_token,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
