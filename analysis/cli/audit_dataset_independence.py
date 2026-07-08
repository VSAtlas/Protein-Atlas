from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.dataset_independence_audit import audit_dataset_independence


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audit train-vs-calibration/benchmark canonical pair and source overlap."
    )
    parser.add_argument("--train", required=True, type=Path)
    parser.add_argument(
        "--candidate",
        nargs="+",
        required=True,
        help="Candidate benchmark/calibration tables as NAME=PATH or PATH.",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--allow-pair-overlap", action="store_true")
    parser.add_argument("--strict-source-overlap", action="store_true")
    args = parser.parse_args(argv)
    audit_dataset_independence(
        args.train,
        args.candidate,
        args.out_dir,
        strict_pair_overlap=not args.allow_pair_overlap,
        strict_source_overlap=args.strict_source_overlap,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
