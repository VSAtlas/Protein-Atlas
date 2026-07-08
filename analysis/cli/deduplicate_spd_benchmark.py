from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.spd import deduplicate_spd_benchmark


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Deduplicate an SPD benchmark to one row per drug-target pair."
    )
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--ml-ready-out", type=Path, default=None)
    args = parser.parse_args(argv)
    deduplicate_spd_benchmark(
        args.benchmark,
        args.out,
        ml_ready_out_path=args.ml_ready_out,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
