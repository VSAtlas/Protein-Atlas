from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.toxcast import build_toxcast_benchmark
from analysis.io import load_config


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build ToxCast activity benchmark.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--toxcast", required=True, type=Path)
    parser.add_argument("--mapping", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--potent-threshold-nm", type=float, default=10000.0)
    args = parser.parse_args(argv)
    _config = load_config(args.config)
    build_toxcast_benchmark(args.pair_table, args.toxcast, args.mapping, args.out, args.potent_threshold_nm)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

