from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.papyrus_chembl import build_bioactivity_benchmark
from analysis.io import load_config


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build Papyrus/ChEMBL-style Atlas bioactivity benchmark.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--bioactivity", required=True, type=Path)
    parser.add_argument("--mapping", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--label-prefix", default=None)
    parser.add_argument("--active-threshold-nm", type=float, default=None)
    args = parser.parse_args(argv)
    config = load_config(args.config)
    default_threshold = config.get("bioactivity", {}).get("active_threshold_nM", 10000.0) if isinstance(config.get("bioactivity"), dict) else 10000.0
    build_bioactivity_benchmark(
        args.pair_table,
        args.bioactivity,
        args.mapping,
        args.out,
        active_threshold_nM=float(args.active_threshold_nm or default_threshold),
        label_prefix=args.label_prefix,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
