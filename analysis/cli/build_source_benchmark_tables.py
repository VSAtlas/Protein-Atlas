from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.source_benchmark_tables import build_source_benchmark_tables


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build separate curated, ToxCast, and SPD Atlas ML benchmark tables.")
    parser.add_argument("--bioactivity-source", required=True, type=Path)
    parser.add_argument("--spd-table", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    build_source_benchmark_tables(
        args.bioactivity_source,
        args.out_dir,
        spd_table_path=args.spd_table,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
