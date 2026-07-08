from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.labels import build_external_label_tables
from analysis.io import load_config


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build joined local-first Atlas external label tables.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args(argv)
    build_external_label_tables(args.pair_table, load_config(args.config), args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
