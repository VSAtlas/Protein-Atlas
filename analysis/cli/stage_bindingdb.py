from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.bindingdb import stage_bindingdb_for_atlas


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stream-filter BindingDB to Atlas drugs and target UniProt IDs.")
    parser.add_argument("--bindingdb", required=True, type=Path)
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--chunksize", type=int, default=200_000)
    args = parser.parse_args(argv)
    stage_bindingdb_for_atlas(
        args.bindingdb,
        args.pair_table,
        args.mapping,
        args.out,
        chunksize=args.chunksize,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
