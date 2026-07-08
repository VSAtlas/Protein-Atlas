from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.excape import stage_excape_for_atlas


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stream-filter ExCAPE/PubChem+ChEMBL activity rows to Atlas pairs.")
    parser.add_argument("--excape", required=True, type=Path)
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--chunksize", type=int, default=500_000)
    args = parser.parse_args(argv)
    stage_excape_for_atlas(
        args.excape,
        args.pair_table,
        args.mapping,
        args.out,
        chunksize=args.chunksize,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
