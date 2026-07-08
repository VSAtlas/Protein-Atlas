from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.banana import build_banana_input_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build BANANA/BigBind inference inputs from Atlas pair rows.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--fda-mapping", type=Path, default=Path("chemdb/data/fda_mapping_from_pdbqt.csv"))
    parser.add_argument(
        "--pocket-map",
        type=Path,
        default=None,
        help="Optional CSV with pdb_id,pocket_pdb columns. Missing pockets remain explicit.",
    )
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    build_banana_input_table(
        args.pair_table,
        args.out,
        fda_mapping_path=args.fda_mapping,
        pocket_map_path=args.pocket_map,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
