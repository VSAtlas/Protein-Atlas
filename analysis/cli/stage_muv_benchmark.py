from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.muv import stage_muv_benchmark


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage MoleculeNet MUV as benchmark-only Atlas evidence.")
    parser.add_argument("--muv-csv", required=True, type=Path)
    parser.add_argument("--mapping", type=Path, default=Path("chemdb/data/fda_mapping_from_pdbqt.csv"))
    parser.add_argument("--no-atlas-mapping", action="store_true")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    stage_muv_benchmark(
        args.muv_csv,
        args.out,
        mapping_path=None if args.no_atlas_mapping else args.mapping,
        map_to_atlas=not args.no_atlas_mapping,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
