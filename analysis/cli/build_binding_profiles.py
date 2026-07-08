from __future__ import annotations

import argparse
from pathlib import Path

from analysis.profiles.binding_profiles import build_binding_profile


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build drug x target Atlas binding profiles.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--value", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--target-space", type=Path, default=None)
    parser.add_argument("--dense", action="store_true")
    parser.add_argument("--protein-class", default=None)
    parser.add_argument("--min-structure-quality", type=float, default=None)
    parser.add_argument("--pdb-panel", nargs="*", default=None)
    args = parser.parse_args(argv)
    build_binding_profile(
        args.pair_table,
        args.out,
        value_col=args.value,
        target_space_path=args.target_space,
        dense=args.dense,
        protein_class=args.protein_class,
        min_structure_quality=args.min_structure_quality,
        pdb_panel=args.pdb_panel,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
