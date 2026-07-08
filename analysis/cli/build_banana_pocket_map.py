from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.banana_pockets import build_banana_pocket_map


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build BANANA pocket-residue PDB files and a pdb_id,pocket_pdb map."
    )
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--raw-pdb-dir",
        type=Path,
        default=None,
        help="Raw PDB cache directory. Defaults to <out parent>/raw_pdbs.",
    )
    parser.add_argument(
        "--pocket-dir",
        type=Path,
        default=None,
        help="Pocket PDB output directory. Defaults to <out parent>/pockets.",
    )
    parser.add_argument(
        "--receptor-search-dir",
        action="append",
        type=Path,
        default=[],
        help="Directory to search for existing receptor PDBs; can be repeated.",
    )
    parser.add_argument("--radius", type=float, default=8.0)
    parser.add_argument("--min-ligand-heavy-atoms", type=int, default=6)
    parser.add_argument("--no-download", action="store_true")
    parser.add_argument("--allow-glycan-centers", action="store_true")
    parser.add_argument("--no-geometric-fallback", action="store_true")
    parser.add_argument("--include-hetatm", action="store_true")
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args(argv)
    raw_pdb_dir = args.raw_pdb_dir or (args.out.parent / "raw_pdbs")
    pocket_dir = args.pocket_dir or (args.out.parent / "pockets")
    build_banana_pocket_map(
        args.pair_table,
        args.out,
        raw_pdb_dir=raw_pdb_dir,
        pocket_dir=pocket_dir,
        receptor_search_dirs=args.receptor_search_dir or None,
        radius_a=args.radius,
        min_ligand_heavy_atoms=args.min_ligand_heavy_atoms,
        allow_download=not args.no_download,
        allow_glycan_centers=args.allow_glycan_centers,
        allow_geometric_fallback=not args.no_geometric_fallback,
        include_hetatm=args.include_hetatm,
        require_complete=not args.allow_incomplete,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
