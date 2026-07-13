"""CLI for identity-safe BindingDB FDA inactive staging."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.bindingdb_fda_inactive_stage import (
    build_target_specs,
    stage_bindingdb_fda_inactives,
)
from analysis.ml.spd_identity_reconciliation import default_fda_mapping_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Stage measured human BindingDB Ki/Kd/IC50 FDA inactives using only "
            "exact full-InChIKey identity matches to checksum-validated canonical "
            "PDBQTs. Untested pairs remain unknown; parent/name-only matching is forbidden."
        )
    )
    parser.add_argument(
        "--bindingdb",
        required=True,
        type=Path,
        help="Normalized BindingDB TSV.",
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=default_fda_mapping_path(),
        help=(
            "FDA identity mapping (default: the active canonical "
            "chemdb/data/fda_mapping_from_pdbqt.csv)."
        ),
    )
    parser.add_argument(
        "--target-uniprot",
        required=True,
        action="append",
        help="Human target UniProt accession; repeat with aligned gene/PDB options.",
    )
    parser.add_argument(
        "--target-gene",
        required=True,
        action="append",
        help="Target gene symbol; repeat in the same order as --target-uniprot.",
    )
    parser.add_argument(
        "--pdb-id",
        required=True,
        action="append",
        help="Atlas PDB ID; repeat in the same order as --target-uniprot.",
    )
    parser.add_argument(
        "--inactive-nm",
        "--inactive-threshold-nm",
        dest="inactive_nm",
        type=float,
        default=10_000.0,
        help="Safe measured-inactive threshold in nM (default: 10000).",
    )
    parser.add_argument(
        "--active-nm",
        "--active-threshold-nm",
        dest="active_nm",
        type=float,
        default=1_000.0,
        help="Same-target active-conflict threshold in nM (default: 1000).",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--library-name", required=True)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.txt"),
        help="Atlas runtime config used to resolve prepped_ligands (default: config.txt).",
    )
    parser.add_argument(
        "--prepped-root",
        type=Path,
        default=None,
        help="Optional explicit prepped_ligands root; otherwise use the path router.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=100_000,
        help="BindingDB TSV streaming chunk size (default: 100000).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    targets = build_target_specs(
        args.target_uniprot, args.target_gene, args.pdb_id
    )
    manifest = stage_bindingdb_fda_inactives(
        bindingdb_path=args.bindingdb,
        mapping_path=args.mapping,
        targets=targets,
        out_dir=args.out_dir,
        library_name=args.library_name,
        inactive_nm=args.inactive_nm,
        active_nm=args.active_nm,
        config_path=args.config,
        prepped_root=args.prepped_root,
        chunksize=args.chunksize,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
