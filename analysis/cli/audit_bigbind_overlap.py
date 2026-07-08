from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.bigbind_overlap_audit import audit_bigbind_overlap


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit Atlas/BANANA model-ready rows against BigBind split membership.")
    parser.add_argument("--atlas-table", required=True, type=Path)
    parser.add_argument("--bigbind-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--fast-raw-ligand",
        action="store_true",
        help="Use cleaned raw SMILES strings for BigBind ligand overlap instead of RDKit canonical/InChIKey generation.",
    )
    args = parser.parse_args(argv)
    summary = audit_bigbind_overlap(
        args.atlas_table,
        args.bigbind_dir,
        args.out_dir,
        canonicalize_bigbind_ligands=not args.fast_raw_ligand,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
