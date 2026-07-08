from __future__ import annotations

import os
import shutil
import sys
from typing import Any

from config.runtime_config import load_config, validate_config
from protein_prep.phenix_tools import run_phenix_pdbtools


def _load_cfg() -> dict[str, Any]:
    cfg = load_config()
    validate_config(cfg)
    return cfg


def clean_pdb(pdb_file: str, output_dir: str | None) -> str:
    cfg = _load_cfg()
    pdb_id = os.path.basename(pdb_file).split(".")[0]

    if not output_dir:
        output_dir = str(
            cfg.get("PHENIX_CLEAN_OUTDIR", os.path.join(str(cfg["OUTPUT_DIR"]), "phenix_clean"))
        )

    out_path = os.path.join(output_dir, pdb_id)
    os.makedirs(out_path, exist_ok=True)

    suffix = str(cfg.get("PHENIX_CLEAN_SUFFIX", "_cleaned.pdb"))
    cleaned_pdb_path = os.path.join(out_path, f"{pdb_id}{suffix}")

    if not run_phenix_pdbtools(pdb_file, cleaned_pdb_path, remove_waters=True):
        shutil.copyfile(pdb_file, cleaned_pdb_path)
        print(f"Phenix unavailable; copied input PDB: {cleaned_pdb_path}")

    cleaned_pdb_path = os.path.abspath(cleaned_pdb_path)
    print(f"Phenix cleaned: {cleaned_pdb_path}")
    return cleaned_pdb_path


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        raise SystemExit("usage: atlas-phenix-clean <pdb_file> [output_dir]")

    pdb_file = args[0]
    output_dir = args[1] if len(args) > 1 else None
    clean_pdb(pdb_file, output_dir)
    return 0
