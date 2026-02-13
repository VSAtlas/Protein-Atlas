"""Canonical output layout and legacy tree migration helpers."""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path
from typing import Dict, Optional, Union

import protein_prep.prep_utils as prep_utils
from protein_prep.os_utils import _as_path
from protein_prep.prep_utils import _resolve_variant_token


def _canon_base(output_root: Path, pdb_id: str) -> Path:
    """Canonical per-protein base: processed_pdbs/<PDB>"""
    output_root = _as_path(output_root)
    return (output_root / pdb_id.upper()).resolve()


def _merge_dir(src: Path, dst: Path) -> None:
    """Merge src directory into dst; remove src after moving."""
    src, dst = src.resolve(), dst.resolve()
    if not src.exists():
        return
    dst.mkdir(parents=True, exist_ok=True)
    for root, dirs, files in os.walk(src):
        r = Path(root)
        rel = r.relative_to(src)
        (dst / rel).mkdir(parents=True, exist_ok=True)
        for d in dirs:
            (dst / rel / d).mkdir(parents=True, exist_ok=True)
        for f in files:
            s = r / f
            t = dst / rel / f
            if t.exists():
                try:
                    if s.stat().st_size == t.stat().st_size:
                        continue
                except Exception:
                    pass
            shutil.move(str(s), str(t))
    try:
        shutil.rmtree(src)
    except Exception:
        pass


def fold_legacy_layout(pdb_id: str, output_root) -> None:
    """
    Extended: migrate uppercase legacy dirs into canonical tree:
      <PDB>_NOLIG            -> processed_pdbs/<PDB>/nolig
      <PDB>_CLEANED_LIGANDS  -> processed_pdbs/<PDB>/ligands_raw
      <PDB>_nolig(.pdb)      -> processed_pdbs/<PDB>/nolig/<PDB>_nolig_phenix_clean.pdb
    """
    try:
        root = _as_path(output_root).resolve()
        pdb_idU = pdb_id.upper()
        base = _canon_base(root, pdb_idU)
        (base / "nolig").mkdir(parents=True, exist_ok=True)
        (base / "ligands_raw").mkdir(parents=True, exist_ok=True)

        legacy_dirs = [
            (root / f"{pdb_id}_nolig", base / "nolig"),
            (root / f"{pdb_id.lower()}_nolig", base / "nolig"),
            (root / f"{pdb_idU}_NOLIG", base / "nolig"),
            (root / f"{pdb_id}_cleaned_ligands", base / "ligands_raw"),
            (root / f"{pdb_id.lower()}_cleaned_ligands", base / "ligands_raw"),
            (root / f"{pdb_idU}_CLEANED_LIGANDS", base / "ligands_raw"),
        ]
        for src, dst in legacy_dirs:
            if src.exists():
                logging.info("Migrating legacy directory %s -> %s", src, dst)
                _merge_dir(src, dst)

        candidates_files = [
            (
                root / f"{pdb_id}_nolig.pdb",
                base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb",
            ),
            (
                root / f"{pdb_id.lower()}_nolig.pdb",
                base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb",
            ),
        ]
        for src, dst in candidates_files:
            if src.exists() and not dst.exists():
                logging.info("Moving legacy file %s -> %s", src, dst)
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src), str(dst))
    except Exception as e:
        logging.warning("fold_legacy_layout (extended) failed for %s: %s", pdb_id, e)


def canon_paths(
    pdb_id: str,
    output_root: Union[str, Path],
    *,
    variant: Optional[str] = None,
) -> Dict[str, Path]:
    root = Path(output_root).resolve()
    base = root / pdb_id.upper()
    variant_token = _resolve_variant_token(prep_utils.config, variant)
    protein_root = base / variant_token if variant_token else base
    receptor_dir = protein_root / "receptor"
    paths = {
        "protein_root": protein_root,
        "raw": protein_root / "raw",
        "work": protein_root / "work",
        "ligands_raw": base / "ligands_raw",
        "nolig": protein_root / "nolig",
        "receptor": receptor_dir,
    }
    if variant_token:
        paths["variant"] = variant_token
        paths["variant_root"] = protein_root
    return paths
