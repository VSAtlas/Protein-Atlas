"""Minimal RDKit SDF serialization (keeps heavy docking pipeline off import paths)."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def write_rdkit_mol_to_sdf(
    mol: Any, sdf_path: Path, *, set_kekulize_false: bool = True
) -> bool:
    from rdkit import Chem

    writer = Chem.SDWriter(str(sdf_path))
    if set_kekulize_false:
        try:
            writer.SetKekulize(False)
        except (AttributeError, TypeError):
            pass
    writer.write(mol)
    writer.close()
    return sdf_path.exists() and sdf_path.stat().st_size > 0
