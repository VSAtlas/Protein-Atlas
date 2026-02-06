"""Hydrogen cleanup utilities for receptor-prep PDB files."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Tuple, Union

from protein_prep.element_guard import _post_write_element_guard
from protein_prep.protonation import conect_coverage


def remove_unbonded_atoms(pdb_path: Union[str, Path]) -> None:
    bonded_atoms = set()
    all_atoms: List[str] = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("CONECT"):
                parts = line.split()
                for atom_serial in parts[1:]:
                    bonded_atoms.add(atom_serial)
            elif line.startswith(("ATOM", "HETATM")):
                all_atoms.append(line)

    filtered: List[str] = []
    for line in all_atoms:
        atom_serial = line[6:11].strip()
        if atom_serial in bonded_atoms or line[76:78].strip() != "H":
            filtered.append(line)
        else:
            logging.info("Removed unbonded hydrogen: %s", line.strip())

    with open(pdb_path, "w", encoding="utf-8") as f:
        f.writelines(filtered)


def remove_implausible_hydrogens_by_distance(pdb_path: Union[str, Path]) -> None:
    atoms: List[Tuple[str, str, Optional[Tuple[float, float, float]]]] = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    atoms.append((line, "UNK", (None, None, None)))
                    continue
                el = line[76:78].strip() or line[12:16].strip()[:1]
                atoms.append((line, el.upper(), (x, y, z)))

    kept: List[str] = []
    heavy_coords = [a[2] for a in atoms if a[1] != "H" and a[2][0] is not None]

    def near_heavy(coord: Optional[Tuple[float, float, float]]) -> bool:
        if coord[0] is None:
            return True
        x, y, z = coord
        for X, Y, Z in heavy_coords:
            dx = x - X
            dy = y - Y
            dz = z - Z
            if (dx * dx + dy * dy + dz * dz) <= (1.35 * 1.35):
                return True
        return False

    for line, el, coord in atoms:
        if el != "H" or near_heavy(coord):
            kept.append(line)
        else:
            logging.info("Removed implausible H: %s", line.strip())

    with open(pdb_path, "w", encoding="utf-8") as out:
        out.writelines(kept)


def clean_hydrogens(
    pdb_path: Union[str, Path],
    use_conect_if_reliable: bool = True,
    conect_min_cov: float = 0.6,
) -> None:
    cov = conect_coverage(pdb_path) if use_conect_if_reliable else 0.0
    if cov >= conect_min_cov:
        remove_unbonded_atoms(pdb_path)
    remove_implausible_hydrogens_by_distance(pdb_path)
    _post_write_element_guard("hydrogen_cleanup", pdb_path)
