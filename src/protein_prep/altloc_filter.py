from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import List, Union


def _maybe_post_write_element_guard(step: str, pdb_path: Union[str, Path]) -> None:
    """
    Keep existing guard behavior when called from automate_protein_prep without
    importing that heavy module in standalone/unit-test contexts.
    """
    automate_mod = sys.modules.get("automate_protein_prep")
    guard = getattr(automate_mod, "_post_write_element_guard", None)
    if callable(guard):
        guard(step, pdb_path)


def filter_altlocs(
    pdb_input_path: Union[str, Path], pdb_output_path: Union[str, Path]
) -> None:
    """Filter alternate locations (altLoc) deterministically and log decisions."""
    lines: List[str] = []
    atoms = {}
    removed_count = 0

    with open(pdb_input_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_name = line[12:16]
                altLoc = line[16]
                chainID = line[21]
                resSeq = line[22:26].strip()
                iCode = line[26]
                key = (chainID, resSeq, iCode, atom_name.strip())
                atoms.setdefault(key, {})[altLoc] = line
            else:
                lines.append(line)

    filtered_atoms: List[str] = []
    for key, altloc_dict in atoms.items():
        altLocs = list(altloc_dict.keys())
        if len(altLocs) > 1:
            logging.info("AltLocs found for %s: %s", key, altLocs)
        if " " in altloc_dict:
            selected = " "
        elif "A" in altloc_dict:
            selected = "A"
        else:
            selected = sorted(altloc_dict.keys())[0]
        if len(altloc_dict) > 1:
            removed = [alt for alt in altLocs if alt != selected]
            removed_count += len(removed)
            logging.info(
                "Keeping altLoc '%s' for atom %s, removed %s", selected, key, removed
            )
        filtered_atoms.append(altloc_dict[selected])

    def _int_safe(s: str) -> int:
        s = s.strip()
        return int(s) if s and s.lstrip("-").isdigit() else 0

    filtered_atoms.sort(
        key=lambda l: (l[21], _int_safe(l[22:26]), l[26], l[12:16].strip())
    )

    with open(pdb_output_path, "w", encoding="utf-8") as f:
        for line in lines:
            f.write(line)
        for atom_line in filtered_atoms:
            f.write(atom_line)
    _maybe_post_write_element_guard("altloc_filter", pdb_output_path)
    logging.info(
        "Filtered altLocs in %s → %s (removed %d alternates)",
        pdb_input_path,
        pdb_output_path,
        removed_count,
    )
