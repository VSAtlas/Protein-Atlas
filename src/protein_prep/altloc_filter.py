from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Union


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
    """Filter alternate locations while preserving record/model order.

    NMR PDBs commonly contain many MODEL/ENDMDL blocks with the same chain,
    residue, and atom identifiers in each model.  A global atom-key sort
    collapses those models and moves all atoms outside their MODEL blocks, so
    this routine selects altLocs per model and then streams records back in the
    original order.
    """
    input_lines = Path(pdb_input_path).read_text(
        encoding="utf-8", errors="ignore"
    ).splitlines(keepends=True)
    choices: dict[tuple[int, str, str, str, str, str], str] = {}
    altlocs_by_key: dict[tuple[int, str, str, str, str, str], set[str]] = {}
    model_index = 0
    for line in input_lines:
        if line.startswith("MODEL"):
            model_index += 1
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = (
            model_index,
            line[0:6],
            line[21],
            line[22:26],
            line[26],
            line[12:16].strip(),
        )
        altlocs_by_key.setdefault(key, set()).add(line[16])

    removed_count = 0
    for key, altlocs in altlocs_by_key.items():
        if " " in altlocs:
            selected = " "
        elif "A" in altlocs:
            selected = "A"
        else:
            selected = sorted(altlocs)[0]
        choices[key] = selected
        if len(altlocs) > 1:
            removed = sorted(alt for alt in altlocs if alt != selected)
            removed_count += len(removed)
            logging.info(
                "Keeping altLoc '%s' for atom %s, removed %s", selected, key, removed
            )

    model_index = 0
    with open(pdb_output_path, "w", encoding="utf-8") as f:
        for line in input_lines:
            if line.startswith("MODEL"):
                model_index += 1
            if not line.startswith(("ATOM  ", "HETATM")):
                f.write(line)
                continue
            key = (
                model_index,
                line[0:6],
                line[21],
                line[22:26],
                line[26],
                line[12:16].strip(),
            )
            selected = choices.get(key, line[16])
            if line[16] != selected:
                continue
            if selected != " ":
                line = line[:16] + " " + line[17:]
            f.write(line)
    _maybe_post_write_element_guard("altloc_filter", pdb_output_path)
    logging.info(
        "Filtered altLocs in %s → %s (removed %d alternates)",
        pdb_input_path,
        pdb_output_path,
        removed_count,
    )
