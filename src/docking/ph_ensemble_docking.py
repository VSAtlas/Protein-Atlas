from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from path_router import load_ph_tags, receptor_file

"""
Helpers for pH-ensemble docking orchestration.

This module centralizes:
- pH-tag enumeration and manifest bookkeeping
- pH-ligand 'context window' microstate prewarm
- Per-ph-label ligand enumeration for PH_LIGAND_MODE=context_window

It intentionally does NOT contain:
- CenterSelector logic
- Multi-stage docking orchestration
- Checkpointing/resume logic
"""


def _parse_ph_values_from_label(label: Optional[str]) -> list[float]:
    """
    Parse one of our PH ensemble labels into a list of numeric pH values.

    Handles simple tags like:
        "pH7_0", "pH8_4"
    and composite tags like:
        "pH7_9+8_4+8_9-dup19"

    Returns a list of floats (e.g. [7.9, 8.4, 8.9]) or an empty list
    if no numeric pH values can be parsed.
    """
    vals: list[float] = []
    if not label:
        return vals

    s = str(label).strip()
    if not s:
        return vals

    # Strip leading "pH" (case-insensitive)
    m = re.search(r"(?i)pH(.+)", s)
    if m:
        s = m.group(1)

    # Drop any suffix after first "-" (e.g. "-dup19")
    if "-" in s:
        s = s.split("-", 1)[0]

    # Split on "+", convert pieces like "7_9" -> 7.9
    for part in s.split("+"):
        part = part.strip()
        if not part:
            continue
        part = part.replace("_", ".")
        try:
            ph = float(part)
        except Exception:
            continue
        # Sanity range for pH values
        if 0.0 < ph < 15.0:
            vals.append(ph)

    return vals


def init_ph_tags_and_manifest(
    cfg: Dict,
    pdb_id: str,
    variant_token: Optional[str],
    legacy_mode: bool,
) -> List[Optional[str]]:
    """
    Initialize pH tags for a given protein/variant and update the
    cfg['_PH_ENSEMBLE_CANONICAL'] manifest, mirroring the existing logic.

    - If PH_ENSEMBLE is true:
        * Log debug map_keys/map_len
        * Use load_ph_tags(...) to get ph_tags
        * If ph_tags is empty: log a warning and return []
        * Else populate cfg['_PH_ENSEMBLE_CANONICAL'][pdb_id]
    - If PH_ENSEMBLE is false:
        * Return [None] to indicate non-ensemble (legacy) mode.
    """
    ph_log = logging.getLogger("ph_ensemble")
    ph_enabled = bool(cfg.get("PH_ENSEMBLE"))

    if ph_enabled:
        try:
            _keys = sorted(list((cfg.get("_PH_ENSEMBLE_CANONICAL") or {}).keys()))
            _len_here = len((cfg.get("_PH_ENSEMBLE_CANONICAL") or {}).get(pdb_id, []))
            ph_log.info("[ph_ensemble.debug] map_keys=%s map_len[%s]=%d", ",".join(_keys), pdb_id, _len_here)
        except Exception:
            pass

        ph_tags = load_ph_tags(pdb_id, variant=variant_token)
        if not ph_tags:
            ph_log.warning(
                "[router.warn] ensemble.json had no members; skipping pdb=%s variant=%s",
                pdb_id,
                variant_token or "legacy",
            )
            return []

        canonical: list[tuple[str, str]] = []
        for tag in ph_tags:
            rec_path = receptor_file(pdb_id, variant=variant_token, ph_tag=tag, legacy=legacy_mode)
            canonical.append((tag, str(rec_path)))
        if canonical:
            cfg.setdefault("_PH_ENSEMBLE_CANONICAL", {})[pdb_id] = canonical

        return ph_tags
    else:
        return [None]


def prewarm_ph_ligand_microstates(
    cfg: Dict,
    ph_tags: List[Optional[str]],
    ph_ligand_root: Optional[Path],
) -> None:
    """
    Prewarm the pH-ligand microstate registry for the union of all pH tags.

    Mirrors the existing 'context window' block in main.py:
    - If PH_LIGAND_MODE != 'context_window' or PH_ENSEMBLE is false, do nothing.
    - Collect numeric pH values from ph_tags via _parse_ph_values_from_label.
    - Build a +-1.0 pH 'window' around each and deduplicate/round to 0.1.
    - Call prep_ligands.enumerate_ligands_for_docking(...) once on the union.
    - Log the same [ph_ligand.*] messages and swallow exceptions gracefully.
    """
    if cfg.get("PH_LIGAND_MODE", "").lower() != "context_window" or not cfg.get("PH_ENSEMBLE"):
        return

    ph_log = logging.getLogger("ph_ensemble")

    try:
        from prep_ligands.prep_ligands_microstates import enumerate_ligands_for_docking

        ph_log.info("[ph_ligand.bridge] active: PH_LIGAND_MODE=context_window")

        context_pHs: list[float] = []
        for tag in ph_tags:
            for pH_num in _parse_ph_values_from_label(tag):
                context_pHs.append(pH_num)

        if not context_pHs:
            ph_log.warning(
                "[ph_ligand.prep] no numeric pH values parsed from tags=%s; skipping ligand prep window",
                ph_tags,
            )
            return

        ligand_window = sorted(
            {
                round(x, 1)
                for x in [p for ph in context_pHs for p in (ph - 1.0, ph, ph + 1.0)]
            }
        )
        ph_log.info(
            "[ph_ligand.prep] preparing ligands for window %s from tags=%s",
            ligand_window,
            ph_tags,
        )

        enumerate_ligands_for_docking(
            requested_ph_values=ligand_window,
            microstate_dedup=True,
            force=False,
            root_dir=ph_ligand_root,
        )
    except Exception as e:
        ph_log.warning(f"[ph_ligand.prep.skip] failed to initialize: {e}")


def enumerate_ligands_for_ph_context(
    cfg: Dict,
    pdb_id: str,
    ph_label: Optional[str],
    ph_ligand_root: Optional[Path],
) -> Optional[List[Path]]:
    """
    For a given ph_label, determine the pH 'window' and call
    prep_ligands.enumerate_ligands_for_docking(...) to get a
    pH-specific set of ligands.

    - If PH_LIGAND_MODE != 'context_window', return None.
    - If ph_label cannot be parsed into numeric pH values, log and return None.
    - If enumeration returns an empty list, log and return None.
    - On success, log [ph_ligand.context] and [ph_ligand.selected] style messages
      and return the list of Path objects.

    The caller (main.py) remains responsible for:
    - Building heavy_atom_counts with _count_heavy_atoms_from_pdbqt
    - Updating pains_flags
    - Falling back to base ligands when this returns None.
    """
    if cfg.get("PH_LIGAND_MODE", "").lower() != "context_window":
        return None

    ph_log = logging.getLogger("ph_ensemble")

    try:
        from prep_ligands.prep_ligands_microstates import enumerate_ligands_for_docking

        ph_values = _parse_ph_values_from_label(ph_label)
        if not ph_values:
            ph_log.warning(
                "[ph_ligand.context] unable to parse numeric pH from ph_label=%s; using base ligands",
                ph_label,
            )
            return None

        ligand_window = sorted(
            {
                round(x, 1)
                for x in [p for ph in ph_values for p in (ph - 1.0, ph, ph + 1.0)]
            }
        )
        ph_log.info(
            "[ph_ligand.context] ph_label=%s values=%s -> ligand window=%s",
            ph_label,
            ph_values,
            ligand_window,
        )

        enumerated = enumerate_ligands_for_docking(
            requested_ph_values=ligand_window,
            microstate_dedup=True,
            force=False,
            root_dir=ph_ligand_root,
        )

        if not enumerated:
            ph_log.warning(
                "[ph_ligand.empty] pdb_id=%s ph=%s window=%s -> no microstates; falling back to base ligands",
                pdb_id,
                ph_label,
                ligand_window,
            )
            return None

        ph_log.info(
            "[ph_ligand.selected] pdb_id=%s ph=%s ligands=%d",
            pdb_id,
            ph_label,
            len(enumerated),
        )
        return list(enumerated)
    except Exception as e:
        ph_log.warning(
            "[ph_ligand.context.skip] failed during pH-specific ligand enumeration: %s",
            e,
        )
        return None
