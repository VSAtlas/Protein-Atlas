from __future__ import annotations

import difflib
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from docking_ligands import _count_heavy_atoms_from_pdbqt, _lib_roots_for_pdb
from input_and_export_functions import _to_bool
from prep_ligands import enumerate_ligands_for_docking
from single_ligand_index import (
    _ensure_single_ligand_index,
    _load_fda_name_map,
    _resolve_single_ligand,
)
from path_router import Paths


def resolve_single_ligand_or_prepare(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    pdb_id: str,
    run_mode: Optional[str],
) -> Optional[Tuple[List[str], Dict[str, int], Dict[str, Any]]]:
    """
    Handles single-ligand mode resolution and optional PH-ligand window priming.
    Returns (ligands, heavy_atom_counts, pains_flags) if single-ligand mode is active and resolved.
    Returns None if single-ligand mode is not active (or fallback was triggered).
    Raises SystemExit(2) if single-ligand is active but not found and fallback is disabled.
    """
    cfg.setdefault("_EFFECTIVE_SINGLE_LIGAND", "")
    single_ligand_hit: Optional[Path] = None
    cfg.pop("_SINGLE_RESOLVED_PATH", None)

    if cfg["_EFFECTIVE_SINGLE_LIGAND"]:
        _ensure_single_ligand_index(cfg, paths, logger)
        cfg.setdefault("paths", {})
        cfg["paths"]["prepped_ligands_dir"] = str(paths.prepped_ligands_dir)

        hit = _resolve_single_ligand(cfg["_EFFECTIVE_SINGLE_LIGAND"], pdb_id, cfg, logger)
        if hit:
            cfg["_SINGLE_RESOLVED_PATH"] = str(hit)
            single_ligand_hit = hit
        else:
            selector_token = cfg["_EFFECTIVE_SINGLE_LIGAND"]
            suggestions: list[str] = []
            try:
                fda_map = _load_fda_name_map(cfg, logger)
                suggestions = difflib.get_close_matches(
                    selector_token,
                    list(fda_map.keys()),
                    n=5,
                    cutoff=0.7,
                )
            except Exception:
                suggestions = []
            if suggestions:
                logger.error("[single.miss.suggest] did_you_mean=%s", ", ".join(suggestions))

            allow_flag = os.environ.get("ALLOW_FDA_FALLBACK")
            if allow_flag is None:
                allow_flag = cfg.get("ALLOW_FDA_FALLBACK", False)
            if not _to_bool(allow_flag):
                logger.error(
                    "[single.block] selector '%s' not found in fda_library via FDA_MAPPING_CSV; aborting instead of fallback.",
                    selector_token,
                )
                raise SystemExit(2)
            logger.warning(
                "[single.block] selector '%s' not found; ALLOW_FDA_FALLBACK enabled, continuing with fallback flow.",
                selector_token,
            )
            cfg["_EFFECTIVE_SINGLE_LIGAND"] = ""

    if cfg.get("_EFFECTIVE_SINGLE_LIGAND"):
        if not single_ligand_hit:
            return None
        ligands = [str(single_ligand_hit)]
        ha = _count_heavy_atoms_from_pdbqt(single_ligand_hit)
        heavy_atom_counts = {str(single_ligand_hit): ha}
        pains_flags = {}
        logger.info(f"[single] Active ? docking only: {single_ligand_hit.name} (heavy={ha})")

        if cfg.get("PH_LIGAND_MODE", "").lower() == "context_window" and cfg.get("PH_ENSEMBLE_IN_PREP"):
            try:
                ph_values = [6.0, 8.0]
                if "_PH_CONTEXT_VALUES" in cfg:
                    ph_values = cfg["_PH_CONTEXT_VALUES"]
                ligand_window = sorted(
                    {round(p, 1) for ph in ph_values for p in (float(ph) - 1.0, float(ph), float(ph) + 1.0)}
                )
                logger.info(f"[single.ph_ligand] Using ligand window {ligand_window}")

                if run_mode in {"dud", "hmdb"}:
                    ph_override_single = run_mode
                elif run_mode == "fda":
                    ph_override_single = "off"
                else:
                    ph_override_single = None
                _, noncontrol_roots_single = _lib_roots_for_pdb(
                    cfg,
                    paths.pdb_id.upper(),
                    paths,
                    logger,
                    test_mode_override=ph_override_single,
                )
                ph_root_path = noncontrol_roots_single[0] if noncontrol_roots_single else None
                ph_root_cfg = str(ph_root_path) if ph_root_path else ""

                logger.info(
                    "[single.ph_ligand.bridge] ph_root_cfg=%s ph_root_path=%s exists=%s",
                    ph_root_cfg,
                    str(ph_root_path) if ph_root_path is not None else "",
                    ph_root_path.exists() if ph_root_path is not None else False,
                )
                if ph_root_path is not None and ph_root_path.exists():
                    enumerate_ligands_for_docking(
                        requested_ph_values=ligand_window,
                        root_dir=ph_root_path,
                        microstate_dedup=True,
                        force=False,
                    )
                else:
                    logger.info(
                        "[single.ph_ligand.bridge.skip] no valid ph_ligand_root; "
                        "skipping microstate priming for single-ligand mode"
                    )
            except Exception as e:
                logger.warning(f"[single.ph_ligand.skip] Could not run PH-ligand window for single mode: {e}")

        return ligands, heavy_atom_counts, pains_flags

    return None
