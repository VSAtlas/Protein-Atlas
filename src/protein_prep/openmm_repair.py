"""OpenMM/PDBFixer structural repair helpers for receptor preparation."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Callable, Mapping


def _cfg_bool(cfg: Mapping[str, Any] | None, key: str, default: bool) -> bool:
    raw = os.environ.get(key)
    if raw is None and cfg is not None:
        raw = cfg.get(key)
    if raw is None:
        return default
    return str(raw).strip().lower() in {"1", "true", "yes", "y", "on"}


def _cfg_int(cfg: Mapping[str, Any] | None, key: str, default: int) -> int:
    raw = os.environ.get(key)
    if raw is None and cfg is not None:
        raw = cfg.get(key)
    try:
        return int(float(str(raw)))
    except Exception:
        return default


def get_target_ph_for_prep(pdb_path: str | Path) -> float:
    """Resolve the pH used by protein-prep protonation helpers."""
    raw_env = os.environ.get("TARGET_PH", "").strip()
    if raw_env:
        try:
            return float(raw_env)
        except Exception:
            logging.warning("[ph] invalid TARGET_PH=%r; falling back to context pH", raw_env)

    try:
        from path_router.context_ph import select_ph_values_for_protonation

        ph_values = select_ph_values_for_protonation(str(pdb_path))
        if ph_values:
            return float(ph_values[0])
    except Exception as exc:
        logging.info("[ph] context pH unavailable for %s: %s", pdb_path, exc)

    return 7.0


def _trim_missing_residues(fixer: Any, max_missing: int, logger: logging.Logger) -> None:
    missing = getattr(fixer, "missingResidues", None)
    if not isinstance(missing, dict) or not missing:
        return

    kept: dict[Any, Any] = {}
    dropped = 0
    for key, residues in missing.items():
        try:
            count = len(residues)
        except Exception:
            count = max_missing + 1
        if count <= max_missing:
            kept[key] = residues
        else:
            dropped += count

    if dropped:
        logger.warning(
            "[pdbfixer] skipped long missing-residue segments total_residues=%d max_segment=%d",
            dropped,
            max_missing,
        )
    fixer.missingResidues = kept


def _load_pdbfixer_api(
    logger: logging.Logger,
) -> tuple[type[Any], Callable[..., None]] | None:
    try:
        from pdbfixer import PDBFixer  # type: ignore[import-untyped]
        from openmm.app import PDBFile  # type: ignore[import-untyped]
    except Exception as exc:
        logger.warning("[pdbfixer] unavailable; skipping openmm repair: %s", exc)
        return None
    return PDBFixer, PDBFile.writeFile


def _configure_missing_residues(
    fixer: Any,
    cfg: Mapping[str, Any] | None,
    logger: logging.Logger,
) -> None:
    if not _cfg_bool(cfg, "PDBFIXER_ADD_MISSING_RESIDUES", False):
        fixer.missingResidues = {}
        return
    fixer.findMissingResidues()
    max_missing = max(0, _cfg_int(cfg, "PDBFIXER_MAX_MISSING_RESIDUES", 12))
    _trim_missing_residues(fixer, max_missing, logger)


def _apply_topology_repairs(
    fixer: Any,
    *,
    target_ph: float,
    cfg: Mapping[str, Any] | None,
    add_hydrogens: bool,
    logger: logging.Logger,
) -> None:
    _configure_missing_residues(fixer, cfg, logger)
    if _cfg_bool(cfg, "PDBFIXER_REPLACE_NONSTANDARD", True):
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(float(target_ph))


def _count_atoms(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        return sum(1 for line in handle if line.startswith(("ATOM  ", "HETATM")))


def repair_with_pdbfixer(
    input_pdb: str | Path,
    output_pdb: str | Path,
    *,
    target_ph: float,
    cfg: Mapping[str, Any] | None = None,
    logger: logging.Logger | None = None,
    add_hydrogens: bool = False,
) -> bool:
    """
    Repair receptor structure with open-source PDBFixer/OpenMM.

    This stage intentionally focuses on protein topology repair. Heterogen policy
    is handled before this stage by Atlas APO/HOLO logic, so heterogens are not
    globally deleted here.
    """
    log = logger or logging.getLogger(__name__)
    if not _cfg_bool(cfg, "PDBFIXER_REPAIR", True):
        log.info("[pdbfixer] disabled by PDBFIXER_REPAIR=0")
        return False

    api = _load_pdbfixer_api(log)
    if api is None:
        return False
    fixer_cls, write_pdb_file = api

    in_path = Path(input_pdb)
    out_path = Path(output_pdb)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        fixer = fixer_cls(filename=str(in_path))
        _apply_topology_repairs(
            fixer,
            target_ph=target_ph,
            cfg=cfg,
            add_hydrogens=add_hydrogens,
            logger=log,
        )
        with out_path.open("w", encoding="utf-8") as handle:
            write_pdb_file(fixer.topology, fixer.positions, handle, keepIds=True)

        atom_count = _count_atoms(out_path)
        if atom_count == 0:
            log.warning("[pdbfixer] wrote no atoms; ignoring output=%s", out_path)
            return False

        log.info(
            "[pdbfixer] repaired input=%s output=%s atoms=%d ph=%.2f add_h=%s",
            in_path,
            out_path,
            atom_count,
            float(target_ph),
            str(bool(add_hydrogens)).lower(),
        )
        return True
    except Exception as exc:
        log.warning("[pdbfixer] repair failed for %s: %s", in_path, exc)
        return False
