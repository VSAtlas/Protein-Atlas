from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from docking.docking_ligands import _count_heavy_atoms_from_pdbqt
from docking.library_mode import compute_allowed_library_roots, parse_test_libraries
from path_router.path_router import Paths, load_ph_tags
from docking.ph_ensemble_docking import (
    enumerate_ligands_for_ph_context,
    init_ph_tags_and_manifest,
    prewarm_ph_ligand_microstates,
)


def resolve_ph_tags_and_root(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    run_mode: Optional[str],
    variant_token: Optional[str],
    legacy_mode: bool,
) -> Tuple[List[str], Optional[Path], Optional[str], Optional[str]]:
    """
    Resolves PH tags, ligand roots, and manifest library info for the subrun.
    Returns (ph_tags, ph_ligand_root, ph_test_mode_override, library_for_manifest).
    """
    ph_enabled = bool(cfg.get("PH_ENSEMBLE"))
    manifest_run_id = cfg.get("RUN_ID")
    library_for_manifest = None
    try:
        tokens_for_manifest = parse_test_libraries(cfg)
        if "dud" in tokens_for_manifest:
            lib_map = cfg.get("_TEST_LIBRARY_CANONICAL", {}) or {}
            library_for_manifest = lib_map.get(paths.pdb_id.upper())
        if not library_for_manifest:
            library_for_manifest = cfg.get("LIBRARY_SUBDIR_DEFAULT")
    except Exception:
        library_for_manifest = cfg.get("LIBRARY_SUBDIR_DEFAULT")

    ph_tags = init_ph_tags_and_manifest(cfg, paths.pdb_id, variant_token, legacy_mode)
    ph_tags = [
        str(tag).strip() for tag in ph_tags if tag is not None and str(tag).strip()
    ]
    if not ph_tags:
        fallback_ph = load_ph_tags(paths.pdb_id, variant=variant_token) or []
        ph_tags = [
            str(tag).strip()
            for tag in fallback_ph
            if tag is not None and str(tag).strip()
        ]
    if not ph_tags:
        ph_tags = ["base"]  # ensure downstream logging/manifest updates occur

    if ph_enabled and not ph_tags:
        logger.info(
            "[subrun.ph] run_mode=%s ph_enabled=True but no ph_tags; skipping PH run",
            run_mode or "None",
        )
        # We return empty ph_tags to signal skip, but caller must handle it.
        # But wait, original code returns early.
        # Here we return empty list, caller checks.
        pass

    tokens_override: Optional[list[str]] = None
    if run_mode is not None:
        token = str(run_mode).strip()
        if token:
            lowered = token.lower()
            if lowered in {"default", "off", "none", "null"}:
                lowered = "fda"
            tokens_override = [lowered]

    noncontrol_roots_ph = compute_allowed_library_roots(
        cfg,
        paths.pdb_id.upper(),
        logger,
        tokens_override=tokens_override,
    )
    ph_ligand_root = noncontrol_roots_ph[0] if noncontrol_roots_ph else None
    ph_test_mode_override = "+".join(tokens_override) if tokens_override else None

    logger.info(
        "[subrun.ph] run_mode=%s ph_tags=%s ph_ligand_root=%s override=%s",
        run_mode or "None",
        ",".join(ph_tags) if ph_tags else "(none)",
        str(ph_ligand_root) if ph_ligand_root else "(none)",
        ph_test_mode_override or "(none)",
    )
    prewarm_ph_ligand_microstates(cfg, ph_tags, ph_ligand_root)

    return ph_tags, ph_ligand_root, ph_test_mode_override, library_for_manifest


def enumerate_ph_ligands_if_needed(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    ph_label: str,
    ph_ligand_root: Optional[Path],
    base_ligands: List[str],
    base_heavy_atoms: Dict[str, int],
    base_pains_flags: Dict[str, Any],
    controls_for_run: List[str],
) -> Tuple[List[str], Dict[str, int], Dict[str, Any]]:
    """
    Enumerates ligands for a specific PH context, merging with base stats.
    Returns (ligands, heavy_atom_counts, pains_flags).
    """
    enumerated = enumerate_ligands_for_ph_context(
        cfg=cfg,
        pdb_id=paths.pdb_id,
        ph_label=ph_label,
        ph_ligand_root=ph_ligand_root,
    )

    if enumerated:
        ligands = [str(p) for p in enumerated]
        heavy_atom_counts = {
            str(p): _count_heavy_atoms_from_pdbqt(p) for p in enumerated
        }
        pains_flags = {
            k: base_pains_flags.get(
                k,
                base_pains_flags.get(Path(k).stem, False),
            )
            for k in ligands
        }
        # Merge controls back in to counts/pains if needed
        for c in controls_for_run:
            if c not in heavy_atom_counts and c in base_heavy_atoms:
                heavy_atom_counts[c] = base_heavy_atoms[c]
            pains_flags.setdefault(
                c,
                base_pains_flags.get(c, base_pains_flags.get(Path(c).stem, False)),
            )
    else:
        ligands = base_ligands[:]
        heavy_atom_counts = dict(base_heavy_atoms)
        pains_flags = dict(base_pains_flags)

    return ligands, heavy_atom_counts, pains_flags
