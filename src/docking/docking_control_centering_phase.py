from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .docking_control_redock import _control_centers_by_ph
from .docking_controls import (
    build_control_lookup,
    extract_ligands_to_nolig,
)
from path_router.path_router import Paths
from prep_ligands.prep_ligands_crystal import prep_ligands_from_pdb
from protein_functions import detect_active_site


def _phase5b_controls_and_control_redock(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    variant_token: Optional[str],
    variant_label: str,
    legacy_mode: bool,
    cleaned_pdb: Optional[str],
    receptor_pdbqt: Optional[str],
    center: Optional[Tuple[float, float, float]],
    box_size: Optional[Tuple[float, float, float]],
) -> Tuple[
    Optional[Tuple[float, float, float]],
    Optional[Tuple[float, float, float]],
    Dict[Optional[str], Tuple[float, float, float]],
    Dict[Optional[str], Tuple[float, float, float]],
    Dict[Optional[str], str],
    List[str],
    Dict[str, Path],
]:
    control_stems: List[str] = []
    control_lookup: Dict[str, Path] = {}
    fallback_center = center
    fallback_box = box_size

    if not receptor_pdbqt:
        logger.warning(
            "[control-redock] skip reason=missing_receptor_pdbqt pdb=%s variant=%s",
            paths.pdb_id,
            variant_label,
        )
        return fallback_center, fallback_box, {}, {}, {}, control_stems, control_lookup

    _lig_count, control_stems = extract_ligands_to_nolig(paths, logger)
    try:
        prep_ligands_from_pdb(
            ligand_output_dir=paths.ligand_output_dir,
            ligands_mol2_dir=paths.ligands_mol2_dir,
            prepped_ligands_dir=paths.prepped_ligands_dir,
        )
        logger.info("[controls] prepped_extracted=true")
    except Exception as e:
        logger.warning("[controls] prepped_extracted=false err=%s", e)

    ctrl_pdbqts: list[Path] = []
    for root in {paths.prepped_ligands_dir, Path(cfg["OUTPUT_LIGANDS_DIR"])}:
        if root.exists():
            ctrl_pdbqts.extend(root.glob("*.pdbqt"))

    logger.info("[controls] prepped_pdbqts=%d", len(ctrl_pdbqts))
    for p in ctrl_pdbqts[:10]:
        logger.info("[controls] prepped_pdbqt name=%s", p.name)

    control_lookup = build_control_lookup(paths)

    if fallback_center is None or fallback_box is None:
        if cleaned_pdb:
            c2, b2, src = detect_active_site(cleaned_pdb)
            if c2 and b2:
                fallback_center = tuple(float(x) for x in c2)
                box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
                fallback_box = tuple(min(box_cap, float(s)) for s in b2)
                logger.info(
                    "[active-site] fallback center=%s box=%s source=%s",
                    fallback_center,
                    fallback_box,
                    src or "activesite",
                )
        if fallback_center is None or fallback_box is None:
            logger.error(
                "[active-site] fallback_missing pdb=%s variant=%s",
                paths.pdb_id,
                variant_label,
            )

    center_by_ph, box_by_ph, source_by_ph = _control_centers_by_ph(
        cfg,
        paths,
        logger,
        variant_token=variant_token,
        legacy_mode=legacy_mode,
        cleaned_pdb=cleaned_pdb or "",
        fallback_center=fallback_center,
        fallback_box=fallback_box,
    )

    return (
        fallback_center,
        fallback_box,
        center_by_ph,
        box_by_ph,
        source_by_ph,
        control_stems,
        control_lookup,
    )
