from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from docking.docking_utils import norm
from docking.docking_ligands import select_ligands_for_next
from docking.docking_ledock import should_run_ledock_for_target
from docking.docking_dock6 import should_run_dock6_for_target


def _apply_force_carry_and_doping(
    cfg: Dict[str, Any],
    docking_mode: str,
    stage_index: int,
    stages_for_run: List[Dict[str, Any]],
    scores: Dict[str, float],
    logger: logging.Logger,
    *,
    stage1_original: List[str],
    forced_extracted_for_stage3: Optional[set] = None,
    invalids: Optional[Dict[str, Tuple[Optional[float], str]]] = None,
    next_stage_is_stage3: bool = False,
    control_norms: Optional[set[str]] = None,
) -> List[str]:
    """
    Mirror the existing Vina selection + doping flow:
      - percentile selection via select_ligands_for_next
      - optional force-carry of extracted controls into stage3
        (filtered back out here; controls are injected at Stage3 run time)
      - optional rescue of self-RMSD near-miss ligands
    """
    use_stage1_base = docking_mode == "polypharmacology" and stage_index == 1

    rescue: List[str] = []
    if invalids and stage_index < len(stages_for_run) - 1:
        for lig, (sc, reason) in invalids.items():
            if (
                sc is not None
                and "self_rmsd_" in str(reason).lower()
                and sc <= float(cfg.get("RESCUE_SELF_RMSD_SCORE_MAX", -8.0))
            ):
                rescue.append((sc, lig))
        rescue = [
            lig
            for _, lig in sorted(rescue)[: int(cfg.get("RESCUE_SELF_RMSD_TOP_N", 10))]
        ]

    selected_raw = select_ligands_for_next(
        docking_mode,
        stage_index,
        stages_for_run,
        scores,
        logger,
        base_pool_n=(len(stage1_original) if use_stage1_base else None),
        force_include=(forced_extracted_for_stage3 if next_stage_is_stage3 else None),
    )

    def _filter_controls(seq: List[str]) -> List[str]:
        if not control_norms:
            return list(seq)
        filt = []
        for lig in seq:
            if norm(lig) in control_norms:
                continue
            filt.append(lig)
        return filt

    selected = _filter_controls(selected_raw)
    removed_selected = len(selected_raw) - len(selected)
    if removed_selected and next_stage_is_stage3:
        logger.info(
            "[Force-carry] Dropped %d control ligand(s) from selection; controls are merged at Stage3 docking time.",
            removed_selected,
        )

    if rescue:
        sel_set = set(selected)
        rescue_unique = [r for r in _filter_controls(rescue) if r not in sel_set]
        ligands = rescue_unique + selected
    else:
        ligands = selected
    return ligands


def _use_ledock(cfg: Dict[str, Any]) -> bool:
    try:
        return should_run_ledock_for_target(cfg)
    except Exception:
        return False


def _use_dock6(cfg: Dict[str, Any]) -> bool:
    try:
        return should_run_dock6_for_target(cfg)
    except Exception:
        return False


def _is_stage3(stage_name: str) -> bool:
    s = str(stage_name).lower()
    return re.search(r"(?:^|_)stage3(?:$|_)", s) is not None


def _split_controls_and_noncontrols(
    ligands: List[str],
    prepped_root: Optional[Path],
    *,
    allowed_subdirs: Optional[set[str]] = None,
    control_stems: Optional[set[str]] = None,
) -> tuple[list[str], list[str]]:
    controls: list[str] = []
    noncontrols: list[str] = []
    allow = allowed_subdirs or {"controls", "reference"}
    ctrl_stems_lower = {s.lower() for s in (control_stems or set())}

    def _is_control_path(p: Path) -> bool:
        if not prepped_root:
            return False
        try:
            rel = p.resolve().relative_to(prepped_root.resolve())
        except Exception:
            return False
        if not rel.parts:
            return False
        return rel.parts[0] in allow

    def _is_control_stem(p: Path) -> bool:
        stem = p.stem.split("_stage")[0].lower()
        return stem in ctrl_stems_lower

    for lig in ligands:
        p = Path(lig)
        if _is_control_path(p) or _is_control_stem(p):
            controls.append(str(p))
        else:
            noncontrols.append(str(p))

    return controls, noncontrols


def _interleave_controls(noncontrols: List[str], controls: List[str]) -> List[str]:
    merged: list[str] = []
    i = 0
    max_len = max(len(noncontrols), len(controls))
    while i < max_len:
        if i < len(noncontrols):
            merged.append(noncontrols[i])
        if i < len(controls):
            merged.append(controls[i])
        i += 1
    return merged
