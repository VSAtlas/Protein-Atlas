# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .docking_stage_runner import RetryManager, run_one_stage
from .docking_subrun_selection import (
    _apply_force_carry_and_doping,
    _interleave_controls,
    _is_stage3,
    _use_ledock,
    _use_dock6,
)
from .docking_centering import CenterSelector
from .docking_utils import (
    _fingerprint_stage,
    early_recenter_decision,
    run_completion_audit,
    norm,
)
from .fallback_recenter import (
    GlobalCenterGuard,
    RecenterParams,
    fallback_recentering_if_empty,
)
from checkpoints import (
    checkpoint_should_skip,
    checkpoint_mark_done,
    checkpoint_invalidate_from,
)
from druggability_orchestrator import decide_engine_policy
from run_manifest import (
    update_manifest_for_druggability_and_engine_plan,
)
from .docking_vina import emit_vina_config
from .run_vina import run_docking_task
from input_and_export_functions import record_score
from record_data import record_le


@dataclass
class MultistageVinaResult:
    scores: Dict[str, float]
    validated: List[str]
    distances: List[float]
    raw_docked: Dict[str, str]
    invalids: Dict[str, Tuple[Optional[float], str]]
    validated_ligands_last: List[str]
    forced_extracted_for_stage3: set
    recenter_attempts: int
    gnina_jobs: List[Dict[str, Any]]
    ledock_jobs: List[Dict[str, Any]]
    dock6_jobs: List[Dict[str, Any]]
    dock6_metrics_by_stage: Dict[str, Dict[Path, Dict[str, Any]]]
    score_history: Dict[str, Dict[str, Dict]]
    docking_mode: str
    switch_history: List[Any]
    global_switches: int


def run_multistage_vina(
    *,
    cfg: Dict[str, Any],
    paths: Any,
    receptor_pdbqt: str,
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    ph_label: Optional[str],
    legacy_mode: bool,
    variant_env: Optional[str],
    stages_for_run: List[Dict[str, Any]],
    stage1_original: List[str],
    controls_for_run: List[str],
    control_norms_for_run: set,
    control_lookup: Dict[str, Path],
    heavy_atom_counts: Dict[str, int],
    pains_flags: Dict[str, Any],
    logger: logging.Logger,
    recenter_params: RecenterParams,
    control_stems: List[str],
) -> MultistageVinaResult:
    
    forced_extracted_for_stage3 = set(controls_for_run)
    logger.info(
        "[Force-carry] Stage3 control pool size=%d noncontrols_stage1=%d",
        len(forced_extracted_for_stage3),
        len(stage1_original),
    )

    score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
    dock6_metrics_by_stage: Dict[str, Dict[Path, Dict[str, Any]]] = defaultdict(dict)
    validated_ligands_last: List[str] = []
    recenter_attempts = 0
    docking_mode = cfg.get("DOCKING_MODE", "discovery").lower()

    retry_mgr = RetryManager()
    gnina_jobs: List[Dict[str, Any]] = []
    ledock_jobs: List[Dict[str, Any]] = []
    dock6_jobs: List[Dict[str, Any]] = []
    ledock_enabled = _use_ledock(cfg)
    dock6_enabled = _use_dock6(cfg)

    ligands = list(stage1_original)
    
    ctrl_stems_lower = {s.lower() for s in control_stems}
    ctrl_blacklist = {t.strip().upper() for t in str(cfg.get("CONTROL_BLACKLIST", "")).split(",") if t.strip()}
    min_ha = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))

    selector = CenterSelector(cfg, logger, control_stems, heavy_atom_counts, center)
    guard = GlobalCenterGuard(
        max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2))
    )
    
    scores: Dict[str, float] = {}
    validated: List[str] = []
    distances: List[float] = []
    raw_docked: Dict[str, str] = {}
    invalids: Dict[str, Tuple[Optional[float], str]] = {}
    
    variant_label = variant_env or "legacy"

    i = 0
    while i < len(stages_for_run):
        guard.reset_stage()
        stage = stages_for_run[i]
        stage_is_stage3 = _is_stage3(stage["name"])
        stage_controls = sorted(controls_for_run, key=lambda p: Path(p).name) if stage_is_stage3 else []
        stage_noncontrols = [
            l
            for l in ligands
            if norm(l) not in control_norms_for_run
            and Path(l).stem.split("_stage")[0].lower() not in ctrl_stems_lower
        ]
        stage_ligands = (
            _interleave_controls(stage_noncontrols, stage_controls)
            if stage_is_stage3
            else stage_noncontrols
        )

        if bool(cfg.get("CHECKPOINT_ENABLE", True)):
            fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
            if checkpoint_should_skip(
                cfg,
                paths.pdb_id,
                stage["name"],
                fp,
                ph_label=ph_label,
                variant=variant_env or None,
                engine="vina",
            ):
                logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                i += 1
                continue

        if not stage_ligands:
            logger.warning(f"No ligands to dock at {stage['name']}; skipping this stage.")
            i += 1
            continue

        stage_ligands_for_audit = [norm(l) for l in stage_ligands]
        stage_dir = paths.docked_stage_dir(variant_env or None, stage["name"], ph_label)

        logger.info(
            "[ligand.stage] stage=%s ph=%s controls=%d noncontrols=%d docking_ligands=%d",
            stage["name"],
            ph_label if ph_label else "base",
            len(stage_controls),
            len(stage_noncontrols),
            len(stage_ligands),
        )
        
        ph_log = logging.getLogger("ph_ensemble")
        if ph_label:
            ph_log.info(
                "[ph_ensemble.dock.stage] pdb_id=%s variant=%s ph=%s stage=%s receptor=%s out=%s",
                paths.pdb_id,
                variant_label,
                ph_label,
                stage["name"],
                receptor_pdbqt,
                str(stage_dir),
            )

        scores, validated, distances, raw_docked, invalids = run_one_stage(
            cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
            stage_ligands, logger, retry_mgr, control_lookup, ph_label=ph_label
        )

        validated_ligands_last = validated

        ctrl_blacklist_set = ctrl_blacklist

        def _is_control(lig: str) -> bool:
            lig_norm = norm(lig)
            stem = Path(lig).stem.split("_stage")[0].lower()
            if stem.upper() in ctrl_blacklist_set:
                return False
            if lig_norm not in control_norms_for_run and stem not in ctrl_stems_lower:
                return False
            ha = heavy_atom_counts.get(lig)
            if ha is not None and ha < min_ha:
                return False
            return True

        control_anchor_hit = any(_is_control(lig) for lig in validated)

        lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", float("inf")))
        lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
        lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))

        qualified_controls = []
        for lig in validated:
            if not _is_control(lig):
                continue
            sc = scores.get(lig)
            if sc is None or not np.isfinite(sc):
                continue
            if sc > lock_score_max:
                continue
            pose_path = raw_docked.get(lig)
            if not pose_path:
                continue
            c = CenterSelector._pdbqt_centroid(pose_path)
            if c is None:
                continue
            if np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                qualified_controls.append(lig)

        if len(qualified_controls) >= lock_min_hits and not guard.locked:
            guard.lock()
            logger.info(
                "[CONTROL-LOCK] Control(s) validated with strong confidence "
                f"(n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Ang); "
                "center is now anchored; future center switches are disabled."
            )

        try:
            processed = {norm(x) for x in stage_ligands}
            valid_set = {norm(x) for x in scores.keys()}
            invalid_set = {norm(x) for x in invalids.keys()}
            both = valid_set & invalid_set
            missing = processed - (valid_set | invalid_set)
            if both or missing:
                logger.error(
                    f"Invariant violation at {stage['name']}: both={len(both)}, missing={len(missing)}"
                )
                if both:
                    logger.error(
                        "Ligands marked both valid & invalid: "
                        + ", ".join(os.path.basename(x) for x in list(both)[:10])
                    )
                if missing:
                    logger.error(
                        "Ligands missing from results: "
                        + ", ".join(os.path.basename(x) for x in list(missing)[:10])
                    )
                    for lig_m in missing:
                        invalids[lig_m] = (None, "not_processed")
        except Exception as _e:
            logger.warning(f"Invariant check failed: {_e}")

        for lig, sc in scores.items():
            record_score(score_history, stage['name'], lig, sc, True)
            record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)
        for lig, (sc, reason) in invalids.items():
            record_score(score_history, stage['name'], lig, sc, False, reason=reason)
            record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)

        promoted_this_stage = False
        try:
            decision = selector.consider_switch(
                stage['name'], scores, validated, raw_docked, receptor_pdbqt, center, guard
            )
            if decision.promoted and decision.new_center is not None:
                old = center
                center = decision.new_center
                promoted_this_stage = True
                guard.mark_switch()
                if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                    checkpoint_invalidate_from(
                        cfg,
                        paths.pdb_id,
                        stages_for_run,
                        start_index=i,
                        ph_label=ph_label,
                        variant=variant_env or None,
                    )
                logger.info(
                    f"[CENTER] Switched from {old} -> {center} ({decision.reason}, "
                    f"SwitchScore={decision.switchscore:.2f}) [global switch]"
                )
        except Exception as e:
            logger.warning(f"CenterSelector failed gracefully: {e}")

        if ph_label:
            ph_log.info(
                "[ph_ensemble.dock.scores] pdb_id=%s variant=%s ph=%s stage=%s valid=%d invalid=%d",
                paths.pdb_id,
                variant_label,
                ph_label,
                stage["name"],
                len(scores),
                len(invalids),
            )

        if not promoted_this_stage:
            restart, center, box_size, redo_ligands, recenter_attempts = early_recenter_decision(
                i, scores, distances, box_size, center, stage1_original, recenter_attempts, recenter_params,
                cfg, paths.pdb_id, receptor_pdbqt, logger, raw_docked, guard, control_anchor_hit
            )
            if restart:
                ligands = redo_ligands
                if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                    checkpoint_invalidate_from(
                        cfg,
                        paths.pdb_id,
                        stages_for_run,
                        start_index=0,
                        ph_label=ph_label,
                        variant=variant_env or None,
                    )
                gnina_jobs.clear()
                i = 0
                continue

        try:
            if cfg.get("ADAPTIVE_SHRINK_ENABLE", True) and validated:
                med = (
                    float(np.median([d for d in distances if isinstance(d, (int, float))]))
                    if distances
                    else None
                )
                if (med is not None) and (med < float(cfg.get("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0))):
                    dec = float(cfg.get("ADAPTIVE_SHRINK_DEC", 4.0))
                    min_box = float(cfg.get("ADAPTIVE_SHRINK_MIN_BOX", 14.0))
                    new_box = tuple(max(min_box, s - dec) for s in box_size)
                    if new_box != box_size:
                        logger.info(
                            f"Adaptive shrink: median dist {med:.2f} A -> box {box_size} -> {new_box}"
                        )
                        box_size = new_box
        except Exception as _e:
            logger.warning(f"Adaptive shrink skipped: {_e}")

        if i < len(stages_for_run) - 1:
            if not scores:
                restart, center, box_size, redo_ligands = fallback_recentering_if_empty(
                    cfg, paths.pdb_id, stage['name'], scores, raw_docked,
                    receptor_pdbqt, center, box_size, stage1_original, logger, guard, control_anchor_hit
                )
                if restart:
                    ligands = redo_ligands
                    if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                        checkpoint_invalidate_from(
                            cfg,
                            paths.pdb_id,
                            stages_for_run,
                            start_index=0,
                            ph_label=ph_label,
                            variant=variant_env or None,
                        )
                    gnina_jobs.clear()
                    i = 0
                    continue
                else:
                    break

            next_stage_is_stage3 = _is_stage3(stages_for_run[i + 1]["name"])
            scores_for_selection = {
                lig: sc for lig, sc in scores.items() if norm(lig) not in control_norms_for_run
            }
            invalids_for_selection = {
                lig: val for lig, val in invalids.items() if norm(lig) not in control_norms_for_run
            }
            ligands = _apply_force_carry_and_doping(
                cfg,
                docking_mode,
                i,
                stages_for_run,
                scores_for_selection,
                logger,
                stage1_original=stage1_original,
                forced_extracted_for_stage3=forced_extracted_for_stage3,
                invalids=invalids_for_selection,
                next_stage_is_stage3=next_stage_is_stage3,
                control_norms=control_norms_for_run,
            )
            if not ligands:
                logger.warning(
                    f"No ligands selected for {stages_for_run[i + 1]['name']}; stopping."
                )
                break

        # Engine follow-ups
        if validated:
            policy = decide_engine_policy(
                cfg=cfg,
                pdb_id=paths.pdb_id,
                variant=variant_env or None,
                ph_label=ph_label,
                center=center,
                ledock_enabled=ledock_enabled,
                dock6_enabled=dock6_enabled,
                logger=logger,
            )

            run_id_token = str(cfg.get("RUN_ID") or "")
            if run_id_token:
                try:
                    update_manifest_for_druggability_and_engine_plan(
                        cfg=cfg,
                        run_id=run_id_token,
                        pdb_id=paths.pdb_id,
                        variant_label=variant_env or None,
                        ph_tag=ph_label,
                        tier=policy.tier,
                        use_gnina=policy.use_gnina,
                        use_ledock=policy.use_ledock,
                        use_dock6=policy.use_dock6,
                    )
                except Exception:
                    logger.warning(
                        "[run-manifest.druggability-plan.skip] run_id=%s pdb=%s variant=%s ph=%s",
                        run_id_token,
                        paths.pdb_id,
                        variant_env,
                        ph_label,
                        exc_info=True,
                    )

            if policy.use_gnina:
                gnina_jobs.append(
                    {
                        "stage_name": stage["name"],
                        "stage_info": dict(stage),
                        "ligands": list(validated),
                        "center": tuple(center) if center is not None else None,
                        "box_size": tuple(box_size) if box_size is not None else None,
                    }
                )
            else:
                logger.info(
                    "[gnina.skip] pdb=%s ph=%s tier=%s reason=%s",
                    paths.pdb_id,
                    ph_label if ph_label else "base",
                    policy.tier,
                    policy.reason,
                )

            if policy.use_ledock:
                ledock_jobs.append(
                    {
                        "stage_name": stage["name"],
                        "stage_info": dict(stage),
                        "ligands": list(validated),
                        "center": tuple(center) if center is not None else None,
                        "box_size": tuple(box_size) if box_size is not None else None,
                    }
                )
            if policy.use_dock6:
                dock6_jobs.append(
                    {
                        "stage_name": stage["name"],
                        "stage_info": dict(stage),
                        "ligands": list(validated),
                        "center": tuple(center) if center is not None else None,
                        "box_size": tuple(box_size) if box_size is not None else None,
                        "variant": variant_env or None,
                        "ph_label": ph_label,
                    }
                )

        threads_per_vina = int(cfg.get("THREADS_PER_VINA", 1))

        def _expected_vina_path(lig: str) -> Path:
            return stage_dir / f"{Path(lig).stem}_{stage['name']}.pdbqt"

        def _rerun_vina_missing(lig: str) -> tuple[bool, Optional[str], Optional[Path]]:
            stage_for_cfg = dict(stage)
            stage_for_cfg["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
            if cfg.get("FAST_MODE"):
                stage_for_cfg["exhaustiveness"] = 1
            conf_path, out_path = emit_vina_config(
                cfg,
                paths.pdb_id,
                receptor_pdbqt,
                center,
                box_size,
                lig,
                stage["name"],
                stage_for_cfg,
                threads_per_vina,
                logger,
                variant=variant_env or None,
                ph_token=ph_label,
                legacy=legacy_mode,
            )
            try:
                Path(conf_path).resolve().relative_to(Path(cfg["CONFIG_RUN_DIR"]).resolve())
            except Exception:
                return False, "rerun_config_outside_run_dir", None
            try:
                run_docking_task(
                    cfg["VINA_EXE"],
                    conf_path,
                    lig,
                    out_path,
                    write_failure_marker_flag=True,
                )
            except Exception as exc:
                return False, f"rerun_error:{exc}", None
            if Path(out_path).exists() and Path(out_path).stat().st_size > 0:
                return True, "rerun_ok", None
            return False, "rerun_no_output", None

        completion_report_vina = run_completion_audit(
            engine="vina",
            pdb_id=paths.pdb_id,
            stage_name=stage["name"],
            ligands=stage_ligands_for_audit,
            expected_output_path=_expected_vina_path,
            rerun_one=_rerun_vina_missing,
            stage_dir=stage_dir,
            cfg=cfg,
            logger=logger,
            retries=1,
            ph_label=ph_label,
            variant=variant_env or None,
        )
        missing_after_vina = completion_report_vina.get("missing_ligands_after") or []
        if missing_after_vina:
            for lig_miss in missing_after_vina:
                if lig_miss not in score_history.get(stage["name"], {}):
                    record_score(
                        score_history,
                        stage["name"],
                        lig_miss,
                        None,
                        False,
                        reason="completion_missing",
                    )
                    record_le(score_history, stage["name"], lig_miss, None, heavy_atom_counts)

        if bool(cfg.get("CHECKPOINT_ENABLE", True)):
            try:
                fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                if completion_report_vina.get("success", False):
                    checkpoint_mark_done(
                        cfg,
                        paths.pdb_id,
                        stage["name"],
                        fp,
                        ph_label=ph_label,
                        variant=variant_env or None,
                    )
                else:
                    logger.warning(
                        "[Checkpoint] defer mark stage=%s pdb=%s reason=completion_missing",
                        stage["name"],
                        paths.pdb_id,
                    )
            except Exception:
                pass

        i += 1

    return MultistageVinaResult(
        scores=scores,
        validated=validated,
        distances=distances,
        raw_docked=raw_docked,
        invalids=invalids,
        validated_ligands_last=validated_ligands_last,
        forced_extracted_for_stage3=forced_extracted_for_stage3,
        recenter_attempts=recenter_attempts,
        gnina_jobs=gnina_jobs,
        ledock_jobs=ledock_jobs,
        dock6_jobs=dock6_jobs,
        dock6_metrics_by_stage=dock6_metrics_by_stage,
        score_history=score_history,
        docking_mode=docking_mode,
        switch_history=getattr(selector, "switch_history", []),
        global_switches=guard.global_switches,
    )
