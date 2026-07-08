from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, cast

from docking.docking_consensus_score import compute_consensus_for_variant_ph
from docking.docking_gnina_support import annotate_gnina_fda_long_csv_with_z_scores_vs_decoys
from docking.docking_ledock import (
    annotate_ledock_fda_long_csv_with_z_scores_vs_decoys,
    write_ledock_scores_csv,
)
from docking.docking_subrun_selection import _use_ledock
from docking.docking_utils import (
    _write_audit_json,
    final_pose_validation,
    norm,
)
from docking.docking_vina import write_scores_csv
from docking.decoy_score_annotation import annotate_fda_long_csv_with_z_scores_vs_decoys
from docking.library_mode import parse_test_libraries
from docking.docking_dock6 import write_dock6_scores_csv
from path_router.path_router import Paths
from prep_docking.prep_for_ledock import ensure_mol2_for_ledock


def run_post_engine_finalize_support(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    variant_env: str,
    variant_label: str,
    variant_token: Optional[str],
    ph_label: str,
    run_mode: Optional[str],
    csv_prefix: str,
    stages_for_run: List[Dict[str, Any]],
    score_history: Dict[str, Dict[str, Dict]],
    validated_ligands_last: List[str],
    ledock_metrics_by_stage: Dict[str, Dict[Path, Dict[str, Any]]],
    dock6_metrics_by_stage: Dict[str, Any],
    receptor_pdbqt: str,
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    docking_mode: str,
    stage1_original: List[str],
    controls_for_run: List[str],
    center_selector_history: List[Any],
    global_switches: int,
) -> None:
    _write_engine_score_csvs(
        cfg=cfg,
        paths=paths,
        ledock_metrics_by_stage=ledock_metrics_by_stage,
        dock6_metrics_by_stage=dock6_metrics_by_stage,
        ph_label=ph_label,
        variant_env=variant_env,
        csv_prefix=csv_prefix,
    )
    _validate_and_write_vina_outputs(
        cfg=cfg,
        paths=paths,
        logger=logger,
        variant_env=variant_env,
        variant_token=variant_token,
        ph_label=ph_label,
        csv_prefix=csv_prefix,
        stages_for_run=stages_for_run,
        score_history=score_history,
        validated_ligands_last=validated_ligands_last,
        receptor_pdbqt=receptor_pdbqt,
        center=center,
        docking_mode=docking_mode,
    )
    has_dud = _annotate_z_scores_if_needed(
        cfg=cfg,
        paths=paths,
        logger=logger,
        ph_label=ph_label,
        run_mode=run_mode,
        csv_prefix=csv_prefix,
    )
    _compute_consensus(
        cfg=cfg,
        paths=paths,
        logger=logger,
        ph_label=ph_label,
        variant_env=variant_env,
        variant_label=variant_label,
        csv_prefix=csv_prefix,
        run_mode=run_mode,
    )
    _prepare_ledock_mol2_if_enabled(
        cfg=cfg,
        logger=logger,
        pdb_id=paths.pdb_id,
        variant_label=variant_label,
        ph_label=ph_label,
        stage1_original=stage1_original,
        controls_for_run=controls_for_run,
    )
    _write_difficulty_audit_summary(
        cfg=cfg,
        paths=paths,
        logger=logger,
        variant_env=variant_env,
        ph_label=ph_label,
        run_mode=run_mode,
        has_dud=has_dud,
        center=center,
        box_size=box_size,
        stage1_original=stage1_original,
        validated_ligands_last=validated_ligands_last,
        center_selector_history=center_selector_history,
        global_switches=global_switches,
        stages_for_run=stages_for_run,
    )


def _write_engine_score_csvs(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    ledock_metrics_by_stage: Dict[str, Dict[Path, Dict[str, Any]]],
    dock6_metrics_by_stage: Dict[str, Any],
    ph_label: str,
    variant_env: str,
    csv_prefix: str,
) -> None:
    if ledock_metrics_by_stage:
        write_ledock_scores_csv(
            cfg,
            paths.pdb_id,
            cast(dict[str, dict[str, dict[str, Any]]], ledock_metrics_by_stage),
            ph_label=ph_label,
            variant=variant_env or None,
            csv_prefix=csv_prefix,
        )

    if dock6_metrics_by_stage:
        write_dock6_scores_csv(
            cfg,
            paths.pdb_id,
            dock6_metrics_by_stage,
            ph_label=ph_label,
            variant=variant_env or None,
            csv_prefix=csv_prefix,
        )


def _validate_and_write_vina_outputs(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    variant_env: str,
    variant_token: Optional[str],
    ph_label: str,
    csv_prefix: str,
    stages_for_run: List[Dict[str, Any]],
    score_history: Dict[str, Dict[str, Dict]],
    validated_ligands_last: List[str],
    receptor_pdbqt: str,
    center: Tuple[float, float, float],
    docking_mode: str,
) -> None:
    final_pose_validation(
        cfg,
        paths.pdb_id,
        stages_for_run,
        receptor_pdbqt,
        center,
        validated_ligands_last,
        score_history,
        docking_mode,
        logger,
        ph_label,
    )

    csv_path = write_scores_csv(
        cfg,
        paths.pdb_id,
        score_history,
        ph_label=ph_label,
        variant=variant_env or None,
        csv_prefix=csv_prefix,
    )
    logger.info(
        "[Scores] ph_label=%s summary=%s",
        ph_label if ph_label else "base",
        csv_path,
    )


def _annotate_z_scores_if_needed(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    ph_label: str,
    run_mode: Optional[str],
    csv_prefix: str,
) -> bool:
    try:
        tokens_now = parse_test_libraries(cfg)
    except Exception:
        tokens_now = []
    has_dud = "dud" in tokens_now

    if run_mode == "fda" and has_dud:
        decoy_csv_prefix = "dud_"
        try:
            from docking.docking_subruns import subruns_for_tokens

            for subrun in subruns_for_tokens(tokens_now):
                if subrun.run_mode == "dud":
                    decoy_csv_prefix = subrun.csv_prefix
                    break
        except Exception as e:
            logger.warning(
                "[z-score.warn] pdb_id=%s ph=%s reason=decoy_prefix_lookup_failed error=%s",
                paths.pdb_id,
                ph_label if ph_label else "base",
                e,
            )
        try:
            annotate_fda_long_csv_with_z_scores_vs_decoys(
                cfg,
                paths.pdb_id,
                ph_label=ph_label,
                csv_prefix=csv_prefix,
                decoy_csv_prefix=decoy_csv_prefix,
                logger=logger,
            )
        except Exception as e:
            logger.warning(
                "[z-score.warn] pdb_id=%s ph=%s reason=%s",
                paths.pdb_id,
                ph_label if ph_label else "base",
                e,
            )
        try:
            annotate_gnina_fda_long_csv_with_z_scores_vs_decoys(
                cfg,
                paths.pdb_id,
                ph_label=ph_label,
                csv_prefix=csv_prefix,
                decoy_csv_prefix=decoy_csv_prefix,
                logger=logger,
            )
        except Exception as e:
            logger.warning(
                "[gnina.z-score.warn] pdb_id=%s ph=%s reason=%s",
                paths.pdb_id,
                ph_label if ph_label else "base",
                e,
            )
        try:
            annotate_ledock_fda_long_csv_with_z_scores_vs_decoys(
                cfg,
                paths.pdb_id,
                ph_label=ph_label,
                csv_prefix=csv_prefix,
                decoy_csv_prefix=decoy_csv_prefix,
                logger=logger,
            )
        except Exception as e:
            logger.warning(
                "[ledock.z-score.warn] pdb_id=%s ph=%s reason=%s",
                paths.pdb_id,
                ph_label if ph_label else "base",
                e,
            )

    return has_dud


def _compute_consensus(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    ph_label: str,
    variant_env: str,
    variant_label: str,
    csv_prefix: str,
    run_mode: Optional[str],
) -> None:
    try:
        compute_consensus_for_variant_ph(
            cfg=cfg,
            paths=paths,
            ph_label=ph_label,
            variant_env=variant_env,
            variant_label=variant_label,
            csv_prefix=csv_prefix,
            run_mode=run_mode,
            logger=logger,
        )
    except Exception:
        logger.exception(
            "[consensus.error] Failed to compute consensus scores for pdb_id=%s variant=%s ph=%s",
            paths.pdb_id,
            variant_label,
            ph_label,
        )


def _prepare_ledock_mol2_if_enabled(
    *,
    cfg: Dict[str, Any],
    logger: logging.Logger,
    pdb_id: str,
    variant_label: str,
    ph_label: str,
    stage1_original: List[str],
    controls_for_run: List[str],
) -> None:
    if _use_ledock(cfg):
        try:
            logger.info(
                "[ledock.mol2] starting_mol2_prep pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_label,
                ph_label if ph_label else "base",
            )
            ligands_for_mol2 = _norm_dedupe(stage1_original + list(controls_for_run))
            ensure_mol2_for_ledock(cfg, ligands_for_mol2, logger)
            logger.info(
                "[ledock.mol2] completed_mol2_prep pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_label,
                ph_label if ph_label else "base",
            )
        except Exception as e:
            logger.warning(
                "[ledock.mol2.warn] pdb=%s variant=%s ph=%s reason=%s",
                pdb_id,
                variant_label,
                ph_label if ph_label else "base",
                e,
                exc_info=True,
            )
    else:
        logger.info(
            "[ledock.mol2] skip_mol2_prep pdb=%s variant=%s ph=%s reason=use_ledock_disabled",
            pdb_id,
            variant_label,
            ph_label if ph_label else "base",
        )


def _norm_dedupe(seq: List[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for item in seq:
        normalized = norm(item)
        if normalized not in seen:
            seen.add(normalized)
            out.append(normalized)
    return out


def _write_difficulty_audit_summary(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    variant_env: str,
    ph_label: str,
    run_mode: Optional[str],
    has_dud: bool,
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    stage1_original: List[str],
    validated_ligands_last: List[str],
    center_selector_history: List[Any],
    global_switches: int,
    stages_for_run: List[Dict[str, Any]],
) -> None:
    difficulty_info = None
    should_eval_difficulty = has_dud
    if should_eval_difficulty and run_mode not in (None, "dud"):
        should_eval_difficulty = False

    if should_eval_difficulty:
        logger.info(
            "[difficulty] pdb=%s action=skip reason=druggability_orchestrator_enabled",
            paths.pdb_id,
        )

    try:
        summary = {
            "pdb_id": paths.pdb_id,
            "center": tuple(map(float, center)) if center else None,
            "box_size": tuple(map(float, box_size)) if box_size else None,
            "n_ligands_stage1": len(stage1_original),
            "n_valid_last_stage": len(validated_ligands_last),
            "switch_history": center_selector_history,
            "global_switches": global_switches,
            "stages": [stage["name"] for stage in stages_for_run],
            "ph_label": ph_label,
        }
        if difficulty_info:
            summary["difficulty"] = difficulty_info.difficulty
            summary["roc_auc"] = difficulty_info.roc_auc
            summary["difficulty_N"] = difficulty_info.N
            summary["difficulty_n_actives"] = difficulty_info.n_actives
        _write_audit_json(
            cfg, paths.pdb_id, summary, ph_label=ph_label, variant=variant_env or None
        )
    except Exception as exc:
        logger.warning("Audit JSON write failed: %s", exc)
