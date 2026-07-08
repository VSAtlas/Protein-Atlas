from __future__ import annotations

import logging
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from docking.docking_dock6 import (
    run_dock6_for_stage,
)
from docking.docking_gnina import run_gnina_for_stage
from docking.docking_gnina_support import write_gnina_scores_csv
from docking.docking_ledock import (
    run_ledock_for_stage,
)
from docking.docking_ligand_selection import compute_stage_membership_from_scores
from docking.docking_subrun_selection import _use_dock6
from docking.docking_subrun_finalize_support import run_post_engine_finalize_support
from docking.docking_utils import (
    _fingerprint_stage,
)
from docking.ligand_metrics import compute_ligand_efficiency, record_le
from docking.score_io import record_score
from path_router.path_router import Paths
from prep_docking.prep_dock6 import ensure_dock6_site
from cli.run_manifest_runtime import (
    update_manifest_for_docking_stage,
    update_manifest_for_protein_success,
)
from docking.checkpoints import checkpoint_mark_done, checkpoint_should_skip


def finalize_ph_subrun(
    *,
    cfg: Dict[str, Any],
    paths: Paths,
    logger: logging.Logger,
    pdb_id: str,
    variant_env: str,
    variant_token: Optional[str],
    variant_label: str,
    legacy_mode: bool,
    ph_label: str,
    ph_start_ts: float,
    receptor_pdbqt: str,
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    run_mode: Optional[str],
    csv_prefix: str,
    stage_name_prefix: str,
    stages_for_run: List[Dict[str, Any]],
    score_history: Dict[str, Dict[str, Dict]],
    validated: List[str],
    validated_ligands_last: List[str],
    gnina_jobs: List[Dict[str, Any]],
    ledock_jobs: List[Dict[str, Any]],
    dock6_jobs: List[Dict[str, Any]],
    dock6_metrics_by_stage: Dict[str, Any],
    heavy_atom_counts: Dict[str, int],
    pains_flags: Dict[str, Any],
    library_for_manifest: Optional[str],
    docking_mode: str,
    stage1_original: List[str],
    controls_for_run: List[str],
    control_lookup: Dict[str, Path],
    recenter_attempts: Any,  # unused but passed for consistency if needed
    center_selector_history: List[Any],
    global_switches: int,
) -> None:
    """
    Orchestrates follow-up docking (GNINA, LeDock, DOCK6), reporting, consensus,
    audit JSON writing, and manifest success/failure updates.
    """
    manifest_run_id = cfg.get("RUN_ID")
    score_history_gnina: defaultdict[str, dict[str, Any]] = defaultdict(dict)
    gnina_records: list[dict[str, Any]] = []
    gnina_metrics_by_stage_best: defaultdict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    gnina_metrics_by_stage_full: defaultdict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    ledock_metrics_by_stage: defaultdict[str, dict[Path, dict[str, Any]]] = defaultdict(dict)
    dock6_enabled = _use_dock6(cfg)

    # --- GNINA Execution ---
    if gnina_jobs:
        logger.info(
            "[gnina.scheduler] pdb=%s variant=%s ph=%s policy=after_all_vina_stages n_jobs=%d",
            paths.pdb_id,
            variant_label,
            ph_label if ph_label else "base",
            len(gnina_jobs),
        )
        for job in gnina_jobs:
            gnina_stage_name = f"gnina_{job['stage_name']}"
            gnina_fp = None
            if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                try:
                    gnina_fp = _fingerprint_stage(
                        cfg, receptor_pdbqt, center, box_size, job["stage_info"]
                    )
                    if checkpoint_should_skip(
                        cfg,
                        paths.pdb_id,
                        gnina_stage_name,
                        gnina_fp,
                        ph_label=ph_label,
                        variant=variant_env or None,
                        engine="gnina",
                    ):
                        logger.info(
                            f"[Checkpoint] Skipping {gnina_stage_name} (fingerprint matched)."
                        )
                        if manifest_run_id:
                            try:
                                update_manifest_for_docking_stage(
                                    cfg,
                                    manifest_run_id,
                                    paths.pdb_id,
                                    variant_label,
                                    gnina_stage_name,
                                    status="completed",
                                    ph_tag=ph_label,
                                    elapsed_sec=0.0,
                                )
                            except Exception:
                                logger.warning(
                                    "[run-manifest.docking-stage] failed to record GNINA checkpoint skip pdb=%s variant=%s ph=%s stage=%s",
                                    paths.pdb_id,
                                    variant_label,
                                    ph_label if ph_label is not None else "base",
                                    gnina_stage_name,
                                    exc_info=True,
                                )
                        continue
                except Exception:
                    gnina_fp = None
            gnina_start_ts = time.time()
            if manifest_run_id:
                try:
                    update_manifest_for_docking_stage(
                        cfg,
                        manifest_run_id,
                        paths.pdb_id,
                        variant_label,
                        gnina_stage_name,
                        status="running",
                        ph_tag=ph_label,
                    )
                except Exception:
                    logger.warning(
                        "[run-manifest.docking-stage] failed to record GNINA start pdb=%s variant=%s ph=%s stage=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label if ph_label is not None else "base",
                        gnina_stage_name,
                        exc_info=True,
                    )

            try:
                scores_gnina, gnina_metrics, gnina_completion = run_gnina_for_stage(
                    cfg=cfg,
                    paths=paths,
                    pdb_id=paths.pdb_id,
                    variant=variant_env or None,
                    ph_label=ph_label,
                    stage_name=job["stage_name"],
                    stage_info=job["stage_info"],
                    ligands=job["ligands"],
                    center=job["center"],
                    box_size=job["box_size"],
                    logger=logger,
                    receptor_pdbqt=receptor_pdbqt,
                    control_lookup=control_lookup,
                )
                gnina_elapsed = time.time() - gnina_start_ts
                if manifest_run_id:
                    try:
                        update_manifest_for_docking_stage(
                            cfg,
                            manifest_run_id,
                            paths.pdb_id,
                            variant_label,
                            gnina_stage_name,
                            status="completed",
                            ph_tag=ph_label,
                            elapsed_sec=gnina_elapsed,
                        )
                    except Exception:
                        logger.warning(
                            "[run-manifest.docking-stage] failed to record GNINA end pdb=%s variant=%s ph=%s stage=%s",
                            paths.pdb_id,
                            variant_label,
                            ph_label if ph_label is not None else "base",
                            gnina_stage_name,
                            exc_info=True,
                        )
                if bool(cfg.get("CHECKPOINT_ENABLE", True)) and gnina_fp is not None:
                    try:
                        if gnina_completion.get("success", False):
                            checkpoint_mark_done(
                                cfg,
                                paths.pdb_id,
                                gnina_stage_name,
                                gnina_fp,
                                ph_label=ph_label,
                                variant=variant_env or None,
                            )
                        else:
                            logger.warning(
                                "[Checkpoint] defer mark gnina_stage=%s pdb=%s reason=completion_missing",
                                gnina_stage_name,
                                paths.pdb_id,
                            )
                    except Exception:
                        pass
            except Exception as e:
                if manifest_run_id:
                    try:
                        update_manifest_for_docking_stage(
                            cfg,
                            manifest_run_id,
                            paths.pdb_id,
                            variant_label,
                            gnina_stage_name,
                            status="failed",
                            ph_tag=ph_label,
                            error=str(e),
                            elapsed_sec=time.time() - gnina_start_ts,
                        )
                    except Exception:
                        logger.warning(
                            "[run-manifest.docking-stage] failed to record GNINA failure pdb=%s variant=%s ph=%s stage=%s",
                            paths.pdb_id,
                            variant_label,
                            ph_label if ph_label is not None else "base",
                            gnina_stage_name,
                            exc_info=True,
                        )
                raise
            gnina_stage_key = f"gnina_{job['stage_name']}"
            for lig, sc in scores_gnina.items():
                metrics_entry = gnina_metrics.get(
                    lig,
                    {
                        "minimized_affinity_kcal": None,
                        "cnn_score": None,
                        "cnn_affinity_pK": None,
                        "gnina_primary_score": sc,
                        "valid": False,
                        "reason": "gnina_failed" if sc is None else "",
                        "self_rmsd": None,
                    },
                )
                primary_score = metrics_entry.get("gnina_primary_score", sc)
                minimized_affinity = metrics_entry.get("minimized_affinity_kcal")
                cnn_score = metrics_entry.get("cnn_score")
                cnn_affinity = metrics_entry.get("cnn_affinity_pK")
                valid_flag = bool(metrics_entry.get("valid", False))
                reason_str = metrics_entry.get("reason", "") or (
                    "gnina_failed" if not valid_flag else ""
                )
                self_rmsd_val = metrics_entry.get("self_rmsd", None)

                ha_val = heavy_atom_counts.get(lig)
                le_val = compute_ligand_efficiency(primary_score, ha_val)
                pains_val = pains_flags.get(lig, pains_flags.get(Path(lig).stem, False))

                gnina_metrics_by_stage_full[gnina_stage_key][lig] = {
                    "minimized_affinity_kcal": minimized_affinity,
                    "cnn_score": cnn_score,
                    "cnn_affinity_pK": cnn_affinity,
                    "gnina_primary_score": primary_score,
                    "valid": bool(valid_flag),
                    "reason": reason_str,
                    "heavy_atoms": int(ha_val)
                    if isinstance(ha_val, (int, float))
                    else None,
                    "le": le_val,
                    "self_rmsd": self_rmsd_val,
                    "pains_flag": bool(pains_val) if pains_val is not None else False,
                }

                gnina_records.append(
                    {
                        "stage_name": gnina_stage_key,
                        "ligand": lig,
                        "primary_score": primary_score,
                        "minimized_affinity_kcal": minimized_affinity,
                        "cnn_score": cnn_score,
                        "cnn_affinity_pK": cnn_affinity,
                        "valid": bool(valid_flag),
                        "reason": reason_str or "",
                        "heavy_atoms": int(ha_val)
                        if isinstance(ha_val, (int, float))
                        else None,
                        "le": le_val,
                        "self_rmsd": self_rmsd_val,
                        "pains_flag": bool(pains_val)
                        if pains_val is not None
                        else False,
                    }
                )
    else:
        logger.info(
            "[gnina.scheduler] pdb=%s variant=%s ph=%s action=skip reason=no_jobs",
            paths.pdb_id,
            variant_label,
            ph_label if ph_label else "base",
        )

    if gnina_records:
        gnina_primary_scores: Dict[str, float] = {}
        for rec in gnina_records:
            lig = rec["ligand"]
            sc = rec.get("primary_score")
            if not isinstance(sc, (int, float)) or not math.isfinite(sc):
                continue
            prev = gnina_primary_scores.get(lig)
            if prev is None or sc > prev:
                gnina_primary_scores[lig] = sc

        gnina_stage_membership: Dict[int, List[str]] = {}
        gnina_stage_index_for_ligand: Dict[str, int] = {}
        if gnina_primary_scores:
            try:
                gnina_stage_membership = compute_stage_membership_from_scores(
                    cfg=cfg,
                    docking_mode=docking_mode,
                    scores=gnina_primary_scores,
                    higher_is_better=True,
                    n_stages=len(stages_for_run),
                )
            except Exception as e:
                logger.warning(
                    "[gnina.staging] failed to compute stage membership: %s", e
                )
                gnina_stage_membership = {}

        for stage_idx, lig_list in (gnina_stage_membership or {}).items():
            for lig in lig_list:
                gnina_stage_index_for_ligand[lig] = stage_idx

        best_rec_for_lig: Dict[str, Dict[str, Any]] = {}
        for rec in gnina_records:
            lig = rec["ligand"]
            sc = rec.get("primary_score")
            prev_record = best_rec_for_lig.get(lig)
            if prev_record is None:
                best_rec_for_lig[lig] = rec
                continue
            prev_sc = prev_record.get("primary_score")
            if (
                isinstance(sc, (int, float))
                and math.isfinite(sc)
                and (
                    not isinstance(prev_sc, (int, float))
                    or not math.isfinite(prev_sc)
                    or sc > prev_sc
                )
            ):
                best_rec_for_lig[lig] = rec

        for lig, rec in best_rec_for_lig.items():
            ligand_stage_idx = gnina_stage_index_for_ligand.get(lig)
            if ligand_stage_idx is not None:
                gnina_stage_name = f"gnina_stage{ligand_stage_idx}"
            else:
                gnina_stage_name = rec.get("stage_name") or "gnina_stage1"

            record_score(
                score_history_gnina,
                gnina_stage_name,
                lig,
                rec["primary_score"],
                rec["valid"],
                reason=rec["reason"] or None,
            )
            record_le(
                score_history_gnina,
                gnina_stage_name,
                lig,
                rec["primary_score"],
                heavy_atom_counts,
            )

            gnina_metrics_by_stage_best[gnina_stage_name][lig] = {
                "minimized_affinity_kcal": rec["minimized_affinity_kcal"],
                "cnn_score": rec["cnn_score"],
                "cnn_affinity_pK": rec["cnn_affinity_pK"],
                "gnina_primary_score": rec["primary_score"],
                "valid": rec["valid"],
                "reason": rec["reason"],
                "heavy_atoms": rec["heavy_atoms"],
                "le": rec["le"],
                "self_rmsd": rec["self_rmsd"],
                "pains_flag": rec["pains_flag"],
            }

    if gnina_metrics_by_stage_full:
        write_gnina_scores_csv(
            cfg,
            paths.pdb_id,
            gnina_metrics_by_stage_full,
            ph_label=ph_label,
            variant=variant_env or None,
            csv_prefix=csv_prefix,
        )

    # --- DOCK6 Execution ---
    if dock6_jobs and dock6_enabled:
        try:
            logger.info(
                "[DOCK6_PREP_ONCE] ensuring DOCK6 site for pdb=%s variant=%s ph=%s",
                paths.pdb_id,
                variant_label,
                ph_label if ph_label else "base",
            )
            ensure_dock6_site(
                cfg=cfg,
                pdb_id=paths.pdb_id,
                variant=variant_env or None,
                ph_label=ph_label,
                logger=logger,
            )
        except Exception as exc:
            logger.warning(
                "[dock6.surface.skip] pdb=%s variant=%s ph=%s reason=%s",
                paths.pdb_id,
                variant_label,
                ph_label if ph_label else "base",
                exc,
            )

    if dock6_jobs:
        logger.info(
            "[dock6.scheduler] pdb=%s variant=%s ph=%s jobs=%d",
            paths.pdb_id,
            variant_label,
            ph_label if ph_label else "base",
            len(dock6_jobs),
        )
        for job in dock6_jobs:
            stage_name = job["stage_name"]
            stage_info = job["stage_info"]
            ligands = [Path(ligand_path) for ligand_path in job.get("ligands", [])]
            dock6_stage_name = f"dock6_{stage_name}"
            dock6_start_ts = time.time()
            if manifest_run_id:
                try:
                    update_manifest_for_docking_stage(
                        cfg,
                        manifest_run_id,
                        paths.pdb_id,
                        variant_label,
                        dock6_stage_name,
                        status="running",
                        elapsed_sec=None,
                        ph_tag=ph_label,
                    )
                except Exception:
                    logger.warning(
                        "[run-manifest.docking-stage] failed to record DOCK6 start pdb=%s variant=%s ph=%s stage=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label if ph_label is not None else "base",
                        dock6_stage_name,
                        exc_info=True,
                    )

            try:
                job_center = job.get("center")
                job_box_size = job.get("box_size")
                if not (
                    isinstance(job_center, tuple)
                    and len(job_center) == 3
                    and isinstance(job_box_size, tuple)
                    and len(job_box_size) == 3
                ):
                    raise ValueError(
                        f"invalid DOCK6 geometry for stage {stage_name}: {job}"
                    )
                dock6_scores, dock6_metrics = run_dock6_for_stage(
                    cfg=cfg,
                    paths=paths,
                    pdb_id=paths.pdb_id,
                    variant=variant_env or None,
                    ph_label=ph_label,
                    stage_name=stage_name,
                    stage_info=stage_info,
                    ligands=ligands,
                    center=job_center,
                    box_size=job_box_size,
                    logger=logger,
                )
                dock6_metrics_by_stage[stage_name] = dock6_metrics
                valid_count = sum(
                    1 for rec in dock6_metrics.values() if rec.get("valid")
                )
                invalid_count = max(len(dock6_metrics) - valid_count, 0)
                dock6_elapsed = time.time() - dock6_start_ts
                if manifest_run_id:
                    try:
                        update_manifest_for_docking_stage(
                            cfg,
                            manifest_run_id,
                            paths.pdb_id,
                            variant_label,
                            dock6_stage_name,
                            status="completed",
                            elapsed_sec=dock6_elapsed,
                            ph_tag=ph_label,
                        )
                    except Exception:
                        logger.warning(
                            "[run-manifest.docking-stage] failed to record DOCK6 completion pdb=%s variant=%s ph=%s stage=%s",
                            paths.pdb_id,
                            variant_label,
                            ph_label if ph_label is not None else "base",
                            dock6_stage_name,
                            exc_info=True,
                        )
                logger.info(
                    "[dock6.done] pdb=%s stage=%s variant=%s ph=%s valid=%d invalid=%d elapsed_sec=%.2f",
                    paths.pdb_id,
                    stage_name,
                    variant_label,
                    ph_label if ph_label else "base",
                    valid_count,
                    invalid_count,
                    dock6_elapsed,
                )
            except Exception as e:
                logger.warning(
                    "[dock6.error] pdb=%s stage=%s variant=%s ph=%s reason=%s",
                    paths.pdb_id,
                    stage_name,
                    variant_label,
                    ph_label if ph_label else "base",
                    e,
                    exc_info=True,
                )
                dock6_metrics_by_stage[stage_name] = {}
                if manifest_run_id:
                    try:
                        update_manifest_for_docking_stage(
                            cfg,
                            manifest_run_id,
                            paths.pdb_id,
                            variant_label,
                            dock6_stage_name,
                            status="failed",
                            ph_tag=ph_label,
                            error=str(e),
                            elapsed_sec=time.time() - dock6_start_ts,
                        )
                    except Exception:
                        logger.warning(
                            "[run-manifest.docking-stage] failed to record DOCK6 failure pdb=%s variant=%s ph=%s stage=%s",
                            paths.pdb_id,
                            variant_label,
                            ph_label if ph_label is not None else "base",
                            dock6_stage_name,
                            exc_info=True,
                        )

    elif dock6_enabled:
        logger.info(
            "[dock6.scheduler] pdb=%s variant=%s ph=%s action=skip reason=no_jobs",
            paths.pdb_id,
            variant_label,
            ph_label if ph_label else "base",
        )

    # --- LeDock Execution ---
    if ledock_jobs:
        logger.info(
            "[ledock.scheduler] pdb=%s variant=%s ph=%s policy=after_gnina n_jobs=%d",
            paths.pdb_id,
            variant_label,
            ph_label if ph_label else "base",
            len(ledock_jobs),
        )
        for job in ledock_jobs:
            stage_name = job["stage_name"]
            ledock_stage_name = f"ledock_{stage_name}"
            ledock_fp = None

            if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                try:
                    ledock_fp = _fingerprint_stage(
                        cfg, receptor_pdbqt, center, box_size, job["stage_info"]
                    )
                    if checkpoint_should_skip(
                        cfg,
                        paths.pdb_id,
                        ledock_stage_name,
                        ledock_fp,
                        ph_label=ph_label,
                        variant=variant_env or None,
                        engine="ledock",
                    ):
                        logger.info(
                            f"[Checkpoint] Skipping {ledock_stage_name} (fingerprint matched)."
                        )
                        if manifest_run_id:
                            try:
                                update_manifest_for_docking_stage(
                                    cfg,
                                    manifest_run_id,
                                    paths.pdb_id,
                                    variant_label,
                                    ledock_stage_name,
                                    status="completed",
                                    ph_tag=ph_label,
                                    elapsed_sec=0.0,
                                )
                            except Exception:
                                logger.warning(
                                    "[run-manifest.docking-stage] failed to record LeDock checkpoint skip pdb=%s variant=%s ph=%s stage=%s",
                                    paths.pdb_id,
                                    variant_label,
                                    ph_label if ph_label is not None else "base",
                                    ledock_stage_name,
                                    exc_info=True,
                                )
                        continue
                except Exception:
                    ledock_fp = None

            ledock_start_ts = time.time()
            if manifest_run_id:
                try:
                    update_manifest_for_docking_stage(
                        cfg,
                        manifest_run_id,
                        paths.pdb_id,
                        variant_label,
                        ledock_stage_name,
                        status="running",
                        elapsed_sec=None,
                        ph_tag=ph_label,
                    )
                except Exception:
                    logger.warning(
                        "[run-manifest.docking-stage] failed to record LeDock start pdb=%s variant=%s ph=%s stage=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label if ph_label is not None else "base",
                        ledock_stage_name,
                        exc_info=True,
                    )
            try:
                scores_ledock, ledock_metrics, ledock_completion = run_ledock_for_stage(
                    cfg=cfg,
                    paths=paths,
                    pdb_id=paths.pdb_id,
                    variant=variant_env or None,
                    ph_label=ph_label,
                    stage_name=stage_name,
                    stage_info=job["stage_info"],
                    ligands=job["ligands"],
                    center=job["center"],
                    box_size=job["box_size"],
                    logger=logger,
                )
                ledock_metrics_by_stage[stage_name] = ledock_metrics
                valid_count = sum(
                    1 for rec in ledock_metrics.values() if rec.get("valid")
                )
                invalid_count = max(len(ledock_metrics) - valid_count, 0)
                ledock_elapsed = time.time() - ledock_start_ts
                if manifest_run_id:
                    try:
                        update_manifest_for_docking_stage(
                            cfg,
                            manifest_run_id,
                            paths.pdb_id,
                            variant_label,
                            ledock_stage_name,
                            status="completed",
                            elapsed_sec=ledock_elapsed,
                            ph_tag=ph_label,
                        )
                    except Exception:
                        logger.warning(
                            "[run-manifest.docking-stage] failed to record LeDock end pdb=%s variant=%s ph=%s stage=%s",
                            paths.pdb_id,
                            variant_label,
                            ph_label if ph_label is not None else "base",
                            ledock_stage_name,
                            exc_info=True,
                        )
                if bool(cfg.get("CHECKPOINT_ENABLE", True)) and ledock_fp is not None:
                    try:
                        if ledock_completion.get("success", False):
                            checkpoint_mark_done(
                                cfg,
                                paths.pdb_id,
                                ledock_stage_name,
                                ledock_fp,
                                ph_label=ph_label,
                                variant=variant_env or None,
                            )
                        else:
                            logger.warning(
                                "[Checkpoint] defer mark ledock_stage=%s pdb=%s reason=completion_missing",
                                ledock_stage_name,
                                paths.pdb_id,
                            )
                    except Exception:
                        pass
                logger.info(
                    "[ledock.done] pdb=%s stage=%s variant=%s ph=%s valid=%d invalid=%d elapsed_sec=%.2f",
                    paths.pdb_id,
                    stage_name,
                    variant_label,
                    ph_label if ph_label else "base",
                    valid_count,
                    invalid_count,
                    time.time() - ledock_start_ts,
                )
            except Exception as e:
                logger.warning(
                    "[ledock.error] pdb=%s stage=%s variant=%s ph=%s reason=%s",
                    paths.pdb_id,
                    stage_name,
                    variant_label,
                    ph_label if ph_label else "base",
                    e,
                    exc_info=True,
                )
                ledock_metrics_by_stage[stage_name] = {}
                if manifest_run_id:
                    try:
                        update_manifest_for_docking_stage(
                            cfg,
                            manifest_run_id,
                            paths.pdb_id,
                            variant_label,
                            ledock_stage_name,
                            status="failed",
                            ph_tag=ph_label,
                            error=str(e),
                            elapsed_sec=time.time() - ledock_start_ts,
                        )
                    except Exception:
                        logger.warning(
                            "[run-manifest.docking-stage] failed to record LeDock failure pdb=%s variant=%s ph=%s stage=%s",
                            paths.pdb_id,
                            variant_label,
                            ph_label if ph_label else "base",
                            ledock_stage_name,
                            exc_info=True,
                        )
    else:
        logger.info(
            "[ledock.scheduler] pdb=%s variant=%s ph=%s action=skip reason=no_jobs",
            paths.pdb_id,
            variant_label,
            ph_label if ph_label else "base",
        )

    run_post_engine_finalize_support(
        cfg=cfg,
        paths=paths,
        logger=logger,
        variant_env=variant_env,
        variant_label=variant_label,
        variant_token=variant_token,
        ph_label=ph_label,
        run_mode=run_mode,
        csv_prefix=csv_prefix,
        stages_for_run=stages_for_run,
        score_history=score_history,
        validated_ligands_last=validated_ligands_last,
        ledock_metrics_by_stage=ledock_metrics_by_stage,
        dock6_metrics_by_stage=dock6_metrics_by_stage,
        receptor_pdbqt=receptor_pdbqt,
        center=center,
        box_size=box_size,
        docking_mode=docking_mode,
        stage1_original=stage1_original,
        controls_for_run=controls_for_run,
        center_selector_history=center_selector_history,
        global_switches=global_switches,
    )

    # --- Manifest Updates ---
    try:
        update_manifest_for_protein_success(
            cfg,
            manifest_run_id or "",
            paths.pdb_id,
            variant_label,
            time.time() - ph_start_ts,
            ph_tag=ph_label,
        )
    except Exception:
        ph_log = logging.getLogger("ph_ensemble")
        ph_log.warning(
            "[run-manifest.skip] pdb=%s variant=%s ph=%s reason=success",
            paths.pdb_id,
            variant_label,
            ph_label if ph_label else "base",
            exc_info=True,
        )
