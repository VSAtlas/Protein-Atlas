from __future__ import annotations

import logging
import os
import subprocess
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

from docking.global_scheduler import acquire_global_cores
from post_docking.rescoring import scorch_selection as _scorch_selection_mod
from post_docking.rescoring import scorch_provisional_cache as _scorch_cache_mod
from post_docking.rescoring.scorch_device import (
    cfg_for_scorch_ligand_count,
    scorch_child_env,
)
from post_docking.rescoring.scorch_lock_fallback import run_with_lock_fallback
from post_docking.rescoring.scorch_types import AnnotateResult, MaterializeResult, StageSpec

EmitTaskEventFn = Callable[..., None]
MaterializeInputsFn = Callable[..., MaterializeResult]
ScorchCommandFn = Callable[[Path, Path, int, Mapping[str, Any]], List[str]]
AnnotateCsvFn = Callable[
    [Path, Dict[str, str], Optional[Dict[str, Dict[str, str]]], logging.Logger],
    AnnotateResult,
]


def _csv_data_row_count(csv_path: Path) -> int:
    try:
        import csv

        with csv_path.open() as handle:
            reader = csv.DictReader(handle)
            return sum(1 for _ in reader)
    except Exception:
        return 0


def _float_cfg_or_env(
    cfg: Mapping[str, object],
    *names: str,
    default: float,
) -> float:
    for name in names:
        raw = os.environ.get(name)
        if raw is None:
            raw = cfg.get(name)
        if raw is None:
            continue
        try:
            value = float(str(raw).strip())
        except Exception:
            continue
        if value >= 0.0:
            return float(value)
    return float(default)


def _optional_float_cfg_or_env(
    cfg: Mapping[str, object],
    *names: str,
) -> Optional[float]:
    for name in names:
        raw = os.environ.get(name)
        if raw is None:
            raw = cfg.get(name)
        if raw is None:
            continue
        try:
            value = float(str(raw).strip())
        except Exception:
            continue
        if value >= 0.0:
            return float(value)
    return None


def _adaptive_scorch_timeout_sec(ligand_count: int) -> float:
    count = max(1, int(ligand_count or 0))
    # SCORCH normally scores 128-ligand CPU shards in a few minutes.  A flat
    # one-hour default lets a single pathological ligand pin a target too long.
    return float(min(1800.0, max(600.0, 60.0 + (5.0 * count))))


def _stage_result_store(cfg: Dict[str, object]) -> Dict[str, Dict[str, object]]:
    raw = cfg.get("_ATLAS_SCORCH_STAGE_RESULTS")
    if not isinstance(raw, dict):
        raw = {}
        cfg["_ATLAS_SCORCH_STAGE_RESULTS"] = raw
    return raw  # type: ignore[return-value]


def _stage_timeout_override(cfg: Mapping[str, object], key: str) -> Optional[float]:
    raw = cfg.get("_ATLAS_SCORCH_SHARD_TIMEOUT_OVERRIDES")
    if not isinstance(raw, Mapping):
        return None
    value = raw.get(key)
    if value is None:
        return None
    try:
        parsed = float(str(value).strip())
    except Exception:
        return None
    return parsed if parsed > 0.0 else None


def score_stage(
    cfg: Dict[str, object],
    spec: StageSpec,
    combo: Tuple[str, str, str],
    run_root: Path,
    post_root: Path,
    receptor: Path,
    threads: int,
    overwrite: bool,
    logger: logging.Logger,
    *,
    decoy_prefix: str,
    component: str,
    emit_task_event: EmitTaskEventFn,
    materialize_inputs: MaterializeInputsFn,
    scorch_command: ScorchCommandFn,
    annotate_csv: AnnotateCsvFn,
    scorch_root: Optional[Path],
    scorch_python: Optional[Path],
    scorch_script: Optional[Path],
    allowed_bases: Optional[Set[str]] = None,
    control_bases: Optional[Set[str]] = None,
    run_mode: str = "fda",
    stage_dirs_override: Optional[Sequence[str]] = None,
    score_csv: Optional[Path] = None,
    selected_stage_by_base: Optional[Dict[str, str]] = None,
    selected_score_by_base: Optional[Dict[str, float]] = None,
    chunk_tag: Optional[str] = None,
    hedge_group_id: Optional[str] = None,
    hedge_role: str = "primary",
    rechunk_generation: int = 0,
    artifact_post_root: Optional[Path] = None,
    cache_read: bool = False,
    cache_write: bool = False,
    cache_phase: str = "provisional",
    run_id: Optional[str] = None,
) -> Tuple[bool, Optional[Path]]:
    pdb_id, variant, ph = combo
    source_combo_post_root = post_root / pdb_id / variant / ph
    combo_post_root = (artifact_post_root or post_root) / pdb_id / variant / ph
    output_name = spec.output_name if run_mode == "fda" else f"{decoy_prefix}_{spec.output_name}"
    if chunk_tag:
        output_name = output_name.replace(".csv", f".{chunk_tag}.csv")
    out_path = combo_post_root / output_name
    combo_post_root.mkdir(parents=True, exist_ok=True)
    rescoring_ligand_count = 0
    run_id_for_cache = str(run_id or run_root.name)

    def _cleanup_stale_part_outputs(*, reuse_only: bool = False) -> None:
        if chunk_tag:
            return
        try:
            pattern = output_name.replace(
                ".csv",
                ".partreuse*.csv" if reuse_only else ".part*.csv",
            )
            for stale in combo_post_root.glob(pattern):
                stale.unlink(missing_ok=True)
        except Exception:
            logger.debug(
                "%s action=stale_part_cleanup status=skip output=%s",
                component,
                out_path,
                exc_info=True,
            )

    _cleanup_stale_part_outputs()

    if out_path.exists() and not overwrite and out_path.stat().st_size > 0:
        existing_rows = _csv_data_row_count(out_path)
        if existing_rows <= 0:
            logger.warning(
                "%s action=score status=rerun source=%s stage=%s reason=existing_output_has_no_rows output=%s",
                component,
                spec.source,
                spec.stage_dir,
                out_path,
            )
        else:
            logger.info(
                "%s action=score status=skip source=%s stage=%s reason=exists rows=%d output=%s",
                component,
                spec.source,
                spec.stage_dir,
                existing_rows,
                out_path,
            )
            return True, out_path
    elif out_path.exists() and not overwrite:
        logger.warning(
            "%s action=score status=rerun source=%s stage=%s reason=existing_output_empty output=%s",
            component,
            spec.source,
            spec.stage_dir,
            out_path,
        )
    def _apply_cached_reuse(ligands: Sequence[Path]) -> Tuple[List[Path], Optional[Path]]:
        nonlocal allowed_bases
        if not cache_read or overwrite or allowed_bases is None:
            return list(ligands), None
        try:
            fields, rows, covered = _scorch_cache_mod.reusable_rows_for_task(
                cfg=cfg,
                run_id=run_id_for_cache,
                combo=combo,
                spec=spec,
                run_mode=run_mode,
                decoy_prefix=decoy_prefix,
                allowed_bases=set(allowed_bases),
                selected_stage_by_base=selected_stage_by_base,
                selected_score_by_base=selected_score_by_base,
                ligands=list(ligands),
                receptor=receptor,
                score_csv=score_csv,
                logger=logger,
            )
        except Exception:
            logger.warning(
                "%s action=cache_reuse status=skip source=%s stage=%s reason=lookup_failed",
                component,
                spec.source,
                spec.stage_dir,
                exc_info=True,
            )
            return list(ligands), None
        requested = {str(x) for x in allowed_bases if str(x).strip()}
        covered_requested = set(covered) & requested
        if not rows or not covered_requested:
            return list(ligands), None
        if covered_requested >= requested:
            try:
                _cleanup_stale_part_outputs()
                _scorch_cache_mod.write_rows_csv(out_path, fields, rows)
            except Exception:
                logger.warning(
                    "%s action=cache_reuse status=skip source=%s stage=%s reason=write_failed output=%s",
                    component,
                    spec.source,
                    spec.stage_dir,
                    out_path,
                    exc_info=True,
                )
                return list(ligands), None
            logger.info(
                "%s action=cache_reuse status=full source=%s stage=%s run_mode=%s ligands=%d output=%s",
                component,
                spec.source,
                spec.stage_dir,
                run_mode,
                len(covered_requested),
                out_path,
            )
            return [], out_path

        logger.info(
            "%s action=cache_reuse status=partial_skip source=%s stage=%s run_mode=%s reusable=%d requested=%d reason=missing_some_final_rows",
            component,
            spec.source,
            spec.stage_dir,
            run_mode,
            len(covered_requested),
            len(requested),
        )
        return list(ligands), None

    lig_path: Optional[Path]
    rescored_stage_by_base: Dict[str, str] = {}
    stage_fallback_reason_by_base: Dict[str, str] = {}
    if spec.source in {"vina", "gnina"}:
        ph_root = run_root / pdb_id / variant / ph
        stage_counts: Optional[Dict[int, int]] = None
        stage_dirs: Sequence[str] = stage_dirs_override or (spec.stage_dir,)
        use_best = spec.stage_dir in {"vina_best", "gnina_best"}
        if use_best:
            if stage_dirs_override is None and spec.stage_dir == "vina_best":
                stage_dirs = ("stage3", "stage2", "stage1")
            elif stage_dirs_override is None and spec.stage_dir == "gnina_best":
                stage_dirs = ("gnina_stage3", "gnina_stage2", "gnina_stage1")
            pose_result = _scorch_selection_mod.collect_best_pose_per_base(
                ph_root,
                stage_dirs,
                allowed_bases,
                logger,
                component=component,
                decoy_prefix=decoy_prefix,
                preferred_stage_by_base=selected_stage_by_base,
            )
            ligands_available = pose_result.ligands
            stage_counts = pose_result.stage_counts
            total_candidates = pose_result.total_candidates
            available_bases = pose_result.available_bases
            rescored_stage_by_base = pose_result.rescored_stage_by_base
            stage_fallback_reason_by_base = pose_result.stage_fallback_reason_by_base
        else:
            pose_result = _scorch_selection_mod.collect_best_pose_per_base(
                ph_root,
                stage_dirs,
                allowed_bases,
                logger,
                component=component,
                decoy_prefix=decoy_prefix,
                preferred_stage_by_base=selected_stage_by_base,
            )
            ligands_available = pose_result.ligands
            stage_counts = pose_result.stage_counts
            total_candidates = pose_result.total_candidates
            available_bases = pose_result.available_bases
            rescored_stage_by_base = pose_result.rescored_stage_by_base
            stage_fallback_reason_by_base = pose_result.stage_fallback_reason_by_base
        rescoring_ligand_count = len(ligands_available)
        if allowed_bases is not None and not ligands_available and total_candidates > 0:
            logger.debug(
                "%s action=select status=debug source=%s stage=%s reason=filtered_empty candidate_bases=%s allowed_sample=%s",
                component,
                spec.source,
                spec.stage_dir,
                sorted(list(available_bases))[:10],
                sorted(list(allowed_bases))[:10],
            )
        if not ligands_available:
            logger.warning(
                "%s action=score status=skip source=%s stage=%s reason=no_pdbqt stage_dirs=%s path=%s",
                component,
                spec.source,
                spec.stage_dir,
                ",".join(stage_dirs),
                ph_root,
            )
            return True, None
        if stage_counts is not None:
            logger.debug(
                "%s action=score source=%s stage=%s candidates=%d ligands_unique=%d stage3=%d stage2=%d stage1=%d stage0=%d",
                component,
                spec.source,
                spec.stage_dir,
                total_candidates,
                len(ligands_available),
                stage_counts.get(3, 0),
                stage_counts.get(2, 0),
                stage_counts.get(1, 0),
                stage_counts.get(0, 0),
            )
        else:
            logger.debug(
                "%s action=score source=%s stage=%s ligands_before=%d ligands_after=%d",
                component,
                spec.source,
                spec.stage_dir,
                total_candidates,
                len(ligands_available),
            )
        ligands_available, cached_full_path = _apply_cached_reuse(ligands_available)
        if cached_full_path is not None:
            return True, cached_full_path
        if not ligands_available:
            logger.info(
                "%s action=score status=skip source=%s stage=%s reason=cache_reuse_exhausted",
                component,
                spec.source,
                spec.stage_dir,
            )
            return True, out_path
        allowed_size = len(allowed_bases) if allowed_bases is not None else 0
        missing_bases = (
            sorted((allowed_bases or set()) - available_bases)[:10]
            if allowed_bases is not None
            else []
        )
        logger.info(
            "%s action=select source=%s stage=%s allowed=%d candidates_bases=%d rescored=%d missing=%d missing_examples=%s",
            component,
            spec.source,
            spec.stage_dir,
            allowed_size,
            len(available_bases),
            len(ligands_available),
            len(missing_bases) if allowed_bases is not None else 0,
            missing_bases,
        )
        if allowed_bases is not None and control_bases is not None:
            noncontrol_allowed = allowed_bases - control_bases
            noncontrol_available = available_bases - control_bases
            if not noncontrol_allowed and noncontrol_available:
                logger.debug(
                    "%s action=select status=debug source=%s stage=%s reason=controls_only_selection noncontrol_available=%d examples=%s",
                    component,
                    spec.source,
                    spec.stage_dir,
                    len(noncontrol_available),
                    sorted(list(noncontrol_available))[:10],
                )
        mat_task_id = (
            f"scorch:materialize:{pdb_id}:{variant}:{ph}:{spec.source}:{run_mode}:{chunk_tag or 'all'}"
        )
        mat_est_sec = max(0.1, rescoring_ligand_count * 0.01)
        emit_task_event(
            "TASK_SUBMITTED",
            task_id=mat_task_id,
            task_type="scorch_materialize",
            want_cores=1,
            min_cores=1,
            est_duration_sec=mat_est_sec,
            pdb_id=pdb_id,
            variant=variant,
            ph=ph,
        )
        lig_path = None
        mat_result = MaterializeResult(
            input_dir=None,
            materialized_count=0,
            failed_count=0,
            failed_examples=(),
        )
        mat_error: Optional[str] = None
        mat_nonfatal_missing = False
        mat_status = "error"
        with acquire_global_cores(
            cfg,
            cores=1,
            min_cores=1,
            priority=1,
            task_id=mat_task_id,
            task_type="scorch_materialize",
            est_duration_sec=mat_est_sec,
        ) as granted:
            emit_task_event(
                "TASK_START",
                task_id=mat_task_id,
                task_type="scorch_materialize",
                want_cores=1,
                min_cores=1,
                granted_cores=int(granted),
                est_duration_sec=mat_est_sec,
                pdb_id=pdb_id,
                variant=variant,
                ph=ph,
            )
            try:
                mat_result = materialize_inputs(
                    ph_root,
                    combo_post_root,
                    spec.stage_dir,
                    ligands_available,
                    overwrite,
                    logger,
                    run_mode=run_mode,
                    chunk_tag=chunk_tag,
                )
                lig_path = mat_result.input_dir
                if lig_path is None:
                    mat_status = "error"
                elif mat_result.failed_count > 0:
                    mat_status = "degraded"
                else:
                    mat_status = "ok"
            except Exception as exc:
                if isinstance(exc, FileNotFoundError):
                    mat_nonfatal_missing = True
                    mat_status = "degraded"
                    logger.warning(
                        "%s action=score status=degraded source=%s stage=%s reason=materialize_exception_nonfatal error=%s",
                        component,
                        spec.source,
                        spec.stage_dir,
                        exc,
                    )
                else:
                    mat_error = str(exc)
                    logger.error(
                        "%s action=score status=failed source=%s stage=%s reason=materialize_exception error=%s",
                        component,
                        spec.source,
                        spec.stage_dir,
                        exc,
                    )
            mat_end_extra: Dict[str, Any] = {}
            if mat_error:
                mat_end_extra["error"] = mat_error
            mat_end_extra["materialized"] = int(mat_result.materialized_count)
            mat_end_extra["failed"] = int(mat_result.failed_count)
            emit_task_event(
                "TASK_END",
                task_id=mat_task_id,
                task_type="scorch_materialize",
                want_cores=1,
                min_cores=1,
                granted_cores=int(granted),
                est_duration_sec=mat_est_sec,
                status=mat_status,
                pdb_id=pdb_id,
                variant=variant,
                ph=ph,
                **mat_end_extra,
            )
        if mat_error:
            return False, None
        if lig_path is None:
            if mat_nonfatal_missing:
                logger.warning(
                    "%s action=score status=skip source=%s stage=%s reason=materialize_exception_nonfatal_no_inputs",
                    component,
                    spec.source,
                    spec.stage_dir,
                )
                return True, None
            logger.error(
                "%s action=score status=failed source=%s stage=%s reason=materialize_failed",
                component,
                spec.source,
                spec.stage_dir,
            )
            return False, None
        rescoring_ligand_count = int(mat_result.materialized_count)
        if mat_result.failed_count > 0:
            logger.warning(
                "%s action=materialize status=degraded source=%s stage=%s materialized=%d failed=%d failed_examples=%s",
                component,
                spec.source,
                spec.stage_dir,
                mat_result.materialized_count,
                mat_result.failed_count,
                list(mat_result.failed_examples),
            )
    else:
        if score_csv is None:
            score_prefix = f"{decoy_prefix}_" if run_mode == "dud" else ""
            score_csv = (
                run_root
                / pdb_id
                / variant
                / ph
                / f"{score_prefix}{spec.source}_docking_score_long.csv"
            )
        if allowed_bases is not None and score_csv is not None:
            _scorch_selection_mod.log_score_csv_coverage(
                spec.source,
                combo,
                score_csv,
                allowed_bases,
                logger,
                component=component,
                decoy_prefix=decoy_prefix,
            )
        ph_root = source_combo_post_root
        stage_dirs = stage_dirs_override or (spec.stage_dir,)
        pose_result = _scorch_selection_mod.collect_best_pose_per_base(
            ph_root,
            stage_dirs,
            allowed_bases,
            logger,
            component=component,
            decoy_prefix=decoy_prefix,
            preferred_stage_by_base=selected_stage_by_base,
        )
        ligands_available = pose_result.ligands
        stage_counts = pose_result.stage_counts
        total_candidates = pose_result.total_candidates
        available_bases = pose_result.available_bases
        rescored_stage_by_base = pose_result.rescored_stage_by_base
        stage_fallback_reason_by_base = pose_result.stage_fallback_reason_by_base
        if allowed_bases is not None and not ligands_available and total_candidates > 0:
            logger.debug(
                "%s action=select status=debug source=%s stage=%s reason=filtered_empty candidate_bases=%s allowed_sample=%s",
                component,
                spec.source,
                spec.stage_dir,
                sorted(list(available_bases))[:10],
                sorted(list(allowed_bases))[:10],
            )
        if not ligands_available:
            logger.warning(
                "%s action=score status=skip source=%s stage=%s reason=no_pdbqt stage_dirs=%s path=%s",
                component,
                spec.source,
                spec.stage_dir,
                ",".join(stage_dirs),
                ph_root,
            )
            return True, None
        rescoring_ligand_count = len(ligands_available)
        logger.debug(
            "%s action=score source=%s stage=%s candidates=%d ligands_unique=%d stage3=%d stage2=%d stage1=%d stage0=%d",
            component,
            spec.source,
            spec.stage_dir,
            total_candidates,
            len(ligands_available),
            stage_counts.get(3, 0),
            stage_counts.get(2, 0),
            stage_counts.get(1, 0),
            stage_counts.get(0, 0),
        )
        ligands_available, cached_full_path = _apply_cached_reuse(ligands_available)
        if cached_full_path is not None:
            return True, cached_full_path
        if not ligands_available:
            logger.info(
                "%s action=score status=skip source=%s stage=%s reason=cache_reuse_exhausted",
                component,
                spec.source,
                spec.stage_dir,
            )
            return True, out_path
        allowed_size = len(allowed_bases) if allowed_bases is not None else 0
        missing_bases = (
            sorted((allowed_bases or set()) - available_bases)[:10]
            if allowed_bases is not None
            else []
        )
        logger.info(
            "%s action=select source=%s stage=%s allowed=%d candidates_bases=%d rescored=%d missing=%d missing_examples=%s",
            component,
            spec.source,
            spec.stage_dir,
            allowed_size,
            len(available_bases),
            len(ligands_available),
            len(missing_bases) if allowed_bases is not None else 0,
            missing_bases,
        )
        if allowed_bases is not None and control_bases is not None:
            noncontrol_allowed = allowed_bases - control_bases
            noncontrol_available = available_bases - control_bases
            if not noncontrol_allowed and noncontrol_available:
                logger.debug(
                    "%s action=select status=debug source=%s stage=%s reason=controls_only_selection noncontrol_available=%d examples=%s",
                    component,
                    spec.source,
                    spec.stage_dir,
                    len(noncontrol_available),
                    sorted(list(noncontrol_available))[:10],
                )
        try:
            mat_result = materialize_inputs(
                ph_root,
                combo_post_root,
                spec.stage_dir,
                ligands_available,
                overwrite,
                logger,
                run_mode=run_mode,
                chunk_tag=chunk_tag,
            )
        except FileNotFoundError as exc:
            logger.warning(
                "%s action=score status=skip source=%s stage=%s reason=materialize_exception_nonfatal error=%s",
                component,
                spec.source,
                spec.stage_dir,
                exc,
            )
            return True, None
        lig_path = mat_result.input_dir
        if lig_path is None:
            logger.error(
                "%s action=score status=failed source=%s stage=%s reason=materialize_failed",
                component,
                spec.source,
                spec.stage_dir,
            )
            return False, None
        rescoring_ligand_count = int(mat_result.materialized_count)
        if mat_result.failed_count > 0:
            logger.warning(
                "%s action=materialize status=degraded source=%s stage=%s materialized=%d failed=%d failed_examples=%s",
                component,
                spec.source,
                spec.stage_dir,
                mat_result.materialized_count,
                mat_result.failed_count,
                list(mat_result.failed_examples),
            )

    score_task_id = f"scorch:score:{pdb_id}:{variant}:{ph}:{spec.source}:{run_mode}:{chunk_tag or 'all'}"
    score_want_cores = max(1, int(threads))
    if rescoring_ligand_count <= max(2, score_want_cores):
        score_want_cores = 1
    score_est_sec = max(0.5, (rescoring_ligand_count * 0.08) / score_want_cores)
    task_extra: Dict[str, Any] = {}
    if hedge_group_id:
        task_extra["hedge_group_id"] = str(hedge_group_id)
    if hedge_role:
        task_extra["hedge_role"] = str(hedge_role)
    if int(rechunk_generation) > 0:
        task_extra["rechunk_generation"] = int(rechunk_generation)
    emit_task_event(
        "TASK_SUBMITTED",
        task_id=score_task_id,
        task_type="scorch_score_chunk",
        want_cores=score_want_cores,
        min_cores=1,
        est_duration_sec=score_est_sec,
        pdb_id=pdb_id,
        variant=variant,
        ph=ph,
        **task_extra,
    )
    proc: Optional[subprocess.CompletedProcess[str]] = None
    score_error: Optional[str] = None
    score_status = "error"
    score_started = time.time()
    configured_timeout_sec = _optional_float_cfg_or_env(
        cfg,
        "ATLAS_SCORCH_SUBPROCESS_TIMEOUT_SEC",
        "SCORCH_SUBPROCESS_TIMEOUT_SEC",
        "ATLAS_SCORCH_SHARD_TIMEOUT_SEC",
        "SCORCH_SHARD_TIMEOUT_SEC",
    )
    stage_result_key = str(hedge_group_id or "").strip()
    if stage_result_key:
        timeout_override = _stage_timeout_override(cfg, stage_result_key)
        if timeout_override is not None:
            configured_timeout_sec = float(timeout_override)
    timeout_sec = (
        float(configured_timeout_sec)
        if configured_timeout_sec is not None
        else _adaptive_scorch_timeout_sec(rescoring_ligand_count)
    )
    heartbeat_sec = _float_cfg_or_env(
        cfg,
        "ATLAS_SCORCH_SUBPROCESS_HEARTBEAT_SEC",
        "SCORCH_SUBPROCESS_HEARTBEAT_SEC",
        default=60.0,
    )

    def _score_heartbeat() -> None:
        logger.info(
            "%s action=score heartbeat=true source=%s stage=%s task_id=%s elapsed_s=%.1f timeout_s=%.1f output=%s",
            component,
            spec.source,
            spec.stage_dir,
            score_task_id,
            max(0.0, time.time() - score_started),
            float(timeout_sec),
            out_path,
        )

    def _record_stage_result(
        *,
        status: str,
        returncode: Optional[int],
        error: str = "",
        stderr: str = "",
    ) -> None:
        if not stage_result_key:
            return
        stderr_tail = str(stderr or "")[-2000:]
        timed_out = bool(
            returncode == 124
            or "SCORCH subprocess timeout after" in stderr_tail
            or "TimeoutExpired" in str(error or "")
        )
        failure_reason = ""
        if str(status) != "ok":
            if timed_out:
                failure_reason = "timeout"
            elif error:
                failure_reason = "exception"
            elif returncode is not None:
                failure_reason = "nonzero_returncode"
            else:
                failure_reason = "unknown"
        _stage_result_store(cfg)[stage_result_key] = {
            "status": str(status),
            "returncode": None if returncode is None else int(returncode),
            "failure_reason": failure_reason,
            "timed_out": bool(timed_out),
            "timeout_sec": float(timeout_sec),
            "elapsed_sec": max(0.0, time.time() - score_started),
            "stderr_tail": stderr_tail,
            "error": str(error or ""),
            "rescoring_ligand_count": int(rescoring_ligand_count),
        }

    with acquire_global_cores(
        cfg,
        cores=score_want_cores,
        min_cores=1,
        priority=2,
        est_duration_sec=score_est_sec,
        task_id=score_task_id,
        task_type="scorch_score_chunk",
    ) as granted:
        score_cfg = cfg_for_scorch_ligand_count(cfg, rescoring_ligand_count, logger)
        cmd = scorch_command(receptor, lig_path, int(granted), score_cfg)
        cmd[cmd.index("{out}")] = str(out_path)
        run_kwargs: Dict[str, Any] = {
            "capture_output": True,
            "text": True,
            "env": scorch_child_env(score_cfg, task_id=score_task_id),
        }
        if timeout_sec > 0:
            run_kwargs["timeout"] = float(timeout_sec)
            run_kwargs["heartbeat_fn"] = _score_heartbeat
            run_kwargs["heartbeat_sec"] = float(max(1.0, heartbeat_sec))
        if scorch_root is not None:
            run_kwargs["cwd"] = str(scorch_root)
        emit_task_event(
            "TASK_START",
            task_id=score_task_id,
            task_type="scorch_score_chunk",
            want_cores=score_want_cores,
            min_cores=1,
            granted_cores=int(granted),
            est_duration_sec=score_est_sec,
            pdb_id=pdb_id,
            variant=variant,
            ph=ph,
            **task_extra,
        )
        try:
            proc, _lock_retry_used = run_with_lock_fallback(
                cmd=cmd,
                run_kwargs=run_kwargs,
                logger=logger,
                component=component,
                source=str(spec.source),
                stage=str(spec.stage_dir),
                scorch_python=scorch_python,
                scorch_script=scorch_script,
            )
            score_status = "ok" if proc.returncode == 0 else "error"
        except Exception as exc:
            score_error = str(exc)
            logger.error(
                "%s action=score status=failed source=%s stage=%s reason=worker_exception error=%s",
                component,
                spec.source,
                spec.stage_dir,
                exc,
            )
        score_end_extra: Dict[str, Any] = {}
        if score_error:
            score_end_extra["error"] = score_error
        emit_task_event(
            "TASK_END",
            task_id=score_task_id,
            task_type="scorch_score_chunk",
            want_cores=score_want_cores,
            min_cores=1,
            granted_cores=int(granted),
            est_duration_sec=score_est_sec,
            status=score_status,
            pdb_id=pdb_id,
            variant=variant,
            ph=ph,
            **score_end_extra,
            **task_extra,
        )
    if score_error is not None or proc is None:
        _record_stage_result(
            status="error",
            returncode=None,
            error=str(score_error or "missing_process"),
        )
        return False, None
    if proc.returncode != 0:
        _record_stage_result(
            status="error",
            returncode=int(proc.returncode),
            stderr=str(proc.stderr or ""),
        )
        logger.error(
            "%s action=score status=failed source=%s stage=%s returncode=%s stderr=%s",
            component,
            spec.source,
            spec.stage_dir,
            proc.returncode,
            proc.stderr.strip(),
        )
        return False, None

    if not out_path.exists() or out_path.stat().st_size == 0:
        _record_stage_result(
            status="error",
            returncode=int(proc.returncode),
            error="empty_output",
            stderr=str(proc.stderr or ""),
        )
        logger.error(
            "%s action=score status=failed source=%s stage=%s reason=empty_output output=%s",
            component,
            spec.source,
            spec.stage_dir,
            out_path,
        )
        return False, None

    metadata = {
        "pdb_id": pdb_id,
        "variant": variant,
        "ph": ph,
        "source": spec.source,
        "stage_dir": spec.stage_dir,
        "run_mode": run_mode,
    }
    row_metadata_by_base: Dict[str, Dict[str, str]] = {}
    if rescored_stage_by_base or selected_stage_by_base or selected_score_by_base:
        for base, rescored_stage in rescored_stage_by_base.items():
            selected_stage = (
                str((selected_stage_by_base or {}).get(base, "")).strip()
                if selected_stage_by_base
                else ""
            )
            selected_score_raw = (
                (selected_score_by_base or {}).get(base)
                if selected_score_by_base
                else None
            )
            fallback_reason = stage_fallback_reason_by_base.get(base, "")
            row_metadata_by_base[base] = {
                "selected_stage": selected_stage,
                "rescored_stage": str(rescored_stage or "").strip(),
                "selected_docking_score": (
                    ""
                    if selected_score_raw is None
                    else f"{float(selected_score_raw):.6g}"
                ),
                "stage_match_flag": (
                    "1"
                    if selected_stage and str(rescored_stage or "").strip() == selected_stage
                    else ("0" if selected_stage else "")
                ),
                "stage_fallback_reason": fallback_reason,
            }
    annotate_result = annotate_csv(out_path, metadata, row_metadata_by_base or None, logger)
    if not annotate_result.ok:
        row_count = _csv_data_row_count(out_path)
        if annotate_result.degraded and row_count > 0:
            logger.warning(
                "%s action=annotate status=degraded source=%s stage=%s reason=%s rows=%d output=%s",
                component,
                spec.source,
                spec.stage_dir,
                annotate_result.reason,
                row_count,
                out_path,
            )
        else:
            _record_stage_result(
                status="error",
                returncode=int(proc.returncode),
                error=f"annotate_failed:{annotate_result.reason}",
                stderr=str(proc.stderr or ""),
            )
            logger.error(
                "%s action=annotate status=failed source=%s stage=%s reason=%s rows=%d output=%s",
                component,
                spec.source,
                spec.stage_dir,
                annotate_result.reason,
                row_count,
                out_path,
            )
            return False, None

    logger.info(
        "%s action=score status=ok source=%s stage=%s ligands=%d output=%s",
        component,
        spec.source,
        spec.stage_dir,
        rescoring_ligand_count,
        out_path,
    )
    _record_stage_result(
        status="ok",
        returncode=int(proc.returncode),
        stderr=str(proc.stderr or ""),
    )
    if cache_write:
        try:
            _scorch_cache_mod.record_task_output(
                cfg=cfg,
                run_id=run_id_for_cache,
                combo=combo,
                spec=spec,
                run_mode=run_mode,
                decoy_prefix=decoy_prefix,
                receptor=receptor,
                score_csv=score_csv,
                allowed_bases=allowed_bases,
                selected_stage_by_base=selected_stage_by_base,
                selected_score_by_base=selected_score_by_base,
                ligands=list(ligands_available),
                output_csv=out_path,
                phase=cache_phase,
                logger=logger,
            )
        except Exception:
            logger.warning(
                "%s action=cache_record status=skip source=%s stage=%s reason=unexpected_error output=%s",
                component,
                spec.source,
                spec.stage_dir,
                out_path,
                exc_info=True,
            )
    return True, out_path
