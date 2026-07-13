from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Set

from docking.global_scheduler import acquire_global_cores
from post_docking.rescoring.scorch_coverage import (
    accept_quarantined_coverage,
    evaluate_scorch_coverage,
    remove_done_sentinel,
    resolve_top_fraction,
    write_coverage_summary,
)
from post_docking.rescoring.scorch_orchestration_build import build_orchestration_tasks
from post_docking.rescoring.scorch_orchestration_execute import execute_orchestration_tasks
from post_docking.rescoring.scorch_orchestration_finalize import finalize_orchestration_combos
from post_docking.rescoring.scorch_shards import execute_sharded_tasks, shard_mode_enabled
from post_docking.rescoring.scorch_orchestration_types import (
    ComboKey,
    DoneSentinelDeps,
    OrchestrationDeps,
    OrchestrationLimits,
    OrchestrationResult,
    receptor_path_for_combo as receptor_path_for_combo,
)
from post_docking.rescoring.scorch_types import StageSpec


def run_combo_task_orchestration(
    *,
    cfg: Dict[str, object],
    include_fda: bool,
    combos: Set[ComboKey],
    specs: list[StageSpec],
    run_root: Path,
    post_run_root: Path,
    processed_root: Path,
    source_post_run_root: Path | None = None,
    top_fraction: float,
    jobs: int,
    threads: int,
    overwrite: bool,
    decoy_prefix: str,
    component: str,
    logger: logging.Logger,
    limits: OrchestrationLimits,
    deps: OrchestrationDeps,
) -> OrchestrationResult:
    source_post_run_root = source_post_run_root or post_run_root
    build_state = build_orchestration_tasks(
        cfg=cfg,
        include_fda=include_fda,
        combos=combos,
        specs=specs,
        run_root=run_root,
        post_run_root=post_run_root,
        processed_root=processed_root,
        source_post_run_root=source_post_run_root,
        top_fraction=top_fraction,
        jobs=jobs,
        threads=threads,
        component=component,
        logger=logger,
        limits=limits,
        deps=deps,
    )
    total_jobs = len(build_state.tasks)
    combos_attempted = len(build_state.combos_with_tasks)
    if shard_mode_enabled(cfg):
        execute_result = execute_sharded_tasks(
            cfg=cfg,
            tasks=build_state.tasks,
            run_root=run_root,
            post_root=source_post_run_root,
            jobs=jobs,
            threads=threads,
            overwrite=overwrite,
            component=component,
            logger=logger,
            deps=deps,
            combo_failed_local=dict(build_state.combo_failed_local),
            decoy_prefix=decoy_prefix,
            top_fraction=top_fraction,
        )
        if execute_result.all_complete:
            local_completed_combos = finalize_orchestration_combos(
                cfg=cfg,
                combos=combos,
                combos_with_tasks=build_state.combos_with_tasks,
                combo_modes=build_state.combo_modes,
                combo_failed_local=execute_result.combo_failed_local,
                specs=specs,
                run_root=run_root,
                post_run_root=post_run_root,
                overwrite=overwrite,
                decoy_prefix=decoy_prefix,
                component=component,
                logger=logger,
                deps=deps,
                total_jobs=total_jobs,
                completed=execute_result.completed,
                failed_jobs=execute_result.failed_jobs,
                combos_attempted=combos_attempted,
                skipped_missing_receptor=build_state.skipped_missing_receptor,
                skipped_missing_consensus=build_state.skipped_missing_consensus,
            )
        else:
            local_completed_combos = set()
            logger.warning(
                "%s action=scorch_shards status=incomplete total=%d completed=%d failed=%d",
                component,
                int(execute_result.total),
                int(execute_result.completed),
                int(execute_result.failed_jobs),
            )
        return OrchestrationResult(
            failed_jobs=execute_result.failed_jobs,
            combos_with_tasks=set(build_state.combos_with_tasks),
            combo_failed_local=dict(execute_result.combo_failed_local),
            local_completed_combos=set(local_completed_combos),
        )
    execute_result = execute_orchestration_tasks(
        cfg=cfg,
        tasks=build_state.tasks,
        run_root=run_root,
        source_post_run_root=source_post_run_root,
        jobs=jobs,
        threads=threads,
        overwrite=overwrite,
        component=component,
        logger=logger,
        limits=limits,
        deps=deps,
        combo_failed_local=dict(build_state.combo_failed_local),
    )
    local_completed_combos = finalize_orchestration_combos(
        cfg=cfg,
        combos=combos,
        combos_with_tasks=build_state.combos_with_tasks,
        combo_modes=build_state.combo_modes,
        combo_failed_local=execute_result.combo_failed_local,
        specs=specs,
        run_root=run_root,
        post_run_root=post_run_root,
        overwrite=overwrite,
        decoy_prefix=decoy_prefix,
        component=component,
        logger=logger,
        deps=deps,
        total_jobs=total_jobs,
        completed=execute_result.completed,
        failed_jobs=execute_result.failed_jobs,
        combos_attempted=combos_attempted,
        skipped_missing_receptor=build_state.skipped_missing_receptor,
        skipped_missing_consensus=build_state.skipped_missing_consensus,
    )
    return OrchestrationResult(
        failed_jobs=execute_result.failed_jobs,
        combos_with_tasks=set(build_state.combos_with_tasks),
        combo_failed_local=dict(execute_result.combo_failed_local),
        local_completed_combos=set(local_completed_combos),
    )


def mark_done_sentinels(
    *,
    cfg: Dict[str, object],
    run_root: Path,
    post_run_root: Path,
    combos_with_tasks: Set[ComboKey],
    combo_failed: Dict[ComboKey, bool],
    decoy_prefix: str = "dud",
    component: str,
    logger: logging.Logger,
    deps: DoneSentinelDeps,
) -> Set[ComboKey]:
    incomplete_combos: Set[ComboKey] = set()
    for combo in sorted(combos_with_tasks):
        if combo_failed.get(combo, True):
            incomplete_combos.add(combo)
            continue
        done_task_id = f"scorch:housekeeping:done:{combo[0]}:{combo[1]}:{combo[2]}"
        deps.emit_task_event(
            "TASK_SUBMITTED",
            task_id=done_task_id,
            task_type="housekeeping",
            want_cores=1,
            min_cores=1,
            est_duration_sec=0.2,
            pdb_id=combo[0],
            variant=combo[1],
            ph=combo[2],
        )
        with acquire_global_cores(
            cfg,
            cores=1,
            min_cores=1,
            priority=0,
            est_duration_sec=0.2,
            task_id=done_task_id,
            task_type="housekeeping",
        ) as granted:
            deps.emit_task_event(
                "TASK_START",
                task_id=done_task_id,
                task_type="housekeeping",
                want_cores=1,
                min_cores=1,
                granted_cores=int(granted),
                est_duration_sec=0.2,
                pdb_id=combo[0],
                variant=combo[1],
                ph=combo[2],
            )
            summary = evaluate_scorch_coverage(
                docked_run_root=run_root,
                post_run_root=post_run_root,
                pdb_id=combo[0],
                variant=combo[1],
                ph=combo[2],
                top_fraction=resolve_top_fraction(cfg),
                decoy_prefix=decoy_prefix,
                accept_quarantined=accept_quarantined_coverage(cfg),
            )
            try:
                write_coverage_summary(post_run_root, summary)
            except Exception:
                logger.debug(
                    "%s action=coverage_summary status=failed combo=%s",
                    component,
                    combo,
                    exc_info=True,
                )
            if not summary.all_complete:
                try:
                    remove_done_sentinel(post_run_root, summary)
                except Exception:
                    logger.debug(
                        "%s action=clear_done status=failed combo=%s",
                        component,
                        combo,
                        exc_info=True,
                    )
                logger.warning(
                    "%s action=done_sentinel status=skip reason=scorch_coverage_incomplete combo=%s details=%s",
                    component,
                    combo,
                    ",".join(summary.reasons) or "unknown",
                )
                incomplete_combos.add(combo)
                continue
            deps.mark_combo_done(post_run_root, combo, logger)
            deps.emit_task_event(
                "TASK_END",
                task_id=done_task_id,
                task_type="housekeeping",
                want_cores=1,
                min_cores=1,
                granted_cores=int(granted),
                est_duration_sec=0.2,
                status="ok",
                pdb_id=combo[0],
                variant=combo[1],
                ph=combo[2],
            )
    return incomplete_combos
