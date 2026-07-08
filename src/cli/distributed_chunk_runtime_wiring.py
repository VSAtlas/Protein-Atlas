from __future__ import annotations

import logging
import os
import time
from typing import Any, Callable, Mapping, Optional

from cli.distributed_chunk_planner import (
    _CHUNK_MAX_ATTEMPTS,
    _distributed_chunk_cpu_per_worker,
    _distributed_chunk_local_workers,
)
from cli.distributed_context import read_chunk_result
from cli.distributed_runtime_claim_types import ClaimLoopContext
from cli.distributed_runtime_completion import DeferredScorchCompletionBookkeeper
from cli.distributed_chunk_runtime_helpers import (
    aggregate_completion_markers_for_owned_scopes,
    aggregate_completion_markers_for_scope,
    collect_combo_chunk_progress,
    collect_terminal_chunk_failures,
)
from cli.run_context import ConfigDict
from docking.global_scheduler import acquire_global_slot
from cli.postrun_hooks_runtime import _maybe_run_scorch_rescore_for_pdb

ScopeKey = tuple[str, str, str]


def build_deferred_scorch_completion_bookkeeper(
    *,
    cfg_v: ConfigDict,
    run_id: str,
    label: str,
    dist_ctx: Any,
    combo_total_chunks: Mapping[ScopeKey, int],
    combo_total_ligands: Mapping[ScopeKey, int],
    combo_owner_task_id: Mapping[ScopeKey, int],
    scorch_queue_service: Any,
    variant_param_for_scorch: Optional[str],
    by_chunk_id: Mapping[str, Mapping[str, Any]],
) -> DeferredScorchCompletionBookkeeper:
    return DeferredScorchCompletionBookkeeper(
        cfg_v=cfg_v,
        run_id=str(run_id),
        label=str(label),
        dist_task_id=int(dist_ctx.task_id),
        combo_total_chunks=combo_total_chunks,
        combo_total_ligands=combo_total_ligands,
        combo_owner_task_id=combo_owner_task_id,
        scorch_queue_service=scorch_queue_service,
        variant_param_for_scorch=variant_param_for_scorch,
        hybrid_start_ts=float(time.time()),
        acquire_global_slot=acquire_global_slot,
        maybe_run_scorch_rescore_for_pdb=_maybe_run_scorch_rescore_for_pdb,
        normalize_ph_tag_token=lambda raw: str(raw).strip().lower(),
        combo_progress_reader=lambda: collect_combo_chunk_progress(
            by_chunk_id=by_chunk_id,
            label=label,
            dist_ctx=dist_ctx,
            max_attempts=int(_CHUNK_MAX_ATTEMPTS),
            read_chunk_result=read_chunk_result,
            cfg_v=cfg_v,
        ),
        before_queue_combo_scorch=lambda combo_scope: aggregate_completion_markers_for_scope(
            combo_scope=combo_scope,
            cfg_v=cfg_v,
            run_id=str(run_id),
            logger=logging.getLogger("completion.aggregate"),
        ),
    )


def build_claim_loop_context(
    *,
    cfg_v: ConfigDict,
    dist_ctx: Any,
    run_id: str,
    label: str,
    variant: Optional[str],
    by_chunk_id: dict[str, dict[str, Any]],
    combo_owner_task_id: Mapping[ScopeKey, int],
    combo_total_chunks: Mapping[ScopeKey, int],
    combo_hybrid_active: bool,
    execution_mode_decision: Any,
    combo_total_weight: float,
    baseline_first_complete_sec: Optional[float],
    baseline_makespan_sec: Optional[float],
    process_one: Callable[..., bool],
    distributed_assigned_pdb_ids: set[str],
    recorded_terminal_chunk_ids: set[str],
    failed_entries: list[tuple[str, str, str, str, str]],
    bar: Any,
) -> ClaimLoopContext:
    local_chunk_workers = _distributed_chunk_local_workers(cfg_v)
    cpu_per_chunk_worker = _distributed_chunk_cpu_per_worker(
        cfg_v, local_workers=local_chunk_workers
    )
    logging.info(
        "[distributed.chunk.local-workers] variant=%s task_id=%d total_cpu=%d local_workers=%d cpu_per_worker=%d planned_vina_slots=%d idle_cpu_budget=%d",
        label.upper(),
        int(dist_ctx.task_id),
        int(cfg_v.get("CPU", os.cpu_count() or 1) or 1),
        int(local_chunk_workers),
        int(cpu_per_chunk_worker),
        int(local_chunk_workers) * int(cpu_per_chunk_worker),
        max(
            0,
            int(cfg_v.get("CPU", os.cpu_count() or 1) or 1)
            - (int(local_chunk_workers) * int(cpu_per_chunk_worker)),
        ),
    )
    return ClaimLoopContext(
        cfg_v=cfg_v,
        dist_ctx=dist_ctx,
        run_id=str(run_id),
        label=str(label),
        variant=None if variant is None else str(variant),
        by_chunk_id=by_chunk_id,
        combo_owner_task_id=combo_owner_task_id,
        combo_total_chunks=combo_total_chunks,
        combo_hybrid_active=bool(combo_hybrid_active),
        execution_mode_decision=execution_mode_decision,
        local_chunk_workers=int(local_chunk_workers),
        cpu_per_chunk_worker=int(cpu_per_chunk_worker),
        combo_total_weight=float(combo_total_weight),
        baseline_first_complete_sec=baseline_first_complete_sec,
        baseline_makespan_sec=baseline_makespan_sec,
        process_one=process_one,
        distributed_assigned_pdb_ids=distributed_assigned_pdb_ids,
        recorded_terminal_chunk_ids=recorded_terminal_chunk_ids,
        failed_entries=failed_entries,
        bar=bar,
    )


def finalize_variant_claim_loop(
    *,
    by_chunk_id: Mapping[str, Mapping[str, Any]],
    label: str,
    dist_ctx: Any,
    cfg_v: ConfigDict,
    run_id: str,
    combo_owner_task_id: Mapping[ScopeKey, int],
    combo_total_chunks: Mapping[ScopeKey, int],
    recorded_terminal_chunk_ids: set[str],
    failed_entries: list[tuple[str, str, str, str, str]],
    completion_bookkeeper: DeferredScorchCompletionBookkeeper,
    claim_loop_rebalance_count: int,
    run_chunk_rebalance_count: int,
) -> int:
    combo_completed_chunks = collect_terminal_chunk_failures(
        by_chunk_id=by_chunk_id,
        label=label,
        dist_ctx=dist_ctx,
        max_attempts=int(_CHUNK_MAX_ATTEMPTS),
        read_chunk_result=read_chunk_result,
        cfg_v=cfg_v,
        recorded_terminal_chunk_ids=recorded_terminal_chunk_ids,
        failed_entries=failed_entries,
    )
    combo_progress = collect_combo_chunk_progress(
        by_chunk_id=by_chunk_id,
        label=label,
        dist_ctx=dist_ctx,
        max_attempts=int(_CHUNK_MAX_ATTEMPTS),
        read_chunk_result=read_chunk_result,
        cfg_v=cfg_v,
    )
    combo_failed_chunks = {
        scope: int(counts[2])
        for scope, counts in combo_progress.items()
        if len(counts) >= 3 and int(counts[2]) > 0
    }
    aggregate_completion_markers_for_owned_scopes(
        combo_owner_task_id=combo_owner_task_id,
        combo_total_chunks=combo_total_chunks,
        dist_ctx=dist_ctx,
        cfg_v=cfg_v,
        run_id=str(run_id),
    )
    run_chunk_rebalance_count += int(max(0, claim_loop_rebalance_count))
    completion_bookkeeper.flush_owned_deferred_submissions(
        combo_completed_chunks=combo_completed_chunks,
        combo_failed_chunks=combo_failed_chunks,
    )
    return int(run_chunk_rebalance_count)
