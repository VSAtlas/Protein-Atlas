from __future__ import annotations

import logging
import time
from typing import Any, Callable, Optional

from cli.distributed_combo_runner import resolve_combo_runtime_mode
from cli.distributed_chunk_planner import (
    _ensure_distributed_manifest_preflight,
    _load_or_create_variant_chunk_plan,
    _distributed_combo_owner_task_id,
    _distributed_chunk_local_workers,
)
from cli.distributed_chunk_runtime_helpers import (
    build_combo_weight_totals,
    build_scope_maps,
)
from cli.distributed_chunk_runtime_wiring import (
    build_claim_loop_context,
    build_deferred_scorch_completion_bookkeeper,
    finalize_variant_claim_loop,
)
from cli.distributed_runtime_claim_loop import run_chunk_claim_loop
from cli.distributed_runtime_completion import DeferredScorchCompletionBookkeeper
from cli.run_context import ConfigDict
from cli.run_process_one import (
    _resume_docking_outputs_complete,
    _resume_manifest_entry,
    _resume_scorch_required,
    _resume_scorch_still_pending,
)


ScopeKey = tuple[str, str, str]

# Compatibility aliases for tests/monkeypatch hooks.
_run_chunk_claim_loop = run_chunk_claim_loop
_DeferredScorchCompletionBookkeeper = DeferredScorchCompletionBookkeeper


def _resume_manifest_scope_entry(
    cfg_v: ConfigDict,
    *,
    scope: ScopeKey,
    fallback_label: str,
) -> Any:
    pdb_scope, variant_scope, ph_scope = scope
    variant_label = str(variant_scope or fallback_label or "legacy").strip() or "legacy"
    ph_token = str(ph_scope or "base").strip() or "base"
    return _resume_manifest_entry(
        cfg_v,
        pdb_id=str(pdb_scope).upper(),
        variant_label=variant_label,
        ph_token=ph_token,
    )


def queue_resume_scorch_only_owned_scopes(
    *,
    cfg_v: ConfigDict,
    dist_ctx: Any,
    label: str,
    combo_owner_task_id: dict[ScopeKey, int],
    combo_total_chunks: dict[ScopeKey, int],
    completion_bookkeeper: DeferredScorchCompletionBookkeeper,
    distributed_assigned_pdb_ids: Optional[set[str]] = None,
) -> tuple[bool, int, int]:
    """Queue SCORCH directly when a resumed distributed plan already finished docking.

    Returns ``(all_docking_complete, owned_pending, queued_or_marked)``.  The
    caller may skip the Vina chunk claim loop only when ``all_docking_complete``
    is true, because otherwise the resume still has real docking work to do.
    """
    manifest_entries = cfg_v.get("_RESUME_MANIFEST_PROTEINS")
    if not isinstance(manifest_entries, dict) or not combo_total_chunks:
        return False, 0, 0

    completed_owned_pending: dict[ScopeKey, int] = {}
    all_complete = True
    owned_complete = 0
    owned_pending = 0
    for scope, total_chunks_raw in sorted(combo_total_chunks.items()):
        total_chunks = int(total_chunks_raw or 0)
        if total_chunks <= 0:
            continue
        entry = _resume_manifest_scope_entry(cfg_v, scope=scope, fallback_label=label)
        if not isinstance(entry, dict) or not _resume_docking_outputs_complete(entry):
            all_complete = False
            continue
        owner_task = int(combo_owner_task_id.get(scope, -1))
        if owner_task != int(dist_ctx.task_id):
            continue
        owned_complete += 1
        if distributed_assigned_pdb_ids is not None:
            distributed_assigned_pdb_ids.add(str(scope[0]).upper())
        if _resume_scorch_required(cfg_v, entry) and _resume_scorch_still_pending(entry):
            completed_owned_pending[scope] = total_chunks
            owned_pending += 1

    queued_before = len(completion_bookkeeper.combo_scorch_submitted)
    if completed_owned_pending:
        if completion_bookkeeper.scorch_queue_service is None:
            logging.warning(
                "[resume.scorch-only.distributed] action=queue_unavailable_continue_chunks variant=%s task_id=%d owned_pending=%d",
                str(label).upper(),
                int(dist_ctx.task_id),
                int(owned_pending),
            )
            return False, int(owned_pending), 0
        completion_bookkeeper.flush_owned_deferred_submissions(
            combo_completed_chunks=completed_owned_pending,
            combo_failed_chunks={},
        )
    queued_after = len(completion_bookkeeper.combo_scorch_submitted)
    queued = max(0, int(queued_after) - int(queued_before))
    logging.info(
        "[resume.scorch-only.distributed] action=%s variant=%s task_id=%d owned_complete=%d owned_pending=%d queued=%d total_scopes=%d",
        "short_circuit_ready" if all_complete else "partial_resume_continue_chunks",
        str(label).upper(),
        int(dist_ctx.task_id),
        int(owned_complete),
        int(owned_pending),
        int(queued),
        int(len(combo_total_chunks)),
    )
    return bool(all_complete), int(owned_pending), int(queued)


def run_distributed_combo_variant(
    *,
    cfg_v: ConfigDict,
    dist_ctx: Any,
    run_id: str,
    label: str,
    combo_items: list[tuple[str, Optional[str]]],
    tokens: list[str],
    bar: Any,
    variant: Optional[str],
    dist_combo_hybrid_mode: bool,
    execution_mode_decision: Any,
    scorch_queue_service: Any,
    process_one: Callable[..., bool],
    recorded_terminal_chunk_ids: set[str],
    failed_entries: list[tuple[str, str, str, str, str]],
    run_scope_completion_times: dict[tuple[str, str, str], float],
    distributed_assigned_pdb_ids: set[str],
    run_chunk_rebalance_count: int,
) -> int:
    target_ligands = int(cfg_v.get("_BENCH_SMALL_TARGET_LIGANDS", 0) or 0)
    sample_seed = int(cfg_v.get("_BENCH_SMALL_SEED", 1337) or 1337)
    _ensure_distributed_manifest_preflight(
        cfg_v,
        run_id=run_id,
        dist_ctx=dist_ctx,
        variant_label=label.upper(),
        combo_items=combo_items,
        run_tokens=tokens,
    )
    plan_started = time.perf_counter()
    chunk_plan = _load_or_create_variant_chunk_plan(
        cfg_v,
        run_id=run_id,
        dist_ctx=dist_ctx,
        variant_label=label.upper(),
        combo_items=combo_items,
        run_tokens=tokens,
        target_ligands=target_ligands,
        sample_seed=sample_seed,
    )
    logging.info(
        "[distributed.chunk-plan] variant=%s chunks=%d elapsed_s=%.2f",
        label.upper(),
        len(chunk_plan),
        time.perf_counter() - plan_started,
    )
    assigned_chunks = [
        c
        for c in chunk_plan
        if int(c.get("assigned_task_id", -1)) == int(dist_ctx.task_id)
    ]
    bar.reset(total=max(1, len(assigned_chunks)))
    logging.info(
        "[distributed.combo-chunk] variant=%s total_chunks=%d assigned=%d task_id=%d/%d",
        label.upper(),
        len(chunk_plan),
        len(assigned_chunks),
        dist_ctx.task_id,
        dist_ctx.task_count,
    )

    by_chunk_id, combo_owner_task_id, combo_total_chunks, combo_total_ligands = (
        build_scope_maps(
            chunk_plan=chunk_plan,
            label=label,
            dist_ctx=dist_ctx,
            owner_task_resolver=_distributed_combo_owner_task_id,
        )
    )
    combo_runtime_mode = resolve_combo_runtime_mode(
        requested_mode=(
            "distributed_combo_hybrid"
            if dist_combo_hybrid_mode
            else "distributed_combo"
        ),
        scope_totals=combo_total_chunks,
    )
    combo_hybrid_active = bool(combo_runtime_mode.hybrid_enabled)
    if dist_combo_hybrid_mode and not combo_hybrid_active:
        logging.info(
            "[execution.mode] requested=distributed_combo_hybrid effective=%s reason=%s variant=%s",
            combo_runtime_mode.effective_mode,
            combo_runtime_mode.reason,
            label.upper(),
        )
    variant_param_for_scorch = None if variant is None else str(variant)

    combo_scope_total_weight, combo_total_weight = build_combo_weight_totals(
        by_chunk_id=by_chunk_id,
        label=label,
    )
    baseline_first_complete_sec: Optional[float] = None
    local_workers_for_baseline = int(_distributed_chunk_local_workers(cfg_v))
    if combo_scope_total_weight:
        baseline_first_complete_sec = float(
            min(combo_scope_total_weight.values())
            / max(1, local_workers_for_baseline)
        )
    baseline_makespan_sec: Optional[float] = None
    if combo_total_weight > 0:
        baseline_makespan_sec = float(
            combo_total_weight / max(1, local_workers_for_baseline)
        )

    completion_bookkeeper = build_deferred_scorch_completion_bookkeeper(
        cfg_v=cfg_v,
        run_id=str(run_id),
        label=str(label),
        dist_ctx=dist_ctx,
        combo_total_chunks=combo_total_chunks,
        combo_total_ligands=combo_total_ligands,
        combo_owner_task_id=combo_owner_task_id,
        scorch_queue_service=scorch_queue_service,
        variant_param_for_scorch=variant_param_for_scorch,
        by_chunk_id=by_chunk_id,
    )
    resume_docking_complete, resume_owned_pending, resume_queued = (
        queue_resume_scorch_only_owned_scopes(
            cfg_v=cfg_v,
            dist_ctx=dist_ctx,
            label=str(label),
            combo_owner_task_id=combo_owner_task_id,
            combo_total_chunks=combo_total_chunks,
            completion_bookkeeper=completion_bookkeeper,
            distributed_assigned_pdb_ids=distributed_assigned_pdb_ids,
        )
    )
    if resume_docking_complete:
        logging.info(
            "[resume.scorch-only.distributed] action=skip_chunk_claim_loop variant=%s task_id=%d owned_pending=%d queued=%d",
            str(label).upper(),
            int(dist_ctx.task_id),
            int(resume_owned_pending),
            int(resume_queued),
        )
        return int(run_chunk_rebalance_count)
    claim_loop_context = build_claim_loop_context(
        cfg_v=cfg_v,
        dist_ctx=dist_ctx,
        run_id=str(run_id),
        label=str(label),
        variant=variant,
        by_chunk_id=by_chunk_id,
        combo_owner_task_id=combo_owner_task_id,
        combo_total_chunks=combo_total_chunks,
        combo_hybrid_active=bool(combo_hybrid_active),
        execution_mode_decision=execution_mode_decision,
        combo_total_weight=float(combo_total_weight),
        baseline_first_complete_sec=baseline_first_complete_sec,
        baseline_makespan_sec=baseline_makespan_sec,
        process_one=process_one,
        distributed_assigned_pdb_ids=distributed_assigned_pdb_ids,
        recorded_terminal_chunk_ids=recorded_terminal_chunk_ids,
        failed_entries=failed_entries,
        bar=bar,
    )
    claim_loop_result = run_chunk_claim_loop(
        context=claim_loop_context,
        completion_bookkeeper=completion_bookkeeper,
    )

    return finalize_variant_claim_loop(
        by_chunk_id=by_chunk_id,
        label=label,
        dist_ctx=dist_ctx,
        cfg_v=cfg_v,
        run_id=str(run_id),
        combo_owner_task_id=combo_owner_task_id,
        combo_total_chunks=combo_total_chunks,
        recorded_terminal_chunk_ids=recorded_terminal_chunk_ids,
        failed_entries=failed_entries,
        completion_bookkeeper=completion_bookkeeper,
        claim_loop_rebalance_count=int(claim_loop_result.rebalance_count),
        run_chunk_rebalance_count=int(run_chunk_rebalance_count),
    )
