from __future__ import annotations

import logging
from typing import Optional

from cli.distributed_chunk_runtime_rebalance import rebuild_candidate_order
from cli.distributed_context import read_chunk_result, try_claim_chunk
from cli.distributed_runtime_claim_prep import (
    chunk_prep_claimable,
    collect_prep_ready_scopes,
)
from cli.distributed_runtime_claim_types import ClaimLoopSession


def rebuild_candidate_order_for_session(
    session: ClaimLoopSession,
    *,
    force: bool,
    reason: str,
) -> None:
    context = session.context
    cfg = session.cfg
    with session.hybrid_control_lock:
        _, completed_snapshot, first_complete_elapsed = (
            session.completion_bookkeeper.snapshot_for_rebalance()
        )
        next_state, changed = rebuild_candidate_order(
            force=bool(force),
            reason=str(reason),
            combo_hybrid_active=bool(context.combo_hybrid_active),
            state=session.hybrid_state,
            by_chunk_id=context.by_chunk_id,
            dist_task_id=int(context.dist_ctx.task_id),
            combo_owner_task_id=context.combo_owner_task_id,
            label=context.label,
            combo_total_chunks=context.combo_total_chunks,
            combo_completed_chunks_live=completed_snapshot,
            local_chunk_workers=int(context.local_chunk_workers),
            prep_ready_scopes=collect_prep_ready_scopes(session),
            hybrid_start_ts=float(session.hybrid_start_ts),
            baseline_makespan_sec=context.baseline_makespan_sec,
            combo_total_weight=float(context.combo_total_weight),
            hybrid_completed_weight=float(session.hybrid_completed_weight),
            hybrid_first_complete_elapsed=first_complete_elapsed,
            baseline_first_complete_sec=context.baseline_first_complete_sec,
            rebalance_sec=int(session.hybrid_rebalance_sec),
            rebalance_chunks=int(session.hybrid_rebalance_chunks),
            use_claim_cache=bool(session.use_claim_cache),
            chunk_max_attempts=int(cfg.chunk_max_attempts),
        )
        if not changed:
            return
        session.hybrid_state = next_state
        with session.claim_cursor_lock:
            session.claim_cursor = 0
        session.chunk_state_tracker.set_claim_cache(session.hybrid_state.claim_cache)
        logging.info(
            "[hybrid.state] variant=%s lane_fraction=%.2f rebalance_count=%d",
            context.label.upper(),
            float(session.hybrid_state.lane_fraction),
            int(session.hybrid_state.rebalance_count),
        )


def claim_next_chunk(session: ClaimLoopSession) -> tuple[Optional[str], int]:
    context = session.context
    cfg = session.cfg
    rebuild_candidate_order_for_session(session, force=False, reason="periodic")
    if session.use_claim_cache and session.hybrid_state.claim_cache is not None:
        total_candidates = session.hybrid_state.claim_cache.candidate_count
        for _ in range(max(0, int(total_candidates))):
            candidate_id, prev_attempt = (
                session.hybrid_state.claim_cache.next_claim_candidate()
            )
            if not candidate_id:
                break
            if not chunk_prep_claimable(session, str(candidate_id)):
                continue
            if try_claim_chunk(
                context.dist_ctx,
                context.cfg_v,
                chunk_id=candidate_id,
                lease_sec=cfg.chunk_claim_lease_sec,
                retry_failed=True,
                max_attempts=cfg.chunk_max_attempts,
            ):
                return candidate_id, int(max(0, prev_attempt))
        return None, 0

    total_candidates = len(session.hybrid_state.candidate_ids)
    if total_candidates <= 0:
        return None, 0
    for _ in range(total_candidates):
        with session.claim_cursor_lock:
            candidate_id = session.hybrid_state.candidate_ids[
                session.claim_cursor % total_candidates
            ]
            session.claim_cursor += 1
        payload = read_chunk_result(context.dist_ctx, context.cfg_v, chunk_id=candidate_id)
        prev_attempt = 0
        if isinstance(payload, dict):
            status = str(payload.get("status") or "").strip().lower()
            try:
                prev_attempt = int(payload.get("attempt") or 1)
            except Exception:
                prev_attempt = 1
            if status == "completed":
                continue
            if status == "terminal_failed":
                continue
            if status == "failed" and int(prev_attempt) >= int(cfg.chunk_max_attempts):
                continue
        if not chunk_prep_claimable(session, str(candidate_id)):
            continue
        if try_claim_chunk(
            context.dist_ctx,
            context.cfg_v,
            chunk_id=candidate_id,
            lease_sec=cfg.chunk_claim_lease_sec,
            retry_failed=True,
            max_attempts=cfg.chunk_max_attempts,
        ):
            return candidate_id, int(max(0, prev_attempt))
    return None, 0
