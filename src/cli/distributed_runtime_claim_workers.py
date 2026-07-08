from __future__ import annotations

import logging
import time
from cli.distributed_runtime_claim_process import process_claimed_chunk
from cli.distributed_runtime_claim_select import claim_next_chunk
from cli.distributed_runtime_claim_types import ClaimLoopSession


def active_workers(session: ClaimLoopSession) -> int:
    with session.active_worker_lock:
        return int(session.active_worker_count)


def set_worker_active(session: ClaimLoopSession, delta: int) -> None:
    with session.active_worker_lock:
        session.active_worker_count = max(
            0, int(session.active_worker_count) + int(delta)
        )


def reconcile_provisional(
    session: ClaimLoopSession,
    reason: str,
    *,
    unresolved_chunks: int,
) -> None:
    context = session.context
    active = active_workers(session)
    idle = max(0, int(context.local_chunk_workers) - int(active))
    session.completion_bookkeeper.reconcile_provisional_owned_submissions(
        reason=reason,
        active_workers=int(active),
        idle_workers=int(idle),
        unresolved_chunks=int(unresolved_chunks),
    )


def chunk_worker_loop(session: ClaimLoopSession, worker_idx: int) -> None:
    context = session.context
    cfg = session.cfg
    idle_since = time.time()
    idle_wait_sec = float(cfg.chunk_idle_wait_sec)
    tail_mode = False
    tail_remaining = 0
    next_tail_probe = 0.0
    while True:
        if session.chunk_state_tracker.all_done():
            return
        claimed_id, claimed_prev_attempt = claim_next_chunk(session)
        if not claimed_id:
            session.completion_bookkeeper.reconcile_ready_owned_deferred_submissions(
                reason="idle_global_progress",
            )
            now = time.time()
            if now >= next_tail_probe:
                tail_remaining = session.chunk_state_tracker.unresolved_chunk_count()
                tail_mode = bool(
                    int(tail_remaining) <= max(2, int(context.local_chunk_workers) * 2)
                )
                next_tail_probe = now + 5.0
            reconcile_provisional(
                session,
                "idle_global_progress",
                unresolved_chunks=int(tail_remaining),
            )
            if time.time() - idle_since >= cfg.chunk_idle_timeout_sec:
                logging.warning(
                    "[distributed.combo-chunk] variant=%s task_id=%d worker=%d/%d waiting_for_chunks_s=%.1f unresolved_chunks=%d tail_mode=%s",
                    context.label.upper(),
                    context.dist_ctx.task_id,
                    int(worker_idx) + 1,
                    int(context.local_chunk_workers),
                    cfg.chunk_idle_timeout_sec,
                    int(tail_remaining),
                    str(bool(tail_mode)).lower(),
                )
                idle_since = time.time()
            sleep_for = 0.25 if tail_mode else idle_wait_sec
            time.sleep(sleep_for)
            idle_wait_cap = 2.0 if tail_mode else 5.0
            idle_wait_sec = min(idle_wait_cap, idle_wait_sec + 0.2)
            continue
        idle_since = time.time()
        idle_wait_sec = float(cfg.chunk_idle_wait_sec)
        set_worker_active(session, 1)
        try:
            process_claimed_chunk(session, claimed_id, int(claimed_prev_attempt))
        finally:
            set_worker_active(session, -1)
        session.completion_bookkeeper.reconcile_ready_owned_deferred_submissions(
            reason="post_chunk_global_progress",
        )
        reconcile_provisional(
            session,
            "post_chunk_global_progress",
            unresolved_chunks=session.chunk_state_tracker.unresolved_chunk_count(),
        )
