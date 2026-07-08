from __future__ import annotations

from typing import Any, Mapping

from cli.distributed_chunk_runtime_rebalance import HybridRebalanceState
from cli.distributed_chunk_runtime_state import chunk_terminal
from cli.distributed_context import write_chunk_result
from cli.distributed_runtime_claim_types import ClaimLoopSession


def write_chunk_result_and_cache(
    session: ClaimLoopSession,
    *,
    chunk_id: str,
    status: str,
    payload: Mapping[str, Any],
) -> None:
    context = session.context
    cfg = session.cfg
    write_chunk_result(
        context.dist_ctx,
        context.cfg_v,
        chunk_id=chunk_id,
        status=status,
        payload=payload,
    )
    try:
        attempt_no = int(payload.get("attempt") or 1)
    except Exception:
        attempt_no = 1
    session.chunk_state_tracker.record_local_result(
        chunk_id=chunk_id,
        status=status,
        attempt=int(attempt_no),
    )
    st = str(status or "").strip().lower()
    terminal = bool(
        chunk_terminal(
            st,
            int(attempt_no),
            max_attempts=int(cfg.chunk_max_attempts),
        )
    )
    if terminal:
        session.hybrid_state = HybridRebalanceState(
            candidate_ids=list(session.hybrid_state.candidate_ids),
            claim_cache=session.hybrid_state.claim_cache,
            lane_fraction=float(session.hybrid_state.lane_fraction),
            last_rebalance_ts=float(session.hybrid_state.last_rebalance_ts),
            chunks_since_rebalance=int(session.hybrid_state.chunks_since_rebalance) + 1,
            rebalance_count=int(session.hybrid_state.rebalance_count),
            consecutive_makespan_violations=int(
                session.hybrid_state.consecutive_makespan_violations
            ),
            promoted_scopes=set(session.hybrid_state.promoted_scopes),
        )
    if st == "completed":
        chunk_payload = context.by_chunk_id.get(str(chunk_id))
        if isinstance(chunk_payload, dict):
            session.hybrid_completed_weight += float(
                chunk_payload.get("weight", 0.0) or 0.0
            )
