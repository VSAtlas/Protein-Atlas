from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor
from typing import Optional

from cli.distributed_chunk_claim_cache import DistributedChunkClaimCache
from cli.distributed_chunk_runtime_rebalance import HybridRebalanceState
from cli.distributed_chunk_runtime_state import ChunkStateTracker
from cli.distributed_context import read_chunk_result
from cli.distributed_runtime_claim_select import rebuild_candidate_order_for_session
from cli.distributed_runtime_claim_types import (
    ClaimLoopConfig,
    ClaimLoopContext,
    ClaimLoopResult,
    ClaimLoopSession,
)

__all__ = [
    "ClaimLoopConfig",
    "ClaimLoopContext",
    "ClaimLoopResult",
    "run_chunk_claim_loop",
]
from cli.distributed_runtime_claim_workers import chunk_worker_loop
from cli.distributed_runtime_completion import DeferredScorchCompletionBookkeeper


def run_chunk_claim_loop(
    *,
    context: ClaimLoopContext,
    completion_bookkeeper: DeferredScorchCompletionBookkeeper,
    config: Optional[ClaimLoopConfig] = None,
) -> ClaimLoopResult:
    cfg = config or ClaimLoopConfig()
    hybrid_lane_fraction = float(
        context.execution_mode_decision.completion_lane_fraction
        if context.combo_hybrid_active
        else 0.30
    )
    hybrid_rebalance_sec = int(max(60, context.execution_mode_decision.rebalance_sec))
    hybrid_rebalance_chunks = int(
        max(16, context.execution_mode_decision.rebalance_chunks)
    )
    use_claim_cache = not bool(context.combo_hybrid_active)
    initial_claim_cache: Optional[DistributedChunkClaimCache] = None
    if use_claim_cache:
        initial_claim_cache = DistributedChunkClaimCache(
            candidate_ids=[],
            max_attempts=int(cfg.chunk_max_attempts),
        )
    hybrid_state = HybridRebalanceState(
        candidate_ids=[],
        claim_cache=initial_claim_cache,
        lane_fraction=float(hybrid_lane_fraction),
        last_rebalance_ts=float(time.time()),
        chunks_since_rebalance=0,
        rebalance_count=0,
        consecutive_makespan_violations=0,
        promoted_scopes=set(),
    )
    chunk_state_tracker = ChunkStateTracker(
        chunk_ids=list(context.by_chunk_id.keys()),
        dist_ctx=context.dist_ctx,
        cfg_v=context.cfg_v,
        read_chunk_result=read_chunk_result,
        max_attempts=int(cfg.chunk_max_attempts),
        reconcile_sec=float(cfg.chunk_state_reconcile_sec),
        use_claim_cache=bool(use_claim_cache),
        claim_cache=hybrid_state.claim_cache,
    )
    session = ClaimLoopSession(
        context=context,
        cfg=cfg,
        completion_bookkeeper=completion_bookkeeper,
        chunk_state_tracker=chunk_state_tracker,
        hybrid_state=hybrid_state,
        claim_cursor=0,
        hybrid_completed_weight=0.0,
        active_worker_count=0,
        use_claim_cache=bool(use_claim_cache),
        hybrid_start_ts=float(completion_bookkeeper.hybrid_start_ts),
        hybrid_rebalance_sec=int(hybrid_rebalance_sec),
        hybrid_rebalance_chunks=int(hybrid_rebalance_chunks),
    )
    rebuild_candidate_order_for_session(session, force=True, reason="initial")
    chunk_state_tracker.refresh(force=True)

    with ThreadPoolExecutor(max_workers=max(1, int(context.local_chunk_workers))) as pool:
        futures = [
            pool.submit(chunk_worker_loop, session, int(worker_idx))
            for worker_idx in range(max(1, int(context.local_chunk_workers)))
        ]
        for fut in futures:
            fut.result()

    return ClaimLoopResult(rebalance_count=int(max(0, session.hybrid_state.rebalance_count)))
