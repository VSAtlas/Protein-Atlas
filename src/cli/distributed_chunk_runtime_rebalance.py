from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Any, Mapping, Optional, Sequence

from cli.chunk_tail_optimizer import build_candidate_ids
from cli.distributed_chunk_claim_cache import DistributedChunkClaimCache
from cli.hybrid_chunk_scheduler import build_hybrid_chunk_order, chunk_scope_key

Scope = tuple[str, str, str]


@dataclass(frozen=True)
class HybridRebalanceState:
    candidate_ids: list[str]
    claim_cache: Optional[DistributedChunkClaimCache]
    lane_fraction: float
    last_rebalance_ts: float
    chunks_since_rebalance: int
    rebalance_count: int
    consecutive_makespan_violations: int
    promoted_scopes: set[Scope]


def scope_progress_ratio_map(
    *,
    combo_total_chunks: Mapping[Scope, int],
    combo_completed_chunks_live: Mapping[Scope, int],
) -> dict[Scope, float]:
    ratios: dict[Scope, float] = {}
    for scope, total in combo_total_chunks.items():
        total_i = int(total)
        if total_i <= 0:
            continue
        done = int(combo_completed_chunks_live.get(scope, 0))
        ratios[scope] = float(max(0.0, min(1.0, float(done) / float(total_i))))
    return ratios


def dedupe_ids(raw_ids: Sequence[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for token in raw_ids:
        cid = str(token or "").strip()
        if not cid or cid in seen:
            continue
        out.append(cid)
        seen.add(cid)
    return out


def rebuild_candidate_order(
    *,
    force: bool,
    reason: str,
    combo_hybrid_active: bool,
    state: HybridRebalanceState,
    by_chunk_id: Mapping[str, Mapping[str, Any]],
    dist_task_id: int,
    combo_owner_task_id: Mapping[Scope, int],
    label: str,
    combo_total_chunks: Mapping[Scope, int],
    combo_completed_chunks_live: Mapping[Scope, int],
    local_chunk_workers: int,
    hybrid_start_ts: float,
    baseline_makespan_sec: Optional[float],
    combo_total_weight: float,
    hybrid_completed_weight: float,
    hybrid_first_complete_elapsed: Optional[float],
    baseline_first_complete_sec: Optional[float],
    rebalance_sec: int,
    rebalance_chunks: int,
    use_claim_cache: bool,
    chunk_max_attempts: int,
    prep_ready_scopes: Optional[set[Scope]] = None,
) -> tuple[HybridRebalanceState, bool]:
    now = time.time()
    if not force and not bool(combo_hybrid_active):
        return state, False
    if (
        not force
        and bool(combo_hybrid_active)
        and (now - float(state.last_rebalance_ts)) < float(rebalance_sec)
        and int(state.chunks_since_rebalance) < int(rebalance_chunks)
    ):
        return state, False

    scope_ratios = scope_progress_ratio_map(
        combo_total_chunks=combo_total_chunks,
        combo_completed_chunks_live=combo_completed_chunks_live,
    )
    assigned_now: list[Mapping[str, Any]] = []
    steal_now: list[Mapping[str, Any]] = []
    variant_token = str(label or "").strip().upper() or "BASE"
    for item in by_chunk_id.values():
        if not isinstance(item, dict):
            continue
        scope = chunk_scope_key(item, variant_label=variant_token)
        assigned_to_task = int(item.get("assigned_task_id", -1)) == int(dist_task_id)
        owned_by_task = int(combo_owner_task_id.get(scope, -1)) == int(dist_task_id)
        if assigned_to_task or owned_by_task:
            assigned_now.append(item)
        else:
            steal_now.append(item)
    order_payload = build_hybrid_chunk_order(
        assigned_chunks=assigned_now,
        steal_chunks=steal_now,
        task_id=int(dist_task_id),
        combo_owner_task_id=combo_owner_task_id,
        variant_label=variant_token,
        scope_progress_ratio=scope_ratios,
        completion_lane_fraction=float(state.lane_fraction),
        use_completion_lane=bool(combo_hybrid_active),
        prep_ready_scopes=set(prep_ready_scopes or set()),
    )
    auto_enable_hedging = bool(
        combo_hybrid_active
        and len(order_payload.preferred_ids)
        <= max(2, int(local_chunk_workers) * 2)
        and len(order_payload.steal_ids) > 0
    )
    ordered_ids, hedging_stats = build_candidate_ids(
        preferred_ids=order_payload.preferred_ids,
        steal_ids=order_payload.steal_ids,
        auto_enable_hedging=auto_enable_hedging,
    )
    next_ids = dedupe_ids(ordered_ids)
    if not next_ids:
        next_ids = dedupe_ids([str(x) for x in by_chunk_id.keys()])

    moved_chunks = int(order_payload.promoted_chunk_count)
    promoted_scopes = set(order_payload.promoted_scopes)
    lane_fraction = float(state.lane_fraction)
    violations = int(state.consecutive_makespan_violations)

    lane_workers = int(
        max(
            1,
            min(
                int(local_chunk_workers),
                round(
                    float(local_chunk_workers) * float(max(0.2, min(0.5, lane_fraction)))
                ),
            ),
        )
    )
    logging.info(
        "[hybrid.lane] variant=%s local_workers=%d lane_workers=%d lane_fraction=%.2f active_scopes=%d",
        label.upper(),
        int(local_chunk_workers),
        int(lane_workers),
        float(lane_fraction),
        int(order_payload.active_scope_count),
    )

    if bool(combo_hybrid_active):
        elapsed = max(1e-6, float(time.time() - hybrid_start_ts))
        if baseline_makespan_sec and hybrid_completed_weight > 0:
            throughput = float(hybrid_completed_weight) / float(elapsed)
            if throughput > 1e-9:
                est_makespan = float(combo_total_weight) / float(throughput)
                makespan_degradation = (
                    float(est_makespan) / float(baseline_makespan_sec)
                ) - 1.0
                if makespan_degradation > 0.15:
                    violations += 1
                else:
                    violations = 0
                if violations >= 2 and float(lane_fraction) > 0.20:
                    lane_fraction = float(max(0.20, float(lane_fraction) - 0.10))
                    violations = 0
        if (
            hybrid_first_complete_elapsed is not None
            and baseline_first_complete_sec
            and int(order_payload.active_scope_count) >= 4
        ):
            target_first_complete = float(baseline_first_complete_sec) * 0.70
            if (
                float(hybrid_first_complete_elapsed) > float(target_first_complete)
                and float(lane_fraction) < 0.50
            ):
                lane_fraction = float(min(0.50, float(lane_fraction) + 0.10))

    next_claim_cache: Optional[DistributedChunkClaimCache] = state.claim_cache
    if use_claim_cache:
        next_claim_cache = DistributedChunkClaimCache(
            candidate_ids=next_ids,
            max_attempts=int(chunk_max_attempts),
        )

    logging.info(
        "[distributed.chunk.hedging] variant=%s task_id=%d enabled=%s preferred=%d steal=%d candidates=%d",
        label.upper(),
        int(dist_task_id),
        int(hedging_stats.get("hedging_enabled", 0)),
        int(hedging_stats.get("preferred", 0)),
        int(hedging_stats.get("steal", 0)),
        len(next_ids),
    )
    logging.info(
        "[hybrid.rebalance] variant=%s reason=%s moved_chunks=%d scopes_promoted=%d lane_fraction=%.2f",
        label.upper(),
        str(reason),
        int(moved_chunks),
        int(len(promoted_scopes)),
        float(lane_fraction),
    )
    for promoted_scope in sorted(promoted_scopes)[:12]:
        ratio = float(scope_ratios.get(promoted_scope, 0.0))
        if ratio <= 0.0 or ratio >= 1.0:
            eta_s = 0.0
        else:
            eta_s = float((time.time() - hybrid_start_ts) * ((1.0 - ratio) / ratio))
        logging.info(
            "[hybrid.progress] scope=%s ratio=%.3f eta_s=%.1f",
            "|".join(str(x) for x in promoted_scope),
            float(ratio),
            float(max(0.0, eta_s)),
        )

    next_state = HybridRebalanceState(
        candidate_ids=list(next_ids),
        claim_cache=next_claim_cache,
        lane_fraction=float(lane_fraction),
        last_rebalance_ts=time.time(),
        chunks_since_rebalance=0,
        rebalance_count=int(state.rebalance_count) + 1,
        consecutive_makespan_violations=int(violations),
        promoted_scopes=set(promoted_scopes),
    )
    return next_state, True
