from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Any, Mapping, Sequence


ScopeKey = tuple[str, str, str]


@dataclass(frozen=True)
class HybridChunkOrder:
    preferred_ids: list[str]
    steal_ids: list[str]
    promoted_scopes: list[ScopeKey]
    promoted_chunk_count: int
    active_scope_count: int


def chunk_scope_key(
    chunk_payload: Mapping[str, Any], *, variant_label: str
) -> ScopeKey:
    pdb_scope = (
        str(chunk_payload.get("pdb_id") or "").strip().upper()
        or os.path.splitext(os.path.basename(str(chunk_payload.get("pdb_file") or "")))[
            0
        ].upper()
    )
    raw_ph = str(chunk_payload.get("ph_tag") or "").strip().lower()
    ph_scope = raw_ph if raw_ph else "base"
    return (
        str(pdb_scope),
        str(variant_label or "").strip().upper() or "BASE",
        str(ph_scope),
    )


def _ordered_chunk_ids(
    chunk_items: Sequence[Mapping[str, Any]],
    *,
    task_id: int,
    combo_owner_task_id: Mapping[ScopeKey, int],
    variant_label: str,
    prefer_local_owner: bool,
    promoted_scopes: set[ScopeKey],
    prep_ready_scopes: set[ScopeKey],
    scope_progress_ratio: Mapping[ScopeKey, float],
    use_completion_lane: bool,
) -> list[str]:
    def _rank_key(item: Mapping[str, Any]) -> tuple[float, float, float, float, float]:
        scope = chunk_scope_key(item, variant_label=variant_label)
        owner_task = int(combo_owner_task_id.get(scope, -1))
        owner_penalty = (
            0.0 if (prefer_local_owner and owner_task == int(task_id)) else 1.0
        )
        lane_penalty = (
            0.0 if (use_completion_lane and scope in promoted_scopes) else 1.0
        )
        prep_ready_penalty = 0.0 if scope in prep_ready_scopes else 1.0
        ratio = float(scope_progress_ratio.get(scope, 0.0))
        return (
            lane_penalty,
            owner_penalty,
            prep_ready_penalty,
            -ratio,
            -float(item.get("weight", 0.0)),
        )

    ranked = sorted(chunk_items, key=_rank_key)
    buckets: dict[ScopeKey, list[str]] = {}
    scope_order: list[ScopeKey] = []
    for item in ranked:
        chunk_id = str(item.get("chunk_id") or "").strip()
        if not chunk_id:
            continue
        scope = chunk_scope_key(item, variant_label=variant_label)
        if scope not in buckets:
            buckets[scope] = []
            scope_order.append(scope)
        buckets[scope].append(chunk_id)

    if scope_order and not promoted_scopes and not prep_ready_scopes:
        shift = int(task_id) % len(scope_order)
        if shift:
            scope_order = scope_order[shift:] + scope_order[:shift]

    ordered: list[str] = []
    if use_completion_lane and promoted_scopes:
        promoted_order = [scope for scope in scope_order if scope in promoted_scopes]
        for scope in promoted_order:
            ordered.extend(buckets.get(scope, []))
        scope_order = [scope for scope in scope_order if scope not in promoted_scopes]

    idx_map = {scope: 0 for scope in scope_order}
    while scope_order:
        emitted = False
        for scope in scope_order:
            idx = int(idx_map.get(scope, 0))
            bucket = buckets.get(scope, [])
            if idx >= len(bucket):
                continue
            ordered.append(bucket[idx])
            idx_map[scope] = idx + 1
            emitted = True
        if not emitted:
            break
    return ordered


def build_hybrid_chunk_order(
    *,
    assigned_chunks: Sequence[Mapping[str, Any]],
    steal_chunks: Sequence[Mapping[str, Any]],
    task_id: int,
    combo_owner_task_id: Mapping[ScopeKey, int],
    variant_label: str,
    scope_progress_ratio: Mapping[ScopeKey, float],
    completion_lane_fraction: float,
    use_completion_lane: bool,
    prep_ready_scopes: set[ScopeKey] | None = None,
) -> HybridChunkOrder:
    prep_ready_scope_set = set(prep_ready_scopes or set())
    active_scopes = sorted(
        {
            chunk_scope_key(item, variant_label=variant_label)
            for item in [*assigned_chunks, *steal_chunks]
            if isinstance(item, Mapping)
        }
    )
    promoted_scopes: list[ScopeKey] = []
    promoted_scope_set: set[ScopeKey] = set()
    if use_completion_lane and active_scopes:
        active_in_progress = [
            scope
            for scope in active_scopes
            if 0.0 < float(scope_progress_ratio.get(scope, 0.0)) < 1.0
        ]
        ranked_scopes = sorted(
            active_in_progress,
            key=lambda scope: float(scope_progress_ratio.get(scope, 0.0)),
            reverse=True,
        )
        lane_scope_count = int(
            max(
                1,
                min(
                    len(ranked_scopes),
                    math.ceil(
                        float(len(active_scopes))
                        * float(max(0.2, min(0.5, completion_lane_fraction)))
                    ),
                ),
            )
        )
        promoted_scopes = ranked_scopes[:lane_scope_count]
        promoted_scope_set = set(promoted_scopes)

    preferred_ids = _ordered_chunk_ids(
        [item for item in assigned_chunks if isinstance(item, Mapping)],
        task_id=int(task_id),
        combo_owner_task_id=combo_owner_task_id,
        variant_label=variant_label,
        prefer_local_owner=True,
        promoted_scopes=promoted_scope_set,
        prep_ready_scopes=prep_ready_scope_set,
        scope_progress_ratio=scope_progress_ratio,
        use_completion_lane=bool(use_completion_lane),
    )
    steal_ids = _ordered_chunk_ids(
        [item for item in steal_chunks if isinstance(item, Mapping)],
        task_id=int(task_id),
        combo_owner_task_id=combo_owner_task_id,
        variant_label=variant_label,
        prefer_local_owner=False,
        promoted_scopes=promoted_scope_set,
        prep_ready_scopes=prep_ready_scope_set,
        scope_progress_ratio=scope_progress_ratio,
        use_completion_lane=bool(use_completion_lane),
    )
    promoted_chunk_count = 0
    if promoted_scope_set:
        promoted_chunk_count = sum(
            1
            for item in assigned_chunks
            if isinstance(item, Mapping)
            and chunk_scope_key(item, variant_label=variant_label) in promoted_scope_set
        )
    return HybridChunkOrder(
        preferred_ids=preferred_ids,
        steal_ids=steal_ids,
        promoted_scopes=promoted_scopes,
        promoted_chunk_count=int(promoted_chunk_count),
        active_scope_count=int(len(active_scopes)),
    )
