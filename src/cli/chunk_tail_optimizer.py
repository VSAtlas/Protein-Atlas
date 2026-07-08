from __future__ import annotations

import hashlib
import math
import os
from typing import Any, Sequence

from config.normalize import _to_bool


def dynamic_rechunk_enabled() -> bool:
    return _to_bool(os.environ.get("ATLAS_DYNAMIC_RECHUNK"), True)


def chunk_hedging_enabled(*, auto_enable: bool = False) -> bool:
    raw = os.environ.get("ATLAS_CHUNK_HEDGING")
    if raw is None or not str(raw).strip():
        return bool(auto_enable)
    token = str(raw).strip().lower()
    if token in {"auto", "adaptive"}:
        return bool(auto_enable)
    return _to_bool(raw, False)


def _rechunk_target_size(default_target: int) -> int:
    raw = os.environ.get("ATLAS_DYNAMIC_RECHUNK_TARGET_SIZE", "")
    try:
        val = int(raw)
    except Exception:
        val = int(default_target)
    return max(16, min(512, int(val)))


def _hash_chunk_id(*parts: str) -> str:
    raw = "|".join(str(p) for p in parts)
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:20]


def _rebalance_assignments(
    chunks: list[dict[str, Any]],
    *,
    task_ids: Sequence[int],
) -> None:
    if not task_ids:
        return
    unique_task_ids = sorted({int(tid) for tid in task_ids})
    loads = {int(tid): 0.0 for tid in unique_task_ids}
    for item in sorted(chunks, key=lambda x: float(x.get("weight", 0.0)), reverse=True):
        target_tid = min(unique_task_ids, key=lambda tid: loads.get(int(tid), 0.0))
        item["assigned_task_id"] = int(target_tid)
        loads[int(target_tid)] = loads.get(int(target_tid), 0.0) + float(
            item.get("weight", 0.0)
        )


def optimize_chunks_for_tail(
    chunks: list[dict[str, Any]],
    *,
    run_id: str,
    variant_label: str,
    task_ids: Sequence[int],
    task_count: int,
    local_workers_per_task: int,
    default_target_size: int,
    min_chunk_size: int,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """
    Adaptive pre-execution rechunk:
      - increases chunk multiplicity toward a tail-friendly target
      - splits only oversized chunks (>= 2 * min_chunk_size)
      - rebalances assignment by LPT after splits
    """
    if not dynamic_rechunk_enabled():
        return list(chunks), {
            "enabled": 0,
            "split_count": 0,
            "chunks_before": int(len(chunks)),
            "chunks_after": int(len(chunks)),
        }
    if not chunks:
        return [], {
            "enabled": 1,
            "split_count": 0,
            "chunks_before": 0,
            "chunks_after": 0,
        }

    target_size = _rechunk_target_size(default_target_size)
    workers_total = max(1, int(task_count) * max(1, int(local_workers_per_task)))
    desired_chunks = max(len(chunks), min(len(chunks) * 3, workers_total * 3))
    out = [dict(item) for item in chunks]
    split_count = 0
    round_idx = 0

    while len(out) < desired_chunks:
        candidates = [
            item
            for item in out
            if int(item.get("ligand_count", 0) or 0) >= max(
                int(target_size * 2), int(min_chunk_size * 2)
            )
        ]
        if not candidates:
            break
        parent = max(candidates, key=lambda item: int(item.get("ligand_count", 0) or 0))
        bases = [str(x) for x in (parent.get("ligand_bases") or []) if str(x).strip()]
        if len(bases) < max(2, int(min_chunk_size * 2)):
            break
        n_parts = max(2, int(math.ceil(float(len(bases)) / float(target_size))))
        part_size = int(math.ceil(float(len(bases)) / float(n_parts)))
        if part_size < max(1, int(min_chunk_size)):
            break

        parent_id = str(parent.get("chunk_id") or "")
        parent_weight = max(0.0, float(parent.get("weight", 0.0)))
        out = [item for item in out if str(item.get("chunk_id") or "") != parent_id]
        for part_idx in range(n_parts):
            start = part_idx * part_size
            chunk_bases = bases[start : start + part_size]
            if not chunk_bases:
                continue
            child = dict(parent)
            child["ligand_bases"] = list(chunk_bases)
            child["ligand_count"] = int(len(chunk_bases))
            child["weight"] = float(parent_weight * (len(chunk_bases) / max(1, len(bases))))
            child["chunk_id"] = _hash_chunk_id(
                run_id,
                variant_label,
                str(parent.get("pdb_id") or ""),
                str(parent.get("ph_tag") or "base"),
                str(parent.get("run_mode") or ""),
                str(parent.get("library_name") or ""),
                str(parent.get("library_root") or ""),
                str(parent_id),
                f"rechunk:{round_idx}:{part_idx}",
                ",".join(chunk_bases),
            )
            child["dynamic_rechunk"] = True
            child["dynamic_rechunk_parent"] = parent_id
            out.append(child)
        split_count += 1
        round_idx += 1
        if round_idx > 1024:
            break

    _rebalance_assignments(out, task_ids=task_ids)
    return out, {
        "enabled": 1,
        "split_count": int(split_count),
        "chunks_before": int(len(chunks)),
        "chunks_after": int(len(out)),
    }


def build_candidate_ids(
    *,
    preferred_ids: Sequence[str],
    steal_ids: Sequence[str],
    auto_enable_hedging: bool = False,
) -> tuple[list[str], dict[str, int]]:
    """
    Build claim candidate order.
    With hedging disabled: preferred first, then steal.
    With hedging enabled: interleave preferred/steal to reduce straggler tails.
    """
    pref = [str(x).strip() for x in preferred_ids if str(x).strip()]
    steal = [str(x).strip() for x in steal_ids if str(x).strip()]
    if not chunk_hedging_enabled(auto_enable=bool(auto_enable_hedging)):
        return [*pref, *steal], {
            "hedging_enabled": 0,
            "preferred": int(len(pref)),
            "steal": int(len(steal)),
            "auto_enabled": int(bool(auto_enable_hedging)),
        }
    out: list[str] = []
    i = 0
    j = 0
    while i < len(pref) or j < len(steal):
        if i < len(pref):
            out.append(pref[i])
            i += 1
        if j < len(steal):
            out.append(steal[j])
            j += 1
    return out, {
        "hedging_enabled": 1,
        "preferred": int(len(pref)),
        "steal": int(len(steal)),
        "auto_enabled": int(bool(auto_enable_hedging)),
    }
