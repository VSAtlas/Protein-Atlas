from __future__ import annotations

import heapq
import random

from cli.distributed_mode_policy import evaluate_balanced_constraints


def _generate_scope_chunks(
    *, scopes: int = 100, chunks_per_scope: int = 24, seed: int = 7
) -> tuple[dict[str, float], dict[str, list[float]]]:
    rng = random.Random(seed)
    totals: dict[str, float] = {}
    chunks: dict[str, list[float]] = {}
    for idx in range(scopes):
        scope = f"P{idx:03d}"
        # Deliberately skewed to emulate heavy-tail protein runtimes.
        total = max(1.0, rng.lognormvariate(0.0, 1.0)) * 200.0
        totals[scope] = float(total)
        chunks[scope] = [float(total) / float(chunks_per_scope)] * int(chunks_per_scope)
    return totals, chunks


def _dispatch_order(order: list[tuple[str, float]], workers: int) -> tuple[float, float]:
    worker_heap = [(0.0, i) for i in range(max(1, int(workers)))]
    remaining: dict[str, int] = {}
    first_done: dict[str, float] = {}
    for scope, _ in order:
        remaining[scope] = int(remaining.get(scope, 0)) + 1
    for scope, duration in order:
        start_t, worker_id = heapq.heappop(worker_heap)
        finish_t = float(start_t) + float(duration)
        heapq.heappush(worker_heap, (finish_t, worker_id))
        remaining[scope] = int(remaining.get(scope, 1)) - 1
        if remaining[scope] <= 0:
            first_done[scope] = finish_t
    makespan = max(float(t) for t, _ in worker_heap)
    first_complete = min(float(t) for t in first_done.values())
    return float(first_complete), float(makespan)


def _simulate_per_protein(
    totals: dict[str, float], *, workers: int
) -> tuple[float, float]:
    worker_heap = [(0.0, i) for i in range(max(1, int(workers)))]
    completions: list[float] = []
    for _scope, total in sorted(totals.items(), key=lambda kv: kv[1], reverse=True):
        start_t, worker_id = heapq.heappop(worker_heap)
        finish_t = float(start_t) + float(total)
        heapq.heappush(worker_heap, (finish_t, worker_id))
        completions.append(finish_t)
    return float(min(completions)), float(max(t for t, _ in worker_heap))


def _simulate_combo_order(chunks: dict[str, list[float]]) -> list[tuple[str, float]]:
    order: list[tuple[str, float]] = []
    mutable = {scope: list(vals) for scope, vals in chunks.items()}
    scope_order = list(mutable.keys())
    while True:
        emitted = False
        for scope in scope_order:
            items = mutable.get(scope) or []
            if not items:
                continue
            order.append((scope, float(items.pop())))
            emitted = True
        if not emitted:
            break
    return order


def _simulate_hybrid_order(
    chunks: dict[str, list[float]], *, lane_fraction: float
) -> list[tuple[str, float]]:
    order: list[tuple[str, float]] = []
    mutable = {scope: list(vals) for scope, vals in chunks.items()}
    total = {scope: len(vals) for scope, vals in mutable.items()}
    done = {scope: 0 for scope in mutable}
    round_robin = list(mutable.keys())
    rr_idx = 0
    lane_every = max(2, int(round(1.0 / max(0.05, min(0.9, float(lane_fraction))))))
    step = 0
    while True:
        active = [scope for scope, vals in mutable.items() if vals]
        if not active:
            break
        picked: str | None = None
        if step % lane_every == 0:
            in_progress = [
                scope
                for scope in active
                if 0 < int(done[scope]) < int(total[scope])
            ]
            if in_progress:
                picked = max(
                    in_progress,
                    key=lambda scope: float(done[scope]) / float(total[scope]),
                )
        if picked is None:
            for _ in range(len(round_robin)):
                candidate = round_robin[rr_idx % len(round_robin)]
                rr_idx += 1
                if mutable[candidate]:
                    picked = candidate
                    break
        if picked is None:
            break
        order.append((picked, float(mutable[picked].pop())))
        done[picked] = int(done[picked]) + 1
        step += 1
    return order


def test_synthetic_ab_constraints_balanced_objective() -> None:
    workers = 32
    totals, chunks = _generate_scope_chunks(scopes=100, chunks_per_scope=24, seed=7)
    per_first, per_make = _simulate_per_protein(totals, workers=workers)
    combo_first, combo_make = _dispatch_order(_simulate_combo_order(chunks), workers)
    hybrid_first, hybrid_make = _dispatch_order(
        _simulate_hybrid_order(chunks, lane_fraction=0.30), workers
    )

    assert per_first <= combo_first
    status = evaluate_balanced_constraints(
        combo_first_complete_s=combo_first,
        hybrid_first_complete_s=hybrid_first,
        combo_makespan_s=combo_make,
        hybrid_makespan_s=hybrid_make,
    )
    assert status.constraint_a_ok is True
    assert status.constraint_b_ok is True
    assert hybrid_first < combo_first
    assert hybrid_make <= combo_make * 1.15
    assert per_make >= 0.0
