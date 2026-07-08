from __future__ import annotations

from cli.hybrid_chunk_scheduler import build_hybrid_chunk_order, chunk_scope_key


def _chunk(
    chunk_id: str,
    pdb_id: str,
    *,
    weight: float,
    assigned_task_id: int,
    ph_tag: str = "base",
) -> dict[str, object]:
    return {
        "chunk_id": chunk_id,
        "pdb_id": pdb_id,
        "pdb_file": f"{pdb_id}.pdb",
        "ph_tag": ph_tag,
        "weight": float(weight),
        "assigned_task_id": int(assigned_task_id),
    }


def test_lane_scope_selection_and_promotion_order() -> None:
    assigned = [
        _chunk("c1", "P001", weight=10.0, assigned_task_id=3),
        _chunk("c2", "P002", weight=9.0, assigned_task_id=3),
        _chunk("c3", "P003", weight=8.0, assigned_task_id=3),
        _chunk("c4", "P004", weight=7.0, assigned_task_id=3),
    ]
    progress = {
        chunk_scope_key(assigned[0], variant_label="HOLO"): 0.80,
        chunk_scope_key(assigned[1], variant_label="HOLO"): 0.65,
        chunk_scope_key(assigned[2], variant_label="HOLO"): 0.30,
        chunk_scope_key(assigned[3], variant_label="HOLO"): 0.10,
    }
    owner = {
        chunk_scope_key(c, variant_label="HOLO"): 3 for c in assigned
    }
    order = build_hybrid_chunk_order(
        assigned_chunks=assigned,
        steal_chunks=[],
        task_id=3,
        combo_owner_task_id=owner,
        variant_label="HOLO",
        scope_progress_ratio=progress,
        completion_lane_fraction=0.50,
        use_completion_lane=True,
    )
    assert len(order.promoted_scopes) == 2
    promoted_ids = {scope[0] for scope in order.promoted_scopes}
    assert promoted_ids == {"P001", "P002"}


def test_completion_lane_drains_promoted_scopes_before_round_robin_rest() -> None:
    assigned = [
        _chunk("p1a", "P001", weight=10.0, assigned_task_id=3),
        _chunk("p1b", "P001", weight=9.0, assigned_task_id=3),
        _chunk("p1c", "P001", weight=8.0, assigned_task_id=3),
        _chunk("p2a", "P002", weight=7.0, assigned_task_id=3),
        _chunk("p2b", "P002", weight=6.0, assigned_task_id=3),
    ]
    progress = {
        chunk_scope_key(assigned[0], variant_label="HOLO"): 0.75,
        chunk_scope_key(assigned[3], variant_label="HOLO"): 0.05,
    }
    owner = {chunk_scope_key(c, variant_label="HOLO"): 3 for c in assigned}
    order = build_hybrid_chunk_order(
        assigned_chunks=assigned,
        steal_chunks=[],
        task_id=3,
        combo_owner_task_id=owner,
        variant_label="HOLO",
        scope_progress_ratio=progress,
        completion_lane_fraction=0.20,
        use_completion_lane=True,
    )

    assert {scope[0] for scope in order.promoted_scopes} == {"P001"}
    assert order.preferred_ids[:3] == ["p1a", "p1b", "p1c"]


def test_prepared_unstarted_scopes_sort_before_unprepared_rest() -> None:
    assigned = [
        _chunk("slow", "P001", weight=100.0, assigned_task_id=3),
        _chunk("ready", "P002", weight=1.0, assigned_task_id=3),
    ]
    ready_scope = chunk_scope_key(assigned[1], variant_label="HOLO")
    owner = {chunk_scope_key(c, variant_label="HOLO"): 3 for c in assigned}
    order = build_hybrid_chunk_order(
        assigned_chunks=assigned,
        steal_chunks=[],
        task_id=3,
        combo_owner_task_id=owner,
        variant_label="HOLO",
        scope_progress_ratio={},
        prep_ready_scopes={ready_scope},
        completion_lane_fraction=0.30,
        use_completion_lane=True,
    )

    assert order.preferred_ids[0] == "ready"


def test_deterministic_claim_order_stability() -> None:
    assigned = [
        _chunk("cA", "PAAA", weight=8.0, assigned_task_id=2),
        _chunk("cB", "PBBB", weight=7.0, assigned_task_id=2),
        _chunk("cC", "PCCC", weight=6.0, assigned_task_id=2),
    ]
    steal = [
        _chunk("sA", "SAAA", weight=5.0, assigned_task_id=1),
        _chunk("sB", "SBBB", weight=4.0, assigned_task_id=1),
    ]
    owner = {chunk_scope_key(c, variant_label="LEGACY"): 2 for c in assigned}
    progress = {
        chunk_scope_key(assigned[0], variant_label="LEGACY"): 0.55,
        chunk_scope_key(assigned[1], variant_label="LEGACY"): 0.15,
        chunk_scope_key(assigned[2], variant_label="LEGACY"): 0.05,
    }
    one = build_hybrid_chunk_order(
        assigned_chunks=assigned,
        steal_chunks=steal,
        task_id=2,
        combo_owner_task_id=owner,
        variant_label="LEGACY",
        scope_progress_ratio=progress,
        completion_lane_fraction=0.30,
        use_completion_lane=True,
    )
    two = build_hybrid_chunk_order(
        assigned_chunks=assigned,
        steal_chunks=steal,
        task_id=2,
        combo_owner_task_id=owner,
        variant_label="LEGACY",
        scope_progress_ratio=progress,
        completion_lane_fraction=0.30,
        use_completion_lane=True,
    )
    assert one.preferred_ids == two.preferred_ids
    assert one.steal_ids == two.steal_ids
