from __future__ import annotations

from cli.distributed_mode_policy import (
    evaluate_balanced_constraints,
    resolve_execution_mode,
)


def test_auto_threshold_boundaries() -> None:
    cfg = {}
    d12 = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=True,
        task_count=8,
        pdb_count=12,
        combo_disabled=False,
        requested_mode_raw="auto",
    )
    d13 = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=True,
        task_count=8,
        pdb_count=13,
        combo_disabled=False,
        requested_mode_raw="auto",
    )
    d60 = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=True,
        task_count=8,
        pdb_count=60,
        combo_disabled=False,
        requested_mode_raw="auto",
    )
    d61 = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=True,
        task_count=8,
        pdb_count=61,
        combo_disabled=False,
        requested_mode_raw="auto",
    )
    assert d12.effective_mode == "distributed_per_protein"
    assert d13.effective_mode == "distributed_combo_hybrid"
    assert d60.effective_mode == "distributed_per_protein"
    assert d61.effective_mode == "distributed_combo_hybrid"
    assert d13.reason == "auto_medium_dense_tasks"


def test_explicit_override_precedence_over_auto() -> None:
    cfg = {}
    decision = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=True,
        task_count=8,
        pdb_count=200,
        combo_disabled=False,
        requested_mode_raw="distributed_combo",
    )
    assert decision.effective_mode == "distributed_combo"
    assert decision.reason == "requested_explicit_mode"


def test_combo_disable_compatibility_wins() -> None:
    cfg = {}
    decision = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=True,
        task_count=8,
        pdb_count=200,
        combo_disabled=True,
        requested_mode_raw="distributed_combo_hybrid",
    )
    assert decision.effective_mode == "distributed_per_protein"
    assert decision.reason == "combo_disabled_env"


def test_balanced_objective_constraint_a_b() -> None:
    status_ok = evaluate_balanced_constraints(
        combo_first_complete_s=1000.0,
        hybrid_first_complete_s=650.0,  # 35% faster
        combo_makespan_s=10000.0,
        hybrid_makespan_s=11200.0,  # 12% slower
    )
    assert status_ok.constraint_a_ok is True
    assert status_ok.constraint_b_ok is True

    status_bad = evaluate_balanced_constraints(
        combo_first_complete_s=1000.0,
        hybrid_first_complete_s=800.0,  # 20% faster
        combo_makespan_s=10000.0,
        hybrid_makespan_s=12500.0,  # 25% slower
    )
    assert status_bad.constraint_a_ok is False
    assert status_bad.constraint_b_ok is False


def test_auto_medium_dense_tasks_prefers_hybrid() -> None:
    cfg = {}
    decision = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=True,
        task_count=16,
        pdb_count=24,
        combo_disabled=False,
        requested_mode_raw="auto",
    )
    assert decision.effective_mode == "distributed_combo_hybrid"
    assert decision.reason == "auto_medium_dense_tasks"


def test_auto_small_dense_tasks_prefers_hybrid() -> None:
    cfg = {}
    decision = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=True,
        task_count=4,
        pdb_count=3,
        combo_disabled=False,
        requested_mode_raw="auto",
    )
    assert decision.effective_mode == "distributed_combo_hybrid"
    assert decision.reason == "auto_small_dense_tasks"


def test_auto_small_full_node_coverage_prefers_hybrid() -> None:
    cfg = {}
    decision = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=True,
        task_count=3,
        pdb_count=3,
        combo_disabled=False,
        requested_mode_raw="auto",
    )
    assert decision.effective_mode == "distributed_combo_hybrid"
    assert decision.reason == "auto_small_dense_tasks"
