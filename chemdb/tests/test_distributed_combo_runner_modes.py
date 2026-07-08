from __future__ import annotations

from cli.distributed_combo_runner import resolve_combo_runtime_mode


def test_combo_hybrid_mode_with_valid_progress() -> None:
    mode = resolve_combo_runtime_mode(
        requested_mode="distributed_combo_hybrid",
        scope_totals={
            ("P001", "HOLO", "base"): 4,
            ("P002", "HOLO", "base"): 2,
        },
    )
    assert mode.effective_mode == "distributed_combo_hybrid"
    assert mode.reason == "requested_hybrid_with_valid_scope_progress"


def test_combo_hybrid_falls_back_when_progress_missing() -> None:
    missing = resolve_combo_runtime_mode(
        requested_mode="distributed_combo_hybrid",
        scope_totals=None,
    )
    assert missing.effective_mode == "distributed_combo"
    assert missing.reason == "fallback_missing_scope_progress"

    corrupt = resolve_combo_runtime_mode(
        requested_mode="distributed_combo_hybrid",
        scope_totals={("P001", "HOLO", "base"): 0},
    )
    assert corrupt.effective_mode == "distributed_combo"
    assert corrupt.reason == "fallback_invalid_scope_progress"
