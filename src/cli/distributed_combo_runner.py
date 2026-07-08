from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class ComboRuntimeMode:
    requested_mode: str
    effective_mode: str
    reason: str

    @property
    def hybrid_enabled(self) -> bool:
        return self.effective_mode == "distributed_combo_hybrid"


def resolve_combo_runtime_mode(
    *,
    requested_mode: str,
    scope_totals: Mapping[Any, Any] | None,
) -> ComboRuntimeMode:
    req = str(requested_mode or "").strip().lower()
    if req not in {"distributed_combo", "distributed_combo_hybrid"}:
        req = "distributed_combo"
    if req == "distributed_combo":
        return ComboRuntimeMode(
            requested_mode=req,
            effective_mode="distributed_combo",
            reason="requested_combo",
        )

    # Hybrid requires valid scope progress bookkeeping. Fallback to combo if missing.
    if not isinstance(scope_totals, Mapping) or not scope_totals:
        return ComboRuntimeMode(
            requested_mode=req,
            effective_mode="distributed_combo",
            reason="fallback_missing_scope_progress",
        )
    valid_totals = 0
    for value in scope_totals.values():
        try:
            if int(value) > 0:
                valid_totals += 1
        except Exception:
            continue
    if valid_totals <= 0:
        return ComboRuntimeMode(
            requested_mode=req,
            effective_mode="distributed_combo",
            reason="fallback_invalid_scope_progress",
        )
    return ComboRuntimeMode(
        requested_mode=req,
        effective_mode="distributed_combo_hybrid",
        reason="requested_hybrid_with_valid_scope_progress",
    )
