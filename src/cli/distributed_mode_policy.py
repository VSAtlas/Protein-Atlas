from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Optional


_VALID_MODES = {
    "auto",
    "single_process",
    "distributed_per_protein",
    "distributed_combo",
    "distributed_combo_hybrid",
}


def _normalize_mode(raw: Any) -> str:
    token = str(raw or "").strip().lower().replace("-", "_")
    if token in {"", "default"}:
        return "auto"
    aliases = {
        "single": "single_process",
        "singleprocess": "single_process",
        "per_protein": "distributed_per_protein",
        "distributed_perprotein": "distributed_per_protein",
        "combo": "distributed_combo",
        "hybrid": "distributed_combo_hybrid",
    }
    token = aliases.get(token, token)
    if token not in _VALID_MODES:
        return "auto"
    return token


@dataclass(frozen=True)
class ExecutionModeDecision:
    requested_mode: str
    effective_mode: str
    reason: str
    pdb_count: int
    imbalance_ratio: float
    completion_lane_fraction: float
    rebalance_sec: int
    rebalance_chunks: int

    @property
    def combo_enabled(self) -> bool:
        return self.effective_mode in {
            "distributed_combo",
            "distributed_combo_hybrid",
        }

    @property
    def hybrid_enabled(self) -> bool:
        return self.effective_mode == "distributed_combo_hybrid"


@dataclass(frozen=True)
class BalancedConstraintStatus:
    first_complete_improvement: float
    makespan_degradation: float
    constraint_a_ok: bool
    constraint_b_ok: bool


def _estimate_imbalance_ratio(
    *,
    pdb_count: int,
    work_hint: Optional[Mapping[str, float]],
) -> float:
    if not work_hint:
        return 1.0
    vals = [
        float(v)
        for v in work_hint.values()
        if isinstance(v, (int, float)) and float(v) > 0.0
    ]
    if len(vals) < 4:
        return 1.0
    vals.sort()
    median = vals[len(vals) // 2]
    p90 = vals[int(max(0, min(len(vals) - 1, math.floor(0.9 * (len(vals) - 1)))))]
    if median <= 0:
        return 1.0
    return float(max(1.0, p90 / median))


def _base_lane_fraction(*, pdb_count: int, imbalance_ratio: float) -> float:
    frac = 0.30
    if pdb_count >= 120:
        frac += 0.05
    if imbalance_ratio >= 3.0:
        frac += 0.05
    if imbalance_ratio >= 5.0:
        frac += 0.05
    return float(max(0.20, min(0.50, frac)))


def resolve_execution_mode(
    *,
    cfg: Mapping[str, Any],
    distributed_enabled: bool,
    task_count: int,
    pdb_count: int,
    combo_disabled: bool,
    requested_mode_raw: Any = None,
    work_hint: Optional[Mapping[str, float]] = None,
) -> ExecutionModeDecision:
    requested = _normalize_mode(requested_mode_raw)
    if requested_mode_raw is None:
        requested = _normalize_mode(cfg.get("EXECUTION_MODE", "auto"))
    imbalance = _estimate_imbalance_ratio(pdb_count=pdb_count, work_hint=work_hint)
    lane_fraction = _base_lane_fraction(
        pdb_count=max(0, int(pdb_count)),
        imbalance_ratio=float(imbalance),
    )

    if not distributed_enabled or int(task_count) <= 1:
        return ExecutionModeDecision(
            requested_mode=requested,
            effective_mode="single_process",
            reason="distributed_disabled_or_single_task",
            pdb_count=int(max(0, pdb_count)),
            imbalance_ratio=float(imbalance),
            completion_lane_fraction=float(lane_fraction),
            rebalance_sec=300,
            rebalance_chunks=64,
        )

    if combo_disabled:
        return ExecutionModeDecision(
            requested_mode=requested,
            effective_mode="distributed_per_protein",
            reason="combo_disabled_env",
            pdb_count=int(max(0, pdb_count)),
            imbalance_ratio=float(imbalance),
            completion_lane_fraction=float(lane_fraction),
            rebalance_sec=300,
            rebalance_chunks=64,
        )

    if requested in {
        "single_process",
        "distributed_per_protein",
        "distributed_combo",
        "distributed_combo_hybrid",
    }:
        return ExecutionModeDecision(
            requested_mode=requested,
            effective_mode=requested,
            reason="requested_explicit_mode",
            pdb_count=int(max(0, pdb_count)),
            imbalance_ratio=float(imbalance),
            completion_lane_fraction=float(lane_fraction),
            rebalance_sec=300,
            rebalance_chunks=64,
        )

    # Auto policy defaults to balanced objective.
    p = int(max(0, pdb_count))
    t = int(max(1, task_count))
    pdbs_per_task = float(p) / float(max(1, t))
    if p <= 12:
        if p <= 6 and t >= p and t > 1:
            effective = "distributed_combo_hybrid"
            reason = "auto_small_dense_tasks"
        else:
            effective = "distributed_per_protein"
            reason = "auto_small_run"
    elif p <= 60:
        if float(imbalance) >= 3.0:
            effective = "distributed_combo_hybrid"
            reason = "auto_medium_imbalanced"
        elif pdbs_per_task <= 2.0:
            effective = "distributed_combo_hybrid"
            reason = "auto_medium_dense_tasks"
        else:
            effective = "distributed_per_protein"
            reason = "auto_medium_balanced"
    else:
        effective = "distributed_combo_hybrid"
        reason = "auto_large_run"

    return ExecutionModeDecision(
        requested_mode=requested,
        effective_mode=effective,
        reason=reason,
        pdb_count=int(max(0, pdb_count)),
        imbalance_ratio=float(imbalance),
        completion_lane_fraction=float(lane_fraction),
        rebalance_sec=300,
        rebalance_chunks=64,
    )


def evaluate_balanced_constraints(
    *,
    combo_first_complete_s: float,
    hybrid_first_complete_s: float,
    combo_makespan_s: float,
    hybrid_makespan_s: float,
) -> BalancedConstraintStatus:
    combo_first = max(1e-9, float(combo_first_complete_s))
    hybrid_first = max(1e-9, float(hybrid_first_complete_s))
    combo_makespan = max(1e-9, float(combo_makespan_s))
    hybrid_makespan = max(1e-9, float(hybrid_makespan_s))
    first_complete_improvement = (combo_first - hybrid_first) / combo_first
    makespan_degradation = (hybrid_makespan - combo_makespan) / combo_makespan
    return BalancedConstraintStatus(
        first_complete_improvement=float(first_complete_improvement),
        makespan_degradation=float(makespan_degradation),
        constraint_a_ok=bool(first_complete_improvement >= 0.30),
        constraint_b_ok=bool(makespan_degradation <= 0.15),
    )
