from __future__ import annotations

import threading
from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Optional

from cli.distributed_chunk_planner import (
    _CHUNK_CLAIM_HEARTBEAT_SEC,
    _CHUNK_CLAIM_LEASE_SEC,
    _CHUNK_IDLE_TIMEOUT_SEC,
    _CHUNK_IDLE_WAIT_SEC,
    _CHUNK_MAX_ATTEMPTS,
    _CHUNK_STATE_RECONCILE_SEC,
    _COMBO_PREP_LEASE_SEC,
    _COMBO_PREP_WAIT_POLL_SEC,
    _COMBO_PREP_WAIT_TIMEOUT_SEC,
)
from cli.distributed_chunk_runtime_rebalance import HybridRebalanceState
from cli.distributed_chunk_runtime_state import ChunkStateTracker
from cli.distributed_runtime_completion import (
    ComboScope,
    DeferredScorchCompletionBookkeeper,
)
from cli.run_context import ConfigDict


@dataclass
class ClaimLoopConfig:
    chunk_claim_heartbeat_sec: float = float(_CHUNK_CLAIM_HEARTBEAT_SEC)
    chunk_claim_lease_sec: float = float(_CHUNK_CLAIM_LEASE_SEC)
    chunk_idle_timeout_sec: float = float(_CHUNK_IDLE_TIMEOUT_SEC)
    chunk_idle_wait_sec: float = float(_CHUNK_IDLE_WAIT_SEC)
    chunk_max_attempts: int = int(_CHUNK_MAX_ATTEMPTS)
    chunk_state_reconcile_sec: float = float(_CHUNK_STATE_RECONCILE_SEC)
    combo_prep_lease_sec: float = float(_COMBO_PREP_LEASE_SEC)
    combo_prep_wait_poll_sec: float = float(_COMBO_PREP_WAIT_POLL_SEC)
    combo_prep_wait_timeout_sec: float = float(_COMBO_PREP_WAIT_TIMEOUT_SEC)


@dataclass
class ClaimLoopContext:
    cfg_v: ConfigDict
    dist_ctx: Any
    run_id: str
    label: str
    variant: Optional[str]
    by_chunk_id: dict[str, dict[str, Any]]
    combo_owner_task_id: Mapping[ComboScope, int]
    combo_total_chunks: Mapping[ComboScope, int]
    combo_hybrid_active: bool
    execution_mode_decision: Any
    local_chunk_workers: int
    cpu_per_chunk_worker: int
    combo_total_weight: float
    baseline_first_complete_sec: Optional[float]
    baseline_makespan_sec: Optional[float]
    process_one: Callable[..., bool]
    distributed_assigned_pdb_ids: set[str]
    recorded_terminal_chunk_ids: set[str]
    failed_entries: list[tuple[str, str, str, str, str]]
    bar: Any


@dataclass
class ClaimLoopResult:
    rebalance_count: int


@dataclass
class ClaimLoopSession:
    context: ClaimLoopContext
    cfg: ClaimLoopConfig
    completion_bookkeeper: DeferredScorchCompletionBookkeeper
    chunk_state_tracker: ChunkStateTracker
    hybrid_state: HybridRebalanceState
    claim_cursor: int
    hybrid_completed_weight: float
    active_worker_count: int
    use_claim_cache: bool
    hybrid_start_ts: float
    hybrid_rebalance_sec: int
    hybrid_rebalance_chunks: int
    claim_cursor_lock: threading.Lock = field(default_factory=threading.Lock)
    hybrid_control_lock: threading.Lock = field(default_factory=threading.Lock)
    progress_lock: threading.Lock = field(default_factory=threading.Lock)
    failure_lock: threading.Lock = field(default_factory=threading.Lock)
    active_worker_lock: threading.Lock = field(default_factory=threading.Lock)
