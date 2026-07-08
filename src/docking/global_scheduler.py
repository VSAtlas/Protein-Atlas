from __future__ import annotations

import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, Iterator

# Keep scheduler tuning internal to avoid many new config lines.
_AGING_STEP_SEC = 8.0
_AGING_MAX_BOOST = 3
_TAIL_MODE_FREE_CORE_FRAC = 0.35
_LARGE_WANT_THRESHOLD = 5
_CRITICAL_TAIL_PENDING_MAX = 6
_CRITICAL_TAIL_FREE_CORE_FRAC = 0.45
_CRITICAL_LARGE_AGE_SEC = 6.0


def fit_class_for_cores(cores: int) -> str:
    c = max(1, int(cores))
    if c == 1:
        return "c1"
    if c == 2:
        return "c2"
    if c <= 4:
        return "c3_4"
    return "c5p"


def _fit_class_index(fit_class: str) -> int:
    order = {"c1": 0, "c2": 1, "c3_4": 2, "c5p": 3}
    return order.get(str(fit_class), 3)


@dataclass
class _PendingRequest:
    req_id: int
    want_cores: int
    min_cores: int
    priority: int
    est_duration_sec: float
    fit_class: str
    cond: threading.Condition
    submitted_ns: int
    task_id: str | None = None
    task_type: str | None = None
    granted: bool = False
    granted_cores: int = 0


class GlobalAdmissionScheduler:
    """
    Core-aware in-process admission scheduler with:
    - best-fit packing
    - elastic min..want core requests
    - aging-based starvation protection
    - pressure-mode tight-fit preference for mixed-size queues
    """

    def __init__(self, total_cores: int) -> None:
        self.total_cores = max(1, int(total_cores))
        self._available = self.total_cores
        self._lock = threading.Lock()
        self._pending: list[_PendingRequest] = []
        self._next_req_id = 1

    def available_cores(self) -> int:
        with self._lock:
            return max(0, int(self._available))

    def queued_requests(self) -> int:
        with self._lock:
            return len(self._pending)

    @staticmethod
    def _age_boost(req: _PendingRequest, now_ns: int) -> int:
        waited_sec = max(0.0, (now_ns - int(req.submitted_ns)) / 1_000_000_000.0)
        boost = int(waited_sec // _AGING_STEP_SEC)
        return max(0, min(_AGING_MAX_BOOST, boost))

    def _under_pressure_locked(self, avail: int) -> bool:
        pending_n = len(self._pending)
        if pending_n > int(self.total_cores) * 2:
            return True
        return avail <= max(1, int(self.total_cores) // 3) and pending_n > int(
            self.total_cores
        )

    def _tail_mode_locked(self, avail: int) -> bool:
        pending_n = len(self._pending)
        if pending_n <= 0:
            return False
        total = max(1, int(self.total_cores))
        # Tail mode is for endgame drain, not general mixed-load operation.
        if pending_n > max(2, total // 3):
            return False
        if avail < max(1, total // 4):
            return False
        free_frac = float(avail) / float(total)
        if pending_n <= max(2, total // 4) and free_frac >= _TAIL_MODE_FREE_CORE_FRAC:
            return True
        fit_now = sum(1 for req in self._pending if int(req.min_cores) <= int(avail))
        low_fit = fit_now <= max(1, pending_n // 2)
        large_pressure = sum(
            1 for req in self._pending if int(req.want_cores) > int(avail)
        )
        mostly_large = large_pressure >= max(1, (pending_n * 2) // 3)
        return low_fit or (free_frac >= _TAIL_MODE_FREE_CORE_FRAC and mostly_large)

    def _critical_tail_mode_locked(self, avail: int) -> bool:
        pending_n = len(self._pending)
        if pending_n <= 0:
            return False
        total = max(1, int(self.total_cores))
        if pending_n > max(_CRITICAL_TAIL_PENDING_MAX, total // 2):
            return False
        free_frac = float(avail) / float(total)
        if free_frac < _CRITICAL_TAIL_FREE_CORE_FRAC:
            return False
        large_pending = [
            req for req in self._pending if int(req.want_cores) >= _LARGE_WANT_THRESHOLD
        ]
        if not large_pending:
            return False
        return any(int(req.min_cores) <= int(avail) for req in large_pending)

    def _effective_want_locked(
        self,
        req: _PendingRequest,
        under_pressure: bool,
        tail_mode: bool,
    ) -> int:
        want = max(1, int(req.want_cores))
        floor = max(1, int(req.min_cores))
        if not under_pressure and not tail_mode:
            return max(floor, want)

        total = int(self.total_cores)
        if total <= 4:
            cap = 2
        elif total <= 8:
            cap = 3
        elif total <= 16:
            cap = 4
        else:
            cap = max(4, total // 4)

        if want >= _LARGE_WANT_THRESHOLD and (under_pressure or tail_mode):
            if total <= 8:
                cap = min(cap, 2)
            elif total <= 16:
                cap = min(cap, 3)
            elif total <= 32:
                cap = min(cap, 4)
            else:
                cap = min(cap, max(4, total // 8))
        return max(floor, min(want, cap))

    def _fit_class_penalty(self, req: _PendingRequest, avail: int) -> int:
        target = _fit_class_index(fit_class_for_cores(avail))
        req_cls = _fit_class_index(req.fit_class)
        return abs(req_cls - target)

    def _pick_candidate_locked(self) -> tuple[int, int] | None:
        if not self._pending:
            return None
        avail = int(self._available)
        if avail <= 0:
            return None

        now_ns = time.monotonic_ns()
        under_pressure = self._under_pressure_locked(avail)
        tail_mode = self._tail_mode_locked(avail)
        critical_tail_mode = self._critical_tail_mode_locked(avail)

        best_idx: int | None = None
        best_grant = 0
        best_key: tuple[float, ...] | None = None

        for idx, req in enumerate(self._pending):
            min_fit = max(1, int(req.min_cores))
            if min_fit > avail:
                continue

            want_eff = max(
                min_fit,
                self._effective_want_locked(req, under_pressure, tail_mode),
            )
            grant = min(want_eff, avail)
            eff_priority = int(req.priority) + self._age_boost(req, now_ns)
            waste = avail - grant
            est = max(0.001, float(req.est_duration_sec))

            if critical_tail_mode:
                waited_sec = max(
                    0.0,
                    (now_ns - int(req.submitted_ns)) / 1_000_000_000.0,
                )
                is_large = int(req.want_cores) >= _LARGE_WANT_THRESHOLD
                # Late-run critical-path mode: prioritize aged large chunks that fit,
                # then remaining large chunks, then fillers.
                if is_large and waited_sec >= _CRITICAL_LARGE_AGE_SEC:
                    critical_rank = 0.0
                elif is_large:
                    critical_rank = 1.0
                else:
                    critical_rank = 2.0
                key = (
                    float(-eff_priority),
                    critical_rank,
                    float(waste),
                    float(est),
                    float(req.req_id),
                )
            elif tail_mode:
                class_penalty = self._fit_class_penalty(req, avail)
                # In tail mode, drain fitting shorter jobs first to reduce stranded end cores.
                key = (
                    float(-eff_priority),
                    float(waste),
                    float(est),
                    float(class_penalty),
                    float(req.req_id),
                )
            elif under_pressure:
                class_penalty = self._fit_class_penalty(req, avail)
                # Fit-first under pressure; duration is a tiebreak.
                key = (
                    float(-eff_priority),
                    float(class_penalty),
                    float(waste),
                    float(est),
                    float(req.req_id),
                )
            else:
                key = (
                    float(-eff_priority),
                    float(waste),
                    float(est),
                    0.0,
                    float(req.req_id),
                )

            if best_key is None or key < best_key:
                best_key = key
                best_idx = idx
                best_grant = grant

        if best_idx is None:
            return None
        return best_idx, best_grant

    def _grant_pending_locked(self) -> None:
        while self._pending:
            picked = self._pick_candidate_locked()
            if picked is None:
                break
            idx, grant = picked
            req = self._pending.pop(idx)
            self._available = max(0, int(self._available) - int(grant))
            req.granted_cores = int(grant)
            req.granted = True
            req.cond.notify()

    def acquire(
        self,
        *,
        cores: int = 1,
        min_cores: int | None = None,
        priority: int = 0,
        est_duration_sec: float | None = None,
        task_id: str | None = None,
        task_type: str | None = None,
        timeout: float | None = None,
    ) -> int:
        want = max(1, min(int(cores), int(self.total_cores)))
        floor = max(1, int(min_cores if min_cores is not None else want))
        floor = min(floor, want)
        est = max(0.001, float(est_duration_sec if est_duration_sec is not None else 1.0))
        deadline = (time.monotonic() + float(timeout)) if timeout is not None else None

        with self._lock:
            req = _PendingRequest(
                req_id=self._next_req_id,
                want_cores=want,
                min_cores=floor,
                priority=int(priority),
                est_duration_sec=est,
                fit_class=fit_class_for_cores(want),
                cond=threading.Condition(self._lock),
                submitted_ns=time.monotonic_ns(),
                task_id=task_id,
                task_type=task_type,
            )
            self._next_req_id += 1
            self._pending.append(req)
            self._grant_pending_locked()

            while not req.granted:
                if deadline is None:
                    req.cond.wait()
                else:
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        self._pending = [p for p in self._pending if p.req_id != req.req_id]
                        raise TimeoutError(
                            f"global admission timeout for task={task_id or req.req_id}"
                        )
                    req.cond.wait(timeout=remaining)
                if not req.granted:
                    self._grant_pending_locked()

            return max(1, int(req.granted_cores))

    def release(self, granted_cores: int = 1) -> None:
        rel = max(1, int(granted_cores))
        with self._lock:
            self._available = min(int(self.total_cores), int(self._available) + rel)
            self._grant_pending_locked()


def get_global_docking_sem(cfg: dict[str, Any] | Any) -> Any:
    if not isinstance(cfg, dict):
        return None
    sem = cfg.get("GLOBAL_DOCKING_SEM")
    if sem is None:
        sem = cfg.get("GLOBAL_LIGAND_SEM")
    if sem is None:
        return None
    if not hasattr(sem, "acquire") or not hasattr(sem, "release"):
        return None
    return sem


def get_global_docking_scheduler(
    cfg: dict[str, Any] | Any,
) -> GlobalAdmissionScheduler | None:
    if not isinstance(cfg, dict):
        return None
    sched = cfg.get("GLOBAL_DOCKING_SCHEDULER")
    if isinstance(sched, GlobalAdmissionScheduler):
        return sched
    return None


def global_scheduler_capacity(cfg: dict[str, Any] | Any, fallback: int) -> int:
    if not isinstance(cfg, dict):
        return max(1, int(fallback))
    raw = cfg.get("GLOBAL_SCHEDULER_CPUS", fallback)
    try:
        cap = int(raw)
    except Exception:
        cap = int(fallback)
    return max(1, cap)


def sem_free_slots(sem: Any, default: int) -> int:
    if hasattr(sem, "available_cores"):
        try:
            return max(0, int(sem.available_cores()))
        except Exception:
            return max(0, int(default))
    try:
        value = int(getattr(sem, "_value"))
    except Exception:
        return max(0, int(default))
    return max(0, value)


@contextmanager
def acquire_global_cores(
    cfg: dict[str, Any] | Any,
    *,
    cores: int = 1,
    min_cores: int | None = None,
    priority: int = 0,
    est_duration_sec: float | None = None,
    task_id: str | None = None,
    task_type: str | None = None,
) -> Iterator[int]:
    scheduler = get_global_docking_scheduler(cfg)
    if scheduler is not None:
        granted = scheduler.acquire(
            cores=max(1, int(cores)),
            min_cores=min_cores,
            priority=priority,
            est_duration_sec=est_duration_sec,
            task_id=task_id,
            task_type=task_type,
        )
        try:
            yield granted
        finally:
            scheduler.release(granted)
        return

    sem = get_global_docking_sem(cfg)
    if sem is None:
        yield max(1, int(cores))
        return
    sem.acquire()
    try:
        yield max(1, int(cores))
    finally:
        sem.release()


@contextmanager
def acquire_global_slot(
    cfg: dict[str, Any] | Any,
    *,
    cores: int = 1,
    min_cores: int | None = None,
    priority: int = 0,
    est_duration_sec: float | None = None,
    task_id: str | None = None,
    task_type: str | None = None,
) -> Iterator[None]:
    with acquire_global_cores(
        cfg,
        cores=cores,
        min_cores=min_cores,
        priority=priority,
        est_duration_sec=est_duration_sec,
        task_id=task_id,
        task_type=task_type,
    ):
        yield
