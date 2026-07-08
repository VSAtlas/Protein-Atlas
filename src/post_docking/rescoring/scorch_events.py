from __future__ import annotations

from typing import Any, Dict, Optional

from docking.global_scheduler import fit_class_for_cores


def emit_bench_event(event: str, **fields: Any) -> None:
    try:
        from chemdb.bench.util_bench.event_log import emit_event

        emit_event(event, **fields)
    except Exception:
        return


def emit_task_event(
    event: str,
    *,
    task_id: str,
    task_type: str,
    want_cores: int = 1,
    min_cores: int = 1,
    granted_cores: Optional[int] = None,
    est_duration_sec: Optional[float] = None,
    **fields: Any,
) -> None:
    want = max(1, int(want_cores))
    floor = max(1, min(int(min_cores), want))
    payload: Dict[str, Any] = {
        "task_id": task_id,
        "task_type": task_type,
        "want_cores": want,
        "min_cores": floor,
        "fit_class": fit_class_for_cores(want),
        "scheduler": "scorch",
    }
    if est_duration_sec is not None:
        payload["est_duration_sec"] = max(0.001, float(est_duration_sec))
    if granted_cores is not None:
        grant = max(1, int(granted_cores))
        payload["granted_cores"] = grant
        payload["cores"] = grant
    else:
        payload["cores"] = want
    payload.update(fields)
    emit_bench_event(event, **payload)
