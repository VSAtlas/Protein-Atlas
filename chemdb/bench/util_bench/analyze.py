from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


_BUCKETS = ("starvation", "dispatch_gap", "fragmentation", "iowait")


@dataclass(frozen=True)
class CpuPoint:
    t_wall_ns: int
    util: float
    idle: float
    iowait: float


@dataclass
class SchedulerState:
    pending: dict[str, "TaskShape"]
    running: dict[str, "TaskShape"]
    idle_workers: set[str]
    task_events_total: int = 0
    task_events_with_min_cores: int = 0
    task_starts_total: int = 0
    task_starts_with_granted_cores: int = 0


@dataclass(frozen=True)
class TaskShape:
    want_cores: int
    min_cores: int
    granted_cores: int
    est_duration_sec: float
    fit_class: str



def _load_events(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(obj, dict):
                rows.append(obj)
    rows.sort(key=lambda x: int(x.get("ts_mono_ns", 0)))
    return rows



def _safe_float(val: Any, default: float = 0.0) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default



def _safe_int(val: Any, default: int = 0) -> int:
    try:
        return int(val)
    except (TypeError, ValueError):
        return default



def _rows_to_cpu_points(rows: list[dict[str, Any]], alloc_cpus: int) -> list[CpuPoint]:
    by_t: dict[int, list[tuple[int, float, float, float]]] = {}
    for row in rows:
        t_ns = _safe_int(row.get("t_wall_ns"), 0)
        if t_ns <= 0:
            continue
        cpu_idx = _safe_int(row.get("cpu_idx"), 0)
        util = max(0.0, min(1.0, _safe_float(row.get("utilization"), 0.0)))
        idle = max(0.0, min(1.0, _safe_float(row.get("idle"), 0.0)))
        iowait = max(0.0, min(1.0, _safe_float(row.get("iowait"), 0.0)))
        by_t.setdefault(t_ns, []).append((cpu_idx, util, idle, iowait))

    points: list[CpuPoint] = []
    for t_ns in sorted(by_t.keys()):
        rows_t = sorted(by_t[t_ns], key=lambda r: r[0])
        rows_t = rows_t[: max(1, min(alloc_cpus, len(rows_t)))]
        n = float(len(rows_t))
        util = sum(r[1] for r in rows_t) / n
        idle = sum(r[2] for r in rows_t) / n
        iowait = sum(r[3] for r in rows_t) / n
        points.append(CpuPoint(t_wall_ns=t_ns, util=util, idle=idle, iowait=iowait))
    return points



def _load_cpu_csv_points(path: Path, alloc_cpus: int) -> list[CpuPoint]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return _rows_to_cpu_points(rows, alloc_cpus)



def _parse_mpstat_line(line: str) -> tuple[str, str, list[float]] | None:
    parts = line.split()
    if len(parts) < 12:
        return None

    cpu_idx = -1
    for idx, tok in enumerate(parts):
        if tok == "all" or tok.isdigit():
            if len(parts) - idx - 1 >= 10:
                cpu_idx = idx
                break
    if cpu_idx < 0:
        return None

    time_key = " ".join(parts[:cpu_idx]).strip()
    cpu = parts[cpu_idx]
    nums = parts[cpu_idx + 1 : cpu_idx + 11]
    vals: list[float] = []
    for tok in nums:
        try:
            vals.append(float(tok.replace(",", ".")))
        except ValueError:
            return None
    return time_key, cpu, vals



def _load_mpstat_points(
    raw_path: Path,
    alloc_cpus: int,
    run_start_wall_ns: int,
    interval_sec: float,
) -> list[CpuPoint]:
    if not raw_path.exists():
        return []

    rows: list[dict[str, Any]] = []
    tick_idx = 0
    time_map: dict[str, int] = {}

    with raw_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parsed = _parse_mpstat_line(line)
            if not parsed:
                continue
            time_key, cpu_tok, vals = parsed
            if cpu_tok == "all":
                continue
            if not cpu_tok.isdigit():
                continue

            if time_key not in time_map:
                time_map[time_key] = tick_idx
                tick_idx += 1
            sample_idx = time_map[time_key]
            t_wall_ns = int(run_start_wall_ns + sample_idx * interval_sec * 1_000_000_000)

            usr, nice, sys_v, iowait, irq, soft, steal, _guest, _gnice, idle = vals
            user = (usr + nice) / 100.0
            system = (sys_v + irq + soft + steal) / 100.0
            idle_f = idle / 100.0
            iowait_f = iowait / 100.0
            util_f = max(0.0, min(1.0, 1.0 - idle_f))
            rows.append(
                {
                    "t_wall_ns": t_wall_ns,
                    "cpu_idx": int(cpu_tok),
                    "user": user,
                    "system": system,
                    "idle": idle_f,
                    "iowait": iowait_f,
                    "utilization": util_f,
                }
            )

    return _rows_to_cpu_points(rows, alloc_cpus)



def _write_normalized_cpu_log(path: Path, points: list[CpuPoint], alloc_cpus: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["t_wall_ns", "cpu_idx", "user", "system", "idle", "iowait", "utilization"])
        for point in points:
            for cpu_idx in range(max(1, alloc_cpus)):
                writer.writerow(
                    [
                        point.t_wall_ns,
                        cpu_idx,
                        "",
                        "",
                        f"{point.idle:.8f}",
                        f"{point.iowait:.8f}",
                        f"{point.util:.8f}",
                    ]
                )



def _task_shape_from_event(event: dict[str, Any], fallback: TaskShape | None = None) -> TaskShape:
    want = _safe_int(event.get("want_cores"), 0)
    if want <= 0:
        want = _safe_int(event.get("cores"), 0)
    if want <= 0 and fallback is not None:
        want = int(fallback.want_cores)
    want = max(1, want)

    min_cores = _safe_int(event.get("min_cores"), 0)
    if min_cores <= 0 and fallback is not None:
        min_cores = int(fallback.min_cores)
    if min_cores <= 0:
        min_cores = want
    min_cores = max(1, min(min_cores, want))

    granted = _safe_int(event.get("granted_cores"), 0)
    if granted <= 0 and fallback is not None:
        granted = int(fallback.granted_cores)
    if granted <= 0:
        granted = _safe_int(event.get("cores"), want)
    granted = max(1, granted)

    est = _safe_float(event.get("est_duration_sec"), 0.0)
    if est <= 0.0 and fallback is not None:
        est = float(fallback.est_duration_sec)
    if est <= 0.0:
        est = 1.0

    fit = str(event.get("fit_class", "")).strip()
    if not fit and fallback is not None:
        fit = str(fallback.fit_class)

    return TaskShape(
        want_cores=want,
        min_cores=min_cores,
        granted_cores=granted,
        est_duration_sec=max(0.001, est),
        fit_class=fit,
    )


def _has_numeric_field(event: dict[str, Any], field: str) -> bool:
    if field not in event:
        return False
    try:
        return float(event[field]) > 0
    except (TypeError, ValueError):
        return False


def _hedge_metrics(events: list[dict[str, Any]]) -> dict[str, float | int]:
    hedge_started = 0
    hedge_canceled = 0
    winner_role_by_group: dict[str, str] = {}
    active_hedge: dict[str, tuple[int, int]] = {}
    hedge_extra_core_seconds = 0.0

    for ev in events:
        et = str(ev.get("event", ""))
        task_id = str(ev.get("task_id", "")).strip()
        group_id = str(ev.get("hedge_group_id", "")).strip()
        role = str(ev.get("hedge_role", "primary")).strip().lower() or "primary"
        ts_mono_ns = _safe_int(ev.get("ts_mono_ns"), 0)

        if et == "TASK_START" and role == "hedge":
            hedge_started += 1
            cores = max(1, _safe_int(ev.get("granted_cores"), _safe_int(ev.get("cores"), 1)))
            if task_id and ts_mono_ns > 0:
                active_hedge[task_id] = (ts_mono_ns, cores)
        elif et == "TASK_END":
            status = str(ev.get("status", "")).strip().lower()
            if role == "hedge" and task_id and task_id in active_hedge and ts_mono_ns > 0:
                start_ns, cores = active_hedge.pop(task_id)
                dt_sec = max(0.0, (ts_mono_ns - start_ns) / 1_000_000_000.0)
                hedge_extra_core_seconds += dt_sec * float(max(1, cores))
            if status == "ok" and group_id and group_id not in winner_role_by_group:
                winner_role_by_group[group_id] = role
        elif et == "HEDGE_CANCELLED":
            hedge_canceled += 1

    hedge_won = sum(1 for role in winner_role_by_group.values() if role == "hedge")
    return {
        "hedge_started": int(hedge_started),
        "hedge_won": int(hedge_won),
        "hedge_canceled": int(hedge_canceled),
        "hedge_extra_core_seconds_est": float(hedge_extra_core_seconds),
    }


def _apply_event(state: SchedulerState, event: dict[str, Any]) -> None:
    et = str(event.get("event", ""))
    task_id = str(event.get("task_id", "")).strip()
    if et in {"TASK_SUBMITTED", "TASK_READY", "TASK_START", "TASK_END"}:
        state.task_events_total += 1
        if _has_numeric_field(event, "min_cores"):
            state.task_events_with_min_cores += 1
    if et == "TASK_START":
        state.task_starts_total += 1
        if _has_numeric_field(event, "granted_cores"):
            state.task_starts_with_granted_cores += 1

    if et in {"TASK_SUBMITTED", "TASK_READY"}:
        if task_id:
            state.pending[task_id] = _task_shape_from_event(
                event,
                fallback=state.pending.get(task_id),
            )
    elif et == "TASK_START":
        if task_id:
            pending_shape = state.pending.pop(task_id, None)
            shape = _task_shape_from_event(event, fallback=pending_shape)
            state.running[task_id] = shape
    elif et == "TASK_END":
        if task_id:
            state.running.pop(task_id, None)
            state.pending.pop(task_id, None)
    elif et == "WORKER_IDLE_START":
        worker_id = str(event.get("worker_id", "")).strip()
        if worker_id:
            state.idle_workers.add(worker_id)
    elif et == "WORKER_IDLE_END":
        worker_id = str(event.get("worker_id", "")).strip()
        if worker_id:
            state.idle_workers.discard(worker_id)



def _classify_idle_bucket(
    *,
    alloc_cpus: int,
    idle_frac: float,
    iowait_frac: float,
    state: SchedulerState,
) -> str:
    runnable = bool(state.pending)
    running_cores = sum(int(shape.granted_cores) for shape in state.running.values())
    free_cores = max(0, int(alloc_cpus) - int(running_cores))
    idle_workers = len(state.idle_workers) if state.idle_workers else free_cores
    pending_shapes = list(state.pending.values())
    min_runnable = (
        min(int(shape.min_cores) for shape in pending_shapes) if pending_shapes else None
    )
    min_want = min(int(shape.want_cores) for shape in pending_shapes) if pending_shapes else None
    want_fit_exists = (
        any(int(shape.want_cores) <= free_cores for shape in pending_shapes)
        if pending_shapes
        else False
    )
    floor_fit_exists = (
        any(int(shape.min_cores) <= free_cores for shape in pending_shapes)
        if pending_shapes
        else False
    )
    any_elastic_pending = (
        any(int(shape.min_cores) < int(shape.want_cores) for shape in pending_shapes)
        if pending_shapes
        else False
    )
    unmet_running_wants = any(
        int(shape.want_cores) > int(shape.granted_cores)
        for shape in state.running.values()
    )

    iowait_dominates = iowait_frac >= max(0.05, idle_frac * 0.5)
    if iowait_dominates:
        return "iowait"
    if not runnable:
        if free_cores > 0 and unmet_running_wants:
            return "fragmentation"
        return "starvation"
    if min_runnable is not None and free_cores < min_runnable:
        return "fragmentation"
    if (
        floor_fit_exists
        and not want_fit_exists
        and min_want is not None
        and free_cores < min_want
        and (any_elastic_pending or unmet_running_wants)
    ):
        return "fragmentation"
    if idle_workers > 0 and free_cores > 0:
        return "dispatch_gap"
    return "starvation"



def analyze_run(run_dir: str | Path, alloc_cpus: int, interval_sec: float) -> dict[str, Any]:
    run_path = Path(run_dir)
    events_path = run_path / "events.jsonl"
    cpu_csv_path = run_path / "cpu.log"
    cpu_mpstat_path = run_path / "cpu_mpstat.log"

    events = _load_events(events_path)

    run_start_event = next((e for e in events if e.get("event") == "RUN_START"), None)
    run_end_event = next((e for e in reversed(events) if e.get("event") == "RUN_END"), None)

    if run_start_event and run_end_event:
        run_start_mono_ns = _safe_int(run_start_event.get("ts_mono_ns"))
        run_end_mono_ns = _safe_int(run_end_event.get("ts_mono_ns"))
        run_start_wall_ns = _safe_int(run_start_event.get("ts_wall_ns"))
    elif events:
        run_start_mono_ns = _safe_int(events[0].get("ts_mono_ns"))
        run_end_mono_ns = _safe_int(events[-1].get("ts_mono_ns"))
        run_start_wall_ns = _safe_int(events[0].get("ts_wall_ns"))
    else:
        now = time_ns_fallback()
        run_start_mono_ns = now
        run_end_mono_ns = now
        run_start_wall_ns = now

    cpu_points = _load_cpu_csv_points(cpu_csv_path, alloc_cpus=alloc_cpus)
    if not cpu_points and cpu_mpstat_path.exists():
        cpu_points = _load_mpstat_points(
            cpu_mpstat_path,
            alloc_cpus=alloc_cpus,
            run_start_wall_ns=run_start_wall_ns,
            interval_sec=interval_sec,
        )
        if cpu_points:
            _write_normalized_cpu_log(cpu_csv_path, cpu_points, alloc_cpus=alloc_cpus)

    if len(cpu_points) < 2:
        wall_sec = max(0.0, (run_end_mono_ns - run_start_mono_ns) / 1_000_000_000.0)
        summary = {
            "AllocCPUs": int(alloc_cpus),
            "WallSeconds": wall_sec,
            "BusyCoreSeconds": 0.0,
            "IdleCoreSeconds": float(alloc_cpus) * wall_sec,
            "IdleFrac": 1.0 if wall_sec > 0 else 0.0,
            "iowait_share": 0.0,
            "attribution": {
                k: {"idle_core_seconds": 0.0, "idle_fraction": 0.0}
                for k in _BUCKETS
            },
            "event_fields": {
                "task_events_total": 0,
                "task_events_with_min_cores": 0,
                "task_starts_total": 0,
                "task_starts_with_granted_cores": 0,
            },
        }
        _write_summaries(run_path, summary)
        return summary

    mono_wall_offset = run_start_mono_ns - run_start_wall_ns

    event_idx = 0
    state = SchedulerState(pending={}, running={}, idle_workers=set())

    busy_core_sec = 0.0
    idle_core_sec = 0.0
    iowait_core_sec = 0.0

    bucket_core_seconds = {k: 0.0 for k in _BUCKETS}

    for i in range(1, len(cpu_points)):
        p0 = cpu_points[i - 1]
        p1 = cpu_points[i]
        dt = max(0.0, (p1.t_wall_ns - p0.t_wall_ns) / 1_000_000_000.0)
        if dt <= 0:
            continue

        t0_mono = p0.t_wall_ns + mono_wall_offset
        while event_idx < len(events) and _safe_int(events[event_idx].get("ts_mono_ns")) <= t0_mono:
            _apply_event(state, events[event_idx])
            event_idx += 1

        util = max(0.0, min(1.0, p0.util))
        idle = max(0.0, min(1.0, p0.idle))
        iowait = max(0.0, min(1.0, p0.iowait))

        busy_cs = util * alloc_cpus * dt
        idle_cs = idle * alloc_cpus * dt
        iowait_cs = iowait * alloc_cpus * dt

        busy_core_sec += busy_cs
        idle_core_sec += idle_cs
        iowait_core_sec += iowait_cs

        if idle_cs > 0:
            bucket = _classify_idle_bucket(
                alloc_cpus=alloc_cpus,
                idle_frac=idle,
                iowait_frac=iowait,
                state=state,
            )
            bucket_core_seconds[bucket] += idle_cs

    wall_sec = max(0.0, (run_end_mono_ns - run_start_mono_ns) / 1_000_000_000.0)
    denom = max(1e-12, alloc_cpus * wall_sec)
    idle_frac = max(0.0, min(1.0, 1.0 - (busy_core_sec / denom)))
    iowait_share = max(0.0, min(1.0, iowait_core_sec / denom))

    attribution = {}
    idle_denom = max(1e-12, idle_core_sec)
    for bucket in _BUCKETS:
        cs = bucket_core_seconds.get(bucket, 0.0)
        attribution[bucket] = {
            "idle_core_seconds": cs,
            "idle_seconds": cs / max(1, alloc_cpus),
            "idle_fraction": cs / idle_denom,
        }

    summary = {
        "AllocCPUs": int(alloc_cpus),
        "WallSeconds": wall_sec,
        "BusyCoreSeconds": busy_core_sec,
        "IdleCoreSeconds": idle_core_sec,
        "IdleFrac": idle_frac,
        "iowait_share": iowait_share,
        "attribution": attribution,
        "event_fields": {
            "task_events_total": int(state.task_events_total),
            "task_events_with_min_cores": int(state.task_events_with_min_cores),
            "task_starts_total": int(state.task_starts_total),
            "task_starts_with_granted_cores": int(state.task_starts_with_granted_cores),
        },
        "hedge": _hedge_metrics(events),
    }
    _write_summaries(run_path, summary)
    return summary



def _write_summaries(run_path: Path, summary: dict[str, Any]) -> None:
    (run_path / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    header = [
        "AllocCPUs",
        "WallSeconds",
        "BusyCoreSeconds",
        "IdleCoreSeconds",
        "IdleFrac",
        "iowait_share",
        "task_events_total",
        "task_events_with_min_cores",
        "task_starts_total",
        "task_starts_with_granted_cores",
        "hedge_started",
        "hedge_won",
        "hedge_canceled",
        "hedge_extra_core_seconds_est",
    ]
    for bucket in _BUCKETS:
        header.append(f"idle_core_seconds_{bucket}")
        header.append(f"idle_fraction_{bucket}")

    row: dict[str, Any] = {
        "AllocCPUs": summary.get("AllocCPUs"),
        "WallSeconds": summary.get("WallSeconds"),
        "BusyCoreSeconds": summary.get("BusyCoreSeconds"),
        "IdleCoreSeconds": summary.get("IdleCoreSeconds"),
        "IdleFrac": summary.get("IdleFrac"),
        "iowait_share": summary.get("iowait_share"),
        "task_events_total": summary.get("event_fields", {}).get("task_events_total", 0),
        "task_events_with_min_cores": summary.get("event_fields", {}).get(
            "task_events_with_min_cores", 0
        ),
        "task_starts_total": summary.get("event_fields", {}).get("task_starts_total", 0),
        "task_starts_with_granted_cores": summary.get("event_fields", {}).get(
            "task_starts_with_granted_cores", 0
        ),
        "hedge_started": summary.get("hedge", {}).get("hedge_started", 0),
        "hedge_won": summary.get("hedge", {}).get("hedge_won", 0),
        "hedge_canceled": summary.get("hedge", {}).get("hedge_canceled", 0),
        "hedge_extra_core_seconds_est": summary.get("hedge", {}).get(
            "hedge_extra_core_seconds_est", 0.0
        ),
    }
    for bucket in _BUCKETS:
        b = summary.get("attribution", {}).get(bucket, {})
        row[f"idle_core_seconds_{bucket}"] = b.get("idle_core_seconds", 0.0)
        row[f"idle_fraction_{bucket}"] = b.get("idle_fraction", 0.0)

    with (run_path / "summary.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header)
        writer.writeheader()
        writer.writerow(row)



def time_ns_fallback() -> int:
    try:
        import time

        return time.time_ns()
    except Exception:
        return 0
