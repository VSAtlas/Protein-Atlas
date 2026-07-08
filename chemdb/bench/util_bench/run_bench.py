from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import socket
import subprocess
import time
from pathlib import Path
from typing import Any

from .analyze import analyze_run
from .cpu_sampler import start_sampler
from .event_log import emit_event, event_run
from .workloads import load_workload


def _parse_cpu_env(raw: str | None) -> int | None:
    if not raw:
        return None
    m = re.search(r"(\d+)", str(raw))
    if not m:
        return None
    try:
        val = int(m.group(1))
    except ValueError:
        return None
    return val if val > 0 else None



def _detect_alloc_cpus(cli_alloc: int | None) -> int:
    if cli_alloc and cli_alloc > 0:
        return int(cli_alloc)

    for key in ("SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE", "SLURM_CPUS_PER_TASK"):
        parsed = _parse_cpu_env(os.environ.get(key))
        if parsed:
            return parsed

    return int(os.cpu_count() or 1)



def _git_sha(repo_root: Path) -> str | None:
    git = shutil.which("git")
    if not git:
        return None
    try:
        proc = subprocess.run(
            [git, "rev-parse", "--short", "HEAD"],
            cwd=str(repo_root),
            check=False,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None
    if proc.returncode != 0:
        return None
    out = (proc.stdout or "").strip()
    return out or None



def _now_run_id(prefix: str = "utilbench") -> str:
    ts = time.strftime("%Y%m%d-%H%M%S")
    return f"{prefix}_{ts}_{os.getpid()}"



def _serialize_task(task: dict[str, Any]) -> dict[str, Any]:
    cores = max(1, int(task.get("cores", 1)))
    expected = task.get("expected_sec")
    est: float | None = None
    if expected is not None:
        try:
            est = float(expected)
        except (TypeError, ValueError):
            est = None
    return {
        "id": str(task.get("id", "")),
        "cores": cores,
        "want_cores": cores,
        "min_cores": 1 if cores >= 5 else cores,
        "elastic": bool(cores >= 5),
        "command": str(task.get("command", "")),
        "expected_sec": est,
        "tags": task.get("tags") or [],
    }


def _fit_class_for_cores(cores: int) -> str:
    c = max(1, int(cores))
    if c == 1:
        return "c1"
    if c == 2:
        return "c2"
    if c <= 4:
        return "c3_4"
    return "c5p"


def _command_with_granted_workers(command: str, workers: int) -> str:
    cmd = str(command or "").strip()
    if not cmd:
        return cmd
    # Stable envelope for elastic tasks; commands may optionally consume this env var.
    return f"UTIL_BENCH_GRANTED_WORKERS={max(1, int(workers))} {cmd}"



def _pick_dispatch_index(
    pending: list[dict[str, Any]],
    idle_workers: int,
    dispatch_policy: str,
) -> int:
    if not pending:
        return -1
    if dispatch_policy == "strict_fifo":
        head = pending[0]
        want = max(1, int(head.get("want_cores", head.get("cores", 1))))
        floor = max(1, int(head.get("min_cores", want)))
        elastic = bool(head.get("elastic", False))
        if want <= idle_workers:
            return 0
        if elastic and floor <= idle_workers:
            return 0
        return -1
    best_idx = -1
    best_key: tuple[float, float, int, int] | None = None
    for idx, task in enumerate(pending):
        want = max(1, int(task.get("want_cores", task.get("cores", 1))))
        floor = max(1, int(task.get("min_cores", want)))
        elastic = bool(task.get("elastic", False))
        if floor > idle_workers:
            continue
        grant = min(idle_workers, want) if elastic else want
        if grant > idle_workers or grant < floor:
            continue
        waste = float(max(0, idle_workers - grant))
        est = float(task.get("expected_sec") or 1.0)
        key = (waste, est, grant, idx)
        if best_key is None or key < best_key:
            best_key = key
            best_idx = idx
    return best_idx


def _run_synthetic_scheduler(
    tasks: list[dict[str, Any]],
    alloc_cpus: int,
    snapshot_sec: float = 0.2,
    dispatch_policy: str = "backfill",
) -> None:
    idle_workers: list[str] = [f"w{i}" for i in range(max(1, alloc_cpus))]
    running: list[dict[str, Any]] = []
    pending = [_serialize_task(t) for t in tasks]

    for task in pending:
        want = max(1, int(task.get("want_cores", task["cores"])))
        floor = max(1, int(task.get("min_cores", want)))
        est = float(task.get("expected_sec") or 1.0)
        common = {
            "task_id": task["id"],
            "cores": want,
            "want_cores": want,
            "min_cores": floor,
            "fit_class": _fit_class_for_cores(want),
            "est_duration_sec": max(0.001, est),
            "task_type": "synthetic",
        }
        emit_event("TASK_SUBMITTED", **common)
        emit_event("TASK_READY", **common)

    for wid in idle_workers:
        emit_event("WORKER_IDLE_START", worker_id=wid)

    last_snapshot = time.monotonic()

    while pending or running:
        dispatched = False

        while pending:
            fit_idx = _pick_dispatch_index(
                pending,
                idle_workers=len(idle_workers),
                dispatch_policy=dispatch_policy,
            )
            if fit_idx < 0:
                break

            task = pending.pop(fit_idx)
            want = max(1, int(task.get("want_cores", task["cores"])))
            floor = max(1, int(task.get("min_cores", want)))
            elastic = bool(task.get("elastic", False))
            granted = min(len(idle_workers), want) if elastic else want
            granted = max(floor, granted)
            if granted > len(idle_workers):
                pending.insert(fit_idx, task)
                break

            worker_slice = sorted(idle_workers[:granted])
            del idle_workers[:granted]

            for wid in worker_slice:
                emit_event("WORKER_IDLE_END", worker_id=wid)

            start_ns = time.monotonic_ns()
            emit_event(
                "TASK_START",
                worker_id=worker_slice[0] if worker_slice else "",
                worker_ids=worker_slice,
                task_id=task["id"],
                cores=granted,
                want_cores=want,
                min_cores=floor,
                granted_cores=granted,
                fit_class=_fit_class_for_cores(want),
                est_duration_sec=max(0.001, float(task.get("expected_sec") or 1.0)),
                task_type="synthetic",
            )

            command = task["command"]
            if elastic and granted != want:
                command = _command_with_granted_workers(command, granted)

            try:
                proc = subprocess.Popen(
                    command,
                    shell=True,
                    executable="/bin/bash",
                )
            except Exception as exc:
                emit_event(
                    "TASK_END",
                    worker_id=worker_slice[0] if worker_slice else "",
                    worker_ids=worker_slice,
                    task_id=task["id"],
                    cores=granted,
                    want_cores=want,
                    min_cores=floor,
                    granted_cores=granted,
                    fit_class=_fit_class_for_cores(want),
                    est_duration_sec=max(0.001, float(task.get("expected_sec") or 1.0)),
                    task_type="synthetic",
                    status="spawn_error",
                    returncode=None,
                    error=str(exc),
                    duration_sec=0.0,
                )
                for wid in worker_slice:
                    idle_workers.append(wid)
                    emit_event("WORKER_IDLE_START", worker_id=wid)
                continue

            running.append(
                {
                    "task": task,
                    "proc": proc,
                    "workers": worker_slice,
                    "start_ns": start_ns,
                    "granted_cores": granted,
                }
            )
            dispatched = True

        for job in list(running):
            proc = job["proc"]
            rc = proc.poll()
            if rc is None:
                continue

            running.remove(job)
            task = job["task"]
            workers = list(job["workers"])
            elapsed = max(0.0, (time.monotonic_ns() - int(job["start_ns"])) / 1_000_000_000.0)

            emit_event(
                "TASK_END",
                worker_id=workers[0] if workers else "",
                worker_ids=workers,
                task_id=task["id"],
                cores=int(job.get("granted_cores", task["cores"])),
                want_cores=int(task.get("want_cores", task["cores"])),
                min_cores=int(task.get("min_cores", task.get("cores", 1))),
                granted_cores=int(job.get("granted_cores", task["cores"])),
                fit_class=_fit_class_for_cores(int(task.get("want_cores", task["cores"]))),
                est_duration_sec=max(0.001, float(task.get("expected_sec") or elapsed)),
                task_type="synthetic",
                status="ok" if rc == 0 else "error",
                returncode=rc,
                duration_sec=elapsed,
            )

            for wid in workers:
                idle_workers.append(wid)
                emit_event("WORKER_IDLE_START", worker_id=wid)

        now = time.monotonic()
        if now - last_snapshot >= max(0.05, snapshot_sec):
            running_cores = sum(int(j.get("granted_cores", j["task"]["cores"])) for j in running)
            emit_event(
                "QUEUE_DEPTH",
                depth=len(pending),
                running=len(running),
                running_cores=running_cores,
                idle_workers=len(idle_workers),
                free_cores=max(0, alloc_cpus - running_cores),
            )
            last_snapshot = now

        if pending or running:
            if not dispatched:
                time.sleep(0.02)



def _build_meta(run_id: str, alloc_cpus: int, interval_sec: float) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parents[3]
    return {
        "run_id": run_id,
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "alloc_cpus": int(alloc_cpus),
        "interval_sec": float(interval_sec),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_cpus_on_node": os.environ.get("SLURM_CPUS_ON_NODE"),
        "git_sha": _git_sha(repo_root),
        "created_wall_ns": time.time_ns(),
    }



def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Atlas utilization benchmark runner")
    ap.add_argument("--workload", required=True, help="Path to workload manifest JSON")
    ap.add_argument("--outdir", default="chemdb/bench/runs", help="Output root directory")
    ap.add_argument("--alloc-cpus", type=int, default=None)
    ap.add_argument("--interval-sec", type=float, default=1.0)
    ap.add_argument(
        "--dispatch-policy",
        choices=("backfill", "strict_fifo"),
        default="backfill",
        help="Synthetic scheduler admission policy",
    )
    ap.add_argument(
        "--force-procstat",
        action="store_true",
        help="Force /proc/stat sampler fallback even if mpstat is available",
    )
    ap.add_argument("--run-id", default="", help="Optional explicit run id")
    args = ap.parse_args(argv)

    workload = load_workload(args.workload)
    alloc_cpus = _detect_alloc_cpus(args.alloc_cpus)

    run_name = str(workload.get("run", {}).get("name") or "utilbench").strip()
    run_id = args.run_id.strip() or _now_run_id(prefix=run_name)

    out_root = Path(args.outdir)
    run_dir = out_root / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    meta = _build_meta(
        run_id=run_id,
        alloc_cpus=alloc_cpus,
        interval_sec=args.interval_sec,
    )
    meta["dispatch_policy"] = str(args.dispatch_policy)
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (run_dir / "workload.json").write_text(
        json.dumps(workload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    sampler = start_sampler(
        run_dir=run_dir,
        interval_sec=float(args.interval_sec),
        alloc_cpus=alloc_cpus,
        force_procstat=bool(args.force_procstat),
    )

    try:
        with event_run(run_dir / "events.jsonl", run_meta=meta):
            _run_synthetic_scheduler(
                tasks=list(workload.get("tasks", [])),
                alloc_cpus=alloc_cpus,
                snapshot_sec=max(0.05, float(args.interval_sec) / 2.0),
                dispatch_policy=str(args.dispatch_policy),
            )
    finally:
        sampler.stop()
        sampler.join(timeout=10)

    summary = analyze_run(run_dir=run_dir, alloc_cpus=alloc_cpus, interval_sec=float(args.interval_sec))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
