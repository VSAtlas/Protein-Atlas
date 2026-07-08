#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


def _build_cmd(repo_root: Path, cmd_args: list[str]) -> list[str]:
    return [sys.executable, str(repo_root / "main.py"), *cmd_args]


def _log_task_cmd(task_id: int, task_max: int, cmd: list[str]) -> None:
    print(
        "[simulate-slurm-array] "
        f"task={task_id}/{task_max} cmd={shlex.join(cmd)}"
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Local simulator for Slurm array execution. "
            "Spawns N local workers with SLURM_ARRAY_* env vars."
        )
    )
    parser.add_argument(
        "--tasks",
        type=int,
        default=2,
        help="Number of simulated array tasks. Ignored when --array is supplied.",
    )
    parser.add_argument(
        "--array",
        default="",
        help="Slurm-style array spec, for example 0-31 or 0-31%%8.",
    )
    parser.add_argument(
        "--cpus-per-task",
        type=int,
        default=4,
        help="Value exposed as SLURM_CPUS_PER_TASK.",
    )
    parser.add_argument(
        "--run-id",
        default="sim_array",
        help="Shared run id used by all local workers.",
    )
    parser.add_argument(
        "--main-args",
        default="",
        help='Arguments passed to main.py, for example "--fast --test-fda --pdbs 1BN1,2OJ9".',
    )
    parser.add_argument(
        "--main-args-file",
        default="",
        help="Read main.py arguments from this file instead of --main-args.",
    )
    parser.add_argument(
        "--summary-json",
        default="",
        help="Optional path for a JSON summary of the simulated array run.",
    )
    parser.add_argument(
        "--serial",
        action="store_true",
        help="Run workers serially instead of concurrently.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned commands and exit.",
    )
    return parser.parse_args()


def _parse_array_spec(raw: str) -> tuple[list[int], int | None]:
    spec = str(raw or "").strip()
    if not spec:
        raise ValueError("array spec is empty")

    concurrency: int | None = None
    if "%" in spec:
        spec, raw_limit = spec.rsplit("%", 1)
        try:
            concurrency = int(raw_limit)
        except ValueError as exc:
            raise ValueError(f"invalid array concurrency limit: {raw_limit!r}") from exc
        if concurrency < 1:
            raise ValueError("--array concurrency limit must be >= 1")

    task_ids: list[int] = []
    for part in spec.split(","):
        item = part.strip()
        if not item:
            continue
        match = re.fullmatch(r"(\d+)(?:-(\d+)(?::(\d+))?)?", item)
        if not match:
            raise ValueError(f"invalid array item: {item!r}")
        start = int(match.group(1))
        end = int(match.group(2) or start)
        step = int(match.group(3) or 1)
        if step < 1:
            raise ValueError(f"invalid array step in {item!r}")
        if end < start:
            raise ValueError(f"array range end precedes start in {item!r}")
        task_ids.extend(range(start, end + 1, step))

    if not task_ids:
        raise ValueError("array spec produced no task IDs")
    return task_ids, concurrency


def _ensure_run_id(args: list[str], run_id: str) -> list[str]:
    tokens = list(args)
    has_run_id = "--run-id" in tokens or "-run-id" in tokens
    if has_run_id:
        return tokens
    return ["--run-id", run_id, *tokens]


def _read_main_args(ns: argparse.Namespace) -> str:
    if ns.main_args and ns.main_args_file:
        raise ValueError("--main-args and --main-args-file are mutually exclusive")
    if ns.main_args_file:
        path = Path(str(ns.main_args_file)).expanduser()
        return path.read_text(encoding="utf-8").strip()
    if ns.main_args:
        return str(ns.main_args)
    raise ValueError("one of --main-args or --main-args-file is required")


def _resolve_task_ids(ns: argparse.Namespace) -> tuple[list[int], int | None]:
    if ns.array:
        return _parse_array_spec(str(ns.array))
    tasks = max(1, int(ns.tasks))
    return list(range(tasks)), None


def _build_env(
    task_id: int,
    task_ids: list[int],
    cpus_per_task: int,
) -> dict[str, str]:
    env = dict(os.environ)
    env["ATLAS_DISTRIBUTED_MODE"] = "slurm_array"
    env["SLURM_ARRAY_TASK_ID"] = str(task_id)
    env["SLURM_ARRAY_TASK_COUNT"] = str(len(task_ids))
    env["SLURM_ARRAY_TASK_MIN"] = str(min(task_ids))
    env["SLURM_ARRAY_TASK_MAX"] = str(max(task_ids))
    env["SLURM_ARRAY_JOB_ID"] = f"local{os.getpid()}"
    env["SLURM_CPUS_PER_TASK"] = str(max(1, int(cpus_per_task)))
    return env


def _run_serial(
    repo_root: Path,
    cmd_args: list[str],
    task_ids: list[int],
    cpus_per_task: int,
) -> dict[int, int | None]:
    exit_codes: dict[int, int | None] = {}
    for task_id in task_ids:
        env = _build_env(task_id, task_ids, cpus_per_task)
        cmd = _build_cmd(repo_root, cmd_args)
        _log_task_cmd(task_id, max(task_ids), cmd)
        proc = subprocess.run(cmd, cwd=str(repo_root), env=env, check=False)
        exit_codes[task_id] = int(proc.returncode)
    return exit_codes


def _run_parallel(
    repo_root: Path,
    cmd_args: list[str],
    task_ids: list[int],
    cpus_per_task: int,
    max_concurrency: int,
) -> dict[int, int | None]:
    procs: list[tuple[int, subprocess.Popen[str]]] = []
    pending = list(task_ids)
    exit_codes: dict[int, int | None] = {}

    try:
        while pending or procs:
            while pending and len(procs) < max(1, int(max_concurrency)):
                task_id = pending.pop(0)
                env = _build_env(task_id, task_ids, cpus_per_task)
                cmd = _build_cmd(repo_root, cmd_args)
                _log_task_cmd(task_id, max(task_ids), cmd)
                proc = subprocess.Popen(
                    cmd,
                    cwd=str(repo_root),
                    env=env,
                    text=True,
                )
                procs.append((task_id, proc))
                time.sleep(0.15)

            for task_id, proc in list(procs):
                code = proc.poll()
                if code is None:
                    continue
                procs.remove((task_id, proc))
                exit_codes[task_id] = int(code)
                print(f"[simulate-slurm-array] task={task_id} exit_code={code}")

            if pending or procs:
                time.sleep(0.05)
    except BaseException:
        for _task_id, proc in procs:
            if proc.poll() is None:
                proc.terminate()
        for _task_id, proc in procs:
            if proc.poll() is None:
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
        raise

    return exit_codes


def _summary_payload(
    *,
    repo_root: Path,
    cmd_args: list[str],
    task_ids: list[int],
    cpus_per_task: int,
    array_spec: str,
    max_concurrency: int,
    elapsed_wall_sec: float,
    exit_codes: dict[int, int | None],
    dry_run: bool,
) -> dict[str, Any]:
    command = _build_cmd(repo_root, cmd_args)
    env0 = _build_env(task_ids[0], task_ids, cpus_per_task)
    return {
        "array": array_spec or f"0-{len(task_ids) - 1}",
        "command": command,
        "command_text": shlex.join(command),
        "cpus_per_task": int(cpus_per_task),
        "dry_run": bool(dry_run),
        "elapsed_wall_sec": round(float(elapsed_wall_sec), 6),
        "exit_codes": {str(k): v for k, v in sorted(exit_codes.items())},
        "max_concurrency": int(max_concurrency),
        "repo_root": str(repo_root),
        "simulated_slurm_env": {
            key: env0[key]
            for key in (
                "ATLAS_DISTRIBUTED_MODE",
                "SLURM_ARRAY_JOB_ID",
                "SLURM_ARRAY_TASK_COUNT",
                "SLURM_ARRAY_TASK_MIN",
                "SLURM_ARRAY_TASK_MAX",
                "SLURM_CPUS_PER_TASK",
            )
        },
        "task_count": len(task_ids),
        "task_ids": task_ids,
    }


def _write_summary(path: str, payload: dict[str, Any]) -> None:
    if not path:
        return
    out = Path(path).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _max_exit_code(exit_codes: dict[int, int | None]) -> int:
    rc = 0
    for code in exit_codes.values():
        if code is None:
            continue
        rc = max(rc, int(code))
    return rc


def _dry_run(
    repo_root: Path,
    cmd_args: list[str],
    task_ids: list[int],
    cpus_per_task: int,
) -> dict[int, int | None]:
    exit_codes: dict[int, int | None] = {}
    for task_id in task_ids:
        env = _build_env(task_id, task_ids, cpus_per_task)
        cmd = _build_cmd(repo_root, cmd_args)
        print(
            "[simulate-slurm-array.dry-run] "
            f"task_id={task_id} cpus_per_task={env['SLURM_CPUS_PER_TASK']} cmd={shlex.join(cmd)}"
        )
        exit_codes[task_id] = None
    return exit_codes


def main() -> int:
    try:
        ns = _parse_args()
        repo_root = Path(__file__).resolve().parents[1]

        task_ids, array_concurrency = _resolve_task_ids(ns)
        cpus_per_task = max(1, int(ns.cpus_per_task))
        raw_main_args = shlex.split(_read_main_args(ns))
        cmd_args = _ensure_run_id(raw_main_args, str(ns.run_id))
        max_concurrency = (
            1
            if ns.serial
            else max(1, int(array_concurrency or len(task_ids)))
        )

        start = time.perf_counter()
        if ns.dry_run:
            exit_codes = _dry_run(repo_root, cmd_args, task_ids, cpus_per_task)
        elif ns.serial:
            exit_codes = _run_serial(repo_root, cmd_args, task_ids, cpus_per_task)
        else:
            exit_codes = _run_parallel(
                repo_root,
                cmd_args,
                task_ids,
                cpus_per_task,
                max_concurrency=max_concurrency,
            )
        elapsed = time.perf_counter() - start

        summary = _summary_payload(
            repo_root=repo_root,
            cmd_args=cmd_args,
            task_ids=task_ids,
            cpus_per_task=cpus_per_task,
            array_spec=str(ns.array or ""),
            max_concurrency=max_concurrency,
            elapsed_wall_sec=elapsed,
            exit_codes=exit_codes,
            dry_run=bool(ns.dry_run),
        )
        _write_summary(str(ns.summary_json or ""), summary)
        return 0 if ns.dry_run else _max_exit_code(exit_codes)
    except Exception as exc:
        print(f"simulate_slurm_array.py: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
