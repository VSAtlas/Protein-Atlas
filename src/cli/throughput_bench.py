from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import re
import shutil
import shlex
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

from config.output_paths import output_root


_PROFILE_MAIN_ARGS = {
    "micro": ["-bench-micro", "-fast"],
    "smoke": ["-bench2", "-fast"],
    "bench2-smoke": ["-bench2", "-fast"],
    "medium": ["-bench2"],
    "bench2-medium": ["-bench2"],
}


def run(argv: Sequence[str], repo_root: Path) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas throughput",
        description="Throughput benchmark harness for local, Slurm-simulated, and Slurm dry-run scheduling.",
    )
    sub = parser.add_subparsers(dest="throughput_cmd", required=True)

    bench = sub.add_parser("bench", help="Run a throughput benchmark profile.")
    bench.add_argument(
        "--profile",
        choices=tuple(_PROFILE_MAIN_ARGS),
        default="smoke",
        help="Workload profile. micro uses -bench-micro -fast; smoke uses -bench2 -fast; medium uses -bench2.",
    )
    bench.add_argument(
        "--mode",
        choices=("local", "slurm-sim", "slurm-dry-run", "all"),
        default="local",
    )
    bench.add_argument("--run-id", default="")
    bench.add_argument("--array", default="", help="Slurm-style array spec, e.g. 0-31%%8.")
    bench.add_argument("--cpus-per-task", type=int, default=4)
    bench.add_argument(
        "--seed",
        type=int,
        default=1337,
        help="Fixed Vina seed for throughput benchmark profiles only.",
    )
    bench.add_argument(
        "--prep-cache-run-id",
        default=os.environ.get("ATLAS_BENCH_MICRO_PREP_CACHE_RUN_ID", ""),
        help="For the micro profile, copy cached receptor prep artifacts from this run id before launching.",
    )
    bench.add_argument(
        "--allow-fresh-prep",
        action="store_true",
        help="For the micro profile, disable cached-prep fail-fast and allow receptor prep.",
    )
    bench.add_argument(
        "--execution-mode",
        choices=("auto", "distributed_per_protein", "distributed_combo", "distributed_combo_hybrid"),
        default="auto",
        help="ATLAS_EXECUTION_MODE used for local Slurm simulation and Slurm dry-run previews.",
    )
    bench.add_argument(
        "--main-args",
        default="",
        help="Override the profile's main.py arguments. Do not include --run-id unless intentional.",
    )
    bench.add_argument("--serial", action="store_true", help="Run Slurm simulation serially.")
    bench.add_argument(
        "--dry-run",
        action="store_true",
        help="Write benchmark config and planned commands without running local or simulated workers.",
    )
    bench.add_argument(
        "--skip-integrity",
        action="store_true",
        help="Do not run atlas analysis throughput integrity after completed local/sim runs.",
    )
    bench.add_argument(
        "--baseline-run-id",
        default="",
        help="Optional baseline run id for throughput acceptance comparison.",
    )
    bench.add_argument(
        "--util-bench",
        action="store_true",
        help="Also run the synthetic util-bench scheduler workload into the artifact bundle.",
    )
    bench.add_argument("--util-alloc-cpus", type=int, default=0)
    bench.add_argument("--util-interval-sec", type=float, default=0.2)
    bench.add_argument(
        "--util-dispatch-policy",
        choices=("backfill", "strict_fifo"),
        default="backfill",
    )
    args = parser.parse_args(list(argv))

    if args.throughput_cmd != "bench":
        return 1
    return _run_bench(args, Path(repo_root))


def _run_bench(args: argparse.Namespace, repo_root: Path) -> int:
    cpus_per_task = int(args.cpus_per_task)
    if cpus_per_task <= 0:
        print("atlas throughput bench: --cpus-per-task must be > 0", file=sys.stderr)
        return 2
    if cpus_per_task > 32:
        print("atlas throughput bench: --cpus-per-task is capped at 32", file=sys.stderr)
        return 2

    run_id = str(args.run_id or _default_run_id(args.profile))
    bench_dir = output_root(repo_root, "data") / run_id / "throughput_bench"
    bench_dir.mkdir(parents=True, exist_ok=True)

    main_args = _main_args_for(args)
    try:
        config = _bench_config(args, repo_root, run_id, main_args)
    except ValueError as exc:
        print(f"atlas throughput bench: {exc}", file=sys.stderr)
        return 2
    _write_json(bench_dir / "bench_config.json", config)

    modes = ["local", "slurm-sim", "slurm-dry-run"] if args.mode == "all" else [str(args.mode)]
    exit_code = 0
    util_summary = ""
    prep_cache = _prepare_micro_prep_cache(args, repo_root, bench_dir, run_id)
    will_execute_main = any(mode in {"local", "slurm-sim"} for mode in modes) and not bool(
        args.dry_run
    )
    if prep_cache.get("required") and will_execute_main and not prep_cache.get("ready"):
        print(
            "atlas throughput bench: micro profile requires cached receptor prep; "
            "use --prep-cache-run-id <RUN_ID> or --allow-fresh-prep",
            file=sys.stderr,
        )
        _write_readiness_report(args, repo_root, bench_dir, run_id, main_args, False)
        print(f"throughput bench artifacts: {_display_path(bench_dir, repo_root)}")
        return 2

    if args.util_bench and not args.dry_run:
        util_rc, util_summary = _run_util_bench(args, bench_dir, run_id)
        exit_code = max(exit_code, util_rc)

    completed_run = False
    for mode in modes:
        if mode == "local":
            rc = _run_local(
                repo_root,
                bench_dir,
                run_id,
                main_args,
                total_cpus=int(args.cpus_per_task),
                seed=int(args.seed),
                dry_run=bool(args.dry_run),
                require_cached_prep=_requires_cached_prep(args),
            )
            completed_run = completed_run or (rc == 0 and not args.dry_run)
        elif mode == "slurm-sim":
            rc = _run_slurm_sim(repo_root, bench_dir, run_id, main_args, args)
            completed_run = completed_run or (rc == 0 and not args.dry_run)
        elif mode == "slurm-dry-run":
            rc = _run_slurm_dry_run(repo_root, bench_dir, run_id, main_args, args)
        else:
            rc = 1
        exit_code = max(exit_code, int(rc))

    if completed_run and not args.skip_integrity:
        integrity_rc = _run_integrity(repo_root, bench_dir, run_id)
        exit_code = max(exit_code, integrity_rc)

    if completed_run and args.baseline_run_id:
        acceptance_rc = _run_acceptance(
            repo_root=repo_root,
            bench_dir=bench_dir,
            baseline_run_id=str(args.baseline_run_id),
            candidate_run_id=run_id,
            util_candidate_summary=util_summary,
        )
        exit_code = max(exit_code, acceptance_rc)

    _write_readiness_report(args, repo_root, bench_dir, run_id, main_args, completed_run)
    print(f"throughput bench artifacts: {_display_path(bench_dir, repo_root)}")
    return exit_code


def _default_run_id(profile: str) -> str:
    stamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    normalized = str(profile).replace("bench2-", "")
    return f"throughput_bench_{normalized}_{stamp}"


def _main_args_for(args: argparse.Namespace) -> list[str]:
    if args.main_args:
        return shlex.split(str(args.main_args))
    return list(_PROFILE_MAIN_ARGS[str(args.profile)])


def _bench_config(
    args: argparse.Namespace,
    repo_root: Path,
    run_id: str,
    main_args: list[str],
) -> dict[str, Any]:
    array_spec = _array_for(args.profile, str(args.array or ""))
    concurrency = _array_concurrency(array_spec)
    return {
        "array": array_spec,
        "baseline_run_id": str(args.baseline_run_id or ""),
        "cpus_per_task": int(args.cpus_per_task),
        "dry_run": bool(args.dry_run),
        "main_args": main_args,
        "main_args_text": shlex.join(main_args),
        "mode": str(args.mode),
        "profile": str(args.profile),
        "allow_fresh_prep": bool(args.allow_fresh_prep),
        "prep_cache_run_id": str(args.prep_cache_run_id or ""),
        "repo_root": str(repo_root),
        "run_id": run_id,
        "seed": int(args.seed),
        "simulated_max_concurrent_cpus": int(concurrency) * int(args.cpus_per_task),
        "workflow": "bench_micro" if str(args.profile) == "micro" else "bench2",
        "execution_mode": str(args.execution_mode),
    }


def _array_for(profile: str, explicit: str) -> str:
    if explicit:
        return explicit
    if str(profile) == "micro":
        return "0-1%2"
    if str(profile).endswith("medium"):
        return "0-31%4"
    return "0-2%3"


def _array_concurrency(array_spec: str) -> int:
    raw = str(array_spec or "").strip()
    limit = None
    if "%" in raw:
        raw, raw_limit = raw.rsplit("%", 1)
        limit = int(raw_limit)
    count = 0
    for item in raw.split(","):
        part = item.strip()
        if not part:
            continue
        match = re.fullmatch(r"(\d+)(?:-(\d+)(?::(\d+))?)?", part)
        if not match:
            raise ValueError(f"invalid --array item: {part!r}")
        start = int(match.group(1))
        end = int(match.group(2) or start)
        step = int(match.group(3) or 1)
        count += len(range(start, end + 1, step))
    if count < 1:
        raise ValueError("--array produced no tasks")
    return max(1, min(count, int(limit or count)))


def _ensure_run_id(main_args: list[str], run_id: str) -> list[str]:
    if "--run-id" in main_args or "-run-id" in main_args:
        return list(main_args)
    return [*main_args, "--run-id", run_id]


def _ensure_local_total_cpus(main_args: list[str], total_cpus: int) -> list[str]:
    if "--total-cpus" in main_args:
        return list(main_args)
    return [*main_args, "--total-cpus", str(max(1, min(32, int(total_cpus))))]


def _requires_cached_prep(args: argparse.Namespace) -> bool:
    return str(args.profile) == "micro" and not bool(args.allow_fresh_prep)


def _micro_pdb_ids() -> list[str]:
    raw = os.environ.get("ATLAS_BENCH_MICRO_PDBS", "bOJG")
    pdb_ids = [token.strip().upper() for token in re.split(r"[,\s]+", raw) if token.strip()]
    return pdb_ids or ["BOJG"]


def _default_prep_cache_run_id(repo_root: Path) -> str:
    default_run = "bench2_fast_cpu_scorch_20260506"
    if (output_root(repo_root, "processed_pdbs") / default_run).exists():
        return default_run
    return ""


def _prepare_micro_prep_cache(
    args: argparse.Namespace, repo_root: Path, bench_dir: Path, run_id: str
) -> dict[str, Any]:
    if str(args.profile) != "micro":
        return {"profile": str(args.profile), "required": False, "ready": True}

    required = _requires_cached_prep(args)
    processed_root = output_root(repo_root, "processed_pdbs")
    source_run_id = str(args.prep_cache_run_id or "").strip() or _default_prep_cache_run_id(
        repo_root
    )
    copied: list[str] = []
    missing: list[str] = []

    for pdb_id in _micro_pdb_ids():
        receptor_rel = Path(pdb_id) / "HOLO" / "receptor"
        dest_dir = processed_root / run_id / receptor_rel
        cleaned = dest_dir / f"{pdb_id}_cleaned.pdb"
        pdbqt = dest_dir / f"{pdb_id}.pdbqt"
        if source_run_id and (not cleaned.exists() or not pdbqt.exists()):
            source_dir = processed_root / source_run_id / receptor_rel
            if source_dir.exists():
                shutil.copytree(source_dir, dest_dir, dirs_exist_ok=True)
                copied.append(str(source_dir))
        if not cleaned.exists():
            missing.append(str(cleaned))
        if not pdbqt.exists():
            missing.append(str(pdbqt))

    payload = {
        "profile": "micro",
        "required": bool(required),
        "ready": not missing,
        "allow_fresh_prep": bool(args.allow_fresh_prep),
        "source_run_id": source_run_id,
        "copied_receptor_dirs": copied,
        "missing": missing,
    }
    _write_json(bench_dir / "prep_cache_summary.json", payload)
    return payload


def _run_local(
    repo_root: Path,
    bench_dir: Path,
    run_id: str,
    main_args: list[str],
    *,
    total_cpus: int,
    seed: int,
    dry_run: bool,
    require_cached_prep: bool,
) -> int:
    local_args = _ensure_local_total_cpus(main_args, int(total_cpus))
    cmd = [sys.executable, str(repo_root / "main.py"), *_ensure_run_id(local_args, run_id)]
    env = _execution_env(
        repo_root,
        "auto",
        seed=int(seed),
        require_cached_prep=bool(require_cached_prep),
    )
    summary_path = bench_dir / "local_summary.json"
    stdout_path = bench_dir / "local_stdout.log"
    stderr_path = bench_dir / "local_stderr.log"
    start = time.perf_counter()

    if dry_run:
        _write_json(
            summary_path,
            _command_summary(cmd=cmd, elapsed=0.0, exit_code=None, dry_run=True),
        )
        return 0

    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        proc = subprocess.run(
            cmd, cwd=repo_root, check=False, stdout=stdout, stderr=stderr, env=env
        )
    elapsed = time.perf_counter() - start
    _write_json(
        summary_path,
        {
            **_command_summary(cmd=cmd, elapsed=elapsed, exit_code=int(proc.returncode), dry_run=False),
            "stderr": str(stderr_path),
            "stdout": str(stdout_path),
        },
    )
    return int(proc.returncode)


def _run_slurm_sim(
    repo_root: Path,
    bench_dir: Path,
    run_id: str,
    main_args: list[str],
    args: argparse.Namespace,
) -> int:
    array_spec = _array_for(str(args.profile), str(args.array or ""))
    cap_error = _concurrent_cpu_cap_error(array_spec, int(args.cpus_per_task))
    if cap_error:
        print(cap_error, file=sys.stderr)
        return 2

    main_args_path = bench_dir / "main_args.txt"
    main_args_path.write_text(shlex.join(main_args) + "\n", encoding="utf-8")
    summary_path = bench_dir / "slurm_sim_summary.json"
    cmd = [
        sys.executable,
        str(repo_root / "tools" / "simulate_slurm_array.py"),
        "--array",
        array_spec,
        "--cpus-per-task",
        str(int(args.cpus_per_task)),
        "--run-id",
        run_id,
        "--main-args-file",
        str(main_args_path),
        "--summary-json",
        str(summary_path),
    ]
    if args.serial:
        cmd.append("--serial")
    if args.dry_run:
        cmd.append("--dry-run")
    env = _execution_env(
        repo_root,
        str(args.execution_mode),
        seed=int(args.seed),
        require_cached_prep=_requires_cached_prep(args),
    )
    return _run_logged(
        cmd,
        repo_root,
        bench_dir / "slurm_sim_stdout.log",
        bench_dir / "slurm_sim_stderr.log",
        env=env,
    )


def _run_slurm_dry_run(
    repo_root: Path,
    bench_dir: Path,
    run_id: str,
    main_args: list[str],
    args: argparse.Namespace,
) -> int:
    array_spec = _array_for(str(args.profile), str(args.array or ""))
    cap_error = _concurrent_cpu_cap_error(array_spec, int(args.cpus_per_task))
    if cap_error:
        print(cap_error, file=sys.stderr)
        return 2
    env = dict(os.environ)
    src_path = str(repo_root / "src")
    env["PYTHONPATH"] = src_path + os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else src_path
    env["ATLAS_THROUGHPUT_BENCH"] = "1"
    env["ATLAS_THROUGHPUT_BENCH_SEED"] = str(int(args.seed))
    if str(args.profile) == "micro":
        env["ATLAS_BENCH_CACHED_PREP_REQUIRED"] = "1" if _requires_cached_prep(args) else "0"
    cmd = [
        sys.executable,
        "-m",
        "cli.atlas_main_cli",
        "slurm",
        "submit",
        "--run-id",
        run_id,
        "--array",
        array_spec,
        "--cpus-per-task",
        str(int(args.cpus_per_task)),
        "--main-args",
        shlex.join(_ensure_run_id(main_args, run_id)),
        "--dry-run",
    ]
    if str(args.execution_mode or "").strip() != "auto":
        cmd.extend(["--execution-mode", str(args.execution_mode)])
    preview_path = bench_dir / "slurm_submit_preview.txt"
    err_path = bench_dir / "slurm_submit_preview.stderr"
    with preview_path.open("w", encoding="utf-8") as stdout, err_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        proc = subprocess.run(cmd, cwd=repo_root, env=env, check=False, stdout=stdout, stderr=stderr)
    return int(proc.returncode)


def _concurrent_cpu_cap_error(array_spec: str, cpus_per_task: int) -> str:
    concurrency = _array_concurrency(array_spec)
    concurrent_cpus = concurrency * int(cpus_per_task)
    if concurrent_cpus <= 32:
        return ""
    return (
        "atlas throughput bench: simulated/planned concurrent CPU budget "
        f"exceeds 32 ({concurrent_cpus})"
    )


def _execution_env(
    repo_root: Path,
    execution_mode: str,
    *,
    seed: int,
    require_cached_prep: bool | None = None,
) -> dict[str, str]:
    env = dict(os.environ)
    src_path = str(repo_root / "src")
    env["PYTHONPATH"] = src_path + os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else src_path
    env["ATLAS_THROUGHPUT_BENCH"] = "1"
    env["ATLAS_THROUGHPUT_BENCH_SEED"] = str(int(seed))
    if require_cached_prep is not None:
        env["ATLAS_BENCH_CACHED_PREP_REQUIRED"] = "1" if require_cached_prep else "0"
    mode = str(execution_mode or "").strip()
    if mode and mode != "auto":
        env["ATLAS_EXECUTION_MODE"] = mode
    return env


def _run_logged(
    cmd: list[str],
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
    *,
    env: dict[str, str] | None = None,
) -> int:
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        proc = subprocess.run(cmd, cwd=cwd, check=False, stdout=stdout, stderr=stderr, env=env)
    return int(proc.returncode)


def _run_integrity(repo_root: Path, bench_dir: Path, run_id: str) -> int:
    from analysis.reporting import throughput_integrity

    stdout = io.StringIO()
    stderr = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        rc = int(
            throughput_integrity.main(
                ["--run-id", run_id, "--repo-root", str(repo_root), "--strict"]
            )
        )
    (bench_dir / "throughput_integrity.stdout").write_text(stdout.getvalue(), encoding="utf-8")
    (bench_dir / "throughput_integrity.stderr").write_text(stderr.getvalue(), encoding="utf-8")
    return rc


def _run_acceptance(
    *,
    repo_root: Path,
    bench_dir: Path,
    baseline_run_id: str,
    candidate_run_id: str,
    util_candidate_summary: str,
) -> int:
    from analysis.reporting import throughput_acceptance

    out_json = bench_dir / "throughput_acceptance.json"
    forwarded = [
        "--repo-root",
        str(repo_root),
        "--baseline-run-id",
        baseline_run_id,
        "--candidate-run-id",
        candidate_run_id,
        "--out-json",
        str(out_json),
    ]
    if util_candidate_summary:
        forwarded.extend(["--util-candidate-summary", util_candidate_summary])

    stdout = io.StringIO()
    stderr = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        rc = int(throughput_acceptance.main(forwarded))
    (bench_dir / "throughput_acceptance.stdout").write_text(stdout.getvalue(), encoding="utf-8")
    (bench_dir / "throughput_acceptance.stderr").write_text(stderr.getvalue(), encoding="utf-8")
    return rc


def _run_util_bench(args: argparse.Namespace, bench_dir: Path, run_id: str) -> tuple[int, str]:
    from chemdb.bench.util_bench import run_bench, workloads

    profile = str(args.profile)
    stage1_tasks = 96 if not profile.endswith("medium") else 384
    workload = workloads.generate_stage_collapse_workload(
        proteins=3,
        stage1_tasks=stage1_tasks,
        seed=1337,
        run_name=f"{run_id}_util",
        stage2_keep_frac=0.25,
        stage3_keep_frac=0.10,
    )
    workload_path = bench_dir / "util_bench_workload.json"
    workloads.write_workload_json(workload_path, workload)
    outdir = bench_dir / "util_bench"
    util_run_id = f"{run_id}_util"
    forwarded = [
        "--workload",
        str(workload_path),
        "--outdir",
        str(outdir),
        "--run-id",
        util_run_id,
        "--interval-sec",
        str(float(args.util_interval_sec)),
        "--dispatch-policy",
        str(args.util_dispatch_policy),
    ]
    if int(args.util_alloc_cpus or 0) > 0:
        forwarded.extend(["--alloc-cpus", str(int(args.util_alloc_cpus))])
    stdout = io.StringIO()
    stderr = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        rc = int(run_bench.main(forwarded))
    (bench_dir / "util_bench.stdout").write_text(stdout.getvalue(), encoding="utf-8")
    (bench_dir / "util_bench.stderr").write_text(stderr.getvalue(), encoding="utf-8")
    summary_path = outdir / util_run_id / "summary.json"
    return rc, str(summary_path) if summary_path.exists() else ""


def _write_readiness_report(
    args: argparse.Namespace,
    repo_root: Path,
    bench_dir: Path,
    run_id: str,
    main_args: list[str],
    completed_run: bool,
) -> None:
    array_spec = _array_for(str(args.profile), str(args.array or ""))
    concurrency = _array_concurrency(array_spec)
    concurrent_cpus = int(concurrency) * int(args.cpus_per_task)
    integrity_path = output_root(repo_root, "data") / run_id / "throughput_integrity.json"
    acceptance_path = bench_dir / "throughput_acceptance.json"
    run_efficiency_path = output_root(repo_root, "manifests") / run_id / "run_efficiency.json"
    submit_preview = bench_dir / "slurm_submit_preview.txt"
    prep_cache_path = bench_dir / "prep_cache_summary.json"

    integrity = _read_json_if_exists(integrity_path)
    acceptance = _read_json_if_exists(acceptance_path)
    efficiency = _read_json_if_exists(run_efficiency_path)
    prep_cache = _read_json_if_exists(prep_cache_path)
    slurm_raw = efficiency.get("slurm")
    slurm_eff: dict[str, Any] = slurm_raw if isinstance(slurm_raw, dict) else {}

    gaps: list[str] = []
    if concurrent_cpus > 32:
        gaps.append("planned concurrent CPU budget exceeds Atlas 32-core docking cap")
    if str(args.mode) != "slurm-dry-run" and not submit_preview.exists():
        gaps.append("Slurm submission preview has not been generated for this run")
    if not completed_run:
        gaps.append("no completed local or simulated run is available for integrity evidence")
    if completed_run and not bool(integrity.get("integrity_pass")):
        gaps.append("strict throughput integrity has not passed")
    if args.baseline_run_id and not bool(acceptance.get("acceptance_pass")):
        gaps.append("baseline acceptance has not passed")
    if not bool(slurm_eff.get("sacct_available", False)):
        gaps.append("real Slurm sacct accounting is unavailable until a cluster run completes")
    gaps.append("local Slurm simulation does not validate queue wait, controller limits, filesystem contention, or accounting")
    if str(args.profile) == "micro":
        gaps.append("micro profile is for scheduler/dev feedback only; use smoke or medium for performance claims")
        if bool(prep_cache.get("required")) and not bool(prep_cache.get("ready")):
            gaps.append("micro cached receptor prep is not ready")

    payload = {
        "run_id": run_id,
        "profile": str(args.profile),
        "main_args": main_args,
        "seed": int(args.seed),
        "array": array_spec,
        "array_concurrency": int(concurrency),
        "cpus_per_task": int(args.cpus_per_task),
        "concurrent_cpu_budget": int(concurrent_cpus),
        "execution_mode": str(args.execution_mode),
        "completed_local_or_simulated_run": bool(completed_run),
        "integrity_pass": integrity.get("integrity_pass"),
        "acceptance_pass": acceptance.get("acceptance_pass"),
        "run_efficiency_json": str(run_efficiency_path) if run_efficiency_path.exists() else "",
        "sacct_available": bool(slurm_eff.get("sacct_available", False)),
        "prep_cache_summary": str(prep_cache_path) if prep_cache_path.exists() else "",
        "prep_cache_ready": prep_cache.get("ready"),
        "slurm_submit_preview": str(submit_preview) if submit_preview.exists() else "",
        "remaining_gaps": gaps,
        "references": {
            "slurm_job_arrays": "https://slurm.schedmd.com/job_array.html",
            "slurm_sacct": "https://slurm.schedmd.com/sacct.html",
            "nextflow_executor_model": "https://docs.seqera.io/nextflow/executor",
            "vina_seed_reproducibility": "https://vina.scripps.edu/manual/",
        },
    }
    _write_json(bench_dir / "slurm_readiness.json", payload)


def _read_json_if_exists(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _command_summary(
    *,
    cmd: list[str],
    elapsed: float,
    exit_code: int | None,
    dry_run: bool,
) -> dict[str, Any]:
    return {
        "command": cmd,
        "command_text": shlex.join(cmd),
        "dry_run": bool(dry_run),
        "elapsed_wall_sec": round(float(elapsed), 6),
        "exit_code": exit_code,
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _display_path(path: Path, repo_root: Path) -> str:
    try:
        return path.relative_to(repo_root).as_posix()
    except ValueError:
        return str(path)
