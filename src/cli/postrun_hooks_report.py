# -*- coding: utf-8 -*-
from __future__ import annotations

import datetime
import json
import logging
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from cli.cli_utils import _norm_pdb_id
from cli.postrun_hooks_common import _is_no_library_docking, _write_json_atomic
from cli.postrun_hooks_support import (
    _as_int,
    _parse_duration_seconds,
    _repo_root,
    _resolve_hook_roots,
    _with_repo_src_on_pythonpath,
)
from cli.run_manifest_runtime import get_manifest_paths, load_run_manifest


def _collect_sacct_cpu_metrics(
    alloc_cpus: int,
    logger: logging.Logger,
    artifact_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "slurm_job_id": None,
        "sacct_available": False,
        "sacct_command": None,
        "sacct_stdout": None,
        "sacct_stderr": None,
        "sacct_returncode": None,
        "elapsed_sec": None,
        "total_cpu_sec": None,
        "alloc_cpus": int(max(1, alloc_cpus)),
        "state": None,
        "max_rss": None,
        "req_mem": None,
        "timelimit_raw": None,
        "exit_code": None,
        "cpu_efficiency": None,
    }
    job_id = str(os.environ.get("SLURM_JOB_ID", "")).strip()
    if not job_id:
        return out
    out["slurm_job_id"] = job_id
    sacct_bin = shutil.which("sacct")
    if not sacct_bin:
        return out
    out["sacct_available"] = True
    cmd = [
        sacct_bin,
        "-X",
        "-j",
        str(job_id),
        "--parsable2",
        "--noheader",
        "--format=JobIDRaw,ElapsedRaw,TotalCPU,AllocCPUS,State,MaxRSS,ReqMem,TimelimitRaw,ExitCode",
    ]
    out["sacct_command"] = " ".join(cmd)
    try:
        proc = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
        )
    except Exception:
        logger.warning("[run-efficiency.sacct] action=invoke_failed", exc_info=True)
        return out
    out["sacct_returncode"] = int(proc.returncode)
    if artifact_dir is not None:
        try:
            artifact_dir.mkdir(parents=True, exist_ok=True)
            command_path = artifact_dir / "sacct_command.txt"
            stdout_path = artifact_dir / "sacct.stdout.psv"
            stderr_path = artifact_dir / "sacct.stderr.txt"
            command_path.write_text(" ".join(cmd) + "\n", encoding="utf-8")
            stdout_path.write_text(proc.stdout or "", encoding="utf-8")
            stderr_path.write_text(proc.stderr or "", encoding="utf-8")
            out["sacct_command"] = str(command_path)
            out["sacct_stdout"] = str(stdout_path)
            out["sacct_stderr"] = str(stderr_path)
        except Exception:
            logger.warning("[run-efficiency.sacct] action=artifact_write_failed", exc_info=True)
    if proc.returncode != 0:
        logger.warning(
            "[run-efficiency.sacct] action=nonzero_returncode rc=%s stderr=%s",
            proc.returncode,
            (proc.stderr or "").strip(),
        )
        return out

    rows: list[list[str]] = []
    for line in (proc.stdout or "").splitlines():
        token = line.strip()
        if not token:
            continue
        rows.append(token.split("|"))
    if not rows:
        return out

    selected = rows[0]
    for row in rows:
        if row and str(row[0]).strip() == job_id:
            selected = row
            break
    if len(selected) < 5:
        return out

    elapsed_sec = _parse_duration_seconds(selected[1])
    total_cpu_sec = _parse_duration_seconds(selected[2])
    alloc_row = _as_int(selected[3], int(max(1, alloc_cpus)))
    state = str(selected[4] or "").strip() or None
    out["elapsed_sec"] = None if elapsed_sec is None else float(elapsed_sec)
    out["total_cpu_sec"] = None if total_cpu_sec is None else float(total_cpu_sec)
    out["alloc_cpus"] = int(max(1, alloc_row))
    out["state"] = state
    if len(selected) > 5:
        out["max_rss"] = str(selected[5] or "").strip() or None
    if len(selected) > 6:
        out["req_mem"] = str(selected[6] or "").strip() or None
    if len(selected) > 7:
        out["timelimit_raw"] = str(selected[7] or "").strip() or None
    if len(selected) > 8:
        out["exit_code"] = str(selected[8] or "").strip() or None
    if (
        elapsed_sec is not None
        and total_cpu_sec is not None
        and elapsed_sec > 0
        and alloc_row > 0
    ):
        denom = float(alloc_row) * float(elapsed_sec)
        out["cpu_efficiency"] = float(total_cpu_sec) / denom if denom > 0 else None
    return out


def _write_run_efficiency_report(cfg: Mapping[str, Any], run_id: str) -> None:
    logger = logging.getLogger("run-efficiency-hook")
    if not run_id:
        logger.info("[run-efficiency.skip] reason=missing_run_id")
        return

    manifest = load_run_manifest(cfg, run_id) or {}
    summary = manifest.get("summary") if isinstance(manifest, dict) else {}
    summary = summary if isinstance(summary, dict) else {}
    timing = manifest.get("timing") if isinstance(manifest, dict) else {}
    timing = timing if isinstance(timing, dict) else {}
    wall_sec = timing.get("wall_time_sec")
    wall_sec_f = float(wall_sec) if isinstance(wall_sec, (int, float)) else None

    alloc_cpus = max(
        1,
        _as_int(
            os.environ.get("SLURM_CPUS_ON_NODE") or cfg.get("CPU") or (os.cpu_count() or 1),
            os.cpu_count() or 1,
        ),
    )
    combo_wall = summary.get("combo_wall_time_sec_sum")
    combo_wall_f = float(combo_wall) if isinstance(combo_wall, (int, float)) else None
    scorch_wall = summary.get("scorch_wall_time_sec_sum")
    scorch_wall_f = (
        float(scorch_wall) if isinstance(scorch_wall, (int, float)) else None
    )
    combo_parallelism = (
        (combo_wall_f / wall_sec_f)
        if combo_wall_f is not None and wall_sec_f is not None and wall_sec_f > 0
        else None
    )
    bench_idle_frac: Optional[float] = None
    bench_summary_path = (
        _repo_root()
        / "chemdb"
        / "bench"
        / "runs"
        / str(run_id)
        / "summary.json"
    )
    if bench_summary_path.exists():
        try:
            payload = json.loads(bench_summary_path.read_text(encoding="utf-8")) or {}
            if isinstance(payload, dict):
                idle = payload.get("IdleFrac")
                if isinstance(idle, (int, float)):
                    bench_idle_frac = float(idle)
        except Exception:
            logger.warning(
                "[run-efficiency.bench] action=summary_read_failed path=%s",
                bench_summary_path,
                exc_info=True,
            )

    manifest_dir, _ = get_manifest_paths(cfg, run_id)
    sacct_metrics = _collect_sacct_cpu_metrics(
        alloc_cpus,
        logger,
        artifact_dir=manifest_dir / "slurm_accounting",
    )
    report = {
        "run_id": str(run_id),
        "generated_at": datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z",
        "alloc_cpus": int(alloc_cpus),
        "wall_time_sec": wall_sec_f,
        "total_combos_scheduled": summary.get("total_combos_scheduled"),
        "total_combos_completed": summary.get("total_combos_completed"),
        "total_combos_failed": summary.get("total_combos_failed"),
        "combo_wall_time_sec_sum": combo_wall_f,
        "combo_parallelism_est": combo_parallelism,
        "scorch_wall_time_sec_sum": scorch_wall_f,
        "scorch_share_of_combo_wall": summary.get("scorch_share_of_combo_wall"),
        "bench_idle_frac": bench_idle_frac,
        "queue_efficiency": (
            (1.0 - bench_idle_frac) if bench_idle_frac is not None else None
        ),
        "slurm": sacct_metrics,
    }
    out_json = manifest_dir / "run_efficiency.json"
    _write_json_atomic(out_json, report)
    logger.info("[run-efficiency.done] run_id=%s output=%s", run_id, out_json)


def _maybe_run_dud_eval(
    cfg: Mapping[str, Any], run_id: str, pdb_files: list[str]
) -> None:
    """
    Best-effort post-run DUD evaluator. Never raises.
    """
    logger = logging.getLogger("dud-eval")
    if not run_id:
        logger.info("[dud-eval.skip] reason=missing_run_id")
        return

    test_mode = str(cfg.get("TEST_MODE_ENABLE", "off")).lower()
    if "dud" not in test_mode and not _manifest_indicates_dud_run(cfg, run_id):
        logger.info("[dud-eval.skip] reason=test_mode_off test_mode=%s", test_mode)
        return
    if _is_no_library_docking(cfg):
        logger.info("[dud-eval.skip] reason=no_library_docking")
        return

    roots = _resolve_hook_roots(cfg)
    out_dir = roots.analysis_root / "analysis" / "dud_eval"
    cmd = [
        sys.executable,
        "-m",
        "analysis.cli.dud_eval",
        "--run-id",
        str(run_id),
        "--docked-root",
        str(roots.docked_root),
        "--post-docked-root",
        str(roots.post_docked_root),
        "--out-dir",
        str(out_dir),
    ]

    if len(pdb_files) == 1:
        pdb_id = _norm_pdb_id(pdb_files[0])
        if pdb_id:
            cmd.extend(["--pdb-id", pdb_id])
            logger.info("[dud-eval.filter] pdb_id=%s", pdb_id)

    logger.info("[dud-eval.run] cmd=%s", shlex.join(cmd))
    env = _with_repo_src_on_pythonpath(dict(os.environ), roots.code_root)
    result = subprocess.run(cmd, cwd=str(roots.code_root), check=False, env=env)
    if result.returncode != 0:
        logger.warning(
            "[dud-eval.fail] run_id=%s returncode=%s", run_id, result.returncode
        )
    else:
        logger.info(
            "[dud-eval.done] run_id=%s returncode=%s out_root=%s docked_root=%s post_docked_root=%s",
            run_id,
            result.returncode,
            str(out_dir),
            str(roots.docked_root),
            str(roots.post_docked_root),
        )


def _manifest_indicates_dud_run(cfg: Mapping[str, Any], run_id: str) -> bool:
    manifest = load_run_manifest(cfg, run_id) or {}
    command = manifest.get("command") if isinstance(manifest, dict) else {}
    if not isinstance(command, Mapping):
        return False
    test_mode = str(command.get("TEST_MODE_ENABLE") or "").strip().lower()
    argv = str(command.get("argv") or "").strip().lower()
    if "dud" in test_mode:
        return True
    try:
        argv_tokens = set(shlex.split(argv))
    except ValueError:
        argv_tokens = set(argv.split())
    return bool(
        argv_tokens & {"-bench", "--bench", "-bench2", "--bench2", "-dude", "--dude"}
    )
def _maybe_run_master_schema_export(cfg: Mapping[str, Any], run_id: str) -> None:
    """
    Best-effort execution of the master schema export. Never raises.
    """
    logger = logging.getLogger("master-schema-export")
    if not run_id:
        logger.info("[master-export.skip] reason=missing_run_id")
        return
    if _is_no_library_docking(cfg):
        logger.info("[master-export.skip] reason=no_library_docking run_id=%s", run_id)
        return

    roots = _resolve_hook_roots(cfg)

    cmd = [
        sys.executable,
        "-m",
        "analysis.cli.master_schema_export",
        "--run-id",
        str(run_id),
        "--repo-root",
        str(roots.analysis_root),
        "--overwrite",
    ]

    logger.info("[master-export.run] cmd=%s", shlex.join(cmd))
    try:
        env = _with_repo_src_on_pythonpath(dict(os.environ), roots.code_root)
        result = subprocess.run(cmd, cwd=str(roots.code_root), check=False, env=env)
        if result.returncode != 0:
            logger.warning(
                "[master-export.fail] run_id=%s returncode=%s",
                run_id,
                result.returncode,
            )
        else:
            logger.info(
                "[master-export.done] run_id=%s returncode=%s",
                run_id,
                result.returncode,
            )
    except Exception:
        logger.warning(
            "[master-export.fail] reason=unexpected_exception", exc_info=True
        )


def _maybe_run_report_generation(cfg: Mapping[str, Any], run_id: str) -> None:
    """
    Best-effort execution of the run report generation. Never raises.
    """
    logger = logging.getLogger("run-report-hook")
    if not run_id:
        logger.info("[report.skip] reason=missing_run_id")
        return
    if _is_no_library_docking(cfg):
        logger.info("[report.skip] reason=no_library_docking run_id=%s", run_id)
        return

    roots = _resolve_hook_roots(cfg)

    cmd = [
        sys.executable,
        "-m",
        "analysis.cli.run_report",
        "--run-id",
        str(run_id),
        "--repo-root",
        str(roots.analysis_root),
        "--emit-html",
        "--vina",
        "--overwrite",
    ]

    logger.info("[report.run] cmd=%s", shlex.join(cmd))
    try:
        env = _with_repo_src_on_pythonpath(dict(os.environ), roots.code_root)
        result = subprocess.run(cmd, cwd=str(roots.code_root), check=False, env=env)
        if result.returncode != 0:
            logger.warning(
                "[report.fail] run_id=%s returncode=%s", run_id, result.returncode
            )
        else:
            logger.info(
                "[report.done] run_id=%s returncode=%s", run_id, result.returncode
            )
            _maybe_run_static_heatmaps(cfg, run_id)
    except Exception:
        logger.warning("[report.fail] reason=unexpected_exception", exc_info=True)


def _maybe_run_static_heatmaps(cfg: Mapping[str, Any], run_id: str) -> None:
    logger = logging.getLogger("heatmap-hook")
    if not run_id:
        logger.info("[heatmap.skip] reason=missing_run_id")
        return

    rscript = shutil.which("Rscript")
    if not rscript:
        logger.info("[heatmap.skip] reason=missing_rscript")
        return

    roots = _resolve_hook_roots(cfg)
    script_path = roots.code_root / "analysis" / "r" / "HeatMap.R"
    if not script_path.exists():
        script_path = roots.code_root / "analysis" / "HeatMap.R"
    if not script_path.exists():
        logger.info("[heatmap.skip] reason=missing_script path=%s", script_path)
        return

    data_dir = roots.analysis_root / "data" / str(run_id)
    inputs = sorted(data_dir.glob("heatmap_input*.csv")) if data_dir.exists() else []
    if not inputs:
        logger.info(
            "[heatmap.skip] reason=missing_input run_id=%s data_dir=%s",
            run_id,
            data_dir,
        )
        return
    logger.info(
        "[heatmap.plan] run_id=%s input_count=%d data_dir=%s",
        run_id,
        len(inputs),
        data_dir,
    )

    for input_path in inputs:
        suffix = input_path.stem[len("heatmap_input") :]
        out_path = data_dir / f"heatmap{suffix}.png"
        cmd = [
            rscript,
            str(script_path),
            "--repo-root",
            str(roots.analysis_root),
            "--run-id",
            str(run_id),
            "--in",
            str(input_path),
            "--out",
            str(out_path),
        ]
        logger.debug("[heatmap.run] cmd=%s", shlex.join(cmd))
        try:
            result = subprocess.run(cmd, cwd=str(roots.code_root), check=False)
        except Exception:
            logger.warning(
                "[heatmap.fail] run_id=%s input=%s reason=unexpected_exception",
                run_id,
                input_path,
                exc_info=True,
            )
            continue
        if result.returncode != 0:
            logger.warning(
                "[heatmap.fail] run_id=%s input=%s returncode=%s",
                run_id,
                input_path,
                result.returncode,
            )
        else:
            logger.info("[heatmap.done] run_id=%s output=%s", run_id, out_path)
