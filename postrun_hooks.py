# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping, Optional

from cli.cli_utils import _norm_pdb_id


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
    if "dud" not in test_mode:
        logger.info("[dud-eval.skip] reason=test_mode_off test_mode=%s", test_mode)
        return
    if cfg.get("NO_LIBRARY_DOCKING"):
        logger.info("[dud-eval.skip] reason=no_library_docking")
        return

    repo_root = Path(__file__).resolve().parent
    preferred = repo_root / "analysis" / "dud_eval.py"
    fallback = repo_root / "dud_eval.py"
    dud_eval_path = preferred if preferred.exists() else fallback
    if not dud_eval_path.exists():
        logger.warning("[dud-eval.skip] reason=script_missing path=%s", dud_eval_path)
        return

    cmd = [
        sys.executable,
        str(dud_eval_path),
        "--run-id",
        str(run_id),
        "--docked-root",
        "docked",
        "--out-dir",
        "analysis/dud_eval",
    ]

    if len(pdb_files) == 1:
        pdb_id = _norm_pdb_id(pdb_files[0])
        if pdb_id:
            cmd.extend(["--pdb-id", pdb_id])
            logger.info("[dud-eval.filter] pdb_id=%s", pdb_id)

    logger.info("[dud-eval.run] cmd=%s", shlex.join(cmd))
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        logger.warning(
            "[dud-eval.fail] run_id=%s returncode=%s", run_id, result.returncode
        )
    else:
        logger.info(
            "[dud-eval.done] run_id=%s returncode=%s out_root=%s",
            run_id,
            result.returncode,
            "analysis/dud_eval",
        )


def _maybe_run_scorch_rescore(
    cfg: Mapping[str, Any], run_id: str, verbose: bool = False
) -> None:
    """
    Best-effort post-run SCORCH rescoring. Never raises.
    """
    logger = logging.getLogger("scorch-rescore-hook")
    if not run_id:
        logger.info("[scorch-rescore.skip] reason=missing_run_id")
        return

    repo_root = Path(__file__).resolve().parent
    script_path = repo_root / "src/post_docking/rescoring/rescoring_scorch.py"
    if not script_path.exists():
        logger.warning(
            "[scorch-rescore.skip] reason=missing_script path=%s", script_path
        )
        return

    def _as_int(val: Any, fallback: int) -> int:
        try:
            return int(val)
        except Exception:
            return fallback

    total_cpu = _as_int(cfg.get("CPU") or (os.cpu_count() or 1), os.cpu_count() or 1)
    total_cpu = max(1, total_cpu)

    jobs = total_cpu
    threads = max(1, total_cpu // jobs)

    scorch_jobs = cfg.get("SCORCH_JOBS")
    scorch_threads = cfg.get("SCORCH_THREADS")
    if scorch_jobs is not None:
        jobs = max(1, min(_as_int(scorch_jobs, jobs), total_cpu))
        threads = max(1, total_cpu // jobs)
    if scorch_threads is not None:
        threads = max(1, min(_as_int(scorch_threads, threads), total_cpu))
        if jobs * threads > total_cpu:
            threads = max(1, total_cpu // jobs)

    cmd = [
        sys.executable,
        str(script_path),
        "--run-id",
        run_id,
        "--repo-root",
        str(repo_root),
        "--threads",
        str(threads),
        "--jobs",
        str(jobs),
    ]
    if verbose or logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        cmd.append("--verbose")

    logger.info(
        "[scorch-rescore.invoke] run_id=%s jobs=%d threads=%d total_cpu=%d cmd=%s",
        run_id,
        jobs,
        threads,
        total_cpu,
        shlex.join(str(c) for c in cmd),
    )
    try:
        proc = subprocess.run(cmd, cwd=str(repo_root), check=False)
    except Exception:
        logger.warning(
            "[scorch-rescore.invoke] action=skip reason=execution_failed", exc_info=True
        )
        return

    if proc.returncode != 0:
        logger.warning(
            "[scorch-rescore.result] status=failed run_id=%s returncode=%d",
            run_id,
            proc.returncode,
        )
    else:
        logger.info("[scorch-rescore.result] status=ok run_id=%s", run_id)


def _maybe_run_scorch_rescore_for_pdb(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    *,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
    verbose: bool = False,
) -> None:
    """
    Best-effort per-protein SCORCH rescoring. Never raises.
    """
    logger = logging.getLogger("scorch-rescore-hook")
    if not run_id:
        logger.info("[scorch-rescore.skip] reason=missing_run_id")
        return
    if not pdb_id:
        logger.info("[scorch-rescore.skip] reason=missing_pdb_id run_id=%s", run_id)
        return

    raw_use = cfg.get("USE_SCORCH", False)
    if isinstance(raw_use, bool):
        use_scorch = raw_use
    else:
        use_scorch = str(raw_use).strip().lower() in {"1", "true", "yes", "on"}
    if not use_scorch:
        logger.info(
            "[scorch-rescore.skip] reason=use_scorch_false run_id=%s pdb_id=%s",
            run_id,
            pdb_id,
        )
        return

    repo_root = Path(__file__).resolve().parent
    script_path = repo_root / "src/post_docking/rescoring/rescoring_scorch.py"
    if not script_path.exists():
        logger.warning(
            "[scorch-rescore.skip] reason=missing_script path=%s", script_path
        )
        return

    def _as_int(val: Any, fallback: int) -> int:
        try:
            return int(val)
        except Exception:
            return fallback

    total_cpu = max(
        1,
        _as_int(cfg.get("CPU") or (os.cpu_count() or 1), os.cpu_count() or 1),
    )
    jobs = 1
    if cfg.get("SCORCH_JOBS") is not None:
        jobs = max(1, min(_as_int(cfg.get("SCORCH_JOBS"), jobs), total_cpu))

    threads_cfg = cfg.get("SCORCH_THREADS")
    threads = total_cpu if threads_cfg is None else _as_int(threads_cfg, total_cpu)
    threads = max(1, min(threads, total_cpu))
    if jobs * threads > total_cpu:
        threads = max(1, total_cpu // jobs)

    cmd = [
        sys.executable,
        str(script_path),
        "--run-id",
        str(run_id),
        "--repo-root",
        str(repo_root),
        "--pdb-id",
        str(pdb_id),
        "--threads",
        str(threads),
        "--jobs",
        str(jobs),
    ]
    if variant:
        cmd.extend(["--variant", str(variant)])
    if ph:
        cmd.extend(["--ph", str(ph)])
    if verbose or logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        cmd.append("--verbose")

    logger.info(
        "[scorch-rescore.invoke] scope=per_pdb run_id=%s pdb_id=%s variant=%s ph=%s jobs=%d threads=%d total_cpu=%d cmd=%s",
        run_id,
        pdb_id,
        variant or "*",
        ph or "*",
        jobs,
        threads,
        total_cpu,
        shlex.join(str(c) for c in cmd),
    )
    try:
        proc = subprocess.run(cmd, cwd=str(repo_root), check=False)
    except Exception:
        logger.warning(
            "[scorch-rescore.invoke] action=skip reason=execution_failed run_id=%s pdb_id=%s",
            run_id,
            pdb_id,
            exc_info=True,
        )
        return

    logger.info(
        "[scorch-rescore.result] scope=per_pdb run_id=%s pdb_id=%s variant=%s ph=%s returncode=%d",
        run_id,
        pdb_id,
        variant or "*",
        ph or "*",
        proc.returncode,
    )


def _log_rescore_verification(run_id: str) -> None:
    """
    Lightweight best-effort check for consensus and SCORCH outputs.
    """
    logger = logging.getLogger("post-run-checks")
    if not run_id:
        logger.info("[post-check.skip] reason=missing_run_id")
        return

    repo_root = Path(__file__).resolve().parent
    docked_root = repo_root / "docked" / run_id
    post_root = repo_root / "post_docked" / run_id

    consensus_paths = [
        p for p in docked_root.rglob("consensus_docking_scores.csv") if p.is_file()
    ]
    scorch_paths = [p for p in post_root.rglob("scorch_scores_all.csv") if p.is_file()]

    if consensus_paths:
        logger.info(
            "[post-check.consensus] found=%d sample=%s",
            len(consensus_paths),
            consensus_paths[0],
        )
    else:
        logger.warning("[post-check.consensus] reason=missing_files run_id=%s", run_id)

    if scorch_paths:
        logger.info(
            "[post-check.scorch] found=%d sample=%s", len(scorch_paths), scorch_paths[0]
        )
    else:
        logger.warning("[post-check.scorch] reason=missing_files run_id=%s", run_id)


def _maybe_run_master_schema_export(cfg: Mapping[str, Any], run_id: str) -> None:
    """
    Best-effort execution of the master schema export. Never raises.
    """
    logger = logging.getLogger("master-schema-export")
    if not run_id:
        logger.info("[master-export.skip] reason=missing_run_id")
        return

    repo_root = Path(__file__).resolve().parent
    script_path = repo_root / "analysis" / "master_schema_export.py"
    if not script_path.exists():
        logger.warning(
            "[master-export.skip] reason=missing_script path=%s", script_path
        )
        return

    cmd = [
        sys.executable,
        "-m",
        "analysis.master_schema_export",
        "--run-id",
        str(run_id),
        "--repo-root",
        str(repo_root),
        "--overwrite",
    ]

    logger.info("[master-export.run] cmd=%s", shlex.join(cmd))
    try:
        result = subprocess.run(cmd, cwd=str(repo_root), check=False)
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

    repo_root = Path(__file__).resolve().parent
    script_path = repo_root / "analysis" / "run_report.py"
    if not script_path.exists():
        logger.warning("[report.skip] reason=missing_script path=%s", script_path)
        return

    cmd = [
        sys.executable,
        "-m",
        "analysis.run_report",
        "--run-id",
        str(run_id),
        "--repo-root",
        str(repo_root),
        "--overwrite",
    ]

    logger.info("[report.run] cmd=%s", shlex.join(cmd))
    try:
        result = subprocess.run(cmd, cwd=str(repo_root), check=False)
        if result.returncode != 0:
            logger.warning(
                "[report.fail] run_id=%s returncode=%s", run_id, result.returncode
            )
        else:
            logger.info(
                "[report.done] run_id=%s returncode=%s", run_id, result.returncode
            )
    except Exception:
        logger.warning("[report.fail] reason=unexpected_exception", exc_info=True)


def _normalize_artifact_retention_mode(raw: Any) -> str:
    token = str(raw or "").strip().lower().replace("-", "_")
    if token in {"", "rerunsafe", "rerun_safe"}:
        return "rerun_safe"
    if token in {"minimaldisk", "minimal_disk"}:
        return "minimal_disk"
    if token in {"off", "none", "false", "0", "disabled", "disable"}:
        return "off"
    return "rerun_safe"


def _as_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    token = str(value).strip().lower()
    if token in {"1", "true", "yes", "on", "y"}:
        return True
    if token in {"0", "false", "no", "off", "n"}:
        return False
    return default


def _to_globs(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(x).strip() for x in value if str(x).strip()]
    token = str(value).strip()
    if not token:
        return []
    if "," in token:
        return [x.strip() for x in token.split(",") if x.strip()]
    return [token]


def _invoke_artifact_retention(
    cfg: Mapping[str, Any],
    run_id: str,
    *,
    pdb_id: Optional[str],
    lock: Optional[Any] = None,
) -> None:
    logger = logging.getLogger("artifact-retention-hook")
    if not run_id:
        logger.info("[artifact-retention.skip] reason=missing_run_id")
        return
    mode = _normalize_artifact_retention_mode(cfg.get("ARTIFACT_RETENTION"))
    if mode == "off":
        logger.info(
            "[artifact-retention.skip] reason=mode_off run_id=%s pdb_id=%s",
            run_id,
            pdb_id or "*",
        )
        return

    dry_run = _as_bool(cfg.get("ARTIFACT_RETENTION_DRY_RUN"), False)
    overwrite = _as_bool(cfg.get("ARTIFACT_RETENTION_OVERWRITE"), False)
    verify_only = _as_bool(cfg.get("ARTIFACT_RETENTION_VERIFY_ONLY"), False)
    include_globs = _to_globs(cfg.get("ARTIFACT_RETENTION_INCLUDE_GLOB"))
    exclude_globs = _to_globs(cfg.get("ARTIFACT_RETENTION_EXCLUDE_GLOB"))

    threads_raw = cfg.get("ARTIFACT_RETENTION_THREADS")
    try:
        threads = (
            max(1, int(threads_raw))
            if threads_raw is not None
            else max(1, os.cpu_count() or 1)
        )
    except Exception:
        threads = max(1, os.cpu_count() or 1)

    repo_root = Path(__file__).resolve().parent
    cmd = [
        sys.executable,
        "-m",
        "post_docking.artifact_retention",
        "--run-id",
        str(run_id),
        "--mode",
        mode,
        "--threads",
        str(threads),
    ]
    if pdb_id:
        cmd.extend(["--pdb-id", str(pdb_id)])
    if dry_run:
        cmd.append("--dry-run")
    if overwrite:
        cmd.append("--overwrite")
    if verify_only:
        cmd.append("--verify-only")
    for pat in include_globs:
        cmd.extend(["--include-glob", pat])
    for pat in exclude_globs:
        cmd.extend(["--exclude-glob", pat])

    logger.info(
        "[artifact-retention.invoke] run_id=%s pdb_id=%s mode=%s dry_run=%s overwrite=%s verify_only=%s cmd=%s",
        run_id,
        pdb_id or "*",
        mode,
        str(dry_run).lower(),
        str(overwrite).lower(),
        str(verify_only).lower(),
        shlex.join(str(c) for c in cmd),
    )

    def _run_cmd() -> None:
        try:
            result = subprocess.run(cmd, cwd=str(repo_root), check=False)
            if result.returncode != 0:
                logger.warning(
                    "[artifact-retention.fail] run_id=%s pdb_id=%s returncode=%s",
                    run_id,
                    pdb_id or "*",
                    result.returncode,
                )
            else:
                logger.info(
                    "[artifact-retention.done] run_id=%s pdb_id=%s returncode=%s",
                    run_id,
                    pdb_id or "*",
                    result.returncode,
                )
        except Exception:
            logger.warning(
                "[artifact-retention.fail] reason=unexpected_exception run_id=%s pdb_id=%s",
                run_id,
                pdb_id or "*",
                exc_info=True,
            )

    if lock is not None:
        with lock:
            _run_cmd()
    else:
        _run_cmd()


def _maybe_run_artifact_retention_for_pdb(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    *,
    lock: Optional[Any] = None,
) -> None:
    _invoke_artifact_retention(cfg, run_id, pdb_id=pdb_id, lock=lock)


def _maybe_run_artifact_retention(cfg: Mapping[str, Any], run_id: str) -> None:
    """
    Optional run-level retention entrypoint (idempotent catch-up path).
    """
    _invoke_artifact_retention(cfg, run_id, pdb_id=None, lock=None)
