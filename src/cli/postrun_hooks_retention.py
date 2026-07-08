# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
import os
import shlex
import subprocess
import sys
from typing import Any, Mapping, Optional

from cli.postrun_hooks_support import _repo_root, _with_repo_src_on_pythonpath


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
    variant: Optional[str] = None,
    ph: Optional[str] = None,
    combos: Optional[list[tuple[str, str, str]]] = None,
    lock: Optional[Any] = None,
) -> None:
    logger = logging.getLogger("artifact-retention-hook")
    if not run_id:
        logger.info("[artifact-retention.skip] reason=missing_run_id")
        return
    mode = _normalize_artifact_retention_mode(cfg.get("ARTIFACT_RETENTION"))
    if mode == "off":
        logger.info(
            "[artifact-retention.skip] reason=mode_off run_id=%s pdb_id=%s variant=%s ph=%s",
            run_id,
            pdb_id or "*",
            variant or "*",
            ph or "*",
        )
        return

    dry_run = _as_bool(cfg.get("ARTIFACT_RETENTION_DRY_RUN"), False)
    overwrite = _as_bool(cfg.get("ARTIFACT_RETENTION_OVERWRITE"), False)
    verify_only = _as_bool(cfg.get("ARTIFACT_RETENTION_VERIFY_ONLY"), False)
    include_globs = _to_globs(cfg.get("ARTIFACT_RETENTION_INCLUDE_GLOB"))
    exclude_globs = _to_globs(cfg.get("ARTIFACT_RETENTION_EXCLUDE_GLOB"))

    threads_raw = cfg.get("CPU")
    try:
        threads = max(1, int(threads_raw)) if threads_raw is not None else 1
    except Exception:
        threads = 1
    group_workers_raw = cfg.get("ARTIFACT_RETENTION_GROUP_WORKERS", 1)
    try:
        group_workers = max(1, int(group_workers_raw))
    except Exception:
        group_workers = 1

    repo_root = _repo_root()
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
        "--group-workers",
        str(group_workers),
    ]
    docked_root = cfg.get("DOCKED_DIR")
    if docked_root:
        cmd.extend(["--docked-root", str(docked_root)])
    post_docked_root = cfg.get("POST_DOCKED_DIR")
    if post_docked_root:
        cmd.extend(["--post-docked-root", str(post_docked_root)])
    combo_scope: list[tuple[str, str, str]] = []
    if combos:
        deduped = {
            (
                str(p or "").strip().upper(),
                str(v or "LEGACY").strip().upper() or "LEGACY",
                str(ph_t or "none").strip() or "none",
            )
            for p, v, ph_t in combos
            if str(p or "").strip()
        }
        combo_scope = sorted(deduped)
        for combo_pdb, combo_variant, combo_ph in combo_scope:
            cmd.extend(
                [
                    "--combo",
                    f"{combo_pdb}:{combo_variant}:{combo_ph}",
                ]
            )
    else:
        if pdb_id:
            cmd.extend(["--pdb-id", str(pdb_id)])
        if variant:
            cmd.extend(["--variant", str(variant)])
        if ph:
            cmd.extend(["--ph", str(ph)])
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
        "[artifact-retention.invoke] run_id=%s pdb_id=%s variant=%s ph=%s combos=%d mode=%s dry_run=%s overwrite=%s verify_only=%s group_workers=%d cmd=%s",
        run_id,
        pdb_id or "*",
        variant or "*",
        ph or "*",
        len(combo_scope),
        mode,
        str(dry_run).lower(),
        str(overwrite).lower(),
        str(verify_only).lower(),
        group_workers,
        shlex.join(str(c) for c in cmd),
    )

    def _run_cmd() -> None:
        try:
            env = _with_repo_src_on_pythonpath(dict(os.environ), repo_root)
            result = subprocess.run(cmd, cwd=str(repo_root), check=False, env=env)
            if result.returncode != 0:
                logger.warning(
                    "[artifact-retention.fail] run_id=%s pdb_id=%s variant=%s ph=%s combos=%d returncode=%s",
                    run_id,
                    pdb_id or "*",
                    variant or "*",
                    ph or "*",
                    len(combo_scope),
                    result.returncode,
                )
            else:
                logger.info(
                    "[artifact-retention.done] run_id=%s pdb_id=%s variant=%s ph=%s combos=%d returncode=%s",
                    run_id,
                    pdb_id or "*",
                    variant or "*",
                    ph or "*",
                    len(combo_scope),
                    result.returncode,
                )
        except Exception:
            logger.warning(
                "[artifact-retention.fail] reason=unexpected_exception run_id=%s pdb_id=%s variant=%s ph=%s combos=%d",
                run_id,
                pdb_id or "*",
                variant or "*",
                ph or "*",
                len(combo_scope),
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


def _maybe_run_artifact_retention_for_combo(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant: Optional[str],
    ph: Optional[str],
    *,
    lock: Optional[Any] = None,
) -> None:
    _invoke_artifact_retention(
        cfg,
        run_id,
        pdb_id=pdb_id,
        variant=variant,
        ph=ph,
        lock=lock,
    )


def _maybe_run_artifact_retention_for_combos(
    cfg: Mapping[str, Any],
    run_id: str,
    combos: list[tuple[str, str, str]],
    *,
    lock: Optional[Any] = None,
) -> None:
    _invoke_artifact_retention(
        cfg,
        run_id,
        pdb_id=None,
        variant=None,
        ph=None,
        combos=combos,
        lock=lock,
    )


def _maybe_run_artifact_retention(cfg: Mapping[str, Any], run_id: str) -> None:
    """
    Optional run-level retention entrypoint (idempotent catch-up path).
    """
    _invoke_artifact_retention(cfg, run_id, pdb_id=None, lock=None)
