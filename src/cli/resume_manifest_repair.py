from __future__ import annotations

import logging
import os
import time
from pathlib import Path
from typing import Any, Iterable, Mapping

from cli.qol.manifest_repair import ManifestRepairSummary, repair_run_manifest
from config.output_paths import runtime_root


_FALSE_TOKENS = {"0", "false", "no", "off", "disable", "disabled"}
_TRUE_TOKENS = {"1", "true", "yes", "on", "enable", "enabled"}


def _bool_token(raw: Any) -> bool | None:
    token = str(raw or "").strip().lower()
    if not token:
        return None
    if token in _TRUE_TOKENS:
        return True
    if token in _FALSE_TOKENS:
        return False
    return None


def _repair_disabled() -> bool:
    explicit = _bool_token(os.environ.get("ATLAS_REPAIR_BEFORE_RESUME"))
    if explicit is not None:
        return not explicit
    return bool(_bool_token(os.environ.get("ATLAS_NO_REPAIR_BEFORE_RESUME")))


def _repair_scorch_selection_requested(cfg: Mapping[str, Any]) -> bool:
    for raw in (
        os.environ.get("ATLAS_REPAIR_SCORCH_SELECTION"),
        os.environ.get("REPAIR_SCORCH_SELECTION"),
        cfg.get("ATLAS_REPAIR_SCORCH_SELECTION"),
        cfg.get("REPAIR_SCORCH_SELECTION"),
    ):
        parsed = _bool_token(raw)
        if parsed is not None:
            return parsed
    return False


def _dedup_paths(paths: Iterable[Path]) -> list[Path]:
    seen: set[str] = set()
    out: list[Path] = []
    for raw in paths:
        try:
            path = Path(raw).expanduser()
        except Exception:
            continue
        key = str(path)
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(path)
    return out


def _root_from_manifest_root(raw: Any) -> Path | None:
    if not raw:
        return None
    path = Path(str(raw)).expanduser()
    if path.name == "manifests" and path.parent.name == "outputs":
        return path.parent.parent
    return path


def _candidate_all_dirs(cfg: Mapping[str, Any], run_id: str) -> list[Path]:
    candidates: list[Path] = []

    for key in ("ATLAS_ALL_DIRS", "ATLAS_RESUME_REPAIR_ROOT"):
        token = str(os.environ.get(key, "")).strip()
        if token:
            candidates.append(Path(token))

    for key in ("ATLAS_ALL_DIRS", "ALL_DIRS"):
        token = str(cfg.get(key) or "").strip()
        if token:
            candidates.append(Path(token))

    for raw in (
        os.environ.get("MANIFESTS_DIR"),
        os.environ.get("ATLAS_MANIFESTS_DIR"),
        cfg.get("MANIFESTS_DIR"),
    ):
        root = _root_from_manifest_root(raw)
        if root is not None:
            candidates.append(root)

    try:
        root = _root_from_manifest_root(runtime_root(cfg, "MANIFESTS_DIR", "manifests"))
        if root is not None:
            candidates.append(root)
    except Exception:
        pass

    candidates.append(Path.cwd())

    runs_root = str(os.environ.get("ATLAS_RUNS_ROOT", "")).strip()
    if runs_root:
        candidates.append(Path(runs_root) / str(run_id))
    for env_name in ("SCRATCH", "WORK2", "WORK"):
        base = str(os.environ.get(env_name, "")).strip()
        if base:
            candidates.append(Path(base) / "atlas" / "runs" / str(run_id))

    return _dedup_paths(candidates)


def _has_distributed_state(all_dirs: Path, run_id: str) -> bool:
    dist_dir = (
        all_dirs
        / "outputs"
        / "manifests"
        / str(run_id)
        / "distributed"
    )
    return dist_dir.exists() and any(dist_dir.glob("combo_chunks_*.json"))


def _select_all_dirs(cfg: Mapping[str, Any], run_id: str) -> Path | None:
    for candidate in _candidate_all_dirs(cfg, run_id):
        if _has_distributed_state(candidate, run_id):
            return candidate
    return None


def _scorch_outputs_exist(all_dirs: Path, run_id: str) -> bool:
    post_root = all_dirs / "outputs" / "post_docked" / str(run_id)
    if not post_root.exists():
        return False
    return any(post_root.rglob("consensus_reranked_scorch.csv")) or any(
        post_root.rglob("scorch_scores_all.csv")
    )


def _scorch_enabled(cfg: Mapping[str, Any]) -> bool | None:
    for raw in (
        os.environ.get("USE_SCORCH"),
        os.environ.get("ATLAS_USE_SCORCH"),
        cfg.get("USE_SCORCH"),
        cfg.get("SCORCH_ENABLED"),
    ):
        parsed = _bool_token(raw)
        if parsed is not None:
            return parsed
    return None


def _require_scorch_for_repair(
    cfg: Mapping[str, Any],
    *,
    all_dirs: Path,
    run_id: str,
) -> bool | None:
    if _scorch_enabled(cfg) is True:
        return True
    if _scorch_outputs_exist(all_dirs, run_id):
        return True
    return None


def _repair_lock_path(all_dirs: Path, run_id: str) -> Path:
    return (
        all_dirs
        / "outputs"
        / "manifests"
        / str(run_id)
        / "run_manifest.yaml.auto_repair.lock"
    )


def _try_acquire_lock(path: Path, *, stale_after_sec: float) -> bool:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(str(path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(f"pid={os.getpid()} started_at={time.time():.3f}\n")
        return True
    except FileExistsError:
        try:
            age = time.time() - path.stat().st_mtime
        except Exception:
            age = 0.0
        if age > stale_after_sec:
            try:
                path.unlink()
            except Exception:
                return False
            return _try_acquire_lock(path, stale_after_sec=stale_after_sec)
        return False


def _wait_for_lock(path: Path, *, wait_sec: float) -> None:
    deadline = time.time() + max(0.0, wait_sec)
    while path.exists() and time.time() < deadline:
        time.sleep(0.5)


def auto_repair_resume_manifest(
    cfg: Mapping[str, Any],
    run_id: str,
    *,
    logger: logging.Logger | None = None,
) -> ManifestRepairSummary | None:
    log = logger or logging.getLogger("resume.manifest-repair")
    if _repair_disabled():
        log.info("[resume.repair] action=skip reason=disabled run_id=%s", run_id)
        return None

    all_dirs = _select_all_dirs(cfg, run_id)
    if all_dirs is None:
        log.info(
            "[resume.repair] action=skip reason=distributed_state_not_found run_id=%s",
            run_id,
        )
        return None

    lock_path = _repair_lock_path(all_dirs, run_id)
    wait_sec = float(os.environ.get("ATLAS_RESUME_REPAIR_WAIT_SEC") or 180.0)
    stale_sec = float(os.environ.get("ATLAS_RESUME_REPAIR_STALE_SEC") or 900.0)
    got_lock = _try_acquire_lock(lock_path, stale_after_sec=stale_sec)
    if not got_lock:
        log.info(
            "[resume.repair] action=wait run_id=%s all_dirs=%s lock=%s",
            run_id,
            all_dirs,
            lock_path,
        )
        _wait_for_lock(lock_path, wait_sec=wait_sec)
        return None

    try:
        require_scorch = _require_scorch_for_repair(
            cfg,
            all_dirs=all_dirs,
            run_id=run_id,
        )
        repair_scorch_selection = _repair_scorch_selection_requested(cfg)
        summary = repair_run_manifest(
            str(run_id),
            root=all_dirs,
            dry_run=False,
            require_scorch=require_scorch,
            repair_scorch_selection=repair_scorch_selection,
            clear_incomplete_scorch_outputs=bool(require_scorch)
            and not bool(repair_scorch_selection),
            trust_existing_scorch_coverage=bool(require_scorch)
            and not bool(repair_scorch_selection),
        )
        log.info(
            "[resume.repair] action=done run_id=%s all_dirs=%s require_scorch=%s repair_scorch_selection=%s chunks=%d/%d targets=%d completed=%d running=%d failed=%d scorch=%d",
            run_id,
            all_dirs,
            summary.require_scorch,
            summary.repair_scorch_selection,
            summary.completed_chunks,
            summary.planned_chunks,
            summary.targets,
            summary.targets_completed,
            summary.targets_running,
            summary.targets_failed,
            summary.scorch_reranked_files,
        )
        return summary
    except Exception:
        log.warning(
            "[resume.repair] action=error run_id=%s all_dirs=%s",
            run_id,
            all_dirs,
            exc_info=True,
        )
        return None
    finally:
        try:
            lock_path.unlink()
        except Exception:
            pass
