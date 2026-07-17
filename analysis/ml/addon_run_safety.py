"""Concurrency and active-run guards for ML add-on scoring."""

from __future__ import annotations

import fcntl
import json
import os
import socket
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator, Mapping

from analysis.reporting.manifest_utils import load_run_manifest
from cli.status_dashboard_slurm import slurm_job_ids, slurm_status
from src.config.output_paths import output_root


ACTIVE_STATUSES = frozenset({"queued", "running", "submitted"})
ACTIVE_SLURM_STATES = frozenset(
    {
        "CONFIGURING",
        "COMPLETING",
        "PENDING",
        "REQUEUED",
        "RESIZING",
        "RUNNING",
        "STAGE_OUT",
        "SUSPENDED",
    }
)
DEFAULT_ACTIVITY_GRACE_SECONDS = 72.0 * 60.0 * 60.0


@dataclass(frozen=True)
class ActiveAtlasRun:
    run_id: str
    status: str
    manifest: str
    started_at: str | None = None
    last_activity_at: str | None = None
    activity_age_seconds: float | None = None
    reason: str = "recent_active_manifest"
    host: str | None = None


def _run_manifest_paths(root: Path) -> list[Path]:
    roots: list[Path] = []
    relocated = str(os.environ.get("MANIFESTS_DIR") or "").strip()
    if relocated:
        roots.append(Path(relocated).expanduser())
    roots.extend((output_root(root, "manifests"), root / "manifests"))

    by_run_id: dict[str, Path] = {}
    for manifest_root in roots:
        if not manifest_root.is_dir():
            continue
        for path in sorted(manifest_root.glob("*/run_manifest.yaml")):
            by_run_id.setdefault(path.parent.name, path)
    return [by_run_id[run_id] for run_id in sorted(by_run_id)]


def active_atlas_runs(
    repo_root: str | Path,
    *,
    exclude_run_ids: tuple[str, ...] = (),
    activity_grace_seconds: float = DEFAULT_ACTIVITY_GRACE_SECONDS,
) -> list[ActiveAtlasRun]:
    """Return runs with process, Slurm, or recent manifest evidence."""

    root = Path(repo_root).resolve()
    excluded = {str(value).strip() for value in exclude_run_ids if str(value).strip()}
    manifest_paths = _run_manifest_paths(root)
    if not manifest_paths:
        return []

    now = time.time()
    grace_seconds = max(0.0, float(activity_grace_seconds))
    process_commands = _atlas_process_commands()
    active: list[ActiveAtlasRun] = []
    for path in manifest_paths:
        candidate_run_id = path.parent.name
        if candidate_run_id in excluded:
            continue
        payload, loaded_path = load_run_manifest(root, candidate_run_id)
        if not isinstance(payload, Mapping):
            continue
        run_id = str(payload.get("run_id") or candidate_run_id).strip()
        status = str(payload.get("status") or "unknown").strip().casefold()
        if not run_id or run_id in excluded or status not in ACTIVE_STATUSES:
            continue

        timing = payload.get("timing")
        started_at = None
        if isinstance(timing, Mapping):
            value = timing.get("started_at") or timing.get("created_at")
            started_at = str(value).strip() if value else None

        resources = payload.get("resources")
        host = None
        if isinstance(resources, Mapping):
            raw_host = resources.get("host")
            host = str(raw_host).strip() if raw_host else None

        manifest_path = loaded_path or path
        latest_activity = _mtime_or_zero(manifest_path) or now
        activity_age = max(0.0, now - latest_activity)
        reason = _active_reason(
            payload=payload,
            run_id=run_id,
            process_commands=process_commands,
            activity_age_seconds=activity_age,
            activity_grace_seconds=grace_seconds,
        )
        if reason is None:
            continue
        active.append(
            ActiveAtlasRun(
                run_id=run_id,
                status=status,
                manifest=str(manifest_path.resolve()),
                started_at=started_at,
                last_activity_at=_timestamp_iso(latest_activity),
                activity_age_seconds=round(activity_age, 3),
                reason=reason,
                host=host,
            )
        )
    return active


def _active_reason(
    *,
    payload: Mapping[str, Any],
    run_id: str,
    process_commands: tuple[str, ...],
    activity_age_seconds: float,
    activity_grace_seconds: float,
) -> str | None:
    if _process_matches_run(run_id, process_commands):
        return "matching_local_process"
    slurm_reason = _active_slurm_reason(payload)
    if slurm_reason:
        return slurm_reason
    if activity_age_seconds <= max(0.0, activity_grace_seconds):
        return "recent_active_manifest"
    return None


def _active_slurm_reason(payload: Mapping[str, Any]) -> str | None:
    if not slurm_job_ids(payload, include_environment=False):
        return None
    status = slurm_status(payload, include_environment=False)
    jobs = status.get("jobs")
    if not isinstance(jobs, list):
        return None
    states = {
        str(job.get("state") or "").strip().upper().split("+")[0]
        for job in jobs
        if isinstance(job, Mapping)
    }
    active_states = sorted(states & ACTIVE_SLURM_STATES)
    if not active_states:
        return None
    return "active_slurm:" + ",".join(active_states)


def _atlas_process_commands() -> tuple[str, ...]:
    current_pids = {os.getpid(), os.getppid()}
    proc_root = Path("/proc")
    if not proc_root.is_dir():
        return ()
    commands: list[str] = []
    for path in proc_root.iterdir():
        if not path.name.isdigit() or int(path.name) in current_pids:
            continue
        try:
            command = (
                (path / "cmdline")
                .read_bytes()
                .replace(b"\0", b" ")
                .decode("utf-8", errors="replace")
            )
        except OSError:
            continue
        command_folded = command.casefold()
        if any(
            token in command_folded
            for token in (
                "atlas_main_cli",
                "main.py",
                "vina",
                "gnina",
                "score_reference_run_vina",
                "rescoring_scorch",
                "scorch.py",
            )
        ):
            commands.append(command)
    return tuple(commands)


def _process_matches_run(run_id: str, commands: tuple[str, ...]) -> bool:
    markers = (
        f"--run-id {run_id}",
        f"--run-id={run_id}",
        f"/{run_id}/",
    )
    return any(marker in command for command in commands for marker in markers)


def _mtime_or_zero(path: Path) -> float:
    try:
        return float(path.stat().st_mtime)
    except OSError:
        return 0.0


def _timestamp_iso(timestamp: float) -> str | None:
    if timestamp <= 0:
        return None
    return datetime.fromtimestamp(timestamp, timezone.utc).isoformat()


def require_idle_atlas(repo_root: str | Path) -> None:
    """Fail closed when any Atlas run has current execution evidence."""

    active = active_atlas_runs(repo_root)
    if not active:
        return
    summary = ", ".join(
        f"{item.run_id} ({item.status}; {item.reason}; "
        f"activity_age={float(item.activity_age_seconds or 0.0) / 3600.0:.1f}h)"
        for item in active
    )
    raise RuntimeError(
        "refusing to launch ML add-on scoring while another Atlas run is active: "
        f"{summary}. Resolve or finalize those runs with `atlas status <RUN_ID> "
        "--errors --explain` before retrying."
    )


@contextmanager
def exclusive_addon_lock(
    repo_root: str | Path,
    *,
    metadata: Mapping[str, object],
) -> Iterator[Path]:
    """Hold a cross-process advisory lock for one add-on scorer at a time."""

    lock_root = output_root(Path(repo_root).resolve(), "manifests")
    lock_root.mkdir(parents=True, exist_ok=True)
    lock_path = lock_root / ".atlas_ml_addon.lock"
    with lock_path.open("a+", encoding="utf-8") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            owner = handle.read().strip() or "unknown owner"
            raise RuntimeError(
                "another Atlas ML add-on scorer already holds the exclusive lock: "
                f"{owner}"
            ) from exc

        owner_payload = {
            "pid": os.getpid(),
            "host": socket.gethostname(),
            "started_at_utc": datetime.now(timezone.utc).isoformat(),
            **{str(key): value for key, value in metadata.items()},
        }
        handle.seek(0)
        handle.truncate()
        handle.write(json.dumps(owner_payload, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
        try:
            yield lock_path
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def active_runs_as_dicts(runs: list[ActiveAtlasRun]) -> list[dict[str, object]]:
    return [asdict(run) for run in runs]


__all__ = [
    "ACTIVE_STATUSES",
    "ACTIVE_SLURM_STATES",
    "ActiveAtlasRun",
    "DEFAULT_ACTIVITY_GRACE_SECONDS",
    "active_atlas_runs",
    "active_runs_as_dicts",
    "exclusive_addon_lock",
    "require_idle_atlas",
]
