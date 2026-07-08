from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from config.output_paths import runtime_root


_INT_RE = re.compile(r"(\d+)")


def _parse_int_env(*names: str) -> int | None:
    for name in names:
        raw = os.environ.get(name)
        if raw is None:
            continue
        match = _INT_RE.search(str(raw))
        if not match:
            continue
        try:
            value = int(match.group(1))
        except Exception:
            continue
        if value > 0:
            return value
    return None


def _parse_float_env(name: str, default: float) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return float(default)
    try:
        value = float(raw)
    except Exception:
        return float(default)
    if value <= 0:
        return float(default)
    return float(value)


def _truthy_env(name: str) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return False
    return str(raw).strip().lower() in {"1", "true", "yes", "on"}


def detect_allocated_cpus(fallback: int) -> int:
    """
    Detect scheduler CPU allocation, falling back to the configured CPU value.
    """
    fallback_cpu = max(1, int(fallback or 1))
    if _truthy_env("ATLAS_USE_NODE_CPUS"):
        scheduler_values = (
            _parse_int_env("SLURM_CPUS_ON_NODE"),
            _parse_int_env("SLURM_JOB_CPUS_PER_NODE"),
        )
        scheduler_cpus = [value for value in scheduler_values if value is not None]
        if scheduler_cpus:
            return max(scheduler_cpus)
        candidates = (
            _parse_int_env("SLURM_CPUS_PER_TASK"),
            _parse_int_env("NSLOTS", "PBS_NP"),
            _parse_int_env("OMP_NUM_THREADS"),
        )
    else:
        candidates = (
            _parse_int_env("SLURM_CPUS_PER_TASK"),
            _parse_int_env("SLURM_CPUS_ON_NODE"),
            _parse_int_env("SLURM_JOB_CPUS_PER_NODE"),
            _parse_int_env("NSLOTS", "PBS_NP"),
            _parse_int_env("OMP_NUM_THREADS"),
        )
    for value in candidates:
        if value is not None and value > 0:
            return int(value)
    return fallback_cpu


def _normalize_mode(raw: Any) -> str:
    token = str(raw or "").strip().lower().replace("-", "_")
    if token in {"", "off", "none", "false", "0", "disabled"}:
        return "off"
    if token in {"slurm", "array", "slurm_array"}:
        return "slurm_array"
    return token


def _manifest_dir(cfg: Mapping[str, Any], run_id: str) -> Path:
    return runtime_root(cfg, "MANIFESTS_DIR", "manifests") / str(run_id)


def _distributed_dir(cfg: Mapping[str, Any], run_id: str) -> Path:
    return _manifest_dir(cfg, run_id) / "distributed"


@dataclass(frozen=True)
class DistributedRunContext:
    mode: str
    enabled: bool
    run_id: str
    task_id: int
    task_count: int
    task_min_id: int
    leader_task_id: int
    barrier_timeout_sec: float
    barrier_poll_sec: float

    @property
    def is_leader(self) -> bool:
        return bool(self.enabled and self.task_id == self.leader_task_id)

    @property
    def task_index(self) -> int:
        return max(0, int(self.task_id) - int(self.task_min_id))

    def expected_task_ids(self) -> list[int]:
        if not self.enabled:
            return [self.task_id]
        start = int(self.task_min_id)
        stop = start + int(self.task_count)
        return list(range(start, stop))

    def marker_dir(self, cfg: Mapping[str, Any]) -> Path:
        return _distributed_dir(cfg, self.run_id) / "markers"

    def marker_path(self, cfg: Mapping[str, Any], task_id: int | None = None) -> Path:
        tid = self.task_id if task_id is None else int(task_id)
        return self.marker_dir(cfg) / f"task_{tid}.json"

    def shard_matches(self, token: str) -> bool:
        if not self.enabled:
            return True
        digest = hashlib.sha1(token.encode("utf-8")).hexdigest()
        bucket = int(digest, 16) % max(1, int(self.task_count))
        return bucket == self.task_index


def resolve_distributed_context(
    cfg: Mapping[str, Any],
    run_id: str,
) -> DistributedRunContext:
    env_mode = os.environ.get("ATLAS_DISTRIBUTED_MODE")
    cfg_mode = cfg.get("DISTRIBUTED_MODE")
    mode = _normalize_mode(env_mode if env_mode is not None else cfg_mode)
    env_mode_explicit_off = env_mode is not None and _normalize_mode(env_mode) == "off"

    slurm_task_count = _parse_int_env("SLURM_ARRAY_TASK_COUNT", "ATLAS_DIST_TASK_COUNT") or 1
    slurm_task_id_raw = os.environ.get("SLURM_ARRAY_TASK_ID")
    atlas_task_id_raw = os.environ.get("ATLAS_DIST_TASK_ID")
    has_task_id = bool(
        (slurm_task_id_raw is not None and str(slurm_task_id_raw).strip() != "")
        or (atlas_task_id_raw is not None and str(atlas_task_id_raw).strip() != "")
    )

    # Resume/relaunch wrappers can drop ATLAS_DISTRIBUTED_MODE while preserving
    # SLURM task metadata. In that case, infer slurm_array unless explicitly off.
    if (
        mode == "off"
        and not env_mode_explicit_off
        and int(slurm_task_count) > 1
        and has_task_id
    ):
        mode = "slurm_array"

    if mode == "off":
        return DistributedRunContext(
            mode="off",
            enabled=False,
            run_id=str(run_id),
            task_id=0,
            task_count=1,
            task_min_id=0,
            leader_task_id=0,
            barrier_timeout_sec=0.0,
            barrier_poll_sec=0.0,
        )

    task_count = int(slurm_task_count)
    task_id = _parse_int_env("SLURM_ARRAY_TASK_ID", "ATLAS_DIST_TASK_ID") or 0
    task_min_id = _parse_int_env("SLURM_ARRAY_TASK_MIN", "ATLAS_DIST_TASK_MIN")
    if task_min_id is None:
        task_min_id = 0
    leader_task_id = (
        _parse_int_env("ATLAS_DIST_LEADER_TASK_ID")
        if os.environ.get("ATLAS_DIST_LEADER_TASK_ID")
        else task_min_id
    )
    if leader_task_id is None:
        leader_task_id = task_min_id

    enabled = bool(task_count > 1)
    if not enabled:
        mode = "off"

    timeout_sec = _parse_float_env("ATLAS_DIST_BARRIER_TIMEOUT_SEC", 10800.0)
    poll_sec = _parse_float_env("ATLAS_DIST_BARRIER_POLL_SEC", 10.0)

    return DistributedRunContext(
        mode=mode,
        enabled=enabled,
        run_id=str(run_id),
        task_id=int(task_id),
        task_count=max(1, int(task_count)),
        task_min_id=int(task_min_id),
        leader_task_id=int(leader_task_id),
        barrier_timeout_sec=float(timeout_sec),
        barrier_poll_sec=float(poll_sec),
    )


def shard_pdb_files_for_context(
    ctx: DistributedRunContext,
    pdb_files: Sequence[str],
) -> list[str]:
    if not ctx.enabled:
        return list(pdb_files)

    ordered = [str(pdb_file) for pdb_file in pdb_files]
    task_ids = ctx.expected_task_ids()
    if task_ids and len(ordered) <= max(1, len(task_ids) * 2):
        try:
            task_pos = task_ids.index(int(ctx.task_id))
        except ValueError:
            task_pos = int(ctx.task_index)
        return [
            pdb_file
            for idx, pdb_file in enumerate(ordered)
            if idx % max(1, len(task_ids)) == task_pos
        ]

    selected: list[str] = []
    for pdb_file in ordered:
        stem = Path(str(pdb_file)).stem.replace("_cleaned", "")
        token = stem.strip().upper() or str(pdb_file).strip().upper()
        if ctx.shard_matches(token):
            selected.append(str(pdb_file))
    return selected


def chunk_plan_path(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    variant_label: str,
) -> Path:
    token = str(variant_label or "legacy").strip().lower().replace(" ", "_")
    return _distributed_dir(cfg, ctx.run_id) / f"combo_chunks_{token}.json"


def chunk_claim_dir(ctx: DistributedRunContext, cfg: Mapping[str, Any]) -> Path:
    return _distributed_dir(cfg, ctx.run_id) / "chunk_claims"


def chunk_result_dir(ctx: DistributedRunContext, cfg: Mapping[str, Any]) -> Path:
    return _distributed_dir(cfg, ctx.run_id) / "chunk_results"


def scorch_scope_claim_dir(ctx: DistributedRunContext, cfg: Mapping[str, Any]) -> Path:
    return _distributed_dir(cfg, ctx.run_id) / "scorch_scope_claims"


def scorch_scope_state_dir(ctx: DistributedRunContext, cfg: Mapping[str, Any]) -> Path:
    return _distributed_dir(cfg, ctx.run_id) / "scorch_scope_state"


def combo_prep_state_dir(ctx: DistributedRunContext, cfg: Mapping[str, Any]) -> Path:
    return _distributed_dir(cfg, ctx.run_id) / "combo_prep_state"


def combo_prep_claim_dir(ctx: DistributedRunContext, cfg: Mapping[str, Any]) -> Path:
    return _distributed_dir(cfg, ctx.run_id) / "combo_prep_claims"


def combo_prep_state_path(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    combo_key: str,
) -> Path:
    return combo_prep_state_dir(ctx, cfg) / f"{str(combo_key)}.json"


def combo_prep_claim_path(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    combo_key: str,
) -> Path:
    return combo_prep_claim_dir(ctx, cfg) / f"{str(combo_key)}.claim.json"


def chunk_result_path(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    chunk_id: str,
) -> Path:
    return chunk_result_dir(ctx, cfg) / f"{str(chunk_id)}.json"


def chunk_claim_path(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    chunk_id: str,
) -> Path:
    return chunk_claim_dir(ctx, cfg) / f"{str(chunk_id)}.claim.json"


def _safe_scope_token(raw: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(raw or "").strip())
    return token.strip("._") or "scope"


def scorch_scope_key(
    *,
    run_id: str,
    pdb_id: str,
    variant: str | None,
    ph: str | None,
    phase: str = "final",
) -> str:
    parts = [
        str(run_id or "").strip(),
        str(phase or "final").strip().lower() or "final",
        str(pdb_id or "").strip().upper(),
        str(variant or "*").strip().upper() or "*",
        str(ph or "base").strip() or "base",
    ]
    digest = hashlib.sha1("|".join(parts).encode("utf-8")).hexdigest()[:16]
    readable = "_".join(_safe_scope_token(part) for part in parts[1:])
    return f"{readable}_{digest}"


def scorch_scope_claim_path(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    scope_key: str,
) -> Path:
    return scorch_scope_claim_dir(ctx, cfg) / f"{str(scope_key)}.claim.json"


def scorch_scope_state_path(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    scope_key: str,
) -> Path:
    return scorch_scope_state_dir(ctx, cfg) / f"{str(scope_key)}.json"


def write_json_atomic(path: Path, payload: Mapping[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=path.name + ".",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            tmp_path = Path(handle.name)
            json.dump(payload, handle, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        tmp_path.replace(path)
        return path
    finally:
        try:
            if tmp_path and tmp_path.exists():
                tmp_path.unlink()
        except Exception:
            pass


def read_scorch_scope_state(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    scope_key: str,
) -> dict[str, Any] | None:
    path = scorch_scope_state_path(ctx, cfg, scope_key=scope_key)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def try_claim_scorch_scope(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    scope_key: str,
    pdb_id: str,
    variant: str | None,
    ph: str | None,
    phase: str = "final",
    lease_sec: float = 21600.0,
    force: bool = False,
) -> bool:
    if not ctx.enabled:
        return True

    claim_path = scorch_scope_claim_path(ctx, cfg, scope_key=scope_key)
    claim_path.parent.mkdir(parents=True, exist_ok=True)
    now = time.time()
    lease_value = float(max(30.0, lease_sec))
    payload = {
        "run_id": ctx.run_id,
        "scope_key": str(scope_key),
        "pdb_id": str(pdb_id or "").upper(),
        "variant": str(variant or "*").upper(),
        "ph": str(ph or "base"),
        "phase": str(phase or "final"),
        "task_id": int(ctx.task_id),
        "task_count": int(ctx.task_count),
        "claimed_at": float(now),
        "updated_at": float(now),
        "lease_sec": float(lease_value),
    }

    def _attempt_create() -> bool:
        try:
            fd = os.open(str(claim_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            return False
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, sort_keys=True)
                handle.flush()
                os.fsync(handle.fileno())
            return True
        except Exception:
            try:
                claim_path.unlink()
            except Exception:
                pass
            return False

    if _attempt_create():
        write_json_atomic(
            scorch_scope_state_path(ctx, cfg, scope_key=scope_key),
            {
                **payload,
                "status": "running",
            },
        )
        return True

    try:
        if claim_path.exists():
            age = max(0.0, now - float(claim_path.stat().st_mtime))
            if age > lease_value:
                try:
                    claim_path.unlink()
                except Exception:
                    pass
                if _attempt_create():
                    write_json_atomic(
                        scorch_scope_state_path(ctx, cfg, scope_key=scope_key),
                        {
                            **payload,
                            "status": "running",
                            "reclaimed_stale": True,
                        },
                    )
                    return True
    except Exception:
        pass
    return False


def renew_scorch_scope_claim(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    scope_key: str,
    lease_sec: float | None = None,
) -> bool:
    if not ctx.enabled:
        return True
    claim_path = scorch_scope_claim_path(ctx, cfg, scope_key=scope_key)
    if not claim_path.exists():
        return False
    now = float(time.time())
    try:
        with claim_path.open("r+", encoding="utf-8") as handle:
            try:
                payload = json.load(handle) or {}
            except Exception:
                payload = {}
            if not isinstance(payload, dict):
                payload = {}
            current_lease = payload.get("lease_sec")
            if lease_sec is None:
                try:
                    lease_value = float(str(current_lease))
                except Exception:
                    lease_value = 21600.0
            else:
                lease_value = float(lease_sec)
            payload["updated_at"] = now
            payload["lease_sec"] = float(max(30.0, lease_value))
            handle.seek(0)
            handle.truncate(0)
            json.dump(payload, handle, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        os.utime(claim_path, (now, now))
        return True
    except FileNotFoundError:
        return False
    except Exception:
        return False


def write_scorch_scope_state(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    scope_key: str,
    status: str,
    payload: Mapping[str, Any] | None = None,
    release_claim: bool = True,
) -> Path:
    now = float(time.time())
    state_payload: dict[str, Any] = {
        "run_id": ctx.run_id,
        "scope_key": str(scope_key),
        "task_id": int(ctx.task_id),
        "task_count": int(ctx.task_count),
        "status": str(status),
        "updated_at": now,
    }
    if payload:
        for key, value in payload.items():
            if key in state_payload:
                continue
            state_payload[str(key)] = value

    path = scorch_scope_state_path(ctx, cfg, scope_key=scope_key)
    written = write_json_atomic(path, state_payload)
    if release_claim:
        try:
            claim = scorch_scope_claim_path(ctx, cfg, scope_key=scope_key)
            if claim.exists():
                claim.unlink()
        except Exception:
            pass
    return written


def chunk_result_exists(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    chunk_id: str,
) -> bool:
    return chunk_result_path(ctx, cfg, chunk_id=chunk_id).exists()


def read_chunk_result(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    chunk_id: str,
) -> dict[str, Any] | None:
    path = chunk_result_path(ctx, cfg, chunk_id=chunk_id)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def read_combo_prep_state(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    combo_key: str,
) -> dict[str, Any] | None:
    path = combo_prep_state_path(ctx, cfg, combo_key=combo_key)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def _prep_state_attempt(payload: Mapping[str, Any] | None) -> int:
    if not isinstance(payload, Mapping):
        return 0
    for key in ("prep_attempt", "attempt"):
        try:
            value = int(payload.get(key) or 0)
        except Exception:
            value = 0
        if value > 0:
            return value
    return 1


def combo_prep_failure_retryable(
    cfg: Mapping[str, Any],
    payload: Mapping[str, Any] | None,
) -> bool:
    if not isinstance(payload, Mapping):
        return False
    status = str(payload.get("status") or "").strip().lower()
    if status and status != "failed":
        return False
    error = str(payload.get("error") or "").strip().lower()
    retryable_errors = {
        "owner_completed_without_receptor_artifact",
    }
    retryable = error in retryable_errors
    retryable = retryable or "openssl_3.3.0" in error or "libcrypto" in error
    if not retryable:
        return False
    try:
        max_attempts = int(
            os.environ.get("ATLAS_COMBO_PREP_MAX_ATTEMPTS")
            or cfg.get("ATLAS_COMBO_PREP_MAX_ATTEMPTS")
            or 2
        )
    except Exception:
        max_attempts = 2
    max_attempts = max(1, int(max_attempts))
    return _prep_state_attempt(payload) < max_attempts


def try_begin_combo_prep(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    combo_key: str,
    lease_sec: float = 14400.0,
) -> bool:
    if not ctx.enabled:
        return True

    state_payload = read_combo_prep_state(ctx, cfg, combo_key=combo_key)
    prep_attempt = 1
    if isinstance(state_payload, dict):
        state_status = str(state_payload.get("status") or "").strip().lower()
        if state_status == "ready":
            return False
        if state_status == "failed":
            if not combo_prep_failure_retryable(cfg, state_payload):
                return False
            prep_attempt = _prep_state_attempt(state_payload) + 1
        elif state_status:
            return False

    claim_path = combo_prep_claim_path(ctx, cfg, combo_key=combo_key)
    claim_path.parent.mkdir(parents=True, exist_ok=True)
    now = time.time()
    payload = {
        "run_id": ctx.run_id,
        "combo_key": str(combo_key),
        "task_id": int(ctx.task_id),
        "task_count": int(ctx.task_count),
        "claimed_at": float(now),
        "lease_sec": float(max(30.0, lease_sec)),
    }

    def _attempt_create() -> bool:
        try:
            fd = os.open(str(claim_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            return False
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, sort_keys=True)
                handle.flush()
                os.fsync(handle.fileno())
            return True
        except Exception:
            try:
                claim_path.unlink()
            except Exception:
                pass
            return False

    if not _attempt_create():
        try:
            if claim_path.exists():
                age = max(0.0, now - float(claim_path.stat().st_mtime))
                if age > float(max(30.0, lease_sec)):
                    try:
                        claim_path.unlink()
                    except Exception:
                        pass
                    if not _attempt_create():
                        return False
                else:
                    return False
            else:
                return False
        except Exception:
            return False

    write_json_atomic(
        combo_prep_state_path(ctx, cfg, combo_key=combo_key),
        {
            "run_id": ctx.run_id,
            "combo_key": str(combo_key),
            "status": "preparing",
            "task_id": int(ctx.task_id),
            "task_count": int(ctx.task_count),
            "updated_at": float(time.time()),
            "lease_sec": float(max(30.0, lease_sec)),
            "prep_attempt": int(prep_attempt),
        },
    )
    return True


def write_combo_prep_state(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    combo_key: str,
    status: str,
    payload: Mapping[str, Any] | None = None,
) -> Path:
    result_payload = {
        "run_id": ctx.run_id,
        "combo_key": str(combo_key),
        "task_id": int(ctx.task_id),
        "task_count": int(ctx.task_count),
        "status": str(status),
        "updated_at": float(time.time()),
    }
    if payload:
        for key, value in payload.items():
            if key in result_payload:
                continue
            result_payload[str(key)] = value

    path = combo_prep_state_path(ctx, cfg, combo_key=combo_key)
    written = write_json_atomic(path, result_payload)
    if str(status).strip().lower() in {"ready", "failed"}:
        try:
            claim = combo_prep_claim_path(ctx, cfg, combo_key=combo_key)
            if claim.exists():
                claim.unlink()
        except Exception:
            pass
    return written


def try_claim_chunk(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    chunk_id: str,
    lease_sec: float = 3600.0,
    retry_failed: bool = False,
    max_attempts: int = 1,
) -> bool:
    if not ctx.enabled:
        return True

    result_path = chunk_result_path(ctx, cfg, chunk_id=chunk_id)
    if result_path.exists():
        if not retry_failed:
            return False
        try:
            payload = json.loads(result_path.read_text(encoding="utf-8")) or {}
        except Exception:
            payload = {}
        status = str((payload or {}).get("status") or "").strip().lower()
        try:
            attempt = int((payload or {}).get("attempt") or 1)
        except Exception:
            attempt = 1
        if status == "failed" and int(attempt) < int(max(1, max_attempts)):
            try:
                result_path.unlink()
            except Exception:
                return False
        else:
            return False

    claim_path = chunk_claim_path(ctx, cfg, chunk_id=chunk_id)
    claim_path.parent.mkdir(parents=True, exist_ok=True)
    now = time.time()
    payload = {
        "run_id": ctx.run_id,
        "chunk_id": str(chunk_id),
        "task_id": int(ctx.task_id),
        "task_count": int(ctx.task_count),
        "claimed_at": float(now),
        "lease_sec": float(max(1.0, lease_sec)),
    }

    def _attempt_create() -> bool:
        try:
            fd = os.open(str(claim_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            return False
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, sort_keys=True)
                handle.flush()
                os.fsync(handle.fileno())
            return True
        except Exception:
            try:
                claim_path.unlink()
            except Exception:
                pass
            return False

    if _attempt_create():
        return True

    try:
        if claim_path.exists():
            age = max(0.0, now - float(claim_path.stat().st_mtime))
            if age > float(max(30.0, lease_sec)):
                try:
                    claim_path.unlink()
                except Exception:
                    pass
                if _attempt_create():
                    return True
    except Exception:
        pass
    return False


def renew_chunk_claim(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    chunk_id: str,
    lease_sec: float | None = None,
) -> bool:
    """
    Best-effort lease heartbeat for an existing chunk claim.
    Returns False when the claim no longer exists or cannot be refreshed.
    """
    if not ctx.enabled:
        return True
    claim_path = chunk_claim_path(ctx, cfg, chunk_id=chunk_id)
    if not claim_path.exists():
        return False
    now = float(time.time())
    try:
        with claim_path.open("r+", encoding="utf-8") as handle:
            try:
                payload = json.load(handle) or {}
            except Exception:
                payload = {}
            if not isinstance(payload, dict):
                payload = {}
            current_lease = payload.get("lease_sec")
            if lease_sec is None:
                try:
                    lease_value = float(str(current_lease))
                except Exception:
                    lease_value = 3600.0
            else:
                lease_value = float(lease_sec)
            try:
                claimed_at = float(payload.get("claimed_at") or now)
            except Exception:
                claimed_at = now
            refreshed = {
                "run_id": ctx.run_id,
                "chunk_id": str(chunk_id),
                "task_id": int(ctx.task_id),
                "task_count": int(ctx.task_count),
                "claimed_at": claimed_at,
                "updated_at": now,
                "lease_sec": float(max(1.0, lease_value)),
            }
            handle.seek(0)
            handle.truncate(0)
            json.dump(refreshed, handle, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        os.utime(claim_path, (now, now))
        return True
    except FileNotFoundError:
        return False
    except Exception:
        return False


def write_chunk_result(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    chunk_id: str,
    status: str,
    payload: Mapping[str, Any] | None = None,
) -> Path:
    result_payload = {
        "run_id": ctx.run_id,
        "chunk_id": str(chunk_id),
        "task_id": int(ctx.task_id),
        "status": str(status),
        "completed_at": float(time.time()),
    }
    if payload:
        for key, value in payload.items():
            if key in result_payload:
                continue
            result_payload[str(key)] = value

    path = chunk_result_path(ctx, cfg, chunk_id=chunk_id)
    written = write_json_atomic(path, result_payload)
    try:
        claim = chunk_claim_path(ctx, cfg, chunk_id=chunk_id)
        if claim.exists():
            claim.unlink()
    except Exception:
        pass
    return written


def remove_own_marker_if_present(ctx: DistributedRunContext, cfg: Mapping[str, Any]) -> None:
    if not ctx.enabled:
        return
    try:
        marker = ctx.marker_path(cfg)
        if marker.exists():
            marker.unlink()
    except Exception:
        logging.getLogger("distributed").debug(
            "[distributed.marker.cleanup] task_id=%s action=skip",
            ctx.task_id,
            exc_info=True,
        )


def write_completion_marker(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    assigned_pdb_ids: Sequence[str],
    failed_entries: Sequence[tuple[str, str, str, str, str]],
) -> Path | None:
    if not ctx.enabled:
        return None
    marker_dir = ctx.marker_dir(cfg)
    marker_dir.mkdir(parents=True, exist_ok=True)
    marker_path = ctx.marker_path(cfg)
    payload = {
        "run_id": ctx.run_id,
        "mode": ctx.mode,
        "task_id": int(ctx.task_id),
        "task_count": int(ctx.task_count),
        "task_min_id": int(ctx.task_min_id),
        "leader_task_id": int(ctx.leader_task_id),
        "assigned_pdb_ids": [str(x) for x in assigned_pdb_ids],
        "failed_entries": [list(entry) for entry in failed_entries],
        "completed_at": time.time(),
    }

    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=marker_path.name + ".",
            suffix=".tmp",
            dir=marker_dir,
            delete=False,
        ) as handle:
            tmp_path = Path(handle.name)
            json.dump(payload, handle, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        tmp_path.replace(marker_path)
        return marker_path
    finally:
        try:
            if tmp_path and tmp_path.exists():
                tmp_path.unlink()
        except Exception:
            pass


def wait_for_all_markers(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
) -> tuple[bool, list[int]]:
    if not ctx.enabled:
        return True, []
    expected_ids = ctx.expected_task_ids()
    if not expected_ids:
        return True, []

    marker_dir = ctx.marker_dir(cfg)
    deadline = time.time() + float(ctx.barrier_timeout_sec)
    while True:
        missing = [
            tid for tid in expected_ids if not (marker_dir / f"task_{tid}.json").exists()
        ]
        if not missing:
            return True, []
        if time.time() >= deadline:
            return False, missing
        time.sleep(float(ctx.barrier_poll_sec))


def collect_failed_entries_from_markers(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
) -> list[tuple[str, str, str, str, str]]:
    if not ctx.enabled:
        return []
    out: list[tuple[str, str, str, str, str]] = []
    for tid in ctx.expected_task_ids():
        marker_path = ctx.marker_path(cfg, task_id=tid)
        if not marker_path.exists():
            continue
        try:
            payload = json.loads(marker_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        entries = payload.get("failed_entries") or []
        if not isinstance(entries, list):
            continue
        for raw in entries:
            if not isinstance(raw, (list, tuple)) or len(raw) < 5:
                continue
            out.append(
                (
                    str(raw[0]),
                    str(raw[1]),
                    str(raw[2]),
                    str(raw[3]),
                    str(raw[4]),
                )
            )
    return out


def _phase_token(raw: str) -> str:
    token = str(raw or "").strip().lower()
    token = re.sub(r"[^a-z0-9_]+", "_", token)
    token = token.strip("_")
    return token or "phase"


def _phase_marker_dir(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    phase: str,
) -> Path:
    return _distributed_dir(cfg, ctx.run_id) / "phase_markers" / _phase_token(phase)


def _phase_marker_path(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    phase: str,
    task_id: int | None = None,
) -> Path:
    tid = int(ctx.task_id) if task_id is None else int(task_id)
    return _phase_marker_dir(ctx, cfg, phase=phase) / f"task_{tid}.json"


def write_phase_marker(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    phase: str,
    payload: Mapping[str, Any] | None = None,
) -> Path | None:
    if not ctx.enabled:
        return None
    marker_path = _phase_marker_path(ctx, cfg, phase=phase)
    marker_payload: dict[str, Any] = {
        "run_id": str(ctx.run_id),
        "phase": _phase_token(phase),
        "task_id": int(ctx.task_id),
        "task_count": int(ctx.task_count),
        "task_min_id": int(ctx.task_min_id),
        "leader_task_id": int(ctx.leader_task_id),
        "updated_at": float(time.time()),
    }
    if payload:
        marker_payload.update({str(k): v for k, v in payload.items()})
    return write_json_atomic(marker_path, marker_payload)


def wait_for_all_phase_markers(
    ctx: DistributedRunContext,
    cfg: Mapping[str, Any],
    *,
    phase: str,
    timeout_sec: float | None = None,
    poll_sec: float | None = None,
) -> tuple[bool, list[int]]:
    if not ctx.enabled:
        return True, []
    expected_ids = ctx.expected_task_ids()
    if not expected_ids:
        return True, []
    marker_dir = _phase_marker_dir(ctx, cfg, phase=phase)
    timeout = (
        float(timeout_sec)
        if timeout_sec is not None
        else float(max(1.0, float(ctx.barrier_timeout_sec)))
    )
    poll = (
        float(poll_sec)
        if poll_sec is not None
        else float(max(0.1, float(ctx.barrier_poll_sec)))
    )
    deadline = time.time() + max(0.1, timeout)
    while True:
        missing = [
            tid for tid in expected_ids if not (marker_dir / f"task_{tid}.json").exists()
        ]
        if not missing:
            return True, []
        if time.time() >= deadline:
            return False, missing
        time.sleep(poll)
