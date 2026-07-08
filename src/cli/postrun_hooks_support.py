from __future__ import annotations

from dataclasses import dataclass
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from config.normalize import _to_bool


def _as_int(val: Any, fallback: int) -> int:
    try:
        return int(val)
    except Exception:
        return fallback


def _parse_duration_seconds(raw: Any) -> Optional[float]:
    token = str(raw or "").strip()
    if not token:
        return None
    if token.isdigit():
        return float(int(token))
    day_part = 0
    hms = token
    if "-" in token:
        left, right = token.split("-", 1)
        try:
            day_part = int(left)
        except Exception:
            return None
        hms = right
    bits = hms.split(":")
    if len(bits) == 3:
        hh_s, mm_s, ss_s = bits
    elif len(bits) == 2:
        hh_s = "0"
        mm_s, ss_s = bits
    else:
        return None
    try:
        hh = int(hh_s)
        mm = int(mm_s)
        ss = float(ss_s)
    except Exception:
        return None
    return float(day_part * 86400 + hh * 3600 + mm * 60 + ss)


def _is_scorch_enabled(cfg: Mapping[str, Any]) -> bool:
    raw_use = cfg.get("USE_SCORCH", False)
    if isinstance(raw_use, bool):
        return raw_use
    return str(raw_use).strip().lower() in {"1", "true", "yes", "on"}


def _resolve_scorch_decoy_prefix(cfg: Mapping[str, Any]) -> str:
    """
    Resolve the effective SCORCH decoy prefix from cfg.
    Prefer rescoring_scorch's own parser to stay aligned with production rules.
    """
    scorch_mod: Any = None
    try:
        from post_docking.rescoring import rescoring_scorch as _scorch_mod_import

        scorch_mod = _scorch_mod_import
    except Exception:
        pass

    if scorch_mod is not None:
        infer = getattr(scorch_mod, "_infer_decoy_prefix_from_test_mode", None)
        if callable(infer):
            return str(infer(dict(cfg)))

    def _prefix_from_env_or_cfg(key: str, fallback: str = "decoys") -> str:
        raw = os.environ.get(key)
        if raw is None:
            raw = cfg.get(key)
        return str(raw or fallback).strip() or fallback

    def _dud_prefix() -> str:
        return _prefix_from_env_or_cfg(
            "DUD_PREFIX",
            _prefix_from_env_or_cfg("DECOY_PREFIX", "decoys"),
        )

    test_mode = str(
        os.environ.get("TEST_MODE_ENABLE")
        if os.environ.get("TEST_MODE_ENABLE") is not None
        else cfg.get("TEST_MODE_ENABLE", "off")
        or "off"
    ).strip().lower()
    if test_mode in {"", "0", "false", "no", "off", "none", "null"}:
        return _dud_prefix()
    if test_mode in {"1", "true", "yes", "on", "dud"}:
        return _dud_prefix()
    first = next((tok for tok in test_mode.replace(",", "+").split("+") if tok.strip()), "")
    if first in {"default", "fda"}:
        return _dud_prefix()
    return first or _dud_prefix()


def _normalize_scorch_parallel_profile(raw: Any) -> str:
    token = str(raw or "").strip().lower()
    if token in {"conservative", "safe", "balanced"}:
        return "conservative"
    return "aggressive"


def _resolve_scorch_parallel_profile(cfg: Mapping[str, Any]) -> str:
    env_raw = os.environ.get("ATLAS_SCORCH_PARALLEL_PROFILE")
    if env_raw:
        return _normalize_scorch_parallel_profile(env_raw)
    return _normalize_scorch_parallel_profile(cfg.get("SCORCH_PARALLEL_PROFILE"))


def _resolve_scorch_execution_mode(cfg: Mapping[str, Any]) -> str:
    def _current_interpreter_matches_env_prefix() -> bool:
        env_prefix = str(
            cfg.get("SCORCH_ENV_PREFIX") or os.environ.get("SCORCH_ENV_PREFIX") or ""
        ).strip()
        if not env_prefix:
            return False
        env_python = Path(env_prefix).expanduser() / "bin" / "python"
        try:
            return env_python.is_file() and env_python.resolve() == Path(sys.executable).resolve()
        except Exception:
            return False

    env_raw = os.environ.get("ATLAS_SCORCH_PERSISTENT_WORKER")
    if env_raw is not None:
        if not _to_bool(env_raw, False):
            return "subprocess"
        return "inprocess" if _current_interpreter_matches_env_prefix() else "subprocess"

    distributed_mode = str(
        os.environ.get("ATLAS_DISTRIBUTED_MODE", cfg.get("ATLAS_DISTRIBUTED_MODE", ""))
    ).strip().lower()
    if distributed_mode != "slurm_array":
        return "subprocess"

    return "inprocess" if _current_interpreter_matches_env_prefix() else "subprocess"


def _scheduler_queue_depth(scheduler: Any) -> int:
    if scheduler is None:
        return 0
    if hasattr(scheduler, "queued_requests"):
        try:
            return max(0, int(scheduler.queued_requests()))
        except Exception:
            return 0
    return 0


def _allocate_scorch_jobs_threads(
    *,
    cpu_budget: int,
    free_cores: int,
    backlog: int,
    allowed_count: int,
    scheduler_queue_depth: int,
    local_queue_depth: int,
    profile: str,
) -> tuple[int, int, float]:
    scorch_mod: Any = None
    try:
        from post_docking.rescoring import rescoring_scorch as _scorch_mod_import

        scorch_mod = _scorch_mod_import
    except Exception:
        pass

    allocator = (
        getattr(scorch_mod, "_allocate_scorch_parallelism", None)
        if scorch_mod is not None
        else None
    )
    if callable(allocator):
        try:
            jobs, threads, util_target = allocator(
                cpu_budget=int(cpu_budget),
                free_cores=int(free_cores),
                backlog=int(backlog),
                allowed_count=int(allowed_count),
                scheduler_queue_depth=int(scheduler_queue_depth),
                local_queue_depth=int(local_queue_depth),
                profile=str(profile),
            )
            return (
                max(1, int(jobs)),
                max(1, int(threads)),
                float(util_target),
            )
        except Exception:
            pass

    budget = max(1, int(cpu_budget))
    free = max(1, min(int(free_cores), budget))
    pending = max(1, int(backlog))
    if profile == "conservative":
        util_target = 0.85
        t_cap = 4
    else:
        util_target = 0.93
        t_cap = 6
    threads = min(t_cap, max(1, free // 32 + 1))
    util_cores = max(1, int(math.floor(util_target * free)))
    jobs = max(1, min(pending, util_cores // max(1, threads)))
    return int(jobs), int(threads), float(util_target)


def _scorch_subprocess_env(cfg: Mapping[str, Any]) -> Dict[str, str]:
    env = dict(os.environ)
    provisional_keys = {
        "SCORCH_PROVISIONAL_ENABLE",
        "SCORCH_PROVISIONAL_REUSE",
        "SCORCH_PROVISIONAL_CACHE_WRITE",
        "SCORCH_PROVISIONAL_TOP_FRACTION",
        "SCORCH_PROVISIONAL_MIN_PROGRESS",
        "SCORCH_PROVISIONAL_RESERVED_WORKERS",
        "SCORCH_PROVISIONAL_CPU_BUDGET",
        "SCORCH_PROVISIONAL_MAX_INFLIGHT",
        "SCORCH_PROVISIONAL_MAX_PER_COMBO",
    }
    forwarded_keys = (
        "USE_SCORCH",
        "TEST_MODE_ENABLE",
        "DUD_PREFIX",
        "DECOY_PREFIX",
        "SCORCH_TOP_FRACTION",
        "SCORCH_PARALLEL_PROFILE",
        "SCORCH_PROVISIONAL_ENABLE",
        "SCORCH_PROVISIONAL_REUSE",
        "SCORCH_PROVISIONAL_CACHE_WRITE",
        "SCORCH_PROVISIONAL_TOP_FRACTION",
        "SCORCH_PROVISIONAL_MIN_PROGRESS",
        "SCORCH_PROVISIONAL_RESERVED_WORKERS",
        "SCORCH_PROVISIONAL_CPU_BUDGET",
        "SCORCH_PROVISIONAL_MAX_INFLIGHT",
        "SCORCH_PROVISIONAL_MAX_PER_COMBO",
        "SCORCH_DEVICE",
        "SCORCH_GPU_BACKEND",
        "SCORCH_GPU_IDS",
        "SCORCH_GPU_BINDING",
        "SCORCH_GPU_WORKERS_PER_GPU",
        "SCORCH_GPU_MIN_LIGANDS",
        "SCORCH_GPU_CHUNK_SIZE",
        "SCORCH",
        "SCORCH_SCRIPT",
        "SCORCH_ENV_PREFIX",
        "SCORCH_ENV_PREFIX_GPU",
        "SCORCH_ENV_PREFIX_AMD",
        "SCORCH_ENV_PREFIX_NVIDIA",
        "SCORCH_ENV",
        "SCORCH_ENV_GPU",
        "SCORCH_ENV_AMD",
        "SCORCH_ENV_NVIDIA",
        "SCORCH_DEVICE_REQUESTED",
        "SCORCH_DEVICE_EFFECTIVE",
        "SCORCH_GPU_BACKEND_EFFECTIVE",
        "SCORCH_GPU_IDS_EFFECTIVE",
        "SCORCH_GPU_WORKERS_PER_GPU_EFFECTIVE",
        "SCORCH_GPU_BINDING_EFFECTIVE",
        "SCORCH_GPU_MIN_LIGANDS_EFFECTIVE",
        "SCORCH_GPU_CHUNK_SIZE_EFFECTIVE",
        "SCORCH_DEVICE_REASON",
        "SCORCH_DEVICE_PROBE",
        "CPU",
        "GLOBAL_SCHEDULER_CPUS",
        "OVERALL_DIR",
        "DOCKED_DIR",
        "POST_DOCKED_DIR",
        "OUTPUT_DIR",
        "DATA_DIR",
        "LOGS_DIR",
        "MANIFESTS_DIR",
    )
    for key in forwarded_keys:
        if key not in cfg:
            if key in provisional_keys:
                env.pop(key, None)
            continue
        value = cfg.get(key)
        if value is None:
            if key in provisional_keys:
                env.pop(key, None)
            continue
        if isinstance(value, bool):
            env[key] = "true" if value else "false"
        else:
            token = str(value).strip()
            if token:
                env[key] = token
    return env


def _with_repo_src_on_pythonpath(env: Dict[str, str], repo_root: Path) -> Dict[str, str]:
    out = dict(env)
    repo_path = str(repo_root.resolve())
    src_path = str((repo_root / "src").resolve())
    existing = out.get("PYTHONPATH", "")
    parts = [p for p in existing.split(os.pathsep) if p]
    if repo_path not in parts:
        parts.insert(0, repo_path)
    if src_path not in parts:
        parts.insert(0, src_path)
    out["PYTHONPATH"] = os.pathsep.join(parts)
    return out


def _scorch_env_python_prefix(cfg: Mapping[str, Any]) -> list[str] | None:
    raw_prefix = str(
        cfg.get("SCORCH_ENV_PREFIX") or os.environ.get("SCORCH_ENV_PREFIX") or ""
    ).strip()
    if raw_prefix:
        env_prefix = Path(raw_prefix).expanduser()
        env_python = env_prefix / "bin" / "python"
        if env_python.is_file() and os.access(env_python, os.X_OK):
            return [str(env_python.resolve())]
        return None

    raw_env = str(
        cfg.get("SCORCH_ENV") or os.environ.get("SCORCH_ENV") or "scorch-env"
    ).strip()
    if raw_env:
        runner = shutil.which("micromamba") or shutil.which("conda")
        if runner:
            return [runner, "run", "-n", raw_env, "python"]
    return None


def _allow_scorch_python_fallback(cfg: Mapping[str, Any]) -> bool:
    raw = os.environ.get("ATLAS_SCORCH_ALLOW_PYTHON_FALLBACK")
    if raw is None:
        raw = cfg.get("ATLAS_SCORCH_ALLOW_PYTHON_FALLBACK")
    return _to_bool(raw, False)


def _scorch_cmd_prefix(script_path: Path, cfg: Mapping[str, Any] | None = None) -> list[str]:
    cfg_map = cfg or {}
    env_prefix = _scorch_env_python_prefix(cfg_map)
    if env_prefix:
        return env_prefix + ["-m", "post_docking.rescoring.rescoring_scorch"]
    if os.environ.get("PYTEST_CURRENT_TEST") or _allow_scorch_python_fallback(cfg_map):
        return [sys.executable, "-m", "post_docking.rescoring.rescoring_scorch"]
    raise RuntimeError(
        "SCORCH is enabled but no usable SCORCH env was resolved. Set "
        "SCORCH_ENV_PREFIX to scorch-env, set SCORCH_ENV to a resolvable env, "
        "or explicitly set ATLAS_SCORCH_ALLOW_PYTHON_FALLBACK=1."
    )


@dataclass(frozen=True)
class _HookRoots:
    code_root: Path
    analysis_root: Path
    docked_root: Path
    post_docked_root: Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _to_resolved_path(raw: Any, *, base: Optional[Path] = None) -> Optional[Path]:
    if raw is None:
        return None
    token = str(raw).strip()
    if not token:
        return None
    path = Path(token).expanduser()
    if not path.is_absolute():
        anchor = base or Path.cwd()
        path = anchor / path
    return path.resolve()


def _infer_analysis_root(
    docked_root: Path,
    post_docked_root: Path,
    fallback: Path,
) -> Path:
    def _candidate(base: Path, leaf: str) -> Optional[Path]:
        if base.name == leaf:
            if base.parent.name == "outputs":
                return base.parent.parent
            return base.parent
        if base.parent.name == leaf:
            if base.parent.parent.name == "outputs":
                return base.parent.parent.parent
            return base.parent.parent
        return None

    docked_candidate = _candidate(docked_root, "docked")
    post_candidate = _candidate(post_docked_root, "post_docked")

    if docked_candidate is not None and post_candidate is not None:
        if docked_candidate == post_candidate:
            return docked_candidate
        return docked_candidate
    if docked_candidate is not None:
        return docked_candidate
    if post_candidate is not None:
        return post_candidate
    return fallback


def _resolve_hook_roots(cfg: Mapping[str, Any]) -> _HookRoots:
    code_root = _repo_root()
    overall_root = (
        _to_resolved_path(os.environ.get("OVERALL_DIR"), base=code_root)
        or _to_resolved_path(cfg.get("OVERALL_DIR"), base=code_root)
    )
    roots_base = overall_root or code_root
    docked_root = (
        _to_resolved_path(os.environ.get("DOCKED_DIR"), base=roots_base)
        or _to_resolved_path(cfg.get("DOCKED_DIR"), base=roots_base)
        or (roots_base / "docked").resolve()
    )
    post_docked_root = (
        _to_resolved_path(os.environ.get("POST_DOCKED_DIR"), base=roots_base)
        or _to_resolved_path(cfg.get("POST_DOCKED_DIR"), base=roots_base)
        or (roots_base / "post_docked").resolve()
    )
    analysis_root = _infer_analysis_root(docked_root, post_docked_root, roots_base)
    return _HookRoots(
        code_root=code_root,
        analysis_root=analysis_root,
        docked_root=docked_root,
        post_docked_root=post_docked_root,
    )
