from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Mapping, Optional


_GPU_EMPTY = {"", "-1", "none", "no", "off", "false", "void", "nodevfiles"}
_NVIDIA_KEYS = ("CUDA_VISIBLE_DEVICES", "NVIDIA_VISIBLE_DEVICES")
_AMD_KEYS = ("ROCR_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES", "GPU_DEVICE_ORDINAL")
_AMD_ENV_NAMES = ("scorch-rocm", "scorch-amd", "scorch-env-rocm", "scorch-env-amd")
_NVIDIA_ENV_NAMES = (
    "scorch-cuda",
    "scorch-nvidia",
    "scorch-env-cuda",
    "scorch-env-nvidia",
)
_AUTO_GPU_MIN_LIGANDS_DEFAULT = 128
_GPU_CHUNK_SIZE_DEFAULT = 1024
_PROBE_CACHE: dict[tuple[str, ...], tuple[bool, str, str]] = {}


@dataclass(frozen=True)
class ScorchDevicePlan:
    requested: str
    effective: str
    backend: str
    gpu_ids: tuple[str, ...]
    env_prefix: Optional[str]
    env_name: Optional[str]
    workers_per_gpu: int
    binding: str
    reason: str
    probe_ok: bool
    probe_message: str

    @property
    def gpu_count(self) -> int:
        return len(self.gpu_ids)

    def manifest_details(self) -> dict[str, Any]:
        return {
            "scorch_device_requested": self.requested,
            "scorch_device_effective": self.effective,
            "scorch_gpu_backend": self.backend,
            "scorch_gpu_count": self.gpu_count,
            "scorch_gpu_ids": list(self.gpu_ids),
            "scorch_gpu_workers_per_gpu": self.workers_per_gpu,
            "scorch_gpu_binding": self.binding,
            "scorch_device_reason": self.reason,
            "scorch_device_probe_ok": self.probe_ok,
            "scorch_device_probe": self.probe_message,
        }


def _token(raw: Any) -> str:
    return str(raw or "").strip()


def _cfg_token(cfg: Mapping[str, Any], key: str, default: str = "") -> str:
    return _token(cfg.get(key, default))


def _cfg_int(
    cfg: Mapping[str, Any], key: str, default: int, *, min_value: int = 0
) -> int:
    try:
        return max(min_value, int(_cfg_token(cfg, key, str(default)) or default))
    except (TypeError, ValueError):
        return max(min_value, int(default))


def _env_int(name: str, *, min_value: int) -> Optional[int]:
    env_value = os.environ.get(name)
    if env_value is None:
        return None
    try:
        return max(min_value, int(env_value.strip()))
    except (TypeError, ValueError):
        return None


def _normalize_device(raw: Any) -> str:
    token = _token(raw).lower()
    if token in {"auto", "gpu", "cuda", "rocm", "hip", "amd", "nvidia"}:
        return "gpu" if token in {"cuda", "rocm", "hip"} else token
    return "cpu"


def _normalize_backend(raw: Any) -> str:
    token = _token(raw).lower()
    if token in {"amd", "rocm", "hip"}:
        return "amd"
    if token in {"nvidia", "cuda"}:
        return "nvidia"
    return "unknown"


def _parse_visible(raw: Any) -> tuple[str, ...]:
    token = _token(raw)
    if token.lower() in _GPU_EMPTY:
        return ()
    parts = [part.strip() for part in token.replace(";", ",").split(",")]
    return tuple(part for part in parts if part and part.lower() not in _GPU_EMPTY)


def _visible_from_env(keys: tuple[str, ...]) -> tuple[str, ...]:
    for key in keys:
        ids = _parse_visible(os.environ.get(key))
        if ids:
            return ids
    return ()


def _ids_from_tool(tool: str, args: list[str]) -> tuple[str, ...]:
    exe = shutil.which(tool)
    if not exe:
        return ()
    try:
        proc = subprocess.run(
            [exe, *args],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return ()
    if proc.returncode != 0:
        return ()
    ids: list[str] = []
    for line in (proc.stdout or "").splitlines():
        token = line.strip()
        if token.isdigit():
            ids.append(token)
            continue
        first_col = token.split(",", 1)[0].strip().lower()
        if first_col.startswith("card") and first_col[4:].isdigit():
            ids.append(first_col[4:])
    return tuple(ids)


def _amd_ids_from_rocm_smi() -> tuple[str, ...]:
    return _ids_from_tool("rocm-smi", ["--showid", "--csv"])


def _nvidia_ids_from_nvidia_smi() -> tuple[str, ...]:
    return _ids_from_tool("nvidia-smi", ["--query-gpu=index", "--format=csv,noheader"])


def _device_backend(cfg: Mapping[str, Any], requested: str) -> str:
    if requested == "amd":
        return "amd"
    if requested == "nvidia":
        return "nvidia"
    return _normalize_backend(
        os.environ.get("ATLAS_SCORCH_GPU_BACKEND")
        or cfg.get("SCORCH_GPU_BACKEND")
        or ""
    )


def _detect_visible_gpus(
    cfg: Mapping[str, Any], requested_backend: str
) -> tuple[str, tuple[str, ...]]:
    forced_backend = requested_backend
    manual_ids = _parse_visible(
        os.environ.get("ATLAS_SCORCH_GPU_IDS") or cfg.get("SCORCH_GPU_IDS") or ""
    )
    amd_ids = _visible_from_env(_AMD_KEYS)
    nvidia_ids = _visible_from_env(_NVIDIA_KEYS)
    slurm_ids = _parse_visible(os.environ.get("SLURM_JOB_GPUS"))
    if forced_backend == "amd":
        if not amd_ids:
            amd_ids = _amd_ids_from_rocm_smi()
    elif forced_backend == "nvidia":
        if not nvidia_ids:
            nvidia_ids = _nvidia_ids_from_nvidia_smi()
    else:
        if not amd_ids:
            amd_ids = _amd_ids_from_rocm_smi()
        if not nvidia_ids:
            nvidia_ids = _nvidia_ids_from_nvidia_smi()
    forced = _forced_backend_visible_ids(
        forced_backend, manual_ids, amd_ids, nvidia_ids, slurm_ids
    )
    if forced is not None:
        return forced
    return _auto_backend_visible_ids(manual_ids, amd_ids, nvidia_ids, slurm_ids)


def _forced_backend_visible_ids(
    forced_backend: str,
    manual_ids: tuple[str, ...],
    amd_ids: tuple[str, ...],
    nvidia_ids: tuple[str, ...],
    slurm_ids: tuple[str, ...],
) -> Optional[tuple[str, tuple[str, ...]]]:
    if forced_backend == "amd":
        return "amd", manual_ids or amd_ids or slurm_ids
    if forced_backend == "nvidia":
        return "nvidia", manual_ids or nvidia_ids or slurm_ids
    return None


def _auto_backend_visible_ids(
    manual_ids: tuple[str, ...],
    amd_ids: tuple[str, ...],
    nvidia_ids: tuple[str, ...],
    slurm_ids: tuple[str, ...],
) -> tuple[str, tuple[str, ...]]:
    if manual_ids:
        return "unknown", manual_ids
    if amd_ids:
        return "amd", amd_ids
    if nvidia_ids:
        return "nvidia", nvidia_ids
    if slurm_ids:
        return "unknown", slurm_ids
    return "unknown", ()


def _generic_env_prefix(cfg: Mapping[str, Any]) -> Optional[str]:
    return _cfg_token(cfg, "SCORCH_ENV_PREFIX") or None


def _generic_env_name(cfg: Mapping[str, Any]) -> str:
    return _cfg_token(cfg, "SCORCH_ENV", "scorch-env") or "scorch-env"


def _env_prefix_candidates(cfg: Mapping[str, Any], backend: str) -> list[Path]:
    candidates: list[Path] = []
    generic = _generic_env_prefix(cfg)
    if generic:
        generic_path = Path(generic).expanduser()
        parent = generic_path.parent if generic_path.name else generic_path
        names = _AMD_ENV_NAMES if backend == "amd" else _NVIDIA_ENV_NAMES
        candidates.extend(parent / name for name in names)

    home = Path.home()
    roots = [
        home / "micromamba" / "envs",
        home / ".mamba" / "envs",
        home / "miniconda3" / "envs",
        home / "anaconda3" / "envs",
    ]
    names = _AMD_ENV_NAMES if backend == "amd" else _NVIDIA_ENV_NAMES
    for root in roots:
        candidates.extend(root / name for name in names)
    return candidates


def _autodetected_env_prefix(cfg: Mapping[str, Any], backend: str) -> str:
    if backend not in {"amd", "nvidia"}:
        return ""
    for candidate in _env_prefix_candidates(cfg, backend):
        python_exe = candidate / "bin" / "python"
        if python_exe.is_file() and os.access(python_exe, os.X_OK):
            return str(candidate.resolve())
    return ""


def _backend_value(
    cfg: Mapping[str, Any],
    backend: str,
    suffix: str,
    generic_key: str,
) -> str:
    keys = []
    if backend == "amd":
        keys.append(f"{generic_key}_AMD")
    elif backend == "nvidia":
        keys.append(f"{generic_key}_NVIDIA")
    keys.append(f"{generic_key}_GPU")
    for key in keys:
        value = _cfg_token(cfg, key)
        if value:
            return value
    env_prefix = (
        _autodetected_env_prefix(cfg, backend)
        if generic_key == "SCORCH_ENV_PREFIX"
        else ""
    )
    if env_prefix:
        return env_prefix
    value = _cfg_token(cfg, generic_key)
    if value:
        return value
    return suffix


def _runner_cmd(env_prefix: Optional[str], env_name: Optional[str]) -> list[str]:
    micromamba = shutil.which("micromamba")
    if micromamba:
        cmd = [micromamba, "run"]
        if env_prefix:
            cmd.extend(["-p", env_prefix])
        elif env_name:
            cmd.extend(["-n", env_name])
        cmd.append("python")
        return cmd
    if env_prefix:
        candidate = Path(env_prefix).expanduser() / "bin" / "python"
        if candidate.is_file():
            return [str(candidate)]
    return [sys.executable]


def _probe_tensorflow_gpu(
    *,
    backend: str,
    gpu_ids: tuple[str, ...],
    env_prefix: Optional[str],
    env_name: Optional[str],
) -> tuple[bool, str]:
    key = (
        backend,
        ",".join(gpu_ids),
        str(env_prefix or ""),
        str(env_name or ""),
    )
    cached = _PROBE_CACHE.get(key)
    if cached is not None:
        return cached[0], cached[1]

    code = (
        "import json, tensorflow as tf; "
        "g=tf.config.list_physical_devices('GPU'); "
        "print(json.dumps({'tensorflow': tf.__version__, "
        "'gpu_count': len(g), 'gpus': [d.name for d in g]}))"
    )
    cmd = _runner_cmd(env_prefix, env_name) + ["-c", code]
    env = scorch_child_env_for(
        effective="gpu",
        backend=backend,
        gpu_ids=gpu_ids,
        binding="all",
        task_id=None,
        base_env=os.environ,
    )
    proc, error = _run_probe_command(cmd, env)
    if error:
        _PROBE_CACHE[key] = (False, error, "")
        return False, error
    ok, message = _probe_result(proc)
    _PROBE_CACHE[key] = (ok, message, "")
    return ok, message


def _run_probe_command(
    cmd: list[str], env: Mapping[str, str]
) -> tuple[Optional[subprocess.CompletedProcess[str]], str]:
    try:
        proc = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
            env=env,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return None, f"probe_exception:{type(exc).__name__}:{exc}"
    return proc, ""


def _probe_result(proc: Optional[subprocess.CompletedProcess[str]]) -> tuple[bool, str]:
    if proc is None:
        return False, "probe_missing_process"
    stdout = (proc.stdout or "").strip().splitlines()
    last = stdout[-1] if stdout else ""
    if proc.returncode != 0:
        stderr = (proc.stderr or "").strip().splitlines()
        message = "probe_failed"
        if stderr:
            message = f"{message}:{stderr[-1][:240]}"
        return False, message
    try:
        payload = json.loads(last)
        count = int(payload.get("gpu_count", 0))
    except (json.JSONDecodeError, TypeError, ValueError):
        count = 0
    if count <= 0:
        return False, "tensorflow_reports_no_gpu"
    return True, last[:500]


def _cpu_plan(
    cfg: Mapping[str, Any],
    *,
    requested: str,
    backend: str,
    gpu_ids: tuple[str, ...],
    workers_per_gpu: int,
    binding: str,
    reason: str,
    probe_message: str,
) -> ScorchDevicePlan:
    return ScorchDevicePlan(
        requested=requested,
        effective="cpu",
        backend=backend,
        gpu_ids=gpu_ids,
        env_prefix=_generic_env_prefix(cfg),
        env_name=_generic_env_name(cfg),
        workers_per_gpu=workers_per_gpu,
        binding=binding,
        reason=reason,
        probe_ok=False,
        probe_message=probe_message,
    )


def _gpu_plan(
    *,
    requested: str,
    backend: str,
    gpu_ids: tuple[str, ...],
    env_prefix: Optional[str],
    env_name: Optional[str],
    workers_per_gpu: int,
    binding: str,
    probe_message: str,
) -> ScorchDevicePlan:
    return ScorchDevicePlan(
        requested=requested,
        effective="gpu",
        backend=backend,
        gpu_ids=gpu_ids,
        env_prefix=env_prefix,
        env_name=env_name,
        workers_per_gpu=workers_per_gpu,
        binding=binding,
        reason="tensorflow_gpu_available",
        probe_ok=True,
        probe_message=probe_message,
    )


def _log_gpu_plan(
    logger: Optional[logging.Logger],
    *,
    requested: str,
    backend: str,
    gpu_count: int,
    env_prefix: Optional[str],
    env_name: Optional[str],
    binding: str,
    workers_per_gpu: int,
) -> None:
    if logger is None:
        return
    logger.info(
        "[scorch-device] requested=%s effective=gpu backend=%s gpu_count=%d env_prefix=%s env_name=%s binding=%s workers_per_gpu=%d",
        requested,
        backend,
        gpu_count,
        env_prefix or "",
        env_name or "",
        binding,
        workers_per_gpu,
    )


def _device_settings(cfg: Mapping[str, Any]) -> tuple[str, str, int]:
    requested = _normalize_device(
        os.environ.get("ATLAS_SCORCH_DEVICE") or cfg.get("SCORCH_DEVICE") or "cpu"
    )
    binding = _cfg_token(cfg, "SCORCH_GPU_BINDING", "single").lower()
    if binding not in {"single", "all"}:
        binding = "single"
    try:
        workers_per_gpu = max(1, int(cfg.get("SCORCH_GPU_WORKERS_PER_GPU", 1) or 1))
    except (TypeError, ValueError):
        workers_per_gpu = 1
    return requested, binding, workers_per_gpu


def _no_visible_gpu_plan(
    cfg: Mapping[str, Any],
    logger: Optional[logging.Logger],
    *,
    requested: str,
    workers_per_gpu: int,
    binding: str,
) -> ScorchDevicePlan:
    if logger is not None and requested == "gpu":
        logger.warning("[scorch-device] requested=gpu effective=cpu reason=no_visible_gpu")
    return _cpu_plan(
        cfg,
        requested=requested,
        backend="none",
        gpu_ids=(),
        workers_per_gpu=workers_per_gpu,
        binding=binding,
        reason="no_visible_gpu",
        probe_message="skipped_no_visible_gpu",
    )


def _probe_failed_plan(
    cfg: Mapping[str, Any],
    logger: Optional[logging.Logger],
    *,
    requested: str,
    backend: str,
    gpu_ids: tuple[str, ...],
    workers_per_gpu: int,
    binding: str,
    probe_message: str,
) -> ScorchDevicePlan:
    if logger is not None:
        logger.warning(
            "[scorch-device] requested=%s effective=cpu backend=%s reason=probe_failed probe=%s",
            requested,
            backend,
            probe_message,
        )
    return _cpu_plan(
        cfg,
        requested=requested,
        backend=backend,
        gpu_ids=gpu_ids,
        workers_per_gpu=workers_per_gpu,
        binding=binding,
        reason="probe_failed",
        probe_message=probe_message,
    )


def build_scorch_device_plan(
    cfg: Mapping[str, Any],
    logger: Optional[logging.Logger] = None,
) -> ScorchDevicePlan:
    requested, binding, workers_per_gpu = _device_settings(cfg)
    if requested == "cpu":
        return _cpu_plan(
            cfg,
            requested=requested,
            backend="none",
            gpu_ids=(),
            workers_per_gpu=workers_per_gpu,
            binding=binding,
            reason="requested_cpu",
            probe_message="skipped_cpu",
        )

    requested_backend = _device_backend(cfg, requested)
    backend, gpu_ids = _detect_visible_gpus(cfg, requested_backend)
    env_prefix = _backend_value(cfg, backend, "", "SCORCH_ENV_PREFIX") or None
    env_name = _backend_value(cfg, backend, "scorch-env", "SCORCH_ENV") or None

    if not gpu_ids:
        return _no_visible_gpu_plan(
            cfg,
            logger,
            requested=requested,
            workers_per_gpu=workers_per_gpu,
            binding=binding,
        )

    probe_ok, probe_message = _probe_tensorflow_gpu(
        backend=backend,
        gpu_ids=gpu_ids,
        env_prefix=env_prefix,
        env_name=env_name,
    )
    if not probe_ok:
        return _probe_failed_plan(
            cfg,
            logger,
            requested=requested,
            backend=backend,
            gpu_ids=gpu_ids,
            workers_per_gpu=workers_per_gpu,
            binding=binding,
            probe_message=probe_message,
        )

    _log_gpu_plan(
        logger,
        requested=requested,
        backend=backend,
        gpu_count=len(gpu_ids),
        env_prefix=env_prefix,
        env_name=env_name,
        binding=binding,
        workers_per_gpu=workers_per_gpu,
    )
    return _gpu_plan(
        requested=requested,
        backend=backend,
        gpu_ids=gpu_ids,
        env_prefix=env_prefix,
        env_name=env_name,
        workers_per_gpu=workers_per_gpu,
        binding=binding,
        probe_message=probe_message,
    )


def apply_scorch_device_plan_to_cfg(
    cfg: Mapping[str, Any], plan: ScorchDevicePlan
) -> dict[str, Any]:
    out = dict(cfg)
    original_env_prefix = _cfg_token(cfg, "SCORCH_ENV_PREFIX")
    original_env_name = _cfg_token(cfg, "SCORCH_ENV")
    if original_env_prefix and not _cfg_token(out, "SCORCH_ENV_PREFIX_AUTO_CPU"):
        out["SCORCH_ENV_PREFIX_AUTO_CPU"] = original_env_prefix
    if original_env_name and not _cfg_token(out, "SCORCH_ENV_AUTO_CPU"):
        out["SCORCH_ENV_AUTO_CPU"] = original_env_name
    if plan.env_prefix:
        out["SCORCH_ENV_PREFIX"] = plan.env_prefix
    if plan.env_name:
        out["SCORCH_ENV"] = plan.env_name
    out["SCORCH_DEVICE_REQUESTED"] = plan.requested
    out["SCORCH_DEVICE_EFFECTIVE"] = plan.effective
    out["SCORCH_GPU_BACKEND_EFFECTIVE"] = plan.backend
    out["SCORCH_GPU_IDS_EFFECTIVE"] = ",".join(plan.gpu_ids)
    out["SCORCH_GPU_WORKERS_PER_GPU_EFFECTIVE"] = plan.workers_per_gpu
    out["SCORCH_GPU_BINDING_EFFECTIVE"] = plan.binding
    out["SCORCH_GPU_MIN_LIGANDS_EFFECTIVE"] = _auto_gpu_min_ligands(out)
    out["SCORCH_GPU_CHUNK_SIZE_EFFECTIVE"] = _gpu_chunk_size(out)
    out["SCORCH_DEVICE_REASON"] = plan.reason
    out["SCORCH_DEVICE_PROBE"] = plan.probe_message
    return out


def _auto_gpu_min_ligands(cfg: Mapping[str, Any]) -> int:
    env_value = _env_int("ATLAS_SCORCH_GPU_MIN_LIGANDS", min_value=0)
    if env_value is not None:
        return env_value
    return _cfg_int(cfg, "SCORCH_GPU_MIN_LIGANDS", _AUTO_GPU_MIN_LIGANDS_DEFAULT)


def _gpu_chunk_size(cfg: Mapping[str, Any]) -> int:
    env_value = _env_int("ATLAS_SCORCH_GPU_CHUNK_SIZE", min_value=1)
    if env_value is not None:
        return env_value
    return _cfg_int(
        cfg,
        "SCORCH_GPU_CHUNK_SIZE",
        _GPU_CHUNK_SIZE_DEFAULT,
        min_value=1,
    )


def scorch_chunk_size_for_cfg(
    cfg: Mapping[str, Any], default_chunk_size: int
) -> int:
    default_size = max(1, int(default_chunk_size))
    env_value = _env_int("ATLAS_SCORCH_CHUNK_SIZE", min_value=1)
    if env_value is not None:
        return env_value
    cfg_value = _cfg_int(cfg, "SCORCH_CHUNK_SIZE", 0, min_value=0)
    if cfg_value > 0:
        return cfg_value
    if _cfg_token(cfg, "SCORCH_DEVICE_EFFECTIVE", "cpu").lower() != "gpu":
        return default_size
    return max(1, _gpu_chunk_size(cfg))


def cfg_for_scorch_ligand_count(
    cfg: Mapping[str, Any],
    ligand_count: int,
    logger: Optional[logging.Logger] = None,
) -> dict[str, Any]:
    out = dict(cfg)
    requested = _cfg_token(out, "SCORCH_DEVICE_REQUESTED", out.get("SCORCH_DEVICE", "cpu"))
    effective = _cfg_token(out, "SCORCH_DEVICE_EFFECTIVE", "cpu")
    min_ligands = _auto_gpu_min_ligands(out)
    if (
        requested.lower() != "auto"
        or effective.lower() != "gpu"
        or min_ligands <= 0
        or int(ligand_count) >= min_ligands
    ):
        return out

    cpu_prefix = (
        _cfg_token(out, "SCORCH_ENV_PREFIX_CPU")
        or _cfg_token(out, "SCORCH_ENV_PREFIX_AUTO_CPU")
    )
    cpu_env = _cfg_token(out, "SCORCH_ENV_CPU") or _cfg_token(out, "SCORCH_ENV_AUTO_CPU")
    if cpu_prefix:
        out["SCORCH_ENV_PREFIX"] = cpu_prefix
    if cpu_env:
        out["SCORCH_ENV"] = cpu_env
    out["SCORCH_DEVICE_EFFECTIVE"] = "cpu"
    out["SCORCH_GPU_BACKEND_EFFECTIVE"] = "none"
    out["SCORCH_GPU_IDS_EFFECTIVE"] = ""
    out["SCORCH_DEVICE_REASON"] = f"auto_small_batch_lt_{min_ligands}"
    out["SCORCH_DEVICE_PROBE"] = f"skipped_small_batch_ligands={ligand_count}"
    if logger is not None:
        logger.info(
            "[scorch-device] requested=auto effective=cpu reason=small_batch ligands=%d min_gpu_ligands=%d",
            int(ligand_count),
            min_ligands,
        )
    return out


def cap_scorch_jobs_for_device(jobs: int, plan: ScorchDevicePlan) -> int:
    requested_jobs = max(1, int(jobs))
    if plan.effective != "gpu":
        return requested_jobs
    cap = max(1, plan.gpu_count * max(1, plan.workers_per_gpu))
    return max(1, min(requested_jobs, cap))


def cap_scorch_queue_workers(worker_cap: int, plan: ScorchDevicePlan) -> int:
    requested_workers = max(1, int(worker_cap))
    if plan.effective != "gpu":
        return requested_workers
    return 1


def _stable_index(task_id: Optional[str], count: int) -> int:
    if count <= 1:
        return 0
    digest = hashlib.sha256(str(task_id or "").encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % count


def scorch_child_env_for(
    *,
    effective: str,
    backend: str,
    gpu_ids: tuple[str, ...],
    binding: str,
    task_id: Optional[str],
    base_env: Mapping[str, str],
) -> dict[str, str]:
    env = dict(base_env)
    env["PYTHONNOUSERSITE"] = "1"
    env["ATLAS_SCORCH_DEVICE_EFFECTIVE"] = effective
    env["ATLAS_SCORCH_GPU_BACKEND"] = backend
    env["ATLAS_SCORCH_GPU_IDS"] = ",".join(gpu_ids)
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "TF_NUM_INTRAOP_THREADS",
        "TF_NUM_INTEROP_THREADS",
    ):
        env.setdefault(key, "1")
    if effective != "gpu":
        env["CUDA_VISIBLE_DEVICES"] = "-1"
        env["NVIDIA_VISIBLE_DEVICES"] = "void"
        env["ROCR_VISIBLE_DEVICES"] = ""
        env["HIP_VISIBLE_DEVICES"] = ""
        env["GPU_DEVICE_ORDINAL"] = ""
        return env

    selected_ids = gpu_ids
    if binding == "single" and gpu_ids:
        selected_ids = (gpu_ids[_stable_index(task_id, len(gpu_ids))],)
    visible = ",".join(selected_ids)
    env["ATLAS_SCORCH_GPU_ASSIGNED"] = visible
    if backend == "amd":
        env["ROCR_VISIBLE_DEVICES"] = visible
        env["HIP_VISIBLE_DEVICES"] = visible
        env["GPU_DEVICE_ORDINAL"] = visible
    elif backend == "nvidia":
        env["CUDA_VISIBLE_DEVICES"] = visible
        env["NVIDIA_VISIBLE_DEVICES"] = visible
    else:
        env["CUDA_VISIBLE_DEVICES"] = visible
        env["ROCR_VISIBLE_DEVICES"] = visible
        env["HIP_VISIBLE_DEVICES"] = visible
    return env


def scorch_child_env(
    cfg: Mapping[str, Any],
    *,
    task_id: Optional[str] = None,
    base_env: Optional[Mapping[str, str]] = None,
) -> dict[str, str]:
    gpu_ids = tuple(
        part
        for part in _cfg_token(cfg, "SCORCH_GPU_IDS_EFFECTIVE").split(",")
        if part.strip()
    )
    return scorch_child_env_for(
        effective=_cfg_token(cfg, "SCORCH_DEVICE_EFFECTIVE", "cpu"),
        backend=_cfg_token(cfg, "SCORCH_GPU_BACKEND_EFFECTIVE", "none"),
        gpu_ids=gpu_ids,
        binding=_cfg_token(cfg, "SCORCH_GPU_BINDING_EFFECTIVE", "single"),
        task_id=task_id,
        base_env=base_env or os.environ,
    )
