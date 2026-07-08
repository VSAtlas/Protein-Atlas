"""Shared MPI launcher resolution and BLAS thread caps for MMGBSA subprocesses."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import MutableMapping

from config.tool_resolver import common_ambertools_prefixes
from config.value_access import cfg_get


def _coerce_int(value: object, default: int) -> int:
    try:
        if isinstance(value, (int, float, str, bytes, bytearray)):
            return int(float(value))
    except (TypeError, ValueError, OverflowError):
        pass
    return default


def resolve_md_mpi_ranks(cfg: object | None) -> int:
    raw = cfg_get(cfg, "MMGBSA_MD_MPI_RANKS", "auto")
    return _resolve_mpi_ranks_core(cfg, raw)


def resolve_mmpbsa_mpi_ranks(cfg: object | None) -> int:
    raw = cfg_get(cfg, "MMGBSA_MMPBSA_MPI_RANKS", None)
    if raw is None or str(raw).strip() == "":
        raw = cfg_get(cfg, "MMGBSA_MD_MPI_RANKS", "auto")
    return _resolve_mpi_ranks_core(cfg, raw)


def _resolve_mpi_ranks_core(cfg: object | None, raw: object) -> int:
    cpu = _coerce_int(cfg_get(cfg, "CPU", 0), 0)
    scheduler_cpu = _coerce_int(cfg_get(cfg, "GLOBAL_SCHEDULER_CPUS", 0), 0)
    env_cpu = _coerce_int(os.environ.get("GLOBAL_SCHEDULER_CPUS"), 0)
    capacity = scheduler_cpu or env_cpu or cpu or (os.cpu_count() or 1)
    capacity = max(1, min(capacity, 32))
    if raw is None or str(raw).strip().lower() == "auto":
        return capacity
    return max(1, min(_coerce_int(raw, capacity), capacity, 32))


def _mpi_exe_in_prefix_bin(prefix: Path) -> str | None:
    for launcher in ("mpirun", "mpiexec"):
        candidate = prefix / "bin" / launcher
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


def _mpi_exe_on_path() -> tuple[str, str] | None:
    for launcher in ("mpirun", "mpiexec"):
        path = shutil.which(launcher)
        if path:
            return path, "PATH"
    return None


def select_mpi_launcher(
    amber_prefix: Path | None,
    *,
    probe_common_prefixes: bool = False,
    cfg_for_common_prefixes: object | None = None,
) -> tuple[str, str]:
    if amber_prefix is not None:
        found = _mpi_exe_in_prefix_bin(amber_prefix)
        if found:
            return found, "prefix"
    if probe_common_prefixes:
        for prefix_candidate in common_ambertools_prefixes(cfg_for_common_prefixes):
            found = _mpi_exe_in_prefix_bin(prefix_candidate)
            if found:
                return found, f"common_prefix:{prefix_candidate}"
    on_path = _mpi_exe_on_path()
    if on_path:
        return on_path
    return "", "missing"


def apply_blas_single_thread_env_defaults(env: MutableMapping[str, str]) -> None:
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")
    env.setdefault("VECLIB_MAXIMUM_THREADS", "1")
