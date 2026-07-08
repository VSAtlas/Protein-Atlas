"""Experimental Vina-GPU adapter.

This module is intentionally not imported by the standard Atlas docking
pipeline. CPU Vina is the supported production backend.
"""

from __future__ import annotations

import fcntl
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Mapping, Sequence

from docking.score_io import extract_best_score


@dataclass(frozen=True)
class VinaGpuResult:
    ligand: str
    output_path: str
    score: float | None
    reason: str | None = None


@dataclass(frozen=True)
class VinaGpuPlan:
    backend: str
    exe: Path
    kernel_dir: Path
    device_kind: str
    gpu_id: str | None
    gpu_ids: tuple[str, ...]
    thread: int
    search_depth: int | None
    chunks_per_gpu: int
    min_chunk_ligands: int


def plan_vina_gpu(cfg: Mapping[str, Any]) -> VinaGpuPlan:
    repo_root = Path(__file__).resolve().parents[2]
    default_root = repo_root / "tools" / "Vina-GPU-2.1" / "AutoDock-Vina-GPU-2.1"

    kind = _detect_device_kind(cfg)
    exe = _first_existing(
        _configured_path(cfg, f"VINA_GPU_EXE_{kind.upper()}"),
        _configured_path(cfg, "VINA_GPU_EXE"),
        default_root / "AutoDock-Vina-GPU-2-1",
    )
    kernel_dir = _first_existing(
        _configured_path(cfg, f"VINA_GPU_OPENCL_BINARY_PATH_{kind.upper()}"),
        _configured_path(cfg, "VINA_GPU_OPENCL_BINARY_PATH"),
        exe.parent,
    )
    gpu_ids = _gpu_ids(cfg)
    gpu_id = gpu_ids[0] if gpu_ids else None
    thread = _positive_int(cfg.get("VINA_GPU_THREAD"), 8000)
    if bool(cfg.get("FAST_MODE")):
        thread = _positive_int(cfg.get("VINA_GPU_FAST_THREAD"), min(thread, 1000))
    search_depth = _optional_positive_int(cfg.get("VINA_GPU_SEARCH_DEPTH"))
    chunks_per_gpu = _positive_int(cfg.get("VINA_GPU_CHUNKS_PER_GPU"), 4)
    min_chunk_ligands = _positive_int(cfg.get("VINA_GPU_MIN_CHUNK_LIGANDS"), 96)

    if not exe.exists():
        raise FileNotFoundError(f"Vina-GPU executable not found: {exe}")
    if not kernel_dir.exists():
        raise FileNotFoundError(f"Vina-GPU OpenCL kernel path not found: {kernel_dir}")
    return VinaGpuPlan(
        backend="vina_gpu",
        exe=exe,
        kernel_dir=kernel_dir,
        device_kind=kind,
        gpu_id=gpu_id,
        gpu_ids=gpu_ids,
        thread=thread,
        search_depth=search_depth,
        chunks_per_gpu=chunks_per_gpu,
        min_chunk_ligands=min_chunk_ligands,
    )


def run_vina_gpu_stage(
    cfg: Mapping[str, Any],
    *,
    pdb_id: str,
    receptor_pdbqt: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    stage_name: str,
    stage_info: Mapping[str, Any],
    ligands: Sequence[str],
    expected_outputs: Mapping[str, str],
    config_run_dir: str,
    variant: str | None,
    ph_label: str | None,
    timeout: int | None = None,
) -> dict[str, VinaGpuResult]:
    plan = plan_vina_gpu(cfg)
    batch_root = _batch_root(
        config_run_dir=config_run_dir,
        pdb_id=pdb_id,
        variant=variant,
        ph_label=ph_label,
        stage_name=stage_name,
    )
    _reset_dir(batch_root)

    results: dict[str, VinaGpuResult] = {}
    seen_stems: set[str] = set()
    batch_ligands: list[str] = []
    for lig in ligands:
        source = Path(lig)
        stem = source.stem
        if stem in seen_stems:
            results[lig] = VinaGpuResult(
                ligand=lig,
                output_path=str(expected_outputs[lig]),
                score=None,
                reason="duplicate_ligand_stem",
            )
            continue
        seen_stems.add(stem)
        batch_ligands.append(lig)

    output_dir_by_lig: dict[str, Path] = {}
    if batch_ligands:
        batch_specs = _prepare_batch_specs(
            plan=plan,
            batch_root=batch_root,
            stage_name=stage_name,
            batch_ligands=batch_ligands,
            receptor_pdbqt=receptor_pdbqt,
            center=center,
            box_size=box_size,
            stage_info=stage_info,
            output_dir_by_lig=output_dir_by_lig,
        )
        failures: list[subprocess.CompletedProcess[str]] = []
        if len(batch_specs) == 1:
            proc = _run_batch_process(
                plan=plan,
                config_run_dir=config_run_dir,
                config_path=batch_specs[0],
                timeout=timeout,
            )
            if proc.returncode != 0:
                failures.append(proc)
        else:
            with ThreadPoolExecutor(max_workers=len(batch_specs)) as pool:
                futures = [
                    pool.submit(
                        _run_batch_process,
                        plan=plan,
                        config_run_dir=config_run_dir,
                        config_path=config_path,
                        timeout=timeout,
                    )
                    for config_path in batch_specs
                ]
                for future in as_completed(futures):
                    proc = future.result()
                    if proc.returncode != 0:
                        failures.append(proc)
        if failures:
            first = failures[0]
            reason = _summarize_failure(first.stderr, first.stdout, first.returncode)
            for lig in batch_ligands:
                results[lig] = VinaGpuResult(
                    ligand=lig,
                    output_path=str(expected_outputs[lig]),
                    score=None,
                    reason=reason,
                )
            return results

    for lig in batch_ligands:
        source = Path(lig)
        output_dir = output_dir_by_lig.get(lig)
        produced = (
            output_dir / f"{source.stem}_out.pdbqt"
            if output_dir is not None
            else Path()
        )
        expected = Path(expected_outputs[lig])
        if not produced.exists() or produced.stat().st_size <= 0:
            results[lig] = VinaGpuResult(
                ligand=lig,
                output_path=str(expected),
                score=None,
                reason="missing_gpu_output",
            )
            continue
        expected.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(produced, expected)
        results[lig] = VinaGpuResult(
            ligand=lig,
            output_path=str(expected),
            score=extract_best_score(str(expected)),
        )
    return results


def _prepare_batch_specs(
    *,
    plan: VinaGpuPlan,
    batch_root: Path,
    stage_name: str,
    batch_ligands: Sequence[str],
    receptor_pdbqt: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    stage_info: Mapping[str, Any],
    output_dir_by_lig: dict[str, Path],
) -> list[Path]:
    chunk_count = _chunk_count(plan, len(batch_ligands))
    chunks = _split_evenly(list(batch_ligands), chunk_count)
    config_paths: list[Path] = []
    for idx, chunk in enumerate(chunks):
        chunk_root = batch_root if len(chunks) == 1 else batch_root / f"chunk_{idx:02d}"
        input_dir = chunk_root / "ligands"
        output_dir = chunk_root / "outputs"
        input_dir.mkdir(parents=True, exist_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)
        for lig in chunk:
            source = Path(lig)
            linked = input_dir / source.name
            try:
                linked.symlink_to(source.resolve())
            except Exception:
                shutil.copy2(source, linked)
            output_dir_by_lig[lig] = output_dir

        config_path = chunk_root / f"{stage_name}_vina_gpu.txt"
        _write_vina_gpu_config(
            config_path,
            receptor_pdbqt=receptor_pdbqt,
            input_dir=input_dir,
            output_dir=output_dir,
            center=center,
            box_size=box_size,
            stage_info=stage_info,
            plan=plan,
        )
        config_paths.append(config_path)
    return config_paths


def _chunk_count(plan: VinaGpuPlan, ligand_count: int) -> int:
    if ligand_count <= 0:
        return 0
    gpu_count = max(1, len(plan.gpu_ids))
    size_limited_chunks = (ligand_count + plan.min_chunk_ligands - 1) // (
        plan.min_chunk_ligands
    )
    target_chunks = max(
        gpu_count,
        min(gpu_count * plan.chunks_per_gpu, size_limited_chunks),
    )
    return min(ligand_count, target_chunks)


def _split_evenly(values: list[str], chunk_count: int) -> list[list[str]]:
    chunks: list[list[str]] = [[] for _ in range(max(1, chunk_count))]
    for idx, value in enumerate(values):
        chunks[idx % len(chunks)].append(value)
    return [chunk for chunk in chunks if chunk]


def _run_batch_process(
    *,
    plan: VinaGpuPlan,
    config_run_dir: str,
    config_path: Path,
    timeout: int | None,
) -> subprocess.CompletedProcess[str]:
    cmd = [str(plan.exe), "--config", str(config_path)]
    with _gpu_lock(config_run_dir, plan) as gpu_id:
        env = _bound_gpu_env(os.environ.copy(), plan, gpu_id)
        return subprocess.run(
            cmd,
            cwd=str(plan.exe.parent),
            env=env,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )


def _configured_path(cfg: Mapping[str, Any], key: str) -> Path | None:
    raw = cfg.get(key)
    if raw is None or str(raw).strip() == "":
        return None
    return Path(str(raw)).expanduser()


def _first_existing(*candidates: Path | None) -> Path:
    for candidate in candidates:
        if candidate is not None and candidate.exists():
            return candidate
    for candidate in candidates:
        if candidate is not None:
            return candidate
    raise FileNotFoundError("No Vina-GPU candidate path configured")


def _positive_int(value: Any, default: int) -> int:
    try:
        parsed = int(value)
        return parsed if parsed > 0 else default
    except Exception:
        return default


def _optional_positive_int(value: Any) -> int | None:
    if value in (None, ""):
        return None
    try:
        parsed = int(value)
        return parsed if parsed > 0 else None
    except Exception:
        return None


def _detect_device_kind(cfg: Mapping[str, Any]) -> str:
    raw = str(cfg.get("VINA_GPU_DEVICE") or cfg.get("VINA_DEVICE") or "auto")
    requested = raw.strip().lower()
    if requested in {"amd", "nvidia"}:
        return requested
    if _visible_env("ROCR_VISIBLE_DEVICES") or _visible_env("HIP_VISIBLE_DEVICES"):
        return "amd"
    if _visible_env("CUDA_VISIBLE_DEVICES"):
        return "nvidia"
    if shutil.which("rocm-smi") or Path("/opt/rocm/bin/rocm-smi").exists():
        return "amd"
    if shutil.which("nvidia-smi"):
        return "nvidia"
    return "amd"


def _visible_env(key: str) -> bool:
    raw = os.environ.get(key)
    return bool(raw and raw.strip() and raw.strip() not in {"-1", "none", "None"})


def _gpu_ids(cfg: Mapping[str, Any]) -> tuple[str, ...]:
    raw = str(
        cfg.get("VINA_GPU_IDS")
        or os.environ.get("VINA_GPU_IDS")
        or os.environ.get("ROCR_VISIBLE_DEVICES")
        or os.environ.get("HIP_VISIBLE_DEVICES")
        or os.environ.get("CUDA_VISIBLE_DEVICES")
        or ""
    ).strip()
    if not raw or raw in {"-1", "none", "None"}:
        return ()
    ids: list[str] = []
    seen: set[str] = set()
    for token in raw.replace(" ", ",").split(","):
        gpu_id = token.strip()
        if not gpu_id or gpu_id in seen:
            continue
        seen.add(gpu_id)
        ids.append(gpu_id)
    return tuple(ids)


def _bound_gpu_env(
    env: dict[str, str], plan: VinaGpuPlan, gpu_id: str | None
) -> dict[str, str]:
    if not gpu_id:
        return env
    if plan.device_kind == "amd":
        env["ROCR_VISIBLE_DEVICES"] = gpu_id
        env["HIP_VISIBLE_DEVICES"] = gpu_id
    elif plan.device_kind == "nvidia":
        env["CUDA_VISIBLE_DEVICES"] = gpu_id
        env["NVIDIA_VISIBLE_DEVICES"] = gpu_id
    return env


@contextmanager
def _gpu_lock(config_run_dir: str, plan: VinaGpuPlan) -> Iterator[str | None]:
    lock_dir = Path(config_run_dir) / "vina_gpu_locks"
    lock_dir.mkdir(parents=True, exist_ok=True)
    gpu_ids = plan.gpu_ids or (None,)
    handles: list[Any] = []
    selected: str | None = None
    acquired = False
    try:
        while not acquired:
            for gpu_id in gpu_ids:
                gpu_token = _safe_token(gpu_id or "default")
                lock_path = lock_dir / f"{_safe_token(plan.device_kind)}_{gpu_token}.lock"
                handle = lock_path.open("w", encoding="utf-8")
                try:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    handle.close()
                    continue
                handles.append(handle)
                selected = gpu_id
                acquired = True
                break
            if not acquired:
                time.sleep(0.25)
        yield selected
    finally:
        for handle in handles:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            finally:
                handle.close()


def _batch_root(
    *,
    config_run_dir: str,
    pdb_id: str,
    variant: str | None,
    ph_label: str | None,
    stage_name: str,
) -> Path:
    safe_variant = _safe_token(variant or "base")
    safe_ph = _safe_token(ph_label or "base")
    return (
        Path(config_run_dir)
        / "vina_gpu_batches"
        / _safe_token(pdb_id)
        / safe_variant
        / safe_ph
        / _safe_token(stage_name)
    )


def _safe_token(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value))


def _reset_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def _write_vina_gpu_config(
    path: Path,
    *,
    receptor_pdbqt: str,
    input_dir: Path,
    output_dir: Path,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    stage_info: Mapping[str, Any],
    plan: VinaGpuPlan,
) -> None:
    lines = [
        f"receptor = {receptor_pdbqt}",
        f"ligand_directory = {input_dir}",
        f"output_directory = {output_dir}",
        f"opencl_binary_path = {plan.kernel_dir}",
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x = {box_size[0]:.3f}",
        f"size_y = {box_size[1]:.3f}",
        f"size_z = {box_size[2]:.3f}",
        f"thread = {plan.thread}",
        f"energy_range = {float(stage_info.get('energy_range', 4))}",
        f"num_modes = {int(stage_info.get('num_modes', 4))}",
    ]
    if plan.search_depth is not None:
        lines.append(f"search_depth = {plan.search_depth}")
    if "seed" in stage_info:
        lines.append(f"seed = {int(stage_info['seed'])}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _summarize_failure(stderr: str, stdout: str, returncode: int) -> str:
    for stream in (stderr, stdout):
        lines = [line.strip() for line in str(stream or "").splitlines() if line.strip()]
        for line in reversed(lines):
            lower = line.lower()
            if "error" in lower or "failed" in lower or "cannot" in lower:
                return line
    return f"vina_gpu_returncode_{returncode}"
