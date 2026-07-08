from __future__ import annotations

import json
import logging
import os
import time
from typing import Any, Mapping, Optional

import yaml  # type: ignore[import-untyped]

from cli.chunk_tail_optimizer import optimize_chunks_for_tail
from cli.distributed_context import chunk_plan_path, write_json_atomic
from cli.planner_chunking import (
    _CHUNK_TARGET_SIZE,
    _distributed_chunk_local_workers,
    _planner_chunk_min_size,
    _split_ligands_into_chunks,
)
from cli.planner_combo_work import _normalize_ph_tag_token
from cli.planner_hash import _hash_token
from cli.planner_manifest_inventory import (
    _resolve_combo_ligand_streams,
    _validate_chunk_stream_resolution,
)
from config.output_paths import runtime_root


def _load_combo_weight_priors(
    cfg: Mapping[str, Any],
    *,
    current_run_id: str,
) -> dict[tuple[str, str, str, str], tuple[float, float]]:
    cached = cfg.get("_COMBO_WEIGHT_PRIORS")
    if isinstance(cached, dict):
        return cached

    priors: dict[tuple[str, str, str, str], tuple[float, float]] = {}
    manifests_root = runtime_root(
        cfg, "MANIFESTS_DIR", "manifests", prefer_existing=True
    )
    data_root = runtime_root(cfg, "DATA_DIR", "data", prefer_existing=True)
    if not manifests_root.exists():
        return priors

    run_dirs = sorted(
        [p for p in manifests_root.iterdir() if p.is_dir()],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    for run_dir in run_dirs[:60]:
        run_id = run_dir.name
        if run_id == str(current_run_id):
            continue
        manifest_path = run_dir / "run_manifest.yaml"
        if not manifest_path.exists():
            continue
        expected_map: dict[tuple[str, str, str, str], int] = {}
        ti_path = data_root / run_id / "throughput_integrity.json"
        if ti_path.exists():
            try:
                ti_payload = json.loads(ti_path.read_text(encoding="utf-8")) or {}
            except Exception:
                ti_payload = {}
            combos = ti_payload.get("combos") or []
            if isinstance(combos, list):
                for combo in combos:
                    if not isinstance(combo, dict):
                        continue
                    key = (
                        str(combo.get("pdb_id") or "").upper(),
                        str(combo.get("variant") or "").upper(),
                        str(combo.get("ph") or ""),
                        str(combo.get("library") or ""),
                    )
                    try:
                        expected = int(combo.get("expected_n") or 0)
                    except Exception:
                        expected = 0
                    if key[0] and key[1] and key[2] and expected > 0:
                        expected_map[key] = expected
        try:
            manifest_payload = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
        except Exception:
            continue
        proteins = manifest_payload.get("proteins") or {}
        if not isinstance(proteins, dict):
            continue
        for protein in proteins.values():
            if not isinstance(protein, dict):
                continue
            status = str(protein.get("status") or "").strip().lower()
            if status not in {"completed", "success", "succeeded"}:
                continue
            pdb_id = str(protein.get("pdb_id") or "").upper()
            variant = str(protein.get("variant") or "").upper()
            ph = str(protein.get("ph") or "")
            library = str(protein.get("library") or "")
            key = (pdb_id, variant, ph, library)
            if key in priors:
                continue
            wall = float(((protein.get("timing") or {}).get("wall_time_sec") or 0.0) or 0.0)
            expected = int(expected_map.get(key, 0))
            if wall <= 0.0 or expected <= 0:
                continue
            prep_sec = max(5.0, 0.15 * wall)
            per_ligand = max(0.05, (wall - prep_sec) / float(expected))
            priors[key] = (prep_sec, per_ligand)

    try:
        cfg["_COMBO_WEIGHT_PRIORS"] = priors  # type: ignore[index]
    except Exception:
        pass
    return priors


def _estimate_combo_weight_seconds(
    priors: Mapping[tuple[str, str, str, str], tuple[float, float]],
    *,
    pdb_id: str,
    variant: str,
    ph: str,
    library: str,
    ligand_count: int,
) -> float:
    key_exact = (pdb_id.upper(), variant.upper(), str(ph), str(library))
    prior = priors.get(key_exact)
    if prior is None:
        for (pdb_k, var_k, ph_k, _lib_k), value in priors.items():
            if pdb_k == key_exact[0] and var_k == key_exact[1] and ph_k == key_exact[2]:
                prior = value
                break
    if prior is None:
        prep_sec, per_ligand = 30.0, 1.0
    else:
        prep_sec, per_ligand = float(prior[0]), float(prior[1])
    return max(1.0, prep_sec + (max(0, int(ligand_count)) * max(0.01, per_ligand)))


_CHUNK_PLAN_SCHEMA_VERSION = 2


def _env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return bool(default)
    token = str(raw).strip().lower()
    if token in {"1", "true", "yes", "on"}:
        return True
    if token in {"0", "false", "no", "off"}:
        return False
    return bool(default)


def _env_int(name: str, default: int, *, min_value: int, max_value: int) -> int:
    raw = os.environ.get(name)
    if raw is None or not str(raw).strip():
        return int(default)
    try:
        parsed = int(str(raw).strip())
    except Exception:
        return int(default)
    return max(int(min_value), min(int(max_value), int(parsed)))


def _env_float(
    name: str,
    default: float,
    *,
    min_value: float,
    max_value: float,
) -> float:
    raw = os.environ.get(name)
    if raw is None or not str(raw).strip():
        return float(default)
    try:
        parsed = float(str(raw).strip())
    except Exception:
        return float(default)
    return max(float(min_value), min(float(max_value), float(parsed)))


def _chunk_plan_compatible(chunks: list[Any]) -> bool:
    if not chunks:
        return True
    for chunk in chunks:
        if not isinstance(chunk, dict):
            return False
        if not str(chunk.get("chunk_id") or "").strip():
            return False
        if not str(chunk.get("run_mode") or "").strip():
            return False
        if not str(chunk.get("library_name") or "").strip():
            return False
        if not str(chunk.get("library_root") or "").strip():
            return False
        if not isinstance(chunk.get("ligand_bases"), list):
            return False
    return True


def _maybe_validate_chunk_plan(
    cfg: Mapping[str, Any],
    chunks: list[dict[str, Any]],
    *,
    loaded_plan: bool = False,
) -> None:
    if not _env_bool("ATLAS_CHUNK_PREFLIGHT_VALIDATE", True):
        return
    if loaded_plan and not _env_bool("ATLAS_CHUNK_PREFLIGHT_VALIDATE_LOADED", False):
        return
    _validate_chunk_stream_resolution(
        cfg,
        chunks=chunks,
        sample_chunks_per_stream=_env_int(
            "ATLAS_CHUNK_PREFLIGHT_SAMPLE_CHUNKS",
            32,
            min_value=1,
            max_value=512,
        ),
        max_loss_fraction=_env_float(
            "ATLAS_CHUNK_PREFLIGHT_MAX_LOSS_FRACTION",
            0.05,
            min_value=0.0,
            max_value=1.0,
        ),
    )


def _distributed_combo_owner_task_id(
    dist_ctx: Any,
    *,
    variant_label: str,
    pdb_id: str,
    ph_tag: Optional[str],
) -> int:
    try:
        task_ids = [int(x) for x in dist_ctx.expected_task_ids()]
    except Exception:
        task_ids = []
    if not task_ids:
        try:
            return int(dist_ctx.task_id)
        except Exception:
            return 1
    ordered = sorted(set(task_ids))
    token = _hash_token(
        str(getattr(dist_ctx, "run_id", "")),
        str(variant_label or "base").upper(),
        str(pdb_id or "").upper(),
        _normalize_ph_tag_token(ph_tag),
    )
    slot = int(token[:8], 16) % len(ordered)
    return int(ordered[slot])


def _build_lpt_combo_chunks(
    cfg: Mapping[str, Any],
    *,
    run_id: str,
    dist_ctx: Any,
    variant_label: str,
    combo_items: list[tuple[str, Optional[str]]],
    run_tokens: list[str],
    target_ligands: int,
    sample_seed: int,
) -> list[dict[str, Any]]:
    priors = _load_combo_weight_priors(cfg, current_run_id=run_id)
    planner_started = time.perf_counter()
    planner_cache_hits_before = int(
        cfg.get("_PLANNER_MANIFEST_CACHE_HITS", 0) or 0
    )
    planner_cache_misses_before = int(
        cfg.get("_PLANNER_MANIFEST_CACHE_MISSES", 0) or 0
    )
    chunks: list[dict[str, Any]] = []
    resolved_combo_cache: dict[
        tuple[str, tuple[str, ...], int, int],
        list[dict[str, Any]],
    ] = {}
    local_workers = _distributed_chunk_local_workers(cfg)
    workers_total = max(1, int(dist_ctx.task_count)) * max(1, int(local_workers))
    production_chunks = int(target_ligands) <= 0
    min_chunk_size = _planner_chunk_min_size(
        1,
        workers_total,
        production=production_chunks,
    )
    scale_with_workers = _env_bool("ATLAS_CHUNK_SCALE_WITH_WORKERS", True)
    for pdb_file, ph_tag in combo_items:
        pdb_id = os.path.splitext(os.path.basename(str(pdb_file)))[0].upper()
        resolved_key = (
            pdb_id,
            tuple(run_tokens),
            int(target_ligands),
            int(sample_seed),
        )
        cached_resolved = resolved_combo_cache.get(resolved_key)
        if cached_resolved is None:
            cached_resolved = _resolve_combo_ligand_streams(
                cfg,
                pdb_id=pdb_id,
                run_tokens=run_tokens,
                target_ligands=target_ligands,
                sample_seed=sample_seed,
            )
            resolved_combo_cache[resolved_key] = [dict(stream) for stream in cached_resolved]
        streams = resolved_combo_cache[resolved_key]
        usable_streams = [
            stream
            for stream in streams
            if isinstance(stream, dict) and stream.get("ligand_bases")
        ]
        if not usable_streams:
            raise RuntimeError(
                "no_plannable_ligands_available:"
                f"pdb={pdb_id}:tokens={'+'.join(run_tokens) or 'fda'}"
            )
        for stream in usable_streams:
            run_mode = str(stream.get("run_mode") or "").strip() or "fda"
            library_name = str(stream.get("library_name") or "").strip()
            library_root = str(stream.get("library_root") or "").strip()
            ligand_bases = [str(x) for x in (stream.get("ligand_bases") or [])]
            stream_min_chunk = _planner_chunk_min_size(
                len(ligand_bases),
                workers_total,
                production=production_chunks,
            )
            stream_min_chunk = max(int(min_chunk_size), int(stream_min_chunk))
            split_chunks = _split_ligands_into_chunks(
                ligand_bases,
                task_count=max(1, int(dist_ctx.task_count)),
                local_workers_per_task=local_workers,
                ligand_weights=(
                    cfg.get("_PLANNER_LIGAND_COMPLEXITY_WEIGHTS")
                    if isinstance(cfg.get("_PLANNER_LIGAND_COMPLEXITY_WEIGHTS"), dict)
                    else None
                ),
                min_chunk_size=stream_min_chunk,
                scale_with_workers=scale_with_workers,
            )
            for idx, chunk_bases in enumerate(split_chunks):
                ph_token = _normalize_ph_tag_token(ph_tag)
                chunk_id = _hash_token(
                    run_id,
                    variant_label,
                    pdb_id,
                    ph_token,
                    run_mode,
                    str(library_name),
                    str(library_root),
                    str(idx),
                    ",".join(chunk_bases),
                )[:20]
                weight = _estimate_combo_weight_seconds(
                    priors,
                    pdb_id=pdb_id,
                    variant=variant_label,
                    ph=ph_token,
                    library=library_name,
                    ligand_count=len(chunk_bases),
                )
                chunks.append(
                    {
                        "chunk_id": chunk_id,
                        "pdb_file": str(pdb_file),
                        "pdb_id": pdb_id,
                        "variant_label": str(variant_label),
                        "ph_tag": None if ph_tag is None else str(ph_tag),
                        "run_mode": run_mode,
                        "library_name": str(library_name),
                        "library_root": str(library_root),
                        "ligand_bases": list(chunk_bases),
                        "ligand_count": int(len(chunk_bases)),
                        "weight": float(weight),
                        "assigned_task_id": int(dist_ctx.task_id),
                    }
                )

    task_ids = dist_ctx.expected_task_ids()
    if not task_ids:
        task_ids = [int(dist_ctx.task_id)]
    chunks, rechunk_stats = optimize_chunks_for_tail(
        chunks,
        run_id=run_id,
        variant_label=str(variant_label),
        task_ids=[int(tid) for tid in task_ids],
        task_count=max(1, int(dist_ctx.task_count)),
        local_workers_per_task=local_workers,
        default_target_size=int(_CHUNK_TARGET_SIZE),
        min_chunk_size=min_chunk_size,
    )
    _maybe_validate_chunk_plan(cfg, chunks)
    logging.getLogger("distributed.chunk").info(
        "[distributed.chunk.rechunk] variant=%s enabled=%s split_count=%d chunks_before=%d chunks_after=%d",
        str(variant_label).upper(),
        int(rechunk_stats.get("enabled", 0)),
        int(rechunk_stats.get("split_count", 0)),
        int(rechunk_stats.get("chunks_before", len(chunks))),
        int(rechunk_stats.get("chunks_after", len(chunks))),
    )
    planner_elapsed = time.perf_counter() - planner_started
    logging.getLogger("distributed.chunk").info(
        "[distributed.chunk.plan-cache] combos=%d unique_pdb=%d manifest_cache_hits=%d manifest_cache_misses=%d elapsed_s=%.2f",
        len(combo_items),
        len(resolved_combo_cache),
        int(cfg.get("_PLANNER_MANIFEST_CACHE_HITS", 0) or 0)
        - planner_cache_hits_before,
        int(cfg.get("_PLANNER_MANIFEST_CACHE_MISSES", 0) or 0)
        - planner_cache_misses_before,
        planner_elapsed,
    )
    return chunks


def _load_or_create_variant_chunk_plan(
    cfg: Mapping[str, Any],
    *,
    run_id: str,
    dist_ctx: Any,
    variant_label: str,
    combo_items: list[tuple[str, Optional[str]]],
    run_tokens: list[str],
    target_ligands: int,
    sample_seed: int,
) -> list[dict[str, Any]]:
    plan_path = chunk_plan_path(dist_ctx, cfg, variant_label=variant_label)
    if plan_path.exists():
        try:
            payload = json.loads(plan_path.read_text(encoding="utf-8")) or {}
            chunks = payload.get("chunks") or []
            if isinstance(chunks, list):
                compatible_chunks = [c for c in chunks if isinstance(c, dict)]
                if (
                    int(payload.get("chunk_plan_schema_version", 0) or 0)
                    >= _CHUNK_PLAN_SCHEMA_VERSION
                    and _chunk_plan_compatible(compatible_chunks)
                ):
                    _maybe_validate_chunk_plan(cfg, compatible_chunks, loaded_plan=True)
                    return compatible_chunks
                logging.warning(
                    "[distributed.chunk-plan] ignoring incompatible plan path=%s schema=%s chunks=%d",
                    plan_path,
                    payload.get("chunk_plan_schema_version"),
                    len(compatible_chunks),
                )
        except Exception:
            logging.warning(
                "[distributed.chunk-plan] failed to read existing plan path=%s",
                plan_path,
                exc_info=True,
            )

    lock_path = plan_path.with_suffix(plan_path.suffix + ".lock")
    got_lock = False
    try:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        os.close(fd)
        got_lock = True
    except FileExistsError:
        got_lock = False

    if got_lock:
        try:
            chunks = _build_lpt_combo_chunks(
                cfg,
                run_id=run_id,
                dist_ctx=dist_ctx,
                variant_label=variant_label,
                combo_items=combo_items,
                run_tokens=run_tokens,
                target_ligands=target_ligands,
                sample_seed=sample_seed,
            )
            payload = {
                "run_id": str(run_id),
                "variant_label": str(variant_label),
                "chunk_plan_schema_version": int(_CHUNK_PLAN_SCHEMA_VERSION),
                "created_at": float(time.time()),
                "task_count": int(dist_ctx.task_count),
                "chunks": chunks,
            }
            write_json_atomic(plan_path, payload)
            return chunks
        finally:
            try:
                lock_path.unlink()
            except Exception:
                pass

    deadline = time.time() + 120.0
    while time.time() < deadline:
        if plan_path.exists():
            try:
                payload = json.loads(plan_path.read_text(encoding="utf-8")) or {}
                chunks = payload.get("chunks") or []
                if isinstance(chunks, list):
                    compatible_chunks = [c for c in chunks if isinstance(c, dict)]
                    if (
                        int(payload.get("chunk_plan_schema_version", 0) or 0)
                        >= _CHUNK_PLAN_SCHEMA_VERSION
                        and _chunk_plan_compatible(compatible_chunks)
                    ):
                        _maybe_validate_chunk_plan(cfg, compatible_chunks, loaded_plan=True)
                        return compatible_chunks
            except Exception:
                pass
        time.sleep(0.5)

    chunks = _build_lpt_combo_chunks(
        cfg,
        run_id=run_id,
        dist_ctx=dist_ctx,
        variant_label=variant_label,
        combo_items=combo_items,
        run_tokens=run_tokens,
        target_ligands=target_ligands,
        sample_seed=sample_seed,
    )
    payload = {
        "run_id": str(run_id),
        "variant_label": str(variant_label),
        "chunk_plan_schema_version": int(_CHUNK_PLAN_SCHEMA_VERSION),
        "created_at": float(time.time()),
        "task_count": int(dist_ctx.task_count),
        "chunks": chunks,
    }
    write_json_atomic(plan_path, payload)
    return chunks
