"""Distributed per-protein manifest state helpers."""

from __future__ import annotations

import json
import logging
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, Mapping, MutableMapping, Optional

from cli.run_manifest_io import (
    get_manifest_paths,
    load_manifest as _load_manifest,
    protein_key as _protein_key,
    write_manifest as _write_manifest,
)
from cli.run_manifest_mutations import ensure_protein as _ensure_protein
from cli.run_manifest_support import (
    distributed_enabled as _distributed_enabled,
    distributed_mode as _distributed_mode,
    distributed_protein_state_relpath,
    distributed_task_count as _distributed_task_count,
    distributed_task_id as _distributed_task_id,
    exclusive_file_lock as _exclusive_file_lock,
    file_sha1 as _file_sha1,
    latest_mtime as _latest_mtime,
    manifest_has_proteins as _manifest_has_proteins,
    parse_int as _parse_int,
    refresh_summary as _refresh_summary,
    suppress_distributed_base_write as _suppress_distributed_base_write,
    utc_now_iso as _utc_now_iso,
)

def distributed_root(cfg: Mapping[str, Any], run_id: str) -> Path:
    manifest_dir, _ = get_manifest_paths(cfg, run_id)
    return manifest_dir / "distributed"


def distributed_protein_state_path(
    cfg: Mapping[str, Any],
    run_id: str,
    protein_key: str,
) -> Path:
    root = distributed_root(cfg, run_id) / "proteins"
    root.mkdir(parents=True, exist_ok=True)
    return root / distributed_protein_state_relpath(protein_key)


def distributed_protein_state_files(cfg: Mapping[str, Any], run_id: str) -> list[Path]:
    if not run_id:
        return []
    protein_dir = distributed_root(cfg, run_id) / "proteins"
    if not protein_dir.exists():
        return []
    return [path for path in sorted(protein_dir.glob("*.json")) if path.is_file()]


def distributed_mutate_protein_entry(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    ph_tag: Optional[str],
    mutator: Any,
) -> bool:
    if not run_id or not _distributed_enabled(cfg):
        return False
    if _suppress_distributed_base_write(cfg, ph_tag=ph_tag):
        logging.info(
            "[run-manifest.distributed.base.skip] run_id=%s pdb=%s variant=%s reason=ph_ensemble_distributed",
            run_id,
            pdb_id,
            variant_label,
        )
        return True
    try:
        protein_key = _protein_key(pdb_id, variant_label, ph_tag)
        state_path = distributed_protein_state_path(cfg, run_id, protein_key)
        lock_path = state_path.with_suffix(state_path.suffix + ".lock")
        with _exclusive_file_lock(lock_path):
            payload: Dict[str, Any] = {}
            if state_path.exists():
                try:
                    payload = json.loads(state_path.read_text(encoding="utf-8")) or {}
                except Exception:
                    payload = {}

            entry = payload.get("entry")
            manifest_stub: Dict[str, Any] = {"proteins": {}}
            if isinstance(entry, Mapping):
                manifest_stub["proteins"][protein_key] = dict(entry)

            protein_entry = _ensure_protein(manifest_stub, pdb_id, variant_label, ph_tag)
            mutator(protein_entry)

            payload = {
                "run_id": str(run_id),
                "protein_key": protein_key,
                "task_id": int(_distributed_task_id()),
                "updated_at": _utc_now_iso(),
                "entry": protein_entry,
            }

            tmp_path: Optional[Path] = None
            try:
                with tempfile.NamedTemporaryFile(
                    mode="w",
                    encoding="utf-8",
                    prefix=state_path.name + ".",
                    suffix=".tmp",
                    dir=state_path.parent,
                    delete=False,
                ) as handle:
                    tmp_path = Path(handle.name)
                    json.dump(payload, handle, sort_keys=True)
                    handle.flush()
                    os.fsync(handle.fileno())
                tmp_path.replace(state_path)
            finally:
                try:
                    if tmp_path and tmp_path.exists():
                        tmp_path.unlink()
                except Exception:
                    pass
        return True
    except Exception:
        logging.warning(
            "[run-manifest.distributed.entry.error] run_id=%s pdb=%s variant=%s ph=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            exc_info=True,
        )
        return False


def reduce_distributed_protein_states(cfg: Mapping[str, Any], run_id: str) -> int:
    """
    Merge per-protein distributed state files into run_manifest.yaml.
    """
    if not run_id:
        return 0
    try:
        state_files = distributed_protein_state_files(cfg, run_id)
        if not state_files:
            return 0

        merged: Dict[str, Dict[str, Any]] = {}
        observed_task_ids: set[int] = set()
        for state_path in state_files:
            try:
                payload = json.loads(state_path.read_text(encoding="utf-8")) or {}
            except Exception:
                continue
            key = str(payload.get("protein_key", "") or "").strip()
            entry = payload.get("entry")
            if not key or not isinstance(entry, Mapping):
                continue
            task_id = _parse_int(payload.get("task_id"))
            if task_id is not None:
                observed_task_ids.add(int(task_id))
            merged[key] = dict(entry)

        if not merged:
            return 0

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            return 0
        before_serialized = json.dumps(manifest, sort_keys=True, default=str)

        proteins = manifest.get("proteins")
        if not isinstance(proteins, MutableMapping):
            proteins = {}
            manifest["proteins"] = proteins

        for key, entry in merged.items():
            proteins[key] = entry

        command = manifest.get("command")
        if not isinstance(command, MutableMapping):
            command = {}
            manifest["command"] = command

        existing_mode = str(command.get("DISTRIBUTED_MODE", "") or "").strip()
        resolved_mode = _distributed_mode(cfg)
        if resolved_mode == "off":
            resolved_mode = existing_mode or "slurm_array"
        command["DISTRIBUTED_MODE"] = resolved_mode

        env_task_count = _distributed_task_count()
        existing_task_count = _parse_int(command.get("DISTRIBUTED_TASK_COUNT"))
        if env_task_count > 1:
            resolved_task_count = env_task_count
        elif existing_task_count is not None and existing_task_count > 1:
            resolved_task_count = existing_task_count
        elif observed_task_ids:
            low = min(observed_task_ids)
            high = max(observed_task_ids)
            resolved_task_count = max(1, high - low + 1, len(observed_task_ids))
        else:
            resolved_task_count = 1
        command["DISTRIBUTED_TASK_COUNT"] = int(resolved_task_count)

        _refresh_summary(manifest)
        after_serialized = json.dumps(manifest, sort_keys=True, default=str)
        if after_serialized != before_serialized:
            _write_manifest(manifest_path, manifest)
        return len(merged)
    except Exception:
        logging.warning(
            "[run-manifest.distributed.reduce.error] run_id=%s",
            run_id,
            exc_info=True,
        )
        return 0


def reconcile_distributed_manifest_state(
    cfg: Mapping[str, Any],
    run_id: str,
    *,
    force: bool = False,
) -> Dict[str, Any]:
    """
    Best-effort reconciliation between distributed per-protein state files and
    canonical run_manifest.yaml.

    Returns a diagnostic dictionary describing action and counts.
    """
    stats: Dict[str, Any] = {
        "run_id": str(run_id or ""),
        "force": bool(force),
        "action": "skip",
        "state_files": 0,
        "merged_entries": 0,
        "manifest_exists": False,
        "manifest_has_proteins": False,
        "manifest_stale": False,
        "manifest_path": None,
    }
    try:
        if not run_id:
            stats["action"] = "skip_missing_run_id"
            return stats

        state_files = distributed_protein_state_files(cfg, run_id)
        stats["state_files"] = len(state_files)
        if not state_files:
            stats["action"] = "skip_no_distributed_state"
            return stats

        _, manifest_path = get_manifest_paths(cfg, run_id)
        stats["manifest_path"] = str(manifest_path)
        stats["manifest_exists"] = bool(manifest_path.exists())
        manifest_mtime = 0.0
        if manifest_path.exists():
            try:
                manifest_mtime = float(manifest_path.stat().st_mtime)
            except Exception:
                manifest_mtime = 0.0

        manifest = _load_manifest(manifest_path)
        has_proteins = _manifest_has_proteins(manifest)
        stats["manifest_has_proteins"] = bool(has_proteins)
        latest_state_mtime = _latest_mtime(state_files)
        stale = (
            (not manifest_path.exists())
            or (latest_state_mtime > (manifest_mtime + 1e-6))
            or (not has_proteins)
        )
        stats["manifest_stale"] = bool(stale)

        if not force and not stale:
            stats["action"] = "skip_up_to_date"
            return stats

        before_digest = _file_sha1(manifest_path)
        merged = reduce_distributed_protein_states(cfg, run_id)
        after_digest = _file_sha1(manifest_path)
        stats["merged_entries"] = int(merged)
        if merged <= 0 and before_digest == after_digest:
            stats["action"] = "reduce_noop"
        elif before_digest == after_digest and before_digest is not None:
            stats["action"] = "reduce_noop"
        else:
            stats["action"] = "reduced"
        return stats
    except Exception:
        logging.warning(
            "[run-manifest.distributed.reconcile.error] run_id=%s",
            run_id,
            exc_info=True,
        )
        stats["action"] = "error"
        return stats
