"""Run-manifest lifecycle update entrypoints (init, config, load, finalize)."""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Dict, Mapping, MutableMapping, Optional, Sequence

from cli.run_manifest_distributed import reconcile_distributed_manifest_state
from cli.run_manifest_io import (
    get_manifest_paths,
    load_manifest as _load_manifest,
    write_manifest as _write_manifest,
)
from cli.run_manifest_support import (
    coerce_bool_token as _coerce_bool_token,
    compute_config_hash,
    distributed_enabled as _distributed_enabled,
    distributed_mode as _distributed_mode,
    distributed_task_count as _distributed_task_count,
    distributed_task_id as _distributed_task_id,
    git_info as _git_info,
    normalize_pdb_id_token as _normalize_pdb_id_token,
    normalize_string_token as _normalize_string_token,
    refresh_summary as _refresh_summary,
    require_yaml as _require_yaml,
    resources_snapshot as _resources_snapshot,
    utc_now_iso as _utc_now_iso,
)

def init_run_manifest(
    cfg: Mapping[str, Any], run_id: str, argv: list[str], log_path: str
) -> None:
    try:
        manifest_dir, manifest_path = get_manifest_paths(cfg, run_id)
        yaml_mod = _require_yaml()
        if yaml_mod is None:
            return

        run_dir = Path(cfg.get("CONFIG_RUN_DIR", ""))
        overall_dir = Path(cfg.get("OVERALL_DIR", "."))
        selection_mode = None
        test_mode = str(cfg.get("TEST_MODE_ENABLE", "off"))

        if test_mode != "off":
            selection_mode = "TEST_LIBRARY_MAP"
        elif cfg.get("_EFFECTIVE_SPECIFIED_PROTEINS"):
            selection_mode = "SPECIFIED_PROTEINS"
        else:
            selection_mode = "AUTO"

        pdb_list = None
        if selection_mode == "SPECIFIED_PROTEINS":
            pdb_list = list(cfg.get("_EFFECTIVE_SPECIFIED_PROTEINS", []) or [])

        manifest: Dict[str, Any] = {
            "run_id": str(run_id),
            "status": "running",
            "command": {
                "argv": " ".join(argv),
                "pdb_selection_mode": selection_mode,
                "TEST_MODE_ENABLE": test_mode,
                "PH_ENSEMBLE_ENABLE": bool(cfg.get("PH_ENSEMBLE", False)),
            },
            "git": _git_info(overall_dir),
            "paths": {
                "run_dir": str(run_dir) if run_dir else None,
                "log_file": str(log_path),
                "docked_dir": str(
                    cfg.get("DOCKED_DIR", overall_dir / "outputs" / "docked")
                ),
                "input_pdb_dir": str(cfg.get("INPUT_DIR", overall_dir / "input_pdbs")),
                "processed_pdb_dir": str(
                    cfg.get(
                        "OUTPUT_DIR",
                        overall_dir / "outputs" / "processed_pdbs" / run_id,
                    )
                ),
                "prepped_ligands_dir": str(
                    cfg.get("PREPPED_LIGANDS_DIR", overall_dir / "prepped_ligands")
                ),
                "config_file": "run_config.yaml",
            },
            "timing": {
                "created_at": _utc_now_iso(),
                "started_at": _utc_now_iso(),
                "finished_at": None,
                "wall_time_sec": None,
            },
            "resources": _resources_snapshot(),
            "summary": {
                "total_proteins_scheduled": 0,
                "total_proteins_completed": 0,
                "total_proteins_failed": 0,
            },
            "proteins": {},
        }

        if pdb_list:
            manifest["command"]["pdb_list"] = pdb_list

        if _distributed_enabled(cfg):
            manifest["command"]["DISTRIBUTED_MODE"] = _distributed_mode(cfg)
            manifest["command"]["DISTRIBUTED_TASK_COUNT"] = _distributed_task_count()
            manifest["command"]["DISTRIBUTED_TASK_ID"] = _distributed_task_id()

        manifest["command"]["config_hash"] = compute_config_hash(cfg)

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to initialize run manifest run_id=%s",
            run_id,
            exc_info=True,
        )

def update_manifest_for_scheduled_proteins(
    cfg: Mapping[str, Any],
    run_id: str,
    scheduled_pdb_ids: Sequence[Any],
) -> None:
    """
    Record the canonical list of PDB IDs queued for this run.
    """
    try:
        if not run_id:
            logging.debug("[run-manifest.scheduled.skip] reason=missing_run_id")
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.scheduled.skip] manifest missing run_id=%s path=%s",
                run_id,
                manifest_path,
            )
            return

        summary = manifest.get("summary")
        if not isinstance(summary, MutableMapping):
            summary = {}
            manifest["summary"] = summary

        norm_ids: set[str] = set()
        for raw in scheduled_pdb_ids:
            nid = _normalize_pdb_id_token(raw)
            if nid:
                norm_ids.add(nid)

        sorted_ids = sorted(norm_ids)
        summary["total_proteins_scheduled"] = len(sorted_ids)
        summary["total_protein_list"] = sorted_ids

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.scheduled.error] run_id=%s",
            run_id,
            exc_info=True,
        )


def update_manifest_for_run_config(
    cfg: Mapping[str, Any],
    run_id: str,
) -> None:
    """
    Record resolved apo/holo mode and water-policy knobs in the manifest.
    """
    try:
        if not run_id:
            logging.debug("[run-manifest.run-config.skip] reason=missing_run_id")
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.run-config.skip] manifest missing run_id=%s path=%s",
                run_id,
                manifest_path,
            )
            return

        cmd = manifest.get("command")
        if not isinstance(cmd, MutableMapping):
            cmd = {}
            manifest["command"] = cmd

        apo_mode = _normalize_string_token(cfg.get("_RESOLVED_APO_HOLO_MODE"))
        if apo_mode is None:
            apo_mode = _normalize_string_token(cfg.get("APO_HOLO_MODE"))
        if apo_mode is not None:
            cmd["APO_HOLO_MODE"] = apo_mode

        if "REMOVE_WATERS" in cfg:
            rw_raw = cfg.get("REMOVE_WATERS")
            rw_val = _coerce_bool_token(rw_raw)
            cmd["REMOVE_WATERS"] = rw_val if rw_val is not None else rw_raw

        if "KEEP_WATERS_WITHIN_A" in cfg:
            val = cfg.get("KEEP_WATERS_WITHIN_A")
            try:
                cmd["KEEP_WATERS_WITHIN_A"] = float(val) if val is not None else None
            except Exception:
                cmd["KEEP_WATERS_WITHIN_A"] = val

        if "WATER_KEEP_POLICY" in cfg:
            val = _normalize_string_token(cfg.get("WATER_KEEP_POLICY"))
            cmd["WATER_KEEP_POLICY"] = (
                val if val is not None else cfg.get("WATER_KEEP_POLICY")
            )

        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.run-config.error] run_id=%s",
            run_id,
            exc_info=True,
        )


def update_manifest_for_config_hash(
    cfg: Mapping[str, Any],
    run_id: str,
) -> None:
    """
    Record the resolved config hash under command.config_hash.
    """
    try:
        if not run_id:
            logging.debug("[run-manifest.config-hash.skip] missing run_id")
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.config-hash.skip] manifest missing run_id=%s path=%s",
                run_id,
                manifest_path,
            )
            return

        cmd = manifest.get("command")
        if not isinstance(cmd, MutableMapping):
            cmd = {}
            manifest["command"] = cmd

        cmd["config_hash"] = compute_config_hash(cfg)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.config-hash.error] run_id=%s",
            run_id,
            exc_info=True,
        )


def load_run_manifest(cfg: Mapping[str, Any], run_id: str) -> Optional[Dict[str, Any]]:
    """
    Best-effort loader for run_manifest.yaml for a given run_id.

    Returns the manifest dict, or None if it can't be loaded.
    """
    try:
        if not run_id:
            logging.debug("[run-manifest.load.skip] reason=missing_run_id")
            return None

        reconcile_stats = reconcile_distributed_manifest_state(cfg, run_id, force=False)
        if reconcile_stats.get("action") in {"reduced", "reduce_noop"}:
            logging.info(
                "[run-manifest.load.reconcile] run_id=%s action=%s state_files=%s merged_entries=%s stale=%s",
                run_id,
                reconcile_stats.get("action"),
                reconcile_stats.get("state_files"),
                reconcile_stats.get("merged_entries"),
                reconcile_stats.get("manifest_stale"),
            )

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.load.missing] run_id=%s path=%s",
                run_id,
                manifest_path,
            )
        return manifest
    except Exception:
        logging.warning(
            "[run-manifest.load.error] run_id=%s",
            run_id,
            exc_info=True,
        )
        return None
def finalize_run_manifest(
    cfg: Mapping[str, Any],
    run_id: str,
    start_time: float,
    failed_entries: list[tuple[str, str, str, str, str]],
) -> None:
    try:
        reconcile_stats = reconcile_distributed_manifest_state(cfg, run_id, force=False)
        if reconcile_stats.get("action") in {"reduced", "reduce_noop"}:
            logging.info(
                "[run-manifest.finalize.reconcile] run_id=%s action=%s state_files=%s merged_entries=%s stale=%s",
                run_id,
                reconcile_stats.get("action"),
                reconcile_stats.get("state_files"),
                reconcile_stats.get("merged_entries"),
                reconcile_stats.get("manifest_stale"),
            )

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            return

        manifest_status = "completed" if len(failed_entries) == 0 else "failed"
        manifest["status"] = manifest_status

        timing = manifest.setdefault("timing", {})
        if isinstance(timing, MutableMapping):
            timing.setdefault("created_at", _utc_now_iso())
            timing.setdefault("started_at", _utc_now_iso())
            timing["finished_at"] = _utc_now_iso()
            timing["wall_time_sec"] = round(time.time() - float(start_time), 3)
        else:
            manifest["timing"] = {
                "created_at": _utc_now_iso(),
                "started_at": _utc_now_iso(),
                "finished_at": _utc_now_iso(),
                "wall_time_sec": round(time.time() - float(start_time), 3),
            }

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to finalize manifest run_id=%s",
            run_id,
            exc_info=True,
        )
