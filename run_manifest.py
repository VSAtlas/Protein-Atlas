"""
Best-effort run manifest writer for Atlas runs.

All helpers must remain non-fatal: any error (missing PyYAML, IO, git)
is logged at WARNING level and silently skipped so scientific behavior
is unchanged.
"""

from __future__ import annotations

import datetime
import logging
import os
import platform
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, Mapping, MutableMapping, Optional, Tuple

try:  # PyYAML is optional; we degrade gracefully if unavailable
    import yaml
except Exception:  # pragma: no cover - import guard
    yaml = None  # type: ignore[assignment]


STAGE_KEYS = ("prep", "pocket_detection", "docking", "postprocessing")
_missing_yaml_logged = False


def _utc_now_iso() -> str:
    """UTC timestamp without microseconds for stable manifest entries."""

    return datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def _require_yaml() -> Optional[Any]:
    """Return the yaml module or log once and skip writes if missing."""

    global _missing_yaml_logged
    if yaml is None:
        if not _missing_yaml_logged:
            logging.warning(
                "[run-manifest] PyYAML not available; manifest writes disabled"
            )
            _missing_yaml_logged = True
        return None
    return yaml


def _default_stage_entry() -> Dict[str, Any]:
    return {
        "status": "pending",
        "error": None,
        "timing": {"started_at": None, "finished_at": None, "wall_time_sec": None},
    }


def _default_protein_entry() -> Dict[str, Any]:
    return {
        "pdb_id": None,
        "library": None,
        "variant": None,
        "ph": None,
        "status": "pending",
        "error": None,
        "timing": {"started_at": None, "finished_at": None, "wall_time_sec": None},
        "stages": {k: _default_stage_entry() for k in STAGE_KEYS},
    }


def _refresh_summary(manifest: MutableMapping[str, Any]) -> None:
    proteins = manifest.get("proteins")
    if not isinstance(proteins, Mapping):
        manifest["summary"] = {
            "total_proteins_scheduled": 0,
            "total_proteins_completed": 0,
            "total_proteins_failed": 0,
        }
        return

    # Best-effort dedupe: if any pH-tagged entry exists for (pdb,variant),
    # drop the corresponding |base entry so summaries align with pH ensembles.
    try:
        ph_pairs: set[tuple[str, str]] = set()
        for key in proteins:
            parts = str(key).split("|", 2)
            if len(parts) != 3:
                continue
            pdb_part, variant_part, ph_token = parts
            if ph_token and ph_token != "base":
                ph_pairs.add((pdb_part, variant_part))

        base_keys_to_drop: list[str] = []
        for k, _ in list(proteins.items()):
            parts = str(k).split("|", 2)
            if len(parts) != 3:
                continue
            pdb_part, variant_part, ph_token = parts
            if ph_token == "base" and (pdb_part, variant_part) in ph_pairs:
                base_keys_to_drop.append(k)

        for k in base_keys_to_drop:
            proteins.pop(k, None)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to dedupe proteins for summary", exc_info=True
        )

    total = len(proteins)
    completed = 0
    failed = 0
    for entry in proteins.values():
        status = None
        if isinstance(entry, Mapping):
            status = entry.get("status")
        if status == "completed":
            completed += 1
        elif status == "failed":
            failed += 1

    manifest["summary"] = {
        "total_proteins_scheduled": total,
        "total_proteins_completed": completed,
        "total_proteins_failed": failed,
    }


def _load_manifest(manifest_path: Path) -> Optional[Dict[str, Any]]:
    yaml_mod = _require_yaml()
    if yaml_mod is None:
        return None

    if not manifest_path.exists():
        return {}

    try:
        with manifest_path.open("r", encoding="utf-8") as fh:
            data = yaml_mod.safe_load(fh) or {}
        if not isinstance(data, dict):
            return {}
        return data
    except Exception:
        logging.warning(
            "[run-manifest] Failed to load manifest path=%s", manifest_path, exc_info=True
        )
        return {}


def _write_manifest(manifest_path: Path, manifest: Dict[str, Any]) -> None:
    yaml_mod = _require_yaml()
    if yaml_mod is None:
        return

    try:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = manifest_path.with_suffix(manifest_path.suffix + ".tmp")
        with tmp_path.open("w", encoding="utf-8") as fh:
            yaml_mod.safe_dump(manifest, fh, default_flow_style=False, sort_keys=False)
        tmp_path.replace(manifest_path)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to write manifest path=%s", manifest_path, exc_info=True
        )


def _protein_key(pdb_id: str, variant_label: Optional[str], ph_tag: Optional[str]) -> str:
    variant_norm = (variant_label or "legacy").strip().upper() or "LEGACY"
    ph_norm = (ph_tag or "base").strip()
    ph_token = ph_norm if ph_norm else "base"
    return f"{str(pdb_id).upper()}|{variant_norm}|{ph_token}"


def get_manifest_paths(cfg: Mapping[str, Any], run_id: str) -> Tuple[Path, Path]:
    """
    Return (manifest_dir, manifest_path) for the current run and ensure the
    directory exists. Never raises.
    """

    try:
        repo_root = Path(cfg.get("OVERALL_DIR", "."))
        manifest_dir = repo_root / "manifests" / str(run_id)
        manifest_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = manifest_dir / "run_manifest.yaml"
        return manifest_dir, manifest_path
    except Exception:
        logging.warning(
            "[run-manifest] Failed to prepare manifest paths run_id=%s", run_id, exc_info=True
        )
        fallback_dir = Path(cfg.get("OVERALL_DIR", ".")) / "manifests"
        return fallback_dir, fallback_dir / "run_manifest.yaml"


def _git_info(repo_root: Path) -> Dict[str, Any]:
    info = {"repo": None, "branch": None, "commit": None, "dirty": None}
    try:
        def _run(cmd: list[str]) -> Optional[str]:
            try:
                res = subprocess.run(
                    cmd,
                    cwd=str(repo_root),
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if res.returncode == 0:
                    out = res.stdout.strip()
                    return out or None
            except Exception:
                return None
            return None

        info["repo"] = _run(["git", "config", "--get", "remote.origin.url"])
        info["branch"] = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
        info["commit"] = _run(["git", "rev-parse", "HEAD"])

        dirty_output = _run(["git", "status", "--porcelain"])
        info["dirty"] = bool(dirty_output) if dirty_output is not None else None
    except Exception:
        logging.warning(
            "[run-manifest] Failed to gather git metadata root=%s", repo_root, exc_info=True
        )
    return info


def _resources_snapshot() -> Dict[str, Any]:
    host = None
    try:
        host = os.uname().nodename
    except Exception:
        try:
            host = platform.node()
        except Exception:
            host = None

    n_cores = None
    try:
        n_cores = os.cpu_count()
    except Exception:
        n_cores = None

    ram_gb = None
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        phys_pages = os.sysconf("SC_PHYS_PAGES")
        ram_bytes = float(page_size) * float(phys_pages)
        ram_gb = round(ram_bytes / (1024 ** 3), 2)
    except Exception:
        ram_gb = None

    gpu = None
    try:
        cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES")
        if cuda_visible:
            gpu = cuda_visible
    except Exception:
        gpu = None

    return {"host": host, "n_cores": n_cores, "ram_gb": ram_gb, "gpu": gpu}


def init_run_manifest(cfg: Mapping[str, Any], run_id: str, argv: list[str], log_path: str) -> None:
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

        manifest = {
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
                "docked_dir": str(cfg.get("DOCKED_DIR", overall_dir / "docked")),
                "input_pdb_dir": str(cfg.get("INPUT_DIR", overall_dir / "input_pdbs")),
                "processed_pdb_dir": str(cfg.get("OUTPUT_DIR", overall_dir / "processed_pdbs")),
                "prepped_ligands_dir": str(cfg.get("PREPPED_LIGANDS_DIR", overall_dir / "prepped_ligands")),
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

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to initialize run manifest run_id=%s", run_id, exc_info=True
        )


def _ensure_protein(
    manifest: MutableMapping[str, Any],
    pdb_id: str,
    variant_label: Optional[str],
    ph_tag: Optional[str],
) -> Dict[str, Any]:
    proteins = manifest.setdefault("proteins", {})
    if not isinstance(proteins, MutableMapping):
        manifest["proteins"] = {}
        proteins = manifest["proteins"]

    key = _protein_key(pdb_id, variant_label, ph_tag)
    entry = proteins.get(key)
    if not isinstance(entry, MutableMapping):
        entry = _default_protein_entry()
        proteins[key] = entry
    else:
        # Ensure stages/timing keys exist
        entry.setdefault("pdb_id", str(pdb_id))
        entry.setdefault("library", None)
        entry.setdefault("variant", None)
        entry.setdefault("ph", ph_tag)
        entry.setdefault("status", "pending")
        entry.setdefault("error", None)
        entry.setdefault(
            "timing", {"started_at": None, "finished_at": None, "wall_time_sec": None}
        )
        stages = entry.setdefault("stages", {})
        for stage_key in STAGE_KEYS:
            stage_entry = stages.get(stage_key)
            if not isinstance(stage_entry, MutableMapping):
                stages[stage_key] = _default_stage_entry()
                continue
            stage_entry.setdefault("status", "pending")
            stage_entry.setdefault("error", None)
            stage_entry.setdefault(
                "timing",
                {"started_at": None, "finished_at": None, "wall_time_sec": None},
            )

    entry["pdb_id"] = str(pdb_id).upper()
    entry["variant"] = (variant_label or "legacy").strip().upper() or "LEGACY"
    entry["ph"] = ph_tag

    return entry  # type: ignore[return-value]


def update_manifest_for_protein_start(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: str,
    library: Optional[str],
    *,
    ph_tag: Optional[str] = None,
) -> None:
    try:
        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
        entry["variant"] = (variant_label or "legacy").upper()
        if library is None:
            library = cfg.get("LIBRARY_SUBDIR_DEFAULT")
        entry["library"] = library
        entry["status"] = "running"
        entry.setdefault("ph", ph_tag)
        timing = entry.get("timing", {})
        if isinstance(timing, MutableMapping):
            if not timing.get("started_at"):
                timing["started_at"] = _utc_now_iso()
        else:
            entry["timing"] = {"started_at": _utc_now_iso(), "finished_at": None, "wall_time_sec": None}

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record start pdb=%s variant=%s ph=%s run_id=%s",
            pdb_id,
            variant_label,
            ph_tag if ph_tag is not None else "base",
            run_id,
            exc_info=True,
        )


def update_manifest_for_protein_success(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: str,
    elapsed_sec: float,
    *,
    ph_tag: Optional[str] = None,
) -> None:
    try:
        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
        entry["status"] = "completed"
        entry["error"] = None

        timing = entry.get("timing", {})
        if isinstance(timing, MutableMapping):
            timing.setdefault("started_at", _utc_now_iso())
            timing["finished_at"] = _utc_now_iso()
            timing["wall_time_sec"] = round(float(elapsed_sec), 3)
        else:
            entry["timing"] = {
                "started_at": _utc_now_iso(),
                "finished_at": _utc_now_iso(),
                "wall_time_sec": round(float(elapsed_sec), 3),
            }

        stages = entry.get("stages", {})
        if isinstance(stages, MutableMapping):
            for stage_key in STAGE_KEYS:
                stage_entry = stages.get(stage_key)
                if not isinstance(stage_entry, MutableMapping):
                    stages[stage_key] = _default_stage_entry()
                    stage_entry = stages[stage_key]
                stage_entry["status"] = "completed"
                stage_entry.setdefault("error", None)
        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record success pdb=%s variant=%s ph=%s run_id=%s",
            pdb_id,
            variant_label,
            ph_tag if ph_tag is not None else "base",
            run_id,
            exc_info=True,
        )


def _extract_error_from_fail_log(fail_log_path: Optional[Path]) -> str:
    try:
        if not fail_log_path or not fail_log_path.exists():
            return "<unknown error>"
        with fail_log_path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                stripped = line.lstrip()
                if stripped.startswith("exception     ="):
                    parts = stripped.split("=", 1)
                    if len(parts) == 2:
                        return parts[1].strip()
        return "<unknown error>"
    except Exception:
        return "<unknown error>"


def update_manifest_for_protein_failure(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: str,
    fail_log_path: Optional[Path],
    *,
    ph_tag: Optional[str] = None,
) -> None:
    try:
        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
        entry["variant"] = (variant_label or "legacy").upper()
        entry["status"] = "failed"
        entry["error"] = _extract_error_from_fail_log(fail_log_path)
        entry["fail_log"] = str(fail_log_path) if fail_log_path is not None else None

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record failure pdb=%s variant=%s ph=%s run_id=%s",
            pdb_id,
            variant_label,
            ph_tag if ph_tag is not None else "base",
            run_id,
            exc_info=True,
        )


def finalize_run_manifest(
    cfg: Mapping[str, Any], run_id: str, start_time: float, failed_entries: list[tuple[str, str, str, str, str]]
) -> None:
    try:
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
            "[run-manifest] Failed to finalize manifest run_id=%s", run_id, exc_info=True
        )
