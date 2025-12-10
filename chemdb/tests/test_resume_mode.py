from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from input_and_export_functions import load_inputs  # noqa: E402
from run_manifest import compute_config_hash, get_manifest_paths  # noqa: E402
from run_manifest import init_run_manifest  # noqa: E402
from path_router import load_ph_tags  # noqa: E402


MAIN_PY = REPO_ROOT / "main.py"
CONFIG_PATH = REPO_ROOT / "config.txt"
RUN_ID = "test_resume_mode"


def _load_manifest(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text())


def _ensure_stage_entry() -> dict[str, Any]:
    return {
        "status": "pending",
        "error": None,
        "timing": {"started_at": None, "finished_at": None, "wall_time_sec": None},
        "details": {},
    }


def _ensure_protein_entry(pdb_id: str, variant: str, ph: str | None = None) -> dict[str, Any]:
    return {
        "pdb_id": pdb_id,
        "library": "test_library_10",
        "variant": variant,
        "ph": ph,
        "status": "pending",
        "error": None,
        "timing": {"started_at": None, "finished_at": None, "wall_time_sec": None},
        "stages": {
            "prep": _ensure_stage_entry(),
            "pocket_detection": _ensure_stage_entry(),
            "docking": _ensure_stage_entry(),
            "postprocessing": _ensure_stage_entry(),
        },
    }


def _run_python(args: list[str], env: Mapping[str, str], *, check: bool = True) -> subprocess.CompletedProcess:
    cmd = ["micromamba", "run", "-n", "docking-env", "python"] + args
    return subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        env=dict(env),
        text=True,
        capture_output=True,
        check=check,
    )


def test_resume_mode(tmp_path: Path) -> None:
    """
    Simulate an interrupted run and verify resume mode reloads the snapshot config.
    """
    assert MAIN_PY.exists(), "main.py not found; cannot run integration test"
    cfg = load_inputs()
    cfg["RUN_ID"] = RUN_ID
    cfg["CONFIG_RUN_DIR"] = str(REPO_ROOT / "configs" / RUN_ID)
    cfg["TEST_MODE_ENABLE"] = "dud+fda"
    cfg["LIBRARY_SUBDIR_DEFAULT"] = cfg.get("TEST_FDA_LIBRARY_SUBDIR", "fda_test_library_10")
    cfg["TEST_LIBRARY_MAP"] = {"TEMP": "test_library_10", "T3MP": "test_library_10"}
    cfg["FAST_MODE"] = True
    cfg["NO_LIBRARY_DOCKING"] = True
    _, manifest_path = get_manifest_paths(cfg, RUN_ID)
    manifest_dir = manifest_path.parent
    manifest_dir.mkdir(parents=True, exist_ok=True)

    # Ensure snapshot exists; if not, run the initial command once.
    resume_env = os.environ.copy()
    resume_env.update(
        {
            "PYTHONUNBUFFERED": "1",
            "CPU": "1",
            "MAX_PARALLEL_JOBS": "1",
            "NO_LIBRARY_DOCKING": "1",
            "TEST_MODE_ENABLE": cfg["TEST_MODE_ENABLE"],
            "FAST_MODE": "1",
            "LIBRARY_SUBDIR_DEFAULT": cfg["LIBRARY_SUBDIR_DEFAULT"],
        }
    )
    initial_cmd = [
        str(MAIN_PY),
        "-test",
        "temp",
        "t3mp",
        "-fast",
        "-test-fda",
        "--run-id",
        RUN_ID,
    ]

    # Avoid expensive reruns if manifest and snapshot already exist.
    baseline_manifest_exists = manifest_path.exists()
    if baseline_manifest_exists:
        manifest = _load_manifest(manifest_path) or {}
        paths = manifest.get("paths") or {}
        run_dir = Path(paths.get("run_dir", REPO_ROOT / "configs" / RUN_ID))
        cfg_file = paths.get("config_file", "run_config.yaml")
        cfg_snap_path = run_dir / cfg_file
    else:
        manifest = {}
        run_dir = REPO_ROOT / "configs" / RUN_ID
        cfg_snap_path = run_dir / "run_config.yaml"
        cfg_file = "run_config.yaml"

    if not cfg_snap_path.exists():
        run_dir.mkdir(parents=True, exist_ok=True)
        if not manifest_path.exists():
            log_stub = REPO_ROOT / "logs" / f"{RUN_ID}_bootstrap.log"
            log_stub.parent.mkdir(parents=True, exist_ok=True)
            init_run_manifest(cfg, RUN_ID, initial_cmd, str(log_stub))
        if not cfg_snap_path.exists():
            cfg_snap_path.parent.mkdir(parents=True, exist_ok=True)
            cfg_snap_path.write_text(yaml.safe_dump(cfg))
        manifest = _load_manifest(manifest_path) if manifest_path.exists() else {}
        paths = manifest.get("paths") or {}
        run_dir = Path(paths.get("run_dir", run_dir))
        cfg_file = paths.get("config_file", "run_config.yaml")
        cfg_snap_path = run_dir / cfg_file

    cfg_snap_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_snap_path.write_text(yaml.safe_dump(cfg))

    # Build interrupted manifest state and save a backup copy of that state.
    manifest = _load_manifest(manifest_path) or {}
    if not manifest:
        log_stub = REPO_ROOT / "logs" / f"{RUN_ID}_bootstrap.log"
        log_stub.parent.mkdir(parents=True, exist_ok=True)
        init_run_manifest(cfg, RUN_ID, initial_cmd, str(log_stub))
        manifest = _load_manifest(manifest_path) or {}
    paths = manifest.get("paths") or {}
    run_dir = Path(paths.get("run_dir", run_dir))
    cfg_file = paths.get("config_file", "run_config.yaml")
    cfg_snap_path = run_dir / cfg_file

    ph_tags = load_ph_tags("T3MP", variant="HOLO")
    ph_token = ph_tags[0] if ph_tags else "base"

    proteins: dict[str, Any] = {}
    for key, entry in (manifest.get("proteins") or {}).items():
        try:
            pdb_part, variant_part, ph_part = str(key).split("|", 2)
        except ValueError:
            continue
        if pdb_part.upper() in {"TEMP", "T3MP"} and variant_part:
            if pdb_part.upper() == "TEMP":
                entry = dict(entry)
                entry["status"] = "completed"
                proteins[f"TEMP|{variant_part}|{ph_part}"] = entry
            elif pdb_part.upper() == "T3MP":
                entry = dict(entry)
                entry["status"] = "pending"
                proteins[f"T3MP|{variant_part}|{ph_part}"] = entry

    if not proteins:
        variant = "HOLO"
        proteins[f"TEMP|{variant}|{ph_token}"] = _ensure_protein_entry("TEMP", variant, ph_token)
        proteins[f"TEMP|{variant}|{ph_token}"]["status"] = "completed"
        proteins[f"T3MP|{variant}|{ph_token}"] = _ensure_protein_entry("T3MP", variant, ph_token)
    else:
        existing_ids = {
            str(k).split("|", 1)[0].upper()
            for k in proteins.keys()
            if isinstance(k, str)
        }
        if "TEMP" not in existing_ids:
            variant = "HOLO"
            proteins[f"TEMP|{variant}|{ph_token}"] = _ensure_protein_entry("TEMP", variant, ph_token)
            proteins[f"TEMP|{variant}|{ph_token}"]["status"] = "completed"
        if "T3MP" not in existing_ids:
            variant = "HOLO"
            proteins[f"T3MP|{variant}|{ph_token}"] = _ensure_protein_entry("T3MP", variant, ph_token)

    manifest["proteins"] = proteins
    cmd_block = manifest.get("command") or {}
    cmd_block["APO_HOLO_MODE"] = "holo"
    cmd_block["PH_ENSEMBLE_ENABLE"] = True
    cmd_block["TEST_MODE_ENABLE"] = "dud+fda"
    cmd_block["config_hash"] = compute_config_hash(cfg)
    manifest["command"] = cmd_block
    protein_ids: list[str] = []
    for key in proteins:
        try:
            pid = str(key).split("|", 1)[0].upper()
        except Exception:
            pid = str(key)
        if pid not in protein_ids:
            protein_ids.append(pid)
    completed_count = sum(
        1
        for entry in proteins.values()
        if isinstance(entry, dict) and str(entry.get("status", "")).lower() == "completed"
    )
    manifest["summary"] = {
        "total_proteins_scheduled": len(protein_ids),
        "total_proteins_completed": completed_count,
        "total_proteins_failed": 0,
        "total_protein_list": protein_ids,
    }
    paths_block = manifest.get("paths") or {}
    paths_block["run_dir"] = str(cfg["CONFIG_RUN_DIR"])
    paths_block["config_file"] = cfg_file
    manifest["paths"] = paths_block
    manifest["status"] = "running"

    backup_path = manifest_dir / "run_manifest.yaml.backup_resume_test"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    # Write interrupted state and back it up so we can restore it after resume.
    manifest_path.write_text(yaml.safe_dump(manifest))
    shutil.copyfile(manifest_path, backup_path)

    # Simulate config drift by altering config.txt (restore later).
    original_config = CONFIG_PATH.read_text()
    config_backup_path = CONFIG_PATH.with_suffix(".backup_resume_test")
    CONFIG_PATH.write_text(original_config)
    shutil.copyfile(CONFIG_PATH, config_backup_path)

    drift_cfg = original_config
    drift_cfg = drift_cfg.replace("APO_HOLO_MODE = holo", "APO_HOLO_MODE = apo")
    drift_cfg = drift_cfg.replace("PH_ENSEMBLE=true", "PH_ENSEMBLE=false")
    drift_cfg = drift_cfg.replace("TEST_MODE_ENABLE = dud+fda", "TEST_MODE_ENABLE = off")
    drift_cfg = drift_cfg.replace('"T3MP": "test_library_10",', "")
    CONFIG_PATH.write_text(drift_cfg)

    try:
        resume_cmd = [str(MAIN_PY), "-resume", "--run-id", RUN_ID]
        _run_python(resume_cmd, env=resume_env)

        manifest_after = _load_manifest(manifest_path)
        assert manifest_after["run_id"] == RUN_ID
        proteins_after = manifest_after.get("proteins") or {}

        def _lookup(pdb: str, variant: str) -> dict[str, Any] | None:
            matches: list[dict[str, Any]] = []
            for key, entry in proteins_after.items():
                parts = str(key).split("|", 2)
                if len(parts) != 3:
                    continue
                pdb_part, variant_part, _ph = parts
                if pdb_part.upper() == pdb and variant_part.lower() == variant.lower():
                    matches.append(entry)
            for entry in matches:
                if str(entry.get("status", "")).lower() == "completed":
                    return entry
            return matches[0] if matches else None

        entry_temp = _lookup("TEMP", "HOLO")
        entry_t3mp = _lookup("T3MP", "HOLO")
        assert entry_temp is not None and entry_temp.get("status") == "completed"
        assert entry_t3mp is not None
        t3mp_status = str(entry_t3mp.get("status", "")).lower()
        assert t3mp_status in {"completed", "pending", "running"}
        if t3mp_status != "completed":
            docked_root = Path((manifest_after.get("paths") or {}).get("docked_dir", REPO_ROOT / "docked"))
            planned = list((docked_root / "T3MP").rglob("planned_ligands_*.txt"))
            assert planned, "Resume run did not emit planned ligands for T3MP"
            assert all(p.stat().st_size > 0 for p in planned)

        summary = manifest_after.get("summary") or {}
        assert summary.get("total_proteins_scheduled") == 2
        completed_count = sum(
            1
            for entry in proteins_after.values()
            if isinstance(entry, dict) and entry.get("status") == "completed"
        )
        summary_completed = summary.get("total_proteins_completed")
        assert summary_completed is not None and int(summary_completed) >= completed_count >= 1
        assert summary.get("total_proteins_failed") == 0

        cmd = manifest_after.get("command") or {}
        assert cmd.get("APO_HOLO_MODE", "").lower() == "holo"
        assert bool(cmd.get("PH_ENSEMBLE_ENABLE")) is True
        assert cmd.get("TEST_MODE_ENABLE") == "dud+fda"
        assert cmd.get("config_hash") == compute_config_hash(cfg)

        snap_cfg = yaml.safe_load(cfg_snap_path.read_text())
        test_map = (snap_cfg or {}).get("TEST_LIBRARY_MAP") or {}
        if isinstance(test_map, str):
            test_map = yaml.safe_load(test_map) or {}
        test_map_norm = {str(k).upper(): v for k, v in test_map.items()}
        assert test_map_norm.get("TEMP") == "test_library_10"
        assert test_map_norm.get("T3MP") == "test_library_10"
    finally:
        # Restore global config drift and interrupted manifest state.
        if config_backup_path.exists():
            shutil.copyfile(config_backup_path, CONFIG_PATH)
            config_backup_path.unlink(missing_ok=True)
        else:
            CONFIG_PATH.write_text(original_config)
        if backup_path.exists():
            shutil.copyfile(backup_path, manifest_path)
