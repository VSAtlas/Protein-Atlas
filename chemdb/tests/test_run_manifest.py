"""Tests for run_manifest summary helpers."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import cast


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cli.run_manifest_runtime import (  # noqa: E402
    _refresh_summary,
    get_manifest_paths,
    init_run_manifest,
    update_manifest_for_docking_overall,
    update_manifest_for_protein_start,
    update_manifest_for_scheduled_proteins,
)

try:
    import yaml  # type: ignore[import-untyped]
except Exception:  # pragma: no cover
    yaml = None  # type: ignore[assignment]


def _load_manifest_dict(manifest_path: Path) -> dict:
    text = manifest_path.read_text(encoding="utf-8")
    if yaml is not None:
        return yaml.safe_load(text) or {}
    return json.loads(text)


def test_refresh_summary_dedupes_base_when_ph_present():
    manifest = {
        "proteins": {
            "TEST|APO|base": {"status": "completed"},
            "TEST|APO|pH7_0": {"status": "completed"},
        },
        "summary": {},
    }

    _refresh_summary(manifest)

    proteins = manifest["proteins"]
    assert "TEST|APO|base" not in proteins
    assert "TEST|APO|pH7_0" in proteins
    assert manifest["summary"]["total_proteins_scheduled"] == 1
    assert manifest["summary"]["total_proteins_completed"] == 1
    assert manifest["summary"]["total_proteins_failed"] == 0


def test_scheduled_protein_summary_preserved(tmp_path: Path) -> None:
    cfg = {
        "CONFIG_RUN_DIR": str(tmp_path / "runs"),
        "OVERALL_DIR": str(tmp_path / "overall"),
        "USE_GNINA": "false",
    }
    run_id = "test_scheduled_manifest"
    argv = ["main.py", "--run-id", run_id]
    log_path = str(tmp_path / "logs" / "dummy.log")

    init_run_manifest(cfg, run_id, argv, log_path)
    _, manifest_path = get_manifest_paths(cfg, run_id)

    update_manifest_for_scheduled_proteins(
        cfg,
        run_id,
        ["test", "TEST ", "-temp", "TEMP.pdb", "zz"],
    )

    manifest = _load_manifest_dict(manifest_path)
    summary = manifest.get("summary") or {}
    assert summary["total_proteins_scheduled"] == 2
    assert set(summary["total_protein_list"]) == {"TEST", "TEMP"}

    manifest["proteins"] = {
        "TEST|LEGACY|base": {"status": "completed", "pdb_id": "TEST"},
        "TEMP|LEGACY|base": {"status": "failed", "pdb_id": "TEMP"},
        "T3MP|LEGACY|base": {"status": "failed", "pdb_id": "T3MP"},
    }

    _refresh_summary(manifest)
    summary = manifest.get("summary") or {}
    assert summary["total_proteins_scheduled"] == 3
    assert set(summary["total_protein_list"]) == {"TEST", "TEMP"}
    assert summary["total_proteins_completed"] == 1
    assert summary["total_proteins_failed"] == 2


def test_distributed_ph_ensemble_suppresses_base_key(monkeypatch, tmp_path: Path) -> None:
    cfg = {
        "CONFIG_RUN_DIR": str(tmp_path / "runs"),
        "OVERALL_DIR": str(tmp_path / "overall"),
        "USE_GNINA": "false",
        "PH_ENSEMBLE": True,
    }
    run_id = "test_dist_ph_base_suppressed"
    argv = ["main.py", "--run-id", run_id]
    log_path = str(tmp_path / "logs" / "dummy.log")

    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "2")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "1")
    monkeypatch.setenv("SLURM_ARRAY_TASK_MIN", "0")

    init_run_manifest(cfg, run_id, argv, log_path)
    update_manifest_for_protein_start(
        cfg,
        run_id,
        "TEST",
        "HOLO",
        "fda_test_library_10",
        ph_tag=None,
    )

    proteins_dir = (
        Path(cast(str, cfg["OVERALL_DIR"]))
        / "manifests"
        / run_id
        / "distributed"
        / "proteins"
    )
    files = sorted(proteins_dir.glob("*.json")) if proteins_dir.exists() else []
    assert files == []

    update_manifest_for_protein_start(
        cfg,
        run_id,
        "TEST",
        "HOLO",
        "fda_test_library_10",
        ph_tag="pH7_0",
    )
    ph_files = sorted(proteins_dir.glob("*.json")) if proteins_dir.exists() else []
    assert len(ph_files) == 1
    payload = json.loads(ph_files[0].read_text(encoding="utf-8"))
    assert payload.get("protein_key") == "TEST|HOLO|pH7_0"


def test_distributed_non_ph_ensemble_keeps_base_key(monkeypatch, tmp_path: Path) -> None:
    cfg = {
        "CONFIG_RUN_DIR": str(tmp_path / "runs"),
        "OVERALL_DIR": str(tmp_path / "overall"),
        "USE_GNINA": "false",
        "PH_ENSEMBLE": False,
    }
    run_id = "test_dist_non_ph_base_kept"
    argv = ["main.py", "--run-id", run_id]
    log_path = str(tmp_path / "logs" / "dummy.log")

    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "2")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "1")
    monkeypatch.setenv("SLURM_ARRAY_TASK_MIN", "0")

    init_run_manifest(cfg, run_id, argv, log_path)
    update_manifest_for_protein_start(
        cfg,
        run_id,
        "TEST",
        "HOLO",
        "fda_test_library_10",
        ph_tag=None,
    )

    proteins_dir = (
        Path(cast(str, cfg["OVERALL_DIR"]))
        / "manifests"
        / run_id
        / "distributed"
        / "proteins"
    )
    files = sorted(proteins_dir.glob("*.json")) if proteins_dir.exists() else []
    assert len(files) == 1
    payload = json.loads(files[0].read_text(encoding="utf-8"))
    assert payload.get("protein_key") == "TEST|HOLO|base"


def test_distributed_ph_ensemble_suppresses_base_for_docking_overall(
    monkeypatch, tmp_path: Path
) -> None:
    cfg = {
        "CONFIG_RUN_DIR": str(tmp_path / "runs"),
        "OVERALL_DIR": str(tmp_path / "overall"),
        "USE_GNINA": "false",
        "PH_ENSEMBLE": True,
    }
    run_id = "test_dist_ph_base_suppressed_docking_overall"
    argv = ["main.py", "--run-id", run_id]
    log_path = str(tmp_path / "logs" / "dummy.log")

    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "2")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "1")
    monkeypatch.setenv("SLURM_ARRAY_TASK_MIN", "0")

    init_run_manifest(cfg, run_id, argv, log_path)
    update_manifest_for_docking_overall(
        cfg,
        run_id,
        "TEST",
        "HOLO",
        ph_tag=None,
        event="start",
    )

    proteins_dir = (
        Path(cast(str, cfg["OVERALL_DIR"]))
        / "manifests"
        / run_id
        / "distributed"
        / "proteins"
    )
    files = sorted(proteins_dir.glob("*.json")) if proteins_dir.exists() else []
    assert files == []

    update_manifest_for_docking_overall(
        cfg,
        run_id,
        "TEST",
        "HOLO",
        ph_tag="pH7_0",
        event="start",
    )
    ph_files = sorted(proteins_dir.glob("*.json")) if proteins_dir.exists() else []
    assert len(ph_files) == 1
    payload = json.loads(ph_files[0].read_text(encoding="utf-8"))
    assert payload.get("protein_key") == "TEST|HOLO|pH7_0"
