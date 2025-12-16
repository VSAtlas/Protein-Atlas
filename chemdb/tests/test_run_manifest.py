"""Tests for run_manifest summary helpers."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run_manifest import (
    _refresh_summary,
    get_manifest_paths,
    init_run_manifest,
    update_manifest_for_scheduled_proteins,
)

try:
    import yaml
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
    assert summary["total_proteins_scheduled"] == 2
    assert set(summary["total_protein_list"]) == {"TEST", "TEMP"}
    assert summary["total_proteins_completed"] == 1
    assert summary["total_proteins_failed"] == 2
