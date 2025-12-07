"""Tests for docking manifest timing + per-stage details."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run_manifest import (
    get_manifest_paths,
    init_run_manifest,
    update_manifest_for_docking_overall,
    update_manifest_for_docking_stage,
    update_manifest_for_protein_start,
    _protein_key,
)

try:  # Optional dependency
    import yaml
except Exception:  # pragma: no cover - optional import guard
    yaml = None  # type: ignore[assignment]


def _load_manifest_dict(manifest_path: Path) -> dict:
    text = manifest_path.read_text(encoding="utf-8")
    if yaml is not None:
        return yaml.safe_load(text) or {}
    return json.loads(text)


def test_docking_stage_details_and_timing(tmp_path: Path) -> None:
    cfg = {
        "CONFIG_RUN_DIR": str(tmp_path / "runs"),
        "OVERALL_DIR": str(tmp_path / "overall"),
    }
    run_id = "test_docking_manifest"
    argv = ["main.py", "--run-id", run_id]
    log_path = str(tmp_path / "logs" / "dummy.log")

    init_run_manifest(cfg, run_id, argv, log_path)
    _, manifest_path = get_manifest_paths(cfg, run_id)
    assert manifest_path.exists()

    pdb_id = "1ABC"
    variant_label = "legacy"
    library = "fda"
    ph_tag = None

    update_manifest_for_protein_start(
        cfg,
        run_id,
        pdb_id,
        variant_label,
        library,
        ph_tag=ph_tag,
    )

    update_manifest_for_docking_overall(
        cfg,
        run_id,
        pdb_id,
        variant_label,
        ph_tag=ph_tag,
        event="start",
    )

    update_manifest_for_docking_stage(
        cfg,
        run_id,
        pdb_id,
        variant_label,
        "stage1",
        status="running",
        ph_tag=ph_tag,
    )
    update_manifest_for_docking_stage(
        cfg,
        run_id,
        pdb_id,
        variant_label,
        "stage1",
        status="completed",
        ph_tag=ph_tag,
        elapsed_sec=1.234,
    )
    update_manifest_for_docking_stage(
        cfg,
        run_id,
        pdb_id,
        variant_label,
        "dud_stage1",
        status="completed",
        ph_tag=ph_tag,
        elapsed_sec=2.0,
    )
    update_manifest_for_docking_stage(
        cfg,
        run_id,
        pdb_id,
        variant_label,
        "hmdb_stage2",
        status="completed",
        ph_tag=ph_tag,
        elapsed_sec=3.5,
    )

    update_manifest_for_docking_overall(
        cfg,
        run_id,
        pdb_id,
        variant_label,
        ph_tag=ph_tag,
        event="end",
        elapsed_sec=10.0,
    )

    manifest = _load_manifest_dict(manifest_path)
    proteins = manifest.get("proteins") or {}
    key = _protein_key(pdb_id, variant_label, ph_tag)
    assert key in proteins

    entry = proteins[key]
    stages = entry.get("stages") or {}
    docking = stages.get("docking") or {}

    timing = docking.get("timing") or {}
    assert timing.get("started_at")
    assert timing.get("finished_at")
    assert isinstance(timing.get("wall_time_sec"), (int, float))

    details = docking.get("details") or {}
    per_stage = details.get("per_stage") or {}
    assert "stage1" in per_stage
    assert "dud_stage1" in per_stage
    assert "hmdb_stage2" in per_stage

    stage_detail = per_stage["stage1"]
    assert stage_detail.get("status") == "completed"
    assert stage_detail.get("subrun") == "primary"
    assert stage_detail.get("stage_base_name") == "stage1"

    dud_stage_detail = per_stage["dud_stage1"]
    assert dud_stage_detail.get("subrun") == "dud"
    assert dud_stage_detail.get("stage_base_name") == "stage1"

    hmdb_stage_detail = per_stage["hmdb_stage2"]
    assert hmdb_stage_detail.get("subrun") == "hmdb"
    assert hmdb_stage_detail.get("stage_base_name") == "stage2"

    stage_timing = stage_detail.get("timing") or {}
    assert stage_timing.get("started_at")
    assert stage_timing.get("finished_at")
    assert isinstance(stage_timing.get("wall_time_sec"), (int, float))

    by_subrun = details.get("by_subrun") or {}
    assert "primary" in by_subrun
    assert "dud" in by_subrun
    assert "hmdb" in by_subrun

    assert "stage1" in by_subrun["primary"]
    assert "dud_stage1" in by_subrun["dud"]
    assert "hmdb_stage2" in by_subrun["hmdb"]

    assert by_subrun["primary"]["stage1"] is stage_detail
    assert by_subrun["dud"]["dud_stage1"] is dud_stage_detail
    assert by_subrun["hmdb"]["hmdb_stage2"] is hmdb_stage_detail
