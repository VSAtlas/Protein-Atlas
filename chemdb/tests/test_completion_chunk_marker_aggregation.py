from __future__ import annotations

import logging
from pathlib import Path

from docking.completion_markers import (
    aggregate_stage_chunk_markers,
    completion_marker_path,
)
from docking.docking_utils import run_completion_audit


def _expected_output(stage_dir: Path, stage_name: str, lig: str) -> Path:
    return stage_dir / f"{Path(lig).stem}_{stage_name}.pdbqt"


def _write_pose(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("MODEL 1\nATOM\nENDMDL\n", encoding="utf-8")


def test_chunk_completion_markers_aggregate_to_canonical(tmp_path: Path) -> None:
    stage_dir = tmp_path / "stage1"
    stage_name = "stage1"
    ligands_chunk_1 = [
        str(tmp_path / "lig_a.pdbqt"),
        str(tmp_path / "lig_b.pdbqt"),
    ]
    ligands_chunk_2 = [
        str(tmp_path / "lig_c.pdbqt"),
        str(tmp_path / "lig_d.pdbqt"),
    ]

    for lig in ligands_chunk_1 + ligands_chunk_2:
        _write_pose(_expected_output(stage_dir, stage_name, lig))

    base_cfg = {"RUN_ID": "run_chunk_marker"}
    run_completion_audit(
        engine="vina",
        pdb_id="BNJS",
        stage_name=stage_name,
        ligands=ligands_chunk_1,
        expected_output_path=lambda lig: _expected_output(stage_dir, stage_name, lig),
        rerun_one=lambda lig: (True, "ok", None),
        stage_dir=stage_dir,
        cfg={**base_cfg, "_CHUNK_ID": "chunk_1"},
        logger=logging.getLogger("test"),
    )
    run_completion_audit(
        engine="vina",
        pdb_id="BNJS",
        stage_name=stage_name,
        ligands=ligands_chunk_2,
        expected_output_path=lambda lig: _expected_output(stage_dir, stage_name, lig),
        rerun_one=lambda lig: (True, "ok", None),
        stage_dir=stage_dir,
        cfg={**base_cfg, "_CHUNK_ID": "chunk_2"},
        logger=logging.getLogger("test"),
    )

    assert completion_marker_path(stage_dir, engine="vina", chunk_id="chunk_1").exists()
    assert completion_marker_path(stage_dir, engine="vina", chunk_id="chunk_2").exists()
    assert not completion_marker_path(stage_dir, engine="vina", chunk_id=None).exists()

    merged = aggregate_stage_chunk_markers(stage_dir, engine="vina")
    assert merged is not None
    assert int(merged.get("expected_count") or 0) == 4
    assert int(merged.get("missing_count_after") or 0) == 0
    assert bool(merged.get("success")) is True
    assert int(merged.get("chunk_marker_count") or 0) == 2
    assert completion_marker_path(stage_dir, engine="vina", chunk_id=None).exists()


def test_non_chunk_completion_still_writes_canonical(tmp_path: Path) -> None:
    stage_dir = tmp_path / "stage2"
    stage_name = "stage2"
    ligands = [str(tmp_path / "lig_e.pdbqt")]
    for lig in ligands:
        _write_pose(_expected_output(stage_dir, stage_name, lig))

    run_completion_audit(
        engine="vina",
        pdb_id="BNJS",
        stage_name=stage_name,
        ligands=ligands,
        expected_output_path=lambda lig: _expected_output(stage_dir, stage_name, lig),
        rerun_one=lambda lig: (True, "ok", None),
        stage_dir=stage_dir,
        cfg={"RUN_ID": "run_non_chunk"},
        logger=logging.getLogger("test"),
    )
    assert completion_marker_path(stage_dir, engine="vina", chunk_id=None).exists()
