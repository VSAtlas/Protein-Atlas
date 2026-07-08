from __future__ import annotations

import csv
import logging
from pathlib import Path

import post_docking.rescoring.rescoring_scorch as rescoring_scorch


def _write_consensus(path: Path, rows: list[tuple[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ligand", "consensus_score"])
        for lig, score in rows:
            writer.writerow([lig, score])


def _touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("POSE\n", encoding="utf-8")


def test_scorch_selection_and_stage_counts_are_exact(tmp_path: Path) -> None:
    logger = logging.getLogger("test_scorch_counts")
    consensus_csv = tmp_path / "consensus_docking_scores.csv"
    rows = [(f"lig_{i:02d}_stage1.pdbqt", float(i)) for i in range(1, 21)]
    _write_consensus(consensus_csv, rows)

    allowed_bases, _, total_rows, selected_rows, *_ = rescoring_scorch._load_consensus_top_bases(
        consensus_csv, 0.4, logger
    )
    assert total_rows == 20
    assert selected_rows == 8
    assert len(allowed_bases) == 8

    ph_root = tmp_path / "combo"
    allowed = sorted(allowed_bases)
    for idx, base in enumerate(allowed):
        if idx < 5:
            _touch(ph_root / "stage3" / f"{base}_stage3.pdbqt")
        elif idx < 7:
            _touch(ph_root / "stage2" / f"{base}_stage2.pdbqt")
        else:
            _touch(ph_root / "stage1" / f"{base}_stage1.pdbqt")
    _touch(ph_root / "stage3" / "extra_not_selected_stage3.pdbqt")

    result = rescoring_scorch._collect_best_pose_per_base(
        ph_root, ("stage3", "stage2", "stage1"), allowed_bases, logger
    )
    ligands = result.ligands
    stage_counts = result.stage_counts
    total_candidates = result.total_candidates
    available_bases = result.available_bases
    assert len(ligands) == 8
    assert total_candidates == 9
    assert stage_counts.get(3, 0) == 5
    assert stage_counts.get(2, 0) == 2
    assert stage_counts.get(1, 0) == 1
    assert allowed_bases.issubset(available_bases)
    assert "extra_not_selected" in available_bases


def test_scorch_excludes_unselected_extras_from_counts(tmp_path: Path) -> None:
    logger = logging.getLogger("test_scorch_counts_extras")
    allowed = {"lig_a", "lig_b", "lig_c"}
    ph_root = tmp_path / "combo"
    _touch(ph_root / "stage3" / "lig_a_stage3.pdbqt")
    _touch(ph_root / "stage2" / "lig_b_stage2.pdbqt")
    _touch(ph_root / "stage1" / "lig_c_stage1.pdbqt")
    _touch(ph_root / "stage3" / "other_01_stage3.pdbqt")
    _touch(ph_root / "stage2" / "other_02_stage2.pdbqt")

    result = rescoring_scorch._collect_best_pose_per_base(
        ph_root, ("stage3", "stage2", "stage1"), allowed, logger
    )
    ligands = result.ligands
    stage_counts = result.stage_counts
    total_candidates = result.total_candidates
    available_bases = result.available_bases
    assert len(ligands) == 3
    assert total_candidates == 5
    assert stage_counts.get(3, 0) == 1
    assert stage_counts.get(2, 0) == 1
    assert stage_counts.get(1, 0) == 1
    assert allowed.issubset(available_bases)
    assert {"other_01", "other_02"}.issubset(available_bases)


def test_stage_preference_falls_back_when_selected_pose_missing(tmp_path: Path) -> None:
    logger = logging.getLogger("test_scorch_counts_fallback")
    ph_root = tmp_path / "combo"
    _touch(ph_root / "stage3" / "lig_x_stage3.pdbqt")

    result = rescoring_scorch._collect_best_pose_per_base(
        ph_root,
        ("stage3", "stage2", "stage1"),
        {"lig_x"},
        logger,
        preferred_stage_by_base={"lig_x": "stage1"},
    )

    assert [p.name for p in result.ligands] == ["lig_x_stage3.pdbqt"]
    assert result.rescored_stage_by_base["lig_x"] == "stage3"
    assert result.stage_fallback_reason_by_base["lig_x"] == "preferred_missing:stage1"
