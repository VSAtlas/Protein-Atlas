# ruff: noqa: E402
from __future__ import annotations

from pathlib import Path
import csv
import logging
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import post_docking.rescoring.rescoring_scorch as rescoring_scorch


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def test_pose_base_from_path_normalizes_dock6_mol2() -> None:
    old_style = Path("actives_final_00001.mol2__dock6_stage1.pdbqt")
    new_style = Path("actives_final_00001__dock6_stage1.pdbqt")

    assert rescoring_scorch._pose_base_from_path(old_style) == "actives_final_00001"
    assert rescoring_scorch._pose_base_from_path(new_style) == "actives_final_00001"


def test_post_engine_stage_preference_prefers_stage3(tmp_path: Path) -> None:
    ph_root = tmp_path / "combo"
    stage_root = ph_root / "dock6_pdbqt"
    stage_root.mkdir(parents=True, exist_ok=True)
    for stage in (1, 2, 3):
        (stage_root / f"lig_a__dock6_stage{stage}.pdbqt").write_text(
            "", encoding="utf-8"
        )

    logger = logging.getLogger("test")
    result = rescoring_scorch._collect_best_pose_per_base(
        ph_root, ("dock6_pdbqt",), None, logger
    )
    ligands = result.ligands
    stage_counts = result.stage_counts
    total = result.total_candidates
    available_bases = result.available_bases

    assert total == 3
    assert len(ligands) == 1
    assert ligands[0].name.endswith("__dock6_stage3.pdbqt")
    assert stage_counts.get(3) == 1
    assert available_bases == {"lig_a"}


def test_stage_preference_uses_selected_stage_when_available(tmp_path: Path) -> None:
    ph_root = tmp_path / "combo"
    stage_root = ph_root / "dock6_pdbqt"
    stage_root.mkdir(parents=True, exist_ok=True)
    stage1 = stage_root / "lig_a__dock6_stage1.pdbqt"
    stage3 = stage_root / "lig_a__dock6_stage3.pdbqt"
    stage1.write_text("", encoding="utf-8")
    stage3.write_text("", encoding="utf-8")

    result = rescoring_scorch._collect_best_pose_per_base(
        ph_root,
        ("dock6_pdbqt",),
        {"lig_a"},
        logging.getLogger("test.stage.preference"),
        preferred_stage_by_base={"lig_a": "dock6_stage1"},
    )

    assert [p.name for p in result.ligands] == ["lig_a__dock6_stage1.pdbqt"]
    assert result.rescored_stage_by_base["lig_a"] == "dock6_stage1"
    assert result.stage_fallback_reason_by_base == {}
