from __future__ import annotations

from pathlib import Path
import csv
import logging
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import rescoring_scorch


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
        (stage_root / f"lig_a__dock6_stage{stage}.pdbqt").write_text("", encoding="utf-8")

    logger = logging.getLogger("test")
    ligands, stage_counts, total = rescoring_scorch._collect_best_pose_per_base(
        ph_root, ("dock6_pdbqt",), None, logger
    )

    assert total == 3
    assert len(ligands) == 1
    assert ligands[0].name.endswith("__dock6_stage3.pdbqt")
    assert stage_counts.get(3) == 1


def test_post_engine_top_fraction_selection(tmp_path: Path) -> None:
    logger = logging.getLogger("test")

    dock6_csv = tmp_path / "dock6_docking_score_long.csv"
    dock6_header = ["run_id", "stage", "ligand", "dock6_grid_score", "dock6_n_poses"]
    dock6_rows = []
    for idx in range(1, 11):
        ligand = f"lig_{idx:02d}.mol2"
        score = -float(idx)
        dock6_rows.append(["RUN1", "dock6_stage1", ligand, f"{score:.2f}", "10"])
    _write_csv(dock6_csv, dock6_header, dock6_rows)

    allowed, candidates, selected, controls, _ = rescoring_scorch._load_post_engine_top_bases(
        dock6_csv, 0.1, logger, "dock6_grid_score", set()
    )
    assert allowed is not None
    assert candidates == 10
    assert selected == 1
    assert controls == 0
    assert allowed == {"lig_10"}

    ledock_csv = tmp_path / "ledock_docking_score_long.csv"
    ledock_header = [
        "run_id",
        "stage",
        "ligand",
        "ledock_best_score_kcal",
        "ledock_cluster_count",
        "ledock_n_poses",
    ]
    ledock_rows = []
    for idx in range(1, 6):
        ligand = f"lig_{idx:02d}.mol2"
        score = -float(idx)
        ledock_rows.append(["RUN1", "ledock_stage1", ligand, f"{score:.2f}", "1", "10"])
    _write_csv(ledock_csv, ledock_header, ledock_rows)

    allowed_ledock, candidates_ledock, selected_ledock, controls_ledock, _ = rescoring_scorch._load_post_engine_top_bases(
        ledock_csv, 0.4, logger, "ledock_best_score_kcal", set()
    )
    assert allowed_ledock is not None
    assert candidates_ledock == 5
    assert selected_ledock == 2
    assert controls_ledock == 0
    assert len(allowed_ledock) == 2
