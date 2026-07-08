# ruff: noqa: E402
from __future__ import annotations

import csv
import logging
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import post_docking.rescoring.rescoring_scorch as rescoring_scorch


def _write_consensus(path: Path, rows: list[tuple[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ligand", "consensus_score"])
        for lig, score in rows:
            writer.writerow([lig, score])


def test_global_selection_shared_across_engines(tmp_path: Path) -> None:
    logger = logging.getLogger("test")
    control_bases = {"ctrl_A", "ctrl_B"}

    consensus_rows: list[tuple[str, float]] = []
    for idx in range(1, 21):
        consensus_rows.append((f"lig_{idx:02d}_stage1.pdbqt", float(idx)))
    consensus_rows.append(("ctrl_A_stage1.pdbqt", 0.5))
    consensus_rows.append(("ctrl_B_stage1.pdbqt", 0.4))

    consensus_csv = tmp_path / "consensus_docking_scores.csv"
    _write_consensus(consensus_csv, consensus_rows)

    (
        allowed_bases,
        allowed_decoy_redundant,
        total_rows,
        selected_rows,
        noncontrol_rows,
        controls_total,
        controls_in_consensus,
    ) = rescoring_scorch._load_consensus_top_bases(
        consensus_csv, 0.10, logger, control_bases
    )

    assert total_rows == 22
    assert noncontrol_rows == 20
    assert controls_total == len(control_bases)
    expected_allowed = (2) + len(control_bases)  # ceil(0.1 * 20) + controls
    assert len(allowed_bases) == expected_allowed
    assert controls_in_consensus == len(control_bases)

    ph_root = tmp_path / "combo"
    ph_root.mkdir(parents=True, exist_ok=True)

    def _touch(folder: Path, name: str) -> Path:
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / name
        path.write_text("", encoding="utf-8")
        return path

    target_base = next(iter(allowed_bases - control_bases))

    for base in allowed_bases:
        _touch(ph_root / "stage3", f"{base}_stage3.pdbqt")
        _touch(ph_root / "gnina_stage3", f"{base}_gnina_stage3.pdbqt")
        _touch(ph_root / "dock6_pdbqt", f"{base}__dock6_stage1.pdbqt")
        _touch(ph_root / "ledock_pdbqt", f"{base}__ledock_stage1.pdbqt")

    _touch(ph_root / "stage1", f"{target_base}_stage1.pdbqt")
    _touch(ph_root / "gnina_stage1", f"{target_base}_gnina_stage1.pdbqt")
    _touch(ph_root / "dock6_pdbqt", f"{target_base}__dock6_stage3.pdbqt")
    _touch(ph_root / "ledock_pdbqt", f"{target_base}__ledock_stage3.pdbqt")

    _touch(ph_root / "stage3", "extra_stage3.pdbqt")
    _touch(ph_root / "dock6_pdbqt", "extra__dock6_stage1.pdbqt")

    def _collect(stage_dirs: tuple[str, ...]) -> tuple[dict[str, str], set[str]]:
        result = rescoring_scorch._collect_best_pose_per_base(
            ph_root, stage_dirs, allowed_bases, logger
        )
        ligs = result.ligands
        available = result.available_bases
        return {
            rescoring_scorch._pose_base_from_path(p): p.name for p in ligs
        }, available

    vina_map, vina_available = _collect(("stage3", "stage2", "stage1"))
    gnina_map, gnina_available = _collect(
        ("gnina_stage3", "gnina_stage2", "gnina_stage1")
    )
    dock6_map, dock6_available = _collect(("dock6_pdbqt",))
    ledock_map, ledock_available = _collect(("ledock_pdbqt",))

    for base_map in (vina_map, gnina_map, dock6_map, ledock_map):
        assert set(base_map.keys()) == allowed_bases

    assert vina_map[target_base].endswith("stage3.pdbqt")
    assert gnina_map[target_base].endswith("stage3.pdbqt")
    assert "dock6_stage3" in dock6_map[target_base]
    assert "ledock_stage3" in ledock_map[target_base]

    for available in (
        vina_available,
        gnina_available,
        dock6_available,
        ledock_available,
    ):
        assert allowed_bases.issubset(available)
