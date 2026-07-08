from __future__ import annotations

import csv
import logging
import math
from pathlib import Path

from post_docking.rescoring.rescoring_scorch import (
    SCORCH_TOP_FRACTION_DEFAULT,
    _chunk_allowed_bases,
    _score_csv_for_spec,
    _select_top_bases_from_score_csv,
    _tail_split_allowed_chunks,
)
from post_docking.rescoring.scorch_selection_invariants import (
    summarize_chunk_partition,
)
from post_docking.rescoring.scorch_types import StageSpec


def _write_score_csv(path: Path, rows: list[tuple[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ligand", "stage3"])
        writer.writeheader()
        for ligand, score in rows:
            writer.writerow({"ligand": ligand, "stage3": f"{score:.6f}"})


def _write_consensus_csv(
    path: Path,
    ligands: list[str],
    scores: list[float] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ligand", "run_mode", "library", "is_decoy", "consensus_score"],
        )
        writer.writeheader()
        for idx, ligand in enumerate(ligands):
            score = scores[idx] if scores is not None else 1.0 - (idx * 0.001)
            writer.writerow(
                {
                    "ligand": ligand,
                    "run_mode": "dud",
                    "library": "dud",
                    "is_decoy": "1",
                    "consensus_score": f"{score:.6f}",
                }
            )


def test_scorch_top_fraction_default_is_ten_percent() -> None:
    assert math.isclose(float(SCORCH_TOP_FRACTION_DEFAULT), 0.10, rel_tol=0.0, abs_tol=1e-12)


def test_top_fraction_selection_is_global_not_per_chunk(tmp_path: Path) -> None:
    score_csv = tmp_path / "docking_score_summary.csv"
    rows = [
        (f"lig_{idx:03d}_stage3.pdbqt", float(idx))
        for idx in range(100)
    ]
    _write_score_csv(score_csv, rows)

    selection = _select_top_bases_from_score_csv(
        score_csv,
        frac=0.10,
        control_bases=set(),
        higher_is_better=False,
        logger=logging.getLogger("test-scorch-global"),
        score_cols=["stage3"],
    )
    allowed = selection.allowed_bases
    n_pool = selection.n_pool
    k = selection.k
    controls = selection.controls_total

    assert n_pool == 100
    assert k == 10
    assert controls == 0
    assert len(allowed) == 10
    assert allowed == {f"lig_{idx:03d}" for idx in range(10)}


def test_chunking_preserves_selected_base_set(tmp_path: Path) -> None:
    score_csv = tmp_path / "docking_score_summary.csv"
    rows = [
        (f"lig_{idx:03d}_stage3.pdbqt", float(idx))
        for idx in range(37)
    ]
    _write_score_csv(score_csv, rows)
    control_bases = {"ctrl_a", "ctrl_b"}

    selection = _select_top_bases_from_score_csv(
        score_csv,
        frac=0.10,
        control_bases=control_bases,
        higher_is_better=False,
        logger=logging.getLogger("test-scorch-partition"),
        score_cols=["stage3"],
    )
    allowed = selection.allowed_bases
    k = selection.k
    controls = selection.controls_total
    assert controls == len(control_bases)
    assert k == 4
    assert len(allowed) == k + len(control_bases)

    chunks = _chunk_allowed_bases(allowed, chunk_size=3)
    chunks = _tail_split_allowed_chunks(chunks, tail_mode=True, min_split=2)
    stats = summarize_chunk_partition(allowed, chunks)
    assert stats.preserved is True
    assert stats.duplicate_count == 0
    assert stats.chunk_unique == len(allowed)


def test_selection_tracks_stage_provenance_for_best_score_columns(tmp_path: Path) -> None:
    score_csv = tmp_path / "docking_score_summary.csv"
    score_csv.parent.mkdir(parents=True, exist_ok=True)
    with score_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ligand", "stage3", "stage2", "stage1"])
        writer.writeheader()
        writer.writerow(
            {
                "ligand": "lig_a_stage3.pdbqt",
                "stage3": "-6.0",
                "stage2": "-9.0",
                "stage1": "-8.0",
            }
        )
        writer.writerow(
            {
                "ligand": "lig_b_stage3.pdbqt",
                "stage3": "-7.0",
                "stage2": "-6.0",
                "stage1": "-5.0",
            }
        )

    selection = _select_top_bases_from_score_csv(
        score_csv,
        frac=1.0,
        control_bases=set(),
        higher_is_better=False,
        logger=logging.getLogger("test-scorch-stage-provenance"),
        score_cols=["stage3", "stage2", "stage1"],
    )

    assert selection.selected_stage_by_base["lig_a"] == "stage2"
    assert math.isclose(selection.selected_score_by_base["lig_a"], -9.0)
    assert selection.selected_stage_by_base["lig_b"] == "stage3"


def test_dud_scorch_selection_uses_consensus_score_not_stream_order(
    tmp_path: Path,
) -> None:
    combo = ("P1", "", "")
    combo_root = tmp_path / "P1"
    ligands = [f"decoys_fda_{idx:05d}.pdbqt" for idx in range(50)]
    scores = [0.001 * idx for idx in range(50)]
    for rank, idx in enumerate((17, 33, 41, 4, 29), start=1):
        scores[idx] = 10.0 - rank
    _write_consensus_csv(
        combo_root / "dud_consensus_docking_scores.csv",
        ligands,
        scores,
    )
    _write_score_csv(
        combo_root / "dud_docking_score_summary.csv",
        [(f"decoys_fda_{idx:05d}_dud_stage1.pdbqt", -10.0 - idx) for idx in range(2)],
    )

    score_csv, score_cols, higher_is_better = _score_csv_for_spec(
        tmp_path,
        combo,
        StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv"),
        "dud",
    )

    assert score_csv == combo_root / "dud_consensus_docking_scores.csv"
    selection = _select_top_bases_from_score_csv(
        score_csv,
        frac=0.10,
        control_bases=set(),
        higher_is_better=higher_is_better,
        logger=logging.getLogger("test-scorch-dud-consensus-source"),
        score_cols=score_cols,
    )

    assert selection.n_pool == 50
    assert selection.k == 5
    assert selection.allowed_bases == {
        "decoys_fda_00017",
        "decoys_fda_00033",
        "decoys_fda_00041",
        "decoys_fda_00004",
        "decoys_fda_00029",
    }
