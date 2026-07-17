from __future__ import annotations

import csv
from pathlib import Path

import pytest

from analysis.ml.consensus_z_repair import _read_reference_run


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_reference_layout(root: Path, relative: Path) -> Path:
    score_dir = root / relative
    _write_csv(
        score_dir / "dud_consensus_docking_scores.csv",
        [
            {"is_decoy": "1", "consensus_score": "0.0"},
            {"is_decoy": "1", "consensus_score": "2.0"},
            {"is_decoy": "0", "consensus_score": "99.0"},
        ],
    )
    _write_csv(
        score_dir / "consensus_docking_scores.csv",
        [{"ligand": "drug.pdbqt", "consensus_score": "1.5"}],
    )
    return score_dir


def test_reference_repair_discovers_flat_layout(tmp_path: Path) -> None:
    score_dir = _write_reference_layout(tmp_path, Path("1abc"))

    nulls, scores = _read_reference_run(tmp_path)

    assert nulls["1ABC"]["reference_layout"] == "flat"
    assert nulls["1ABC"]["reference_score_dir"] == str(score_dir)
    assert nulls["1ABC"]["mu"] == pytest.approx(1.0)
    assert nulls["1ABC"]["sigma"] == pytest.approx(1.0)
    assert scores[("1ABC", "drug")] == pytest.approx(1.5)


def test_reference_repair_discovers_nested_variant_ph_layout(
    tmp_path: Path,
) -> None:
    score_dir = _write_reference_layout(
        tmp_path,
        Path("1abc") / "HOLO" / "pH7_0",
    )

    nulls, scores = _read_reference_run(tmp_path)

    assert nulls["1ABC"]["reference_layout"] == "HOLO/pH7_0"
    assert nulls["1ABC"]["reference_score_dir"] == str(score_dir)
    assert scores[("1ABC", "drug")] == pytest.approx(1.5)


def test_reference_repair_rejects_ambiguous_strata(tmp_path: Path) -> None:
    _write_reference_layout(tmp_path, Path("1abc") / "HOLO" / "pH7_0")
    _write_reference_layout(tmp_path, Path("1abc") / "APO" / "pH7_0")

    with pytest.raises(
        ValueError,
        match="ambiguous DUD consensus layouts.*1ABC",
    ):
        _read_reference_run(tmp_path)
