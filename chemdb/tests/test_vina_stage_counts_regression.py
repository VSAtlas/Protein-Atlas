from __future__ import annotations

import csv
from pathlib import Path

from docking.docking_vina import write_scores_csv


def _base_cfg(tmp_path: Path) -> dict[str, object]:
    cfg: dict[str, object] = {
        "OVERALL_DIR": str(tmp_path),
        "INPUT_DIR": str(tmp_path / "input_pdbs"),
        "OUTPUT_DIR": str(tmp_path / "processed_pdbs"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "PREPPED_LIGANDS_DIR": str(tmp_path / "prepped_ligands"),
        "RUN_ID": "stage_count_regress",
    }
    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR", "PREPPED_LIGANDS_DIR"):
        Path(str(cfg[key])).mkdir(parents=True, exist_ok=True)
    return cfg


def _score_record(score: float) -> dict[str, object]:
    return {
        "score": score,
        "valid": True,
        "reason": "",
        "heavy_atoms": 12,
        "le": None,
        "pains_flag": False,
    }


def _read_stage_counts(long_csv: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with long_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            stage = str(row.get("stage", "")).strip()
            if not stage:
                continue
            out[stage] = out.get(stage, 0) + 1
    return out


def test_vina_long_csv_preserves_exact_stage_counts(tmp_path: Path) -> None:
    cfg = _base_cfg(tmp_path)
    score_history = {
        "stage1": {f"lig_{i:02d}.pdbqt": _score_record(-7.0 - i) for i in range(10)},
        "stage2": {f"lig_{i:02d}.pdbqt": _score_record(-8.0 - i) for i in range(4)},
        "stage3": {f"lig_{i:02d}.pdbqt": _score_record(-9.0 - i) for i in range(2)},
    }
    summary_csv = Path(write_scores_csv(cfg, "TEST", score_history))
    long_csv = summary_csv.with_name("docking_score_long.csv")
    assert long_csv.exists()

    counts = _read_stage_counts(long_csv)
    assert counts == {"stage1": 10, "stage2": 4, "stage3": 2}


def test_vina_stage_counts_stable_after_chunk_merge(tmp_path: Path) -> None:
    cfg = _base_cfg(tmp_path)
    cfg["_CHUNK_LIGAND_KEYS"] = ["a_00"]
    first = {
        "stage1": {"a_00.pdbqt": _score_record(-7.0)},
        "stage2": {"a_00.pdbqt": _score_record(-8.0)},
    }
    cfg["_CHUNK_LIGAND_KEYS"] = ["b_00"]
    second = {
        "stage1": {"b_00.pdbqt": _score_record(-7.5)},
        "stage2": {"b_00.pdbqt": _score_record(-8.5)},
        "stage3": {"b_00.pdbqt": _score_record(-9.5)},
    }
    write_scores_csv(cfg, "TEST", first)
    summary_csv = Path(write_scores_csv(cfg, "TEST", second))
    long_csv = summary_csv.with_name("docking_score_long.csv")
    counts = _read_stage_counts(long_csv)
    assert counts == {"stage1": 2, "stage2": 2, "stage3": 1}
