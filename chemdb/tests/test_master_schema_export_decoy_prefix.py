from __future__ import annotations

import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.reporting import master_schema_export as master_schema_export  # noqa: E402


def _write_decoy_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "consensus_score",
        "consensus_score_pre",
        "ligand_file",
        "library",
        "is_decoy",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def test_decoy_prefix_consensus_used(tmp_path: Path) -> None:
    decoy_prefix = "fda_dud"
    run_id = "unitrun"
    pdb_id = "1ABC"
    variant = "HOLO"
    ph = "pH7_0"
    combo_dir = tmp_path / "post_docked" / run_id / pdb_id / variant / ph
    combo_dir.mkdir(parents=True, exist_ok=True)

    rows = [
        {
            "consensus_score": "1.0",
            "consensus_score_pre": "0.9",
            "ligand_file": "fda_dud_lig1.pdbqt",
            "library": "DECOY",
            "is_decoy": "1",
        },
        {
            "consensus_score": "0.8",
            "consensus_score_pre": "0.7",
            "ligand_file": "fda_dud_lig2.pdbqt",
            "library": "DECOY",
            "is_decoy": "1",
        },
    ]
    _write_decoy_csv(
        combo_dir / f"{decoy_prefix}_consensus_reranked_scorch.csv", rows
    )

    field, scores, fallback_used, _scores_by_field, _source = (
        master_schema_export._get_decoy_scores_for_combo(
            combo_dir, rows, min_decoys=1, decoy_prefix=decoy_prefix
        )
    )

    assert field is not None
    assert fallback_used is False
    assert len(scores) == len(rows)


def test_legacy_dud_consensus_used(tmp_path: Path) -> None:
    decoy_prefix = "fda_dud"
    run_id = "unitrun"
    pdb_id = "2XYZ"
    variant = "APO"
    ph = "pH5_0"
    combo_dir = tmp_path / "post_docked" / run_id / pdb_id / variant / ph
    combo_dir.mkdir(parents=True, exist_ok=True)

    rows = [
        {
            "consensus_score": "2.0",
            "consensus_score_pre": "1.9",
            "ligand_file": "dud_lig1.pdbqt",
            "library": "DECOY",
            "is_decoy": "1",
        }
    ]
    _write_decoy_csv(combo_dir / "dud_consensus_reranked_scorch.csv", rows)

    field, scores, fallback_used, _scores_by_field, _source = (
        master_schema_export._get_decoy_scores_for_combo(
            combo_dir, rows, min_decoys=1, decoy_prefix=decoy_prefix
        )
    )

    assert field is not None
    assert fallback_used is False
    assert len(scores) == len(rows)
