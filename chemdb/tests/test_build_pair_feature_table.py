from __future__ import annotations

import csv
from pathlib import Path

from analysis.build_pair_table import build_pair_feature_table


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_build_pair_feature_table_preserves_required_columns(tmp_path: Path) -> None:
    screening = tmp_path / "ranked_pairs.csv"
    annotations = tmp_path / "annotations.csv"
    _write_rows(
        screening,
        [
            {
                "ligand_base": "drug_a.sanitized",
                "target_id": "P12345",
                "pdb_id": "1ABC",
                "z_selected": "2.5",
                "free_cmax": "0.5",
                "Ki": "4",
            }
        ],
    )
    _write_rows(
        annotations,
        [
            {
                "drug_id": "drug_a",
                "target_id": "P12345",
                "target_adr_evidence": "1",
                "pathway_evidence": "1",
                "ligand_chemotype": "Kinase-like",
                "literature_supported_label": "1",
            }
        ],
    )

    rows = build_pair_feature_table(screening, annotations)

    assert rows[0]["drug_id"] == "drug_a"
    assert rows[0]["atlas_score"] == "2.5"
    assert rows[0]["exposure_plausibility"] == "1"
    assert rows[0]["target_adr_evidence"] == "1"
    assert rows[0]["pathway_evidence"] == "1"
    assert rows[0]["literature_supported_label"] == "1"


def test_build_pair_feature_table_converts_lower_is_better_score(tmp_path: Path) -> None:
    screening = tmp_path / "ranked_pairs.csv"
    _write_rows(
        screening,
        [{"drug_id": "drug_b", "target_id": "target_b", "pdb_id": "2DEF", "selected_docking_score": "-7.2"}],
    )

    rows = build_pair_feature_table(screening, None)

    assert rows[0]["atlas_score"] == "7.2"

