from __future__ import annotations

import csv
import json
from pathlib import Path

from chemdb import pocket_eval_dataset


def _build_synthetic_inputs():
    pockets_plan = [
        {
            "pocket_id": "P1",
            "center": (0.0, 0.0, 0.0),
            "box_size": (10.0, 10.0, 10.0),
            "pocket": {"rank": 1, "score": 12.3, "residue_summary": "A:10"},
        },
        {
            "pocket_id": "P2",
            "center": (1.0, 1.0, 1.0),
            "box_size": (8.0, 8.0, 8.0),
            "pocket": {"rank": 2, "score": 9.1, "residue_summary": "B:20"},
        },
    ]

    calibrators_used = [
        {"ligand_id": "lig_strong_1", "label": "strong", "smiles": "C"},
        {"ligand_id": "lig_strong_2", "label": "strong", "smiles": "CC"},
        {"ligand_id": "lig_non_1", "label": "non", "smiles": "CCC"},
        {"ligand_id": "lig_non_2", "label": "non", "smiles": "CCCC"},
    ]

    scores_by_pocket = {
        "P1": [
            {"ligand_id": "lig_strong_1", "label": "strong", "score": -8.0},
            {"ligand_id": "lig_non_1", "label": "non", "score": -6.0},
        ],
        "P2": [{"ligand_id": "lig_strong_1", "label": "strong", "score": -7.5}],
    }

    prep_report = {
        "missing_ligand_ids": ["lig_non_2"],
        "invalid_ligands": [
            {"ligand_id": "lig_strong_2", "pdbqt_path": "bad.pdbqt", "bad_type": "ZN"}
        ],
        "prepared_ligand_ids": ["lig_strong_1", "lig_non_1"],
    }

    return pockets_plan, calibrators_used, scores_by_pocket, prep_report


def test_pocket_eval_dataset_export_writes_required_artifacts(tmp_path: Path) -> None:
    (
        pockets_plan,
        calibrators_used,
        scores_by_pocket,
        prep_report,
    ) = _build_synthetic_inputs()

    rows = pocket_eval_dataset.build_pocket_eval_dataset_long(
        run_id="RUN123",
        pdb_id="TEST",
        variant="HOLO",
        calibrators_used=calibrators_used,
        pockets_plan=pockets_plan,
        scores_by_pocket=scores_by_pocket,
        prep_report=prep_report,
    )
    rows, splits_payload = pocket_eval_dataset.make_pocket_eval_folds(
        rows,
        n_splits=2,
        seed=7,
        group_key="ligand_id",
        stratify_labels=False,
        strategy="hash_group",
    )

    dock_root = tmp_path / "pocket_eval"
    pocket_eval_dataset.write_pocket_eval_artifacts(
        dock_root,
        rows,
        config_snapshot={"run_id": "RUN123", "pdb_id": "TEST"},
        pockets_used=pocket_eval_dataset.build_pockets_used(pockets_plan),
        splits_payload=splits_payload,
        formats="csv",
    )

    dataset_csv = dock_root / "pocket_eval_dataset_long.csv"
    splits_path = dock_root / "pocket_eval_splits.json"
    pockets_path = dock_root / "pocket_eval_pockets_used.json"
    meta_path = dock_root / "pocket_eval_dataset_meta.json"

    assert dataset_csv.is_file()
    assert splits_path.is_file()
    assert pockets_path.is_file()
    assert meta_path.is_file()

    with dataset_csv.open("r", encoding="utf-8", newline="") as handle:
        rows_out = list(csv.DictReader(handle))
    assert len(rows_out) == 8

    statuses = {row["dock_status"] for row in rows_out}
    assert "prep_missing" in statuses
    assert "prep_invalid" in statuses
    assert "missing_score" in statuses

    splits_payload_disk = json.loads(splits_path.read_text(encoding="utf-8"))
    rows2 = pocket_eval_dataset.build_pocket_eval_dataset_long(
        run_id="RUN123",
        pdb_id="TEST",
        variant="HOLO",
        calibrators_used=calibrators_used,
        pockets_plan=pockets_plan,
        scores_by_pocket=scores_by_pocket,
        prep_report=prep_report,
    )
    _, splits_payload_2 = pocket_eval_dataset.make_pocket_eval_folds(
        rows2,
        n_splits=2,
        seed=7,
        group_key="ligand_id",
        stratify_labels=False,
        strategy="hash_group",
    )
    assert splits_payload_disk["group_to_fold"] == splits_payload_2["group_to_fold"]


def test_pocket_eval_dataset_cv_outputs(tmp_path: Path) -> None:
    (
        pockets_plan,
        calibrators_used,
        scores_by_pocket,
        prep_report,
    ) = _build_synthetic_inputs()
    rows = pocket_eval_dataset.build_pocket_eval_dataset_long(
        run_id="RUN123",
        pdb_id="TEST",
        variant="HOLO",
        calibrators_used=calibrators_used,
        pockets_plan=pockets_plan,
        scores_by_pocket=scores_by_pocket,
        prep_report=prep_report,
    )
    rows, splits_payload = pocket_eval_dataset.make_pocket_eval_folds(
        rows,
        n_splits=2,
        seed=3,
        group_key="ligand_id",
        stratify_labels=False,
        strategy="hash_group",
    )

    dock_root = tmp_path / "pocket_eval_cv"
    pocket_eval_dataset.write_pocket_eval_artifacts(
        dock_root,
        rows,
        config_snapshot={"run_id": "RUN123", "pdb_id": "TEST"},
        pockets_used=pocket_eval_dataset.build_pockets_used(pockets_plan),
        splits_payload=splits_payload,
        formats="csv",
        cv_enable=True,
        cv_top_m=1,
        cv_require_scores=False,
    )

    cv_path = dock_root / "pocket_eval_cv_assigned.json"
    cv_rows_path = dock_root / "pocket_eval_cv_assigned_per_ligand.csv"
    assert cv_path.is_file()
    payload = json.loads(cv_path.read_text(encoding="utf-8"))
    assert "cv_assigned_auc_mean" in payload
    assert "per_fold" in payload
    assert cv_rows_path.is_file()
