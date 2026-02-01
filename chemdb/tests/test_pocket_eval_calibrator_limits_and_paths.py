from __future__ import annotations

import logging
from pathlib import Path

import pocket_eval
import pytest


def test_select_calibrator_rows_respects_maximum(monkeypatch):
    cal_rows = []
    for i in range(50):
        cal_rows.append(
            {"ligand_id": f"chembl:strong_{i}", "smiles": "C", "label": "strong"}
        )
        cal_rows.append({"ligand_id": f"chembl:non_{i}", "smiles": "C", "label": "non"})
        cal_rows.append(
            {"ligand_id": f"chembl:weak_{i}", "smiles": "C", "label": "weak"}
        )

    rows_used, weak_ignored = pocket_eval._select_calibrator_rows_for_eval_and_prep(
        cal_rows, max_total=10, seed=0, logger=logging.getLogger("test")
    )

    assert len(rows_used) == 10
    assert all(r["label"] in {"strong", "non"} for r in rows_used)
    assert weak_ignored == 50


def test_prepare_calibrator_inputs_override_writes_to_extracted(tmp_path, monkeypatch):
    prepped_dir = tmp_path / "prepped_ligands" / "6LU7_calibrator"
    inputs_dir = tmp_path / "extracted_ligands" / "6LU7_calibrator"

    rows = [
        {"ligand_id": "chembl:strong_0", "smiles": "C", "label": "strong"},
        {"ligand_id": "chembl:non_0", "smiles": "CC", "label": "non"},
    ]

    def fake_convert_smi_to_sdf(smi_path, sdf_path, require_output=True):
        Path(smi_path).parent.mkdir(parents=True, exist_ok=True)
        Path(smi_path).write_text("SMI", encoding="utf-8")
        Path(sdf_path).write_text("SDF", encoding="utf-8")

    def fake_prep_ligands_with_mgltools(
        force=False, microstate_dedup=False, ph_values=None, root_dir=None
    ):
        for row in rows:
            sanitized = pocket_eval._sanitize_ligand_name_for_filename(row["ligand_id"])
            Path(root_dir).mkdir(parents=True, exist_ok=True)
            Path(root_dir, f"{sanitized}.pdbqt").write_text("PDBQT", encoding="utf-8")

    monkeypatch.setattr(pocket_eval, "convert_smi_to_sdf", fake_convert_smi_to_sdf)
    monkeypatch.setattr(
        pocket_eval, "prep_ligands_with_mgltools", fake_prep_ligands_with_mgltools
    )

    cfg = {"FORCE_REPROCESS": True}
    prepared, prep_report = pocket_eval.prepare_calibrator_ligands(
        rows, prepped_dir, cfg, logger=None, inputs_dir_override=inputs_dir
    )

    assert (inputs_dir / "calibrators.smi").is_file()
    assert (inputs_dir / "calibrators.sdf").is_file()
    assert not (prepped_dir / "inputs" / "calibrators.sdf").exists()
    assert len(prepared) == len(rows)
    assert prep_report["missing_ligand_ids"] == []


def test_filter_supported_calibrators_skips_unsupported_atoms(monkeypatch):
    if pocket_eval.Chem is None:
        pytest.skip("RDKit not available; filter is no-op without Chem")

    rows = [
        {"ligand_id": "ok1", "smiles": "CCO", "label": "strong"},
        {"ligand_id": "boron", "smiles": "B", "label": "non"},
    ]
    kept = pocket_eval._filter_supported_calibrators(rows, logging.getLogger("test"))
    kept_ids = {r["ligand_id"] for r in kept}
    assert "ok1" in kept_ids
    assert "boron" not in kept_ids
