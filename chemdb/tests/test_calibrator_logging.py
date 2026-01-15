import json
from pathlib import Path

from calibrator import chemdbl_calibrator


def test_calibrator_logging_and_audit(monkeypatch):
    def fake_resolve(pdb_id, deepcoy_root=chemdbl_calibrator.DEFAULT_DEEPCOY_ROOT):
        return "P0DTD1", [], {"source": "rcsb_polymer_entity"}

    def fake_fetch(
        uniprot_id,
        pdb_id,
        cache_dir,
        timeout,
        retries,
        activity_types,
        chembl_max_phase,
        label_thresholds,
        session=None,
        debug=False,
        debug_max_ids=25,
        debug_rejection_samples_max=25,
        unp_start=None,
        unp_end=None,
        target_keywords=None,
    ):
        assert debug is True
        return {"strong": [], "weak": [], "non": []}, {
            "cached": False,
            "counts": {"targets": 1, "assays": 1, "activities": 3, "molecules": 2},
            "debug": {
                "target_ids_sample": ["T1"],
                "assay_ids_sample": ["A1"],
                "skips": {"skip_relation_not_equal": 2},
                "activity_types_seen_top": [{"type": "Ki", "count": 3}],
            },
        }

    monkeypatch.setattr(chemdbl_calibrator, "resolve_target", fake_resolve)
    monkeypatch.setattr(chemdbl_calibrator, "fetch_chembl_labeled_smiles", fake_fetch)

    run_tag = "pytest_run"
    log_file = chemdbl_calibrator.DEFAULT_LOG_ROOT / run_tag / "6LU7.log"
    previous_mtime = log_file.stat().st_mtime if log_file.exists() else 0

    chemdbl_calibrator.main(["--pdbs", "6LU7", "--run-tag", run_tag, "--debug-chembl"])

    assert log_file.is_file()
    assert log_file.stat().st_mtime >= previous_mtime
    log_text = log_file.read_text()
    assert "[calibrator.source_audit]" in log_text

    audit_path = Path("extracted_ligands/6LU7_calibrator/calibrator_audit.json")
    assert audit_path.is_file()
    audit = json.loads(audit_path.read_text())
    assert audit["run_tag"] == run_tag
    assert audit["chembl_meta"]["debug"]["target_ids_sample"] == ["T1"]
