import json
from pathlib import Path

from calibrator import chemdbl_calibrator


def test_calibrator_telemetry_printed_and_urls_full(monkeypatch):
    sample_url = (
        "https://www.ebi.ac.uk/chembl/api/data/target.json?"
        "target_components__accession=P0DTD1&limit=20&offset=0"
    )

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
        debug=False,
        debug_max_ids=25,
        debug_rejection_samples_max=25,
    ):
        telemetry = {
            "activity_sanity": {
                "n_activities_total": 2,
                "n_with_molecule_chembl_id": 2,
                "n_with_pchembl_value": 1,
                "relation_histogram": {"=": 2},
                "units_histogram": {"nM": 2},
                "n_standard_value_parseable": 2,
            },
            "labeling_sanity": {
                "n_labeled_strong": 0,
                "n_labeled_weak": 0,
                "n_labeled_non": 0,
                "n_rejected_by_relation": 0,
                "n_rejected_by_units": 0,
                "n_rejected_by_missing_pchembl": 1,
                "n_rejected_by_missing_standard_value": 1,
                "n_rejected_by_value_threshold": 0,
                "n_rejected_by_type_not_allowed": 0,
                "n_rejected_by_phase_filtered": 0,
            },
            "molecule_sanity": {
                "n_molecule_fetch_attempted": 1,
                "n_molecule_fetch_failed_http": 0,
                "n_molecule_missing_smiles": 0,
                "n_molecule_parsed_smiles_ok": 0,
            },
            "request_urls": [sample_url],
            "request_urls_truncated": False,
        }
        meta = {
            "cached": False,
            "counts": {"targets": 1, "assays": 1, "activities": 2, "molecules": 0},
            "telemetry": telemetry,
        }
        return {"strong": [], "weak": [], "non": []}, meta

    monkeypatch.setattr(chemdbl_calibrator, "resolve_target", fake_resolve)
    monkeypatch.setattr(chemdbl_calibrator, "fetch_chembl_labeled_smiles", fake_fetch)

    run_tag = "pytest_run_telemetry"
    log_file = chemdbl_calibrator.DEFAULT_LOG_ROOT / run_tag / "6LU7.log"
    chemdbl_calibrator.main(["--pdbs", "6LU7", "--run-tag", run_tag])

    assert log_file.is_file()
    log_text = log_file.read_text()
    assert "activity_sanity" in log_text
    assert "molecule_sanity" in log_text
    assert f"url={sample_url}" in log_text

    audit_path = Path("extracted_ligands/6LU7_calibrator/calibrator_audit.json")
    meta_path = Path("extracted_ligands/6LU7_calibrator/calibrator_meta.json")
    assert audit_path.is_file()
    assert meta_path.is_file()
    audit = json.loads(audit_path.read_text())
    stored_meta = json.loads(meta_path.read_text())
    assert audit["telemetry"]["request_urls"] == [sample_url]
    assert stored_meta.get("chembl_request_urls_sample") == [sample_url]
