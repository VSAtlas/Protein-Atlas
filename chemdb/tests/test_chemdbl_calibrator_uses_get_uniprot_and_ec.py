import json

from calibrator import chemdbl_calibrator


def test_chemdbl_calibrator_uses_get_uniprot_and_ec(tmp_path, monkeypatch):
    calls = {}

    def fake_get_uniprot(pdb_id):
        calls["pdb_id"] = pdb_id
        return "P12345", ["1.1.1.1"]

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
        cache_dir.mkdir(parents=True, exist_ok=True)
        return {"strong": ["CCO"], "weak": [], "non": []}, {
            "cached": False,
            "counts": {"targets": 1, "assays": 1, "activities": 1, "molecules": 1},
        }

    monkeypatch.setattr(
        "calibrator.chemdbl_calibrator.get_uniprot_and_ec", fake_get_uniprot
    )
    monkeypatch.setattr(
        "calibrator.chemdbl_calibrator.fetch_chembl_labeled_smiles", fake_fetch
    )

    out_root = tmp_path / "outputs"
    meta = chemdbl_calibrator.run_calibrator_for_pdb(
        "test",
        out_root=out_root,
        cache_dir=tmp_path / "cache",
        activity_types=["Ki"],
        label_thresholds=chemdbl_calibrator.DEFAULT_LABEL_THRESHOLDS,
        deepcoy_root=tmp_path / "deepcoy",
        log_dir=tmp_path / "logs",
        run_tag="unit_test",
    )

    assert calls["pdb_id"] == "TEST"
    output_dir = out_root / "TEST_calibrator"
    assert (output_dir / "strong_binders.smi").is_file()
    stored = json.loads((output_dir / "calibrator_meta.json").read_text())
    assert stored["mapping"]["source"] == "rcsb_polymer_entity"
    assert stored["uniprot_id"] == "P12345"
    assert meta["counts"]["strong"] == 1
