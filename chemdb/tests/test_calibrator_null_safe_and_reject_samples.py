from DeepCoy_duds import external_sources as es


def test_fetch_chembl_labeled_smiles_null_safe_and_reject_samples(
    monkeypatch, tmp_path
):
    activities = [
        {
            "activity_chembl_id": "ACT1",
            "assay_chembl_id": "ASSAY1",
            "molecule_chembl_id": "M1",
            "standard_relation": ">",
            "standard_units": "nM",
            "standard_value": "50",
            "standard_type": "Ki",
            "pchembl_value": None,
        },
        {
            "activity_chembl_id": "ACT2",
            "assay_chembl_id": "ASSAY1",
            "molecule_chembl_id": "M2",
            "standard_relation": "=",
            "standard_units": "nM",
            "standard_value": "10",
            "standard_type": "Ki",
            "pchembl_value": None,
        },
        {
            "activity_chembl_id": "ACT3",
            "assay_chembl_id": "ASSAY1",
            "molecule_chembl_id": "M3",
            "standard_relation": "=",
            "standard_units": "nM",
            "standard_value": "200",
            "standard_type": "Ki",
            "pchembl_value": None,
        },
    ]
    molecule_payloads = {
        "M1": {"molecule_structures": None, "max_phase": None},
        "M2": {
            "molecule_structures": {"canonical_smiles": "CCC"},
            "max_phase": "unknown",
        },
        "M3": None,
    }

    def fake_paginated(
        url,
        params,
        headers,
        timeout,
        retries,
        result_key,
        session=None,
        limit=200,
        max_pages=50,
        request_logger=None,
    ):
        if "target" in url:
            yield {"target_chembl_id": "T1"}
            return
        if "assay" in url:
            yield {"assay_chembl_id": "A1"}
            return
        if "activity" in url:
            for activity in activities:
                yield activity

    def fake_get_json(
        url, params, headers, timeout, retries, session=None, request_log=None
    ):
        if "molecule" in url:
            mol_id = url.split("/")[-1].split(".")[0]
            return molecule_payloads.get(mol_id)
        return {}

    def fake_select(
        uniprot,
        unp_start,
        unp_end,
        keywords,
        prefer_single_protein=True,
        audit=None,
        **kwargs,
    ):
        return ["T1"], {
            "selected_ids": ["T1"],
            "returned_count": 1,
            "mode": "metadata_fallback",
            "reason": "test",
        }

    monkeypatch.setattr(es, "_chembl_paginated", fake_paginated)
    monkeypatch.setattr(es, "_get_json", fake_get_json)
    monkeypatch.setattr(es, "select_chembl_targets_for_chain", fake_select)

    labels, meta = es.fetch_chembl_labeled_smiles(
        "P12345",
        "TEST",
        cache_dir=tmp_path,
        timeout=1,
        retries=0,
        activity_types=["Ki"],
        chembl_max_phase=1,
        label_thresholds={
            "pchembl_strong": 7.0,
            "pchembl_weak": 5.0,
            "standard_value_nm_strong": 100.0,
            "standard_value_nm_weak": 10000.0,
        },
        debug=True,
        debug_max_ids=2,
        debug_rejection_samples_max=2,
    )

    assert labels["strong"] == ["CCC"]
    skips = meta["debug"]["skips"]
    assert skips.get("skip_molecule_structures_null", 0) >= 1
    assert skips.get("skip_payload_not_dict", 0) >= 1
    assert skips.get("skip_phase_unknown", 0) >= 1
    assert skips.get("skip_molecule_no_smiles", 0) >= 1

    telemetry = meta["telemetry"]["labeling_sanity"]
    assert telemetry.get("n_rejected_by_value_threshold") == 1
    samples = telemetry.get("rejected_value_threshold_samples") or []
    assert len(samples) == 1
    sample = samples[0]
    assert sample["standard_relation"] == ">"
    assert sample["molecule_chembl_id"] == "M1"
    assert sample["standard_value"] == "50"
