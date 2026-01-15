from DeepCoy_duds import external_sources as es


def test_chembl_labeled_fetcher_conflict_priority(monkeypatch, tmp_path):
    activity_payload = [
        {
            "pchembl_value": 4.5,
            "standard_type": "Ki",
            "standard_units": "nM",
            "standard_value": "12000",
            "standard_relation": "=",
            "molecule_chembl_id": "CHEMBL_NON",
        },
        {
            "pchembl_value": 7.2,
            "standard_type": "Ki",
            "standard_units": "nM",
            "standard_value": "10",
            "standard_relation": "=",
            "molecule_chembl_id": "CHEMBL_STRONG",
        },
        {
            "standard_type": "Ki",
            "standard_units": "nM",
            "standard_value": "500",
            "standard_relation": "=",
            "molecule_chembl_id": "CHEMBL_WEAK",
        },
        {
            "pchembl_value": 4.8,
            "standard_type": "Ki",
            "standard_units": "nM",
            "standard_value": "20000",
            "standard_relation": "=",
            "molecule_chembl_id": "CHEMBL_SHARED_NON",
        },
        {
            "standard_type": "Ki",
            "standard_units": "nM",
            "standard_value": "50",
            "standard_relation": "=",
            "molecule_chembl_id": "CHEMBL_SHARED_NON",
        },
    ]

    molecule_payload = {
        "CHEMBL_NON": {
            "molecule_structures": {"canonical_smiles": "CNon"},
            "max_phase": 3,
        },
        "CHEMBL_STRONG": {
            "molecule_structures": {"canonical_smiles": "CShared"},
            "max_phase": 3,
        },
        "CHEMBL_WEAK": {
            "molecule_structures": {"canonical_smiles": "CWeak"},
            "max_phase": 3,
        },
        "CHEMBL_SHARED_NON": {
            "molecule_structures": {"canonical_smiles": "CShared"},
            "max_phase": 3,
        },
    }

    def fake_get_json(
        url, params, headers, timeout, retries, session=None, request_log=None
    ):
        if "target" in url:
            return {
                "targets": [{"target_chembl_id": "T1"}],
                "page_meta": {
                    "total_count": 1,
                    "offset": params.get("offset", 0),
                    "limit": params.get("limit", 200),
                },
            }
        if "assay" in url:
            return {
                "assays": [{"assay_chembl_id": "A1"}],
                "page_meta": {
                    "total_count": 1,
                    "offset": params.get("offset", 0),
                    "limit": params.get("limit", 200),
                },
            }
        if "activity" in url:
            return {
                "activities": activity_payload,
                "page_meta": {
                    "total_count": len(activity_payload),
                    "offset": params.get("offset", 0),
                    "limit": params.get("limit", 500),
                },
            }
        if "molecule" in url:
            chembl_id = url.split("/")[-1].split(".")[0]
            return molecule_payload[chembl_id]
        raise AssertionError(f"Unexpected URL {url}")

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

    monkeypatch.setattr(es, "_get_json", fake_get_json)
    monkeypatch.setattr(es, "select_chembl_targets_for_chain", fake_select)

    labels, meta = es.fetch_chembl_labeled_smiles(
        "P12345",
        "TEST",
        tmp_path,
        timeout=5,
        retries=0,
        activity_types=["Ki", "Kd", "IC50", "EC50"],
        chembl_max_phase=None,
        label_thresholds={
            "pchembl_strong": 7.0,
            "pchembl_weak": 5.0,
            "standard_value_nm_strong": 100,
            "standard_value_nm_weak": 10000,
        },
    )

    assert "CShared" in labels["strong"]
    assert "CShared" not in labels["non"]
    assert "CWeak" in labels["weak"]
    assert "CNon" in labels["non"]
    assert meta["counts"]["targets"] == 1
    assert meta["counts"]["assays"] == 1
    assert meta["counts"]["activities"] == len(activity_payload)
