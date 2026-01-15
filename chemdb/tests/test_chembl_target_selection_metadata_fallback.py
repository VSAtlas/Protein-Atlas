from DeepCoy_duds import external_sources as es


def test_chembl_target_selection_metadata_fallback(monkeypatch):
    targets = [
        {
            "target_chembl_id": "CHEMBL_MAIN",
            "target_type": "SINGLE PROTEIN",
            "pref_name": "SARS-CoV-2 main protease (3CLpro)",
        },
        {
            "target_chembl_id": "CHEMBL_POLY",
            "target_type": "POLYPROTEIN",
            "pref_name": "ORF1ab polyprotein",
        },
        {
            "target_chembl_id": "CHEMBL_OTHER",
            "target_type": "SINGLE PROTEIN",
            "pref_name": "Unrelated protease",
        },
    ]

    def fake_get_json_paged(
        url,
        params,
        headers,
        timeout,
        retries,
        session=None,
        item_key=None,
        max_pages=None,
        max_items=None,
        audit=None,
        source=None,
        purpose=None,
        request_logger=None,
    ):
        return targets, []

    monkeypatch.setattr(es, "_get_json_paged", fake_get_json_paged)
    selected, meta = es.select_chembl_targets_for_chain(
        "P0DTD1",
        None,
        None,
        ["3clpro", "main protease"],
        return_meta=True,
    )

    assert selected == ["CHEMBL_MAIN"]
    assert meta.get("mode") == "metadata_fallback"
    assert "keyword" in meta.get("reason", "") or meta.get("reason") == "keyword_match"
