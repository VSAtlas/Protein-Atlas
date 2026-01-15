from calibrator import uniprot_resolver


def test_uniprot_resolver_dbref_segments_for_6lu7():
    mapping = uniprot_resolver.resolve_chain_uniprot_segments(
        "6lu7", write_file=False
    )
    assert "A" in mapping
    entry = mapping["A"]
    assert entry["uniprot"] == "P0DTD1"
    assert entry["pdb_start"] == 1
    assert entry["pdb_end"] == 306
    assert entry["unp_start"] == 3264
    assert entry["unp_end"] == 3569
    assert entry["source"] == "dbref"
