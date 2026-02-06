from protein_prep import aliases_policy


def test_aliases_policy_normalization_and_constants() -> None:
    assert aliases_policy._normalize_resname("zn2") == "ZN"
    assert aliases_policy._normalize_resname("cl-") == "CL"
    assert aliases_policy._normalize_resname("hoh") == "HOH"

    assert "HOH" in aliases_policy._WATER_NAMES
    assert "ZN" in aliases_policy._ION_AUDIT_METALS
    assert "NA" in aliases_policy._ION_AUDIT_SIMPLE_IONS
    assert "MG" in aliases_policy._ELEM_CANON


def test_flatten_semicolons_is_stable() -> None:
    assert aliases_policy._flatten_semicolons(["A; b", " c "]) == ["A", "B", "C"]
