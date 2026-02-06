from chemdb.tools import fill_fda_mapping_names as tool


def test_name_quality_helpers_flag_low_quality_patterns():
    assert tool.is_cas_only_name("924644-56-4")
    assert tool.is_formula_like_name("C17H12Cl2F3N7O2S")
    assert tool.is_unii_like_name("N3O8L7R66H")
    assert tool.is_spaced_iupacish_name("3 hydroxy 2 phenyl 1h inden 1 one")


def test_name_quality_helpers_allow_common_drug_names():
    good_names = [
        "Aspirin",
        "Phenindione",
        "Acetaminophen",
        "Cinnarizine",
    ]
    for name in good_names:
        assert not tool.is_cas_only_name(name)
        assert not tool.is_formula_like_name(name)
        assert not tool.is_unii_like_name(name)
        assert not tool.is_spaced_iupacish_name(name)
        assert not tool.is_low_quality_display_name(name, strict=True)
