from chemdb.tools import fill_fda_mapping_names as tool


def test_pubchem_candidate_scoring_prefers_drug_like_synonym():
    candidates = [
        "924644-56-4",
        "C17H12Cl2F3N7O2S",
        "N3O8L7R66H",
        "3 hydroxy 2 phenyl 1h inden 1 one",
        "Cinnarizine",
    ]
    source_headings = {
        "924644-56-4": "Other Identifiers",
        "C17H12Cl2F3N7O2S": "Computed Descriptors",
        "N3O8L7R66H": "Other Identifiers",
        "3 hydroxy 2 phenyl 1h inden 1 one": "Synonyms",
        "Cinnarizine": "Drug and Medication Information",
    }

    selected = tool._select_pubchem_candidate(
        candidates,
        pubchem_iupac_name="",
        source_headings={key.lower(): value for key, value in source_headings.items()},
    )
    assert selected == "Cinnarizine"
