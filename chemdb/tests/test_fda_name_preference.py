import csv

from prep_ligands.metabolite_resolver import load_library_index


def _write_mapping_csv(path, rows):
    fieldnames = [
        "scheme",
        "file_num",
        "path",
        "pubchem_name",
        "pubchem_iupac_name",
        "pubchem_record_title",
        "display_name",
        "rxnorm_generic_name",
        "pubchem_synonyms",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_prefers_display_over_pubchem_name(tmp_path):
    mapping_path = tmp_path / "mapping.csv"
    long_iupac = (
        "4,4-dimethyl-2-(1H-1,2,4-triazol-1-ylmethyl)-"
        "1-(4-chlorophenyl)pentan-3-ol"
    )
    _write_mapping_csv(
        mapping_path,
        [
            {
                "scheme": "rdk",
                "file_num": "704",
                "path": "/x/rdk_0000704.pdbqt",
                "pubchem_name": long_iupac,
                "pubchem_iupac_name": long_iupac,
                "pubchem_record_title": "Itraconazole, (R)-(-)-",
                "display_name": "itraconazole (r)",
                "rxnorm_generic_name": "",
                "pubchem_synonyms": "",
            }
        ],
    )

    idx = load_library_index(str(mapping_path))
    rec = idx.id_to_rec["rdk_0000704"]

    assert rec.name == "itraconazole (r)"


def test_synonym_fallback_beats_iupac(tmp_path):
    mapping_path = tmp_path / "mapping.csv"
    long_iupac = (
        "6,7-dimethoxy-2-[4-(phenylmethyl)piperazin-1-yl]"
        "quinazolin-4-amine"
    )
    _write_mapping_csv(
        mapping_path,
        [
            {
                "scheme": "rdk",
                "file_num": "705",
                "path": "/x/rdk_0000705.pdbqt",
                "pubchem_name": long_iupac,
                "pubchem_iupac_name": long_iupac,
                "pubchem_record_title": "",
                "display_name": "",
                "rxnorm_generic_name": "",
                "pubchem_synonyms": "CHEMBL123; Itraconazole; 154003-19-7",
            }
        ],
    )

    idx = load_library_index(str(mapping_path))
    rec = idx.id_to_rec["rdk_0000705"]

    assert rec.name == "Itraconazole"
