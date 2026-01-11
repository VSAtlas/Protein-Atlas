import json
from pathlib import Path

import pytest

from external_sources import (
    _parse_bindingdb_entries,
    filter_chembl_activity,
    parse_chembl_molecule,
    parse_iuphar_structure,
)


def load_fixture(name: str):
    here = Path(__file__).resolve().parent
    with open(here / "fixtures" / name, "r") as f:
        return json.load(f)


def test_bindingdb_parsing_filters_by_cutoff():
    entries = load_fixture("bindingdb_entries.json")
    smiles = _parse_bindingdb_entries(entries, affinity_cutoff_nm=100)
    assert "CCO" in smiles
    assert "CCC" not in smiles


def test_chembl_activity_and_molecule_parsing():
    activities = load_fixture("chembl_activity.json")
    good = [
        a
        for a in activities
        if filter_chembl_activity(a, cutoff_nm=100, allowed_types=["Ki", "Kd", "IC50"])
    ]
    assert len(good) == 1
    mol_json = load_fixture("chembl_molecule.json")
    smi = parse_chembl_molecule(mol_json, max_phase=4)
    assert smi == "C1=CC=CC=C1"
    mol_high = load_fixture("chembl_molecule_high_phase.json")
    smi_high = parse_chembl_molecule(mol_high, max_phase=4)
    assert smi_high is None


def test_iuphar_structure_parsing_smiles():
    struct = load_fixture("iuphar_structure.json")
    smi = parse_iuphar_structure(struct)
    assert smi == "CN1CCCC1"


@pytest.mark.skip(reason="RDKit-based InChI parsing may be unavailable in minimal env")
def test_iuphar_structure_parsing_inchi():
    struct = load_fixture("iuphar_structure_inchi.json")
    smi = parse_iuphar_structure(struct)
    assert smi
