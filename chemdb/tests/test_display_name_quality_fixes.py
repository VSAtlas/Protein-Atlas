from __future__ import annotations

import pytest

from chemdb.tools import fill_fda_mapping_names as tool


def _row_template(**overrides: str) -> dict[str, str]:
    row = {
        "file_num": "1",
        "display_name": "UNK_1",
        "generic_name": "",
        "inchikey": "",
        "remark_inchikey": "",
        "smiles_neutral": "",
        "smiles": "",
        "remark_smiles": "",
        "inchi": "",
        "remark_inchi": "",
        "pubchem_iupac_name": "",
        "pubchem_record_title": "",
        "pubchem_name": "",
        "pubchem_synonyms": "",
        "pubchem_cid_resolved": "",
        "pubchem_parent_cid": "",
        "pubchem_unii_list": "",
        "rxnorm_generic_name": "",
        "rxnorm_rxcui": "",
    }
    row.update(overrides)
    return row


def _run_row(row: dict[str, str], *, rxnorm_client: object | None = None) -> tool.RowResult:
    return tool._process_row(
        row_index=1,
        row=row,
        only_fix_bad_display_names=True,
        pubchem_client=None,
        rxnorm_client=rxnorm_client,  # type: ignore[arg-type]
        unichem_client=None,
        chembl_client=None,
        gsrs_client=None,
        drugcentral_index=tool.DrugCentralIndex(),
        enable_unichem=False,
        enable_pubchem_rerank=False,
        enable_gsrs=False,
    )


def test_inchikey_like_is_bad_and_overwritten():
    row = _row_template(
        display_name="VWMJHAFYPMOMGF-UHFFFAOYSA-N",
        rxnorm_generic_name="metformin",
    )
    result = _run_row(row)
    assert result.rxnorm_promoted
    assert row["display_name"] == "metformin"
    assert not tool.is_inchikey_string(row["display_name"])


def test_tokenized_inchikey_is_bad():
    row = _row_template(
        display_name="vwmjhafypmomgf uhfffaoysa n",
        rxnorm_generic_name="metformin",
    )
    result = _run_row(row)
    assert result.rxnorm_promoted
    assert row["display_name"] == "metformin"
    assert not tool.is_inchikey_string(row["display_name"])


def test_tautomer_canonicalization_allows_drugcentral_match_for_phenindione():
    if not tool._rdkit_available():
        pytest.skip("RDKit is required for tautomer canonicalization test.")

    smiles_enol = "O=C1C(c2ccccc2)=C(O)c2ccccc21"
    expected_inchikey = "NFBAXHOPROOJAW-UHFFFAOYSA-N"
    derived = tool.rdkit_inchikey_from_smiles(smiles_enol).upper()
    assert derived == expected_inchikey

    row = _row_template(display_name="UNK_4688", smiles=smiles_enol)
    index = tool.DrugCentralIndex(
        full={expected_inchikey: "Phenindione"},
        two_d={expected_inchikey[:14]: "Phenindione"},
    )
    result = tool._process_row(
        row_index=4688,
        row=row,
        only_fix_bad_display_names=True,
        pubchem_client=None,
        rxnorm_client=None,
        unichem_client=None,
        chembl_client=None,
        gsrs_client=None,
        drugcentral_index=index,
        enable_unichem=False,
        enable_pubchem_rerank=False,
        enable_gsrs=False,
    )

    assert result.drugcentral_resolved
    assert row["display_name"].lower() == "phenindione"


def test_salvage_semicolon_name_hint_resolves_acetaminophen(monkeypatch: pytest.MonkeyPatch):
    row = _row_template(
        display_name="UNK_Pharmace",
        smiles="powder; Acetaminophen",
    )

    def _fake_resolve_rxnorm(term: str, _client: object) -> tool.RxNormResult:
        if term.lower() == "acetaminophen":
            return tool.RxNormResult(rxcui="161", name="acetaminophen")
        return tool.RxNormResult()

    monkeypatch.setattr(tool, "resolve_rxnorm", _fake_resolve_rxnorm)
    result = _run_row(row, rxnorm_client=object())

    assert result.rxnorm_promoted
    assert row["display_name"] == "acetaminophen"
    assert not row["display_name"].upper().startswith("UNK_")
