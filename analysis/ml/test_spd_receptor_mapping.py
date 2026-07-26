from __future__ import annotations

import copy
import json

import pandas as pd
import pytest

from analysis.ml.spd_receptor_mapping import (
    DEFAULT_CONTRACT_PATH,
    OUTPUT_COLUMNS,
    POLICY_VERSION,
    ReceptorMappingContractError,
    apply_receptor_target_contract,
    apply_spd_receptor_mapping_contract,
    load_spd_receptor_mapping_contract,
    validate_spd_receptor_mapping_contract,
)


def _payload() -> dict[str, object]:
    return json.loads(DEFAULT_CONTRACT_PATH.read_text(encoding="utf-8"))


def _mapping(payload: dict[str, object], pdb_id: str) -> dict[str, object]:
    mappings = payload["mappings"]
    assert isinstance(mappings, list)
    return next(row for row in mappings if row["pdb_id"] == pdb_id)


def test_default_contract_encodes_locked_decisions_and_provenance() -> None:
    contract = load_spd_receptor_mapping_contract()
    mappings = contract.by_pdb()

    assert contract.policy_version == POLICY_VERSION
    assert set(mappings) == {"9I52", "6HUJ", "8HCQ"}
    revised = mappings["9I52"]
    assert (
        revised.legacy_target_gene,
        revised.legacy_target_uniprot,
        revised.revised_target_gene,
        revised.revised_target_uniprot,
    ) == ("ADRB2", "P07550", "DRD1", "P21728")
    assert revised.selected_chain == revised.site_chain == "R"
    assert revised.prepared_chains == ("A", "B", "G", "R")
    assert revised.site_relevant_chains == ("R",)
    assert revised.ligand_site_evidence == {
        "ligand": "A1I",
        "ligand_atom_count": 26,
        "contact_cutoff_angstrom": 5.0,
        "target_gene": "DRD1",
        "target_uniprot": "P21728",
        "chain": "R",
        "ligand_atoms_within_cutoff": 26,
        "atom_pairs_within_cutoff": 242,
        "minimum_distance_angstrom": 2.9937418,
    }
    assert revised.strict_eligible is True
    assert revised.requires_target_derived_rebuild is True

    for pdb_id in ("6HUJ", "8HCQ"):
        unresolved = mappings[pdb_id]
        assert unresolved.status == "unresolved"
        assert unresolved.strict_eligible is False
        assert unresolved.revised_target_id is None
        assert unresolved.revised_target_gene is None
        assert unresolved.revised_target_uniprot is None
        assert unresolved.ligand_site_evidence is None

    required_evidence = {
        "audit_mapping_review",
        "pdbe_sifts_snapshot",
        "uniprot_snapshot",
        "local_structure",
        "prepared_receptor",
    }
    assert all(
        {item.kind for item in mapping.evidence} == required_evidence
        for mapping in mappings.values()
    )


def test_apply_adds_only_namespaced_contract_columns_and_preserves_input() -> None:
    frame = pd.DataFrame(
        {
            "pdb_id": ["9i52", "6HUJ", "8HCQ", "1ABC", pd.NA],
            "target_gene": ["ADRB2", "GABRA1", "EDNRA", "OTHER", "MISSING"],
            "spd_binding_label": [1, 0, 1, 0, 1],
        }
    )

    result = apply_spd_receptor_mapping_contract(frame)

    assert result["target_gene"].tolist() == frame["target_gene"].tolist()
    assert result["spd_binding_label"].tolist() == frame["spd_binding_label"].tolist()
    assert set(result.columns) - set(frame.columns) == set(OUTPUT_COLUMNS)
    assert result.loc[0, "spd_receptor_mapping_revised_target_gene"] == "DRD1"
    assert result.loc[0, "spd_receptor_mapping_revised_target_uniprot"] == "P21728"
    assert result.loc[0, "spd_receptor_mapping_strict_eligible"] == True  # noqa: E712
    assert result.loc[0, "spd_receptor_mapping_requires_target_derived_rebuild"] == True  # noqa: E712
    assert result.loc[1, "spd_receptor_mapping_strict_eligible"] == False  # noqa: E712
    assert pd.isna(result.loc[1, "spd_receptor_mapping_revised_target_gene"])
    assert pd.isna(result.loc[2, "spd_receptor_mapping_revised_target_uniprot"])
    assert result.loc[1, "spd_receptor_mapping_exclusion_reason"]
    assert result.loc[3, "spd_receptor_mapping_status"] == "unaffected"
    assert result.loc[3, "spd_receptor_mapping_strict_eligible"] == True  # noqa: E712
    assert result.loc[3, "spd_receptor_mapping_requires_target_derived_rebuild"] == False  # noqa: E712
    assert result.loc[4, "spd_receptor_mapping_policy_version"] == POLICY_VERSION


def test_reconciled_wrapper_preserves_legacy_and_applies_strict_identity() -> None:
    source = pd.DataFrame(
        {
            "pdb_id": ["9I52", "6HUJ", "8HCQ", "1ABC"],
            "target_id": ["ADRB2", "GABRA1", "EDNRA", "OTHER"],
            "target_gene": ["ADRB2", "GABRA1", "EDNRA", "OTHER"],
            "target_uniprot": ["P07550", "P14867", "P25101", "P00001"],
            "spd_binding_label": [1, 0, 1, 0],
        }
    )

    legacy, legacy_summary = apply_receptor_target_contract(source, mode="legacy")
    strict, strict_summary = apply_receptor_target_contract(source, mode="strict")

    pd.testing.assert_frame_equal(legacy, source)
    assert strict.loc[0, ["target_id", "target_gene", "target_uniprot"]].tolist() == [
        "DRD1",
        "DRD1",
        "P21728",
    ]
    assert strict.loc[1, ["target_id", "target_gene", "target_uniprot"]].tolist() == [
        "GABRA1",
        "GABRA1",
        "P14867",
    ]
    assert strict.loc[2, ["target_id", "target_gene", "target_uniprot"]].tolist() == [
        "EDNRA",
        "EDNRA",
        "P25101",
    ]
    assert strict.loc[3, ["target_id", "target_gene", "target_uniprot"]].tolist() == [
        "OTHER",
        "OTHER",
        "P00001",
    ]
    assert strict["spd_binding_label"].tolist() == source["spd_binding_label"].tolist()
    assert legacy_summary["rows_in_contract"] == strict_summary["rows_in_contract"] == 3
    assert strict_summary["rows_strict_eligible"] == 2
    assert strict_summary["rows_strict_ineligible"] == 2
    assert strict_summary["rows_requiring_target_derived_rebuild"] == 1


def test_validator_rejects_duplicate_pdb() -> None:
    payload = _payload()
    mappings = payload["mappings"]
    assert isinstance(mappings, list)
    mappings.append(copy.deepcopy(mappings[0]))

    with pytest.raises(ReceptorMappingContractError, match="duplicate PDB"):
        validate_spd_receptor_mapping_contract(payload)


def test_validator_rejects_missing_provenance() -> None:
    payload = _payload()
    row = _mapping(payload, "9I52")
    evidence = row["evidence"]
    assert isinstance(evidence, list)
    evidence[:] = [item for item in evidence if item["kind"] != "prepared_receptor"]

    with pytest.raises(ReceptorMappingContractError, match="invalid evidence kinds"):
        validate_spd_receptor_mapping_contract(payload)


@pytest.mark.parametrize(
    ("pdb_id", "field", "value", "message"),
    [
        ("9I52", "revised_target_uniprot", "P07550", "unexpected revised_target_uniprot"),
        (
            "6HUJ",
            "revised_target_gene",
            "GABRA1",
            "unresolved mapping cannot declare a revised identity",
        ),
        (
            "8HCQ",
            "strict_eligible",
            True,
            "unresolved mapping cannot be strict eligible",
        ),
    ],
)
def test_validator_rejects_identity_or_eligibility_drift(
    pdb_id: str,
    field: str,
    value: object,
    message: str,
) -> None:
    payload = _payload()
    _mapping(payload, pdb_id)[field] = value

    with pytest.raises(ReceptorMappingContractError, match=message):
        validate_spd_receptor_mapping_contract(payload)


def test_validator_rejects_unexpected_pdb_identity() -> None:
    payload = _payload()
    _mapping(payload, "9I52")["pdb_id"] = "9ZZZ"

    with pytest.raises(ReceptorMappingContractError, match="unexpected PDB identity"):
        validate_spd_receptor_mapping_contract(payload)


def test_strict_application_rejects_conflicting_scoped_input_identity() -> None:
    source = pd.DataFrame(
        {
            "pdb_id": ["9I52"],
            "target_id": ["UNRELATED"],
            "target_gene": ["UNRELATED"],
            "target_uniprot": ["P00000"],
        }
    )

    with pytest.raises(ReceptorMappingContractError, match="input target_id conflicts"):
        apply_receptor_target_contract(source, mode="strict")

    legacy, _ = apply_receptor_target_contract(source, mode="legacy")
    pd.testing.assert_frame_equal(legacy, source)
