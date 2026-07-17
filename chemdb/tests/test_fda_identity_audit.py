from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, inchi

from prep_ligands import fda_identity_audit as audit
from prep_ligands.fda_identity_audit import (
    FDAIdentityAuditConfig,
    run_fda_identity_audit,
)


FIELDS = [
    "scheme",
    "file_num",
    "path",
    "match_method",
    "sdf_index",
    "sdf_title",
    "smiles",
    "inchikey",
    "id",
    "cas",
    "display_name",
    "pubchem_record_title",
    "pubchem_cid_resolved",
    "pubchem_synonyms",
    "drugcentral_id",
    "pdbqt_heavy_atoms",
    "sdf_heavy_atoms",
]


def _write_mapping(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def _row(index: int, smiles: str, name: str, **overrides: str) -> dict[str, str]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    exact_key = inchi.MolToInchiKey(mol)
    row = {
        "scheme": "rdk",
        "file_num": str(index),
        "path": f"/stale/prepped_ligands/rdk_{index:07d}.pdbqt",
        "match_method": "index-fallback",
        "sdf_index": str(index + 1),
        "sdf_title": f"source_{index + 1}",
        "smiles": smiles,
        "inchikey": exact_key,
        "id": f"DC{index + 1}",
        "cas": "",
        "display_name": name,
        "pubchem_record_title": name,
        "pubchem_cid_resolved": str(1000 + index),
        "pubchem_synonyms": name,
        "drugcentral_id": "",
        "pdbqt_heavy_atoms": str(mol.GetNumHeavyAtoms()),
        "sdf_heavy_atoms": str(mol.GetNumHeavyAtoms()),
    }
    row.update(overrides)
    return row


def _write_pdbqt(path: Path, smiles: str) -> None:
    path.write_text(f"REMARK SMILES {smiles}\nTORSDOF 0\n", encoding="utf-8")


def _write_coordinate_pdbqt(path: Path, smiles: str) -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=17) == 0
    path.write_text(Chem.MolToPDBBlock(mol), encoding="utf-8")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_source_sdf(
    path: Path,
    rows: list[dict[str, str]],
    *,
    source_names: list[str] | None = None,
) -> None:
    writer = Chem.SDWriter(str(path))
    for index, row in enumerate(rows, start=1):
        mol = Chem.MolFromSmiles(row["smiles"])
        assert mol is not None
        source_name = (
            source_names[index - 1]
            if source_names is not None
            else row.get("display_name") or f"source_{index}"
        )
        mol.SetProp("_Name", source_name)
        mol.SetProp("ID", row["id"])
        mol.SetProp("PREFERRED_NAME", source_name)
        writer.write(mol)
    writer.close()
    metadata = path.parent / "FDA_Approved.csv"
    metadata.write_text(
        "".join(
            f"{row['id']},{(source_names or [item['display_name'] for item in rows])[index]}\n"
            for index, row in enumerate(rows)
        ),
        encoding="utf-8",
    )
    manifest = {
        "schema_version": 1,
        "filter_version": "drugcentral-fda-exact-v1",
        "filtered_sdf": str(path),
        "filtered_sdf_sha256": _sha256(path),
        "kept_records": len(rows),
        "approval_metadata": [{"path": str(metadata), "sha256": _sha256(metadata)}],
        "approval_key_counts": {
            "drugcentral_ids": len(rows),
            "names": len(rows),
            "cas_numbers": 0,
            "inchi_keys": 0,
        },
    }
    (path.parent / "fda_filter_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_structure_audit_separates_forms_nonmedication_and_mismatches(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    library = tmp_path / "library"
    output = tmp_path / "audit"
    classified = tmp_path / "classified.csv"
    source_sdf = tmp_path / "source.sdf"
    library.mkdir()

    exact = "CCO"
    salt = "CNCCC(c1ccccc1)Oc1ccc(C(F)(F)F)cc1.[Cl-]"
    salt_parent = "CNCCC(c1ccccc1)Oc1ccc(C(F)(F)F)cc1"
    additive = "O=[Ti]=O"
    mismatch_mapping = "c1ccccc1"
    mismatch_actual = "c1ccncc1"
    unverified = "CC(=O)O"
    combination = "c1ccccc1CCCCCCCCCC.c1ccccc1CCCCCCCCCC"
    rows = [
        _row(0, exact, "ethanol"),
        _row(1, salt, "example hydrochloride"),
        _row(2, additive, "titanium dioxide"),
        _row(3, mismatch_mapping, "benzene"),
        _row(4, unverified, "acetic acid"),
        _row(5, combination, "combination example"),
    ]
    _write_mapping(mapping, rows)
    original_mapping = mapping.read_bytes()

    _write_source_sdf(source_sdf, rows)

    _write_pdbqt(library / "rdk_0000000.pdbqt", exact)
    _write_pdbqt(library / "rdk_0000001.pdbqt", salt_parent)
    _write_pdbqt(library / "rdk_0000002.pdbqt", additive)
    _write_pdbqt(library / "rdk_0000003.pdbqt", mismatch_actual)
    _write_pdbqt(library / "rdk_0000005.pdbqt", combination)
    classified.write_text(
        "rdk_id,display_name,substance_class,status,note\n"
        "rdk_0000002,titanium dioxide,color_additive,valid,fixture\n",
        encoding="utf-8",
    )

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source_sdf,
            library_dir=library,
            output_dir=output,
            classified_nonmed_csv=classified,
        )
    )

    assert mapping.read_bytes() == original_mapping
    assert result.total_rows == 6
    assert result.blocker_count == 1
    by_id = {row["rdk_id"]: row for row in _read_csv(result.row_audit_csv)}
    assert by_id["rdk_0000000"]["identity_verdict"] == "verified_exact"
    assert by_id["rdk_0000000"]["fda_drug_eligible"] == "true"
    assert by_id["rdk_0000001"]["identity_verdict"] == "verified_parent"
    assert by_id["rdk_0000001"]["substance_class"] == "salt_or_solvate"
    assert by_id["rdk_0000001"]["fda_drug_eligible"] == "true"
    assert by_id["rdk_0000002"]["substance_class"] == "additive_or_nonmedication"
    assert by_id["rdk_0000002"]["fda_drug_eligible"] == "false"
    assert by_id["rdk_0000003"]["identity_verdict"] == "probable_mismatch"
    assert by_id["rdk_0000004"]["identity_verdict"] == "unverifiable"
    assert "ordinal_index_fallback_unverified" in by_id["rdk_0000004"]["reason_codes"]
    assert by_id["rdk_0000005"]["substance_class"] == "multi_active_combination"
    assert by_id["rdk_0000005"]["fda_drug_eligible"] == "false"

    quarantined = _read_csv(result.quarantine_manifest_csv)
    assert [row["rdk_id"] for row in quarantined] == ["rdk_0000003"]
    eligible = _read_csv(result.eligible_manifest_csv)
    assert {row["source_rdk_ids"] for row in eligible} == {
        "rdk_0000000",
        "rdk_0000001",
    }
    catalog = _read_csv(result.source_catalog_csv)
    assert len(catalog) == 6
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["schema_version"] == 4
    assert summary["policy"]["canonical_mapping_mutated"] is False
    assert summary["policy"]["index_fallback_is_verified"] is False
    assert summary["substance_class_counts"]["salt_or_solvate"] == 1
    assert len(result.evidence_jsonl.read_text(encoding="utf-8").splitlines()) == 6


def test_mapping_requires_a_ligand_identifier_column(tmp_path: Path) -> None:
    mapping = tmp_path / "invalid.csv"
    mapping.write_text("smiles,display_name\nCCO,ethanol\n", encoding="utf-8")

    with pytest.raises(ValueError, match="no ligand identifier column"):
        run_fda_identity_audit(
            FDAIdentityAuditConfig(mapping_csv=mapping, output_dir=tmp_path / "out")
        )


def test_structure_linked_wrong_name_is_quarantined(tmp_path: Path) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    row = _row(
        0,
        "CCO",
        "aspirin",
    )
    aspirin = _row(1, "CC(=O)Oc1ccccc1C(=O)O", "aspirin")
    _write_mapping(mapping, [row])
    _write_source_sdf(source, [row, aspirin], source_names=["ethanol", "aspirin"])
    _write_pdbqt(library / "rdk_0000000.pdbqt", "CCO")

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )

    audited = _read_csv(result.row_audit_csv)[0]
    assert audited["identity_verdict"] == "probable_mismatch"
    assert audited["name_corroborated"] == "false"
    assert "probable_name_structure_mismatch" in audited["reason_codes"]
    assert _read_csv(result.quarantine_manifest_csv)[0]["rdk_id"] == "rdk_0000000"


def test_stereo_conflicts_and_combination_parent_loss_are_not_eligible(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    specified = "C[C@H](O)C(=O)O"
    opposite = "C[C@@H](O)C(=O)O"
    unspecified = "CC(O)C(=O)O"
    combination = "CC(=O)NC1=CC=C(C=C1)O.Cn1c(=O)c2c(ncn2C)n(C)c1=O"
    caffeine = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
    rows = [
        _row(0, specified, "lactic acid"),
        _row(1, specified, "lactic acid"),
        _row(2, combination, "acetaminophen caffeine combination"),
    ]
    _write_mapping(mapping, rows)
    _write_source_sdf(source, rows)
    _write_pdbqt(library / "rdk_0000000.pdbqt", opposite)
    _write_pdbqt(library / "rdk_0000001.pdbqt", unspecified)
    _write_pdbqt(library / "rdk_0000002.pdbqt", caffeine)

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    by_id = {row["rdk_id"]: row for row in _read_csv(result.row_audit_csv)}

    assert by_id["rdk_0000000"]["identity_verdict"] == "probable_mismatch"
    assert by_id["rdk_0000001"]["identity_verdict"] == "ambiguous"
    assert by_id["rdk_0000002"]["substance_class"] == "multi_active_combination"
    assert by_id["rdk_0000002"]["identity_verdict"] == "ambiguous"
    assert all(row["fda_drug_eligible"] == "false" for row in by_id.values())


def test_stable_source_id_beats_shifted_ordinal_and_pdbqt_key_conflict_blocks_verification(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    target = _row(0, "CCO", "ethanol", sdf_index="1", id="TARGET")
    decoy = _row(9, "c1ccccc1", "benzene", sdf_index="1", id="DECOY")
    _write_mapping(mapping, [target])
    _write_source_sdf(source, [decoy, target])
    wrong_key = inchi.MolToInchiKey(Chem.MolFromSmiles("c1ccccc1"))
    (library / "rdk_0000000.pdbqt").write_text(
        f"REMARK SMILES CCO\nREMARK INCHIKEY {wrong_key}\nTORSDOF 0\n",
        encoding="utf-8",
    )

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["source_match_method"] == "stable_id_inchikey"
    assert audited["approval_verified"] == "true"
    assert audited["identity_verdict"] == "unverifiable"
    assert audited["fda_drug_eligible"] == "false"
    assert "pdbqt_stated_inchikey_conflict" in audited["reason_codes"]


def test_unmanifested_source_cannot_prove_fda_approval(tmp_path: Path) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    row = _row(0, "CCO", "ethanol")
    _write_mapping(mapping, [row])
    _write_source_sdf(source, [row])
    (tmp_path / "fda_filter_manifest.json").unlink()
    _write_pdbqt(library / "rdk_0000000.pdbqt", "CCO")

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["identity_verdict"] == "verified_exact"
    assert audited["approval_verified"] == "false"
    assert audited["name_corroborated"] == "false"
    assert audited["fda_drug_eligible"] == "false"
    assert "approval_manifest_missing" in audited["reason_codes"]


def test_conflicting_stable_identifiers_are_ambiguous_not_quarantined(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    wrong_id = _row(8, "c1ccccc1", "benzene", id="CROSSWIRED")
    target = _row(0, "CCO", "ethanol", id="TARGET")
    mapped = dict(target, id="CROSSWIRED")
    _write_mapping(mapping, [mapped])
    _write_source_sdf(source, [wrong_id, target])
    _write_pdbqt(library / "rdk_0000000.pdbqt", "CCO")

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["source_match_method"] == "conflicting_stable_identifiers"
    assert audited["identity_verdict"] == "ambiguous"
    assert audited["approval_verified"] == "false"
    assert _read_csv(result.quarantine_manifest_csv) == []


def test_missing_name_is_review_only_not_probable_mismatch(tmp_path: Path) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    row = _row(0, "CCO", "")
    _write_mapping(mapping, [row])
    _write_source_sdf(source, [row], source_names=["ethanol"])
    _write_pdbqt(library / "rdk_0000000.pdbqt", "CCO")

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["identity_verdict"] == "verified_exact"
    assert audited["fda_drug_eligible"] == "false"
    assert "display_name_missing" in audited["reason_codes"]
    assert _read_csv(result.quarantine_manifest_csv) == []


def test_salt_parent_pdbqt_and_exact_source_sidecar_are_consistent(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    salt = "CNCCC(c1ccccc1)Oc1ccc(C(F)(F)F)cc1.[Cl-]"
    parent = "CNCCC(c1ccccc1)Oc1ccc(C(F)(F)F)cc1"
    row = _row(0, salt, "example hydrochloride")
    _write_mapping(mapping, [row])
    _write_source_sdf(source, [row])
    salt_key = inchi.MolToInchiKey(Chem.MolFromSmiles(salt))
    pdbqt = library / "rdk_0000000.pdbqt"
    pdbqt.write_text(
        f"REMARK SMILES {parent}\nREMARK INCHIKEY {salt_key}\nTORSDOF 0\n",
        encoding="utf-8",
    )
    pdbqt.with_suffix(".ligprep_source.json").write_text(
        json.dumps(
            {
                "chemistry_authoritative": True,
                "source_sdf": str(source),
                "source_record_index": 1,
                "source_record_name": "example hydrochloride",
            }
        ),
        encoding="utf-8",
    )

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["identity_verdict"] == "verified_parent"
    assert audited["approval_verified"] == "true"
    assert audited["fda_drug_eligible"] == "true"
    assert "pdbqt_stated_inchikey_conflict" not in audited["reason_codes"]
    assert "prepared_candidate_structure_disagreement" not in audited["reason_codes"]


@pytest.mark.parametrize(
    ("display_name", "source_name", "smiles"),
    [
        (
            "hydrochloride",
            "metformin hydrochloride",
            "CN(C)C(=[NH2+])NC(=N)N.[Cl-]",
        ),
        (
            "aspirin and caffeine",
            "aspirin",
            "CC(=O)Oc1ccccc1C(=O)O",
        ),
        (
            "hydrochloride",
            "hydrochloride",
            "CN(C)C(=[NH2+])NC(=N)N.[Cl-]",
        ),
        (
            "aspirin and caffeine",
            "aspirin and caffeine",
            "CC(=O)Oc1ccccc1C(=O)O",
        ),
        (
            "D lactic acid",
            "L lactic acid",
            "CC(O)C(=O)O",
        ),
    ],
)
def test_component_or_formulation_subphrases_do_not_corroborate_drug_names(
    tmp_path: Path,
    display_name: str,
    source_name: str,
    smiles: str,
) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    row = _row(0, smiles, display_name)
    _write_mapping(mapping, [row])
    _write_source_sdf(source, [row], source_names=[source_name])
    _write_pdbqt(library / "rdk_0000000.pdbqt", smiles)

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["name_corroborated"] == "false"
    assert audited["fda_drug_eligible"] == "false"
    assert audited["identity_verdict"] != "probable_mismatch"
    assert _read_csv(result.quarantine_manifest_csv) == []


def test_unlisted_valid_alias_is_review_only_not_quarantined(tmp_path: Path) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    row = _row(0, "CC(=O)NC1=CC=C(C=C1)O", "paracetamol")
    _write_mapping(mapping, [row])
    _write_source_sdf(source, [row], source_names=["acetaminophen"])
    _write_pdbqt(library / "rdk_0000000.pdbqt", row["smiles"])

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["identity_verdict"] == "verified_exact"
    assert audited["name_corroborated"] == "false"
    assert audited["fda_drug_eligible"] == "false"
    assert _read_csv(result.quarantine_manifest_csv) == []


def test_multi_active_source_cannot_approve_one_parent_component(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    caffeine = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
    combination = f"CC(=O)NC1=CC=C(C=C1)O.{caffeine}"
    mapped = _row(0, caffeine, "acetaminophen and caffeine", id="COMBO")
    source_row = _row(9, combination, "acetaminophen and caffeine", id="COMBO")
    _write_mapping(mapping, [mapped])
    _write_source_sdf(source, [source_row])
    _write_pdbqt(library / "rdk_0000000.pdbqt", caffeine)

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["identity_verdict"] == "ambiguous"
    assert audited["approval_verified"] == "false"
    assert audited["fda_drug_eligible"] == "false"
    assert "approved_source_parent_is_mixture" in audited["reason_codes"]


def test_false_string_sidecar_authority_is_rejected(tmp_path: Path) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    row = _row(0, "CCO", "ethanol")
    _write_mapping(mapping, [row])
    _write_source_sdf(source, [row])
    pdbqt = library / "rdk_0000000.pdbqt"
    pdbqt.write_text("TORSDOF 0\n", encoding="utf-8")
    pdbqt.with_suffix(".ligprep_source.json").write_text(
        json.dumps(
            {
                "chemistry_authoritative": "false",
                "source_sdf": str(source),
                "source_record_index": 1,
                "source_record_name": "ethanol",
            }
        ),
        encoding="utf-8",
    )

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["identity_verdict"] == "unverifiable"
    assert "source_sidecar_not_authoritative" in audited["reason_codes"]


def test_unmatched_claimed_source_id_blocks_approval(tmp_path: Path) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    source_row = _row(0, "CCO", "ethanol", id="SOURCE_ID")
    mapped = dict(source_row, id="MISSING_ID")
    _write_mapping(mapping, [mapped])
    _write_source_sdf(source, [source_row])
    _write_pdbqt(library / "rdk_0000000.pdbqt", "CCO")

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["source_match_method"] == "stable_inchikey_unmatched_id"
    assert audited["identity_verdict"] == "ambiguous"
    assert audited["approval_verified"] == "false"
    assert "source_identifier_unmatched" in audited["reason_codes"]


def test_coordinate_only_pdbqts_report_exact_parent_and_disagreement(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    library = tmp_path / "library"
    library.mkdir()
    rows = [
        _row(0, "CCO", "ethanol"),
        _row(1, "CCN.[Cl-]", "ethylamine hydrochloride"),
        _row(2, "c1ccccc1", "benzene"),
        _row(3, "CC(=O)O", "acetic acid"),
    ]
    _write_mapping(mapping, rows)
    _write_coordinate_pdbqt(library / "rdk_0000000.pdbqt", "CCO")
    _write_coordinate_pdbqt(library / "rdk_0000001.pdbqt", "CCN")
    _write_coordinate_pdbqt(library / "rdk_0000002.pdbqt", "c1ccncc1")

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    by_id = {row["rdk_id"]: row for row in _read_csv(result.row_audit_csv)}

    assert by_id["rdk_0000000"]["identity_verdict"] == ("connectivity_consistent_exact")
    assert by_id["rdk_0000001"]["identity_verdict"] == (
        "connectivity_consistent_parent"
    )
    assert by_id["rdk_0000002"]["identity_verdict"] == "ambiguous"
    assert (
        "pdbqt_connectivity_topology_disagreement"
        in (by_id["rdk_0000002"]["reason_codes"])
    )
    assert by_id["rdk_0000003"]["identity_verdict"] == "unverifiable"
    assert all(row["fda_drug_eligible"] == "false" for row in by_id.values())

    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))
    assert summary["prepared_pdbqt_connectivity_counts"] == {
        "connectivity_consistent_exact": 1,
        "connectivity_consistent_parent": 1,
        "connectivity_disagrees": 1,
        "not_supplied": 1,
    }
    assert summary["prepared_connectivity_consistent_rows"] == 2
    assert summary["prepared_connectivity_evaluable_rows"] == 3
    assert summary["mapping_structure_parsed_rows"] == 4
    assert summary["high_fidelity_prepared_structure_rows"] == 0


def test_unmatched_id_prevents_wrong_source_structure_quarantine(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    library.mkdir()
    mapped = _row(
        0,
        "CCO",
        "ethanol",
        id="WRONG_MATCH",
        drugcentral_id="MISSING_CLAIM",
    )
    wrong_source = _row(9, "c1ccccc1", "benzene", id="WRONG_MATCH")
    _write_mapping(mapping, [mapped])
    _write_source_sdf(source, [wrong_source])
    _write_pdbqt(library / "rdk_0000000.pdbqt", "CCO")

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["source_match_method"] == "stable_id_unmatched_id"
    assert audited["identity_verdict"] == "ambiguous"
    assert audited["approval_verified"] == "false"
    assert "source_identifier_unmatched" in audited["reason_codes"]
    assert _read_csv(result.quarantine_manifest_csv) == []


def test_library_index_prefers_canonical_pdbqt_over_temporary_artifacts(
    tmp_path: Path,
) -> None:
    canonical = tmp_path / "rdk_0000000.pdbqt"
    temporary = tmp_path / "rdk_0000000.obabel_tmp.pdbqt"
    canonical.write_text("TORSDOF 0\n", encoding="utf-8")
    temporary.write_text("", encoding="utf-8")

    indexed = audit._library_index(tmp_path)

    assert indexed["rdk_0000000"] == [canonical]


def test_connectivity_uses_openbabel_fallback_when_coordinate_parse_fails(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    pdbqt = tmp_path / "rdk_0000000.pdbqt"
    pdbqt.write_text("ATOM malformed\n", encoding="utf-8")
    mapping = audit._structure_from_smiles("CCO", "mapping")
    monkeypatch.setattr(audit, "_pdbqt_coordinate_mol", lambda _lines: None)
    monkeypatch.setattr(
        audit,
        "_obabel_pdbqt_mol",
        lambda _path: Chem.MolFromSmiles("CCO"),
    )

    evidence = audit._pdbqt_connectivity_evidence(pdbqt, mapping)

    assert evidence.verdict == "connectivity_consistent_exact"
    assert evidence.method == "openbabel_pdbqt_graph"


def test_mapping_inchikey_layer_difference_is_ambiguous_but_graph_change_blocks(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    library = tmp_path / "library"
    library.mkdir()
    specified = "C[C@H](O)C(=O)O"
    opposite_key = inchi.MolToInchiKey(Chem.MolFromSmiles("C[C@@H](O)C(=O)O"))
    benzene_key = inchi.MolToInchiKey(Chem.MolFromSmiles("c1ccccc1"))
    rows = [
        _row(0, specified, "lactic acid", inchikey=opposite_key),
        _row(1, "CCO", "ethanol", inchikey=benzene_key),
    ]
    _write_mapping(mapping, rows)
    _write_pdbqt(library / "rdk_0000000.pdbqt", specified)
    _write_pdbqt(library / "rdk_0000001.pdbqt", "CCO")

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    by_id = {row["rdk_id"]: row for row in _read_csv(result.row_audit_csv)}

    assert by_id["rdk_0000000"]["identity_verdict"] == "ambiguous"
    assert (
        "mapping_stored_inchikey_stereo_or_protonation_mismatch"
        in (by_id["rdk_0000000"]["reason_codes"])
    )
    assert by_id["rdk_0000001"]["identity_verdict"] == "probable_mismatch"
    assert (
        "mapping_stored_inchikey_connectivity_mismatch"
        in (by_id["rdk_0000001"]["reason_codes"])
    )
    assert [row["rdk_id"] for row in _read_csv(result.quarantine_manifest_csv)] == [
        "rdk_0000001"
    ]


def test_disagreeing_duplicate_connectivity_downgrades_high_fidelity_match(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "mapping.csv"
    source = tmp_path / "source.sdf"
    library = tmp_path / "library"
    high_fidelity = library / "trusted"
    low_fidelity = library / "reconstructed"
    high_fidelity.mkdir(parents=True)
    low_fidelity.mkdir()
    row = _row(0, "CCO", "ethanol")
    _write_mapping(mapping, [row])
    _write_source_sdf(source, [row])
    _write_pdbqt(high_fidelity / "rdk_0000000.pdbqt", "CCO")
    _write_coordinate_pdbqt(low_fidelity / "rdk_0000000.pdbqt", "CCN")

    result = run_fda_identity_audit(
        FDAIdentityAuditConfig(
            mapping_csv=mapping,
            source_sdf=source,
            library_dir=library,
            output_dir=tmp_path / "audit",
        )
    )
    audited = _read_csv(result.row_audit_csv)[0]

    assert audited["identity_verdict"] == "ambiguous"
    assert audited["approval_verified"] == "true"
    assert audited["name_corroborated"] == "true"
    assert audited["fda_drug_eligible"] == "false"
    assert "prepared_high_fidelity_connectivity_conflict" in (audited["reason_codes"])
    assert _read_csv(result.quarantine_manifest_csv) == []


def test_audit_refuses_output_file_collision_with_mapping(tmp_path: Path) -> None:
    output_dir = tmp_path / "audit"
    output_dir.mkdir()
    mapping = output_dir / "fda_identity_audit.csv"
    _write_mapping(mapping, [_row(0, "CCO", "ethanol")])

    with pytest.raises(ValueError, match="collides with an input"):
        run_fda_identity_audit(
            FDAIdentityAuditConfig(mapping_csv=mapping, output_dir=output_dir)
        )
