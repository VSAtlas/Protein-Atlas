from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from analysis.atlas_database.receptor_evidence import (
    CHEMISTRY_CLASSIFICATION_METHOD,
    ReceptorEvidenceError,
    audit_prepared_receptor_chemistry,
    build_receptor_evidence_records,
    write_receptor_evidence,
)


def _atom_line(
    serial: int,
    atom: str,
    residue: str,
    resseq: int,
    element: str,
    *,
    het: bool = True,
) -> str:
    record = "HETATM" if het else "ATOM  "
    return (
        f"{record}{serial:5d} {atom:^4s} {residue:>3s} A{resseq:4d}    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}"
        f"          {element:>2s}\n"
    )


def test_exact_prepared_content_controls_observed_apo_holo(tmp_path: Path) -> None:
    holo = tmp_path / "holo.pdbqt"
    holo.write_text(
        "".join(
            (
                _atom_line(1, "ZN", "ZN", 1, "ZN"),
                _atom_line(2, "C1", "FAD", 2, "C"),
                _atom_line(3, "C2", "FAD", 2, "C"),
                _atom_line(4, "CA", "ALA", 3, "C", het=False),
            )
        ),
        encoding="utf-8",
    )
    apo = tmp_path / "apo.pdbqt"
    apo.write_text(_atom_line(1, "CA", "ALA", 1, "C", het=False), encoding="utf-8")

    holo_result = audit_prepared_receptor_chemistry(
        holo,
        canonical_metals={"ZN"},
        canonical_cofactors={"FAD"},
    )
    apo_result = audit_prepared_receptor_chemistry(
        apo,
        canonical_metals={"ZN"},
        canonical_cofactors={"FAD"},
    )

    assert holo_result["classification_method"] == CHEMISTRY_CLASSIFICATION_METHOD
    assert holo_result["receptor_classification"] == "HOLO"
    assert holo_result["retained_metal_atom_count"] == 1
    assert holo_result["retained_cofactor_residue_count"] == 1
    assert holo_result["retained_cofactor_tokens"] == ["FAD"]
    assert len(holo_result["prepared_receptor_sha256"]) == 64
    assert apo_result["receptor_classification"] == "APO"
    assert apo_result["retained_metal_atom_count"] == 0
    assert apo_result["retained_cofactor_residue_count"] == 0


def test_invalid_and_ambiguous_content_never_becomes_apo(tmp_path: Path) -> None:
    empty = tmp_path / "empty.pdbqt"
    empty.write_text("", encoding="utf-8")
    heterogen_only = tmp_path / "heterogen_only.pdbqt"
    heterogen_only.write_text(
        _atom_line(1, "ZN", "ZN", 1, "ZN"),
        encoding="utf-8",
    )
    malformed = tmp_path / "malformed.pdbqt"
    malformed.write_text(
        "ATOM      1  CA  ALA A   1    not-coordinates\n",
        encoding="utf-8",
    )
    ambiguous = tmp_path / "ambiguous.pdbqt"
    ambiguous.write_text(
        "".join(
            (
                _atom_line(1, "CA", "ALA", 1, "C", het=False),
                _atom_line(2, "NA", "NA", 2, "NA"),
                _atom_line(3, "CL", "CL", 3, "CL"),
                _atom_line(4, "P", "ATP", 4, "P"),
            )
        ),
        encoding="utf-8",
    )

    empty_result = audit_prepared_receptor_chemistry(empty)
    heterogen_result = audit_prepared_receptor_chemistry(
        heterogen_only,
        canonical_metals={"ZN"},
        canonical_cofactors=set(),
    )
    malformed_result = audit_prepared_receptor_chemistry(malformed)
    ambiguous_result = audit_prepared_receptor_chemistry(
        ambiguous,
        canonical_metals=set(),
        canonical_cofactors=set(),
    )

    assert empty_result["receptor_classification"] == ""
    assert empty_result["chemistry_evidence_status"] == "unresolved"
    assert {"empty_file", "no_atom_records", "no_protein_coordinates"}.issubset(
        empty_result["validation_errors"]
    )
    assert heterogen_result["receptor_classification"] == ""
    assert "no_protein_coordinates" in heterogen_result["validation_errors"]
    assert malformed_result["receptor_classification"] == ""
    assert "malformed_atom_coordinates" in malformed_result["validation_errors"]
    assert ambiguous_result["receptor_classification"] == ""
    assert ambiguous_result["chemistry_evidence_status"] == "ambiguous"
    assert ambiguous_result["retained_ambiguous_cation_tokens"] == ["NA"]
    assert ambiguous_result["retained_halide_tokens"] == ["CL"]
    assert ambiguous_result["retained_nucleotide_like_tokens"] == ["ATP"]


def _release_fixture(tmp_path: Path) -> tuple[Path, Path]:
    run_id = "run-1"
    run_manifest = (
        tmp_path
        / "outputs"
        / "manifests"
        / run_id
        / "run_manifest.yaml"
    )
    processed = tmp_path / "outputs" / "processed_pdbs" / run_id
    holo_receptor = processed / "1AAA" / "receptor" / "1AAA.pdbqt"
    apo_receptor = processed / "2BBB" / "receptor" / "2BBB.pdbqt"
    holo_receptor.parent.mkdir(parents=True)
    apo_receptor.parent.mkdir(parents=True)
    holo_receptor.write_text(
        "".join(
            (
                _atom_line(1, "ZN", "ZN", 1, "ZN"),
                _atom_line(2, "CA", "ALA", 2, "C", het=False),
            )
        ),
        encoding="utf-8",
    )
    apo_receptor.write_text(
        _atom_line(1, "CA", "ALA", 1, "C", het=False), encoding="utf-8"
    )
    recorded_root = Path("/scratch/relocated/run-1")
    audit = (
        tmp_path
        / "outputs"
        / "docked"
        / run_id
        / "1AAA"
        / "HOLO"
        / "metal_site_audit.json"
    )
    audit.parent.mkdir(parents=True)
    audit.write_text(
        json.dumps(
            {
                "pdb_id": "1AAA",
                "variant": "HOLO",
                "ph_label": None,
                "source_files": {
                    "receptor_pdbqt": str(
                        recorded_root
                        / "outputs"
                        / "processed_pdbs"
                        / run_id
                        / "1AAA"
                        / "receptor"
                        / "1AAA.pdbqt"
                    )
                },
                "coordination_summary": {"retained_site_count": 1},
            }
        ),
        encoding="utf-8",
    )
    stale_processed_audit = (
        processed
        / "1AAA"
        / "work"
        / "1AAA_reduced.post_reduce.retained_hets_audit.json"
    )
    stale_processed_audit.parent.mkdir(parents=True)
    stale_processed_audit.write_text(
        json.dumps(
            {
                "pdb_id": "1AAA",
                "variant": "APO",
                "target_pdb": str(holo_receptor),
                "after": {"metals": 1, "cofactors": 0},
            }
        ),
        encoding="utf-8",
    )
    run_manifest.parent.mkdir(parents=True)
    run_manifest.write_text(
        yaml.safe_dump(
            {
                "run_id": run_id,
                "paths": {"all_dirs": str(recorded_root)},
                "proteins": {
                    "1AAA|HOLO|base": {
                        "pdb_id": "1AAA",
                        "variant": "HOLO",
                        "stages": {
                            "prep": {
                                "details": {
                                    "receptor_pdbqt": str(
                                        recorded_root
                                        / "outputs"
                                        / "processed_pdbs"
                                        / run_id
                                        / "1AAA"
                                        / "receptor"
                                        / "1AAA.pdbqt"
                                    )
                                }
                            }
                        },
                    },
                    "2BBB|HOLO|base": {
                        "pdb_id": "2BBB",
                        "variant": "HOLO",
                        "stages": {
                            "prep": {
                                "details": {
                                    "receptor_pdbqt": str(
                                        recorded_root
                                        / "outputs"
                                        / "processed_pdbs"
                                        / run_id
                                        / "2BBB"
                                        / "receptor"
                                        / "2BBB.pdbqt"
                                    )
                                }
                            }
                        },
                    },
                    "3CCC|HOLO|base": {
                        "pdb_id": "3CCC",
                        "variant": "HOLO",
                        "stages": {"prep": {"details": {}}},
                    },
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    release_manifest = tmp_path / "release.yaml"
    release_manifest.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "release_id": "receptor-evidence-fixture",
                "runs": [run_id],
                "scientific_policies": {
                    "receptor_selection": {"policy": "frozen"},
                    "native_redocking": {"policy": "required"},
                    "failure_handling": {"policy": "complete"},
                    "normalization": {"policy": "declared"},
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return release_manifest, run_manifest


def test_batch_evidence_preserves_conflicts_audits_and_missing_files(
    tmp_path: Path,
) -> None:
    release_manifest, _ = _release_fixture(tmp_path)

    records, summary = build_receptor_evidence_records(
        release_manifest,
        tmp_path,
        expected_context_count=3,
    )

    by_pdb = {row["pdb_id"]: row for row in records}
    assert by_pdb["1AAA"]["receptor_classification"] == "HOLO"
    assert by_pdb["1AAA"]["provenance"]["requested_observed_conflict"] is False
    audits = by_pdb["1AAA"]["provenance"]["structured_receptor_audits"]
    assert {item["audit_kind"] for item in audits} == {"metal_site_interaction"}
    assert audits[0]["binding_status"] == "exact_context_path_and_receptor"
    assert audits[0]["payload"]["coordination_summary"]["retained_site_count"] == 1
    assert len(audits[0]["source_sha256"]) == 64
    rejections = by_pdb["1AAA"]["provenance"][
        "structured_receptor_audit_rejections"
    ]
    assert [item["rejection_reason"] for item in rejections] == [
        "declared_variant_mismatch"
    ]
    assert by_pdb["2BBB"]["receptor_classification"] == "APO"
    assert by_pdb["2BBB"]["variant"] == "HOLO"
    assert by_pdb["2BBB"]["provenance"]["requested_observed_conflict"] is True
    assert by_pdb["3CCC"]["receptor_classification"] == ""
    assert by_pdb["3CCC"]["provenance"]["chemistry_evidence_status"] == "unresolved"
    assert summary["context_count"] == 3
    assert summary["chemistry_observed_count"] == 2
    assert summary["chemistry_unresolved_count"] == 1
    assert summary["requested_observed_conflict_count"] == 1
    assert summary["structured_receptor_audit_count"] == 1
    assert summary["structured_receptor_audit_rejection_count"] == 1
    assert summary["structured_receptor_audit_kind_counts"] == {
        "metal_site_interaction": 1
    }
    assert len(summary["context_inventory_sha256"]) == 64

    with pytest.raises(ReceptorEvidenceError, match="requires 90 contexts"):
        build_receptor_evidence_records(
            release_manifest,
            tmp_path,
            expected_context_count=90,
        )


def test_cleaned_and_docking_disagreement_stays_unresolved(tmp_path: Path) -> None:
    release_manifest, run_manifest_path = _release_fixture(tmp_path)
    processed = tmp_path / "outputs" / "processed_pdbs" / "run-1"
    cleaned = processed / "2BBB" / "receptor" / "2BBB_cleaned.pdb"
    cleaned.write_text(
        "".join(
            (
                _atom_line(1, "ZN", "ZN", 1, "ZN"),
                _atom_line(2, "CA", "ALA", 2, "C", het=False),
            )
        ),
        encoding="utf-8",
    )
    run_manifest = yaml.safe_load(run_manifest_path.read_text(encoding="utf-8"))
    run_manifest["proteins"]["2BBB|HOLO|base"]["stages"]["prep"]["details"][
        "cleaned_pdb"
    ] = str(cleaned)
    run_manifest_path.write_text(
        yaml.safe_dump(run_manifest, sort_keys=False),
        encoding="utf-8",
    )

    records, _ = build_receptor_evidence_records(release_manifest, tmp_path)
    record = next(row for row in records if row["pdb_id"] == "2BBB")

    assert record["receptor_classification"] == ""
    assert record["chemistry_evidence_status"] == "unresolved"
    assert (
        record["provenance"]["preparation_comparison_status"]
        == "classification_disagreement"
    )
    cleaned_evidence = record["provenance"]["cleaned_receptor_chemistry"]
    docking_evidence = record["provenance"]["exact_docking_receptor_chemistry"]
    assert cleaned_evidence["receptor_classification"] == "HOLO"
    assert docking_evidence["receptor_classification"] == "APO"
    assert (
        cleaned_evidence["prepared_receptor_sha256"]
        != docking_evidence["prepared_receptor_sha256"]
    )
    assert record["prepared_receptor_sha256"] == docking_evidence[
        "prepared_receptor_sha256"
    ]
    assert record["qualification_status"] == "pending_receptor_quality_review"


def test_production_audits_are_bound_to_exact_variant_and_ph(
    tmp_path: Path,
) -> None:
    release_manifest, run_manifest_path = _release_fixture(tmp_path)
    processed = tmp_path / "outputs" / "processed_pdbs" / "run-1"
    ph_receptor = (
        processed
        / "1AAA"
        / "HOLO"
        / "receptor"
        / "ph_ensemble"
        / "1AAA_pH7_0.pdbqt"
    )
    ph_receptor.parent.mkdir(parents=True)
    ph_receptor.write_text(
        "".join(
            (
                _atom_line(1, "ZN", "ZN", 1, "ZN"),
                _atom_line(2, "CA", "ALA", 2, "C", het=False),
            )
        ),
        encoding="utf-8",
    )
    ph_audit = (
        tmp_path
        / "outputs"
        / "docked"
        / "run-1"
        / "1AAA"
        / "HOLO"
        / "pH7_0"
        / "metal_site_audit.json"
    )
    ph_audit.parent.mkdir(parents=True)
    ph_audit.write_text(
        json.dumps(
            {
                "pdb_id": "1AAA",
                "variant": "HOLO",
                "ph_label": "pH7_0",
                "source_files": {"receptor_pdbqt": str(ph_receptor)},
                "coordination_summary": {"context": "pH7_0"},
            }
        ),
        encoding="utf-8",
    )
    unrelated_audit = (
        tmp_path
        / "outputs"
        / "docked"
        / "run-1"
        / "1AAA"
        / "APO"
        / "metal_site_audit.json"
    )
    unrelated_audit.parent.mkdir(parents=True)
    unrelated_audit.write_text(
        json.dumps(
            {
                "pdb_id": "1AAA",
                "variant": "APO",
                "ph_label": None,
                "source_files": {"receptor_pdbqt": str(ph_receptor)},
            }
        ),
        encoding="utf-8",
    )
    run_manifest = yaml.safe_load(run_manifest_path.read_text(encoding="utf-8"))
    run_manifest["proteins"]["1AAA|HOLO|pH7_0"] = {
        "pdb_id": "1AAA",
        "variant": "HOLO",
        "ph": "pH7_0",
        "stages": {"prep": {"details": {}}},
    }
    run_manifest_path.write_text(
        yaml.safe_dump(run_manifest, sort_keys=False),
        encoding="utf-8",
    )

    records, summary = build_receptor_evidence_records(
        release_manifest,
        tmp_path,
        expected_context_count=4,
    )
    one_aaa = [row for row in records if row["pdb_id"] == "1AAA"]
    by_ph = {row["ph_label"]: row for row in one_aaa}
    base_audits = by_ph[""]["provenance"]["structured_receptor_audits"]
    ph_audits = by_ph["pH7_0"]["provenance"]["structured_receptor_audits"]

    assert len(base_audits) == 1
    assert base_audits[0]["payload"]["ph_label"] is None
    assert len(ph_audits) == 1
    assert ph_audits[0]["payload"]["ph_label"] == "pH7_0"
    assert ph_audits[0]["binding_status"] == "exact_context_path_and_receptor"
    assert by_ph["pH7_0"]["provenance"][
        "exact_docking_receptor_resolution"
    ] == "canonical_output_path"
    assert summary["structured_receptor_audit_count"] == 2
    assert summary["structured_receptor_audit_rejection_count"] == 2


def test_exact_inventory_keys_and_hash_fail_closed(tmp_path: Path) -> None:
    release_manifest, _ = _release_fixture(tmp_path)
    _, first_summary = build_receptor_evidence_records(release_manifest, tmp_path)

    _, pinned_summary = build_receptor_evidence_records(
        release_manifest,
        tmp_path,
        expected_context_keys=first_summary["context_keys"],
        expected_inventory_sha256=first_summary["context_inventory_sha256"],
    )

    assert (
        pinned_summary["context_inventory_sha256"]
        == first_summary["context_inventory_sha256"]
    )
    wrong_keys = [*first_summary["context_keys"]]
    wrong_keys[-1] = "run-1|9ZZZ|HOLO|base"
    with pytest.raises(ReceptorEvidenceError, match="inventory key mismatch"):
        build_receptor_evidence_records(
            release_manifest,
            tmp_path,
            expected_context_keys=wrong_keys,
        )
    with pytest.raises(ReceptorEvidenceError, match="SHA-256 mismatch"):
        build_receptor_evidence_records(
            release_manifest,
            tmp_path,
            expected_inventory_sha256="0" * 64,
        )


def test_writes_hashed_annotation_and_refuses_implicit_overwrite(tmp_path: Path) -> None:
    release_manifest, _ = _release_fixture(tmp_path)
    output = tmp_path / "receptors.yaml"

    summary = write_receptor_evidence(
        release_manifest,
        tmp_path,
        output,
        expected_context_count=3,
    )

    assert output.is_file()
    assert output.with_suffix(".summary.json").is_file()
    assert summary["annotation_sha256"]
    payload = yaml.safe_load(output.read_text(encoding="utf-8"))
    assert len(payload["records"]) == 3
    with pytest.raises(ReceptorEvidenceError, match="output already exists"):
        write_receptor_evidence(
            release_manifest,
            tmp_path,
            output,
            expected_context_count=3,
        )
