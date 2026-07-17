from __future__ import annotations

import csv
import hashlib
import json
import sqlite3
from copy import deepcopy
from pathlib import Path

import pytest

from analysis.atlas_database.exports import _selected_ligand_stereo_projection
from analysis.atlas_database.ligand_stereo_evidence import (
    LigandStereoEvidenceError,
    ingest_ligand_stereo_evidence,
)
from analysis.atlas_database.schema import create_schema
from analysis.cli.build_atlas_database import _sanitize_public_database
from prep_ligands.ligand_stereo_audit import (
    LigandStereoAuditConfig,
    LigandStereoSource,
    audit_ligand_stereochemistry,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _audit_rows(
    tmp_path: Path,
    rows: list[dict[str, str]],
    *,
    stem: str,
) -> tuple[Path, Path, list[dict[str, str]]]:
    source = tmp_path / f"{stem}.csv"
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    result = audit_ligand_stereochemistry(
        LigandStereoAuditConfig(
            library_id="FDA",
            output_dir=tmp_path / f"{stem}-audit",
            sources=(
                LigandStereoSource(
                    path=source,
                    kind="csv",
                    id_field="ligand_id",
                    structure_field="smiles",
                    structure_format="smiles",
                ),
            ),
        )
    )
    with result.records_csv.open(encoding="utf-8", newline="") as handle:
        audited = list(csv.DictReader(handle))
    return source, result.records_csv, audited


def _database(
    database: str | Path,
    ligand_ids: list[str],
) -> sqlite3.Connection:
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    create_schema(connection)
    connection.execute(
        """INSERT INTO releases(
        release_id, schema_version, manifest_sha256, manifest_json, imported_at
        ) VALUES ('release', '1', ?, '{}', '2026-07-17T00:00:00Z')""",
        ("a" * 64,),
    )
    connection.execute(
        """INSERT INTO runs(
        run_id, release_id, manifest_path, manifest_sha256, manifest_json
        ) VALUES ('run', 'release', 'manifest.yaml', ?, '{}')""",
        ("b" * 64,),
    )
    connection.execute(
        """INSERT INTO receptor_contexts(
        receptor_context_id, run_id, pdb_id, variant, ph_label, library
        ) VALUES (1, 'run', '1ABC', 'HOLO', '7.4', 'FDA')"""
    )
    for ligand_id, canonical_id in enumerate(ligand_ids, start=1):
        connection.execute(
            "INSERT INTO ligands(ligand_id, canonical_id) VALUES (?, ?)",
            (ligand_id, canonical_id),
        )
        connection.execute(
            """INSERT INTO pair_cells(
            receptor_context_id, ligand_id, final_status
            ) VALUES (1, ?, 'valid')""",
            (ligand_id,),
        )
    return connection


def test_csv_audit_reports_stereo_without_assigning_missing_configuration(
    tmp_path: Path,
) -> None:
    source, _annotation, audited = _audit_rows(
        tmp_path,
        [
            {"ligand_id": "specified", "smiles": "C[C@H](O)F"},
            {"ligand_id": "unresolved", "smiles": "CC(O)F"},
            {"ligand_id": "achiral", "smiles": "CCO"},
        ],
        stem="source",
    )
    by_id = {row["ligand_canonical_id"]: row for row in audited}

    assert by_id["specified"]["audit_status"] == "fully_specified"
    assert by_id["specified"]["specified_stereo_count"] == "1"

    unresolved = by_id["unresolved"]
    assert unresolved["audit_status"] == "all_unspecified"
    assert unresolved["unspecified_stereo_count"] == "1"
    assert unresolved["has_unresolved_potential_stereo"] == "true"
    assert "@" not in unresolved["canonical_isomeric_smiles"]

    assert by_id["achiral"]["audit_status"] == "no_potential_stereo"
    assert by_id["achiral"]["potential_stereo_count"] == "0"
    assert {row["source_file_sha256"] for row in audited} == {_sha256(source)}


def test_ingestion_accepts_zero_count_and_rejects_tampered_record_hash(
    tmp_path: Path,
) -> None:
    _source, annotation, audited = _audit_rows(
        tmp_path,
        [
            {"ligand_id": "achiral", "smiles": "CCO"},
            {"ligand_id": "unresolved", "smiles": "CC(O)F"},
        ],
        stem="roundtrip",
    )

    connection = _database(":memory:", ["achiral", "unresolved"])
    try:
        assert (
            ingest_ligand_stereo_evidence(
                connection,
                audited,
                annotation,
                _sha256(annotation),
            )
            == 2
        )
        stored = connection.execute(
            """SELECT audit_status, selected_for_release, selection_method
            FROM ligand_stereo_evidence
            WHERE source_record_id='achiral'"""
        ).fetchone()
        assert dict(stored) == {
            "audit_status": "no_potential_stereo",
            "selected_for_release": 1,
            "selection_method": "only_evidence_record",
        }
    finally:
        connection.close()

    tampered = deepcopy(audited)
    tampered[0]["source_record_sha256"] = "0" * 64
    connection = _database(":memory:", ["achiral", "unresolved"])
    try:
        with pytest.raises(
            LigandStereoEvidenceError,
            match="source_record_sha256 does not match source record",
        ):
            ingest_ligand_stereo_evidence(
                connection,
                tampered,
                annotation,
                _sha256(annotation),
            )
    finally:
        connection.close()


def test_duplicate_evidence_requires_exactly_one_explicit_selection(
    tmp_path: Path,
) -> None:
    _source_a, _annotation_a, first = _audit_rows(
        tmp_path,
        [{"ligand_id": "drug", "smiles": "CC(O)F"}],
        stem="candidate-a",
    )
    _source_b, _annotation_b, second = _audit_rows(
        tmp_path,
        [{"ligand_id": "drug", "smiles": "C[C@H](O)F"}],
        stem="candidate-b",
    )
    combined = first + second
    annotation = tmp_path / "combined.json"
    annotation.write_text(json.dumps(combined), encoding="utf-8")

    connection = _database(":memory:", ["drug"])
    try:
        with pytest.raises(
            LigandStereoEvidenceError,
            match="requires an explicit selected_for_release boolean",
        ):
            ingest_ligand_stereo_evidence(
                connection,
                combined,
                annotation,
                _sha256(annotation),
            )

        combined[0]["selected_for_release"] = "true"
        combined[1]["selected_for_release"] = "false"
        annotation.write_text(json.dumps(combined), encoding="utf-8")

        assert (
            ingest_ligand_stereo_evidence(
                connection,
                combined,
                annotation,
                _sha256(annotation),
            )
            == 2
        )
        selected = connection.execute(
            """SELECT source_record_sha256
            FROM ligand_stereo_evidence
            WHERE selected_for_release=1"""
        ).fetchall()
        assert [row[0] for row in selected] == [
            combined[0]["source_record_sha256"]
        ]
    finally:
        connection.close()


def test_public_projection_redacts_raw_evidence_and_legacy_db_is_empty(
    tmp_path: Path,
) -> None:
    _source, annotation, audited = _audit_rows(
        tmp_path,
        [{"ligand_id": "drug", "smiles": "CC(O)F"}],
        stem="public",
    )
    database = tmp_path / "public.sqlite"
    connection = _database(database, ["drug"])
    try:
        ingest_ligand_stereo_evidence(
            connection,
            audited,
            annotation,
            _sha256(annotation),
        )
        connection.commit()
    finally:
        connection.close()

    _sanitize_public_database(database, tmp_path)

    with sqlite3.connect(database) as connection:
        connection.row_factory = sqlite3.Row
        stored = connection.execute(
            """SELECT ligand_source_path, annotation_source_path,
                      source_record_evidence_json, source_record_json
            FROM ligand_stereo_evidence"""
        ).fetchone()
        assert stored is not None
        assert not Path(stored["ligand_source_path"]).is_absolute()
        assert not Path(stored["annotation_source_path"]).is_absolute()
        assert json.loads(stored["source_record_evidence_json"])[
            "public_projection"
        ] == "omitted"
        assert json.loads(stored["source_record_json"])["public_projection"] == "omitted"

        projection = _selected_ligand_stereo_projection(connection)
        assert projection[1]["selected_stereo_evidence_count"] == 1
        assert projection[1]["selected_stereo_evidence"][0][
            "source_record_evidence"
        ]["public_projection"] == "omitted"

    with sqlite3.connect(":memory:") as legacy:
        legacy.execute(
            "CREATE TABLE ligands(ligand_id INTEGER PRIMARY KEY, canonical_id TEXT)"
        )
        assert _selected_ligand_stereo_projection(legacy) == {}
