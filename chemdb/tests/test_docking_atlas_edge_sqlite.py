from __future__ import annotations

import hashlib
import json
import sqlite3
from pathlib import Path

import pytest

from analysis.atlas_database.receptor_evidence import CHEMISTRY_CLASSIFICATION_METHOD
from analysis.atlas_database.schema import create_schema
from analysis.reporting.docking_atlas_delivery import build_deployment_preflight
from analysis.reporting.docking_atlas_edge import main as edge_main
from analysis.reporting.docking_atlas_edge_sqlite import build_streaming_edge_bundle


def _release_database(
    site_dir: Path,
    *,
    target_count: int = 60,
    drug_count: int = 100,
    legacy_score_source: bool = False,
    apo_target_count: int = 0,
) -> Path:
    downloads = site_dir / "downloads"
    downloads.mkdir(parents=True)
    database = downloads / "docking_atlas.sqlite"
    with sqlite3.connect(database) as connection:
        create_schema(connection)
        connection.execute(
            """INSERT INTO releases(
            release_id, schema_version, title, manifest_sha256, manifest_json,
            imported_at) VALUES (?, ?, ?, ?, ?, ?)""",
            (
                "atlas-stream-test",
                "0.1",
                "Streaming Atlas test",
                "a" * 64,
                json.dumps({"summary": "Failure-complete streaming test"}),
                "2026-07-13T00:00:00Z",
            ),
        )
        connection.execute(
            """INSERT INTO runs(
            run_id, release_id, status, manifest_path, manifest_sha256,
            manifest_json) VALUES (?, ?, ?, ?, ?, ?)""",
            (
                "stream-run",
                "atlas-stream-test",
                "completed",
                "manifests/stream-run.yaml",
                "b" * 64,
                "{}",
            ),
        )
        contexts = [
            (
                index + 1,
                "stream-run",
                f"P{index:04d}",
                "HOLO",
                "pH7.4",
                "completed",
                "qualified",
            )
            for index in range(target_count)
        ]
        connection.executemany(
            """INSERT INTO receptor_contexts(
            receptor_context_id, run_id, pdb_id, variant, ph_label, status,
            native_redock_status) VALUES (?, ?, ?, ?, ?, ?, ?)""",
            contexts,
        )
        connection.executemany(
            """INSERT INTO receptor_annotations(
            receptor_annotation_id, receptor_context_id, receptor_classification,
            qualification_status, native_redock_status, classification_method,
            chemistry_evidence_status, prepared_receptor_sha256,
            canonical_chemistry_policy_sha256, source_path, source_sha256,
            source_record_index, source_record_json)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            [
                (
                    index + 1,
                    index + 1,
                    "APO" if index < apo_target_count else "HOLO",
                    "qualified",
                    "qualified",
                    CHEMISTRY_CLASSIFICATION_METHOD,
                    "observed",
                    "d" * 64,
                    "e" * 64,
                    "annotations/receptors.csv",
                    "c" * 64,
                    index,
                    "{}",
                )
                for index in range(target_count)
            ],
        )
        connection.executemany(
            "INSERT INTO ligands(ligand_id, canonical_id, display_name) VALUES (?, ?, ?)",
            [
                (index + 1, f"DRUG-{index:04d}", f"Drug {index:04d}")
                for index in range(drug_count)
            ],
        )
        pairs: list[tuple[object, ...]] = []
        pair_id = 0
        for target_index in range(target_count):
            for drug_index in range(drug_count):
                pair_id += 1
                result_json = json.dumps({"final_score_source": "consensus_vs_decoy_z"})
                pairs.append(
                    (
                        pair_id,
                        target_index + 1,
                        drug_index + 1,
                        "valid",
                        1,
                        1,
                        1,
                        "selected_final_score_pose",
                        float(drug_index),
                        "consensus_vs_decoy_z",
                        result_json,
                    )
                )
        connection.executemany(
            """INSERT INTO pair_cells(
            pair_cell_id, receptor_context_id, ligand_id, final_status,
            expected, has_result, pose_valid, pose_validation_scope, final_score,
            final_score_source, result_json)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            pairs,
        )
        connection.executemany(
            """INSERT INTO artifacts(
            artifact_id, run_id, receptor_context_id, pair_cell_id, ligand_id,
            artifact_role, artifact_scope, member_name, sha256, size_bytes,
            file_type, verified, artifact_json)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            [
                (
                    1,
                    "stream-run",
                    1,
                    1,
                    1,
                    "docking_pose",
                    "pair",
                    "pose-1.sdf",
                    "d" * 64,
                    123,
                    "chemical/x-mdl-sdfile",
                    1,
                    '{"selected_for_release":true}',
                ),
                (
                    2,
                    "stream-run",
                    1,
                    1,
                    1,
                    "docking_pose",
                    "pair",
                    "unselected-pose.sdf",
                    "e" * 64,
                    122,
                    "chemical/x-mdl-sdfile",
                    1,
                    '{"selected_for_release":false}',
                ),
                (
                    3,
                    "stream-run",
                    1,
                    1,
                    1,
                    "raw_stdout",
                    "pair",
                    "stdout.txt",
                    "f" * 64,
                    121,
                    "text/plain",
                    1,
                    '{"selected_for_release":true}',
                ),
                (
                    4,
                    "stream-run",
                    1,
                    1,
                    1,
                    "pose_image",
                    "pair",
                    "unverified.png",
                    "a" * 64,
                    120,
                    "image/png",
                    0,
                    "{}",
                ),
                (
                    5,
                    "stream-run",
                    1,
                    1,
                    1,
                    "pose_validation_report",
                    "pair",
                    "invalid-hash.json",
                    "invalid-sha256",
                    119,
                    "application/json",
                    1,
                    "{}",
                ),
            ],
        )
        if legacy_score_source:
            connection.executescript(
                """
                ALTER TABLE pair_cells DROP COLUMN final_score_source;
                ALTER TABLE pair_cells DROP COLUMN pose_validation_method;
                ALTER TABLE pair_cells DROP COLUMN pose_validation_scope;
                ALTER TABLE pair_cells DROP COLUMN pose_validation_thresholds_json;
                ALTER TABLE receptor_annotations DROP COLUMN classification_method;
                ALTER TABLE receptor_annotations DROP COLUMN chemistry_evidence_status;
                ALTER TABLE receptor_annotations DROP COLUMN prepared_receptor_sha256;
                ALTER TABLE receptor_annotations
                    DROP COLUMN canonical_chemistry_policy_sha256;
                ALTER TABLE receptor_annotations DROP COLUMN retained_metal_atom_count;
                ALTER TABLE receptor_annotations
                    DROP COLUMN retained_cofactor_residue_count;
                ALTER TABLE receptor_annotations DROP COLUMN requested_observed_conflict;
                ALTER TABLE receptor_annotations DROP COLUMN native_redock_reason;
                PRAGMA user_version=2;
                """
            )
        connection.commit()
    digest = hashlib.sha256(database.read_bytes()).hexdigest()
    (site_dir / "site_manifest.json").write_text(
        json.dumps(
            {
                "release_id": "atlas-stream-test",
                "downloads": [
                    {
                        "label": "SQLite release database",
                        "url": "downloads/docking_atlas.sqlite",
                        "content_hash": f"sha256:{digest}",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return database


def _read_object(output_dir: Path, key: str) -> dict[str, object]:
    return json.loads((output_dir / "objects" / key).read_text(encoding="utf-8"))


def test_sqlite_streaming_export_uses_coarse_shards_not_pair_objects(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    database = _release_database(site_dir)
    output_dir = tmp_path / "edge"

    summary = build_streaming_edge_bundle(
        database,
        output_dir,
        source_site_dir=site_dir,
        download_base_url="https://downloads.example.org",
        batch_rows=200,
        coarse_shard_rows=1_000,
        max_pairs=7_000,
    )

    assert summary["input_mode"] == "sqlite_streaming"
    assert summary["publication_scale_supported"] is True
    assert summary["counts"] == {
        "targets": 60,
        "drugs": 100,
        "pairs": 6_000,
        "pair_detail_shards": 6,
    }
    assert summary["bounded_memory"] == {
        "source_fetch_batch_rows": 200,
        "coarse_pair_shard_rows": 1_000,
        "max_buffered_pair_rows": 1_000,
        "work_database": "on_disk_disposable",
        "work_database_peak_bytes": summary["bounded_memory"][
            "work_database_peak_bytes"
        ],
        "pair_collection_materialized": False,
    }
    assert summary["object_shape"]["pair_object_per_cell"] is False
    assert summary["object_shape"]["pair_objects_emitted"] == 0
    assert summary["r2_object_count"] < 300
    assert not (output_dir / ".atlas-edge-stream-work.sqlite").exists()

    object_manifest = json.loads(
        (output_dir / "object_manifest.json").read_text(encoding="utf-8")
    )
    keys = [str(row["key"]) for row in object_manifest["objects"]]
    assert not any("/records/pairs/" in key for key in keys)
    shard_keys = [key for key in keys if "/records/pair-shards/" in key]
    assert len(shard_keys) == 6
    assert all(
        len(_read_object(output_dir, key)["records"]) <= 1_000 for key in shard_keys
    )

    release_token = summary["release_token"]
    target_index = _read_object(
        output_dir, f"releases/{release_token}/indexes/targets.json"
    )
    first_target = target_index["records"][0]
    target_record = _read_object(
        output_dir,
        f"releases/{release_token}/records/targets/{first_target['route_id']}.json",
    )
    assert len(target_record["pairs"]) == 100
    assert target_record["pairs"][0]["rank"] is None
    assert target_record["pairs"][0]["rank_track_label"] == "Unqualified"
    assert target_record["pairs"][0]["pair_shard_id"]

    pair_index = _read_object(
        output_dir, f"releases/{release_token}/indexes/pairs.json"
    )
    assert pair_index["kind"] == "coarse_pair_shards"
    assert sum(int(row["count"]) for row in pair_index["shards"]) == 6_000
    first_shard = _read_object(output_dir, shard_keys[0])
    first_pair_record = first_shard["records"][0]
    public_artifacts = first_pair_record["artifacts"]
    assert first_pair_record["pair"]["artifact_count"] == 1
    assert first_pair_record["pair"]["verified_artifact_count"] == 1
    assert len(public_artifacts) == 1
    selected_artifact = public_artifacts[0]
    assert selected_artifact["artifact_role"] == "docking_pose"
    assert selected_artifact["publication_policy"] == "public_if_selected_for_release"

    report = build_deployment_preflight(
        output_dir, site_dir=site_dir, provider="generic"
    )
    assert report["status"] == "ready"
    assert report["edge_objects"]["count"] == summary["r2_object_count"]


def test_sqlite_streaming_export_never_invents_ranks_from_unbound_scores(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "mixed-site"
    database = _release_database(
        site_dir,
        target_count=2,
        drug_count=3,
        apo_target_count=1,
    )
    output_dir = tmp_path / "mixed-edge"

    summary = build_streaming_edge_bundle(
        database,
        output_dir,
        batch_rows=100,
        coarse_shard_rows=100,
    )

    assert summary["coverage"]["rank_eligible_count"] == 0
    assert summary["coverage"]["apo_exploratory_rank_eligible_count"] == 0
    release_token = summary["release_token"]
    target_index = _read_object(
        output_dir, f"releases/{release_token}/indexes/targets.json"
    )
    target = _read_object(
        output_dir,
        f"releases/{release_token}/records/targets/"
        f"{target_index['records'][0]['route_id']}.json",
    )
    assert target["entity"]["receptor_classification"] == "APO"
    assert target["pairs"][0]["rank"] is None
    assert target["pairs"][0]["rank_track_label"] == "Unqualified"
    assert target["pairs"][0]["ranking_track"] is None
    assert target["pairs"][0]["ranking_eligibility_reason"] == "apo_receptor"
    assert target["pairs"][0]["rank_within_receptor"] is None
    assert target["pairs"][0]["apo_rank_within_receptor"] is None

    drug_index = _read_object(
        output_dir, f"releases/{release_token}/indexes/drugs.json"
    )
    drug = _read_object(
        output_dir,
        f"releases/{release_token}/records/drugs/"
        f"{drug_index['records'][0]['route_id']}.json",
    )
    assert {row["rank_track_label"] for row in drug["pairs"]} == {"Unqualified"}
    assert {row["rank"] for row in drug["pairs"]} == {None}
    assert {row["ranking_track"] for row in drug["pairs"]} == {None}
    assert all(row["rank_across_receptors"] is None for row in drug["pairs"])
    assert all(row["apo_rank_across_receptors"] is None for row in drug["pairs"])


def test_sqlite_streaming_export_supports_legacy_schema_and_score_source(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "legacy-site"
    database = _release_database(
        site_dir,
        target_count=2,
        drug_count=3,
        legacy_score_source=True,
    )
    output_dir = tmp_path / "legacy-edge"

    summary = build_streaming_edge_bundle(
        database,
        output_dir,
        batch_rows=100,
        coarse_shard_rows=100,
    )

    assert summary["source_sqlite_user_version"] == 2
    release_token = summary["release_token"]
    pair_index = _read_object(
        output_dir, f"releases/{release_token}/indexes/pairs.json"
    )
    shard = _read_object(
        output_dir,
        f"releases/{release_token}/records/pair-shards/{pair_index['shards'][0]['id']}.json",
    )
    assert {row["pair"]["final_score_source"] for row in shard["records"]} == {
        "consensus_vs_decoy_z"
    }
    assert {row["pair"]["pose_validation_scope"] for row in shard["records"]} == {None}
    assert {row["pair"]["rank_eligible"] for row in shard["records"]} == {0}
    assert {row["pair"]["ranking_eligibility_reason"] for row in shard["records"]} == {
        "receptor_chemistry_method_unqualified"
    }
    target_index = _read_object(
        output_dir, f"releases/{release_token}/indexes/targets.json"
    )
    target = _read_object(
        output_dir,
        f"releases/{release_token}/records/targets/"
        f"{target_index['records'][0]['route_id']}.json",
    )
    assert target["entity"]["receptor_classification_method"] is None


def test_streaming_cli_default_output_is_sibling_of_site(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    site_dir = tmp_path / "site"
    database = _release_database(site_dir, target_count=2, drug_count=3)

    with pytest.raises(ValueError, match="must be disjoint"):
        build_streaming_edge_bundle(
            database,
            site_dir / "edge",
            batch_rows=100,
            coarse_shard_rows=100,
        )
    assert not (site_dir / "edge").exists()

    assert (
        edge_main(
            [
                "--database",
                str(database),
                "--batch-rows",
                "100",
                "--coarse-shard-rows",
                "100",
            ]
        )
        == 0
    )
    summary = json.loads(capsys.readouterr().out)
    assert summary["input_mode"] == "sqlite_streaming"
    assert (tmp_path / "edge" / "edge_bundle_summary.json").is_file()
    assert not (site_dir / "edge").exists()


def test_sqlite_streaming_export_fails_before_output_above_max_pairs(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    database = _release_database(site_dir, target_count=2, drug_count=3)
    output_dir = tmp_path / "edge"

    with pytest.raises(ValueError, match="above --max-pairs"):
        build_streaming_edge_bundle(database, output_dir, max_pairs=5)

    assert not output_dir.exists()


def test_streaming_worker_and_spa_resolve_coarse_pair_shards(tmp_path: Path) -> None:
    site_dir = tmp_path / "site"
    database = _release_database(site_dir, target_count=2, drug_count=3)
    output_dir = tmp_path / "edge"
    build_streaming_edge_bundle(
        database,
        output_dir,
        batch_rows=100,
        coarse_shard_rows=100,
    )

    worker = (output_dir / "src" / "index.mjs").read_text(encoding="utf-8")
    app = (output_dir / "public" / "assets" / "app.js").read_text(encoding="utf-8")
    assert 'parts[3] === "pair-shards"' in worker
    assert "records/pair-shards" in worker
    assert 'config.pair_storage !== "coarse_shards"' in app
    assert "pair_shard_id" in app
    assert "pair is absent from its declared shard" in app
