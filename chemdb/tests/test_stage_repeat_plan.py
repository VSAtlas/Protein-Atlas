from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

import pytest

from analysis.atlas_database.score_source_contract import (
    LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z,
)
from analysis.atlas_database.schema import create_schema
from cli import stage_repeat_plan
from cli.qol import stages
from cli.qol.dispatch import dispatch
from cli.stage_repeat_plan import (
    PlanRequest,
    StageRepeatPlanError,
    build_stage_repeat_plan,
)


def _write_file(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _write_rows(path: Path, rows: list[dict[str, object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _ledger_rows(result: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for raw_path in result["ledger_paths"]:  # type: ignore[union-attr]
        ledger_path = Path(raw_path)
        rows.extend(
            json.loads(line)
            for line in ledger_path.read_text(encoding="utf-8").splitlines()
        )
    return rows


def _csv_pair_row(
    root: Path,
    ligand: str,
    *,
    variant: str = "HOLO",
    ph_label: str = "base",
) -> dict[str, object]:
    receptor = _write_file(root / "artifacts" / "receptor.pdbqt", "RECEPTOR\n")
    ligand_path = _write_file(root / "artifacts" / f"{ligand}.pdbqt", f"{ligand}\n")
    return {
        "run_id": "run_a",
        "pdb_id": "1ABC",
        "variant": variant,
        "ph_label": ph_label,
        "ligand_base": ligand,
        "prepared_receptor_path": str(receptor),
        "prepared_ligand_path": str(ligand_path),
        "center_x": 1.0,
        "center_y": 2.0,
        "center_z": 3.0,
        "box_x": 20.0,
        "box_y": 21.0,
        "box_z": 22.0,
    }


def _request(**overrides: object) -> PlanRequest:
    values: dict[str, object] = {
        "run_id": "run_a",
        "stage": "vina",
        "receptors": (),
        "all_receptors": True,
        "ligands": (),
        "all_ligands": True,
        "all_contexts": True,
        "fraction": 1.0,
        "count": None,
        "seed": 17,
        "cpus": 8,
        "source_stage": "stage1",
    }
    values.update(overrides)
    return PlanRequest(**values)  # type: ignore[arg-type]


def test_hash_selection_is_deterministic_and_persisted(tmp_path: Path) -> None:
    rows = [_csv_pair_row(tmp_path, f"lig_{index}") for index in range(6)]
    source = _write_rows(tmp_path / "master_rows.csv", rows)
    request = _request(fraction=0.5, count=None, seed=991, selection_strategy="hash")

    first = build_stage_repeat_plan(
        tmp_path,
        request,
        master_rows=source,
        output=tmp_path / "plans" / "one.json",
    )
    second = build_stage_repeat_plan(
        tmp_path,
        request,
        master_rows=source,
        output=tmp_path / "plans" / "two.json",
    )
    stable_paths = [
        Path(first["plan_path"]),
        Path(first["pairs_path"]),
        *(Path(path) for path in first["ledger_paths"]),
    ]
    stable_snapshots = {
        path: (path.read_bytes(), path.stat().st_ino, path.stat().st_mtime_ns)
        for path in stable_paths
    }
    reused = build_stage_repeat_plan(
        tmp_path,
        request,
        master_rows=source,
        output=tmp_path / "plans" / "one.json",
    )

    first_plan = first["plan"]
    second_plan = second["plan"]
    assert first["reused_existing"] is False
    assert second["reused_existing"] is False
    assert reused["reused_existing"] is True
    assert first_plan["plan_id"] == second_plan["plan_id"]
    assert reused["plan"] == first_plan
    for path, (payload, inode, modified_ns) in stable_snapshots.items():
        assert path.read_bytes() == payload
        assert path.stat().st_ino == inode
        assert path.stat().st_mtime_ns == modified_ns
    assert first_plan["summary"] == {"selected": 3, "ready": 3, "blocked": 0}
    first_rows = _ledger_rows(first)
    second_rows = _ledger_rows(second)
    assert [row["ligand_id"] for row in first_rows] == [
        row["ligand_id"] for row in second_rows
    ]
    assert "rows" not in first_plan
    assert first_plan["selection"]["algorithm"] == "sha256_seeded_pair_order_v1"
    assert first_plan["selection"]["strategy"] == "hash"
    assert first_plan["selection"]["seed"] == 991
    assert first_plan["scalability"]["json_embeds_pair_rows"] is False
    assert first_plan["execution"]["executor_status"] == "planner_only"
    assert first_plan["execution"]["can_execute"] is False
    assert first_plan["execution"]["gap"] == (
        "exact_pair_plan_executor_not_implemented"
    )
    assert Path(first["pairs_path"]).is_file()
    assert all(Path(path).is_file() for path in first["ledger_paths"])


def test_plan_refuses_to_overwrite_conflicting_existing_output(
    tmp_path: Path,
) -> None:
    source = _write_rows(
        tmp_path / "master_rows.csv", [_csv_pair_row(tmp_path, "lig_a")]
    )
    output = tmp_path / "plan.json"
    result = build_stage_repeat_plan(
        tmp_path,
        _request(),
        master_rows=source,
        output=output,
    )
    pairs_path = Path(result["pairs_path"])
    conflicting_payload = pairs_path.read_bytes() + b"conflict\n"
    pairs_path.write_bytes(conflicting_payload)

    with pytest.raises(
        StageRepeatPlanError,
        match="conflicting stage-repeat output already exists",
    ):
        build_stage_repeat_plan(
            tmp_path,
            _request(),
            master_rows=source,
            output=output,
        )

    assert pairs_path.read_bytes() == conflicting_payload


def test_plan_refuses_incomplete_existing_output_set(tmp_path: Path) -> None:
    source = _write_rows(
        tmp_path / "master_rows.csv", [_csv_pair_row(tmp_path, "lig_a")]
    )
    output = tmp_path / "plan.json"
    result = build_stage_repeat_plan(
        tmp_path,
        _request(),
        master_rows=source,
        output=output,
    )
    missing_ledger = Path(result["ledger_paths"][0])
    missing_ledger.unlink()
    plan_payload = output.read_bytes()

    with pytest.raises(
        StageRepeatPlanError,
        match="incomplete stage-repeat output set already exists",
    ):
        build_stage_repeat_plan(
            tmp_path,
            _request(),
            master_rows=source,
            output=output,
        )

    assert output.read_bytes() == plan_payload
    assert not missing_ledger.exists()


def test_plan_rejects_oversubscription_and_ambiguous_contexts(tmp_path: Path) -> None:
    rows = [
        _csv_pair_row(tmp_path, "lig_a", variant="APO"),
        _csv_pair_row(tmp_path, "lig_b", variant="HOLO"),
    ]
    source = _write_rows(tmp_path / "master_rows.csv", rows)

    with pytest.raises(StageRepeatPlanError, match="cpus must be between 1 and 32"):
        build_stage_repeat_plan(
            tmp_path,
            _request(cpus=33),
            master_rows=source,
            output=tmp_path / "bad-cpu.json",
        )

    with pytest.raises(StageRepeatPlanError, match="multiple receptor contexts"):
        build_stage_repeat_plan(
            tmp_path,
            _request(all_contexts=False),
            master_rows=source,
            output=tmp_path / "ambiguous.json",
        )


def test_pose_validation_rejects_unknown_score_pose_semantics(tmp_path: Path) -> None:
    row = _csv_pair_row(tmp_path, "lig_a")
    generic_pose = _write_file(tmp_path / "artifacts" / "generic_pose.pdbqt", "POSE\n")
    row["pose_path"] = str(generic_pose)
    source = _write_rows(tmp_path / "master_rows.csv", [row])
    request = _request(stage="pose-validation", source_stage="")

    blocked_result = build_stage_repeat_plan(
        tmp_path,
        request,
        master_rows=source,
        output=tmp_path / "blocked.json",
    )
    blocked = blocked_result["plan"]
    assert blocked["summary"] == {"selected": 1, "ready": 0, "blocked": 1}
    blocked_row = _ledger_rows(blocked_result)[0]
    assert "selected_docking_pose_missing_or_unverified" in blocked_row["blockers"]

    row["selected_pose_path"] = str(generic_pose)
    row["selected_pose_stage"] = "stage3"
    source = _write_rows(tmp_path / "master_rows_explicit.csv", [row])
    ready_result = build_stage_repeat_plan(
        tmp_path,
        request,
        master_rows=source,
        output=tmp_path / "ready.json",
    )
    ready = ready_result["plan"]
    ready_row = _ledger_rows(ready_result)[0]
    assert ready["summary"] == {"selected": 1, "ready": 0, "blocked": 1}
    assert "score_source_pose_linkage_unknown" in ready_row["blockers"]
    assert "selected_docking_pose" not in ready_row["artifacts"]


def _insert_sqlite_fixture(
    path: Path,
    *,
    linked_pose: bool,
    score_source: str = "vina_pose_score",
    result_json: dict[str, object] | None = None,
    pose_metadata: dict[str, object] | None = None,
) -> None:
    digest_a = "a" * 64
    digest_b = "b" * 64
    digest_c = "c" * 64
    digest_d = "d" * 64
    connection = sqlite3.connect(path)
    create_schema(connection)
    connection.execute(
        "INSERT INTO releases VALUES (?, ?, ?, ?, ?, ?, ?)",
        ("release", "5", "test", "now", digest_a, "{}", "now"),
    )
    connection.execute(
        """INSERT INTO runs
        (run_id, release_id, status, manifest_path, manifest_sha256,
         command_json, git_json, paths_json, timing_json, resources_json, manifest_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (
            "run_a",
            "release",
            "complete",
            "manifest",
            digest_a,
            "{}",
            "{}",
            "{}",
            "{}",
            "{}",
            "{}",
        ),
    )
    context_id = connection.execute(
        """INSERT INTO receptor_contexts
        (run_id, pdb_id, variant, ph_label, center_json, box_json)
        VALUES (?, ?, ?, ?, ?, ?)""",
        ("run_a", "1ABC", "HOLO", "", "[1,2,3]", "[20,20,20]"),
    ).lastrowid
    ligand_id = connection.execute(
        "INSERT INTO ligands (canonical_id) VALUES (?)", ("lig_a",)
    ).lastrowid
    pair_id = connection.execute(
        """INSERT INTO pair_cells
        (receptor_context_id, ligand_id, final_status, expected, has_result)
        VALUES (?, ?, ?, 1, 1)""",
        (context_id, ligand_id, "calculated_unvalidated"),
    ).lastrowid
    completion_id = connection.execute(
        """INSERT INTO completion_records
        (run_id, receptor_context_id, engine, stage, chunk_id,
         completion_path, completion_sha256, completion_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
        (
            "run_a",
            context_id,
            "vina",
            "stage3",
            "chunk",
            "completion.json",
            digest_b,
            "{}",
        ),
    ).lastrowid
    connection.execute(
        """INSERT INTO docking_attempts
        (pair_cell_id, completion_record_id, engine, stage, chunk_id,
         status, selected_for_release, selection_method)
        VALUES (?, ?, ?, ?, ?, ?, 1, ?)""",
        (pair_id, completion_id, "vina", "stage3", "chunk", "completed", "unique"),
    )
    result_id = connection.execute(
        """INSERT INTO result_attempts
        (pair_cell_id, completion_record_id, completion_link_method,
         completion_link_evidence_json, input_csv_path, input_csv_sha256,
         source_row_number, result_sha256, final_status, final_score_source,
         result_json, selected_for_release, selection_method)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1, ?)""",
        (
            pair_id,
            completion_id,
            "manifest_exact",
            json.dumps({"completion_sha256": digest_b}),
            "master.csv",
            digest_c,
            2,
            digest_d,
            "calculated_unvalidated",
            score_source,
            json.dumps(result_json or {}),
            "unique",
        ),
    ).lastrowid
    assert result_id is not None
    entries = [
        (None, context_id, None, "prepared_receptor", "receptor", "", digest_a, {}),
        (None, None, ligand_id, "prepared_ligand", "ligand", "", digest_b, {}),
        (
            pair_id,
            context_id,
            ligand_id,
            "docking_pose",
            "pair",
            "stage3",
            digest_c,
            (
                {
                    **({"result_sha256": digest_d} if linked_pose else {}),
                    **(pose_metadata or {}),
                }
            ),
        ),
    ]
    for index, (
        artifact_pair,
        artifact_context,
        artifact_ligand,
        role,
        scope,
        stage,
        digest,
        metadata,
    ) in enumerate(entries):
        connection.execute(
            """INSERT INTO artifacts
            (run_id, receptor_context_id, pair_cell_id, ligand_id, stage,
             archive_path, member_name, artifact_role, artifact_scope, sha256,
             size_bytes, file_type, verified, artifact_json)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1, ?)""",
            (
                "run_a",
                artifact_context,
                artifact_pair,
                artifact_ligand,
                stage,
                f"archive_{index}.tar",
                f"member_{index}.pdbqt",
                role,
                scope,
                digest,
                10,
                "pdbqt",
                json.dumps(metadata),
            ),
        )
    connection.commit()
    connection.close()


def test_sqlite_contributor_pose_plan_requires_result_lineage(
    tmp_path: Path,
) -> None:
    database = tmp_path / "atlas.sqlite"
    _insert_sqlite_fixture(
        database,
        linked_pose=False,
        score_source="consensus_vs_decoy_z",
        result_json={"contributing_pose_sha256s": ["c" * 64]},
        pose_metadata={"score_pose_role": "contributing_pose"},
    )
    request = _request(stage="pose-validation", source_stage="")

    blocked_result = build_stage_repeat_plan(
        tmp_path,
        request,
        database=database,
        output=tmp_path / "db-blocked.json",
    )
    blocked_row = _ledger_rows(blocked_result)[0]
    assert "score_source_pose_identity_unresolved" in blocked_row["blockers"]

    connection = sqlite3.connect(database)
    connection.execute(
        "UPDATE artifacts SET artifact_json=? WHERE artifact_role='docking_pose'",
        (
            json.dumps(
                {
                    "result_sha256": "d" * 64,
                    "score_pose_role": "contributing_pose",
                }
            ),
        ),
    )
    connection.commit()
    connection.close()
    ready_result = build_stage_repeat_plan(
        tmp_path,
        request,
        database=database,
        output=tmp_path / "db-ready.json",
    )
    assert ready_result["plan"]["summary"] == {
        "selected": 1,
        "ready": 1,
        "blocked": 0,
    }


def test_partial_plan_requires_explicit_strategy(tmp_path: Path) -> None:
    source = _write_rows(
        tmp_path / "master_rows.csv",
        [_csv_pair_row(tmp_path, f"lig_{index}") for index in range(4)],
    )
    with pytest.raises(StageRepeatPlanError, match="partial plans require"):
        build_stage_repeat_plan(
            tmp_path,
            _request(fraction=0.5, count=None),
            master_rows=source,
            output=tmp_path / "strategy-missing.json",
        )
    with pytest.raises(StageRepeatPlanError, match="top-score selection"):
        build_stage_repeat_plan(
            tmp_path,
            _request(
                fraction=0.5,
                count=None,
                selection_strategy="top-score",
            ),
            master_rows=source,
            output=tmp_path / "top-score-pending.json",
        )


def test_plan_json_uses_bounded_jsonl_shards(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(stage_repeat_plan, "PAIR_LEDGER_SHARD_ROWS", 2)
    source = _write_rows(
        tmp_path / "master_rows.csv",
        [_csv_pair_row(tmp_path, f"lig_{index}") for index in range(5)],
    )
    result = build_stage_repeat_plan(
        tmp_path,
        _request(),
        master_rows=source,
        output=tmp_path / "sharded.json",
    )
    persisted = json.loads(Path(result["plan_path"]).read_text(encoding="utf-8"))
    assert "rows" not in persisted
    assert persisted["pair_ledger"]["row_count"] == 5
    assert [shard["row_count"] for shard in persisted["pair_ledger"]["shards"]] == [
        2,
        2,
        1,
    ]
    assert len(result["ledger_paths"]) == 3
    assert len(_ledger_rows(result)) == 5


def test_sqlite_source_refuses_above_bounded_cap_before_materialization(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    database = tmp_path / "atlas.sqlite"
    _insert_sqlite_fixture(database, linked_pose=True)
    monkeypatch.setattr(stage_repeat_plan, "MAX_IN_MEMORY_PLAN_PAIRS", 0)
    with pytest.raises(
        StageRepeatPlanError, match="streaming_stage_repeat_plan_required"
    ):
        build_stage_repeat_plan(
            tmp_path,
            _request(stage="pose-validation", source_stage=""),
            database=database,
            output=tmp_path / "too-large.json",
        )


def test_aggregate_score_requires_complete_contributor_pose_lineage(
    tmp_path: Path,
) -> None:
    database = tmp_path / "atlas.sqlite"
    _insert_sqlite_fixture(
        database, linked_pose=True, score_source="consensus_vs_decoy_z"
    )
    request = _request(stage="pose-validation", source_stage="")
    blocked_result = build_stage_repeat_plan(
        tmp_path,
        request,
        database=database,
        output=tmp_path / "aggregate-blocked.json",
    )
    assert (
        "score_source_pose_identity_unresolved"
        in _ledger_rows(blocked_result)[0]["blockers"]
    )

    connection = sqlite3.connect(database)
    connection.execute(
        "UPDATE result_attempts SET result_json=?",
        (json.dumps({"contributing_pose_sha256s": ["c" * 64]}),),
    )
    connection.execute(
        "UPDATE artifacts SET artifact_json=? WHERE artifact_role='docking_pose'",
        (
            json.dumps(
                {
                    "result_sha256": "d" * 64,
                    "score_pose_role": "contributing_pose",
                }
            ),
        ),
    )
    connection.commit()
    connection.close()

    ready_result = build_stage_repeat_plan(
        tmp_path,
        request,
        database=database,
        output=tmp_path / "aggregate-ready.json",
    )
    ready_row = _ledger_rows(ready_result)[0]
    assert ready_row["status"] == "ready"
    assert "selected_docking_pose" not in ready_row["artifacts"]
    contributors = ready_row["artifacts"]["contributing_docking_poses"]
    assert [pose["sha256"] for pose in contributors] == ["c" * 64]


def test_reconstructed_legacy_scorch_source_fails_closed_without_pose_id(
    tmp_path: Path,
) -> None:
    database = tmp_path / "atlas.sqlite"
    _insert_sqlite_fixture(database, linked_pose=True, score_source="")
    effective_source = f"reconstructed:{LEGACY_SCORCH_NEUTRAL_CNN_VS_DECOY_Z}"
    with sqlite3.connect(database) as connection:
        connection.execute(
            """UPDATE result_attempts
               SET final_score_source=NULL,
                   final_score_source_reconstructed=?""",
            (effective_source,),
        )

    result = build_stage_repeat_plan(
        tmp_path,
        _request(stage="pose-validation", source_stage=""),
        database=database,
        output=tmp_path / "legacy-scorch-blocked.json",
    )
    row = _ledger_rows(result)[0]
    assert row["selected_result"]["score_source"] == effective_source
    assert "score_source_pose_identity_unresolved" in row["blockers"]


def test_vina_requires_exact_source_stage(tmp_path: Path) -> None:
    source = _write_rows(
        tmp_path / "master_rows.csv", [_csv_pair_row(tmp_path, "lig_a")]
    )
    with pytest.raises(StageRepeatPlanError, match="source-stage is required"):
        build_stage_repeat_plan(
            tmp_path,
            _request(source_stage=""),
            master_rows=source,
            output=tmp_path / "missing-stage.json",
        )


def test_dispatch_stages_repeat_emits_json_plan(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    source = _write_rows(
        tmp_path / "master_rows.csv", [_csv_pair_row(tmp_path, "lig_a")]
    )
    output = tmp_path / "cli-plan.json"
    monkeypatch.setattr(stages._bindings, "repo_root", lambda: tmp_path)

    return_code = dispatch(
        [
            "stages",
            "repeat",
            "--run-id",
            "run_a",
            "--stage",
            "vina",
            "--all-receptors",
            "--all-ligands",
            "--all-contexts",
            "--source-stage",
            "stage1",
            "--master-rows",
            str(source),
            "--output",
            str(output),
            "--cpus",
            "8",
            "--json",
        ]
    )

    assert return_code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["plan_path"] == str(output)
    assert payload["selected"] == 1
    assert payload["ready"] == 1
    assert payload["blocked"] == 0
    assert payload["dry_run"] is True
    assert payload["execution_status"] == "not_executed"
    assert payload["executor_status"] == "planner_only"
    assert payload["reused_existing"] is False
    assert output.is_file()
