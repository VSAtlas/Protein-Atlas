from __future__ import annotations

import csv
import json
import logging
import os
import threading
import time
from pathlib import Path

from post_docking.rescoring.scorch_orchestration_types import OrchestrationDeps
from post_docking.rescoring.scorch_orchestration_finalize import (
    finalize_orchestration_combos,
)
from post_docking.rescoring.scorch_postprocess import aggregate_combo
from post_docking.rescoring.scorch_shards import (
    _adaptive_scorch_timeout_sec,
    _claim_next_record,
    _result_path,
    _split_timeout_record,
    ensure_shard_plan,
    execute_sharded_tasks,
    planned_completed_output_csvs,
    read_all_shard_plans,
    read_all_shard_records,
    read_all_shard_results,
    repair_missing_completed_shard_results,
    repair_missing_split_child_plan_records,
    summarize_shard_ledger,
    task_to_shard_record,
    try_claim_shard,
    write_shard_result,
)
from post_docking.rescoring.scorch_types import SelectionResult, StageSpec


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else ["Ligand_ID", "SCORCH_score"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _deps(score_stage):
    return OrchestrationDeps(
        scheduler_runtime_snapshot=lambda *_args, **_kwargs: (8, 0),
        discover_mode_dirs=lambda *_args, **_kwargs: {"dud": [], "fda": []},
        load_control_bases=lambda *_args, **_kwargs: set(),
        stage_dir_candidates=lambda *_args, **_kwargs: [],
        score_csv_for_spec=lambda *_args, **_kwargs: (None, [], False),
        select_top_bases_from_score_csv=lambda *_args, **_kwargs: SelectionResult(
            allowed_bases=set(),
            n_pool=0,
            k=0,
            controls_total=0,
            selected_stage_by_base={},
            selected_score_by_base={},
        ),
        adaptive_chunk_size=lambda **_kwargs: 1,
        chunk_allowed_bases=lambda bases, _size: [set(bases)],
        tail_split_allowed_chunks=lambda chunks, **_kwargs: list(chunks),
        score_stage=score_stage,
        task_estimate_seconds=lambda *_args, **_kwargs: 1.0,
        task_allowed_count=lambda task, *_args, **_kwargs: len(task[3]),
        split_task_for_rechunk=lambda task, *_args, **_kwargs: [task],
        elastic_scorch_threads=lambda **_kwargs: 1,
        aggregate_combo=lambda *_args, **_kwargs: None,
        annotate_scorch_z_scores=lambda *_args, **_kwargs: None,
        emit_task_event=lambda *_args, **_kwargs: None,
        emit_bench_event=lambda *_args, **_kwargs: None,
        summarize_chunk_partition=lambda *_args, **_kwargs: None,
    )


def test_adaptive_scorch_timeout_splits_small_shards_promptly() -> None:
    assert _adaptive_scorch_timeout_sec(1) == 300.0
    assert _adaptive_scorch_timeout_sec(5) == 600.0
    assert _adaptive_scorch_timeout_sec(48) >= 696.0
    assert _adaptive_scorch_timeout_sec(96) >= 1272.0
    assert _adaptive_scorch_timeout_sec(1000) == 2400.0


def test_shard_plan_reuses_equivalent_and_replaces_incompatible(tmp_path: Path) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(tmp_path / "manifests")}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task_a = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    task_b = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_b", "lig_c"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part999",
    )

    first = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task_a],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    same = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task_a],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    second = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task_b],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert same["shards"][0]["shard_id"] == first["shards"][0]["shard_id"]
    assert same["shards"][0]["chunk_tag"] == "part000"
    assert second["shards"][0]["shard_id"] != first["shards"][0]["shard_id"]
    assert second["shards"][0]["chunk_tag"] == "part999"


def test_global_ledger_repairs_missing_results_and_summarizes_all_plans(
    tmp_path: Path,
) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(tmp_path / "manifests")}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    combos = [("PDB1", "HOLO", "ph_7_0"), ("PDB2", "APO", "ph_8_0")]

    for idx, combo in enumerate(combos):
        task = (
            spec,
            combo,
            tmp_path / f"receptor_{idx}.pdbqt",
            {f"lig_{idx}"},
            set(),
            "fda",
            None,
            None,
            {},
            {},
            f"part{idx:03d}",
        )
        plan = ensure_shard_plan(
            cfg=cfg,
            run_id=run_root.name,
            combo=combo,
            tasks=[task],
            post_root=post_root,
            decoy_prefix="dud",
            top_fraction=0.10,
        )
        record = plan["shards"][0]
        _write_csv(
            Path(record["output_csv"]),
            [{"Ligand_ID": f"lig_{idx}", "source": "vina", "SCORCH_score": "1.0"}],
        )

    assert len(read_all_shard_plans(cfg, run_root.name)) == 2
    assert len(read_all_shard_records(cfg, run_root.name)) == 2
    assert read_all_shard_results(cfg, run_root.name) == {}

    repaired = repair_missing_completed_shard_results(
        cfg,
        run_root.name,
        logger=logging.getLogger("test.scorch.shards.ledger_repair"),
    )

    assert repaired == 2
    results = read_all_shard_results(cfg, run_root.name)
    assert len(results) == 2
    assert {payload["status"] for payload in results.values()} == {"completed"}
    assert all(payload["repaired_from_output"] is True for payload in results.values())
    assert all(payload["final"] is True for payload in results.values())
    assert all(payload["provisional"] is False for payload in results.values())

    summary = summarize_shard_ledger(cfg, run_root.name, write=True)

    assert summary["total"] == 2
    assert summary["completed"] == 2
    assert summary["row_total"] == 2
    assert summary["all_complete"] is True
    assert summary["premature"] == 0
    assert summary["provisional"] == 0
    assert {entry["ledger_state"] for entry in summary["shards"]} == {"completed"}
    assert all(entry["result_kind"] == "final" for entry in summary["shards"])


def test_planned_completed_output_csvs_requires_complete_combo_mode(
    tmp_path: Path,
) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(tmp_path / "manifests")}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    tasks = [
        (
            spec,
            combo,
            tmp_path / "receptor.pdbqt",
            {f"lig_{idx}"},
            set(),
            "fda",
            None,
            None,
            {},
            {},
            f"part{idx:03d}",
        )
        for idx in range(2)
    ]
    plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=tasks,
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    first, second = plan["shards"]
    _write_csv(
        Path(first["output_csv"]),
        [{"Ligand_ID": "lig_0", "source": "vina", "SCORCH_score": "1.0"}],
    )

    partial = planned_completed_output_csvs(
        cfg,
        run_root.name,
        combo,
        [spec],
        "fda",
        "dud",
    )

    assert partial == []
    assert len(read_all_shard_results(cfg, run_root.name)) == 1

    _write_csv(
        Path(second["output_csv"]),
        [{"Ligand_ID": "lig_1", "source": "vina", "SCORCH_score": "1.0"}],
    )
    complete = planned_completed_output_csvs(
        cfg,
        run_root.name,
        combo,
        [spec],
        "fda",
        "dud",
    )

    assert complete == [Path(first["output_csv"]), Path(second["output_csv"])]
    assert len(read_all_shard_results(cfg, run_root.name)) == 2


def test_plan_bound_aggregation_repairs_lost_split_children(
    tmp_path: Path,
) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(tmp_path / "manifests")}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "", "")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a", "lig_b"},
        set(),
        "dud",
        None,
        None,
        {},
        {},
        "part000",
    )
    plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task],
        post_root=post_root,
        decoy_prefix="decoys",
        top_fraction=0.10,
    )
    parent = plan["shards"][0]
    children = _split_timeout_record(
        run_id=run_root.name,
        record=parent,
        post_root=post_root,
        decoy_prefix="decoys",
        top_fraction=0.10,
        attempt=2,
    )
    write_shard_result(
        cfg,
        run_root.name,
        parent,
        status="split",
        payload={
            "attempt": 1,
            "failure_reason": "timeout_split",
            "child_shard_ids": [child["shard_id"] for child in children],
            "child_count": len(children),
        },
    )
    for child in children:
        ligands = list(child["allowed_bases"])
        assert len(ligands) == 1
        _write_csv(
            Path(child["output_csv"]),
            [{"Ligand_ID": ligands[0], "source": "vina", "SCORCH_score": "1.0"}],
        )
        write_shard_result(
            cfg,
            run_root.name,
            child,
            status="completed",
            payload={
                "attempt": 2,
                "output_csv": child["output_csv"],
                "row_count": 1,
                "expected_row_count": 1,
            },
        )

    assert len(read_all_shard_records(cfg, run_root.name)) == 1
    repaired = repair_missing_split_child_plan_records(cfg, run_root.name)
    assert repaired == 2
    assert len(read_all_shard_records(cfg, run_root.name)) == 3

    complete = planned_completed_output_csvs(
        cfg,
        run_root.name,
        combo,
        [spec],
        "dud",
        "decoys",
    )

    assert sorted(str(path) for path in complete or []) == sorted(
        str(child["output_csv"]) for child in children
    )

    _result_path(cfg, run_root.name, str(parent["shard_id"])).unlink()
    complete_without_parent_result = planned_completed_output_csvs(
        cfg,
        run_root.name,
        combo,
        [spec],
        "dud",
        "decoys",
    )

    assert sorted(str(path) for path in complete_without_parent_result or []) == sorted(
        str(child["output_csv"]) for child in children
    )


def test_planned_completed_output_csvs_keeps_valid_control_only_output(
    tmp_path: Path,
) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(tmp_path / "manifests")}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"ctrl_lig"},
        {"ctrl_lig"},
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    record = plan["shards"][0]
    assert record["expected_output_rows"] == 0
    _write_csv(
        Path(record["output_csv"]),
        [{"Ligand_ID": "ctrl_lig", "source": "vina", "SCORCH_score": "1.0"}],
    )

    complete = planned_completed_output_csvs(
        cfg,
        run_root.name,
        combo,
        [spec],
        "fda",
        "dud",
    )

    assert complete == [Path(record["output_csv"])]


def test_claim_next_record_skips_process_local_active_shard(tmp_path: Path) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(tmp_path / "manifests")}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "", "")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    tasks = [
        (
            spec,
            combo,
            tmp_path / "receptor.pdbqt",
            {"lig_a", "lig_b"},
            set(),
            "fda",
            None,
            None,
            {},
            {},
            "part000",
        ),
        (
            spec,
            combo,
            tmp_path / "receptor.pdbqt",
            {"lig_c"},
            set(),
            "fda",
            None,
            None,
            {},
            {},
            "part001",
        ),
    ]
    plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=tasks,
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    records = list(plan["shards"])
    first = max(records, key=lambda rec: int(rec.get("allowed_count") or 0))

    claimed = _claim_next_record(
        cfg,
        run_root.name,
        records,
        lease_sec=600.0,
        max_attempts=2,
        overwrite=True,
        logger=logging.getLogger("test.scorch.claim.skip"),
        skip_shard_ids={str(first["shard_id"])},
    )

    assert claimed is not None
    assert claimed[0]["shard_id"] != first["shard_id"]


def test_shard_rejects_stale_output_with_wrong_ligands(tmp_path: Path) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"expected_lig"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    stale_output = (
        post_root
        / combo[0]
        / combo[1]
        / combo[2]
        / "scorch_scores_vina_best.part000.csv"
    )
    _write_csv(
        stale_output,
        [{"Ligand_ID": "wrong_lig_stage1", "SCORCH_score": "0.1"}],
    )
    calls = {"score_stage": 0}

    def fake_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *_args,
    ):
        del (
            cfg_arg,
            run_root_arg,
            receptor_arg,
            threads_arg,
            logger_arg,
            control_bases,
            run_mode,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        calls["score_stage"] += 1
        output = (
            post_root_arg
            / combo_arg[0]
            / combo_arg[1]
            / combo_arg[2]
            / spec_arg.output_name.replace(".csv", f".{chunk_tag}.csv")
        )
        assert not output.exists()
        _write_csv(
            output,
            [{"Ligand_ID": next(iter(allowed_bases)), "SCORCH_score": "1.0"}],
        )
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.stale"),
        deps=_deps(fake_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete
    assert calls["score_stage"] == 1
    with stale_output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["Ligand_ID"] for row in rows] == ["expected_lig"]


def test_execute_sharded_tasks_writes_results_and_reuses_claims(tmp_path: Path) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    tasks = [
        (
            spec,
            combo,
            tmp_path / "receptor.pdbqt",
            {f"lig_{idx}"},
            set(),
            "fda",
            None,
            None,
            {},
            {},
            f"part{idx:03d}",
        )
        for idx in range(2)
    ]

    def fake_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *_args,
    ):
        del (
            cfg_arg,
            run_root_arg,
            receptor_arg,
            threads_arg,
            logger_arg,
            control_bases,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        output = (
            post_root_arg
            / combo_arg[0]
            / combo_arg[1]
            / combo_arg[2]
            / spec_arg.output_name.replace(".csv", f".{chunk_tag}.csv")
        )
        _write_csv(
            output,
            [
                {
                    "Ligand_ID": next(iter(allowed_bases)),
                    "source": run_mode,
                    "SCORCH_score": "1.0",
                }
            ],
        )
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=tasks,
        run_root=run_root,
        post_root=post_root,
        jobs=2,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards"),
        deps=_deps(fake_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.failed_jobs == 0
    assert result.completed == 2
    result_files = sorted((tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_results").glob("*.json"))
    assert len(result_files) == 2
    summary_path = tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["productive"] == 2
    assert summary["zero_output"] == 0
    assert summary["row_total"] == 2


def test_execute_sharded_tasks_ignores_unrelated_failed_ledger_scope(
    tmp_path: Path,
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    stale_combo = ("STALE", "HOLO", "ph_7_0")
    stale_task = (
        spec,
        stale_combo,
        tmp_path / "receptor.pdbqt",
        {"stale_lig"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    stale_plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=stale_combo,
        tasks=[stale_task],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    stale_record = stale_plan["shards"][0]
    result_path = (
        tmp_path
        / "manifests"
        / run_root.name
        / "distributed"
        / "scorch_shard_results"
        / f"{stale_record['shard_id']}.json"
    )
    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.write_text(
        json.dumps(
            {
                "status": "failed",
                "attempt": 2,
                "shard_id": stale_record["shard_id"],
                "failure_reason": "unrelated",
            }
        ),
        encoding="utf-8",
    )
    combo = ("PDB1", "HOLO", "ph_7_0")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"fresh_lig"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )

    def fake_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *_args,
    ):
        del (
            cfg_arg,
            run_root_arg,
            receptor_arg,
            threads_arg,
            logger_arg,
            control_bases,
            run_mode,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        output = (
            post_root_arg
            / combo_arg[0]
            / combo_arg[1]
            / combo_arg[2]
            / spec_arg.output_name.replace(".csv", f".{chunk_tag}.csv")
        )
        _write_csv(
            output,
            [{"Ligand_ID": next(iter(allowed_bases)), "source": "vina", "SCORCH_score": "1.0"}],
        )
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.scoped"),
        deps=_deps(fake_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.failed_jobs == 0
    assert result.completed == 1
    global_summary = summarize_shard_ledger(cfg, run_root.name)
    assert global_summary["failed"] == 1


def test_try_claim_shard_retries_completed_record_with_missing_output(
    tmp_path: Path,
) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(tmp_path / "manifests")}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    record = plan["shards"][0]
    result_path = (
        tmp_path
        / "manifests"
        / run_root.name
        / "distributed"
        / "scorch_shard_results"
        / f"{record['shard_id']}.json"
    )
    _write_csv(post_root / "placeholder.csv", [{"Ligand_ID": "x", "SCORCH_score": "1"}])
    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.write_text(
        json.dumps(
            {
                "status": "completed",
                "attempt": 1,
                "shard_id": record["shard_id"],
                "output_csv": str(post_root / "missing.csv"),
            }
        ),
        encoding="utf-8",
    )

    claimed, attempt = try_claim_shard(
        cfg,
        run_root.name,
        record,
        lease_sec=60.0,
        max_attempts=2,
    )

    assert claimed is True
    assert attempt == 2
    assert not result_path.exists()


def test_try_claim_shard_preserves_quarantined_record_by_default(
    tmp_path: Path,
) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(tmp_path / "manifests")}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    record = plan["shards"][0]
    output_path = Path(str(record["output_csv"]))
    _write_csv(
        output_path,
        [
            {
                "Ligand_ID": "lig_a",
                "SCORCH_score": "",
                "scorch_status": "quarantined",
            }
        ],
    )
    result_path = (
        tmp_path
        / "manifests"
        / run_root.name
        / "distributed"
        / "scorch_shard_results"
        / f"{record['shard_id']}.json"
    )
    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.write_text(
        json.dumps(
            {
                "status": "quarantined",
                "attempt": 2,
                "shard_id": record["shard_id"],
                "output_csv": str(output_path),
                "quarantine_output": True,
            }
        ),
        encoding="utf-8",
    )

    claimed, attempt = try_claim_shard(
        cfg,
        run_root.name,
        record,
        lease_sec=60.0,
        max_attempts=2,
    )

    assert claimed is False
    assert attempt == 2
    assert result_path.exists()
    assert output_path.exists()


def test_try_claim_shard_retries_quarantined_record_when_enabled(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setenv("ATLAS_SCORCH_RETRY_QUARANTINED_ON_RESUME", "1")
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(tmp_path / "manifests")}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    record = plan["shards"][0]
    output_path = Path(str(record["output_csv"]))
    _write_csv(
        output_path,
        [
            {
                "Ligand_ID": "lig_a",
                "SCORCH_score": "",
                "scorch_status": "quarantined",
            }
        ],
    )
    result_path = (
        tmp_path
        / "manifests"
        / run_root.name
        / "distributed"
        / "scorch_shard_results"
        / f"{record['shard_id']}.json"
    )
    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.write_text(
        json.dumps(
            {
                "status": "quarantined",
                "attempt": 2,
                "shard_id": record["shard_id"],
                "output_csv": str(output_path),
                "quarantine_output": True,
            }
        ),
        encoding="utf-8",
    )

    claimed, attempt = try_claim_shard(
        cfg,
        run_root.name,
        record,
        lease_sec=60.0,
        max_attempts=2,
    )

    assert claimed is True
    assert attempt == 3
    assert not result_path.exists()
    assert not output_path.exists()


def test_execute_sharded_tasks_accounts_explicit_quarantine_without_failed_job(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setenv("ATLAS_SCORCH_LIGAND_QUARANTINE_BASES", "lig_bad")
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_bad"},
        set(),
        "fda",
        None,
        None,
        {},
        {"lig_bad": -12.3},
        "part000",
    )
    calls = {"score_stage": 0}

    def should_not_score(*_args, **_kwargs):
        calls["score_stage"] += 1
        raise AssertionError("explicit quarantine should not invoke SCORCH")

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.explicit_quarantine"),
        deps=_deps(should_not_score),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.failed_jobs == 0
    assert result.completed == 1
    assert calls["score_stage"] == 0

    summary = json.loads(
        (tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["all_complete"] is True
    assert summary["all_scored"] is False
    assert summary["accounted"] == 1
    assert summary["scored_completed"] == 0
    assert summary["failed"] == 0
    assert summary["quarantined"] == 1
    assert summary["quarantined_rows"] == 1

    paths = planned_completed_output_csvs(
        cfg,
        run_root.name,
        combo,
        [spec],
        "fda",
        "dud",
    )
    assert paths == [Path(task_to_shard_record(
        run_id=run_root.name,
        task=task,
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )["output_csv"])]
    with paths[0].open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["Ligand_ID"] == "lig_bad"
    assert rows[0]["scorch_status"] == "quarantined"
    assert rows[0]["scorch_unscorable_flag"] == "1"


def test_try_claim_shard_reclaims_old_foreign_slurm_claim(
    tmp_path: Path,
    monkeypatch,
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_ORPHAN_CLAIM_SEC": 60,
    }
    monkeypatch.setenv("SLURM_JOB_ID", "new_job")
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    record = plan["shards"][0]
    claim_path = (
        tmp_path
        / "manifests"
        / run_root.name
        / "distributed"
        / "scorch_shard_claims"
        / f"{record['shard_id']}.claim.json"
    )
    claim_path.parent.mkdir(parents=True, exist_ok=True)
    claim_path.write_text(
        json.dumps(
            {
                "run_id": run_root.name,
                "shard_id": record["shard_id"],
                "hostname": "old-node",
                "pid": 12345,
                "slurm_job_id": "old_job",
                "lease_sec": 7200,
            }
        ),
        encoding="utf-8",
    )
    old_mtime = time.time() - 120
    os.utime(claim_path, (old_mtime, old_mtime))

    claimed, attempt = try_claim_shard(
        cfg,
        run_root.name,
        record,
        lease_sec=7200.0,
        max_attempts=2,
    )

    assert claimed is True
    assert attempt == 1
    payload = json.loads(claim_path.read_text(encoding="utf-8"))
    assert payload["slurm_job_id"] == "new_job"


def test_try_claim_shard_preserves_same_slurm_job_claim(
    tmp_path: Path,
    monkeypatch,
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_ORPHAN_CLAIM_SEC": 60,
    }
    monkeypatch.setenv("SLURM_JOB_ID", "same_job")
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    plan = ensure_shard_plan(
        cfg=cfg,
        run_id=run_root.name,
        combo=combo,
        tasks=[task],
        post_root=post_root,
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    record = plan["shards"][0]
    claim_path = (
        tmp_path
        / "manifests"
        / run_root.name
        / "distributed"
        / "scorch_shard_claims"
        / f"{record['shard_id']}.claim.json"
    )
    claim_path.parent.mkdir(parents=True, exist_ok=True)
    claim_path.write_text(
        json.dumps(
            {
                "run_id": run_root.name,
                "shard_id": record["shard_id"],
                "hostname": "other-node",
                "pid": 12345,
                "slurm_job_id": "same_job",
                "lease_sec": 7200,
            }
        ),
        encoding="utf-8",
    )
    old_mtime = time.time() - 120
    os.utime(claim_path, (old_mtime, old_mtime))

    claimed, attempt = try_claim_shard(
        cfg,
        run_root.name,
        record,
        lease_sec=7200.0,
        max_attempts=2,
    )

    assert claimed is False
    assert attempt == 1
    payload = json.loads(claim_path.read_text(encoding="utf-8"))
    assert payload["slurm_job_id"] == "same_job"


def test_execute_sharded_tasks_retries_failed_shard_and_completes_tail(
    tmp_path: Path,
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    tasks = [
        (
            spec,
            combo,
            tmp_path / "receptor.pdbqt",
            {f"lig_{idx}"},
            set(),
            "fda",
            None,
            None,
            {},
            {},
            f"part{idx:03d}",
        )
        for idx in range(2)
    ]
    calls: dict[str, int] = {}

    def flaky_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *_args,
    ):
        del (
            cfg_arg,
            run_root_arg,
            receptor_arg,
            threads_arg,
            logger_arg,
            control_bases,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        calls[chunk_tag] = calls.get(chunk_tag, 0) + 1
        if chunk_tag == "part000" and calls[chunk_tag] == 1:
            return False, None
        output = post_root_arg / combo_arg[0] / combo_arg[1] / combo_arg[2] / spec_arg.output_name.replace(
            ".csv", f".{chunk_tag}.csv"
        )
        _write_csv(
            output,
            [{"Ligand_ID": next(iter(allowed_bases)), "source": run_mode, "SCORCH_score": "1.0"}],
        )
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=tasks,
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.retry"),
        deps=_deps(flaky_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.failed_jobs == 0
    assert result.completed == 2
    assert calls["part000"] == 2
    assert calls["part001"] == 1


def test_execute_sharded_tasks_splits_large_timeout_shard(tmp_path: Path) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
        "ATLAS_SCORCH_SHARD_TIMEOUT_SPLIT_MIN": 3,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a", "lig_b", "lig_c", "lig_d"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    calls: list[str] = []

    def timeout_then_children_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *args,
    ):
        del (
            run_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            control_bases,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        shard_id = str(args[0])
        calls.append(str(chunk_tag))
        if chunk_tag == "part000":
            cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
                "status": "error",
                "returncode": 124,
                "failure_reason": "timeout",
                "timed_out": True,
                "timeout_sec": 600.0,
            }
            return False, None
        output = post_root_arg / combo_arg[0] / combo_arg[1] / combo_arg[2] / spec_arg.output_name.replace(
            ".csv", f".{chunk_tag}.csv"
        )
        _write_csv(
            output,
            [
                {"Ligand_ID": base, "source": run_mode, "SCORCH_score": "1.0"}
                for base in sorted(allowed_bases)
            ],
        )
        cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
            "status": "ok",
            "returncode": 0,
            "failure_reason": "",
            "timed_out": False,
            "timeout_sec": 600.0,
        }
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.timeout_split"),
        deps=_deps(timeout_then_children_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.failed_jobs == 0
    assert calls[0] == "part000"
    assert sorted(calls[1:]) == [
        "part000t02r00",
        "part000t02r01",
        "part000t02r02",
        "part000t02r03",
    ]
    result_payloads = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted((tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_results").glob("*.json"))
    ]
    assert sorted(payload["status"] for payload in result_payloads) == [
        "completed",
        "completed",
        "completed",
        "completed",
        "split",
    ]
    summary = json.loads(
        (tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["split"] == 1
    assert summary["productive"] == 4
    assert summary["all_complete"] is True


def test_overwrite_preserves_split_parent_and_summary_supersedes_failed_parent(
    tmp_path: Path,
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
        "ATLAS_SCORCH_SHARD_TIMEOUT_SPLIT_MIN": 3,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a", "lig_b", "lig_c", "lig_d"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    calls: list[str] = []

    def split_then_score_children(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *args,
    ):
        del (
            run_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            control_bases,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        shard_id = str(args[0])
        calls.append(str(chunk_tag))
        if chunk_tag == "part000":
            cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
                "status": "error",
                "returncode": 124,
                "failure_reason": "timeout",
                "timed_out": True,
                "timeout_sec": 600.0,
            }
            return False, None
        output = post_root_arg / combo_arg[0] / combo_arg[1] / combo_arg[2] / spec_arg.output_name.replace(
            ".csv", f".{chunk_tag}.csv"
        )
        _write_csv(
            output,
            [
                {"Ligand_ID": base, "source": run_mode, "SCORCH_score": "1.0"}
                for base in sorted(allowed_bases)
            ],
        )
        cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
            "status": "ok",
            "returncode": 0,
            "failure_reason": "",
            "timed_out": False,
            "timeout_sec": 600.0,
        }
        return True, output

    first = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.overwrite_split.initial"),
        deps=_deps(split_then_score_children),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    assert first.all_complete is True
    assert first.failed_jobs == 0
    assert calls[0] == "part000"

    calls.clear()

    def overwrite_reruns_children_only(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *args,
    ):
        if chunk_tag == "part000":
            raise AssertionError("overwrite must not re-run split parent shards")
        return split_then_score_children(
            cfg_arg,
            spec_arg,
            combo_arg,
            run_root_arg,
            post_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            allowed_bases,
            control_bases,
            run_mode,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
            chunk_tag,
            *args,
        )

    second = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=True,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.overwrite_split.second"),
        deps=_deps(overwrite_reruns_children_only),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )
    assert second.all_complete is True
    assert second.failed_jobs == 0
    assert "part000" not in calls

    records = read_all_shard_records(cfg, "rid")
    parent = next(record for record in records if record.get("chunk_tag") == "part000")
    write_shard_result(
        cfg,
        "rid",
        parent,
        status="failed",
        payload={"attempt": 2, "failure_reason": "timeout", "timed_out": True},
    )
    summary = summarize_shard_ledger(cfg, "rid", records=records, write=True)
    assert summary["failed"] == 0
    assert summary["superseded"] == 1
    assert summary["all_complete"] is True
    parent_entry = next(item for item in summary["shards"] if item["shard_id"] == parent["shard_id"])
    assert parent_entry["ledger_state"] == "superseded"
    assert parent_entry["superseded_by_children"] is True


def test_execute_sharded_tasks_splits_small_timeout_shard_by_default(tmp_path: Path) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a", "lig_b", "lig_c", "lig_d", "lig_e"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    calls: list[str] = []

    def timeout_then_children_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *args,
    ):
        del (
            run_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            control_bases,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        shard_id = str(args[0])
        calls.append(str(chunk_tag))
        if chunk_tag == "part000":
            cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
                "status": "error",
                "returncode": 124,
                "failure_reason": "timeout",
                "timed_out": True,
                "timeout_sec": 600.0,
            }
            return False, None
        output = (
            post_root_arg
            / combo_arg[0]
            / combo_arg[1]
            / combo_arg[2]
            / spec_arg.output_name.replace(".csv", f".{chunk_tag}.csv")
        )
        _write_csv(
            output,
            [
                {"Ligand_ID": base, "source": run_mode, "SCORCH_score": "1.0"}
                for base in sorted(allowed_bases)
            ],
        )
        cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
            "status": "ok",
            "returncode": 0,
            "failure_reason": "",
            "timed_out": False,
            "timeout_sec": 600.0,
        }
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.small_timeout_split"),
        deps=_deps(timeout_then_children_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert calls[0] == "part000"
    assert sorted(calls[1:]) == [
        "part000t02r00",
        "part000t02r01",
        "part000t02r02",
        "part000t02r03",
        "part000t02r04",
    ]
    summary = json.loads(
        (tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["split"] == 1
    assert summary["all_complete"] is True


def test_execute_sharded_tasks_fans_out_timeout_tail_to_keep_workers_busy(
    tmp_path: Path,
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {f"lig_{idx:02d}" for idx in range(16)},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "parent",
    )
    lock = threading.Lock()
    enough_children_started = threading.Event()
    active_children = 0
    max_active_children = 0
    child_calls: list[str] = []
    min_parallel_children = 4

    def timeout_then_fanout_children_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *args,
    ):
        nonlocal active_children, max_active_children
        del (
            run_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            control_bases,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        shard_id = str(args[0])
        if chunk_tag == "parent":
            cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
                "status": "error",
                "returncode": 124,
                "failure_reason": "timeout",
                "timed_out": True,
                "timeout_sec": 300.0,
            }
            return False, None
        if str(chunk_tag).startswith("parentt02"):
            with lock:
                child_calls.append(str(chunk_tag))
                active_children += 1
                max_active_children = max(max_active_children, active_children)
                if active_children >= min_parallel_children:
                    enough_children_started.set()
            enough_children_started.wait(timeout=2.0)
        output = (
            post_root_arg
            / combo_arg[0]
            / combo_arg[1]
            / combo_arg[2]
            / spec_arg.output_name.replace(".csv", f".{chunk_tag}.csv")
        )
        _write_csv(
            output,
            [
                {"Ligand_ID": base, "source": run_mode, "SCORCH_score": "1.0"}
                for base in sorted(allowed_bases)
            ],
        )
        if str(chunk_tag).startswith("parentt02"):
            with lock:
                active_children -= 1
        cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
            "status": "ok",
            "returncode": 0,
            "failure_reason": "",
            "timed_out": False,
            "timeout_sec": 300.0,
        }
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=16,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.timeout_fanout_workers"),
        deps=_deps(timeout_then_fanout_children_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert len(child_calls) == 16
    assert max_active_children >= min_parallel_children
    summary = json.loads(
        (tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["split"] == 1
    assert summary["productive"] == 16
    assert summary["all_complete"] is True


def test_execute_sharded_tasks_keeps_idle_workers_for_runtime_splits(
    tmp_path: Path,
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    tasks = [
        (
            spec,
            combo,
            tmp_path / "receptor.pdbqt",
            {"lig_fast"},
            set(),
            "fda",
            None,
            None,
            {},
            {},
            "fast",
        ),
        (
            spec,
            combo,
            tmp_path / "receptor.pdbqt",
            {"lig_a", "lig_b"},
            set(),
            "fda",
            None,
            None,
            {},
            {},
            "parent",
        ),
    ]
    lock = threading.Lock()
    second_child_started = threading.Event()
    active_children = 0
    max_active_children = 0

    def timeout_then_child_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *args,
    ):
        nonlocal active_children, max_active_children
        del (
            run_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            control_bases,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        shard_id = str(args[0])
        if chunk_tag == "parent":
            cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
                "status": "error",
                "returncode": 124,
                "failure_reason": "timeout",
                "timed_out": True,
                "timeout_sec": 300.0,
            }
            return False, None
        if str(chunk_tag).startswith("parentt02"):
            with lock:
                active_children += 1
                max_active_children = max(max_active_children, active_children)
                if active_children >= 2:
                    second_child_started.set()
            second_child_started.wait(timeout=2.0)
        output = (
            post_root_arg
            / combo_arg[0]
            / combo_arg[1]
            / combo_arg[2]
            / spec_arg.output_name.replace(".csv", f".{chunk_tag}.csv")
        )
        _write_csv(
            output,
            [
                {"Ligand_ID": base, "source": run_mode, "SCORCH_score": "1.0"}
                for base in sorted(allowed_bases)
            ],
        )
        if str(chunk_tag).startswith("parentt02"):
            with lock:
                active_children -= 1
        cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
            "status": "ok",
            "returncode": 0,
            "failure_reason": "",
            "timed_out": False,
            "timeout_sec": 300.0,
        }
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=tasks,
        run_root=run_root,
        post_root=post_root,
        jobs=2,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.idle_split_workers"),
        deps=_deps(timeout_then_child_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert max_active_children == 2


def test_execute_sharded_tasks_does_not_split_single_expected_timeout(
    tmp_path: Path,
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
        "ATLAS_SCORCH_SHARD_MAX_ATTEMPTS": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_timeout", "ctrl_a", "ctrl_b", "ctrl_c"},
        {"ctrl_a", "ctrl_b", "ctrl_c"},
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )

    def timeout_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *args,
    ):
        del (
            spec_arg,
            combo_arg,
            run_root_arg,
            post_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            allowed_bases,
            control_bases,
            run_mode,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
            chunk_tag,
        )
        shard_id = str(args[0])
        cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
            "status": "error",
            "returncode": 124,
            "failure_reason": "timeout",
            "timed_out": True,
            "timeout_sec": 600.0,
        }
        return False, None

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.single_expected_timeout"),
        deps=_deps(timeout_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.failed_jobs == 0
    summary = json.loads(
        (tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["split"] == 0
    assert summary["completed"] == 1
    assert summary["scored_completed"] == 0
    assert summary["accounted"] == 1
    assert summary["failed"] == 0
    assert summary["quarantined"] == 1
    assert summary["quarantined_rows"] == 1
    assert summary["all_complete"] is True
    assert summary["all_scored"] is False
    output = (
        post_root
        / combo[0]
        / combo[1]
        / combo[2]
        / "scorch_scores_vina_best.part000.csv"
    )
    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["Ligand_ID"] == "lig_timeout"
    assert rows[0]["SCORCH_score"] == "0"
    assert rows[0]["scorch_unscorable_flag"] == "1"


def test_execute_sharded_tasks_quarantines_repeat_tiny_timeout(tmp_path: Path) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
        "ATLAS_SCORCH_SHARD_TIMEOUT_SPLIT_MIN": 2,
        "ATLAS_SCORCH_SHARD_MAX_ATTEMPTS": 2,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_timeout"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    calls = 0

    def repeated_timeout_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *args,
    ):
        nonlocal calls
        del (
            spec_arg,
            combo_arg,
            run_root_arg,
            post_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            allowed_bases,
            control_bases,
            run_mode,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
            chunk_tag,
        )
        calls += 1
        shard_id = str(args[0])
        cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
            "status": "error",
            "returncode": 124,
            "failure_reason": "timeout",
            "timed_out": True,
            "timeout_sec": 600.0 * calls,
        }
        return False, None

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.quarantine"),
        deps=_deps(repeated_timeout_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.failed_jobs == 0
    assert calls == 2
    result_files = sorted((tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_results").glob("*.json"))
    result_payload = json.loads(result_files[0].read_text(encoding="utf-8"))
    assert result_payload["status"] == "quarantined"
    assert result_payload["failure_reason"] == "timeout"
    assert result_payload["returncode"] == 124
    assert result_payload["quarantine_output"] is True
    summary = json.loads(
        (tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["completed"] == 1
    assert summary["scored_completed"] == 0
    assert summary["accounted"] == 1
    assert summary["failed"] == 0
    assert summary["quarantined"] == 1
    assert summary["quarantined_rows"] == 1
    assert summary["all_complete"] is True
    assert summary["all_scored"] is False
    assert summary["shards"][0]["ledger_state"] == "quarantined"
    assert summary["shards"][0]["final"] is False


def test_execute_sharded_tasks_globally_quarantines_repeated_ligand_timeout(
    tmp_path: Path,
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
        "ATLAS_SCORCH_SHARD_TIMEOUT_SPLIT_MIN": 2,
        "ATLAS_SCORCH_SHARD_MAX_ATTEMPTS": 4,
        "ATLAS_SCORCH_LIGAND_TIMEOUT_QUARANTINE_THRESHOLD": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    tasks = [
        (
            spec,
            ("PDB1", "HOLO", "ph_7_0"),
            tmp_path / "receptor1.pdbqt",
            {"lig_timeout"},
            set(),
            "fda",
            None,
            None,
            {},
            {"lig_timeout": 1.0},
            "part000",
        ),
        (
            spec,
            ("PDB2", "HOLO", "ph_7_0"),
            tmp_path / "receptor2.pdbqt",
            {"lig_timeout"},
            set(),
            "fda",
            None,
            None,
            {},
            {"lig_timeout": 0.9},
            "part000",
        ),
    ]
    calls = 0

    def repeated_ligand_timeout_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *args,
    ):
        nonlocal calls
        del (
            spec_arg,
            combo_arg,
            run_root_arg,
            post_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            allowed_bases,
            control_bases,
            run_mode,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
            chunk_tag,
        )
        calls += 1
        shard_id = str(args[0])
        cfg_arg.setdefault("_ATLAS_SCORCH_STAGE_RESULTS", {})[shard_id] = {
            "status": "error",
            "returncode": 124,
            "failure_reason": "timeout",
            "timed_out": True,
            "timeout_sec": 300.0,
        }
        return False, None

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=tasks,
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.global_quarantine"),
        deps=_deps(repeated_ligand_timeout_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.failed_jobs == 0
    assert calls == 1
    results = read_all_shard_results(cfg, "rid")
    assert len(results) == 2
    assert {payload["status"] for payload in results.values()} == {"quarantined"}
    assert any(payload.get("global_ligand_quarantine") for payload in results.values())
    summary = json.loads(
        (tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["quarantined"] == 2
    assert summary["quarantined_rows"] == 2
    assert summary["completed"] == 2
    assert summary["scored_completed"] == 0
    assert summary["accounted"] == 2
    assert summary["failed"] == 0
    assert summary["all_complete"] is True
    assert summary["all_scored"] is False


def test_execute_sharded_tasks_retries_partial_success_output(tmp_path: Path) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"lig_a", "lig_b"},
        set(),
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )
    calls = 0

    def partial_then_complete_score_stage(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *_args,
    ):
        nonlocal calls
        del (
            cfg_arg,
            run_root_arg,
            receptor_arg,
            threads_arg,
            logger_arg,
            control_bases,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        calls += 1
        output = post_root_arg / combo_arg[0] / combo_arg[1] / combo_arg[2] / spec_arg.output_name.replace(
            ".csv", f".{chunk_tag}.csv"
        )
        if calls == 2:
            assert not output.exists()
            assert overwrite_arg is True
        else:
            assert overwrite_arg is False
        bases = sorted(allowed_bases)
        rows = bases[:1] if calls == 1 else bases
        _write_csv(
            output,
            [
                {"Ligand_ID": base, "source": run_mode, "SCORCH_score": "1.0"}
                for base in rows
            ],
        )
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.partial"),
        deps=_deps(partial_then_complete_score_stage),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.completed == 1
    assert calls == 2
    result_files = sorted((tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_results").glob("*.json"))
    result_payload = json.loads(result_files[0].read_text(encoding="utf-8"))
    assert result_payload["status"] == "completed"
    assert result_payload["row_count"] == 2
    assert result_payload["expected_row_count"] == 2


def test_execute_sharded_tasks_excludes_controls_from_expected_rows(tmp_path: Path) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC": 1,
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    task = (
        spec,
        combo,
        tmp_path / "receptor.pdbqt",
        {"native_control", "lig_a"},
        {"native_control"},
        "fda",
        None,
        None,
        {},
        {},
        "part000",
    )

    def score_with_control_rows(
        cfg_arg,
        spec_arg,
        combo_arg,
        run_root_arg,
        post_root_arg,
        receptor_arg,
        threads_arg,
        overwrite_arg,
        logger_arg,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
        *_args,
    ):
        del (
            cfg_arg,
            run_root_arg,
            receptor_arg,
            threads_arg,
            overwrite_arg,
            logger_arg,
            stage_dirs_override,
            score_csv,
            selected_stage_by_base,
            selected_score_by_base,
        )
        output = post_root_arg / combo_arg[0] / combo_arg[1] / combo_arg[2] / spec_arg.output_name.replace(
            ".csv", f".{chunk_tag}.csv"
        )
        rows = [
            {
                "Ligand_ID": (
                    f"{base}.sanitized" if base in set(control_bases) else base
                ),
                "source": run_mode,
                "SCORCH_score": "1.0",
            }
            for base in sorted(set(allowed_bases))
        ]
        _write_csv(output, rows)
        return True, output

    result = execute_sharded_tasks(
        cfg=cfg,
        tasks=[task],
        run_root=run_root,
        post_root=post_root,
        jobs=1,
        threads=1,
        overwrite=False,
        component="[test]",
        logger=logging.getLogger("test.scorch.shards.controls"),
        deps=_deps(score_with_control_rows),
        combo_failed_local={},
        decoy_prefix="dud",
        top_fraction=0.10,
    )

    assert result.all_complete is True
    assert result.completed == 1
    result_files = sorted((tmp_path / "manifests" / "rid" / "distributed" / "scorch_shard_results").glob("*.json"))
    result_payload = json.loads(result_files[0].read_text(encoding="utf-8"))
    assert result_payload["row_count"] == 2
    assert result_payload["expected_row_count"] == 1


def test_aggregate_combo_includes_tail_and_runtime_rechunk_outputs(tmp_path: Path) -> None:
    combo = ("PDB1", "HOLO", "ph_7_0")
    combo_dir = tmp_path / "post" / "rid" / combo[0] / combo[1] / combo[2]
    _write_csv(
        combo_dir / "scorch_scores_vina_best.tailt00.csv",
        [{"Ligand_ID": "a.pdbqt", "source": "vina", "SCORCH_score": "1.0"}],
    )
    _write_csv(
        combo_dir / "scorch_scores_vina_best.rt001r00.csv",
        [{"Ligand_ID": "b.pdbqt", "source": "vina", "SCORCH_score": "2.0"}],
    )

    out = aggregate_combo(
        tmp_path / "post" / "rid",
        [StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")],
        combo,
        logging.getLogger("test.scorch.aggregate"),
        component="[test]",
    )

    assert out is not None
    with out.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["Ligand_ID"] for row in rows} == {"a.pdbqt", "b.pdbqt"}


def test_aggregate_combo_plan_bound_excludes_unplanned_stale_parts(tmp_path: Path) -> None:
    combo = ("PDB1", "HOLO", "ph_7_0")
    combo_dir = tmp_path / "post" / "rid" / combo[0] / combo[1] / combo[2]
    planned = combo_dir / "scorch_scores_vina_best.part000.csv"
    stale = combo_dir / "scorch_scores_vina_best.part999.csv"
    _write_csv(
        planned,
        [{"Ligand_ID": "planned.pdbqt", "source": "vina", "SCORCH_score": "1.0"}],
    )
    _write_csv(
        stale,
        [{"Ligand_ID": "stale.pdbqt", "source": "vina", "SCORCH_score": "9.0"}],
    )

    out = aggregate_combo(
        tmp_path / "post" / "rid",
        [StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")],
        combo,
        logging.getLogger("test.scorch.aggregate.plan_bound"),
        component="[test]",
        input_csvs=[planned],
    )

    assert out is not None
    with out.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["Ligand_ID"] for row in rows] == ["planned.pdbqt"]


def test_aggregate_combo_plan_bound_filters_unplanned_rows(tmp_path: Path) -> None:
    combo = ("PDB1", "HOLO", "ph_7_0")
    combo_dir = tmp_path / "post" / "rid" / combo[0] / combo[1] / combo[2]
    planned = combo_dir / "scorch_scores_vina_best.part000.csv"
    _write_csv(
        planned,
        [
            {"Ligand_ID": "planned.pdbqt", "source": "vina", "SCORCH_score": "1.0"},
            {
                "Ligand_ID": "CLR_P401.sanitized_stage1",
                "source": "vina",
                "SCORCH_score": "9.0",
            },
        ],
    )

    out = aggregate_combo(
        tmp_path / "post" / "rid",
        [StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")],
        combo,
        logging.getLogger("test.scorch.aggregate.plan_filter"),
        component="[test]",
        input_csvs=[planned],
        allowed_ligand_bases={"planned"},
    )

    assert out is not None
    with out.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["Ligand_ID"] for row in rows] == ["planned.pdbqt"]


def test_finalize_regenerates_stale_reranked_outputs(tmp_path: Path) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path)}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "", "")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    consensus = run_root / "PDB1" / "consensus_docking_scores.csv"
    raw_scorch = post_root / "PDB1" / "scorch_scores_all.csv"
    _write_csv(
        consensus,
        [
            {
                "ligand": "lig_a.pdbqt",
                "run_mode": "fda",
                "consensus_score": "1.0",
                "best_engine": "vina",
            }
        ],
    )
    _write_csv(
        raw_scorch,
        [
            {
                "Ligand_ID": "lig_a_stage1",
                "run_mode": "fda",
                "source": "vina",
                "SCORCH_score": "0.8",
                "SCORCH_certainty": "0.9",
            }
        ],
    )
    stale_out = post_root / "PDB1" / "consensus_reranked_scorch.csv"
    _write_csv(
        stale_out,
        [{"ligand": "lig_a.pdbqt", "run_mode": "fda", "rescored_flag": "0"}],
    )

    calls: list[dict[str, object]] = []
    deps = OrchestrationDeps(
        scheduler_runtime_snapshot=lambda *_args, **_kwargs: (8, 0),
        discover_mode_dirs=lambda *_args, **_kwargs: {"dud": [], "fda": []},
        load_control_bases=lambda *_args, **_kwargs: set(),
        stage_dir_candidates=lambda *_args, **_kwargs: [],
        score_csv_for_spec=lambda *_args, **_kwargs: (None, [], False),
        select_top_bases_from_score_csv=lambda *_args, **_kwargs: SelectionResult(
            allowed_bases=set(),
            n_pool=0,
            k=0,
            controls_total=0,
            selected_stage_by_base={},
            selected_score_by_base={},
        ),
        adaptive_chunk_size=lambda **_kwargs: 1,
        chunk_allowed_bases=lambda bases, _size: [set(bases)],
        tail_split_allowed_chunks=lambda chunks, **_kwargs: list(chunks),
        score_stage=lambda *_args, **_kwargs: (True, raw_scorch),
        task_estimate_seconds=lambda *_args, **_kwargs: 1.0,
        task_allowed_count=lambda task, *_args, **_kwargs: len(task[3]),
        split_task_for_rechunk=lambda task, *_args, **_kwargs: [task],
        elastic_scorch_threads=lambda **_kwargs: 1,
        aggregate_combo=lambda *_args, **_kwargs: raw_scorch,
        annotate_scorch_z_scores=lambda *_args, **_kwargs: None,
        emit_task_event=lambda *_args, **_kwargs: None,
        emit_bench_event=lambda *_args, **_kwargs: None,
        summarize_chunk_partition=lambda *_args, **_kwargs: None,
        find_consensus_csv=lambda _dock_dir: consensus,
        rerank_consensus_with_scorch=lambda *args, **kwargs: calls.append(
            {"args": args, "kwargs": kwargs}
        )
        or True,
    )

    finalize_orchestration_combos(
        cfg=cfg,
        combos={combo},
        combos_with_tasks={combo},
        combo_modes={combo: {"fda"}},
        combo_failed_local={combo: False},
        specs=[spec],
        run_root=run_root,
        post_run_root=post_root,
        overwrite=False,
        decoy_prefix="dud",
        component="[test]",
        logger=logging.getLogger("test.scorch.finalize"),
        deps=deps,
        total_jobs=1,
        completed=1,
        failed_jobs=0,
        combos_attempted=1,
        skipped_missing_receptor=0,
        skipped_missing_consensus=0,
    )

    assert len(calls) == 1
    assert calls[0]["kwargs"]["overwrite"] is True


def test_finalize_does_not_rerank_with_incomplete_plan_bound_decoy_mode(
    tmp_path: Path, monkeypatch
) -> None:
    import post_docking.rescoring.scorch_orchestration_finalize as finalize_mod

    cfg = {"OVERALL_DIR": str(tmp_path)}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "", "")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    consensus = run_root / "PDB1" / "consensus_docking_scores.csv"
    raw_fda = post_root / "PDB1" / "scorch_scores_all.csv"
    stale_dud = post_root / "PDB1" / "dud_scorch_scores_all.csv"
    _write_csv(
        consensus,
        [
            {
                "ligand": "lig_a.pdbqt",
                "run_mode": "fda",
                "consensus_score": "1.0",
                "best_engine": "vina",
            }
        ],
    )
    _write_csv(
        raw_fda,
        [
            {
                "Ligand_ID": "lig_a_stage1",
                "run_mode": "fda",
                "source": "vina",
                "SCORCH_score": "0.8",
            }
        ],
    )
    _write_csv(
        stale_dud,
        [
            {
                "Ligand_ID": "stale_decoy_stage1",
                "run_mode": "dud",
                "source": "vina",
                "SCORCH_score": "0.1",
            }
        ],
    )

    def fake_planned_csvs(_cfg, _run_id, _combo, _specs, run_mode, _decoy_prefix):
        return [raw_fda] if run_mode == "fda" else []

    def fake_planned_bases(_cfg, _run_id, _combo, _specs, run_mode, _decoy_prefix):
        return {"lig_a"} if run_mode == "fda" else set()

    monkeypatch.setattr(finalize_mod, "planned_completed_output_csvs", fake_planned_csvs)
    monkeypatch.setattr(finalize_mod, "planned_completed_output_bases", fake_planned_bases)
    aggregate_calls: list[dict[str, object]] = []
    rerank_calls: list[dict[str, object]] = []
    deps = OrchestrationDeps(
        scheduler_runtime_snapshot=lambda *_args, **_kwargs: (8, 0),
        discover_mode_dirs=lambda *_args, **_kwargs: {"dud": [], "fda": []},
        load_control_bases=lambda *_args, **_kwargs: set(),
        stage_dir_candidates=lambda *_args, **_kwargs: [],
        score_csv_for_spec=lambda *_args, **_kwargs: (None, [], False),
        select_top_bases_from_score_csv=lambda *_args, **_kwargs: SelectionResult(
            allowed_bases=set(),
            n_pool=0,
            k=0,
            controls_total=0,
            selected_stage_by_base={},
            selected_score_by_base={},
        ),
        adaptive_chunk_size=lambda **_kwargs: 1,
        chunk_allowed_bases=lambda bases, _size: [set(bases)],
        tail_split_allowed_chunks=lambda chunks, **_kwargs: list(chunks),
        score_stage=lambda *_args, **_kwargs: (True, raw_fda),
        task_estimate_seconds=lambda *_args, **_kwargs: 1.0,
        task_allowed_count=lambda task, *_args, **_kwargs: len(task[3]),
        split_task_for_rechunk=lambda task, *_args, **_kwargs: [task],
        elastic_scorch_threads=lambda **_kwargs: 1,
        aggregate_combo=lambda *_args, **kwargs: aggregate_calls.append(kwargs)
        or (raw_fda if kwargs.get("run_mode") == "fda" else None),
        annotate_scorch_z_scores=lambda *_args, **_kwargs: None,
        emit_task_event=lambda *_args, **_kwargs: None,
        emit_bench_event=lambda *_args, **_kwargs: None,
        summarize_chunk_partition=lambda *_args, **_kwargs: None,
        find_consensus_csv=lambda _dock_dir: consensus,
        rerank_consensus_with_scorch=lambda *args, **kwargs: rerank_calls.append(
            {"args": args, "kwargs": kwargs}
        )
        or True,
    )

    completed = finalize_orchestration_combos(
        cfg=cfg,
        combos={combo},
        combos_with_tasks={combo},
        combo_modes={combo: {"fda", "dud"}},
        combo_failed_local={combo: False},
        specs=[spec],
        run_root=run_root,
        post_run_root=post_root,
        overwrite=False,
        decoy_prefix="dud",
        component="[test]",
        logger=logging.getLogger("test.scorch.finalize.incomplete_plan"),
        deps=deps,
        total_jobs=2,
        completed=1,
        failed_jobs=0,
        combos_attempted=1,
        skipped_missing_receptor=0,
        skipped_missing_consensus=0,
    )

    assert completed == set()
    assert rerank_calls == []
    assert any(call.get("run_mode") == "dud" for call in aggregate_calls)
    assert any(
        call.get("run_mode") == "fda"
        and call.get("allowed_ligand_bases") == {"lig_a"}
        for call in aggregate_calls
    )


def test_finalize_reranks_dud_only_plan_bound_mode(
    tmp_path: Path, monkeypatch
) -> None:
    import post_docking.rescoring.scorch_orchestration_finalize as finalize_mod

    cfg = {"OVERALL_DIR": str(tmp_path)}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "", "")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    consensus = run_root / "PDB1" / "consensus_docking_scores.csv"
    raw_dud = post_root / "PDB1" / "dud_scorch_scores_all.csv"
    _write_csv(
        consensus,
        [
            {
                "ligand": "decoys_final_00001.pdbqt",
                "run_mode": "dud",
                "consensus_score": "1.0",
                "best_engine": "vina",
            }
        ],
    )
    _write_csv(
        raw_dud,
        [
            {
                "Ligand_ID": "decoys_final_00001_stage1",
                "run_mode": "dud",
                "source": "vina",
                "SCORCH_score": "0.8",
                "SCORCH_certainty": "0.9",
            }
        ],
    )

    def fake_planned_csvs(_cfg, _run_id, _combo, _specs, run_mode, _decoy_prefix):
        return [raw_dud] if run_mode == "dud" else []

    def fake_planned_bases(_cfg, _run_id, _combo, _specs, run_mode, _decoy_prefix):
        return {"decoys_final_00001"} if run_mode == "dud" else set()

    monkeypatch.setattr(finalize_mod, "planned_completed_output_csvs", fake_planned_csvs)
    monkeypatch.setattr(finalize_mod, "planned_completed_output_bases", fake_planned_bases)
    aggregate_calls: list[dict[str, object]] = []
    rerank_calls: list[dict[str, object]] = []

    def fake_rerank(*args, **kwargs):
        rerank_calls.append({"args": args, "kwargs": kwargs})
        out_csv = args[2]
        row = {
            "ligand": "decoys_final_00001.pdbqt",
            "run_mode": "dud",
            "rescored_flag": "1",
            "consensus_score": "1.0",
            "SCORCH_score_used": "0.8",
        }
        _write_csv(out_csv, [row])
        _write_csv(out_csv.with_name("dud_consensus_reranked_scorch.csv"), [row])
        return True

    deps = OrchestrationDeps(
        scheduler_runtime_snapshot=lambda *_args, **_kwargs: (8, 0),
        discover_mode_dirs=lambda *_args, **_kwargs: {"dud": [], "fda": []},
        load_control_bases=lambda *_args, **_kwargs: set(),
        stage_dir_candidates=lambda *_args, **_kwargs: [],
        score_csv_for_spec=lambda *_args, **_kwargs: (None, [], False),
        select_top_bases_from_score_csv=lambda *_args, **_kwargs: SelectionResult(
            allowed_bases=set(),
            n_pool=0,
            k=0,
            controls_total=0,
            selected_stage_by_base={},
            selected_score_by_base={},
        ),
        adaptive_chunk_size=lambda **_kwargs: 1,
        chunk_allowed_bases=lambda bases, _size: [set(bases)],
        tail_split_allowed_chunks=lambda chunks, **_kwargs: list(chunks),
        score_stage=lambda *_args, **_kwargs: (True, raw_dud),
        task_estimate_seconds=lambda *_args, **_kwargs: 1.0,
        task_allowed_count=lambda task, *_args, **_kwargs: len(task[3]),
        split_task_for_rechunk=lambda task, *_args, **_kwargs: [task],
        elastic_scorch_threads=lambda **_kwargs: 1,
        aggregate_combo=lambda *_args, **kwargs: aggregate_calls.append(kwargs)
        or (raw_dud if kwargs.get("run_mode") == "dud" else None),
        annotate_scorch_z_scores=lambda *_args, **_kwargs: None,
        emit_task_event=lambda *_args, **_kwargs: None,
        emit_bench_event=lambda *_args, **_kwargs: None,
        summarize_chunk_partition=lambda *_args, **_kwargs: None,
        find_consensus_csv=lambda _dock_dir: consensus,
        rerank_consensus_with_scorch=fake_rerank,
    )

    completed = finalize_orchestration_combos(
        cfg=cfg,
        combos={combo},
        combos_with_tasks={combo},
        combo_modes={combo: {"dud"}},
        combo_failed_local={combo: False},
        specs=[spec],
        run_root=run_root,
        post_run_root=post_root,
        overwrite=False,
        decoy_prefix="dud",
        component="[test]",
        logger=logging.getLogger("test.scorch.finalize.dud_only"),
        deps=deps,
        total_jobs=1,
        completed=1,
        failed_jobs=0,
        combos_attempted=1,
        skipped_missing_receptor=0,
        skipped_missing_consensus=0,
    )

    assert completed == {combo}
    assert len(rerank_calls) == 1
    assert list(rerank_calls[0]["args"][1]) == [raw_dud]
    assert any(call.get("run_mode") == "dud" for call in aggregate_calls)
    assert (post_root / "PDB1" / "consensus_reranked_scorch.csv").exists()
    assert (post_root / "PDB1" / "dud_consensus_reranked_scorch.csv").exists()
