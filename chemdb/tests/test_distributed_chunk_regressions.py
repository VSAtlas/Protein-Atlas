from __future__ import annotations

import csv
import json
import logging
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import main
from cli import distributed_chunk_runtime
from cli import planner_chunking
from cli import distributed_chunk_planner
from cli.distributed_chunk_runtime_helpers import (
    build_scope_maps,
    collect_combo_chunk_progress,
    collect_terminal_chunk_failures,
)
from cli.distributed_runtime_completion import DeferredScorchCompletionBookkeeper
from cli.chunk_tail_optimizer import build_candidate_ids, optimize_chunks_for_tail
from cli.run_process_one import (
    ProcessOneContext,
    ProcessOneDeps,
    ProcessOneSharedState,
    build_process_one_runner,
)
from docking.docking_runtime_state import _signal_distributed_prep_ready
from docking.docking_vina import write_scores_csv
from docking.pose_validation import validate_pose_pdbqt
from cli.run_bootstrap_helpers import init_config_run_dir
from config.output_paths import output_root


def _write_summary_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["run_id", "Ligand", "stage1"])
        writer.writeheader()
        writer.writerows(rows)


def _base_cfg(tmp_path: Path, *, chunk_keys: list[str]) -> dict[str, object]:
    cfg: dict[str, object] = {
        "OVERALL_DIR": str(tmp_path),
        "INPUT_DIR": str(tmp_path / "input_pdbs"),
        "OUTPUT_DIR": str(tmp_path / "processed_pdbs"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "PREPPED_LIGANDS_DIR": str(tmp_path / "prepped_ligands"),
        "RUN_ID": "regress_run",
        "_CHUNK_LIGAND_KEYS": list(chunk_keys),
    }
    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR", "PREPPED_LIGANDS_DIR"):
        Path(str(cfg[key])).mkdir(parents=True, exist_ok=True)
    return cfg


def _score_history(ligand_name: str, score: float) -> dict[str, dict[str, dict[str, object]]]:
    return {
        "stage1": {
            ligand_name: {
                "score": score,
                "valid": True,
                "reason": "",
                "heavy_atoms": 12,
                "le": None,
                "pains_flag": False,
            }
        }
    }


def test_verify_chunk_outputs_handles_run_scoped_docked_dir(tmp_path: Path) -> None:
    run_id = "rid1"
    pdb_id = "T123"
    docked_run_root = tmp_path / "docked" / run_id
    summary_path = docked_run_root / pdb_id / "HOLO" / "docking_score_summary.csv"
    _write_summary_csv(
        summary_path,
        [
            {
                "run_id": run_id,
                "Ligand": "lig_a.pdbqt",
                "stage1": "-7.10",
            }
        ],
    )

    ok, reason = main._verify_chunk_combo_outputs(
        {"DOCKED_DIR": str(docked_run_root)},
        run_id=run_id,
        pdb_id=pdb_id,
        variant_label="HOLO",
        ph_tag=None,
        library_name="",
        chunk_ligand_bases=["lig_a"],
    )
    assert ok, reason


def test_distributed_manifest_dir_honors_relocated_manifest_root(tmp_path: Path) -> None:
    manifest_root = tmp_path / "scratch" / "outputs" / "manifests"
    cfg = {
        "OVERALL_DIR": str(tmp_path / "work_repo"),
        "MANIFESTS_DIR": str(manifest_root),
    }

    assert distributed_chunk_planner._distributed_manifest_dir(
        cfg, run_id="reloc_chunks"
    ) == (manifest_root / "reloc_chunks" / "distributed")


def test_verify_chunk_outputs_degrades_when_summary_missing(tmp_path: Path) -> None:
    ok, reason = main._verify_chunk_combo_outputs(
        {"DOCKED_DIR": str(tmp_path / "docked")},
        run_id="rid1",
        pdb_id="T123",
        variant_label="HOLO",
        ph_tag=None,
        library_name="",
        chunk_ligand_bases=["lig_a"],
    )
    assert ok is True
    assert reason.startswith("degraded_missing_docking_summary:")


def test_verify_chunk_outputs_degrades_when_summary_has_nul(tmp_path: Path) -> None:
    run_id = "rid2"
    pdb_id = "T999"
    summary_path = tmp_path / "docked" / run_id / pdb_id / "HOLO" / "docking_score_summary.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_bytes(b"run_id,Ligand,stage1\nrid2,lig_a.pdbqt,\x00\n")

    ok, reason = main._verify_chunk_combo_outputs(
        {"DOCKED_DIR": str(tmp_path / "docked")},
        run_id=run_id,
        pdb_id=pdb_id,
        variant_label="HOLO",
        ph_tag=None,
        library_name="",
        chunk_ligand_bases=["lig_a"],
    )
    assert ok is True
    assert reason.startswith("degraded_missing_chunk_scores:")


def test_verify_chunk_outputs_chunk_mode_requires_stage_marker(tmp_path: Path) -> None:
    ok, reason = main._verify_chunk_combo_outputs(
        {"DOCKED_DIR": str(tmp_path / "docked")},
        run_id="rid_missing_marker",
        chunk_id="chunk_a",
        pdb_id="T321",
        variant_label="HOLO",
        ph_tag="pH7_0",
        library_name="",
        chunk_ligand_bases=["lig_a"],
    )
    assert ok is False
    assert reason.startswith("retryable_missing_stage1_chunk_marker:")
    assert ":missing_docking_summary:" in reason


def test_verify_chunk_outputs_chunk_mode_missing_all_scores_is_fatal(
    tmp_path: Path,
) -> None:
    run_id = "rid_missing_scores"
    pdb_id = "T654"
    chunk_id = "chunk_b"
    stage1_dir = tmp_path / "docked" / run_id / pdb_id / "HOLO" / "pH7_0" / "stage1"
    stage1_dir.mkdir(parents=True, exist_ok=True)
    (stage1_dir / f"completion_vina__chunk_{chunk_id}.json").write_text(
        "{\"ok\": true}\n", encoding="utf-8"
    )
    summary_path = (
        tmp_path
        / "docked"
        / run_id
        / pdb_id
        / "HOLO"
        / "pH7_0"
        / "docking_score_summary.csv"
    )
    _write_summary_csv(
        summary_path,
        [
            {
                "run_id": run_id,
                "Ligand": "lig_other.pdbqt",
                "stage1": "-6.10",
            }
        ],
    )

    ok, reason = main._verify_chunk_combo_outputs(
        {"DOCKED_DIR": str(tmp_path / "docked")},
        run_id=run_id,
        chunk_id=chunk_id,
        pdb_id=pdb_id,
        variant_label="HOLO",
        ph_tag="pH7_0",
        library_name="",
        chunk_ligand_bases=["lig_a"],
    )
    assert ok is False
    assert reason.startswith("missing_chunk_scores:missing=1:expected=1:scored=1")


def test_verify_chunk_outputs_chunk_mode_partial_scores_are_retryable(
    tmp_path: Path,
) -> None:
    run_id = "rid_partial_scores"
    pdb_id = "T655"
    chunk_id = "chunk_c"
    stage1_dir = tmp_path / "docked" / run_id / pdb_id / "HOLO" / "pH7_0" / "stage1"
    stage1_dir.mkdir(parents=True, exist_ok=True)
    (stage1_dir / f"completion_vina__chunk_{chunk_id}.json").write_text(
        "{\"ok\": true}\n", encoding="utf-8"
    )
    summary_path = (
        tmp_path
        / "docked"
        / run_id
        / pdb_id
        / "HOLO"
        / "pH7_0"
        / "docking_score_summary.csv"
    )
    _write_summary_csv(
        summary_path,
        [
            {
                "run_id": run_id,
                "Ligand": "lig_a.pdbqt",
                "stage1": "-6.10",
            }
        ],
    )

    ok, reason = main._verify_chunk_combo_outputs(
        {"DOCKED_DIR": str(tmp_path / "docked")},
        run_id=run_id,
        chunk_id=chunk_id,
        pdb_id=pdb_id,
        variant_label="HOLO",
        ph_tag="pH7_0",
        library_name="",
        chunk_ligand_bases=["lig_a", "lig_b"],
    )
    assert ok is False
    assert reason.startswith("missing_chunk_scores:missing=1:expected=2:scored=1")


def test_chunk_mode_retryable_verification_failure_is_chunk_local(tmp_path: Path) -> None:
    failed_root = tmp_path / "failed"
    failed_root.mkdir(parents=True, exist_ok=True)
    cfg: dict[str, object] = {"CPU": 1}
    context = ProcessOneContext(
        run_id="rid_retryable",
        label="HOLO",
        mode="distributed",
        variant="holo",
        is_resume=False,
        completed_combo_lookup=set(),
        dist_combo_chunk_mode=True,
        cfg_v=cfg,
        stages=None,
        params=None,
        tokens=[],
        global_start=0.0,
    )
    shared_state = ProcessOneSharedState(
        failed_root=str(failed_root),
        failed_entries=[],
        run_scope_completion_times={},
        retention_lock=threading.Lock(),
        pending_retention_combos=set(),
        pending_coverage_refresh={},
        pending_retention_lock=threading.Lock(),
        scorch_queue_service=None,
    )
    deps = ProcessOneDeps(
        normalize_ph_tag_token=lambda ph: str(ph or "").strip() or "base",
        chunk_ligand_key=lambda name: str(name),
        process_one_protein=lambda _cfg, _pdb_file, _stages, _params: None,
        update_manifest_for_protein_start=lambda *_args, **_kwargs: None,
        update_manifest_for_protein_success=lambda *_args, **_kwargs: None,
        update_manifest_for_protein_failure=lambda *_args, **_kwargs: None,
        verify_chunk_combo_outputs=lambda *_args, **_kwargs: (
            False,
            "retryable_missing_stage1_chunk_marker:/m:missing_docking_summary:/s",
        ),
        resolve_combo_output_dir=lambda *_args, **_kwargs: tmp_path,
        resolve_combo_docking_summary_csv=lambda *_args, **_kwargs: tmp_path / "missing.csv",
        load_scored_ligand_keys_from_summary=lambda *_args, **_kwargs: set(),
        update_combo_coverage_snapshot=lambda *_args, **_kwargs: None,
        resolve_combo_post_consensus_csv=lambda *_args, **_kwargs: tmp_path / "missing_post.csv",
        load_post_scored_ligand_keys=lambda *_args, **_kwargs: set(),
        acquire_global_slot=lambda *_args, **_kwargs: threading.Lock(),
        maybe_run_scorch_rescore_for_pdb=lambda *_args, **_kwargs: None,
        maybe_run_artifact_retention_for_combo=lambda *_args, **_kwargs: None,
    )
    process_one = build_process_one_runner(
        context=context,
        shared_state=shared_state,
        deps=deps,
    )

    ok = process_one(
        "BOJG.pdb",
        cfg,
        "pH7_2",
        chunk_ligand_bases=["lig_a"],
        chunk_id="c1",
        chunk_library_name="bench2",
    )
    assert ok is False
    assert shared_state.failed_entries == []
    assert cfg["_LAST_CHUNK_FAILURE_RETRYABLE"] is True
    assert str(cfg["_LAST_CHUNK_FAILURE_REASON"]).startswith(
        "chunk_output_verification_failed:"
    )
    assert list(failed_root.glob("*.log")) == []


def test_collect_terminal_chunk_failures_respects_retry_window() -> None:
    by_chunk_id = {
        "chunk_a": {
            "pdb_id": "BOJG",
            "pdb_file": "BOJG.pdb",
            "ph_tag": "pH7_2",
        }
    }
    recorded_terminal_chunk_ids: set[str] = set()
    failed_entries: list[tuple[str, str, str, str, str]] = []
    dist_ctx = type("Ctx", (), {"task_id": 1})()

    collect_terminal_chunk_failures(
        by_chunk_id=by_chunk_id,
        label="HOLO",
        dist_ctx=dist_ctx,
        max_attempts=3,
        read_chunk_result=lambda *_args, **_kwargs: {
            "status": "failed",
            "attempt": 1,
            "task_id": 1,
            "chunk_reason": "chunk_output_verification_failed:retryable_missing_stage1_chunk_marker:/m:missing_docking_summary:/s",
        },
        cfg_v={},
        recorded_terminal_chunk_ids=recorded_terminal_chunk_ids,
        failed_entries=failed_entries,
    )
    assert failed_entries == []

    collect_terminal_chunk_failures(
        by_chunk_id=by_chunk_id,
        label="HOLO",
        dist_ctx=dist_ctx,
        max_attempts=3,
        read_chunk_result=lambda *_args, **_kwargs: {
            "status": "failed",
            "attempt": 3,
            "task_id": 1,
            "chunk_reason": "chunk_output_verification_failed:retryable_missing_stage1_chunk_marker:/m:missing_docking_summary:/s",
        },
        cfg_v={},
        recorded_terminal_chunk_ids=recorded_terminal_chunk_ids,
        failed_entries=failed_entries,
    )
    assert len(failed_entries) == 1
    assert failed_entries[0][3] == "ChunkTerminalFailed"


def test_collect_combo_chunk_progress_counts_global_terminal_by_scope() -> None:
    by_chunk_id = {
        "a1": {"pdb_id": "BNNQ", "pdb_file": "BNNQ.pdb", "ph_tag": "pH5_6"},
        "a2": {"pdb_id": "BNNQ", "pdb_file": "BNNQ.pdb", "ph_tag": "pH5_6"},
        "b1": {"pdb_id": "BOJG", "pdb_file": "BOJG.pdb", "ph_tag": "pH7_2"},
    }
    results = {
        "a1": {"status": "completed", "attempt": 1},
        "a2": {"status": "failed", "attempt": 3},
        "b1": {"status": "failed", "attempt": 1},
    }

    progress = collect_combo_chunk_progress(
        by_chunk_id=by_chunk_id,
        label="HOLO",
        dist_ctx=object(),
        max_attempts=3,
        read_chunk_result=lambda *_args, **kwargs: results.get(str(kwargs["chunk_id"])),
        cfg_v={},
    )

    assert progress[("BNNQ", "HOLO", "pH5_6")] == (2, 1, 1)
    assert ("BOJG", "HOLO", "pH7_2") not in progress


def test_deferred_scorch_reconcile_queues_owned_scope_from_global_progress() -> None:
    scope_owned = ("BNNQ", "HOLO", "pH5_6")
    scope_other = ("BOJG", "HOLO", "pH7_2")
    submitted: list[tuple[str, str | None, str | None, int | None]] = []
    aggregated: list[tuple[str, str, str]] = []

    class Queue:
        def submit(
            self,
            *,
            cfg,
            pdb_id,
            variant=None,
            ph=None,
            allowed_count_hint=None,
        ) -> bool:
            submitted.append((pdb_id, variant, ph, allowed_count_hint))
            return True

    bookkeeper = DeferredScorchCompletionBookkeeper(
        cfg_v={},
        run_id="rid",
        label="HOLO",
        dist_task_id=0,
        combo_total_chunks={scope_owned: 2, scope_other: 2},
        combo_total_ligands={scope_owned: 32, scope_other: 32},
        combo_owner_task_id={scope_owned: 0, scope_other: 1},
        scorch_queue_service=Queue(),
        variant_param_for_scorch="HOLO",
        hybrid_start_ts=0.0,
        acquire_global_slot=lambda *_args, **_kwargs: threading.Lock(),
        maybe_run_scorch_rescore_for_pdb=lambda *_args, **_kwargs: 0,
        normalize_ph_tag_token=lambda raw: str(raw or "base"),
        combo_progress_reader=lambda: {
            scope_owned: (2, 2),
            scope_other: (2, 2),
        },
        before_queue_combo_scorch=lambda scope: aggregated.append(scope),
        global_reconcile_interval_sec=0.0,
    )

    assert (
        bookkeeper.reconcile_ready_owned_deferred_submissions(
            reason="test",
            force=True,
        )
        == 1
    )
    assert submitted == [("BNNQ", "HOLO", "pH5_6", 32)]
    assert aggregated == [scope_owned]
    assert scope_owned in bookkeeper.combo_scorch_submitted

    assert (
        bookkeeper.reconcile_ready_owned_deferred_submissions(
            reason="test_second",
            force=True,
        )
        == 0
    )
    assert len(submitted) == 1


def test_deferred_scorch_reconcile_blocks_final_when_chunks_failed() -> None:
    scope_owned = ("BNNQ", "HOLO", "pH5_6")
    submitted: list[tuple[str, str | None, str | None, int | None]] = []

    class Queue:
        def submit(
            self,
            *,
            cfg,
            pdb_id,
            variant=None,
            ph=None,
            allowed_count_hint=None,
        ) -> bool:
            submitted.append((pdb_id, variant, ph, allowed_count_hint))
            return True

    bookkeeper = DeferredScorchCompletionBookkeeper(
        cfg_v={},
        run_id="rid",
        label="HOLO",
        dist_task_id=0,
        combo_total_chunks={scope_owned: 2},
        combo_total_ligands={scope_owned: 32},
        combo_owner_task_id={scope_owned: 0},
        scorch_queue_service=Queue(),
        variant_param_for_scorch="HOLO",
        hybrid_start_ts=0.0,
        acquire_global_slot=lambda *_args, **_kwargs: threading.Lock(),
        maybe_run_scorch_rescore_for_pdb=lambda *_args, **_kwargs: 0,
        normalize_ph_tag_token=lambda raw: str(raw or "base"),
        combo_progress_reader=lambda: {scope_owned: (2, 1, 1)},
        global_reconcile_interval_sec=0.0,
    )

    assert (
        bookkeeper.reconcile_ready_owned_deferred_submissions(
            reason="failed_scope",
            force=True,
        )
        == 0
    )
    assert submitted == []
    assert scope_owned not in bookkeeper.combo_scorch_submitted
    assert bookkeeper.combo_failed_chunks_live[scope_owned] == 1


def test_deferred_scorch_live_terminal_failure_does_not_queue_final() -> None:
    scope_owned = ("BNNQ", "HOLO", "pH5_6")
    submitted: list[tuple[str, str | None, str | None, int | None]] = []

    class Queue:
        def submit(
            self,
            *,
            cfg,
            pdb_id,
            variant=None,
            ph=None,
            allowed_count_hint=None,
        ) -> bool:
            submitted.append((pdb_id, variant, ph, allowed_count_hint))
            return True

    bookkeeper = DeferredScorchCompletionBookkeeper(
        cfg_v={},
        run_id="rid",
        label="HOLO",
        dist_task_id=0,
        combo_total_chunks={scope_owned: 2},
        combo_total_ligands={scope_owned: 32},
        combo_owner_task_id={scope_owned: 0},
        scorch_queue_service=Queue(),
        variant_param_for_scorch="HOLO",
        hybrid_start_ts=0.0,
        acquire_global_slot=lambda *_args, **_kwargs: threading.Lock(),
        maybe_run_scorch_rescore_for_pdb=lambda *_args, **_kwargs: 0,
        normalize_ph_tag_token=lambda raw: str(raw or "base"),
    )

    bookkeeper.mark_chunk_completed(scope_owned, reason="ok")
    bookkeeper.mark_chunk_terminal_failure(scope_owned, reason="terminal_failed")

    assert submitted == []
    assert scope_owned not in bookkeeper.combo_scorch_submitted
    assert bookkeeper.combo_completed_chunks_live[scope_owned] == 1
    assert bookkeeper.combo_failed_chunks_live[scope_owned] == 1


def test_deferred_scorch_shard_mode_keeps_global_progress_owner_scoped() -> None:
    scope_owned = ("BNNQ", "HOLO", "pH5_6")
    scope_other = ("BOJG", "HOLO", "pH7_2")
    submitted: list[tuple[str, str | None, str | None, int | None]] = []

    class Queue:
        def submit(
            self,
            *,
            cfg,
            pdb_id,
            variant=None,
            ph=None,
            allowed_count_hint=None,
        ) -> bool:
            submitted.append((pdb_id, variant, ph, allowed_count_hint))
            return True

    bookkeeper = DeferredScorchCompletionBookkeeper(
        cfg_v={"SCORCH_SHARDS_ENABLE": True},
        run_id="rid",
        label="HOLO",
        dist_task_id=0,
        combo_total_chunks={scope_owned: 2, scope_other: 2},
        combo_total_ligands={scope_owned: 32, scope_other: 32},
        combo_owner_task_id={scope_owned: 0, scope_other: 1},
        scorch_queue_service=Queue(),
        variant_param_for_scorch="HOLO",
        hybrid_start_ts=0.0,
        acquire_global_slot=lambda *_args, **_kwargs: threading.Lock(),
        maybe_run_scorch_rescore_for_pdb=lambda *_args, **_kwargs: 0,
        normalize_ph_tag_token=lambda raw: str(raw or "base"),
        combo_progress_reader=lambda: {
            scope_owned: (2, 2),
            scope_other: (2, 2),
        },
        global_reconcile_interval_sec=0.0,
    )

    assert (
        bookkeeper.reconcile_ready_owned_deferred_submissions(
            reason="test",
            force=True,
        )
        == 1
    )
    assert submitted == [("BNNQ", "HOLO", "pH5_6", 32)]
    assert scope_owned in bookkeeper.combo_scorch_submitted
    assert scope_other not in bookkeeper.combo_scorch_submitted


def test_deferred_scorch_shard_mode_final_flush_keeps_owner_scoped() -> None:
    scope_owned = ("BNNQ", "HOLO", "pH5_6")
    scope_other = ("BOJG", "HOLO", "pH7_2")
    submitted: list[tuple[str, str | None, str | None, int | None]] = []

    class Queue:
        def submit(
            self,
            *,
            cfg,
            pdb_id,
            variant=None,
            ph=None,
            allowed_count_hint=None,
        ) -> bool:
            submitted.append((pdb_id, variant, ph, allowed_count_hint))
            return True

    bookkeeper = DeferredScorchCompletionBookkeeper(
        cfg_v={"SCORCH_SHARDS_ENABLE": True},
        run_id="rid",
        label="HOLO",
        dist_task_id=0,
        combo_total_chunks={scope_owned: 2, scope_other: 2},
        combo_total_ligands={scope_owned: 32, scope_other: 32},
        combo_owner_task_id={scope_owned: 0, scope_other: 1},
        scorch_queue_service=Queue(),
        variant_param_for_scorch="HOLO",
        hybrid_start_ts=0.0,
        acquire_global_slot=lambda *_args, **_kwargs: threading.Lock(),
        maybe_run_scorch_rescore_for_pdb=lambda *_args, **_kwargs: 0,
        normalize_ph_tag_token=lambda raw: str(raw or "base"),
    )

    bookkeeper.flush_owned_deferred_submissions(
        combo_completed_chunks={scope_owned: 2, scope_other: 2}
    )

    assert submitted == [("BNNQ", "HOLO", "pH5_6", 32)]
    assert scope_owned in bookkeeper.combo_scorch_submitted
    assert scope_other not in bookkeeper.combo_scorch_submitted


def test_deferred_scorch_final_flush_requires_full_clean_completion() -> None:
    scope_partial = ("BNNQ", "HOLO", "pH5_6")
    scope_failed = ("BOJG", "HOLO", "pH7_2")
    scope_clean = ("BNJS", "HOLO", "pH7_0")
    submitted: list[tuple[str, str | None, str | None, int | None]] = []

    class Queue:
        def submit(
            self,
            *,
            cfg,
            pdb_id,
            variant=None,
            ph=None,
            allowed_count_hint=None,
        ) -> bool:
            submitted.append((pdb_id, variant, ph, allowed_count_hint))
            return True

    bookkeeper = DeferredScorchCompletionBookkeeper(
        cfg_v={"SCORCH_SHARDS_ENABLE": True},
        run_id="rid",
        label="HOLO",
        dist_task_id=0,
        combo_total_chunks={scope_partial: 2, scope_failed: 2, scope_clean: 2},
        combo_total_ligands={scope_partial: 32, scope_failed: 32, scope_clean: 32},
        combo_owner_task_id={scope_partial: 0, scope_failed: 0, scope_clean: 0},
        scorch_queue_service=Queue(),
        variant_param_for_scorch="HOLO",
        hybrid_start_ts=0.0,
        acquire_global_slot=lambda *_args, **_kwargs: threading.Lock(),
        maybe_run_scorch_rescore_for_pdb=lambda *_args, **_kwargs: 0,
        normalize_ph_tag_token=lambda raw: str(raw or "base"),
    )

    bookkeeper.flush_owned_deferred_submissions(
        combo_completed_chunks={
            scope_partial: 1,
            scope_failed: 1,
            scope_clean: 2,
        },
        combo_failed_chunks={scope_failed: 1},
    )

    assert submitted == [("BNJS", "HOLO", "pH7_0", 32)]
    assert scope_clean in bookkeeper.combo_scorch_submitted
    assert scope_partial not in bookkeeper.combo_scorch_submitted
    assert scope_failed not in bookkeeper.combo_scorch_submitted


def test_completed_distributed_resume_queues_scorch_on_every_task(monkeypatch) -> None:
    pdb_ids = [f"P{i:03d}" for i in range(90)]
    chunk_plan = [
        {
            "chunk_id": f"chunk_{pdb_id}",
            "pdb_id": pdb_id,
            "pdb_file": f"{pdb_id}.pdb",
            "ph_tag": None,
            "ligand_bases": [f"{pdb_id}_lig_{j}" for j in range(10)],
            "assigned_task_id": i % 10,
            "weight": 1.0,
        }
        for i, pdb_id in enumerate(pdb_ids)
    ]
    cfg = {
        "CPU": 112,
        "USE_SCORCH": True,
        "SCORCH_SHARDS_ENABLE": True,
        "_RESUME_MANIFEST_PROTEINS": {
            f"{pdb_id}|LEGACY|base": {
                "status": "running",
                "stages": {
                    "distributed_chunk_coverage": {
                        "status": "completed",
                        "details": {
                            "planned_chunks": 1,
                            "completed_chunks": 1,
                            "failed_chunks": 0,
                            "missing_result_chunks": 0,
                            "planned_ligands": 10,
                            "completed_ligands": 10,
                            "has_consensus_csv": True,
                            "require_scorch": True,
                        },
                    },
                    "postprocessing": {
                        "status": "running",
                        "details": {"has_scorch_reranked_csv": False},
                    },
                },
            }
            for pdb_id in pdb_ids
        },
    }
    submissions: list[tuple[int, str, int | None]] = []

    class DistCtx:
        enabled = True
        task_count = 10

        def __init__(self, task_id: int) -> None:
            self.task_id = task_id

        def expected_task_ids(self):
            return list(range(10))

    class Queue:
        def __init__(self, task_id: int) -> None:
            self.task_id = task_id

        def submit(
            self,
            *,
            cfg,
            pdb_id,
            variant=None,
            ph=None,
            allowed_count_hint=None,
        ) -> bool:
            submissions.append((self.task_id, str(pdb_id), allowed_count_hint))
            return True

    class Bar:
        def reset(self, total=None) -> None:
            self.total = total

    mode_decision = type(
        "Decision",
        (),
        {
            "completion_lane_fraction": 0.30,
            "rebalance_sec": 60,
            "rebalance_chunks": 16,
        },
    )()

    monkeypatch.setattr(
        distributed_chunk_runtime,
        "_ensure_distributed_manifest_preflight",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        distributed_chunk_runtime,
        "_load_or_create_variant_chunk_plan",
        lambda *_args, **_kwargs: list(chunk_plan),
    )
    monkeypatch.setattr(
        distributed_chunk_runtime,
        "run_chunk_claim_loop",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("completed resume should not enter the Vina chunk loop")
        ),
    )

    assigned_ids: set[str] = set()
    for task_id in range(10):
        distributed_chunk_runtime.run_distributed_combo_variant(
            cfg_v=dict(cfg),
            dist_ctx=DistCtx(task_id),
            run_id="rid_resume",
            label="LEGACY",
            combo_items=[(f"{pdb_id}.pdb", None) for pdb_id in pdb_ids],
            tokens=["dud", "fda"],
            bar=Bar(),
            variant=None,
            dist_combo_hybrid_mode=True,
            execution_mode_decision=mode_decision,
            scorch_queue_service=Queue(task_id),
            process_one=lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("completed resume should not rerun docking chunks")
            ),
            recorded_terminal_chunk_ids=set(),
            failed_entries=[],
            run_scope_completion_times={},
            distributed_assigned_pdb_ids=assigned_ids,
            run_chunk_rebalance_count=0,
        )

    by_task = {task_id: 0 for task_id in range(10)}
    for task_id, _pdb_id, allowed_count_hint in submissions:
        by_task[task_id] += 1
        assert allowed_count_hint == 10
    assert by_task == {task_id: 9 for task_id in range(10)}
    assert len({pdb_id for _task_id, pdb_id, _hint in submissions}) == 90
    assert assigned_ids == set(pdb_ids)


def test_update_combo_coverage_snapshot_atomic_temp_is_unique(tmp_path: Path) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path)}
    run_id = "rid_cov"
    pdb_id = "ABCD"
    variant = "HOLO"
    ph = "pH7_0"

    def _write(idx: int) -> None:
        main._update_combo_coverage_snapshot(
            cfg,
            run_id=run_id,
            pdb_id=pdb_id,
            variant_label=variant,
            ph_tag=ph,
            library_name="libA",
            expected_keys={f"lig_{idx:03d}"},
        )

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(_write, range(40)))

    snapshot = (
        output_root(tmp_path, "manifests")
        / run_id
        / "coverage"
        / f"{pdb_id}__{variant}__{ph}.json"
    )
    assert snapshot.exists()
    payload = json.loads(snapshot.read_text(encoding="utf-8"))
    assert payload["pdb_id"] == pdb_id
    assert payload["variant"] == variant
    assert payload["ph"] == ph
    assert isinstance(payload.get("expected_ligands"), list)


def test_update_combo_coverage_snapshot_honors_relocated_manifests_dir(
    tmp_path: Path,
) -> None:
    manifest_root = tmp_path / "scratch" / "outputs" / "manifests"
    cfg = {"OVERALL_DIR": str(tmp_path), "MANIFESTS_DIR": str(manifest_root)}

    main._update_combo_coverage_snapshot(
        cfg,
        run_id="rid_reloc_cov",
        pdb_id="ABCD",
        variant_label="HOLO",
        ph_tag="pH7_0",
        library_name="libA",
        expected_keys={"lig_a"},
    )

    assert (
        manifest_root / "rid_reloc_cov" / "coverage" / "ABCD__HOLO__pH7_0.json"
    ).exists()


def test_chunk_vina_writes_merge_instead_of_overwrite(tmp_path: Path) -> None:
    cfg = _base_cfg(tmp_path, chunk_keys=["lig_a"])
    summary_path = Path(write_scores_csv(cfg, "TEST", _score_history("lig_a.pdbqt", -7.10)))

    cfg["_CHUNK_LIGAND_KEYS"] = ["lig_b"]
    write_scores_csv(cfg, "TEST", _score_history("lig_b.pdbqt", -8.20))

    with summary_path.open("r", newline="", encoding="utf-8") as handle:
        summary_rows = list(csv.DictReader(handle))
    assert {row["Ligand"] for row in summary_rows} == {"lig_a.pdbqt", "lig_b.pdbqt"}

    long_path = summary_path.with_name("docking_score_long.csv")
    with long_path.open("r", newline="", encoding="utf-8") as handle:
        long_rows = list(csv.DictReader(handle))
    assert {row["ligand"] for row in long_rows} == {"lig_a.pdbqt", "lig_b.pdbqt"}


def test_init_config_run_dir_ignores_concurrent_missing_rmtree(
    tmp_path: Path, monkeypatch
) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "CONFIGS_DIR": str(tmp_path / "configs"),
        "RUN_ID": "run_cfg",
        "RESET_CONFIGS": True,
    }
    run_dir = Path(str(cfg["CONFIGS_DIR"])) / str(cfg["RUN_ID"])
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "old.txt").write_text("old", encoding="utf-8")

    def _raise_missing(*_args, **_kwargs) -> None:
        raise FileNotFoundError("simulated concurrent delete")

    monkeypatch.setattr("cli.run_bootstrap_helpers.shutil.rmtree", _raise_missing)
    init_config_run_dir(cfg, run_id=cfg["RUN_ID"], reset=True, logger=None)
    assert Path(str(cfg["CONFIG_RUN_DIR"])).exists()


def test_validate_pose_empty_receptor_returns_invalid(tmp_path: Path) -> None:
    receptor = tmp_path / "empty_receptor.pdbqt"
    ligand = tmp_path / "ligand.pdbqt"
    receptor.write_text("REMARK EMPTY\n", encoding="utf-8")
    ligand.write_text(
        "ATOM      1  C   LIG A   1      10.000  10.000  10.000  0.00  0.00      C \n",
        encoding="utf-8",
    )

    result = validate_pose_pdbqt(
        protein_pdbqt=str(receptor),
        ligand_pdbqt=str(ligand),
        pocket_center=(10.0, 10.0, 10.0),
    )
    assert result["valid"] is False
    assert result["reason"] == "Protein has no heavy atoms"


def test_distributed_chunk_local_workers_cpu_budget() -> None:
    assert planner_chunking._distributed_chunk_local_workers({"CPU": 1}) == 1
    assert planner_chunking._distributed_chunk_local_workers({"CPU": 16}) == 16
    assert planner_chunking._distributed_chunk_local_workers({"CPU": 48}) == 48
    assert planner_chunking._distributed_chunk_local_workers({"CPU": 109}) == 109
    assert planner_chunking._distributed_chunk_local_workers({"CPU": 110}) == 110
    assert planner_chunking._distributed_chunk_local_workers({"CPU": 112}) == 112
    assert planner_chunking._distributed_chunk_local_workers({"CPU": 113}) == 113
    assert planner_chunking._distributed_chunk_local_workers({"CPU": 400}) == 400
    assert (
        planner_chunking._distributed_chunk_cpu_per_worker(
            {"CPU": 64},
            local_workers=64,
        )
        == 1
    )
    assert (
        planner_chunking._distributed_chunk_cpu_per_worker(
            {"CPU": 110},
            local_workers=110,
        )
        == 1
    )


class _FakeDistContext:
    task_id = 0

    def expected_task_ids(self) -> list[int]:
        return [0, 1, 2]


def test_scope_owner_assignment_balances_small_scope_sets() -> None:
    chunk_plan = []
    for pdb_id, ph_tag in [
        ("BNNQ", "pH5_6"),
        ("BNJS", "pH7_0"),
        ("BOJG", "pH7_2"),
        ("BOJG", "pH7_7"),
    ]:
        for idx in range(3):
            chunk_plan.append(
                {
                    "chunk_id": f"{pdb_id}_{ph_tag}_{idx}",
                    "pdb_file": f"{pdb_id}.pdb",
                    "pdb_id": pdb_id,
                    "ph_tag": ph_tag,
                    "ligand_bases": [f"lig_{idx}_a", f"lig_{idx}_b"],
                }
            )

    _, owner_map, total_chunks, _ = build_scope_maps(
        chunk_plan=chunk_plan,
        label="HOLO",
        dist_ctx=_FakeDistContext(),
        owner_task_resolver=lambda *args, **kwargs: 1,
    )

    counts = {task_id: list(owner_map.values()).count(task_id) for task_id in [0, 1, 2]}
    bojg_owners = {owner for scope, owner in owner_map.items() if scope[0] == "BOJG"}
    assert set(owner_map.values()) == {0, 1, 2}
    assert len(bojg_owners) == 1
    assert max(counts.values()) - min(counts.values()) <= 1
    assert sorted(total_chunks.values()) == [3, 3, 3, 3]


def test_split_ligands_into_chunks_scales_with_local_workers() -> None:
    ligands = [f"lig_{idx:05d}" for idx in range(2048)]
    chunks_low = planner_chunking._split_ligands_into_chunks(
        ligands, task_count=3, local_workers_per_task=1
    )
    chunks_high = planner_chunking._split_ligands_into_chunks(
        ligands, task_count=3, local_workers_per_task=8
    )
    assert len(chunks_high) >= len(chunks_low)
    assert len(chunks_high) > len(chunks_low)


def test_split_ligands_into_chunks_adapts_min_size_for_large_workers(
    monkeypatch,
) -> None:
    monkeypatch.delenv("ATLAS_CHUNK_MIN_SIZE", raising=False)
    monkeypatch.setenv("ATLAS_ADAPTIVE_CHUNK_MIN_SIZE", "1")
    ligands = [f"lig_{idx:05d}" for idx in range(512)]

    chunks_low = planner_chunking._split_ligands_into_chunks(
        ligands, task_count=1, local_workers_per_task=8
    )
    chunks_high = planner_chunking._split_ligands_into_chunks(
        ligands, task_count=3, local_workers_per_task=112
    )

    assert len(chunks_high) > len(chunks_low)
    assert max(len(chunk) for chunk in chunks_high) <= 4
    assert sorted(item for chunk in chunks_high for item in chunk) == sorted(ligands)


def test_split_ligands_into_chunks_uses_complexity_weights() -> None:
    ligands = [f"lig_{idx:02d}" for idx in range(16)]
    weights = {lig: 1.0 for lig in ligands}
    weights["lig_00"] = 50.0
    weights["lig_01"] = 40.0
    chunks = planner_chunking._split_ligands_into_chunks(
        ligands,
        task_count=2,
        local_workers_per_task=2,
        ligand_weights=weights,
    )
    weighted_loads = [
        sum(float(weights.get(item, 1.0)) for item in chunk) for chunk in chunks
    ]
    assert sorted(item for chunk in chunks for item in chunk) == sorted(ligands)
    assert max(weighted_loads) - min(weighted_loads) < 50.0


def test_distributed_prep_ready_signal_writes_reuse_state(tmp_path: Path) -> None:
    state_path = tmp_path / "prep_state.json"
    claim_path = tmp_path / "prep_claim.lock"
    claim_path.write_text("claimed\n", encoding="utf-8")
    cfg = {
        "RUN_ID": "rid",
        "_DISTRIBUTED_PREP_READY_SIGNAL": {
            "state_path": str(state_path),
            "claim_path": str(claim_path),
            "combo_key": "scope_a",
            "task_id": 2,
            "task_count": 3,
            "owner_task_id": 2,
            "pdb_id": "BNJS",
            "variant_label": "HOLO",
            "ph_tag": "pH7_0",
        },
        "_COMBO_PREP_CENTER_BY_PH": {"pH7_0": [1.0, 2.0, 3.0]},
        "_COMBO_PREP_BOX_BY_PH": {"pH7_0": [24.0, 24.0, 24.0]},
        "_COMBO_PREP_SOURCE_BY_PH": {"pH7_0": "control_redock"},
        "_COMBO_PREP_CONTROL_STEMS": ["lig_a"],
        "_COMBO_PREP_CONTROL_LOOKUP": {"lig_a": "/tmp/lig_a.pdb"},
    }

    _signal_distributed_prep_ready(
        cfg,
        logging.getLogger("test"),
        pdb_id="BNJS",
        variant_label="HOLO",
        active_ph_label="pH7_0",
    )

    payload = json.loads(state_path.read_text(encoding="utf-8"))
    assert payload["status"] == "ready"
    assert payload["note"] == "early_prep_ready"
    assert payload["center_by_ph"] == {"pH7_0": [1.0, 2.0, 3.0]}
    assert payload["box_by_ph"] == {"pH7_0": [24.0, 24.0, 24.0]}
    assert not claim_path.exists()


def test_optimize_chunks_for_tail_splits_heavy_chunks(monkeypatch) -> None:
    monkeypatch.setenv("ATLAS_DYNAMIC_RECHUNK", "1")
    chunks = [
        {
            "chunk_id": "a",
            "pdb_id": "BNJS",
            "ph_tag": "pH7_0",
            "library_name": "bench",
            "ligand_bases": [f"lig_{i:03d}" for i in range(300)],
            "ligand_count": 300,
            "weight": 300.0,
            "assigned_task_id": 0,
        },
        {
            "chunk_id": "b",
            "pdb_id": "BOJG",
            "ph_tag": "pH7_2",
            "library_name": "bench",
            "ligand_bases": [f"lig_b_{i:03d}" for i in range(64)],
            "ligand_count": 64,
            "weight": 64.0,
            "assigned_task_id": 1,
        },
    ]
    out, stats = optimize_chunks_for_tail(
        chunks,
        run_id="rid",
        variant_label="HOLO",
        task_ids=[0, 1, 2, 3],
        task_count=4,
        local_workers_per_task=4,
        default_target_size=96,
        min_chunk_size=32,
    )
    assert int(stats.get("enabled", 0)) == 1
    assert int(stats.get("split_count", 0)) >= 1
    assert int(stats.get("chunks_after", 0)) > int(stats.get("chunks_before", 0))
    assert sum(int(x.get("ligand_count", 0) or 0) for x in out) == 364


def test_build_candidate_ids_hedging_interleaves(monkeypatch) -> None:
    monkeypatch.setenv("ATLAS_CHUNK_HEDGING", "1")
    ordered, stats = build_candidate_ids(
        preferred_ids=["p1", "p2", "p3"],
        steal_ids=["s1", "s2"],
    )
    assert int(stats.get("hedging_enabled", 0)) == 1
    assert ordered[:4] == ["p1", "s1", "p2", "s2"]


def test_build_candidate_ids_auto_hedging_interleaves(monkeypatch) -> None:
    monkeypatch.delenv("ATLAS_CHUNK_HEDGING", raising=False)
    ordered, stats = build_candidate_ids(
        preferred_ids=["p1", "p2", "p3"],
        steal_ids=["s1", "s2"],
        auto_enable_hedging=True,
    )
    assert int(stats.get("hedging_enabled", 0)) == 1
    assert int(stats.get("auto_enabled", 0)) == 1
    assert ordered[:4] == ["p1", "s1", "p2", "s2"]
