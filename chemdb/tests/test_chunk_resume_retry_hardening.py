from __future__ import annotations

import threading
from pathlib import Path

from cli.distributed_chunk_runtime_helpers import collect_terminal_chunk_failures
from cli.planner_validation import _verify_chunk_combo_outputs
from cli.run_process_one import (
    ProcessOneContext,
    ProcessOneDeps,
    ProcessOneSharedState,
    build_process_one_runner,
)


def test_verify_chunk_outputs_missing_marker_and_summary_is_retryable(
    tmp_path: Path,
) -> None:
    ok, reason = _verify_chunk_combo_outputs(
        {"DOCKED_DIR": str(tmp_path / "docked")},
        run_id="rid_retryable_verify",
        chunk_id="chunk_a",
        pdb_id="BOJG",
        variant_label="HOLO",
        ph_tag="pH7_2",
        library_name="bench2",
        chunk_ligand_bases=["lig_a"],
    )
    assert ok is False
    assert reason.startswith("retryable_missing_stage1_chunk_marker:")
    assert ":missing_docking_summary:" in reason


def test_chunk_mode_retryable_verification_failure_stays_chunk_local(
    tmp_path: Path,
) -> None:
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


def test_chunk_mode_success_does_not_mark_whole_target_completed(
    tmp_path: Path,
) -> None:
    failed_root = tmp_path / "failed"
    failed_root.mkdir(parents=True, exist_ok=True)
    cfg: dict[str, object] = {"CPU": 1}
    success_updates: list[tuple[tuple[object, ...], dict[str, object]]] = []
    context = ProcessOneContext(
        run_id="rid_chunk_success",
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
        update_manifest_for_protein_success=lambda *args, **kwargs: success_updates.append(
            (args, kwargs)
        ),
        update_manifest_for_protein_failure=lambda *_args, **_kwargs: None,
        verify_chunk_combo_outputs=lambda *_args, **_kwargs: (True, "ok"),
        resolve_combo_output_dir=lambda *_args, **_kwargs: tmp_path,
        resolve_combo_docking_summary_csv=lambda *_args, **_kwargs: tmp_path / "summary.csv",
        load_scored_ligand_keys_from_summary=lambda *_args, **_kwargs: {"lig_a"},
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

    assert ok is True
    assert success_updates == []
    assert shared_state.failed_entries == []


def test_resume_scorch_only_does_not_rerun_target_docking(tmp_path: Path) -> None:
    failed_root = tmp_path / "failed"
    failed_root.mkdir(parents=True, exist_ok=True)
    cfg: dict[str, object] = {
        "CPU": 1,
        "USE_SCORCH": True,
        "_RESUME_MANIFEST_PROTEINS": {
            "BOJG|HOLO|pH7_2": {
                "status": "running",
                "stages": {
                    "distributed_chunk_coverage": {
                        "status": "running",
                        "details": {
                            "planned_chunks": 2,
                            "completed_chunks": 2,
                            "failed_chunks": 0,
                            "missing_result_chunks": 0,
                            "planned_ligands": 10,
                            "completed_ligands": 10,
                            "has_consensus_csv": True,
                            "require_scorch": True,
                            "has_scorch_reranked_csv": False,
                        },
                    },
                    "postprocessing": {
                        "status": "running",
                        "details": {"has_scorch_reranked_csv": False},
                    },
                },
            }
        },
    }
    context = ProcessOneContext(
        run_id="rid_resume_scorch_only",
        label="HOLO",
        mode="distributed",
        variant="holo",
        is_resume=True,
        completed_combo_lookup=set(),
        dist_combo_chunk_mode=False,
        cfg_v=cfg,
        stages=None,
        params=None,
        tokens=[],
        global_start=0.0,
    )
    submitted: list[dict[str, object]] = []

    class FakeScorchQueue:
        def submit(self, **kwargs):
            submitted.append(dict(kwargs))
            return True

    shared_state = ProcessOneSharedState(
        failed_root=str(failed_root),
        failed_entries=[],
        run_scope_completion_times={},
        retention_lock=threading.Lock(),
        pending_retention_combos=set(),
        pending_coverage_refresh={},
        pending_retention_lock=threading.Lock(),
        scorch_queue_service=FakeScorchQueue(),
    )
    start_updates: list[object] = []
    deps = ProcessOneDeps(
        normalize_ph_tag_token=lambda ph: str(ph or "").strip() or "base",
        chunk_ligand_key=lambda name: str(name),
        process_one_protein=lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("target docking should not run")
        ),
        update_manifest_for_protein_start=lambda *args, **_kwargs: start_updates.append(
            args
        ),
        update_manifest_for_protein_success=lambda *_args, **_kwargs: None,
        update_manifest_for_protein_failure=lambda *_args, **_kwargs: None,
        verify_chunk_combo_outputs=lambda *_args, **_kwargs: (True, "ok"),
        resolve_combo_output_dir=lambda *_args, **_kwargs: tmp_path,
        resolve_combo_docking_summary_csv=lambda *_args, **_kwargs: tmp_path / "summary.csv",
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

    ok = process_one("BOJG.pdb", cfg, "pH7_2")

    assert ok is True
    assert start_updates == []
    assert len(submitted) == 1
    assert submitted[0]["pdb_id"] == "BOJG"
    assert submitted[0]["allowed_count_hint"] == 10
    assert ("BOJG", "HOLO", "pH7_2") in shared_state.pending_retention_combos


def test_collect_terminal_chunk_failures_only_records_terminal_attempts() -> None:
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
