from __future__ import annotations

import json
from types import SimpleNamespace

from tools import finalize_distributed_run


def test_finalize_collects_nonterminal_manifest_entries_without_consensus(
    monkeypatch, tmp_path
) -> None:
    docked_root = tmp_path / "docked"
    done_root = docked_root / "rid" / "RUNN"
    done_root.mkdir(parents=True)
    (done_root / "consensus_docking_scores.csv").write_text(
        "ligand,best_score\nlig,1\n",
        encoding="utf-8",
    )
    manifest = {
        "proteins": {
            "DONE|LEGACY|base": {
                "pdb_id": "DONE",
                "variant": "LEGACY",
                "ph": "base",
                "status": "completed",
            },
            "RUNN|LEGACY|base": {
                "pdb_id": "RUNN",
                "variant": "LEGACY",
                "ph": "base",
                "status": "running",
            },
            "PEND|LEGACY|base": {
                "pdb_id": "PEND",
                "variant": "LEGACY",
                "ph": "base",
                "status": "pending",
            },
        }
    }
    monkeypatch.setattr(
        finalize_distributed_run,
        "load_run_manifest",
        lambda _cfg, _run_id: manifest,
    )

    entries = finalize_distributed_run._collect_incomplete_entries_from_manifest(
        {"DOCKED_DIR": str(docked_root)}, "rid"
    )

    assert (
        "RUNN",
        "LEGACY",
        "base",
        "DistributedIncomplete",
        "nonterminal_status=running;coverage_unresolved_nonterminal",
    ) in entries
    assert (
        "PEND",
        "LEGACY",
        "base",
        "DistributedIncomplete",
        "nonterminal_status=pending;missing_consensus_csv",
    ) in entries
    assert (
        "DONE",
        "LEGACY",
        "base",
        "DistributedIncomplete",
        "missing_consensus_csv",
    ) in entries


def test_finalize_requires_chunk_plan_coverage_for_consensus(monkeypatch, tmp_path) -> None:
    docked_root = tmp_path / "docked"
    manifests_root = tmp_path / "manifests"
    for pdb_id, rows in {"FULL": ["a", "b"], "PART": ["a"]}.items():
        out_dir = docked_root / "rid" / pdb_id
        out_dir.mkdir(parents=True)
        (out_dir / "consensus_docking_scores.csv").write_text(
            "ligand,best_score\n" + "\n".join(f"{lig},1" for lig in rows) + "\n",
            encoding="utf-8",
        )

    dist_dir = manifests_root / "rid" / "distributed"
    (dist_dir / "chunk_results").mkdir(parents=True)
    plan = {
        "chunks": [
            {
                "chunk_id": "full_a",
                "pdb_id": "FULL",
                "variant_label": "LEGACY",
                "ph_tag": "base",
                "ligand_bases": ["a"],
            },
            {
                "chunk_id": "full_b",
                "pdb_id": "FULL",
                "variant_label": "LEGACY",
                "ph_tag": "base",
                "ligand_bases": ["b"],
            },
            {
                "chunk_id": "part_a",
                "pdb_id": "PART",
                "variant_label": "LEGACY",
                "ph_tag": "base",
                "ligand_bases": ["a"],
            },
            {
                "chunk_id": "part_b",
                "pdb_id": "PART",
                "variant_label": "LEGACY",
                "ph_tag": "base",
                "ligand_bases": ["b"],
            },
        ]
    }
    (dist_dir / "combo_chunks_0.json").write_text(
        finalize_distributed_run.json.dumps(plan),
        encoding="utf-8",
    )
    for chunk_id in ("full_a", "full_b", "part_a"):
        (dist_dir / "chunk_results" / f"{chunk_id}.json").write_text(
            finalize_distributed_run.json.dumps(
                {"chunk_id": chunk_id, "status": "completed"}
            ),
            encoding="utf-8",
        )
    manifest = {
        "proteins": {
            "FULL|LEGACY|base": {
                "pdb_id": "FULL",
                "variant": "LEGACY",
                "ph": "base",
                "status": "completed",
            },
            "PART|LEGACY|base": {
                "pdb_id": "PART",
                "variant": "LEGACY",
                "ph": "base",
                "status": "completed",
            },
        }
    }
    monkeypatch.setattr(
        finalize_distributed_run,
        "load_run_manifest",
        lambda _cfg, _run_id: manifest,
    )

    entries = finalize_distributed_run._collect_incomplete_entries_from_manifest(
        {
            "DOCKED_DIR": str(docked_root),
            "MANIFESTS_DIR": str(manifests_root),
        },
        "rid",
    )

    assert all(entry[0] != "FULL" for entry in entries)
    assert (
        "PART",
        "LEGACY",
        "base",
        "DistributedIncomplete",
        "chunk_plan_incomplete=1/2",
    ) in entries


def test_finalize_treats_required_scorch_repair_incomplete_as_failed() -> None:
    complete = SimpleNamespace(
        targets=2,
        targets_completed=2,
        targets_running=0,
        targets_failed=0,
    )
    incomplete = SimpleNamespace(
        targets=2,
        targets_completed=1,
        targets_running=1,
        targets_failed=0,
    )

    assert finalize_distributed_run._scorch_repair_summary_incomplete(True, complete) is False
    assert finalize_distributed_run._scorch_repair_summary_incomplete(True, incomplete) is True
    assert finalize_distributed_run._scorch_repair_summary_incomplete(False, incomplete) is False


def test_finalize_collects_scorch_queue_marker_failures(tmp_path) -> None:
    marker = (
        tmp_path
        / "manifests"
        / "rid"
        / "distributed"
        / "phase_markers"
        / "scorch_drain"
        / "task_0.json"
    )
    marker.parent.mkdir(parents=True)
    marker.write_text(
        json.dumps({"scorch_queue_failed": 1, "scorch_reconcile_failed": 2}),
        encoding="utf-8",
    )

    failures = finalize_distributed_run._collect_scorch_queue_marker_failures(
        {"MANIFESTS_DIR": str(tmp_path / "manifests")}, "rid"
    )

    assert ("scorch_drain", "task_0.json", 1) in failures
    assert ("scorch_drain", "task_0.json", 2) in failures


def test_finalize_fails_when_scorch_queue_markers_report_failures(
    monkeypatch,
    tmp_path,
) -> None:
    marker = (
        tmp_path
        / "manifests"
        / "rid"
        / "distributed"
        / "phase_markers"
        / "scorch_queue_closed"
        / "task_0.json"
    )
    marker.parent.mkdir(parents=True)
    marker.write_text(json.dumps({"scorch_queue_failed": 1}), encoding="utf-8")

    incomplete = SimpleNamespace(
        targets=1,
        targets_completed=0,
        targets_running=1,
        targets_failed=0,
        completed_chunks=0,
        planned_chunks=1,
        scorch_reranked_files=0,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_parse_args",
        lambda: SimpleNamespace(run_id="rid", no_hooks=True, reconcile_only=False),
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_load_effective_cfg",
        lambda _run_id: {"USE_SCORCH": False, "MANIFESTS_DIR": str(tmp_path / "manifests")},
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "reconcile_distributed_manifest_state",
        lambda *_args, **_kwargs: {
            "action": "ok",
            "state_files": 0,
            "merged_entries": 0,
            "manifest_stale": False,
        },
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_backfill_manifest_libraries_from_chunks",
        lambda *_args, **_kwargs: 0,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_collect_failed_entries_from_markers",
        lambda *_args, **_kwargs: [],
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_collect_incomplete_entries_from_manifest",
        lambda *_args, **_kwargs: [],
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_resolve_manifest_started_epoch",
        lambda *_args, **_kwargs: 0.0,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "finalize_run_manifest",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "load_run_manifest",
        lambda *_args, **_kwargs: {"status": "completed", "summary": {"total_protein_list": []}},
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "repair_run_manifest",
        lambda *_args, **_kwargs: incomplete,
    )

    assert finalize_distributed_run.main() == 4


def test_finalize_drops_stale_marker_failures_after_manifest_repair(
    monkeypatch,
    tmp_path,
) -> None:
    complete = SimpleNamespace(
        targets=1,
        targets_completed=1,
        targets_running=0,
        targets_failed=0,
        completed_chunks=1,
        planned_chunks=1,
        scorch_reranked_files=1,
    )
    manifest = {
        "status": "completed",
        "proteins": {
            "PDB1|LEGACY|base": {
                "pdb_id": "PDB1",
                "variant": "LEGACY",
                "ph": "base",
                "status": "completed",
            }
        },
        "summary": {"total_protein_list": []},
    }
    finalized: dict[str, object] = {}

    monkeypatch.setattr(
        finalize_distributed_run,
        "_parse_args",
        lambda: SimpleNamespace(run_id="rid", no_hooks=True, reconcile_only=False),
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_load_effective_cfg",
        lambda _run_id: {"USE_SCORCH": True, "MANIFESTS_DIR": str(tmp_path / "manifests")},
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "reconcile_distributed_manifest_state",
        lambda *_args, **_kwargs: {
            "action": "ok",
            "state_files": 0,
            "merged_entries": 0,
            "manifest_stale": False,
        },
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_backfill_manifest_libraries_from_chunks",
        lambda *_args, **_kwargs: 0,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_collect_failed_entries_from_markers",
        lambda *_args, **_kwargs: [
            ("PDB1", "LEGACY", "base", "StaleWorkerFailure", "old_timeout")
        ],
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_collect_incomplete_entries_from_manifest",
        lambda *_args, **_kwargs: [],
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_resolve_manifest_started_epoch",
        lambda *_args, **_kwargs: 0.0,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "finalize_run_manifest",
        lambda *_args, **kwargs: finalized.update(
            {"failed_entries": kwargs.get("failed_entries")}
        ),
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "load_run_manifest",
        lambda *_args, **_kwargs: manifest,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "repair_run_manifest",
        lambda *_args, **_kwargs: complete,
    )

    assert finalize_distributed_run.main() == 0
    assert finalized["failed_entries"] == []


def test_finalize_ignores_stale_state_when_manifest_repair_is_complete(
    monkeypatch,
    tmp_path,
) -> None:
    marker = (
        tmp_path
        / "manifests"
        / "rid"
        / "distributed"
        / "phase_markers"
        / "scorch_drain"
        / "task_0.json"
    )
    marker.parent.mkdir(parents=True)
    marker.write_text(json.dumps({"scorch_queue_failed": 1}), encoding="utf-8")
    complete = SimpleNamespace(
        targets=1,
        targets_completed=1,
        targets_running=0,
        targets_failed=0,
        completed_chunks=1,
        planned_chunks=1,
        missing_chunks=0,
        failed_chunks=0,
        scorch_reranked_files=1,
    )
    manifest = {
        "status": "completed",
        "proteins": {
            "PDB1|LEGACY|base": {
                "pdb_id": "PDB1",
                "variant": "LEGACY",
                "ph": "base",
                "status": "completed",
            }
        },
        "summary": {"total_protein_list": []},
    }
    finalized: dict[str, object] = {}

    monkeypatch.setattr(
        finalize_distributed_run,
        "_parse_args",
        lambda: SimpleNamespace(run_id="rid", no_hooks=True, reconcile_only=False),
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_load_effective_cfg",
        lambda _run_id: {"USE_SCORCH": True, "MANIFESTS_DIR": str(tmp_path / "manifests")},
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "reconcile_distributed_manifest_state",
        lambda *_args, **_kwargs: {
            "action": "ok",
            "state_files": 0,
            "merged_entries": 0,
            "manifest_stale": False,
        },
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_backfill_manifest_libraries_from_chunks",
        lambda *_args, **_kwargs: 0,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_collect_failed_entries_from_markers",
        lambda *_args, **_kwargs: [],
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_collect_incomplete_entries_from_manifest",
        lambda *_args, **_kwargs: [
            ("PDB1", "LEGACY", "base", "DistributedIncomplete", "stale")
        ],
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "_resolve_manifest_started_epoch",
        lambda *_args, **_kwargs: 0.0,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "finalize_run_manifest",
        lambda *_args, **kwargs: finalized.update(
            {"failed_entries": kwargs.get("failed_entries")}
        ),
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "load_run_manifest",
        lambda *_args, **_kwargs: manifest,
    )
    monkeypatch.setattr(
        finalize_distributed_run,
        "repair_run_manifest",
        lambda *_args, **_kwargs: complete,
    )

    assert finalize_distributed_run.main() == 0
    assert finalized["failed_entries"] == []
