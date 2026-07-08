from __future__ import annotations

from pathlib import Path

from cli.distributed_context import (
    DistributedRunContext,
    shard_pdb_files_for_context,
    wait_for_all_phase_markers,
    write_phase_marker,
)


def _ctx(tmp_path: Path, *, task_id: int, task_count: int = 2) -> DistributedRunContext:
    return DistributedRunContext(
        mode="slurm_array",
        enabled=True,
        run_id="run_phase_barrier",
        task_id=int(task_id),
        task_count=int(task_count),
        task_min_id=0,
        leader_task_id=0,
        barrier_timeout_sec=1.0,
        barrier_poll_sec=0.01,
    )


def test_phase_marker_wait_requires_all_tasks(tmp_path: Path) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path)}
    ctx0 = _ctx(tmp_path, task_id=0)
    ctx1 = _ctx(tmp_path, task_id=1)

    write_phase_marker(ctx0, cfg, phase="scorch_drain")

    ok_before, missing_before = wait_for_all_phase_markers(
        ctx0,
        cfg,
        phase="scorch_drain",
        timeout_sec=0.05,
        poll_sec=0.01,
    )
    assert ok_before is False
    assert missing_before == [1]

    write_phase_marker(ctx1, cfg, phase="scorch_drain")

    ok_after, missing_after = wait_for_all_phase_markers(
        ctx0,
        cfg,
        phase="scorch_drain",
        timeout_sec=0.2,
        poll_sec=0.01,
    )
    assert ok_after is True
    assert missing_after == []


def test_small_per_protein_shards_round_robin_without_empty_tasks(tmp_path: Path) -> None:
    pdbs = ["a.pdb", "b.pdb", "c.pdb"]
    assigned = [
        shard_pdb_files_for_context(_ctx(tmp_path, task_id=task_id, task_count=3), pdbs)
        for task_id in range(3)
    ]
    assert assigned == [["a.pdb"], ["b.pdb"], ["c.pdb"]]
