from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemdb.bench.util_bench.run_bench import main as run_bench_main
from chemdb.bench.util_bench.workloads import (
    generate_stage_collapse_workload,
    generate_synthetic_workload,
    prepare_synthetic_ligand_fixture,
    write_workload_json,
)



def test_util_bench_smoke(tmp_path: Path) -> None:
    run_root = tmp_path / "bench_runs"
    run_root.mkdir(parents=True, exist_ok=True)

    workload = generate_synthetic_workload(
        n_tasks=8,
        seed=7,
        mean_sec=0.25,
        sigma=0.25,
        outlier_fraction=0.0,
        core_weights={1: 1.0},
        run_name="smoke_util",
    )
    workload_path = write_workload_json(tmp_path / "workload.json", workload)

    rc = run_bench_main(
        [
            "--workload",
            str(workload_path),
            "--outdir",
            str(run_root),
            "--alloc-cpus",
            "2",
            "--interval-sec",
            "0.2",
            "--force-procstat",
            "--run-id",
            "smoke_run",
        ]
    )
    assert rc == 0

    out_dir = run_root / "smoke_run"
    summary_path = out_dir / "summary.json"
    assert summary_path.exists()

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    required = {
        "IdleFrac",
        "IdleCoreSeconds",
        "BusyCoreSeconds",
        "WallSeconds",
        "AllocCPUs",
        "iowait_share",
    }
    assert required.issubset(summary.keys())
    assert 0.0 <= float(summary["IdleFrac"]) <= 1.0
    assert int(summary["AllocCPUs"]) == 2
    event_fields = summary.get("event_fields", {})
    assert int(event_fields.get("task_events_with_min_cores", 0)) > 0
    assert int(event_fields.get("task_starts_with_granted_cores", 0)) > 0
    hedge = summary.get("hedge", {})
    assert "hedge_started" in hedge
    assert "hedge_won" in hedge



def test_prepare_synthetic_ligand_fixture_from_test_library10(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[2]
    fixture_subdir = "fda_test_library_10"
    fixture_root = repo_root / "chemdb" / "tests" / "fixtures" / "prepped_ligands"
    input_fixture_root = repo_root / "chemdb" / "tests" / "fixtures" / "input_pdbs"
    if not (input_fixture_root / "TEST.pdb").exists():
        pytest.skip("TEST.pdb missing")
    if not (fixture_root / fixture_subdir).is_dir():
        pytest.skip("bundled prepped ligand fixture missing")

    payload = prepare_synthetic_ligand_fixture(
        repo_root=repo_root,
        out_root=tmp_path,
        new_pdb_id="TB01",
        source_pdb_id="TEST",
        source_library_subdir=fixture_subdir,
        n_ligands=4,
    )

    mapped = payload["test_library_map"]
    assert mapped["TEST"] == fixture_subdir
    assert mapped["TB01"] == payload["clone_library_subdir"]
    assert int(payload["copied_ligands"]) >= 1

    cloned_pdb = Path(payload["input_dir"]) / "TB01.pdb"
    assert cloned_pdb.exists()


def test_stage_collapse_workload_shape() -> None:
    workload = generate_stage_collapse_workload(
        proteins=2,
        stage1_tasks=10,
        stage2_keep_frac=0.3,
        stage3_keep_frac=0.2,
        seed=11,
    )
    tasks = workload["tasks"]
    assert len(tasks) == 30
    stage_counts: dict[str, int] = {"stage1": 0, "stage2": 0, "stage3": 0}
    for row in tasks:
        tags = row.get("tags") or []
        for stage in stage_counts:
            if stage in tags:
                stage_counts[stage] += 1
                break
    assert stage_counts["stage1"] == 20
    assert stage_counts["stage2"] == 6
    assert stage_counts["stage3"] == 4
