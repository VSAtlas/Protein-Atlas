from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

from cli import distributed_chunk_planner as chunk_planner
from cli.distributed_context import DistributedRunContext


def _combo_dock_dir(cfg: dict[str, str], run_id: str) -> Path:
    return Path(str(cfg["DOCKED_DIR"])) / run_id / "BNJS" / "HOLO" / "pH7_0"


def _combo_post_dir(cfg: dict[str, str], run_id: str) -> Path:
    return Path(str(cfg["POST_DOCKED_DIR"])) / run_id / "BNJS" / "HOLO" / "pH7_0"


def _write_docking_consensus(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ligand", "run_mode", "is_decoy", "consensus_score"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "ligand": "actives_final_00001.pdbqt",
                "run_mode": "fda",
                "is_decoy": "0",
                "consensus_score": "-9.0",
            }
        )


def _write_raw_scorch_csv(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Ligand_ID", "run_mode", "score"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "Ligand_ID": "actives_final_00001_stage1",
                "run_mode": "fda",
                "score": "0.9",
            }
        )


def _write_post_csv(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ligand_base",
                "ligand",
                "run_mode",
                "rescored_flag",
                "final_score",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "ligand_base": "actives_final_00001",
                "ligand": "actives_final_00001.pdbqt",
                "run_mode": "fda",
                "rescored_flag": "1",
                "final_score": "1.0",
            }
        )


def _write_complete_scorch_evidence(cfg: dict[str, str], run_id: str) -> None:
    _write_docking_consensus(
        _combo_dock_dir(cfg, run_id) / "consensus_docking_scores.csv"
    )
    combo_post = _combo_post_dir(cfg, run_id)
    _write_raw_scorch_csv(combo_post / "scorch_scores_all.csv")
    _write_post_csv(combo_post / "consensus_reranked_scorch.csv")


def test_reconcile_missing_post_outputs_reruns_scorch(
    tmp_path: Path, monkeypatch
) -> None:
    run_id = "reconcile_smoke"
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "POST_DOCKED_DIR": str(tmp_path / "post_docked"),
        "DISTRIBUTED_MODE": "off",
    }
    _write_docking_consensus(_combo_dock_dir(cfg, run_id) / "consensus_docking_scores.csv")
    calls: list[dict[str, object]] = []

    def _fake_rerun(
        cfg_in,
        run_id_in,
        pdb_id,
        *,
        variant=None,
        ph=None,
        allowed_count_hint=None,
        execution_mode="subprocess",
        overwrite=False,
        verbose=False,
        **_kwargs,
    ):
        calls.append(
            {
                "run_id": run_id_in,
                "pdb_id": pdb_id,
                "variant": variant,
                "ph": ph,
                "allowed_count_hint": allowed_count_hint,
                "execution_mode": execution_mode,
                "overwrite": overwrite,
                "verbose": verbose,
            }
        )
        combo_post = (
            Path(str(cfg_in["POST_DOCKED_DIR"]))
            / str(run_id_in)
            / str(pdb_id)
            / str(variant)
            / str(ph)
        )
        _write_raw_scorch_csv(combo_post / "scorch_scores_all.csv")
        _write_post_csv(combo_post / "consensus_reranked_scorch.csv")
        (combo_post / "scorch").mkdir(parents=True, exist_ok=True)
        (combo_post / "scorch" / "_DONE").write_text("ok\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(chunk_planner, "_maybe_run_scorch_rescore_for_pdb", _fake_rerun)
    monkeypatch.setattr(
        chunk_planner,
        "_resolve_scorch_execution_mode",
        lambda _cfg: "inprocess",
    )

    summary = chunk_planner._reconcile_missing_post_scorch_outputs(
        cfg,
        run_id=run_id,
        dist_ctx=SimpleNamespace(enabled=False, task_id=0),
        coverage_refresh={
            ("BNJS", "HOLO", "pH7_0"): {
                "library": "bench_pur2",
                "expected": {"actives_final_00001"},
            }
        },
    )

    assert summary == {"attempted": 1, "succeeded": 1, "failed": 0, "skipped_owner": 0}
    assert len(calls) == 1
    assert calls[0]["overwrite"] is False
    assert calls[0]["execution_mode"] == "inprocess"
    assert calls[0]["allowed_count_hint"] == 1


def test_reconcile_marks_done_for_complete_post_csv_without_rerun(
    tmp_path: Path, monkeypatch
) -> None:
    run_id = "reconcile_mark_done"
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "POST_DOCKED_DIR": str(tmp_path / "post_docked"),
        "DISTRIBUTED_MODE": "off",
    }
    _write_complete_scorch_evidence(cfg, run_id)
    combo_post = _combo_post_dir(cfg, run_id)

    def _fail_if_called(*_args, **_kwargs):
        raise AssertionError("complete post-SCORCH output should only be marked done")

    monkeypatch.setattr(
        chunk_planner, "_maybe_run_scorch_rescore_for_pdb", _fail_if_called
    )

    summary = chunk_planner._reconcile_missing_post_scorch_outputs(
        cfg,
        run_id=run_id,
        dist_ctx=SimpleNamespace(enabled=False, task_id=0),
        coverage_refresh={
            ("BNJS", "HOLO", "pH7_0"): {
                "library": "bench_pur2",
                "expected": {"actives_final_00001"},
            }
        },
    )

    assert summary == {"attempted": 0, "succeeded": 1, "failed": 0, "skipped_owner": 0}
    assert (combo_post / "scorch" / "_DONE").read_text(encoding="utf-8") == "ok\n"


def test_reconcile_missing_post_outputs_skips_non_owner(
    tmp_path: Path, monkeypatch
) -> None:
    run_id = "reconcile_owner"
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "POST_DOCKED_DIR": str(tmp_path / "post_docked"),
        "DISTRIBUTED_MODE": "slurm_array",
    }
    ctx = DistributedRunContext(
        mode="slurm_array",
        enabled=True,
        run_id=run_id,
        task_id=0,
        task_count=2,
        task_min_id=0,
        leader_task_id=0,
        barrier_timeout_sec=1.0,
        barrier_poll_sec=0.01,
    )

    monkeypatch.setattr(
        chunk_planner,
        "_distributed_combo_owner_task_id",
        lambda *a, **k: 1,
    )

    def _fail_if_called(*_args, **_kwargs):
        raise AssertionError("scorch rerun should not be called by non-owner task")

    monkeypatch.setattr(
        chunk_planner,
        "_maybe_run_scorch_rescore_for_pdb",
        _fail_if_called,
    )

    summary = chunk_planner._reconcile_missing_post_scorch_outputs(
        cfg,
        run_id=run_id,
        dist_ctx=ctx,
        coverage_refresh={
            ("BNJS", "HOLO", "pH7_0"): {
                "library": "bench_pur2",
                "expected": {"actives_final_00001"},
            }
        },
    )
    assert summary == {"attempted": 0, "succeeded": 0, "failed": 0, "skipped_owner": 1}
