from __future__ import annotations

import csv
import logging
from pathlib import Path

from analysis.reporting import master_schema_export
from cli.distributed_runtime_completion import DeferredScorchCompletionBookkeeper
from post_docking.rescoring import scorch_provisional_cache as cache
from post_docking.rescoring.scorch_stage_runner import score_stage
from post_docking.rescoring.scorch_types import AnnotateResult, StageSpec


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def test_provisional_cache_reuses_only_exact_row_keys(tmp_path: Path) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path), "DATA_DIR": str(tmp_path / "data")}
    run_id = "run_exact"
    combo = ("1ABC", "HOLO", "pH7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    ligand = tmp_path / "lig1_stage1.pdbqt"
    ligand.write_text("LIGAND A\n", encoding="utf-8")
    output_csv = tmp_path / "scores.csv"
    _write_csv(
        output_csv,
        [
            {
                "ligand_file": ligand.name,
                "Ligand_ID": ligand.name,
                "SCORCH_score": "0.75",
                "SCORCH_certainty": "0.9",
            }
        ],
    )

    cache.record_task_output(
        cfg=cfg,
        run_id=run_id,
        combo=combo,
        spec=spec,
        run_mode="fda",
        decoy_prefix="dud",
        receptor=receptor,
        score_csv=None,
        allowed_bases={"lig1"},
        selected_stage_by_base={"lig1": "stage1"},
        selected_score_by_base={"lig1": -9.1},
        ligands=[ligand],
        output_csv=output_csv,
        phase="provisional",
        logger=logging.getLogger("test"),
    )

    fields, rows, covered = cache.reusable_rows_for_task(
        cfg=cfg,
        run_id=run_id,
        combo=combo,
        spec=spec,
        run_mode="fda",
        decoy_prefix="dud",
        allowed_bases={"lig1"},
        selected_stage_by_base={"lig1": "stage1"},
        selected_score_by_base={"lig1": -9.1},
        ligands=[ligand],
        receptor=receptor,
        logger=logging.getLogger("test"),
    )
    assert covered == {"lig1"}
    assert rows[0]["SCORCH_score"] == "0.75"
    assert "SCORCH_score" in fields

    _fields, rows, covered = cache.reusable_rows_for_task(
        cfg=cfg,
        run_id=run_id,
        combo=combo,
        spec=spec,
        run_mode="fda",
        decoy_prefix="dud",
        allowed_bases={"lig1"},
        selected_stage_by_base={"lig1": "stage1"},
        selected_score_by_base={"lig1": -8.0},
        ligands=[ligand],
        receptor=receptor,
        logger=logging.getLogger("test"),
    )
    assert covered == set()
    assert rows == []

    receptor.write_text("RECEPTOR CHANGED\n", encoding="utf-8")
    _fields, rows, covered = cache.reusable_rows_for_task(
        cfg=cfg,
        run_id=run_id,
        combo=combo,
        spec=spec,
        run_mode="fda",
        decoy_prefix="dud",
        allowed_bases={"lig1"},
        selected_stage_by_base={"lig1": "stage1"},
        selected_score_by_base={"lig1": -9.1},
        ligands=[ligand],
        receptor=receptor,
        logger=logging.getLogger("test"),
    )
    assert covered == set()
    assert rows == []


def test_score_stage_full_cache_reuse_skips_scorch_execution(tmp_path: Path) -> None:
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "DATA_DIR": str(tmp_path / "data"),
    }
    run_id = "run_reuse"
    combo = ("1ABC", "HOLO", "pH7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    run_root = tmp_path / "docked" / run_id
    stage1 = run_root / "1ABC" / "HOLO" / "pH7_0" / "stage1"
    stage1.mkdir(parents=True)
    ligand = stage1 / "lig1_stage1.pdbqt"
    ligand.write_text("LIGAND A\n", encoding="utf-8")
    receptor = tmp_path / "processed" / "1ABC" / "HOLO" / "receptor" / "ph_ensemble" / "1ABC_pH7_0.pdbqt"
    receptor.parent.mkdir(parents=True)
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    seed_csv = tmp_path / "seed_scores.csv"
    _write_csv(
        seed_csv,
        [
            {
                "pdb_id": "1ABC",
                "variant": "HOLO",
                "ph": "pH7_0",
                "source": "vina",
                "stage_dir": "vina_best",
                "run_mode": "fda",
                "ligand_file": ligand.name,
                "Ligand_ID": ligand.name,
                "SCORCH_score": "0.88",
            }
        ],
    )
    cache.record_task_output(
        cfg=cfg,
        run_id=run_id,
        combo=combo,
        spec=spec,
        run_mode="fda",
        decoy_prefix="dud",
        receptor=receptor,
        score_csv=None,
        allowed_bases={"lig1"},
        selected_stage_by_base={"lig1": "stage1"},
        selected_score_by_base={"lig1": -9.1},
        ligands=[ligand],
        output_csv=seed_csv,
        phase="provisional",
        logger=logging.getLogger("test"),
    )
    stale_part = (
        tmp_path
        / "post_docked"
        / run_id
        / "1ABC"
        / "HOLO"
        / "pH7_0"
        / "scorch_scores_vina_best.part000.csv"
    )
    _write_csv(stale_part, [{"ligand_file": "stale.pdbqt", "SCORCH_score": "99"}])

    def _unexpected(*_args, **_kwargs):
        raise AssertionError("SCORCH execution should not run on a full cache hit")

    ok, out_path = score_stage(
        cfg,
        spec,
        combo,
        run_root,
        tmp_path / "post_docked" / run_id,
        receptor,
        1,
        False,
        logging.getLogger("test"),
        decoy_prefix="dud",
        component="[test]",
        emit_task_event=lambda *_args, **_kwargs: None,
        materialize_inputs=_unexpected,
        scorch_command=_unexpected,
        annotate_csv=lambda *_args, **_kwargs: AnnotateResult(True, False),
        scorch_root=None,
        scorch_python=None,
        scorch_script=None,
        allowed_bases={"lig1"},
        control_bases=set(),
        selected_stage_by_base={"lig1": "stage1"},
        selected_score_by_base={"lig1": -9.1},
        cache_read=True,
        run_id=run_id,
    )

    assert ok is True
    assert out_path is not None
    with out_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["SCORCH_score"] == "0.88"
    assert not stale_part.exists()


def test_score_stage_existing_output_cleans_stale_part_files(tmp_path: Path) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path)}
    run_id = "run_existing"
    combo = ("1ABC", "HOLO", "pH7_0")
    spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")
    post_combo = tmp_path / "post_docked" / run_id / "1ABC" / "HOLO" / "pH7_0"
    final_csv = post_combo / "scorch_scores_vina_best.csv"
    stale_part = post_combo / "scorch_scores_vina_best.part000.csv"
    _write_csv(final_csv, [{"ligand_file": "fresh.pdbqt", "SCORCH_score": "1"}])
    _write_csv(stale_part, [{"ligand_file": "stale.pdbqt", "SCORCH_score": "99"}])
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")

    def _unexpected(*_args, **_kwargs):
        raise AssertionError("existing output should skip SCORCH execution")

    ok, out_path = score_stage(
        cfg,
        spec,
        combo,
        tmp_path / "docked" / run_id,
        tmp_path / "post_docked" / run_id,
        receptor,
        1,
        False,
        logging.getLogger("test"),
        decoy_prefix="dud",
        component="[test]",
        emit_task_event=lambda *_args, **_kwargs: None,
        materialize_inputs=_unexpected,
        scorch_command=_unexpected,
        annotate_csv=lambda *_args, **_kwargs: AnnotateResult(True, False),
        scorch_root=None,
        scorch_python=None,
        scorch_script=None,
        allowed_bases={"fresh"},
        control_bases=set(),
        run_id=run_id,
    )

    assert ok is True
    assert out_path == final_csv
    assert not stale_part.exists()


def test_master_schema_discovery_ignores_hidden_provisional_cache(tmp_path: Path) -> None:
    run_root = tmp_path / "repo"
    visible = run_root / "post_docked" / "run1" / "1ABC" / "HOLO" / "pH7_0"
    hidden = (
        run_root
        / "outputs"
        / "data"
        / "run1"
        / ".scorch_provisional_cache"
        / "artifacts"
        / "post_docked"
        / "run1"
        / "1ABC"
        / "HOLO"
        / "pH7_0"
    )
    _write_csv(visible / "consensus_reranked_scorch.csv", [{"ligand": "visible"}])
    _write_csv(hidden / "consensus_reranked_scorch.csv", [{"ligand": "hidden"}])

    found = master_schema_export._discover_consensus_files(
        run_root,
        ["fda"],
        "dud",
        logging.getLogger("test"),
    )

    assert found == [visible / "consensus_reranked_scorch.csv"]


def test_provisional_submission_does_not_mark_final_submitted() -> None:
    class FakeQueue:
        def __init__(self) -> None:
            self.calls: list[dict[str, object]] = []

        def submit_provisional(self, **kwargs):
            self.calls.append(kwargs)
            return True

    scope = ("1ABC", "HOLO", "pH7_0")
    queue = FakeQueue()
    bookkeeper = DeferredScorchCompletionBookkeeper(
        cfg_v={
            "SCORCH_PROVISIONAL_ENABLE": True,
            "SCORCH_PROVISIONAL_MIN_PROGRESS": 0.10,
        },
        run_id="run1",
        label="HOLO",
        dist_task_id=0,
        combo_total_chunks={scope: 10},
        combo_total_ligands={scope: 100},
        combo_owner_task_id={scope: 0},
        scorch_queue_service=queue,
        variant_param_for_scorch="HOLO",
        hybrid_start_ts=0.0,
        acquire_global_slot=lambda *_args, **_kwargs: None,
        maybe_run_scorch_rescore_for_pdb=lambda *_args, **_kwargs: None,
        normalize_ph_tag_token=lambda value: str(value or "base"),
    )

    bookkeeper.mark_chunk_completed(scope, reason="test")
    queued = bookkeeper.reconcile_provisional_owned_submissions(
        reason="idle",
        active_workers=3,
        idle_workers=1,
        unresolved_chunks=9,
    )

    assert queued == 1
    assert len(queue.calls) == 1
    assert bookkeeper.combo_scorch_submitted == set()
    assert bookkeeper.combo_provisional_submitted
