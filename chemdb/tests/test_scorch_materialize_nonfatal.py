from __future__ import annotations

import csv
from contextlib import contextmanager
from pathlib import Path

import post_docking.rescoring.rescoring_scorch as rescoring_scorch


def test_score_stage_file_not_found_in_materialize_is_nonfatal(
    monkeypatch, tmp_path: Path
) -> None:
    spec = rescoring_scorch.StageSpec(
        source="vina",
        stage_dir="vina_best",
        output_name="scorch_scores_vina_best.csv",
    )
    combo = ("BNJS", "HOLO", "pH7_0")
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("REMARK receptor\n", encoding="utf-8")

    run_root = tmp_path / "docked" / "run_x"
    post_root = tmp_path / "post_docked" / "run_x"
    ph_root = run_root / combo[0] / combo[1] / combo[2] / "stage3"
    ph_root.mkdir(parents=True, exist_ok=True)
    ligand = ph_root / "actives_final_00001_stage3.pdbqt"
    ligand.write_text("REMARK ligand\n", encoding="utf-8")

    monkeypatch.setattr(rescoring_scorch, "_emit_task_event", lambda *a, **k: None)

    @contextmanager
    def _fake_acquire(*args, **kwargs):
        del args, kwargs
        yield 1

    monkeypatch.setattr(rescoring_scorch, "acquire_global_cores", _fake_acquire)
    monkeypatch.setattr(
        rescoring_scorch,
        "_collect_best_pose_per_base",
        lambda *a, **k: ([ligand], {3: 1}, 1, {"actives_final_00001"}),
    )

    def _raise_missing(*args, **kwargs):
        del args, kwargs
        raise FileNotFoundError("actives_final_00001_stage3.pdbqt")

    monkeypatch.setattr(rescoring_scorch, "_materialize_inputs", _raise_missing)
    monkeypatch.setattr(
        rescoring_scorch.subprocess,
        "run",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("subprocess.run should not be called after nonfatal materialize miss")
        ),
    )

    ok, out_path = rescoring_scorch._score_stage(
        cfg={},
        spec=spec,
        combo=combo,
        run_root=run_root,
        post_root=post_root,
        receptor=receptor,
        threads=2,
        overwrite=True,
        logger=rescoring_scorch.logging.getLogger("test.scorch.materialize.nonfatal"),
        allowed_bases=None,
        control_bases=None,
        run_mode="fda",
    )
    assert ok is True
    assert out_path is None


def test_score_stage_passes_chunk_tag_to_materialize(monkeypatch, tmp_path: Path) -> None:
    spec = rescoring_scorch.StageSpec(
        source="vina",
        stage_dir="vina_best",
        output_name="scorch_scores_vina_best.csv",
    )
    combo = ("BNJS", "HOLO", "pH7_0")
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("REMARK receptor\n", encoding="utf-8")

    run_root = tmp_path / "docked" / "run_x"
    post_root = tmp_path / "post_docked" / "run_x"
    ph_root = run_root / combo[0] / combo[1] / combo[2] / "stage3"
    ph_root.mkdir(parents=True, exist_ok=True)
    ligand = ph_root / "actives_final_00001_stage3.pdbqt"
    ligand.write_text("REMARK ligand\n", encoding="utf-8")

    monkeypatch.setattr(rescoring_scorch, "_emit_task_event", lambda *a, **k: None)

    @contextmanager
    def _fake_acquire(*args, **kwargs):
        del args, kwargs
        yield 1

    monkeypatch.setattr(rescoring_scorch, "acquire_global_cores", _fake_acquire)

    seen: dict[str, object] = {}

    def _fake_materialize(
        ph_root_arg,
        combo_post_root_arg,
        stage_dir_arg,
        ligands_arg,
        overwrite_arg,
        logger_arg,
        *,
        run_mode,
        chunk_tag=None,
    ):
        del ph_root_arg, overwrite_arg, logger_arg
        seen["chunk_tag"] = chunk_tag
        seen["ligands"] = [p.name for p in ligands_arg]
        input_dir = combo_post_root_arg / ".scorch_inputs" / run_mode / stage_dir_arg
        if chunk_tag:
            input_dir = input_dir / str(chunk_tag)
        input_dir.mkdir(parents=True, exist_ok=True)
        staged = input_dir / ligands_arg[0].name
        staged.write_text("POSE\n", encoding="utf-8")
        return rescoring_scorch.MaterializeResult(
            input_dir=input_dir,
            materialized_count=1,
            failed_count=0,
            failed_examples=(),
        )

    monkeypatch.setattr(rescoring_scorch, "_materialize_inputs", _fake_materialize)

    def _fake_run_with_lock_fallback(**kwargs):
        cmd = list(kwargs["cmd"])
        out_path = Path(cmd[cmd.index("--out") + 1])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["Ligand_ID", "SCORCH_score", "SCORCH_certainty"],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "Ligand_ID": "actives_final_00001_stage3",
                    "SCORCH_score": "0.5",
                    "SCORCH_certainty": "0.8",
                }
            )

        class _Proc:
            returncode = 0
            stderr = ""

        return _Proc(), False

    monkeypatch.setattr(
        rescoring_scorch._scorch_stage_runner_mod,
        "run_with_lock_fallback",
        _fake_run_with_lock_fallback,
    )

    ok, out_path = rescoring_scorch._score_stage(
        cfg={},
        spec=spec,
        combo=combo,
        run_root=run_root,
        post_root=post_root,
        receptor=receptor,
        threads=2,
        overwrite=True,
        logger=rescoring_scorch.logging.getLogger(
            "test.scorch.materialize.chunk_forward"
        ),
        allowed_bases={"actives_final_00001"},
        control_bases=None,
        run_mode="fda",
        chunk_tag="part001",
    )

    assert ok is True
    assert out_path is not None
    assert seen["chunk_tag"] == "part001"
    assert seen["ligands"] == [ligand.name]
