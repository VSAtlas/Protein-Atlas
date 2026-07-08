from __future__ import annotations

import csv
import logging
import subprocess
from contextlib import contextmanager
from pathlib import Path

from post_docking.rescoring import scorch_stage_runner
from post_docking.rescoring.scorch_types import AnnotateResult, MaterializeResult, StageSpec


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else ["Ligand_ID", "SCORCH_score"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_score_stage_uses_adaptive_timeout_for_small_shards(
    tmp_path: Path,
    monkeypatch,
) -> None:
    cfg: dict[str, object] = {"OVERALL_DIR": str(tmp_path)}
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    stage_root = run_root / combo[0] / combo[1] / combo[2] / "stage1"
    stage_root.mkdir(parents=True)
    for idx in range(7):
        (stage_root / f"lig_{idx:03d}_stage1.pdbqt").write_text("MODEL\nENDMDL\n")
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("RECEPTOR\n")

    captured: dict[str, float] = {}

    @contextmanager
    def fake_acquire_global_cores(*_args, cores=1, **_kwargs):
        yield int(cores)

    def fake_materialize_inputs(
        _ph_root,
        combo_post_root,
        _stage_dir,
        ligands_available,
        _overwrite,
        _logger,
        **_kwargs,
    ):
        input_dir = combo_post_root / ".inputs"
        input_dir.mkdir(parents=True, exist_ok=True)
        for ligand in ligands_available:
            (input_dir / Path(ligand).name).write_text(Path(ligand).read_text())
        return MaterializeResult(
            input_dir=input_dir,
            materialized_count=len(ligands_available),
            failed_count=0,
            failed_examples=(),
        )

    def fake_scorch_command(_receptor, _ligands, _threads, _cfg):
        return ["python", "scorch.py", "--out", "{out}"]

    def fake_run_with_lock_fallback(*, cmd, run_kwargs, **_kwargs):
        captured["timeout"] = float(run_kwargs["timeout"])
        out_path = Path(cmd[cmd.index("--out") + 1])
        _write_csv(out_path, [{"Ligand_ID": "lig_000", "SCORCH_score": "1.0"}])
        return subprocess.CompletedProcess(cmd, 0, "", ""), False

    monkeypatch.setattr(scorch_stage_runner, "acquire_global_cores", fake_acquire_global_cores)
    monkeypatch.setattr(scorch_stage_runner, "run_with_lock_fallback", fake_run_with_lock_fallback)

    ok, out_path = scorch_stage_runner.score_stage(
        cfg,
        StageSpec("vina", "stage1", "scorch_scores.csv"),
        combo,
        run_root,
        post_root,
        receptor,
        threads=1,
        overwrite=False,
        logger=logging.getLogger("test.scorch.stage.timeout"),
        decoy_prefix="decoys",
        component="[test]",
        emit_task_event=lambda *_args, **_kwargs: None,
        materialize_inputs=fake_materialize_inputs,
        scorch_command=fake_scorch_command,
        annotate_csv=lambda *_args, **_kwargs: AnnotateResult(True, False),
        scorch_root=None,
        scorch_python=None,
        scorch_script=None,
        allowed_bases=None,
        control_bases=set(),
        run_mode="fda",
    )

    assert ok is True
    assert out_path is not None
    assert captured["timeout"] == 600.0


def test_score_stage_honors_explicit_timeout(tmp_path: Path, monkeypatch) -> None:
    cfg: dict[str, object] = {
        "OVERALL_DIR": str(tmp_path),
        "SCORCH_SHARD_TIMEOUT_SEC": "42",
    }
    run_root = tmp_path / "docked" / "rid"
    post_root = tmp_path / "post" / "rid"
    combo = ("PDB1", "HOLO", "ph_7_0")
    stage_root = run_root / combo[0] / combo[1] / combo[2] / "stage1"
    stage_root.mkdir(parents=True)
    (stage_root / "lig_stage1.pdbqt").write_text("MODEL\nENDMDL\n")
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("RECEPTOR\n")
    captured: dict[str, float] = {}

    @contextmanager
    def fake_acquire_global_cores(*_args, cores=1, **_kwargs):
        yield int(cores)

    def fake_run_with_lock_fallback(*, cmd, run_kwargs, **_kwargs):
        captured["timeout"] = float(run_kwargs["timeout"])
        out_path = Path(cmd[cmd.index("--out") + 1])
        _write_csv(out_path, [{"Ligand_ID": "lig", "SCORCH_score": "1.0"}])
        return subprocess.CompletedProcess(cmd, 0, "", ""), False

    monkeypatch.setattr(scorch_stage_runner, "acquire_global_cores", fake_acquire_global_cores)
    monkeypatch.setattr(scorch_stage_runner, "run_with_lock_fallback", fake_run_with_lock_fallback)

    ok, _out_path = scorch_stage_runner.score_stage(
        cfg,
        StageSpec("vina", "stage1", "scorch_scores.csv"),
        combo,
        run_root,
        post_root,
        receptor,
        threads=1,
        overwrite=False,
        logger=logging.getLogger("test.scorch.stage.timeout.override"),
        decoy_prefix="decoys",
        component="[test]",
        emit_task_event=lambda *_args, **_kwargs: None,
        materialize_inputs=lambda _ph_root, combo_post_root, _stage_dir, ligands, *_args, **_kwargs: MaterializeResult(
            input_dir=combo_post_root / ".inputs",
            materialized_count=len(ligands),
            failed_count=0,
            failed_examples=(),
        ),
        scorch_command=lambda *_args, **_kwargs: ["python", "scorch.py", "--out", "{out}"],
        annotate_csv=lambda *_args, **_kwargs: AnnotateResult(True, False),
        scorch_root=None,
        scorch_python=None,
        scorch_script=None,
        allowed_bases=None,
        control_bases=set(),
        run_mode="fda",
    )

    assert ok is True
    assert captured["timeout"] == 42.0
