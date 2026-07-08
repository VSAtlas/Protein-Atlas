from __future__ import annotations

import argparse
from pathlib import Path

from cli import qol_cli


def _submit_args(**overrides: object) -> argparse.Namespace:
    values: dict[str, object] = {
        "bench2_canary": True,
        "array": qol_cli.SLURM_SUBMIT_DEFAULT_ARRAY,
        "cpus_per_task": qol_cli.SLURM_SUBMIT_DEFAULT_CPUS_PER_TASK,
        "script": qol_cli.SLURM_SUBMIT_DEFAULT_SCRIPT,
        "finalize_script": "tools/slurm/run_spr_finalize.sh",
        "partition": None,
        "account": None,
        "time": None,
        "job_name": qol_cli.SLURM_SUBMIT_DEFAULT_JOB_NAME,
        "run_id": "",
        "main_args": "",
        "all_dirs": "",
        "prepped_root": "",
        "dry_run": False,
        "no_finalize": False,
        "mail_user": "",
        "mail_type": "BEGIN,END,FAIL",
        "exclusive_node": False,
        "execution_mode": "",
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_slurm_submit_bench2_canary_sets_low_su_defaults() -> None:
    args = _submit_args()

    qol_cli._apply_slurm_submit_presets(args)

    assert args.array == qol_cli.SLURM_BENCH2_CANARY_ARRAY
    assert args.cpus_per_task == qol_cli.SLURM_BENCH2_CANARY_CPUS_PER_TASK
    assert args.script == qol_cli.SLURM_BENCH2_CANARY_SCRIPT
    assert args.time == qol_cli.SLURM_BENCH2_CANARY_TIME
    assert args.job_name == qol_cli.SLURM_BENCH2_CANARY_JOB_NAME
    assert args.main_args == qol_cli.SLURM_BENCH2_CANARY_MAIN_ARGS
    assert args.exclusive_node is True


def test_slurm_submit_bench2_canary_preserves_explicit_overrides() -> None:
    args = _submit_args(
        array="0-5%2",
        cpus_per_task=24,
        script="tools/slurm/custom.sh",
        time="00:10:00",
        job_name="custom_canary",
        main_args="-bench2 -fast --run-id {run_id} --custom",
        exclusive_node=True,
    )

    qol_cli._apply_slurm_submit_presets(args)

    assert args.array == "0-5%2"
    assert args.cpus_per_task == 24
    assert args.script == "tools/slurm/custom.sh"
    assert args.time == "00:10:00"
    assert args.job_name == "custom_canary"
    assert args.main_args == "-bench2 -fast --run-id {run_id} --custom"
    assert args.exclusive_node is True


def test_slurm_submit_presets_skip_when_canary_disabled() -> None:
    args = _submit_args(bench2_canary=False)

    qol_cli._apply_slurm_submit_presets(args)

    assert args.array == qol_cli.SLURM_SUBMIT_DEFAULT_ARRAY
    assert args.cpus_per_task == qol_cli.SLURM_SUBMIT_DEFAULT_CPUS_PER_TASK
    assert args.script == qol_cli.SLURM_SUBMIT_DEFAULT_SCRIPT
    assert args.time is None
    assert args.job_name == qol_cli.SLURM_SUBMIT_DEFAULT_JOB_NAME
    assert args.main_args == ""
    assert args.exclusive_node is False


def test_slurm_cpus_per_task_accepts_auto_for_whole_node() -> None:
    assert qol_cli._slurm_cpus_per_task_value("auto") is None
    assert qol_cli._slurm_cpus_per_task_value("whole-node") is None
    assert qol_cli._slurm_cpus_per_task_value("112") == 112


def test_slurm_submit_relocated_dry_run_exports_stable_repo_root(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    repo = tmp_path / "repo"
    slurm_dir = repo / "tools" / "slurm"
    slurm_dir.mkdir(parents=True)
    (slurm_dir / "run_relocated_spr.sh").write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    (slurm_dir / "run_relocated_finalize.sh").write_text(
        "#!/usr/bin/env bash\n", encoding="utf-8"
    )

    monkeypatch.setattr(qol_cli, "_repo_root", lambda: repo)
    args = _submit_args(
        bench2_canary=False,
        run_id="spd_resume",
        script="tools/slurm/run_relocated_spr.sh",
        finalize_script="tools/slurm/run_spr_finalize.sh",
        all_dirs=str(tmp_path / "scratch" / "spd_resume"),
        prepped_root=str(tmp_path / "scratch" / "prepped_ligands"),
        main_args="--resume --run-id {run_id}",
        array="0-9%10",
        cpus_per_task="auto",
        partition="spr",
        account="MCB24041",
        time="20:00:00",
        job_name="spd90_resume_scorch",
        dry_run=True,
        exclusive_node=True,
    )

    assert qol_cli._cmd_slurm_submit(args) == 0
    out = capsys.readouterr().out

    assert f"--chdir={repo}" in out
    assert f"ATLAS_REPO_ROOT={repo}" in out
    assert f"ATLAS_PREPPED_ROOT={tmp_path / 'scratch' / 'prepped_ligands'}" in out
    assert "ATLAS_USE_NODE_CPUS=1" in out
    finalizer_line = next(line for line in out.splitlines() if line.startswith("finalizer: "))
    assert f"ATLAS_PREPPED_ROOT={tmp_path / 'scratch' / 'prepped_ligands'}" in finalizer_line
    assert "tools/slurm/run_relocated_finalize.sh" in finalizer_line
    assert "tools/slurm/run_spr_finalize.sh" not in finalizer_line


def test_slurm_submit_relocated_auto_exports_scratch_prepped_root(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    repo = tmp_path / "repo"
    slurm_dir = repo / "tools" / "slurm"
    prepped_root = repo / "prepped_ligands"
    slurm_dir.mkdir(parents=True)
    prepped_root.mkdir(parents=True)
    (prepped_root / "_manifest.json").write_text("{}", encoding="utf-8")
    (slurm_dir / "run_relocated_spr.sh").write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    (slurm_dir / "run_relocated_finalize.sh").write_text(
        "#!/usr/bin/env bash\n", encoding="utf-8"
    )

    monkeypatch.setattr(qol_cli, "_repo_root", lambda: repo)
    args = _submit_args(
        bench2_canary=False,
        run_id="spd_resume",
        script="tools/slurm/run_relocated_spr.sh",
        finalize_script="tools/slurm/run_relocated_finalize.sh",
        all_dirs=str(tmp_path / "scratch" / "spd_resume"),
        main_args="--resume --run-id {run_id}",
        dry_run=True,
    )

    assert qol_cli._cmd_slurm_submit(args) == 0
    out = capsys.readouterr().out

    assert f"ATLAS_PREPPED_ROOT={prepped_root}" in out
