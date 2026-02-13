from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from types import SimpleNamespace

import postrun_hooks
import post_docking.rescoring.rescoring_scorch as rescoring_scorch


def _reset_scorch_globals() -> None:
    rescoring_scorch.SCORCH_SCRIPT = None
    rescoring_scorch.SCORCH_ENV = "scorch-env"
    rescoring_scorch.SCORCH_ENV_PREFIX = None
    rescoring_scorch.SCORCH_ROOT = None


def test_preflight_uses_env_prefix_when_configured(
    monkeypatch, tmp_path: Path
) -> None:
    _reset_scorch_globals()
    script_path = tmp_path / "scorch.py"
    script_path.write_text("print('ok')\n", encoding="utf-8")
    env_prefix = tmp_path / "scorch_env"
    env_prefix.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(
        rescoring_scorch.shutil,
        "which",
        lambda name: "/usr/bin/micromamba" if name == "micromamba" else None,
    )

    cfg = {
        "SCORCH": str(script_path),
        "SCORCH_ENV": "cfg-name-env",
        "SCORCH_ENV_PREFIX": str(env_prefix),
    }
    logger = logging.getLogger("test_scorch_env_prefix")

    assert rescoring_scorch._preflight(cfg, logger) is True

    cmd = rescoring_scorch._scorch_command(
        Path("receptor.pdbqt"), Path("ligands.sdf"), threads=1
    )
    assert cmd[:4] == ["micromamba", "run", "-p", str(env_prefix.resolve())]
    assert "-n" not in cmd


def test_preflight_uses_env_name_by_default(monkeypatch, tmp_path: Path) -> None:
    _reset_scorch_globals()
    script_path = tmp_path / "scorch.py"
    script_path.write_text("print('ok')\n", encoding="utf-8")

    monkeypatch.setattr(
        rescoring_scorch.shutil,
        "which",
        lambda name: "/usr/bin/micromamba" if name == "micromamba" else None,
    )

    cfg = {"SCORCH": str(script_path)}
    logger = logging.getLogger("test_scorch_env_name")

    assert rescoring_scorch._preflight(cfg, logger) is True

    cmd = rescoring_scorch._scorch_command(
        Path("receptor.pdbqt"), Path("ligands.sdf"), threads=1
    )
    assert cmd[:4] == ["micromamba", "run", "-n", "scorch-env"]
    assert "-p" not in cmd


def test_preflight_fails_when_env_prefix_path_missing(
    monkeypatch, tmp_path: Path
) -> None:
    _reset_scorch_globals()
    script_path = tmp_path / "scorch.py"
    script_path.write_text("print('ok')\n", encoding="utf-8")
    missing_prefix = tmp_path / "does_not_exist"

    monkeypatch.setattr(
        rescoring_scorch.shutil,
        "which",
        lambda name: "/usr/bin/micromamba" if name == "micromamba" else None,
    )

    cfg = {"SCORCH": str(script_path), "SCORCH_ENV_PREFIX": str(missing_prefix)}
    logger = logging.getLogger("test_scorch_env_prefix_missing")

    assert rescoring_scorch._preflight(cfg, logger) is False


def test_per_pdb_hook_uses_cpu_for_jobs_and_threads_one(monkeypatch) -> None:
    captured: list[list[str]] = []

    def _fake_run(cmd, cwd=None, check=False):
        del cwd, check
        cmd_list = [str(x) for x in cmd]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)

    cfg = {
        "USE_SCORCH": True,
        "CPU": 12,
        "SCORCH_THREADS": 99,
        "DOCKED_DIR": "/tmp/custom_docked",
        "POST_DOCKED_DIR": "/tmp/custom_post_docked",
    }
    postrun_hooks._maybe_run_scorch_rescore_for_pdb(cfg, "run_x", "TEST")

    assert captured
    cmd = captured[0]
    assert cmd[cmd.index("--jobs") + 1] == "12"
    assert cmd[cmd.index("--threads") + 1] == "1"
    assert cmd[cmd.index("--docked-root") + 1] == "/tmp/custom_docked"
    assert cmd[cmd.index("--post-docked-root") + 1] == "/tmp/custom_post_docked"


def test_resolve_roots_prefers_cli_over_cfg(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    cfg = {
        "DOCKED_DIR": str(tmp_path / "cfg_docked"),
        "POST_DOCKED_DIR": str(tmp_path / "cfg_post_docked"),
        "OUTPUT_DIR": str(tmp_path / "cfg_processed"),
    }

    args_cfg = SimpleNamespace(
        repo_root=str(repo_root),
        docked_root=None,
        post_docked_root=None,
    )
    resolved_repo, resolved_docked, resolved_post, resolved_processed = (
        rescoring_scorch._resolve_roots(args_cfg, cfg)
    )
    assert resolved_repo == repo_root.resolve()
    assert resolved_docked == Path(cfg["DOCKED_DIR"]).resolve()
    assert resolved_post == Path(cfg["POST_DOCKED_DIR"]).resolve()
    assert resolved_processed == Path(cfg["OUTPUT_DIR"]).resolve()

    args_cli = SimpleNamespace(
        repo_root=str(repo_root),
        docked_root=str(tmp_path / "cli_docked"),
        post_docked_root=str(tmp_path / "cli_post_docked"),
    )
    _, cli_docked, cli_post, cli_processed = rescoring_scorch._resolve_roots(
        args_cli, cfg
    )
    assert cli_docked == (tmp_path / "cli_docked").resolve()
    assert cli_post == (tmp_path / "cli_post_docked").resolve()
    assert cli_processed == Path(cfg["OUTPUT_DIR"]).resolve()
