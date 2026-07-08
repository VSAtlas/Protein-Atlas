from __future__ import annotations

import logging
import subprocess
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

import cli.postrun_hooks_scorch as scorch_hook
import cli.postrun_hooks_scorch_queue as scorch_queue
import cli.postrun_hooks_runtime as postrun_hooks
import cli.postrun_hooks_support as postrun_support
from cli.distributed_context import (
    resolve_distributed_context,
    scorch_scope_key,
    try_claim_scorch_scope,
)
import post_docking.rescoring.rescoring_scorch as rescoring_scorch
import post_docking.rescoring.scorch_housekeeping as scorch_housekeeping
import post_docking.rescoring.scorch_runtime_state as scorch_runtime_state
from post_docking.rescoring.scorch_device import scorch_child_env_for


def _reset_scorch_globals() -> None:
    scorch_runtime_state.SCORCH_SCRIPT = None
    scorch_runtime_state.SCORCH_ENV = "scorch-env"
    scorch_runtime_state.SCORCH_ENV_PREFIX = None
    scorch_runtime_state.SCORCH_ROOT = None
    scorch_runtime_state.SCORCH_PYTHON = None
    scorch_runtime_state.SCORCH_USE_MICROMAMBA = False


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
    assert str(script_path.resolve()) in cmd
    assert "None" not in cmd
    assert "-n" not in cmd


def test_preflight_uses_prefix_python_when_available(
    monkeypatch, tmp_path: Path
) -> None:
    _reset_scorch_globals()
    script_path = tmp_path / "scorch.py"
    script_path.write_text("print('ok')\n", encoding="utf-8")
    env_prefix = tmp_path / "scorch_env"
    env_python = env_prefix / "bin" / "python"
    env_python.parent.mkdir(parents=True, exist_ok=True)
    env_python.write_text("#!/bin/sh\n", encoding="utf-8")
    env_python.chmod(0o755)

    monkeypatch.setattr(
        rescoring_scorch.shutil,
        "which",
        lambda name: "/usr/bin/micromamba" if name == "micromamba" else None,
    )

    cfg = {"SCORCH": str(script_path), "SCORCH_ENV_PREFIX": str(env_prefix)}
    logger = logging.getLogger("test_scorch_env_prefix_python")

    assert rescoring_scorch._preflight(cfg, logger) is True

    cmd = rescoring_scorch._scorch_command(
        Path("receptor.pdbqt"), Path("ligands.sdf"), threads=1
    )
    assert cmd[:2] == [str(env_python.resolve()), str(script_path.resolve())]
    assert "micromamba" not in cmd


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
    assert str(script_path.resolve()) in cmd
    assert "None" not in cmd
    assert "-p" not in cmd


def test_preflight_uses_current_python_when_already_in_scorch_env(
    monkeypatch, tmp_path: Path
) -> None:
    _reset_scorch_globals()
    script_path = tmp_path / "scorch.py"
    script_path.write_text("print('ok')\n", encoding="utf-8")
    env_prefix = tmp_path / "scorch-env"
    env_python = env_prefix / "bin" / "python"
    env_python.parent.mkdir(parents=True, exist_ok=True)
    env_python.write_text("#!/bin/sh\n", encoding="utf-8")
    env_python.chmod(0o755)

    monkeypatch.setattr(
        rescoring_scorch.shutil,
        "which",
        lambda name: "/usr/bin/micromamba" if name == "micromamba" else None,
    )
    monkeypatch.setattr(scorch_runtime_state.sys, "executable", str(env_python))

    cfg = {"SCORCH": str(script_path), "SCORCH_ENV_PREFIX": str(env_prefix)}
    logger = logging.getLogger("test_scorch_current_env")

    assert rescoring_scorch._preflight(cfg, logger) is True

    cmd = rescoring_scorch._scorch_command(
        Path("receptor.pdbqt"), Path("ligands.sdf"), threads=1
    )
    assert cmd[:2] == [str(env_python.resolve()), str(script_path.resolve())]
    assert "micromamba" not in cmd
    assert "None" not in cmd


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


def test_main_applies_repo_scorch_default_to_cfg_override(
    monkeypatch, tmp_path: Path
) -> None:
    repo_root = tmp_path / "atlas" / "code" / "protein_automation"
    script_path = tmp_path / "atlas" / "tools" / "SCORCH" / "scorch.py"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text("print('ok')\n", encoding="utf-8")
    repo_root.mkdir(parents=True, exist_ok=True)
    captured_cfg: dict[str, object] = {}

    def _capture_preflight(cfg, logger):
        del logger
        captured_cfg.update(cfg)
        return False

    monkeypatch.delenv("SCORCH", raising=False)
    monkeypatch.delenv("SCORCH_SCRIPT", raising=False)
    monkeypatch.delenv("SCORCH_ENV", raising=False)
    monkeypatch.setattr(rescoring_scorch, "_preflight", _capture_preflight)

    rc = rescoring_scorch.main(
        ["--run-id", "run_x", "--repo-root", str(repo_root)],
        cfg_override={"USE_SCORCH": True},
    )

    assert rc == 1
    assert captured_cfg["SCORCH"] == str(script_path.resolve())
    assert captured_cfg["SCORCH_ENV"] == "scorch-env"


def test_per_pdb_hook_uses_elastic_jobs_and_threads(monkeypatch) -> None:
    captured: list[list[str]] = []

    def _fake_run(cmd, cwd=None, check=False, env=None):
        del cwd, check, env
        cmd_list = [str(x) for x in cmd]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )

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
    jobs = int(cmd[cmd.index("--jobs") + 1])
    threads = int(cmd[cmd.index("--threads") + 1])
    assert jobs >= 1
    assert threads >= 1
    assert jobs * threads <= 12
    assert threads > 1
    assert cmd[cmd.index("--docked-root") + 1] == "/tmp/custom_docked"
    assert cmd[cmd.index("--post-docked-root") + 1] == "/tmp/custom_post_docked"


def test_per_pdb_hook_forwards_runtime_scorch_config_to_subprocess(
    monkeypatch,
) -> None:
    captured_env: list[dict[str, str]] = []

    def _fake_run(cmd, cwd=None, check=False, env=None):
        del cwd, check
        cmd_list = [str(x) for x in cmd]
        captured_env.append(dict(env or {}))
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )

    cfg = {
        "USE_SCORCH": True,
        "TEST_MODE_ENABLE": "dud",
        "SCORCH_TOP_FRACTION": 0.1,
        "CPU": 32,
        "DOCKED_DIR": "/work/docked",
        "POST_DOCKED_DIR": "/work/post_docked",
    }

    rc = postrun_hooks._maybe_run_scorch_rescore_for_pdb(
        cfg,
        "run_x",
        "TEST",
        variant="HOLO",
        ph="pH7_0",
    )

    assert rc == 0
    assert captured_env
    env = captured_env[0]
    assert env["USE_SCORCH"] == "true"
    assert env["TEST_MODE_ENABLE"] == "dud"
    assert env["SCORCH_TOP_FRACTION"] == "0.1"
    assert env["CPU"] == "32"
    assert env["DOCKED_DIR"] == "/work/docked"
    assert env["POST_DOCKED_DIR"] == "/work/post_docked"


def test_scorch_subprocess_env_drops_stale_provisional_flags(monkeypatch) -> None:
    monkeypatch.setenv("SCORCH_PROVISIONAL_ENABLE", "true")
    monkeypatch.setenv("SCORCH_PROVISIONAL_CACHE_WRITE", "true")
    monkeypatch.setenv("SCORCH_PROVISIONAL_REUSE", "true")

    env = postrun_hooks._scorch_subprocess_env({"USE_SCORCH": True})

    assert env["USE_SCORCH"] == "true"
    assert "SCORCH_PROVISIONAL_ENABLE" not in env
    assert "SCORCH_PROVISIONAL_CACHE_WRITE" not in env
    assert "SCORCH_PROVISIONAL_REUSE" not in env


def test_prep_for_scorch_forwards_worker_count(tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def _fake_run(cmd, capture_output=False, text=False):
        del capture_output, text
        cmd_list = [str(part) for part in cmd]
        captured["cmd"] = cmd_list
        return subprocess.CompletedProcess(cmd_list, 0, "", "")

    ok = scorch_housekeeping.run_prep_for_scorch(
        "run_workers",
        tmp_path,
        tmp_path / "docked",
        tmp_path / "post_docked",
        False,
        logging.getLogger("test_prep_workers"),
        component="[test]",
        sys_executable="python",
        run_subprocess=_fake_run,
        max_workers=5,
    )

    assert ok
    cmd = captured["cmd"]
    assert cmd[cmd.index("--workers") + 1] == "5"


def test_scorch_queue_max_workers_env_caps_default(monkeypatch) -> None:
    monkeypatch.setenv("ATLAS_SCORCH_QUEUE_MAX_WORKERS", "1")

    service = postrun_hooks._start_scorch_mixed_queue_service(
        {"USE_SCORCH": True, "CPU": 112},
        "run_x",
    )

    assert service is not None
    summary = service.close_and_wait()
    assert summary["worker_cap"] == 1


def test_scorch_queue_cpu_shard_mode_scales_outer_cap_on_large_nodes(monkeypatch) -> None:
    monkeypatch.delenv("ATLAS_SCORCH_QUEUE_MAX_WORKERS", raising=False)
    monkeypatch.delenv("SCORCH_QUEUE_MAX_WORKERS", raising=False)
    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")

    service = postrun_hooks._start_scorch_mixed_queue_service(
        {"USE_SCORCH": True, "CPU": 112, "ATLAS_SCORCH_SHARDS_ENABLE": True},
        "run_x",
    )

    assert service is not None
    summary = service.close_and_wait()
    assert summary["worker_cap"] == 8


def test_scorch_queue_cpu_shard_mode_avoids_full_node_tail_underutilization(
    monkeypatch,
) -> None:
    monkeypatch.delenv("ATLAS_SCORCH_QUEUE_MAX_WORKERS", raising=False)
    monkeypatch.delenv("SCORCH_QUEUE_MAX_WORKERS", raising=False)
    monkeypatch.delenv("ATLAS_SCORCH_QUEUE_MIN_CPU_PER_WORKER", raising=False)
    monkeypatch.delenv("SCORCH_QUEUE_MIN_CPU_PER_WORKER", raising=False)
    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    captured: list[dict[str, object]] = []
    release = threading.Event()

    def _fake_rescore(*_args, **kwargs):
        captured.append(dict(kwargs))
        release.wait(timeout=2.0)
        return 0

    monkeypatch.setattr(scorch_queue, "_maybe_run_scorch_rescore_for_pdb", _fake_rescore)
    cfg = {"USE_SCORCH": True, "CPU": 112, "ATLAS_SCORCH_SHARDS_ENABLE": True}
    service = postrun_hooks._start_scorch_mixed_queue_service(cfg, "run_util")

    assert service is not None
    for idx in range(12):
        assert service.submit(
            cfg=cfg,
            pdb_id=f"P{idx:03d}",
            variant="HOLO",
            ph="pH7_0",
        )

    deadline = time.time() + 2.0
    while time.time() < deadline and len(captured) < 8:
        time.sleep(0.01)
    release.set()
    summary = service.close_and_wait()

    assert summary["failed"] == 0
    assert summary["worker_cap"] == 8
    assert summary["workers_peak"] >= 8
    assert len(captured) >= 8
    assert max(int(item["cpu_budget_hint"]) for item in captured) >= 56
    assert max(int(item["free_cores_hint"]) for item in captured) >= 56
    assert min(int(item["cpu_budget_hint"]) for item in captured) <= 14


def test_scorch_queue_cpu_shard_mode_prevents_node_oversubscription(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.delenv("ATLAS_SCORCH_QUEUE_MAX_WORKERS", raising=False)
    monkeypatch.delenv("SCORCH_QUEUE_MAX_WORKERS", raising=False)
    monkeypatch.delenv("ATLAS_SCORCH_QUEUE_MIN_CPU_PER_WORKER", raising=False)
    monkeypatch.delenv("SCORCH_QUEUE_MIN_CPU_PER_WORKER", raising=False)
    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    captured: list[dict[str, object]] = []
    release = threading.Event()

    def _fake_rescore(*_args, **kwargs):
        captured.append(dict(kwargs))
        release.wait(timeout=2.0)
        return 0

    monkeypatch.setattr(scorch_queue, "_maybe_run_scorch_rescore_for_pdb", _fake_rescore)
    cfg = {"USE_SCORCH": True, "CPU": 112, "ATLAS_SCORCH_SHARDS_ENABLE": True}
    service = postrun_hooks._start_scorch_mixed_queue_service(cfg, "run_no_oversub")

    assert service is not None
    for idx in range(12):
        assert service.submit(
            cfg=cfg,
            pdb_id=f"P{idx:03d}",
            variant="HOLO",
            ph="pH7_0",
            allowed_count_hint=512,
        )

    deadline = time.time() + 2.0
    while time.time() < deadline and len(captured) < 8:
        time.sleep(0.01)
    release.set()
    summary = service.close_and_wait()

    assert summary["failed"] == 0
    assert captured
    max_wrappers = int(summary["worker_cap"])
    budgets = [int(item["cpu_budget_hint"]) for item in captured]
    max_budget = max(budgets)
    assert max_wrappers >= 4
    assert max_budget >= 56
    assert min(budgets) <= 14

    launched: list[list[str]] = []
    docked_dir = tmp_path / "docked"
    (docked_dir / "run_no_oversub").mkdir(parents=True)

    def _fake_run_scorch_hook_command(**kwargs):
        cmd_list = [str(x) for x in kwargs["cmd"]]
        launched.append(cmd_list)
        return subprocess.CompletedProcess(kwargs["cmd"], 0, "", "")

    monkeypatch.setattr(
        scorch_hook,
        "_run_scorch_hook_command",
        _fake_run_scorch_hook_command,
    )
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )

    rc = postrun_hooks._maybe_run_scorch_rescore_for_pdb(
        {
            "USE_SCORCH": True,
            "CPU": 112,
            "ATLAS_SCORCH_SHARDS_ENABLE": True,
            "DOCKED_DIR": str(docked_dir),
        },
        "run_no_oversub",
        "TEST",
        allowed_count_hint=512,
        free_cores_hint=max_budget,
        cpu_budget_hint=max_budget,
        scheduler_queue_depth_hint=0,
        backlog_hint=32,
        profile_hint="aggressive",
    )

    assert rc == 0
    assert launched
    cmd = launched[0]
    jobs = int(cmd[cmd.index("--jobs") + 1])
    threads = int(cmd[cmd.index("--threads") + 1])
    assert jobs * threads <= max(1, max_budget // 2)
    assert jobs * threads >= max(16, max_budget // 4)


def test_scorch_queue_reaps_done_futures_before_dedupe(monkeypatch) -> None:
    monkeypatch.setattr(
        scorch_queue,
        "_maybe_run_scorch_rescore_for_pdb",
        lambda *_args, **_kwargs: 0,
    )
    service = postrun_hooks.ScorchMixedQueueService(
        run_id="run_reap",
        worker_cap=1,
        total_cpu=4,
        execution_mode="subprocess",
        logger=logging.getLogger("test_scorch_queue_reap"),
        verbose=False,
    )

    assert service.submit(cfg={}, pdb_id="ABCD", variant="HOLO", ph="pH7_0")
    deadline = time.time() + 2.0
    while time.time() < deadline:
        if all(future.done() for future in service._futures):
            break
        time.sleep(0.01)
    assert service.submit(cfg={}, pdb_id="ABCD", variant="HOLO", ph="pH7_0")
    summary = service.close_and_wait()

    assert summary["submitted"] == 2
    assert summary["failed"] == 0


def test_scorch_queue_counts_none_return_as_failure(monkeypatch) -> None:
    monkeypatch.setattr(
        scorch_queue,
        "_maybe_run_scorch_rescore_for_pdb",
        lambda *_args, **_kwargs: None,
    )
    service = postrun_hooks.ScorchMixedQueueService(
        run_id="run_none_failure",
        worker_cap=1,
        total_cpu=4,
        execution_mode="subprocess",
        logger=logging.getLogger("test_scorch_queue_none_failure"),
        verbose=False,
    )

    assert service.submit(cfg={}, pdb_id="ABCD", variant="HOLO", ph="pH7_0")
    summary = service.close_and_wait()

    assert summary["submitted"] == 1
    assert summary["failed"] == 1


def test_scorch_queue_slices_cpu_budget_across_target_workers(monkeypatch) -> None:
    monkeypatch.delenv("ATLAS_SCORCH_QUEUE_MIN_CPU_PER_WORKER", raising=False)
    captured: list[dict[str, object]] = []
    release = threading.Event()

    def _fake_rescore(*_args, **kwargs):
        captured.append(dict(kwargs))
        release.wait(timeout=2.0)
        return 0

    monkeypatch.setattr(scorch_queue, "_maybe_run_scorch_rescore_for_pdb", _fake_rescore)
    service = postrun_hooks.ScorchMixedQueueService(
        run_id="run_budget",
        worker_cap=8,
        total_cpu=48,
        execution_mode="subprocess",
        logger=logging.getLogger("test_scorch_queue_budget"),
        verbose=False,
    )

    cfg = {"USE_SCORCH": True, "CPU": 48}
    assert service.submit(cfg=cfg, pdb_id="AAAA", variant="HOLO", ph="pH7_0")
    assert service.submit(cfg=cfg, pdb_id="BBBB", variant="HOLO", ph="pH7_0")
    assert service.submit(cfg=cfg, pdb_id="CCCC", variant="HOLO", ph="pH7_0")
    release.set()
    summary = service.close_and_wait()

    assert summary["failed"] == 0
    assert summary["submitted"] == 3
    assert captured
    assert max(int(item["cpu_budget_hint"]) for item in captured) <= 24
    assert max(int(item["free_cores_hint"]) for item in captured) <= 24


def test_scorch_hook_treats_cpu_budget_hint_as_cap(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: list[list[str]] = []
    docked_dir = tmp_path / "docked"
    (docked_dir / "run_budget_cap").mkdir(parents=True)

    def _fake_run_scorch_hook_command(**kwargs):
        cmd_list = [str(x) for x in kwargs["cmd"]]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(kwargs["cmd"], 0, "", "")

    monkeypatch.setattr(
        scorch_hook,
        "_run_scorch_hook_command",
        _fake_run_scorch_hook_command,
    )
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )

    rc = postrun_hooks._maybe_run_scorch_rescore_for_pdb(
        {
            "USE_SCORCH": True,
            "CPU": 48,
            "ATLAS_SCORCH_SHARDS_ENABLE": True,
            "DOCKED_DIR": str(docked_dir),
        },
        "run_budget_cap",
        "TEST",
        allowed_count_hint=512,
        free_cores_hint=24,
        cpu_budget_hint=24,
        scheduler_queue_depth_hint=0,
        backlog_hint=32,
        profile_hint="aggressive",
    )

    assert rc == 0
    assert captured
    cmd = captured[0]
    jobs = int(cmd[cmd.index("--jobs") + 1])
    threads = int(cmd[cmd.index("--threads") + 1])
    assert jobs * threads <= 24


def test_distributed_scorch_scope_claim_blocks_duplicate_hook_launch(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cfg = {
        "USE_SCORCH": True,
        "CPU": 4,
        "ATLAS_SCORCH_SHARDS_ENABLE": False,
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "POST_DOCKED_DIR": str(tmp_path / "post_docked"),
    }
    run_id = "run_claim"
    scope_key = scorch_scope_key(
        run_id=run_id,
        pdb_id="BNNQ",
        variant="HOLO",
        ph="pH5_6",
        phase="final",
    )
    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    monkeypatch.setenv("ATLAS_DIST_TASK_COUNT", "2")
    monkeypatch.setenv("ATLAS_DIST_TASK_MIN", "0")
    monkeypatch.setenv("ATLAS_DIST_TASK_ID", "0")
    dist_ctx0 = resolve_distributed_context(cfg, run_id)

    assert try_claim_scorch_scope(
        dist_ctx0,
        cfg,
        scope_key=scope_key,
        pdb_id="BNNQ",
        variant="HOLO",
        ph="pH5_6",
        lease_sec=600.0,
    )

    launched: list[list[str]] = []

    def _fake_run_scorch_hook_command(**kwargs):
        launched.append([str(x) for x in kwargs["cmd"]])
        return subprocess.CompletedProcess(kwargs["cmd"], 0, "", "")

    monkeypatch.setattr(
        scorch_hook,
        "_run_scorch_hook_command",
        _fake_run_scorch_hook_command,
    )
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )
    monkeypatch.setenv("ATLAS_DIST_TASK_ID", "1")

    rc = scorch_hook._maybe_run_scorch_rescore_for_pdb(
        cfg,
        run_id,
        "BNNQ",
        variant="HOLO",
        ph="pH5_6",
        allowed_count_hint=10,
    )

    assert rc is None
    assert launched == []


def test_distributed_shard_mode_skips_coarse_scorch_scope_claim(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cfg = {
        "USE_SCORCH": True,
        "CPU": 4,
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "POST_DOCKED_DIR": str(tmp_path / "post_docked"),
    }
    run_id = "run_sharded_claim"
    scope_key = scorch_scope_key(
        run_id=run_id,
        pdb_id="BNNQ",
        variant="HOLO",
        ph="pH5_6",
        phase="final",
    )
    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    monkeypatch.setenv("ATLAS_DIST_TASK_COUNT", "2")
    monkeypatch.setenv("ATLAS_DIST_TASK_MIN", "0")
    monkeypatch.setenv("ATLAS_DIST_TASK_ID", "0")
    dist_ctx0 = resolve_distributed_context(cfg, run_id)
    assert try_claim_scorch_scope(
        dist_ctx0,
        cfg,
        scope_key=scope_key,
        pdb_id="BNNQ",
        variant="HOLO",
        ph="pH5_6",
        lease_sec=600.0,
    )

    launched: list[list[str]] = []

    def _fake_run_scorch_hook_command(**kwargs):
        launched.append([str(x) for x in kwargs["cmd"]])
        return subprocess.CompletedProcess(kwargs["cmd"], 0, "", "")

    monkeypatch.setattr(
        scorch_hook,
        "_run_scorch_hook_command",
        _fake_run_scorch_hook_command,
    )
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )
    monkeypatch.setenv("ATLAS_DIST_TASK_ID", "1")

    rc = scorch_hook._maybe_run_scorch_rescore_for_pdb(
        cfg,
        run_id,
        "BNNQ",
        variant="HOLO",
        ph="pH5_6",
        allowed_count_hint=10,
    )

    assert rc == 0
    assert launched
    assert "--pdb-id" in launched[0]


def test_distributed_cpu_shard_mode_defaults_to_productive_scorch_process_count(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cfg = {
        "USE_SCORCH": True,
        "CPU": 112,
        "MANIFESTS_DIR": str(tmp_path / "manifests"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "POST_DOCKED_DIR": str(tmp_path / "post_docked"),
    }
    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    monkeypatch.setenv("ATLAS_DIST_TASK_COUNT", "3")
    monkeypatch.delenv("ATLAS_SCORCH_SHARD_LOCAL_WORKERS", raising=False)
    monkeypatch.delenv("SCORCH_SHARD_LOCAL_WORKERS", raising=False)
    monkeypatch.delenv("ATLAS_SCORCH_SHARD_THREADS", raising=False)
    monkeypatch.delenv("SCORCH_SHARD_THREADS", raising=False)
    monkeypatch.delenv("ATLAS_SCORCH_HOOK_TIMEOUT_SEC", raising=False)
    monkeypatch.delenv("SCORCH_HOOK_TIMEOUT_SEC", raising=False)

    launched: list[tuple[list[str], dict[str, str], dict[str, object]]] = []

    def _fake_run_scorch_hook_command(**kwargs):
        launched.append(
            (
                [str(x) for x in kwargs["cmd"]],
                dict(kwargs.get("env") or {}),
                dict(kwargs.get("cfg") or {}),
            )
        )
        return subprocess.CompletedProcess(kwargs["cmd"], 0, "", "")

    monkeypatch.setattr(
        scorch_hook,
        "_run_scorch_hook_command",
        _fake_run_scorch_hook_command,
    )
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )

    rc = scorch_hook._maybe_run_scorch_rescore_for_pdb(
        cfg,
        "run_safe_shards",
        "BNNQ",
        allowed_count_hint=512,
        backlog_hint=30,
        free_cores_hint=112,
        cpu_budget_hint=112,
        profile_hint="aggressive",
    )

    assert rc == 0
    assert launched
    cmd, env, hook_cfg = launched[0]
    jobs = int(cmd[cmd.index("--jobs") + 1])
    threads = int(cmd[cmd.index("--threads") + 1])
    assert jobs == 32
    assert threads == 1
    assert env["TF_NUM_INTRAOP_THREADS"] == "1"
    assert env["TF_NUM_INTEROP_THREADS"] == "1"
    assert float(hook_cfg["ATLAS_SCORCH_HOOK_TIMEOUT_SEC"]) == 5400.0


def test_scorch_child_env_caps_thread_pools_by_default() -> None:
    env = scorch_child_env_for(
        effective="cpu",
        backend="none",
        gpu_ids=(),
        binding="single",
        task_id=None,
        base_env={},
    )

    assert env["OMP_NUM_THREADS"] == "1"
    assert env["OPENBLAS_NUM_THREADS"] == "1"
    assert env["MKL_NUM_THREADS"] == "1"
    assert env["NUMEXPR_NUM_THREADS"] == "1"
    assert env["TF_NUM_INTRAOP_THREADS"] == "1"
    assert env["TF_NUM_INTEROP_THREADS"] == "1"
    assert env["CUDA_VISIBLE_DEVICES"] == "-1"


def test_scorch_child_env_preserves_explicit_thread_pool_overrides() -> None:
    env = scorch_child_env_for(
        effective="cpu",
        backend="none",
        gpu_ids=(),
        binding="single",
        task_id=None,
        base_env={"OMP_NUM_THREADS": "8", "TF_NUM_INTRAOP_THREADS": "2"},
    )

    assert env["OMP_NUM_THREADS"] == "8"
    assert env["TF_NUM_INTRAOP_THREADS"] == "2"
    assert env["TF_NUM_INTEROP_THREADS"] == "1"


def test_resolve_scorch_execution_mode_uses_subprocess_when_scorch_env_differs(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.delenv("ATLAS_SCORCH_PERSISTENT_WORKER", raising=False)
    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    env_prefix = tmp_path / "scorch-env"
    (env_prefix / "bin").mkdir(parents=True)
    (env_prefix / "bin" / "python").write_text("#!/bin/sh\n", encoding="utf-8")
    assert (
        postrun_hooks._resolve_scorch_execution_mode(
            {"SCORCH_ENV_PREFIX": str(env_prefix)}
        )
        == "subprocess"
    )


def test_resolve_scorch_execution_mode_env_can_disable(monkeypatch) -> None:
    monkeypatch.setenv("ATLAS_SCORCH_PERSISTENT_WORKER", "0")
    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    assert postrun_hooks._resolve_scorch_execution_mode({}) == "subprocess"


def test_resolve_scorch_execution_mode_env_true_still_requires_same_env(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setenv("ATLAS_SCORCH_PERSISTENT_WORKER", "1")
    monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
    env_prefix = tmp_path / "other-scorch-env"
    (env_prefix / "bin").mkdir(parents=True)
    (env_prefix / "bin" / "python").write_text("#!/bin/sh\n", encoding="utf-8")
    assert (
        postrun_hooks._resolve_scorch_execution_mode(
            {"SCORCH_ENV_PREFIX": str(env_prefix)}
        )
        == "subprocess"
    )


def test_scorch_cmd_prefix_uses_configured_env_python(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.delenv("PYTEST_CURRENT_TEST", raising=False)
    env_prefix = tmp_path / "scorch-env"
    env_python = env_prefix / "bin" / "python"
    env_python.parent.mkdir(parents=True)
    env_python.write_text("#!/bin/sh\n", encoding="utf-8")
    env_python.chmod(0o755)

    cmd = postrun_hooks._scorch_cmd_prefix(
        tmp_path / "rescoring_scorch.py",
        {"SCORCH_ENV_PREFIX": str(env_prefix)},
    )

    assert cmd[:3] == [
        str(env_python.resolve()),
        "-m",
        "post_docking.rescoring.rescoring_scorch",
    ]


def test_scorch_cmd_prefix_requires_env_outside_pytest(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.delenv("PYTEST_CURRENT_TEST", raising=False)
    monkeypatch.delenv("SCORCH_ENV_PREFIX", raising=False)
    monkeypatch.delenv("SCORCH_ENV", raising=False)
    monkeypatch.delenv("ATLAS_SCORCH_ALLOW_PYTHON_FALLBACK", raising=False)
    monkeypatch.setattr(postrun_support.shutil, "which", lambda _name: None)

    with pytest.raises(RuntimeError, match="no usable SCORCH env"):
        postrun_hooks._scorch_cmd_prefix(tmp_path / "rescoring_scorch.py", {})


def test_scorch_cmd_prefix_defaults_to_named_scorch_env(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.delenv("SCORCH_ENV_PREFIX", raising=False)
    monkeypatch.delenv("SCORCH_ENV", raising=False)
    monkeypatch.setattr(
        postrun_support.shutil,
        "which",
        lambda name: "/usr/bin/micromamba" if name == "micromamba" else None,
    )

    cmd = postrun_hooks._scorch_cmd_prefix(tmp_path / "rescoring_scorch.py", {})

    assert cmd[:4] == ["/usr/bin/micromamba", "run", "-n", "scorch-env"]


def test_per_pdb_hook_uses_adaptive_shard_caps(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: list[list[str]] = []
    docked_dir = tmp_path / "docked"
    (docked_dir / "run_x").mkdir(parents=True)

    def _fake_run_scorch_hook_command(**kwargs):
        cmd_list = [str(x) for x in kwargs["cmd"]]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(kwargs["cmd"], 0, "", "")

    monkeypatch.setattr(
        scorch_hook,
        "_run_scorch_hook_command",
        _fake_run_scorch_hook_command,
    )
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )

    rc = postrun_hooks._maybe_run_scorch_rescore_for_pdb(
        {
            "USE_SCORCH": True,
            "CPU": 112,
            "ATLAS_SCORCH_SHARDS_ENABLE": True,
            "DOCKED_DIR": str(docked_dir),
        },
        "run_x",
        "TEST",
        allowed_count_hint=128,
        free_cores_hint=112,
        scheduler_queue_depth_hint=0,
        backlog_hint=1,
        profile_hint="aggressive",
    )

    assert rc == 0
    assert captured
    cmd = captured[0]
    assert int(cmd[cmd.index("--jobs") + 1]) == 8
    assert int(cmd[cmd.index("--threads") + 1]) == 1


def test_per_pdb_hook_inprocess_mode_skips_subprocess(monkeypatch) -> None:
    invoked: list[list[str]] = []

    def _fake_main(argv=None, cfg_override=None):
        del cfg_override
        invoked.append([str(x) for x in list(argv or [])])
        return 0

    monkeypatch.setattr(rescoring_scorch, "main", _fake_main)
    monkeypatch.setattr(
        postrun_hooks.subprocess,
        "run",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("unexpected subprocess")),
    )
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )
    cfg = {
        "USE_SCORCH": True,
        "CPU": 8,
    }
    rc = postrun_hooks._maybe_run_scorch_rescore_for_pdb(
        cfg,
        "run_inprocess",
        "TEST",
        execution_mode="inprocess",
    )
    assert rc == 0
    assert invoked
    assert invoked[0][0] == "--run-id"


def test_per_pdb_hook_inprocess_argv_ignores_launcher_prefix(monkeypatch) -> None:
    invoked: list[list[str]] = []

    def _fake_main(argv=None, cfg_override=None):
        del cfg_override
        invoked.append([str(x) for x in list(argv or [])])
        return 0

    monkeypatch.setattr(rescoring_scorch, "main", _fake_main)
    monkeypatch.setattr(
        postrun_hooks,
        "_scorch_cmd_prefix",
        lambda _script, _cfg=None: [
            "micromamba",
            "run",
            "-p",
            "/tmp/scorch-env",
            "python",
            "-m",
            "post_docking.rescoring.rescoring_scorch",
        ],
    )
    monkeypatch.setattr(
        scorch_hook,
        "_validate_scorch_scope_coverage",
        lambda **_kwargs: True,
    )
    rc = postrun_hooks._maybe_run_scorch_rescore_for_pdb(
        {"USE_SCORCH": True, "CPU": 8},
        "run_inprocess",
        "TEST",
        execution_mode="inprocess",
    )
    assert rc == 0
    assert invoked
    assert invoked[0][0] == "--run-id"
    assert "-m" not in invoked[0]


def test_provisional_per_pdb_hook_uses_hidden_post_root_and_skips_manifest(
    monkeypatch,
    tmp_path: Path,
) -> None:
    captured_cmd: list[list[str]] = []
    manifest_calls: list[dict[str, object]] = []

    def _fake_run(cmd, cwd=None, check=False, env=None):
        del cwd, check, env
        cmd_list = [str(x) for x in cmd]
        captured_cmd.append(cmd_list)
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        postrun_hooks,
        "update_manifest_for_postprocessing",
        lambda *args, **kwargs: manifest_calls.append({"args": args, "kwargs": kwargs}),
    )

    cfg = {
        "USE_SCORCH": True,
        "CPU": 8,
        "OVERALL_DIR": str(tmp_path),
        "DATA_DIR": str(tmp_path / "outputs" / "data"),
        "DOCKED_DIR": str(tmp_path / "outputs" / "docked"),
        "POST_DOCKED_DIR": str(tmp_path / "outputs" / "post_docked"),
        "SCORCH_PROVISIONAL_ENABLE": True,
        "SCORCH_PROVISIONAL_CPU_BUDGET": 1,
    }

    rc = postrun_hooks._maybe_run_scorch_rescore_for_pdb(
        cfg,
        "run_prov",
        "TEST",
        execution_mode="subprocess",
        provisional=True,
    )

    assert rc == 0
    assert captured_cmd
    cmd = captured_cmd[0]
    post_root = Path(cmd[cmd.index("--post-docked-root") + 1])
    assert post_root == Path(cfg["POST_DOCKED_DIR"])
    assert "--skip-autofix" in cmd
    assert manifest_calls == []


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
