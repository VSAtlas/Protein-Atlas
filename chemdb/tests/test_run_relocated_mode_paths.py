from __future__ import annotations

import os
import shutil
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import run_relocated_mode  # noqa: E402


def _write_prepped_payload(prepped_root: Path) -> None:
    lib_dir = prepped_root / "demo_lib"
    lib_dir.mkdir(parents=True, exist_ok=True)
    (lib_dir / "demo_1.pdbqt").write_text("MODEL 1\nENDMDL\n", encoding="utf-8")


def _require_tar_tools() -> None:
    missing = [tool for tool in ("tar", "unzstd") if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"requires tools on PATH: {', '.join(missing)}")


def _create_tar_zst(source_root: Path, archive_path: Path) -> None:
    commands = [
        [
            "tar",
            "--zstd",
            "-cf",
            str(archive_path),
            "-C",
            str(source_root),
            "prepped_ligands",
        ],
        [
            "tar",
            "--use-compress-program=zstd",
            "-cf",
            str(archive_path),
            "-C",
            str(source_root),
            "prepped_ligands",
        ],
    ]
    failures: list[str] = []
    for cmd in commands:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode == 0:
            return
        failures.append((proc.stderr or proc.stdout or "").strip())
    pytest.skip(
        "unable to create .tar.zst with local tar/zstd tooling: "
        + " | ".join(msg for msg in failures if msg)
    )


def _make_archive_with_intermediates_link(
    tmp_path: Path, link_target: str
) -> tuple[Path, Path]:
    source_root = tmp_path / "source_bundle"
    lib_dir = source_root / "prepped_ligands" / "3KFA"
    lib_dir.mkdir(parents=True, exist_ok=True)
    (lib_dir / "demo_1.pdbqt").write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    os.symlink(link_target, lib_dir / "intermediates_src")

    archive_root = tmp_path / "archive_bundle"
    archive_root.mkdir(parents=True, exist_ok=True)
    archive_path = archive_root / "prepped_ligands.tar.zst"
    _create_tar_zst(source_root, archive_path)
    return archive_root, archive_path


def test_extracted_root_optional_uses_config_resolution(
    monkeypatch, tmp_path: Path
) -> None:
    prepped_root = tmp_path / "prepped"
    _write_prepped_payload(prepped_root)

    monkeypatch.delenv("EXTRACTED_LIGANDS_DIR", raising=False)
    monkeypatch.delenv("DOCKED_DIR", raising=False)
    monkeypatch.delenv("POST_DOCKED_DIR", raising=False)

    captured: dict[str, object] = {}

    def _fake_run(cmd, cwd=None, env=None, check=False, **kwargs):
        del check, kwargs
        captured["cmd"] = [str(x) for x in cmd]
        captured["cwd"] = cwd
        captured["env"] = dict(env or {})
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(run_relocated_mode.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        run_relocated_mode,
        "_build_main_command",
        lambda repo_root, raw: [sys.executable, str(repo_root / "main.py"), "--noop"],
    )

    argv = [
        "run_relocated_mode.py",
        "--prepped-root",
        str(prepped_root),
        "--skip-extract",
        "--main-args",
        "--run-id reloc_test --pdb TEST --fast --test-fda",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    rc = run_relocated_mode.main()
    assert rc == 0

    env = captured["env"]
    assert env["PREPPED_LIGANDS_DIR"] == str(prepped_root.resolve())
    assert "EXTRACTED_LIGANDS_DIR" not in env
    assert "DOCKED_DIR" not in env
    assert "POST_DOCKED_DIR" not in env


def test_wrapper_sets_optional_roots_when_provided(
    monkeypatch, tmp_path: Path
) -> None:
    prepped_root = tmp_path / "prepped"
    extracted_root = tmp_path / "extracted"
    docked_root = tmp_path / "docked_out"
    post_docked_root = tmp_path / "post_docked_out"
    _write_prepped_payload(prepped_root)

    captured: dict[str, object] = {}

    def _fake_run(cmd, cwd=None, env=None, check=False, **kwargs):
        del check, kwargs
        captured["cmd"] = [str(x) for x in cmd]
        captured["cwd"] = cwd
        captured["env"] = dict(env or {})
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(run_relocated_mode.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        run_relocated_mode,
        "_build_main_command",
        lambda repo_root, raw: [sys.executable, str(repo_root / "main.py"), "--noop"],
    )

    argv = [
        "run_relocated_mode.py",
        "--prepped-root",
        str(prepped_root),
        "--extracted-root",
        str(extracted_root),
        "--docked-root",
        str(docked_root),
        "--post-docked-root",
        str(post_docked_root),
        "--skip-extract",
        "--main-args",
        "--run-id reloc_test --pdb TEST --fast --test-fda",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    rc = run_relocated_mode.main()
    assert rc == 0

    env = captured["env"]
    assert env["PREPPED_LIGANDS_DIR"] == str(prepped_root.resolve())
    assert env["EXTRACTED_LIGANDS_DIR"] == str(extracted_root.resolve())
    assert env["DOCKED_DIR"] == str(docked_root.resolve())
    assert env["POST_DOCKED_DIR"] == str(post_docked_root.resolve())

    assert extracted_root.exists()
    assert docked_root.exists()
    assert post_docked_root.exists()


def test_prepped_archive_extract_targets_prepped_root(
    monkeypatch, tmp_path: Path
) -> None:
    prepped_root = tmp_path / "prepped"
    archive_path = tmp_path / "prepped_ligands.tar.zst"
    archive_path.write_text("placeholder", encoding="utf-8")

    extracted_call: dict[str, object] = {}
    captured: dict[str, object] = {}

    def _fake_extract(archive: Path, dest: Path, force_extract: bool) -> bool:
        extracted_call["archive"] = archive
        extracted_call["dest"] = dest
        extracted_call["force_extract"] = force_extract
        _write_prepped_payload(dest)
        return True

    def _fake_run(cmd, cwd=None, env=None, check=False, **kwargs):
        del check, kwargs
        captured["cmd"] = [str(x) for x in cmd]
        captured["cwd"] = cwd
        captured["env"] = dict(env or {})
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(run_relocated_mode, "_extract_prepped_archive", _fake_extract)
    monkeypatch.setattr(run_relocated_mode.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        run_relocated_mode,
        "_build_main_command",
        lambda repo_root, raw: [sys.executable, str(repo_root / "main.py"), "--noop"],
    )

    argv = [
        "run_relocated_mode.py",
        "--prepped-root",
        str(prepped_root),
        "--prepped-archive",
        str(archive_path),
        "--main-args",
        "--run-id reloc_test --pdb TEST --fast --test-fda",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    rc = run_relocated_mode.main()
    assert rc == 0
    assert extracted_call["archive"] == archive_path.resolve()
    assert extracted_call["dest"] == prepped_root.resolve()
    assert extracted_call["force_extract"] is False
    assert captured["env"]["PREPPED_LIGANDS_DIR"] == str(prepped_root.resolve())


def test_all_dirs_sets_bundle_and_precreates_dirs(
    monkeypatch, tmp_path: Path
) -> None:
    all_root = tmp_path / "atlas2" / "protein_automation"
    captured: dict[str, object] = {}

    def _fake_run(cmd, cwd=None, env=None, check=False, **kwargs):
        del check, kwargs
        captured["cmd"] = [str(x) for x in cmd]
        captured["cwd"] = cwd
        captured["env"] = dict(env or {})
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(run_relocated_mode, "_has_prepped_payload", lambda _: True)
    monkeypatch.setattr(run_relocated_mode.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        run_relocated_mode,
        "_build_main_command",
        lambda repo_root, raw: [sys.executable, str(repo_root / "main.py"), "--noop"],
    )

    argv = [
        "run_relocated_mode.py",
        "--all-dirs",
        str(all_root),
        "--skip-extract",
        "--main-args",
        "--run-id reloc_test --pdb TEST --fast --test-fda",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    rc = run_relocated_mode.main()
    assert rc == 0

    env = captured["env"]
    assert env["PREPPED_LIGANDS_DIR"] == str((all_root / "prepped_ligands").resolve())
    assert env["EXTRACTED_LIGANDS_DIR"] == str(
        (all_root / "extracted_ligands").resolve()
    )
    assert env["DOCKED_DIR"] == str((all_root / "outputs" / "docked").resolve())
    assert env["POST_DOCKED_DIR"] == str(
        (all_root / "outputs" / "post_docked").resolve()
    )
    assert env["OUTPUT_DIR"] == str(
        (all_root / "outputs" / "processed_pdbs").resolve()
    )
    assert env["CONFIGS_DIR"] == str((all_root / "outputs" / "configs").resolve())
    assert env["MANIFESTS_DIR"] == str(
        (all_root / "outputs" / "manifests").resolve()
    )
    assert env["LOGS_DIR"] == str((all_root / "outputs" / "logs").resolve())
    assert env["DATA_DIR"] == str((all_root / "outputs" / "data").resolve())

    assert all_root.exists()
    assert (all_root / "prepped_ligands").exists()
    assert (all_root / "extracted_ligands").exists()
    assert (all_root / "outputs" / "docked").exists()
    assert (all_root / "outputs" / "post_docked").exists()
    assert (all_root / "outputs" / "processed_pdbs").exists()
    assert (all_root / "outputs" / "configs").exists()
    assert (all_root / "outputs" / "manifests").exists()
    assert (all_root / "outputs" / "logs").exists()
    assert (all_root / "outputs" / "data").exists()


def test_finalizer_manifest_dir_honors_relocated_manifest_root(tmp_path: Path) -> None:
    from tools import finalize_distributed_run

    scratch_manifest_root = tmp_path / "scratch" / "outputs" / "manifests"
    cfg = {
        "OVERALL_DIR": str(tmp_path / "work_repo"),
        "MANIFESTS_DIR": str(scratch_manifest_root),
    }

    assert finalize_distributed_run._manifest_run_dir(cfg, "reloc_final") == (
        scratch_manifest_root / "reloc_final"
    )


def test_explicit_roots_override_all_dirs_for_overlap(
    monkeypatch, tmp_path: Path
) -> None:
    all_root = tmp_path / "all"
    explicit_prepped = tmp_path / "explicit_prepped"
    explicit_docked = tmp_path / "explicit_docked"
    explicit_post = tmp_path / "explicit_post"
    _write_prepped_payload(explicit_prepped)

    captured: dict[str, object] = {}

    def _fake_run(cmd, cwd=None, env=None, check=False, **kwargs):
        del check, kwargs
        captured["cmd"] = [str(x) for x in cmd]
        captured["cwd"] = cwd
        captured["env"] = dict(env or {})
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(run_relocated_mode.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        run_relocated_mode,
        "_build_main_command",
        lambda repo_root, raw: [sys.executable, str(repo_root / "main.py"), "--noop"],
    )

    argv = [
        "run_relocated_mode.py",
        "--all-dirs",
        str(all_root),
        "--prepped-root",
        str(explicit_prepped),
        "--docked-root",
        str(explicit_docked),
        "--post-docked-root",
        str(explicit_post),
        "--skip-extract",
        "--main-args",
        "--run-id reloc_test --pdb TEST --fast --test-fda",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    rc = run_relocated_mode.main()
    assert rc == 0

    env = captured["env"]
    assert env["PREPPED_LIGANDS_DIR"] == str(explicit_prepped.resolve())
    assert env["DOCKED_DIR"] == str(explicit_docked.resolve())
    assert env["POST_DOCKED_DIR"] == str(explicit_post.resolve())
    assert env["EXTRACTED_LIGANDS_DIR"] == str((all_root / "extracted_ligands").resolve())
    assert env["OUTPUT_DIR"] == str(
        (all_root / "outputs" / "processed_pdbs").resolve()
    )
    assert env["CONFIGS_DIR"] == str((all_root / "outputs" / "configs").resolve())
    assert env["MANIFESTS_DIR"] == str(
        (all_root / "outputs" / "manifests").resolve()
    )
    assert env["LOGS_DIR"] == str((all_root / "outputs" / "logs").resolve())
    assert env["DATA_DIR"] == str((all_root / "outputs" / "data").resolve())


def test_default_archive_used_when_not_specified(
    monkeypatch, tmp_path: Path
) -> None:
    prepped_root = tmp_path / "prepped"
    default_archive = tmp_path / "repo_default_prepped.tar.zst"
    default_archive.write_text("placeholder", encoding="utf-8")

    extracted_call: dict[str, object] = {}
    captured: dict[str, object] = {}

    def _fake_extract(archive: Path, dest: Path, force_extract: bool) -> bool:
        extracted_call["archive"] = archive
        extracted_call["dest"] = dest
        extracted_call["force_extract"] = force_extract
        _write_prepped_payload(dest)
        return True

    def _fake_run(cmd, cwd=None, env=None, check=False, **kwargs):
        del check, kwargs
        captured["cmd"] = [str(x) for x in cmd]
        captured["cwd"] = cwd
        captured["env"] = dict(env or {})
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(
        run_relocated_mode, "_default_prepped_archive_path", lambda _repo: default_archive
    )
    monkeypatch.setattr(run_relocated_mode, "_extract_prepped_archive", _fake_extract)
    monkeypatch.setattr(run_relocated_mode.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        run_relocated_mode,
        "_build_main_command",
        lambda repo_root, raw: [sys.executable, str(repo_root / "main.py"), "--noop"],
    )

    argv = [
        "run_relocated_mode.py",
        "--prepped-root",
        str(prepped_root),
        "--main-args",
        "--run-id reloc_test --pdb TEST --fast --test-fda",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    rc = run_relocated_mode.main()
    assert rc == 0
    assert extracted_call["archive"] == default_archive.resolve()
    assert extracted_call["dest"] == prepped_root.resolve()
    assert extracted_call["force_extract"] is False


def test_dry_run_does_not_require_prepped_payload_or_extract(
    monkeypatch, tmp_path: Path
) -> None:
    prepped_root = tmp_path / "prepped"
    default_archive = tmp_path / "missing_default.tar.zst"
    extract_called = {"value": False}

    def _fake_extract(archive: Path, dest: Path, force_extract: bool) -> bool:
        del archive, dest, force_extract
        extract_called["value"] = True
        return False

    monkeypatch.setattr(
        run_relocated_mode, "_default_prepped_archive_path", lambda _repo: default_archive
    )
    monkeypatch.setattr(run_relocated_mode, "_extract_prepped_archive", _fake_extract)
    monkeypatch.setattr(run_relocated_mode, "_has_prepped_payload", lambda _: False)

    argv = [
        "run_relocated_mode.py",
        "--prepped-root",
        str(prepped_root),
        "--main-args",
        "--run-id reloc_test --pdb TEST --fast --test-fda",
        "--dry-run",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    rc = run_relocated_mode.main()
    assert rc == 0
    assert extract_called["value"] is False


def test_permission_error_during_mkdir_is_friendly(monkeypatch, tmp_path: Path) -> None:
    denied_root = tmp_path / "denied_root"

    def _deny_mkdir(self, parents=False, exist_ok=False):
        del self, parents, exist_ok
        raise PermissionError("nope")

    monkeypatch.setattr(run_relocated_mode.Path, "mkdir", _deny_mkdir)

    argv = [
        "run_relocated_mode.py",
        "--all-dirs",
        str(denied_root),
        "--main-args",
        "--run-id reloc_test --pdb TEST --fast --test-fda",
        "--dry-run",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(RuntimeError, match="Unable to create directory"):
        run_relocated_mode.main()


def test_extract_prepped_archive_real_tar_dangling_symlink_does_not_fail(
    tmp_path: Path,
) -> None:
    _require_tar_tools()
    dangling = "/nonexistent/intermediates/3KFA"
    _, archive_path = _make_archive_with_intermediates_link(tmp_path, dangling)

    prepped_root = tmp_path / "relocated_bundle" / "prepped_ligands"
    extracted = run_relocated_mode._extract_prepped_archive(
        archive_path, prepped_root, force_extract=False
    )
    assert extracted is True
    link_path = prepped_root / "3KFA" / "intermediates_src"
    assert link_path.is_symlink()


def test_extract_prepped_archive_real_tar_rewrites_absolute_symlink_target(
    tmp_path: Path,
) -> None:
    _require_tar_tools()
    old_root = tmp_path / "archive_bundle"
    old_target = old_root / "intermediates" / "3KFA"
    _, archive_path = _make_archive_with_intermediates_link(tmp_path, str(old_target))

    prepped_root = tmp_path / "relocated_bundle" / "prepped_ligands"
    run_relocated_mode._extract_prepped_archive(
        archive_path, prepped_root, force_extract=False
    )
    link_path = prepped_root / "3KFA" / "intermediates_src"
    expected = prepped_root.parent / "intermediates" / "3KFA"
    assert link_path.is_symlink()
    assert os.readlink(link_path) == str(expected)


def test_extract_prepped_archive_real_tar_keeps_relative_symlink_target(
    tmp_path: Path,
) -> None:
    _require_tar_tools()
    relative_target = "../shared/intermediates"
    _, archive_path = _make_archive_with_intermediates_link(tmp_path, relative_target)

    prepped_root = tmp_path / "relocated_bundle" / "prepped_ligands"
    run_relocated_mode._extract_prepped_archive(
        archive_path, prepped_root, force_extract=False
    )
    link_path = prepped_root / "3KFA" / "intermediates_src"
    assert link_path.is_symlink()
    assert os.readlink(link_path) == relative_target


def test_extract_prepped_archive_real_tar_keeps_external_absolute_symlink_target(
    tmp_path: Path,
) -> None:
    _require_tar_tools()
    external_target = "/opt/shared/intermediates/3KFA"
    _, archive_path = _make_archive_with_intermediates_link(tmp_path, external_target)

    prepped_root = tmp_path / "relocated_bundle" / "prepped_ligands"
    run_relocated_mode._extract_prepped_archive(
        archive_path, prepped_root, force_extract=False
    )
    link_path = prepped_root / "3KFA" / "intermediates_src"
    assert link_path.is_symlink()
    assert os.readlink(link_path) == external_target


def test_extract_prepped_archive_uses_lock_and_parent_tmpdir(
    tmp_path: Path, monkeypatch
) -> None:
    archive_path = tmp_path / "prepped_ligands.tar.zst"
    archive_path.write_text("placeholder", encoding="utf-8")
    prepped_root = tmp_path / "relocated_bundle" / "prepped_ligands"

    lock_calls: list[Path] = []

    @contextmanager
    def _fake_lock(lock_path: Path):
        lock_calls.append(lock_path)
        yield

    tmp_targets: list[Path] = []

    def _fake_run(cmd, capture_output=False, text=False):
        del capture_output, text
        cmd_list = [str(x) for x in cmd]
        extract_root = Path(cmd_list[cmd_list.index("-C") + 1])
        tmp_targets.append(extract_root)
        payload = extract_root / "prepped_ligands" / "demo_lib"
        payload.mkdir(parents=True, exist_ok=True)
        (payload / "demo_1.pdbqt").write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(run_relocated_mode, "_ensure_tool_available", lambda _: None)
    monkeypatch.setattr(run_relocated_mode, "_exclusive_extract_lock", _fake_lock)
    monkeypatch.setattr(run_relocated_mode.subprocess, "run", _fake_run)

    extracted = run_relocated_mode._extract_prepped_archive(
        archive_path, prepped_root, force_extract=False
    )
    assert extracted is True
    assert lock_calls == [prepped_root.parent / ".prepped_extract.lock"]
    assert tmp_targets
    assert all(path.parent == prepped_root.parent for path in tmp_targets)
    assert (prepped_root / "demo_lib" / "demo_1.pdbqt").exists()


def test_relocated_slurm_sim_keeps_docked_outputs_run_scoped(
    monkeypatch, tmp_path: Path
) -> None:
    from src.path_router import path_router as router

    all_root = tmp_path / "stor" / "home" / "mpg2352" / "relocated_sim"
    run_id = "sim_reloc_run"
    pdb_id = "BNJS"
    observed_by_task: dict[str, str] = {}

    monkeypatch.setattr(run_relocated_mode, "_has_prepped_payload", lambda _: True)

    def _fake_run(cmd, cwd=None, env=None, check=False, **kwargs):
        del cwd, check, kwargs
        cmd_list = [str(x) for x in cmd]
        parsed_run_id = run_id
        if "--run-id" in cmd_list:
            parsed_run_id = cmd_list[cmd_list.index("--run-id") + 1]
        task_id = str((env or {}).get("SLURM_ARRAY_TASK_ID", "0"))

        with pytest.MonkeyPatch.context() as local_mp:
            local_mp.setenv("DOCKED_DIR", str((env or {})["DOCKED_DIR"]))
            local_mp.delenv("ATLAS_RUN_ID", raising=False)
            local_mp.setattr(router, "_ROUTER_ROOTS", None)

            # Simulate a pre-main router touch before RUN_ID is exported.
            _ = router.docked_dir(pdb_id)

            local_mp.setenv("ATLAS_RUN_ID", parsed_run_id)
            resolved = router.docked_dir(pdb_id)
            resolved.mkdir(parents=True, exist_ok=True)

        observed_by_task[task_id] = str(resolved)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(run_relocated_mode.subprocess, "run", _fake_run)

    for task_id in ("0", "1"):
        monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
        monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "2")
        monkeypatch.setenv("SLURM_ARRAY_TASK_MIN", "0")
        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", task_id)
        monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "localsim")

        argv = [
            "run_relocated_mode.py",
            "--all-dirs",
            str(all_root),
            "--skip-extract",
            "--main-args",
            f"--run-id {run_id} --pdb {pdb_id} --fast --test-fda",
        ]
        monkeypatch.setattr(sys, "argv", argv)
        rc = run_relocated_mode.main()
        assert rc == 0

    expected_scoped = all_root / "outputs" / "docked" / run_id / pdb_id
    unexpected_flat = all_root / "outputs" / "docked" / pdb_id

    assert expected_scoped.exists()
    assert not unexpected_flat.exists()
    assert observed_by_task == {"0": str(expected_scoped), "1": str(expected_scoped)}


def test_relocated_slurm_sim_keeps_post_docked_outputs_run_scoped(
    monkeypatch, tmp_path: Path
) -> None:
    all_root = tmp_path / "stor" / "home" / "mpg2352" / "relocated_post_sim"
    run_id = "sim_post_run"
    pdb_id = "BNJS"
    observed_by_task: dict[str, str] = {}

    monkeypatch.setattr(run_relocated_mode, "_has_prepped_payload", lambda _: True)

    def _fake_run(cmd, cwd=None, env=None, check=False, **kwargs):
        del cwd, check, kwargs
        cmd_list = [str(x) for x in cmd]
        parsed_run_id = run_id
        if "--run-id" in cmd_list:
            parsed_run_id = cmd_list[cmd_list.index("--run-id") + 1]
        task_id = str((env or {}).get("SLURM_ARRAY_TASK_ID", "0"))

        post_root = Path(str((env or {})["POST_DOCKED_DIR"]))
        scoped = post_root / parsed_run_id / pdb_id / "HOLO" / "pH7_0"
        scoped.mkdir(parents=True, exist_ok=True)
        (scoped / "scorch_scores_all.csv").write_text(
            "scorch_mu_decoy,scorch_sigma_decoy,scorch_n_decoys\n0.5,0.2,10\n",
            encoding="utf-8",
        )

        observed_by_task[task_id] = str(scoped)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(run_relocated_mode.subprocess, "run", _fake_run)

    for task_id in ("0", "1"):
        monkeypatch.setenv("ATLAS_DISTRIBUTED_MODE", "slurm_array")
        monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "2")
        monkeypatch.setenv("SLURM_ARRAY_TASK_MIN", "0")
        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", task_id)
        monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "localsim")

        argv = [
            "run_relocated_mode.py",
            "--all-dirs",
            str(all_root),
            "--skip-extract",
            "--main-args",
            f"--run-id {run_id} --pdb {pdb_id} --fast --test-fda",
        ]
        monkeypatch.setattr(sys, "argv", argv)
        rc = run_relocated_mode.main()
        assert rc == 0

    expected_scoped = all_root / "outputs" / "post_docked" / run_id / pdb_id
    unexpected_flat = all_root / "outputs" / "post_docked" / pdb_id

    assert expected_scoped.exists()
    assert not unexpected_flat.exists()
    assert all(
        path.startswith(str(expected_scoped)) for path in observed_by_task.values()
    )
