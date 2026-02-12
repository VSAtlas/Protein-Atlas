from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import run_relocated_mode  # noqa: E402


def _write_prepped_payload(prepped_root: Path) -> None:
    lib_dir = prepped_root / "demo_lib"
    lib_dir.mkdir(parents=True, exist_ok=True)
    (lib_dir / "demo_1.pdbqt").write_text("MODEL 1\nENDMDL\n", encoding="utf-8")


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
    assert env["DOCKED_DIR"] == str((all_root / "docked").resolve())
    assert env["POST_DOCKED_DIR"] == str((all_root / "post_docked").resolve())
    assert env["OUTPUT_DIR"] == str((all_root / "processed_pdbs").resolve())
    assert env["CONFIGS_DIR"] == str((all_root / "configs").resolve())

    assert all_root.exists()
    assert (all_root / "prepped_ligands").exists()
    assert (all_root / "extracted_ligands").exists()
    assert (all_root / "docked").exists()
    assert (all_root / "post_docked").exists()
    assert (all_root / "processed_pdbs").exists()
    assert (all_root / "configs").exists()


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
    assert env["OUTPUT_DIR"] == str((all_root / "processed_pdbs").resolve())
    assert env["CONFIGS_DIR"] == str((all_root / "configs").resolve())


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
