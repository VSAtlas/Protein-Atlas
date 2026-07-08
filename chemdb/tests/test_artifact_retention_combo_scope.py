from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import cli.postrun_hooks_runtime as postrun_hooks
import pytest
from post_docking import artifact_retention


def _touch(path: Path, content: str = "X") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def test_postrun_artifact_retention_combo_scope_builds_expected_args(monkeypatch) -> None:
    captured: list[list[str]] = []

    def _fake_run(cmd, cwd=None, check=False, env=None):
        del cwd, check, env
        cmd_list = [str(x) for x in cmd]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)

    cfg = {
        "ARTIFACT_RETENTION": "rerun_safe",
        "CPU": 4,
        "DOCKED_DIR": "/tmp/docked_root",
        "POST_DOCKED_DIR": "/tmp/post_docked_root",
    }
    postrun_hooks._maybe_run_artifact_retention_for_combo(
        cfg,
        "run_combo_scope",
        "BOJG",
        "HOLO",
        "pH7_7",
    )

    assert captured
    cmd = captured[0]
    assert cmd[cmd.index("--pdb-id") + 1] == "BOJG"
    assert cmd[cmd.index("--variant") + 1] == "HOLO"
    assert cmd[cmd.index("--ph") + 1] == "pH7_7"


def test_postrun_artifact_retention_multi_combo_scope_builds_expected_args(
    monkeypatch,
) -> None:
    captured: list[list[str]] = []

    def _fake_run(cmd, cwd=None, check=False, env=None):
        del cwd, check, env
        cmd_list = [str(x) for x in cmd]
        captured.append(cmd_list)
        return subprocess.CompletedProcess(cmd_list, 0)

    monkeypatch.setattr(postrun_hooks.subprocess, "run", _fake_run)

    cfg = {
        "ARTIFACT_RETENTION": "rerun_safe",
        "CPU": 8,
        "ARTIFACT_RETENTION_GROUP_WORKERS": 2,
    }
    postrun_hooks._maybe_run_artifact_retention_for_combos(
        cfg,
        "run_combo_scope",
        [
            ("BOJG", "HOLO", "pH7_7"),
            ("BNJS", "HOLO", "pH7_0"),
            ("BNJS", "HOLO", "pH7_0"),
        ],
    )

    assert captured
    cmd = captured[0]
    assert cmd[cmd.index("--group-workers") + 1] == "2"
    combo_values = [
        cmd[idx + 1] for idx, token in enumerate(cmd) if token == "--combo"
    ]
    assert combo_values == ["BNJS:HOLO:pH7_0", "BOJG:HOLO:pH7_7"]


def test_artifact_retention_cli_rejects_variant_or_ph_without_pdb() -> None:
    with pytest.raises(SystemExit) as exc_info:
        artifact_retention.main(
            [
                "--run-id",
                "run_combo_scope",
                "--mode",
                "rerun_safe",
                "--variant",
                "HOLO",
            ]
        )
    assert "--variant/--ph filters require --pdb-id" in str(exc_info.value)

    with pytest.raises(SystemExit) as exc_info:
        artifact_retention.main(
            [
                "--run-id",
                "run_combo_scope",
                "--mode",
                "rerun_safe",
                "--ph",
                "pH7_7",
            ]
        )
    assert "--variant/--ph filters require --pdb-id" in str(exc_info.value)


def test_artifact_retention_cli_rejects_combo_with_legacy_scope() -> None:
    with pytest.raises(SystemExit) as exc_info:
        artifact_retention.main(
            [
                "--run-id",
                "run_combo_scope",
                "--mode",
                "rerun_safe",
                "--combo",
                "BOJG:HOLO:pH7_7",
                "--pdb-id",
                "BOJG",
            ]
        )
    assert "mutually exclusive" in str(exc_info.value)


def test_artifact_retention_discover_candidates_filters_combo(tmp_path: Path) -> None:
    run_id = "combo_filter_run"
    run_docked = tmp_path / "docked" / run_id
    run_post = tmp_path / "post_docked" / run_id
    _touch(run_docked / "BOJG" / "HOLO" / "pH7_7" / "stage1" / "lig_a.pdbqt")
    _touch(run_docked / "BOJG" / "HOLO" / "pH7_2" / "stage1" / "lig_b.pdbqt")
    _touch(run_post / "BOJG" / "HOLO" / "pH7_7" / "dock6_pdbqt" / "lig_c.pdbqt")
    _touch(run_docked / "BNJS" / "HOLO" / "pH7_0" / "stage1" / "lig_d.pdbqt")

    groups = artifact_retention._discover_candidates(
        run_id=run_id,
        run_docked=run_docked,
        run_post_docked=run_post,
        pdb_id_filter="BOJG",
        variant_filter="HOLO",
        ph_filter="ph7_7",
        include_globs=[],
        exclude_globs=[],
        logger=logging.getLogger("test_artifact_retention_discover"),
    )

    assert len(groups) == 2
    keys = {(key.source_root, key.pdb_id, key.variant, key.ph) for key in groups.keys()}
    assert ("docked", "BOJG", "HOLO", "pH7_7") in keys
    assert ("post_docked", "BOJG", "HOLO", "pH7_7") in keys
