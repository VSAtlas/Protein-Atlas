from __future__ import annotations

import json
import shutil
from pathlib import Path

from post_docking import artifact_retention


def _touch(path: Path, content: str = "X") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _patch_fast_retention(monkeypatch) -> None:
    monkeypatch.setattr(artifact_retention, "_tooling_preflight", lambda: (True, "ok"))
    monkeypatch.setattr(
        artifact_retention,
        "_compress_tar_to_zst",
        lambda tar_path, archive_path, threads: shutil.move(str(tar_path), str(archive_path)),
    )
    monkeypatch.setattr(artifact_retention, "_verify_archive", lambda archive_path, entries: (True, []))


def _load_index(run_docked: Path) -> dict:
    index_path = run_docked / "_artifact_archives" / "archive_index.json"
    assert index_path.exists()
    return json.loads(index_path.read_text(encoding="utf-8"))


def test_retention_archives_and_deletes_exact_counts(monkeypatch, tmp_path: Path) -> None:
    _patch_fast_retention(monkeypatch)
    run_id = "ret_count_run"
    run_docked = tmp_path / "docked" / run_id
    run_post = tmp_path / "post_docked" / run_id
    _touch(run_docked / "BOJG" / "HOLO" / "pH7_7" / "stage1" / "lig_a.pdbqt")
    _touch(run_docked / "BOJG" / "HOLO" / "pH7_7" / "stage1" / "lig_b.pdbqt")
    _touch(run_post / "BOJG" / "HOLO" / "pH7_7" / "dock6_pdbqt" / "lig_c.pdbqt")
    keep_csv = run_docked / "BOJG" / "HOLO" / "pH7_7" / "docking_score_summary.csv"
    _touch(keep_csv, "csv\n")

    rc = artifact_retention.main(
        [
            "--run-id",
            run_id,
            "--mode",
            "rerun_safe",
            "--docked-root",
            str(tmp_path / "docked"),
            "--post-docked-root",
            str(tmp_path / "post_docked"),
            "--threads",
            "4",
            "--group-workers",
            "2",
        ]
    )
    assert rc == 0
    payload = _load_index(run_docked)
    groups = payload.get("groups") or []
    assert sum(int(g.get("entry_count", 0)) for g in groups) == 3
    assert not (run_docked / "BOJG" / "HOLO" / "pH7_7" / "stage1" / "lig_a.pdbqt").exists()
    assert not (run_docked / "BOJG" / "HOLO" / "pH7_7" / "stage1" / "lig_b.pdbqt").exists()
    assert not (run_post / "BOJG" / "HOLO" / "pH7_7" / "dock6_pdbqt" / "lig_c.pdbqt").exists()
    assert keep_csv.exists()


def test_retention_minimal_disk_deletes_extra_suffixes(monkeypatch, tmp_path: Path) -> None:
    _patch_fast_retention(monkeypatch)
    run_id = "ret_min_disk_run"
    run_docked = tmp_path / "docked" / run_id
    _touch(run_docked / "BNJS" / "HOLO" / "pH7_0" / "stage1" / "lig_a.pdbqt")
    extra_tmp = run_docked / "BNJS" / "HOLO" / "pH7_0" / "stage1" / "left.tmp"
    _touch(extra_tmp)

    rc = artifact_retention.main(
        [
            "--run-id",
            run_id,
            "--mode",
            "minimal_disk",
            "--docked-root",
            str(tmp_path / "docked"),
            "--post-docked-root",
            str(tmp_path / "post_docked"),
        ]
    )
    assert rc == 0
    payload = _load_index(run_docked)
    groups = payload.get("groups") or []
    assert sum(int(g.get("entry_count", 0)) for g in groups) == 1
    assert not extra_tmp.exists()


def test_retention_verify_only_does_not_mutate(monkeypatch, tmp_path: Path) -> None:
    _patch_fast_retention(monkeypatch)
    run_id = "ret_verify_run"
    run_docked = tmp_path / "docked" / run_id
    target = run_docked / "BNNQ" / "HOLO" / "pH5_6" / "stage1" / "lig_keep.pdbqt"
    _touch(target)

    rc = artifact_retention.main(
        [
            "--run-id",
            run_id,
            "--mode",
            "rerun_safe",
            "--verify-only",
            "--docked-root",
            str(tmp_path / "docked"),
            "--post-docked-root",
            str(tmp_path / "post_docked"),
        ]
    )
    assert rc == 0
    assert target.exists()


def test_retention_combo_scope_filters_exact_groups(monkeypatch, tmp_path: Path) -> None:
    _patch_fast_retention(monkeypatch)
    run_id = "ret_combo_run"
    run_docked = tmp_path / "docked" / run_id
    keep = run_docked / "BOJG" / "HOLO" / "pH7_7" / "stage1" / "lig_keep.pdbqt"
    skip = run_docked / "BOJG" / "HOLO" / "pH7_2" / "stage1" / "lig_skip.pdbqt"
    _touch(keep)
    _touch(skip)

    rc = artifact_retention.main(
        [
            "--run-id",
            run_id,
            "--mode",
            "rerun_safe",
            "--combo",
            "BOJG:HOLO:pH7_7",
            "--docked-root",
            str(tmp_path / "docked"),
            "--post-docked-root",
            str(tmp_path / "post_docked"),
        ]
    )
    assert rc == 0
    payload = _load_index(run_docked)
    groups = payload.get("groups") or []
    assert len(groups) == 1
    assert sum(int(g.get("entry_count", 0)) for g in groups) == 1
    assert not keep.exists()
    assert skip.exists()
