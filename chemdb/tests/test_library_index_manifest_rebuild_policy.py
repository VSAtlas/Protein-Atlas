from __future__ import annotations

import json
import logging
import time
from pathlib import Path

from prep_ligands.library_index import LibraryIndex


def _write_pdbqt(path: Path, text: str = "MODEL\nATOM\nENDMDL\n") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_manifest(root: Path, payload: dict) -> Path:
    path = root / "_manifest.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path


def test_manifest_reused_when_mtime_drift_but_references_valid(tmp_path: Path) -> None:
    root = tmp_path / "lib"
    lig = root / "actives" / "lig_a.pdbqt"
    _write_pdbqt(lig)
    manifest = _write_manifest(
        root,
        {
            "mtime": 0.0,
            "child_dir_mtimes": {},
            "entries": {"lig_a": "actives/lig_a.pdbqt"},
            "filenames": {"lig_a.pdbqt": "actives/lig_a.pdbqt"},
        },
    )
    pre = manifest.stat().st_mtime_ns
    time.sleep(0.01)
    # Change root mtime without breaking manifest references.
    (root / "note.txt").write_text("touch\n", encoding="utf-8")

    idx = LibraryIndex(logger=logging.getLogger("test.libidx.reuse"))
    idx.load([root])
    post = manifest.stat().st_mtime_ns
    assert pre == post
    assert idx.lookup("lig_a.pdbqt", [root]) == lig


def test_manifest_rebuilds_when_referenced_path_missing(tmp_path: Path) -> None:
    root = tmp_path / "lib"
    _write_pdbqt(root / "actives" / "lig_a.pdbqt")
    manifest = _write_manifest(
        root,
        {
            "mtime": 0.0,
            "entries": {"missing": "actives/lig_missing.pdbqt"},
            "filenames": {"lig_missing.pdbqt": "actives/lig_missing.pdbqt"},
        },
    )
    pre = manifest.stat().st_mtime_ns
    time.sleep(0.01)

    idx = LibraryIndex(logger=logging.getLogger("test.libidx.missing"))
    idx.load([root])
    post = manifest.stat().st_mtime_ns
    assert post > pre
    assert idx.lookup("lig_a.pdbqt", [root]) == (root / "actives" / "lig_a.pdbqt")


def test_manifest_rebuilds_when_referenced_path_empty(tmp_path: Path) -> None:
    root = tmp_path / "lib"
    _write_pdbqt(root / "actives" / "lig_a.pdbqt")
    _write_pdbqt(root / "actives" / "lig_empty.pdbqt", text="")
    manifest = _write_manifest(
        root,
        {
            "mtime": 0.0,
            "entries": {"lig_empty": "actives/lig_empty.pdbqt"},
            "filenames": {"lig_empty.pdbqt": "actives/lig_empty.pdbqt"},
        },
    )
    pre = manifest.stat().st_mtime_ns
    time.sleep(0.01)

    idx = LibraryIndex(logger=logging.getLogger("test.libidx.empty"))
    idx.load([root])
    post = manifest.stat().st_mtime_ns
    assert post > pre
    assert idx.lookup("lig_a.pdbqt", [root]) == (root / "actives" / "lig_a.pdbqt")


def test_manifest_rebuilds_when_path_traversal_reference_present(tmp_path: Path) -> None:
    root = tmp_path / "lib"
    _write_pdbqt(root / "actives" / "lig_a.pdbqt")
    manifest = _write_manifest(
        root,
        {
            "mtime": 0.0,
            "entries": {"evil": "../outside.pdbqt"},
            "filenames": {"outside.pdbqt": "../outside.pdbqt"},
        },
    )
    pre = manifest.stat().st_mtime_ns
    time.sleep(0.01)

    idx = LibraryIndex(logger=logging.getLogger("test.libidx.traversal"))
    idx.load([root])
    post = manifest.stat().st_mtime_ns
    assert post > pre
    assert idx.lookup("lig_a.pdbqt", [root]) == (root / "actives" / "lig_a.pdbqt")
