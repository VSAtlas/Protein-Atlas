from __future__ import annotations

from pathlib import Path

import pytest

import main
from cli import run_profiles


def _touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("HEADER TEST\n", encoding="utf-8")


def test_dude_overrides_test_mode_enable_with_full_map_coverage(
    tmp_path: Path,
) -> None:
    input_dir = tmp_path / "input_pdbs"
    _touch(input_dir / "1ABC.pdb")
    _touch(input_dir / "2DEF.pdb")
    _touch(input_dir / "3GHI.pdb")

    cfg = {
        "INPUT_DIR": str(input_dir),
        "TEST_LIBRARY_MAP": {"1ABC": "lib1", "2DEF": "lib2", "3GHI": "lib3"},
        "SPECIFIED_PROTEINS": "",
        "TEST_MODE_ENABLE": "off",
    }

    run_profiles._apply_dude_overrides(cfg)
    assert cfg["TEST_MODE_ENABLE"] == "dud"

    selected = main.select_pdb_files_for_run(cfg, ["main.py"])
    assert sorted(selected) == ["1ABC.pdb", "2DEF.pdb", "3GHI.pdb"]


def test_dude_filter_requires_mappings_for_explicitly_specified_proteins(
    tmp_path: Path,
) -> None:
    input_dir = tmp_path / "input_pdbs"
    _touch(input_dir / "1ABC.pdb")
    _touch(input_dir / "2DEF.pdb")
    _touch(input_dir / "3GHI.pdb")

    cfg = {
        "INPUT_DIR": str(input_dir),
        "TEST_LIBRARY_MAP": {"1ABC": "lib1", "2DEF": "lib2"},
        "SPECIFIED_PROTEINS": "",
        "TEST_MODE_ENABLE": "off",
    }

    run_profiles._apply_dude_overrides(cfg)
    with pytest.raises(SystemExit) as excinfo:
        main.select_pdb_files_for_run(cfg, ["main.py", "--pdb", "3GHI"])
    assert excinfo.value.code == 2


def test_dude_flag_is_not_parsed_as_specified_protein(tmp_path: Path) -> None:
    input_dir = tmp_path / "input_pdbs"
    _touch(input_dir / "1ABC.pdb")
    _touch(input_dir / "2DEF.pdb")
    _touch(input_dir / "3GHI.pdb")

    cfg = {
        "INPUT_DIR": str(input_dir),
        "TEST_LIBRARY_MAP": {"1ABC": "lib1", "2DEF": "lib2", "3GHI": "lib3"},
        "SPECIFIED_PROTEINS": "",
        "TEST_MODE_ENABLE": "off",
    }

    run_profiles._apply_dude_overrides(cfg)
    selected = main.select_pdb_files_for_run(cfg, ["main.py", "-dude", "-fast"])
    assert sorted(selected) == ["1ABC.pdb", "2DEF.pdb", "3GHI.pdb"]
    assert cfg["_EFFECTIVE_SPECIFIED_PROTEINS"] == []


def test_dude_requires_non_empty_test_library_map(tmp_path: Path) -> None:
    input_dir = tmp_path / "input_pdbs"
    _touch(input_dir / "1ABC.pdb")
    _touch(input_dir / "2DEF.pdb")

    cfg = {
        "INPUT_DIR": str(input_dir),
        "TEST_LIBRARY_MAP": {},
        "SPECIFIED_PROTEINS": "",
        "TEST_MODE_ENABLE": "off",
    }

    run_profiles._apply_dude_overrides(cfg)
    with pytest.raises(SystemExit) as excinfo:
        main.select_pdb_files_for_run(cfg, ["main.py"])
    assert excinfo.value.code == 2


def test_mixed_mode_dud_plus_fda_allows_unmapped_proteins(tmp_path: Path) -> None:
    input_dir = tmp_path / "input_pdbs"
    _touch(input_dir / "1ABC.pdb")
    _touch(input_dir / "2DEF.pdb")
    _touch(input_dir / "3GHI.pdb")

    cfg = {
        "INPUT_DIR": str(input_dir),
        "TEST_LIBRARY_MAP": {"1ABC": "lib1"},
        "SPECIFIED_PROTEINS": "",
        "TEST_MODE_ENABLE": "fda+dud",
    }

    selected = main.select_pdb_files_for_run(cfg, ["main.py"])
    assert sorted(selected) == ["1ABC.pdb", "2DEF.pdb", "3GHI.pdb"]


def test_mixed_mode_dud_plus_fda_with_empty_map_does_not_error(tmp_path: Path) -> None:
    input_dir = tmp_path / "input_pdbs"
    _touch(input_dir / "1ABC.pdb")
    _touch(input_dir / "2DEF.pdb")

    cfg = {
        "INPUT_DIR": str(input_dir),
        "TEST_LIBRARY_MAP": {},
        "SPECIFIED_PROTEINS": "",
        "TEST_MODE_ENABLE": "fda+dud",
    }

    selected = main.select_pdb_files_for_run(cfg, ["main.py"])
    assert sorted(selected) == ["1ABC.pdb", "2DEF.pdb"]


def test_pdb_selection_splits_comma_values_without_treating_flags_as_targets(
    tmp_path: Path,
) -> None:
    input_dir = tmp_path / "input_pdbs"
    _touch(input_dir / "1BN1.pdb")
    _touch(input_dir / "2OJ9.pdb")
    _touch(input_dir / "PDBS.pdb")

    base_cfg = {"INPUT_DIR": str(input_dir), "SPECIFIED_PROTEINS": ""}

    cfg = dict(base_cfg)
    selected = main.select_pdb_files_for_run(
        cfg, ["main.py", "--pdbs", "1BN1,2OJ9", "--fast"]
    )
    assert sorted(selected) == ["1BN1.pdb", "2OJ9.pdb"]
    assert cfg["_EFFECTIVE_SPECIFIED_PROTEINS"] == ["1BN1", "2OJ9"]

    cfg = dict(base_cfg)
    selected = main.select_pdb_files_for_run(
        cfg, ["main.py", "--pdb", "1BN1,2OJ9", "--fast"]
    )
    assert sorted(selected) == ["1BN1.pdb", "2OJ9.pdb"]
    assert cfg["_EFFECTIVE_SPECIFIED_PROTEINS"] == ["1BN1", "2OJ9"]
