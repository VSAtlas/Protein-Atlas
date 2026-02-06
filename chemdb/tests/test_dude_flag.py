from __future__ import annotations

from pathlib import Path

import main


def _touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("HEADER TEST\n", encoding="utf-8")


def test_dude_overrides_test_mode_enable_and_filters_to_test_library_map(
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

    main._apply_dude_overrides(cfg)
    assert cfg["TEST_MODE_ENABLE"] == "dud"

    selected = main.select_pdb_files_for_run(cfg, ["main.py"])
    assert sorted(selected) == ["1ABC.pdb", "2DEF.pdb"]


def test_dude_filter_keeps_explicitly_specified_proteins(tmp_path: Path) -> None:
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

    main._apply_dude_overrides(cfg)
    selected = main.select_pdb_files_for_run(cfg, ["main.py", "--pdb", "3GHI"])
    assert selected == ["3GHI.pdb"]


def test_dude_flag_is_not_parsed_as_specified_protein(tmp_path: Path) -> None:
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

    main._apply_dude_overrides(cfg)
    selected = main.select_pdb_files_for_run(cfg, ["main.py", "-dude", "-fast"])
    assert sorted(selected) == ["1ABC.pdb", "2DEF.pdb"]
