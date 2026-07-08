from pathlib import Path

import pytest

from config.runtime_config import load_config
from config.runtime_config import validate_config
from config.normalize import normalize_config


def test_inline_comments_are_stripped_and_typed(tmp_path: Path, monkeypatch) -> None:
    for key in (
        "CPU",
        "USE_GNINA",
        "EARLY_RECENTER_RATIO",
    ):
        monkeypatch.delenv(key, raising=False)

    cfg_path = tmp_path / "config.txt"
    cfg_path.write_text(
        "\n".join(
            [
                "CPU=1 # core count",
                "USE_GNINA=true # engine flag",
                "EARLY_RECENTER_RATIO=0.70 # ratio",
                'SOME_PATH="C:/path/with#hash" # comment',
                "",
            ]
        )
    )

    cfg = load_config(config_path=str(cfg_path), base_dir=tmp_path)

    assert cfg["CPU"] == 1
    assert cfg["USE_GNINA"] is True
    assert isinstance(cfg["EARLY_RECENTER_RATIO"], float)
    assert cfg["EARLY_RECENTER_RATIO"] == 0.70
    assert "with#hash" in cfg["SOME_PATH"]
    assert cfg["SOME_PATH"].strip('"') == "C:/path/with#hash"


def test_validate_config_missing_input_dir_error_has_doctor_hint(tmp_path: Path) -> None:
    overall_dir = tmp_path / "overall"
    overall_dir.mkdir()
    missing_input_dir = tmp_path / "input_pdbs_missing"

    cfg = {
        "OUTPUT_DIR": str(tmp_path / "processed"),
        "INPUT_DIR": str(missing_input_dir),
        "PDBQT_DIR": str(tmp_path / "pdbqt"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "OVERALL_DIR": str(overall_dir),
        "MGLTOOLS_PYTHON": "/bin/echo",
        "PREPARE_RECEPTOR_SCRIPT": "/bin/echo",
        "VINA_EXE": "/bin/echo",
        "CPU": 1,
        "PYMOL_EXE": "/bin/echo",
        "PREPPED_LIGANDS_DIR": str(tmp_path / "prepped_ligands"),
        "EXTRACTED_LIGANDS_DIR": str(tmp_path / "extracted_ligands"),
        "LIGANDS_MOL2_DIR": str(tmp_path / "ligands_mol2"),
        "P2RANK_OUTPUT_DIR": str(tmp_path / "p2rank"),
        "CONFIGS_DIR": str(tmp_path / "configs"),
    }

    with pytest.raises(FileNotFoundError) as excinfo:
        validate_config(cfg)

    msg = str(excinfo.value)
    assert "INPUT_DIR" in msg
    assert "input_pdbs" in msg
    assert "atlas --doctor" in msg


def test_normalize_config_derives_adt_scripts_from_pythonsh(tmp_path: Path) -> None:
    mgl_root = tmp_path / "mgltools"
    bin_dir = mgl_root / "bin"
    utilities_dir = (
        mgl_root / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24"
    )
    bin_dir.mkdir(parents=True)
    utilities_dir.mkdir(parents=True)

    pythonsh = bin_dir / "pythonsh"
    prepare_receptor = utilities_dir / "prepare_receptor4.py"
    prepare_ligand = utilities_dir / "prepare_ligand4.py"
    pythonsh.write_text("#!/bin/sh\n")
    prepare_receptor.write_text("# receptor\n")
    prepare_ligand.write_text("# ligand\n")

    cfg = normalize_config(
        {
            "OVERALL_DIR": str(tmp_path),
            "MGLTOOLS_PYTHON": str(pythonsh),
            "PREPARE_RECEPTOR_SCRIPT": "",
            "PREPARE_LIGAND_SCRIPT": "",
        }
    )

    assert cfg["MGLTOOLS_PATH"] == str(mgl_root)
    assert cfg["PREPARE_RECEPTOR_SCRIPT"] == str(prepare_receptor)
    assert cfg["PREPARE_LIGAND_SCRIPT"] == str(prepare_ligand)
