import logging
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from path_router import make_paths
from docking_ligands import prepare_and_filter_ligands


def _make_fake_lig(folder: Path, stem: str) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    lig_path = folder / f"{stem}.pdbqt"
    lines = ["REMARK fake ligand for DUD/FDA mode test"]
    for i in range(1, 6):
        lines.append(
            f"ATOM  {i:5d}  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C"
        )
    lines.append("END")
    lig_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return lig_path


def test_prepare_and_filter_ligands_dud_fda_modes(tmp_path):
    root = tmp_path / "fda_dud_mode_test"
    root.mkdir()

    cfg = {
        "OVERALL_DIR": str(root),
        "INPUT_DIR": str(root / "input_pdbs"),
        "OUTPUT_DIR": str(root / "processed_pdbs"),
        "DOCKED_DIR": str(root / "docked"),
        "OUTPUT_LIGANDS_DIR": str(root / "prepped_ligands"),
        "LIBRARY_SUBDIR_DEFAULT": "fda_test",
        "TEST_LIBRARY_MAP": {"TESTPDB": "dud_test"},
        "TEST_MODE_ENABLE": "off",
        "USE_GNINA": "false",
        "RUN_ID": "test_run",
    }

    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR", "OUTPUT_LIGANDS_DIR"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)

    lig_root = Path(cfg["OUTPUT_LIGANDS_DIR"])
    dud_root = lig_root / "dud_test"
    fda_root = lig_root / "fda_test"

    _make_fake_lig(dud_root, "dud_test_00001")
    _make_fake_lig(fda_root, "fda_test_00001")

    paths = make_paths(cfg, base_id="TESTPDB", pdb_file="TESTPDB.pdb")

    logger = logging.getLogger("fda_dud_mode_test")
    logger.setLevel(logging.INFO)

    def libs_for_mode(mode: str) -> Counter:
        cfg["TEST_MODE_ENABLE"] = mode
        ligands, _, _ = prepare_and_filter_ligands(cfg, paths, logger)
        return Counter(Path(path).parent.name for path in ligands)

    libs_off = libs_for_mode("off")
    libs_dud = libs_for_mode("dud")
    libs_both = libs_for_mode("fda+dud")

    assert libs_off == Counter({"fda_test": 1})
    assert libs_dud == Counter({"dud_test": 1})
    assert libs_both == Counter({"dud_test": 1, "fda_test": 1})


def test_prepare_and_filter_ligands_filters_ph_dirs_when_off(tmp_path):
    root = tmp_path / "ph_filter"
    root.mkdir()

    cfg = {
        "OVERALL_DIR": str(root),
        "INPUT_DIR": str(root / "input_pdbs"),
        "OUTPUT_DIR": str(root / "processed_pdbs"),
        "DOCKED_DIR": str(root / "docked"),
        "OUTPUT_LIGANDS_DIR": str(root / "prepped_ligands"),
        "LIBRARY_SUBDIR_DEFAULT": "ph_lib",
        "TEST_MODE_ENABLE": "off",
        "USE_GNINA": "false",
        "RUN_ID": "test_run",
    }

    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR", "OUTPUT_LIGANDS_DIR"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)

    lib_root = Path(cfg["OUTPUT_LIGANDS_DIR"]) / cfg["LIBRARY_SUBDIR_DEFAULT"]
    base_lig = _make_fake_lig(lib_root, "base_a")
    ph_lig = _make_fake_lig(lib_root / "7.4", "ph_dir_a")
    micro_lig = _make_fake_lig(lib_root / "microstates", "micro_a")

    paths = make_paths(cfg, base_id="TESTPDB", pdb_file="TESTPDB.pdb")

    logger = logging.getLogger("ph_filter_test")
    logger.setLevel(logging.INFO)

    cfg["PH_LIGAND_MODE"] = "off"
    ligs_off, _, _ = prepare_and_filter_ligands(cfg, paths, logger)
    assert {Path(p) for p in ligs_off} == {base_lig}

    cfg["PH_LIGAND_MODE"] = "context_window"
    ligs_on, _, _ = prepare_and_filter_ligands(cfg, paths, logger)
    assert {Path(p) for p in ligs_on} == {base_lig, ph_lig, micro_lig}
