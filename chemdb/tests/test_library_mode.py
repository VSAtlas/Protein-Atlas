import logging
from pathlib import Path

from docking.library_mode import compute_allowed_library_roots, parse_test_libraries
from docking.docking_subruns import subruns_for_tokens
from docking.docking_ligands import prepare_and_filter_ligands
from path_router import make_paths


def _make_fake_lig(folder: Path, stem: str) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    lig_path = folder / f"{stem}.pdbqt"
    lines = ["REMARK fake ligand for test mode"]
    for i in range(1, 6):
        lines.append(
            f"ATOM  {i:5d}  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C"
        )
    lines.append("END")
    lig_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return lig_path


def test_parse_test_libraries_custom_token_with_underscore(monkeypatch):
    monkeypatch.setenv("TEST_MODE_ENABLE", "fda_dud+fda")
    assert parse_test_libraries({}) == ["fda_dud", "fda"]

    monkeypatch.setenv("TEST_MODE_ENABLE", "dud+fda")
    assert parse_test_libraries({}) == ["dud", "fda"]


def test_compute_allowed_library_roots_custom_and_fda(tmp_path):
    base = tmp_path / "ligands"
    (base / "fda_library").mkdir(parents=True, exist_ok=True)
    (base / "test_library_10").mkdir(parents=True, exist_ok=True)

    cfg = {
        "OUTPUT_LIGANDS_DIR": str(base),
        "LIBRARY_SUBDIR_DEFAULT": "fda_library",
        "TEST_MODE_ENABLE": "test_library_10+fda",
    }
    logger = logging.getLogger("test.library_roots")

    roots = compute_allowed_library_roots(cfg, "1ABC", logger)
    assert set(roots) == {base / "fda_library", base / "test_library_10"}


def test_dud_uses_test_library_map(tmp_path):
    base = tmp_path / "ligands"
    (base / "bench_lib").mkdir(parents=True, exist_ok=True)

    cfg = {
        "OUTPUT_LIGANDS_DIR": str(base),
        "TEST_LIBRARY_MAP": {"1ABC": "bench_lib"},
        "TEST_MODE_ENABLE": "dud",
    }
    logger = logging.getLogger("test.dud_map")

    roots = compute_allowed_library_roots(cfg, "1ABC", logger)
    assert base / "bench_lib" in roots


def test_use_deepcoy_false_skips_autogen_in_prepare_and_filter(tmp_path, monkeypatch):
    root = tmp_path / "deepcoy_off"
    cfg = {
        "OVERALL_DIR": str(root),
        "INPUT_DIR": str(root / "input_pdbs"),
        "OUTPUT_DIR": str(root / "processed_pdbs"),
        "DOCKED_DIR": str(root / "docked"),
        "OUTPUT_LIGANDS_DIR": str(root / "prepped_ligands"),
        "LIBRARY_SUBDIR_DEFAULT": "fda_test",
        "TEST_MODE_ENABLE": "dud",
        "TEST_LIBRARY_MAP": {},
        "USE_DEEPCOY": "off",
        "DEEPCOY_ENABLE_AUTOGEN_SDF": "on",
        "DEEPCOY_ENABLE_AUTOGEN_PDBQT": "on",
        "RUN_ID": "test_run",
        "USE_GNINA": "false",
    }

    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR", "OUTPUT_LIGANDS_DIR"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)

    fda_root = Path(cfg["OUTPUT_LIGANDS_DIR"]) / cfg["LIBRARY_SUBDIR_DEFAULT"]
    _make_fake_lig(fda_root, "fda_test_00001")

    paths = make_paths(cfg, base_id="TESTPDB", pdb_file="TESTPDB.pdb")
    logger = logging.getLogger("deepcoy_off_test")
    logger.setLevel(logging.INFO)

    def _boom(*_args, **_kwargs):
        raise AssertionError("DeepCoy autogen should not run when USE_DEEPCOY=off")

    monkeypatch.setattr(
        "docking.docking_ligands.ensure_deepcoy_decoy_pdbqts",
        _boom,
    )
    monkeypatch.delenv("TEST_MODE_ENABLE", raising=False)

    prepare_and_filter_ligands(cfg, paths, logger)


def test_subruns_include_custom_token(monkeypatch):
    monkeypatch.setenv("TEST_MODE_ENABLE", "test_library_10+fda")
    tokens = parse_test_libraries({})
    subruns = subruns_for_tokens(tokens)

    assert subruns[0].run_mode == "test_library_10"
    assert subruns[0].csv_prefix == "test_library_10_"
    assert subruns[0].stage_name_prefix == "test_library_10_"
    assert subruns[1].run_mode == "fda"
    assert subruns[1].csv_prefix == ""
    assert subruns[1].stage_name_prefix == ""
