import logging
from pathlib import Path

from config.runtime_config import load_inputs
from docking.library_mode import compute_allowed_library_roots, parse_test_libraries
from docking.docking_subruns import _run_scoped_receptor_fallback, subruns_for_tokens
from cli.run_bootstrap import _retarget_default_test_library_map_for_test_fda


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


def test_test_mode_enable_env_override_feeds_runtime_config(monkeypatch):
    monkeypatch.setenv("TEST_MODE_ENABLE", "dud")
    cfg = load_inputs()

    assert cfg["TEST_MODE_ENABLE"] == "dud"
    assert parse_test_libraries(cfg) == ["dud"]


def test_test_fda_library_subdir_env_override_feeds_runtime_config(monkeypatch):
    monkeypatch.setenv("TEST_FDA_LIBRARY_SUBDIR", "dry_bench_64")
    cfg = load_inputs()

    assert cfg["TEST_FDA_LIBRARY_SUBDIR"] == "dry_bench_64"


def test_test_fda_override_retargets_default_test_map_for_dud_stream(tmp_path):
    base = tmp_path / "ligands"
    (base / "dry_bench_64").mkdir(parents=True)
    (base / "fda_test_library_10").mkdir(parents=True)
    cfg = {
        "OUTPUT_LIGANDS_DIR": str(base),
        "LIBRARY_SUBDIR_DEFAULT": "dry_bench_64",
        "TEST_FDA_LIBRARY_SUBDIR": "dry_bench_64",
        "TEST_MODE_ENABLE": "dud+fda",
        "TEST_LIBRARY_MAP": {"TEST": "fda_test_library_10"},
        "_TEST_LIBRARY_MAP_CACHE_KEY": ("dict", "stale"),
        "_TEST_LIBRARY_MAP_CACHE_VAL": {"TEST": "fda_test_library_10"},
        "_ALLOWED_LIBRARY_ROOTS_CACHE": {"stale": [str(base / "fda_test_library_10")]},
        "_LIB_ROOTS_MAP_LOG_FINGERPRINT": "dict:1",
    }

    assert _retarget_default_test_library_map_for_test_fda(cfg, "dry_bench_64") is True
    assert cfg["TEST_LIBRARY_MAP"]["TEST"] == "dry_bench_64"
    assert "_TEST_LIBRARY_MAP_CACHE_KEY" not in cfg
    assert "_TEST_LIBRARY_MAP_CACHE_VAL" not in cfg
    assert "_ALLOWED_LIBRARY_ROOTS_CACHE" not in cfg

    roots = compute_allowed_library_roots(
        cfg,
        "TEST",
        logging.getLogger("test.library_roots"),
        tokens_override=["dud"],
    )
    assert roots == [base / "dry_bench_64"]


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


def test_run_scoped_receptor_fallback_resolves_unscoped_output_dir(tmp_path):
    run_id = "run123"
    receptor = (
        tmp_path
        / "outputs"
        / "processed_pdbs"
        / run_id
        / "TEST"
        / "HOLO"
        / "receptor"
        / "TEST.pdbqt"
    )
    receptor.parent.mkdir(parents=True)
    receptor.write_text("REMARK receptor\n", encoding="utf-8")

    resolved = _run_scoped_receptor_fallback(
        {
            "RUN_ID": run_id,
            "OUTPUT_DIR": str(tmp_path / "outputs" / "processed_pdbs"),
        },
        "TEST",
        "HOLO",
        None,
    )

    assert resolved == receptor
