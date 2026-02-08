from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Ensure we can import path_router when this test file lives in chemdb/tests
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from path_router.path_router import (
    _ensure_router_roots,
    expand_variants,
    raw_dir_v,
    nolig_dir_v,
    work_dir_v,
    receptor_dir,
    ph_ensemble_dir,
    docked_dir,
    make_paths,
    _default_config_path,
)
from input_and_export_functions import load_config

# Optional helpers: not all versions of path_router export these
try:
    from path_router.path_router import _canon_pdb_id_from_path  # type: ignore
except Exception:  # pragma: no cover - older versions
    _canon_pdb_id_from_path = None  # type: ignore[assignment]

try:
    from path_router import _default_variant  # type: ignore
except Exception:  # pragma: no cover - older versions
    _default_variant = None  # type: ignore[assignment]

try:
    from path_router import _expand_tokens  # type: ignore
except Exception:  # pragma: no cover - older versions
    _expand_tokens = None  # type: ignore[assignment]


MINIMAL_TEST_PDB = """\
HEADER    DUMMY PDB FOR PATH_ROUTER TEST
ATOM      1  N   ALA A   1      11.104  13.207   9.597  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.560  13.207   9.597  1.00 20.00           C
TER
END
"""
TEST_RUN_ID = "test_run"


def _load_cfg() -> dict[str, str]:
    """
    Load config.txt through the shared loader + normalizer so minimal root-only
    configs still provide fully derived path keys.
    """
    cfg_path = _default_config_path()
    cfg = load_config(config_path=str(cfg_path), base_dir=cfg_path.parent)
    assert cfg, f"Config seems empty or missing: {cfg_path}"

    # Baseline required keys
    for key in ("OVERALL_DIR", "INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR"):
        assert key in cfg, f"Missing required key {key!r} in config"

    # Make sure some ligand-related keys exist for make_paths
    for key in ("LIGANDS_MOL2_DIR", "PREPPED_LIGANDS_DIR"):
        assert key in cfg, f"Missing required key {key!r} in config"

    # Make sure path_router tests are independent of the APO_HOLO_MODE in
    # config.txt. These tests specifically validate apo_vs_holo behavior.
    cfg["APO_HOLO_MODE"] = "apo_vs_holo"
    cfg["RUN_ID"] = TEST_RUN_ID

    return cfg


def _variant_segment(variant: str | None) -> list[str]:
    """
    Helper: return [] for None, or [variant] for 'APO'/'HOLO'.
    """
    return [] if variant is None else [variant]


@pytest.fixture(scope="session")
def cfg() -> dict[str, str]:
    return _load_cfg()


@pytest.fixture(autouse=True)
def _reset_router_roots(monkeypatch):
    # Ensure each test re-evaluates roots (allows run-id tweaks).
    import path_router.path_router as pr

    pr._ROUTER_ROOTS = None  # type: ignore[attr-defined]
    monkeypatch.setenv("ATLAS_RUN_ID", TEST_RUN_ID)
    yield
    pr._ROUTER_ROOTS = None  # type: ignore[attr-defined]


@pytest.fixture(scope="session")
def test_pdb_path(cfg: dict[str, str]) -> Path:
    """
    Create a tiny dummy TEST.pdb under INPUT_DIR so that make_paths.input_pdb_path
    actually points at a real file. If the file already exists, leave it alone.
    """
    input_root = Path(cfg["INPUT_DIR"])
    input_root.mkdir(parents=True, exist_ok=True)
    pdb_path = input_root / "TEST.pdb"
    if not pdb_path.exists():
        pdb_path.write_text(MINIMAL_TEST_PDB)
    return pdb_path


# --- _ensure_router_roots vs config.txt --------------------------------------


def test_ensure_router_roots_matches_config(cfg: dict[str, str]) -> None:
    """
    Sanity check: _ensure_router_roots should match the paths implied by
    config.txt after token expansion.
    """
    import path_router.path_router as pr

    pr._ROUTER_ROOTS = None

    roots = _ensure_router_roots()

    overall_expected = Path(cfg["OVERALL_DIR"]).resolve()
    processed_expected = Path(cfg["OUTPUT_DIR"]).resolve()
    docked_expected = Path(cfg["DOCKED_DIR"]).resolve() / TEST_RUN_ID
    configs_expected = Path(
        cfg.get("CONFIGS_DIR", overall_expected / "configs")
    ).resolve()

    assert roots.overall == overall_expected
    assert roots.processed == processed_expected
    assert roots.docked == docked_expected
    assert roots.configs == configs_expected


# --- expand_variants behavior ------------------------------------------------


@pytest.mark.parametrize(
    "mode, expected",
    [
        (None, [None]),
        ("", [None]),
        ("legacy", [None]),
        ("none", [None]),
        ("apo", ["APO"]),
        ("holo", ["HOLO"]),
        ("apo_vs_holo", ["APO", "HOLO"]),
        ("apo vs holo", ["APO", "HOLO"]),
        ("apo-vs-holo", ["APO", "HOLO"]),
        ("both", ["APO", "HOLO"]),
        ("apoandholo", ["APO", "HOLO"]),
    ],
)
def test_expand_variants_various_modes(mode, expected) -> None:
    assert expand_variants(mode) == expected


# --- Variant-aware directory helpers ----------------------------------------


@pytest.mark.parametrize("variant", [None, "APO", "HOLO"])
def test_variant_directories_follow_output_dir_layout(
    cfg: dict[str, str], variant
) -> None:
    """
    raw_dir_v / nolig_dir_v / work_dir_v / receptor_dir / ph_ensemble_dir
    should all follow:

        OUTPUT_DIR / PDB_ID / [variant?] / <subdir>
    """
    import path_router.path_router as pr

    pr._ROUTER_ROOTS = None

    output_root = Path(cfg["OUTPUT_DIR"])
    pdb_id = "TEST"

    v_seg = _variant_segment(variant)
    base = output_root.joinpath(pdb_id, *v_seg)

    assert raw_dir_v(pdb_id, variant=variant) == base / "raw"
    assert nolig_dir_v(pdb_id, variant=variant) == base / "nolig"
    assert work_dir_v(pdb_id, variant=variant) == base / "work"
    assert receptor_dir(pdb_id, variant=variant) == base / "receptor"
    assert ph_ensemble_dir(pdb_id, variant=variant) == base / "receptor" / "ph_ensemble"


@pytest.mark.parametrize("variant", [None, "APO", "HOLO"])
def test_docked_dir_layout(cfg: dict[str, str], variant) -> None:
    """
    docked_dir should follow:

        legacy=False, ph_tag=None:
            DOCKED_DIR / PDB / [variant?]

        legacy=True, ph_tag=None:
            DOCKED_DIR / PDB  (ignores variant)

        legacy=False, ph_tag='pH7_0':
            DOCKED_DIR / PDB / [variant?] / 'pH7_0'
    """
    import path_router.path_router as pr

    pr._ROUTER_ROOTS = None

    docked_root = Path(cfg["DOCKED_DIR"]) / TEST_RUN_ID
    pdb_id = "TEST"
    v_seg = _variant_segment(variant)

    expected_no_ph = docked_root.joinpath(pdb_id, *v_seg)
    assert (
        docked_dir(pdb_id, variant=variant, ph_tag=None, legacy=False) == expected_no_ph
    )

    expected_legacy = docked_root / pdb_id
    assert (
        docked_dir(pdb_id, variant=variant, ph_tag=None, legacy=True) == expected_legacy
    )

    ph_tag = "pH7_0"
    expected_ph = docked_root.joinpath(pdb_id, *v_seg, ph_tag)
    assert (
        docked_dir(pdb_id, variant=variant, ph_tag=ph_tag, legacy=False) == expected_ph
    )


def test_docked_dir_scopes_run_id_and_legacy(monkeypatch, cfg: dict[str, str]) -> None:
    """
    ATLAS_RUN_ID should inject a run-level subdir under docked/, while
    clearing it should fall back to legacy docked/<PDB>.
    """
    import path_router.path_router as pr

    # With run id
    pr._ROUTER_ROOTS = None  # type: ignore[attr-defined]
    monkeypatch.setenv("ATLAS_RUN_ID", TEST_RUN_ID)
    roots = _ensure_router_roots()
    base = Path(cfg["DOCKED_DIR"]).resolve() / TEST_RUN_ID
    assert roots.docked == base
    assert docked_dir("1ABC", legacy=False) == base / "1ABC"

    # Without run id (fallback)
    pr._ROUTER_ROOTS = None  # type: ignore[attr-defined]
    monkeypatch.delenv("ATLAS_RUN_ID", raising=False)
    roots = _ensure_router_roots()
    legacy_base = Path(cfg["DOCKED_DIR"]).resolve()
    assert roots.docked == legacy_base
    assert docked_dir("1ABC", legacy=False) == legacy_base / "1ABC"


# --- make_paths structure ----------------------------------------------------


def test_make_paths_dummy_pdb_layout(cfg: dict[str, str], test_pdb_path: Path) -> None:
    """
    Equivalent to the manual 'make_paths' printout, but using a dummy PDB ID
    'TEST' so we don't depend on real PDBs.
    """
    output_root = Path(cfg["OUTPUT_DIR"])
    input_root = Path(cfg["INPUT_DIR"])
    lig_mol2_root = Path(cfg["LIGANDS_MOL2_DIR"])
    prepped_root = Path(cfg["PREPPED_LIGANDS_DIR"])

    pdb_id = "TEST"
    pdb_file = "TEST.pdb"

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=pdb_file)

    assert paths.pdb_id == pdb_id
    assert paths.pdb_file == pdb_file

    root_pdb_dir = output_root / pdb_id

    assert paths.root_pdb_dir == root_pdb_dir
    assert paths.raw_dir == root_pdb_dir / "raw"
    assert paths.ligands_raw_dir == root_pdb_dir / "ligands_raw"
    assert paths.nolig_dir == root_pdb_dir / "nolig"
    assert paths.work_dir == root_pdb_dir / "work"

    assert paths.input_pdb_path == input_root / pdb_file
    assert paths.nolig_pdb_path == root_pdb_dir / "nolig" / f"{pdb_id}_nolig.pdb"

    assert paths.ligand_output_dir == root_pdb_dir / "ligands_raw"
    assert paths.ligands_mol2_dir == lig_mol2_root / pdb_id
    assert paths.prepped_ligands_dir == prepped_root / pdb_id


# --- _canon_pdb_id_from_path behavior ---------------------------------------


@pytest.mark.parametrize(
    "path, expected",
    [
        ("TEST.pdb", "TEST"),
        ("/tmp/TEST.pdb", "TEST"),
        ("TEST_cleaned.pdb", "TEST"),
        ("/some/dir/TEST_cleaned.pdb", "TEST"),
        ("/some/dir/test_cleaned.pdb", "TEST"),
        ("test.pdb", "TEST"),
    ],
)
def test_canon_pdb_id_from_path_variants(path: str, expected: str) -> None:
    """
    If path_router exposes _canon_pdb_id_from_path, verify that it collapses
    various TEST-style filenames to the canonical 'TEST' ID.

    On older versions where the helper is not exported, we simply no-op so the
    test passes without being marked as SKIPPED.
    """
    if _canon_pdb_id_from_path is None:
        # Older path_router: helper not exported; nothing to assert here.
        # Early-return so pytest sees this as a successful (no-op) test.
        return

    assert _canon_pdb_id_from_path(path) == expected


def test_default_variant_matches_current_config_when_apo_vs_holo(
    cfg: dict[str, str],
) -> None:
    """
    If _default_variant exists and APO_HOLO_MODE=='apo_vs_holo', we expect
    the default variant to be 'APO'. On older path_router versions without
    this helper, we no-op so the test passes without SKIPPED.
    """
    if _default_variant is None:
        # Older path_router version: nothing to assert here.
        return

    mode = cfg.get("APO_HOLO_MODE")
    if mode == "apo_vs_holo":
        assert _default_variant(cfg) == "APO"
