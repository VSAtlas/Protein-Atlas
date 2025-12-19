from __future__ import annotations

import csv
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Sequence
from unittest.mock import patch

import pytest

# Repository root (chemdb/tests -> chemdb -> repo root)
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
TESTS_ROOT = Path(__file__).resolve().parent
if str(TESTS_ROOT) not in sys.path:
    sys.path.insert(0, str(TESTS_ROOT))

MAIN_PY = REPO_ROOT / "main.py"
RUN_ID = "test_run"

# Ensure router helpers see a stable run id in this process.
os.environ.setdefault("ATLAS_RUN_ID", RUN_ID)

from apo_holo_mode import resolve_apo_holo_mode  # noqa: E402
from docking_ligands import (  # noqa: E402
    _coerce_test_map,
    _is_under_ph_subdir,
    _lib_roots_for_pdb,
    _resolve_test_mode,
)
from input_and_export_functions import load_inputs  # noqa: E402
from path_router import docked_dir, load_ph_tags, make_paths  # noqa: E402
from ph_ensemble_docking import enumerate_ligands_for_ph_context  # noqa: E402
from testdata_hmdb_test_library import ensure_hmdb_test_library, infer_ligand_roots  # noqa: E402


def _build_env(tmp_path: Path, extra_env: dict[str, str] | None = None) -> dict[str, str]:
    env = os.environ.copy()
    base = {
        "PYTHONUNBUFFERED": "1",
        "DOCKED_DIR": tmp_path / "docked",
        "CONFIGS_DIR": tmp_path / "configs",
        "P2RANK_OUTPUT_DIR": tmp_path / "p2rank_out",
        "LIGAND_EXTRACTED_DIR": tmp_path / "extracted_ligands",
        "RUN_ID": RUN_ID,
        "CPU": "1",
        "MAX_PARALLEL_JOBS": "1",
        "FORCE_REPROCESS": "0",
        # Keep integration runs small by forcing the tiny test FDA library.
        "TEST_FDA_LIBRARY_SUBDIR": "fda_test_library_10",
        "LIBRARY_SUBDIR_DEFAULT": "fda_test_library_10",
        # Ensure test mode stays on (allows dud/fda small libraries).
        "TEST_MODE_ENABLE": "dud+fda",
        # Let downstream logic know we're in pytest-driven TEST_LIBRARY_MAP flows.
        "PDB_SELECTION_MODE": "TEST_LIBRARY_MAP",
        # Force Vina-only runs during integration tests.
        "USE_GNINA": "false",
    }
    env.update({k: str(v) for k, v in base.items()})
    env["ATLAS_RUN_ID"] = RUN_ID
    if extra_env:
        env.update({k: str(v) for k, v in extra_env.items()})

    for key in ("DOCKED_DIR", "CONFIGS_DIR", "P2RANK_OUTPUT_DIR", "LIGAND_EXTRACTED_DIR"):
        Path(env[key]).mkdir(parents=True, exist_ok=True)
    return env


def _load_cfg(env: dict[str, str]) -> dict:
    with patch.dict(os.environ, env, clear=False):
        return load_inputs()


def _pick_test_pdb(cfg: dict) -> str:
    input_dir = Path(cfg["INPUT_DIR"])
    test_map = _coerce_test_map(cfg.get("TEST_LIBRARY_MAP", {}))
    lig_root = Path(
        cfg.get("OUTPUT_LIGANDS_DIR")
        or cfg.get("PREPPED_LIGANDS_ROOT")
        or REPO_ROOT / "prepped_ligands"
    )

    preferred: list[str] = ["TEST"]
    preferred.extend([pid for pid in test_map.keys() if pid not in preferred])
    preferred.extend(sorted(p.stem.upper() for p in input_dir.glob("*.pdb")))

    for pdb_id in preferred:
        pdb_path = input_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            continue
        lib = test_map.get(pdb_id)
        if lib and not (lig_root / lib).exists():
            continue
        return pdb_id

    pytest.skip("No suitable input PDB found for integration tests.")


def _score_csvs(base_dir: Path, pattern: str = "*docking_score_long.csv") -> list[Path]:
    return [p for p in base_dir.rglob(pattern) if p.is_file()]


def _csv_has_rows(csv_path: Path) -> bool:
    try:
        lines = [ln for ln in csv_path.read_text().splitlines() if ln.strip()]
    except FileNotFoundError:
        return False
    return len(lines) > 1


def _ligands_from_csv(csv_path: Path, stage: str | None = None) -> set[str]:
    ligands: set[str] = set()
    with csv_path.open(newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            return set()

        col_idx = {name: idx for idx, name in enumerate(header)}
        lig_idx = col_idx.get("ligand")
        stage_idx = col_idx.get("stage")
        if lig_idx is None:
            return set()

        for row in reader:
            if stage and stage_idx is not None and row[stage_idx] != stage:
                continue
            ligands.add(Path(row[lig_idx]).stem)
    return ligands


def _write_fake_pdbqt(path: Path) -> Path:
    """Minimal PDBQT writer to satisfy filtering and planned-ligand checks."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["REMARK fake ligand for pH filter test"]
    for i in range(1, 6):
        lines.append(
            f"ATOM  {i:5d}  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run_main_cli(
    tmp_path: Path,
    extra_args: Sequence[str] | None = None,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    env_vars = dict(env) if env is not None else _build_env(tmp_path)
    config_dir = Path(env_vars.get("CONFIGS_DIR", tmp_path / "configs"))
    config_dir.mkdir(parents=True, exist_ok=True)

    # Force fast + test FDA runs to keep integration tests lightweight.
    args = [sys.executable, str(MAIN_PY), "-test", "-fast", "-test-fda", "--run-id", RUN_ID]
    if extra_args:
        cleaned: list[str] = []
        skip = False
        for tok in extra_args:
            if skip:
                skip = False
                continue
            if tok in {"--run-id", "-run-id"}:
                skip = True
                continue
            cleaned.append(tok)
        args.extend(cleaned)
    if "--configs-dir" not in args:
        args.extend(["--configs-dir", str(config_dir)])

    # Tests reuse fixed run_ids; clear any stale main log before invoking the CLI.
    run_id = None
    for idx, token in enumerate(args):
        if token in {"--run-id", "-run-id"} and (idx + 1) < len(args):
            run_id = args[idx + 1]
            break
    if run_id:
        log_path = REPO_ROOT / "logs" / f"main_{run_id}.log"
        try:
            log_path.unlink()
        except FileNotFoundError:
            pass

    cp = subprocess.run(
        args,
        env=env_vars,
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if cp.returncode != 0:
        print("=== main.py stdout/stderr ===")
        print(cp.stdout)
    return cp


def enumerate_ph_ligands_for_test_pdb(
    env: dict[str, str], pdb_id: str | None = None
) -> set[str]:
    with patch.dict(os.environ, env, clear=False):
        cfg = load_inputs()
        # Mirror the CLI -test-fda override used by run_main_cli so enumeration
        # matches the runtime library selection.
        cfg["LIBRARY_SUBDIR_DEFAULT"] = cfg.get("TEST_FDA_LIBRARY_SUBDIR", "fda_test_library_10")
        cfg["PH_LIGAND_MODE"] = "context_window"
        chosen = pdb_id or _pick_test_pdb(cfg)
        ph_tags = load_ph_tags(chosen, variant=None)
        ph_label = ph_tags[0] if ph_tags else None
        paths = make_paths(cfg, base_id=chosen, pdb_file=f"{chosen}.pdb")
        logger = logging.getLogger("test.ph_ligands")
        _, noncontrol_roots = _lib_roots_for_pdb(
            cfg,
            chosen.upper(),
            paths,
            logger,
            test_mode_override=_resolve_test_mode(cfg),
        )
        ph_root = None
        if noncontrol_roots:
            default_root_name = Path(str(cfg.get("LIBRARY_SUBDIR_DEFAULT", ""))).name

            # Prefer the configured default library (e.g., fda_test_library_10) when present.
            for root in noncontrol_roots:
                if default_root_name and Path(root).name == default_root_name:
                    ph_root = root
                    break

            # Otherwise, prefer non-decoy roots (e.g., test_library_10) over decoys_* entries.
            if ph_root is None:
                for root in noncontrol_roots:
                    name = Path(root).name
                    if not name.startswith("decoys_"):
                        ph_root = root
                        break

            if ph_root is None:
                ph_root = noncontrol_roots[0]

        logger.info(
            "[test.ph_ligands] chosen ph_ligand_root=%s (noncontrol_roots=%s)",
            ph_root,
            noncontrol_roots,
        )
        ligs = enumerate_ligands_for_ph_context(
            cfg=cfg,
            pdb_id=chosen,
            ph_label=ph_label,
            ph_ligand_root=ph_root,
        )
    return {Path(p).stem for p in ligs or []}


def _prepare_run(tmp_path: Path, extra_env: dict[str, str] | None = None) -> tuple[dict[str, str], str]:
    env = _build_env(tmp_path, extra_env)
    cfg = _load_cfg(env)
    test_mode_raw = str(env.get("TEST_MODE_ENABLE", "")).lower()
    if "hmdb" in test_mode_raw:
        extracted_root, prepped_root = infer_ligand_roots(cfg)
        ensure_hmdb_test_library(extracted_root, prepped_root)
    pdb_id = _pick_test_pdb(cfg)
    return env, pdb_id


def test_apo_mode_generates_csv(tmp_path: Path) -> None:
    env, pdb_id = _prepare_run(tmp_path, {"APO_HOLO_MODE": "apo"})
    cp = run_main_cli(tmp_path, extra_args=["-fast", "-test-fda","-pdb", pdb_id, "--run-id", "apo_mode"], env=env)
    assert cp.returncode == 0

    with patch.dict(os.environ, env, clear=False):
        cfg = load_inputs()
        mode, variants = resolve_apo_holo_mode(cfg)
        legacy = mode == "legacy"
        assert variants, "No variants resolved for APO mode"
        variant_root = docked_dir(pdb_id, variant=variants[0], ph_tag=None, legacy=legacy)
        csv_paths = _score_csvs(variant_root)
        assert csv_paths, f"No docking_score_long.csv found under {variant_root}"
        assert any(_csv_has_rows(p) for p in csv_paths)
        assert any(_ligands_from_csv(p) for p in csv_paths)


def test_holo_mode_generates_csv(tmp_path: Path) -> None:
    env, pdb_id = _prepare_run(tmp_path, {"APO_HOLO_MODE": "holo"})
    cp = run_main_cli(tmp_path, extra_args=["-pdb", pdb_id, "--run-id", "holo_mode"], env=env)
    assert cp.returncode == 0

    with patch.dict(os.environ, env, clear=False):
        cfg = load_inputs()
        mode, variants = resolve_apo_holo_mode(cfg)
        legacy = mode == "legacy"
        assert variants, "No variants resolved for HOLO mode"
        variant_root = docked_dir(pdb_id, variant=variants[0], ph_tag=None, legacy=legacy)
        csv_paths = _score_csvs(variant_root)
        assert csv_paths, f"No docking_score_long.csv found under {variant_root}"
        assert any(_csv_has_rows(p) for p in csv_paths)
        assert any(_ligands_from_csv(p) for p in csv_paths)


def test_apo_vs_holo_mode_generates_csv_per_variant(tmp_path: Path) -> None:
    env, pdb_id = _prepare_run(tmp_path, {"APO_HOLO_MODE": "apo_vs_holo"})
    cp = run_main_cli(
        tmp_path,
        extra_args=["-pdb", pdb_id, "--run-id", "apo_vs_holo_mode"],
        env=env,
    )
    assert cp.returncode == 0

    with patch.dict(os.environ, env, clear=False):
        cfg = load_inputs()
        mode, variants = resolve_apo_holo_mode(cfg)
        legacy = mode == "legacy"
        assert variants, "No variants resolved for apo_vs_holo mode"

        for variant in variants:
            variant_root = docked_dir(pdb_id, variant=variant, ph_tag=None, legacy=legacy)
            csv_paths = _score_csvs(variant_root)
            assert csv_paths, f"Missing docking_score_long.csv for variant {variant}"
            assert any(_csv_has_rows(p) for p in csv_paths)
            assert any(_ligands_from_csv(p) for p in csv_paths)


def test_ph_ligand_mode_context_window_ligands_match_enumeration(tmp_path: Path) -> None:
    env, pdb_id = _prepare_run(
        tmp_path,
        {
            "PH_LIGAND_MODE": "context_window",
            "APO_HOLO_MODE": "apo",
        },
    )
    cp = run_main_cli(
        tmp_path,
        extra_args=[
            "-pdb",
            pdb_id,
            "--ph-ligand-mode",
            "context_window",
            "--run-id",
            "ph_context_window",
        ],
        env=env,
    )
    assert cp.returncode == 0

    expected = enumerate_ph_ligands_for_test_pdb(env, pdb_id=pdb_id)
    assert expected, "PH ligand enumeration returned no ligands for context_window mode"

    with patch.dict(os.environ, env, clear=False):
        cfg = load_inputs()
        mode, variants = resolve_apo_holo_mode(cfg)
        legacy = mode == "legacy"
        variant_root = docked_dir(pdb_id, variant=variants[0], ph_tag=None, legacy=legacy)
        csv_paths = _score_csvs(variant_root)
        assert csv_paths, "No docking_score_long.csv produced for PH ligand mode"

        actual: set[str] = set()
        for csv_path in csv_paths:
            actual.update(_ligands_from_csv(csv_path, stage="stage1"))

    assert actual == expected


def test_no_docking_mode_creates_planned_ligands_only(tmp_path: Path) -> None:
    env, pdb_id = _prepare_run(tmp_path, {"NO_LIBRARY_DOCKING": "1"})
    cp = run_main_cli(
        tmp_path,
        extra_args=["-pdb", pdb_id, "--no-docking", "--run-id", "no_docking_mode"],
        env=env,
    )
    assert cp.returncode == 0

    with patch.dict(os.environ, env, clear=False):
        cfg = load_inputs()
        mode, variants = resolve_apo_holo_mode(cfg)
        legacy = mode == "legacy"

        for variant in variants:
            variant_root = docked_dir(pdb_id, variant=variant, ph_tag=None, legacy=legacy)
            planned = list(variant_root.rglob("planned_ligands_*.txt"))
            assert planned, f"No planned_ligands_*.txt files for variant {variant}"
            assert all(p.stat().st_size > 0 for p in planned)

            for stage_dir in variant_root.rglob("stage*"):
                assert not any(stage_dir.glob("*_stage*.pdbqt")), (
                    f"Found docking outputs in no-docking mode under {stage_dir}"
                )


def test_main_fast_no_docking_excludes_ph_and_microstate_dirs(tmp_path: Path) -> None:
    pdb_id = "1BCD"
    env = _build_env(
        tmp_path,
        {
            "OVERALL_DIR": tmp_path,
            "INPUT_DIR": REPO_ROOT / "input_pdbs",
            "OUTPUT_DIR": tmp_path / "processed_pdbs",
            "PDBQT_DIR": tmp_path / "pdbqts",
            "OUTPUT_LIGANDS_DIR": tmp_path / "prepped_ligands",
            "ONLY_PDBS": pdb_id,
            "PH_LIGAND_MODE": "off",
        },
    )

    with patch.dict(os.environ, env, clear=False):
        cfg_for_roots = load_inputs()
    test_map = _coerce_test_map(cfg_for_roots.get("TEST_LIBRARY_MAP", {}))
    mapped = test_map.get(pdb_id)
    base_root = Path(cfg_for_roots["OUTPUT_LIGANDS_DIR"])
    library_root = base_root / (mapped or cfg_for_roots.get("LIBRARY_SUBDIR_DEFAULT", "fda_library"))
    library_root.mkdir(parents=True, exist_ok=True)

    lig_root = library_root
    base_lig = _write_fake_pdbqt(lig_root / "base_a.pdbqt")
    _write_fake_pdbqt(lig_root / "7.4" / "ph_dir_a.pdbqt")
    _write_fake_pdbqt(lig_root / "microstates" / "micro_a.pdbqt")

    cp = run_main_cli(
        tmp_path,
        extra_args=["-pdb", pdb_id, "--no-docking"],
        env=env,
    )
    assert cp.returncode == 0

    with patch.dict(os.environ, env, clear=False):
        cfg = load_inputs()
        mode, variants = resolve_apo_holo_mode(cfg)
        legacy = mode == "legacy"
        assert variants, "No variants resolved for PH-off filtering check"
        variant_root = docked_dir(pdb_id, variant=variants[0], ph_tag=None, legacy=legacy)

    planned_lists = list(variant_root.rglob("planned_ligands_*.txt"))
    assert planned_lists, "No planned ligands written for PH-off filter check"

    planned_paths: list[Path] = []
    for txt in planned_lists:
        planned_paths.extend(
            [Path(ln.strip()) for ln in txt.read_text().splitlines() if ln.strip()]
        )

    assert base_lig.resolve() in {p.resolve() for p in planned_paths}
    ph_roots = [lig_root]
    is_micro = lambda p: any(part.lower() == "microstates" for part in p.parts)
    assert not any(is_micro(p) for p in planned_paths)
    assert not any(_is_under_ph_subdir(p, ph_roots) for p in planned_paths)


TEST_MODES_AND_PREFIXES = [
    ("off", [""]),
    ("dud", [""]),
    ("hmdb", ["hmdb_"]),
    ("hmdb+dud", ["hmdb_", "dud_"]),
    ("hmdb+fda", ["hmdb_", ""]),
    ("fda+dud", ["dud_", ""]),
    ("fda+dud+hmdb", ["dud_", "hmdb_", ""]),
]


@pytest.mark.parametrize("mode,expected_prefixes", TEST_MODES_AND_PREFIXES)
def test_test_mode_enable_generates_prefixed_csvs(
    tmp_path: Path, mode: str, expected_prefixes: Iterable[str]
) -> None:
    env, pdb_id = _prepare_run(
        tmp_path,
        {
            "TEST_MODE_ENABLE": mode,
            "NO_LIBRARY_DOCKING": "1",
        },
    )
    cp = run_main_cli(
        tmp_path,
        extra_args=["-pdb", pdb_id, "--no-docking", "--run-id", f"test_mode_{mode}"],
        env=env,
    )
    assert cp.returncode == 0

    with patch.dict(os.environ, env, clear=False):
        cfg = load_inputs()
        mode_token, variants = resolve_apo_holo_mode(cfg)
        legacy = mode_token == "legacy"

        for variant in variants:
            variant_root = docked_dir(pdb_id, variant=variant, ph_tag=None, legacy=legacy)
            for prefix in expected_prefixes:
                long_csv = _score_csvs(variant_root, pattern=f"{prefix}docking_score_long.csv")
                summary_csv = _score_csvs(variant_root, pattern=f"{prefix}docking_score_summary.csv")
                if long_csv or summary_csv:
                    assert all(_csv_has_rows(p) for p in long_csv), "Empty docking_score_long.csv detected"
                    assert all(p.stat().st_size > 0 for p in summary_csv), "Empty docking_score_summary.csv detected"
                    continue

                suffix = "run"
                if prefix == "dud_":
                    suffix = "dud"
                elif prefix == "hmdb_":
                    suffix = "hmdb"
                elif "fda" in mode.lower():
                    suffix = "fda"

                planned = list(variant_root.rglob(f"planned_ligands_{suffix}.txt"))
                assert planned, f"Missing planned ligands list for mode={mode} prefix={prefix} variant={variant}"
                assert all(p.stat().st_size > 0 for p in planned), (
                    f"Planned ligands list empty for mode={mode} prefix={prefix} variant={variant}"
                )
