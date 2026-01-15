from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import pytest
import yaml

from ligand_pocket import compute_box_from_ligand_coords
from pocket_eval import pocket_id_from_entry

REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN_PY = REPO_ROOT / "main.py"
PDB_ID = "6LU7"
INPUT_PDB = REPO_ROOT / "input_pdbs" / f"{PDB_ID}.pdb"
POCKETS_DIR = REPO_ROOT / "processed_pdbs" / PDB_ID / "pockets"
POCKETS_JSON = POCKETS_DIR / "pockets.json"
MANIFESTS_DIR = REPO_ROOT / "manifests"
FIXTURE_CACHE = REPO_ROOT / "chemdb" / "tests" / "data" / "calibrator" / "6LU7_calibrator_cache.json"

MMGBSA_DISABLE_ENV = {
    "MMGBSA_ENABLED": "false",
    "MMGBSA_MD_ENABLED": "false",
    "MMGBSA_MMPBSA_ENABLED": "false",
    "MMGBSA_MMPBSA_RUN": "false",
    "MD_FIVE_REPLICATE": "false",
}
ENGINE_DISABLE_ENV = {
    "USE_GNINA": "false",
    "USE_LEDOCK": "false",
    "USE_DOCK6": "false",
    "USE_SCORCH": "false",
}


def _round_center(center):
    return [round(float(x), 3) for x in center]


def _round_box(box):
    return [round(float(x), 1) for x in box]


def _center_and_box_from_pocket(entry):
    coords = entry.get("coords") or []
    if coords:
        parsed = []
        for xyz in coords:
            if not isinstance(xyz, (list, tuple)) or len(xyz) < 3:
                continue
            parsed.append((float(xyz[0]), float(xyz[1]), float(xyz[2])))
        if parsed:
            return compute_box_from_ligand_coords(parsed)

    center = entry.get("center")
    bounds = entry.get("bounds") or {}
    min_b = bounds.get("min") if isinstance(bounds, dict) else None
    max_b = bounds.get("max") if isinstance(bounds, dict) else None
    if center and min_b and max_b:
        cx, cy, cz = (float(center[0]), float(center[1]), float(center[2]))
        sx = float(max_b[0]) - float(min_b[0]) + 5.0
        sy = float(max_b[1]) - float(min_b[1]) + 5.0
        sz = float(max_b[2]) - float(min_b[2]) + 5.0
        return (cx, cy, cz), (sx, sy, sz)
    return None, None


def _read_pocket_details(run_id: str):
    manifest_path = MANIFESTS_DIR / run_id / "run_manifest.yaml"
    assert manifest_path.is_file(), f"Manifest not found at {manifest_path}"
    with manifest_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    proteins = data.get("proteins", {}) or {}
    for entry in proteins.values():
        if entry.get("pdb_id") != PDB_ID:
            continue
        stages = entry.get("stages", {}) or {}
        pocket = stages.get("pocket_detection", {}) or {}
        return pocket.get("details", {}) or {}
    return {}


def _run_main(run_id: str, env_overrides: dict[str, str]) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env.update(env_overrides)
    env.setdefault("PYTHONUNBUFFERED", "1")
    cmd = [
        sys.executable,
        "main.py",
        "-6lu7",
        "-fast",
        "--test-fda",
        "--run-id",
        run_id,
    ]
    return subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )


def _best_pocket_id(perf_pockets):
    def _metric(value, default):
        try:
            val = float(value)
        except Exception:
            return default
        if math.isnan(val):
            return default
        return val

    def _key(entry):
        metrics = entry.get("metrics") or {}
        auc = _metric(metrics.get("auc_mean"), float("-inf"))
        ef = _metric(metrics.get("ef1_mean"), float("-inf"))
        std = _metric(metrics.get("auc_std"), float("inf"))
        return (-auc, -ef, std, entry.get("pocket_id", ""))

    ranked = sorted(perf_pockets, key=_key)
    return ranked[0].get("pocket_id") if ranked else None


def test_pocket_eval_flag_on_selects_best_pocket(tmp_path: Path) -> None:
    if not MAIN_PY.exists() or not INPUT_PDB.exists():
        pytest.skip("Missing main.py or input PDB for pocket eval on test.")
    if not FIXTURE_CACHE.exists():
        pytest.skip("Missing calibrator cache fixture for pocket eval on test.")

    run_id = f"pocket_eval_on_{uuid.uuid4().hex[:8]}"
    dock_root = REPO_ROOT / "docked" / run_id / PDB_ID / "pocket_eval"
    dock_root.mkdir(parents=True, exist_ok=True)
    cache_target = dock_root / "calibrator_cache.json"
    shutil.copyfile(FIXTURE_CACHE, cache_target)
    perf_path = dock_root / "pocket_performance.json"
    env = {
        "POCKET_EVAL": "true",
        "POCKET_EVAL_MAX_CALIBRATORS": "6",
        "POCKET_EVAL_FOLDS": "2",
        "POCKET_EVAL_SEED": "0",
        **MMGBSA_DISABLE_ENV,
        **ENGINE_DISABLE_ENV,
    }

    _run_main(run_id, env)

    assert perf_path.is_file(), "pocket_performance.json missing when POCKET_EVAL=true"
    perf = json.loads(perf_path.read_text(encoding="utf-8"))
    perf_pockets = perf.get("pockets") or []
    assert perf_pockets, "No pockets recorded in pocket_performance.json"

    assert POCKETS_JSON.is_file(), "pockets.json missing for pocket eval test"
    pockets_payload = json.loads(POCKETS_JSON.read_text(encoding="utf-8"))
    pockets_list = pockets_payload.get("pockets") or []
    assert len(perf_pockets) <= len(pockets_list)
    assert len(perf_pockets) >= 2

    for entry in perf_pockets:
        metrics = entry.get("metrics") or {}
        auc_mean = metrics.get("auc_mean")
        assert isinstance(auc_mean, (int, float))
        assert not math.isnan(float(auc_mean))
        assert int(metrics.get("n_strong", 0)) > 0
        assert int(metrics.get("n_non", 0)) > 0
        artifacts = entry.get("artifacts") or {}
        dock_dir = Path(artifacts.get("dock_dir", ""))
        scores_path = Path(artifacts.get("scores_path", ""))
        assert dock_dir.is_dir()
        assert scores_path.is_file()
        assert scores_path.stat().st_size > 0

    selected_id = perf.get("selected_pocket_id")
    assert selected_id, "selected_pocket_id missing in pocket_performance.json"
    assert perf.get("selected_reason"), "selected_reason missing in pocket_performance.json"

    expected_best = _best_pocket_id(perf_pockets)
    assert expected_best == selected_id

    pocket_map = {}
    for entry in pockets_list:
        pid = pocket_id_from_entry(entry)
        center, box = _center_and_box_from_pocket(entry)
        if center and box:
            pocket_map[pid] = (center, box)

    expected_center, expected_box = pocket_map.get(expected_best, (None, None))
    assert expected_center and expected_box

    details = _read_pocket_details(run_id)
    assert details, "Missing pocket_detection details in run manifest"
    assert details.get("center") == _round_center(expected_center)
    assert details.get("box_size") == _round_box(expected_box)
