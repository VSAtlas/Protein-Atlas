from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest
import yaml

import pocket_policy
from ligand_pocket import extract_and_remove_ligands

REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN_PY = REPO_ROOT / "main.py"
PDB_ID = "6LU7"
INPUT_PDB = REPO_ROOT / "input_pdbs" / f"{PDB_ID}.pdb"
POCKETS_DIR = REPO_ROOT / "processed_pdbs" / PDB_ID / "pockets"
MANIFESTS_DIR = REPO_ROOT / "manifests"

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


def _expected_center_box(tmp_path: Path):
    cleaned_path = tmp_path / f"{PDB_ID}_cleaned.pdb"
    ligands_dir = tmp_path / "ligands"
    ligands, _ = extract_and_remove_ligands(
        str(INPUT_PDB), str(cleaned_path), str(ligands_dir)
    )
    center, box_size, source, _ = pocket_policy.select_pocket(
        str(cleaned_path),
        str(INPUT_PDB),
        PDB_ID,
        ligands,
        ligands_dir,
        logging.getLogger("pocket-eval-test"),
    )
    return center, box_size, source


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


def test_pocket_eval_flag_off_preserves_selection(tmp_path: Path) -> None:
    if not MAIN_PY.exists() or not INPUT_PDB.exists():
        pytest.skip("Missing main.py or input PDB for pocket eval off test.")

    run_id = f"pocket_eval_off_{uuid.uuid4().hex[:8]}"
    env = {"POCKET_EVAL": "false", **MMGBSA_DISABLE_ENV, **ENGINE_DISABLE_ENV}

    _run_main(run_id, env)

    perf_path = (
        REPO_ROOT
        / "docked"
        / run_id
        / PDB_ID
        / "pocket_eval"
        / "pocket_performance.json"
    )
    assert (
        not perf_path.exists()
    ), "pocket_performance.json should not be created when POCKET_EVAL=false"
    pockets_json_path = POCKETS_DIR / "pockets.json"
    assert pockets_json_path.is_file(), "pockets.json missing for POCKET_EVAL=false test"
    _ = json.loads(pockets_json_path.read_text(encoding="utf-8"))

    expected_center, expected_box, _source = _expected_center_box(tmp_path)
    details = _read_pocket_details(run_id)
    assert details, "Missing pocket_detection details in run manifest"

    actual_center = details.get("center")
    actual_box = details.get("box_size")
    assert actual_center == _round_center(expected_center)
    assert actual_box == _round_box(expected_box)
    assert details.get("method") != "pocket_eval_override"
