import os
import sys
import subprocess
from pathlib import Path

import pytest
import yaml

PDBS_OF_INTEREST = [
    "1B9V",
    "1BCD",
    "1C8K",
    "1D3G",
    "1E66",
    "1H00",
    "1J4H",
    "1KVO",
    "1L2S",
    "1LRU",
    "1MV9",
    "1NJS",
    "1Q4X",
    "1S3B",
    "1SJ0",
    "1SQT",
    "1SYN",
    "1UYG",
    "1XL2",
    "1YPE",
    "2AA2",
    "2AM9",
    "2AYW",
    "2B8T",
    "2CNK",
    "2ETR",
    "2H7L",
    "2I0E",
    "2I78",
    "2OF2",
    "2OWB",
    "TEST",
]


@pytest.mark.slow
def test_control_redock_sets_method_control_for_dud_targets(tmp_path):
    root = Path(__file__).resolve().parents[2]
    manifests_dir = root / "manifests"
    manifests_dir.mkdir(exist_ok=True)

    before = {d.name for d in manifests_dir.iterdir() if d.is_dir()}

    pdbs_arg = " ".join(PDBS_OF_INTEREST)
    cmd = [
        sys.executable,
        os.fspath(root / "main.py"),
        "-fast",
        "-pdbs",
        pdbs_arg,
        "--single",
        "dexamethasone",
    ]

    result = subprocess.run(
        cmd,
        cwd=root,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, (
        "main.py failed with non-zero exit code\n"
        f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
    )

    after = {d.name for d in manifests_dir.iterdir() if d.is_dir()}
    new_ids = sorted(after - before)
    assert len(new_ids) == 1, f"Expected exactly one new manifest folder, got: {new_ids}"
    run_id = new_ids[0]

    manifest_path = manifests_dir / run_id / "run_manifest.yaml"
    assert manifest_path.is_file(), f"Manifest not found at {manifest_path}"

    with manifest_path.open("r") as f:
        data = yaml.safe_load(f) or {}
    proteins = data.get("proteins", {}) or {}

    required = set(PDBS_OF_INTEREST)
    seen_control = {pdb_id: False for pdb_id in required}

    for key, entry in proteins.items():
        pdb_id = entry.get("pdb_id")
        if pdb_id not in required:
            continue

        status = entry.get("status")
        stages = entry.get("stages", {}) or {}
        pocket = stages.get("pocket_detection", {}) or {}
        details = pocket.get("details", {}) or {}
        method = details.get("method")

        if status == "completed" and method == "control":
            seen_control[pdb_id] = True

    missing = [p for p, ok in seen_control.items() if not ok]

    assert not missing, (
        "Expected control redocking to be used for pocket_detection "
        "for all DUD targets, but these PDBs had no completed entry "
        "with method='control': "
        + ", ".join(sorted(missing))
        + f"\nRun ID: {run_id}"
    )
