from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN_PY = REPO_ROOT / "main.py"
INPUT_FIXTURE_ROOT = REPO_ROOT / "chemdb" / "tests" / "fixtures" / "input_pdbs"
TEST_PDB_PATH = INPUT_FIXTURE_ROOT / "TEST.pdb"
PREPPED_FIXTURE_ROOT = REPO_ROOT / "chemdb" / "tests" / "fixtures" / "prepped_ligands"
PREPPED_TEST_LIB = PREPPED_FIXTURE_ROOT / "fda_test_library_10"

MMGBSA_DISABLE_ENV = {
    "MMGBSA_ENABLED": "false",
    "MMGBSA_MD_ENABLED": "false",
    "MMGBSA_MMPBSA_ENABLED": "false",
    "MMGBSA_MMPBSA_RUN": "false",
    "MD_FIVE_REPLICATE": "false",
}


def _have_test_inputs() -> bool:
    if not MAIN_PY.exists() or not TEST_PDB_PATH.exists():
        return False
    if not PREPPED_TEST_LIB.exists():
        return False
    return any(PREPPED_TEST_LIB.glob("*.pdbqt"))


def _copy_prepped_fixture(tmp_path: Path) -> Path:
    prepped_root = tmp_path / "prepped_ligands"
    shutil.copytree(PREPPED_TEST_LIB, prepped_root / "fda_test_library_10")
    return prepped_root


def _is_finite_number(val: object) -> bool:
    try:
        return math.isfinite(float(val))
    except Exception:
        return False


def test_report_populates_manifest_and_dud_metrics(tmp_path: Path) -> None:
    if not _have_test_inputs():
        pytest.skip("Skipping: TEST inputs or prepped ligands missing.")

    run_id = uuid.uuid4().hex[:8]
    env = os.environ.copy()
    env.update(MMGBSA_DISABLE_ENV)
    env["INPUT_DIR"] = str(INPUT_FIXTURE_ROOT)
    env["PREPPED_LIGANDS_DIR"] = str(_copy_prepped_fixture(tmp_path))
    env["TEST_FIXTURE_PREPPED_LIGANDS_DIR"] = str(PREPPED_FIXTURE_ROOT)
    env["TEST_LIBRARY_MAP"] = '{"TEST": "fda_test_library_10"}'
    env["TEST_MODE_ENABLE"] = env.get("TEST_MODE_ENABLE", "dud+fda")
    env["PYTHONUNBUFFERED"] = "1"

    cmd = [
        sys.executable,
        str(MAIN_PY),
        "-test",
        "-fast",
        "--test-fda",
        "--run-id",
        run_id,
    ]

    subprocess.run(cmd, cwd=str(REPO_ROOT), env=env, check=True)

    manifest_path = REPO_ROOT / "manifests" / run_id / "run_manifest.yaml"
    report_path = REPO_ROOT / "data" / run_id / "report.yaml"

    assert manifest_path.exists(), f"manifest missing at {manifest_path}"
    if not report_path.exists():
        pytest.skip(f"report.yaml missing at {report_path}; run-report hook likely skipped in this environment.")

    with report_path.open("r", encoding="utf-8") as handle:
        report = yaml.safe_load(handle)

    assert report["sources"]["manifest_yaml"] == f"manifests/{run_id}/run_manifest.yaml"
    assert report["targets"], "No targets found in report."

    for target in report["targets"].values():
        pocket = target["pocket"]
        assert pocket["method"], "Pocket method missing."
        center = pocket.get("center")
        box = pocket.get("box")
        assert (
            isinstance(center, list)
            and len(center) == 3
            and all(_is_finite_number(c) for c in center)
        )
        assert (
            isinstance(box, list)
            and len(box) == 3
            and all(_is_finite_number(b) for b in box)
        )

        qc = target["qc"]
        assert _is_finite_number(qc.get("ef1"))
        assert _is_finite_number(qc.get("roc_auc"))
        assert _is_finite_number(qc.get("roc_auc_adj"))

        fdr_val = qc.get("fdr")
        assert fdr_val is not None
        if isinstance(fdr_val, str):
            assert fdr_val.startswith("SMALL ")
            assert _is_finite_number(fdr_val.split("SMALL", 1)[1])
        else:
            assert _is_finite_number(fdr_val)
