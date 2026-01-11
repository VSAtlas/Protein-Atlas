from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DUD_EVAL = REPO_ROOT / "analysis" / "dud_eval.py"


def _write_csv(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "ligand,score",
        "foo_active_1.pdbqt,-8.0",
        "foo_decoy_1.pdbqt,-6.0",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _make_env(tmp_path: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.update(
        {
            "OVERALL_DIR": str(tmp_path),
            "INPUT_DIR": str(tmp_path / "input_pdbs"),
            "OUTPUT_DIR": str(tmp_path / "processed_pdbs"),
            "PDBQT_DIR": str(tmp_path / "pdbqts"),
            "DOCKED_DIR": str(tmp_path / "docked"),
            "OUTPUT_LIGANDS_DIR": str(tmp_path / "prepped_ligands"),
            "LIGAND_EXTRACTED_DIR": str(tmp_path / "extracted_ligands"),
            "LIGANDS_MOL2_DIR": str(tmp_path / "ligands_mol2"),
            "CONFIGS_DIR": str(tmp_path / "configs"),
            "PYTHONUNBUFFERED": "1",
        }
    )
    for key in (
        "INPUT_DIR",
        "OUTPUT_DIR",
        "PDBQT_DIR",
        "DOCKED_DIR",
        "OUTPUT_LIGANDS_DIR",
        "LIGAND_EXTRACTED_DIR",
        "LIGANDS_MOL2_DIR",
        "CONFIGS_DIR",
    ):
        Path(env[key]).mkdir(parents=True, exist_ok=True)
    (Path(env["INPUT_DIR"]) / "1ABC.pdb").write_text(
        "HEADER 1ABC\nEND\n", encoding="utf-8"
    )
    return env


def _run_dud_eval(
    tmp_path: Path, env: dict[str, str], extra_args: list[str]
) -> subprocess.CompletedProcess:
    cmd = [
        sys.executable,
        str(DUD_EVAL),
        "--docked-root",
        str(Path(env["DOCKED_DIR"])),
        "--out-dir",
        str(tmp_path / "out"),
    ]
    cmd.extend(extra_args)
    return subprocess.run(
        cmd,
        env=env,
        cwd=str(tmp_path),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def test_dud_eval_scans_run_and_legacy_roots(tmp_path: Path) -> None:
    env = _make_env(tmp_path)
    run_id = "20250101_000000"

    run_root = Path(env["DOCKED_DIR"]) / run_id / "1ABC"
    legacy_root = Path(env["DOCKED_DIR"]) / "1ABC"
    _write_csv(run_root / "docking_score_long.csv")
    _write_csv(legacy_root / "docking_score_long.csv")

    # Run with explicit run-id: should scan docked/<run_id>/...
    result = _run_dud_eval(tmp_path, env, ["--run-id", run_id])
    assert result.returncode == 0, result.stdout
    analysis_root = Path(tmp_path / "out")
    assert (analysis_root / f"summary{run_id}.tsv").exists()
    target_dir = analysis_root / "1ABC"
    assert (target_dir / "metrics.tsv").exists()

    # Run without run-id: should scan run folders then legacy without crashing.
    result2 = _run_dud_eval(tmp_path, env, [])
    assert result2.returncode == 0, result2.stdout
    assert (analysis_root / "summary.tsv").exists()
