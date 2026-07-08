#!/usr/bin/env python3
"""
Reproduce DUD-eval discovery behavior when manifest data is missing.

Scenario A:
  - run-scoped docked outputs exist
  - manifest is absent
  - input_pdbs/<PDB>.pdb is absent
  Expectation: targets are discovered and evaluation runs.

Scenario B:
  - same as A, but input_pdbs/<PDB>.pdb exists
  Expectation: targets are discovered and evaluation runs.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def _write_docking_csv(path: Path, run_id: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "ligand_file,docking_score,run_id",
                f"actives_demo_0001.pdbqt,-9.1,{run_id}",
                f"decoys_demo_0001.pdbqt,-6.2,{run_id}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _prepare_fixture(root: Path, run_id: str, pdb_id: str) -> tuple[Path, Path]:
    docked_root = root / "docked"
    post_docked_root = root / "post_docked"

    dock_csv = (
        docked_root
        / run_id
        / pdb_id
        / "HOLO"
        / "pH7_0"
        / "docking_score_long.csv"
    )
    cons_csv = (
        docked_root
        / run_id
        / pdb_id
        / "HOLO"
        / "pH7_0"
        / "consensus_docking_scores.csv"
    )
    reranked_csv = (
        post_docked_root
        / run_id
        / pdb_id
        / "HOLO"
        / "pH7_0"
        / "consensus_reranked_scorch.csv"
    )

    _write_docking_csv(dock_csv, run_id)
    _write_docking_csv(cons_csv, run_id)
    _write_docking_csv(reranked_csv, run_id)
    return docked_root, post_docked_root


def _run_eval(
    *,
    repo_root: Path,
    fixture_root: Path,
    run_id: str,
    mode_tag: str,
) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    py_path = env.get("PYTHONPATH", "").strip()
    env["PYTHONPATH"] = (
        str(repo_root) if not py_path else f"{str(repo_root)}:{py_path}"
    )

    cmd = [
        sys.executable,
        "-m",
        "analysis.cli.dud_eval",
        "--docked-root",
        str(fixture_root / "docked"),
        "--post-docked-root",
        str(fixture_root / "post_docked"),
        "--run-id",
        run_id,
        "--out-dir",
        str(fixture_root / f"analysis_out_{mode_tag}"),
        "--log-level",
        "INFO",
    ]
    return subprocess.run(
        cmd,
        cwd=str(fixture_root),
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def _print_result(tag: str, cp: subprocess.CompletedProcess[str]) -> None:
    print(f"\n[{tag}] returncode={cp.returncode}")
    tail = (cp.stdout or "").splitlines()[-20:]
    if tail:
        print(f"[{tag}] stdout.tail:")
        for line in tail:
            print(line)
    err_tail = (cp.stderr or "").splitlines()[-20:]
    if err_tail:
        print(f"[{tag}] stderr.tail:")
        for line in err_tail:
            print(line)


def _reranked_summary_path(fixture_root: Path, mode_tag: str, run_id: str) -> Path:
    return (
        fixture_root
        / f"analysis_out_{mode_tag}"
        / run_id
        / "post_docked"
        / f"consensus_reranked_scorch_summary_{run_id}.tsv"
    )


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Reproduce dud_eval discovery behavior for missing-manifest runs."
    )
    ap.add_argument(
        "--repo-root",
        default=str(Path(__file__).resolve().parents[1]),
        help="Path to protein_automation repo root (for PYTHONPATH).",
    )
    ap.add_argument("--run-id", default="spr_repro_missing_manifest")
    ap.add_argument("--pdb-id", default="BNJS")
    ap.add_argument(
        "--work-root",
        default="/stor/home/mpg2352",
        help="Parent dir for temporary fixture (avoid /tmp on this environment).",
    )
    ap.add_argument(
        "--keep",
        action="store_true",
        help="Keep fixture directory instead of deleting it.",
    )
    args = ap.parse_args()

    repo_root = Path(args.repo_root).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    fixture = Path(
        tempfile.mkdtemp(prefix="dud_eval_repro_", dir=str(work_root))
    ).resolve()

    print(f"fixture={fixture}")
    try:
        _prepare_fixture(fixture, args.run_id, args.pdb_id)

        # Scenario A: no input_pdbs/<PDB>.pdb
        cp_a = _run_eval(
            repo_root=repo_root,
            fixture_root=fixture,
            run_id=args.run_id,
            mode_tag="missing_input_pdb",
        )
        _print_result("missing_input_pdb", cp_a)

        # Scenario B: add input_pdbs/<PDB>.pdb
        pdb_path = fixture / "input_pdbs" / f"{args.pdb_id}.pdb"
        pdb_path.parent.mkdir(parents=True, exist_ok=True)
        pdb_path.write_text(
            "\n".join(
                [
                    f"HEADER    {args.pdb_id} REPRO",
                    "ATOM      1  N   MET A   1      11.104  13.207   4.599  1.00 20.00           N",
                    "END",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

        cp_b = _run_eval(
            repo_root=repo_root,
            fixture_root=fixture,
            run_id=args.run_id,
            mode_tag="with_input_pdb",
        )
        _print_result("with_input_pdb", cp_b)

        # Expectations for reproducer:
        #   Both scenarios should proceed now, and post-docked reranked summary
        #   should be emitted in each mode.
        a_succeeded = cp_a.returncode == 0
        b_succeeded = cp_b.returncode == 0
        a_reranked_summary = _reranked_summary_path(
            fixture, "missing_input_pdb", args.run_id
        ).is_file()
        b_reranked_summary = _reranked_summary_path(
            fixture, "with_input_pdb", args.run_id
        ).is_file()

        print(
            f"\nsummary: scenarioA_succeeded={a_succeeded} scenarioB_succeeded={b_succeeded} "
            f"scenarioA_reranked_summary={a_reranked_summary} scenarioB_reranked_summary={b_reranked_summary}"
        )

        if not (a_succeeded and b_succeeded and a_reranked_summary and b_reranked_summary):
            print("repro_result=UNEXPECTED")
            return 1
        print("repro_result=EXPECTED")
        return 0
    finally:
        if args.keep:
            print(f"fixture_kept={fixture}")
        else:
            shutil.rmtree(fixture, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
