"""Small DUD-E water-policy enrichment benchmark manifest and staging helpers."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, cast

from path_router.path_router import _ensure_router_roots


@dataclass(frozen=True)
class WaterDudeTarget:
    pdb_id: str
    dude_target: str
    library: str
    conservation_refs: tuple[str, ...]
    water_rationale: str


WATER_DUDE_TARGETS: tuple[WaterDudeTarget, ...] = (
    WaterDudeTarget(
        pdb_id="3EML",
        dude_target="aa2ar",
        library="water_aa2ar",
        conservation_refs=("3PWH", "4EIY", "5IU4"),
        water_rationale="A2A receptor has a conserved ligand-site water network.",
    ),
    WaterDudeTarget(
        pdb_id="1UYG",
        dude_target="hs90a",
        library="water_hs90a",
        conservation_refs=("1UY6", "1YER", "1YC4"),
        water_rationale="HSP90 N-terminal ligands frequently retain or displace ordered pocket waters.",
    ),
    WaterDudeTarget(
        pdb_id="1XL2",
        dude_target="hivpr",
        library="water_hivpr",
        conservation_refs=("1HXW", "1HSG", "2Q3K"),
        water_rationale="HIV protease has a conserved flap/ligand bridging-water motif.",
    ),
)
WATER_POLICIES = ("remove_all", "site_only", "keep_all")


def write_plan(out_path: Path, *, max_ligands: int) -> dict[str, object]:
    """Write a reproducible plan for remove-all/site-only/keep-all water runs."""

    plan = _plan_payload(max_ligands=max_ligands)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(plan, indent=2, sort_keys=True), encoding="utf-8")
    return plan


def stage_libraries(
    *,
    source_prepped_root: Path,
    output_prepped_root: Path,
    max_ligands: int,
    max_actives: int,
) -> list[dict[str, object]]:
    """Copy small active/decoy PDBQT subsets into water benchmark libraries."""

    reports: list[dict[str, object]] = []
    for target in WATER_DUDE_TARGETS:
        src = source_prepped_root / target.dude_target
        dst = output_prepped_root / target.library
        reports.append(
            _stage_one_library(
                src,
                dst,
                max_ligands=max_ligands,
                max_actives=max_actives,
                target=target,
            )
        )
    return reports


def _plan_payload(*, max_ligands: int) -> dict[str, object]:
    refs = {
        target.pdb_id: list(target.conservation_refs)
        for target in WATER_DUDE_TARGETS
    }
    pdb_ids = ",".join(target.pdb_id for target in WATER_DUDE_TARGETS)
    return {
        "benchmark": "water_dude_mini",
        "max_ligands_per_target": int(max_ligands),
        "targets": [_target_payload(target) for target in WATER_DUDE_TARGETS],
        "water_policies": WATER_POLICIES,
        "conservation_refs_env": json.dumps(refs, sort_keys=True),
        "library_map": {
            target.pdb_id: target.library for target in WATER_DUDE_TARGETS
        },
        "commands": _run_commands(pdb_ids, refs),
        "evaluation": {
            "metrics": ["ROC_AUC", "PR_AUC", "logAUC", "BEDROC_alpha_20", "EF@1%", "EF@5%"],
            "compare": "Run analysis.cli.dud_eval for each run_id, then compare selected/site_only against remove_all and keep_all.",
        },
    }


def _run_commands(pdb_ids: str, refs: dict[str, list[str]]) -> list[str]:
    commands: list[str] = []
    refs_env = json.dumps(refs, sort_keys=True)
    for policy in WATER_POLICIES:
        run_id = f"water_dude_{policy}"
        commands.append(
            "ATLAS_WATER_BENCH_POLICY={policy} "
            "ATLAS_WATER_CONSERVATION_REFS='{refs}' "
            "python main.py --waterbench -fast --run-id {run_id} --pdbs {pdbs}".format(
                policy=policy,
                refs=refs_env,
                run_id=run_id,
                pdbs=pdb_ids,
            )
        )
        commands.append(
            "python -m analysis.cli.dud_eval --run-id {run_id} "
            "--docked-root docked --post-docked-root post_docked "
            "--out-dir analysis/dud_eval".format(run_id=run_id)
        )
    return commands


def run_policy_benchmark(
    *,
    run_prefix: str,
    policies: Sequence[str],
    out_report: Path,
    eval_out_dir: Path,
    skip_docking: bool,
    skip_eval: bool,
    cpu: int,
    max_ligands: int,
    fast: bool = True,
) -> dict[str, object]:
    """Run remove-all/site-only/keep-all water policies and summarize metrics."""

    repo_root = _default_repo_root()
    refs = {
        target.pdb_id: list(target.conservation_refs)
        for target in WATER_DUDE_TARGETS
    }
    pdb_ids = ",".join(target.pdb_id for target in WATER_DUDE_TARGETS)
    env_base = _waterbench_env(refs, cpu=cpu, max_ligands=max_ligands)
    runs: list[dict[str, object]] = []
    for policy in policies:
        if policy not in WATER_POLICIES:
            raise ValueError(f"unknown water policy: {policy}")
        run_id = f"{run_prefix}_{policy}"
        env = dict(env_base)
        env["ATLAS_WATER_BENCH_POLICY"] = policy
        dock_cmd = [
            sys.executable,
            "main.py",
            "--waterbench",
            "--run-id",
            run_id,
            "--pdbs",
            pdb_ids,
        ]
        if fast:
            dock_cmd.insert(3, "-fast")
        eval_cmd = [
            sys.executable,
            "-m",
            "analysis.cli.dud_eval",
            "--run-id",
            run_id,
            "--docked-root",
            str(repo_root / "outputs" / "docked"),
            "--post-docked-root",
            str(repo_root / "outputs" / "post_docked"),
            "--out-dir",
            str(eval_out_dir),
            "--prepped-root",
            str(repo_root / "prepped_ligands"),
            "--report-library",
        ]
        dock_rc = 0 if skip_docking else _run_logged(dock_cmd, repo_root, env)
        eval_rc = 0
        if not skip_eval and dock_rc == 0:
            eval_rc = _run_logged(eval_cmd, repo_root, env)
        runs.append(
            {
                "policy": policy,
                "run_id": run_id,
                "docking_returncode": dock_rc,
                "dud_eval_returncode": eval_rc,
                "docking_command": dock_cmd,
                "dud_eval_command": eval_cmd,
            }
        )
    summary = summarize_policy_metrics(
        run_ids=[str(row["run_id"]) for row in runs],
        eval_out_dir=eval_out_dir,
    )
    payload: dict[str, object] = {
        "benchmark": "water_dude_mini",
        "fast_mode": bool(fast),
        "max_ligands_per_target": int(max_ligands),
        "policies": list(policies),
        "targets": [_target_payload(target) for target in WATER_DUDE_TARGETS],
        "runs": runs,
        "summary": summary,
    }
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def summarize_policy_metrics(
    *,
    run_ids: Sequence[str],
    eval_out_dir: Path,
) -> list[dict[str, object]]:
    """Collect DUD-E metric TSV rows for water-policy comparison."""

    rows: list[dict[str, object]] = []
    for run_id in run_ids:
        policy = _policy_from_run_id(run_id)
        for metrics_path in sorted((eval_out_dir / run_id).rglob("metrics.tsv")):
            with metrics_path.open("r", encoding="utf-8", errors="ignore", newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    row["run_id"] = run_id
                    row["policy"] = policy
                    row["metrics_path"] = str(metrics_path)
                    rows.append(row)
    summary_tsv = eval_out_dir / "water_policy_summary.tsv"
    summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        fields = sorted({key for row in rows for key in row})
        with summary_tsv.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    return rows


def _waterbench_env(
    refs: dict[str, list[str]],
    *,
    cpu: int,
    max_ligands: int,
) -> dict[str, str]:
    env = dict(os.environ)
    env["ATLAS_WATER_CONSERVATION_REFS"] = json.dumps(refs, sort_keys=True)
    env["ATLAS_WATER_BENCH_MAX_LIGANDS"] = str(max(1, int(max_ligands)))
    env["CPU"] = str(max(1, min(int(cpu), 32)))
    env["QUIET_CONSOLE_OVERRIDE"] = "1"
    env["FILTER_VINA_STDOUT"] = "1"
    tmp_dir = Path.home() / ".atlas_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    env["TMPDIR"] = str(tmp_dir)
    env_bin = str(Path(sys.executable).resolve().parent)
    env["PATH"] = os.pathsep.join([env_bin, env.get("PATH", "")])
    return env


def _run_logged(cmd: Sequence[str], cwd: Path, env: dict[str, str]) -> int:
    print("[water_dude.run]", " ".join(cmd), flush=True)
    return subprocess.run(list(cmd), cwd=cwd, env=env, check=False).returncode


def _policy_from_run_id(run_id: str) -> str:
    for policy in WATER_POLICIES:
        if run_id.endswith(f"_{policy}") or f"_{policy}_" in run_id:
            return policy
    return ""


def _target_payload(target: WaterDudeTarget) -> dict[str, object]:
    return {
        "pdb_id": target.pdb_id,
        "dude_target": target.dude_target,
        "library": target.library,
        "conservation_refs": list(target.conservation_refs),
        "water_rationale": target.water_rationale,
    }


def _stage_one_library(
    src: Path,
    dst: Path,
    *,
    max_ligands: int,
    max_actives: int,
    target: WaterDudeTarget,
) -> dict[str, object]:
    actives = _sorted_pdbqts(src, "actives")
    decoys = _sorted_pdbqts(src, "decoys")
    if not actives or not decoys:
        return _stage_report(target, "missing_source_pdbqt", src, dst, 0, 0)
    take_actives = min(max_actives, len(actives), max(1, max_ligands // 20))
    take_decoys = min(len(decoys), max(0, max_ligands - take_actives))
    _copy_subset(actives[:take_actives], dst / "actives")
    _copy_subset(decoys[:take_decoys], dst / "decoys")
    _write_manifest(dst, target, take_actives, take_decoys, max_ligands)
    return _stage_report(target, "staged", src, dst, take_actives, take_decoys)


def _sorted_pdbqts(root: Path, label: str) -> list[Path]:
    direct = sorted((root / label).glob("*.pdbqt"))
    if direct:
        return direct
    prefix = "actives_final" if label == "actives" else "decoys_final"
    return sorted(root.glob(f"{prefix}*.pdbqt"))


def _copy_subset(paths: Sequence[Path], dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    for path in paths:
        shutil.copy2(path, dst / path.name)


def _write_manifest(
    dst: Path,
    target: WaterDudeTarget,
    actives: int,
    decoys: int,
    max_ligands: int,
) -> None:
    payload = {
        "target": _target_payload(target),
        "actives": int(actives),
        "decoys": int(decoys),
        "total": int(actives + decoys),
        "max_ligands": int(max_ligands),
    }
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "_water_benchmark_manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def _stage_report(
    target: WaterDudeTarget,
    status: str,
    src: Path,
    dst: Path,
    actives: int,
    decoys: int,
) -> dict[str, object]:
    return {
        "pdb_id": target.pdb_id,
        "dude_target": target.dude_target,
        "library": target.library,
        "status": status,
        "source": str(src),
        "destination": str(dst),
        "actives": int(actives),
        "decoys": int(decoys),
        "total": int(actives + decoys),
    }


def _default_repo_root() -> Path:
    return _ensure_router_roots().overall


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    plan = sub.add_parser("plan", help="Write the water benchmark run plan JSON.")
    plan.add_argument("--out", type=Path, default=Path("docs/water_dude_benchmark_plan.json"))
    plan.add_argument("--max-ligands", type=int, default=999)
    stage = sub.add_parser("stage-libraries", help="Copy <=999 ligand PDBQT subsets.")
    stage.add_argument("--source-prepped-root", type=Path, default=None)
    stage.add_argument("--output-prepped-root", type=Path, default=None)
    stage.add_argument("--max-ligands", type=int, default=999)
    stage.add_argument("--max-actives", type=int, default=25)
    stage.add_argument("--report", type=Path, default=Path("docs/water_dude_stage_report.json"))
    run = sub.add_parser("run-policies", help="Run and summarize water-policy benchmarks.")
    run.add_argument("--run-prefix", default="water_dude")
    run.add_argument("--policies", default=",".join(WATER_POLICIES))
    run.add_argument("--report", type=Path, default=Path("docs/water_dude_run_report.json"))
    run.add_argument("--eval-out-dir", type=Path, default=Path("outputs/data/water_dude_eval"))
    run.add_argument("--skip-docking", action="store_true")
    run.add_argument("--skip-eval", action="store_true")
    run.add_argument("--cpu", type=int, default=16)
    run.add_argument("--max-ligands", type=int, default=96)
    run.add_argument(
        "--no-fast",
        action="store_true",
        help="Run main.py --waterbench without -fast for higher-search docking.",
    )
    summarize = sub.add_parser("summarize", help="Summarize existing water-policy metrics.")
    summarize.add_argument("--run-ids", required=True)
    summarize.add_argument("--eval-out-dir", type=Path, default=Path("outputs/data/water_dude_eval"))
    summarize.add_argument("--report", type=Path, default=Path("docs/water_dude_summary.json"))
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    repo_root = _default_repo_root()
    if args.cmd == "plan":
        write_plan(args.out, max_ligands=args.max_ligands)
        return 0
    if args.cmd == "stage-libraries":
        source = args.source_prepped_root or repo_root / "prepped_ligands"
        output = args.output_prepped_root or repo_root / "prepped_ligands"
        reports = stage_libraries(
            source_prepped_root=source,
            output_prepped_root=output,
            max_ligands=int(args.max_ligands),
            max_actives=int(args.max_actives),
        )
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(reports, indent=2, sort_keys=True), encoding="utf-8")
        return 0
    if args.cmd == "run-policies":
        policies = [token.strip() for token in args.policies.split(",") if token.strip()]
        payload = run_policy_benchmark(
            run_prefix=str(args.run_prefix),
            policies=policies,
            out_report=args.report,
            eval_out_dir=args.eval_out_dir,
            skip_docking=bool(args.skip_docking),
            skip_eval=bool(args.skip_eval),
            cpu=int(args.cpu),
            max_ligands=int(args.max_ligands),
            fast=not bool(args.no_fast),
        )
        run_rows = cast(list[dict[str, object]], payload["runs"])
        failed = [
            row
            for row in run_rows
            if row["docking_returncode"] != 0 or row["dud_eval_returncode"] != 0
        ]
        return 1 if failed else 0
    if args.cmd == "summarize":
        rows = summarize_policy_metrics(
            run_ids=[token.strip() for token in args.run_ids.split(",") if token.strip()],
            eval_out_dir=args.eval_out_dir,
        )
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(rows, indent=2, sort_keys=True), encoding="utf-8")
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
