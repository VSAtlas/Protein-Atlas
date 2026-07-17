from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any


def add_analysis_subcommands(subparsers: argparse._SubParsersAction[Any]) -> None:
    dud = subparsers.add_parser("dud-eval", help="Run DUD/decoy benchmark evaluation.")
    dud.add_argument("run_id")
    dud.add_argument("--pdb-id", action="append", default=[])
    dud.add_argument("--valid", action="store_true")
    dud.add_argument("--score-col")
    dud.add_argument("--lig-col")
    dud.add_argument("--out-dir", default="analysis/dud_eval")
    dud.add_argument("--docked-root")
    dud.add_argument("--post-docked-root")
    dud.add_argument("--log-level", default="INFO")

    throughput = subparsers.add_parser(
        "throughput", help="Run throughput integrity or acceptance checks."
    )
    throughput_sub = throughput.add_subparsers(dest="throughput_cmd", required=True)
    integrity = throughput_sub.add_parser(
        "integrity", help="Validate prep/docking/post-docking throughput integrity."
    )
    integrity.add_argument("run_id")
    integrity.add_argument("--strict", dest="strict", action="store_true", default=True)
    integrity.add_argument("--no-strict", dest="strict", action="store_false")
    integrity.add_argument("--verbose", action="store_true")

    acceptance = throughput_sub.add_parser(
        "acceptance", help="Compare throughput KPIs for a baseline and candidate run."
    )
    acceptance.add_argument("baseline_run_id")
    acceptance.add_argument("candidate_run_id")
    acceptance.add_argument("--util-baseline-summary", default="")
    acceptance.add_argument("--util-candidate-summary", default="")
    acceptance.add_argument("--target-wall-improvement-pct", type=float, default=10.0)
    acceptance.add_argument("--target-idle-core-improvement-pct", type=float, default=15.0)
    acceptance.add_argument("--target-idle-frac-abs", type=float, default=0.03)
    acceptance.add_argument("--target-frag-dispatch-improvement-pct", type=float, default=10.0)
    acceptance.add_argument("--out-json", default="")

    interactions = subparsers.add_parser(
        "interactions", help="Export report/master rows to analysis interaction datasets."
    )
    interactions_sub = interactions.add_subparsers(dest="interactions_cmd", required=True)
    export = interactions_sub.add_parser(
        "export", help="Convert master_rows/heatmap input to partitioned Parquet."
    )
    export.add_argument("run_id")
    export.add_argument("--input")
    export.add_argument("--out-dir")
    export.add_argument("--partition-cols", default="run_id,target_id")
    export.add_argument("--targets")
    export.add_argument("--overwrite", action="store_true")
    export.add_argument("--verbose", action="store_true")

    holo = subparsers.add_parser(
        "holo-integrity",
        help=(
            "Audit frozen-SPD integrity from one manifest-bound artifact lineage."
        ),
    )
    holo.add_argument("run_id")
    holo.add_argument("--repo-root")
    holo.add_argument("--selection")
    holo.add_argument("--strict-manifest")
    holo.add_argument("--out-dir")
    holo.add_argument("--strict", dest="strict", action="store_true", default=True)
    holo.add_argument("--no-strict", dest="strict", action="store_false")

    dataset_eda = subparsers.add_parser(
        "dataset-eda", help="Run read-only EDA reports for an Atlas ML/tabular dataset."
    )
    dataset_eda.add_argument("--dataset")
    dataset_eda.add_argument("--dataset-name")
    dataset_eda.add_argument("--list-matches", action="store_true")
    dataset_eda.add_argument("--out-dir", default="outputs/data/dataset_eda")
    dataset_eda.add_argument("--label-col")
    dataset_eda.add_argument("--tool", action="append", default=[])
    dataset_eda.add_argument("--tools")
    dataset_eda.add_argument("--read-rows", type=int)
    dataset_eda.add_argument("--sample-rows", type=int, default=10_000)
    dataset_eda.add_argument("--random-state", type=int, default=13)
    dataset_eda.add_argument("--max-association-columns", type=int, default=80)
    dataset_eda.add_argument("--max-category-levels", type=int, default=200)
    dataset_eda.add_argument("--max-mi-features", type=int, default=300)
    dataset_eda.add_argument("--network-threshold", type=float, default=0.35)
    dataset_eda.add_argument("--network-max-edges", type=int, default=300)
    dataset_eda.add_argument("--full-ydata", action="store_true")
    dataset_eda.add_argument("--dython-nominal-assoc", default="cramer")
    dataset_eda.add_argument("--include-id-like", action="store_true")
    dataset_eda.add_argument("--include-high-cardinality", action="store_true")
    dataset_eda.add_argument("--fail-on-missing", action="store_true")


def run_analysis_subcommand(args: argparse.Namespace, repo_root: Path) -> int | None:
    command = getattr(args, "analysis_cmd", "")
    if command == "dud-eval":
        return _run_dud_eval(args, repo_root)
    if command == "throughput":
        return _run_throughput(args, repo_root)
    if command == "interactions":
        return _run_interactions(args, repo_root)
    if command == "holo-integrity":
        return _run_holo_integrity(args, repo_root)
    if command == "dataset-eda":
        return _run_dataset_eda(args, repo_root)
    return None


def _run_dud_eval(args: argparse.Namespace, repo_root: Path) -> int:
    forwarded = [
        "--run-id",
        str(args.run_id),
        "--out-dir",
        str(args.out_dir),
        "--log-level",
        str(args.log_level),
    ]
    if args.docked_root:
        forwarded.extend(["--docked-root", str(args.docked_root)])
    else:
        forwarded.extend(["--docked-root", str(repo_root / "outputs" / "docked")])
    if args.post_docked_root:
        forwarded.extend(["--post-docked-root", str(args.post_docked_root)])
    else:
        forwarded.extend(["--post-docked-root", str(repo_root / "outputs" / "post_docked")])
    for pdb_id in args.pdb_id or []:
        forwarded.extend(["--pdb-id", str(pdb_id)])
    if args.valid:
        forwarded.append("--valid")
    if args.score_col:
        forwarded.extend(["--score-col", str(args.score_col)])
    if args.lig_col:
        forwarded.extend(["--lig-col", str(args.lig_col)])
    return _run_noargv_main("atlas analysis dud-eval", forwarded, "analysis.dud_eval_core.orchestrate")


def _run_throughput(args: argparse.Namespace, repo_root: Path) -> int:
    if args.throughput_cmd == "integrity":
        from analysis.reporting import throughput_integrity

        forwarded = ["--run-id", str(args.run_id), "--repo-root", str(repo_root)]
        forwarded.append("--strict" if args.strict else "--no-strict")
        if args.verbose:
            forwarded.append("--verbose")
        return int(throughput_integrity.main(forwarded))

    if args.throughput_cmd == "acceptance":
        from analysis.reporting import throughput_acceptance

        forwarded = [
            "--repo-root",
            str(repo_root),
            "--baseline-run-id",
            str(args.baseline_run_id),
            "--candidate-run-id",
            str(args.candidate_run_id),
            "--target-wall-improvement-pct",
            str(args.target_wall_improvement_pct),
            "--target-idle-core-improvement-pct",
            str(args.target_idle_core_improvement_pct),
            "--target-idle-frac-abs",
            str(args.target_idle_frac_abs),
            "--target-frag-dispatch-improvement-pct",
            str(args.target_frag_dispatch_improvement_pct),
        ]
        if args.util_baseline_summary:
            forwarded.extend(["--util-baseline-summary", str(args.util_baseline_summary)])
        if args.util_candidate_summary:
            forwarded.extend(["--util-candidate-summary", str(args.util_candidate_summary)])
        if args.out_json:
            forwarded.extend(["--out-json", str(args.out_json)])
        return int(throughput_acceptance.main(forwarded))
    return 1


def _run_interactions(args: argparse.Namespace, repo_root: Path) -> int:
    if args.interactions_cmd != "export":
        return 1
    from analysis.cli import convert_interactions_to_parquet

    forwarded = ["--run-id", str(args.run_id), "--repo-root", str(repo_root)]
    if args.input:
        forwarded.extend(["--input", str(args.input)])
    if args.out_dir:
        forwarded.extend(["--out-dir", str(args.out_dir)])
    if args.partition_cols:
        forwarded.extend(["--partition-cols", str(args.partition_cols)])
    if args.targets:
        forwarded.extend(["--targets", str(args.targets)])
    if args.overwrite:
        forwarded.append("--overwrite")
    if args.verbose:
        forwarded.append("--verbose")
    return int(convert_interactions_to_parquet.main(forwarded))


def _run_holo_integrity(args: argparse.Namespace, repo_root: Path) -> int:
    from analysis.cli import spd_holo_integrity

    effective_root = Path(args.repo_root).expanduser() if args.repo_root else repo_root
    forwarded = [
        "--run-id",
        str(args.run_id),
        "--repo-root",
        str(effective_root),
        "--strict" if args.strict else "--no-strict",
    ]
    if args.selection:
        forwarded.extend(["--selection", str(args.selection)])
    if args.strict_manifest:
        forwarded.extend(["--strict-manifest", str(args.strict_manifest)])
    if args.out_dir:
        forwarded.extend(["--out-dir", str(args.out_dir)])
    return int(spd_holo_integrity.main(forwarded))


def _run_dataset_eda(args: argparse.Namespace, repo_root: Path) -> int:
    from analysis.cli import run_dataset_eda

    forwarded = ["--repo-root", str(repo_root), "--out-dir", str(args.out_dir)]
    if args.dataset:
        forwarded.extend(["--dataset", str(args.dataset)])
    if args.dataset_name:
        forwarded.extend(["--dataset-name", str(args.dataset_name)])
    if args.list_matches:
        forwarded.append("--list-matches")
    if args.label_col:
        forwarded.extend(["--label-col", str(args.label_col)])
    if args.tools:
        forwarded.extend(["--tools", str(args.tools)])
    for tool in args.tool or []:
        forwarded.extend(["--tool", str(tool)])
    if args.read_rows is not None:
        forwarded.extend(["--read-rows", str(args.read_rows)])
    forwarded.extend(
        [
            "--sample-rows",
            str(args.sample_rows),
            "--random-state",
            str(args.random_state),
            "--max-association-columns",
            str(args.max_association_columns),
            "--max-category-levels",
            str(args.max_category_levels),
            "--max-mi-features",
            str(args.max_mi_features),
            "--network-threshold",
            str(args.network_threshold),
            "--network-max-edges",
            str(args.network_max_edges),
            "--dython-nominal-assoc",
            str(args.dython_nominal_assoc),
        ]
    )
    if args.full_ydata:
        forwarded.append("--full-ydata")
    if args.include_id_like:
        forwarded.append("--include-id-like")
    if args.include_high_cardinality:
        forwarded.append("--include-high-cardinality")
    if args.fail_on_missing:
        forwarded.append("--fail-on-missing")
    return int(run_dataset_eda.main(forwarded))


def _run_noargv_main(prog: str, forwarded: list[str], module_name: str) -> int:
    old_argv = sys.argv[:]
    try:
        sys.argv = [prog, *forwarded]
        module = __import__(module_name, fromlist=["main"])
        rc = module.main()
        return int(rc or 0)
    finally:
        sys.argv = old_argv
