from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
from typing import Optional

from config.runtime_config import load_config
from post_docking.mmgbsa.mmgbsa_batch import (
    extract_artifact,
    group_hits,
    load_report_hits,
    run_jobs,
    write_plan,
)
from post_docking.mmgbsa.mmgbsa_resp import finalize_resp_workflow, plan_resp_workflow
from post_docking.mmgbsa.mmgbsa_staged import (
    run_staged_aggregate_job,
    run_staged_md_job,
    stage_from_report,
)
from post_docking.mmgbsa.mmgbsa_wang_lite import (
    build_wang_lite_combined_report,
    build_wang_lite_report,
    write_wang_crystal_pose_report,
)
from post_docking.mmgbsa.mmgbsa_wang2016_benchmark import build_wang2016_benchmark
from post_docking.mmgbsa.mmgbsa_wang2016_production import (
    write_wang2016_production_package,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="atlas-mmgbsa",
        description="Run MM/GBSA as a report-driven post-processing workflow.",
    )
    subparsers = parser.add_subparsers(dest="command")
    batch = subparsers.add_parser("batch-from-report")
    batch.add_argument("--report", required=True, type=Path)
    batch.add_argument("--run-id", required=True)
    batch.add_argument("--config", default="config.txt")
    batch.add_argument("--default-stage-dir", default="stage1")
    batch.add_argument("--default-variant", default="")
    batch.add_argument("--default-ph-label", default="")
    batch.add_argument("--require-significant", action="store_true")
    batch.add_argument("--max-q-value", type=float, default=None)
    batch.add_argument("--score-column", default="")
    batch.add_argument("--score-direction", default="lower")
    batch.add_argument("--top-n", type=int, default=0)
    batch.add_argument("--artifact", action="append", default=[], type=Path)
    batch.add_argument("--artifact-dest", type=Path, default=None)
    batch.add_argument("--overall-dir", type=Path, default=None)
    batch.add_argument("--post-docked-dir", type=Path, default=None)
    batch.add_argument("--docked-dir", type=Path, default=None)
    batch.add_argument("--configs-dir", type=Path, default=None)
    batch.add_argument("--plan-out", type=Path, default=None)
    batch.add_argument("--test-mode", default="off")
    batch.add_argument("--dry-run", action="store_true")
    batch.add_argument("--log-level", default="INFO")

    staged = subparsers.add_parser("staged-from-report")
    staged.add_argument("--report", required=True, type=Path)
    staged.add_argument("--run-id", required=True)
    staged.add_argument("--out-dir", type=Path, default=None)
    staged.add_argument("--config", default="config.txt")
    staged.add_argument(
        "--validation-tier",
        choices=("custom", "contract", "topology", "tiny-md", "canary", "wang-lite", "production"),
        default="custom",
    )
    staged.add_argument("--fixture-out-dir", type=Path, default=None)
    staged.add_argument("--default-stage-dir", default="stage1")
    staged.add_argument("--default-variant", default="")
    staged.add_argument("--default-ph-label", default="")
    staged.add_argument("--require-significant", action="store_true")
    staged.add_argument("--max-q-value", type=float, default=None)
    staged.add_argument("--score-column", default="")
    staged.add_argument("--score-direction", default="lower")
    staged.add_argument("--top-n", type=int, default=0)
    staged.add_argument("--artifact", action="append", default=[], type=Path)
    staged.add_argument("--artifact-dest", type=Path, default=None)
    staged.add_argument("--overall-dir", type=Path, default=None)
    staged.add_argument("--post-docked-dir", type=Path, default=None)
    staged.add_argument("--docked-dir", type=Path, default=None)
    staged.add_argument("--configs-dir", type=Path, default=None)
    staged.add_argument("--no-prepare", action="store_true")
    staged.add_argument("--dry-run", action="store_true")
    staged.add_argument("--test-mode", default="off")
    staged.add_argument("--replicates", type=int, default=1)
    staged.add_argument("--prod-ps", type=float, default=4000.0)
    staged.add_argument("--frame-stride-ps", type=float, default=10.0)
    staged.add_argument("--md-engine", choices=("sander", "openmm"), default="openmm")
    staged.add_argument("--gpus", default="0")
    staged.add_argument("--cpus", type=int, default=32)
    staged.add_argument("--md-ranks", type=int, default=1)
    staged.add_argument("--mmpbsa-ranks", type=int, default=32)
    staged.add_argument("--openmm-start-stage", default="heat")
    staged.add_argument("--python-exe", default="")
    staged.add_argument("--log-level", default="INFO")

    staged_md = subparsers.add_parser("staged-run-md")
    staged_md.add_argument("--job", required=True, type=Path)
    staged_md.add_argument("--config", default="config.txt")
    staged_md.add_argument("--force", action="store_true")
    staged_md.add_argument("--log-level", default="INFO")

    staged_agg = subparsers.add_parser("staged-aggregate")
    staged_agg.add_argument("--job", required=True, type=Path)
    staged_agg.add_argument("--config", default="config.txt")
    staged_agg.add_argument("--force", action="store_true")
    staged_agg.add_argument("--log-level", default="INFO")

    wang_lite = subparsers.add_parser("wang-lite-report")
    wang_lite.add_argument("--staged-dir", required=True, type=Path)
    wang_lite.add_argument("--out-dir", required=True, type=Path)
    wang_lite.add_argument("--log-level", default="INFO")

    wang_lite_combined = subparsers.add_parser("wang-lite-combined-report")
    wang_lite_combined.add_argument("--staged-dir", required=True, action="append", type=Path)
    wang_lite_combined.add_argument("--out-dir", required=True, type=Path)
    wang_lite_combined.add_argument("--log-level", default="INFO")

    wang_crystal = subparsers.add_parser("wang2016-crystal-report")
    wang_crystal.add_argument("--report", required=True, type=Path)
    wang_crystal.add_argument("--out-csv", required=True, type=Path)
    wang_crystal.add_argument("--run-id", default="")
    wang_crystal.add_argument("--log-level", default="INFO")

    resp_plan = subparsers.add_parser("resp-plan")
    resp_plan.add_argument("--ligand", required=True, type=Path)
    resp_plan.add_argument("--out-dir", required=True, type=Path)
    resp_plan.add_argument("--net-charge", required=True, type=int)
    resp_plan.add_argument("--config", default="config.txt")
    resp_plan.add_argument("--residue-name", default="LIG")
    resp_plan.add_argument("--atom-type", default="gaff2")
    resp_plan.add_argument("--qm-engine", default="gaussian")
    resp_plan.add_argument("--qm-route", default="")
    resp_plan.add_argument("--run-antechamber", action="store_true")
    resp_plan.add_argument("--force", action="store_true")
    resp_plan.add_argument("--log-level", default="INFO")

    resp_finalize = subparsers.add_parser("resp-finalize")
    resp_finalize.add_argument("--qm-output", required=True, type=Path)
    resp_finalize.add_argument("--out-dir", required=True, type=Path)
    resp_finalize.add_argument("--net-charge", required=True, type=int)
    resp_finalize.add_argument("--config", default="config.txt")
    resp_finalize.add_argument("--residue-name", default="LIG")
    resp_finalize.add_argument("--atom-type", default="gaff2")
    resp_finalize.add_argument("--qm-engine", default="gaussian")
    resp_finalize.add_argument("--force", action="store_true")
    resp_finalize.add_argument("--log-level", default="INFO")

    bench = subparsers.add_parser("benchmark-wang2016")
    bench.add_argument("--out-dir", required=True, type=Path)
    bench.add_argument("--atlas-results", type=Path, default=None)
    bench.add_argument("--reference-csv", type=Path, default=None)
    bench.add_argument("--atlas-energy-column", default="delta_total")
    bench.add_argument("--reference-energy-column", default="delta_g_binding")
    bench.add_argument("--supplement-url", default="")
    bench.add_argument("--no-download", action="store_true")
    bench.add_argument("--download-pdbs", action="store_true")
    bench.add_argument("--pdbbind-root", type=Path, default=None)
    bench.add_argument("--pdbbind-zip", type=Path, default=None)
    bench.add_argument("--log-level", default="INFO")

    wang_prod = subparsers.add_parser("wang2016-production-package")
    wang_prod.add_argument("--report", required=True, type=Path)
    wang_prod.add_argument("--run-id", required=True)
    wang_prod.add_argument("--out-dir", required=True, type=Path)
    wang_prod.add_argument("--artifact", action="append", default=[], type=Path)
    wang_prod.add_argument("--artifact-dest", type=Path, default=None)
    wang_prod.add_argument("--overall-dir", type=Path, default=None)
    wang_prod.add_argument("--post-docked-dir", type=Path, default=None)
    wang_prod.add_argument("--docked-dir", type=Path, default=None)
    wang_prod.add_argument("--configs-dir", type=Path, default=None)
    wang_prod.add_argument("--cpus", type=int, default=32)
    wang_prod.add_argument("--partition", default="")
    wang_prod.add_argument("--hours", type=int, default=72)
    wang_prod.add_argument("--python-exe", default="")
    wang_prod.add_argument("--md-prod-ps", type=float, default=10000.0)
    wang_prod.add_argument("--md-replicates", type=int, default=3)
    wang_prod.add_argument("--frame-stride-ps", type=float, default=10.0)
    wang_prod.add_argument("--analysis-start-ps", type=float, default=0.0)
    wang_prod.add_argument("--analysis-end-ps", type=float, default=0.0)
    wang_prod.add_argument("--analysis-interval", type=int, default=1)
    wang_prod.add_argument("--md-engine", choices=("sander", "openmm"), default="sander")
    wang_prod.add_argument("--gpus", default="")
    wang_prod.add_argument("--openmm-start-stage", default="prod")
    wang_prod.add_argument(
        "--validation-mode",
        action="store_true",
        help="Allow shortened Wang smoke/validation runs without production replicate gates.",
    )
    wang_prod.add_argument("--log-level", default="INFO")
    return parser


def _default_staged_out_dir(run_id: str, validation_tier: str) -> Path:
    tier = str(validation_tier or "custom").strip().replace("-", "_") or "custom"
    return Path("outputs") / "mmgbsa_benchmarks" / run_id / f"staged_{tier}"


def _print_result(result: dict[str, object]) -> None:
    print(json.dumps(result, indent=2, sort_keys=True))


def _apply_path_overrides(cfg: dict[str, object], args: argparse.Namespace) -> None:
    if args.overall_dir is not None:
        cfg["OVERALL_DIR"] = str(args.overall_dir.expanduser())
    if args.post_docked_dir is not None:
        cfg["POST_DOCKED_DIR"] = str(args.post_docked_dir.expanduser())
    if args.docked_dir is not None:
        cfg["DOCKED_DIR"] = str(args.docked_dir.expanduser())
    if args.configs_dir is not None:
        cfg["CONFIGS_DIR"] = str(args.configs_dir.expanduser())
        os.environ["CONFIGS_DIR"] = str(cfg["CONFIGS_DIR"])


def _run_staged_from_report(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)
    _apply_path_overrides(cfg, args)
    result = stage_from_report(
        report=args.report,
        run_id=args.run_id,
        out_dir=args.out_dir or _default_staged_out_dir(args.run_id, args.validation_tier),
        cfg=cfg,
        config_path=args.config,
        artifacts=args.artifact,
        artifact_dest=args.artifact_dest,
        default_stage_dir=args.default_stage_dir,
        default_variant=args.default_variant,
        default_ph_label=args.default_ph_label,
        require_significant=args.require_significant,
        max_q_value=args.max_q_value,
        score_column=args.score_column,
        score_direction=args.score_direction,
        top_n=args.top_n,
        run_prepare=not args.no_prepare,
        test_mode=args.test_mode,
        dry_run=args.dry_run,
        replicates=args.replicates,
        prod_ps=args.prod_ps,
        frame_stride_ps=args.frame_stride_ps,
        md_engine=args.md_engine,
        gpus=args.gpus,
        cpus=args.cpus,
        md_ranks=args.md_ranks,
        mmpbsa_ranks=args.mmpbsa_ranks,
        openmm_start_stage=args.openmm_start_stage,
        python_exe=args.python_exe,
        validation_tier=args.validation_tier,
        fixture_out_dir=args.fixture_out_dir,
    )
    _print_result(result)
    return 0


def _run_staged_md(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)
    result = run_staged_md_job(job_json=args.job, cfg=cfg, force=args.force)
    _print_result(result)
    return 0 if result.get("ok") else 1


def _run_staged_aggregate(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)
    result = run_staged_aggregate_job(job_json=args.job, cfg=cfg, force=args.force)
    _print_result(result)
    return 0 if result.get("ok") else 1


def _run_resp_plan(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)
    result = plan_resp_workflow(
        ligand_path=args.ligand,
        out_dir=args.out_dir,
        net_charge=args.net_charge,
        cfg=cfg,
        residue_name=args.residue_name,
        atom_type=args.atom_type,
        qm_engine=args.qm_engine,
        qm_route=args.qm_route or None,
        run_antechamber=args.run_antechamber,
        force=args.force,
    )
    _print_result(result)
    return 0


def _run_resp_finalize(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)
    result = finalize_resp_workflow(
        qm_output=args.qm_output,
        out_dir=args.out_dir,
        net_charge=args.net_charge,
        cfg=cfg,
        residue_name=args.residue_name,
        atom_type=args.atom_type,
        qm_engine=args.qm_engine,
        force=args.force,
    )
    _print_result(result)
    return 0


def _run_benchmark_wang2016(args: argparse.Namespace) -> int:
    result = build_wang2016_benchmark(
        out_dir=args.out_dir,
        atlas_results=args.atlas_results,
        reference_csv=args.reference_csv,
        atlas_energy_column=args.atlas_energy_column,
        reference_energy_column=args.reference_energy_column,
        download=not args.no_download,
        download_pdbs=args.download_pdbs,
        pdbbind_root=args.pdbbind_root,
        pdbbind_zip=args.pdbbind_zip,
        supplement_url=args.supplement_url or None,
    )
    _print_result(result)
    return 0


def _run_wang2016_production_package(args: argparse.Namespace) -> int:
    result = write_wang2016_production_package(
        report=args.report,
        run_id=args.run_id,
        out_dir=args.out_dir,
        artifacts=args.artifact,
        artifact_dest=args.artifact_dest,
        overall_dir=args.overall_dir,
        post_docked_dir=args.post_docked_dir,
        docked_dir=args.docked_dir,
        configs_dir=args.configs_dir,
        cpus=args.cpus,
        partition=args.partition,
        hours=args.hours,
        python_exe=args.python_exe or None,
        md_prod_ps=args.md_prod_ps,
        md_replicates=args.md_replicates,
        frame_stride_ps=args.frame_stride_ps,
        analysis_start_ps=args.analysis_start_ps,
        analysis_end_ps=args.analysis_end_ps,
        analysis_interval=args.analysis_interval,
        production_strict=not args.validation_mode,
        md_engine=args.md_engine,
        gpus=args.gpus,
        openmm_start_stage=args.openmm_start_stage,
    )
    _print_result(result)
    return 0


def _run_wang_lite_report(args: argparse.Namespace) -> int:
    result = build_wang_lite_report(staged_dir=args.staged_dir, out_dir=args.out_dir)
    _print_result(result)
    return 0 if result.get("ok") else 1


def _run_wang_lite_combined_report(args: argparse.Namespace) -> int:
    result = build_wang_lite_combined_report(
        staged_dirs=args.staged_dir,
        out_dir=args.out_dir,
    )
    _print_result(result)
    return 0 if result.get("ok") else 1


def _run_wang_crystal_report(args: argparse.Namespace) -> int:
    result = write_wang_crystal_pose_report(
        report=args.report,
        out_csv=args.out_csv,
        run_id=args.run_id,
    )
    _print_result(result)
    return 0 if result.get("ok") else 1


def _run_batch_from_report(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)
    _apply_path_overrides(cfg, args)

    artifact_dest = args.artifact_dest
    if artifact_dest is None:
        artifact_dest = Path(str(cfg.get("POST_DOCKED_DIR", "post_docked"))) / args.run_id
    for archive in args.artifact:
        extract_artifact(archive, artifact_dest)

    hits = load_report_hits(
        args.report,
        run_id=args.run_id,
        default_stage_dir=args.default_stage_dir,
        default_variant=args.default_variant,
        default_ph_label=args.default_ph_label,
        require_significant=args.require_significant,
        max_q_value=args.max_q_value,
        score_column=args.score_column,
        score_direction=args.score_direction,
        top_n=max(0, args.top_n),
    )
    jobs = group_hits(hits)
    results = run_jobs(jobs, cfg=cfg, test_mode=args.test_mode, dry_run=args.dry_run)
    if args.plan_out:
        write_plan(args.plan_out, results)
    _print_result(results)
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 2

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    handlers = {
        "batch-from-report": _run_batch_from_report,
        "benchmark-wang2016": _run_benchmark_wang2016,
        "resp-finalize": _run_resp_finalize,
        "resp-plan": _run_resp_plan,
        "staged-aggregate": _run_staged_aggregate,
        "staged-from-report": _run_staged_from_report,
        "staged-run-md": _run_staged_md,
        "wang-lite-combined-report": _run_wang_lite_combined_report,
        "wang-lite-report": _run_wang_lite_report,
        "wang2016-crystal-report": _run_wang_crystal_report,
        "wang2016-production-package": _run_wang2016_production_package,
    }
    handler = handlers.get(args.command)
    if handler is None:
        parser.print_help()
        return 2
    return handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
