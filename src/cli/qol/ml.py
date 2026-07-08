from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
from typing import Any, Sequence

from config.output_paths import output_root

from cli.qol import _bindings


def _default_run_dir(repo_root: Path, run_id: str) -> Path:
    return output_root(repo_root, "data") / str(run_id)


def _default_ml_out_dir(repo_root: Path, run_id: str, name: str) -> Path:
    return _default_run_dir(repo_root, run_id) / "ml" / name


def _run_module_main(module_name: str, forwarded: list[str]) -> int:
    module = importlib.import_module(module_name)
    return int(module.main(forwarded) or 0)


def _cmd_ml(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas ml",
        description="Discoverable wrappers for Atlas ML training, audits, and source/PU sensitivity workflows.",
    )
    sub = parser.add_subparsers(dest="ml_cmd", required=True)

    train = sub.add_parser("train", help="Run the one-command Atlas ML training pass.")
    train.add_argument("--run-id", default="pilotstudy")
    train.add_argument("--run-dir")
    train.add_argument("--out-dir")
    train.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    train.add_argument("--sample-rows", type=int, default=0)
    train.add_argument("--sample-seed", type=int, default=42)
    train.add_argument("--sample-strategy", choices=["stratified", "random"], default="stratified")
    train.add_argument("--split-manifest")
    train.add_argument("--expert-split-manifest", nargs="*", default=None)
    train.add_argument("--no-mlflow", action="store_true")
    train.add_argument("--no-feature-cache", action="store_true")
    train.add_argument("--refresh-feature-cache", action="store_true")

    audit = sub.add_parser("audit", help="Run the reusable ML audit suite.")
    audit.add_argument("--run-id", default="pilotstudy")
    audit.add_argument("--dataset", required=True)
    audit.add_argument("--label", required=True)
    audit.add_argument("--feature-set")
    audit.add_argument("--out-dir")
    audit.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")

    source_pu = sub.add_parser("source-pu", help="Run source-transfer and PU sensitivity suites.")
    source_pu.add_argument("--run-id", default="pilotstudy")
    source_pu.add_argument("--run-dir")
    source_pu.add_argument("--out-dir")
    source_pu.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")

    hpo = sub.add_parser("hpo", help="Run an optional Optuna-backed hyperparameter sweep.")
    hpo.add_argument("--run-id", default="pilotstudy")
    hpo.add_argument("--dataset", required=True)
    hpo.add_argument("--label", required=True)
    hpo.add_argument("--out-dir")
    hpo.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    hpo.add_argument("--pruner", choices=["median", "successive_halving", "hyperband", "none"], default="median")
    hpo.add_argument("--split-manifest")
    hpo.add_argument("--no-mlflow", action="store_true")

    chemprop = sub.add_parser("chemprop", help="Stage or run an external Chemprop molecular baseline.")
    chemprop.add_argument("--run-id", default="pilotstudy")
    chemprop.add_argument("--dataset", required=True)
    chemprop.add_argument("--label", required=True)
    chemprop.add_argument("--out-dir")

    tdc = sub.add_parser("tdc", help="Stage a PyTDC or local TDC-style benchmark dataset.")
    tdc.add_argument("--run-id", default="pilotstudy")
    tdc_source = tdc.add_mutually_exclusive_group(required=True)
    tdc_source.add_argument("--input-csv")
    tdc_source.add_argument("--name")
    tdc.add_argument("--out")
    tdc.add_argument("--dataset-name")

    reinvent = sub.add_parser("reinvent", help="Capture or launch REINVENT4 generated molecule outputs.")
    reinvent.add_argument("--run-id", default="pilotstudy")
    reinvent.add_argument("--generated-smiles")
    reinvent.add_argument("--out-dir")

    gen_bench = sub.add_parser("gen-bench", help="Prepare/capture GuacaMol or MOSES generated-SMILES benchmarks.")
    gen_bench.add_argument("--run-id", default="pilotstudy")
    gen_bench.add_argument("--suite", required=True, choices=["guacamol", "moses"])
    gen_bench.add_argument("--generated-smiles", required=True)
    gen_bench.add_argument("--out-dir")

    doctor = sub.add_parser("doctor", help="Preflight ML inputs and metadata before training.")
    doctor.add_argument("--run-id", default="pilotstudy")
    doctor.add_argument("--run-dir")
    doctor.add_argument("--dataset")
    doctor.add_argument("--label")
    doctor.add_argument("--feature-set")
    doctor.add_argument("--strict", action="store_true")
    doctor.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    doctor.add_argument("--write-report", action="store_true")
    doctor.add_argument("--json", action="store_true", dest="json_out")


    prepare = sub.add_parser("prepare-run", help="Validate an imported run and build fast ML iteration artifacts.")
    prepare.add_argument("--run-id", default="pilotstudy")
    prepare.add_argument("--run-dir")
    prepare.add_argument("--ml-root")
    prepare.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    prepare.add_argument("--sample-rows", type=int, default=5000)
    prepare.add_argument("--sample-seed", type=int, default=42)
    prepare.add_argument("--skip-splits", action="store_true")
    prepare.add_argument("--skip-sample-train", action="store_true")

    ready = sub.add_parser("ready", help="Alias for prepare-run: doctor, split lock, sample train, leaderboard.")
    ready.add_argument("--run-id", default="pilotstudy")
    ready.add_argument("--run-dir")
    ready.add_argument("--ml-root")
    ready.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    ready.add_argument("--sample-rows", type=int, default=5000)
    ready.add_argument("--sample-seed", type=int, default=42)
    ready.add_argument("--skip-splits", action="store_true")
    ready.add_argument("--skip-sample-train", action="store_true")

    leaderboard = sub.add_parser("leaderboard", help="Build a run-scoped ML leaderboard from manifests.")
    leaderboard.add_argument("--run-id", default="pilotstudy")
    leaderboard.add_argument("--ml-root")
    leaderboard.add_argument("--out-dir")

    splits = sub.add_parser("splits", help="Lock or validate run-scoped ML split manifests.")
    splits_sub = splits.add_subparsers(dest="splits_cmd", required=True)
    splits_lock = splits_sub.add_parser("lock")
    splits_lock.add_argument("--run-id", default="pilotstudy")
    splits_lock.add_argument("--run-dir")
    splits_lock.add_argument("--out-dir")
    splits_lock.add_argument("--experts", nargs="*", default=None)
    splits_lock.add_argument("--splits", nargs="*", default=None)
    splits_lock.add_argument("--seed", type=int, default=42)
    splits_lock.add_argument("--validation-fraction", type=float, default=0.15)
    splits_validate = splits_sub.add_parser("validate")
    splits_validate.add_argument("--run-id", default="pilotstudy")
    splits_validate.add_argument("--run-dir")
    splits_validate.add_argument("--split-root")

    external_eval = sub.add_parser("external-eval", help="Evaluate a trained Atlas model on an external labeled table.")
    external_eval.add_argument("--run-id", default="pilotstudy")
    external_eval.add_argument("--model-dir", required=True)
    external_eval.add_argument("--dataset", required=True)
    external_eval.add_argument("--label")
    external_eval.add_argument("--name", required=True)
    external_eval.add_argument("--out-dir")

    status = sub.add_parser("status", help="Summarize run-scoped ML artifacts.")
    status.add_argument("--run-id", default="pilotstudy")
    status.add_argument("--out-dir")
    status.add_argument("--json", action="store_true", dest="json_out")

    args, extra = parser.parse_known_args(list(argv))
    repo_root = _bindings.repo_root()
    if args.ml_cmd == "train":
        return _cmd_ml_train(args, extra, repo_root)
    if args.ml_cmd == "audit":
        return _cmd_ml_audit(args, extra, repo_root)
    if args.ml_cmd == "source-pu":
        return _cmd_ml_source_pu(args, extra, repo_root)
    if args.ml_cmd == "hpo":
        return _cmd_ml_hpo(args, extra, repo_root)
    if args.ml_cmd == "chemprop":
        return _cmd_ml_chemprop(args, extra, repo_root)
    if args.ml_cmd == "tdc":
        return _cmd_ml_tdc(args, extra, repo_root)
    if args.ml_cmd == "reinvent":
        return _cmd_ml_reinvent(args, extra, repo_root)
    if args.ml_cmd == "gen-bench":
        return _cmd_ml_gen_bench(args, extra, repo_root)
    if args.ml_cmd in {"prepare-run", "ready"}:
        return _cmd_ml_prepare_run(args, extra, repo_root)
    if args.ml_cmd == "leaderboard":
        return _cmd_ml_leaderboard(args, extra, repo_root)
    if args.ml_cmd == "splits":
        return _cmd_ml_splits(args, extra, repo_root)
    if args.ml_cmd == "external-eval":
        return _cmd_ml_external_eval(args, extra, repo_root)
    if args.ml_cmd == "doctor":
        return _cmd_ml_doctor(args, repo_root)
    if args.ml_cmd == "status":
        return _cmd_ml_status(args, repo_root)
    return 1


def _cmd_ml_train(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    forwarded = [
        "--repo-root",
        str(repo_root),
        "--run-id",
        run_id,
        "--out-dir",
        str(Path(args.out_dir) if args.out_dir else _default_ml_out_dir(repo_root, run_id, "training_pass")),
        "--claim-mode",
        str(args.claim_mode),
    ]
    if args.run_dir:
        forwarded.extend(["--run-dir", str(args.run_dir)])
    else:
        forwarded.extend(["--run-dir", str(_default_run_dir(repo_root, run_id))])
    if args.sample_rows:
        forwarded.extend(["--sample-rows", str(args.sample_rows)])
    forwarded.extend(["--sample-seed", str(args.sample_seed), "--sample-strategy", str(args.sample_strategy)])
    if args.split_manifest:
        forwarded.extend(["--split-manifest", str(args.split_manifest)])
    for item in args.expert_split_manifest or []:
        forwarded.extend(["--expert-split-manifest", str(item)])
    if args.no_mlflow:
        forwarded.append("--no-mlflow")
    if args.no_feature_cache:
        forwarded.append("--no-feature-cache")
    if args.refresh_feature_cache:
        forwarded.append("--refresh-feature-cache")
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.run_ml_training_pass", forwarded)


def _cmd_ml_audit(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    forwarded = [
        "--dataset",
        str(args.dataset),
        "--label",
        str(args.label),
        "--out-dir",
        str(Path(args.out_dir) if args.out_dir else _default_ml_out_dir(repo_root, run_id, "audit_suite")),
        "--claim-mode",
        str(args.claim_mode),
    ]
    if args.feature_set:
        forwarded.extend(["--feature-set", str(args.feature_set)])
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.run_ml_audit_suite", forwarded)


def _cmd_ml_source_pu(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    forwarded = [
        "--run-dir",
        str(Path(args.run_dir) if args.run_dir else _default_run_dir(repo_root, run_id)),
        "--out-dir",
        str(Path(args.out_dir) if args.out_dir else _default_ml_out_dir(repo_root, run_id, "source_pu_suite")),
        "--claim-mode",
        str(args.claim_mode),
    ]
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.run_ml_source_pu_suite", forwarded)


def _cmd_ml_hpo(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    forwarded = [
        "--dataset",
        str(args.dataset),
        "--label",
        str(args.label),
        "--out-dir",
        str(Path(args.out_dir) if args.out_dir else _default_ml_out_dir(repo_root, run_id, "optuna_sweep")),
        "--claim-mode",
        str(args.claim_mode),
        "--repo-root",
        str(repo_root),
        "--pruner",
        str(args.pruner),
    ]
    if args.split_manifest:
        forwarded.extend(["--split-manifest", str(args.split_manifest)])
    if args.no_mlflow:
        forwarded.append("--no-mlflow")
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.run_ml_optuna_sweep", forwarded)


def _cmd_ml_chemprop(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    forwarded = [
        "--dataset",
        str(args.dataset),
        "--label",
        str(args.label),
        "--out-dir",
        str(Path(args.out_dir) if args.out_dir else _default_ml_out_dir(repo_root, run_id, "chemprop")),
    ]
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.run_chemprop_training", forwarded)


def _cmd_ml_tdc(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    tdc_root = _default_ml_out_dir(repo_root, run_id, "tdc")
    dataset_name = str(args.dataset_name or args.name or (Path(args.input_csv).stem if args.input_csv else "tdc_dataset"))
    out_path = Path(args.out) if args.out else tdc_root / f"{dataset_name}.csv"
    forwarded = ["--out", str(out_path)]
    if args.input_csv:
        forwarded.extend(["--input-csv", str(args.input_csv)])
    if args.name:
        forwarded.extend(["--name", str(args.name)])
    if args.dataset_name:
        forwarded.extend(["--dataset-name", str(args.dataset_name)])
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.stage_tdc_dataset", forwarded)


def _cmd_ml_reinvent(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    forwarded = [
        "--out-dir",
        str(Path(args.out_dir) if args.out_dir else _default_ml_out_dir(repo_root, run_id, "reinvent_generation")),
    ]
    if args.generated_smiles:
        forwarded.extend(["--generated-smiles", str(args.generated_smiles)])
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.run_reinvent_generation", forwarded)


def _cmd_ml_gen_bench(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    forwarded = [
        "--suite",
        str(args.suite),
        "--generated-smiles",
        str(args.generated_smiles),
        "--out-dir",
        str(Path(args.out_dir) if args.out_dir else _default_ml_out_dir(repo_root, run_id, f"{args.suite}_benchmark")),
    ]
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.benchmark_generated_smiles", forwarded)




def _cmd_ml_prepare_run(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    ml_root = Path(args.ml_root) if args.ml_root else _default_run_dir(repo_root, run_id) / "ml"
    run_dir = Path(args.run_dir) if args.run_dir else _default_run_dir(repo_root, run_id)
    forwarded = [
        "--run-id", run_id,
        "--run-dir", str(run_dir),
        "--ml-root", str(ml_root),
        "--repo-root", str(repo_root),
        "--claim-mode", str(args.claim_mode),
        "--sample-rows", str(args.sample_rows),
        "--sample-seed", str(args.sample_seed),
    ]
    if args.skip_splits:
        forwarded.append("--skip-splits")
    if args.skip_sample_train:
        forwarded.append("--skip-sample-train")
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.run_ml_readiness", forwarded)


def _cmd_ml_leaderboard(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    ml_root = Path(args.ml_root) if args.ml_root else _default_run_dir(repo_root, run_id) / "ml"
    out_dir = Path(args.out_dir) if args.out_dir else ml_root / "leaderboard"
    forwarded = ["--ml-root", str(ml_root), "--out-dir", str(out_dir)]
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.build_ml_leaderboard", forwarded)


def _cmd_ml_splits(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    run_dir = Path(args.run_dir) if args.run_dir else _default_run_dir(repo_root, run_id)
    if args.splits_cmd == "lock":
        out_dir = Path(args.out_dir) if args.out_dir else _default_run_dir(repo_root, run_id) / "ml" / "splits"
        forwarded = ["lock", "--run-id", run_id, "--run-dir", str(run_dir), "--out-dir", str(out_dir), "--seed", str(args.seed), "--validation-fraction", str(args.validation_fraction)]
        if args.experts:
            forwarded.append("--experts")
            forwarded.extend(args.experts)
        if args.splits:
            forwarded.append("--splits")
            forwarded.extend(args.splits)
    else:
        split_root = Path(args.split_root) if args.split_root else _default_run_dir(repo_root, run_id) / "ml" / "splits"
        forwarded = ["validate", "--run-dir", str(run_dir), "--split-root", str(split_root)]
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.run_ml_split_lock", forwarded)


def _cmd_ml_external_eval(args: argparse.Namespace, extra: list[str], repo_root: Path) -> int:
    run_id = str(args.run_id)
    out_dir = Path(args.out_dir) if args.out_dir else _default_run_dir(repo_root, run_id) / "ml" / "external_eval" / str(args.name)
    forwarded = ["--model-dir", str(args.model_dir), "--dataset", str(args.dataset), "--name", str(args.name), "--out-dir", str(out_dir)]
    if args.label:
        forwarded.extend(["--label", str(args.label)])
    forwarded.extend(extra)
    return _run_module_main("analysis.cli.run_ml_external_eval", forwarded)

def _cmd_ml_doctor(args: argparse.Namespace, repo_root: Path) -> int:
    from analysis.ml.preflight import build_ml_doctor_report, write_doctor_report

    run_id = str(args.run_id)
    run_dir = Path(args.run_dir) if args.run_dir else _default_run_dir(repo_root, run_id)
    report = build_ml_doctor_report(
        run_id=run_id,
        run_dir=run_dir,
        dataset=Path(args.dataset) if args.dataset else None,
        label_col=args.label,
        feature_set=args.feature_set,
        claim_mode=args.claim_mode,
        strict=bool(args.strict),
    )
    if args.write_report:
        report_root = (Path(args.run_dir) / "ml" / "readiness") if args.run_dir else _default_ml_out_dir(repo_root, run_id, "readiness")
        report_path = report_root / "ml_doctor_report.json"
        write_doctor_report(report, report_path)
        report["report_path"] = str(report_path)
    if args.json_out:
        print(json.dumps(report, indent=2, sort_keys=True, default=str))
    else:
        _print_doctor_report(report)
        if args.write_report:
            print(f"Report: {report['report_path']}")
    return 1 if report.get("blockers") else 0


def _print_doctor_report(report: dict[str, Any]) -> None:
    print(f"ML doctor status: {report['status']}")
    print(f"Run dir: {report['run_dir']}")
    for blocker in report.get("blockers", []):
        print(f"BLOCKER: {blocker}")
    for warning in report.get("warnings", []):
        print(f"WARNING: {warning}")
    if not report.get("blockers") and not report.get("warnings"):
        print("No ML preflight blockers or warnings detected.")

def _cmd_ml_status(args: argparse.Namespace, repo_root: Path) -> int:
    run_id = str(args.run_id)
    root = Path(args.out_dir) if args.out_dir else _default_run_dir(repo_root, run_id) / "ml"
    manifest_paths = [
        root / "training_pass" / "ml_training_pass_manifest.json",
        root / "source_pu_suite" / "ml_source_pu_suite_manifest.json",
        root / "audit_suite" / "ml_audit_suite_manifest.json",
        root / "optuna_sweep" / "ml_optuna_sweep_manifest.json",
        root / "optuna_sweep" / "hpo" / "optuna" / "optuna_sweep_manifest.json",
        root / "chemprop" / "chemprop_training_manifest.json",
        root / "reinvent_generation" / "reinvent_generation_manifest.json",
        root / "ml_feature_metadata_manifest.json",
        root / "readiness" / "ml_readiness_manifest.json",
        root / "readiness" / "ml_doctor_report.json",
        root / "leaderboard" / "leaderboard_manifest.json",
        root / "splits" / "split_lock_manifest.json",
    ]
    manifest_paths.extend(sorted(root.glob("*/ml_optuna_sweep_manifest.json")))
    manifest_paths.extend(sorted(root.glob("*/hpo/optuna/optuna_sweep_manifest.json")))
    manifest_paths.extend(sorted(root.glob("*/chemprop_training_manifest.json")))
    manifest_paths.extend(sorted(root.glob("tdc/*.manifest.json")))
    manifest_paths.extend(sorted(root.glob("*/reinvent_generation_manifest.json")))
    manifest_paths.extend(sorted(root.glob("*/*_benchmark_manifest.json")))
    manifest_paths.extend(sorted(root.glob("external_eval/*/external_eval_manifest.json")))
    summaries = []
    seen_paths: set[Path] = set()
    for path in manifest_paths:
        if path in seen_paths:
            continue
        seen_paths.add(path)
        if not path.exists():
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            summaries.append({"path": str(path), "status": "unreadable", "error": str(exc)})
            continue
        summaries.append(_summarize_manifest(path, payload))
    report = {"run_id": run_id, "ml_root": str(root), "manifests": summaries}
    if args.json_out:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(f"ML root: {root}")
        if not summaries:
            print("No ML manifests found.")
        for summary in summaries:
            print(f"{summary['kind']}: {summary['path']}")
            if summary.get("status"):
                print(f"  status: {summary['status']}")
            if summary.get("adapter"):
                print(f"  adapter: {summary['adapter']}")
            if summary.get("rows") is not None:
                print(f"  rows: {summary['rows']}")
            if summary.get("failed"):
                print(f"  failed: {summary['failed']}")
            if summary.get("skipped"):
                print(f"  skipped: {summary['skipped']}")
            if summary.get("trained"):
                print(f"  trained: {summary['trained']}")
    return 0


def _summarize_manifest(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    entries = payload.get("experts") or payload.get("results") or []
    failed = [str(item.get("expert") or item.get("name") or item.get("suite")) for item in entries if item.get("status") == "failed"]
    skipped = [str(item.get("expert") or item.get("name") or item.get("suite")) for item in entries if item.get("status") == "skipped"]
    trained = [str(item.get("expert") or item.get("name") or item.get("suite")) for item in entries if item.get("status") == "trained"]
    return {
        "kind": path.stem,
        "path": str(path),
        "status": payload.get("status") or payload.get("overall_status"),
        "adapter": payload.get("adapter"),
        "rows": payload.get("rows") or payload.get("staged_rows") or payload.get("n_trials"),
        "failed": failed,
        "skipped": skipped,
        "trained": trained,
        "claim_status": (payload.get("claim_readiness") or {}).get("overall_status"),
    }
