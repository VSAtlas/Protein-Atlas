from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_MODELS = [
    "logistic_regression",
    "elastic_net",
    "random_forest",
    "xgboost",
    "lightgbm",
    "catboost",
    "ebm",
    "svm",
    "knn",
]

EXPERT_DEFAULTS = {
    "binding": {
        "label": "spd_binding_label",
        "feature_set": "spd_binding_nonleaky",
    },
    "exposure": {
        "label": "spd_exposure_relevant",
        "feature_set": "spd_exposure_nonleaky",
    },
    "tissue": {
        "label": "tissue_site_label",
        "feature_set": "spd_tissue_site_nonleaky",
    },
    "mechanism": {
        "label": "mechanism_ml_label",
        "feature_set": "flat_mechanism_baseline",
    },
}


def _read_record(model_dir: Path) -> dict[str, Any]:
    path = model_dir / "model_run_record.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _run_training(
    command: list[str],
    *,
    log_path: Path,
    timeout: int | None,
    thread_limit: int,
) -> tuple[int, str | None]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    bounded_threads = str(max(1, min(int(thread_limit), 32)))
    for name in [
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ]:
        env[name] = bounded_threads
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(command) + "\n\n")
        try:
            completed = subprocess.run(
                command,
                check=False,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as exc:
            handle.write(f"\nTIMEOUT after {timeout} seconds\n")
            return 124, str(exc)
    return int(completed.returncode), None


def _summary_row(
    *,
    expert: str,
    model: str,
    status: str,
    model_dir: Path,
    command: list[str],
    record: dict[str, Any],
    error: str | None = None,
) -> dict[str, Any]:
    return {
        "expert": expert,
        "model_type": model,
        "status": status,
        "label_used": record.get("label_used"),
        "feature_set": record.get("feature_set"),
        "split_method": record.get("split_method"),
        "PU_strategy": record.get("PU_strategy"),
        "calibration_method": record.get("calibration_method"),
        "number_of_rows": record.get("number_of_rows"),
        "number_of_positives": record.get("number_of_positives"),
        "n_train": record.get("n_train"),
        "n_test": record.get("n_test"),
        "AUROC": record.get("AUROC"),
        "PR_AUC": record.get("PR_AUC"),
        "precision_at_K": record.get("precision_at_K"),
        "enrichment_at_K": record.get("enrichment_at_K"),
        "Brier_score": record.get("Brier_score"),
        "selected_model": record.get("selected_model"),
        "dataset_version": record.get("dataset_version"),
        "git_commit": record.get("git_commit"),
        "model_dir": str(model_dir),
        "command": " ".join(command),
        "error": error,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run a model-family sweep on cached Atlas ML expert tables.")
    parser.add_argument("--binding-dataset")
    parser.add_argument("--exposure-dataset")
    parser.add_argument("--tissue-dataset")
    parser.add_argument("--mechanism-dataset")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--experts", nargs="+", choices=sorted(EXPERT_DEFAULTS), default=sorted(EXPERT_DEFAULTS))
    parser.add_argument("--models", nargs="+", default=DEFAULT_MODELS)
    parser.add_argument("--split", default="drug_holdout")
    parser.add_argument("--calibration", default="sigmoid", choices=["none", "raw", "sigmoid", "isotonic"])
    parser.add_argument("--pu-mode", default="standard_binary")
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--model-n-jobs", type=int, default=16)
    parser.add_argument("--per-run-timeout", type=int, default=3600)
    parser.add_argument("--claim-mode", default="exploratory", choices=["exploratory", "publication"])
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--strict-feature-set", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--mechanism-label", default=EXPERT_DEFAULTS["mechanism"]["label"])
    parser.add_argument("--mechanism-feature-set", default=EXPERT_DEFAULTS["mechanism"]["feature_set"])
    parser.add_argument(
        "--group-reweight-cols",
        nargs="*",
        default=None,
        help="Optional capped inverse-frequency reweighting over source/family strata.",
    )
    parser.add_argument(
        "--group-reweight-include-label",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include the supervised label in group-reweight strata.",
    )
    parser.add_argument("--group-reweight-max-factor", type=float, default=5.0)
    parser.add_argument(
        "--tissue-pu-mode",
        choices=["bagging_pu", "stratified_bagging_pu"],
        default="stratified_bagging_pu",
    )
    parser.add_argument("--tissue-pu-bags", type=int, default=50)
    parser.add_argument("--tissue-pu-unlabeled-ratio", type=float, default=1.0)
    parser.add_argument("--tissue-top-k", nargs="+", type=int, default=[10, 20, 50])
    parser.add_argument("--tissue-max-eval-background", type=int, default=10000)
    args = parser.parse_args(argv)

    datasets = {
        "binding": args.binding_dataset,
        "exposure": args.exposure_dataset,
        "tissue": args.tissue_dataset,
        "mechanism": args.mechanism_dataset,
    }
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []

    for expert in args.experts:
        dataset = datasets.get(expert)
        if not dataset:
            rows.append(
                {
                    "expert": expert,
                    "status": "skipped_missing_dataset",
                    "model_type": None,
                    "error": f"no --{expert}-dataset provided",
                }
            )
            continue
        label = args.mechanism_label if expert == "mechanism" else EXPERT_DEFAULTS[expert]["label"]
        feature_set = args.mechanism_feature_set if expert == "mechanism" else EXPERT_DEFAULTS[expert]["feature_set"]
        for model in args.models:
            model_dir = out / expert / model / "model"
            log_path = out / expert / model / "train.log"
            if expert == "tissue":
                command = [
                    sys.executable,
                    "-m",
                    "analysis.cli.run_tissue_pu_recovery",
                    "--dataset",
                    str(dataset),
                    "--label",
                    label,
                    "--feature-set",
                    feature_set,
                    "--model",
                    model,
                    "--split",
                    args.split,
                    "--pu-mode",
                    args.tissue_pu_mode,
                    "--pu-bags",
                    str(args.tissue_pu_bags),
                    "--pu-unlabeled-ratio",
                    str(args.tissue_pu_unlabeled_ratio),
                    "--top-k",
                    *[str(k) for k in args.tissue_top_k],
                    "--max-eval-background",
                    str(args.tissue_max_eval_background),
                    "--repo-root",
                    args.repo_root,
                    "--out-dir",
                    str(model_dir),
                ]
            else:
                command = [
                    sys.executable,
                    "-m",
                    "analysis.cli.train_ml_model",
                    "--dataset",
                    str(dataset),
                    "--label",
                    label,
                    "--feature-set",
                    feature_set,
                    "--model",
                    model,
                    "--model-n-jobs",
                    str(args.model_n_jobs),
                    "--split",
                    args.split,
                    "--pu-mode",
                    args.pu_mode,
                    "--bootstraps",
                    str(args.bootstraps),
                    "--calibration",
                    args.calibration,
                    "--claim-mode",
                    args.claim_mode,
                    "--repo-root",
                    args.repo_root,
                    "--out-dir",
                    str(model_dir),
                ]
            if expert != "tissue" and args.group_reweight_cols:
                command.extend(["--group-reweight-cols", *args.group_reweight_cols])
                command.extend(["--group-reweight-max-factor", str(args.group_reweight_max_factor)])
                if not args.group_reweight_include_label:
                    command.append("--no-group-reweight-include-label")
            if args.strict_feature_set:
                command.append("--strict-feature-set")
            if args.skip_existing and (model_dir / "model_run_record.json").exists():
                record = _read_record(model_dir)
                rows.append(
                    _summary_row(
                        expert=expert,
                        model=model,
                        status="skipped_existing",
                        model_dir=model_dir,
                        command=command,
                        record=record,
                    )
                )
                records.append({"expert": expert, "model_type": model, **record})
                continue
            returncode, error = _run_training(
                command,
                log_path=log_path,
                timeout=args.per_run_timeout,
                thread_limit=args.model_n_jobs,
            )
            record = _read_record(model_dir)
            status = "ok" if returncode == 0 and record else f"failed_returncode_{returncode}"
            row = _summary_row(
                expert=expert,
                model=model,
                status=status,
                model_dir=model_dir,
                command=command,
                record=record,
                error=error,
            )
            rows.append(row)
            if status == "ok":
                records.append({"expert": expert, "model_type": model, **record})
            else:
                failures.append(row)
            pd.DataFrame(rows).to_csv(out / "model_family_sweep_summary.csv", index=False)
            pd.DataFrame(failures).to_csv(out / "model_family_sweep_failures.csv", index=False)
            with (out / "model_family_sweep_records.jsonl").open("w", encoding="utf-8") as handle:
                for record_row in records:
                    handle.write(json.dumps(record_row, sort_keys=True, default=str) + "\n")

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "experts": args.experts,
        "models": args.models,
        "split": args.split,
        "calibration": args.calibration,
        "pu_mode": args.pu_mode,
        "bootstraps": args.bootstraps,
        "model_n_jobs": args.model_n_jobs,
        "per_run_timeout": args.per_run_timeout,
        "strict_feature_set": bool(args.strict_feature_set),
        "group_reweight_cols": args.group_reweight_cols or [],
        "group_reweight_include_label": bool(args.group_reweight_include_label),
        "group_reweight_max_factor": float(args.group_reweight_max_factor),
        "tissue_pu_mode": args.tissue_pu_mode,
        "tissue_pu_bags": int(args.tissue_pu_bags),
        "tissue_pu_unlabeled_ratio": float(args.tissue_pu_unlabeled_ratio),
        "tissue_top_k": [int(k) for k in args.tissue_top_k],
        "tissue_max_eval_background": int(args.tissue_max_eval_background),
        "summary": str(out / "model_family_sweep_summary.csv"),
        "records_jsonl": str(out / "model_family_sweep_records.jsonl"),
        "failures": str(out / "model_family_sweep_failures.csv"),
    }
    (out / "model_family_sweep_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
