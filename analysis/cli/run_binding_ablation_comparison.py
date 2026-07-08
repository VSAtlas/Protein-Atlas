from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.feature_sets import RETAINED_CONTEXT_AUDIT_COLUMNS, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.splits import make_split, split_overlap_summary
from analysis.ml.train_classifier_core import _fit_design_matrix, _fit_model_with_optional_weights, _model, _model_probabilities, _transform_design_matrix
from analysis.statistics import auroc, average_precision, brier_score, ranked_binary_metrics

RAW_BASELINES = [
    {
        "ablation": "banana_only_raw",
        "score_candidates": ["banana_binding_probability", "banana_score_normalized", "banana_score"],
        "probability_like": True,
    },
    {
        "ablation": "consensus_score_raw",
        "score_candidates": ["consensus_score"],
        "probability_like": False,
    },
]

TRAINED_ABLATIONS = [
    ("banana_only_trained", "banana_only"),
    ("consensus_only_trained", "consensus_only"),
    ("atlas_consensus_without_banana", "spd_binding_consensus_context_no_banana"),
    ("atlas_consensus_without_banana_no_target_class", "spd_binding_consensus_context_no_banana_no_target_class"),
    ("rdkit_target_metadata_only", "spd_binding_rdkit_target_metadata"),
    ("rdkit_chemistry_only_no_target_class", "spd_binding_rdkit_chemistry_only"),
    ("atlas_consensus_plus_banana_scores", "spd_binding_consensus_banana_scores_only"),
    ("full_clean_binding_model", "spd_binding_nonleaky"),
    ("full_clean_binding_model_no_target_class", "spd_binding_nonleaky_no_target_class"),
]


def _load_record(model_dir: Path) -> dict[str, Any]:
    path = model_dir / "model_run_record.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _topk(labels: list[int], scores: list[float], *, k: int = 20) -> dict[str, float | int]:
    if not labels:
        return {"precision_at_K": math.nan, "enrichment_at_K": math.nan, "K": int(k)}
    pairs = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    top_n = min(int(k), len(pairs))
    prevalence = sum(labels) / len(labels) if labels else 0.0
    precision = sum(label for _score, label in pairs[:top_n]) / top_n if top_n else 0.0
    return {
        "precision_at_K": float(precision),
        "enrichment_at_K": float(precision / prevalence) if prevalence else 0.0,
        "K": int(top_n),
    }


def _raw_baseline_rows(
    dataset: Path,
    *,
    label_col: str,
    split: str,
    seed: int,
    top_k: int,
) -> list[dict[str, Any]]:
    frame = pd.read_csv(dataset, low_memory=False)
    frame["_atlas_observed_label"] = binary_label_series(frame[label_col])
    data = frame.loc[frame["_atlas_observed_label"].notna()].copy()
    data[label_col] = data["_atlas_observed_label"].astype(int)
    train_idx, test_idx = make_split(data, split_mode=split, seed=seed)
    train = data.loc[train_idx]
    test = data.loc[test_idx]
    split_summary = split_overlap_summary(train, test, split)
    rows: list[dict[str, Any]] = []
    for spec in RAW_BASELINES:
        score_col = next((col for col in spec["score_candidates"] if col in test.columns and test[col].notna().any()), None)
        if score_col is None:
            rows.append(
                {
                    "ablation": spec["ablation"],
                    "kind": "raw_score",
                    "split_method": split,
                    "status": "skipped_missing_score",
                    "score_col": None,
                }
            )
            continue
        work = test[[label_col, score_col]].dropna().copy()
        labels = pd.to_numeric(work[label_col], errors="coerce").astype(int).tolist()
        scores = pd.to_numeric(work[score_col], errors="coerce").tolist()
        status = "ok" if len(set(labels)) == 2 else "insufficient_labels"
        row: dict[str, Any] = {
            "ablation": spec["ablation"],
            "kind": "raw_score",
            "model_type": "none",
            "feature_set": None,
            "label_used": label_col,
            "score_col": score_col,
            "split_method": split,
            "status": status,
            "number_of_rows": int(len(data)),
            "number_of_positives": int(data[label_col].sum()),
            "n_train": int(len(train)),
            "n_test": int(len(test)),
            "n_scored_test": int(len(work)),
            "test_positive_rate": float(pd.Series(labels).mean()) if labels else math.nan,
            "split_passes_holdout": bool(split_summary.get("passes_holdout", False)),
        }
        if status == "ok":
            row.update(ranked_binary_metrics(scores, labels, [0.01, 0.05, 0.10]))
            row.update(_topk(labels, scores, k=top_k))
            row["Brier_score"] = float(brier_score(scores, labels)) if spec.get("probability_like") else math.nan
            row["AUROC"] = float(auroc(scores, labels))
            row["PR_AUC"] = float(average_precision(scores, labels))
        rows.append(row)
    return rows


def _fast_trained_row(
    *,
    frame: pd.DataFrame,
    label_col: str,
    ablation: str,
    feature_set: str,
    model: str,
    split: str,
    seed: int,
    model_n_jobs: int,
    top_k: int,
    strict_feature_set: bool,
) -> dict[str, Any]:
    data = frame.loc[frame["_atlas_observed_label"].notna()].copy()
    data[label_col] = data["_atlas_observed_label"].astype(int)
    requested = get_feature_set(feature_set)
    missing = [feature for feature in requested if feature not in data.columns]
    all_missing = [feature for feature in requested if feature in data.columns and data[feature].notna().sum() == 0]
    if strict_feature_set and (missing or all_missing):
        return {
            "ablation": ablation,
            "kind": "trained_model_metrics_only",
            "model_type": model,
            "feature_set": feature_set,
            "split_method": split,
            "status": "failed_strict_feature_set",
            "error": f"missing={missing[:10]} all_missing={all_missing[:10]}",
        }
    features = [feature for feature in requested if feature in data.columns and feature not in all_missing]
    retained_not_trained = sorted(col for col in RETAINED_CONTEXT_AUDIT_COLUMNS if col in data.columns and col not in features)
    if not features:
        return {
            "ablation": ablation,
            "kind": "trained_model_metrics_only",
            "model_type": model,
            "feature_set": feature_set,
            "split_method": split,
            "status": "failed_no_features",
        }
    try:
        train_idx, test_idx = make_split(data, split_mode=split, seed=seed)
        train = data.loc[train_idx].copy()
        test = data.loc[test_idx].copy()
        split_summary = split_overlap_summary(train, test, split)
        x_train, preprocessing = _fit_design_matrix(train, features)
        x_test = _transform_design_matrix(test, preprocessing)
        clf = _model(
            model,
            seed,
            class_weight="balanced",
            model_params={"n_jobs": model_n_jobs, "thread_count": model_n_jobs} if model_n_jobs > 0 else None,
        )
        clf = _fit_model_with_optional_weights(clf, x_train, train[label_col])
        scores = _model_probabilities(clf, x_test)
        labels = pd.to_numeric(test[label_col], errors="coerce").astype(int).tolist()
    except Exception as exc:  # noqa: BLE001 - metrics-only comparison records failures.
        return {
            "ablation": ablation,
            "kind": "trained_model_metrics_only",
            "model_type": model,
            "feature_set": feature_set,
            "split_method": split,
            "status": "failed",
            "error": str(exc),
        }
    status = "ok" if len(set(labels)) == 2 else "insufficient_labels"
    row: dict[str, Any] = {
        "ablation": ablation,
        "kind": "trained_model_metrics_only",
        "model_type": model,
        "feature_set": feature_set,
        "label_used": label_col,
        "score_col": "ml_prediction_score",
        "split_method": split,
        "status": status,
        "number_of_rows": int(len(data)),
        "number_of_positives": int(data[label_col].sum()),
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "n_scored_test": int(len(test)),
        "test_positive_rate": float(pd.Series(labels).mean()) if labels else math.nan,
        "split_passes_holdout": bool(split_summary.get("passes_holdout", False)),
        "retained_not_trained_features": ";".join(retained_not_trained),
    }
    if status == "ok":
        row.update(ranked_binary_metrics(scores, labels, [0.01, 0.05, 0.10]))
        row.update(_topk(labels, scores, k=top_k))
        row["Brier_score"] = float(brier_score(scores, labels))
        row["AUROC"] = float(auroc(scores, labels))
        row["PR_AUC"] = float(average_precision(scores, labels))
    return row


def _trained_row(
    *,
    dataset: Path,
    label_col: str,
    ablation: str,
    feature_set: str,
    model: str,
    split: str,
    out_dir: Path,
    seed: int,
    bootstraps: int,
    calibration: str,
    model_n_jobs: int,
    claim_mode: str,
    repo_root: Path,
    strict_feature_set: bool,
    skip_existing: bool,
) -> dict[str, Any]:
    model_dir = out_dir / split / ablation / model / "model"
    log_path = model_dir.parent / "train.log"
    if skip_existing and (model_dir / "model_run_record.json").exists():
        record = _load_record(model_dir)
        status = "skipped_existing" if record else "skipped_existing_missing_record"
    else:
        command = [
            sys.executable,
            "-m",
            "analysis.cli.train_ml_model",
            "--dataset",
            str(dataset),
            "--label",
            label_col,
            "--feature-set",
            feature_set,
            "--model",
            model,
            "--model-n-jobs",
            str(model_n_jobs),
            "--split",
            split,
            "--bootstraps",
            str(bootstraps),
            "--calibration",
            calibration,
            "--claim-mode",
            claim_mode,
            "--repo-root",
            str(repo_root),
            "--out-dir",
            str(model_dir),
        ]
        if strict_feature_set:
            command.append("--strict-feature-set")
        model_dir.parent.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        bounded_threads = str(max(1, min(int(model_n_jobs or 1), 32)))
        for name in [
            "OPENBLAS_NUM_THREADS",
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
        ]:
            env[name] = bounded_threads
        env.setdefault("RDKIT_LOG_LEVEL", "ERROR")
        with log_path.open("w", encoding="utf-8") as handle:
            handle.write(" ".join(command) + "\n\n")
            completed = subprocess.run(
                command,
                check=False,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
                text=True,
            )
        record = _load_record(model_dir)
        if completed.returncode != 0 or not record:
            return {
                "ablation": ablation,
                "kind": "trained_model",
                "model_type": model,
                "feature_set": feature_set,
                "split_method": split,
                "status": f"failed_returncode_{completed.returncode}",
                "error": f"see {log_path}",
                "model_dir": str(model_dir),
            }
        status = "ok"
    return {
        "ablation": ablation,
        "kind": "trained_model",
        "model_type": model,
        "feature_set": feature_set,
        "label_used": record.get("label_used", label_col),
        "split_method": record.get("split_method", split),
        "status": status,
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
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare BANANA, Atlas consensus, descriptor/context, and full clean binding models on the same splits."
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", default="spd_binding_label")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--splits",
        nargs="+",
        default=["target_holdout", "target_family_holdout", "chemical_cluster_holdout"],
    )
    parser.add_argument("--models", nargs="+", default=["logistic_regression", "random_forest", "lightgbm"])
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--calibration", choices=["none", "raw", "sigmoid", "isotonic"], default="sigmoid")
    parser.add_argument("--model-n-jobs", type=int, default=8)
    parser.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--strict-feature-set", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--metrics-only", action="store_true", help="Run fast comparison metrics without full model artifacts.")
    args = parser.parse_args(argv)

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for split in args.splits:
        try:
            rows.extend(_raw_baseline_rows(args.dataset, label_col=args.label, split=split, seed=args.seed, top_k=args.top_k))
        except Exception as exc:  # noqa: BLE001 - keep sweep running across splits.
            rows.append(
                {
                    "kind": "raw_score",
                    "split_method": split,
                    "status": "failed",
                    "error": str(exc),
                }
            )
        pd.DataFrame(rows).to_csv(out / "binding_ablation_comparison.csv", index=False)
        fast_frame = None
        if args.metrics_only:
            fast_frame = pd.read_csv(args.dataset, low_memory=False)
            fast_frame["_atlas_observed_label"] = binary_label_series(fast_frame[args.label])
        for ablation, feature_set in TRAINED_ABLATIONS:
            for model in args.models:
                if args.metrics_only:
                    row = _fast_trained_row(
                        frame=fast_frame,
                        label_col=args.label,
                        ablation=ablation,
                        feature_set=feature_set,
                        model=model,
                        split=split,
                        seed=args.seed,
                        model_n_jobs=args.model_n_jobs,
                        top_k=args.top_k,
                        strict_feature_set=bool(args.strict_feature_set),
                    )
                else:
                    row = _trained_row(
                        dataset=args.dataset,
                        label_col=args.label,
                        ablation=ablation,
                        feature_set=feature_set,
                        model=model,
                        split=split,
                        out_dir=out,
                        seed=args.seed,
                        bootstraps=args.bootstraps,
                        calibration=args.calibration,
                        model_n_jobs=args.model_n_jobs,
                        claim_mode=args.claim_mode,
                        repo_root=args.repo_root,
                        strict_feature_set=bool(args.strict_feature_set),
                        skip_existing=bool(args.skip_existing),
                    )
                rows.append(row)
                pd.DataFrame(rows).to_csv(out / "binding_ablation_comparison.csv", index=False)

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": str(args.dataset),
        "label": args.label,
        "splits": args.splits,
        "models": args.models,
        "raw_baselines": RAW_BASELINES,
        "trained_ablations": TRAINED_ABLATIONS,
        "metrics_only": bool(args.metrics_only),
        "summary": str(out / "binding_ablation_comparison.csv"),
        "notes": [
            "Raw consensus scores are ranking scores, not calibrated probabilities; Brier is omitted for those rows.",
            "BANANA-only raw rows evaluate frozen BANANA outputs without training an Atlas classifier.",
            "Trained rows use Atlas train_ml_model with identical split mode and seed for each ablation.",
            "No-target-class ablations keep target_family and protein_class in the dataset but exclude them from the model feature set.",
        ],
    }
    (out / "binding_ablation_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
