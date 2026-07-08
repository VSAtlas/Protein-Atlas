from __future__ import annotations

import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

LEDGER_COLUMNS = [
    "run_id",
    "date",
    "git_commit",
    "dataset_version",
    "dataset_path",
    "label_used",
    "feature_set",
    "excluded_columns",
    "split_method",
    "PU_strategy",
    "model_type",
    "selected_model",
    "hyperparameters",
    "calibration_method",
    "number_of_rows",
    "number_of_positives",
    "n_train",
    "n_test",
    "AUROC",
    "PR_AUC",
    "precision_at_K",
    "precision_K",
    "enrichment_at_K",
    "enrichment_K",
    "Brier_score",
    "notes",
    "model_dir",
]


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _read_metric_csv(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        frame = pd.read_csv(path)
    except Exception:
        return {}
    if {"metric", "value"}.issubset(frame.columns):
        return {str(row.metric): row.value for row in frame.itertuples(index=False)}
    return {}


def _git_commit(repo_root: str | Path | None) -> str:
    if repo_root is None:
        return "unknown"
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(repo_root),
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception:
        return "unknown"
    value = result.stdout.strip()
    return value if result.returncode == 0 and value else "unknown"


def _label_count(manifest: dict[str, Any], positive_key: str = "1") -> int:
    counts = manifest.get("label_counts") or {}
    if not isinstance(counts, dict):
        return 0
    for key in [positive_key, "1.0", "True", "true"]:
        if key in counts:
            try:
                return int(float(counts[key]))
            except Exception:
                return 0
    return 0


def _topk_metrics(pred_path: Path, label_col: str, *, score_col: str = "ml_prediction_score", k: int = 20) -> dict[str, Any]:
    if not pred_path.exists():
        return {"precision_at_K": None, "enrichment_at_K": None, "precision_K": int(k), "enrichment_K": int(k)}
    try:
        pred = pd.read_csv(pred_path)
    except Exception:
        return {"precision_at_K": None, "enrichment_at_K": None, "precision_K": int(k), "enrichment_K": int(k)}
    if label_col not in pred.columns or score_col not in pred.columns:
        return {"precision_at_K": None, "enrichment_at_K": None, "precision_K": int(k), "enrichment_K": int(k)}
    work = pred[[label_col, score_col]].copy()
    work[label_col] = pd.to_numeric(work[label_col], errors="coerce")
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    work = work.dropna(subset=[label_col, score_col])
    if work.empty:
        return {"precision_at_K": None, "enrichment_at_K": None, "precision_K": int(k), "enrichment_K": int(k)}
    top_n = min(max(1, int(k)), len(work))
    ranked = work.sort_values(score_col, ascending=False)
    precision = float(ranked.head(top_n)[label_col].mean())
    prevalence = float(work[label_col].mean())
    enrichment = precision / prevalence if prevalence else 0.0
    return {
        "precision_at_K": precision,
        "enrichment_at_K": enrichment,
        "precision_K": int(top_n),
        "enrichment_K": int(top_n),
    }


def build_model_run_record(
    model_dir: str | Path,
    *,
    repo_root: str | Path | None = None,
    run_id: str | None = None,
    task: str | None = None,
    notes: str = "",
    precision_k: int = 20,
    extra_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    root = Path(model_dir)
    card = _read_json(root / "model_card.json")
    dataset = _read_json(root / "dataset_version_manifest.json")
    split = _read_json(root / "split_manifest.json")
    hpo = _read_json(root / "hpo_manifest.json")
    calibration = _read_json(root / "probability_calibration_manifest.json")
    metrics = dict(card.get("metrics") or {})
    metrics.update(_read_metric_csv(root / "model_metrics.csv"))
    label_col = str(card.get("label_col") or dataset.get("label_col") or "")
    topk = _topk_metrics(root / "model_predictions.csv", label_col, k=precision_k) if label_col else {}
    registry = _read_json(root / "registry.json")
    params = registry.get("params") or {}
    hyperparameters = hpo.get("best_params") or hpo.get("effective_model_params") or params.get("model_params") or {}
    split_summary = split.get("split_summary") or {}
    if not isinstance(split_summary, dict):
        split_summary = {}
    resolved_run_id = run_id
    if resolved_run_id is None and extra_config:
        raw_run_id = extra_config.get("run_id")
        resolved_run_id = str(raw_run_id) if raw_run_id is not None else None
    record: dict[str, Any] = {
        "run_id": resolved_run_id,
        "date": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(repo_root),
        "dataset_version": dataset.get("dataset_version_id") or card.get("dataset_version_id"),
        "dataset_path": dataset.get("dataset_path") or card.get("dataset_path"),
        "label_used": label_col,
        "feature_set": card.get("feature_set") or dataset.get("feature_set"),
        "excluded_columns": ";".join(str(v) for v in dataset.get("exclude_features", [])),
        "split_method": split.get("split_mode") or card.get("split_mode"),
        "PU_strategy": card.get("pu_mode") or params.get("pu_mode"),
        "model_type": card.get("model_type") or params.get("model_type"),
        "selected_model": card.get("selected_model") or params.get("selected_model"),
        "hyperparameters": json.dumps(hyperparameters, sort_keys=True, default=str),
        "calibration_method": calibration.get("method", "none"),
        "number_of_rows": dataset.get("n_rows"),
        "number_of_positives": _label_count(dataset),
        "n_train": split.get("n_train") or split_summary.get("n_train"),
        "n_test": split.get("n_test") or split_summary.get("n_test"),
        "AUROC": metrics.get("AUROC"),
        "PR_AUC": metrics.get("AUPRC"),
        "precision_at_K": topk.get("precision_at_K"),
        "precision_K": topk.get("precision_K", precision_k),
        "enrichment_at_K": topk.get("enrichment_at_K"),
        "enrichment_K": topk.get("enrichment_K", precision_k),
        "Brier_score": metrics.get("Brier"),
        "notes": notes,
        "model_dir": str(root),
    }
    if task:
        record["task"] = task
    if extra_config:
        record["config"] = extra_config
    return record


def write_model_run_record(
    model_dir: str | Path,
    *,
    repo_root: str | Path | None = None,
    run_id: str | None = None,
    task: str | None = None,
    notes: str = "",
    precision_k: int = 20,
    ledger_path: str | Path | None = None,
    extra_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    root = Path(model_dir)
    record = build_model_run_record(
        root,
        repo_root=repo_root,
        run_id=run_id,
        task=task,
        notes=notes,
        precision_k=precision_k,
        extra_config=extra_config,
    )
    root.mkdir(parents=True, exist_ok=True)
    (root / "model_run_record.json").write_text(json.dumps(record, indent=2, sort_keys=True, default=str), encoding="utf-8")
    pd.DataFrame([record]).to_csv(root / "model_run_record.csv", index=False)
    if ledger_path is not None:
        append_model_run_ledger(ledger_path, record)
    return record


def append_model_run_ledger(path: str | Path, record: dict[str, Any]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    exists = out.exists()
    columns = [*LEDGER_COLUMNS, *[key for key in record if key not in LEDGER_COLUMNS]]
    with out.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerow(record)
