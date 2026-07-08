from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"status": "unreadable", "error": str(exc), "path": str(path)}


def _metrics_from_csv(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        frame = pd.read_csv(path)
    except Exception:
        return {}
    if {"metric", "value"}.issubset(frame.columns):
        return {str(row.metric): row.value for row in frame.itertuples(index=False)}
    return {}


def _mlflow_run_id(registry: dict[str, Any]) -> str | None:
    tracking = registry.get("external_tracking") or {}
    mlflow = tracking.get("mlflow") or {}
    return mlflow.get("run_id")


def collect_ml_leaderboard(ml_root: str | Path) -> pd.DataFrame:
    root = Path(ml_root)
    rows: list[dict[str, Any]] = []
    for card_path in sorted(root.glob("**/model_card.json")):
        model_dir = card_path.parent
        card = _read_json(card_path)
        dataset = _read_json(model_dir / "dataset_version_manifest.json")
        split = _read_json(model_dir / "split_manifest.json")
        readiness = _read_json(model_dir / "model_claim_readiness.json")
        registry = _read_json(model_dir / "registry.json")
        metrics = dict(card.get("metrics") or {})
        metrics.update(_metrics_from_csv(model_dir / "model_metrics.csv"))
        sample_manifest = None
        dataset_path = str(dataset.get("dataset_path") or "")
        sample_candidate = Path(dataset_path).parent / "sample_manifest.json" if dataset_path else None
        if sample_candidate and sample_candidate.exists():
            sample_manifest = _read_json(sample_candidate)
        rows.append(
            {
                "kind": "model",
                "expert": model_dir.parent.name if model_dir.name == "model" else model_dir.name,
                "model_dir": str(model_dir),
                "model_type": card.get("model_type"),
                "selected_model": card.get("selected_model"),
                "label_col": card.get("label_col"),
                "feature_set": card.get("feature_set"),
                "dataset_hash": dataset.get("dataset_version_id") or card.get("dataset_version_id"),
                "split_hash": split.get("split_manifest_hash") or split.get("dataset_hash"),
                "split_mode": split.get("split_mode") or card.get("split_mode"),
                "sampled": bool(sample_manifest and sample_manifest.get("sampled")),
                "claim_status": readiness.get("overall_status") or (card.get("claim_readiness") or {}).get("overall_status"),
                "AUROC": metrics.get("AUROC"),
                "AUPRC": metrics.get("AUPRC"),
                "Brier": metrics.get("Brier"),
                "ECE": metrics.get("ECE"),
                "adaptive_ECE": metrics.get("adaptive_ECE"),
                "conformal_status": _read_json(model_dir / "conformal_summary.json").get("status"),
                "artifact_path": str(model_dir / "model_artifacts.json"),
                "mlflow_run_id": _mlflow_run_id(registry),
            }
        )
    for manifest_path in sorted(root.glob("**/external_eval_manifest.json")):
        manifest = _read_json(manifest_path)
        metrics = manifest.get("metrics") or {}
        rows.append(
            {
                "kind": "external_eval",
                "expert": manifest.get("name"),
                "model_dir": manifest.get("model_dir"),
                "label_col": manifest.get("label_col"),
                "feature_set": manifest.get("feature_set"),
                "dataset_hash": manifest.get("dataset_hash"),
                "split_hash": "external_eval",
                "split_mode": "external_eval",
                "sampled": False,
                "claim_status": manifest.get("status"),
                "AUROC": metrics.get("AUROC"),
                "AUPRC": metrics.get("AUPRC"),
                "Brier": metrics.get("Brier"),
                "ECE": metrics.get("ECE"),
                "adaptive_ECE": metrics.get("adaptive_ECE"),
                "artifact_path": str(manifest_path),
                "mlflow_run_id": None,
            }
        )
    for manifest_path in sorted(root.glob("**/optuna_sweep_manifest.json")):
        manifest = _read_json(manifest_path)
        if manifest.get("status") == "trained":
            continue
        rows.append(
            {
                "kind": "hpo",
                "expert": manifest_path.parent.parent.parent.name if len(manifest_path.parts) > 3 else "optuna",
                "model_dir": manifest.get("model_dir"),
                "selected_model": None,
                "dataset_hash": None,
                "split_hash": None,
                "split_mode": None,
                "sampled": False,
                "claim_status": (manifest.get("claim_readiness") or {}).get("overall_status"),
                "AUROC": (manifest.get("metrics") or {}).get("AUROC"),
                "AUPRC": (manifest.get("metrics") or {}).get("AUPRC"),
                "artifact_path": str(manifest_path),
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(columns=["kind", "expert", "model_dir", "AUROC", "AUPRC", "claim_status", "artifact_path"])
    return frame


def write_ml_leaderboard(ml_root: str | Path, out_dir: str | Path | None = None) -> dict[str, Any]:
    root = Path(ml_root)
    out = Path(out_dir) if out_dir else root / "leaderboard"
    out.mkdir(parents=True, exist_ok=True)
    frame = collect_ml_leaderboard(root)
    csv_path = out / "model_leaderboard.csv"
    json_path = out / "model_leaderboard.json"
    frame.to_csv(csv_path, index=False)
    json_path.write_text(frame.to_json(orient="records", indent=2), encoding="utf-8")
    manifest = {"status": "written", "ml_root": str(root), "rows": int(len(frame)), "csv": str(csv_path), "json": str(json_path)}
    (out / "leaderboard_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
