from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.calibration.metrics import calibration_metrics
from analysis.ml.dataset_manifest import write_dataset_version_manifest
from analysis.ml.labels import binary_label_series
from analysis.ml.split_manifest import dataframe_content_hash
from analysis.ml.train_classifier_core import (
    _model_probabilities,
    _transform_design_matrix,
    load_preprocessing_artifact,
)
from analysis.ml.train_classifier_outputs import (
    _append_applicability_domain,
    _write_grouped_calibration_outputs,
    _write_reliability_outputs,
    write_subgroup_metrics,
)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def _write_status_csv(path: Path, *, status: str, reason: str) -> None:
    pd.DataFrame([{"status": status, "reason": reason}]).to_csv(path, index=False)


def _write_blocker_manifest(
    *,
    out: Path,
    reason: str,
    details: dict[str, Any],
    model_root: Path,
    dataset: str | Path,
    name: str,
    label_col: str | None,
) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    _write_json(out / "conformal_summary.json", {"status": "blocker", "reason": reason, **details})
    _write_json(out / "external_calibration_summary.json", {"status": "blocker", "reason": reason, **details})
    _write_status_csv(out / "model_applicability_domain_summary.csv", status="blocker", reason=reason)
    manifest = {
        "status": "blocked",
        "reason": reason,
        "details": details,
        "name": name,
        "model_dir": str(model_root),
        "dataset": str(dataset),
        "label_col": label_col,
        "outputs": {
            "calibration_summary": str(out / "external_calibration_summary.json"),
            "conformal_summary": str(out / "conformal_summary.json"),
            "applicability_domain": str(out / "model_applicability_domain_summary.csv"),
        },
    }
    _write_json(out / "external_eval_manifest.json", manifest)
    return manifest


def _validate_model_columns(model: Any, x_eval: pd.DataFrame, preprocessing: dict[str, Any]) -> None:
    design_columns = [str(col) for col in preprocessing.get("design_columns", [])]
    fitted_columns = [str(col) for col in getattr(model, "feature_names_in_", [])]
    if fitted_columns and fitted_columns != design_columns:
        raise ValueError("model fitted columns do not match preprocessing artifact design_columns")
    if list(x_eval.columns) != (fitted_columns or design_columns):
        raise ValueError("external design matrix columns do not match fitted model columns")


def _write_external_calibration_summary(
    *,
    pred: pd.DataFrame,
    label_col: str,
    metrics: dict[str, Any],
    out: Path,
) -> dict[str, Any]:
    labels = pd.to_numeric(pred[label_col], errors="coerce")
    if labels.nunique() < 2:
        summary = {
            "status": "skipped",
            "reason": "external_labels_single_class",
            "n_labeled": int(labels.notna().sum()),
            "n_classes": int(labels.nunique()),
        }
    else:
        summary = {
            "status": "ok",
            "n_labeled": int(labels.notna().sum()),
            "metrics": metrics,
            "outputs": {
                "reliability": str(out / "model_reliability_table.csv"),
                "grouped_calibration": str(out / "model_grouped_calibration.csv"),
            },
        }
    _write_json(out / "external_calibration_summary.json", summary)
    return summary


def _append_external_conformal_sets(
    pred: pd.DataFrame,
    *,
    model_root: Path,
    label_col: str,
    out: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    training_summary_path = model_root / "conformal_summary.json"
    if not training_summary_path.exists():
        summary = {"status": "skipped", "reason": "training_conformal_summary_missing"}
        _write_json(out / "conformal_summary.json", summary)
        return pred, summary
    training_summary = json.loads(training_summary_path.read_text(encoding="utf-8"))
    qhat = training_summary.get("qhat")
    alpha = training_summary.get("alpha", 0.10)
    if training_summary.get("status") != "ok" or qhat is None:
        summary = {
            "status": "skipped",
            "reason": "training_conformal_unavailable",
            "training_status": training_summary.get("status"),
            "training_reason": training_summary.get("reason"),
        }
        _write_json(out / "conformal_summary.json", summary)
        return pred, summary

    enriched = pred.copy()
    sets: list[str] = []
    sizes: list[int] = []
    for prob in pd.to_numeric(enriched["ml_prediction_score"], errors="coerce"):
        labels_in_set: list[str] = []
        if pd.notna(prob):
            if float(prob) <= float(qhat):
                labels_in_set.append("0")
            if 1.0 - float(prob) <= float(qhat):
                labels_in_set.append("1")
        sets.append(";".join(labels_in_set))
        sizes.append(len(labels_in_set))
    enriched["conformal_prediction_set"] = sets
    enriched["conformal_set_size"] = sizes
    enriched["conformal_alpha"] = float(alpha)
    enriched["conformal_qhat"] = float(qhat)
    true_labels = pd.to_numeric(enriched[label_col], errors="coerce")
    enriched["conformal_contains_true"] = [
        (str(int(label)) in set_value.split(";")) if pd.notna(label) else None
        for label, set_value in zip(true_labels, sets)
    ]
    coverage_values = pd.Series(enriched["conformal_contains_true"]).dropna()
    summary = {
        "status": "ok",
        "source": "training_conformal_summary",
        "alpha": float(alpha),
        "target_coverage": float(1.0 - float(alpha)),
        "qhat": float(qhat),
        "n_prediction": int(len(enriched)),
        "empirical_coverage": float(coverage_values.astype(bool).mean()) if len(coverage_values) else None,
        "average_set_size": float(pd.Series(sizes).mean()) if sizes else None,
        "caveat": "External conformal sets reuse the training calibration qhat; validity depends on exchangeability with that calibration split.",
    }
    _write_json(out / "conformal_summary.json", summary)
    enriched[
        [
            col
            for col in [
                "drug_id",
                "target_id",
                "pdb_id",
                label_col,
                "ml_prediction_score",
                "conformal_prediction_set",
                "conformal_set_size",
                "conformal_contains_true",
            ]
            if col in enriched.columns
        ]
    ].to_csv(out / "conformal_predictions.csv", index=False)
    return enriched, summary


def _append_external_applicability_domain(
    pred: pd.DataFrame,
    *,
    x_eval: pd.DataFrame,
    preprocessing: dict[str, Any],
    model_root: Path,
    label_col: str,
    out: Path,
) -> tuple[pd.DataFrame, str]:
    matrix_ref = preprocessing.get("training_design_matrix")
    if not matrix_ref:
        _write_status_csv(
            out / "model_applicability_domain_summary.csv",
            status="skipped",
            reason="training_design_matrix_reference_missing",
        )
        return pred, "skipped"
    matrix_path = Path(str(matrix_ref))
    if not matrix_path.exists():
        relocated = model_root / matrix_path.name
        matrix_path = relocated if relocated.exists() else matrix_path
    if not matrix_path.exists():
        _write_status_csv(
            out / "model_applicability_domain_summary.csv",
            status="skipped",
            reason="training_design_matrix_missing",
        )
        return pred, "skipped"
    x_train = pd.read_csv(matrix_path, low_memory=False)
    design_columns = [str(col) for col in preprocessing.get("design_columns", [])]
    if list(x_train.columns) != design_columns:
        _write_status_csv(
            out / "model_applicability_domain_summary.csv",
            status="blocker",
            reason="training_design_matrix_columns_mismatch",
        )
        return pred, "blocker"
    enriched = _append_applicability_domain(pred, x_train, x_eval, label_col, out_path=out)
    if not (out / "model_applicability_domain_summary.csv").exists():
        _write_status_csv(
            out / "model_applicability_domain_summary.csv",
            status="skipped",
            reason="applicability_domain_backend_unavailable",
        )
        return enriched, "skipped"
    return enriched, "ok"


def run_external_model_eval(
    *,
    model_dir: str | Path,
    dataset: str | Path,
    label_col: str | None,
    name: str,
    out_dir: str | Path,
) -> dict[str, Any]:
    model_root = Path(model_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    card = json.loads((model_root / "model_card.json").read_text(encoding="utf-8"))
    features = list(card.get("features") or [])
    if not features:
        raise ValueError(f"model card has no features: {model_root / 'model_card.json'}")
    effective_label = label_col or str(card.get("label_col") or "")

    try:
        preprocessing = load_preprocessing_artifact(model_root)
    except ValueError as exc:
        _write_blocker_manifest(
            out=out,
            reason="missing_preprocessing_artifact",
            details={"error": str(exc)},
            model_root=model_root,
            dataset=dataset,
            name=name,
            label_col=effective_label,
        )
        raise

    frame = pd.read_csv(dataset, low_memory=False)
    raw_features = [str(feature) for feature in preprocessing.get("raw_features", [])]
    missing_raw = [feature for feature in raw_features if feature not in frame.columns]
    if missing_raw:
        _write_blocker_manifest(
            out=out,
            reason="missing_raw_features",
            details={"missing_raw_features": missing_raw},
            model_root=model_root,
            dataset=dataset,
            name=name,
            label_col=effective_label,
        )
        raise ValueError(f"external evaluation dataset is missing raw feature columns: {', '.join(missing_raw)}")

    data = frame.copy()
    labels = binary_label_series(data[effective_label]) if effective_label in data.columns else pd.Series([pd.NA] * len(data), index=data.index)
    labeled = labels.notna()
    if not labeled.any():
        raise ValueError("external evaluation requires at least one labeled row")
    data = data.loc[labeled].copy()
    data[effective_label] = labels.loc[labeled].astype(int)
    with (model_root / "trained_model.pkl").open("rb") as handle:
        model = pickle.load(handle)
    try:
        x_eval = _transform_design_matrix(data, preprocessing, strict_categories=True)
        _validate_model_columns(model, x_eval, preprocessing)
    except ValueError as exc:
        _write_blocker_manifest(
            out=out,
            reason="preprocessing_replay_failed",
            details={"error": str(exc)},
            model_root=model_root,
            dataset=dataset,
            name=name,
            label_col=effective_label,
        )
        raise

    probs = _model_probabilities(model, x_eval)
    metadata = [
        "drug_id", "target_id", "pdb_id", "scaffold_key", "chemical_cluster", "target_family",
        "protein_class", "label_source", "source_family", "upstream_source", "assay_type", "endpoint_type",
        "activity_type", "smiles", "canonical_smiles", effective_label,
    ]
    pred = data[[col for col in metadata if col in data.columns]].copy()
    pred["ml_prediction_score"] = probs
    metrics = calibration_metrics([float(v) for v in probs], data[effective_label].astype(int).tolist()) if data[effective_label].nunique() > 1 else {}
    pd.DataFrame([{"metric": key, "value": value} for key, value in metrics.items()]).to_csv(out / "external_metrics.csv", index=False)
    _write_reliability_outputs(pred, effective_label, out)
    _write_grouped_calibration_outputs(pred, effective_label, out)
    calibration_summary = _write_external_calibration_summary(pred=pred, label_col=effective_label, metrics=metrics, out=out)
    pred, conformal_summary = _append_external_conformal_sets(pred, model_root=model_root, label_col=effective_label, out=out)
    pred, applicability_status = _append_external_applicability_domain(
        pred,
        x_eval=x_eval,
        preprocessing=preprocessing,
        model_root=model_root,
        label_col=effective_label,
        out=out,
    )
    pred.to_csv(out / "external_predictions.csv", index=False)
    subgroup = write_subgroup_metrics(pred, effective_label, out)
    dataset_manifest = write_dataset_version_manifest(
        dataset_path=dataset,
        frame=frame,
        out_path=out / "external_dataset_manifest.json",
        label_col=effective_label,
        feature_set=str(card.get("feature_set") or ""),
        features=features,
        exclude_features=[],
        provenance={"stage": "external_eval", "model_dir": str(model_root), "name": name},
    )
    manifest = {
        "status": "evaluated",
        "name": name,
        "model_dir": str(model_root),
        "dataset": str(dataset),
        "dataset_hash": dataframe_content_hash(frame),
        "dataset_manifest": dataset_manifest,
        "label_col": effective_label,
        "feature_set": card.get("feature_set"),
        "n_rows": int(len(data)),
        "metrics": metrics,
        "calibration_status": calibration_summary.get("status"),
        "conformal_status": conformal_summary.get("status"),
        "applicability_domain_status": applicability_status,
        "preprocessing_artifact": str(model_root / "model_preprocessing.json"),
        "subgroup_metrics_rows": int(len(subgroup)),
        "outputs": {
            "predictions": str(out / "external_predictions.csv"),
            "metrics": str(out / "external_metrics.csv"),
            "calibration_summary": str(out / "external_calibration_summary.json"),
            "reliability": str(out / "model_reliability_table.csv"),
            "grouped_calibration": str(out / "model_grouped_calibration.csv"),
            "subgroup_metrics": str(out / "model_subgroup_metrics.csv"),
            "conformal_summary": str(out / "conformal_summary.json"),
            "applicability_domain": str(out / "model_applicability_domain_summary.csv"),
        },
    }
    _write_json(out / "external_eval_manifest.json", manifest)
    return manifest
