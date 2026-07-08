from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.calibration.metrics import calibration_metrics
from analysis.ml.labels import binary_label_series


def _source_mask(series: pd.Series, tokens: list[str]) -> pd.Series:
    text = series.fillna("").astype(str)
    mask = pd.Series(False, index=series.index)
    for token in tokens:
        value = str(token).strip()
        if not value:
            continue
        mask |= text.eq(value) | text.str.split(";").map(lambda parts: value in parts)
    return mask


def _balance(df: pd.DataFrame, label_col: str, group_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group_key, group in df.groupby(group_cols, dropna=False):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        row = {col: value for col, value in zip(group_cols, group_key, strict=False)}
        row.update(
            {
                "n": int(len(group)),
                "positive_rate": float(group[label_col].mean()) if len(group) else float("nan"),
                "n_drugs": int(group["drug_id"].nunique()) if "drug_id" in group else 0,
                "n_targets": int(group["target_id"].nunique()) if "target_id" in group else 0,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def _feature_shift(
    train: pd.DataFrame,
    test: pd.DataFrame,
    feature_cols: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for col in feature_cols:
        if col not in train.columns or col not in test.columns:
            continue
        train_values = pd.to_numeric(train[col], errors="coerce")
        test_values = pd.to_numeric(test[col], errors="coerce")
        pooled_iqr = pd.concat([train_values, test_values]).quantile(0.75) - pd.concat([train_values, test_values]).quantile(0.25)
        if pd.isna(pooled_iqr) or pooled_iqr == 0:
            pooled_iqr = 1.0
        rows.append(
            {
                "feature": col,
                "train_nonmissing": int(train_values.notna().sum()),
                "test_nonmissing": int(test_values.notna().sum()),
                "train_missing_fraction": float(train_values.isna().mean()),
                "test_missing_fraction": float(test_values.isna().mean()),
                "train_median": float(train_values.median()) if train_values.notna().any() else float("nan"),
                "test_median": float(test_values.median()) if test_values.notna().any() else float("nan"),
                "median_shift_iqr_units": float((test_values.median() - train_values.median()) / pooled_iqr)
                if train_values.notna().any() and test_values.notna().any()
                else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def _overlap(train: pd.DataFrame, test: pd.DataFrame, col: str) -> dict[str, object]:
    if col not in train.columns or col not in test.columns:
        return {"field": col, "available": False}
    train_values = set(train[col].dropna().astype(str))
    test_values = set(test[col].dropna().astype(str))
    overlap = train_values & test_values
    return {
        "field": col,
        "available": True,
        "train_unique": len(train_values),
        "test_unique": len(test_values),
        "overlap_unique": len(overlap),
        "test_overlap_fraction": len(overlap) / len(test_values) if test_values else 0.0,
    }


def _prediction_audit(predictions_path: str | Path | None, label_col: str) -> tuple[pd.DataFrame, dict[str, Any]]:
    if predictions_path is None or not Path(predictions_path).exists():
        return pd.DataFrame(), {"status": "missing_predictions"}
    pred = pd.read_csv(predictions_path, low_memory=False)
    if label_col not in pred.columns or "ml_prediction_score" not in pred.columns:
        return pred, {"status": "missing_label_or_prediction_column"}
    labels = binary_label_series(pred[label_col])
    pred = pred.loc[labels.notna()].copy()
    pred[label_col] = labels.dropna().astype(int)
    rows = []
    group_col = "label_source" if "label_source" in pred.columns else None
    groups = pred.groupby(group_col, dropna=False) if group_col else [("all", pred)]
    for source, group in groups:
        scores = pd.to_numeric(group["ml_prediction_score"], errors="coerce")
        valid = group[label_col].notna() & scores.notna()
        metric_values = (
            calibration_metrics(scores.loc[valid].tolist(), group.loc[valid, label_col].astype(int).tolist())
            if valid.sum() and group.loc[valid, label_col].nunique() > 1
            else {}
        )
        rows.append(
            {
                "source": source,
                "n": int(len(group)),
                "positive_rate": float(group[label_col].mean()),
                "mean_prediction": float(scores.mean()),
                "median_prediction": float(scores.median()),
                **metric_values,
            }
        )
    return pd.DataFrame(rows), {"status": "written", "n_predictions": int(len(pred))}


def audit_source_transfer(
    dataset_path: str | Path,
    label_col: str,
    out_dir: str | Path,
    *,
    train_sources: list[str],
    test_sources: list[str],
    predictions_path: str | Path | None = None,
    feature_cols: list[str] | None = None,
) -> dict[str, Any]:
    df = pd.read_csv(dataset_path, low_memory=False)
    if "label_source" not in df.columns:
        raise ValueError("source-transfer audit requires label_source")
    labels = binary_label_series(df[label_col])
    data = df.loc[labels.notna()].copy()
    data[label_col] = labels.dropna().astype(int)
    train = data.loc[_source_mask(data["label_source"], train_sources)].copy()
    test = data.loc[_source_mask(data["label_source"], test_sources) & ~data.index.isin(train.index)].copy()
    if train.empty or test.empty:
        raise ValueError(f"empty train/test source split: train={len(train)} test={len(test)}")

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    _balance(data, label_col, ["label_source"]).to_csv(out / "source_label_balance.csv", index=False)
    _balance(data, label_col, ["target_id", "label_source"]).to_csv(out / "target_source_label_balance.csv", index=False)
    context_cols = [col for col in ["label_source", "assay_type", "endpoint_type", "activity_type", "standard_type"] if col in data.columns]
    if context_cols:
        _balance(data, label_col, context_cols).to_csv(out / "assay_context_label_balance.csv", index=False)
    _feature_shift(train, test, feature_cols or ["atlas_score", "consensus_score", "free_cmax_um", "cmax_um", "fraction_unbound_plasma"]).to_csv(
        out / "feature_shift_train_vs_test.csv",
        index=False,
    )
    pred_audit, pred_manifest = _prediction_audit(predictions_path, label_col)
    if not pred_audit.empty:
        pred_audit.to_csv(out / "prediction_calibration_by_source.csv", index=False)

    train_rate = float(train[label_col].mean())
    test_rate = float(test[label_col].mean())
    findings = []
    if train_rate > 0 and test_rate > 0 and max(train_rate, test_rate) / min(train_rate, test_rate) >= 5:
        findings.append("large_label_prevalence_shift")
    target_overlap = _overlap(train, test, "target_id")
    drug_overlap = _overlap(train, test, "drug_id")
    target_test_overlap = pd.to_numeric(
        pd.Series([target_overlap.get("test_overlap_fraction", 0.0)]),
        errors="coerce",
    ).fillna(0.0).iloc[0]
    drug_test_overlap = pd.to_numeric(
        pd.Series([drug_overlap.get("test_overlap_fraction", 0.0)]),
        errors="coerce",
    ).fillna(0.0).iloc[0]
    if target_overlap.get("available") and float(target_test_overlap) < 0.5:
        findings.append("low_target_overlap")
    if drug_overlap.get("available") and float(drug_test_overlap) < 0.2:
        findings.append("low_drug_overlap")
    if pred_manifest.get("status") == "written" and not pred_audit.empty:
        weighted_pred = float((pred_audit["mean_prediction"] * pred_audit["n"]).sum() / pred_audit["n"].sum())
        if test_rate > 0 and weighted_pred / test_rate >= 5:
            findings.append("prediction_prevalence_miscalibration")
        if "ECE" in pred_audit.columns and pd.to_numeric(pred_audit["ECE"], errors="coerce").dropna().gt(0.20).any():
            findings.append("source_specific_calibration_error_high")

    manifest: dict[str, Any] = {
        "dataset_path": str(dataset_path),
        "label_col": label_col,
        "train_sources": train_sources,
        "test_sources": test_sources,
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "train_positive_rate": train_rate,
        "test_positive_rate": test_rate,
        "target_overlap": target_overlap,
        "drug_overlap": drug_overlap,
        "prediction_audit": pred_manifest,
        "findings": findings,
        "recommended_strategy": [
            "Do not claim ChEMBL/Papyrus to ToxCast transfer as a positive result until source-calibrated and target-matched performance improves.",
            "Report ToxCast as an independent in-vitro benchmark with its own prevalence and calibration table.",
            "Train source-specific or source-calibrated experts rather than one pooled bioactivity classifier when label definitions differ.",
            "Use BigBind/ChEMBL/Papyrus for docking-vs-bioactivity and keep ToxCast as orthogonal HTS validation.",
        ],
    }
    (out / "source_transfer_audit_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
