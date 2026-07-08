from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.calibration.metrics import calibration_metrics
from analysis.ml.baseline_panel import write_standard_baseline_panel
from analysis.ml.decision_metrics import score_metric_row
from analysis.ml.feature_sets import effective_exclude_features, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.train_classifier import _design_matrix, _model
from analysis.ml.train_classifier_core import (
    _fit_model_with_optional_weights,
    combine_sample_weights,
    group_reweighting_weights,
)


def _source_mask(df: pd.DataFrame, sources: list[str], source_col: str = "label_source") -> pd.Series:
    source_text = df.get(source_col, pd.Series("", index=df.index)).fillna("").astype(str)
    mask = pd.Series(False, index=df.index)
    for source in sources:
        token = str(source).strip()
        if token:
            mask |= source_text.eq(token) | source_text.str.split(";").map(lambda parts: token in parts)
    return mask


def _stratified_half_split(df: pd.DataFrame, label_col: str, seed: int) -> tuple[pd.Index, pd.Index]:
    cal_parts = []
    test_parts = []
    for _label, group in df.groupby(label_col):
        shuffled = group.sample(frac=1.0, random_state=seed)
        cut = max(1, len(shuffled) // 2)
        if len(shuffled) - cut == 0 and len(shuffled) > 1:
            cut = len(shuffled) - 1
        cal_parts.append(shuffled.index[:cut])
        test_parts.append(shuffled.index[cut:])
    cal_idx = cal_parts[0].append(cal_parts[1:]) if cal_parts else pd.Index([])
    test_idx = test_parts[0].append(test_parts[1:]) if test_parts else pd.Index([])
    return cal_idx, test_idx


def _fit_calibrator(method: str, scores: pd.Series, labels: pd.Series) -> Any:
    if method == "logistic":
        from sklearn.linear_model import LogisticRegression

        model = LogisticRegression(max_iter=1000)
        model.fit(scores.to_numpy().reshape(-1, 1), labels.astype(int).to_numpy())
        return model
    if method == "isotonic":
        from sklearn.isotonic import IsotonicRegression

        model = IsotonicRegression(out_of_bounds="clip")
        model.fit(scores.to_numpy(), labels.astype(int).to_numpy())
        return model
    raise ValueError(f"unsupported calibration method: {method}")


def _apply_calibrator(model: Any, method: str, scores: pd.Series) -> list[float]:
    if method == "logistic":
        return [float(v) for v in model.predict_proba(scores.to_numpy().reshape(-1, 1))[:, 1]]
    return [float(v) for v in model.predict(scores.to_numpy())]


def _grouped_metrics(pred: pd.DataFrame, label_col: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group_col in ["label_source", "source_family", "scaffold_key", "ligand_chemotype", "target_family", "protein_class"]:
        if group_col not in pred.columns:
            continue
        for value, group in pred.groupby(group_col, dropna=False):
            for score_col, method in [
                ("raw_prediction_score", "raw_transfer"),
                ("source_calibrated_prediction_score", "source_calibrated"),
            ]:
                if score_col not in group.columns:
                    continue
                row = score_metric_row(group, label_col=label_col, score_col=score_col, method=method)
                row["group_col"] = group_col
                row["group_value"] = value
                rows.append(row)
    return pd.DataFrame(rows)


def run_source_calibration(
    dataset_path: str | Path,
    label_col: str,
    feature_set: str,
    model_type: str,
    train_sources: list[str],
    calibration_source: str,
    out_dir: str | Path,
    *,
    calibration_method: str = "logistic",
    seed: int = 42,
    class_weight: str | None = "balanced",
    group_reweight_cols: list[str] | None = None,
    group_reweight_include_label: bool = True,
    group_reweight_max_factor: float = 5.0,
    source_col: str = "label_source",
) -> dict[str, object]:
    df = pd.read_csv(dataset_path, low_memory=False)
    excluded = effective_exclude_features(label_col, None)
    features = [feature for feature in get_feature_set(feature_set) if feature in df.columns and feature not in excluded]
    assert_no_leakage(features)
    data = df.dropna(subset=[label_col]).copy()
    data[label_col] = binary_label_series(data[label_col])
    data = data.dropna(subset=[label_col]).copy()
    data[label_col] = data[label_col].astype(int)
    train = data.loc[_source_mask(data, train_sources, source_col=source_col)].copy()
    calibration_source_rows = data.loc[
        _source_mask(data, [calibration_source], source_col=source_col) & ~data.index.isin(train.index)
    ].copy()
    if train[label_col].nunique() < 2:
        raise ValueError("source-calibration train rows do not contain both classes")
    if calibration_source_rows[label_col].nunique() < 2:
        raise ValueError("calibration source rows do not contain both classes")
    cal_idx, test_idx = _stratified_half_split(calibration_source_rows, label_col, seed)
    cal = calibration_source_rows.loc[cal_idx]
    test = calibration_source_rows.loc[test_idx]
    x_train = _design_matrix(train, features)
    model = _model(model_type, seed, class_weight=class_weight)
    group_weights, group_manifest = group_reweighting_weights(
        train,
        label_col,
        group_cols=group_reweight_cols,
        include_label=group_reweight_include_label,
        max_factor=group_reweight_max_factor,
    )
    weights = combine_sample_weights(group_weights)
    _fit_model_with_optional_weights(model, x_train, train[label_col], weights)
    x_cal = _design_matrix(cal, features).reindex(columns=x_train.columns, fill_value=0)
    x_test = _design_matrix(test, features).reindex(columns=x_train.columns, fill_value=0)
    cal_raw = pd.Series(model.predict_proba(x_cal)[:, 1], index=cal.index)
    test_raw = pd.Series(model.predict_proba(x_test)[:, 1], index=test.index)
    calibrator = _fit_calibrator(calibration_method, cal_raw, cal[label_col])
    test_calibrated = _apply_calibrator(calibrator, calibration_method, test_raw)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    pred = test[
        [
            col
            for col in [
                "drug_id",
                "target_id",
                "pdb_id",
                "label_source",
                "source_family",
                "external_source_family",
                "external_upstream_source",
                "scaffold_key",
                "ligand_chemotype",
                "target_family",
                "protein_class",
                "assay_type",
                "endpoint_type",
                "activity_type",
                label_col,
            ]
            if col in test.columns
        ]
    ].copy()
    pred["raw_prediction_score"] = test_raw.values
    pred["source_calibrated_prediction_score"] = test_calibrated
    pred.to_csv(out_path / "source_calibrated_predictions.csv", index=False)
    (out_path / "group_reweighting_manifest.json").write_text(
        json.dumps(group_manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    if weights is not None:
        pd.DataFrame({"row_index": [str(idx) for idx in weights.index], "sample_weight": weights.astype(float).tolist()}).to_csv(
            out_path / "group_reweighting_weights.csv",
            index=False,
        )
    baseline_frame = test.copy()
    baseline_frame["ml_prediction_score"] = test_calibrated
    write_standard_baseline_panel(baseline_frame, label_col=label_col, out_path=out_path / "standard_baseline_panel.csv", seed=seed)
    _grouped_metrics(pred, label_col).to_csv(out_path / "source_calibration_grouped_metrics.csv", index=False)
    metric_rows = []
    raw_metrics = calibration_metrics(test_raw.tolist(), test[label_col].tolist())
    calibrated_metrics = calibration_metrics(test_calibrated, test[label_col].tolist())
    for metric, value in raw_metrics.items():
        metric_rows.append({"method": "raw_transfer", "metric": metric, "value": value})
    for metric, value in calibrated_metrics.items():
        metric_rows.append({"method": f"{calibration_source}_{calibration_method}_calibrated", "metric": metric, "value": value})
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(out_path / "source_calibration_metrics.csv", index=False)
    balance = pd.DataFrame(
        [
            {"split": "train", "source": ";".join(train_sources), "n": len(train), "positive_rate": float(train[label_col].mean())},
            {"split": "calibration", "source": calibration_source, "n": len(cal), "positive_rate": float(cal[label_col].mean())},
            {"split": "test", "source": calibration_source, "n": len(test), "positive_rate": float(test[label_col].mean())},
        ]
    )
    balance.to_csv(out_path / "source_calibration_label_balance.csv", index=False)
    manifest = {
        "dataset_path": str(dataset_path),
        "feature_set": feature_set,
        "model_type": model_type,
        "source_col": source_col,
        "train_sources": train_sources,
        "calibration_source": calibration_source,
        "calibration_method": calibration_method,
        "group_reweighting": group_manifest,
        "interpretation": "Uses a held-out source calibration split to assess whether source-specific calibration improves transfer; not a pure no-target-source external validation.",
    }
    (out_path / "source_calibration_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "metrics": metric_rows,
        "n_train": len(train),
        "n_calibration": len(cal),
        "n_test": len(test),
        "group_reweighting": group_manifest,
    }
