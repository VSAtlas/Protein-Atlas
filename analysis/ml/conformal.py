from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd


def _binary_probabilities(model: Any, x: pd.DataFrame) -> list[float]:
    if x.empty:
        return []
    if hasattr(model, "predict_proba"):
        return [float(v) for v in model.predict_proba(x)[:, 1]]
    raw = model.decision_function(x)
    return [float(1 / (1 + math.exp(-float(v)))) for v in raw]


def _conformal_quantile(scores: list[float], alpha: float) -> float | None:
    if not scores:
        return None
    ordered = sorted(float(score) for score in scores)
    rank = min(len(ordered), math.ceil((len(ordered) + 1) * (1.0 - alpha)))
    return float(ordered[max(0, rank - 1)])


def append_split_conformal_sets(
    pred: pd.DataFrame,
    *,
    model: Any,
    x_calibration: pd.DataFrame,
    y_calibration: pd.Series,
    alpha: float = 0.10,
    label_col: str,
    out_path: str | Path,
    calibration_role: str,
) -> pd.DataFrame:
    """Append split-conformal binary prediction sets when calibration rows exist.

    The nonconformity score is 1 minus the model probability assigned to the
    observed class. This gives prediction sets with finite-sample coverage under
    the usual conformal exchangeability assumption; the manifest records when
    the calibration split also served model selection so readers can treat it as
    a caveated uncertainty layer.
    """

    out = Path(out_path)
    summary_path = out / "conformal_summary.json"
    if pred.empty or x_calibration.empty or y_calibration.empty:
        summary = {
            "status": "skipped",
            "reason": "no_calibration_rows",
            "alpha": alpha,
            "calibration_role": calibration_role,
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
        return pred

    cal_probs = _binary_probabilities(model, x_calibration)
    labels = pd.to_numeric(y_calibration, errors="coerce").dropna().astype(int)
    if len(labels) != len(cal_probs) or labels.empty:
        summary = {
            "status": "skipped",
            "reason": "invalid_calibration_labels",
            "alpha": alpha,
            "n_calibration": int(len(labels)),
            "calibration_role": calibration_role,
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
        return pred

    nonconformity = [
        (1.0 - prob) if label == 1 else prob
        for prob, label in zip(cal_probs, labels.tolist())
    ]
    qhat = _conformal_quantile(nonconformity, alpha)
    if qhat is None:
        summary = {
            "status": "skipped",
            "reason": "empty_nonconformity_scores",
            "alpha": alpha,
            "calibration_role": calibration_role,
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
        return pred

    enriched = pred.copy()
    sets: list[str] = []
    sizes: list[int] = []
    for prob in pd.to_numeric(enriched["ml_prediction_score"], errors="coerce"):
        if pd.isna(prob):
            labels_in_set: list[str] = []
        else:
            labels_in_set = []
            if float(prob) <= qhat:
                labels_in_set.append("0")
            if 1.0 - float(prob) <= qhat:
                labels_in_set.append("1")
        sets.append(";".join(labels_in_set))
        sizes.append(len(labels_in_set))
    enriched["conformal_prediction_set"] = sets
    enriched["conformal_set_size"] = sizes
    enriched["conformal_alpha"] = float(alpha)
    enriched["conformal_qhat"] = float(qhat)
    if label_col in enriched.columns:
        true_labels = pd.to_numeric(enriched[label_col], errors="coerce")
        enriched["conformal_contains_true"] = [
            (str(int(label)) in set_value.split(";")) if pd.notna(label) else None
            for label, set_value in zip(true_labels, sets)
        ]
    coverage = (
        float(pd.Series(enriched["conformal_contains_true"]).dropna().astype(bool).mean())
        if "conformal_contains_true" in enriched.columns
        and pd.Series(enriched["conformal_contains_true"]).dropna().shape[0]
        else None
    )
    summary = {
        "status": "ok",
        "alpha": float(alpha),
        "target_coverage": float(1.0 - alpha),
        "empirical_coverage": coverage,
        "average_set_size": float(pd.Series(sizes).mean()) if sizes else None,
        "n_calibration": int(len(labels)),
        "n_prediction": int(len(enriched)),
        "qhat": float(qhat),
        "calibration_role": calibration_role,
        "caveat": (
            "Use a calibration split disjoint from model selection for strict "
            "claim-grade conformal coverage."
            if calibration_role != "pure_calibration"
            else ""
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
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
    return enriched
