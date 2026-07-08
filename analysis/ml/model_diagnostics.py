from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.feature_sets import effective_exclude_features, get_feature_set
from analysis.ml.labels import binary_label_series
from analysis.ml.leakage_checks import assert_no_leakage
from analysis.ml.splits import make_split, split_overlap_summary
from analysis.ml.train_classifier import _design_matrix, _model


def _metrics(y_true: pd.Series, y_score: pd.Series) -> dict[str, float]:
    from sklearn.metrics import (
        average_precision_score,
        brier_score_loss,
        log_loss,
        roc_auc_score,
    )

    labels = y_true.astype(int)
    scores = y_score.astype(float)
    out: dict[str, float] = {
        "log_loss": float(log_loss(labels, scores.clip(1e-6, 1 - 1e-6), labels=[0, 1])),
        "brier": float(brier_score_loss(labels, scores)),
        "auprc": float(average_precision_score(labels, scores)),
    }
    out["auroc"] = float(roc_auc_score(labels, scores)) if labels.nunique() == 2 else float("nan")
    return out


def _plot_learning_curve(metrics: pd.DataFrame, out_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(metrics["train_fraction"], metrics["train_log_loss"], marker="o", label="train")
        ax.plot(metrics["train_fraction"], metrics["test_log_loss"], marker="o", label="test")
        ax.set_xlabel("Training fraction")
        ax.set_ylabel("Log loss")
        ax.set_title("Learning curve")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / "learning_curve_log_loss.png", dpi=160)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(metrics["train_fraction"], metrics["test_auprc"], marker="o", label="test AUPRC")
        ax.plot(metrics["train_fraction"], metrics["test_auroc"], marker="o", label="test AUROC")
        ax.set_xlabel("Training fraction")
        ax.set_ylabel("Metric")
        ax.set_title("Held-out performance by training size")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / "learning_curve_auc.png", dpi=160)
        plt.close(fig)
    except Exception:
        return


def run_ml_diagnostics(
    dataset_path: str | Path,
    label_col: str,
    feature_set: str,
    model_type: str,
    split_mode: str,
    out_dir: str | Path,
    *,
    seed: int = 42,
    exclude_features: list[str] | None = None,
    train_fractions: list[float] | None = None,
    class_weight: str | None = "balanced",
) -> dict[str, Any]:
    df = pd.read_csv(dataset_path, low_memory=False)
    excluded = effective_exclude_features(label_col, exclude_features)
    features = [feature for feature in get_feature_set(feature_set) if feature in df.columns and feature not in excluded]
    assert_no_leakage(features)
    data = df.dropna(subset=[label_col]).copy()
    data[label_col] = binary_label_series(data[label_col])
    data = data.dropna(subset=[label_col]).copy()
    data[label_col] = data[label_col].astype(int)
    train_idx, test_idx = make_split(data, split_mode=split_mode, seed=seed)
    train = data.loc[train_idx].copy()
    test = data.loc[test_idx].copy()
    split_summary = split_overlap_summary(train, test, split_mode)
    if not split_summary.get("passes_holdout", False):
        raise ValueError(f"holdout split leakage detected: {split_summary}")
    fractions = train_fractions or [0.1, 0.25, 0.5, 0.75, 1.0]
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    learning_rows: list[dict[str, Any]] = []
    for fraction in fractions:
        subset = (
            train.groupby(label_col, group_keys=False)
            .sample(frac=min(1.0, max(float(fraction), 0.0)), random_state=seed)
            .copy()
        )
        if subset[label_col].nunique() < 2:
            continue
        x_train = _design_matrix(subset, features)
        x_test = _design_matrix(test, features).reindex(columns=x_train.columns, fill_value=0)
        model = _model(model_type, seed, class_weight=class_weight)
        model.fit(x_train, subset[label_col])
        train_score = pd.Series(model.predict_proba(x_train)[:, 1], index=subset.index)
        test_score = pd.Series(model.predict_proba(x_test)[:, 1], index=test.index)
        train_metrics = _metrics(subset[label_col], train_score)
        test_metrics = _metrics(test[label_col], test_score)
        learning_rows.append(
            {
                "train_fraction": float(fraction),
                "n_train": int(len(subset)),
                "n_test": int(len(test)),
                **{f"train_{key}": val for key, val in train_metrics.items()},
                **{f"test_{key}": val for key, val in test_metrics.items()},
            }
        )
    learning = pd.DataFrame(learning_rows)
    learning.to_csv(out_path / "learning_curve_metrics.csv", index=False)
    _plot_learning_curve(learning, out_path)

    ablation_rows: list[dict[str, Any]] = []
    full_x_train = _design_matrix(train, features)
    full_x_test = _design_matrix(test, features).reindex(columns=full_x_train.columns, fill_value=0)
    full_model = _model(model_type, seed, class_weight=class_weight)
    full_model.fit(full_x_train, train[label_col])
    full_score = pd.Series(full_model.predict_proba(full_x_test)[:, 1], index=test.index)
    full_metrics = _metrics(test[label_col], full_score)
    for feature in features:
        reduced = [col for col in features if col != feature]
        if not reduced:
            continue
        x_train = _design_matrix(train, reduced)
        x_test = _design_matrix(test, reduced).reindex(columns=x_train.columns, fill_value=0)
        model = _model(model_type, seed, class_weight=class_weight)
        model.fit(x_train, train[label_col])
        score = pd.Series(model.predict_proba(x_test)[:, 1], index=test.index)
        metrics = _metrics(test[label_col], score)
        ablation_rows.append(
            {
                "removed_feature": feature,
                **metrics,
                "delta_auprc_vs_full": metrics["auprc"] - full_metrics["auprc"],
                "delta_auroc_vs_full": metrics["auroc"] - full_metrics["auroc"],
                "delta_log_loss_vs_full": metrics["log_loss"] - full_metrics["log_loss"],
            }
        )
    pd.DataFrame(ablation_rows).to_csv(out_path / "feature_ablation_metrics.csv", index=False)

    baseline_rows: list[dict[str, Any]] = []
    for feature in features:
        vals = pd.to_numeric(test[feature], errors="coerce")
        if vals.notna().sum() < 2:
            continue
        lo = vals.min()
        hi = vals.max()
        if pd.isna(lo) or pd.isna(hi) or lo == hi:
            continue
        score = ((vals - lo) / (hi - lo)).fillna(0.5)
        baseline_rows.append({"score_feature": feature, **_metrics(test[label_col], score)})
    pd.DataFrame(baseline_rows).to_csv(out_path / "score_baseline_metrics.csv", index=False)
    manifest = {
        "dataset": str(dataset_path),
        "label_col": label_col,
        "feature_set": feature_set,
        "features": features,
        "split_summary": split_summary,
        "full_model_metrics": full_metrics,
        "outputs": [
            "learning_curve_metrics.csv",
            "learning_curve_log_loss.png",
            "learning_curve_auc.png",
            "feature_ablation_metrics.csv",
            "score_baseline_metrics.csv",
        ],
    }
    (out_path / "ml_diagnostics_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
