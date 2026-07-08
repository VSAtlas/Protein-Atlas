from __future__ import annotations

import random
from pathlib import Path

import pandas as pd

from analysis.calibration.metrics import calibration_metrics


def _bootstrap_metric_ci(
    pred: pd.DataFrame,
    label_col: str,
    *,
    n_bootstraps: int,
    seed: int,
) -> pd.DataFrame:
    if n_bootstraps <= 0 or pred.empty:
        return pd.DataFrame(columns=["metric", "ci_low", "ci_high", "n_bootstraps"])
    rng = random.Random(seed)
    rows: dict[str, list[float]] = {}
    labels = pred[label_col].astype(int).tolist()
    probs = pred["ml_prediction_score"].astype(float).tolist()
    n = len(pred)
    for _idx in range(n_bootstraps):
        sample_idx = [rng.randrange(n) for _ in range(n)]
        sample_labels = [labels[i] for i in sample_idx]
        if len(set(sample_labels)) < 2:
            continue
        sample_probs = [probs[i] for i in sample_idx]
        metrics = calibration_metrics(sample_probs, sample_labels)
        for metric, value in metrics.items():
            if pd.notna(value):
                rows.setdefault(metric, []).append(float(value))
    out = []
    for metric, values in rows.items():
        if not values:
            continue
        series = pd.Series(values)
        out.append(
            {
                "metric": metric,
                "ci_low": float(series.quantile(0.025)),
                "ci_high": float(series.quantile(0.975)),
                "n_bootstraps": int(len(values)),
            }
        )
    return pd.DataFrame(out)


def _write_reliability_outputs(pred: pd.DataFrame, label_col: str, out_path: Path) -> None:
    work = pred[[label_col, "ml_prediction_score"]].dropna().copy()
    if work.empty:
        return
    work[label_col] = work[label_col].astype(int)
    work["bin"] = pd.cut(
        work["ml_prediction_score"].astype(float),
        bins=[i / 10 for i in range(11)],
        include_lowest=True,
    )
    table = (
        work.groupby("bin", observed=False)
        .agg(
            n=(label_col, "size"),
            mean_predicted_probability=("ml_prediction_score", "mean"),
            observed_positive_rate=(label_col, "mean"),
        )
        .reset_index()
    )
    table["bin"] = table["bin"].astype(str)
    table.to_csv(out_path / "model_reliability_table.csv", index=False)
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.plot([0, 1], [0, 1], color="0.5", linestyle="--", linewidth=1)
        ax.plot(
            table["mean_predicted_probability"],
            table["observed_positive_rate"],
            marker="o",
            linewidth=1.5,
        )
        ax.set_xlabel("Mean predicted probability")
        ax.set_ylabel("Observed positive rate")
        ax.set_title("ML reliability")
        fig.tight_layout()
        fig.savefig(out_path / "ml_calibration_curve.png", dpi=160)
        plt.close(fig)
    except ImportError:
        return
    except Exception:
        return


def _write_grouped_calibration_outputs(pred: pd.DataFrame, label_col: str, out_path: Path) -> None:
    group_cols = [
        "label_source",
        "source_family",
        "external_source_family",
        "external_upstream_source",
        "external_evidence_sources",
        "mechanism_label_source",
        "scaffold_key",
        "chemical_cluster",
        "target_family",
        "protein_class",
        "assay_type",
        "endpoint_type",
        "activity_type",
    ]
    rows: list[dict[str, object]] = []
    for col in group_cols:
        if col not in pred.columns:
            continue
        for value, group in pred.groupby(col, dropna=False):
            labels = pd.to_numeric(group[label_col], errors="coerce")
            scores = pd.to_numeric(group["ml_prediction_score"], errors="coerce")
            valid = labels.notna() & scores.notna()
            metric_values = (
                calibration_metrics(scores.loc[valid].tolist(), labels.loc[valid].astype(int).tolist())
                if valid.sum() and labels.loc[valid].nunique() > 1
                else {}
            )
            rows.append(
                {
                    "group_col": col,
                    "group_value": value,
                    "n": int(len(group)),
                    "positive_rate": float(labels.mean()) if labels.notna().any() else None,
                    **metric_values,
                }
            )
    pd.DataFrame(rows).to_csv(out_path / "model_grouped_calibration.csv", index=False)


def write_subgroup_metrics(pred: pd.DataFrame, label_col: str, out_path: Path) -> pd.DataFrame:
    group_cols = [
        "target_id",
        "target_family",
        "protein_class",
        "label_source",
        "source_family",
        "upstream_source",
        "external_source_family",
        "external_upstream_source",
        "external_evidence_sources",
        "mechanism_label_source",
        "scaffold_key",
        "chemical_cluster",
        "ligand_chemotype",
        "assay_type",
        "endpoint_type",
        "activity_type",
        "applicability_domain",
    ]
    rows: list[dict[str, object]] = []
    for col in group_cols:
        if col not in pred.columns:
            continue
        for value, group in pred.groupby(col, dropna=False):
            labels = pd.to_numeric(group[label_col], errors="coerce")
            scores = pd.to_numeric(group["ml_prediction_score"], errors="coerce")
            valid = labels.notna() & scores.notna()
            metric_values = (
                calibration_metrics(scores.loc[valid].tolist(), labels.loc[valid].astype(int).tolist())
                if valid.sum() and labels.loc[valid].nunique() > 1
                else {}
            )
            rows.append(
                {
                    "group_col": col,
                    "group_value": value,
                    "n": int(len(group)),
                    "n_positive": int((labels == 1).sum()),
                    "n_negative": int((labels == 0).sum()),
                    "positive_rate": float(labels.mean()) if labels.notna().any() else None,
                    "mean_prediction": float(scores.mean()) if scores.notna().any() else None,
                    "median_prediction": float(scores.median()) if scores.notna().any() else None,
                    **metric_values,
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(out_path / "model_subgroup_metrics.csv", index=False)
    return frame


def _write_source_transfer_diagnostics(
    train: pd.DataFrame,
    test: pd.DataFrame,
    pred: pd.DataFrame,
    label_col: str,
    split_column: str | None,
    out_path: Path,
) -> None:
    if not split_column or split_column not in train.columns or split_column not in test.columns:
        return
    rows = []
    for split_name, frame in (("train", train), ("test", test)):
        grouped = frame.groupby(split_column, dropna=False)
        for source, group in grouped:
            rows.append(
                {
                    "split": split_name,
                    "source": source,
                    "n": int(len(group)),
                    "positive_rate": float(group[label_col].mean()) if len(group) else float("nan"),
                }
            )
    pd.DataFrame(rows).to_csv(out_path / "source_transfer_label_balance.csv", index=False)
    if split_column in pred.columns:
        score_rows = []
        grouped_pred = pred.groupby(split_column, dropna=False)
        for source, group in grouped_pred:
            labels = pd.to_numeric(group[label_col], errors="coerce")
            scores = pd.to_numeric(group["ml_prediction_score"], errors="coerce")
            valid = labels.notna() & scores.notna()
            metric_values = (
                calibration_metrics(scores.loc[valid].tolist(), labels.loc[valid].astype(int).tolist())
                if valid.sum() and labels.loc[valid].nunique() > 1
                else {}
            )
            score_rows.append(
                {
                    "source": source,
                    "n": int(len(group)),
                    "positive_rate": float(labels.mean()) if labels.notna().any() else float("nan"),
                    "mean_prediction": float(group["ml_prediction_score"].mean()),
                    "median_prediction": float(group["ml_prediction_score"].median()),
                    **metric_values,
                }
            )
        pd.DataFrame(score_rows).to_csv(out_path / "source_transfer_prediction_balance.csv", index=False)


def _append_applicability_domain(
    pred: pd.DataFrame,
    x_train: pd.DataFrame,
    x_test: pd.DataFrame,
    label_col: str,
    out_path: Path,
) -> pd.DataFrame:
    if x_train.empty or x_test.empty:
        return pred
    try:
        from sklearn.neighbors import NearestNeighbors
    except ImportError:
        return pred

    means = x_train.mean(axis=0)
    stds = x_train.std(axis=0).replace(0, 1.0).fillna(1.0)
    train_scaled = ((x_train - means) / stds).fillna(0.0)
    test_scaled = ((x_test - means) / stds).fillna(0.0)
    n_neighbors = 2 if len(train_scaled) > 1 else 1
    train_nn = NearestNeighbors(n_neighbors=n_neighbors)
    train_nn.fit(train_scaled)
    train_distances = train_nn.kneighbors(train_scaled, return_distance=True)[0]
    train_reference_distance = train_distances[:, -1]
    threshold = float(pd.Series(train_reference_distance).quantile(0.95))
    test_nn = NearestNeighbors(n_neighbors=1)
    test_nn.fit(train_scaled)
    test_distance = test_nn.kneighbors(test_scaled, return_distance=True)[0][:, 0]
    out = pred.copy()
    out["applicability_distance"] = test_distance
    out["applicability_domain"] = [
        "in_domain" if float(distance) <= threshold else "out_of_domain" for distance in test_distance
    ]
    distance_by_row = pd.Series(test_distance, index=out.index)

    rows: list[dict[str, object]] = []
    for domain, group in out.groupby("applicability_domain", dropna=False):
        labels = pd.to_numeric(group[label_col], errors="coerce")
        scores = pd.to_numeric(group["ml_prediction_score"], errors="coerce")
        valid = labels.notna() & scores.notna()
        metric_values = (
            calibration_metrics(scores.loc[valid].tolist(), labels.loc[valid].astype(int).tolist())
            if valid.sum() and labels.loc[valid].nunique() > 1
            else {}
        )
        rows.append(
            {
                "applicability_domain": domain,
                "n": int(len(group)),
                "positive_rate": float(labels.mean()) if labels.notna().any() else None,
                "distance_threshold_train_q95": threshold,
                "median_distance": float(distance_by_row.loc[group.index].median()) if len(group) else None,
                **metric_values,
            }
        )
    pd.DataFrame(rows).to_csv(out_path / "model_applicability_domain_summary.csv", index=False)
    return out
