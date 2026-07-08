from __future__ import annotations

import pickle
from pathlib import Path
from typing import Sequence

import pandas as pd

from analysis.calibration.metrics import calibration_metrics
from analysis.ml.splits import make_split


def _fit_predict(method: str, train_x: pd.Series, train_y: pd.Series, test_x: pd.Series) -> tuple[object, list[float]]:
    from sklearn.isotonic import IsotonicRegression
    from sklearn.linear_model import LogisticRegression

    if method == "isotonic":
        model = IsotonicRegression(out_of_bounds="clip")
        model.fit(train_x.to_numpy(), train_y.to_numpy())
        return model, [float(v) for v in model.predict(test_x.to_numpy())]
    if method == "logistic":
        model = LogisticRegression(max_iter=1000)
        model.fit(train_x.to_numpy().reshape(-1, 1), train_y.to_numpy())
        return model, [float(v) for v in model.predict_proba(test_x.to_numpy().reshape(-1, 1))[:, 1]]
    raise ValueError(f"unsupported calibration method: {method}")


def run_atlas_calibration(
    benchmark_table_path: str | Path,
    label_col: str,
    score_col: str = "atlas_score",
    out_dir: str | Path = "outputs/calibration",
    methods: Sequence[str] = ("logistic", "isotonic"),
    split_mode: str = "random",
    seed: int = 42,
) -> dict[str, object]:
    df = pd.read_csv(benchmark_table_path)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    data = df.dropna(subset=[score_col, label_col]).copy()
    data[label_col] = data[label_col].astype(bool).astype(int)
    train_idx, test_idx = make_split(data, split_mode=split_mode, seed=seed)
    train = data.loc[train_idx]
    test = data.loc[test_idx]
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    predictions: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []
    fitted: dict[str, object] = {}
    for method in methods:
        model, probs = _fit_predict(method, train[score_col], train[label_col], test[score_col])
        fitted[method] = model
        pred = test[["drug_id", "target_id", score_col, label_col]].copy() if {"drug_id", "target_id"}.issubset(test.columns) else test[[score_col, label_col]].copy()
        pred["calibration_method"] = method
        pred["calibrated_activity_probability"] = probs
        predictions.append(pred)
        metrics = calibration_metrics(probs, test[label_col].tolist(), test[score_col].tolist())
        metric_rows.extend({"method": method, "metric": key, "value": value} for key, value in metrics.items())
    pred_df = pd.concat(predictions, ignore_index=True) if predictions else pd.DataFrame()
    pred_df.to_csv(out_path / "calibration_predictions.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(out_path / "calibration_metrics.csv", index=False)
    with (out_path / "calibration_model.pkl").open("wb") as handle:
        pickle.dump(fitted, handle)
    write_calibration_plots(pred_df, label_col, out_path.parent / "figures")
    return {"metrics": metric_rows, "n_train": len(train), "n_test": len(test)}


def write_calibration_plots(predictions: pd.DataFrame, label_col: str, figures_dir: str | Path) -> None:
    if predictions.empty:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(figures_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 4), dpi=160)
    for method, group in predictions.groupby("calibration_method"):
        bins = pd.cut(group["calibrated_activity_probability"], bins=10, include_lowest=True)
        rel = group.groupby(bins, observed=False).agg(prob=("calibrated_activity_probability", "mean"), obs=(label_col, "mean")).dropna()
        ax.plot(rel["prob"], rel["obs"], marker="o", label=method)
    ax.plot([0, 1], [0, 1], color="black", linewidth=1)
    ax.set_xlabel("Predicted probability")
    ax.set_ylabel("Observed fraction")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "calibration_curve.png")
    fig.savefig(out_dir / "reliability_diagram.png")
    plt.close(fig)

