from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from analysis.statistics import auroc, average_precision


def _read_profile(path: str | Path) -> pd.DataFrame:
    table_path = Path(path)
    if table_path.suffix.lower() == ".parquet":
        return pd.read_parquet(table_path)
    return pd.read_csv(table_path, index_col=0)


def _read_labels(path: str | Path) -> pd.DataFrame:
    labels = pd.read_csv(path)
    if {"drug_id", "side_effect", "label"}.issubset(labels.columns):
        return labels.pivot_table(index="drug_id", columns="side_effect", values="label", aggfunc="max").fillna(0)
    if "drug_id" in labels.columns:
        return labels.set_index("drug_id")
    return pd.read_csv(path, index_col=0)


def _split_drugs(drugs: list[str], seed: int) -> tuple[list[str], list[str]]:
    import random

    rng = random.Random(seed)
    shuffled = list(drugs)
    rng.shuffle(shuffled)
    cut = max(1, int(len(shuffled) * 0.8)) if len(shuffled) > 1 else len(shuffled)
    return shuffled[:cut], shuffled[cut:]


def run_sawada_side_effect_baseline(
    profile_path: str | Path,
    labels_path: str | Path,
    out_dir: str | Path,
    *,
    model_type: str = "logistic_regression",
    split_mode: str = "drug_holdout",
    seed: int = 42,
) -> dict[str, Any]:
    from sklearn.ensemble import RandomForestClassifier  # type: ignore[import-untyped]
    from sklearn.impute import SimpleImputer  # type: ignore[import-untyped]
    from sklearn.linear_model import LogisticRegression  # type: ignore[import-untyped]
    from sklearn.metrics import brier_score_loss  # type: ignore[import-untyped]
    from sklearn.pipeline import make_pipeline  # type: ignore[import-untyped]

    profile = _read_profile(profile_path)
    labels = _read_labels(labels_path)
    common = sorted(set(profile.index.astype(str)) & set(labels.index.astype(str)))
    profile = profile.loc[common]
    labels = labels.loc[common].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    train_drugs, test_drugs = _split_drugs(common, seed)
    if split_mode != "drug_holdout":
        raise ValueError("Sawada-mode baseline currently supports split_mode='drug_holdout'")
    x_train = profile.loc[train_drugs]
    x_test = profile.loc[test_drugs] if test_drugs else profile.loc[train_drugs]
    y_train = labels.loc[train_drugs]
    y_test = labels.loc[test_drugs] if test_drugs else labels.loc[train_drugs]
    predictions: list[dict[str, Any]] = []
    per_label: list[dict[str, Any]] = []
    def model_factory() -> object:
        if model_type == "random_forest":
            return RandomForestClassifier(n_estimators=200, random_state=seed, n_jobs=1)
        return LogisticRegression(max_iter=1000, class_weight="balanced")
    for side_effect in y_train.columns:
        y_col = y_train[side_effect]
        if y_col.nunique() < 2:
            continue
        model = make_pipeline(SimpleImputer(strategy="median"), model_factory())
        model.fit(x_train, y_col)
        probs = model.predict_proba(x_test)[:, 1]
        truth = y_test[side_effect].astype(int).tolist()
        for drug_id, prob, label in zip(x_test.index, probs, truth):
            predictions.append({"drug_id": drug_id, "side_effect": side_effect, "prediction": prob, "label": label})
        per_label.append(
            {
                "side_effect": side_effect,
                "n_test": len(truth),
                "positives": sum(truth),
                "AUPRC": average_precision(probs.tolist(), truth),
                "AUROC": auroc(probs.tolist(), truth),
                "Brier": brier_score_loss(truth, probs) if len(set(truth)) > 1 else "",
            }
        )
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pred_df = pd.DataFrame(predictions)
    metric_df = pd.DataFrame(per_label)
    pred_df.to_csv(out / "predictions.csv", index=False)
    metric_df.to_csv(out / "per_side_effect_metrics.csv", index=False)
    summary = {
        "model_type": model_type,
        "split_mode": split_mode,
        "n_drugs": len(common),
        "n_train_drugs": len(train_drugs),
        "n_test_drugs": len(test_drugs) if test_drugs else len(train_drugs),
        "macro_AUPRC": float(pd.to_numeric(metric_df.get("AUPRC", pd.Series(dtype=float)), errors="coerce").mean()) if not metric_df.empty else 0.0,
        "macro_AUROC": float(pd.to_numeric(metric_df.get("AUROC", pd.Series(dtype=float)), errors="coerce").mean()) if not metric_df.empty else 0.0,
    }
    pd.DataFrame([summary]).to_csv(out / "metrics.csv", index=False)
    return summary
