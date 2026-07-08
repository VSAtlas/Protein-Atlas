from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.ml.labels import binary_label_series
from analysis.ml.splits import make_split


DEFAULT_FORBIDDEN_SELECTION_FEATURES = {
    "atlas_score",
    "SCORCH_score_used",
    "final_score",
    "consensus_score",
    "z_score",
    "z_selected",
    "fdr_q_value",
    "mechanism_graph_score",
    "heatmap_rank",
}


def _distance_from_positive_centroid(frame: pd.DataFrame, features: list[str], label_col: str) -> pd.Series:
    labels = binary_label_series(frame[label_col])
    positives = frame.loc[labels.eq(1), features].apply(pd.to_numeric, errors="coerce")
    candidates = frame.loc[labels.isna(), features].apply(pd.to_numeric, errors="coerce")
    if positives.empty or candidates.empty:
        return pd.Series(dtype=float)
    center = positives.median(axis=0, skipna=True)
    scale = positives.subtract(center).abs().median(axis=0, skipna=True).replace(0, 1.0).fillna(1.0)
    distances = candidates.subtract(center).divide(scale).pow(2).sum(axis=1).pow(0.5)
    return distances.sort_values(ascending=False)


def select_fold_reliable_negatives(
    df: pd.DataFrame,
    *,
    label_col: str,
    selection_features: list[str],
    split_mode: str = "drug_holdout",
    seed: int = 42,
    test_fraction: float = 0.2,
    n_per_positive: float = 1.0,
    max_negatives: int | None = None,
    fold_name: str = "train",
) -> pd.DataFrame:
    """Mine reliable negatives inside the training fold using external features only."""

    forbidden = DEFAULT_FORBIDDEN_SELECTION_FEATURES & set(selection_features)
    if forbidden:
        raise ValueError(
            "reliable-negative selection features cannot include Atlas/model score fields: "
            + ", ".join(sorted(forbidden))
        )
    missing = [feature for feature in selection_features if feature not in df.columns]
    if missing:
        raise ValueError(f"reliable-negative selection features not found: {missing}")
    train_idx, _ = make_split(df, split_mode=split_mode, seed=seed, test_fraction=test_fraction)
    train = df.loc[train_idx].copy()
    labels = binary_label_series(train[label_col])
    n_pos = int(labels.eq(1).sum())
    if n_pos == 0:
        return pd.DataFrame(columns=[*df.columns, "negative_evidence_type"])
    distances = _distance_from_positive_centroid(train, selection_features, label_col)
    n_select = int(round(n_pos * float(n_per_positive)))
    if max_negatives is not None:
        n_select = min(n_select, int(max_negatives))
    selected_idx = distances.head(max(0, n_select)).index
    out = train.loc[selected_idx].copy()
    out[label_col] = 0
    out["mechanism_label"] = 0
    out["negative_evidence_type"] = "fold_specific_reliable_negative"
    out["negative_source"] = "external_similarity_mining"
    out["negative_confidence"] = 0.35
    out["negative_selection_fold"] = fold_name
    out["negative_selection_features"] = ";".join(selection_features)
    out["reliable_negative_distance"] = distances.loc[selected_idx]
    return out


def build_reliable_negative_table(
    dataset_path: str | Path,
    out_path: str | Path,
    *,
    label_col: str,
    selection_features: list[str],
    split_mode: str = "drug_holdout",
    seed: int = 42,
    test_fraction: float = 0.2,
    n_per_positive: float = 1.0,
    max_negatives: int | None = None,
) -> pd.DataFrame:
    df = pd.read_csv(dataset_path, low_memory=False)
    out = select_fold_reliable_negatives(
        df,
        label_col=label_col,
        selection_features=selection_features,
        split_mode=split_mode,
        seed=seed,
        test_fraction=test_fraction,
        n_per_positive=n_per_positive,
        max_negatives=max_negatives,
    )
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    return out

