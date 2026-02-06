from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


def _target_to_family_id(target: str) -> str:
    text = str(target or "").strip()
    if not text:
        return "unknown"
    if "_" in text:
        return text.split("_", 1)[0]
    return text


def cluster_targets_by_sequence(
    df: pd.DataFrame,
    *,
    method: str = "heuristic",
    identity_threshold: float = 0.3,
    seed: int = 42,
) -> dict[str, str]:
    del identity_threshold, seed
    if method not in {"heuristic", "mmseqs", "blast"}:
        raise ValueError("method must be one of heuristic, mmseqs, blast")
    target_col = "ex_rec_pdb" if "ex_rec_pdb" in df.columns else "target"
    targets = sorted(set(df.get(target_col, pd.Series(dtype=object)).astype(str).tolist()))
    mapping: dict[str, str] = {}
    for target in targets:
        if method == "heuristic":
            mapping[target] = _target_to_family_id(target)
        else:
            # Placeholder until mmseqs/blast integration is added.
            mapping[target] = _target_to_family_id(target)
    return mapping


def assign_family_split_labels(
    df: pd.DataFrame,
    *,
    method: str = "heuristic",
    identity_threshold: float = 0.3,
    seed: int = 42,
) -> pd.Series:
    mapping = cluster_targets_by_sequence(
        df,
        method=method,
        identity_threshold=identity_threshold,
        seed=seed,
    )
    values = df.get("ex_rec_pdb", pd.Series(dtype=object)).astype(str)
    return values.map(mapping).fillna("unknown")


def assign_pocket_similarity_clusters(
    df: pd.DataFrame,
    *,
    feature_columns: list[str],
    n_clusters: int = 5,
    seed: int = 42,
) -> pd.Series:
    if df.empty:
        return pd.Series([], dtype=int, index=df.index)
    cols = [col for col in feature_columns if col in df.columns]
    if not cols:
        return pd.Series(np.zeros(len(df), dtype=int), index=df.index)

    matrix = df.loc[:, cols].astype(float).to_numpy()
    matrix[~np.isfinite(matrix)] = 0.0
    n_used = max(2, min(int(n_clusters), len(df)))
    if len(df) < 2:
        return pd.Series(np.zeros(len(df), dtype=int), index=df.index)
    model = KMeans(n_clusters=n_used, random_state=int(seed), n_init=10)
    labels = model.fit_predict(matrix)
    return pd.Series(labels.astype(int), index=df.index)
