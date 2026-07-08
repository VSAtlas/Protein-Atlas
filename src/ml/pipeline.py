from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import sparse

from ml.models import build_model


def _group_values(df, group_key: str) -> np.ndarray:
    key = str(group_key).strip().lower()
    if key == "target":
        return df["ex_rec_pdb"].astype(str).to_numpy()
    if key == "pocket":
        return df["pocket"].astype(str).to_numpy()
    return (df["ex_rec_pdb"].astype(str) + "::" + df["pocket"].astype(str)).to_numpy()


def select_stage1_keep_mask(
    scores: np.ndarray,
    *,
    keep_top_pct: float | None = None,
    keep_prob_ge: float | None = None,
    groups: np.ndarray | None = None,
    keep_within_group: bool = True,
) -> np.ndarray:
    n = int(len(scores))
    mask = np.zeros(n, dtype=bool)
    if n == 0:
        return mask
    values = np.asarray(scores, dtype=float).reshape(-1)
    if keep_prob_ge is not None:
        mask |= values >= float(keep_prob_ge)

    if keep_top_pct is None:
        return mask if keep_prob_ge is not None else np.ones(n, dtype=bool)

    pct = float(np.clip(keep_top_pct, 0.0, 100.0))
    if pct <= 0.0:
        return mask

    if keep_within_group and groups is not None:
        group_arr = np.asarray(groups, dtype=object)
        for group in np.unique(group_arr):
            idx = np.where(group_arr == group)[0]
            if idx.size == 0:
                continue
            k = max(1, int(np.ceil((pct / 100.0) * idx.size)))
            order = np.argsort(values[idx])[::-1]
            keep_idx = idx[order[:k]]
            mask[keep_idx] = True
    else:
        k = max(1, int(np.ceil((pct / 100.0) * n)))
        order = np.argsort(values)[::-1]
        mask[order[:k]] = True
    return mask


@dataclass
class TwoStageOutputs:
    stage1_prob_holdout: np.ndarray
    stage2_rank_score_holdout: np.ndarray
    stage1_keep_holdout: np.ndarray
    stage1_keep_train: np.ndarray


def run_two_stage_pipeline(
    *,
    train_X: sparse.csr_matrix,
    train_y: np.ndarray,
    holdout_X: sparse.csr_matrix,
    train_df,
    holdout_df,
    stage1_model_family: str,
    stage1_model_params: dict[str, Any],
    stage2_model_family: str,
    stage2_model_params: dict[str, Any],
    group_key: str,
    keep_top_pct: float | None,
    keep_prob_ge: float | None,
    keep_within_group: bool,
    random_seed: int,
) -> TwoStageOutputs:
    stage1 = build_model(
        model_family=stage1_model_family,
        model_params=stage1_model_params,
        random_seed=random_seed,
    )
    stage1.fit(train_X, train_y)
    p_train = np.asarray(stage1.predict_proba(train_X)[:, 1], dtype=float)
    p_holdout = np.asarray(stage1.predict_proba(holdout_X)[:, 1], dtype=float)

    train_groups = _group_values(train_df, group_key=group_key)
    holdout_groups = _group_values(holdout_df, group_key=group_key)
    keep_train = select_stage1_keep_mask(
        p_train,
        keep_top_pct=keep_top_pct,
        keep_prob_ge=keep_prob_ge,
        groups=train_groups,
        keep_within_group=keep_within_group,
    )
    keep_holdout = select_stage1_keep_mask(
        p_holdout,
        keep_top_pct=keep_top_pct,
        keep_prob_ge=keep_prob_ge,
        groups=holdout_groups,
        keep_within_group=keep_within_group,
    )

    # If Stage 1 keeps too little, fallback to all rows.
    if int(keep_train.sum()) < 2:
        keep_train = np.ones_like(keep_train, dtype=bool)
    if int(keep_holdout.sum()) < 1:
        keep_holdout = np.ones_like(keep_holdout, dtype=bool)

    stage2 = build_model(
        model_family=stage2_model_family,
        model_params=stage2_model_params,
        random_seed=random_seed,
    )
    stage2.fit(train_X[keep_train], train_y[keep_train])
    rank_holdout_kept = np.asarray(stage2.raw_score(holdout_X[keep_holdout]), dtype=float)

    # Filtered-out rows are ranked as worst.
    min_score = float(np.min(rank_holdout_kept)) if rank_holdout_kept.size else 0.0
    worst = min_score - 1.0
    full_rank = np.full(len(holdout_df), worst, dtype=float)
    full_rank[np.where(keep_holdout)[0]] = rank_holdout_kept

    return TwoStageOutputs(
        stage1_prob_holdout=p_holdout,
        stage2_rank_score_holdout=full_rank,
        stage1_keep_holdout=keep_holdout,
        stage1_keep_train=keep_train,
    )
