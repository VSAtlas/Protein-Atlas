from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Any

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold

try:  # pragma: no cover - optional sklearn versions
    from sklearn.model_selection import StratifiedGroupKFold
except Exception:  # pragma: no cover - optional sklearn versions
    StratifiedGroupKFold = None  # type: ignore[assignment]


def _normalize_groups(groups: np.ndarray | list[str] | tuple[str, ...]) -> np.ndarray:
    values = np.asarray(groups, dtype=object).reshape(-1)
    normalized: list[str] = []
    for idx, raw in enumerate(values):
        text = str(raw or "").strip()
        normalized.append(text if text else f"__empty_scaffold_{idx}")
    return np.asarray(normalized, dtype=object)


def _clip_probabilities(p: np.ndarray, *, eps: float = 1e-8) -> np.ndarray:
    return np.clip(np.asarray(p, dtype=float), eps, 1.0 - eps)


def _logit_probabilities(p: np.ndarray) -> np.ndarray:
    p_used = _clip_probabilities(p)
    return np.log(p_used / (1.0 - p_used))


def _build_grouped_splits(
    *,
    y: np.ndarray,
    groups: np.ndarray,
    cv_folds: int,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    n_rows = int(y.shape[0])
    n_groups = int(np.unique(groups).shape[0])
    n_splits = min(max(2, int(cv_folds)), n_rows, n_groups)
    if n_splits < 2:
        return []

    if StratifiedGroupKFold is not None:
        splitter = StratifiedGroupKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=int(seed),
        )
        return list(splitter.split(np.zeros((n_rows, 1)), y, groups))

    splitter = GroupKFold(n_splits=n_splits)
    return list(splitter.split(np.zeros((n_rows, 1)), y, groups))


@dataclass
class ProbabilityCalibrator:
    method: str
    model: Any | None = None
    fitted: bool = False
    eps: float = 1e-8
    debug_fold_indices: list[dict[str, Any]] = field(default_factory=list)

    def transform(self, p: np.ndarray) -> np.ndarray:
        p_arr = np.asarray(p, dtype=float).reshape(-1)
        clipped = _clip_probabilities(p_arr, eps=self.eps)
        if not self.fitted or self.model is None:
            return clipped

        method = self.method.strip().lower()
        if method == "sigmoid":
            x = _logit_probabilities(clipped).reshape(-1, 1)
            calibrated = np.asarray(self.model.predict_proba(x)[:, 1], dtype=float)
            return _clip_probabilities(calibrated, eps=self.eps)
        if method == "isotonic":
            calibrated = np.asarray(self.model.predict(clipped), dtype=float)
            return _clip_probabilities(calibrated, eps=self.eps)
        return clipped


def fit_oof_calibrator(
    X: sparse.csr_matrix,
    y: np.ndarray,
    groups: np.ndarray | list[str] | tuple[str, ...],
    base_model_builder: Callable[[int], Any],
    method: str,
    cv_folds: int = 5,
    seed: int = 42,
    fit_query_groups: np.ndarray | list[str] | tuple[str, ...] | None = None,
    fit_sample_weights: np.ndarray | None = None,
    *,
    test_mode: bool = False,
) -> tuple[ProbabilityCalibrator, pd.DataFrame]:
    y_arr = np.asarray(y, dtype=int).reshape(-1)
    if y_arr.size == 0:
        raise ValueError("Cannot calibrate with an empty y array.")
    if X.shape[0] != y_arr.shape[0]:
        raise ValueError("X and y must have matching row counts.")

    group_arr = _normalize_groups(groups)
    if group_arr.shape[0] != y_arr.shape[0]:
        raise ValueError("groups must have the same number of rows as y.")

    oof_p = np.full(y_arr.shape[0], np.nan, dtype=float)
    oof_fold = np.full(y_arr.shape[0], -1, dtype=int)
    debug_splits: list[dict[str, Any]] = []

    splits = _build_grouped_splits(y=y_arr, groups=group_arr, cv_folds=cv_folds, seed=seed)
    for fold_id, (train_idx, val_idx) in enumerate(splits):
        y_train = y_arr[train_idx]
        y_val = y_arr[val_idx]
        if np.unique(y_train).size < 2 or np.unique(y_val).size < 2:
            continue
        model = base_model_builder(int(seed) + fold_id)
        fit_kwargs: dict[str, Any] = {}
        if fit_query_groups is not None:
            q = np.asarray(fit_query_groups, dtype=object).reshape(-1)
            fit_kwargs["query_groups"] = q[train_idx]
        if fit_sample_weights is not None:
            w = np.asarray(fit_sample_weights, dtype=float).reshape(-1)
            fit_kwargs["sample_weight"] = w[train_idx]
        model.fit(X[train_idx], y_train, **fit_kwargs)
        p_val = np.asarray(model.predict_proba(X[val_idx])[:, 1], dtype=float)
        oof_p[val_idx] = p_val
        oof_fold[val_idx] = int(fold_id)

        if test_mode:
            debug_splits.append(
                {
                    "fold_id": int(fold_id),
                    "train_idx": np.asarray(train_idx, dtype=int).tolist(),
                    "val_idx": np.asarray(val_idx, dtype=int).tolist(),
                }
            )

    valid_mask = np.isfinite(oof_p) & (oof_fold >= 0)
    method_norm = str(method).strip().lower()
    calibrator = ProbabilityCalibrator(method=method_norm)

    if valid_mask.any() and np.unique(y_arr[valid_mask]).size >= 2:
        if method_norm == "sigmoid":
            x_fit = _logit_probabilities(oof_p[valid_mask]).reshape(-1, 1)
            lr = LogisticRegression(
                solver="lbfgs",
                max_iter=1000,
                random_state=int(seed),
            )
            lr.fit(x_fit, y_arr[valid_mask])
            calibrator.model = lr
            calibrator.fitted = True
        elif method_norm == "isotonic":
            iso = IsotonicRegression(y_min=0.0, y_max=1.0, out_of_bounds="clip")
            iso.fit(oof_p[valid_mask], y_arr[valid_mask])
            calibrator.model = iso
            calibrator.fitted = True

    if test_mode:
        calibrator.debug_fold_indices = debug_splits

    oof_df = pd.DataFrame(
        {
            "row_id": np.arange(y_arr.shape[0], dtype=int),
            "fold_id": oof_fold.astype(int),
            "group": group_arr.astype(str),
            "y_true": y_arr.astype(int),
            "p_uncalibrated": oof_p.astype(float),
            "used_for_calibration": valid_mask.astype(int),
        }
    )
    return calibrator, oof_df


def apply_calibration(calibrator: ProbabilityCalibrator, p: np.ndarray) -> np.ndarray:
    return calibrator.transform(np.asarray(p, dtype=float))
