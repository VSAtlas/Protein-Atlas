from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import sparse
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression

try:  # pragma: no cover - optional dependency
    import pandas as pd  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - optional dependency
    pd = None  # type: ignore[assignment]

try:  # pragma: no cover - optional dependency
    from lightgbm import LGBMClassifier, LGBMRanker  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - optional dependency
    LGBMClassifier = None  # type: ignore[assignment]
    LGBMRanker = None  # type: ignore[assignment]

try:  # pragma: no cover - optional dependency
    from xgboost import XGBClassifier, XGBRanker  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - optional dependency
    XGBClassifier = None  # type: ignore[assignment]
    XGBRanker = None  # type: ignore[assignment]


def _as_model_input(
    X: sparse.csr_matrix,
    *,
    requires_dense: bool,
    with_column_names: bool = False,
) -> sparse.csr_matrix | np.ndarray | Any:
    base_input: sparse.csr_matrix | np.ndarray = X.toarray() if requires_dense else X
    if not with_column_names or pd is None:
        return base_input

    if sparse.issparse(base_input):
        n_features = int(base_input.shape[1])
        columns = [f"Column_{idx}" for idx in range(n_features)]
        return pd.DataFrame.sparse.from_spmatrix(base_input, columns=columns)

    arr = np.asarray(base_input)
    n_features = int(arr.shape[1]) if arr.ndim == 2 else 1
    columns = [f"Column_{idx}" for idx in range(n_features)]
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return pd.DataFrame(arr, columns=columns)


def _clip_probabilities(p: np.ndarray, *, eps: float = 1e-8) -> np.ndarray:
    return np.clip(np.asarray(p, dtype=float), eps, 1.0 - eps)


def _logit_from_probabilities(p: np.ndarray) -> np.ndarray:
    p_used = _clip_probabilities(p)
    return np.log(p_used / (1.0 - p_used))


@dataclass
class _BaseWrapper:
    estimator: Any
    requires_dense: bool = False
    is_ranker: bool = False
    use_column_named_frame: bool = False

    @staticmethod
    def _group_sizes(query_groups: np.ndarray | None) -> list[int] | None:
        if query_groups is None:
            return None
        groups = np.asarray(query_groups, dtype=object).reshape(-1)
        if groups.size == 0:
            return None
        sizes: list[int] = []
        last = groups[0]
        count = 1
        for item in groups[1:]:
            if item == last:
                count += 1
            else:
                sizes.append(int(count))
                last = item
                count = 1
        sizes.append(int(count))
        return sizes

    def fit(
        self,
        X: sparse.csr_matrix,
        y: np.ndarray,
        *,
        sample_weight: np.ndarray | None = None,
        query_groups: np.ndarray | None = None,
    ) -> "_BaseWrapper":
        fit_kwargs: dict[str, Any] = {}
        if sample_weight is not None:
            fit_kwargs["sample_weight"] = np.asarray(sample_weight, dtype=float)
        if self.is_ranker and query_groups is not None:
            group_sizes = self._group_sizes(query_groups)
            if group_sizes:
                fit_kwargs["group"] = group_sizes
        if self.is_ranker and "group" not in fit_kwargs:
            raise ValueError("Ranker models require query_groups for fit().")
        self.estimator.fit(
            _as_model_input(
                X,
                requires_dense=self.requires_dense,
                with_column_names=self.use_column_named_frame,
            ),
            y,
            **fit_kwargs,
        )
        return self

    def predict_proba(self, X: sparse.csr_matrix) -> np.ndarray:
        if self.is_ranker:
            raw = self.raw_score(X)
            p = 1.0 / (1.0 + np.exp(-np.asarray(raw, dtype=float)))
            return np.vstack([1.0 - p, p]).T
        probs = self.estimator.predict_proba(
            _as_model_input(
                X,
                requires_dense=self.requires_dense,
                with_column_names=self.use_column_named_frame,
            )
        )
        probs_arr = np.asarray(probs, dtype=float)
        if probs_arr.ndim != 2 or probs_arr.shape[1] < 2:
            raise ValueError("Expected predict_proba output with shape (n_rows, >=2).")
        return probs_arr

    def raw_score(self, X: sparse.csr_matrix) -> np.ndarray:
        if hasattr(self.estimator, "decision_function"):
            raw = np.asarray(
                self.estimator.decision_function(
                    _as_model_input(
                        X,
                        requires_dense=self.requires_dense,
                        with_column_names=self.use_column_named_frame,
                    )
                ),
                dtype=float,
            )
            if raw.ndim == 2 and raw.shape[1] > 1:
                raw = raw[:, 1]
            return raw.reshape(-1)

        if hasattr(self.estimator, "predict"):
            try:
                raw_margin = np.asarray(
                    self.estimator.predict(
                        _as_model_input(
                            X,
                            requires_dense=self.requires_dense,
                            with_column_names=self.use_column_named_frame,
                        ),
                        output_margin=True,
                    ),
                    dtype=float,
                )
                return raw_margin.reshape(-1)
            except Exception:
                pass

        p_active = self.predict_proba(X)[:, 1]
        return _logit_from_probabilities(p_active)


class SklearnLogRegWrapper(_BaseWrapper):
    pass


class HistGradientBoostingWrapper(_BaseWrapper):
    pass


class LightGBMWrapper(_BaseWrapper):
    pass


class XGBoostWrapper(_BaseWrapper):
    pass


class LightGBMRankerWrapper(_BaseWrapper):
    pass


class XGBoostRankerWrapper(_BaseWrapper):
    pass


def build_model(
    *,
    model_family: str,
    model_params: dict[str, Any] | None,
    random_seed: int,
) -> _BaseWrapper:
    family = str(model_family).strip().lower()
    params = dict(model_params or {})

    if family in {"logreg", "logistic", "logistic_regression"}:
        estimator = LogisticRegression(
            C=float(params.get("C", 1.0)),
            class_weight=params.get("class_weight", "balanced"),
            max_iter=int(params.get("max_iter", 1000)),
            solver=str(params.get("solver", "liblinear")),
            penalty=str(params.get("penalty", "l2")),
            random_state=int(random_seed),
        )
        return SklearnLogRegWrapper(estimator=estimator, requires_dense=False)

    if family in {"hgb", "histgb", "hist_gradient_boosting"}:
        estimator = HistGradientBoostingClassifier(
            max_depth=params.get("max_depth"),
            learning_rate=float(params.get("learning_rate", 0.1)),
            max_leaf_nodes=int(params.get("max_leaf_nodes", 31)),
            min_samples_leaf=int(params.get("min_samples_leaf", 20)),
            random_state=int(random_seed),
        )
        return HistGradientBoostingWrapper(estimator=estimator, requires_dense=True)

    if family in {"lightgbm", "lgbm"}:
        if LGBMClassifier is None:
            raise RuntimeError(
                "model_family='lightgbm' requested but lightgbm is not installed."
            )
        estimator = LGBMClassifier(
            learning_rate=float(params.get("learning_rate", 0.1)),
            num_leaves=int(params.get("num_leaves", 31)),
            min_child_samples=int(params.get("min_child_samples", 20)),
            n_estimators=int(params.get("n_estimators", 300)),
            max_depth=int(params["max_depth"]) if params.get("max_depth") is not None else -1,
            subsample=float(params.get("subsample", 1.0)),
            colsample_bytree=float(params.get("colsample_bytree", 1.0)),
            class_weight=params.get("class_weight", "balanced"),
            n_jobs=int(params.get("n_jobs", 1)),
            random_state=int(random_seed),
            verbosity=-1,
        )
        return LightGBMWrapper(
            estimator=estimator,
            requires_dense=False,
            use_column_named_frame=True,
        )

    if family in {"xgboost", "xgb"}:
        if XGBClassifier is None:
            raise RuntimeError(
                "model_family='xgboost' requested but xgboost is not installed."
            )
        estimator = XGBClassifier(
            max_depth=int(params.get("max_depth", 6)),
            learning_rate=float(params.get("eta", params.get("learning_rate", 0.1))),
            n_estimators=int(params.get("n_estimators", 300)),
            subsample=float(params.get("subsample", 1.0)),
            colsample_bytree=float(params.get("colsample_bytree", params.get("colsample", 1.0))),
            reg_lambda=float(params.get("reg_lambda", 1.0)),
            reg_alpha=float(params.get("reg_alpha", 0.0)),
            objective="binary:logistic",
            eval_metric="logloss",
            random_state=int(random_seed),
            n_jobs=int(params.get("n_jobs", 1)),
            verbosity=0,
        )
        return XGBoostWrapper(estimator=estimator, requires_dense=False)

    if family in {"lightgbm_rank", "lgbm_rank", "lightgbm_ranker"}:
        if LGBMRanker is None:
            raise RuntimeError(
                "model_family='lightgbm_rank' requested but lightgbm is not installed."
            )
        estimator = LGBMRanker(
            objective=str(params.get("objective", "lambdarank")),
            learning_rate=float(params.get("learning_rate", 0.1)),
            num_leaves=int(params.get("num_leaves", 31)),
            min_child_samples=int(params.get("min_child_samples", 20)),
            n_estimators=int(params.get("n_estimators", 300)),
            max_depth=int(params["max_depth"]) if params.get("max_depth") is not None else -1,
            subsample=float(params.get("subsample", 1.0)),
            colsample_bytree=float(params.get("colsample_bytree", 1.0)),
            n_jobs=int(params.get("n_jobs", 1)),
            random_state=int(random_seed),
            verbosity=-1,
        )
        return LightGBMRankerWrapper(
            estimator=estimator,
            requires_dense=False,
            is_ranker=True,
            use_column_named_frame=True,
        )

    if family in {"xgboost_rank", "xgb_rank", "xgboost_ranker"}:
        if XGBRanker is None:
            raise RuntimeError(
                "model_family='xgboost_rank' requested but xgboost is not installed."
            )
        estimator = XGBRanker(
            objective=str(params.get("objective", "rank:pairwise")),
            max_depth=int(params.get("max_depth", 6)),
            learning_rate=float(params.get("eta", params.get("learning_rate", 0.1))),
            n_estimators=int(params.get("n_estimators", 300)),
            subsample=float(params.get("subsample", 1.0)),
            colsample_bytree=float(params.get("colsample_bytree", params.get("colsample", 1.0))),
            reg_lambda=float(params.get("reg_lambda", 1.0)),
            reg_alpha=float(params.get("reg_alpha", 0.0)),
            eval_metric=str(params.get("eval_metric", "ndcg")),
            random_state=int(random_seed),
            n_jobs=int(params.get("n_jobs", 1)),
            verbosity=0,
        )
        return XGBoostRankerWrapper(estimator=estimator, requires_dense=False, is_ranker=True)

    raise ValueError(f"Unsupported model_family={model_family!r}")
