from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from analysis.dud_eval import ef_at_fractions
from ml.models import build_model


def _get_ranker():
    for family in ("lightgbm_rank", "xgboost_rank"):
        try:
            model = build_model(model_family=family, model_params={"n_estimators": 80}, random_seed=7)
            return model
        except RuntimeError:
            continue
    return None


def test_ranker_training_and_scoring_direction():
    model = _get_ranker()
    if model is None:
        pytest.skip("No ranker backend available (lightgbm/xgboost).")

    n_groups = 8
    rows_per_group = 8
    n_rows = n_groups * rows_per_group
    groups = np.asarray([f"G{idx // rows_per_group}" for idx in range(n_rows)], dtype=object)
    y = np.asarray(
        [1 if (idx % rows_per_group) < 2 else 0 for idx in range(n_rows)],
        dtype=int,
    )
    x = np.asarray(
        [
            [float(y[idx]), float((idx % rows_per_group) / rows_per_group), float(idx % 3)]
            for idx in range(n_rows)
        ],
        dtype=float,
    )
    X = sparse.csr_matrix(x)
    model.fit(X, y, query_groups=groups)
    score = np.asarray(model.raw_score(X), dtype=float)
    assert score.shape == (n_rows,)

    order = np.argsort(score)[::-1]
    top_rate = float(np.mean(y[order[: max(4, n_rows // 4)]]))
    bottom_rate = float(np.mean(y[order[-max(4, n_rows // 4) :]]))
    assert top_rate > bottom_rate

    ef = ef_at_fractions(y, scores_low_is_better=-score, fractions=(0.01, 0.1))
    assert float(ef["EF@10%"]) >= 1.0
