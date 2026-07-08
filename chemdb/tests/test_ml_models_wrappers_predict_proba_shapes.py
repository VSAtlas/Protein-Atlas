from __future__ import annotations

import numpy as np
from scipy import sparse

from ml.models import build_model


def test_ml_models_wrappers_predict_proba_shapes():
    y = np.asarray(([0, 1] * 20), dtype=int)
    x_dense = np.asarray(
        [[float(idx % 5), float((idx * 7) % 9), float(y[idx])] for idx in range(len(y))],
        dtype=float,
    )
    X = sparse.csr_matrix(x_dense)

    tested = 0
    for family, params in (
        (
            "logreg",
            {
                "C": 1.0,
                "class_weight": "balanced",
                "max_iter": 500,
                "solver": "liblinear",
                "penalty": "l2",
            },
        ),
        ("lightgbm", {"n_estimators": 50, "num_leaves": 31}),
        ("xgboost", {"n_estimators": 50, "max_depth": 4, "eta": 0.1}),
    ):
        try:
            model = build_model(model_family=family, model_params=params, random_seed=7)
        except RuntimeError:
            continue
        model.fit(X, y)
        prob = model.predict_proba(X)
        assert prob.shape == (len(y), 2)
        assert np.all(prob >= 0.0)
        assert np.all(prob <= 1.0)
        assert np.allclose(np.sum(prob, axis=1), 1.0, atol=1e-6)

        raw = model.raw_score(X)
        assert raw.shape == (len(y),)
        assert np.all(np.isfinite(raw))
        tested += 1

    assert tested >= 1
