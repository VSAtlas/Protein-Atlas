from __future__ import annotations

import numpy as np
from scipy import sparse

from ml.calibration import fit_oof_calibrator
from ml.models import build_model


def test_ml_calibration_oof_is_leakage_safe():
    n_groups = 12
    rows_per_group = 6
    n_rows = n_groups * rows_per_group

    y = np.asarray([(idx % 2) for idx in range(n_rows)], dtype=int)
    groups = np.asarray(
        [f"scaffold_{idx // rows_per_group}" for idx in range(n_rows)],
        dtype=object,
    )
    x_dense = np.asarray(
        [
            [float(idx % 7), float((idx * 3) % 11), float(y[idx])]
            for idx in range(n_rows)
        ],
        dtype=float,
    )
    X = sparse.csr_matrix(x_dense)

    def _builder(seed: int):
        return build_model(
            model_family="logreg",
            model_params={
                "C": 1.0,
                "class_weight": "balanced",
                "max_iter": 500,
                "solver": "liblinear",
                "penalty": "l2",
            },
            random_seed=seed,
        )

    calibrator, oof_df = fit_oof_calibrator(
        X=X,
        y=y,
        groups=groups,
        base_model_builder=_builder,
        method="sigmoid",
        cv_folds=4,
        seed=123,
        test_mode=True,
    )

    assert len(oof_df) == n_rows
    assert int(oof_df["used_for_calibration"].sum()) == n_rows
    assert np.all(np.isfinite(oof_df["p_uncalibrated"].to_numpy(dtype=float)))

    for split in calibrator.debug_fold_indices:
        train_set = set(split["train_idx"])
        val_set = set(split["val_idx"])
        assert train_set.isdisjoint(val_set)

    _, oof_df_2 = fit_oof_calibrator(
        X=X,
        y=y,
        groups=groups,
        base_model_builder=_builder,
        method="sigmoid",
        cv_folds=4,
        seed=123,
        test_mode=True,
    )
    assert np.allclose(
        oof_df["p_uncalibrated"].to_numpy(dtype=float),
        oof_df_2["p_uncalibrated"].to_numpy(dtype=float),
    )
    assert np.array_equal(
        oof_df["fold_id"].to_numpy(dtype=int),
        oof_df_2["fold_id"].to_numpy(dtype=int),
    )
