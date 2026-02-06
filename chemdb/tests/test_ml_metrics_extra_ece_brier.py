from __future__ import annotations

import numpy as np

from ml.metrics_extra import brier_score, expected_calibration_error, reliability_curve


def test_ml_metrics_extra_ece_brier():
    y_true = np.asarray([0, 1, 0, 1], dtype=int)
    p = np.asarray([0.1, 0.9, 0.4, 0.6], dtype=float)

    score = brier_score(y_true, p)
    assert np.isclose(score, 0.085, atol=1e-12)

    ece = expected_calibration_error(y_true, p, n_bins=2)
    assert np.isclose(ece, 0.25, atol=1e-12)

    centers, frac_pos, mean_pred, counts = reliability_curve(y_true, p, n_bins=2)
    assert centers.shape == (2,)
    assert counts.tolist() == [2, 2]
    assert np.isclose(frac_pos[0], 0.0)
    assert np.isclose(frac_pos[1], 1.0)
    assert np.isclose(mean_pred[0], 0.25)
    assert np.isclose(mean_pred[1], 0.75)
