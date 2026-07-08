from __future__ import annotations

import math

import numpy as np
import pytest

from analysis.dud_eval_core.metrics import bedroc, ef_at_fractions


def _bedroc_expected(y_true: np.ndarray, y_score: np.ndarray, alpha: float) -> float:
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    n_total = len(y_true)
    n_actives = int(y_true.sum())
    if n_total == 0 or n_actives == 0 or n_actives == n_total:
        return float("nan")

    order = np.argsort(-y_score)
    ranks = np.nonzero(y_true[order] == 1)[0] + 1
    s = float(np.sum(np.exp(-alpha * ranks / n_total)))
    ra = n_actives / n_total
    k1 = (ra * (1.0 - math.exp(-alpha))) / (math.exp(alpha / n_total) - 1.0)
    k2 = (ra * math.sinh(alpha / 2.0)) / (
        math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra)
    )
    return float((s / k1) * k2)


def _bedroc_zero_based_reference(
    y_true: np.ndarray, y_score: np.ndarray, alpha: float
) -> float:
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    n_total = len(y_true)
    n_actives = int(y_true.sum())
    if n_total == 0 or n_actives == 0 or n_actives == n_total:
        return float("nan")

    order = np.argsort(-y_score)
    ranks = np.nonzero(y_true[order] == 1)[0]
    s = float(np.sum(np.exp(-alpha * ranks / n_total)))
    ra = n_actives / n_total
    k1 = (ra * (1.0 - math.exp(-alpha))) / (math.exp(alpha / n_total) - 1.0)
    k2 = (ra * math.sinh(alpha / 2.0)) / (
        math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra)
    )
    return float((s / k1) * k2)


def test_bedroc_uses_one_based_ranks() -> None:
    y_true = np.array([1, 0, 0, 1, 0, 0], dtype=int)
    y_score = np.array([0.95, 0.8, 0.7, 0.6, 0.2, 0.1], dtype=float)
    alpha = 20.0

    expected = _bedroc_expected(y_true, y_score, alpha)
    got = bedroc(y_true, y_score, alpha)
    old_zero_based = _bedroc_zero_based_reference(y_true, y_score, alpha)

    assert got == pytest.approx(expected, rel=1e-12, abs=1e-12)
    assert got != pytest.approx(old_zero_based)


def test_ef_uses_ceil_cutoff_of_total_ranked_list_fraction() -> None:
    n_total = 149
    y_true = np.zeros(n_total, dtype=int)
    y_true[1] = 1
    scores_low = np.arange(n_total, dtype=float)

    ef = ef_at_fractions(y_true, scores_low, fractions=(0.01,))

    # top 1% means top ceil(0.01 * 149) = 2 molecules in ranked list
    expected = (1.0 / 2.0) / (1.0 / 149.0)
    assert ef["EF@1%"] == pytest.approx(expected)
