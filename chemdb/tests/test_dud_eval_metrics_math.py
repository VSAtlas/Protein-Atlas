from __future__ import annotations

import numpy as np
import pytest
from sklearn.metrics import roc_curve

from analysis.dud_eval_metrics import bedroc, ef_at_fractions, log_auc_from_roc, pr_auc


def test_ef_at_fractions_simple_case() -> None:
    y_true = np.array([1, 1, 0, 0])
    scores_low = np.array([-9.0, -8.0, -4.0, -3.0])

    ef = ef_at_fractions(y_true, scores_low)

    assert ef["EF@1%"] == pytest.approx(2.0)
    assert ef["EF@2%"] == pytest.approx(2.0)
    assert ef["EF@5%"] == pytest.approx(2.0)
    assert ef["EF@10%"] == pytest.approx(2.0)


def test_pr_auc_and_log_auc_and_bedroc_monotonic_sanity() -> None:
    y_true = np.array([1, 0, 1, 0])
    good_scores = np.array([0.9, 0.2, 0.8, 0.1])
    bad_scores = np.array([0.6, 0.5, 0.4, 0.3])

    pr_good, _, _ = pr_auc(y_true, good_scores)
    pr_bad, _, _ = pr_auc(y_true, bad_scores)
    assert pr_good > pr_bad

    fpr_good, tpr_good, _ = roc_curve(y_true, good_scores)
    fpr_bad, tpr_bad, _ = roc_curve(y_true, bad_scores)
    log_good = log_auc_from_roc(fpr_good, tpr_good, lam=1e-3)
    log_bad = log_auc_from_roc(fpr_bad, tpr_bad, lam=1e-3)
    assert log_good > log_bad

    bed_good = bedroc(y_true, good_scores, alpha=20.0)
    bed_bad = bedroc(y_true, bad_scores, alpha=20.0)
    assert np.isfinite(bed_good)
    assert np.isfinite(bed_bad)
    assert bed_good > bed_bad


def test_bedroc_degenerate_set_returns_nan() -> None:
    y_true_all_actives = np.array([1, 1, 1])
    scores = np.array([1.0, 0.5, 0.2])
    assert np.isnan(bedroc(y_true_all_actives, scores, alpha=20.0))
