import numpy as np

from analysis.dud_eval import ef_at_fractions, pr_auc
from ml.evaluate import evaluate_holdout_metrics


def test_evaluate_holdout_metrics_uses_expected_keys_and_score_direction():
    y_true = np.asarray([1, 0, 1, 0, 0, 1], dtype=int)
    p_active = np.asarray([0.98, 0.10, 0.80, 0.05, 0.20, 0.70], dtype=float)

    report = evaluate_holdout_metrics(y_true, p_active)

    expected_metric_keys = {
        "N",
        "n_actives",
        "actives_fraction",
        "PR_AUC",
        "EF@1%",
        "EF@2%",
        "EF@5%",
        "EF@10%",
        "precision@1%",
        "precision@2%",
        "precision@5%",
        "precision@10%",
    }
    assert expected_metric_keys.issubset(report.keys())

    direct_pr_auc, _, _ = pr_auc(y_true, p_active)
    assert np.isclose(float(report["PR_AUC"]), float(direct_pr_auc))
    assert 0.0 <= float(report["PR_AUC"]) <= 1.0

    expected_ef = ef_at_fractions(y_true, scores_low_is_better=-p_active)
    for key in ("EF@1%", "EF@2%", "EF@5%", "EF@10%"):
        assert np.isclose(float(report[key]), float(expected_ef[key]))

