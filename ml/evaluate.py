from __future__ import annotations

import math
from typing import Iterable

import numpy as np

from analysis.dud_eval import ef_at_fractions, pr_auc


DEFAULT_FRACTIONS = (0.01, 0.02, 0.05, 0.10)


def evaluate_holdout_metrics(
    y_true: np.ndarray | Iterable[int],
    p_active: np.ndarray | Iterable[float],
    *,
    fractions: tuple[float, ...] = DEFAULT_FRACTIONS,
) -> dict[str, float | int]:
    y_arr = np.asarray(list(y_true), dtype=int)
    p_arr = np.asarray(list(p_active), dtype=float)
    if y_arr.shape[0] != p_arr.shape[0]:
        raise ValueError("y_true and p_active must have the same number of elements.")

    n_rows = int(y_arr.shape[0])
    n_actives = int(y_arr.sum())
    base_rate = float(n_actives / n_rows) if n_rows else float("nan")

    if n_rows and n_actives > 0:
        pr_auc_value, _, _ = pr_auc(y_arr, p_arr)
    else:
        pr_auc_value = float("nan")

    ef_values = ef_at_fractions(y_arr, scores_low_is_better=-p_arr, fractions=fractions)

    report: dict[str, float | int] = {
        "N": n_rows,
        "n_actives": n_actives,
        "actives_fraction": base_rate,
        "PR_AUC": float(pr_auc_value),
    }
    report.update(ef_values)

    for fraction in fractions:
        pct = int(round(100 * fraction))
        ef_key = f"EF@{pct}%"
        precision_key = f"precision@{pct}%"
        ef_value = float(report.get(ef_key, float("nan")))
        if math.isnan(ef_value) or math.isnan(base_rate):
            report[precision_key] = float("nan")
        else:
            report[precision_key] = float(ef_value * base_rate)

    return report

