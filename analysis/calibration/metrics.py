from __future__ import annotations

import math
from typing import Sequence

from analysis.statistics import (
    average_precision,
    auroc,
    brier_score,
    expected_calibration_error,
    ranked_binary_metrics,
    threshold_binary_metrics,
)


def _clip_probability(value: float) -> float:
    return min(1.0 - 1e-15, max(1e-15, float(value)))


def _log_loss(probabilities: Sequence[float], labels: Sequence[int]) -> float:
    if not probabilities:
        return 0.0
    total = 0.0
    for prob, label in zip(probabilities, labels):
        p = _clip_probability(float(prob))
        total += -(int(label) * math.log(p) + (1 - int(label)) * math.log(1.0 - p))
    return total / len(probabilities)


def _max_calibration_error(probabilities: Sequence[float], labels: Sequence[int], *, n_bins: int = 10) -> float:
    if not probabilities:
        return 0.0
    worst = 0.0
    for idx in range(n_bins):
        lo = idx / n_bins
        hi = (idx + 1) / n_bins
        members = [
            (p, y)
            for p, y in zip(probabilities, labels)
            if lo <= p < hi or (idx == n_bins - 1 and p == 1.0)
        ]
        if not members:
            continue
        conf = sum(float(p) for p, _y in members) / len(members)
        acc = sum(int(y) for _p, y in members) / len(members)
        worst = max(worst, abs(acc - conf))
    return worst


def _adaptive_ece(probabilities: Sequence[float], labels: Sequence[int], *, n_bins: int = 10) -> float:
    if not probabilities:
        return 0.0
    pairs = sorted((float(p), int(y)) for p, y in zip(probabilities, labels))
    n = len(pairs)
    if n == 0:
        return 0.0
    ece = 0.0
    bin_size = max(1, math.ceil(n / n_bins))
    for start in range(0, n, bin_size):
        members = pairs[start : start + bin_size]
        conf = sum(p for p, _y in members) / len(members)
        acc = sum(y for _p, y in members) / len(members)
        ece += (len(members) / n) * abs(acc - conf)
    return ece


def calibration_metrics(probabilities: Sequence[float], labels: Sequence[int], scores: Sequence[float] | None = None) -> dict[str, float]:
    ranking_scores = list(scores or probabilities)
    out = ranked_binary_metrics(ranking_scores, labels, [0.01, 0.05, 0.10])
    out["Brier"] = brier_score(probabilities, labels)
    out["ECE"] = expected_calibration_error(probabilities, labels)
    out["adaptive_ECE"] = _adaptive_ece(probabilities, labels)
    out["MCE"] = _max_calibration_error(probabilities, labels)
    out["log_loss"] = _log_loss(probabilities, labels)
    out.update(threshold_binary_metrics(probabilities, labels, threshold=0.5))
    out["AUROC"] = auroc(ranking_scores, labels)
    out["AUPRC"] = average_precision(ranking_scores, labels)
    return out
