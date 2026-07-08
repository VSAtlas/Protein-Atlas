from __future__ import annotations

import math
import random
from typing import Sequence


STATISTICAL_FORMULAS = (
    {
        "title": "Enrichment factor",
        "formula": "EF@f = (positives_in_top_k / k) / (total_positives / N), with k = ceil(f * N)",
        "notes": "Implemented by ranked_binary_metrics for top fractions such as 1%, 5%, and 10%.",
        "source": "analysis.statistics.ranked_binary_metrics",
    },
    {
        "title": "Average precision / AUPRC",
        "formula": "AP = (1 / P) * sum_i precision_at_i * I(row_i is positive)",
        "notes": "Rows are sorted by score descending; P is the number of labeled positives.",
        "source": "analysis.statistics.average_precision",
    },
    {
        "title": "AUROC",
        "formula": "AUROC = Pr(score_positive > score_negative) + 0.5 * Pr(tie)",
        "notes": "Implemented as pairwise positive-vs-negative score comparison.",
        "source": "analysis.statistics.auroc",
    },
    {
        "title": "F1 score",
        "formula": "F1 = 2 * precision * recall / (precision + recall), after thresholding probabilities into predicted positives/negatives",
        "notes": "Atlas reports F1@0.5 for classifier outputs. Use AUPRC/enrichment/top-K for ranking claims and F1 only when a decision threshold is meaningful.",
        "source": "analysis.statistics.threshold_binary_metrics",
    },
    {
        "title": "Expected calibration error",
        "formula": "ECE = sum_b (n_b / N) * |observed_positive_rate_b - mean_predicted_probability_b|",
        "notes": "Used for ML calibration summaries; adaptive-bin ECE and maximum calibration error are also emitted by analysis.calibration.metrics.",
        "source": "analysis.statistics.expected_calibration_error",
    },
    {
        "title": "Split-conformal prediction set",
        "formula": "S(x) = {y: nonconformity_y(x) <= q_(1-alpha)}, where q is the finite-sample calibration quantile",
        "notes": "For binary classifiers Atlas uses nonconformity 1 minus the probability assigned to the candidate class.",
        "source": "analysis.ml.conformal.append_split_conformal_sets",
    },
    {
        "title": "Permutation / label-shuffle p-value",
        "formula": "p_perm = (1 + #{null_metric >= observed_metric}) / (1 + n_permutations)",
        "notes": "Used by permutation_p_value and the enrichment suite's label-shuffle controls.",
        "source": "analysis.statistics.permutation_p_value",
    },
)


def empirical_p_value(observed: float, null_scores: Sequence[float], *, higher_is_better: bool = True) -> float:
    if higher_is_better:
        extreme = sum(1 for score in null_scores if score >= observed)
    else:
        extreme = sum(1 for score in null_scores if score <= observed)
    return (extreme + 1.0) / (len(null_scores) + 1.0)


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    q_values = [1.0] * len(p_values)
    running = 1.0
    total = len(p_values)
    for rank_from_end, (idx, p_value) in enumerate(reversed(indexed), start=1):
        rank = total - rank_from_end + 1
        running = min(running, p_value * total / rank)
        q_values[idx] = min(1.0, running)
    return q_values


def ranked_binary_metrics(scores: Sequence[float], labels: Sequence[int], top_fractions: Sequence[float]) -> dict[str, float]:
    pairs = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    n = len(pairs)
    positives = sum(label for _score, label in pairs)
    overall_rate = positives / n if n else 0.0
    out: dict[str, float] = {}
    for fraction in top_fractions:
        top_n = max(1, math.ceil(n * fraction)) if n else 0
        top_rate = sum(label for _score, label in pairs[:top_n]) / top_n if top_n else 0.0
        out[f"EF@{int(fraction * 100)}%"] = top_rate / overall_rate if overall_rate else 0.0
    out["AUPRC"] = average_precision(scores, labels)
    out["AUROC"] = auroc(scores, labels)
    return out


def average_precision(scores: Sequence[float], labels: Sequence[int]) -> float:
    pairs = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    positives = sum(labels)
    if positives <= 0:
        return 0.0
    hits = 0
    total = 0.0
    for rank, (_score, label) in enumerate(pairs, start=1):
        if label:
            hits += 1
            total += hits / rank
    return total / positives


def auroc(scores: Sequence[float], labels: Sequence[int]) -> float:
    pos = [score for score, label in zip(scores, labels) if label]
    neg = [score for score, label in zip(scores, labels) if not label]
    if not pos or not neg:
        return 0.0
    wins = 0.0
    for p_score in pos:
        for n_score in neg:
            if p_score > n_score:
                wins += 1.0
            elif p_score == n_score:
                wins += 0.5
    return wins / (len(pos) * len(neg))



def threshold_binary_metrics(
    probabilities: Sequence[float],
    labels: Sequence[int],
    *,
    threshold: float = 0.5,
) -> dict[str, float]:
    if not probabilities:
        return {
            f"precision@{threshold:g}": 0.0,
            f"recall@{threshold:g}": 0.0,
            f"F1@{threshold:g}": 0.0,
            f"specificity@{threshold:g}": 0.0,
            f"accuracy@{threshold:g}": 0.0,
            f"balanced_accuracy@{threshold:g}": 0.0,
        }
    predicted = [1 if float(prob) >= threshold else 0 for prob in probabilities]
    truth = [1 if int(label) else 0 for label in labels]
    tp = sum(1 for pred, label in zip(predicted, truth) if pred == 1 and label == 1)
    fp = sum(1 for pred, label in zip(predicted, truth) if pred == 1 and label == 0)
    fn = sum(1 for pred, label in zip(predicted, truth) if pred == 0 and label == 1)
    tn = sum(1 for pred, label in zip(predicted, truth) if pred == 0 and label == 0)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    specificity = tn / (tn + fp) if tn + fp else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    accuracy = (tp + tn) / len(truth) if truth else 0.0
    balanced_accuracy = (recall + specificity) / 2.0
    suffix = f"@{threshold:g}"
    return {
        f"precision{suffix}": precision,
        f"recall{suffix}": recall,
        f"F1{suffix}": f1,
        f"specificity{suffix}": specificity,
        f"accuracy{suffix}": accuracy,
        f"balanced_accuracy{suffix}": balanced_accuracy,
        f"tp{suffix}": float(tp),
        f"fp{suffix}": float(fp),
        f"fn{suffix}": float(fn),
        f"tn{suffix}": float(tn),
    }

def brier_score(probabilities: Sequence[float], labels: Sequence[int]) -> float:
    if not probabilities:
        return 0.0
    return sum((prob - label) ** 2 for prob, label in zip(probabilities, labels)) / len(probabilities)


def expected_calibration_error(probabilities: Sequence[float], labels: Sequence[int], *, n_bins: int = 10) -> float:
    if not probabilities:
        return 0.0
    ece = 0.0
    total = len(probabilities)
    for idx in range(n_bins):
        lo = idx / n_bins
        hi = (idx + 1) / n_bins
        members = [(p, y) for p, y in zip(probabilities, labels) if lo <= p < hi or (idx == n_bins - 1 and p == 1.0)]
        if not members:
            continue
        conf = sum(p for p, _y in members) / len(members)
        acc = sum(y for _p, y in members) / len(members)
        ece += (len(members) / total) * abs(acc - conf)
    return ece


def permutation_p_value(
    scores: Sequence[float],
    labels: Sequence[int],
    metric_name: str,
    *,
    n_permutations: int,
    seed: int,
) -> float:
    observed = ranked_binary_metrics(scores, labels, [0.01, 0.05, 0.10]).get(metric_name, 0.0)
    rng = random.Random(seed)
    null = list(labels)
    extreme = 0
    for _idx in range(n_permutations):
        rng.shuffle(null)
        value = ranked_binary_metrics(scores, null, [0.01, 0.05, 0.10]).get(metric_name, 0.0)
        if value >= observed:
            extreme += 1
    return (extreme + 1.0) / (n_permutations + 1.0)
