from __future__ import annotations

import math
from typing import Iterable, List, Sequence

try:  # optional dependency
    from sklearn.metrics import roc_auc_score  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - optional dependency
    roc_auc_score = None


def _rankdata(scores: Sequence[float]) -> List[float]:
    indexed = list(enumerate(scores))
    indexed.sort(key=lambda item: item[1])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + j + 2) / 2.0  # 1-based ranks
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg_rank
        i = j + 1
    return ranks


def compute_auc(y_true: Sequence[int], y_score: Sequence[float]) -> float:
    if not y_true or len(y_true) != len(y_score):
        return float("nan")
    n_pos = int(sum(1 for y in y_true if y == 1))
    n_neg = int(len(y_true) - n_pos)
    if n_pos <= 0 or n_neg <= 0:
        return float("nan")
    if roc_auc_score is not None:
        try:
            return float(roc_auc_score(y_true, y_score))
        except Exception:
            return float("nan")

    ranks = _rankdata(list(map(float, y_score)))
    sum_pos = sum(r for r, y in zip(ranks, y_true) if y == 1)
    auc = (sum_pos - (n_pos * (n_pos + 1) / 2.0)) / (n_pos * n_neg)
    return float(auc)


def compute_enrichment_factor(
    y_true: Sequence[int],
    y_score: Sequence[float],
    top_frac: float,
    *,
    lower_is_better: bool = True,
) -> float:
    if not y_true or len(y_true) != len(y_score):
        return float("nan")
    if top_frac is None:
        return float("nan")
    try:
        frac = float(top_frac)
    except Exception:
        return float("nan")
    if frac <= 0 or frac > 1:
        return float("nan")

    n_total = len(y_true)
    n_top = max(1, int(math.ceil(frac * n_total)))
    order = sorted(
        range(n_total),
        key=lambda i: float(y_score[i]),
        reverse=not lower_is_better,
    )
    total_hits = int(sum(1 for y in y_true if y == 1))
    if total_hits <= 0:
        return float("nan")
    hits_top = int(sum(1 for i in order[:n_top] if y_true[i] == 1))
    ef = (hits_top / n_top) / (total_hits / n_total)
    return float(ef)


def parse_top_fracs(
    value: object, default: Iterable[float] = (0.01, 0.005)
) -> List[float]:
    if value is None:
        return list(default)
    if isinstance(value, (list, tuple)):
        out: List[float] = []
        for item in value:
            try:
                out.append(float(item))
            except Exception:
                continue
        return out or list(default)
    if isinstance(value, str):
        parts = [p.strip() for p in value.split(",") if p.strip()]
        out = []
        for part in parts:
            try:
                out.append(float(part))
            except Exception:
                continue
        return out or list(default)
    try:
        return [float(str(value))]
    except Exception:
        return list(default)
