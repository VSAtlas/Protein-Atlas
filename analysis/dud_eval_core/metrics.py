from __future__ import annotations

# Metric equations used by DUD eval:
# ROC_AUC:
#   ROC_AUC = AUC(FPR, TPR), computed from y_true and y_score_high_is_better.
# PR_AUC:
#   PR_AUC = AUC(recall, precision), with precision/recall from
#   precision_recall_curve(y_true, y_score_high_is_better).
# logAUC (lambda = lam, default lam=1e-3):
#   logAUC = [ integral_{FPR=lam..1} TPR(FPR) d(log10(FPR)) ] / log10(1/lam).
#   Implemented by sorting ROC points by FPR, inserting/interpolating a point at
#   FPR=lam if needed, then trapezoidal integration in log10(FPR) space.
# logAUC_adj:
#   logAUC_adj = logAUC - 0.14462.
# BEDROC_alpha_20 (general alpha shown; alpha=20 in common usage):
#   N = total compounds, n = number of actives, ra = n/N.
#   Rank by descending y_score_high_is_better, active ranks r_i are 1-based.
#   s  = sum_i exp(-alpha * r_i / N)
#   k1 = (ra * (1 - exp(-alpha))) / (exp(alpha / N) - 1)
#   k2 = (ra * sinh(alpha / 2)) /
#        (cosh(alpha / 2) - cosh(alpha / 2 - alpha * ra))
#   BEDROC = (s / k1) * k2.
# EF at fraction f of the total ranked list (EF@1%, EF@2%, EF@5%, EF@10%):
#   Sort low-is-better scores ascending, k = max(1, ceil(f * N)),
#   found = actives in top-k, hit_rate = found/k, base_rate = n_act/N.
#   EF@f = hit_rate / base_rate = (found * N) / (k * n_act), when n_act > 0.
#   This is rank-list-fraction EF (not ROC/FPR-based EF).

import math
from typing import Dict, Tuple

import numpy as np
from sklearn.metrics import auc, precision_recall_curve

def ef_at_fractions(
    y_true: np.ndarray,
    scores_low_is_better: np.ndarray,
    fractions=(0.01, 0.02, 0.05, 0.10),
) -> Dict[str, float]:
    """Compute EF@f using top ceil(f*N) compounds of the ranked total list."""
    N = len(y_true)
    n_act = int(y_true.sum())
    out = {}
    if N == 0 or n_act == 0:
        return {f"EF@{int(fr * 100)}%": float("nan") for fr in fractions}
    order = np.argsort(scores_low_is_better)  # lowest score first
    y_sorted = y_true[order]
    cum_act = np.cumsum(y_sorted)
    for fr in fractions:
        k = max(1, int(math.ceil(fr * N)))
        found = int(cum_act[k - 1])
        hit_rate = found / k
        base_rate = n_act / N
        out[f"EF@{int(fr * 100)}%"] = (
            float(hit_rate / base_rate) if base_rate > 0 else float("nan")
        )
    return out

def pr_auc(
    y_true: np.ndarray, y_score_high_is_better: np.ndarray
) -> Tuple[float, np.ndarray, np.ndarray]:
    precision, recall, _ = precision_recall_curve(y_true, y_score_high_is_better)
    return float(auc(recall, precision)), precision, recall

def log_auc_from_roc(fpr: np.ndarray, tpr: np.ndarray, lam: float = 1e-3) -> float:
    order = np.argsort(fpr)
    fpr = fpr[order]
    tpr = tpr[order]
    mask = fpr >= lam
    if not np.any(mask):
        return 0.0
    # insert point at lam if needed
    if not np.isclose(fpr[mask][0], lam):
        i = np.searchsorted(fpr, lam)
        x0, x1 = fpr[i - 1], fpr[i]
        y0, y1 = tpr[i - 1], tpr[i]
        ylam = y0 + (y1 - y0) * (lam - x0) / (x1 - x0)
        fpr = np.insert(fpr, i, lam)
        tpr = np.insert(tpr, i, ylam)
        mask = fpr >= lam
    xf = fpr[mask]
    yf = tpr[mask]
    logx = np.log10(xf)
    num: float = float(np.sum((logx[1:] - logx[:-1]) * (yf[1:] + yf[:-1]) / 2.0))
    denom = math.log10(1.0 / lam)
    return float(num / denom) if denom > 0 else float("nan")

def bedroc(
    y_true: np.ndarray, y_score_high_is_better: np.ndarray, alpha: float = 20.0
) -> float:
    """Compute BEDROC with 1-based active ranks in the exponential term."""
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score_high_is_better, dtype=float)
    N = len(y_true)
    n = int(y_true.sum())
    if N == 0 or n == 0 or n == N:
        return float("nan")
    order = np.argsort(-y_score)
    # Use 1-based ranks for BEDROC/RIE consistency.
    ranks = np.nonzero(y_true[order] == 1)[0] + 1
    s = float(np.sum(np.exp(-alpha * ranks / N)))
    ra = n / N
    # constants per Truchon & Bayly (2007), matching common implementations
    k1 = (ra * (1.0 - math.exp(-alpha))) / (math.exp(alpha / N) - 1.0)
    k2 = (ra * math.sinh(alpha / 2.0)) / (
        math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra)
    )
    return float((s / k1) * k2)
