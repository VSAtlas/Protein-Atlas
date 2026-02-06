from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _as_arrays(y_true: np.ndarray, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    y_arr = np.asarray(y_true, dtype=float).reshape(-1)
    p_arr = np.asarray(p, dtype=float).reshape(-1)
    if y_arr.shape[0] != p_arr.shape[0]:
        raise ValueError("y_true and p must have the same number of rows.")
    if y_arr.size == 0:
        raise ValueError("Cannot compute metrics for empty arrays.")
    p_arr = np.clip(p_arr, 0.0, 1.0)
    return y_arr, p_arr


def brier_score(y_true: np.ndarray, p: np.ndarray) -> float:
    y_arr, p_arr = _as_arrays(y_true, p)
    return float(np.mean((p_arr - y_arr) ** 2))


def reliability_curve(
    y_true: np.ndarray,
    p: np.ndarray,
    n_bins: int = 15,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    y_arr, p_arr = _as_arrays(y_true, p)
    n_bins_used = max(2, int(n_bins))
    bin_edges = np.linspace(0.0, 1.0, n_bins_used + 1)
    bin_ids = np.clip(np.digitize(p_arr, bin_edges, right=True) - 1, 0, n_bins_used - 1)

    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    frac_pos = np.full(n_bins_used, np.nan, dtype=float)
    mean_pred = np.full(n_bins_used, np.nan, dtype=float)
    counts = np.zeros(n_bins_used, dtype=int)

    for idx in range(n_bins_used):
        mask = bin_ids == idx
        count = int(mask.sum())
        counts[idx] = count
        if count <= 0:
            continue
        frac_pos[idx] = float(np.mean(y_arr[mask]))
        mean_pred[idx] = float(np.mean(p_arr[mask]))

    return centers, frac_pos, mean_pred, counts


def expected_calibration_error(
    y_true: np.ndarray,
    p: np.ndarray,
    n_bins: int = 15,
) -> float:
    _centers, frac_pos, mean_pred, counts = reliability_curve(y_true, p, n_bins=n_bins)
    total = float(np.sum(counts))
    if total <= 0:
        return float("nan")

    ece = 0.0
    for idx, count in enumerate(counts):
        if count <= 0:
            continue
        if not np.isfinite(frac_pos[idx]) or not np.isfinite(mean_pred[idx]):
            continue
        ece += (count / total) * abs(float(frac_pos[idx]) - float(mean_pred[idx]))
    return float(ece)


def save_reliability_plot(
    path: Path,
    curve_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> None:
    centers, frac_pos, mean_pred, counts = curve_data
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    valid = (counts > 0) & np.isfinite(frac_pos) & np.isfinite(mean_pred)

    plt.figure(figsize=(5.6, 4.8))
    plt.plot([0.0, 1.0], [0.0, 1.0], linestyle="--", color="gray", linewidth=1.0, label="Perfect")
    if np.any(valid):
        plt.plot(
            mean_pred[valid],
            frac_pos[valid],
            marker="o",
            linewidth=1.5,
            label="Model",
        )
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    plt.xlabel("Mean predicted probability")
    plt.ylabel("Observed positive fraction")
    plt.title("Reliability Diagram")
    plt.grid(alpha=0.25)
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()
