from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


_STRENGTH_COLUMNS = (
    "pchembl_value",
    "pchembl",
    "activity_value",
    "standard_value",
)


@dataclass(frozen=True)
class LabelStrength:
    values: np.ndarray
    source_column: str | None


def _to_float_or_nan(value: Any) -> float:
    try:
        parsed = float(value)
        if np.isfinite(parsed):
            return float(parsed)
    except Exception:
        pass
    return float("nan")


def parse_activity_strength(df: pd.DataFrame) -> LabelStrength:
    for col in _STRENGTH_COLUMNS:
        if col not in df.columns:
            continue
        values = np.asarray([_to_float_or_nan(v) for v in df[col].tolist()], dtype=float)
        if np.isfinite(values).any():
            return LabelStrength(values=values, source_column=col)
    return LabelStrength(values=np.full(len(df), np.nan, dtype=float), source_column=None)


def derive_sample_weight(df: pd.DataFrame, *, weight_cap: float = 5.0) -> np.ndarray:
    strength = parse_activity_strength(df)
    weights = np.ones(len(df), dtype=float)
    if len(df) == 0 or strength.source_column is None:
        return weights

    values = strength.values
    finite = np.isfinite(values)
    if not finite.any():
        return weights

    min_v = float(np.nanmin(values[finite]))
    max_v = float(np.nanmax(values[finite]))
    if max_v <= min_v:
        return weights

    scaled = (values - min_v) / (max_v - min_v + 1e-12)
    scaled[~finite] = 0.0
    active = np.asarray(df.get("active", pd.Series([0] * len(df))), dtype=int)
    weights = 1.0 + scaled * np.clip(active, 0, 1)
    return np.clip(weights, 1.0, max(1.0, float(weight_cap)))


def apply_label_smoothing(y: np.ndarray, eps: float) -> np.ndarray:
    y_arr = np.asarray(y, dtype=float).reshape(-1)
    eps_used = float(np.clip(eps, 0.0, 0.49))
    return y_arr * (1.0 - eps_used) + 0.5 * eps_used
