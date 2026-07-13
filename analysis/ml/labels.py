from __future__ import annotations

import pandas as pd


TRUE_LABELS = {"1", "1.0", "true", "yes", "y", "positive", "active"}
FALSE_LABELS = {"0", "0.0", "false", "no", "n", "negative", "inactive"}


def truthy_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(TRUE_LABELS)


def binary_label_series(series: pd.Series) -> pd.Series:
    missing = series.isna()
    text = series.astype("object").where(~missing, "").astype(str).str.strip().str.lower()
    out = pd.Series(pd.NA, index=series.index, dtype="Int64")
    out = out.mask(text.isin(TRUE_LABELS), 1)
    out = out.mask(text.isin(FALSE_LABELS), 0)
    numeric = pd.to_numeric(series, errors="coerce")
    out = out.mask((~missing) & numeric.eq(1), 1)
    out = out.mask((~missing) & numeric.eq(0), 0)
    return out


def training_eligibility_mask(frame: pd.DataFrame) -> tuple[pd.Series, dict[str, int]]:
    """Apply fail-closed identity and explicit source/benchmark training policy."""

    eligible = pd.Series(True, index=frame.index, dtype=bool)
    reasons: dict[str, int] = {}

    if "identity_training_allowed" in frame.columns:
        allowed = truthy_series(frame["identity_training_allowed"])
        blocked = ~allowed
        eligible &= allowed
        reasons["identity_training_disallowed"] = int(blocked.sum())
    elif "training_allowed" in frame.columns:
        value = frame["training_allowed"]
        text = value.astype("object").where(value.notna(), "").astype(str)
        explicitly_disallowed = text.str.strip().str.lower().isin(FALSE_LABELS)
        eligible &= ~explicitly_disallowed
        reasons["source_training_disallowed"] = int(explicitly_disallowed.sum())

    benchmark = pd.Series(False, index=frame.index, dtype=bool)
    for column in ("benchmark_only", "is_benchmark_only"):
        if column in frame.columns:
            benchmark |= truthy_series(frame[column])
    eligible &= ~benchmark
    reasons["benchmark_only_excluded"] = int(benchmark.sum())
    reasons["input_rows"] = int(len(frame))
    reasons["eligible_rows"] = int(eligible.sum())
    reasons["excluded_rows"] = int((~eligible).sum())
    return eligible, reasons
