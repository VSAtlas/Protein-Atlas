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
