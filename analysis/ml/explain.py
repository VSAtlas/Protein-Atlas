from __future__ import annotations

from pathlib import Path

import pandas as pd


def write_feature_importance(model: object, feature_names: list[str], out_path: str | Path) -> pd.DataFrame:
    values = None
    if hasattr(model, "coef_"):
        values = getattr(model, "coef_")[0]
    elif hasattr(model, "feature_importances_"):
        values = getattr(model, "feature_importances_")
    if values is None:
        values = [0.0] * len(feature_names)
    if len(values) != len(feature_names):
        values = list(values)[: len(feature_names)] + [0.0] * max(0, len(feature_names) - len(values))
    out = pd.DataFrame({"feature": feature_names, "importance": values}).sort_values("importance", key=lambda s: s.abs(), ascending=False)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return out
