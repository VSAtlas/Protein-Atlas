from __future__ import annotations

import json
from typing import Mapping

import pandas as pd

from analysis.normalization import bounds, normalize as normalize_numeric, numeric_value


DEFAULT_WEIGHTS = {
    "atlas_score": 0.40,
    "exposure_plausibility": 0.15,
    "tissue_expression": 0.10,
    "target_adr_evidence": 0.15,
    "pathway_evidence": 0.10,
    "structure_quality": 0.10,
}
NORMALIZED_COLUMNS = {"atlas_score", "tissue_expression", "structure_quality"}


def compute_priority_score(
    df: pd.DataFrame,
    weights: Mapping[str, float] | None = None,
    normalize_values: bool = True,
    missing_strategy: str = "ignore_and_rescale",
    normalize: bool | None = None,
) -> pd.DataFrame:
    if normalize is not None:
        normalize_values = normalize
    if missing_strategy != "ignore_and_rescale":
        raise ValueError("only missing_strategy='ignore_and_rescale' is supported")
    out = df.copy()
    active_weights = dict(weights or DEFAULT_WEIGHTS)
    column_bounds = {
        name: bounds(numeric_value(value) for value in out[name])
        for name in active_weights
        if name in out.columns and name in NORMALIZED_COLUMNS
    }
    scores: list[float | None] = []
    component_json: list[str] = []
    missing_components: list[str] = []
    for _idx, row in out.iterrows():
        components: dict[str, float] = {}
        missing: list[str] = []
        for name, weight in active_weights.items():
            if name not in out.columns:
                missing.append(name)
                continue
            value = numeric_value(row.get(name))
            if value is None:
                missing.append(name)
                continue
            if normalize_values and name in NORMALIZED_COLUMNS:
                value = normalize_numeric(value, column_bounds.get(name))
            if value is None:
                missing.append(name)
                continue
            components[name] = float(value) * float(weight)
        available_weight = sum(float(active_weights[name]) for name in components)
        if not components or available_weight <= 0:
            scores.append(None)
        else:
            scores.append(sum(components.values()) / available_weight)
        component_json.append(json.dumps(components, sort_keys=True))
        missing_components.append(";".join(missing))
    out["priority_score"] = scores
    out["priority_score_components_json"] = component_json
    out["priority_score_missing_components"] = missing_components
    return out
