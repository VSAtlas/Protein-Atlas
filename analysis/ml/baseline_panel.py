from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.decision_metrics import random_baseline_rows, score_metric_row


BASELINE_SCORE_SPECS = [
    ("atlas_score_only", ["atlas_score", "z_score", "z_selected"]),
    ("vina_consensus_only", ["consensus_score", "final_score"]),
    ("scorch_final_only", ["pilot_SCORCH_score_used", "SCORCH_score_used", "pilot_final_score"]),
    ("banana_bigbind_only", ["banana_binding_probability", "banana_score_normalized", "banana_score"]),
    ("atlas_scorch_available", ["pilot_final_score", "SCORCH_score_used", "consensus_score"]),
    ("atlas_scorch_banana_available", ["banana_atlas_blend_score", "binding_expert_score"]),
    ("trained_model", ["ml_prediction_score"]),
]


def _first_existing_numeric(df: pd.DataFrame, columns: list[str]) -> str | None:
    for col in columns:
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any():
            return col
    return None


def write_standard_baseline_panel(
    df: pd.DataFrame,
    *,
    label_col: str,
    out_path: str | Path,
    seed: int = 42,
) -> pd.DataFrame:
    """Write a standard baseline table for every trained Atlas ML run."""

    rows: list[dict[str, Any]] = []
    for method, columns in BASELINE_SCORE_SPECS:
        score_col = _first_existing_numeric(df, columns)
        if score_col is None:
            rows.append({"method": method, "score_col": "", "status": "missing_score_column", "n": 0})
            continue
        rows.append(score_metric_row(df, label_col=label_col, score_col=score_col, method=method))
    rows.extend(random_baseline_rows(df, label_col=label_col, seed=seed))
    table = pd.DataFrame(rows)
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    return table
