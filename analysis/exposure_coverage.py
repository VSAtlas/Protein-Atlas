from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


EXPOSURE_COLUMNS = ("free_cmax_um", "cmax_um", "fraction_unbound_plasma")


def _coverage_row(df: pd.DataFrame, group: str, group_value: str) -> dict[str, Any]:
    row: dict[str, Any] = {"group": group, "group_value": group_value, "n_rows": int(len(df))}
    for col in EXPOSURE_COLUMNS:
        if col in df.columns:
            values = pd.to_numeric(df[col], errors="coerce")
            row[f"{col}_nonmissing"] = int(values.notna().sum())
            row[f"{col}_coverage"] = float(values.notna().mean()) if len(values) else 0.0
        else:
            row[f"{col}_nonmissing"] = 0
            row[f"{col}_coverage"] = 0.0
    return row


def write_exposure_coverage_report(
    input_path: str | Path,
    out_dir: str | Path,
    *,
    group_cols: list[str] | None = None,
) -> pd.DataFrame:
    df = pd.read_csv(input_path, low_memory=False)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    rows = [_coverage_row(df, "overall", "all")]
    for group_col in group_cols or ["label_source", "pdb_id", "target_id"]:
        if group_col not in df.columns:
            continue
        for group_value, group in df.groupby(group_col, dropna=False):
            rows.append(_coverage_row(group, group_col, str(group_value)))
    report = pd.DataFrame(rows)
    report.to_csv(out_path / "exposure_coverage_summary.csv", index=False)
    manifest = {
        "input_path": str(input_path),
        "exposure_columns": list(EXPOSURE_COLUMNS),
        "policy": "Missing exposure is unknown and is not converted to negative evidence.",
    }
    (out_path / "exposure_coverage_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return report
