from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.split_manifest import dataframe_content_hash


def _label_counts(frame: pd.DataFrame, label_col: str) -> dict[str, int]:
    if label_col not in frame.columns:
        return {}
    values = pd.to_numeric(frame[label_col], errors="coerce")
    return {str(key): int(value) for key, value in values.value_counts(dropna=False).items()}


def write_dataset_version_manifest(
    *,
    dataset_path: str | Path,
    frame: pd.DataFrame,
    out_path: str | Path,
    label_col: str,
    feature_set: str | None,
    features: list[str] | None = None,
    requested_features: list[str] | None = None,
    missing_features: list[str] | None = None,
    all_missing_features: list[str] | None = None,
    exclude_features: list[str] | None = None,
    retained_excluded_features: list[str] | None = None,
    retained_context_not_trained: list[str] | None = None,
    provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(frame.columns)
    manifest: dict[str, Any] = {
        "dataset_path": str(dataset_path),
        "dataset_version_id": dataframe_content_hash(frame),
        "n_rows": int(len(frame)),
        "n_columns": int(len(columns)),
        "columns": columns,
        "label_col": label_col,
        "label_counts": _label_counts(frame, label_col),
        "feature_set": feature_set,
        "features": list(features or []),
        "requested_features": list(requested_features or []),
        "missing_features": list(missing_features or []),
        "all_missing_features": list(all_missing_features or []),
        "feature_coverage_fraction": (
            len(features or []) / len(requested_features or []) if requested_features else None
        ),
        "exclude_features": list(exclude_features or []),
        "retained_excluded_features": list(retained_excluded_features or []),
        "retained_excluded_policy": (
            "Columns listed here were present in the dataset and retained for audit/reporting, "
            "but excluded from the model design matrix."
            if retained_excluded_features
            else "No present dataset columns were explicitly excluded from the model design matrix."
        ),
        "retained_context_not_trained": list(retained_context_not_trained or []),
        "retained_context_not_trained_policy": (
            "High-risk/context columns listed here were present in the dataset but were not selected "
            "as model features by the active feature set."
            if retained_context_not_trained
            else "No audited high-risk/context columns were present outside the active feature set."
        ),
        "source_columns_present": [
            col
            for col in ["label_source", "source_family", "upstream_source", "assay_type", "endpoint_type"]
            if col in frame.columns
        ],
        "split_metadata_columns_present": [
            col
            for col in [
                "drug_id",
                "target_id",
                "scaffold_key",
                "chemical_cluster",
                "target_family",
                "protein_class",
                "activity_publication_year",
                "database_release_year",
            ]
            if col in frame.columns
        ],
        "provenance": provenance or {},
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
