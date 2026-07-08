from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.audit_utils import load_table
from analysis.ml.feature_sets import effective_exclude_features, get_feature_set


def _maybe_features(dataset: pd.DataFrame, label_col: str, feature_set: str | None, exclude_features: list[str] | None) -> list[str]:
    if not feature_set:
        return []
    excluded = effective_exclude_features(label_col, exclude_features)
    return [feature for feature in get_feature_set(feature_set) if feature in dataset.columns and feature not in excluded]


def write_ml_data_card(
    dataset_path: str | Path,
    label_col: str,
    feature_set: str | None,
    exclude_features: list[str] | None,
    out_dir: Path,
) -> dict[str, Any]:
    data = load_table(dataset_path)
    features = _maybe_features(data, label_col, feature_set, exclude_features)
    label = pd.to_numeric(data[label_col], errors="coerce") if label_col in data.columns else pd.Series(dtype=float)
    observed_label = label.dropna()
    numeric_features = [feature for feature in features if feature in data.columns and pd.api.types.is_numeric_dtype(data[feature])]
    year_cols = [
        col
        for col in [
            "source_available_date",
            "source_release_date",
            "document_date",
            "document_year",
            "activity_publication_year",
            "database_release_year",
            "label_publication_year",
            "evidence_publication_year",
        ]
        if col in data.columns
    ]
    card = {
        "dataset_path": str(dataset_path),
        "n_rows": int(len(data)),
        "n_columns": int(len(data.columns)),
        "label_col": label_col,
        "n_label_observed": int(observed_label.shape[0]) if len(label) else 0,
        "n_positive": int((label == 1).sum()) if len(label) else 0,
        "n_negative": int((label == 0).sum()) if len(label) else 0,
        "positive_rate": float((observed_label == 1).mean()) if len(observed_label) else None,
        "n_drugs": int(data["drug_id"].nunique()) if "drug_id" in data.columns else None,
        "n_targets": int(data["target_id"].nunique()) if "target_id" in data.columns else None,
        "n_scaffolds": int(data["scaffold_key"].nunique()) if "scaffold_key" in data.columns else None,
        "n_chemical_clusters": int(data["chemical_cluster"].nunique()) if "chemical_cluster" in data.columns else None,
        "n_target_families": int(data["target_family"].nunique()) if "target_family" in data.columns else None,
        "n_source_families": int(data["source_family"].nunique()) if "source_family" in data.columns else None,
        "temporal_metadata_columns": year_cols,
        "has_temporal_metadata": bool(year_cols),
        "feature_set": feature_set,
        "features_available": features,
        "numeric_feature_count": len(numeric_features),
        "recommended_reporting": [
            "Report AUROC together with AUPRC, calibration, enrichment@k, and label prevalence.",
            "Treat random split as an optimistic smoke test; prefer drug, target, scaffold, source, and temporal holdouts.",
            "Use chemical-cluster and target-family holdouts for OOD claims when metadata are available.",
            "Run train-vs-benchmark independence audits before claiming external validation.",
            "Keep unlabeled rows unknown unless a source provides measured or reliable negative evidence.",
            "Report applicability-domain and conformal/uncertainty outputs for predictions that will guide follow-up experiments.",
        ],
    }
    path = out_dir / "ml_data_card.json"
    path.write_text(json.dumps(card, indent=2, sort_keys=True), encoding="utf-8")
    return card
