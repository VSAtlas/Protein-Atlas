from __future__ import annotations

from typing import Sequence

from analysis.ml.feature_sets import FORBIDDEN_FEATURES
from analysis.ml.feature_sets import SPD_EXPOSURE_DEFINITION_FEATURES

FORBIDDEN_AVAILABILITY_FEATURES = {
    "rescored_flag",
    "rescored_stage",
    "cnn_rescored_flag",
    "rescore_available",
    "scorch_available",
    "SCORCH_available",
}

FORBIDDEN_PARTIAL_RESCORING_FEATURES = {
    "SCORCH_score_used",
    "SCORCH_certainty_used",
    "scorch_composite",
    "scorch_pct",
    "scorch_source_used",
    "final_rank",
    "pilot_SCORCH_score_used",
    "pilot_final_score",
}

FORBIDDEN_EXPOSURE_LABEL_FEATURES = set(SPD_EXPOSURE_DEFINITION_FEATURES)

FORBIDDEN_SOURCE_LABEL_FEATURES = {
    "source_objective",
    "label_source",
    "source_label_policy",
    "source_conflict",
    "label_source_count",
    "source_family",
    "upstream_source",
    "parent_sources",
    "source_release",
    "source_version",
    "assay_type",
    "assay_mode",
    "endpoint_type",
    "activity_type",
    "standard_type",
}


def find_leaky_features(
    features: Sequence[str],
    *,
    allow_label_definition_features: bool = False,
    allow_partial_rescoring_features: bool = False,
) -> list[str]:
    exposure_definition_features = (
        set()
        if allow_label_definition_features
        else FORBIDDEN_EXPOSURE_LABEL_FEATURES
    )
    partial_rescoring_features = (
        set()
        if allow_partial_rescoring_features
        else FORBIDDEN_PARTIAL_RESCORING_FEATURES
    )
    leaky = [
        feature
        for feature in features
        if feature in FORBIDDEN_FEATURES
        or feature in FORBIDDEN_AVAILABILITY_FEATURES
        or feature in partial_rescoring_features
        or feature in exposure_definition_features
        or feature in FORBIDDEN_SOURCE_LABEL_FEATURES
        or feature.lower().startswith("uniprot_")
    ]
    leaky.extend(feature for feature in features if feature.endswith("_label") and feature != "predicted_label")
    leaky.extend(feature for feature in features if feature.endswith("_label_status"))
    leaky.extend(feature for feature in features if feature.endswith("_threshold_nM"))
    leaky.extend(feature for feature in features if feature.endswith("_assay_ids"))
    leaky.extend(feature for feature in features if feature.endswith("_assay_count"))
    leaky.extend(feature for feature in features if feature.endswith("_activity_nM"))
    leaky.extend(
        feature
        for feature in features
        if "rescored" in feature.lower() or feature.lower().endswith("_available")
    )
    return sorted(set(leaky))


def assert_no_leakage(
    features: Sequence[str],
    *,
    allow_label_definition_features: bool = False,
    allow_partial_rescoring_features: bool = False,
) -> None:
    leaky = find_leaky_features(
        features,
        allow_label_definition_features=allow_label_definition_features,
        allow_partial_rescoring_features=allow_partial_rescoring_features,
    )
    if leaky:
        raise ValueError(f"leaky predictive features are not allowed: {', '.join(leaky)}")
