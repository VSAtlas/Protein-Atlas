"""Explicit provenance contract for Atlas pose-validity evidence."""

from __future__ import annotations

import json
from typing import Any, Mapping

LEGACY_POSE_VALIDATION_METHOD = "atlas_posebusters_any_stage_legacy_v1"
LEGACY_POSE_VALIDATION_SCOPE = "any_stage_for_ligand"
QUALIFYING_POSE_VALIDATION_SCOPE = "selected_final_score_pose"
POSEBUSTERS_RELATIVE_DISTANCE_CLI_DEFAULT = 0.92
POSEBUSTERS_REQUIRED_CHECKS = (
    "mol_pred_loaded",
    "sanitization",
    "passes_valence_checks",
    "internal_steric_clash",
    "bond_lengths",
    "bond_angles",
)


def legacy_pose_validation_thresholds() -> dict[str, Any]:
    """Describe the legacy master-row predicate without claiming run config."""
    return {
        "required_boolean_checks": list(POSEBUSTERS_REQUIRED_CHECKS),
        "protein_clash_rule": {
            "either": [
                "most_extreme_clash_protein",
                "most_extreme_relative_distance_protein >= relative_distance_cutoff",
            ]
        },
        "relative_distance_cutoff": {
            "cli_default": POSEBUSTERS_RELATIVE_DISTANCE_CLI_DEFAULT,
            "exact_run_value": None,
            "provenance": "tools/pose_bust.py default; absent from legacy master rows",
        },
    }


def legacy_pose_validation_thresholds_json() -> str:
    """Return deterministic JSON for the legacy PoseBusters predicate."""
    return json.dumps(
        legacy_pose_validation_thresholds(),
        sort_keys=True,
        separators=(",", ":"),
    )


def pose_validation_browser_contract(
    scientific_policies: Mapping[str, Any],
) -> dict[str, Any]:
    """Expose what the legacy flag proves and the decision it cannot prove."""
    explicit_policy = scientific_policies.get("pose_validation")
    return {
        "stored_pass_field": "pose_valid",
        "legacy_source_field": "pose_valid_any",
        "legacy_method": LEGACY_POSE_VALIDATION_METHOD,
        "legacy_scope": LEGACY_POSE_VALIDATION_SCOPE,
        "ranking_required_scope": QUALIFYING_POSE_VALIDATION_SCOPE,
        "aggregation": "true when any PoseBusters stage row passes for a ligand",
        "thresholds": legacy_pose_validation_thresholds(),
        "exact_selected_final_score_pose_alignment": "unverified",
        "release_policy": explicit_policy
        if isinstance(explicit_policy, Mapping)
        else None,
    }


__all__ = [
    "LEGACY_POSE_VALIDATION_METHOD",
    "LEGACY_POSE_VALIDATION_SCOPE",
    "POSEBUSTERS_RELATIVE_DISTANCE_CLI_DEFAULT",
    "POSEBUSTERS_REQUIRED_CHECKS",
    "QUALIFYING_POSE_VALIDATION_SCOPE",
    "legacy_pose_validation_thresholds",
    "legacy_pose_validation_thresholds_json",
    "pose_validation_browser_contract",
]
