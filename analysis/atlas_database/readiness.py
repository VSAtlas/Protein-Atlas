"""Descriptive release-readiness metrics for Atlas SQLite snapshots."""

from __future__ import annotations

import json
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from analysis.atlas_database.annotations import is_explicitly_qualified
from analysis.atlas_database.exports import (
    ranking_eligibility,
    receptor_ranking_eligibility,
)

READINESS_SCHEMA_VERSION = 3
_FAILURE_STATUSES = {
    "expected",
    "failed",
    "invalid",
    "missing_result",
    "unsuccessful",
}


def _fraction(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def _status_counts(values: Iterable[Any]) -> dict[str, int]:
    counts = Counter(str(value or "missing") for value in values)
    return dict(sorted(counts.items()))


def _final_score_source(
    result_json: Any, final_score: Any, stored_source: Any = None
) -> str:
    if final_score is None:
        return "missing_final_score"
    source = str(stored_source or "").strip()
    if source:
        return source
    try:
        payload = json.loads(str(result_json or "{}"))
    except (TypeError, ValueError):
        payload = {}
    source = str(payload.get("final_score_source") or "").strip()
    return source or "legacy_unclassified_final_score"


def _check(
    name: str,
    *,
    complete: bool,
    observed: int,
    total: int,
    missing: int,
) -> dict[str, Any]:
    return {
        "check": name,
        "status": "complete" if complete else "missing_coverage",
        "observed_count": observed,
        "total_count": total,
        "missing_count": missing,
        "coverage_fraction": _fraction(observed, total),
    }


def audit_release_readiness(
    database_path: Path,
    *,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """Query one release database and return descriptive readiness evidence."""
    connection = sqlite3.connect(database_path)
    connection.row_factory = sqlite3.Row
    try:
        release_row = connection.execute(
            "SELECT release_id FROM releases LIMIT 1"
        ).fetchone()
        if release_row is None:
            raise ValueError("database contains no release record")
        release_id = str(release_row["release_id"])
        pairs = [
            dict(row)
            for row in connection.execute(
                """SELECT p.pair_cell_id, p.final_status, p.failure_code,
                p.failure_reason, p.expected, p.has_result, p.pose_valid,
                p.is_control, p.is_decoy, p.final_score,
                p.final_score_source, p.result_json, r.variant,
                a.receptor_classification,
                a.qualification_status AS receptor_qualification_status,
                COALESCE(a.native_redock_status, r.native_redock_status)
                    AS native_redock_status
                FROM pair_cells p
                JOIN receptor_contexts r
                  ON r.receptor_context_id=p.receptor_context_id
                LEFT JOIN receptor_annotations a
                  ON a.receptor_context_id=r.receptor_context_id"""
            )
        ]
        receptors = [
            dict(row)
            for row in connection.execute(
                """SELECT r.receptor_context_id, r.variant,
                a.receptor_classification,
                a.qualification_status AS receptor_qualification_status,
                COALESCE(a.native_redock_status, r.native_redock_status)
                    AS native_redock_status,
                COALESCE(a.native_redock_reason, r.native_redock_reason)
                    AS native_redock_reason
                FROM receptor_contexts r
                LEFT JOIN receptor_annotations a
                  ON a.receptor_context_id=r.receptor_context_id"""
            )
        ]
        artifacts = [
            dict(row)
            for row in connection.execute(
                """SELECT artifact_id, pair_cell_id, sha256, verified
                FROM artifacts"""
            )
        ]
    finally:
        connection.close()

    pair_count = len(pairs)
    expected_pairs = [row for row in pairs if row["expected"] == 1]
    nonterminal_expected = [
        row for row in expected_pairs if row["final_status"] == "expected"
    ]
    failure_rows = [
        row
        for row in pairs
        if row["failure_code"]
        or str(row["final_status"] or "").lower() in _FAILURE_STATUSES
    ]
    failures_with_reason = [row for row in failure_rows if row["failure_reason"]]
    scored = [row for row in pairs if row["final_score"] is not None]
    valid_drug_candidates = [
        row
        for row in pairs
        if row["pose_valid"] == 1 and receptor_ranking_eligibility(row)[0] == 1
    ]
    scored_valid_drug_candidates = [
        row for row in valid_drug_candidates if row["final_score"] is not None
    ]
    ranking_results = [ranking_eligibility(row) for row in pairs]
    rank_eligible = [
        row for row, (eligible, _) in zip(pairs, ranking_results) if eligible == 1
    ]
    ranking_reason_counts = _status_counts(
        reason or "eligible" for _, reason in ranking_results
    )
    pose_valid = [row for row in pairs if row["pose_valid"] == 1]
    pose_invalid = [row for row in pairs if row["pose_valid"] == 0]
    pose_unvalidated = [row for row in pairs if row["pose_valid"] is None]

    quality_status_counts = _status_counts(
        row["receptor_qualification_status"] for row in receptors
    )
    redock_status_counts = _status_counts(
        row["native_redock_status"] for row in receptors
    )
    classification_policy_pending = [
        row for row in receptors if str(row["receptor_classification"] or "").strip()
    ]
    headline_receptors = [
        row
        for row in receptors
        if receptor_ranking_eligibility(row)[1] != "apo_receptor"
    ]
    explicitly_quality_qualified = [
        row
        for row in headline_receptors
        if is_explicitly_qualified(row["receptor_qualification_status"])
    ]
    explicitly_qualified = [
        row for row in receptors if is_explicitly_qualified(row["native_redock_status"])
    ]
    variant_counts = {"HOLO": 0, "APO": 0, "other": 0}
    for row in receptors:
        variant = str(row["variant"] or "").strip().upper()
        variant_counts[variant if variant in {"HOLO", "APO"} else "other"] += 1

    artifact_count = len(artifacts)
    artifacts_hashed = [row for row in artifacts if str(row["sha256"] or "").strip()]
    artifacts_verified = [row for row in artifacts if row["verified"] == 1]
    pairs_with_artifacts = {
        int(row["pair_cell_id"]) for row in artifacts if row["pair_cell_id"] is not None
    }
    source_counts = _status_counts(
        _final_score_source(
            row["result_json"], row["final_score"], row["final_score_source"]
        )
        for row in pairs
    )
    classified_sources = sum(
        count
        for source, count in source_counts.items()
        if source not in {"missing_final_score", "legacy_unclassified_final_score"}
    )

    checklist = [
        _check(
            "expected_pair_cells_terminal",
            complete=not nonterminal_expected,
            observed=len(expected_pairs) - len(nonterminal_expected),
            total=len(expected_pairs),
            missing=len(nonterminal_expected),
        ),
        _check(
            "failure_reason_coverage",
            complete=len(failures_with_reason) == len(failure_rows),
            observed=len(failures_with_reason),
            total=len(failure_rows),
            missing=len(failure_rows) - len(failures_with_reason),
        ),
        _check(
            "valid_drug_final_score_coverage",
            complete=(len(scored_valid_drug_candidates) == len(valid_drug_candidates)),
            observed=len(scored_valid_drug_candidates),
            total=len(valid_drug_candidates),
            missing=(len(valid_drug_candidates) - len(scored_valid_drug_candidates)),
        ),
        _check(
            "pose_validation_coverage",
            complete=not pose_unvalidated,
            observed=len(pose_valid) + len(pose_invalid),
            total=pair_count,
            missing=len(pose_unvalidated),
        ),
        _check(
            "receptor_classification_policy_approval",
            complete=not classification_policy_pending,
            observed=len(receptors) - len(classification_policy_pending),
            total=len(receptors),
            missing=len(classification_policy_pending),
        ),
        _check(
            "receptor_quality_explicit_qualification",
            complete=len(explicitly_quality_qualified) == len(headline_receptors),
            observed=len(explicitly_quality_qualified),
            total=len(headline_receptors),
            missing=len(headline_receptors) - len(explicitly_quality_qualified),
        ),
        _check(
            "native_redock_explicit_qualification",
            complete=len(explicitly_qualified) == len(receptors),
            observed=len(explicitly_qualified),
            total=len(receptors),
            missing=len(receptors) - len(explicitly_qualified),
        ),
        _check(
            "pair_artifact_coverage",
            complete=len(pairs_with_artifacts) == pair_count,
            observed=len(pairs_with_artifacts),
            total=pair_count,
            missing=pair_count - len(pairs_with_artifacts),
        ),
        _check(
            "artifact_hash_coverage",
            complete=len(artifacts_hashed) == artifact_count,
            observed=len(artifacts_hashed),
            total=artifact_count,
            missing=artifact_count - len(artifacts_hashed),
        ),
        _check(
            "artifact_verification_coverage",
            complete=len(artifacts_verified) == artifact_count,
            observed=len(artifacts_verified),
            total=artifact_count,
            missing=artifact_count - len(artifacts_verified),
        ),
        _check(
            "final_score_source_classification",
            complete=classified_sources == len(scored),
            observed=classified_sources,
            total=len(scored),
            missing=len(scored) - classified_sources,
        ),
    ]
    blockers = [
        {
            "code": item["check"],
            "missing_count": item["missing_count"],
            "message": (
                f"{item['missing_count']} of {item['total_count']} records lack "
                f"coverage for {item['check']}"
            ),
        }
        for item in checklist
        if item["missing_count"] > 0
    ]
    result = {
        "readiness_schema_version": READINESS_SCHEMA_VERSION,
        "release_id": release_id,
        "database_file": database_path.name,
        "metrics": {
            "pair_cells": {
                "count": pair_count,
                "expected_count": len(expected_pairs),
                "nonterminal_expected_count": len(nonterminal_expected),
                "status_counts": _status_counts(row["final_status"] for row in pairs),
                "failure_row_count": len(failure_rows),
                "failure_reason_present_count": len(failures_with_reason),
                "failure_reason_missing_count": (
                    len(failure_rows) - len(failures_with_reason)
                ),
            },
            "scores": {
                "final_score_present_count": len(scored),
                "final_score_missing_count": pair_count - len(scored),
                "final_score_coverage_fraction": _fraction(len(scored), pair_count),
                "rank_eligible_count": len(rank_eligible),
                "rank_eligible_fraction": _fraction(len(rank_eligible), pair_count),
                "ranking_eligibility_reason_counts": ranking_reason_counts,
                "final_score_source_counts": source_counts,
                "classified_final_score_source_count": classified_sources,
                "valid_drug_candidate_count": len(valid_drug_candidates),
                "valid_drug_candidate_with_final_score_count": len(
                    scored_valid_drug_candidates
                ),
                "valid_drug_final_score_coverage_fraction": _fraction(
                    len(scored_valid_drug_candidates),
                    len(valid_drug_candidates),
                ),
            },
            "pose_validation": {
                "valid_count": len(pose_valid),
                "invalid_count": len(pose_invalid),
                "unvalidated_count": len(pose_unvalidated),
            },
            "receptor_quality": {
                "receptor_classification": {
                    "policy": "pending_user_approval",
                    "recorded_but_unapproved_count": len(classification_policy_pending),
                    "ranking_behavior": "excluded_without_free_text_inference",
                },
                "receptor_count": len(receptors),
                "headline_non_apo_receptor_count": len(headline_receptors),
                "status_counts": quality_status_counts,
                "explicitly_qualified_count": len(explicitly_quality_qualified),
                "not_explicitly_qualified_count": (
                    len(headline_receptors) - len(explicitly_quality_qualified)
                ),
            },
            "native_redocking": {
                "receptor_count": len(receptors),
                "status_counts": redock_status_counts,
                "explicitly_qualified_count": len(explicitly_qualified),
                "not_explicitly_qualified_count": (
                    len(receptors) - len(explicitly_qualified)
                ),
            },
            "receptor_variants": variant_counts,
            "artifacts": {
                "artifact_count": artifact_count,
                "hashed_count": len(artifacts_hashed),
                "verified_count": len(artifacts_verified),
                "pair_with_artifact_count": len(pairs_with_artifacts),
                "pair_artifact_coverage_fraction": _fraction(
                    len(pairs_with_artifacts), pair_count
                ),
                "hash_coverage_fraction": _fraction(
                    len(artifacts_hashed), artifact_count
                ),
                "verified_coverage_fraction": _fraction(
                    len(artifacts_verified), artifact_count
                ),
            },
        },
        "checklist": checklist,
        "blockers": blockers,
        "blocker_count": len(blockers),
    }
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(result, indent=2, sort_keys=True), encoding="utf-8"
        )
    return result
