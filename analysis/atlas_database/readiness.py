"""Descriptive release-readiness metrics for Atlas SQLite snapshots."""

from __future__ import annotations

import json
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from analysis.atlas_database.annotations import is_explicitly_qualified
from analysis.atlas_database.exports import (
    apo_exploratory_ranking_eligibility,
    ranking_eligibility,
    receptor_chemistry_evidence_eligibility,
    receptor_ranking_eligibility,
)
from analysis.atlas_database.pose_validation_contract import (
    QUALIFYING_POSE_VALIDATION_SCOPE,
)
from analysis.atlas_database.score_source_contract import (
    classify_materialized_source,
    score_source_browser_contract,
)

READINESS_SCHEMA_VERSION = 6
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


def _has_structured_evidence(value: Any) -> bool:
    try:
        payload = json.loads(str(value or ""))
    except (TypeError, ValueError):
        return False
    return isinstance(payload, dict) and bool(payload)


def _declared_final_score_source(
    result_json: Any, final_score: Any, stored_source: Any = None
) -> str | None:
    if final_score is None:
        return None
    source = str(stored_source or "").strip()
    if source:
        return source
    try:
        payload = json.loads(str(result_json or "{}"))
    except (TypeError, ValueError):
        payload = {}
    source = str(payload.get("final_score_source") or "").strip()
    return source or None


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
                p.pose_validation_scope, p.is_control, p.is_decoy,
                p.final_score,
                p.final_score_source, p.final_score_source_reconstructed,
                p.final_score_source_classification,
                p.final_score_source_evidence_json,
                sra.input_csv_sha256 AS selected_result_input_csv_sha256,
                sra.source_row_number AS selected_result_source_row_number,
                sra.result_sha256 AS selected_result_sha256,
                p.result_json, r.variant,
                a.receptor_classification,
                a.classification_method AS receptor_classification_method,
                a.chemistry_evidence_status,
                a.prepared_receptor_sha256,
                a.canonical_chemistry_policy_sha256,
                a.requested_observed_conflict,
                a.qualification_status AS receptor_qualification_status,
                COALESCE(a.native_redock_status, r.native_redock_status)
                    AS native_redock_status
                FROM pair_cells p
                JOIN receptor_contexts r
                  ON r.receptor_context_id=p.receptor_context_id
                LEFT JOIN receptor_annotations a
                  ON a.receptor_context_id=r.receptor_context_id
                LEFT JOIN result_attempts sra
                  ON sra.pair_cell_id=p.pair_cell_id
                 AND sra.selected_for_release=1"""
            )
        ]
        receptors = [
            dict(row)
            for row in connection.execute(
                """SELECT r.receptor_context_id, r.variant,
                a.receptor_classification,
                a.chemistry_evidence_status,
                a.classification_method AS receptor_classification_method,
                a.prepared_receptor_sha256,
                a.canonical_chemistry_policy_sha256,
                a.requested_observed_conflict,
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
        receptor_audits = [
            dict(row)
            for row in connection.execute(
                """SELECT audit_kind, parse_status FROM receptor_audits
                ORDER BY receptor_audit_id"""
            )
        ]
        selected_result_attempts = [
            dict(row)
            for row in connection.execute(
                """SELECT result_attempt_id, completion_record_id,
                completion_link_method, completion_link_evidence_json,
                CASE WHEN EXISTS (
                    SELECT 1 FROM docking_attempts d
                    WHERE d.pair_cell_id=result_attempts.pair_cell_id
                      AND d.completion_record_id=
                          result_attempts.completion_record_id
                      AND d.selected_for_release=1
                ) THEN 1 ELSE 0 END AS selected_completion_match
                FROM result_attempts
                WHERE selected_for_release=1"""
            )
        ]
    finally:
        connection.close()

    pair_count = len(pairs)
    for row in pairs:
        row["final_score_source"] = _declared_final_score_source(
            row["result_json"], row["final_score"], row["final_score_source"]
        )
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
        if row["pose_valid"] == 1
        and row["pose_validation_scope"] == QUALIFYING_POSE_VALIDATION_SCOPE
        and receptor_ranking_eligibility(row)[0] == 1
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
    apo_ranking_results = [apo_exploratory_ranking_eligibility(row) for row in pairs]
    apo_rank_eligible = [
        row for row, (eligible, _) in zip(pairs, apo_ranking_results) if eligible == 1
    ]
    apo_ranking_reason_counts = _status_counts(
        reason or "eligible" for _, reason in apo_ranking_results
    )
    pose_valid = [row for row in pairs if row["pose_valid"] == 1]
    pose_invalid = [row for row in pairs if row["pose_valid"] == 0]
    pose_unvalidated = [row for row in pairs if row["pose_valid"] is None]
    pose_selected_scope = [
        row
        for row in pose_valid
        if row["pose_validation_scope"] == QUALIFYING_POSE_VALIDATION_SCOPE
    ]
    pose_scope_unqualified = [
        row
        for row in pose_valid
        if row["pose_validation_scope"] != QUALIFYING_POSE_VALIDATION_SCOPE
    ]

    quality_status_counts = _status_counts(
        row["receptor_qualification_status"] for row in receptors
    )
    redock_status_counts = _status_counts(
        row["native_redock_status"] for row in receptors
    )
    controlled_classifications = [
        row
        for row in receptors
        if str(row["receptor_classification"] or "").strip().upper() in {"APO", "HOLO"}
    ]
    classification_missing_or_uncontrolled = [
        row for row in receptors if row not in controlled_classifications
    ]
    classification_counts = _status_counts(
        str(row["receptor_classification"] or "").strip().upper() or "missing"
        for row in receptors
    )
    chemistry_evidence_status_counts = _status_counts(
        row["chemistry_evidence_status"] for row in receptors
    )
    chemistry_conflicts = [
        row for row in receptors if row["requested_observed_conflict"] == 1
    ]
    chemistry_evidence_results = [
        receptor_chemistry_evidence_eligibility(row) for row in receptors
    ]
    chemistry_evidence_qualified = [
        row
        for row, (eligible, _) in zip(receptors, chemistry_evidence_results)
        if eligible == 1
    ]
    headline_receptors = [
        row
        for row in receptors
        if str(row["receptor_classification"] or "").strip().upper() == "HOLO"
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
    scored_source_classifications = [
        classify_materialized_source(row) for row in scored
    ]
    source_counts = _status_counts(row["final_score_source"] for row in scored)
    effective_source_counts = _status_counts(
        item.effective_source for item in scored_source_classifications
    )
    source_classification_counts = _status_counts(
        item.classification for item in scored_source_classifications
    )
    classified_sources = sum(item.eligible for item in scored_source_classifications)
    pose_linkage_requirement_counts = _status_counts(
        item.pose_linkage_requirement for item in scored_source_classifications
    )
    score_pose_linkage_unresolved = sum(
        item.eligible
        and bool(item.pose_linkage_requirement)
        and (
            "no_unique" in str(item.pose_linkage_requirement)
            or "no_immutable_pose" in str(item.pose_linkage_requirement)
            or "requires_immutable_pose_identifier"
            in str(item.pose_linkage_requirement)
        )
        for item in scored_source_classifications
    )
    result_attempts_with_completion_lineage = [
        row
        for row in selected_result_attempts
        if row["completion_record_id"] is not None
        and str(row["completion_link_method"] or "").strip()
        and _has_structured_evidence(row["completion_link_evidence_json"])
        and row["selected_completion_match"] == 1
    ]

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
            "pose_validation_selected_pose_alignment",
            complete=not pose_scope_unqualified,
            observed=len(pose_selected_scope),
            total=len(pose_valid),
            missing=len(pose_scope_unqualified),
        ),
        _check(
            "receptor_classification_controlled_coverage",
            complete=not classification_missing_or_uncontrolled,
            observed=len(controlled_classifications),
            total=len(receptors),
            missing=len(classification_missing_or_uncontrolled),
        ),
        _check(
            "receptor_chemistry_exact_evidence",
            complete=len(chemistry_evidence_qualified) == len(receptors),
            observed=len(chemistry_evidence_qualified),
            total=len(receptors),
            missing=len(receptors) - len(chemistry_evidence_qualified),
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
        _check(
            "final_score_pose_linkage_resolution",
            complete=score_pose_linkage_unresolved == 0,
            observed=len(scored) - score_pose_linkage_unresolved,
            total=len(scored),
            missing=score_pose_linkage_unresolved,
        ),
        _check(
            "selected_result_completion_lineage",
            complete=(
                len(result_attempts_with_completion_lineage)
                == len(selected_result_attempts)
            ),
            observed=len(result_attempts_with_completion_lineage),
            total=len(selected_result_attempts),
            missing=(
                len(selected_result_attempts)
                - len(result_attempts_with_completion_lineage)
            ),
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
                "apo_exploratory_rank_eligible_count": len(apo_rank_eligible),
                "apo_exploratory_rank_eligible_fraction": _fraction(
                    len(apo_rank_eligible), pair_count
                ),
                "apo_exploratory_ranking_eligibility_reason_counts": (
                    apo_ranking_reason_counts
                ),
                "final_score_source_counts": source_counts,
                "effective_final_score_source_counts": effective_source_counts,
                "final_score_source_classification_counts": (
                    source_classification_counts
                ),
                "classified_final_score_source_count": classified_sources,
                "final_score_pose_linkage_requirement_counts": (
                    pose_linkage_requirement_counts
                ),
                "final_score_pose_linkage_unresolved_count": (
                    score_pose_linkage_unresolved
                ),
                "final_score_source_contract": score_source_browser_contract(),
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
                "selected_pose_scope_count": len(pose_selected_scope),
                "scope_unqualified_count": len(pose_scope_unqualified),
            },
            "receptor_quality": {
                "receptor_classification": {
                    "vocabulary": ["APO", "HOLO"],
                    "classification_counts": classification_counts,
                    "controlled_count": len(controlled_classifications),
                    "missing_or_uncontrolled_count": len(
                        classification_missing_or_uncontrolled
                    ),
                    "chemistry_evidence_status_counts": (
                        chemistry_evidence_status_counts
                    ),
                    "requested_observed_conflict_count": len(chemistry_conflicts),
                    "exact_evidence_qualified_count": len(
                        chemistry_evidence_qualified
                    ),
                    "exact_evidence_missing_count": (
                        len(receptors) - len(chemistry_evidence_qualified)
                    ),
                    "ranking_behavior": "separate exact-evidence HOLO and APO tracks",
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
            "receptor_audits": {
                "count": len(receptor_audits),
                "kind_counts": _status_counts(
                    row["audit_kind"] for row in receptor_audits
                ),
                "status_counts": _status_counts(
                    row["parse_status"] for row in receptor_audits
                ),
                "accepted_count": sum(
                    row["parse_status"] == "parsed" for row in receptor_audits
                ),
                "rejected_count": sum(
                    str(row["parse_status"]).startswith("rejected:")
                    for row in receptor_audits
                ),
                "qualification_behavior": (
                    "descriptive evidence only; audit presence never qualifies "
                    "receptor quality or native redocking"
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
            "attempt_lineage": {
                "selected_result_attempt_count": len(selected_result_attempts),
                "completion_linked_count": len(result_attempts_with_completion_lineage),
                "completion_link_missing_count": (
                    len(selected_result_attempts)
                    - len(result_attempts_with_completion_lineage)
                ),
                "policy": (
                    "a selected result must name its causal completion record, "
                    "link method, structured evidence, and matching selected "
                    "completion attempt"
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
