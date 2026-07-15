"""Publication-safe browser and tabular exports for Atlas release databases."""

from __future__ import annotations

import csv
import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

from analysis.atlas_database.annotations import is_explicitly_qualified
from analysis.atlas_database.artifact_contract import (
    ARTIFACT_PUBLICATION_POLICIES,
    ARTIFACT_ROLES,
    artifact_public_record_allowed,
    artifact_publication_policy,
)
from analysis.atlas_database.pose_validation_contract import (
    QUALIFYING_POSE_VALIDATION_SCOPE,
    pose_validation_browser_contract,
)
from analysis.atlas_database.receptor_evidence import (
    CHEMISTRY_CLASSIFICATION_METHOD,
)
from analysis.atlas_database.score_source_contract import (
    classify_materialized_source,
    score_source_browser_contract,
)

_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")


PAIR_COLUMNS = (
    "pair_cell_id",
    "run_id",
    "pdb_id",
    "variant",
    "ph_label",
    "target_key",
    "target_id",
    "receptor_label",
    "receptor_classification",
    "receptor_classification_method",
    "chemistry_evidence_status",
    "prepared_receptor_sha256",
    "canonical_chemistry_policy_sha256",
    "retained_metal_atom_count",
    "retained_cofactor_residue_count",
    "requested_observed_conflict",
    "receptor_qualification_status",
    "native_redock_status",
    "ligand_id",
    "ligand_canonical_id",
    "drug_id",
    "ligand_display_name",
    "final_status",
    "failure_code",
    "failure_reason",
    "expected",
    "has_result",
    "pose_valid",
    "pose_validation_method",
    "pose_validation_scope",
    "pose_validation_thresholds_json",
    "is_control",
    "is_decoy",
    "final_score",
    "final_score_source",
    "final_score_source_reconstructed",
    "final_score_source_effective",
    "final_score_source_canonical",
    "final_score_source_family",
    "final_score_pose_linkage_requirement",
    "final_score_source_classification",
    "final_score_source_evidence_json",
    "final_score_source_eligible",
    "primary_score_present",
    "rank_eligible",
    "ranking_eligibility_reason",
    "apo_exploratory_rank_eligible",
    "apo_exploratory_ranking_eligibility_reason",
    "ranking_track",
    "final_rank",
    "rank_within_receptor",
    "rank_across_receptors",
    "apo_rank_within_receptor",
    "apo_rank_across_receptors",
    "atlas_score",
    "atlas_score_source",
    "selected_docking_score",
    "consensus_score",
    "score_component_coverage",
    "artifact_count",
    "verified_artifact_count",
)


def _receptor_label(row: Mapping[str, Any]) -> str:
    return "|".join(
        value
        for value in (
            str(row.get("pdb_id") or ""),
            str(row.get("variant") or ""),
            str(row.get("ph_label") or ""),
        )
        if value
    )


def _target_key(row: Mapping[str, Any]) -> str:
    return "|".join(
        value for value in (str(row.get("run_id") or ""), _receptor_label(row)) if value
    )


def _component_coverage(row: Mapping[str, Any]) -> str:
    names = (
        "final_score",
        "atlas_score",
        "selected_docking_score",
        "consensus_score",
    )
    return ",".join(name for name in names if row.get(name) is not None)


def _declared_final_score_source(
    result_json: Any, final_score: Any, stored_source: Any = None
) -> str | None:
    if final_score is None:
        return None
    source = str(stored_source or "").strip()
    if source:
        return source
    try:
        result = json.loads(str(result_json or "{}"))
    except (TypeError, ValueError):
        result = {}
    source = str(result.get("final_score_source") or "").strip()
    return source or None


def receptor_chemistry_evidence_eligibility(
    row: Mapping[str, Any],
) -> tuple[int, str | None]:
    """Require exact, policy-bound prepared-receptor chemistry evidence."""
    if str(row.get("receptor_classification_method") or "").strip() != (
        CHEMISTRY_CLASSIFICATION_METHOD
    ):
        return 0, "receptor_chemistry_method_unqualified"
    if str(row.get("chemistry_evidence_status") or "").strip().lower() != "observed":
        return 0, "receptor_chemistry_evidence_unresolved"
    if not _SHA256.fullmatch(str(row.get("prepared_receptor_sha256") or "").strip()):
        return 0, "prepared_receptor_hash_missing"
    if not _SHA256.fullmatch(
        str(row.get("canonical_chemistry_policy_sha256") or "").strip()
    ):
        return 0, "receptor_chemistry_policy_hash_missing"
    return 1, None


def _score_pose_linkage_is_unresolved(requirement: Any) -> bool:
    token = str(requirement or "").strip().lower()
    return bool(token) and (
        "no_unique" in token
        or "no_immutable_pose" in token
        or "requires_immutable_pose_identifier" in token
    )


def receptor_ranking_eligibility(
    row: Mapping[str, Any],
) -> tuple[int, str | None]:
    """Apply only the explicit receptor-level headline ranking gates."""
    if row.get("is_control") == 1:
        return 0, "native_control"
    if row.get("is_decoy") == 1:
        return 0, "decoy"
    classification = str(row.get("receptor_classification") or "").strip().upper()
    if classification == "APO":
        return 0, "apo_receptor"
    if not classification:
        return 0, "receptor_classification_missing"
    if classification != "HOLO":
        return 0, "receptor_classification_uncontrolled"
    chemistry_eligible, reason = receptor_chemistry_evidence_eligibility(row)
    if not chemistry_eligible:
        return chemistry_eligible, reason
    if not is_explicitly_qualified(row.get("receptor_qualification_status")):
        return 0, "receptor_quality_not_qualified"
    if not is_explicitly_qualified(row.get("native_redock_status")):
        return 0, "native_redock_not_qualified"
    return 1, None


def apo_receptor_ranking_eligibility(
    row: Mapping[str, Any],
) -> tuple[int, str | None]:
    """Apply receptor gates for the separate exploratory APO ranking track."""
    if row.get("is_control") == 1:
        return 0, "native_control"
    if row.get("is_decoy") == 1:
        return 0, "decoy"
    classification = str(row.get("receptor_classification") or "").strip().upper()
    if not classification:
        return 0, "receptor_classification_missing"
    if classification == "HOLO":
        return 0, "holo_receptor"
    if classification != "APO":
        return 0, "receptor_classification_uncontrolled"
    chemistry_eligible, reason = receptor_chemistry_evidence_eligibility(row)
    if not chemistry_eligible:
        return chemistry_eligible, reason
    if not is_explicitly_qualified(row.get("receptor_qualification_status")):
        return 0, "receptor_quality_not_qualified"
    if not is_explicitly_qualified(row.get("native_redock_status")):
        return 0, "native_redock_not_qualified"
    return 1, None


def _pair_ranking_eligibility(row: Mapping[str, Any]) -> tuple[int, str | None]:
    """Apply pair gates shared by HOLO and separate APO ranking tracks."""
    if row.get("final_score") is None:
        return 0, "missing_final_score"
    if row.get("pose_valid") == 0:
        return 0, "pose_invalid"
    if row.get("pose_valid") is None:
        return 0, "pose_validation_missing"
    if row.get("pose_validation_scope") != QUALIFYING_POSE_VALIDATION_SCOPE:
        return 0, "pose_validation_scope_unqualified"
    source = classify_materialized_source(row)
    if not source.eligible:
        return 0, "final_score_source_unclassified"
    if _score_pose_linkage_is_unresolved(source.pose_linkage_requirement):
        return 0, "final_score_pose_linkage_unresolved"
    return 1, None


def ranking_eligibility(row: Mapping[str, Any]) -> tuple[int, str | None]:
    """Apply receptor qualification and pair-level HOLO headline gates."""
    receptor_eligible, reason = receptor_ranking_eligibility(row)
    if not receptor_eligible:
        return receptor_eligible, reason
    return _pair_ranking_eligibility(row)


def apo_exploratory_ranking_eligibility(
    row: Mapping[str, Any],
) -> tuple[int, str | None]:
    """Apply the same pair gates to APO without pooling APO with HOLO ranks."""
    receptor_eligible, reason = apo_receptor_ranking_eligibility(row)
    if not receptor_eligible:
        return receptor_eligible, reason
    return _pair_ranking_eligibility(row)


def _rank(
    rows: list[dict[str, Any]],
    group_name: str,
    output_name: str,
    *,
    eligibility_name: str,
) -> None:
    groups: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[row[group_name]].append(row)
    for values in groups.values():
        scored = [row for row in values if row[eligibility_name]]
        scored.sort(
            key=lambda row: (
                -float(row["final_score"]),
                str(row["ligand_canonical_id"]),
                str(row["target_key"]),
            )
        )
        previous: float | None = None
        rank = 0
        for index, row in enumerate(scored, 1):
            score = float(row["final_score"])
            if previous is None or score != previous:
                rank, previous = index, score
            row[output_name] = rank


def _pair_rows(connection: sqlite3.Connection) -> list[dict[str, Any]]:
    connection.row_factory = sqlite3.Row
    public_artifact_counts: dict[int, int] = defaultdict(int)
    for artifact in connection.execute(
        """SELECT pair_cell_id, artifact_role, artifact_json, verified, sha256
        FROM artifacts
        WHERE pair_cell_id IS NOT NULL
        ORDER BY artifact_id"""
    ):
        try:
            artifact_entry = json.loads(str(artifact["artifact_json"] or "{}"))
        except (TypeError, ValueError):
            artifact_entry = {}
        if artifact_public_record_allowed(
            str(artifact["artifact_role"]),
            artifact_entry,
            verified=artifact["verified"],
            sha256=artifact["sha256"],
        ):
            public_artifact_counts[int(artifact["pair_cell_id"])] += 1
    query = """
        SELECT p.pair_cell_id, r.run_id, r.pdb_id, r.variant, r.ph_label,
            ra.receptor_classification,
            ra.classification_method AS receptor_classification_method,
            ra.chemistry_evidence_status,
            ra.prepared_receptor_sha256,
            ra.canonical_chemistry_policy_sha256,
            ra.retained_metal_atom_count,
            ra.retained_cofactor_residue_count,
            ra.requested_observed_conflict,
            ra.qualification_status AS receptor_qualification_status,
            COALESCE(ra.native_redock_status, r.native_redock_status)
                AS native_redock_status,
            l.ligand_id, l.canonical_id AS ligand_canonical_id,
            COALESCE(l.display_name, l.canonical_id) AS ligand_display_name,
            p.final_status, p.failure_code, p.failure_reason, p.expected,
            p.has_result, p.pose_valid, p.pose_validation_method,
            p.pose_validation_scope, p.pose_validation_thresholds_json,
            p.is_control, p.is_decoy, p.final_score,
            p.final_score_source, p.final_score_source_reconstructed,
            p.final_score_source_classification,
            p.final_score_source_evidence_json,
            sra.input_csv_sha256 AS selected_result_input_csv_sha256,
            sra.source_row_number AS selected_result_source_row_number,
            sra.result_sha256 AS selected_result_sha256,
            p.final_rank, p.atlas_score, p.atlas_score_source,
            p.selected_docking_score, p.consensus_score, p.result_json
        FROM pair_cells p
        JOIN receptor_contexts r
          ON r.receptor_context_id = p.receptor_context_id
        JOIN ligands l ON l.ligand_id = p.ligand_id
        LEFT JOIN receptor_annotations ra
          ON ra.receptor_context_id = r.receptor_context_id
        LEFT JOIN result_attempts sra
          ON sra.pair_cell_id = p.pair_cell_id
         AND sra.selected_for_release = 1
        GROUP BY p.pair_cell_id
        ORDER BY r.run_id, r.pdb_id, r.variant, r.ph_label, l.canonical_id
    """
    rows: list[dict[str, Any]] = []
    for raw in connection.execute(query):
        row = dict(raw)
        row["receptor_label"] = _receptor_label(row)
        row["target_key"] = _target_key(row)
        row["target_id"] = row["target_key"]
        row["drug_id"] = row["ligand_canonical_id"]
        row["final_score_source"] = _declared_final_score_source(
            row.pop("result_json", None),
            row["final_score"],
            row.get("final_score_source"),
        )
        source = classify_materialized_source(row)
        row["final_score_source_reconstructed"] = source.reconstructed_source
        row["final_score_source_effective"] = source.effective_source
        row["final_score_source_canonical"] = source.canonical_source
        row["final_score_source_family"] = source.family
        row["final_score_pose_linkage_requirement"] = source.pose_linkage_requirement
        row["final_score_source_classification"] = source.classification
        row["final_score_source_eligible"] = int(source.eligible)
        row["primary_score_present"] = int(row["final_score"] is not None)
        (
            row["rank_eligible"],
            row["ranking_eligibility_reason"],
        ) = ranking_eligibility(row)
        (
            row["apo_exploratory_rank_eligible"],
            row["apo_exploratory_ranking_eligibility_reason"],
        ) = apo_exploratory_ranking_eligibility(row)
        row["ranking_track"] = (
            "qualified_holo"
            if row["rank_eligible"]
            else "exploratory_apo"
            if row["apo_exploratory_rank_eligible"]
            else None
        )
        row["score_component_coverage"] = _component_coverage(row)
        artifact_count = public_artifact_counts.get(int(row["pair_cell_id"]), 0)
        row["artifact_count"] = artifact_count
        row["verified_artifact_count"] = artifact_count
        row["rank_within_receptor"] = None
        row["rank_across_receptors"] = None
        row["apo_rank_within_receptor"] = None
        row["apo_rank_across_receptors"] = None
        row.pop("selected_result_input_csv_sha256", None)
        row.pop("selected_result_source_row_number", None)
        row.pop("selected_result_sha256", None)
        rows.append(row)
    _rank(
        rows,
        "target_key",
        "rank_within_receptor",
        eligibility_name="rank_eligible",
    )
    _rank(
        rows,
        "ligand_id",
        "rank_across_receptors",
        eligibility_name="rank_eligible",
    )
    _rank(
        rows,
        "target_key",
        "apo_rank_within_receptor",
        eligibility_name="apo_exploratory_rank_eligible",
    )
    _rank(
        rows,
        "ligand_id",
        "apo_rank_across_receptors",
        eligibility_name="apo_exploratory_rank_eligible",
    )
    return rows


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(PAIR_COLUMNS))
        writer.writeheader()
        writer.writerows({name: row.get(name) for name in PAIR_COLUMNS} for row in rows)


def _write_parquet(path: Path, rows: list[dict[str, Any]]) -> None:
    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "Parquet export requires pandas and a Parquet engine"
        ) from exc
    try:
        pd.DataFrame(rows, columns=list(PAIR_COLUMNS)).to_parquet(path, index=False)
    except (ImportError, ValueError) as exc:
        raise RuntimeError(
            "Parquet export requires pyarrow or fastparquet in the active environment"
        ) from exc


def _counts(rows: Iterable[Mapping[str, Any]], name: str) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        counts[str(row.get(name) or "unknown")] += 1
    return dict(sorted(counts.items()))


def _browser_payload(
    connection: sqlite3.Connection, rows: list[dict[str, Any]]
) -> dict[str, Any]:
    connection.row_factory = sqlite3.Row
    release_row = connection.execute("SELECT * FROM releases LIMIT 1").fetchone()
    if release_row is None:
        raise ValueError("database contains no release record")
    release = dict(release_row)
    raw_manifest = release.pop("manifest_json", None)
    try:
        manifest = json.loads(str(raw_manifest or "{}"))
    except (TypeError, ValueError):
        manifest = {}
    release["id"] = release["release_id"]
    release["version"] = release["release_id"]
    description = manifest.get("description") or manifest.get("summary")
    if isinstance(description, str) and description.strip():
        release["summary"] = description.strip()
        release["description"] = description.strip()
    policies = {
        row["policy_name"]: json.loads(row["policy_json"])
        for row in connection.execute(
            "SELECT policy_name, policy_json FROM scientific_policies "
            "ORDER BY policy_name"
        )
    }
    pairs_by_target: dict[str, list[dict[str, Any]]] = defaultdict(list)
    pairs_by_ligand: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        pairs_by_target[str(row["target_key"])].append(row)
        pairs_by_ligand[int(row["ligand_id"])].append(row)
    target_rows = [
        dict(row)
        for row in connection.execute(
            """SELECT r.receptor_context_id, r.run_id, r.pdb_id, r.variant,
            r.ph_label, r.library, r.status, r.failure_reason, r.pocket_method,
            r.center_json, r.box_json,
            COALESCE(a.native_redock_status, r.native_redock_status)
                AS native_redock_status,
            COALESCE(a.native_redock_reason, r.native_redock_reason)
                AS native_redock_reason,
            a.receptor_classification, a.qualification_status,
            a.classification_method, a.chemistry_evidence_status,
            a.prepared_receptor_sha256, a.canonical_chemistry_policy_sha256,
            a.retained_metal_atom_count, a.retained_cofactor_residue_count,
            a.requested_observed_conflict,
            a.qualification_reason, a.native_redock_rmsd,
            p.protein_key, p.uniprot_id, p.gene_symbol,
            p.display_name
            FROM receptor_contexts r
            LEFT JOIN receptor_annotations a
              ON a.receptor_context_id=r.receptor_context_id
            LEFT JOIN protein_identities p
              ON p.protein_identity_id=a.protein_identity_id
            ORDER BY r.run_id, r.pdb_id, r.variant, r.ph_label"""
        )
    ]
    receptor_audits = [
        dict(row)
        for row in connection.execute(
            """SELECT receptor_audit_id, receptor_context_id, audit_kind,
            source_sha256, parse_status
            FROM receptor_audits
            ORDER BY receptor_context_id, audit_kind, receptor_audit_id"""
        )
    ]
    audits_by_context: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for audit in receptor_audits:
        audits_by_context[int(audit["receptor_context_id"])].append(audit)
    targets: list[dict[str, Any]] = []
    for target in target_rows:
        target["receptor_label"] = _receptor_label(target)
        target["target_key"] = _target_key(target)
        target["id"] = target["target_key"]
        target["target_id"] = target["target_key"]
        target_audits = audits_by_context[int(target["receptor_context_id"])]
        target["receptor_audit_count"] = len(target_audits)
        target["receptor_audit_kinds"] = sorted(
            {str(audit["audit_kind"]) for audit in target_audits}
        )
        target["accepted_receptor_audit_kinds"] = sorted(
            {
                str(audit["audit_kind"])
                for audit in target_audits
                if audit["parse_status"] == "parsed"
            }
        )
        target["receptor_audit_status_counts"] = _counts(target_audits, "parse_status")
        target["rejected_receptor_audit_count"] = sum(
            str(audit["parse_status"]).startswith("rejected:")
            for audit in target_audits
        )
        for name in ("center_json", "box_json"):
            value = target.pop(name, None)
            target[name.removesuffix("_json")] = json.loads(value) if value else None
        target_pairs = pairs_by_target[str(target["target_key"])]
        target["pair_count"] = len(target_pairs)
        target["primary_score_count"] = sum(
            int(row["primary_score_present"]) for row in target_pairs
        )
        target["status_counts"] = _counts(target_pairs, "final_status")
        targets.append(target)

    ligand_rows = [
        dict(row)
        for row in connection.execute(
            """SELECT ligand_id, canonical_id,
            COALESCE(display_name, canonical_id) AS display_name
            FROM ligands ORDER BY canonical_id"""
        )
    ]
    ligands: list[dict[str, Any]] = []
    for ligand in ligand_rows:
        ligand["id"] = ligand["canonical_id"]
        ligand["drug_id"] = ligand["canonical_id"]
        ligand_pairs = pairs_by_ligand[int(ligand["ligand_id"])]
        ligand["pair_count"] = len(ligand_pairs)
        ligand["primary_score_count"] = sum(
            int(row["primary_score_present"]) for row in ligand_pairs
        )
        ligand["status_counts"] = _counts(ligand_pairs, "final_status")
        ligands.append(ligand)

    artifacts = []
    for row in connection.execute(
        """SELECT a.artifact_id, a.pair_cell_id, a.receptor_context_id,
            a.ligand_id, l.canonical_id AS ligand_canonical_id,
            a.artifact_role, a.artifact_scope, a.stage, a.mode,
            a.member_name, a.sha256, a.size_bytes, a.file_type, a.verified,
            a.artifact_json
            FROM artifacts a LEFT JOIN ligands l ON l.ligand_id=a.ligand_id
            ORDER BY a.artifact_id"""
    ):
        artifact = dict(row)
        try:
            artifact_entry = json.loads(str(artifact.pop("artifact_json") or "{}"))
        except (TypeError, ValueError):
            artifact_entry = {}
        role = str(artifact["artifact_role"])
        if not artifact_public_record_allowed(
            role,
            artifact_entry,
            verified=artifact.get("verified"),
            sha256=artifact.get("sha256"),
        ):
            continue
        artifact["publication_policy"] = artifact_publication_policy(role)
        artifacts.append(artifact)
    scored = [row for row in rows if row["primary_score_present"]]
    rank_eligible = [row for row in rows if row["rank_eligible"]]
    apo_rank_eligible = [row for row in rows if row["apo_exploratory_rank_eligible"]]
    source_counts = _counts(scored, "final_score_source")
    effective_source_counts = _counts(scored, "final_score_source_effective")
    source_classification_counts = _counts(scored, "final_score_source_classification")
    return {
        "payload_schema_version": 1,
        "release": release,
        "scientific_policies": policies,
        "pose_validation_contract": pose_validation_browser_contract(policies),
        "ranking_contract": {
            "headline_track": {
                "name": "qualified_holo",
                "receptor_classification": "HOLO",
                "eligible_field": "rank_eligible",
                "within_receptor_rank_field": "rank_within_receptor",
                "across_receptors_rank_field": "rank_across_receptors",
            },
            "exploratory_track": {
                "name": "exploratory_apo",
                "receptor_classification": "APO",
                "eligible_field": "apo_exploratory_rank_eligible",
                "within_receptor_rank_field": "apo_rank_within_receptor",
                "across_receptors_rank_field": "apo_rank_across_receptors",
                "label": "APO exploratory",
            },
            "apo_holo_pooling": "forbidden",
            "score_direction": "higher_is_better",
        },
        "artifact_contract": {
            "embedding": "top_level_collection",
            "pair_join_field": "pair_cell_id",
            "public_urls_included": False,
            "absolute_paths_included": False,
            "role_field": "artifact_role",
            "controlled_roles": list(ARTIFACT_ROLES),
            "publication_policies": dict(sorted(ARTIFACT_PUBLICATION_POLICIES.items())),
            "private_roles_are_never_projected": True,
            "private_role_metadata_included": False,
            "conditional_pose_rule": (
                "docking_pose is public only when its exact artifact is selected "
                "for the frozen release"
            ),
            "sanitization_required_before_publication": True,
        },
        "score_contract": {
            "primary_field": "final_score",
            "declared_source_field": "final_score_source",
            "reconstructed_source_field": "final_score_source_reconstructed",
            "effective_source_field": "final_score_source_effective",
            "source_classification_field": "final_score_source_classification",
            "direction": "higher_is_better",
            "no_fallback": True,
            "source_classification": score_source_browser_contract(),
            "normalized_comparison_field": "atlas_score",
            "normalized_comparison_source_field": "atlas_score_source",
        },
        "coverage": {
            "pair_count": len(rows),
            "primary_score_count": len(scored),
            "primary_score_missing_count": len(rows) - len(scored),
            "primary_score_fraction": (len(scored) / len(rows)) if rows else 0.0,
            "rank_eligible_count": len(rank_eligible),
            "rank_eligible_fraction": (len(rank_eligible) / len(rows) if rows else 0.0),
            "apo_exploratory_rank_eligible_count": len(apo_rank_eligible),
            "apo_exploratory_rank_eligible_fraction": (
                len(apo_rank_eligible) / len(rows) if rows else 0.0
            ),
            "final_score_source_counts": source_counts,
            "effective_final_score_source_counts": effective_source_counts,
            "final_score_source_classification_counts": (source_classification_counts),
            "classified_final_score_source_count": sum(
                int(row["final_score_source_eligible"]) for row in scored
            ),
            "pair_status_counts": _counts(rows, "final_status"),
            "receptor_classification_counts": _counts(
                targets, "receptor_classification"
            ),
            "chemistry_evidence_status_counts": _counts(
                targets, "chemistry_evidence_status"
            ),
            "receptor_audit_count": len(receptor_audits),
            "receptor_audit_kind_counts": _counts(receptor_audits, "audit_kind"),
            "receptor_audit_status_counts": _counts(receptor_audits, "parse_status"),
        },
        "targets": targets,
        "ligands": ligands,
        "pairs": [{name: row.get(name) for name in PAIR_COLUMNS} for row in rows],
        "artifacts": artifacts,
        "receptor_audits": receptor_audits,
    }


def export_release_database(
    database_path: Path,
    output_dir: Path,
    *,
    include_parquet: bool = True,
) -> dict[str, Any]:
    """Export one SQLite snapshot for browser use and offline analysis."""
    output_dir.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(database_path)
    try:
        rows = _pair_rows(connection)
        payload = _browser_payload(connection, rows)
    finally:
        connection.close()

    browser_path = output_dir / "release_browser.json"
    csv_path = output_dir / "pairs.csv"
    browser_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    _write_csv(csv_path, rows)
    parquet_path: Path | None = None
    parquet_error: str | None = None
    if include_parquet:
        candidate = output_dir / "pairs.parquet"
        try:
            _write_parquet(candidate, rows)
            parquet_path = candidate
        except RuntimeError as exc:
            parquet_error = str(exc)

    coverage = payload["coverage"]
    summary = {
        "database_file": database_path.name,
        "browser_payload": browser_path.name,
        "pairs_csv": csv_path.name,
        "pairs_parquet": parquet_path.name if parquet_path else None,
        "parquet_error": parquet_error,
        "pair_count": len(rows),
        "primary_score_count": coverage["primary_score_count"],
        "primary_score_missing_count": coverage["primary_score_missing_count"],
        "primary_score_fraction": coverage["primary_score_fraction"],
        "rank_eligible_count": coverage["rank_eligible_count"],
        "rank_eligible_fraction": coverage["rank_eligible_fraction"],
        "apo_exploratory_rank_eligible_count": coverage[
            "apo_exploratory_rank_eligible_count"
        ],
        "apo_exploratory_rank_eligible_fraction": coverage[
            "apo_exploratory_rank_eligible_fraction"
        ],
        "final_score_source_counts": coverage["final_score_source_counts"],
        "effective_final_score_source_counts": coverage[
            "effective_final_score_source_counts"
        ],
        "final_score_source_classification_counts": coverage[
            "final_score_source_classification_counts"
        ],
        "classified_final_score_source_count": coverage[
            "classified_final_score_source_count"
        ],
    }
    summary_path = output_dir / "export_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    summary["summary_file"] = summary_path.name
    return summary
