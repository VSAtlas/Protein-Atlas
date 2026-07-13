"""Publication-safe browser and tabular exports for Atlas release databases."""

from __future__ import annotations

import csv
import json
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

from analysis.atlas_database.annotations import is_explicitly_qualified


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
    "is_control",
    "is_decoy",
    "final_score",
    "final_score_source",
    "primary_score_present",
    "rank_eligible",
    "ranking_eligibility_reason",
    "final_rank",
    "rank_within_receptor",
    "rank_across_receptors",
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


def _final_score_source(
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
    return source or "legacy_unclassified_final_score"


def receptor_ranking_eligibility(
    row: Mapping[str, Any],
) -> tuple[int, str | None]:
    """Apply only the explicit receptor-level headline ranking gates."""
    if row.get("is_control") == 1:
        return 0, "native_control"
    if row.get("is_decoy") == 1:
        return 0, "decoy"
    variant = str(row.get("variant") or "").strip().upper()
    if variant == "APO":
        return 0, "apo_receptor"
    if variant != "HOLO":
        return 0, "receptor_variant_not_holo"
    # Classification vocabulary and run/source conflict handling are scientific
    # policies awaiting user approval. Preserve the annotation for display and
    # provenance, but never interpret its free text in a headline rank.
    if str(row.get("receptor_classification") or "").strip():
        return 0, "receptor_classification_policy_pending"
    if not is_explicitly_qualified(row.get("receptor_qualification_status")):
        return 0, "receptor_quality_not_qualified"
    if not is_explicitly_qualified(row.get("native_redock_status")):
        return 0, "native_redock_not_qualified"
    return 1, None


def ranking_eligibility(row: Mapping[str, Any]) -> tuple[int, str | None]:
    """Apply receptor qualification and pair-level primary ranking gates."""
    receptor_eligible, reason = receptor_ranking_eligibility(row)
    if not receptor_eligible:
        return receptor_eligible, reason
    if row.get("final_score") is None:
        return 0, "missing_final_score"
    if row.get("pose_valid") == 0:
        return 0, "pose_invalid"
    if row.get("pose_valid") is None:
        return 0, "pose_validation_missing"
    return 1, None


def _rank(rows: list[dict[str, Any]], group_name: str, output_name: str) -> None:
    groups: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[row[group_name]].append(row)
    for values in groups.values():
        scored = [row for row in values if row["rank_eligible"]]
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
    query = """
        SELECT p.pair_cell_id, r.run_id, r.pdb_id, r.variant, r.ph_label,
            ra.receptor_classification,
            ra.qualification_status AS receptor_qualification_status,
            COALESCE(ra.native_redock_status, r.native_redock_status)
                AS native_redock_status,
            l.ligand_id, l.canonical_id AS ligand_canonical_id,
            COALESCE(l.display_name, l.canonical_id) AS ligand_display_name,
            p.final_status, p.failure_code, p.failure_reason, p.expected,
            p.has_result, p.pose_valid, p.is_control, p.is_decoy, p.final_score,
            p.final_score_source, p.final_rank, p.atlas_score, p.atlas_score_source,
            p.selected_docking_score, p.consensus_score, p.result_json,
            COUNT(a.artifact_id) AS artifact_count,
            COALESCE(SUM(CASE WHEN a.verified = 1 THEN 1 ELSE 0 END), 0)
                AS verified_artifact_count
        FROM pair_cells p
        JOIN receptor_contexts r
          ON r.receptor_context_id = p.receptor_context_id
        JOIN ligands l ON l.ligand_id = p.ligand_id
        LEFT JOIN receptor_annotations ra
          ON ra.receptor_context_id = r.receptor_context_id
        LEFT JOIN artifacts a ON a.pair_cell_id = p.pair_cell_id
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
        row["final_score_source"] = _final_score_source(
            row.pop("result_json", None),
            row["final_score"],
            row.get("final_score_source"),
        )
        row["primary_score_present"] = int(row["final_score"] is not None)
        (
            row["rank_eligible"],
            row["ranking_eligibility_reason"],
        ) = ranking_eligibility(row)
        row["score_component_coverage"] = _component_coverage(row)
        row["rank_within_receptor"] = None
        row["rank_across_receptors"] = None
        rows.append(row)
    _rank(rows, "target_key", "rank_within_receptor")
    _rank(rows, "ligand_id", "rank_across_receptors")
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
    targets: list[dict[str, Any]] = []
    for target in target_rows:
        target["receptor_label"] = _receptor_label(target)
        target["target_key"] = _target_key(target)
        target["id"] = target["target_key"]
        target["target_id"] = target["target_key"]
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

    artifacts = [
        dict(row)
        for row in connection.execute(
            """SELECT artifact_id, pair_cell_id, receptor_context_id, stage, mode,
            member_name, sha256, size_bytes, file_type, verified
            FROM artifacts ORDER BY artifact_id"""
        )
    ]
    scored = [row for row in rows if row["primary_score_present"]]
    rank_eligible = [row for row in rows if row["rank_eligible"]]
    source_counts = _counts(scored, "final_score_source")
    return {
        "payload_schema_version": 1,
        "release": release,
        "scientific_policies": policies,
        "artifact_contract": {
            "embedding": "top_level_collection",
            "pair_join_field": "pair_cell_id",
            "public_urls_included": False,
            "absolute_paths_included": False,
        },
        "score_contract": {
            "primary_field": "final_score",
            "source_field": "final_score_source",
            "direction": "higher_is_better",
            "no_fallback": True,
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
            "final_score_source_counts": source_counts,
            "classified_final_score_source_count": sum(
                count
                for source, count in source_counts.items()
                if source != "legacy_unclassified_final_score"
            ),
            "pair_status_counts": _counts(rows, "final_status"),
        },
        "targets": targets,
        "ligands": ligands,
        "pairs": [{name: row.get(name) for name in PAIR_COLUMNS} for row in rows],
        "artifacts": artifacts,
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
        "final_score_source_counts": coverage["final_score_source_counts"],
    }
    summary_path = output_dir / "export_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    summary["summary_file"] = summary_path.name
    return summary
