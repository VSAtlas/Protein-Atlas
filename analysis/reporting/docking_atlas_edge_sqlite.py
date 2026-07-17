"""Bounded-memory SQLite-to-edge projection with coarse pair-detail shards."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sqlite3
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from analysis.atlas_database.artifact_contract import (
    artifact_public_record_allowed,
    artifact_publication_policy,
)
from analysis.atlas_database.exports import (
    PAIR_COLUMNS,
    _component_coverage,
    _selected_ligand_stereo_projection,
    apo_exploratory_ranking_eligibility,
    ranking_eligibility,
)
from analysis.atlas_database.score_source_contract import (
    classify_materialized_source,
    score_source_browser_contract,
)
from analysis.reporting.docking_atlas_delivery import (
    DOWNLOAD_MANIFEST_NAME,
    build_download_projection,
)
from analysis.reporting.docking_atlas_edge import (
    _BUNDLE_MARKER,
    _counterpart_summary,
    _display_label,
    _pair_summary,
    _prepare_output,
    _release_token,
    _route_token,
    _validate_locations,
    _write_bytes,
    _write_json,
)
from analysis.reporting.docking_atlas_edge_assets import (
    APP_CSS,
    APP_JS,
    INDEX_HTML,
    PUBLIC_HEADERS,
    WORKER_SOURCE,
    WRANGLER_TEMPLATE,
)

STREAM_LAYOUT = "sqlite_streaming_coarse_v1"
DEFAULT_BATCH_ROWS = 1_000
DEFAULT_COARSE_SHARD_ROWS = 2_000
DEFAULT_MAX_PAIRS = 2_000_000
MAX_ENTITY_COUNT = 100_000
_PRIVATE_PATH_MARKERS = ("/stor/", "/home/", "/tmp/")
_SUPPORTED_SCHEMA_VERSIONS = {2, 3, 4, 5, 6, 7}
_WORK_NAME = ".atlas-edge-stream-work.sqlite"

_INTEGER_COLUMNS = {
    "pair_cell_id",
    "ligand_id",
    "expected",
    "has_result",
    "pose_valid",
    "is_control",
    "is_decoy",
    "primary_score_present",
    "final_score_source_eligible",
    "rank_eligible",
    "apo_exploratory_rank_eligible",
    "final_rank",
    "rank_within_receptor",
    "rank_across_receptors",
    "apo_rank_within_receptor",
    "apo_rank_across_receptors",
    "artifact_count",
    "verified_artifact_count",
}
_REAL_COLUMNS = {
    "final_score",
    "atlas_score",
    "selected_docking_score",
    "consensus_score",
}
_EXTRA_WORK_COLUMNS = (
    "ordinal",
    "target_route_id",
    "target_label",
    "drug_route_id",
    "drug_label",
    "pair_route_id",
    "pair_shard_id",
)


def _json_line(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")


def _safe_value(value: Any, label: str) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            _safe_value(item, f"{label}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _safe_value(item, f"{label}[{index}]")
    elif isinstance(value, str) and any(
        marker in value for marker in _PRIVATE_PATH_MARKERS
    ):
        raise ValueError(f"public edge value contains a machine-local path: {label}")


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_entry(objects_root: Path, path: Path) -> dict[str, Any]:
    return {
        "key": path.relative_to(objects_root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": _file_sha256(path),
        "content_type": "application/json; charset=utf-8",
        "cache_control": "public, max-age=31536000, immutable",
    }


def _write_small_object(
    objects_root: Path,
    relative_key: str,
    value: Any,
    entries: list[dict[str, Any]],
) -> dict[str, Any]:
    _safe_value(value, relative_key)
    path = objects_root / relative_key
    _write_bytes(path, _json_line(value) + b"\n")
    entry = _file_entry(objects_root, path)
    entries.append(entry)
    return entry


def _write_stream_object(
    objects_root: Path,
    relative_key: str,
    metadata: Mapping[str, Any],
    array_name: str,
    rows: Iterable[Mapping[str, Any]],
    entries: list[dict[str, Any]],
) -> dict[str, Any]:
    """Write one JSON array incrementally without retaining its rows."""
    _safe_value(metadata, relative_key)
    path = objects_root / relative_key
    path.parent.mkdir(parents=True, exist_ok=True)
    prefix = _json_line(dict(metadata))
    if prefix == b"{}":
        opening = b'{"' + array_name.encode("ascii") + b'":['
    else:
        opening = prefix[:-1] + b',"' + array_name.encode("ascii") + b'":['
    with path.open("wb") as handle:
        handle.write(opening)
        first = True
        for index, row in enumerate(rows):
            _safe_value(row, f"{relative_key}.{array_name}[{index}]")
            if not first:
                handle.write(b",")
            handle.write(_json_line(row))
            first = False
        handle.write(b"]}\n")
    entry = _file_entry(objects_root, path)
    entries.append(entry)
    return entry


def _readonly_connection(database_path: Path) -> sqlite3.Connection:
    uri = f"{database_path.resolve().as_uri()}?mode=ro&immutable=1"
    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA query_only=ON")
    return connection


def _table_columns(connection: sqlite3.Connection, table: str) -> set[str]:
    return {str(row[1]) for row in connection.execute(f"PRAGMA table_info({table})")}


def _sql_column(table_alias: str, columns: set[str], name: str) -> str:
    if name in columns:
        return f"{table_alias}.{name}"
    return "NULL"


def _validate_source_database(
    connection: sqlite3.Connection, database_path: Path
) -> tuple[int, int]:
    version = int(connection.execute("PRAGMA user_version").fetchone()[0])
    if version not in _SUPPORTED_SCHEMA_VERSIONS:
        raise ValueError(
            f"unsupported Atlas database schema version {version}; expected one of "
            f"{sorted(_SUPPORTED_SCHEMA_VERSIONS)}"
        )
    required = {
        "releases",
        "scientific_policies",
        "runs",
        "receptor_contexts",
        "ligands",
        "pair_cells",
        "artifacts",
        "protein_identities",
        "receptor_annotations",
    }
    present = {
        str(row[0])
        for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )
    }
    missing = sorted(required - present)
    if missing:
        raise ValueError(f"Atlas database is missing required tables: {missing}")
    release_count = int(
        connection.execute("SELECT COUNT(*) FROM releases").fetchone()[0]
    )
    if release_count != 1:
        raise ValueError(
            f"streaming edge export requires exactly one release row; found {release_count}"
        )
    pair_count = int(
        connection.execute("SELECT COUNT(*) FROM pair_cells").fetchone()[0]
    )
    for suffix in ("-wal", "-shm"):
        sidecar = database_path.with_name(database_path.name + suffix)
        if sidecar.exists():
            raise ValueError(
                f"streaming edge export refuses a live SQLite sidecar: {sidecar.name}"
            )
    return version, pair_count


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


def _declared_final_score_source(
    result_json: Any, final_score: Any, stored_source: Any
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


def _release_payload(
    connection: sqlite3.Connection,
) -> tuple[dict[str, Any], dict[str, Any]]:
    row = connection.execute("SELECT * FROM releases").fetchone()
    if row is None:  # guarded by validation
        raise ValueError("database contains no release row")
    release = dict(row)
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
    policies: dict[str, Any] = {}
    for policy in connection.execute(
        "SELECT policy_name, policy_json FROM scientific_policies ORDER BY policy_name"
    ):
        try:
            policies[str(policy["policy_name"])] = json.loads(policy["policy_json"])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"invalid scientific policy JSON: {policy['policy_name']}"
            ) from exc
    _safe_value(release, "release")
    _safe_value(policies, "scientific_policies")
    return release, policies


def _entities(
    connection: sqlite3.Connection,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, str],
    dict[str, str],
    dict[str, str],
    dict[str, str],
]:
    targets: list[dict[str, Any]] = []
    target_routes: dict[str, str] = {}
    target_labels: dict[str, str] = {}
    context_columns = _table_columns(connection, "receptor_contexts")
    annotation_columns = _table_columns(connection, "receptor_annotations")
    query = f"""
        SELECT r.receptor_context_id, r.run_id, r.pdb_id, r.variant,
        r.ph_label, r.library, r.status, r.failure_reason, r.pocket_method,
        r.center_json, r.box_json,
        COALESCE(
            {_sql_column("a", annotation_columns, "native_redock_status")},
            {_sql_column("r", context_columns, "native_redock_status")}
        ) AS native_redock_status,
        COALESCE(
            {_sql_column("a", annotation_columns, "native_redock_reason")},
            {_sql_column("r", context_columns, "native_redock_reason")}
        ) AS native_redock_reason,
        {_sql_column("a", annotation_columns, "receptor_classification")}
            AS receptor_classification,
        {_sql_column("a", annotation_columns, "classification_method")}
            AS receptor_classification_method,
        {_sql_column("a", annotation_columns, "chemistry_evidence_status")}
            AS chemistry_evidence_status,
        {_sql_column("a", annotation_columns, "prepared_receptor_sha256")}
            AS prepared_receptor_sha256,
        {_sql_column("a", annotation_columns, "canonical_chemistry_policy_sha256")}
            AS canonical_chemistry_policy_sha256,
        {_sql_column("a", annotation_columns, "retained_metal_atom_count")}
            AS retained_metal_atom_count,
        {_sql_column("a", annotation_columns, "retained_cofactor_residue_count")}
            AS retained_cofactor_residue_count,
        {_sql_column("a", annotation_columns, "requested_observed_conflict")}
            AS requested_observed_conflict,
        {_sql_column("a", annotation_columns, "qualification_status")}
            AS qualification_status,
        {_sql_column("a", annotation_columns, "qualification_reason")}
            AS qualification_reason,
        {_sql_column("a", annotation_columns, "native_redock_rmsd")}
            AS native_redock_rmsd,
        p.protein_key, p.uniprot_id, p.gene_symbol, p.display_name
        FROM receptor_contexts r
        LEFT JOIN receptor_annotations a
          ON a.receptor_context_id=r.receptor_context_id
        LEFT JOIN protein_identities p
          ON p.protein_identity_id=a.protein_identity_id
        ORDER BY r.run_id, r.pdb_id, r.variant, r.ph_label
    """
    for raw in connection.execute(query):
        target = dict(raw)
        target["receptor_label"] = _receptor_label(target)
        target_id = _target_key(target)
        target["target_key"] = target_id
        target["id"] = target_id
        target["target_id"] = target_id
        for name in ("center_json", "box_json"):
            value = target.pop(name, None)
            try:
                target[name.removesuffix("_json")] = (
                    json.loads(value) if value else None
                )
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid receptor {name}: {target_id}") from exc
        route_id = _route_token(target_id, "target identifier")
        target["route_id"] = route_id
        target["pair_count"] = 0
        target["primary_score_count"] = 0
        target["status_counts"] = {}
        _safe_value(target, f"target {target_id}")
        targets.append(target)
        target_routes[target_id] = route_id
        target_labels[target_id] = _display_label(target, target_id)

    stereo_by_ligand = _selected_ligand_stereo_projection(connection)
    drugs: list[dict[str, Any]] = []
    drug_routes: dict[str, str] = {}
    drug_labels: dict[str, str] = {}
    for raw in connection.execute(
        """SELECT ligand_id, canonical_id,
        COALESCE(display_name, canonical_id) AS display_name
        FROM ligands ORDER BY canonical_id"""
    ):
        drug = dict(raw)
        drug_id = str(drug["canonical_id"])
        drug["id"] = drug_id
        drug["drug_id"] = drug_id
        route_id = _route_token(drug_id, "drug identifier")
        drug["route_id"] = route_id
        drug["pair_count"] = 0
        drug["primary_score_count"] = 0
        drug["status_counts"] = {}
        drug.update(
            stereo_by_ligand.get(
                int(drug["ligand_id"]),
                {
                    "selected_stereo_evidence_count": 0,
                    "selected_stereo_evidence": [],
                },
            )
        )
        _safe_value(drug, f"drug {drug_id}")
        drugs.append(drug)
        drug_routes[drug_id] = route_id
        drug_labels[drug_id] = _display_label(drug, drug_id)
    if len(targets) > MAX_ENTITY_COUNT or len(drugs) > MAX_ENTITY_COUNT:
        raise ValueError(
            "streaming edge export refuses more than "
            f"{MAX_ENTITY_COUNT} targets or drugs"
        )
    return (
        targets,
        drugs,
        target_routes,
        target_labels,
        drug_routes,
        drug_labels,
    )


def _work_schema(connection: sqlite3.Connection) -> None:
    definitions: list[str] = []
    for name in PAIR_COLUMNS:
        column_type = (
            "INTEGER"
            if name in _INTEGER_COLUMNS
            else "REAL"
            if name in _REAL_COLUMNS
            else "TEXT"
        )
        definitions.append(f"{name} {column_type}")
    definitions.extend(
        (
            "ordinal INTEGER NOT NULL UNIQUE",
            "target_route_id TEXT NOT NULL",
            "target_label TEXT NOT NULL",
            "drug_route_id TEXT NOT NULL",
            "drug_label TEXT NOT NULL",
            "pair_route_id TEXT NOT NULL UNIQUE",
            "pair_shard_id INTEGER NOT NULL",
        )
    )
    connection.execute(
        "CREATE TABLE stream_pairs ("
        + ",".join(definitions)
        + ", PRIMARY KEY(pair_cell_id))"
    )


def _source_pair_query(connection: sqlite3.Connection) -> str:
    pair_columns = _table_columns(connection, "pair_cells")
    annotation_columns = _table_columns(connection, "receptor_annotations")
    attempt_columns = _table_columns(connection, "result_attempts")
    if attempt_columns:
        attempt_projection = """
            sra.input_csv_sha256 AS selected_result_input_csv_sha256,
            sra.source_row_number AS selected_result_source_row_number,
            sra.result_sha256 AS selected_result_sha256,"""
        attempt_join = """
        LEFT JOIN result_attempts sra
          ON sra.pair_cell_id=p.pair_cell_id
         AND sra.selected_for_release=1"""
    else:
        attempt_projection = """
            NULL AS selected_result_input_csv_sha256,
            NULL AS selected_result_source_row_number,
            NULL AS selected_result_sha256,"""
        attempt_join = ""
    return f"""
        SELECT p.pair_cell_id, r.run_id, r.pdb_id, r.variant, r.ph_label,
            {_sql_column("a", annotation_columns, "receptor_classification")}
                AS receptor_classification,
            {_sql_column("a", annotation_columns, "classification_method")}
                AS receptor_classification_method,
            {_sql_column("a", annotation_columns, "chemistry_evidence_status")}
                AS chemistry_evidence_status,
            {_sql_column("a", annotation_columns, "prepared_receptor_sha256")}
                AS prepared_receptor_sha256,
            {_sql_column("a", annotation_columns, "canonical_chemistry_policy_sha256")}
                AS canonical_chemistry_policy_sha256,
            {_sql_column("a", annotation_columns, "retained_metal_atom_count")}
                AS retained_metal_atom_count,
            {_sql_column("a", annotation_columns, "retained_cofactor_residue_count")}
                AS retained_cofactor_residue_count,
            {_sql_column("a", annotation_columns, "requested_observed_conflict")}
                AS requested_observed_conflict,
            {_sql_column("a", annotation_columns, "qualification_status")}
                AS receptor_qualification_status,
            COALESCE(
                {_sql_column("a", annotation_columns, "native_redock_status")},
                r.native_redock_status
            ) AS native_redock_status,
            l.ligand_id, l.canonical_id AS ligand_canonical_id,
            COALESCE(l.display_name, l.canonical_id) AS ligand_display_name,
            p.final_status, p.failure_code, p.failure_reason, p.expected,
            p.has_result, p.pose_valid,
            {_sql_column("p", pair_columns, "pose_validation_method")}
                AS pose_validation_method,
            {_sql_column("p", pair_columns, "pose_validation_scope")}
                AS pose_validation_scope,
            {_sql_column("p", pair_columns, "pose_validation_thresholds_json")}
                AS pose_validation_thresholds_json,
            p.is_control, p.is_decoy, p.final_score,
            {_sql_column("p", pair_columns, "final_score_source")}
                AS final_score_source,
            {_sql_column("p", pair_columns, "final_score_source_reconstructed")}
                AS final_score_source_reconstructed,
            {_sql_column("p", pair_columns, "final_score_source_classification")}
                AS final_score_source_classification,
            {_sql_column("p", pair_columns, "final_score_source_evidence_json")}
                AS final_score_source_evidence_json,
            {attempt_projection}
            p.final_rank, p.atlas_score, p.atlas_score_source,
            p.selected_docking_score, p.consensus_score, p.result_json,
            0 AS artifact_count,
            0 AS verified_artifact_count
        FROM pair_cells p
        JOIN receptor_contexts r
          ON r.receptor_context_id = p.receptor_context_id
        JOIN ligands l ON l.ligand_id = p.ligand_id
        LEFT JOIN receptor_annotations a
          ON a.receptor_context_id = r.receptor_context_id
        {attempt_join}
        ORDER BY p.pair_cell_id
    """


def _normalize_pair(
    raw: Mapping[str, Any],
    *,
    ordinal: int,
    coarse_shard_rows: int,
    target_routes: Mapping[str, str],
    target_labels: Mapping[str, str],
    drug_routes: Mapping[str, str],
    drug_labels: Mapping[str, str],
) -> dict[str, Any]:
    row = dict(raw)
    row["receptor_label"] = _receptor_label(row)
    target_id = _target_key(row)
    drug_id = str(row["ligand_canonical_id"])
    if target_id not in target_routes:
        raise ValueError(f"pair references unknown target: {target_id}")
    if drug_id not in drug_routes:
        raise ValueError(f"pair references unknown drug: {drug_id}")
    row["target_key"] = target_id
    row["target_id"] = target_id
    row["drug_id"] = drug_id
    row["final_score_source"] = _declared_final_score_source(
        row.pop("result_json", None),
        row.get("final_score"),
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
    row["primary_score_present"] = int(row.get("final_score") is not None)
    row["rank_eligible"], row["ranking_eligibility_reason"] = ranking_eligibility(row)
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
    row["rank_within_receptor"] = None
    row["rank_across_receptors"] = None
    row["apo_rank_within_receptor"] = None
    row["apo_rank_across_receptors"] = None
    pair_id = str(row["pair_cell_id"])
    normalized = {name: row.get(name) for name in PAIR_COLUMNS}
    normalized.update(
        {
            "ordinal": ordinal,
            "target_route_id": target_routes[target_id],
            "target_label": target_labels[target_id],
            "drug_route_id": drug_routes[drug_id],
            "drug_label": drug_labels[drug_id],
            "pair_route_id": _route_token(pair_id, "pair cell identifier"),
            "pair_shard_id": (ordinal - 1) // coarse_shard_rows,
        }
    )
    _safe_value(normalized, f"pair {pair_id}")
    return normalized


def _populate_work_database(
    source: sqlite3.Connection,
    work: sqlite3.Connection,
    *,
    pair_count: int,
    batch_rows: int,
    coarse_shard_rows: int,
    target_routes: Mapping[str, str],
    target_labels: Mapping[str, str],
    drug_routes: Mapping[str, str],
    drug_labels: Mapping[str, str],
) -> dict[str, Any]:
    columns = (*PAIR_COLUMNS, *_EXTRA_WORK_COLUMNS)
    placeholders = ",".join("?" for _ in columns)
    insert_sql = (
        f"INSERT INTO stream_pairs ({','.join(columns)}) VALUES ({placeholders})"
    )
    status_counts: dict[str, int] = defaultdict(int)
    score_source_counts: dict[str, int] = defaultdict(int)
    effective_score_source_counts: dict[str, int] = defaultdict(int)
    score_source_classification_counts: dict[str, int] = defaultdict(int)
    primary_score_count = 0
    classified_score_source_count = 0
    rank_eligible_count = 0
    apo_rank_eligible_count = 0
    source_cursor = source.execute(_source_pair_query(source))
    ordinal = 0
    while True:
        batch = source_cursor.fetchmany(batch_rows)
        if not batch:
            break
        values: list[tuple[Any, ...]] = []
        for raw in batch:
            ordinal += 1
            pair = _normalize_pair(
                raw,
                ordinal=ordinal,
                coarse_shard_rows=coarse_shard_rows,
                target_routes=target_routes,
                target_labels=target_labels,
                drug_routes=drug_routes,
                drug_labels=drug_labels,
            )
            values.append(tuple(pair.get(name) for name in columns))
            status_counts[str(pair.get("final_status") or "unknown")] += 1
            if pair["primary_score_present"]:
                primary_score_count += 1
                score_source_counts[
                    str(pair.get("final_score_source") or "unknown")
                ] += 1
                effective_score_source_counts[
                    str(pair.get("final_score_source_effective") or "unknown")
                ] += 1
                score_source_classification_counts[
                    str(pair.get("final_score_source_classification") or "unknown")
                ] += 1
                classified_score_source_count += int(
                    pair.get("final_score_source_eligible") or 0
                )
            rank_eligible_count += int(pair["rank_eligible"] or 0)
            apo_rank_eligible_count += int(pair["apo_exploratory_rank_eligible"] or 0)
        work.executemany(insert_sql, values)
        work.commit()
    if ordinal != pair_count:
        raise ValueError(
            f"pair count changed during streaming export: expected {pair_count}, read {ordinal}"
        )
    work.executescript(
        """
        CREATE INDEX idx_stream_target
          ON stream_pairs(target_id, rank_eligible, final_score);
        CREATE INDEX idx_stream_drug
          ON stream_pairs(drug_id, rank_eligible, final_score);
        CREATE INDEX idx_stream_apo_target
          ON stream_pairs(
            target_id, apo_exploratory_rank_eligible, final_score
          );
        CREATE INDEX idx_stream_apo_drug
          ON stream_pairs(
            drug_id, apo_exploratory_rank_eligible, final_score
          );
        CREATE INDEX idx_stream_shard ON stream_pairs(pair_shard_id, ordinal);
        CREATE TABLE target_ranks AS
          SELECT pair_cell_id,
          RANK() OVER (PARTITION BY target_id ORDER BY final_score DESC) AS rank_value
          FROM stream_pairs WHERE rank_eligible=1;
        CREATE UNIQUE INDEX idx_target_rank_pair ON target_ranks(pair_cell_id);
        CREATE TABLE drug_ranks AS
          SELECT pair_cell_id,
          RANK() OVER (PARTITION BY ligand_id ORDER BY final_score DESC) AS rank_value
          FROM stream_pairs WHERE rank_eligible=1;
        CREATE UNIQUE INDEX idx_drug_rank_pair ON drug_ranks(pair_cell_id);
        CREATE TABLE apo_target_ranks AS
          SELECT pair_cell_id,
          RANK() OVER (PARTITION BY target_id ORDER BY final_score DESC) AS rank_value
          FROM stream_pairs WHERE apo_exploratory_rank_eligible=1;
        CREATE UNIQUE INDEX idx_apo_target_rank_pair
          ON apo_target_ranks(pair_cell_id);
        CREATE TABLE apo_drug_ranks AS
          SELECT pair_cell_id,
          RANK() OVER (PARTITION BY ligand_id ORDER BY final_score DESC) AS rank_value
          FROM stream_pairs WHERE apo_exploratory_rank_eligible=1;
        CREATE UNIQUE INDEX idx_apo_drug_rank_pair
          ON apo_drug_ranks(pair_cell_id);
        UPDATE stream_pairs SET rank_within_receptor=(
          SELECT rank_value FROM target_ranks
          WHERE target_ranks.pair_cell_id=stream_pairs.pair_cell_id
        );
        UPDATE stream_pairs SET rank_across_receptors=(
          SELECT rank_value FROM drug_ranks
          WHERE drug_ranks.pair_cell_id=stream_pairs.pair_cell_id
        );
        UPDATE stream_pairs SET apo_rank_within_receptor=(
          SELECT rank_value FROM apo_target_ranks
          WHERE apo_target_ranks.pair_cell_id=stream_pairs.pair_cell_id
        );
        UPDATE stream_pairs SET apo_rank_across_receptors=(
          SELECT rank_value FROM apo_drug_ranks
          WHERE apo_drug_ranks.pair_cell_id=stream_pairs.pair_cell_id
        );
        DROP TABLE target_ranks;
        DROP TABLE drug_ranks;
        DROP TABLE apo_target_ranks;
        DROP TABLE apo_drug_ranks;
        """
    )
    work.commit()
    return {
        "pair_count": pair_count,
        "primary_score_count": primary_score_count,
        "primary_score_missing_count": pair_count - primary_score_count,
        "primary_score_fraction": primary_score_count / pair_count
        if pair_count
        else 0.0,
        "rank_eligible_count": rank_eligible_count,
        "rank_eligible_fraction": rank_eligible_count / pair_count
        if pair_count
        else 0.0,
        "apo_exploratory_rank_eligible_count": apo_rank_eligible_count,
        "apo_exploratory_rank_eligible_fraction": (
            apo_rank_eligible_count / pair_count if pair_count else 0.0
        ),
        "classified_final_score_source_count": classified_score_source_count,
        "effective_final_score_source_counts": dict(
            sorted(effective_score_source_counts.items())
        ),
        "final_score_source_classification_counts": dict(
            sorted(score_source_classification_counts.items())
        ),
        "final_score_source_counts": dict(sorted(score_source_counts.items())),
        "pair_status_counts": dict(sorted(status_counts.items())),
    }


def _apply_entity_counts(
    work: sqlite3.Connection,
    entities: list[dict[str, Any]],
    *,
    id_name: str,
    work_id_name: str,
) -> None:
    by_id = {str(entity[id_name]): entity for entity in entities}
    query = f"""
        SELECT {work_id_name} AS entity_id, final_status, COUNT(*) AS pair_count,
        SUM(primary_score_present) AS primary_score_count
        FROM stream_pairs GROUP BY {work_id_name}, final_status
    """
    for raw in work.execute(query):
        entity = by_id[str(raw["entity_id"])]
        count = int(raw["pair_count"] or 0)
        entity["pair_count"] += count
        entity["primary_score_count"] += int(raw["primary_score_count"] or 0)
        entity["status_counts"][str(raw["final_status"] or "unknown")] = count


def _work_rows(
    work: sqlite3.Connection,
    where: str,
    parameters: Sequence[Any],
    order: str,
    *,
    batch_rows: int,
) -> Iterator[dict[str, Any]]:
    query = f"SELECT * FROM stream_pairs WHERE {where} ORDER BY {order}"
    cursor = work.execute(query, tuple(parameters))
    while True:
        rows = cursor.fetchmany(batch_rows)
        if not rows:
            return
        for row in rows:
            yield dict(row)


def _summary(row: Mapping[str, Any], entity_kind: str) -> dict[str, Any]:
    summary = _pair_summary(
        row,
        pair_id=str(row["pair_cell_id"]),
        pair_route_id=str(row["pair_route_id"]),
        target_id=str(row["target_id"]),
        target_route_id=str(row["target_route_id"]),
        target_label=str(row["target_label"]),
        drug_id=str(row["drug_id"]),
        drug_route_id=str(row["drug_route_id"]),
        drug_label=str(row["drug_label"]),
    )
    summary["pair_shard_id"] = f"{int(row['pair_shard_id']):06d}"
    return _counterpart_summary(summary, entity_kind)


def _artifact_select(connection: sqlite3.Connection) -> list[str]:
    columns = _table_columns(connection, "artifacts")
    allowed = (
        "artifact_id",
        "pair_cell_id",
        "receptor_context_id",
        "ligand_id",
        "artifact_role",
        "artifact_scope",
        "stage",
        "mode",
        "member_name",
        "sha256",
        "size_bytes",
        "file_type",
        "verified",
        "artifact_json",
    )
    return [name for name in allowed if name in columns]


def _artifacts_for_pairs(
    source: sqlite3.Connection, pair_ids: Sequence[int]
) -> dict[int, list[dict[str, Any]]]:
    by_pair: dict[int, list[dict[str, Any]]] = defaultdict(list)
    if not pair_ids:
        return by_pair
    selected = _artifact_select(source)
    for offset in range(0, len(pair_ids), 400):
        chunk = pair_ids[offset : offset + 400]
        placeholders = ",".join("?" for _ in chunk)
        query = (
            f"SELECT {','.join(selected)} FROM artifacts "
            f"WHERE pair_cell_id IN ({placeholders}) ORDER BY artifact_id"
        )
        for raw in source.execute(query, tuple(chunk)):
            artifact = dict(raw)
            pair_id = int(artifact["pair_cell_id"])
            try:
                artifact_entry = json.loads(
                    str(artifact.pop("artifact_json", None) or "{}")
                )
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
            _safe_value(artifact, f"artifact {artifact.get('artifact_id')}")
            by_pair[pair_id].append(artifact)
    return by_pair


def _pair_record(
    row: Mapping[str, Any],
    *,
    release_id: str,
    release_token: str,
    artifacts: Mapping[int, list[dict[str, Any]]],
) -> dict[str, Any]:
    pair_id = int(row["pair_cell_id"])
    public_artifacts = artifacts.get(pair_id, [])
    pair = {name: row.get(name) for name in PAIR_COLUMNS}
    pair["artifact_count"] = len(public_artifacts)
    pair["verified_artifact_count"] = len(public_artifacts)
    return {
        "schema_version": 2,
        "release_id": release_id,
        "release_token": release_token,
        "id": str(pair_id),
        "route_id": row["pair_route_id"],
        "pair_shard_id": f"{int(row['pair_shard_id']):06d}",
        "target": {
            "id": row["target_id"],
            "route_id": row["target_route_id"],
            "label": row["target_label"],
        },
        "drug": {
            "id": row["drug_id"],
            "route_id": row["drug_route_id"],
            "label": row["drug_label"],
        },
        "pair": pair,
        "artifacts": public_artifacts,
    }


def _write_entity_objects(
    work: sqlite3.Connection,
    *,
    entities: Sequence[Mapping[str, Any]],
    kind: str,
    id_field: str,
    route_field: str,
    release_id: str,
    release_token: str,
    objects_root: Path,
    prefix: str,
    entries: list[dict[str, Any]],
    batch_rows: int,
) -> list[dict[str, Any]]:
    index: list[dict[str, Any]] = []
    entity_kind = "targets" if kind == "targets" else "drugs"
    work_id = "target_id" if kind == "targets" else "drug_id"
    rank_name = "rank_within_receptor" if kind == "targets" else "rank_across_receptors"
    apo_rank_name = (
        "apo_rank_within_receptor" if kind == "targets" else "apo_rank_across_receptors"
    )
    label_name = "drug_label" if kind == "targets" else "target_label"
    for entity in entities:
        entity_id = str(entity[id_field])
        route_id = str(entity[route_field])
        rows = _work_rows(
            work,
            f"{work_id}=?",
            (entity_id,),
            (
                f"ranking_track IS NULL, "
                f"COALESCE({rank_name}, {apo_rank_name}), "
                f"final_score IS NULL, final_score DESC, "
                f"lower({label_name}), pair_cell_id"
            ),
            batch_rows=batch_rows,
        )
        summaries = (_summary(row, entity_kind) for row in rows)
        _write_stream_object(
            objects_root,
            f"{prefix}/records/{kind}/{route_id}.json",
            {
                "schema_version": 2,
                "release_id": release_id,
                "release_token": release_token,
                "id": entity_id,
                "route_id": route_id,
                "entity": dict(entity),
            },
            "pairs",
            summaries,
            entries,
        )
        index.append(
            {
                "id": entity_id,
                "route_id": route_id,
                "label": _display_label(entity, entity_id),
                "pair_count": int(entity.get("pair_count") or 0),
                "primary_score_count": int(entity.get("primary_score_count") or 0),
            }
        )
    return index


def _write_pair_shards(
    source: sqlite3.Connection,
    work: sqlite3.Connection,
    *,
    pair_count: int,
    coarse_shard_rows: int,
    release_id: str,
    release_token: str,
    objects_root: Path,
    prefix: str,
    entries: list[dict[str, Any]],
    batch_rows: int,
) -> list[dict[str, Any]]:
    descriptors: list[dict[str, Any]] = []
    shard_count = math.ceil(pair_count / coarse_shard_rows) if pair_count else 0
    for shard_number in range(shard_count):
        pair_ids = [
            int(row[0])
            for row in work.execute(
                "SELECT pair_cell_id FROM stream_pairs WHERE pair_shard_id=? "
                "ORDER BY ordinal",
                (shard_number,),
            )
        ]
        artifacts = _artifacts_for_pairs(source, pair_ids)
        rows = _work_rows(
            work,
            "pair_shard_id=?",
            (shard_number,),
            "ordinal",
            batch_rows=batch_rows,
        )
        records = (
            _pair_record(
                row,
                release_id=release_id,
                release_token=release_token,
                artifacts=artifacts,
            )
            for row in rows
        )
        shard_id = f"{shard_number:06d}"
        entry = _write_stream_object(
            objects_root,
            f"{prefix}/records/pair-shards/{shard_id}.json",
            {
                "schema_version": 2,
                "release_id": release_id,
                "release_token": release_token,
                "kind": "pair_records",
                "shard_id": shard_id,
                "first_ordinal": shard_number * coarse_shard_rows + 1,
                "count": len(pair_ids),
            },
            "records",
            records,
            entries,
        )
        descriptors.append(
            {
                "id": shard_id,
                "count": len(pair_ids),
                "first_ordinal": shard_number * coarse_shard_rows + 1,
                "sha256": entry["sha256"],
                "size_bytes": entry["size_bytes"],
                "api_path": f"/api/releases/{release_token}/pair-shards/{shard_id}",
            }
        )
    return descriptors


def _source_unchanged(database_path: Path, initial: tuple[int, int, int]) -> None:
    current = database_path.stat()
    observed = (current.st_size, current.st_mtime_ns, current.st_ino)
    if observed != initial:
        raise ValueError("source SQLite file changed during streaming edge export")


def build_streaming_edge_bundle(
    database_path: Path,
    output_dir: Path,
    *,
    source_site_dir: Path | None = None,
    overwrite: bool = False,
    download_base_url: str | None = None,
    batch_rows: int = DEFAULT_BATCH_ROWS,
    coarse_shard_rows: int = DEFAULT_COARSE_SHARD_ROWS,
    max_pairs: int = DEFAULT_MAX_PAIRS,
) -> dict[str, Any]:
    """Stream a public Atlas SQLite snapshot into coarse edge objects."""
    database_path = Path(database_path)
    output_dir = Path(output_dir)
    if database_path.is_symlink() or not database_path.is_file():
        raise ValueError(
            f"source SQLite database is missing or symlinked: {database_path}"
        )
    if not 100 <= batch_rows <= 10_000:
        raise ValueError("batch_rows must be between 100 and 10,000")
    if not 100 <= coarse_shard_rows <= 10_000:
        raise ValueError("coarse_shard_rows must be between 100 and 10,000")
    if max_pairs < 1:
        raise ValueError("max_pairs must be positive")
    resolved_database = database_path.resolve()
    resolved_output = output_dir.resolve()
    if resolved_database == resolved_output or resolved_database.is_relative_to(
        resolved_output
    ):
        raise ValueError("edge output must not contain the source SQLite database")
    if source_site_dir is not None:
        _validate_locations(Path(source_site_dir), output_dir)
    initial_stat = database_path.stat()
    stat_identity = (
        initial_stat.st_size,
        initial_stat.st_mtime_ns,
        initial_stat.st_ino,
    )
    source_sha256 = _file_sha256(database_path)
    source = _readonly_connection(database_path)
    try:
        schema_version, pair_count = _validate_source_database(source, database_path)
        if pair_count > max_pairs:
            raise ValueError(
                f"streaming edge export refuses {pair_count} pairs above --max-pairs "
                f"{max_pairs}"
            )
        release, policies = _release_payload(source)
        release_id = str(release["release_id"])
        release_token = _release_token(release_id, source_sha256)
        (
            targets,
            drugs,
            target_routes,
            target_labels,
            drug_routes,
            drug_labels,
        ) = _entities(source)
        inferred_site = source_site_dir
        if inferred_site is None:
            candidate = database_path.parent.parent
            if (candidate / "site_manifest.json").is_file():
                inferred_site = candidate
        if inferred_site is not None and source_site_dir is None:
            _validate_locations(Path(inferred_site), output_dir)
        download_projection = (
            build_download_projection(
                Path(inferred_site),
                release_id=release_id,
                release_token=release_token,
                download_base_url=download_base_url,
            )
            if download_base_url and inferred_site is not None
            else None
        )
        if download_base_url and inferred_site is None:
            raise ValueError(
                "--download-base-url requires --source-site-dir or a database under "
                "site/downloads/"
            )

        _prepare_output(output_dir, overwrite)
        _write_json(
            output_dir / _BUNDLE_MARKER,
            {
                "schema_version": 1,
                "generator": "atlas publish edge-bundle",
                "state": "building",
            },
        )
        work_path = output_dir / _WORK_NAME
        work = sqlite3.connect(work_path)
        work.row_factory = sqlite3.Row
        try:
            work.executescript(
                "PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; "
                "PRAGMA temp_store=FILE; PRAGMA cache_size=-16384;"
            )
            _work_schema(work)
            coverage = _populate_work_database(
                source,
                work,
                pair_count=pair_count,
                batch_rows=batch_rows,
                coarse_shard_rows=coarse_shard_rows,
                target_routes=target_routes,
                target_labels=target_labels,
                drug_routes=drug_routes,
                drug_labels=drug_labels,
            )
            _apply_entity_counts(
                work, targets, id_name="target_id", work_id_name="target_id"
            )
            _apply_entity_counts(work, drugs, id_name="drug_id", work_id_name="drug_id")
            work_bytes = work_path.stat().st_size
            objects_root = output_dir / "objects"
            public_root = output_dir / "public"
            prefix = f"releases/{release_token}"
            object_entries: list[dict[str, Any]] = []

            target_index = _write_entity_objects(
                work,
                entities=targets,
                kind="targets",
                id_field="target_id",
                route_field="route_id",
                release_id=release_id,
                release_token=release_token,
                objects_root=objects_root,
                prefix=prefix,
                entries=object_entries,
                batch_rows=batch_rows,
            )
            drug_index = _write_entity_objects(
                work,
                entities=drugs,
                kind="drugs",
                id_field="drug_id",
                route_field="route_id",
                release_id=release_id,
                release_token=release_token,
                objects_root=objects_root,
                prefix=prefix,
                entries=object_entries,
                batch_rows=batch_rows,
            )
            _write_small_object(
                objects_root,
                f"{prefix}/indexes/targets.json",
                {
                    "schema_version": 2,
                    "release_token": release_token,
                    "kind": "targets",
                    "records": target_index,
                },
                object_entries,
            )
            _write_small_object(
                objects_root,
                f"{prefix}/indexes/drugs.json",
                {
                    "schema_version": 2,
                    "release_token": release_token,
                    "kind": "drugs",
                    "records": drug_index,
                },
                object_entries,
            )
            pair_shards = _write_pair_shards(
                source,
                work,
                pair_count=pair_count,
                coarse_shard_rows=coarse_shard_rows,
                release_id=release_id,
                release_token=release_token,
                objects_root=objects_root,
                prefix=prefix,
                entries=object_entries,
                batch_rows=batch_rows,
            )
            _write_small_object(
                objects_root,
                f"{prefix}/indexes/pairs.json",
                {
                    "schema_version": 2,
                    "release_token": release_token,
                    "kind": "coarse_pair_shards",
                    "count": pair_count,
                    "coarse_shard_rows": coarse_shard_rows,
                    "shards": pair_shards,
                },
                object_entries,
            )
            release_manifest = {
                "schema_version": 2,
                "release": release,
                "release_token": release_token,
                "source_database_sha256": source_sha256,
                "storage_layout": STREAM_LAYOUT,
                "counts": {
                    "targets": len(targets),
                    "drugs": len(drugs),
                    "pairs": pair_count,
                    "pair_detail_shards": len(pair_shards),
                },
                "coverage": coverage,
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
                "score_contract": {
                    "primary_field": "final_score",
                    "declared_source_field": "final_score_source",
                    "reconstructed_source_field": ("final_score_source_reconstructed"),
                    "effective_source_field": "final_score_source_effective",
                    "source_classification_field": (
                        "final_score_source_classification"
                    ),
                    "direction": "higher_is_better",
                    "no_fallback": True,
                    "source_classification": score_source_browser_contract(),
                    "normalized_comparison_field": "atlas_score",
                    "normalized_comparison_source_field": "atlas_score_source",
                },
                "scientific_policies": policies,
                "downloads": (
                    download_projection["public_entries"] if download_projection else []
                ),
                "route_templates": {
                    "release": f"/releases/{release_token}",
                    "target": f"/releases/{release_token}/targets/{{target_route_id}}",
                    "drug": f"/releases/{release_token}/drugs/{{drug_route_id}}",
                    "pair": f"/releases/{release_token}/pairs/{{pair_route_id}}?shard={{pair_shard_id}}",
                },
                "api_templates": {
                    "target": f"/api/releases/{release_token}/targets/{{target_route_id}}",
                    "drug": f"/api/releases/{release_token}/drugs/{{drug_route_id}}",
                    "pair_shard": f"/api/releases/{release_token}/pair-shards/{{pair_shard_id}}",
                },
            }
            _write_small_object(
                objects_root,
                f"{prefix}/manifest.json",
                release_manifest,
                object_entries,
            )

            _write_bytes(public_root / "index.html", INDEX_HTML.encode("utf-8"))
            _write_bytes(public_root / "_headers", PUBLIC_HEADERS.encode("utf-8"))
            _write_bytes(public_root / "assets" / "app.css", APP_CSS.encode("utf-8"))
            _write_bytes(public_root / "assets" / "app.js", APP_JS.encode("utf-8"))
            _write_json(
                public_root / "release-config.json",
                {
                    "release_id": release_id,
                    "release_token": release_token,
                    "storage_layout": STREAM_LAYOUT,
                    "pair_storage": "coarse_shards",
                },
            )
            _write_bytes(
                output_dir / "src" / "index.mjs", WORKER_SOURCE.encode("utf-8")
            )
            _write_bytes(
                output_dir / "wrangler.toml", WRANGLER_TEMPLATE.encode("utf-8")
            )

            object_entries.sort(key=lambda row: str(row["key"]))
            object_manifest = {
                "schema_version": 2,
                "release_id": release_id,
                "release_token": release_token,
                "storage_layout": STREAM_LAYOUT,
                "source": {
                    "file": database_path.name,
                    "sha256": source_sha256,
                    "sqlite_user_version": schema_version,
                },
                "object_prefix": prefix,
                "immutability": (
                    "content-qualified release prefix; never overwrite object keys"
                ),
                "objects": object_entries,
            }
            _write_json(output_dir / "object_manifest.json", object_manifest)
            if download_projection is not None:
                _write_json(output_dir / DOWNLOAD_MANIFEST_NAME, download_projection)
            summary = {
                "schema_version": 2,
                "release_id": release_id,
                "release_token": release_token,
                "source_database_sha256": source_sha256,
                "source_sqlite_user_version": schema_version,
                "input_mode": "sqlite_streaming",
                "storage_layout": STREAM_LAYOUT,
                "publication_scale_supported": True,
                "bounded_memory": {
                    "source_fetch_batch_rows": batch_rows,
                    "coarse_pair_shard_rows": coarse_shard_rows,
                    "max_buffered_pair_rows": max(batch_rows, coarse_shard_rows),
                    "work_database": "on_disk_disposable",
                    "work_database_peak_bytes": work_bytes,
                    "pair_collection_materialized": False,
                },
                "object_shape": {
                    "pair_object_per_cell": False,
                    "pair_objects_emitted": 0,
                    "pair_detail_shards": len(pair_shards),
                    "target_records": len(targets),
                    "drug_records": len(drugs),
                },
                "static_fallback_modified": False,
                "deploy_performed": False,
                "counts": release_manifest["counts"],
                "coverage": coverage,
                "r2_object_count": len(object_entries),
                "r2_bytes": sum(int(row["size_bytes"]) for row in object_entries),
                "entrypoint": f"/releases/{release_token}",
                "object_manifest": "object_manifest.json",
                "download_projection": (
                    DOWNLOAD_MANIFEST_NAME if download_projection is not None else None
                ),
                "download_count": (
                    int(download_projection["upload_count"])
                    if download_projection
                    else 0
                ),
                "download_bytes": (
                    int(download_projection["upload_bytes"])
                    if download_projection
                    else 0
                ),
                "worker": "src/index.mjs",
                "wrangler_template": "wrangler.toml",
            }
            _source_unchanged(database_path, stat_identity)
            _write_json(output_dir / "edge_bundle_summary.json", summary)
            _write_json(
                output_dir / _BUNDLE_MARKER,
                {
                    "schema_version": 1,
                    "generator": "atlas publish edge-bundle",
                    "release_token": release_token,
                    "storage_layout": STREAM_LAYOUT,
                },
            )
            return summary
        finally:
            work.close()
            work_path.unlink(missing_ok=True)
    finally:
        source.close()


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas publish edge-bundle --database",
        description=(
            "Stream a public Atlas SQLite snapshot into coarse read-only edge "
            "objects without emitting one object per pair cell."
        ),
    )
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--source-site-dir", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--download-base-url")
    parser.add_argument("--batch-rows", type=int, default=DEFAULT_BATCH_ROWS)
    parser.add_argument(
        "--coarse-shard-rows", type=int, default=DEFAULT_COARSE_SHARD_ROWS
    )
    parser.add_argument("--max-pairs", type=int, default=DEFAULT_MAX_PAIRS)
    args = parser.parse_args(list(argv) if argv is not None else None)
    try:
        summary = build_streaming_edge_bundle(
            args.database,
            args.out_dir,
            source_site_dir=args.source_site_dir,
            overwrite=args.overwrite,
            download_base_url=args.download_base_url,
            batch_rows=args.batch_rows,
            coarse_shard_rows=args.coarse_shard_rows,
            max_pairs=args.max_pairs,
        )
    except (OSError, ValueError, sqlite3.DatabaseError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["STREAM_LAYOUT", "build_streaming_edge_bundle", "main"]
