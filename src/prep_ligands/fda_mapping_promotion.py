"""Controlled promotion of an independently verified FDA v3 mapping."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import stat
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence, cast


PROMOTION_SCHEMA = "atlas.fda-mapping-promotion.v1"
VERIFIER_SCHEMA = "atlas.fda-terminal-bundle-verification.v1"
V3_RESOLVER_VERSION = "fda-terminal-resolution-v3"
CONFIG_KEY = "FDA_MAPPING_CSV"
_TERMINAL_DISPOSITIONS = frozenset(
    {
        "confirmed_fda_active_ingredient",
        "confirmed_fda_salt_parent_docked",
        "retained_salt_counterion",
        "collapsed_combination_product",
        "additive_excipient",
        "identified_non_fda",
        "incorrect_name_structure_mapping_resolved",
        "prepared_file_unusable",
    }
)
_REQUIRED_VERIFIER_CHECK_NAMES = frozenset(
    {
        "mapping_row_and_rdk_uniqueness",
        "prepared_path_sha_parity",
        "prepared_connectivity_and_score_reuse_partition",
        "deserpidine_regression_sentinel",
        "terminal_disposition_partition",
        "zero_nonterminal_usable_rows",
        "fda_rows_have_structure_and_approval_evidence",
        "fda_canonical_parent_readiness",
        "unsafe_parent_collapse_repair_quarantine",
        "named_library_manifest_counts_and_hashes",
        "fda_source_record_partition",
        "pubchem_manifest_batch_hashes_and_offline_coverage",
    }
)

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_RDK_RE = re.compile(r"rdk[_-]?(\d+)", re.IGNORECASE)
_CANONICAL_RDK_RE = re.compile(r"rdk_\d{7}")
_REUSE_STATUS_BY_CONNECTIVITY = {
    "reuse_supported_coordinate_exact": "coordinate_connectivity_exact",
    "reuse_supported_coordinate_parent": "coordinate_connectivity_parent",
}
_SCORE_LIGAND_FIELDS = (
    "rdk_id",
    "ligand_base",
    "ligand",
    "ligand_file",
    "Ligand_ID",
)
_SCORE_HASH_FIELDS = (
    "ligand_pdbqt_sha256",
    "pdbqt_sha256",
    "ligand_sha256",
    "input_ligand_sha256",
    "historical_pdbqt_sha256",
)
_SCORE_PROVENANCE_FIELDS = (
    "score_source_csv",
    "score_source_sha256",
    "score_source_row_number",
    "fda_mapping_v3_csv",
    "fda_mapping_v3_sha256",
)
_SCORE_JOIN_FIELDS = (
    "candidate_rdk_id",
    "joined_rdk_id",
    "score_pdbqt_sha256",
    "mapping_pdbqt_sha256",
    "exact_identity_join_status",
    "identity_resolution_status",
    "resolved_preferred_name",
    "terminal_disposition",
    "regulatory_status",
    "legacy_score_reuse_status",
    "prepared_connectivity_state",
    "score_reuse_authorized",
)
_REQUIRED_MAPPING_FIELDS = frozenset(
    {
        "rdk_id",
        "terminal_resolution_version",
        "identity_resolution_status",
        "resolved_preferred_name",
        "terminal_disposition",
        "selected_pdbqt_sha256",
    }
)


@dataclass(frozen=True)
class FDAMappingPromotionOutputs:
    promotion_manifest_json: Path
    score_identity_joins_csv: Path | None


@dataclass(frozen=True)
class _PromotionInputs:
    mapping: Path
    verifier: Path
    scores: tuple[Path, ...]
    config: Path | None


@dataclass(frozen=True)
class _PromotionEvidence:
    mapping_fields: list[str]
    mapping_rows: list[dict[str, str]]
    mapping_sha256: str
    verifier: dict[str, Any]
    verifier_sha256: str


@dataclass(frozen=True)
class _ScoreArtifact:
    summary: dict[str, Any]
    payload: bytes | None


@dataclass(frozen=True)
class _ConfigArtifact:
    change: dict[str, Any]
    payload: bytes | None
    original_payload: bytes | None


@dataclass(frozen=True)
class _StagedPromotion:
    score: Path | None
    manifest: Path
    config: Path | None


@dataclass(frozen=True)
class _ScoreSource:
    path: Path
    sha256: str
    fields: list[str]
    rows: list[dict[str, str]]


@dataclass(frozen=True)
class _ScoreJoinResolution:
    candidate: str
    observed_checksum: str
    status: str
    mapping_row: Mapping[str, str] | None


def promote_verified_fda_mapping_v3(
    *,
    mapping_v3_csv: Path,
    verifier_json: Path,
    output_dir: Path,
    score_csvs: Sequence[Path] = (),
    update_config: Path | None = None,
) -> tuple[FDAMappingPromotionOutputs, dict[str, Any]]:
    """Publish a hash-bound pointer manifest without copying the mapping.

    No configuration is changed unless ``update_config`` names an existing file.
    Score identities are attached only when an exact RDK identifier and an exact
    PDBQT checksum both match the verified mapping row.
    """

    inputs = _resolve_promotion_inputs(
        mapping_v3_csv, verifier_json, score_csvs, update_config
    )
    resolved_output_dir = Path(output_dir).expanduser().resolve()
    outputs = _promotion_outputs(resolved_output_dir, bool(inputs.scores))
    _validate_promotion_scope(inputs, outputs)
    evidence = _load_promotion_evidence(inputs)
    score_artifact = _prepare_score_artifact(inputs, outputs, evidence)
    config_artifact = _prepare_config_artifact(inputs.config, inputs.mapping)
    manifest = _build_promotion_manifest(
        inputs, evidence, score_artifact.summary, config_artifact.change
    )
    manifest_payload = _json_bytes(manifest)
    _assert_promotion_inputs_stable(inputs, evidence)
    _publish_promotion(
        resolved_output_dir,
        inputs,
        outputs,
        score_artifact,
        config_artifact,
        manifest_payload,
    )
    return outputs, manifest


def _resolve_promotion_inputs(
    mapping_v3_csv: Path,
    verifier_json: Path,
    score_csvs: Sequence[Path],
    update_config: Path | None,
) -> _PromotionInputs:
    scores = tuple(
        _required_file(path, f"score CSV {index}")
        for index, path in enumerate(score_csvs, start=1)
    )
    config = (
        _required_file(update_config, "explicit config")
        if update_config is not None
        else None
    )
    inputs = _PromotionInputs(
        mapping=_required_file(mapping_v3_csv, "v3 mapping"),
        verifier=_required_file(verifier_json, "independent verifier report"),
        scores=scores,
        config=config,
    )
    _validate_distinct_inputs(inputs)
    return inputs


def _validate_distinct_inputs(inputs: _PromotionInputs) -> None:
    paths = [inputs.mapping, inputs.verifier, *inputs.scores]
    if inputs.config is not None:
        paths.append(inputs.config)
    if len(set(paths)) != len(paths):
        raise ValueError("mapping, verifier, score, and config paths must be distinct")


def _promotion_outputs(
    output_dir: Path, has_scores: bool
) -> FDAMappingPromotionOutputs:
    return FDAMappingPromotionOutputs(
        promotion_manifest_json=output_dir
        / "fda_mapping_v3_promotion_manifest.json",
        score_identity_joins_csv=(
            output_dir / "fda_mapping_v3_score_identity_joins.csv"
            if has_scores
            else None
        ),
    )


def _promotion_output_paths(
    outputs: FDAMappingPromotionOutputs,
) -> tuple[Path, ...]:
    paths = [outputs.promotion_manifest_json]
    if outputs.score_identity_joins_csv is not None:
        paths.append(outputs.score_identity_joins_csv)
    return tuple(paths)


def _validate_promotion_scope(
    inputs: _PromotionInputs, outputs: FDAMappingPromotionOutputs
) -> None:
    input_paths = {inputs.mapping, inputs.verifier, *inputs.scores}
    if inputs.config is not None:
        input_paths.add(inputs.config)
    output_paths = _promotion_output_paths(outputs)
    if input_paths.intersection(output_paths):
        raise ValueError("promotion outputs cannot overwrite an input")
    existing = sorted(str(path) for path in output_paths if path.exists())
    if existing:
        raise FileExistsError(
            "immutable promotion output already exists: " + ", ".join(existing)
        )


def _load_promotion_evidence(inputs: _PromotionInputs) -> _PromotionEvidence:
    mapping_sha256 = _sha256(inputs.mapping)
    mapping_fields, mapping_rows = _load_and_validate_mapping(inputs.mapping)
    verifier = _load_and_validate_verifier(
        inputs.verifier, mapping_sha256, expected_rows=len(mapping_rows)
    )
    return _PromotionEvidence(
        mapping_fields=mapping_fields,
        mapping_rows=mapping_rows,
        mapping_sha256=mapping_sha256,
        verifier=verifier,
        verifier_sha256=_sha256(inputs.verifier),
    )


def _prepare_score_artifact(
    inputs: _PromotionInputs,
    outputs: FDAMappingPromotionOutputs,
    evidence: _PromotionEvidence,
) -> _ScoreArtifact:
    summary = _empty_score_summary(inputs.mapping, evidence.mapping_sha256)
    if not inputs.scores:
        return _ScoreArtifact(summary=summary, payload=None)
    fields, rows, summary = exact_score_identity_joins(
        mapping_path=inputs.mapping,
        mapping_sha256=evidence.mapping_sha256,
        mapping_rows=evidence.mapping_rows,
        score_csvs=inputs.scores,
    )
    payload = _csv_bytes(fields, rows)
    output = outputs.score_identity_joins_csv
    if output is None:
        raise RuntimeError("score output path is unavailable")
    summary["output"] = {
        "path": str(output),
        "sha256": _sha256_bytes(payload),
        "bytes": len(payload),
    }
    return _ScoreArtifact(summary=summary, payload=payload)


def _prepare_config_artifact(
    config_path: Path | None, mapping_path: Path
) -> _ConfigArtifact:
    change, payload, original = _prepare_config_change(config_path, mapping_path)
    return _ConfigArtifact(
        change=change,
        payload=payload,
        original_payload=original,
    )


def _build_promotion_manifest(
    inputs: _PromotionInputs,
    evidence: _PromotionEvidence,
    score_summary: Mapping[str, Any],
    config_change: Mapping[str, Any],
) -> dict[str, Any]:
    mapping_artifact = _verifier_mapping_artifact(evidence.verifier)
    return {
        "schema": PROMOTION_SCHEMA,
        "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "promoted_verified_mapping_pointer",
        "policy": {
            "mapping_copied": False,
            "automatic_config_mutation": False,
            "config_update_requires_explicit_path": True,
            "score_join_requires_exact_rdk_id_and_checksum": True,
            "score_reuse_requires_explicit_status_and_connectivity": True,
        },
        "mapping": {
            "path": str(inputs.mapping),
            "sha256": evidence.mapping_sha256,
            "bytes": inputs.mapping.stat().st_size,
            "rows": len(evidence.mapping_rows),
            "columns": len(evidence.mapping_fields),
            "terminal_resolution_version": V3_RESOLVER_VERSION,
        },
        "independent_verifier": {
            "path": str(inputs.verifier),
            "sha256": evidence.verifier_sha256,
            "schema": evidence.verifier["schema"],
            "passed": True,
            "reported_mapping_path": mapping_artifact.get("path", ""),
            "reported_mapping_sha256": evidence.mapping_sha256,
            "check_count": len(evidence.verifier["checks"]),
        },
        "score_join": score_summary,
        "config_update": dict(config_change),
    }


def _json_bytes(payload: Mapping[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _assert_promotion_inputs_stable(
    inputs: _PromotionInputs, evidence: _PromotionEvidence
) -> None:
    if _sha256(inputs.mapping) != evidence.mapping_sha256:
        raise RuntimeError("v3 mapping changed while preparing promotion")
    if _sha256(inputs.verifier) != evidence.verifier_sha256:
        raise RuntimeError("verifier report changed while preparing promotion")


def _publish_promotion(
    output_dir: Path,
    inputs: _PromotionInputs,
    outputs: FDAMappingPromotionOutputs,
    score: _ScoreArtifact,
    config: _ConfigArtifact,
    manifest_payload: bytes,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    staged_paths: list[Path] = []
    published: list[Path] = []
    try:
        staged = _stage_promotion(
            inputs, outputs, score, config, manifest_payload, staged_paths
        )
        _commit_promotion(inputs, outputs, config, staged, published)
    except Exception:
        _rollback_promotion(inputs, config, published)
        raise
    finally:
        _cleanup_staged(staged_paths)


def _stage_promotion(
    inputs: _PromotionInputs,
    outputs: FDAMappingPromotionOutputs,
    score: _ScoreArtifact,
    config: _ConfigArtifact,
    manifest_payload: bytes,
    staged_paths: list[Path],
) -> _StagedPromotion:
    staged_score = None
    if outputs.score_identity_joins_csv is not None and score.payload is not None:
        staged_score = _stage_bytes(
            outputs.score_identity_joins_csv, score.payload
        )
        staged_paths.append(staged_score)
    staged_manifest = _stage_bytes(outputs.promotion_manifest_json, manifest_payload)
    staged_paths.append(staged_manifest)
    staged_config = _stage_config(inputs.config, config, staged_paths)
    return _StagedPromotion(
        score=staged_score,
        manifest=staged_manifest,
        config=staged_config,
    )


def _stage_config(
    config_path: Path | None,
    config: _ConfigArtifact,
    staged_paths: list[Path],
) -> Path | None:
    if not config.change["changed"] or config_path is None or config.payload is None:
        return None
    staged = _stage_bytes(
        config_path,
        config.payload,
        mode=stat.S_IMODE(config_path.stat().st_mode),
    )
    staged_paths.append(staged)
    return staged


def _commit_promotion(
    inputs: _PromotionInputs,
    outputs: FDAMappingPromotionOutputs,
    config: _ConfigArtifact,
    staged: _StagedPromotion,
    published: list[Path],
) -> None:
    if staged.score is not None and outputs.score_identity_joins_csv is not None:
        _publish_new_file(staged.score, outputs.score_identity_joins_csv)
        published.append(outputs.score_identity_joins_csv)
    if staged.config is not None and inputs.config is not None:
        _assert_config_stable(inputs.config, config.change)
        os.replace(staged.config, inputs.config)
    _publish_new_file(staged.manifest, outputs.promotion_manifest_json)
    published.append(outputs.promotion_manifest_json)


def _assert_config_stable(config_path: Path, change: Mapping[str, Any]) -> None:
    if _sha256(config_path) != change["previous_sha256"]:
        raise RuntimeError("explicit config changed while preparing update")


def _rollback_promotion(
    inputs: _PromotionInputs,
    config: _ConfigArtifact,
    published: Sequence[Path],
) -> None:
    for path in reversed(published):
        path.unlink(missing_ok=True)
    if config.change["changed"] and inputs.config is not None:
        _restore_config(inputs.config, config.original_payload)


def _restore_config(config_path: Path, original_payload: bytes | None) -> None:
    if original_payload is None:
        return
    rollback = _stage_bytes(
        config_path,
        original_payload,
        mode=stat.S_IMODE(config_path.stat().st_mode),
    )
    os.replace(rollback, config_path)


def _cleanup_staged(paths: Sequence[Path]) -> None:
    for path in paths:
        path.unlink(missing_ok=True)


def exact_score_identity_joins(
    *,
    mapping_path: Path,
    mapping_sha256: str,
    mapping_rows: Sequence[Mapping[str, str]],
    score_csvs: Sequence[Path],
) -> tuple[list[str], list[dict[str, Any]], dict[str, Any]]:
    """Join scores only through a matching exact RDK ID and PDBQT SHA-256."""

    _validate_score_mapping_reference(mapping_path, mapping_sha256, mapping_rows)
    mapping_by_id = {row["rdk_id"]: row for row in mapping_rows}
    original_fields, output, sources, statuses = _collect_score_join_rows(
        mapping_path, mapping_sha256, mapping_by_id, score_csvs
    )
    fields = _score_join_output_fields(original_fields)
    summary = _score_join_summary(
        mapping_path, mapping_sha256, sources, output, statuses
    )
    return fields, output, summary


def _validate_score_mapping_reference(
    mapping_path: Path,
    mapping_sha256: str,
    mapping_rows: Sequence[Mapping[str, str]],
) -> None:
    if not _SHA256_RE.fullmatch(mapping_sha256):
        raise ValueError("mapping SHA-256 is malformed")
    if _sha256(mapping_path) != mapping_sha256:
        raise ValueError("mapping SHA-256 does not match mapping_path")
    _validate_score_mapping_ids(mapping_rows)


def _validate_score_mapping_ids(
    mapping_rows: Sequence[Mapping[str, str]],
) -> None:
    mapping_ids = [_clean(row.get("rdk_id")) for row in mapping_rows]
    if (
        any(not _CANONICAL_RDK_RE.fullmatch(rdk_id) for rdk_id in mapping_ids)
        or len(set(mapping_ids)) != len(mapping_ids)
    ):
        raise ValueError("score join requires unique canonical mapping rdk_id values")


def _collect_score_join_rows(
    mapping_path: Path,
    mapping_sha256: str,
    mapping_by_id: Mapping[str, Mapping[str, str]],
    score_csvs: Sequence[Path],
) -> tuple[
    list[str],
    list[dict[str, Any]],
    list[dict[str, Any]],
    Counter[str],
]:
    original_fields: list[str] = []
    output: list[dict[str, Any]] = []
    sources: list[dict[str, Any]] = []
    statuses: Counter[str] = Counter()
    for score_path in score_csvs:
        source = _read_stable_score_source(score_path)
        _extend_unique_fields(original_fields, source.fields)
        sources.append(_score_source_summary(source))
        joined_rows = _join_score_source_rows(
            source, mapping_path, mapping_sha256, mapping_by_id
        )
        output.extend(joined_rows)
        statuses.update(row["exact_identity_join_status"] for row in joined_rows)
    return original_fields, output, sources, statuses


def _read_stable_score_source(path: Path) -> _ScoreSource:
    source_sha = _sha256(path)
    fields, rows = _read_csv(path)
    if _sha256(path) != source_sha:
        raise RuntimeError(f"score CSV changed while reading: {path}")
    return _ScoreSource(path=path, sha256=source_sha, fields=fields, rows=rows)


def _extend_unique_fields(existing: list[str], incoming: Sequence[str]) -> None:
    for field in incoming:
        if field not in existing:
            existing.append(field)


def _score_source_summary(source: _ScoreSource) -> dict[str, Any]:
    return {
        "path": str(source.path),
        "sha256": source.sha256,
        "bytes": source.path.stat().st_size,
        "rows": len(source.rows),
    }


def _join_score_source_rows(
    source: _ScoreSource,
    mapping_path: Path,
    mapping_sha256: str,
    mapping_by_id: Mapping[str, Mapping[str, str]],
) -> list[dict[str, Any]]:
    return [
        _joined_score_source_row(
            source, row_number, row, mapping_path, mapping_sha256, mapping_by_id
        )
        for row_number, row in enumerate(source.rows, start=1)
    ]


def _joined_score_source_row(
    source: _ScoreSource,
    row_number: int,
    row: Mapping[str, str],
    mapping_path: Path,
    mapping_sha256: str,
    mapping_by_id: Mapping[str, Mapping[str, str]],
) -> dict[str, Any]:
    return {
        **row,
        "score_source_csv": str(source.path),
        "score_source_sha256": source.sha256,
        "score_source_row_number": row_number,
        "fda_mapping_v3_csv": str(mapping_path),
        "fda_mapping_v3_sha256": mapping_sha256,
        **_exact_score_join(row, mapping_by_id),
    }


def _score_join_output_fields(original_fields: Sequence[str]) -> list[str]:
    reserved = [*_SCORE_PROVENANCE_FIELDS, *_SCORE_JOIN_FIELDS]
    return [
        *_SCORE_PROVENANCE_FIELDS,
        *(field for field in original_fields if field not in reserved),
        *_SCORE_JOIN_FIELDS,
    ]


def _score_join_summary(
    mapping_path: Path,
    mapping_sha256: str,
    sources: Sequence[Mapping[str, Any]],
    output: Sequence[Mapping[str, Any]],
    statuses: Mapping[str, int],
) -> dict[str, Any]:
    return {
        "requested": True,
        "mapping_sha256_stamped_on_every_row": True,
        "join_keys": ["rdk_id", "historical_ligand_pdbqt_sha256"],
        "reuse_authorization_pairs": dict(_REUSE_STATUS_BY_CONNECTIVITY),
        "mapping": {"path": str(mapping_path), "sha256": mapping_sha256},
        "sources": sources,
        "rows": len(output),
        "status_counts": dict(sorted(statuses.items())),
    }


def _exact_score_join(
    score_row: Mapping[str, str],
    mapping_by_id: Mapping[str, Mapping[str, str]],
) -> dict[str, str]:
    resolution = _resolve_score_join(score_row, mapping_by_id)
    return {
        "candidate_rdk_id": resolution.candidate,
        "joined_rdk_id": _joined_rdk_id(resolution),
        "score_pdbqt_sha256": resolution.observed_checksum,
        "mapping_pdbqt_sha256": _candidate_mapping_checksum(
            resolution.candidate, mapping_by_id
        ),
        "exact_identity_join_status": resolution.status,
        **_joined_identity_fields(resolution.mapping_row),
    }


def _resolve_score_join(
    score_row: Mapping[str, str],
    mapping_by_id: Mapping[str, Mapping[str, str]],
) -> _ScoreJoinResolution:
    ids, id_invalid = _score_rdk_ids(score_row)
    hashes, hash_invalid = _score_hashes(score_row)
    candidate = _single_value(ids)
    observed = _single_value(hashes)
    id_status = _identifier_block_status(ids, id_invalid, candidate, mapping_by_id)
    if id_status:
        return _blocked_score_resolution(candidate, observed, id_status)
    hash_status = _checksum_block_status(hashes, hash_invalid)
    if hash_status:
        return _blocked_score_resolution(candidate, observed, hash_status)
    return _mapping_checksum_resolution(candidate, observed, mapping_by_id)


def _single_value(values: set[str]) -> str:
    return next(iter(values)) if len(values) == 1 else ""


def _identifier_block_status(
    ids: set[str],
    invalid: bool,
    candidate: str,
    mapping_by_id: Mapping[str, Mapping[str, str]],
) -> str:
    if len(ids) > 1:
        return "blocked_conflicting_rdk_ids"
    if not ids:
        return "blocked_invalid_rdk_id" if invalid else "blocked_missing_rdk_id"
    if candidate not in mapping_by_id:
        return "blocked_rdk_id_not_in_mapping"
    return ""


def _checksum_block_status(hashes: set[str], invalid: bool) -> str:
    if len(hashes) > 1:
        return "blocked_conflicting_score_checksums"
    if invalid:
        return "blocked_invalid_score_checksum"
    if not hashes:
        return "blocked_missing_score_checksum"
    return ""


def _mapping_checksum_resolution(
    candidate: str,
    observed: str,
    mapping_by_id: Mapping[str, Mapping[str, str]],
) -> _ScoreJoinResolution:
    mapping_row = mapping_by_id[candidate]
    expected = _clean(mapping_row.get("selected_pdbqt_sha256")).lower()
    if not _SHA256_RE.fullmatch(expected):
        return _blocked_score_resolution(
            candidate, observed, "blocked_mapping_checksum_unavailable"
        )
    if expected != observed:
        return _blocked_score_resolution(candidate, observed, "blocked_checksum_mismatch")
    return _ScoreJoinResolution(
        candidate=candidate,
        observed_checksum=observed,
        status="joined_exact_rdk_id_and_checksum",
        mapping_row=mapping_row,
    )


def _blocked_score_resolution(
    candidate: str, observed: str, status: str
) -> _ScoreJoinResolution:
    return _ScoreJoinResolution(
        candidate=candidate,
        observed_checksum=observed,
        status=status,
        mapping_row=None,
    )


def _joined_rdk_id(resolution: _ScoreJoinResolution) -> str:
    return resolution.candidate if resolution.mapping_row is not None else ""


def _candidate_mapping_checksum(
    candidate: str,
    mapping_by_id: Mapping[str, Mapping[str, str]],
) -> str:
    row = mapping_by_id.get(candidate)
    return _clean(row.get("selected_pdbqt_sha256")).lower() if row else ""


def _joined_identity_fields(
    mapping_row: Mapping[str, str] | None,
) -> dict[str, str]:
    if mapping_row is None:
        return _empty_identity_fields()
    reuse_status = _clean(mapping_row.get("legacy_score_reuse_status"))
    connectivity_state = _clean(mapping_row.get("prepared_connectivity_state"))
    reuse_authorized = bool(
        _REUSE_STATUS_BY_CONNECTIVITY.get(reuse_status) == connectivity_state
    )
    return {
        "identity_resolution_status": _clean(
            mapping_row.get("identity_resolution_status")
        ),
        "resolved_preferred_name": _clean(mapping_row.get("resolved_preferred_name")),
        "terminal_disposition": _clean(mapping_row.get("terminal_disposition")),
        "regulatory_status": _clean(mapping_row.get("regulatory_status")),
        "legacy_score_reuse_status": reuse_status,
        "prepared_connectivity_state": connectivity_state,
        "score_reuse_authorized": str(reuse_authorized).lower(),
    }


def _empty_identity_fields() -> dict[str, str]:
    return {
        "identity_resolution_status": "",
        "resolved_preferred_name": "",
        "terminal_disposition": "",
        "regulatory_status": "",
        "legacy_score_reuse_status": "",
        "prepared_connectivity_state": "",
        "score_reuse_authorized": "false",
    }


def _score_rdk_ids(row: Mapping[str, str]) -> tuple[set[str], bool]:
    values = [_clean(row.get(field)) for field in _SCORE_LIGAND_FIELDS]
    populated = [value for value in values if value]
    parsed = {_parse_rdk_id(value) for value in populated}
    parsed.discard("")
    return parsed, bool(populated and not parsed)


def _score_hashes(row: Mapping[str, str]) -> tuple[set[str], bool]:
    values = [_clean(row.get(field)).lower() for field in _SCORE_HASH_FIELDS]
    populated = [value for value in values if value]
    valid = {value for value in populated if _SHA256_RE.fullmatch(value)}
    invalid = any(not _SHA256_RE.fullmatch(value) for value in populated)
    return valid, invalid


def _parse_rdk_id(value: str) -> str:
    stem = Path(value).name
    for suffix in (".pdbqt", ".mol2", ".sdf", ".pdb"):
        if stem.casefold().endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    match = _RDK_RE.fullmatch(stem)
    return f"rdk_{int(match.group(1)):07d}" if match else ""


def _load_and_validate_mapping(
    path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    fields, rows = _read_csv(path)
    _require_mapping_fields(fields)
    _require_mapping_rows(rows)
    _validate_mapping_ids(rows)
    _validate_mapping_versions(rows)
    _validate_mapping_terminal_statuses(rows)
    _validate_mapping_dispositions(rows)
    return fields, rows


def _require_mapping_fields(fields: Sequence[str]) -> None:
    missing = sorted(_REQUIRED_MAPPING_FIELDS - set(fields))
    if missing:
        raise ValueError("v3 mapping is missing columns: " + ", ".join(missing))


def _require_mapping_rows(rows: Sequence[Mapping[str, str]]) -> None:
    if not rows:
        raise ValueError("v3 mapping has no rows")


def _validate_mapping_ids(rows: Sequence[Mapping[str, str]]) -> None:
    ids = [_clean(row.get("rdk_id")) for row in rows]
    if any(not _CANONICAL_RDK_RE.fullmatch(value) for value in ids):
        raise ValueError("v3 mapping contains a missing or malformed rdk_id")
    if len(set(ids)) != len(ids):
        raise ValueError("v3 mapping contains duplicate rdk_id values")


def _validate_mapping_versions(rows: Sequence[Mapping[str, str]]) -> None:
    if any(
        _clean(row.get("terminal_resolution_version")) != V3_RESOLVER_VERSION
        for row in rows
    ):
        raise ValueError("mapping is not uniformly terminal resolver v3 output")


def _validate_mapping_terminal_statuses(
    rows: Sequence[Mapping[str, str]],
) -> None:
    if any(
        not _clean(row.get("identity_resolution_status")).startswith("terminal_")
        for row in rows
    ):
        raise ValueError("v3 mapping contains a nonterminal identity status")


def _validate_mapping_dispositions(rows: Sequence[Mapping[str, str]]) -> None:
    observed = {_clean(row.get("terminal_disposition")) for row in rows}
    unexpected = sorted(observed - _TERMINAL_DISPOSITIONS)
    if unexpected:
        raise ValueError(
            "v3 mapping contains unsupported terminal dispositions: "
            + ", ".join(value or "<blank>" for value in unexpected)
        )


def _load_and_validate_verifier(
    path: Path,
    mapping_sha256: str,
    *,
    expected_rows: int,
) -> dict[str, Any]:
    report = _read_verifier_json(path)
    _validate_verifier_status(report)
    _validate_verifier_checks(report, expected_rows=expected_rows)
    mapping_artifact = _verifier_mapping_artifact(report)
    _validate_verifier_mapping_sha(mapping_artifact, mapping_sha256)
    return report


def _read_verifier_json(path: Path) -> dict[str, Any]:
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid verifier JSON: {exc}") from exc
    if not isinstance(report, dict):
        raise ValueError("verifier report must be a JSON object")
    return report


def _validate_verifier_status(report: Mapping[str, Any]) -> None:
    if report.get("schema") != VERIFIER_SCHEMA:
        raise ValueError("verifier report has an unsupported schema")
    if report.get("passed") is not True:
        raise ValueError("verifier report did not pass")


def _validate_verifier_checks(
    report: Mapping[str, Any],
    *,
    expected_rows: int,
) -> None:
    checks = _require_verifier_checks(report)
    _validate_verifier_check_entries(checks)
    typed_checks = cast(list[dict[str, Any]], checks)
    _validate_verifier_check_names(typed_checks)
    _validate_verifier_mapping_row_count(typed_checks, expected_rows)


def _require_verifier_checks(report: Mapping[str, Any]) -> list[Any]:
    checks = report.get("checks")
    if not isinstance(checks, list):
        raise ValueError("verifier report has no independent checks")
    if len(checks) != len(_REQUIRED_VERIFIER_CHECK_NAMES):
        raise ValueError(
            "verifier report must contain exactly "
            f"{len(_REQUIRED_VERIFIER_CHECK_NAMES)} independent checks"
        )
    return checks


def _valid_verifier_check_entry(check: object) -> bool:
    if not isinstance(check, dict):
        return False
    return check.get("passed") is True and bool(_clean(check.get("name")))


def _validate_verifier_check_entries(checks: Sequence[object]) -> None:
    if not all(_valid_verifier_check_entry(check) for check in checks):
        raise ValueError("verifier report contains a failed or malformed check")


def _validate_verifier_check_names(checks: Sequence[Mapping[str, Any]]) -> None:
    raw_names = [check.get("name") for check in checks]
    names = [_clean(name) for name in raw_names]
    noncanonical_names = any(
        raw_name != clean_name for raw_name, clean_name in zip(raw_names, names)
    )
    duplicate_names = sorted(
        name for name, count in Counter(names).items() if count != 1
    )
    missing = sorted(_REQUIRED_VERIFIER_CHECK_NAMES - set(names))
    unexpected = sorted(set(names) - _REQUIRED_VERIFIER_CHECK_NAMES)
    if noncanonical_names or duplicate_names or missing or unexpected:
        raise ValueError(
            "verifier check-name contract mismatch: "
            f"duplicates={duplicate_names}, missing={missing}, unexpected={unexpected}"
        )


def _validate_verifier_mapping_row_count(
    checks: Sequence[Mapping[str, Any]], expected_rows: int
) -> None:
    uniqueness = next(
        check
        for check in checks
        if check["name"] == "mapping_row_and_rdk_uniqueness"
    )
    metrics = uniqueness.get("metrics")
    if not isinstance(metrics, dict) or metrics.get("mapping_rows") != expected_rows:
        raise ValueError(
            "verifier mapping-row count does not match the supplied v3 mapping"
        )


def _verifier_mapping_artifact(report: Mapping[str, Any]) -> dict[str, Any]:
    artifacts = report.get("artifacts")
    if not isinstance(artifacts, dict):
        return {}
    mapping = artifacts.get("mapping")
    return mapping if isinstance(mapping, dict) else {}


def _validate_verifier_mapping_sha(
    mapping_artifact: Mapping[str, Any], mapping_sha256: str
) -> None:
    reported_sha = _clean(mapping_artifact.get("sha256")).lower()
    if not _SHA256_RE.fullmatch(reported_sha):
        raise ValueError("verifier report has no valid mapping SHA-256")
    if reported_sha != mapping_sha256:
        raise ValueError(
            "verifier mapping SHA-256 does not match the supplied v3 mapping"
        )


def _prepare_config_change(
    config_path: Path | None, mapping_path: Path
) -> tuple[dict[str, Any], bytes | None, bytes | None]:
    if config_path is None:
        return {
            "requested": False,
            "changed": False,
            "path": "",
            "key": CONFIG_KEY,
        }, None, None
    original = config_path.read_bytes()
    try:
        text = original.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("explicit config must be UTF-8 text") from exc
    lines = text.splitlines(keepends=True)
    key_re = re.compile(rf"^\s*{re.escape(CONFIG_KEY)}\s*=")
    positions = [index for index, line in enumerate(lines) if key_re.match(line)]
    if len(positions) > 1:
        raise ValueError(f"explicit config contains duplicate {CONFIG_KEY} entries")
    replacement = f"{CONFIG_KEY}={mapping_path}\n"
    if positions:
        ending = "\r\n" if lines[positions[0]].endswith("\r\n") else "\n"
        replacement = replacement.rstrip("\n") + ending
        lines[positions[0]] = replacement
    else:
        if lines and not lines[-1].endswith(("\n", "\r")):
            lines[-1] += "\n"
        lines.append(replacement)
    updated = "".join(lines).encode("utf-8")
    return {
        "requested": True,
        "changed": updated != original,
        "path": str(config_path),
        "key": CONFIG_KEY,
        "value": str(mapping_path),
        "previous_sha256": _sha256_bytes(original),
        "updated_sha256": _sha256_bytes(updated),
    }, updated, original


def _empty_score_summary(mapping_path: Path, mapping_sha256: str) -> dict[str, Any]:
    return {
        "requested": False,
        "mapping_sha256_stamped_on_every_row": True,
        "join_keys": ["rdk_id", "historical_ligand_pdbqt_sha256"],
        "reuse_authorization_pairs": dict(_REUSE_STATUS_BY_CONNECTIVITY),
        "mapping": {"path": str(mapping_path), "sha256": mapping_sha256},
        "sources": [],
        "rows": 0,
        "status_counts": {},
    }


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                raise ValueError(f"CSV has no header: {path}")
            if len(set(reader.fieldnames)) != len(reader.fieldnames):
                raise ValueError(f"CSV has duplicate header fields: {path}")
            return list(reader.fieldnames), [
                {key: _clean(value) for key, value in row.items()}
                for row in reader
            ]
    except csv.Error as exc:
        raise ValueError(f"invalid CSV {path}: {exc}") from exc


def _csv_bytes(
    fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(fields), extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def _stage_bytes(path: Path, payload: bytes, mode: int | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    staged = path.with_name(f".{path.name}.{os.getpid()}.part")
    if staged.exists():
        raise FileExistsError(f"staging path already exists: {staged}")
    try:
        with staged.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        if mode is not None:
            os.chmod(staged, mode)
    except Exception:
        staged.unlink(missing_ok=True)
        raise
    return staged


def _publish_new_file(staged: Path, destination: Path) -> None:
    os.link(staged, destination)
    try:
        staged.unlink()
    except Exception:
        destination.unlink(missing_ok=True)
        raise


def _required_file(path: Path, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} is not a file: {resolved}")
    return resolved


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _clean(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.casefold() in {"nan", "none", "null"} else text


__all__ = [
    "FDAMappingPromotionOutputs",
    "exact_score_identity_joins",
    "promote_verified_fda_mapping_v3",
]
