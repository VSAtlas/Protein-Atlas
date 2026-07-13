"""Fail-closed identity reconciliation for Phase 1 SPD model tables."""

from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd

from src import path_router


RECONCILIATION_VERSION = "atlas-spd-phase1-identity-v1"
SUMMARY_SCHEMA = "atlas.spd-phase1-identity-reconciliation-summary.v1"
AUDIT_SCHEMA = "atlas.spd-phase1-identity-reconciliation-row-audit.v1"

DEFAULT_IDENTITY_KEY_COLUMNS = (
    "rdk_id",
    "ligand_base",
    "_join_rdk",
    "ligand_rdk_id",
    "rdk",
)
FALLBACK_IDENTITY_KEY_COLUMNS = ("ligand", "ligand_file", "Ligand_ID")
DEFAULT_SPD_INCHIKEY_COLUMNS = (
    "spd_inchikey",
    "spd_drug_inchikey",
    "spd_identity_inchikey",
)
NAME_COLUMNS = (
    "drug_id",
    "drug",
    "drug_name",
    "display_name",
    "generic_name",
    "compound_name",
    "ligand_name",
)

# These columns are owned by SPD assay/exposure evidence. External PK context,
# source-specific external labels, combined labels, and model scores are omitted
# deliberately.
SPD_OWNED_VALUE_COLUMNS = frozenset(
    {
        "spd_binding_label",
        "spd_binding_ml_label",
        "spd_activity_label",
        "spd_activity_ml_label",
        "spd_exposure_label",
        "spd_exposure_ml_label",
        "spd_exposure_relevant",
        "spd_exposure_weak",
        "spd_exposure_unlikely",
        "ml_binary_label",
        "spd_ac50",
        "spd_ac50_nm",
        "spd_ac50_um",
        "ac50_nm",
        "ac50_um",
        "free_cmax_nm",
        "free_cmax_um",
        "free_cmax",
        "total_cmax_nm",
        "total_cmax_um",
        "combined_free_cmax_um",
        "combined_cmax_um",
        "cmax_um",
        "exposure_margin",
        "spd_exposure_margin",
        "spd_activity_relation",
        "activity_relation",
    }
)
SPD_EVIDENCE_PRESENCE_COLUMNS = frozenset(
    column for column in SPD_OWNED_VALUE_COLUMNS if column.startswith("spd_")
)

_RDK_RE = re.compile(
    r"rdk[_-]?(\d+)(?:_(?:dud_)?stage\d+)?",
    flags=re.IGNORECASE,
)
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_TRUE_VALUES = frozenset({"1", "1.0", "true", "yes", "y"})
_REUSE_STATUS_BY_CONNECTIVITY = {
    "reuse_supported_coordinate_exact": "coordinate_connectivity_exact",
    "reuse_supported_coordinate_parent": "coordinate_connectivity_parent",
}

_REQUIRED_MAPPING_COLUMNS = frozenset(
    {
        "rdk_id",
        "mapping_row_number",
        "display_name",
        "generic_name",
        "preferred_identity",
        "preferred_identity_source",
        "preferred_identity_source_id",
        "identity_resolution_status",
        "identity_structure_validated",
        "identity_exact_inchikey",
        "identity_parent_inchikey",
        "terminal_disposition",
        "terminal_reason_codes",
        "regulatory_status",
        "canonical_parent_id",
        "selected_pdbqt_path",
        "selected_pdbqt_sha256",
        "canonical_pdbqt_path",
        "canonical_pdbqt_sha256",
        "prepared_connectivity_state",
        "legacy_score_reuse_status",
    }
)

_MAPPING_ATTACHMENTS = {
    "rdk_id": "canonical_identity_rdk_id",
    "mapping_row_number": "identity_mapping_row_number",
    "display_name": "canonical_identity_display_name",
    "generic_name": "canonical_identity_generic_name",
    "preferred_identity": "canonical_identity_name",
    "preferred_identity_source": "canonical_identity_source",
    "preferred_identity_source_id": "canonical_identity_source_id",
    "identity_resolution_status": "canonical_identity_resolution_status",
    "identity_structure_validated": "canonical_identity_structure_validated",
    "identity_exact_inchikey": "canonical_identity_exact_inchikey",
    "identity_parent_inchikey": "canonical_identity_parent_inchikey",
    "identity_structure_relation": "canonical_identity_structure_relation",
    "identity_change_applied": "canonical_identity_change_applied",
    "identity_change_reason_codes": "canonical_identity_change_reason_codes",
    "identity_evidence_codes": "canonical_identity_evidence_codes",
    "approval_evidence_codes": "canonical_identity_approval_evidence_codes",
    "terminal_disposition": "canonical_identity_terminal_disposition",
    "terminal_reason_codes": "canonical_identity_terminal_reason_codes",
    "terminal_resolution_version": "canonical_identity_mapping_version",
    "regulatory_status": "canonical_identity_regulatory_status",
    "regulatory_evidence": "canonical_identity_regulatory_evidence",
    "canonical_parent_id": "canonical_identity_parent_id",
    "canonical_materialization_status": "canonical_identity_materialization_status",
    "selected_pdbqt_path": "identity_mapping_selected_pdbqt_path",
    "selected_pdbqt_sha256": "identity_mapping_selected_pdbqt_sha256",
    "canonical_pdbqt_path": "identity_mapping_canonical_pdbqt_path",
    "canonical_pdbqt_sha256": "identity_mapping_canonical_pdbqt_sha256",
    "prepared_connectivity_state": "canonical_identity_prepared_connectivity_state",
    "legacy_score_reuse_status": "canonical_identity_legacy_score_reuse_status",
    "selected_pdbqt_checksum_status_v3": "canonical_identity_selected_checksum_status",
    "legacy_reconciliation_version": "canonical_identity_legacy_mapping_version",
    "legacy_display_name": "canonical_identity_mapping_legacy_display_name",
    "legacy_generic_name": "canonical_identity_mapping_legacy_generic_name",
}


class IdentityReconciliationError(ValueError):
    """Raised when reconciliation inputs or output paths are unsafe."""


@dataclass(frozen=True)
class ReconciliationOutputs:
    reconciled_csv: Path
    row_audit_csv: Path
    summary_json: Path


@dataclass(frozen=True)
class _EvidenceSpec:
    slot: str
    expected_hash_column: str
    path_columns: tuple[str, ...]
    hash_columns: tuple[str, ...]


_EVIDENCE_SPECS = (
    _EvidenceSpec(
        slot="historical",
        expected_hash_column="selected_pdbqt_sha256",
        path_columns=(
            "historical_pdbqt_path",
            "score_pdbqt_path",
            "ligand_pdbqt_path",
            "pdbqt_path",
            "ligand_path",
        ),
        hash_columns=(
            "historical_pdbqt_sha256",
            "score_pdbqt_sha256",
            "ligand_pdbqt_sha256",
            "pdbqt_sha256",
            "ligand_sha256",
            "input_ligand_sha256",
        ),
    ),
    _EvidenceSpec(
        slot="selected",
        expected_hash_column="selected_pdbqt_sha256",
        path_columns=("selected_pdbqt_path", "current_pdbqt_path"),
        hash_columns=(
            "selected_pdbqt_sha256",
            "mapping_pdbqt_sha256",
            "current_pdbqt_sha256",
            "observed_current_pdbqt_sha256",
            "verified_pdbqt_sha256",
        ),
    ),
    _EvidenceSpec(
        slot="canonical",
        expected_hash_column="canonical_pdbqt_sha256",
        path_columns=("canonical_pdbqt_path", "named_parent_pdbqt_path"),
        hash_columns=("canonical_pdbqt_sha256", "named_parent_pdbqt_sha256"),
    ),
)


def default_fda_mapping_path() -> Path:
    """Return the canonical FDA mapping path through the repository router."""

    repo_root = Path(path_router.__file__).resolve().parents[2]
    return repo_root / "chemdb" / "data" / "fda_mapping_from_pdbqt.csv"


def default_output_paths(output_csv: str | Path) -> ReconciliationOutputs:
    """Derive the row-audit and summary paths beside a reconciled CSV."""

    reconciled = Path(output_csv)
    stem = reconciled.stem
    return ReconciliationOutputs(
        reconciled_csv=reconciled,
        row_audit_csv=reconciled.with_name(f"{stem}_identity_row_audit.csv"),
        summary_json=reconciled.with_name(f"{stem}_identity_summary.json"),
    )


def normalize_rdk_identity(value: object) -> str:
    """Normalize an exact RDK filename/base identifier without name fallback."""

    text = _clean(value)
    if not text:
        return ""
    stem = Path(text).name
    for suffix in (".pdbqt", ".mol2", ".sdf", ".pdb"):
        if stem.casefold().endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    match = _RDK_RE.fullmatch(stem)
    return f"rdk_{int(match.group(1)):07d}" if match else ""


def reconcile_spd_phase1_identity(
    input_csv: str | Path,
    output_csv: str | Path,
    *,
    mapping_csv: str | Path | None = None,
    row_audit_csv: str | Path | None = None,
    summary_json: str | Path | None = None,
    identity_key_columns: Sequence[str] | None = None,
    spd_inchikey_columns: Sequence[str] | None = None,
    additional_spd_null_columns: Sequence[str] = (),
    fail_closed: bool = True,
) -> dict[str, Any]:
    """Reconcile a Phase 1 SPD table against the canonical FDA/RDK mapping.

    Rows are joined only by normalized RDK identifiers. Name fields are updated
    in place only after that join, while their original values are retained as
    ``legacy_*`` columns. A connectivity mismatch between a validated mapping
    InChIKey and an SPD InChIKey is the only condition that nulls SPD-owned
    measurements or labels. Scores, external labels, and combined labels are
    never used as join keys and are never modified.
    """

    source_path = _required_file(Path(input_csv), "Phase 1 SPD input CSV")
    map_path = _required_file(
        Path(mapping_csv) if mapping_csv is not None else default_fda_mapping_path(),
        "canonical FDA mapping CSV",
    )
    outputs = _resolve_outputs(output_csv, row_audit_csv, summary_json)
    _validate_distinct_paths(source_path, map_path, outputs)

    input_sha256 = _sha256_file(source_path)
    mapping_sha256 = _sha256_file(map_path)
    source = pd.read_csv(source_path, dtype=object, low_memory=False)
    mapping = _load_mapping(map_path)
    mapping_by_key = {
        str(row["_normalized_rdk_id"]): row for row in mapping.to_dict(orient="records")
    }

    requested_keys = tuple(identity_key_columns or DEFAULT_IDENTITY_KEY_COLUMNS)
    key_columns = _available_key_columns(source, requested_keys)
    if not key_columns:
        available = ", ".join(str(column) for column in source.columns)
        raise IdentityReconciliationError(
            "input has no usable exact RDK identity column; expected one of "
            f"{requested_keys!r}. Available columns: {available}"
        )

    requested_spd_keys = tuple(spd_inchikey_columns or DEFAULT_SPD_INCHIKEY_COLUMNS)
    available_spd_keys = tuple(
        column for column in requested_spd_keys if column in source.columns
    )
    spd_owned_columns = _select_spd_owned_columns(
        source.columns, additional_spd_null_columns
    )
    name_columns = tuple(column for column in NAME_COLUMNS if column in source.columns)
    reserved = _reserved_output_columns(name_columns)
    collisions = sorted(reserved.intersection(source.columns))
    if collisions:
        raise IdentityReconciliationError(
            "input already contains reconciliation output columns: "
            + ", ".join(collisions)
        )

    protected_columns = _protected_columns(source.columns)
    protected_snapshot = source[list(protected_columns)].copy(deep=True)
    reconciled = source.copy(deep=True)
    for column in name_columns:
        reconciled[f"legacy_{column}"] = source[column]

    row_contexts = _build_row_contexts(
        source,
        mapping_by_key=mapping_by_key,
        key_columns=key_columns,
        spd_inchikey_columns=available_spd_keys,
        spd_owned_columns=spd_owned_columns,
        source_dir=source_path.parent,
    )
    context = pd.DataFrame(row_contexts, index=source.index)
    for column in context.columns:
        reconciled[column] = context[column]

    _attach_mapping_fields(
        reconciled,
        mapping_by_key,
        mapping_path=map_path,
        mapping_sha256=mapping_sha256,
    )
    _reconcile_name_fields(reconciled, name_columns)
    original_spd_values = _null_mismatched_spd_values(
        reconciled, source, spd_owned_columns
    )
    reconciled["identity_original_spd_values_json"] = original_spd_values
    _apply_identity_training_gate(reconciled, source)
    _assert_protected_columns_unchanged(
        reconciled, protected_snapshot, protected_columns
    )
    if len(reconciled) != len(source) or not reconciled.index.equals(source.index):
        raise IdentityReconciliationError(
            "row count or order changed during identity reconciliation"
        )

    audit = _build_row_audit(reconciled, source.columns, name_columns)
    _assert_inputs_stable(source_path, input_sha256, map_path, mapping_sha256)
    _atomic_to_csv(reconciled, outputs.reconciled_csv)
    _atomic_to_csv(audit, outputs.row_audit_csv)

    summary = _build_summary(
        source_path=source_path,
        input_sha256=input_sha256,
        source=source,
        map_path=map_path,
        mapping_sha256=mapping_sha256,
        mapping=mapping,
        outputs=outputs,
        reconciled=reconciled,
        key_columns=key_columns,
        spd_inchikey_columns=available_spd_keys,
        spd_owned_columns=spd_owned_columns,
        name_columns=name_columns,
        protected_columns=protected_columns,
        fail_closed=fail_closed,
    )
    summary["outputs"]["reconciled_csv"]["sha256"] = _sha256_file(
        outputs.reconciled_csv
    )
    summary["outputs"]["row_audit_csv"]["sha256"] = _sha256_file(outputs.row_audit_csv)
    _atomic_write_json(summary, outputs.summary_json)
    return summary


def _clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _truthy(value: object) -> bool:
    return _clean(value).casefold() in _TRUE_VALUES


def _explicit_non_fda_probe(row: pd.Series) -> bool:
    """Return true only for rows explicitly quarantined from FDA/SPD truth."""

    domains = " ".join(
        _clean(row.get(column)).casefold()
        for column in (
            "external_addon_ligand_domain",
            "ligand_domain",
            "training_domain",
        )
        if column in row.index
    )
    domain_is_probe = "non_fda" in domains or "probe" in domains
    production_values = [
        row.get(column)
        for column in (
            "external_addon_production_truth_allowed",
            "production_truth_allowed",
            "primary_fda_claim_allowed",
        )
        if column in row.index and _clean(row.get(column))
    ]
    explicitly_not_production = bool(production_values) and not any(
        _truthy(value) for value in production_values
    )
    return domain_is_probe and explicitly_not_production


def _required_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise IdentityReconciliationError(f"{label} does not exist: {resolved}")
    return resolved


def _resolve_outputs(
    output_csv: str | Path,
    row_audit_csv: str | Path | None,
    summary_json: str | Path | None,
) -> ReconciliationOutputs:
    defaults = default_output_paths(output_csv)
    return ReconciliationOutputs(
        reconciled_csv=Path(output_csv).expanduser().resolve(),
        row_audit_csv=(
            Path(row_audit_csv).expanduser().resolve()
            if row_audit_csv is not None
            else defaults.row_audit_csv.expanduser().resolve()
        ),
        summary_json=(
            Path(summary_json).expanduser().resolve()
            if summary_json is not None
            else defaults.summary_json.expanduser().resolve()
        ),
    )


def _validate_distinct_paths(
    source_path: Path,
    mapping_path: Path,
    outputs: ReconciliationOutputs,
) -> None:
    paths = {
        "input": source_path,
        "mapping": mapping_path,
        "reconciled": outputs.reconciled_csv,
        "row audit": outputs.row_audit_csv,
        "summary": outputs.summary_json,
    }
    by_path: dict[Path, list[str]] = {}
    for label, path in paths.items():
        by_path.setdefault(path, []).append(label)
    collisions = [labels for labels in by_path.values() if len(labels) > 1]
    if collisions:
        rendered = "; ".join("/".join(labels) for labels in collisions)
        raise IdentityReconciliationError(
            f"input, mapping, and output paths must be distinct: {rendered}"
        )


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_mapping(path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False)
    missing = sorted(_REQUIRED_MAPPING_COLUMNS.difference(mapping.columns))
    if missing:
        raise IdentityReconciliationError(
            f"canonical mapping is missing required columns: {', '.join(missing)}"
        )
    mapping["_normalized_rdk_id"] = mapping["rdk_id"].map(normalize_rdk_identity)
    invalid = mapping["_normalized_rdk_id"].eq("")
    if invalid.any():
        rows = ", ".join(str(index + 2) for index in mapping.index[invalid][:10])
        raise IdentityReconciliationError(
            f"canonical mapping has invalid RDK identifiers at CSV rows: {rows}"
        )
    duplicated = mapping["_normalized_rdk_id"].duplicated(keep=False)
    if duplicated.any():
        values = sorted(mapping.loc[duplicated, "_normalized_rdk_id"].unique())
        raise IdentityReconciliationError(
            "canonical mapping has ambiguous normalized RDK identifiers: "
            + ", ".join(values[:10])
        )
    for column in ("selected_pdbqt_sha256", "canonical_pdbqt_sha256"):
        values = mapping[column].map(_clean).str.lower()
        invalid_hash = values.ne("") & ~values.str.fullmatch(_SHA256_RE)
        if invalid_hash.any():
            rows = ", ".join(
                str(index + 2) for index in mapping.index[invalid_hash][:10]
            )
            raise IdentityReconciliationError(
                f"canonical mapping has invalid {column} values at CSV rows: {rows}"
            )
        mapping[column] = values
    return mapping


def _available_key_columns(
    source: pd.DataFrame, requested: Sequence[str]
) -> tuple[str, ...]:
    available = tuple(column for column in requested if column in source.columns)
    if available:
        return available
    return tuple(
        column for column in FALLBACK_IDENTITY_KEY_COLUMNS if column in source.columns
    )


def _reserved_output_columns(name_columns: Sequence[str]) -> set[str]:
    reserved = {
        "identity_reconciliation_version",
        "identity_reconciliation_policy",
        "identity_input_row_number",
        "identity_join_key",
        "identity_join_key_source",
        "identity_join_status",
        "identity_spd_inchikey_columns",
        "identity_spd_inchikey_values",
        "identity_spd_structure_relation",
        "identity_spd_labels_nulled",
        "identity_spd_nulled_columns",
        "identity_original_spd_values_json",
        "identity_verification_state",
        "identity_blocking",
        "identity_blocking_reason_codes",
        "identity_name_reconciled",
        "identity_name_difference_fields",
        "identity_mapping_path",
        "identity_mapping_sha256",
        "identity_row_audit_schema",
        "identity_combined_label_review_required",
        "identity_source_training_allowed",
        "identity_training_allowed",
        "identity_training_exclusion_reason",
    }
    for spec in _EVIDENCE_SPECS:
        reserved.update(
            {
                f"identity_{spec.slot}_pdbqt_status",
                f"identity_{spec.slot}_pdbqt_expected_sha256",
                f"identity_{spec.slot}_pdbqt_declared_sha256",
                f"identity_{spec.slot}_pdbqt_observed_sha256",
                f"identity_{spec.slot}_pdbqt_path_status",
                f"identity_{spec.slot}_pdbqt_evidence_columns",
            }
        )
    reserved.update(_MAPPING_ATTACHMENTS.values())
    reserved.update(f"legacy_{column}" for column in name_columns)
    return reserved


def _select_spd_owned_columns(
    columns: Iterable[str], additional: Sequence[str]
) -> tuple[str, ...]:
    by_casefold = {str(column).casefold(): str(column) for column in columns}
    selected = {
        original
        for normalized, original in by_casefold.items()
        if normalized in SPD_OWNED_VALUE_COLUMNS
    }
    for column in additional:
        if column not in columns:
            raise IdentityReconciliationError(
                f"requested SPD null column is absent: {column}"
            )
        selected.add(column)
    return tuple(column for column in columns if column in selected)


def _protected_columns(columns: Iterable[str]) -> tuple[str, ...]:
    prefixes = (
        "external_",
        "bindingdb_",
        "chembl_",
        "toxcast_",
        "pubchem_",
        "scenario_",
        "combined_activity_",
    )
    protected = []
    for column in columns:
        normalized = str(column).casefold()
        if "score" in normalized or normalized.startswith(prefixes):
            protected.append(str(column))
    return tuple(protected)


def _build_row_contexts(
    source: pd.DataFrame,
    *,
    mapping_by_key: Mapping[str, Mapping[str, object]],
    key_columns: Sequence[str],
    spd_inchikey_columns: Sequence[str],
    spd_owned_columns: Sequence[str],
    source_dir: Path,
) -> list[dict[str, object]]:
    path_hash_cache: dict[Path, tuple[str, str]] = {}
    rows: list[dict[str, object]] = []
    for offset, (_, row) in enumerate(source.iterrows(), start=2):
        key, key_source, join_status = _resolve_row_key(
            row, key_columns, mapping_by_key
        )
        non_fda_probe = _explicit_non_fda_probe(row)
        if non_fda_probe and join_status in {
            "missing_rdk_identity",
            "invalid_rdk_identity",
        }:
            key = _clean(row.get("ligand_base"))
            key_source = "ligand_base"
            join_status = "not_applicable_non_fda_probe"
        mapping_row = mapping_by_key.get(key)
        spd_evidence = _compare_spd_identity(row, mapping_row, spd_inchikey_columns)
        pdbqt_evidence: dict[str, object] = {}
        pdbqt_blocking: list[str] = []
        for spec in _EVIDENCE_SPECS:
            evidence = _verify_pdbqt_evidence(
                row,
                mapping_row,
                spec,
                source_dir=source_dir,
                path_hash_cache=path_hash_cache,
            )
            if (
                non_fda_probe
                and evidence[f"identity_{spec.slot}_pdbqt_status"]
                == "mapping_hash_unavailable"
            ):
                declared = evidence[f"identity_{spec.slot}_pdbqt_declared_sha256"]
                observed = evidence[f"identity_{spec.slot}_pdbqt_observed_sha256"]
                if declared and observed and declared == observed:
                    evidence[f"identity_{spec.slot}_pdbqt_status"] = (
                        "verified_non_fda_file_and_declared_hash"
                    )
                elif observed:
                    evidence[f"identity_{spec.slot}_pdbqt_status"] = (
                        "verified_non_fda_file_hash"
                    )
                elif declared:
                    evidence[f"identity_{spec.slot}_pdbqt_status"] = (
                        "verified_non_fda_declared_hash"
                    )
            pdbqt_evidence.update(evidence)
            status = str(evidence[f"identity_{spec.slot}_pdbqt_status"])
            if status in {
                "invalid_declared_hash",
                "conflicting_declared_hashes",
                "conflicting_file_hashes",
                "declared_hash_path_mismatch",
                "mapping_hash_unavailable",
                "hash_mismatch",
                "path_unreadable",
            }:
                pdbqt_blocking.append(f"{spec.slot}_pdbqt_{status}")

        # Generic activity/PK columns may be populated by external sources.
        # Only explicit SPD-prefixed values establish that this row carries SPD
        # truth, while the broader owned set is still nulled on a proven
        # structure-versus-SPD mismatch.
        spd_values_present = any(
            column in row.index and _clean(row.get(column))
            for column in SPD_EVIDENCE_PRESENCE_COLUMNS
        )
        score_values_present = any(
            "score" in str(column).casefold() and _clean(row.get(column))
            for column in source.columns
        )
        reuse_authorized = _mapping_reuse_authorized(mapping_row)
        blocking_reasons = _blocking_reasons(
            join_status=join_status,
            spd_relation=str(spd_evidence["identity_spd_structure_relation"]),
            spd_values_present=spd_values_present,
            score_values_present=score_values_present,
            reuse_authorized=reuse_authorized,
            pdbqt_reasons=pdbqt_blocking,
        )
        verification_state = _verification_state(
            join_status,
            str(spd_evidence["identity_spd_structure_relation"]),
            pdbqt_evidence,
            spd_values_present,
            score_values_present,
            reuse_authorized,
        )
        rows.append(
            {
                "identity_reconciliation_version": RECONCILIATION_VERSION,
                "identity_reconciliation_policy": (
                    "exact_normalized_rdk_join;no_name_score_join;missing_is_unknown"
                ),
                "identity_input_row_number": offset,
                "identity_join_key": key,
                "identity_join_key_source": key_source,
                "identity_join_status": join_status,
                **spd_evidence,
                **pdbqt_evidence,
                "identity_verification_state": verification_state,
                "identity_blocking": bool(blocking_reasons),
                "identity_blocking_reason_codes": ";".join(blocking_reasons),
                "identity_spd_labels_nulled": False,
                "identity_spd_nulled_columns": "",
                "identity_name_reconciled": False,
                "identity_name_difference_fields": "",
                "identity_combined_label_review_required": False,
            }
        )
    return rows


def _resolve_row_key(
    row: pd.Series,
    key_columns: Sequence[str],
    mapping_by_key: Mapping[str, Mapping[str, object]],
) -> tuple[str, str, str]:
    populated = [(column, _clean(row.get(column))) for column in key_columns]
    populated = [(column, value) for column, value in populated if value]
    parsed = [(column, normalize_rdk_identity(value)) for column, value in populated]
    valid = [(column, value) for column, value in parsed if value]
    values = {value for _, value in valid}
    if not populated:
        return "", "", "missing_rdk_identity"
    if not valid:
        return "", ";".join(column for column, _ in populated), "invalid_rdk_identity"
    if len(values) > 1:
        return (
            "",
            ";".join(column for column, _ in valid),
            "conflicting_rdk_identities",
        )
    key = next(iter(values))
    source_columns = ";".join(column for column, value in valid if value == key)
    if key not in mapping_by_key:
        return key, source_columns, "rdk_identity_not_in_mapping"
    return key, source_columns, "joined_exact_normalized_rdk"


def _normalize_inchikey(value: object) -> str:
    normalized = _clean(value).upper()
    return normalized if _INCHIKEY_RE.fullmatch(normalized) else ""


def _compare_spd_identity(
    row: pd.Series,
    mapping_row: Mapping[str, object] | None,
    columns: Sequence[str],
) -> dict[str, str]:
    populated = [(column, _clean(row.get(column))) for column in columns]
    populated = [(column, value) for column, value in populated if value]
    valid = [(column, _normalize_inchikey(value)) for column, value in populated]
    valid = [(column, value) for column, value in valid if value]
    values = {value for _, value in valid}
    rendered_columns = ";".join(column for column, _ in populated)
    rendered_values = ";".join(sorted(values))
    if mapping_row is None:
        relation = "unresolved_mapping"
    elif not populated:
        relation = "missing_spd_identity"
    elif len(valid) != len(populated):
        relation = "invalid_spd_identity"
    elif len({value[:14] for value in values}) > 1:
        relation = "conflicting_spd_identities"
    elif not _truthy(mapping_row.get("identity_structure_validated")):
        relation = "mapping_identity_unverifiable"
    else:
        expected = {
            _normalize_inchikey(mapping_row.get(column))
            for column in ("identity_exact_inchikey", "identity_parent_inchikey")
        }
        expected.discard("")
        expected_connectivity = {value[:14] for value in expected}
        observed_connectivity = {value[:14] for value in values}
        if not expected_connectivity:
            relation = "mapping_identity_unverifiable"
        elif observed_connectivity.intersection(expected_connectivity):
            relation = (
                "exact_match"
                if bool(values.intersection(expected))
                else "connectivity_match"
            )
        else:
            relation = "demonstrated_structure_spd_mismatch"
    return {
        "identity_spd_inchikey_columns": rendered_columns,
        "identity_spd_inchikey_values": rendered_values,
        "identity_spd_structure_relation": relation,
    }


def _verify_pdbqt_evidence(
    row: pd.Series,
    mapping_row: Mapping[str, object] | None,
    spec: _EvidenceSpec,
    *,
    source_dir: Path,
    path_hash_cache: dict[Path, tuple[str, str]],
) -> dict[str, str]:
    prefix = f"identity_{spec.slot}_pdbqt"
    expected = (
        _clean(mapping_row.get(spec.expected_hash_column)).lower()
        if mapping_row is not None
        else ""
    )
    declared_items = [
        (column, _clean(row.get(column)).lower())
        for column in spec.hash_columns
        if column in row.index and _clean(row.get(column))
    ]
    path_items = [
        (column, _clean(row.get(column)))
        for column in spec.path_columns
        if column in row.index and _clean(row.get(column))
    ]
    evidence_columns = [column for column, _ in declared_items + path_items]
    declared_valid = {
        value for _, value in declared_items if _SHA256_RE.fullmatch(value)
    }
    declared_invalid = any(
        not _SHA256_RE.fullmatch(value) for _, value in declared_items
    )

    observed_hashes: set[str] = set()
    path_states: list[str] = []
    for _, path_text in path_items:
        path = _resolve_evidence_path(path_text, source_dir)
        file_hash, path_state = _hash_evidence_path(path, path_hash_cache)
        path_states.append(path_state)
        if file_hash:
            observed_hashes.add(file_hash)

    path_status = _collapse_path_states(path_states)
    if not evidence_columns:
        status = "not_provided"
    elif declared_invalid:
        status = "invalid_declared_hash"
    elif len(declared_valid) > 1:
        status = "conflicting_declared_hashes"
    elif len(observed_hashes) > 1:
        status = "conflicting_file_hashes"
    elif declared_valid and observed_hashes and declared_valid != observed_hashes:
        status = "declared_hash_path_mismatch"
    elif not expected or not _SHA256_RE.fullmatch(expected):
        status = "mapping_hash_unavailable"
    else:
        available = declared_valid | observed_hashes
        if not available:
            status = "path_unreadable"
        elif available != {expected}:
            status = "hash_mismatch"
        elif declared_valid and observed_hashes:
            status = "verified_file_and_declared_hash"
        elif observed_hashes:
            status = "verified_file_hash"
        else:
            status = "verified_declared_hash"
    return {
        f"{prefix}_status": status,
        f"{prefix}_expected_sha256": expected,
        f"{prefix}_declared_sha256": _single_value(declared_valid),
        f"{prefix}_observed_sha256": _single_value(observed_hashes),
        f"{prefix}_path_status": path_status,
        f"{prefix}_evidence_columns": ";".join(evidence_columns),
    }


def _resolve_evidence_path(value: str, source_dir: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    source_candidate = (source_dir / path).resolve()
    if source_candidate.exists():
        return source_candidate
    repo_root = Path(path_router.__file__).resolve().parents[2]
    repo_candidate = (repo_root / path).resolve()
    return repo_candidate if repo_candidate.exists() else source_candidate


def _hash_evidence_path(
    path: Path, cache: dict[Path, tuple[str, str]]
) -> tuple[str, str]:
    if path in cache:
        return cache[path]
    if not path.exists():
        result = ("", "missing")
    elif not path.is_file():
        result = ("", "not_file")
    else:
        try:
            result = (_sha256_file(path), "verified_file")
        except OSError:
            result = ("", "read_error")
    cache[path] = result
    return result


def _collapse_path_states(states: Sequence[str]) -> str:
    if not states:
        return "not_provided"
    unique = sorted(set(states))
    return unique[0] if len(unique) == 1 else "conflicting:" + ";".join(unique)


def _single_value(values: set[str]) -> str:
    return next(iter(values)) if len(values) == 1 else ""


def _mapping_reuse_authorized(mapping_row: Mapping[str, object] | None) -> bool:
    if mapping_row is None:
        return False
    reuse_status = _clean(mapping_row.get("legacy_score_reuse_status"))
    connectivity = _clean(mapping_row.get("prepared_connectivity_state"))
    return _REUSE_STATUS_BY_CONNECTIVITY.get(reuse_status) == connectivity


def _blocking_reasons(
    *,
    join_status: str,
    spd_relation: str,
    spd_values_present: bool,
    score_values_present: bool,
    reuse_authorized: bool,
    pdbqt_reasons: Sequence[str],
) -> list[str]:
    reasons = list(pdbqt_reasons)
    if join_status == "not_applicable_non_fda_probe":
        if spd_values_present:
            reasons.append("non_fda_probe_has_spd_values")
    elif join_status != "joined_exact_normalized_rdk":
        reasons.append(join_status)
    if spd_relation == "demonstrated_structure_spd_mismatch":
        reasons.append("demonstrated_structure_spd_mismatch")
    elif spd_values_present and spd_relation not in {
        "exact_match",
        "connectivity_match",
    }:
        reasons.append(f"spd_labels_{spd_relation}")
    if (
        score_values_present
        and join_status == "joined_exact_normalized_rdk"
        and not reuse_authorized
    ):
        reasons.append("mapping_score_reuse_not_authorized")
    return list(dict.fromkeys(reasons))


def _verification_state(
    join_status: str,
    spd_relation: str,
    pdbqt_evidence: Mapping[str, object],
    spd_values_present: bool,
    score_values_present: bool,
    reuse_authorized: bool,
) -> str:
    if join_status == "not_applicable_non_fda_probe":
        return (
            "non_fda_probe_quarantined"
            if not spd_values_present
            else "non_fda_probe_has_spd_values"
        )
    if join_status != "joined_exact_normalized_rdk":
        return "identity_unresolved"
    if spd_relation == "demonstrated_structure_spd_mismatch":
        return "structure_vs_spd_identity_mismatch"
    statuses = [
        str(pdbqt_evidence[f"identity_{spec.slot}_pdbqt_status"])
        for spec in _EVIDENCE_SPECS
    ]
    if any(
        status.endswith("mismatch") or status.startswith("conflicting")
        for status in statuses
    ):
        return "pdbqt_identity_mismatch"
    valid_pdbqt_statuses = {
        "not_provided",
        "verified_declared_hash",
        "verified_file_hash",
        "verified_file_and_declared_hash",
    }
    if any(status not in valid_pdbqt_statuses for status in statuses):
        return "pdbqt_identity_unverifiable"
    if spd_relation in {"exact_match", "connectivity_match"}:
        return (
            "verified_structure_and_spd_identity"
            if reuse_authorized or not score_values_present
            else "verified_spd_identity_score_reuse_unverifiable"
        )
    if spd_values_present:
        return "structure_vs_spd_identity_unverifiable"
    if any(status.startswith("verified_") for status in statuses):
        return "verified_pdbqt_identity_no_spd_labels"
    return "canonical_identity_mapped_no_spd_labels"


def _attach_mapping_fields(
    reconciled: pd.DataFrame,
    mapping_by_key: Mapping[str, Mapping[str, object]],
    *,
    mapping_path: Path,
    mapping_sha256: str,
) -> None:
    records = reconciled["identity_join_key"].map(mapping_by_key.get)
    for mapping_column, output_column in _MAPPING_ATTACHMENTS.items():
        reconciled[output_column] = records.map(
            lambda record: _clean(record.get(mapping_column)) if record else ""
        )
    reconciled["identity_mapping_path"] = str(mapping_path)
    reconciled["identity_mapping_sha256"] = mapping_sha256


def _reconcile_name_fields(
    reconciled: pd.DataFrame, name_columns: Sequence[str]
) -> None:
    preferred = reconciled["canonical_identity_name"].fillna("").astype(str)
    display = reconciled["canonical_identity_display_name"].fillna("").astype(str)
    generic = reconciled["canonical_identity_generic_name"].fillna("").astype(str)
    canonical_by_column = {
        "display_name": display.where(display.str.strip().ne(""), preferred),
        "generic_name": generic.where(generic.str.strip().ne(""), preferred),
    }
    differences: list[list[str]] = [[] for _ in range(len(reconciled))]
    positions = {index: position for position, index in enumerate(reconciled.index)}
    joined = reconciled["identity_join_status"].eq("joined_exact_normalized_rdk")
    validated = reconciled["canonical_identity_structure_validated"].map(_truthy)
    for column in name_columns:
        canonical = canonical_by_column.get(column, preferred)
        update = joined & validated & canonical.str.strip().ne("")
        changed = update & (
            reconciled[f"legacy_{column}"].fillna("").astype(str) != canonical
        )
        reconciled.loc[update, column] = canonical.loc[update]
        for index in reconciled.index[changed]:
            differences[positions[index]].append(column)
    rendered = [";".join(columns) for columns in differences]
    reconciled["identity_name_difference_fields"] = rendered
    reconciled["identity_name_reconciled"] = [bool(value) for value in rendered]


def _null_mismatched_spd_values(
    reconciled: pd.DataFrame,
    source: pd.DataFrame,
    spd_owned_columns: Sequence[str],
) -> pd.Series:
    for column in ("spd_label_status", "spd_missing_reason"):
        if column in reconciled.columns:
            reconciled[column] = reconciled[column].astype("object")
    mismatch = reconciled["identity_spd_structure_relation"].eq(
        "demonstrated_structure_spd_mismatch"
    )
    original_json = pd.Series("{}", index=reconciled.index, dtype="object")
    for index in reconciled.index[mismatch]:
        original: dict[str, object] = {}
        nulled: list[str] = []
        for column in spd_owned_columns:
            value = source.at[index, column]
            if pd.isna(value) or _clean(value) == "":
                continue
            original[column] = _json_value(value)
            reconciled.at[index, column] = pd.NA
            nulled.append(column)
        reconciled.at[index, "identity_spd_labels_nulled"] = bool(nulled)
        reconciled.at[index, "identity_spd_nulled_columns"] = ";".join(nulled)
        original_json.at[index] = json.dumps(
            original, sort_keys=True, separators=(",", ":")
        )
        if "spd_label_status" in reconciled.columns:
            reconciled.at[index, "spd_label_status"] = "identity_mismatch"
        if "spd_missing_reason" in reconciled.columns:
            reconciled.at[index, "spd_missing_reason"] = (
                "demonstrated_structure_spd_identity_mismatch"
            )
        combined_columns = [
            column
            for column in reconciled.columns
            if str(column).casefold().startswith("combined_activity_")
            and "label" in str(column).casefold()
            and _clean(reconciled.at[index, column])
        ]
        reconciled.at[index, "identity_combined_label_review_required"] = bool(
            combined_columns
        )
    return original_json


def _apply_identity_training_gate(
    reconciled: pd.DataFrame, source: pd.DataFrame
) -> None:
    """Prevent unresolved identity rows from silently entering model fitting."""

    if "training_allowed" in source.columns:
        source_allowed = source["training_allowed"].map(_truthy)
        reconciled["identity_source_training_allowed"] = source[
            "training_allowed"
        ]
    else:
        source_allowed = pd.Series(True, index=source.index, dtype=bool)
        reconciled["identity_source_training_allowed"] = True
    identity_allowed = ~reconciled["identity_blocking"].fillna(True).astype(bool)
    final_allowed = source_allowed & identity_allowed
    reconciled["identity_training_allowed"] = final_allowed
    if "training_allowed" in reconciled.columns:
        reconciled["training_allowed"] = final_allowed
    reconciled["identity_training_exclusion_reason"] = ""
    source_blocked = ~source_allowed
    reconciled.loc[
        source_blocked, "identity_training_exclusion_reason"
    ] = "source_policy_disallows_training"
    identity_blocked = ~identity_allowed
    reconciled.loc[
        identity_blocked, "identity_training_exclusion_reason"
    ] = reconciled.loc[identity_blocked, "identity_blocking_reason_codes"]


def _json_value(value: object) -> object:
    if value is None or pd.isna(value):
        return None
    if hasattr(value, "item"):
        try:
            return value.item()  # type: ignore[no-any-return, union-attr]
        except (AttributeError, ValueError):
            pass
    return value if isinstance(value, (str, int, float, bool)) else str(value)


def _assert_protected_columns_unchanged(
    reconciled: pd.DataFrame,
    snapshot: pd.DataFrame,
    columns: Sequence[str],
) -> None:
    changed = [
        column for column in columns if not reconciled[column].equals(snapshot[column])
    ]
    if changed:
        raise IdentityReconciliationError(
            "score or external evidence columns changed unexpectedly: "
            + ", ".join(changed)
        )


def _build_row_audit(
    reconciled: pd.DataFrame,
    source_columns: Iterable[str],
    name_columns: Sequence[str],
) -> pd.DataFrame:
    source_identifiers = [
        column
        for column in (
            "ligand_base",
            "rdk_id",
            "_join_rdk",
            "target_id",
            "target_gene",
            "target_uniprot",
            "pdb_id",
            "spd_drug_id",
            "spd_inchikey",
        )
        if column in source_columns
    ]
    legacy_names = [f"legacy_{column}" for column in name_columns]
    identity_columns = [
        column
        for column in reconciled.columns
        if column.startswith("identity_") or column.startswith("canonical_identity_")
    ]
    columns = list(
        dict.fromkeys([*source_identifiers, *legacy_names, *identity_columns])
    )
    audit = reconciled[columns].copy()
    audit.insert(0, "identity_row_audit_schema", AUDIT_SCHEMA)
    return audit


def _value_counts(series: pd.Series) -> dict[str, int]:
    counts = series.fillna("<missing>").astype(str).value_counts(dropna=False)
    return {str(key): int(value) for key, value in counts.sort_index().items()}


def _build_summary(
    *,
    source_path: Path,
    input_sha256: str,
    source: pd.DataFrame,
    map_path: Path,
    mapping_sha256: str,
    mapping: pd.DataFrame,
    outputs: ReconciliationOutputs,
    reconciled: pd.DataFrame,
    key_columns: Sequence[str],
    spd_inchikey_columns: Sequence[str],
    spd_owned_columns: Sequence[str],
    name_columns: Sequence[str],
    protected_columns: Sequence[str],
    fail_closed: bool,
) -> dict[str, Any]:
    blocking_rows = int(reconciled["identity_blocking"].fillna(False).sum())
    mismatch = reconciled["identity_spd_structure_relation"].eq(
        "demonstrated_structure_spd_mismatch"
    )
    nulled_rows = reconciled["identity_spd_labels_nulled"].fillna(False).astype(bool)
    nulled_cells = int(
        reconciled.loc[nulled_rows, "identity_spd_nulled_columns"]
        .fillna("")
        .astype(str)
        .map(lambda value: len([item for item in value.split(";") if item]))
        .sum()
    )
    summary: dict[str, Any] = {
        "schema_version": SUMMARY_SCHEMA,
        "reconciliation_version": RECONCILIATION_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "release_status": "blocked" if blocking_rows else "ready",
        "fail_closed_enforced": bool(fail_closed),
        "recommended_exit_code": 2 if fail_closed and blocking_rows else 0,
        "input": {
            "path": str(source_path),
            "sha256": input_sha256,
            "rows": int(len(source)),
            "columns": int(len(source.columns)),
        },
        "mapping": {
            "path": str(map_path),
            "sha256": mapping_sha256,
            "rows": int(len(mapping)),
            "key": "rdk_id normalized to rdk_NNNNNNN",
        },
        "outputs": {
            "reconciled_csv": {"path": str(outputs.reconciled_csv)},
            "row_audit_csv": {"path": str(outputs.row_audit_csv)},
            "summary_json": {"path": str(outputs.summary_json)},
        },
        "policy": {
            "identity_join": "exact normalized RDK/ligand_base only; names are never join keys",
            "score_policy": "rows and score columns are preserved; scores are never reassigned by name",
            "missing_policy": "missing identity or label evidence remains unknown, never negative",
            "mismatch_policy": (
                "only a validated canonical-vs-SPD InChIKey connectivity mismatch nulls "
                "SPD-owned binding/exposure/AC50/free-Cmax fields"
            ),
            "external_label_policy": "external and combined label columns are preserved unchanged",
            "training_policy": (
                "training_allowed is fail-closed by source policy and identity_blocking; "
                "the original source value is preserved in identity_source_training_allowed"
            ),
            "historical_hash_policy": (
                "historical score/PDBQT hashes are checked against the mapping selected_pdbqt_sha256"
            ),
        },
        "columns": {
            "identity_key_columns": list(key_columns),
            "spd_inchikey_columns": list(spd_inchikey_columns),
            "legacy_name_columns": [f"legacy_{column}" for column in name_columns],
            "spd_owned_null_columns": list(spd_owned_columns),
            "protected_score_external_columns": list(protected_columns),
        },
        "counts": {
            "rows": int(len(reconciled)),
            "blocking_rows": blocking_rows,
            "demonstrated_structure_spd_mismatch_rows": int(mismatch.sum()),
            "spd_values_nulled_rows": int(nulled_rows.sum()),
            "spd_values_nulled_cells": nulled_cells,
            "name_reconciled_rows": int(
                reconciled["identity_name_reconciled"].fillna(False).sum()
            ),
            "combined_label_review_required_rows": int(
                reconciled["identity_combined_label_review_required"]
                .fillna(False)
                .sum()
            ),
            "identity_training_allowed_rows": int(
                reconciled["identity_training_allowed"].fillna(False).sum()
            ),
            "identity_training_excluded_rows": int(
                (~reconciled["identity_training_allowed"].fillna(False)).sum()
            ),
            "join_status": _value_counts(reconciled["identity_join_status"]),
            "verification_state": _value_counts(
                reconciled["identity_verification_state"]
            ),
            "spd_structure_relation": _value_counts(
                reconciled["identity_spd_structure_relation"]
            ),
            "blocking_reason_codes": _reason_counts(
                reconciled["identity_blocking_reason_codes"]
            ),
            "pdbqt_status": {
                spec.slot: _value_counts(
                    reconciled[f"identity_{spec.slot}_pdbqt_status"]
                )
                for spec in _EVIDENCE_SPECS
            },
        },
    }
    return summary


def _reason_counts(series: pd.Series) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in series.fillna("").astype(str):
        for reason in (item for item in value.split(";") if item):
            counts[reason] = counts.get(reason, 0) + 1
    return dict(sorted(counts.items()))


def _assert_inputs_stable(
    source_path: Path,
    input_sha256: str,
    mapping_path: Path,
    mapping_sha256: str,
) -> None:
    if _sha256_file(source_path) != input_sha256:
        raise IdentityReconciliationError(
            f"input changed during reconciliation: {source_path}"
        )
    if _sha256_file(mapping_path) != mapping_sha256:
        raise IdentityReconciliationError(
            f"mapping changed during reconciliation: {mapping_path}"
        )


def _atomic_to_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        frame.to_csv(temporary, index=False)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_write_json(payload: Mapping[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
