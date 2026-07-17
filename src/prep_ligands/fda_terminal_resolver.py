"""Terminal, structure-first v3 disposition for every legacy FDA-map row."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from rdkit import Chem

from prep_ligands.fda_identity_audit import (
    _structure_from_mol,
    _structure_from_smiles,
)


SCHEMA_VERSION = 3
RESOLVER_VERSION = "fda-terminal-resolution-v3"

TERMINAL_DISPOSITIONS = frozenset(
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

_EXTRA_MAPPING_FIELDS = [
    "terminal_resolution_version",
    "v2_display_name",
    "preferred_identity",
    "preferred_identity_source",
    "preferred_identity_source_id",
    "identity_structure_relation",
    "identity_structure_validated",
    "identity_resolution_status",
    "resolved_preferred_name",
    "resolved_drugcentral_id",
    "identity_exact_inchikey",
    "identity_parent_inchikey",
    "regulatory_status",
    "regulatory_evidence",
    "identity_evidence_codes",
    "approval_evidence_codes",
    "terminal_disposition",
    "terminal_reason_codes",
    "canonical_parent_id",
    "canonical_materialization_status",
    "canonical_pdbqt_path",
    "canonical_pdbqt_sha256",
    "selected_pdbqt_checksum_status_v3",
    "prepared_connectivity_state",
    "legacy_score_reuse_status",
]

_DISPOSITION_FIELDS = [
    "mapping_row_number",
    "rdk_id",
    "v2_display_name",
    "preferred_identity",
    "preferred_identity_source",
    "preferred_identity_source_id",
    "identity_structure_relation",
    "identity_structure_validated",
    "identity_resolution_status",
    "mapping_exact_inchikey",
    "mapping_parent_inchikey",
    "identity_exact_inchikey",
    "identity_parent_inchikey",
    "regulatory_status",
    "regulatory_evidence",
    "identity_evidence_codes",
    "approval_evidence_codes",
    "terminal_disposition",
    "terminal_reason_codes",
    "legacy_file_usable",
    "selected_pdbqt_path",
    "selected_pdbqt_sha256",
    "selected_pdbqt_checksum_status_v3",
    "prepared_connectivity_state",
    "legacy_score_reuse_status",
    "canonical_parent_id",
    "canonical_materialization_status",
    "canonical_pdbqt_path",
    "canonical_pdbqt_sha256",
]

_COMMON_FORM_SUFFIXES = (
    "monohydrate",
    "dihydrate",
    "trihydrate",
    "hydrate",
    "hydrochloride",
    "hydrobromide",
    "mesylate",
    "besylate",
    "maleate",
    "fumarate",
    "tartrate",
    "succinate",
    "phosphate",
    "sulfate",
    "acetate",
    "sodium",
    "potassium",
    "calcium",
)

_FORMULA_LIKE_NAME_RE = re.compile(r"^(?:[A-Z][a-z]?\d*){2,}$")
_CAS_NAME_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")


@dataclass(frozen=True)
class FDATerminalV3Outputs:
    mapping_v3_csv: Path
    terminal_dispositions_csv: Path
    fda_subset_csv: Path
    excluded_subset_csv: Path
    summary_json: Path


@dataclass(frozen=True)
class _DCRecord:
    index: int
    drugcentral_id: str
    preferred_name: str
    names: tuple[str, ...]
    exact_inchikey: str
    parent_inchikey: str
    connectivity_key: str
    approved_snapshot: bool


@dataclass(frozen=True)
class _DCIndex:
    records: tuple[_DCRecord, ...]
    by_exact: Mapping[str, tuple[_DCRecord, ...]]
    by_parent: Mapping[str, tuple[_DCRecord, ...]]
    by_connectivity: Mapping[str, tuple[_DCRecord, ...]]


@dataclass(frozen=True)
class _PubChemRecord:
    cid: str
    title: str
    iupac_name: str
    exact_inchikey: str
    parent_inchikey: str
    connectivity_key: str
    source_content_sha256: str


@dataclass(frozen=True)
class _ExactKeyCache:
    by_exact: Mapping[str, _PubChemRecord]
    outcomes: Mapping[str, str]
    query_by_rdk_id: Mapping[str, str]


@dataclass(frozen=True)
class _UnsafeParentQuarantine:
    by_parent: Mapping[str, Mapping[str, str]]
    source_indices: frozenset[int]
    row_count: int


@dataclass(frozen=True)
class _IdentityResolution:
    preferred_name: str
    source: str
    source_id: str
    relation: str
    exact_inchikey: str
    parent_inchikey: str
    structure_validated: bool
    regulatory_status: str
    regulatory_evidence: str
    reason_codes: tuple[str, ...]
    canonical_parent_key: str


def terminal_v3_paths(output_dir: Path) -> FDATerminalV3Outputs:
    return FDATerminalV3Outputs(
        mapping_v3_csv=output_dir / "fda_legacy_mapping_v3.csv",
        terminal_dispositions_csv=output_dir / "fda_terminal_dispositions_v3.csv",
        fda_subset_csv=output_dir / "fda_verified_fda_subset_v3.csv",
        excluded_subset_csv=output_dir / "fda_excluded_or_nonfda_subset_v3.csv",
        summary_json=output_dir / "fda_terminal_v3_summary.json",
    )


def resolve_fda_terminal_v3(
    *,
    mapping_v2_csv: Path,
    unresolved_v2_csv: Path,
    repaired_mapping_csv: Path,
    full_drugcentral_sdf: Path,
    approved_source_sdf: Path,
    fda_approved_csv: Path,
    pubchem_records_jsonl: Path,
    pubchem_manifest_json: Path,
    pubchem_inchikey_records_jsonl: Path,
    pubchem_inchikey_manifest_json: Path,
    drugsfda_products_txt: Path,
    named_library_manifest_csv: Path,
    repair_quarantine_csv: Path,
    output_dir: Path,
) -> tuple[FDATerminalV3Outputs, dict[str, Any]]:
    """Resolve preferred identities and emit one terminal disposition per row."""

    inputs = _resolved_inputs(
        mapping_v2_csv=mapping_v2_csv,
        unresolved_v2_csv=unresolved_v2_csv,
        repaired_mapping_csv=repaired_mapping_csv,
        full_drugcentral_sdf=full_drugcentral_sdf,
        approved_source_sdf=approved_source_sdf,
        fda_approved_csv=fda_approved_csv,
        pubchem_records_jsonl=pubchem_records_jsonl,
        pubchem_manifest_json=pubchem_manifest_json,
        pubchem_inchikey_records_jsonl=pubchem_inchikey_records_jsonl,
        pubchem_inchikey_manifest_json=pubchem_inchikey_manifest_json,
        drugsfda_products_txt=drugsfda_products_txt,
        named_library_manifest_csv=named_library_manifest_csv,
        repair_quarantine_csv=repair_quarantine_csv,
    )
    named_summary = inputs["named_library_manifest_csv"].with_name(
        "fda_named_library_summary.json"
    )
    if not named_summary.is_file():
        raise FileNotFoundError(
            f"canonical named-library summary is not a file: {named_summary}"
        )
    inputs["named_library_summary_json"] = named_summary
    output_dir = Path(output_dir).expanduser().resolve()
    outputs = terminal_v3_paths(output_dir)
    _validate_output_scope(outputs, inputs.values())

    mapping_fields, mapping_rows = _read_csv(inputs["mapping_v2_csv"])
    _, unresolved_rows = _read_csv(inputs["unresolved_v2_csv"])
    unresolved = _validate_mapping_partition(mapping_rows, unresolved_rows)
    _, repaired_rows = _read_csv(inputs["repaired_mapping_csv"])
    repaired = _validate_repaired_mapping(mapping_rows, repaired_rows)
    approved_ids, approved_names = _read_approved_snapshot(
        inputs["fda_approved_csv"]
    )
    official_ingredients = _read_official_ingredients(
        inputs["drugsfda_products_txt"]
    )
    dc_index = _read_full_drugcentral(
        inputs["full_drugcentral_sdf"], approved_ids
    )
    pubchem = _read_pubchem_cache(
        inputs["pubchem_records_jsonl"], inputs["pubchem_manifest_json"]
    )
    exact_key_cache = _read_exact_key_cache(
        inputs["pubchem_inchikey_records_jsonl"],
        inputs["pubchem_inchikey_manifest_json"],
    )
    canonical = _read_and_validate_named_library(
        inputs["named_library_manifest_csv"],
        inputs["named_library_summary_json"],
    )
    unsafe_parents = _read_unsafe_parent_quarantine(
        inputs["repair_quarantine_csv"]
    )
    approved_source_count = _validate_approved_source_partition(
        inputs["approved_source_sdf"], canonical, unsafe_parents
    )

    v3_rows: list[dict[str, Any]] = []
    dispositions: list[dict[str, Any]] = []
    for row in mapping_rows:
        disposition = _resolve_terminal_row(
            row=row,
            unresolved=unresolved.get(row["rdk_id"]),
            repaired_row=repaired[row["rdk_id"]],
            dc_index=dc_index,
            pubchem=pubchem,
            exact_key_cache=exact_key_cache,
            approved_names=approved_names,
            official_ingredients=official_ingredients,
            canonical=canonical,
            unsafe_parents=unsafe_parents.by_parent,
        )
        dispositions.append(disposition)
        v3_rows.append(_mapping_v3_row(row, disposition, mapping_fields))

    _validate_terminal_rows(mapping_rows, dispositions)
    _validate_regression_rows(mapping_rows, dispositions)
    fda_ids = {
        row["rdk_id"]
        for row in dispositions
        if row["regulatory_status"]
        in {"drugcentral_fda_approved", "fda_approved_current_or_historical"}
        and row["canonical_materialization_status"] == "ready"
    }
    fda_rows = [row for row in v3_rows if row["rdk_id"] in fda_ids]
    excluded_rows = [row for row in v3_rows if row["rdk_id"] not in fda_ids]

    output_dir.mkdir(parents=True, exist_ok=True)
    v3_fields = [*mapping_fields, *_new_fields(mapping_fields)]
    _write_csv(outputs.mapping_v3_csv, v3_rows, v3_fields)
    _write_csv(
        outputs.terminal_dispositions_csv, dispositions, _DISPOSITION_FIELDS
    )
    _write_csv(outputs.fda_subset_csv, fda_rows, v3_fields)
    _write_csv(outputs.excluded_subset_csv, excluded_rows, v3_fields)
    summary = _summary(
        inputs=inputs,
        outputs=outputs,
        mapping_rows=mapping_rows,
        unresolved_rows=unresolved_rows,
        dispositions=dispositions,
        dc_index=dc_index,
        pubchem=pubchem,
        exact_key_cache=exact_key_cache,
        canonical=canonical,
        unsafe_parents=unsafe_parents,
        approved_source_count=approved_source_count,
    )
    _write_json(outputs.summary_json, summary)
    return outputs, summary


def _resolve_terminal_row(
    *,
    row: Mapping[str, str],
    unresolved: Mapping[str, str] | None,
    repaired_row: Mapping[str, str],
    dc_index: _DCIndex,
    pubchem: Mapping[str, _PubChemRecord],
    exact_key_cache: _ExactKeyCache,
    approved_names: set[str],
    official_ingredients: set[str],
    canonical: Mapping[str, Mapping[str, str]],
    unsafe_parents: Mapping[str, Mapping[str, str]],
) -> dict[str, Any]:
    mapping_exact, mapping_parent = _mapping_keys(row, unresolved)
    checksum_status, usable = _verify_selected_pdbqt(row)
    connectivity_state = repaired_row["prepared_connectivity_state"]
    identity = _identity_for_row(
        row=row,
        mapping_exact=mapping_exact,
        mapping_parent=mapping_parent,
        dc_index=dc_index,
        pubchem=pubchem,
        exact_key_cache=exact_key_cache,
        approved_names=approved_names,
        official_ingredients=official_ingredients,
    )
    canonical_row = canonical.get(
        identity.canonical_parent_key,
        unsafe_parents.get(identity.canonical_parent_key, {}),
    )
    disposition, reasons = _terminal_disposition(
        row=row,
        unresolved=unresolved,
        identity=identity,
        usable=usable,
        canonical_row=canonical_row,
    )
    return {
        "mapping_row_number": row["mapping_row_number"],
        "rdk_id": row["rdk_id"],
        "v2_display_name": row.get("display_name", ""),
        "preferred_identity": identity.preferred_name,
        "preferred_identity_source": identity.source,
        "preferred_identity_source_id": identity.source_id,
        "identity_structure_relation": identity.relation,
        "identity_structure_validated": _bool(identity.structure_validated),
        "identity_resolution_status": (
            "terminal_structure_validated"
            if identity.structure_validated and identity.preferred_name
            else "terminal_unresolved"
        ),
        "mapping_exact_inchikey": mapping_exact,
        "mapping_parent_inchikey": mapping_parent,
        "identity_exact_inchikey": identity.exact_inchikey,
        "identity_parent_inchikey": identity.parent_inchikey,
        "regulatory_status": identity.regulatory_status,
        "regulatory_evidence": identity.regulatory_evidence,
        "identity_evidence_codes": _identity_evidence_codes(identity),
        "approval_evidence_codes": _approval_evidence_codes(identity),
        "terminal_disposition": disposition,
        "terminal_reason_codes": _joined([*identity.reason_codes, *reasons]),
        "legacy_file_usable": _bool(usable),
        "selected_pdbqt_path": row.get("selected_pdbqt_path", ""),
        "selected_pdbqt_sha256": row.get("selected_pdbqt_sha256", ""),
        "selected_pdbqt_checksum_status_v3": checksum_status,
        "prepared_connectivity_state": connectivity_state,
        "legacy_score_reuse_status": _score_reuse_status(
            connectivity_state, canonical_row
        ),
        "canonical_parent_id": canonical_row.get("canonical_parent_id", ""),
        "canonical_materialization_status": canonical_row.get(
            "materialization_status", ""
        ),
        "canonical_pdbqt_path": canonical_row.get("output_pdbqt_path", ""),
        "canonical_pdbqt_sha256": canonical_row.get("output_pdbqt_sha256", ""),
    }


def _identity_for_row(
    *,
    row: Mapping[str, str],
    mapping_exact: str,
    mapping_parent: str,
    dc_index: _DCIndex,
    pubchem: Mapping[str, _PubChemRecord],
    exact_key_cache: _ExactKeyCache,
    approved_names: set[str],
    official_ingredients: set[str],
) -> _IdentityResolution:
    existing = _existing_manifest_identity(row)
    if existing is not None:
        return existing
    drugcentral = _drugcentral_identity_for_keys(
        mapping_exact, mapping_parent, dc_index
    )
    if drugcentral is not None:
        return drugcentral
    return _pubchem_identity_for_row(
        row=row,
        mapping_exact=mapping_exact,
        mapping_parent=mapping_parent,
        dc_index=dc_index,
        pubchem=pubchem,
        exact_key_cache=exact_key_cache,
        approved_names=approved_names,
        official_ingredients=official_ingredients,
    )


def _drugcentral_identity_for_keys(
    mapping_exact: str,
    mapping_parent: str,
    dc_index: _DCIndex,
) -> _IdentityResolution | None:
    exact = _unique_record(dc_index.by_exact.get(mapping_exact, ()))
    if exact is not None:
        return _drugcentral_identity(exact, "exact")
    parent = _unique_record(dc_index.by_parent.get(mapping_parent, ()))
    return _drugcentral_identity(parent, "parent") if parent is not None else None


def _pubchem_identity_for_row(
    *,
    row: Mapping[str, str],
    mapping_exact: str,
    mapping_parent: str,
    dc_index: _DCIndex,
    pubchem: Mapping[str, _PubChemRecord],
    exact_key_cache: _ExactKeyCache,
    approved_names: set[str],
    official_ingredients: set[str],
) -> _IdentityResolution:
    cid = _clean(row.get("pubchem_cid_resolved"))
    pubchem_record = pubchem.get(cid)
    pubchem_relation = _structure_relation(
        mapping_exact,
        mapping_parent,
        pubchem_record.exact_inchikey if pubchem_record else "",
        pubchem_record.parent_inchikey if pubchem_record else "",
    )
    corroborated_dc = _corroborated_connectivity_record(
        row,
        mapping_parent,
        dc_index,
        approved_names,
        official_ingredients,
    )
    if pubchem_record and pubchem_relation in {"exact", "parent"}:
        return _pubchem_identity(
            pubchem_record, pubchem_relation, corroborated_dc
        )
    if (
        pubchem_record
        and pubchem_relation == "connectivity"
        and corroborated_dc is not None
    ):
        return _connectivity_corroborated_identity(
            pubchem_record, corroborated_dc
        )
    exact_query = mapping_exact or exact_key_cache.query_by_rdk_id.get(
        _clean(row.get("rdk_id")), ""
    )
    exact_record = exact_key_cache.by_exact.get(exact_query)
    if exact_record is not None:
        return _pubchem_identity(exact_record, "exact", None)
    return _structure_only_identity(
        row=row,
        mapping_exact=exact_query,
        mapping_parent=mapping_parent,
        pubchem_relation=pubchem_relation,
    )


def _existing_manifest_identity(
    row: Mapping[str, str],
) -> _IdentityResolution | None:
    source_id = _clean(row.get("resolved_drugcentral_id"))
    preferred = _clean(row.get("resolved_preferred_name"))
    exact = _clean(row.get("approved_full_form_exact_inchikey"))
    parent = _clean(row.get("approved_parent_inchikey"))
    if not (source_id and preferred and exact and parent):
        return None
    return _IdentityResolution(
        preferred_name=preferred,
        source="manifest_approved_v2",
        source_id=f"drugcentral:{source_id}",
        relation="approved_exact_or_parent",
        exact_inchikey=exact,
        parent_inchikey=parent,
        structure_validated=True,
        regulatory_status="drugcentral_fda_approved",
        regulatory_evidence="validated_drugcentral_fda_manifest_v2",
        reason_codes=("existing_manifest_approved_identity",),
        canonical_parent_key=parent,
    )


def _drugcentral_identity(record: _DCRecord, relation: str) -> _IdentityResolution:
    approved = record.approved_snapshot
    return _IdentityResolution(
        preferred_name=record.preferred_name,
        source="full_drugcentral_structure",
        source_id=f"drugcentral:{record.drugcentral_id}",
        relation=relation,
        exact_inchikey=record.exact_inchikey,
        parent_inchikey=record.parent_inchikey,
        structure_validated=True,
        regulatory_status=(
            "drugcentral_fda_approved"
            if approved
            else "not_fda_approved_drugcentral_snapshot"
        ),
        regulatory_evidence=(
            "drugcentral_id_in_fda_approved_snapshot"
            if approved
            else "drugcentral_id_absent_from_fda_approved_snapshot"
        ),
        reason_codes=(f"unique_drugcentral_{relation}_structure_match",),
        canonical_parent_key=record.parent_inchikey if approved else "",
    )


def _pubchem_identity(
    record: _PubChemRecord,
    relation: str,
    corroborated_dc: _DCRecord | None,
) -> _IdentityResolution:
    if corroborated_dc is not None:
        return _IdentityResolution(
            preferred_name=corroborated_dc.preferred_name,
            source="pubchem_structure_plus_unique_fda_corroboration",
            source_id=(
                f"pubchem:{record.cid};drugcentral:{corroborated_dc.drugcentral_id}"
            ),
            relation=relation,
            exact_inchikey=record.exact_inchikey,
            parent_inchikey=record.parent_inchikey,
            structure_validated=True,
            regulatory_status="fda_approved_current_or_historical",
            regulatory_evidence="unique_approved_dc_name_and_drugsfda_ingredient",
            reason_codes=(
                f"pubchem_{relation}_structure_validated",
                "approved_connectivity_name_officially_corroborated",
            ),
            canonical_parent_key=corroborated_dc.parent_inchikey,
        )
    preferred = _preferred_pubchem_name(record)
    return _IdentityResolution(
        preferred_name=preferred,
        source="pubchem_structure",
        source_id=f"pubchem:{record.cid}",
        relation=relation,
        exact_inchikey=record.exact_inchikey,
        parent_inchikey=record.parent_inchikey,
        structure_validated=bool(preferred),
        regulatory_status="fda_not_verified",
        regulatory_evidence="pubchem_identity_is_not_fda_approval_evidence",
        reason_codes=(f"pubchem_{relation}_structure_validated",),
        canonical_parent_key="",
    )


def _preferred_pubchem_name(record: _PubChemRecord) -> str:
    title = record.title.strip()
    if title and not (
        _FORMULA_LIKE_NAME_RE.fullmatch(title)
        or _CAS_NAME_RE.fullmatch(title)
    ):
        return title
    return record.iupac_name or title


def _connectivity_corroborated_identity(
    record: _PubChemRecord, dc_record: _DCRecord
) -> _IdentityResolution:
    return _IdentityResolution(
        preferred_name=dc_record.preferred_name,
        source="pubchem_connectivity_plus_unique_fda_corroboration",
        source_id=f"pubchem:{record.cid};drugcentral:{dc_record.drugcentral_id}",
        relation="connectivity",
        exact_inchikey=record.exact_inchikey,
        parent_inchikey=record.parent_inchikey,
        structure_validated=True,
        regulatory_status="fda_approved_current_or_historical",
        regulatory_evidence="unique_approved_dc_name_and_drugsfda_ingredient",
        reason_codes=(
            "pubchem_connectivity_validated",
            "approved_connectivity_name_officially_corroborated",
        ),
        canonical_parent_key=dc_record.parent_inchikey,
    )


def _blocking_identity(
    pubchem_record: _PubChemRecord | None, relation: str
) -> _IdentityResolution:
    reasons = ["usable_row_lacks_structure_validated_preferred_identity"]
    if pubchem_record is None:
        reasons.append("pubchem_record_missing")
    else:
        reasons.append(f"pubchem_structure_{relation or 'unresolved'}")
    return _IdentityResolution(
        preferred_name="",
        source="unresolved",
        source_id=(f"pubchem:{pubchem_record.cid}" if pubchem_record else ""),
        relation=relation or "unresolved",
        exact_inchikey=(pubchem_record.exact_inchikey if pubchem_record else ""),
        parent_inchikey=(pubchem_record.parent_inchikey if pubchem_record else ""),
        structure_validated=False,
        regulatory_status="fda_not_verified",
        regulatory_evidence="insufficient_structure_validated_identity",
        reason_codes=tuple(reasons),
        canonical_parent_key="",
    )


def _structure_only_identity(
    *,
    row: Mapping[str, str],
    mapping_exact: str,
    mapping_parent: str,
    pubchem_relation: str,
) -> _IdentityResolution:
    if not mapping_exact or not mapping_parent:
        return _blocking_identity(None, pubchem_relation)
    if pubchem_relation == "connectivity":
        legacy_name = _clean(row.get("display_name")) or "legacy ligand"
        preferred = f"{legacy_name} stereoisomer [{mapping_exact}]"
        source = "mapping_exact_stereoisomer_pdbqt_lineage"
    else:
        preferred = f"structure-only identity [{mapping_exact}]"
        source = "mapping_exact_structure_pdbqt_lineage"
    return _IdentityResolution(
        preferred_name=preferred,
        source=source,
        source_id=f"inchikey:{mapping_exact}",
        relation="exact",
        exact_inchikey=mapping_exact,
        parent_inchikey=mapping_parent,
        structure_validated=True,
        regulatory_status="not_fda_approved_structure_only",
        regulatory_evidence="structure_lineage_is_not_fda_approval_evidence",
        reason_codes=(source,),
        canonical_parent_key="",
    )


def _corroborated_connectivity_record(
    row: Mapping[str, str],
    mapping_parent: str,
    dc_index: _DCIndex,
    approved_names: set[str],
    official_ingredients: set[str],
) -> _DCRecord | None:
    connectivity = _connectivity(mapping_parent)
    candidate = _unique_record(dc_index.by_connectivity.get(connectivity, ()))
    if candidate is None or not candidate.approved_snapshot:
        return None
    row_names = _row_name_keys(row)
    dc_names = {_name_key(name) for name in candidate.names if _name_key(name)}
    if not row_names.intersection(dc_names):
        return None
    if _name_key(candidate.preferred_name) not in approved_names:
        return None
    if not any(_official_name_match(name, official_ingredients) for name in candidate.names):
        return None
    return candidate


def _terminal_disposition(
    *,
    row: Mapping[str, str],
    unresolved: Mapping[str, str] | None,
    identity: _IdentityResolution,
    usable: bool,
    canonical_row: Mapping[str, str],
) -> tuple[str, tuple[str, ...]]:
    early = _early_terminal_disposition(row, unresolved, identity, usable)
    if early is not None:
        return early
    if identity.regulatory_status in {
        "drugcentral_fda_approved",
        "fda_approved_current_or_historical",
    }:
        return _fda_disposition(canonical_row)
    evidence = (
        "positive_full_drugcentral_nonapproval_evidence"
        if identity.regulatory_status == "not_fda_approved_drugcentral_snapshot"
        else "structure_identity_resolved_without_fda_evidence"
    )
    return "identified_non_fda", (evidence,)


def _early_terminal_disposition(
    row: Mapping[str, str],
    unresolved: Mapping[str, str] | None,
    identity: _IdentityResolution,
    usable: bool,
) -> tuple[str, tuple[str, ...]] | None:
    if not usable:
        return "prepared_file_unusable", ("legacy_pdbqt_unusable",)
    category = _clean((unresolved or {}).get("audit_category"))
    curated = _curated_disposition(category)
    if curated is not None:
        return curated
    if not identity.structure_validated or not identity.preferred_name:
        return "incorrect_name_structure_mapping_resolved", (
            "completion_blocker_usable_identity_unresolved",
        )
    structure_only = identity.source.startswith("mapping_exact_")
    if category == "incorrect name–structure mapping" or (
        not structure_only and _name_disagrees(row, identity.preferred_name)
    ):
        return "incorrect_name_structure_mapping_resolved", (
            "legacy_name_corrected_to_structure_validated_identity",
        )
    return None


def _curated_disposition(
    category: str,
) -> tuple[str, tuple[str, ...]] | None:
    if category == "additive/excipient":
        return "additive_excipient", ("curated_nonmedication_disposition",)
    if category == "collapsed combination product":
        return "collapsed_combination_product", (
            "combination_identity_not_single_approved_parent",
        )
    category_key = _name_key(category)
    if "salt" in category_key or "counterion" in category_key:
        return "retained_salt_counterion", (
            "curated_salt_or_counterion_disposition",
        )
    return None


def _fda_disposition(
    canonical_row: Mapping[str, str],
) -> tuple[str, tuple[str, ...]]:
    if not canonical_row:
        return "incorrect_name_structure_mapping_resolved", (
            "completion_blocker_fda_parent_not_in_named_library",
        )
    if canonical_row.get("materialization_status") == "quarantined_unsafe_parent_collapse":
        return "incorrect_name_structure_mapping_resolved", (
            "repair_quarantine_source_record_linked",
            "unsafe_parent_collapse",
        )
    source_classes = canonical_row.get("approved_source_classes", "")
    form_relation = canonical_row.get("approved_form_relation", "")
    if "salt_or_solvate" in source_classes or form_relation != "full_equals_parent":
        return "confirmed_fda_salt_parent_docked", (
            "canonical_fda_salt_parent_linked",
        )
    return "confirmed_fda_active_ingredient", (
        "canonical_fda_parent_linked",
    )


def _mapping_v3_row(
    row: Mapping[str, str],
    disposition: Mapping[str, Any],
    mapping_fields: Sequence[str],
) -> dict[str, Any]:
    output = {field: row.get(field, "") for field in mapping_fields}
    preferred = str(disposition["preferred_identity"])
    output["v2_display_name"] = row.get("display_name", "")
    if preferred:
        output["display_name"] = preferred
        output["generic_name"] = preferred
    output.update(
        {
            field: disposition.get(field, "")
            for field in _EXTRA_MAPPING_FIELDS
            if field != "v2_display_name"
        }
    )
    output["terminal_resolution_version"] = RESOLVER_VERSION
    output["identity_resolution_status"] = disposition.get(
        "identity_resolution_status", ""
    )
    output["resolved_preferred_name"] = preferred
    drugcentral_match = re.search(
        r"(?:^|;)drugcentral:([^;]+)",
        str(disposition.get("preferred_identity_source_id", "")),
    )
    output["resolved_drugcentral_id"] = (
        drugcentral_match.group(1)
        if drugcentral_match
        else row.get("resolved_drugcentral_id", "")
    )
    return output


def _mapping_keys(
    row: Mapping[str, str], unresolved: Mapping[str, str] | None
) -> tuple[str, str]:
    if unresolved is not None:
        return (
            _clean(unresolved.get("legacy_exact_inchikey")).upper(),
            _clean(unresolved.get("legacy_parent_inchikey")).upper(),
        )
    return (
        _clean(row.get("approved_full_form_exact_inchikey") or row.get("inchikey")).upper(),
        _clean(row.get("approved_parent_inchikey")).upper(),
    )


def _structure_relation(
    mapping_exact: str,
    mapping_parent: str,
    candidate_exact: str,
    candidate_parent: str,
) -> str:
    if mapping_exact and mapping_exact == candidate_exact:
        return "exact"
    if mapping_parent and mapping_parent == candidate_parent:
        return "parent"
    if _connectivity(mapping_parent) and _connectivity(mapping_parent) == _connectivity(candidate_parent):
        return "connectivity"
    return "mismatch" if candidate_exact or candidate_parent else "missing"


def _verify_selected_pdbqt(row: Mapping[str, str]) -> tuple[str, bool]:
    path_text = _clean(row.get("selected_pdbqt_path"))
    expected = _clean(row.get("selected_pdbqt_sha256")).lower()
    if not path_text or not expected:
        return "unusable_missing_path_or_hash", False
    path = Path(path_text).expanduser().resolve()
    if not path.is_file() or path.stat().st_size <= 0:
        return "unusable_missing_or_empty", False
    observed = _sha256(path)
    if observed != expected:
        raise ValueError(f"legacy PDBQT checksum mismatch: {path}")
    return "verified_match", True


def _read_full_drugcentral(path: Path, approved_ids: set[str]) -> _DCIndex:
    records: list[_DCRecord] = []
    by_exact: dict[str, list[_DCRecord]] = defaultdict(list)
    by_parent: dict[str, list[_DCRecord]] = defaultdict(list)
    by_connectivity: dict[str, list[_DCRecord]] = defaultdict(list)
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    for index, mol in enumerate(supplier, start=1):
        record = _drugcentral_record(mol, index, approved_ids)
        records.append(record)
        _append_dc_index(record, by_exact, by_parent, by_connectivity)
    return _DCIndex(
        records=tuple(records),
        by_exact={key: tuple(value) for key, value in by_exact.items()},
        by_parent={key: tuple(value) for key, value in by_parent.items()},
        by_connectivity={key: tuple(value) for key, value in by_connectivity.items()},
    )


def _drugcentral_record(
    mol: Chem.Mol | None, index: int, approved_ids: set[str]
) -> _DCRecord:
    if mol is None:
        raise ValueError(f"full DrugCentral SDF record failed to parse: {index}")
    structure = _structure_from_mol(mol, f"full_drugcentral:{index}")
    dc_id = _mol_prop(mol, "ID")
    preferred = _mol_prop(mol, "PREFERRED_NAME") or _mol_prop(mol, "_Name")
    if not (dc_id and preferred and structure.exact_inchikey and structure.parent_inchikey):
        raise ValueError(f"incomplete full DrugCentral record: {index}")
    return _DCRecord(
        index=index,
        drugcentral_id=dc_id,
        preferred_name=preferred,
        names=_dc_names(mol, preferred),
        exact_inchikey=structure.exact_inchikey,
        parent_inchikey=structure.parent_inchikey,
        connectivity_key=structure.parent_connectivity_key,
        approved_snapshot=dc_id in approved_ids,
    )


def _append_dc_index(
    record: _DCRecord,
    by_exact: dict[str, list[_DCRecord]],
    by_parent: dict[str, list[_DCRecord]],
    by_connectivity: dict[str, list[_DCRecord]],
) -> None:
    by_exact[record.exact_inchikey].append(record)
    by_parent[record.parent_inchikey].append(record)
    by_connectivity[record.connectivity_key].append(record)


def _read_pubchem_cache(
    records_path: Path, manifest_path: Path
) -> dict[str, _PubChemRecord]:
    manifest = _read_json(manifest_path)
    _validate_pubchem_manifest(records_path, manifest)
    output: dict[str, _PubChemRecord] = {}
    with records_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            record = _pubchem_record(json.loads(line), line_number)
            if record.cid in output:
                raise ValueError(f"duplicate PubChem CID {record.cid}")
            output[record.cid] = record
    if len(output) != _as_int(manifest.get("returned_count")):
        raise ValueError("PubChem returned_count differs from records")
    return output


def _read_exact_key_cache(
    records_path: Path, manifest_path: Path
) -> _ExactKeyCache:
    manifest = _read_json(manifest_path)
    _validate_exact_key_manifest(records_path, manifest)
    requested = {
        _clean(value).upper() for value in manifest.get("requested_inchikeys", [])
    }
    selected: dict[str, _PubChemRecord] = {}
    outcomes: dict[str, str] = {}
    query_by_rdk_id: dict[str, str] = {}
    with records_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            raw = json.loads(line)
            query = _clean(raw.get("query_inchikey")).upper()
            if not query or query in outcomes or query not in requested:
                raise ValueError(f"invalid exact-key cache query at line {line_number}")
            outcomes[query] = _clean(raw.get("outcome"))
            _index_exact_query_rows(raw, query, query_by_rdk_id)
            record = _select_exact_key_record(raw, query, line_number)
            if record is not None:
                selected[query] = record
    if set(outcomes) != requested:
        raise ValueError("exact-key cache coverage differs from requested keys")
    if Counter(outcomes.values()) != {
        str(key): _as_int(value)
        for key, value in dict(manifest.get("outcome_counts", {})).items()
    }:
        raise ValueError("exact-key cache outcome counts differ from manifest")
    return _ExactKeyCache(
        by_exact=selected,
        outcomes=outcomes,
        query_by_rdk_id=query_by_rdk_id,
    )


def _validate_exact_key_manifest(
    records_path: Path, manifest: Mapping[str, Any]
) -> None:
    if manifest.get("schema") != "atlas.pubchem-inchikey-cache.v1":
        raise ValueError("unsupported exact-key cache manifest schema")
    if _sha256(records_path) != _clean(manifest.get("records_sha256")).lower():
        raise ValueError("exact-key cache records hash differs from manifest")
    if records_path.stat().st_size != _as_int(manifest.get("records_bytes")):
        raise ValueError("exact-key cache records byte count differs from manifest")
    if _as_int(manifest.get("invalid_structure_query_count")):
        raise ValueError("exact-key cache contains invalid structures")
    requested = manifest.get("requested_inchikeys", [])
    if len(requested) != len(set(requested)):
        raise ValueError("exact-key cache request list contains duplicates")
    if _as_int(manifest.get("covered_unique_inchikeys")) != len(requested):
        raise ValueError("exact-key cache manifest coverage is incomplete")
    if float(str(manifest.get("coverage_percent", 0))) != 100.0:
        raise ValueError("exact-key cache coverage is not 100 percent")
    if sum(
        _as_int(value)
        for value in dict(manifest.get("outcome_counts", {})).values()
    ) != len(requested):
        raise ValueError("exact-key cache outcome total differs from requests")


def _index_exact_query_rows(
    raw: Mapping[str, Any], query: str, query_by_rdk_id: dict[str, str]
) -> None:
    for input_row in raw.get("input_rows", []):
        rdk_id = _clean(input_row.get("rdk_id"))
        existing = query_by_rdk_id.get(rdk_id)
        if rdk_id and existing not in {None, query}:
            raise ValueError(f"exact-key cache has conflicting keys for {rdk_id}")
        if rdk_id:
            query_by_rdk_id[rdk_id] = query


def _select_exact_key_record(
    raw: Mapping[str, Any], query: str, line_number: int
) -> _PubChemRecord | None:
    raw_records = list(raw.get("records", []))
    outcome = _clean(raw.get("outcome"))
    if outcome == "no_hit":
        if raw_records:
            raise ValueError(f"no-hit exact-key query has records: {query}")
        return None
    parsed = [_pubchem_record(value, line_number) for value in raw_records]
    if not parsed or any(record.exact_inchikey != query for record in parsed):
        raise ValueError(f"exact-key hit does not match query: {query}")
    if not _exact_key_names_agree(parsed):
        return None
    return min(parsed, key=lambda record: _as_int(record.cid))


def _exact_key_names_agree(records: Sequence[_PubChemRecord]) -> bool:
    if len(records) <= 1:
        return True
    names = {_name_key(record.title or record.iupac_name) for record in records}
    return len(names) == 1 and "" not in names


def _validate_pubchem_manifest(
    records_path: Path, manifest: Mapping[str, Any]
) -> None:
    if manifest.get("schema") != "atlas.pubchem-cid-cache.v1":
        raise ValueError("unsupported PubChem cache manifest schema")
    if _sha256(records_path) != _clean(manifest.get("records_sha256")).lower():
        raise ValueError("PubChem records hash differs from manifest")
    if records_path.stat().st_size != _as_int(manifest.get("records_bytes")):
        raise ValueError("PubChem records byte count differs from manifest")
    if _as_int(manifest.get("missing_count")) or _as_int(
        manifest.get("invalid_structure_count")
    ):
        raise ValueError("PubChem cache is incomplete or has invalid structures")


def _pubchem_record(raw: Mapping[str, Any], line_number: int) -> _PubChemRecord:
    cid = str(raw.get("CID", "")).strip()
    structure = _structure_from_smiles(
        _clean(raw.get("SMILES")), f"pubchem:{cid}"
    )
    if not cid or not structure.exact_inchikey or not structure.parent_inchikey:
        raise ValueError(f"invalid PubChem record at line {line_number}")
    if _clean(raw.get("InChIKey")).upper() != structure.exact_inchikey:
        raise ValueError(f"PubChem InChIKey mismatch for CID {cid}")
    return _PubChemRecord(
        cid=cid,
        title=_clean(raw.get("Title")),
        iupac_name=_clean(raw.get("IUPACName")),
        exact_inchikey=structure.exact_inchikey,
        parent_inchikey=structure.parent_inchikey,
        connectivity_key=structure.parent_connectivity_key,
        source_content_sha256=_clean(raw.get("_source_content_sha256")),
    )


def _read_and_validate_named_library(
    path: Path,
    summary_path: Path,
) -> dict[str, Mapping[str, str]]:
    _, rows = _read_csv(path)
    counts = Counter(row["materialization_status"] for row in rows)
    summary = _read_json(summary_path)
    _validate_named_library_summary(path, rows, counts, summary)
    output: dict[str, Mapping[str, str]] = {}
    for row in rows:
        _append_named_library_row(output, row)
    return output


def _validate_named_library_summary(
    path: Path,
    rows: Sequence[Mapping[str, str]],
    counts: Counter[str],
    summary: Mapping[str, Any],
) -> None:
    expected_counts = {
        str(key): _as_int(value)
        for key, value in dict(summary.get("status_counts", {})).items()
    }
    if not expected_counts or dict(counts) != expected_counts:
        raise ValueError(f"canonical named-library parity mismatch: {dict(counts)}")
    if _clean(summary.get("manifest_sha256")).lower() != _sha256(path):
        raise ValueError("canonical named-library manifest hash differs from summary")
    if _as_int(summary.get("canonical_parent_rows")) != len(rows):
        raise ValueError("canonical parent row count differs from summary")
    ready_count = counts.get("ready", 0)
    if ready_count != _as_int(summary.get("ready_pdbqt_count")):
        raise ValueError("canonical ready count differs from summary")
    inventory = dict(summary.get("library_inventory", {}))
    if ready_count != _as_int(inventory.get("count")):
        raise ValueError("canonical inventory count differs from summary")


def _append_named_library_row(
    output: dict[str, Mapping[str, str]], row: Mapping[str, str]
) -> None:
    parent = row["parent_inchikey"]
    if not parent or parent in output:
        raise ValueError("canonical named library requires unique parent keys")
    if row["materialization_status"] == "ready":
        _validate_manifest_file(
            row["output_pdbqt_path"],
            row["output_pdbqt_sha256"],
            row["output_pdbqt_bytes"],
        )
    output[parent] = row


def _read_unsafe_parent_quarantine(
    path: Path,
) -> _UnsafeParentQuarantine:
    _, rows = _read_csv(path)
    output: dict[str, Mapping[str, str]] = {}
    source_indices: set[int] = set()
    source_row_count = 0
    for row in rows:
        if row.get("scope") != "approved_source_record":
            continue
        source_row_count += 1
        if row.get("reason_codes") != "unsafe_parent_collapse":
            raise ValueError("unsupported approved-source quarantine reason")
        parent = row.get("parent_inchikey", "")
        if not parent:
            raise ValueError("approved-source quarantine row lacks parent InChIKey")
        source_indices.add(_as_int(row.get("source_record_indices")))
        output.setdefault(parent, _unsafe_parent_manifest_row(row, parent))
    if not output or 0 in source_indices or len(source_indices) != source_row_count:
        raise ValueError("repair quarantine contains no approved-source rows")
    return _UnsafeParentQuarantine(
        by_parent=output,
        source_indices=frozenset(source_indices),
        row_count=source_row_count,
    )


def _unsafe_parent_manifest_row(
    row: Mapping[str, str], parent: str
) -> Mapping[str, str]:
    return {
        "canonical_parent_id": f"fda_unsafe_parent_{parent}",
        "preferred_name": row.get("display_or_preferred_names", ""),
        "parent_inchikey": parent,
        "source_record_indices": row.get("source_record_indices", ""),
        "approved_source_classes": "unsafe_approved_source",
        "approved_form_relation": "unsafe_parent_collapse",
        "materialization_status": "quarantined_unsafe_parent_collapse",
        "status_reason_codes": "repair_quarantine_source_record_linked;unsafe_parent_collapse",
        "output_pdbqt_path": "",
        "output_pdbqt_sha256": "",
    }


def _validate_approved_source_partition(
    source_sdf: Path,
    canonical: Mapping[str, Mapping[str, str]],
    unsafe: _UnsafeParentQuarantine,
) -> int:
    supplier = Chem.SDMolSupplier(
        str(source_sdf), removeHs=False, sanitize=False
    )
    source_count = len(supplier)
    if source_count == 0 or any(mol is None for mol in supplier):
        raise ValueError("filtered approved source SDF is empty or unparsable")
    named_indices = {
        index
        for row in canonical.values()
        for index in _source_record_indices(row.get("source_record_indices", ""))
    }
    expected = set(range(1, source_count + 1))
    if named_indices.intersection(unsafe.source_indices):
        raise ValueError("approved source partition overlaps named and unsafe rows")
    if named_indices.union(unsafe.source_indices) != expected:
        raise ValueError("approved source partition is incomplete")
    return source_count


def _source_record_indices(value: str) -> set[int]:
    indices = {
        _as_int(token)
        for token in re.split(r"[;,|\s]+", value)
        if token.strip()
    }
    if not indices or 0 in indices:
        raise ValueError(f"invalid source record index set: {value}")
    return indices


def _validate_manifest_file(path_text: str, digest: str, size: str) -> None:
    path = Path(path_text).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"canonical PDBQT missing: {path}")
    if path.stat().st_size != _as_int(size) or _sha256(path) != digest.lower():
        raise ValueError(f"canonical PDBQT checksum/size mismatch: {path}")


def _read_approved_snapshot(path: Path) -> tuple[set[str], set[str]]:
    ids: set[str] = set()
    names: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.reader(handle):
            if row and _clean(row[0]):
                ids.add(_clean(row[0]))
            if len(row) > 1 and _name_key(row[1]):
                names.add(_name_key(row[1]))
    if not ids or not names:
        raise ValueError("FDA Approved snapshot has no usable IDs or names")
    return ids, names


def _read_official_ingredients(path: Path) -> set[str]:
    output: set[str] = set()
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "ActiveIngredient" not in (reader.fieldnames or []):
            raise ValueError("Drugs@FDA Products.txt lacks ActiveIngredient")
        for row in reader:
            for ingredient in re.split(r"[;,/]", row.get("ActiveIngredient", "")):
                key = _name_key(ingredient)
                if key:
                    output.add(key)
                    output.add(_strip_form_suffix(key))
    if not output:
        raise ValueError("Drugs@FDA Products.txt contains no active ingredients")
    return output


def _official_name_match(name: str, ingredients: set[str]) -> bool:
    key = _name_key(name)
    return bool(key and (key in ingredients or _strip_form_suffix(key) in ingredients))


def _strip_form_suffix(key: str) -> str:
    value = key
    changed = True
    while changed:
        changed = False
        for suffix in _COMMON_FORM_SUFFIXES:
            token = _name_key(suffix)
            if value.endswith(token) and len(value) > len(token) + 2:
                value = value[: -len(token)]
                changed = True
                break
    return value


def _row_name_keys(row: Mapping[str, str]) -> set[str]:
    fields = (
        "display_name",
        "generic_name",
        "pubchem_name",
        "pubchem_record_title",
        "rxnorm_generic_name",
    )
    return {_name_key(row.get(field, "")) for field in fields if _name_key(row.get(field, ""))}


def _name_disagrees(row: Mapping[str, str], preferred: str) -> bool:
    preferred_key = _strip_form_suffix(_name_key(preferred))
    legacy_keys = {
        _strip_form_suffix(key) for key in _row_name_keys(row) if key
    }
    return bool(preferred_key and legacy_keys and preferred_key not in legacy_keys)


def _identity_evidence_codes(identity: _IdentityResolution) -> str:
    if identity.source in {
        "mapping_exact_stereoisomer_pdbqt_lineage",
        "mapping_exact_structure_pdbqt_lineage",
    }:
        return identity.source
    if identity.source == "manifest_approved_v2":
        return "manifest_approved_structure"
    if identity.source == "full_drugcentral_structure":
        return f"drugcentral_{identity.relation}"
    if identity.relation in {"exact", "parent"}:
        return f"pubchem_{identity.relation}"
    if identity.relation == "connectivity" and identity.structure_validated:
        return "pubchem_validated_connectivity_name"
    return ""


def _approval_evidence_codes(identity: _IdentityResolution) -> str:
    if identity.regulatory_status not in {
        "drugcentral_fda_approved",
        "fda_approved_current_or_historical",
    }:
        return ""
    if identity.source == "manifest_approved_v2":
        return "drugcentral_approved_manifest"
    if identity.source == "full_drugcentral_structure":
        return "drugcentral_approved_id"
    return "drugsatfda_active_ingredient"


def _dc_names(mol: Chem.Mol, preferred: str) -> tuple[str, ...]:
    names = {preferred, _mol_prop(mol, "_Name")}
    for field in ("PREFERRED_NAME", "SYNONYMS"):
        value = _mol_prop(mol, field)
        names.update(part.strip() for part in re.split(r"[;|\r\n]+", value) if part.strip())
    return tuple(sorted(name for name in names if name))


def _mol_prop(mol: Chem.Mol, field: str) -> str:
    return _clean(mol.GetProp(field)) if mol.HasProp(field) else ""


def _unique_record(records: Sequence[_DCRecord]) -> _DCRecord | None:
    by_id = {record.drugcentral_id: record for record in records}
    return next(iter(by_id.values())) if len(by_id) == 1 else None


def _connectivity(inchikey: str) -> str:
    return inchikey.split("-", 1)[0] if inchikey else ""


def _name_key(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", _clean(value).casefold())


def _validate_mapping_partition(
    mapping_rows: Sequence[Mapping[str, str]],
    unresolved_rows: Sequence[Mapping[str, str]],
) -> dict[str, Mapping[str, str]]:
    mapping_ids = [row.get("rdk_id", "") for row in mapping_rows]
    if not all(mapping_ids) or len(set(mapping_ids)) != len(mapping_ids):
        raise ValueError("v2 mapping requires unique nonempty rdk_id values")
    unresolved = {row.get("rdk_id", ""): row for row in unresolved_rows}
    expected = {
        row["rdk_id"]
        for row in mapping_rows
        if row.get("identity_resolution_status") == "unresolved_requires_second_pass"
    }
    if set(unresolved) != expected or "" in unresolved:
        raise ValueError("v2 unresolved queue does not match v2 mapping statuses")
    return unresolved


def _validate_repaired_mapping(
    mapping_rows: Sequence[Mapping[str, str]],
    repaired_rows: Sequence[Mapping[str, str]],
) -> dict[str, Mapping[str, str]]:
    repaired = {row.get("rdk_id", ""): row for row in repaired_rows}
    mapping_ids = {row["rdk_id"] for row in mapping_rows}
    if len(repaired) != len(repaired_rows) or set(repaired) != mapping_ids:
        raise ValueError("repaired mapping is not a 1:1 RDK partition")
    mapping_by_id = {row["rdk_id"]: row for row in mapping_rows}
    for rdk_id, row in repaired.items():
        _validate_repaired_row(mapping_by_id[rdk_id], row)
    return repaired


def _validate_repaired_row(
    mapping_row: Mapping[str, str], repaired_row: Mapping[str, str]
) -> None:
    parity_fields = (
        "mapping_row_number",
        "selected_pdbqt_path",
        "selected_pdbqt_sha256",
    )
    if any(mapping_row.get(field) != repaired_row.get(field) for field in parity_fields):
        raise ValueError(f"repaired mapping parity mismatch: {mapping_row['rdk_id']}")
    if repaired_row.get("prepared_connectivity_state") not in {
        "coordinate_connectivity_exact",
        "coordinate_connectivity_parent",
        "ambiguous",
        "unverifiable",
    }:
        raise ValueError(f"unsupported prepared connectivity state: {mapping_row['rdk_id']}")


def _score_reuse_status(
    connectivity_state: str, canonical_row: Mapping[str, str]
) -> str:
    if canonical_row.get("materialization_status") == "quarantined_unsafe_parent_collapse":
        return "no_reuse_unsafe_parent_collapse"
    if connectivity_state == "coordinate_connectivity_exact":
        return "reuse_supported_coordinate_exact"
    if connectivity_state == "coordinate_connectivity_parent":
        return "reuse_supported_coordinate_parent"
    if connectivity_state == "ambiguous":
        return "no_reuse_ambiguous_connectivity"
    return "no_reuse_missing_file"


def _validate_terminal_rows(
    mapping_rows: Sequence[Mapping[str, str]],
    dispositions: Sequence[Mapping[str, Any]],
) -> None:
    if len(dispositions) != len(mapping_rows):
        raise RuntimeError("v3 terminal disposition cardinality changed")
    ids = [str(row["rdk_id"]) for row in dispositions]
    if len(set(ids)) != len(ids):
        raise RuntimeError("v3 terminal dispositions contain duplicate rdk_id")
    if any(not row["terminal_disposition"] for row in dispositions):
        raise RuntimeError("v3 row lacks terminal disposition")


def _validate_regression_rows(
    mapping_rows: Sequence[Mapping[str, str]],
    dispositions: Sequence[Mapping[str, Any]],
) -> None:
    _validate_row_4714(dispositions)
    _validate_cisapride_rows(mapping_rows, dispositions)
    _validate_lobaplatin_row(dispositions)


def _validate_row_4714(dispositions: Sequence[Mapping[str, Any]]) -> None:
    by_row = {row["mapping_row_number"]: row for row in dispositions}
    row_4714 = by_row.get("4714")
    if row_4714 is None or _name_key(row_4714["preferred_identity"]) != "deserpidine":
        raise RuntimeError("row 4714 preferred identity regression")
    if row_4714["identity_structure_validated"] != "true":
        raise RuntimeError("row 4714 structure validation regression")


def _validate_cisapride_rows(
    mapping_rows: Sequence[Mapping[str, str]],
    dispositions: Sequence[Mapping[str, Any]],
) -> None:
    cisapride_ids = {
        row["rdk_id"]
        for row in mapping_rows
        if "cisapride" in _row_name_keys(row)
    }
    if not cisapride_ids:
        raise RuntimeError("cisapride regression cohort is absent")
    fda_statuses = {
        "drugcentral_fda_approved",
        "fda_approved_current_or_historical",
    }
    cisapride_rows = [
        row for row in dispositions if row["rdk_id"] in cisapride_ids
    ]
    if any(
        row["identity_structure_validated"] != "true"
        or row["regulatory_status"] not in fda_statuses
        for row in cisapride_rows
    ):
        raise RuntimeError("cisapride FDA false-negative regression")


def _validate_lobaplatin_row(
    dispositions: Sequence[Mapping[str, Any]],
) -> None:
    rows = [row for row in dispositions if row["rdk_id"] == "rdk_0004569"]
    if len(rows) != 1:
        raise RuntimeError("lobaplatin regression row is absent or duplicated")
    row = rows[0]
    if _name_key(row["preferred_identity"]) != "lobaplatin":
        raise RuntimeError("lobaplatin exact-key preferred identity regression")
    if row["identity_evidence_codes"] != "pubchem_exact":
        raise RuntimeError("lobaplatin exact-key evidence regression")


def _summary(
    *,
    inputs: Mapping[str, Path],
    outputs: FDATerminalV3Outputs,
    mapping_rows: Sequence[Mapping[str, str]],
    unresolved_rows: Sequence[Mapping[str, str]],
    dispositions: Sequence[Mapping[str, Any]],
    dc_index: _DCIndex,
    pubchem: Mapping[str, _PubChemRecord],
    exact_key_cache: _ExactKeyCache,
    canonical: Mapping[str, Mapping[str, str]],
    unsafe_parents: _UnsafeParentQuarantine,
    approved_source_count: int,
) -> dict[str, Any]:
    usable_blockers, canonical_blockers, unexpected_dispositions = (
        _completion_blockers(dispositions)
    )
    blockers = [
        *usable_blockers,
        *canonical_blockers,
        *unexpected_dispositions,
    ]
    complete = not blockers
    output_map = _summary_output_map(outputs)
    return {
        "schema_version": SCHEMA_VERSION,
        "resolver_version": RESOLVER_VERSION,
        "completion_status": "complete" if complete else "blocked",
        "complete": complete,
        "input_paths": {key: str(path) for key, path in inputs.items()},
        "input_sha256": {key: _sha256(path) for key, path in inputs.items()},
        "input_bytes": {key: path.stat().st_size for key, path in inputs.items()},
        "output_paths": {key: str(path) for key, path in output_map.items()},
        "output_sha256": {key: _sha256(path) for key, path in output_map.items()},
        "counts": _summary_counts(
            mapping_rows=mapping_rows,
            unresolved_rows=unresolved_rows,
            dispositions=dispositions,
            dc_index=dc_index,
            pubchem=pubchem,
            exact_key_cache=exact_key_cache,
            canonical=canonical,
            unsafe_parents=unsafe_parents,
            approved_source_count=approved_source_count,
            usable_blockers=usable_blockers,
            canonical_blockers=canonical_blockers,
            unexpected_dispositions=unexpected_dispositions,
        ),
        "blocking_rdk_ids": sorted({row["rdk_id"] for row in blockers}),
        "policy": _summary_policy(),
    }


def _completion_blockers(
    dispositions: Sequence[Mapping[str, Any]],
) -> tuple[list[Mapping[str, Any]], list[Mapping[str, Any]], list[Mapping[str, Any]]]:
    usable_blockers = [
        row
        for row in dispositions
        if row["legacy_file_usable"] == "true"
        and row["identity_structure_validated"] != "true"
    ]
    canonical_blockers = [
        row
        for row in dispositions
        if row["regulatory_status"]
        in {"drugcentral_fda_approved", "fda_approved_current_or_historical"}
        and not row["canonical_parent_id"]
    ]
    unexpected_dispositions = [
        row
        for row in dispositions
        if row["terminal_disposition"] not in TERMINAL_DISPOSITIONS
    ]
    return usable_blockers, canonical_blockers, unexpected_dispositions


def _summary_output_map(outputs: FDATerminalV3Outputs) -> dict[str, Path]:
    return {
        "mapping_v3_csv": outputs.mapping_v3_csv,
        "terminal_dispositions_csv": outputs.terminal_dispositions_csv,
        "fda_subset_csv": outputs.fda_subset_csv,
        "excluded_subset_csv": outputs.excluded_subset_csv,
    }


def _summary_counts(
    *,
    mapping_rows: Sequence[Mapping[str, str]],
    unresolved_rows: Sequence[Mapping[str, str]],
    dispositions: Sequence[Mapping[str, Any]],
    dc_index: _DCIndex,
    pubchem: Mapping[str, _PubChemRecord],
    exact_key_cache: _ExactKeyCache,
    canonical: Mapping[str, Mapping[str, str]],
    unsafe_parents: _UnsafeParentQuarantine,
    approved_source_count: int,
    usable_blockers: Sequence[Mapping[str, Any]],
    canonical_blockers: Sequence[Mapping[str, Any]],
    unexpected_dispositions: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    return {
        "mapping_rows": len(mapping_rows),
        "v2_unresolved_rows": len(unresolved_rows),
        "full_drugcentral_records": len(dc_index.records),
        "pubchem_records": len(pubchem),
        "pubchem_exact_key_queries": len(exact_key_cache.outcomes),
        "pubchem_exact_key_selected_records": len(exact_key_cache.by_exact),
        "pubchem_exact_key_outcomes": dict(
            sorted(Counter(exact_key_cache.outcomes.values()).items())
        ),
        "canonical_parent_rows": len(canonical),
        "approved_source_records": approved_source_count,
        "named_approved_source_records": (
            approved_source_count - unsafe_parents.row_count
        ),
        "unsafe_approved_source_records": unsafe_parents.row_count,
        "unsafe_approved_source_parents": len(unsafe_parents.by_parent),
        "canonical_ready_rows": sum(
            row["materialization_status"] == "ready" for row in canonical.values()
        ),
        "canonical_quarantine_rows": sum(
            row["materialization_status"] == "quarantined_non_dockable"
            for row in canonical.values()
        ),
        "usable_rows": sum(
            row["legacy_file_usable"] == "true" for row in dispositions
        ),
        "usable_rows_without_structure_validated_identity": len(usable_blockers),
        "fda_rows_without_canonical_parent": len(canonical_blockers),
        "unexpected_terminal_dispositions": len(unexpected_dispositions),
        "identity_source": _field_counts(dispositions, "preferred_identity_source"),
        "identity_structure_relation": _field_counts(
            dispositions, "identity_structure_relation"
        ),
        "regulatory_status": _field_counts(dispositions, "regulatory_status"),
        "terminal_disposition": _field_counts(dispositions, "terminal_disposition"),
        "prepared_connectivity_state": _field_counts(
            dispositions, "prepared_connectivity_state"
        ),
        "legacy_score_reuse_status": _field_counts(
            dispositions, "legacy_score_reuse_status"
        ),
    }


def _summary_policy() -> dict[str, bool]:
    return {
        "identity_resolution_is_structure_first": True,
        "pubchem_identity_implies_fda": False,
        "unknown_identity_forced_non_fda": False,
        "every_row_has_terminal_disposition": True,
        "usable_row_requires_structure_validated_preferred_identity": True,
        "fda_row_requires_canonical_parent_link": True,
        "legacy_mapping_v2_mutated": False,
        "legacy_pdbqt_files_mutated": False,
    }


def _resolved_inputs(**values: Path) -> dict[str, Path]:
    output = {key: Path(path).expanduser().resolve() for key, path in values.items()}
    for key, path in output.items():
        if not path.is_file():
            raise FileNotFoundError(f"{key} is not a file: {path}")
    return output


def _validate_output_scope(
    outputs: FDATerminalV3Outputs, input_paths: Iterable[Path]
) -> None:
    output_paths = set(_output_map(outputs).values())
    if len(output_paths) != 5 or output_paths.intersection(set(input_paths)):
        raise ValueError("v3 outputs collide with each other or inputs")


def _output_map(outputs: FDATerminalV3Outputs) -> dict[str, Path]:
    return {
        "mapping_v3_csv": outputs.mapping_v3_csv,
        "terminal_dispositions_csv": outputs.terminal_dispositions_csv,
        "fda_subset_csv": outputs.fda_subset_csv,
        "excluded_subset_csv": outputs.excluded_subset_csv,
        "summary_json": outputs.summary_json,
    }


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError(f"CSV has no header: {path}")
        return fields, [
            {key: _clean(value) for key, value in row.items()} for row in reader
        ]


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"JSON root is not an object: {path}")
    return payload


def _write_csv(
    path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _new_fields(existing: Sequence[str]) -> list[str]:
    seen = set(existing)
    return [field for field in _EXTRA_MAPPING_FIELDS if field not in seen]


def _field_counts(
    rows: Sequence[Mapping[str, Any]], field: str
) -> dict[str, int]:
    return dict(sorted(Counter(str(row.get(field, "")) for row in rows).items()))


def _joined(values: Iterable[object]) -> str:
    return ";".join(dict.fromkeys(str(value).strip() for value in values if str(value).strip()))


def _clean(value: object) -> str:
    text = str(value or "").strip()
    return "" if text.casefold() == "nan" else text


def _as_int(value: object) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError):
        return 0


def _bool(value: bool) -> str:
    return "true" if value else "false"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = [
    "FDATerminalV3Outputs",
    "RESOLVER_VERSION",
    "resolve_fda_terminal_v3",
    "terminal_v3_paths",
]
