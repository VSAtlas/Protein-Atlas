from __future__ import annotations

import csv
import hashlib
import json
import math
import numbers
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, cast

import pandas as pd


CONSENSUS_PREFIX = "z_vs_compare_run_consensus"
SCORCH_PREFIX = "z_vs_compare_run_scorch_composite"
PAIR_COLUMNS = ("pdb_id", "ligand_base")
LABEL_STATUS = "strict_measured_matched_inactive"
MISSING_TEXT = frozenset({"", "na", "n/a", "nan", "none", "null", "<na>"})
PROVENANCE_SUFFIXES = (
    "_source",
    "_run_id",
    "_run_sha256",
    "_null_sha256",
    "_decoy_n",
    "_decoy_mu",
    "_decoy_sigma",
)

_KNOWN_SUFFIX = re.compile(r"\.(?:pdbqt|mol2|sdf)(?:\.gz)?$", re.IGNORECASE)
_LEGACY_STAGE_SUFFIX = re.compile(
    r"(?:_(?:dud_)?(?:gnina|dock6|ledock)?_?stage\d+)+$", re.IGNORECASE
)

STANDARD_PROVENANCE_COLUMNS = (
    "source_objective",
    "label_source",
    "source_family",
    "upstream_source",
    "source_label_policy",
    "assay_type",
    "endpoint_type",
    "training_allowed",
    "matched_inactive_pair_key",
    "matched_inactive_score_source",
    "matched_inactive_selected_row",
    "matched_inactive_score_row",
    "matched_inactive_selected_target_uniprot",
    "matched_inactive_selected_target_gene",
    "matched_inactive_selected_generic_name",
    "matched_inactive_selected_display_name",
    "matched_inactive_target_mapping_status",
    "matched_inactive_target_mapping_key",
)

TARGET_AUTHORITATIVE_METADATA = (
    "target_family",
    "protein_class",
    "structure_quality",
)
LIGAND_AUTHORITATIVE_METADATA = (
    "chemical_cluster",
    "butina_cluster",
    "scaffold_key",
    "ligand_chemotype",
)
LIGAND_DESCRIPTOR_COLUMNS = frozenset(
    {
        "molecular_weight",
        "mol_wt",
        "logp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "formal_charge",
        "aromatic_rings",
        "fraction_csp3",
        "qed",
    }
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        if bool(pd.isna(value)):
            return True
    except (TypeError, ValueError):
        pass
    return str(value).strip().casefold() in MISSING_TEXT


def _finite_float(value: Any) -> float | None:
    if _is_missing(value):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _clean_text(value: Any) -> str:
    return "" if _is_missing(value) else str(value).strip()


def _first_nonmissing(row: Mapping[str, Any], *columns: str) -> Any:
    for column in columns:
        value = row.get(column)
        if not _is_missing(value):
            return value
    return ""


def _normalize_pdb_id(value: Any) -> str:
    return "" if _is_missing(value) else str(value).strip().upper()


def _normalize_target_identity(value: Any) -> str:
    return "" if _is_missing(value) else str(value).strip().casefold()


def _normalize_ligand_base(value: Any) -> str:
    if _is_missing(value):
        return ""
    text = str(value).strip().replace("\\", "/").rsplit("/", 1)[-1]
    text = text.replace(".sanitized", "")
    text = _KNOWN_SUFFIX.sub("", text)
    text = _LEGACY_STAGE_SUFFIX.sub("", text)
    text = re.sub(r"__+", "_", text).rstrip("_")
    return text.casefold()


def _pair_key(pdb_id: Any, ligand_base: Any) -> tuple[str, str] | None:
    key = (_normalize_pdb_id(pdb_id), _normalize_ligand_base(ligand_base))
    return key if all(key) else None


def _key_text(key: tuple[str, str] | None) -> str:
    return "" if key is None else f"{key[0]}||{key[1]}"


def _read_csv_strict(path: Path, *, name: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{name} CSV does not exist: {path}")
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{name} CSV is empty: {path}") from exc
    duplicates = sorted(
        column for column, count in Counter(header).items() if count > 1
    )
    if duplicates:
        raise ValueError(f"{name} CSV has duplicate columns {duplicates}: {path}")
    if not header or any(not str(column).strip() for column in header):
        raise ValueError(f"{name} CSV has a blank or missing header: {path}")
    return pd.read_csv(path, low_memory=False)


def _require_columns(
    frame: pd.DataFrame, required: Iterable[str], *, name: str
) -> None:
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{name} is missing required columns: {', '.join(missing)}")


def _comparison_columns(frame: pd.DataFrame) -> list[str]:
    return [
        column
        for column in frame.columns
        if column == CONSENSUS_PREFIX
        or column.startswith(f"{CONSENSUS_PREFIX}_")
        or column == SCORCH_PREFIX
        or column.startswith(f"{SCORCH_PREFIX}_")
    ]


def _canonical_value(value: Any) -> str:
    if _is_missing(value):
        return "missing:"
    if isinstance(value, bool):
        return f"bool:{str(value).lower()}"
    if isinstance(value, numbers.Number):
        number = float(cast(Any, value))
        if math.isfinite(number):
            return f"number:{format(number, '.17g')}"
    return f"text:{str(value).strip()}"


def _row_payload(row: Mapping[str, Any], columns: Iterable[str]) -> tuple[str, ...]:
    return tuple(_canonical_value(row.get(column)) for column in columns)


def _frame_hash(frame: pd.DataFrame, columns: Iterable[str]) -> str:
    selected_columns = list(columns)
    digest = hashlib.sha256()
    digest.update(json.dumps(selected_columns, separators=(",", ":")).encode())
    for values in frame.loc[:, selected_columns].itertuples(index=False, name=None):
        digest.update(
            json.dumps(
                [_canonical_value(value) for value in values],
                ensure_ascii=True,
                separators=(",", ":"),
            ).encode("utf-8")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _nonmissing_hash(frame: pd.DataFrame, columns: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for row_number, values in enumerate(
        frame.loc[:, list(columns)].itertuples(index=False, name=None)
    ):
        for column, value in zip(columns, values, strict=False):
            if _is_missing(value):
                continue
            digest.update(
                f"{row_number}\0{column}\0{_canonical_value(value)}\n".encode("utf-8")
            )
    return digest.hexdigest()


def _label_contract_valid(row: Mapping[str, Any]) -> bool:
    label = _finite_float(row.get("projected_label"))
    status = _clean_text(row.get("projected_label_status")).casefold()
    return label == 0.0 and status == LABEL_STATUS


def _truthy(value: Any) -> bool:
    return _clean_text(value).casefold() in {"1", "true", "t", "yes", "y"}


def _training_allowed(row: Mapping[str, Any]) -> bool:
    production = _truthy(row.get("production_truth_allowed"))
    benchmark = _truthy(row.get("benchmark_only"))
    probe_only = _truthy(row.get("probe_sensitivity_only"))
    return production and not benchmark and not probe_only


def _provenance_errors(row: Mapping[str, Any], prefix: str) -> list[str]:
    errors: list[str] = []
    for suffix in PROVENANCE_SUFFIXES:
        column = f"{prefix}{suffix}"
        value = row.get(column)
        if suffix in {"_decoy_n", "_decoy_mu", "_decoy_sigma"}:
            if _finite_float(value) is None:
                errors.append(column)
        elif _is_missing(value):
            errors.append(column)
    decoy_n = _finite_float(row.get(f"{prefix}_decoy_n"))
    decoy_sigma = _finite_float(row.get(f"{prefix}_decoy_sigma"))
    if decoy_n is not None and decoy_n <= 0:
        errors.append(f"{prefix}_decoy_n_nonpositive")
    if decoy_sigma is not None and decoy_sigma <= 0:
        errors.append(f"{prefix}_decoy_sigma_nonpositive")
    return errors


def _score_choice(row: Mapping[str, Any]) -> tuple[float | None, str, list[str]]:
    scorch = _finite_float(row.get(SCORCH_PREFIX))
    consensus = _finite_float(row.get(CONSENSUS_PREFIX))
    if scorch is not None:
        return (
            scorch,
            "comparison_run_scorch_composite",
            _provenance_errors(row, SCORCH_PREFIX),
        )
    if consensus is not None:
        return (
            consensus,
            "comparison_run_consensus",
            _provenance_errors(row, CONSENSUS_PREFIX),
        )
    return None, "", []


def _duplicate_score_audit(
    scores: pd.DataFrame,
    keys: list[tuple[str, str] | None],
    comparison_columns: list[str],
) -> tuple[set[tuple[str, str]], list[dict[str, Any]]]:
    valid_positions: dict[tuple[str, str], list[int]] = {}
    rows: list[dict[str, Any]] = []
    for position, key in enumerate(keys):
        if key is None:
            rows.append(
                {
                    "score_row": position,
                    "normalized_pair_key": "",
                    "conflict_type": "invalid_score_pair_key",
                    "duplicate_count": 0,
                    "payload_hash": "",
                }
            )
            continue
        valid_positions.setdefault(key, []).append(position)

    rejected: set[tuple[str, str]] = set()
    for key, positions in sorted(valid_positions.items()):
        if len(positions) == 1:
            continue
        payloads = [
            _row_payload(scores.iloc[position].to_dict(), comparison_columns)
            for position in positions
        ]
        conflict_type = (
            "duplicate_identical_score_key"
            if len(set(payloads)) == 1
            else "duplicate_conflicting_score_key"
        )
        rejected.add(key)
        for position, payload in zip(positions, payloads, strict=False):
            rows.append(
                {
                    "score_row": position,
                    "pdb_id": key[0],
                    "ligand_base": key[1],
                    "normalized_pair_key": _key_text(key),
                    "conflict_type": conflict_type,
                    "duplicate_count": len(positions),
                    "payload_hash": _json_sha256(payload),
                }
            )
    return rejected, rows


def _selected_duplicate_keys(
    keys: list[tuple[str, str] | None],
) -> set[tuple[str, str]]:
    counts = Counter(key for key in keys if key is not None)
    return {key for key, count in counts.items() if count > 1}


def _ligand_authoritative_columns(base: pd.DataFrame) -> list[str]:
    columns = ["drug_id"]
    for column in base.columns:
        normalized = column.casefold()
        if (
            column in LIGAND_AUTHORITATIVE_METADATA
            or normalized.startswith("rdkit_")
            or normalized in LIGAND_DESCRIPTOR_COLUMNS
        ):
            columns.append(column)
    return list(dict.fromkeys(columns))


def _target_authoritative_columns(base: pd.DataFrame) -> list[str]:
    return [
        column
        for column in ("target_id", "target_gene", *TARGET_AUTHORITATIVE_METADATA)
        if column in base.columns
    ]


def _build_base_identity_map(
    base: pd.DataFrame,
    *,
    key_column: str,
    identity_column: str,
    value_columns: list[str],
    normalize_key: Callable[[Any], str],
    map_type: str,
) -> tuple[dict[str, dict[str, Any]], dict[str, str], list[dict[str, Any]]]:
    positions_by_key: dict[str, list[int]] = {}
    audit_rows: list[dict[str, Any]] = []
    for position, value in enumerate(base[key_column]):
        key = normalize_key(value)
        if not key:
            audit_rows.append(
                {
                    "map_type": map_type,
                    "normalized_key": "",
                    "base_row_count": 1,
                    "mapping_status": "invalid_base_key",
                    "conflict_columns": "",
                    "base_row": position,
                }
            )
            continue
        positions_by_key.setdefault(key, []).append(position)

    records: dict[str, dict[str, Any]] = {}
    statuses: dict[str, str] = {}
    for key, positions in sorted(positions_by_key.items()):
        values: dict[str, Any] = {}
        conflicts: list[str] = []
        for column in value_columns:
            unique: dict[str, Any] = {}
            for value in base.iloc[positions][column]:
                if _is_missing(value):
                    continue
                unique.setdefault(_canonical_value(value), value)
            if len(unique) > 1:
                conflicts.append(column)
                values[column] = pd.NA
            elif unique:
                values[column] = next(iter(unique.values()))
            else:
                values[column] = pd.NA

        if conflicts:
            status = "ambiguous_authoritative_metadata"
        elif _is_missing(values.get(identity_column)):
            status = "missing_authoritative_identity"
        else:
            status = "unique"
            records[key] = values
        statuses[key] = status
        audit_rows.append(
            {
                "map_type": map_type,
                "normalized_key": key,
                "base_row_count": len(positions),
                "mapping_status": status,
                "conflict_columns": ";".join(conflicts),
                "base_row": "",
                **{
                    f"authoritative_{column}": values.get(column)
                    for column in value_columns
                },
            }
        )
    return records, statuses, audit_rows


def _build_base_target_alias_map(
    base: pd.DataFrame,
    *,
    value_columns: list[str],
) -> tuple[dict[str, dict[str, Any]], dict[str, str], list[dict[str, Any]]]:
    alias_rows: list[dict[str, Any]] = []
    alias_columns = ["target_id"]
    if "target_gene" in base.columns:
        alias_columns.append("target_gene")
    for row in base.to_dict(orient="records"):
        aliases = {
            _normalize_target_identity(row.get(column)) for column in alias_columns
        }
        aliases.discard("")
        for alias in sorted(aliases):
            alias_rows.append(
                {
                    "_target_alias": alias,
                    **{column: row.get(column) for column in value_columns},
                }
            )
    alias_frame = pd.DataFrame(alias_rows, columns=["_target_alias", *value_columns])
    if alias_frame.empty:
        return {}, {}, []
    return _build_base_identity_map(
        alias_frame,
        key_column="_target_alias",
        identity_column="target_id",
        value_columns=value_columns,
        normalize_key=_normalize_target_identity,
        map_type="target_id_or_gene_to_target",
    )


def _candidate_target_record(
    row: Mapping[str, Any],
    *,
    pdb_map: Mapping[str, Mapping[str, Any]],
    alias_map: Mapping[str, Mapping[str, Any]],
) -> Mapping[str, Any] | None:
    status = _clean_text(row.get("matched_inactive_target_mapping_status"))
    mapping_key = _clean_text(row.get("matched_inactive_target_mapping_key"))
    if status == "authoritative_pdb_verified":
        return pdb_map.get(_normalize_pdb_id(mapping_key))
    if status == "authoritative_target_gene_verified":
        return alias_map.get(_normalize_target_identity(mapping_key))
    return None


def _identity_comparison(selected_value: Any, authoritative_value: Any) -> str:
    if _is_missing(selected_value):
        return "selected_missing"
    if _is_missing(authoritative_value):
        return "authoritative_missing"
    selected_text = str(selected_value).strip()
    authoritative_text = str(authoritative_value).strip()
    if selected_text == authoritative_text:
        return "exact_match"
    if selected_text.casefold() == authoritative_text.casefold():
        return "case_only_difference"
    return "different_identity_or_namespace"


def _candidate_row(
    selected_row: Mapping[str, Any],
    score_row: Mapping[str, Any],
    *,
    selected_row_number: int,
    score_row_number: int,
    key: tuple[str, str],
    final_score: float,
    final_score_source: str,
    comparison_columns: Iterable[str],
    authoritative_drug: Mapping[str, Any],
    authoritative_target: Mapping[str, Any],
    target_mapping_status: str,
    target_mapping_key: str,
) -> dict[str, Any]:
    candidate = dict(selected_row)
    selected_snapshot_columns = (
        "target_uniprot",
        "target_gene",
        "generic_name",
        "display_name",
        *authoritative_drug.keys(),
        *authoritative_target.keys(),
    )
    for column in dict.fromkeys(selected_snapshot_columns):
        candidate[f"matched_inactive_selected_{column}"] = selected_row.get(column)
    candidate.update(authoritative_drug)
    candidate.update(authoritative_target)
    for column in comparison_columns:
        candidate[column] = score_row.get(column)
    candidate.update(
        {
            "pdb_id": key[0],
            "ligand_base": _clean_text(selected_row.get("ligand_base")),
            "final_score": final_score,
            "final_score_source": final_score_source,
            "spd_binding_label": 0,
            "spd_binding_label_policy_version": "spd_binding_censor_aware_v1",
            "combined_activity_ml_label": 0,
            "source_objective": "spd_binding_activity",
            "label_source": _first_nonmissing(
                selected_row,
                "evidence_source_names",
                "evidence_sources",
                "representative_source",
            ),
            "source_family": _first_nonmissing(
                selected_row, "evidence_source_names", "representative_source"
            ),
            "upstream_source": _first_nonmissing(
                selected_row, "evidence_sources", "representative_source"
            ),
            "source_label_policy": LABEL_STATUS,
            "assay_type": selected_row.get("activity_type"),
            "endpoint_type": selected_row.get("activity_type"),
            "training_allowed": _training_allowed(selected_row),
            "matched_inactive_pair_key": _key_text(key),
            "matched_inactive_score_source": final_score_source,
            "matched_inactive_selected_row": selected_row_number,
            "matched_inactive_score_row": score_row_number,
            "matched_inactive_selected_target_uniprot": selected_row.get(
                "target_uniprot"
            ),
            "matched_inactive_selected_target_gene": selected_row.get("target_gene"),
            "matched_inactive_selected_generic_name": selected_row.get("generic_name"),
            "matched_inactive_selected_display_name": selected_row.get("display_name"),
            "matched_inactive_target_mapping_status": target_mapping_status,
            "matched_inactive_target_mapping_key": target_mapping_key,
        }
    )
    if _is_missing(candidate.get("activity_value_nm")):
        candidate["activity_value_nm"] = selected_row.get("activity_nM")
    if _is_missing(candidate.get("activity_units")):
        candidate["activity_units"] = "nM"
    return candidate


def _write_csv(
    rows: pd.DataFrame,
    path: Path,
    *,
    columns: Iterable[str] | None = None,
) -> None:
    if columns is not None:
        rows = rows.reindex(columns=list(columns))
    rows.to_csv(path, index=False)


def merge_matched_inactive_scores(
    base_table_path: str | Path,
    selected_pairs_path: str | Path,
    candidate_score_path: str | Path,
    out_path: str | Path,
    *,
    audit_dir: str | Path | None = None,
    allow_base_pair_collisions: bool = False,
) -> dict[str, Any]:
    """Append strictly measured, fully scored matched inactives to an ML table."""

    base_path = Path(base_table_path)
    selected_path = Path(selected_pairs_path)
    score_path = Path(candidate_score_path)
    output_path = Path(out_path)
    audit_path = (
        Path(audit_dir)
        if audit_dir is not None
        else output_path.with_name(f"{output_path.stem}_audit")
    )
    audit_path.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    base = _read_csv_strict(base_path, name="base table")
    selected = _read_csv_strict(selected_path, name="selected matched-inactive pairs")
    scores = _read_csv_strict(score_path, name="candidate score table")
    _require_columns(base, (*PAIR_COLUMNS, "drug_id", "target_id"), name="base table")
    _require_columns(
        selected,
        (*PAIR_COLUMNS, "projected_label", "projected_label_status"),
        name="selected matched-inactive pairs",
    )
    required_score_columns = {
        *PAIR_COLUMNS,
        CONSENSUS_PREFIX,
        SCORCH_PREFIX,
        *(
            f"{prefix}{suffix}"
            for prefix in (CONSENSUS_PREFIX, SCORCH_PREFIX)
            for suffix in PROVENANCE_SUFFIXES
        ),
    }
    _require_columns(scores, required_score_columns, name="candidate score table")
    comparison_columns = _comparison_columns(scores)
    ligand_authoritative_columns = _ligand_authoritative_columns(base)
    target_authoritative_columns = _target_authoritative_columns(base)
    target_alias_authoritative_columns = [
        column
        for column in target_authoritative_columns
        if column != "structure_quality"
    ]
    (
        authoritative_drug_map,
        authoritative_drug_status,
        base_drug_map_audit,
    ) = _build_base_identity_map(
        base,
        key_column="ligand_base",
        identity_column="drug_id",
        value_columns=ligand_authoritative_columns,
        normalize_key=_normalize_ligand_base,
        map_type="ligand_to_drug",
    )
    (
        authoritative_target_map,
        authoritative_target_status,
        base_target_map_audit,
    ) = _build_base_identity_map(
        base,
        key_column="pdb_id",
        identity_column="target_id",
        value_columns=target_authoritative_columns,
        normalize_key=_normalize_pdb_id,
        map_type="pdb_to_target",
    )
    (
        authoritative_target_alias_map,
        authoritative_target_alias_status,
        base_target_alias_map_audit,
    ) = _build_base_target_alias_map(
        base,
        value_columns=target_alias_authoritative_columns,
    )

    base_keys = [
        _pair_key(pdb_id, ligand_base)
        for pdb_id, ligand_base in base.loc[:, list(PAIR_COLUMNS)].itertuples(
            index=False, name=None
        )
    ]
    selected_keys = [
        _pair_key(pdb_id, ligand_base)
        for pdb_id, ligand_base in selected.loc[:, list(PAIR_COLUMNS)].itertuples(
            index=False, name=None
        )
    ]
    score_keys = [
        _pair_key(pdb_id, ligand_base)
        for pdb_id, ligand_base in scores.loc[:, list(PAIR_COLUMNS)].itertuples(
            index=False, name=None
        )
    ]
    base_key_counts = Counter(key for key in base_keys if key is not None)
    duplicate_selected = _selected_duplicate_keys(selected_keys)
    rejected_score_keys, score_conflicts = _duplicate_score_audit(
        scores, score_keys, comparison_columns
    )
    score_positions = {
        key: position
        for position, key in enumerate(score_keys)
        if key is not None and key not in rejected_score_keys
    }

    candidate_rows: list[dict[str, Any]] = []
    coverage_rows: list[dict[str, Any]] = []
    collision_rows: list[dict[str, Any]] = []
    identity_discrepancy_rows: list[dict[str, Any]] = []
    for selected_position, (key, selected_values) in enumerate(
        zip(selected_keys, selected.to_dict(orient="records"), strict=False)
    ):
        reasons: list[str] = []
        if key is None:
            reasons.append("invalid_selected_pair_key")
        elif key in duplicate_selected:
            reasons.append("duplicate_selected_pair_key")
        label_valid = _label_contract_valid(selected_values)
        if not label_valid:
            reasons.append("invalid_selected_label_contract")

        drug_key = key[1] if key else ""
        target_key = key[0] if key else ""
        drug_mapping_status = authoritative_drug_status.get(
            drug_key, "missing_base_key"
        )
        pdb_target_mapping_status = authoritative_target_status.get(
            target_key, "missing_base_key"
        )
        target_gene_key = _normalize_target_identity(selected_values.get("target_gene"))
        target_gene_mapping_status = "not_attempted_pdb_mapping_available"
        target_mapping_key = target_key
        authoritative_drug = authoritative_drug_map.get(drug_key)
        authoritative_target = authoritative_target_map.get(target_key)
        if authoritative_target is not None:
            target_mapping_status = "authoritative_pdb_verified"
        else:
            target_mapping_status = pdb_target_mapping_status
            if pdb_target_mapping_status == "missing_base_key":
                target_gene_mapping_status = authoritative_target_alias_status.get(
                    target_gene_key, "missing_base_key"
                )
                verified_gene_target = authoritative_target_alias_map.get(
                    target_gene_key
                )
                if verified_gene_target is not None:
                    authoritative_target = {
                        column: verified_gene_target.get(column, pd.NA)
                        for column in target_authoritative_columns
                    }
                    target_mapping_status = "authoritative_target_gene_verified"
                    target_mapping_key = target_gene_key

        if authoritative_drug is None:
            reason = (
                "ambiguous_authoritative_base_drug_mapping"
                if drug_mapping_status == "ambiguous_authoritative_metadata"
                else "missing_authoritative_base_drug_mapping"
            )
            reasons.append(reason)
        if authoritative_target is None:
            if pdb_target_mapping_status == "ambiguous_authoritative_metadata":
                reason = "ambiguous_authoritative_base_target_mapping"
            elif target_gene_mapping_status == "ambiguous_authoritative_metadata":
                reason = "ambiguous_authoritative_base_target_gene_mapping"
            else:
                reason = "missing_authoritative_base_target_mapping"
            reasons.append(reason)

        base_count = int(base_key_counts.get(key, 0)) if key is not None else 0
        if base_count and not allow_base_pair_collisions:
            reasons.append("base_pair_collision")
            collision_rows.append(
                {
                    "selected_row": selected_position,
                    "pdb_id": key[0] if key else "",
                    "ligand_base": key[1] if key else "",
                    "normalized_pair_key": _key_text(key),
                    "base_pair_count": base_count,
                    "action": "skipped",
                }
            )

        score_position: int | None = None
        if key is not None and key in rejected_score_keys:
            reasons.append("duplicate_or_conflicting_score_key")
        elif key is not None:
            score_position = score_positions.get(key)
            if score_position is None:
                reasons.append("missing_candidate_score_row")

        final_score: float | None = None
        final_score_source = ""
        provenance_errors: list[str] = []
        score_values: dict[str, Any] = {}
        if score_position is not None:
            score_values = scores.iloc[score_position].to_dict()
            final_score, final_score_source, provenance_errors = _score_choice(
                score_values
            )
            if final_score is None:
                reasons.append("no_finite_comparison_score")
            if provenance_errors:
                reasons.append("incomplete_selected_score_provenance")

        appended = not reasons
        if (
            appended
            and key is not None
            and final_score is not None
            and score_position is not None
            and authoritative_drug is not None
            and authoritative_target is not None
        ):
            candidate_rows.append(
                _candidate_row(
                    selected_values,
                    score_values,
                    selected_row_number=selected_position,
                    score_row_number=score_position,
                    key=key,
                    final_score=final_score,
                    final_score_source=final_score_source,
                    comparison_columns=comparison_columns,
                    authoritative_drug=authoritative_drug,
                    authoritative_target=authoritative_target,
                    target_mapping_status=target_mapping_status,
                    target_mapping_key=target_mapping_key,
                )
            )

        authoritative_drug_id = (
            authoritative_drug.get("drug_id") if authoritative_drug else pd.NA
        )
        authoritative_target_id = (
            authoritative_target.get("target_id") if authoritative_target else pd.NA
        )
        identity_row: dict[str, Any] = {
            "selected_row": selected_position,
            "normalized_pair_key": _key_text(key),
            "normalized_ligand_base": drug_key,
            "normalized_pdb_id": target_key,
            "normalized_selected_target_gene": target_gene_key,
            "drug_mapping_status": drug_mapping_status,
            "pdb_target_mapping_status": pdb_target_mapping_status,
            "target_gene_mapping_status": target_gene_mapping_status,
            "target_mapping_status": target_mapping_status,
            "target_mapping_key": target_mapping_key,
            "selected_generic_name": selected_values.get("generic_name"),
            "selected_target_uniprot": selected_values.get("target_uniprot"),
            "selected_target_gene": selected_values.get("target_gene"),
            "authoritative_drug_id": authoritative_drug_id,
            "authoritative_target_id": authoritative_target_id,
            "generic_name_vs_drug_id": _identity_comparison(
                selected_values.get("generic_name"), authoritative_drug_id
            ),
            "target_gene_vs_target_id": _identity_comparison(
                selected_values.get("target_gene"), authoritative_target_id
            ),
            "target_uniprot_vs_target_id": _identity_comparison(
                selected_values.get("target_uniprot"), authoritative_target_id
            ),
            "appended": appended,
            "skip_reasons": ";".join(dict.fromkeys(reasons)),
        }
        for column in ligand_authoritative_columns:
            selected_value = selected_values.get(column)
            authoritative_value = (
                authoritative_drug.get(column) if authoritative_drug else pd.NA
            )
            identity_row[f"selected_{column}"] = selected_value
            identity_row[f"authoritative_{column}"] = authoritative_value
            identity_row[f"{column}_comparison"] = _identity_comparison(
                selected_value, authoritative_value
            )
        for column in target_authoritative_columns:
            selected_value = selected_values.get(column)
            authoritative_value = (
                authoritative_target.get(column) if authoritative_target else pd.NA
            )
            identity_row[f"selected_{column}"] = selected_value
            identity_row[f"authoritative_{column}"] = authoritative_value
            identity_row[f"{column}_comparison"] = _identity_comparison(
                selected_value, authoritative_value
            )
        identity_discrepancy_rows.append(identity_row)
        coverage_rows.append(
            {
                "selected_row": selected_position,
                "pdb_id": key[0]
                if key
                else _normalize_pdb_id(selected_values.get("pdb_id")),
                "ligand_base": key[1]
                if key
                else _normalize_ligand_base(selected_values.get("ligand_base")),
                "normalized_pair_key": _key_text(key),
                "label_contract_valid": label_valid,
                "base_pair_count": base_count,
                "drug_mapping_status": drug_mapping_status,
                "pdb_target_mapping_status": pdb_target_mapping_status,
                "target_gene_mapping_status": target_gene_mapping_status,
                "target_mapping_status": target_mapping_status,
                "target_mapping_key": target_mapping_key,
                "normalized_selected_target_gene": target_gene_key,
                "authoritative_drug_id": authoritative_drug_id,
                "authoritative_target_id": authoritative_target_id,
                "score_row": score_position if score_position is not None else "",
                CONSENSUS_PREFIX: score_values.get(CONSENSUS_PREFIX),
                SCORCH_PREFIX: score_values.get(SCORCH_PREFIX),
                "selected_final_score": final_score,
                "selected_final_score_source": final_score_source,
                "score_provenance_errors": ";".join(provenance_errors),
                "appended": appended,
                "skip_reasons": ";".join(dict.fromkeys(reasons)),
            }
        )

    candidate_frame = pd.DataFrame(candidate_rows)
    mandatory_columns = [
        "final_score",
        "final_score_source",
        "spd_binding_label",
        "combined_activity_ml_label",
        "activity_value_nm",
        "activity_units",
        *STANDARD_PROVENANCE_COLUMNS,
    ]
    selected_snapshot_columns = [
        f"matched_inactive_selected_{column}"
        for column in dict.fromkeys(
            [
                "target_uniprot",
                "target_gene",
                "generic_name",
                "display_name",
                *ligand_authoritative_columns,
                *target_authoritative_columns,
            ]
        )
    ]
    output_columns = list(
        dict.fromkeys(
            [
                *base.columns,
                *selected.columns,
                *comparison_columns,
                *mandatory_columns,
                *selected_snapshot_columns,
            ]
        )
    )
    base_extended = base.reindex(columns=output_columns)
    candidate_frame = candidate_frame.reindex(columns=output_columns)
    merged = pd.concat([base_extended, candidate_frame], ignore_index=True, sort=False)

    base_columns = list(base.columns)
    base_prefix = merged.iloc[: len(base)].loc[:, base_columns]
    base_frame_hash = _frame_hash(base, base_columns)
    prefix_frame_hash = _frame_hash(base_prefix, base_columns)
    base_nonmissing_hash = _nonmissing_hash(base, base_columns)
    prefix_nonmissing_hash = _nonmissing_hash(base_prefix, base_columns)
    appended_labels = pd.to_numeric(
        candidate_frame.get("spd_binding_label", pd.Series(dtype=float)),
        errors="coerce",
    )
    appended_combined_labels = pd.to_numeric(
        candidate_frame.get("combined_activity_ml_label", pd.Series(dtype=float)),
        errors="coerce",
    )
    candidate_keys = {
        _pair_key(row.get("pdb_id"), row.get("ligand_base")) for row in candidate_rows
    }
    base_key_set = {key for key in base_keys if key is not None}
    appended_drug_ids_authoritative = all(
        (
            (key := _pair_key(row.get("pdb_id"), row.get("ligand_base"))) is not None
            and key[1] in authoritative_drug_map
            and _canonical_value(row.get("drug_id"))
            == _canonical_value(authoritative_drug_map[key[1]]["drug_id"])
        )
        for row in candidate_rows
    )
    expected_candidate_targets = [
        _candidate_target_record(
            row,
            pdb_map=authoritative_target_map,
            alias_map=authoritative_target_alias_map,
        )
        for row in candidate_rows
    ]
    appended_target_ids_authoritative = all(
        expected is not None
        and _canonical_value(row.get("target_id"))
        == _canonical_value(expected.get("target_id"))
        for row, expected in zip(
            candidate_rows, expected_candidate_targets, strict=False
        )
    )
    appended_target_metadata_authoritative = all(
        expected is not None
        and all(
            _canonical_value(row.get(column)) == _canonical_value(expected.get(column))
            for column in target_authoritative_columns
        )
        for row, expected in zip(
            candidate_rows, expected_candidate_targets, strict=False
        )
    )
    target_gene_fallbacks_verified = all(
        (
            _clean_text(row.get("matched_inactive_target_mapping_status"))
            != "authoritative_target_gene_verified"
        )
        or (
            (
                mapping_key := _normalize_target_identity(
                    row.get("matched_inactive_target_mapping_key")
                )
            )
            != ""
            and mapping_key
            == _normalize_target_identity(
                row.get("matched_inactive_selected_target_gene")
            )
            and mapping_key in authoritative_target_alias_map
        )
        for row in candidate_rows
    )
    pdb_primary_preference_preserved = all(
        (
            _clean_text(row.get("matched_inactive_target_mapping_status"))
            != "authoritative_target_gene_verified"
        )
        or _normalize_pdb_id(row.get("pdb_id")) not in authoritative_target_map
        for row in candidate_rows
    )
    invariants = [
        ("base_rows_are_output_prefix", base_frame_hash == prefix_frame_hash),
        (
            "base_nonmissing_values_preserved",
            base_nonmissing_hash == prefix_nonmissing_hash,
        ),
        ("base_row_count_preserved", len(base_prefix) == len(base)),
        ("output_row_count_expected", len(merged) == len(base) + len(candidate_frame)),
        ("output_columns_unique", merged.columns.is_unique),
        ("all_appended_spd_labels_zero", appended_labels.eq(0).all()),
        ("all_appended_combined_labels_zero", appended_combined_labels.eq(0).all()),
        (
            "all_appended_final_scores_finite",
            pd.to_numeric(candidate_frame.get("final_score"), errors="coerce")
            .map(lambda value: pd.notna(value) and math.isfinite(float(value)))
            .all(),
        ),
        (
            "base_pair_collisions_absent_when_rejected",
            allow_base_pair_collisions or not bool(candidate_keys & base_key_set),
        ),
        ("candidate_pair_keys_unique", len(candidate_keys) == len(candidate_frame)),
        (
            "all_appended_drug_ids_authoritative",
            appended_drug_ids_authoritative,
        ),
        (
            "all_appended_target_ids_authoritative",
            appended_target_ids_authoritative,
        ),
        (
            "all_appended_target_metadata_authoritative",
            appended_target_metadata_authoritative,
        ),
        (
            "all_target_gene_fallbacks_verified",
            target_gene_fallbacks_verified,
        ),
        (
            "pdb_primary_preference_preserved",
            pdb_primary_preference_preserved,
        ),
        (
            "candidate_count_matches_coverage",
            len(candidate_frame) == sum(bool(row["appended"]) for row in coverage_rows),
        ),
    ]
    invariants_frame = pd.DataFrame(
        [{"invariant": name, "passed": bool(passed)} for name, passed in invariants]
    )
    failed = invariants_frame.loc[~invariants_frame["passed"], "invariant"].tolist()
    if failed:
        raise RuntimeError("merge invariants failed: " + ", ".join(failed))

    merged.to_csv(output_path, index=False)
    coverage = pd.DataFrame(coverage_rows)
    if coverage.empty:
        coverage = pd.DataFrame(
            columns=[
                "selected_row",
                "pdb_id",
                "ligand_base",
                "normalized_pair_key",
                "label_contract_valid",
                "base_pair_count",
                "drug_mapping_status",
                "pdb_target_mapping_status",
                "target_gene_mapping_status",
                "target_mapping_status",
                "target_mapping_key",
                "normalized_selected_target_gene",
                "authoritative_drug_id",
                "authoritative_target_id",
                "score_row",
                CONSENSUS_PREFIX,
                SCORCH_PREFIX,
                "selected_final_score",
                "selected_final_score_source",
                "score_provenance_errors",
                "appended",
                "skip_reasons",
            ]
        )
    skipped = coverage.loc[~coverage["appended"]].copy()
    conflicts = pd.DataFrame(score_conflicts)
    collisions = pd.DataFrame(collision_rows)
    identity_discrepancies = pd.DataFrame(identity_discrepancy_rows)
    base_drug_identity_map = pd.DataFrame(base_drug_map_audit)
    base_target_identity_map = pd.DataFrame(base_target_map_audit)
    base_target_alias_identity_map = pd.DataFrame(base_target_alias_map_audit)
    score_source_counts = (
        candidate_frame.get("final_score_source", pd.Series(dtype="object"))
        .fillna("missing")
        .astype(str)
        .value_counts(dropna=False)
        .rename_axis("final_score_source")
        .reset_index(name="candidate_rows")
    )
    label_contract = pd.DataFrame(
        [
            {
                "contract": "selected_projected_label_is_zero",
                "required_value": "0",
                "selected_rows_passing": int(
                    sum(
                        _finite_float(row.get("projected_label")) == 0.0
                        for row in selected.to_dict(orient="records")
                    )
                ),
                "appended_rows_passing": int(appended_labels.eq(0).sum()),
                "appended_rows": int(len(candidate_frame)),
                "passed": bool(appended_labels.eq(0).all()),
            },
            {
                "contract": "selected_projected_label_status",
                "required_value": LABEL_STATUS,
                "selected_rows_passing": int(
                    selected["projected_label_status"]
                    .fillna("")
                    .astype(str)
                    .str.strip()
                    .str.casefold()
                    .eq(LABEL_STATUS)
                    .sum()
                ),
                "appended_rows_passing": int(len(candidate_frame)),
                "appended_rows": int(len(candidate_frame)),
                "passed": True,
            },
        ]
    )

    audit_outputs = {
        "candidate_coverage": audit_path / "candidate_coverage.csv",
        "score_key_conflicts": audit_path / "score_key_conflicts.csv",
        "base_pair_collisions": audit_path / "base_pair_collisions.csv",
        "base_drug_identity_map": audit_path / "base_drug_identity_map.csv",
        "base_target_identity_map": audit_path / "base_target_identity_map.csv",
        "base_target_alias_identity_map": (
            audit_path / "base_target_alias_identity_map.csv"
        ),
        "selected_identity_discrepancies": (
            audit_path / "selected_identity_discrepancies.csv"
        ),
        "skipped_rows": audit_path / "skipped_rows.csv",
        "score_source_counts": audit_path / "score_source_counts.csv",
        "label_contract": audit_path / "label_contract.csv",
        "invariants": audit_path / "invariants.csv",
        "appended_candidate_rows": audit_path / "appended_candidate_rows.csv",
    }
    _write_csv(coverage, audit_outputs["candidate_coverage"])
    _write_csv(
        conflicts,
        audit_outputs["score_key_conflicts"],
        columns=(
            "score_row",
            "pdb_id",
            "ligand_base",
            "normalized_pair_key",
            "conflict_type",
            "duplicate_count",
            "payload_hash",
        ),
    )
    _write_csv(
        collisions,
        audit_outputs["base_pair_collisions"],
        columns=(
            "selected_row",
            "pdb_id",
            "ligand_base",
            "normalized_pair_key",
            "base_pair_count",
            "action",
        ),
    )
    _write_csv(base_drug_identity_map, audit_outputs["base_drug_identity_map"])
    _write_csv(base_target_identity_map, audit_outputs["base_target_identity_map"])
    _write_csv(
        base_target_alias_identity_map,
        audit_outputs["base_target_alias_identity_map"],
    )
    _write_csv(
        identity_discrepancies,
        audit_outputs["selected_identity_discrepancies"],
    )
    _write_csv(skipped, audit_outputs["skipped_rows"])
    _write_csv(score_source_counts, audit_outputs["score_source_counts"])
    _write_csv(label_contract, audit_outputs["label_contract"])
    _write_csv(invariants_frame, audit_outputs["invariants"])
    _write_csv(candidate_frame, audit_outputs["appended_candidate_rows"])

    source_counts = {
        str(row["final_score_source"]): int(row["candidate_rows"])
        for row in score_source_counts.to_dict(orient="records")
    }
    skip_counts = Counter(
        reason
        for value in skipped.get("skip_reasons", pd.Series(dtype="object")).fillna("")
        for reason in str(value).split(";")
        if reason
    )
    manifest: dict[str, Any] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "base_table": str(base_path),
        "selected_pairs": str(selected_path),
        "candidate_scores": str(score_path),
        "output": str(output_path),
        "audit_dir": str(audit_path),
        "join_key": ["normalized pdb_id", "normalized ligand_base"],
        "base_collision_policy": (
            "append_duplicate_without_overwrite"
            if allow_base_pair_collisions
            else "reject_candidate"
        ),
        "score_policy": (
            f"finite {SCORCH_PREFIX} first, else finite {CONSENSUS_PREFIX}; "
            "chosen comparison provenance is mandatory"
        ),
        "label_policy": (
            "append only selected rows with projected_label=0 and "
            f"projected_label_status={LABEL_STATUS}; write both model labels as 0"
        ),
        "identity_policy": (
            "drug_id and ligand metadata come only from a strict unique normalized "
            "ligand_base base map. Target identity uses the strict PDB map first; "
            "only a missing PDB key may use selected target_gene after that alias "
            "resolves uniquely against canonical base target_id/target_gene metadata. "
            "Selected target_uniprot is never an identity fallback."
        ),
        "counts": {
            "base_rows": int(len(base)),
            "selected_rows": int(len(selected)),
            "candidate_score_rows": int(len(scores)),
            "appended_rows": int(len(candidate_frame)),
            "output_rows": int(len(merged)),
            "skipped_selected_rows": int(len(skipped)),
            "invalid_score_key_rows": int(
                sum(
                    row.get("conflict_type") == "invalid_score_pair_key"
                    for row in score_conflicts
                )
            ),
            "duplicate_selected_keys": int(len(duplicate_selected)),
            "rejected_duplicate_score_keys": int(len(rejected_score_keys)),
            "base_pair_collision_rows": int(len(collisions)),
            "authoritative_target_gene_verified_rows": int(
                candidate_frame.get(
                    "matched_inactive_target_mapping_status",
                    pd.Series(dtype="object"),
                )
                .eq("authoritative_target_gene_verified")
                .sum()
            ),
            "authoritative_pdb_verified_rows": int(
                candidate_frame.get(
                    "matched_inactive_target_mapping_status",
                    pd.Series(dtype="object"),
                )
                .eq("authoritative_pdb_verified")
                .sum()
            ),
            "score_source_counts": source_counts,
            "skip_reason_counts": dict(sorted(skip_counts.items())),
        },
        "hashes": {
            "base_table_sha256": _sha256(base_path),
            "selected_pairs_sha256": _sha256(selected_path),
            "candidate_scores_sha256": _sha256(score_path),
            "base_frame_sha256": base_frame_hash,
            "base_nonmissing_sha256": base_nonmissing_hash,
            "appended_rows_sha256": _frame_hash(candidate_frame, output_columns),
            "output_sha256": _sha256(output_path),
            "audit_sha256": {
                name: _sha256(path) for name, path in audit_outputs.items()
            },
        },
        "invariants": {
            str(row["invariant"]): bool(row["passed"])
            for row in invariants_frame.to_dict(orient="records")
        },
        "outputs": {
            "merged_table": str(output_path),
            **{name: str(path) for name, path in audit_outputs.items()},
        },
    }
    manifest_path = audit_path / "merge_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return manifest


__all__ = ["merge_matched_inactive_scores"]
