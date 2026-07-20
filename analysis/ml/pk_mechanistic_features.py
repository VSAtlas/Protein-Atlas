from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
from typing import Any, Literal

import pandas as pd


SCHEMA_VERSION = "atlas_pk_mechanistic_features_v1"
SPD_FREE_CMAX_TARGET = "spd_free_cmax_um"
PK_CONTEXT_FREE_CMAX_TARGET = "pk_context_free_cmax_um"
TARGET_ENDPOINTS = (SPD_FREE_CMAX_TARGET, PK_CONTEXT_FREE_CMAX_TARGET)
TARGET_VALUE_COLUMNS = {
    SPD_FREE_CMAX_TARGET: "free_cmax_um",
    PK_CONTEXT_FREE_CMAX_TARGET: "pk_context_free_cmax_um",
}

EXCLUDED_GENERIC_ALIASES = (
    "value",
    "unit",
    "normalized_value",
    "normalized_unit",
    "solubility",
    "permeability",
    "clearance",
    "cl",
    "ka",
)

FieldKind = Literal["numeric", "text", "boolean"]


@dataclass(frozen=True)
class Candidate:
    column: str
    context_column: str | None = None
    required_contexts: tuple[str, ...] = ()
    route_column: str | None = None
    reference_route_column: str | None = None
    require_extravascular_iv_reference: bool = False


@dataclass(frozen=True)
class FieldSpec:
    output: str
    kind: FieldKind
    candidates: tuple[Candidate, ...]
    value_feature: bool = False


def _plain(*columns: str) -> tuple[Candidate, ...]:
    candidates: list[Candidate] = []
    seen: set[str] = set()
    for column in columns:
        for name in (column, f"pk_context_{column}"):
            if name not in seen:
                candidates.append(Candidate(name))
                seen.add(name)
    return tuple(candidates)


def _contextual(
    column: str,
    context_column: str,
    required_context: str,
) -> Candidate:
    return Candidate(
        column=column,
        context_column=context_column,
        required_contexts=(required_context,),
    )


def _field(
    suffix: str,
    kind: FieldKind,
    candidates: tuple[Candidate, ...],
    *,
    value_feature: bool = False,
) -> FieldSpec:
    return FieldSpec(
        output=f"pk_mech_{suffix}",
        kind=kind,
        candidates=candidates,
        value_feature=value_feature,
    )


def _clearance_specs(
    endpoint: str,
    *,
    value_aliases: tuple[str, ...],
    unit_aliases: tuple[str, ...],
) -> tuple[FieldSpec, ...]:
    context = f"{endpoint}_clearance"
    joined_prefix = f"pk_clearance_{endpoint}_"
    value_candidates = (
        *_plain(*value_aliases),
        Candidate(f"{joined_prefix}clearance_value"),
        _contextual("clearance_value", "measurement_context", context),
        _contextual(
            "pk_context_clearance_value",
            "pk_context_measurement_context",
            context,
        ),
    )
    unit_candidates = (
        *_plain(*unit_aliases),
        Candidate(f"{joined_prefix}clearance_unit"),
        _contextual("clearance_unit", "measurement_context", context),
        _contextual(
            "pk_context_clearance_unit",
            "pk_context_measurement_context",
            context,
        ),
    )
    context_candidates = (
        *_plain(f"{endpoint}_clearance_context"),
        Candidate(f"{joined_prefix}measurement_context"),
        _contextual("measurement_context", "measurement_context", context),
        _contextual(
            "pk_context_measurement_context",
            "pk_context_measurement_context",
            context,
        ),
    )
    training_candidates = (
        Candidate(f"{endpoint}_clearance_training_allowed"),
        Candidate(f"pk_context_{endpoint}_clearance_training_allowed"),
        Candidate(f"{joined_prefix}training_allowed"),
        _contextual(
            "pk_context_training_allowed",
            "pk_context_measurement_context",
            context,
        ),
    )
    return (
        _field(
            f"{endpoint}_clearance_value",
            "numeric",
            value_candidates,
            value_feature=True,
        ),
        _field(f"{endpoint}_clearance_unit", "text", unit_candidates),
        _field(f"{endpoint}_clearance_context", "text", context_candidates),
        _field(
            f"{endpoint}_clearance_training_allowed",
            "boolean",
            training_candidates,
        ),
    )


_BASE_FIELD_SPECS = (
    _field(
        "dose_normalized_value",
        "numeric",
        _plain("dose_normalized_value", "normalized_dose_value", "dose_value"),
        value_feature=True,
    ),
    _field(
        "dose_normalized_unit",
        "text",
        _plain("dose_normalized_unit", "normalized_dose_unit", "dose_unit"),
    ),
    _field(
        "dose_context",
        "text",
        _plain(
            "dose_normalization_context",
            "dose_context",
            "dose_context_type",
            "dose_provenance_type",
        ),
    ),
    _field(
        "route",
        "text",
        _plain("route", "route_of_administration", "dosage_route"),
    ),
    _field("formulation", "text", _plain("formulation", "dosage_form")),
    _field("regimen", "text", _plain("regimen", "dosing_regimen")),
    _field("steady_state", "boolean", _plain("steady_state", "is_steady_state")),
    _field(
        "pka_value",
        "numeric",
        _plain("pka_value", "pka"),
        value_feature=True,
    ),
    _field("pka_type", "text", _plain("pka_type", "pka_class", "ionization_type")),
    _field("pka_context", "text", _plain("pka_context", "ionization_context")),
    _field(
        "logd_7_4",
        "numeric",
        _plain("logd_7_4", "logd7_4", "logd_ph_7_4", "logd_ph7_4"),
        value_feature=True,
    ),
    _field(
        "logd_7_4_context",
        "text",
        _plain("logd_7_4_context", "logd_context"),
    ),
    _field(
        "aqueous_solubility_value",
        "numeric",
        _plain("aqueous_solubility_value", "water_solubility_value"),
        value_feature=True,
    ),
    _field(
        "aqueous_solubility_unit",
        "text",
        _plain("aqueous_solubility_unit", "water_solubility_unit"),
    ),
    _field(
        "aqueous_solubility_context",
        "text",
        _plain("aqueous_solubility_context", "water_solubility_context"),
    ),
    _field(
        "biorelevant_solubility_value",
        "numeric",
        _plain("biorelevant_solubility_value"),
        value_feature=True,
    ),
    _field(
        "biorelevant_solubility_unit",
        "text",
        _plain("biorelevant_solubility_unit"),
    ),
    _field(
        "biorelevant_solubility_medium",
        "text",
        _plain("biorelevant_solubility_medium", "biorelevant_medium"),
    ),
    _field(
        "caco2_papp_value",
        "numeric",
        _plain(
            "caco2_papp_value",
            "caco_2_papp_value",
            "caco2_apparent_permeability_value",
        ),
        value_feature=True,
    ),
    _field(
        "caco2_papp_unit",
        "text",
        _plain(
            "caco2_papp_unit",
            "caco_2_papp_unit",
            "caco2_apparent_permeability_unit",
        ),
    ),
    _field(
        "caco2_papp_direction",
        "text",
        _plain("caco2_papp_direction", "caco_2_papp_direction", "caco2_direction"),
    ),
    _field(
        "peff_value",
        "numeric",
        _plain("peff_value", "effective_permeability_value"),
        value_feature=True,
    ),
    _field(
        "peff_unit",
        "text",
        _plain("peff_unit", "effective_permeability_unit"),
    ),
    _field(
        "peff_context",
        "text",
        _plain("peff_context", "effective_permeability_context"),
    ),
    _field(
        "absorption_rate_ka_value",
        "numeric",
        _plain(
            "absorption_rate_ka_value",
            "absorption_rate_constant_value",
            "ka_absorption_value",
        ),
        value_feature=True,
    ),
    _field(
        "absorption_rate_ka_unit",
        "text",
        _plain(
            "absorption_rate_ka_unit",
            "absorption_rate_constant_unit",
            "ka_absorption_unit",
        ),
    ),
    _field(
        "absorption_rate_ka_context",
        "text",
        _plain("absorption_rate_ka_context"),
    ),
    _field(
        "tmax_value",
        "numeric",
        _plain("tmax_value", "t_max_value", "tmax"),
        value_feature=True,
    ),
    _field("tmax_unit", "text", _plain("tmax_unit", "t_max_unit")),
    _field("tmax_context", "text", _plain("tmax_context", "t_max_context")),
    _field(
        "fraction_unbound_plasma",
        "numeric",
        _plain("fraction_unbound_plasma", "fu_plasma", "fup"),
        value_feature=True,
    ),
    _field(
        "fraction_unbound_plasma_context",
        "text",
        _plain("fraction_unbound_plasma_context", "fu_plasma_context"),
    ),
    _field(
        "blood_to_plasma_ratio",
        "numeric",
        _plain("blood_to_plasma_ratio", "blood_plasma_ratio", "b_p_ratio"),
        value_feature=True,
    ),
    _field(
        "blood_to_plasma_context",
        "text",
        _plain("blood_to_plasma_context", "blood_plasma_ratio_context"),
    ),
)

_CLEARANCE_FIELD_SPECS = (
    *_clearance_specs(
        "intrinsic",
        value_aliases=("intrinsic_clearance_value", "clint_value", "cl_int_value"),
        unit_aliases=("intrinsic_clearance_unit", "clint_unit", "cl_int_unit"),
    ),
    *_clearance_specs(
        "hepatic",
        value_aliases=("hepatic_clearance_value", "clh_value"),
        unit_aliases=("hepatic_clearance_unit", "clh_unit"),
    ),
    *_clearance_specs(
        "renal",
        value_aliases=("renal_clearance_value", "clr_value"),
        unit_aliases=("renal_clearance_unit", "clr_unit"),
    ),
    *_clearance_specs(
        "systemic",
        value_aliases=("systemic_clearance_value",),
        unit_aliases=("systemic_clearance_unit",),
    ),
    *_clearance_specs(
        "apparent_oral",
        value_aliases=("apparent_oral_clearance_value",),
        unit_aliases=("apparent_oral_clearance_unit",),
    ),
    *_clearance_specs(
        "apparent_plasma",
        value_aliases=("apparent_plasma_clearance_value",),
        unit_aliases=("apparent_plasma_clearance_unit",),
    ),
    *_clearance_specs(
        "plasma",
        value_aliases=("plasma_clearance_value",),
        unit_aliases=("plasma_clearance_unit",),
    ),
    *_clearance_specs(
        "blood",
        value_aliases=("blood_clearance_value",),
        unit_aliases=("blood_clearance_unit",),
    ),
    *_clearance_specs(
        "total_body",
        value_aliases=("total_body_clearance_value",),
        unit_aliases=("total_body_clearance_unit",),
    ),
)


def _absolute_bioavailability_candidate(
    column: str,
    *,
    prefix: str = "",
) -> Candidate:
    return Candidate(
        column=f"{prefix}{column}",
        context_column=f"{prefix}measurement_context",
        required_contexts=("absolute_bioavailability",),
        route_column=f"{prefix}route",
        reference_route_column=f"{prefix}reference_route",
        require_extravascular_iv_reference=True,
    )


_BIOAVAILABILITY_FIELD_SPECS = (
    _field(
        "bioavailability_f_value",
        "numeric",
        (
            *_plain("bioavailability_f_value", "absolute_bioavailability_value"),
            Candidate("pk_bioavailability_absolute_bioavailability_value"),
            _absolute_bioavailability_candidate("bioavailability_value"),
            _absolute_bioavailability_candidate(
                "bioavailability_value",
                prefix="pk_context_",
            ),
        ),
        value_feature=True,
    ),
    _field(
        "bioavailability_f_unit",
        "text",
        (
            *_plain("bioavailability_f_unit", "absolute_bioavailability_unit"),
            Candidate("pk_bioavailability_absolute_bioavailability_unit"),
            _absolute_bioavailability_candidate("bioavailability_unit"),
            _absolute_bioavailability_candidate(
                "bioavailability_unit",
                prefix="pk_context_",
            ),
        ),
    ),
    _field(
        "bioavailability_f_context",
        "text",
        (
            *_plain("bioavailability_f_context", "absolute_bioavailability_context"),
            Candidate("pk_bioavailability_absolute_measurement_context"),
            _absolute_bioavailability_candidate("measurement_context"),
            _absolute_bioavailability_candidate(
                "measurement_context",
                prefix="pk_context_",
            ),
        ),
    ),
    _field(
        "bioavailability_f_reference_route",
        "text",
        (
            *_plain(
                "bioavailability_f_reference_route",
                "absolute_bioavailability_reference_route",
            ),
            Candidate("pk_bioavailability_absolute_reference_route"),
            _absolute_bioavailability_candidate("reference_route"),
            _absolute_bioavailability_candidate(
                "reference_route",
                prefix="pk_context_",
            ),
        ),
    ),
    _field(
        "bioavailability_f_training_allowed",
        "boolean",
        (
            Candidate("bioavailability_f_training_allowed"),
            Candidate("absolute_bioavailability_training_allowed"),
            Candidate("pk_context_bioavailability_f_training_allowed"),
            Candidate("pk_context_absolute_bioavailability_training_allowed"),
            Candidate("pk_bioavailability_absolute_training_allowed"),
            _absolute_bioavailability_candidate(
                "training_allowed",
                prefix="pk_context_",
            ),
        ),
    ),
)

_TAIL_FIELD_SPECS = (
    _field(
        "vss_value",
        "numeric",
        _plain("vss_value", "steady_state_volume_of_distribution_value"),
        value_feature=True,
    ),
    _field(
        "vss_unit",
        "text",
        _plain("vss_unit", "steady_state_volume_of_distribution_unit"),
    ),
    _field("vss_context", "text", _plain("vss_context", "steady_state_volume_context")),
    _field(
        "vd_value",
        "numeric",
        _plain("vd_value", "volume_of_distribution_value"),
        value_feature=True,
    ),
    _field("vd_unit", "text", _plain("vd_unit", "volume_of_distribution_unit")),
    _field(
        "vd_context",
        "text",
        _plain("vd_context", "volume_of_distribution_context"),
    ),
    _field(
        "half_life_value",
        "numeric",
        _plain(
            "half_life_value",
            "terminal_half_life_value",
            "elimination_half_life_value",
        ),
        value_feature=True,
    ),
    _field(
        "half_life_unit",
        "text",
        _plain(
            "half_life_unit",
            "terminal_half_life_unit",
            "elimination_half_life_unit",
        ),
    ),
    _field(
        "half_life_context",
        "text",
        _plain("half_life_context", "half_life_type"),
    ),
    _field("species", "text", _plain("species", "organism")),
    _field(
        "population",
        "text",
        _plain("population", "subject_population", "cohort"),
    ),
    _field(
        "body_weight_value",
        "numeric",
        _plain("body_weight_value", "subject_body_weight_value"),
        value_feature=True,
    ),
    _field(
        "body_weight_unit",
        "text",
        _plain("body_weight_unit", "subject_body_weight_unit"),
    ),
    _field(
        "body_weight_context",
        "text",
        _plain("body_weight_context", "subject_body_weight_context"),
    ),
    _field("source_name", "text", _plain("source_name")),
    _field("source_version", "text", _plain("source_version")),
    _field("source_record_id", "text", _plain("source_record_id")),
    _field("source_url", "text", _plain("source_url")),
    _field("study_id", "text", _plain("study_id")),
    _field("reference", "text", _plain("reference")),
    _field("measurement_context", "text", _plain("measurement_context")),
    _field("reference_route", "text", _plain("reference_route")),
    _field("extraction_method", "text", _plain("extraction_method")),
    _field("source_confidence", "text", _plain("source_confidence")),
    _field("context_status", "text", _plain("context_status")),
    _field(
        "training_allowed",
        "boolean",
        (
            Candidate("pk_source_training_allowed"),
            Candidate("pk_context_training_allowed"),
        ),
    ),
    _field("license_note", "text", _plain("license_note")),
    _field("source_missing_reason", "text", _plain("missing_reason")),
)

PK_MECHANISTIC_FIELD_SPECS = (
    *_BASE_FIELD_SPECS,
    *_CLEARANCE_FIELD_SPECS,
    *_BIOAVAILABILITY_FIELD_SPECS,
    *_TAIL_FIELD_SPECS,
)
PK_MECHANISTIC_COLUMNS = tuple(spec.output for spec in PK_MECHANISTIC_FIELD_SPECS)
PK_MECHANISTIC_VALUE_COLUMNS = tuple(
    spec.output for spec in PK_MECHANISTIC_FIELD_SPECS if spec.value_feature
)

ADMINISTRATION_ALIGNMENT_FIELDS = {
    "dose": "pk_mech_dose_normalized_value",
    "route": "pk_mech_route",
    "formulation": "pk_mech_formulation",
    "regimen": "pk_mech_regimen",
    "steady_state": "pk_mech_steady_state",
}
PK_ENDPOINT_ALIGNMENT_COLUMNS = (
    "pk_mech_target_endpoint",
    "pk_mech_target_value_column",
    "pk_mech_target_source_name",
    "pk_mech_target_study_id",
    "pk_mech_target_scenario_id",
    *tuple(
        column
        for field in ADMINISTRATION_ALIGNMENT_FIELDS
        for column in (
            f"pk_mech_{field}_endpoint_alignment_status",
            f"pk_mech_{field}_endpoint_alignment_reason",
            f"pk_mech_{field}_training_allowed",
        )
    ),
)
PK_MECHANISTIC_AUDIT_COLUMNS = (
    "pk_mech_materialization_status",
    "pk_mech_missing_reasons",
)
PK_MECHANISTIC_OUTPUT_COLUMNS = (
    *PK_MECHANISTIC_COLUMNS,
    *PK_ENDPOINT_ALIGNMENT_COLUMNS,
    *PK_MECHANISTIC_AUDIT_COLUMNS,
)

_MISSING_TEXT = {"", "nan", "none", "null", "unknown", "n/a", "na", "not available"}
_TRUE_TEXT = {"1", "true", "yes", "y"}
_FALSE_TEXT = {"0", "false", "no", "n"}


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False
    try:
        if bool(missing):
            return True
    except (TypeError, ValueError):
        pass
    return isinstance(value, str) and value.strip().casefold() in _MISSING_TEXT


def _coerce(value: Any, kind: FieldKind) -> Any:
    if _is_missing(value):
        return pd.NA
    if kind == "numeric":
        numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        return pd.NA if pd.isna(numeric) else numeric
    if kind == "boolean":
        if isinstance(value, bool):
            return value
        normalized = str(value).strip().casefold()
        if normalized in _TRUE_TEXT:
            return True
        if normalized in _FALSE_TEXT:
            return False
        return pd.NA
    return str(value).strip()


def _canonical(value: Any, kind: FieldKind) -> tuple[str, Any]:
    if kind == "numeric":
        return ("numeric", float(value))
    if kind == "boolean":
        return ("boolean", bool(value))
    return ("text", str(value).strip().casefold())


def _normalized_text(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series("", index=frame.index, dtype="object")
    return frame[column].fillna("").astype(str).str.strip().str.casefold()


def _candidate_condition(frame: pd.DataFrame, candidate: Candidate) -> pd.Series:
    condition = pd.Series(True, index=frame.index, dtype=bool)
    if candidate.required_contexts:
        if candidate.context_column is None or candidate.context_column not in frame:
            return pd.Series(False, index=frame.index, dtype=bool)
        condition &= _normalized_text(frame, candidate.context_column).isin(
            candidate.required_contexts
        )
    if candidate.require_extravascular_iv_reference:
        if (
            candidate.route_column is None
            or candidate.reference_route_column is None
            or candidate.route_column not in frame
            or candidate.reference_route_column not in frame
        ):
            return pd.Series(False, index=frame.index, dtype=bool)
        route = _normalized_text(frame, candidate.route_column)
        reference_route = _normalized_text(frame, candidate.reference_route_column)
        condition &= route.ne("") & route.ne("intravenous")
        condition &= reference_route.eq("intravenous")
    return condition


def _merge_output_columns(
    frame: pd.DataFrame,
    columns: dict[str, Any],
) -> pd.DataFrame:
    positional: dict[str, Any] = {}
    for column, values in columns.items():
        if isinstance(values, pd.Series):
            positional[column] = values.to_numpy(copy=False)
        elif isinstance(values, (str, bytes)) or not hasattr(values, "__len__"):
            positional[column] = [values] * len(frame)
        else:
            positional[column] = list(values)
    updates = pd.DataFrame(positional, index=frame.index)
    existing = [column for column in columns if column in frame]
    if existing:
        replacement = set(existing)
        base_columns = {
            column: (
                updates[column].array if column in replacement else frame[column].array
            )
            for column in frame
        }
        out = pd.DataFrame(base_columns, index=frame.index)
    else:
        out = frame.copy()
    added = [column for column in columns if column not in out]
    if added:
        out = pd.concat([out, updates[added]], axis=1)
    return out


def _resolve_field(
    frame: pd.DataFrame,
    spec: FieldSpec,
) -> tuple[pd.Series, list[str], list[str | None], dict[str, Any]]:
    existing = (
        frame[spec.output]
        if spec.output in frame
        else pd.Series(pd.NA, index=frame.index, dtype="object")
    )
    existing_coerced = existing.map(lambda value: _coerce(value, spec.kind))
    output = existing_coerced.copy()
    candidate_rows: list[tuple[Candidate, pd.Series, pd.Series, pd.Series]] = []
    present_columns: list[str] = []
    for candidate in spec.candidates:
        if candidate.column not in frame:
            continue
        present_columns.append(candidate.column)
        raw = frame[candidate.column]
        raw_present = raw.map(lambda value: not _is_missing(value))
        condition = _candidate_condition(frame, candidate)
        coerced = raw.map(lambda value: _coerce(value, spec.kind))
        parsed = coerced.map(lambda value: not _is_missing(value))
        candidate_rows.append(
            (
                candidate,
                coerced,
                raw_present & condition & parsed,
                raw_present & ~condition,
            )
        )

    reasons: list[str] = []
    selected_sources: list[str | None] = []
    selected_counts: Counter[str] = Counter()
    conflicts = 0
    context_rejected = 0
    invalid_values = 0
    for position in range(len(frame)):
        values: list[tuple[Candidate, Any]] = []
        any_raw = False
        any_invalid = False
        any_context_rejected = False
        for candidate, coerced, valid, rejected in candidate_rows:
            raw_value = frame[candidate.column].iloc[position]
            raw_present = not _is_missing(raw_value)
            any_raw |= raw_present
            any_context_rejected |= bool(rejected.iloc[position])
            if raw_present and not rejected.iloc[position] and not valid.iloc[position]:
                any_invalid = True
            if valid.iloc[position]:
                values.append((candidate, coerced.iloc[position]))

        unique: dict[tuple[str, Any], tuple[Candidate, Any]] = {}
        for candidate, value in values:
            unique.setdefault(_canonical(value, spec.kind), (candidate, value))
        resolved_source: str | None = None
        resolved_value: Any = pd.NA
        if len(unique) == 1:
            resolved_source = next(iter(unique.values()))[0].column
            resolved_value = next(iter(unique.values()))[1]
            selected_counts[resolved_source] += 1
        elif len(unique) > 1:
            conflicts += 1

        if _is_missing(output.iloc[position]) and not _is_missing(resolved_value):
            output.iloc[position] = resolved_value
        selected_sources.append(resolved_source)
        if not _is_missing(output.iloc[position]):
            reasons.append("")
        elif len(unique) > 1:
            reasons.append("conflicting_recognized_values")
        elif any_context_rejected:
            reasons.append("incompatible_or_missing_endpoint_context")
            context_rejected += 1
        elif any_invalid:
            reasons.append("invalid_source_value")
            invalid_values += 1
        elif present_columns and any_raw:
            reasons.append("recognized_source_value_unusable")
        elif present_columns:
            reasons.append("recognized_source_value_missing")
        else:
            reasons.append("no_recognized_source_column")

    if spec.kind == "numeric":
        output = pd.to_numeric(output, errors="coerce").astype("Float64")
    elif spec.kind == "boolean":
        output = output.astype("boolean")
    else:
        output = output.astype("string")

    report = {
        "column": spec.output,
        "kind": spec.kind,
        "recognized_source_columns_present": present_columns,
        "selected_source_counts": dict(sorted(selected_counts.items())),
        "conflicting_rows": conflicts,
        "context_rejected_rows": context_rejected,
        "invalid_source_value_rows": invalid_values,
    }
    return output, reasons, selected_sources, report


def _text_at(frame: pd.DataFrame, position: int, *columns: str) -> str:
    for column in columns:
        if column not in frame:
            continue
        value = frame[column].iloc[position]
        if not _is_missing(value):
            return str(value).strip()
    return ""


def _truthy_at(frame: pd.DataFrame, position: int, *columns: str) -> bool:
    value = _text_at(frame, position, *columns).casefold()
    return value in _TRUE_TEXT


def _exact_pk_context_identity(frame: pd.DataFrame, position: int) -> bool:
    model_inchikey = _text_at(
        frame,
        position,
        "inchikey",
        "drug_inchikey",
        "ligand_inchikey",
    ).upper()
    context_inchikey = _text_at(frame, position, "pk_context_inchikey").upper()
    return bool(
        model_inchikey
        and context_inchikey
        and model_inchikey == context_inchikey
    )


def _value_present(frame: pd.DataFrame, position: int, column: str) -> bool:
    return column in frame and not _is_missing(frame[column].iloc[position])


def _target_provenance(
    frame: pd.DataFrame,
    position: int,
    target_endpoint: str,
) -> tuple[str, str, str]:
    if target_endpoint == SPD_FREE_CMAX_TARGET:
        source = "SPD"
        context_source = _text_at(frame, position, "pk_context_source_name")
        if context_source.casefold() == "spd":
            return (
                source,
                _text_at(frame, position, "pk_context_study_id"),
                _text_at(
                    frame,
                    position,
                    "pk_context_pk_context_id",
                    "pk_context_context_id",
                ),
            )
        row_source = _text_at(frame, position, "source_name")
        scenario = ""
        study = ""
        if row_source.casefold() == "spd":
            scenario = _text_at(frame, position, "pk_context_id", "source_record_id")
            study = _text_at(frame, position, "study_id")
        return source, study, scenario
    return (
        _text_at(frame, position, "pk_context_source_name"),
        _text_at(frame, position, "pk_context_study_id"),
        _text_at(
            frame,
            position,
            "pk_context_pk_context_id",
            "pk_context_context_id",
        ),
    )


def _alignment(
    frame: pd.DataFrame,
    position: int,
    *,
    field: str,
    output_column: str,
    selected_source: str | None,
    target_endpoint: str,
    require_exact_identity: bool = False,
) -> tuple[str, str, bool]:
    if not _value_present(frame, position, output_column):
        return (
            "feature_missing",
            "mechanistic administration field is unavailable",
            False,
        )
    target_column = TARGET_VALUE_COLUMNS[target_endpoint]
    if not _value_present(frame, position, target_column):
        return "target_missing", f"{target_column} is unavailable", False
    if selected_source is None:
        return (
            "source_scenario_unresolved",
            "materialized value lacks a unique recognized source column",
            False,
        )

    from_pk_context = selected_source.startswith("pk_context_")
    if target_endpoint == SPD_FREE_CMAX_TARGET:
        if from_pk_context:
            source = _text_at(frame, position, "pk_context_source_name")
            scenario = _text_at(
                frame,
                position,
                "pk_context_pk_context_id",
                "pk_context_context_id",
            )
            record = _text_at(
                frame,
                position,
                "pk_context_study_id",
                "pk_context_source_record_id",
            )
            target_value = pd.to_numeric(
                pd.Series([frame[target_column].iloc[position]]), errors="coerce"
            ).iloc[0]
            context_value = pd.to_numeric(
                pd.Series([frame["pk_context_free_cmax_um"].iloc[position]]),
                errors="coerce",
            ).iloc[0]
            same_value = (
                pd.notna(target_value)
                and pd.notna(context_value)
                and abs(float(target_value) - float(context_value))
                <= 1e-9 * max(abs(float(target_value)), 1.0)
            )
            if (
                source.casefold() != "spd"
                or not scenario
                or not record
                or not same_value
            ):
                return (
                    "mismatched_external_context_for_spd_target",
                    "pk_context source/scenario does not exactly match the SPD free-Cmax target",
                    False,
                )
            if require_exact_identity and not _exact_pk_context_identity(
                frame, position
            ):
                return (
                    "non_exact_pk_context_join",
                    "training requires a full-InChIKey PK context join",
                    False,
                )
            dose_context = _text_at(
                frame, position, "pk_context_dose_context_type"
            ).casefold()
            if field == "dose" and dose_context != "cmax_study_matched":
                return (
                    "dose_not_cmax_scenario_matched",
                    "dose_context_type is not cmax_study_matched",
                    False,
                )
            rights = _truthy_at(frame, position, "pk_context_training_allowed")
            if not rights:
                return (
                    "aligned_same_cmax_scenario",
                    "SPD endpoint context aligns, but source training_allowed is not affirmative",
                    False,
                )
            return (
                "aligned_same_cmax_scenario",
                "recovered SPD source, record, and Cmax scenario align",
                True,
            )
        source = _text_at(frame, position, "source_name")
        scenario = _text_at(frame, position, "pk_context_id", "source_record_id")
        measurement = _text_at(frame, position, "measurement_context").casefold()
        dose_context = _text_at(frame, position, "dose_context_type").casefold()
        cmax_scenario = (
            measurement == "observed_cmax" or dose_context == "cmax_study_matched"
        )
        if source.casefold() != "spd" or not scenario or not cmax_scenario:
            return (
                "source_study_or_cmax_scenario_unverified",
                "SPD source and same-row Cmax scenario identity are not jointly verified",
                False,
            )
        # Backward-compatible direct calls may identify an explicitly
        # SPD-namespaced field. The production resolver does not emit these
        # names, so ordinary dose/route columns cannot borrow generic rights.
        if not selected_source.casefold().startswith("spd_"):
            return (
                "pk_specific_training_rights_missing",
                "generic training_allowed cannot authorize an unnamespaced PK feature",
                False,
            )
        rights = _truthy_at(frame, position, "training_allowed")
    else:
        if not from_pk_context:
            return (
                "mismatched_noncontext_feature_for_pk_context_target",
                "administration value was not selected from the pk_context target scenario",
                False,
            )
        source = _text_at(frame, position, "pk_context_source_name")
        scenario = _text_at(
            frame,
            position,
            "pk_context_pk_context_id",
            "pk_context_context_id",
        )
        record = _text_at(
            frame,
            position,
            "pk_context_study_id",
            "pk_context_source_record_id",
        )
        dose_context = _text_at(
            frame,
            position,
            "pk_context_dose_context_type",
        ).casefold()
        if not source or not scenario or not record:
            return (
                "source_study_or_cmax_scenario_unverified",
                "pk_context source, study/record, and scenario identity are required",
                False,
            )
        if require_exact_identity and not _exact_pk_context_identity(
            frame, position
        ):
            return (
                "non_exact_pk_context_join",
                "training requires a full-InChIKey PK context join",
                False,
            )
        if field == "dose" and dose_context != "cmax_study_matched":
            return (
                "dose_not_cmax_scenario_matched",
                "dose_context_type is not cmax_study_matched",
                False,
            )
        rights = _truthy_at(frame, position, "pk_context_training_allowed")

    if not rights:
        return (
            "aligned_same_cmax_scenario",
            "endpoint context aligns, but source training_allowed is not affirmative",
            False,
        )
    return (
        "aligned_same_cmax_scenario",
        "source, study/record, and Cmax scenario align with the selected target",
        True,
    )


def _add_endpoint_alignment(
    frame: pd.DataFrame,
    *,
    target_endpoint: str,
    selection_sources: dict[str, list[str | None]],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    target_sources: list[str | Any] = []
    target_studies: list[str | Any] = []
    target_scenarios: list[str | Any] = []
    alignment_counts: dict[str, Counter[str]] = {
        field: Counter() for field in ADMINISTRATION_ALIGNMENT_FIELDS
    }
    training_counts: dict[str, int] = {
        field: 0 for field in ADMINISTRATION_ALIGNMENT_FIELDS
    }
    results: dict[str, dict[str, list[Any]]] = {
        field: {"status": [], "reason": [], "allowed": []}
        for field in ADMINISTRATION_ALIGNMENT_FIELDS
    }
    for position in range(len(frame)):
        source, study, scenario = _target_provenance(
            frame,
            position,
            target_endpoint,
        )
        target_sources.append(source or pd.NA)
        target_studies.append(study or pd.NA)
        target_scenarios.append(scenario or pd.NA)
        for field, output_column in ADMINISTRATION_ALIGNMENT_FIELDS.items():
            selected = selection_sources.get(output_column, [None] * len(frame))[
                position
            ]
            status, reason, allowed = _alignment(
                frame,
                position,
                field=field,
                output_column=output_column,
                selected_source=selected,
                target_endpoint=target_endpoint,
                require_exact_identity=True,
            )
            results[field]["status"].append(status)
            results[field]["reason"].append(reason)
            results[field]["allowed"].append(allowed)
            alignment_counts[field][status] += 1
            training_counts[field] += int(allowed)

    alignment_columns: dict[str, Any] = {
        "pk_mech_target_endpoint": [target_endpoint] * len(frame),
        "pk_mech_target_value_column": [TARGET_VALUE_COLUMNS[target_endpoint]]
        * len(frame),
        "pk_mech_target_source_name": target_sources,
        "pk_mech_target_study_id": target_studies,
        "pk_mech_target_scenario_id": target_scenarios,
    }
    for field in ADMINISTRATION_ALIGNMENT_FIELDS:
        alignment_columns[f"pk_mech_{field}_endpoint_alignment_status"] = results[
            field
        ]["status"]
        alignment_columns[f"pk_mech_{field}_endpoint_alignment_reason"] = results[
            field
        ]["reason"]
        alignment_columns[f"pk_mech_{field}_training_allowed"] = results[field][
            "allowed"
        ]
    aligned = _merge_output_columns(frame, alignment_columns)
    return aligned, {
        "target_endpoint": target_endpoint,
        "target_value_column": TARGET_VALUE_COLUMNS[target_endpoint],
        "alignment_status_counts": {
            field: dict(sorted(counts.items()))
            for field, counts in alignment_counts.items()
        },
        "training_allowed_rows": training_counts,
    }


def materialize_pk_mechanistic_features(
    frame: pd.DataFrame,
    *,
    target_endpoint: str = SPD_FREE_CMAX_TARGET,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Copy recognized PK evidence into explicit, context-preserving columns."""

    if target_endpoint not in TARGET_ENDPOINTS:
        raise ValueError(
            f"target_endpoint must be one of {TARGET_ENDPOINTS}; got {target_endpoint!r}"
        )
    if not frame.columns.is_unique:
        raise ValueError("Input table has duplicate column names")

    reason_by_column: dict[str, list[str]] = {}
    selection_sources: dict[str, list[str | None]] = {}
    field_reports: list[dict[str, Any]] = []
    materialized_values: dict[str, pd.Series] = {}
    for spec in PK_MECHANISTIC_FIELD_SPECS:
        values, reasons, sources, report = _resolve_field(frame, spec)
        materialized_values[spec.output] = values
        reason_by_column[spec.output] = reasons
        selection_sources[spec.output] = sources
        field_reports.append(report)

    out = _merge_output_columns(frame, materialized_values)
    value_coverage = out[list(PK_MECHANISTIC_VALUE_COLUMNS)].apply(
        lambda column: column.map(lambda value: not _is_missing(value))
    )
    value_counts = value_coverage.sum(axis=1)
    materialization_status = pd.Series("partial", index=out.index, dtype="object")
    materialization_status.loc[value_counts.eq(0)] = "no_mechanistic_values"
    materialization_status.loc[value_counts.eq(len(PK_MECHANISTIC_VALUE_COLUMNS)),] = (
        "complete"
    )

    compact_missing: list[Any] = []
    for position in range(len(out)):
        issues = {
            column: reason_by_column[column][position]
            for column in PK_MECHANISTIC_VALUE_COLUMNS
            if reason_by_column[column][position]
            not in {"", "no_recognized_source_column"}
        }
        if not issues and value_counts.iloc[position] == 0:
            issues = {"all": "no_recognized_mechanistic_source_columns"}
        compact_missing.append(
            json.dumps(issues, sort_keys=True, separators=(",", ":"))
            if issues
            else pd.NA
        )
    out = _merge_output_columns(
        out,
        {
            "pk_mech_materialization_status": materialization_status,
            "pk_mech_missing_reasons": compact_missing,
        },
    )
    out, alignment_manifest = _add_endpoint_alignment(
        out,
        target_endpoint=target_endpoint,
        selection_sources=selection_sources,
    )
    for report in field_reports:
        column = report["column"]
        non_missing = int(out[column].map(lambda value: not _is_missing(value)).sum())
        reason_counts = Counter(reason for reason in reason_by_column[column] if reason)
        report.update(
            {
                "non_missing_rows": non_missing,
                "missing_rows": int(len(out) - non_missing),
                "coverage_fraction": (
                    round(non_missing / len(out), 8) if len(out) else 0.0
                ),
                "missing_reason_counts": dict(sorted(reason_counts.items())),
            }
        )

    status_counts = (
        out["pk_mech_materialization_status"].value_counts(dropna=False).to_dict()
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "row_count": int(len(out)),
        "materialized_columns": list(PK_MECHANISTIC_COLUMNS),
        "endpoint_alignment_columns": list(PK_ENDPOINT_ALIGNMENT_COLUMNS),
        "audit_columns": list(PK_MECHANISTIC_AUDIT_COLUMNS),
        "policy": {
            "operation": "copy_only",
            "unit_conversion": "none",
            "imputation": "none",
            "conflicting_recognized_values": "leave_missing",
            "generic_aliases_excluded": list(EXCLUDED_GENERIC_ALIASES),
            "source_training_allowed": (
                "pk_mech_training_allowed is copied from source and never inferred"
            ),
            "administration_training_allowed": (
                "derived from endpoint alignment and affirmative source rights"
            ),
            "default_target": SPD_FREE_CMAX_TARGET,
            "administration_alignment": (
                "training requires the same source/study-or-record/Cmax scenario "
                "as the selected target plus affirmative source training_allowed"
            ),
            "mismatched_context": "audit_or_sensitivity_only",
        },
        "field_coverage": field_reports,
        "fully_missing_columns": [
            report["column"]
            for report in field_reports
            if report["non_missing_rows"] == 0
        ],
        "materialization_status_counts": {
            str(key): int(value) for key, value in sorted(status_counts.items())
        },
        "endpoint_alignment": alignment_manifest,
        "invariants": {
            "row_count_preserved": len(out) == len(frame),
            "input_columns_preserved_as_prefix": list(out.columns[: len(frame.columns)])
            == list(frame.columns),
            "output_columns_unique": out.columns.is_unique,
            "no_unit_conversion_or_imputation": True,
        },
    }
    return out, manifest


__all__ = [
    "ADMINISTRATION_ALIGNMENT_FIELDS",
    "PK_CONTEXT_FREE_CMAX_TARGET",
    "PK_ENDPOINT_ALIGNMENT_COLUMNS",
    "PK_MECHANISTIC_AUDIT_COLUMNS",
    "PK_MECHANISTIC_COLUMNS",
    "PK_MECHANISTIC_FIELD_SPECS",
    "PK_MECHANISTIC_OUTPUT_COLUMNS",
    "PK_MECHANISTIC_VALUE_COLUMNS",
    "SCHEMA_VERSION",
    "SPD_FREE_CMAX_TARGET",
    "TARGET_ENDPOINTS",
    "TARGET_VALUE_COLUMNS",
    "materialize_pk_mechanistic_features",
]
