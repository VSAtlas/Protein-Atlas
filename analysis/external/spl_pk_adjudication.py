from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Collection, Mapping, Sequence
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

from analysis.external.spl_pk_semantics import (
    SEMANTIC_COLUMNS,
    enrich_spl_semantic_evidence,
)


SCHEMA_VERSION = "atlas_spl_pk_contextual_adjudication_v4"
MANIFEST_NAME = "spl_pk_adjudication_manifest.json"

ALL_ADJUDICATED_NAME = "all_adjudicated.csv"
ACCEPTED_CLEARANCE_NAME = "accepted_clearance_context.csv"
ACCEPTED_BIOAVAILABILITY_NAME = "accepted_absolute_bioavailability_context.csv"
MAX_DOSE_SENSITIVITY_NAME = "max_recommended_dose_sensitivity.csv"
CONFLICT_SUMMARY_NAME = "conflict_summary.csv"
REVIEW_TEMPLATE_NAME = "review_decisions_template.csv"

ADJUDICATION_STATUSES = ("accepted", "ambiguous", "rejected")
REVIEW_DECISIONS = ("approved", "rejected", "deferred")
_CLEARANCE_ENDPOINTS = {
    "renal_clearance",
    "apparent_oral_clearance",
    "apparent_plasma_clearance",
    "total_body_clearance",
    "systemic_clearance",
    "plasma_clearance",
    "blood_clearance",
}
_MAX_DOSE_ENDPOINTS = {
    "maximum_daily_dose",
    "maximum_single_dose",
    "not_to_exceed_dose",
}
_PARENT_ANALYTE_STATUSES = {"context_label_match", "single_label_analyte"}
_AMBIGUOUS_ANALYTE_STATUSES = {"cache_drug_id_fallback", "not_detected"}
_TARGET_CANDIDATE_TYPES = {
    "clearance",
    "absolute_bioavailability",
    "maximum_labeled_dose",
}

_REQUIRED_INPUT_COLUMNS = {
    "candidate_id",
    "candidate_type",
    "drug_id",
    "normalized_value",
    "normalized_unit",
    "normalization_status",
    "measurement_context",
    "analyte",
    "analyte_status",
    "route",
    "population",
    "section",
    "source_excerpt",
}

_REVIEW_DECISION_ALIASES = {
    "accept": "approved",
    "accepted": "approved",
    "approve": "approved",
    "approved": "approved",
    "approve_for_training": "approved",
    "defer": "deferred",
    "deferred": "deferred",
    "needs_review": "deferred",
    "reject": "rejected",
    "rejected": "rejected",
}

_ADJUDICATION_COLUMNS = (
    *SEMANTIC_COLUMNS,
    "input_row_number",
    "adjudication_status",
    "adjudication_reasons",
    "evidence_tier",
    "machine_adjudication_status",
    "machine_adjudication_reasons",
    "machine_evidence_tier",
    "endpoint_semantics",
    "unit_basis",
    "human_context",
    "parent_analyte_confirmed",
    "extravascular_route",
    "iv_reference_confirmed",
    "semantic_evidence_key",
    "comparison_key",
    "canonical_candidate_id",
    "duplicate_group_size",
    "collapsed_duplicate_count",
    "duplicate_candidate_ids",
    "source_record_ids",
    "spl_versions",
    "cache_file_sha256s",
    "section_text_sha256s",
    "source_lineage_json",
    "conflict_id",
    "conflicting_normalized_values",
    "review_decision",
    "review_application_status",
    "reviewer",
    "reviewed_at",
    "review_notes",
    "reviewed_value",
    "reviewed_unit",
    "reviewed_value_qualifier",
    "reviewed_analyte",
    "reviewed_endpoint_semantics",
    "reviewed_reference_basis",
    "reviewed_extravascular_route",
    "reviewed_formulation",
    "reviewed_population",
    "review_decision_source",
    "review_decision_file_sha256",
    "review_decision_row_json",
    "review_approved",
    "sensitivity_only",
    "training_allowed",
    "adjudication_schema_version",
)

_REVIEW_TEMPLATE_COLUMNS = (
    "canonical_candidate_id",
    "candidate_type",
    "machine_adjudication_status",
    "eligible_for_approval",
    "sensitivity_only",
    "drug_id",
    "normalized_value",
    "normalized_unit",
    "endpoint_semantics",
    "semantic_clause",
    "route_evidence_text",
    "route_evidence_json",
    "extravascular_route",
    "iv_reference_confirmed",
    "source_record_ids",
    "spl_versions",
    "section_text_sha256s",
    "source_lineage_json",
    "review_decision",
    "reviewer",
    "reviewed_at",
    "review_notes",
    "reviewed_value",
    "reviewed_unit",
    "reviewed_value_qualifier",
    "reviewed_analyte",
    "reviewed_endpoint_semantics",
    "reviewed_reference_basis",
    "reviewed_extravascular_route",
    "reviewed_formulation",
    "reviewed_population",
)

_CONFLICT_COLUMNS = (
    "conflict_id",
    "comparison_key",
    "candidate_type",
    "drug_id",
    "endpoint_semantics",
    "analyte",
    "normalized_unit",
    "unit_basis",
    "route",
    "formulation",
    "population",
    "regimen",
    "steady_state",
    "normalized_values",
    "semantic_evidence_keys",
    "candidate_ids",
    "candidate_count",
    "source_lineage_json",
    "adjudication_reason",
)

_LINEAGE_FIELDS = (
    "candidate_id",
    "source_name",
    "source_record_id",
    "source_set_id",
    "source_document_id",
    "spl_version",
    "effective_time",
    "source_url",
    "cache_file",
    "cache_file_sha256",
    "section",
    "section_item_index",
    "section_text_sha256",
    "match_start",
    "match_end",
    "source_excerpt",
)

_NON_HUMAN_RE = re.compile(
    r"\b(?:animal|animals|canine|dog|dogs|mice|mouse|monkey|monkeys|rabbit|"
    r"rabbits|rat|rats|rodent|rodents)\b",
    re.IGNORECASE,
)
_HUMAN_RE = re.compile(
    r"\b(?:human|humans|adult|adults|subject|subjects|patient|patients|"
    r"volunteer|volunteers|participant|participants)\b",
    re.IGNORECASE,
)
_PEDIATRIC_RE = re.compile(
    r"\b(?:pediatric|paediatric|child|children|adolescent|adolescents|infant|"
    r"infants|neonate|neonates)\b",
    re.IGNORECASE,
)
_ADULT_RE = re.compile(r"\b(?:adult|adults|geriatric|elderly)\b", re.IGNORECASE)
_ABSOLUTE_BA_RE = re.compile(
    r"\babsolute\s+(?:oral\s+)?bioavailability\b|\bAbsBio\b",
    re.IGNORECASE,
)
_RELATIVE_BA_RE = re.compile(
    r"\brelative\s+(?:oral\s+)?bioavailability\b|"
    r"\bbioavailability\s+(?:was\s+)?relative\b",
    re.IGNORECASE,
)
_IV_RE = re.compile(r"\b(?:intravenous(?:ly)?|IV)\b", re.IGNORECASE)
_RANGE_RE = re.compile(
    r"\d\s*(?:-|\u2013|\u2014|to)\s*\d|\+/-|\u00b1",
    re.IGNORECASE,
)
_MAXIMUM_RE = re.compile(
    r"\b(?:maximum|maximal|max\.)\b|\bnot\s+(?:to\s+)?exceed\b",
    re.IGNORECASE,
)
_RECOMMENDED_DOSE_RE = re.compile(
    r"\b(?:recommended|recommendation|should\s+not\s+exceed|"
    r"must\s+not\s+exceed|do\s+not\s+exceed|not\s+(?:to\s+)?exceed)\b",
    re.IGNORECASE,
)
_DAILY_RE = re.compile(r"\b(?:daily|per\s+day|each\s+day)\b", re.IGNORECASE)
_SINGLE_DOSE_RE = re.compile(
    r"\b(?:single|per)\s+(?:treatment\s+)?dose\b",
    re.IGNORECASE,
)
_NOT_TO_EXCEED_RE = re.compile(r"\bnot\s+(?:to\s+)?exceed\b", re.IGNORECASE)
_STUDY_DOSE_RE = re.compile(
    r"\b(?:maximum\s+tolerated\s+dose|MTD|highest\s+(?:dose\s+)?studied|"
    r"highest\s+studied\s+dose|pharmacokinetic\s+(?:study|trial)|"
    r"PK\s+(?:study|trial)|dose[- ]escalation)\b",
    re.IGNORECASE,
)

_ROUTE_PATTERNS = (
    ("intravenous", _IV_RE),
    ("oral", re.compile(r"\b(?:oral(?:ly)?|PO)\b", re.IGNORECASE)),
    ("subcutaneous", re.compile(r"\bsubcutaneous(?:ly)?\b", re.IGNORECASE)),
    ("intramuscular", re.compile(r"\bintramuscular(?:ly)?\b", re.IGNORECASE)),
    ("inhalation", re.compile(r"\b(?:inhaled|inhalation)\b", re.IGNORECASE)),
    ("transdermal", re.compile(r"\btransdermal(?:ly)?\b", re.IGNORECASE)),
    ("sublingual", re.compile(r"\bsublingual(?:ly)?\b", re.IGNORECASE)),
    ("buccal", re.compile(r"\bbuccal(?:ly)?\b", re.IGNORECASE)),
    ("intranasal", re.compile(r"\b(?:intranasal(?:ly)?|nasal)\b", re.IGNORECASE)),
    ("vaginal", re.compile(r"\bvaginal(?:ly)?\b", re.IGNORECASE)),
    ("rectal", re.compile(r"\brectal(?:ly)?\b", re.IGNORECASE)),
    ("topical", re.compile(r"\btopical(?:ly)?\b", re.IGNORECASE)),
)


def _clean(value: Any) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def _canonical_text(value: Any) -> str:
    return _clean(value).casefold()


def _unique(values: Collection[str]) -> list[str]:
    return sorted({value for value in values if value})


def _append_reason(reasons: list[str], reason: str) -> None:
    if reason and reason not in reasons:
        reasons.append(reason)


def _decimal_text(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    try:
        number = Decimal(text)
    except InvalidOperation:
        return ""
    if not number.is_finite():
        return ""
    normalized = format(number.normalize(), "f")
    return "0" if normalized in {"-0", ""} else normalized


def _numeric_value(row: Mapping[str, Any]) -> float | None:
    text = _clean(row.get("normalized_value"))
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def _unit_basis(unit: Any) -> str:
    normalized = _canonical_text(unit).replace("\u00b2", "2")
    if re.search(r"/(?:kg|kilograms?)(?:/|$)", normalized):
        return "body_weight_normalized"
    if re.search(r"/(?:m2|square\s*meters?)(?:/|$)", normalized):
        return "body_surface_area_normalized"
    return "absolute"


def _same_clause_route_tokens(row: Mapping[str, Any]) -> list[str]:
    clause = _clean(row.get("semantic_clause"))
    return [name for name, pattern in _ROUTE_PATTERNS if pattern.search(clause)]


def _route_tokens(row: Mapping[str, Any]) -> list[str]:
    if _clean(row.get("candidate_type")) in {"clearance", "absolute_bioavailability"}:
        return _same_clause_route_tokens(row)
    text = " ".join((_clean(row.get("route")), _clean(row.get("source_excerpt"))))
    return [name for name, pattern in _ROUTE_PATTERNS if pattern.search(text)]


def _normalization_reasons(
    row: Mapping[str, Any],
    *,
    allow_zero: bool = False,
) -> list[str]:
    reasons: list[str] = []
    raw = _clean(row.get("raw_value"))
    exclusion = _canonical_text(row.get("exclusion_reason"))
    if (
        _RANGE_RE.search(raw)
        or "numeric_range_requires_adjudication" in exclusion
        or "summary_statistic_requires_adjudication" in exclusion
    ):
        reasons.append("range_or_variability_not_allowed")
    if _clean(row.get("normalization_status")) != "normalized_unambiguous":
        reasons.append("normalization_not_unambiguous")
    value = _numeric_value(row)
    if value is None:
        reasons.append("normalized_value_missing_or_invalid")
    elif value < 0 or (value == 0 and not allow_zero):
        reasons.append("normalized_value_not_positive")
    if not _clean(row.get("normalized_unit")):
        reasons.append("normalized_unit_missing")
    return reasons


def _common_context(
    row: Mapping[str, Any],
    *,
    require_parent_analyte: bool,
    allow_zero: bool = False,
) -> tuple[list[str], list[str], str, bool]:
    rejected = _normalization_reasons(row, allow_zero=allow_zero)
    ambiguous: list[str] = []
    if _clean(row.get("semantic_context_status")) != "verified":
        reason = (
            _clean(row.get("semantic_context_reasons")) or "semantic_context_missing"
        )
        _append_reason(ambiguous, f"same_clause_semantics_not_verified:{reason}")
    population_context = " ".join(
        (
            _clean(row.get("population")),
            _clean(row.get("source_excerpt")),
            _clean(row.get("exclusion_reason")),
        )
    )
    if _NON_HUMAN_RE.search(
        population_context
    ) or "non_human_context" in _canonical_text(row.get("exclusion_reason")):
        _append_reason(rejected, "non_human_context")
        human_context = "non_human"
    elif _HUMAN_RE.search(population_context):
        human_context = "explicit_human"
    else:
        human_context = "human_spl_label_context"

    analyte_status = _clean(row.get("analyte_status"))
    parent_confirmed = bool(
        _clean(row.get("adjudicated_analyte"))
        and _canonical_text(row.get("adjudicated_analyte"))
        == _canonical_text(row.get("drug_id"))
    )
    if require_parent_analyte and not parent_confirmed:
        if analyte_status in _AMBIGUOUS_ANALYTE_STATUSES or not parent_confirmed:
            _append_reason(ambiguous, "parent_analyte_not_confirmed")
        else:
            _append_reason(rejected, "parent_or_single_label_analyte_required")
    return rejected, ambiguous, human_context, parent_confirmed


def _clearance_decision(row: Mapping[str, Any]) -> dict[str, Any]:
    rejected, ambiguous, human_context, parent_confirmed = _common_context(
        row,
        require_parent_analyte=True,
    )
    endpoint = _clean(row.get("semantic_measurement_context")) or _clean(
        row.get("measurement_context")
    )
    if endpoint not in _CLEARANCE_ENDPOINTS:
        _append_reason(rejected, "explicit_clearance_endpoint_required")
    unit = _clean(row.get("normalized_unit"))
    if unit not in {"L/h", "L/h/kg", "L/h/m2"}:
        _append_reason(rejected, "unsupported_clearance_unit")
    if _clean(row.get("section")) not in {"pharmacokinetics", "clinical_pharmacology"}:
        _append_reason(rejected, "clearance_requires_pk_label_section")

    routes = _same_clause_route_tokens(row)
    if not routes:
        _append_reason(ambiguous, "clearance_route_not_confirmed_in_value_clause")
    if endpoint == "apparent_oral_clearance" and "oral" not in routes:
        if routes:
            _append_reason(rejected, "apparent_oral_clearance_route_not_oral")
        else:
            _append_reason(
                ambiguous,
                "apparent_oral_clearance_route_not_confirmed_in_value_clause",
            )
    return {
        "rejected": rejected,
        "ambiguous": ambiguous,
        "endpoint_semantics": endpoint or "clearance_unspecified",
        "unit_basis": _unit_basis(unit),
        "human_context": human_context,
        "parent_analyte_confirmed": parent_confirmed,
        "extravascular_route": "",
        "iv_reference_confirmed": "",
    }


def _bioavailability_decision(row: Mapping[str, Any]) -> dict[str, Any]:
    rejected, ambiguous, human_context, parent_confirmed = _common_context(
        row,
        require_parent_analyte=True,
        allow_zero=True,
    )
    clause = _clean(row.get("semantic_clause"))
    absolute_reported = bool(_ABSOLUTE_BA_RE.search(clause))
    if not absolute_reported:
        _append_reason(
            ambiguous,
            "explicit_absolute_bioavailability_not_confirmed_in_value_clause",
        )
    if _RELATIVE_BA_RE.search(clause):
        if absolute_reported:
            _append_reason(
                ambiguous,
                "relative_bioavailability_also_present_requires_review",
            )
        else:
            _append_reason(rejected, "relative_bioavailability_not_allowed")

    value = _numeric_value(row)
    if _clean(row.get("normalized_unit")) != "%":
        _append_reason(rejected, "bioavailability_percent_unit_required")
    if value is not None and not 0 <= value <= 100:
        _append_reason(rejected, "bioavailability_outside_zero_to_100_percent")
    if _clean(row.get("measurement_context")) != "absolute_bioavailability":
        _append_reason(rejected, "absolute_bioavailability_endpoint_required")
    if _clean(row.get("section")) not in {"pharmacokinetics", "clinical_pharmacology"}:
        _append_reason(rejected, "bioavailability_requires_pk_label_section")

    routes = _same_clause_route_tokens(row)
    extravascular = [route for route in routes if route != "intravenous"]
    iv_confirmed = "intravenous" in routes
    if not extravascular:
        _append_reason(ambiguous, "extravascular_route_unconfirmed_in_value_clause")
    if not iv_confirmed:
        _append_reason(ambiguous, "intravenous_reference_unconfirmed_in_value_clause")
    route_semantics = (
        "+".join(extravascular) if extravascular else "extravascular_unconfirmed"
    )
    return {
        "rejected": rejected,
        "ambiguous": ambiguous,
        "endpoint_semantics": f"absolute_bioavailability:{route_semantics}_vs_intravenous",
        "unit_basis": "absolute",
        "human_context": human_context,
        "parent_analyte_confirmed": parent_confirmed,
        "extravascular_route": ";".join(extravascular),
        "iv_reference_confirmed": str(iv_confirmed).lower(),
    }


def _maximum_dose_decision(row: Mapping[str, Any]) -> dict[str, Any]:
    rejected, ambiguous, human_context, parent_confirmed = _common_context(
        row,
        require_parent_analyte=False,
    )
    endpoint = _clean(row.get("semantic_measurement_context")) or _clean(
        row.get("measurement_context")
    )
    clause = _clean(row.get("semantic_clause"))
    if _clean(row.get("section")) != "dosage_and_administration":
        _append_reason(rejected, "maximum_dose_requires_dosage_and_administration")
    if endpoint not in _MAX_DOSE_ENDPOINTS:
        _append_reason(rejected, "maximum_dose_interval_unspecified")
    pediatric_context = bool(_PEDIATRIC_RE.search(clause))
    adult_context = bool(_ADULT_RE.search(clause))
    if pediatric_context and adult_context:
        _append_reason(
            ambiguous,
            "mixed_adult_pediatric_maximum_requires_structured_review",
        )
        human_context = "mixed_adult_pediatric"
    elif pediatric_context:
        _append_reason(rejected, "pediatric_maximum_dose_not_allowed")
        human_context = "pediatric"
    elif adult_context:
        human_context = "adult"
    else:
        _append_reason(
            ambiguous, "explicit_adult_context_not_confirmed_in_value_clause"
        )
        human_context = "adult_context_unconfirmed"
    if _STUDY_DOSE_RE.search(clause):
        _append_reason(rejected, "pk_study_highest_studied_or_mtd_not_allowed")
    if not _MAXIMUM_RE.search(clause):
        _append_reason(
            ambiguous,
            "explicit_maximum_dose_language_not_confirmed_in_value_clause",
        )
    if not _RECOMMENDED_DOSE_RE.search(clause):
        _append_reason(
            ambiguous,
            "explicit_recommended_dose_context_not_confirmed_in_value_clause",
        )
    if endpoint == "maximum_daily_dose" and not _DAILY_RE.search(clause):
        _append_reason(ambiguous, "daily_interval_not_confirmed")
    elif endpoint == "maximum_single_dose" and not _SINGLE_DOSE_RE.search(clause):
        _append_reason(ambiguous, "single_dose_interval_not_confirmed")
    elif endpoint == "not_to_exceed_dose" and not _NOT_TO_EXCEED_RE.search(clause):
        _append_reason(ambiguous, "not_to_exceed_language_not_confirmed")
    unit = _clean(row.get("normalized_unit"))
    if not re.fullmatch(r"mg(?:/(?:kg|m2|day|h|dose)){0,2}", unit):
        _append_reason(rejected, "unsupported_normalized_maximum_dose_unit")
    return {
        "rejected": rejected,
        "ambiguous": ambiguous,
        "endpoint_semantics": endpoint or "maximum_labeled_dose_unspecified_interval",
        "unit_basis": _unit_basis(unit),
        "human_context": human_context,
        "parent_analyte_confirmed": parent_confirmed,
        "extravascular_route": "",
        "iv_reference_confirmed": "",
    }


def _decision(row: Mapping[str, Any]) -> dict[str, Any]:
    candidate_type = _clean(row.get("candidate_type"))
    if candidate_type == "clearance":
        decision = _clearance_decision(row)
    elif candidate_type == "absolute_bioavailability":
        decision = _bioavailability_decision(row)
    elif candidate_type == "maximum_labeled_dose":
        decision = _maximum_dose_decision(row)
    else:
        decision = {
            "rejected": ["candidate_type_not_in_contextual_adjudication_scope"],
            "ambiguous": [],
            "endpoint_semantics": _clean(row.get("measurement_context"))
            or "unsupported",
            "unit_basis": _unit_basis(row.get("normalized_unit")),
            "human_context": "not_adjudicated",
            "parent_analyte_confirmed": False,
            "extravascular_route": "",
            "iv_reference_confirmed": "",
        }

    rejected = list(decision.pop("rejected"))
    ambiguous = list(decision.pop("ambiguous"))
    if rejected:
        status = "rejected"
        reasons = rejected + [reason for reason in ambiguous if reason not in rejected]
        tier = "excluded"
    elif ambiguous:
        status = "ambiguous"
        reasons = ambiguous
        tier = "tier_2_context_incomplete"
    else:
        status = "accepted"
        reasons = ["strict_contextual_criteria_met"]
        tier = (
            "sensitivity_only_explicit_label_maximum"
            if candidate_type == "maximum_labeled_dose"
            else "tier_1_explicit_label_context"
        )
    decision.update(
        {
            "adjudication_status": status,
            "adjudication_reasons": ";".join(reasons),
            "evidence_tier": tier,
        }
    )
    return decision


def _hash_payload(prefix: str, payload: Mapping[str, Any]) -> str:
    serialized = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    )
    return f"{prefix}_{hashlib.sha256(serialized.encode('utf-8')).hexdigest()[:24]}"


def _key_context(row: Mapping[str, Any]) -> dict[str, str]:
    routes = ";".join(_unique(_route_tokens(row)))
    return {
        "drug_id": _canonical_text(row.get("drug_id")),
        "candidate_type": _clean(row.get("candidate_type")),
        "endpoint_semantics": _clean(row.get("endpoint_semantics")),
        "analyte": _canonical_text(row.get("analyte")),
        "normalized_unit": _clean(row.get("normalized_unit")),
        "unit_basis": _clean(row.get("unit_basis")),
        "route": routes,
        "formulation": _canonical_text(row.get("formulation")),
        "population": _canonical_text(row.get("population")),
        "regimen": _canonical_text(row.get("regimen")),
        "steady_state": _canonical_text(row.get("steady_state")),
        "extravascular_route": _clean(row.get("extravascular_route")),
        "iv_reference_confirmed": _clean(row.get("iv_reference_confirmed")),
    }


def _assign_evidence_keys(row: dict[str, Any]) -> None:
    comparison = _key_context(row)
    semantic = dict(comparison)
    semantic["normalized_value"] = _decimal_text(row.get("normalized_value"))
    if not semantic["normalized_value"]:
        semantic["raw_value"] = _canonical_text(row.get("raw_value"))
        semantic["raw_unit"] = _canonical_text(row.get("raw_unit"))
    row["comparison_key"] = _hash_payload("cmp", comparison)
    row["semantic_evidence_key"] = _hash_payload("sem", semantic)


def _row_identity(row: Mapping[str, Any]) -> tuple[str, int]:
    candidate_id = _clean(row.get("candidate_id"))
    try:
        row_number = int(_clean(row.get("input_row_number")))
    except ValueError:
        row_number = 0
    return candidate_id, row_number


def _lineage(rows: Sequence[Mapping[str, Any]]) -> str:
    lineage = [
        {field: _clean(row.get(field)) for field in _LINEAGE_FIELDS}
        | {"input_row_number": _clean(row.get("input_row_number"))}
        for row in sorted(rows, key=_row_identity)
    ]
    return json.dumps(lineage, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _joined_values(rows: Sequence[Mapping[str, Any]], field: str) -> str:
    return ";".join(_unique([_clean(row.get(field)) for row in rows]))


def _apply_duplicate_metadata(rows: list[dict[str, Any]]) -> None:
    groups: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[_clean(row.get("semantic_evidence_key"))].append(row)
    for group in groups.values():
        ordered = sorted(group, key=_row_identity)
        candidate_ids = _joined_values(ordered, "candidate_id")
        canonical_id = _clean(ordered[0].get("candidate_id")) or (
            f"input_row_{ordered[0]['input_row_number']}"
        )
        lineage = _lineage(ordered)
        aggregate = {
            "canonical_candidate_id": canonical_id,
            "duplicate_group_size": len(ordered),
            "collapsed_duplicate_count": max(0, len(ordered) - 1),
            "duplicate_candidate_ids": candidate_ids,
            "source_record_ids": _joined_values(ordered, "source_record_id"),
            "spl_versions": _joined_values(ordered, "spl_version"),
            "cache_file_sha256s": _joined_values(ordered, "cache_file_sha256"),
            "section_text_sha256s": _joined_values(ordered, "section_text_sha256"),
            "source_lineage_json": lineage,
        }
        for row in group:
            row.update(aggregate)


def _mark_conflicts(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if (
            _clean(row.get("candidate_type")) in _TARGET_CANDIDATE_TYPES
            and _clean(row.get("adjudication_status")) != "rejected"
            and _decimal_text(row.get("normalized_value"))
        ):
            groups[_clean(row.get("comparison_key"))].append(row)

    summaries: list[dict[str, Any]] = []
    for comparison_key, group in sorted(groups.items()):
        values = _unique([_decimal_text(row.get("normalized_value")) for row in group])
        if len(values) <= 1:
            continue
        conflict_id = _hash_payload(
            "conflict",
            {"comparison_key": comparison_key, "normalized_values": values},
        )
        for row in group:
            reasons = [
                reason
                for reason in _clean(row.get("adjudication_reasons")).split(";")
                if reason and reason != "strict_contextual_criteria_met"
            ]
            _append_reason(reasons, "comparable_normalized_value_conflict")
            row.update(
                {
                    "adjudication_status": "ambiguous",
                    "adjudication_reasons": ";".join(reasons),
                    "evidence_tier": "tier_2_conflicting_label_evidence",
                    "conflict_id": conflict_id,
                    "conflicting_normalized_values": ";".join(values),
                }
            )
        first = sorted(group, key=_row_identity)[0]
        summaries.append(
            {
                "conflict_id": conflict_id,
                "comparison_key": comparison_key,
                "candidate_type": _clean(first.get("candidate_type")),
                "drug_id": _clean(first.get("drug_id")),
                "endpoint_semantics": _clean(first.get("endpoint_semantics")),
                "analyte": _clean(first.get("analyte")),
                "normalized_unit": _clean(first.get("normalized_unit")),
                "unit_basis": _clean(first.get("unit_basis")),
                "route": _joined_values(group, "route"),
                "formulation": _joined_values(group, "formulation"),
                "population": _joined_values(group, "population"),
                "regimen": _joined_values(group, "regimen"),
                "steady_state": _joined_values(group, "steady_state"),
                "normalized_values": ";".join(values),
                "semantic_evidence_keys": _joined_values(
                    group, "semantic_evidence_key"
                ),
                "candidate_ids": _joined_values(group, "candidate_id"),
                "candidate_count": len(group),
                "source_lineage_json": _lineage(group),
                "adjudication_reason": "comparable_normalized_value_conflict",
            }
        )
    return summaries


def _initialize_review_metadata(rows: Sequence[dict[str, Any]]) -> None:
    for row in rows:
        sensitivity_only = _clean(row.get("candidate_type")) == "maximum_labeled_dose"
        row.update(
            {
                "machine_adjudication_status": _clean(row.get("adjudication_status")),
                "machine_adjudication_reasons": _clean(row.get("adjudication_reasons")),
                "machine_evidence_tier": _clean(row.get("evidence_tier")),
                "review_decision": "",
                "review_application_status": "review_not_provided",
                "reviewer": "",
                "reviewed_at": "",
                "review_notes": "",
                "reviewed_value": "",
                "reviewed_unit": "",
                "reviewed_value_qualifier": "",
                "reviewed_analyte": "",
                "reviewed_endpoint_semantics": "",
                "reviewed_reference_basis": "",
                "reviewed_extravascular_route": "",
                "reviewed_formulation": "",
                "reviewed_population": "",
                "review_decision_source": "",
                "review_decision_file_sha256": "",
                "review_decision_row_json": "",
                "review_approved": "false",
                "sensitivity_only": _bool_text(sensitivity_only),
                "training_allowed": "false",
            }
        )


def _normalize_review_decision(value: Any) -> str:
    decision = re.sub(r"[^a-z0-9]+", "_", _canonical_text(value)).strip("_")
    if not decision:
        return ""
    normalized = _REVIEW_DECISION_ALIASES.get(decision)
    if normalized is None:
        allowed = ", ".join(REVIEW_DECISIONS)
        raise ValueError(
            f"unsupported review_decision {value!r}; expected one of: {allowed}"
        )
    return normalized


def _read_review_decisions(
    path: Path,
) -> tuple[dict[str, dict[str, str]], str]:
    if not path.is_file():
        raise ValueError(f"review decision CSV not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"review decision CSV has no header: {path}")
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError("review decision CSV contains duplicate column names")
        missing = {"canonical_candidate_id", "review_decision"}.difference(fieldnames)
        if missing:
            raise ValueError(
                "review decision CSV missing required columns: "
                + ", ".join(sorted(missing))
            )

        decisions: dict[str, dict[str, str]] = {}
        for row_number, raw_row in enumerate(reader, start=2):
            if None in raw_row:
                raise ValueError(
                    f"review decision CSV row {row_number} has more values than columns"
                )
            cleaned = {str(key): _clean(value) for key, value in raw_row.items()}
            canonical_id = cleaned.get("canonical_candidate_id", "")
            if not canonical_id:
                if any(cleaned.values()):
                    raise ValueError(
                        f"review decision CSV row {row_number} has no "
                        "canonical_candidate_id"
                    )
                continue
            if canonical_id in decisions:
                raise ValueError(
                    "review decision CSV repeats canonical_candidate_id "
                    f"{canonical_id!r}"
                )
            try:
                decision = _normalize_review_decision(
                    cleaned.get("review_decision", "")
                )
            except ValueError as exc:
                raise ValueError(
                    f"review decision CSV row {row_number}: {exc}"
                ) from exc
            decisions[canonical_id] = {
                "review_decision": decision,
                "reviewer": cleaned.get("reviewer", ""),
                "reviewed_at": (
                    cleaned.get("reviewed_at", "") or cleaned.get("review_date", "")
                ),
                "review_notes": (
                    cleaned.get("review_notes", "") or cleaned.get("notes", "")
                ),
                "reviewed_value": cleaned.get("reviewed_value", ""),
                "reviewed_unit": cleaned.get("reviewed_unit", ""),
                "reviewed_value_qualifier": cleaned.get("reviewed_value_qualifier", ""),
                "reviewed_analyte": cleaned.get("reviewed_analyte", ""),
                "reviewed_endpoint_semantics": cleaned.get(
                    "reviewed_endpoint_semantics", ""
                ),
                "reviewed_reference_basis": cleaned.get("reviewed_reference_basis", ""),
                "reviewed_extravascular_route": cleaned.get(
                    "reviewed_extravascular_route", ""
                ),
                "reviewed_formulation": cleaned.get("reviewed_formulation", ""),
                "reviewed_population": cleaned.get("reviewed_population", ""),
                "review_decision_row_json": json.dumps(
                    cleaned,
                    sort_keys=True,
                    separators=(",", ":"),
                    ensure_ascii=True,
                ),
            }
    return decisions, _sha256(path)


_REVIEWED_BA_QUALIFIERS = {"point", "mean", "central_estimate"}
_REVIEWED_EXTRAVASCULAR_ROUTES = {
    "oral",
    "subcutaneous",
    "intramuscular",
    "inhalation",
    "intranasal",
    "transdermal",
    "sublingual",
    "buccal",
    "vaginal",
    "rectal",
    "topical",
}
_REVIEWED_BA_REFERENCE_BASES = {
    "explicit_absolute_bioavailability",
    "explicit_iv_comparator",
}
_REVIEWED_MAX_DOSE_QUALIFIERS = {
    "maximum",
    "recommended_maximum",
}
_REVIEWED_MAX_DOSE_REFERENCE_BASES = {
    "explicit_recommended_adult_maximum",
}
_REVIEWED_MAX_DOSE_ROUTES = _REVIEWED_EXTRAVASCULAR_ROUTES | {
    "intravenous",
    "submucosal",
}


def _maximum_dose_structured_review_errors(
    row: Mapping[str, Any],
) -> list[str]:
    required = (
        "reviewed_value",
        "reviewed_unit",
        "reviewed_value_qualifier",
        "reviewed_analyte",
        "reviewed_endpoint_semantics",
        "reviewed_reference_basis",
        "reviewed_extravascular_route",
        "reviewed_population",
    )
    errors = [f"missing_{field}" for field in required if not _clean(row.get(field))]
    try:
        reviewed_value = float(_clean(row.get("reviewed_value")))
    except ValueError:
        errors.append("reviewed_value_not_numeric")
    else:
        if reviewed_value <= 0:
            errors.append("reviewed_value_not_positive")
    unit = _canonical_text(row.get("reviewed_unit"))
    if not re.fullmatch(r"mg(?:/(?:kg|m2|day|h|dose)){0,2}", unit):
        errors.append("reviewed_maximum_dose_unit_not_supported")
    if (
        _canonical_text(row.get("reviewed_value_qualifier"))
        not in _REVIEWED_MAX_DOSE_QUALIFIERS
    ):
        errors.append("reviewed_value_qualifier_not_supported")
    if _canonical_text(row.get("reviewed_analyte")) != _canonical_text(
        row.get("drug_id")
    ):
        errors.append("reviewed_analyte_does_not_match_drug")
    if (
        _canonical_text(row.get("reviewed_endpoint_semantics"))
        not in _MAX_DOSE_ENDPOINTS
    ):
        errors.append("reviewed_endpoint_not_supported_maximum_dose")
    if (
        _canonical_text(row.get("reviewed_reference_basis"))
        not in _REVIEWED_MAX_DOSE_REFERENCE_BASES
    ):
        errors.append("reviewed_reference_basis_not_supported")
    routes = {
        _canonical_text(route)
        for route in re.split(
            r"[;,]",
            _clean(row.get("reviewed_extravascular_route")),
        )
        if _canonical_text(route)
    }
    if not routes or not routes.issubset(_REVIEWED_MAX_DOSE_ROUTES):
        errors.append("reviewed_route_not_supported")
    population = _canonical_text(row.get("reviewed_population"))
    if "adult" not in population or "pediatric" in population:
        errors.append("reviewed_population_not_explicitly_adult")
    return list(dict.fromkeys(errors))


def _structured_review_errors(row: Mapping[str, Any]) -> list[str]:
    candidate_type = _clean(row.get("candidate_type"))
    if candidate_type == "maximum_labeled_dose":
        return _maximum_dose_structured_review_errors(row)
    if candidate_type != "absolute_bioavailability":
        return ["structured_ambiguous_review_not_supported_for_candidate_type"]

    required = (
        "reviewed_value",
        "reviewed_unit",
        "reviewed_value_qualifier",
        "reviewed_analyte",
        "reviewed_endpoint_semantics",
        "reviewed_reference_basis",
        "reviewed_extravascular_route",
    )
    errors = [f"missing_{field}" for field in required if not _clean(row.get(field))]
    try:
        reviewed_value = float(_clean(row.get("reviewed_value")))
    except ValueError:
        errors.append("reviewed_value_not_numeric")
    else:
        if not 0 <= reviewed_value <= 100:
            errors.append("reviewed_value_outside_zero_to_100_percent")
    if _canonical_text(row.get("reviewed_unit")) not in {"%", "percent"}:
        errors.append("reviewed_unit_not_percent")
    if (
        _canonical_text(row.get("reviewed_value_qualifier"))
        not in _REVIEWED_BA_QUALIFIERS
    ):
        errors.append("reviewed_value_qualifier_not_supported")
    if _canonical_text(row.get("reviewed_analyte")) != _canonical_text(
        row.get("drug_id")
    ):
        errors.append("reviewed_analyte_does_not_match_drug")
    if (
        _canonical_text(row.get("reviewed_endpoint_semantics"))
        != "absolute_bioavailability"
    ):
        errors.append("reviewed_endpoint_not_absolute_bioavailability")
    if (
        _canonical_text(row.get("reviewed_reference_basis"))
        not in _REVIEWED_BA_REFERENCE_BASES
    ):
        errors.append("reviewed_reference_basis_not_supported")
    if (
        _canonical_text(row.get("reviewed_extravascular_route"))
        not in _REVIEWED_EXTRAVASCULAR_ROUTES
    ):
        errors.append("reviewed_extravascular_route_not_supported")
    return list(dict.fromkeys(errors))


def _apply_structured_review(row: dict[str, Any]) -> None:
    reviewed_value = float(_clean(row.get("reviewed_value")))
    route = _canonical_text(row.get("reviewed_extravascular_route"))
    if _clean(row.get("candidate_type")) == "maximum_labeled_dose":
        reviewed_unit = _canonical_text(row.get("reviewed_unit"))
        row.update(
            {
                "normalized_value": f"{reviewed_value:g}",
                "normalized_unit": reviewed_unit,
                "unit_basis": _unit_basis(reviewed_unit),
                "adjudicated_analyte": _clean(row.get("reviewed_analyte")),
                "parent_analyte_confirmed": "true",
                "route": route,
                "population": _clean(row.get("reviewed_population")),
                "endpoint_semantics": _canonical_text(
                    row.get("reviewed_endpoint_semantics")
                ),
                "semantic_context_status": "review_verified_sensitivity_only",
            }
        )
        return
    row.update(
        {
            "normalized_value": f"{reviewed_value:g}",
            "normalized_unit": "%",
            "unit_basis": "absolute_bioavailability_percent",
            "adjudicated_analyte": _clean(row.get("reviewed_analyte")),
            "parent_analyte_confirmed": "true",
            "extravascular_route": route,
            "iv_reference_confirmed": "true",
            "endpoint_semantics": (f"absolute_bioavailability:{route}_vs_intravenous"),
            "semantic_context_status": "review_verified",
        }
    )


_REVIEWABLE_MAXIMUM_DOSE_REJECTIONS = {
    "maximum_dose_interval_unspecified",
    "parent_analyte_not_in_value_clause",
    "combination_product_maximum_requires_manual_ingredient_binding",
    "maximum_dose_route_not_in_source_context",
}


def _maximum_dose_machine_rejection_reviewable(
    row: Mapping[str, Any],
) -> bool:
    if _clean(row.get("candidate_type")) != "maximum_labeled_dose":
        return False
    reasons = [
        reason
        for reason in _clean(row.get("machine_adjudication_reasons")).split(";")
        if reason
    ]
    if not reasons:
        return False
    return all(
        reason in _REVIEWABLE_MAXIMUM_DOSE_REJECTIONS
        or reason.startswith("same_clause_semantics_not_verified:")
        for reason in reasons
    )


def _apply_review_decisions(
    rows: Sequence[dict[str, Any]],
    decisions: Mapping[str, Mapping[str, str]],
    *,
    source_path: Path | None,
    source_sha256: str,
) -> None:
    canonical_ids = {_clean(row.get("canonical_candidate_id")) for row in rows}
    unknown_decisions = sorted(
        canonical_id
        for canonical_id, decision in decisions.items()
        if canonical_id not in canonical_ids and _clean(decision.get("review_decision"))
    )
    if unknown_decisions:
        raise ValueError(
            "review decision CSV references unknown canonical_candidate_id values: "
            + ", ".join(unknown_decisions)
        )

    for row in rows:
        canonical_id = _clean(row.get("canonical_candidate_id"))
        decision_row = decisions.get(canonical_id)
        if decision_row is None:
            if source_path is not None:
                row["review_application_status"] = "review_decision_missing"
            continue

        decision = _clean(decision_row.get("review_decision"))
        row.update(
            {
                "review_decision": decision,
                "reviewer": _clean(decision_row.get("reviewer")),
                "reviewed_at": _clean(decision_row.get("reviewed_at")),
                "review_notes": _clean(decision_row.get("review_notes")),
                "reviewed_value": _clean(decision_row.get("reviewed_value")),
                "reviewed_unit": _clean(decision_row.get("reviewed_unit")),
                "reviewed_value_qualifier": _clean(
                    decision_row.get("reviewed_value_qualifier")
                ),
                "reviewed_analyte": _clean(decision_row.get("reviewed_analyte")),
                "reviewed_endpoint_semantics": _clean(
                    decision_row.get("reviewed_endpoint_semantics")
                ),
                "reviewed_reference_basis": _clean(
                    decision_row.get("reviewed_reference_basis")
                ),
                "reviewed_extravascular_route": _clean(
                    decision_row.get("reviewed_extravascular_route")
                ),
                "reviewed_formulation": _clean(
                    decision_row.get("reviewed_formulation")
                ),
                "reviewed_population": _clean(decision_row.get("reviewed_population")),
                "review_decision_source": str(source_path) if source_path else "",
                "review_decision_file_sha256": source_sha256,
                "review_decision_row_json": _clean(
                    decision_row.get("review_decision_row_json")
                ),
            }
        )
        if not decision:
            row["review_application_status"] = "review_decision_pending"
        elif decision == "rejected":
            row["review_application_status"] = "review_rejected"
        elif decision == "deferred":
            row["review_application_status"] = "review_deferred"
        elif _clean(
            row.get("machine_adjudication_status")
        ) == "rejected" and not _maximum_dose_machine_rejection_reviewable(row):
            row["review_application_status"] = "approval_ineligible_machine_" + _clean(
                row.get("machine_adjudication_status")
            )
        elif not all(
            _clean(row.get(field))
            for field in ("reviewer", "reviewed_at", "review_notes")
        ):
            row["review_application_status"] = (
                "approval_ineligible_missing_review_metadata"
            )
        else:
            machine_status = _clean(row.get("machine_adjudication_status"))
            structured_errors = (
                _structured_review_errors(row)
                if machine_status == "ambiguous"
                or _maximum_dose_machine_rejection_reviewable(row)
                else []
            )
            if structured_errors:
                row["review_application_status"] = (
                    "approval_ineligible_structured_review:"
                    + ";".join(structured_errors)
                )
            else:
                if machine_status == "ambiguous":
                    _apply_structured_review(row)
                sensitivity_only = _clean(row.get("sensitivity_only")) == "true"
                row["review_approved"] = "true"
                row["review_application_status"] = (
                    "approved_sensitivity_only"
                    if sensitivity_only
                    else "approved_for_context_training"
                )
                row["training_allowed"] = _bool_text(not sensitivity_only)


def _review_template_rows(
    rows: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    groups: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[_clean(row.get("canonical_candidate_id"))].append(row)

    template: list[dict[str, Any]] = []
    for canonical_id in sorted(groups):
        group = sorted(groups[canonical_id], key=_row_identity)
        eligible_rows = [
            row
            for row in group
            if _clean(row.get("machine_adjudication_status"))
            in {"accepted", "ambiguous"}
            or _maximum_dose_machine_rejection_reviewable(row)
        ]
        canonical = eligible_rows[0] if eligible_rows else group[0]
        machine_status = _clean(canonical.get("machine_adjudication_status"))
        if (
            machine_status == "rejected"
            and not _maximum_dose_machine_rejection_reviewable(canonical)
        ):
            continue
        template.append(
            {
                column: _clean(canonical.get(column))
                for column in _REVIEW_TEMPLATE_COLUMNS
            }
            | {
                "eligible_for_approval": _bool_text(
                    machine_status in {"accepted", "ambiguous"}
                    or _maximum_dose_machine_rejection_reviewable(canonical)
                )
            }
        )
    return template


def adjudicate_spl_pk_rows(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Adjudicate candidate rows without dropping source records or averaging values."""
    adjudicated: list[dict[str, Any]] = []
    for row_number, source_row in enumerate(rows, start=1):
        row: dict[str, Any] = {
            str(key): _clean(value) for key, value in source_row.items()
        }
        row["input_row_number"] = row_number
        row.update(_decision(row))
        row.update(
            {
                "canonical_candidate_id": "",
                "duplicate_group_size": 1,
                "collapsed_duplicate_count": 0,
                "duplicate_candidate_ids": _clean(row.get("candidate_id")),
                "source_record_ids": _clean(row.get("source_record_id")),
                "spl_versions": _clean(row.get("spl_version")),
                "cache_file_sha256s": _clean(row.get("cache_file_sha256")),
                "section_text_sha256s": _clean(row.get("section_text_sha256")),
                "source_lineage_json": "",
                "conflict_id": "",
                "conflicting_normalized_values": "",
                "adjudication_schema_version": SCHEMA_VERSION,
            }
        )
        _assign_evidence_keys(row)
        adjudicated.append(row)

    conflicts = _mark_conflicts(adjudicated)
    _apply_duplicate_metadata(adjudicated)
    _initialize_review_metadata(adjudicated)
    return adjudicated, conflicts


def _collapsed_accepted(
    rows: Sequence[dict[str, Any]],
    *,
    candidate_type: str,
) -> list[dict[str, Any]]:
    groups: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if (
            _clean(row.get("candidate_type")) == candidate_type
            and (
                _clean(row.get("adjudication_status")) == "accepted"
                or _clean(row.get("review_approved")) == "true"
            )
            and _clean(row.get("review_decision")) != "rejected"
        ):
            groups[_clean(row.get("semantic_evidence_key"))].append(row)
    collapsed: list[dict[str, Any]] = []
    for key in sorted(groups):
        group = sorted(groups[key], key=_row_identity)
        canonical = dict(group[0])
        canonical.update(
            {
                "duplicate_group_size": len(group),
                "collapsed_duplicate_count": max(0, len(group) - 1),
                "duplicate_candidate_ids": _joined_values(group, "candidate_id"),
                "source_record_ids": _joined_values(group, "source_record_id"),
                "spl_versions": _joined_values(group, "spl_version"),
                "cache_file_sha256s": _joined_values(group, "cache_file_sha256"),
                "section_text_sha256s": _joined_values(group, "section_text_sha256"),
                "source_lineage_json": _lineage(group),
            }
        )
        collapsed.append(canonical)
    return collapsed


def _read_candidates(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"candidate CSV has no header: {path}")
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError("candidate CSV contains duplicate column names")
        missing = sorted(_REQUIRED_INPUT_COLUMNS.difference(fieldnames))
        if missing:
            raise ValueError(
                f"candidate CSV missing required columns: {', '.join(missing)}"
            )
        rows = []
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(
                    f"candidate CSV row {row_number} has more values than columns"
                )
            rows.append({key: value or "" for key, value in row.items()})
    return fieldnames, rows


def _output_columns(input_columns: Sequence[str]) -> list[str]:
    return list(dict.fromkeys([*input_columns, *_ADJUDICATION_COLUMNS]))


def _write_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    columns: Sequence[str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def adjudicate_spl_pk_candidates(
    *,
    candidates: str | Path,
    out_dir: str | Path,
    review_decisions: str | Path | None = None,
) -> dict[str, Any]:
    """Adjudicate one SPL PK candidate CSV and write contextual evidence outputs."""
    candidate_path = Path(candidates)
    output = Path(out_dir)
    if not candidate_path.is_file():
        raise ValueError(f"candidate CSV not found: {candidate_path}")
    input_columns, source_rows = _read_candidates(candidate_path)
    source_rows = enrich_spl_semantic_evidence(source_rows)
    output.mkdir(parents=True, exist_ok=True)

    review_path = Path(review_decisions) if review_decisions is not None else None
    review_rows: dict[str, dict[str, str]] = {}
    review_sha256 = ""
    if review_path is not None:
        review_rows, review_sha256 = _read_review_decisions(review_path)

    adjudicated, conflicts = adjudicate_spl_pk_rows(source_rows)
    _apply_review_decisions(
        adjudicated,
        review_rows,
        source_path=review_path,
        source_sha256=review_sha256,
    )
    review_template = _review_template_rows(adjudicated)
    clearance = _collapsed_accepted(adjudicated, candidate_type="clearance")
    bioavailability = _collapsed_accepted(
        adjudicated,
        candidate_type="absolute_bioavailability",
    )
    maximum_doses = _collapsed_accepted(
        adjudicated,
        candidate_type="maximum_labeled_dose",
    )

    columns = _output_columns(input_columns)
    outputs: dict[str, tuple[Path, Sequence[Mapping[str, Any]], Sequence[str]]] = {
        "all_adjudicated": (output / ALL_ADJUDICATED_NAME, adjudicated, columns),
        "accepted_clearance_context": (
            output / ACCEPTED_CLEARANCE_NAME,
            clearance,
            columns,
        ),
        "accepted_absolute_bioavailability_context": (
            output / ACCEPTED_BIOAVAILABILITY_NAME,
            bioavailability,
            columns,
        ),
        "max_recommended_dose_sensitivity": (
            output / MAX_DOSE_SENSITIVITY_NAME,
            maximum_doses,
            columns,
        ),
        "conflict_summary": (
            output / CONFLICT_SUMMARY_NAME,
            conflicts,
            _CONFLICT_COLUMNS,
        ),
        "review_decisions_template": (
            output / REVIEW_TEMPLATE_NAME,
            review_template,
            _REVIEW_TEMPLATE_COLUMNS,
        ),
    }
    for path, rows, output_columns in outputs.values():
        _write_csv(path, rows, output_columns)

    status_counts = Counter(
        _clean(row.get("adjudication_status")) for row in adjudicated
    )
    review_counts = Counter(
        _clean(row.get("review_decision")) or "pending" for row in review_template
    )
    type_status_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    for row in adjudicated:
        type_status_counts[_clean(row.get("candidate_type"))][
            _clean(row.get("adjudication_status"))
        ] += 1
    output_manifest = {
        name: {
            "path": str(path),
            "rows": len(rows),
            "sha256": _sha256(path),
        }
        for name, (path, rows, _) in outputs.items()
    }
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "candidate_input": str(candidate_path),
        "candidate_input_sha256": _sha256(candidate_path),
        "review_decisions_input": str(review_path) if review_path else "",
        "review_decisions_input_sha256": review_sha256,
        "review_decision_rows": len(review_rows),
        "review_decision_counts": dict(sorted(review_counts.items())),
        "candidate_rows": len(source_rows),
        "candidate_columns": input_columns,
        "candidate_extraction_versions": _unique(
            [_clean(row.get("extraction_method")) for row in source_rows]
        ),
        "output_dir": str(output),
        "outputs": output_manifest,
        "adjudication_status_counts": {
            status: status_counts.get(status, 0) for status in ADJUDICATION_STATUSES
        },
        "candidate_type_status_counts": {
            candidate_type: {
                status: counts.get(status, 0) for status in ADJUDICATION_STATUSES
            }
            for candidate_type, counts in sorted(type_status_counts.items())
        },
        "accepted_semantic_evidence_rows": (
            len(clearance) + len(bioavailability) + len(maximum_doses)
        ),
        "final_accepted_context_rows_by_type": {
            "clearance": len(clearance),
            "absolute_bioavailability": len(bioavailability),
            "maximum_labeled_dose_sensitivity": len(maximum_doses),
        },
        "training_allowed_context_rows": sum(
            _clean(row.get("training_allowed")) == "true"
            for row in [*clearance, *bioavailability]
        ),
        "accepted_duplicate_rows_collapsed": sum(
            int(row["collapsed_duplicate_count"])
            for row in [*clearance, *bioavailability, *maximum_doses]
        ),
        "conflict_groups": len(conflicts),
        "policy": {
            "candidate_preservation": "all input rows retained in all_adjudicated.csv",
            "status_vocabulary": list(ADJUDICATION_STATUSES),
            "duplicate_handling": (
                "accepted endpoint views collapse stable semantic evidence keys and retain JSON "
                "source lineage"
            ),
            "conflict_handling": (
                "comparable conflicting normalized values are ambiguous and are never averaged"
            ),
            "clearance_semantics": sorted(_CLEARANCE_ENDPOINTS),
            "route_evidence": (
                "machine acceptance requires same-clause route evidence; an explicit "
                "review decision may supply a source-confirmed route with provenance"
            ),
            "review_authority": (
                "machine acceptance is non-training unless an explicit canonical-ID "
                "review decision approves it; ambiguous bioavailability and a fail-closed "
                "whitelist of parser-limited adult maximum-dose rows require structured "
                "reviewer-confirmed value, qualifier, analyte, route, endpoint, and "
                "reference basis; maximum-dose approvals remain sensitivity-only"
            ),
            "clearance_unit_pooling": (
                "absolute, body-weight-normalized, and body-surface-normalized units remain separate"
            ),
            "bioavailability": (
                "explicit source-reported absolute bioavailability may be human-reviewed "
                "without a locally repeated intravenous comparator; relative-only evidence "
                "is rejected and mixed absolute/relative clauses require review"
            ),
            "maximum_dose": (
                "same-clause explicit recommended adult dosage-and-administration "
                "maxima are sensitivity-only and never training-allowed"
            ),
        },
    }
    manifest_path = output / MANIFEST_NAME
    manifest["manifest_output"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


__all__ = [
    "ADJUDICATION_STATUSES",
    "MANIFEST_NAME",
    "REVIEW_DECISIONS",
    "REVIEW_TEMPLATE_NAME",
    "SCHEMA_VERSION",
    "adjudicate_spl_pk_candidates",
    "adjudicate_spl_pk_rows",
]
