from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from analysis.external.openfda_pk import _values


SEMANTIC_COLUMNS = (
    "source_span_verified",
    "semantic_context_status",
    "semantic_measurement_context",
    "semantic_context_reasons",
    "semantic_clause",
    "semantic_clause_start",
    "semantic_clause_end",
    "endpoint_evidence_text",
    "endpoint_evidence_start",
    "endpoint_evidence_end",
    "analyte_evidence_text",
    "analyte_evidence_start",
    "analyte_evidence_end",
    "route_evidence_text",
    "route_evidence_start",
    "route_evidence_end",
    "route_evidence_json",
    "population_evidence_text",
    "population_evidence_start",
    "population_evidence_end",
    "adjudicated_analyte",
)

_ENDPOINT_PATTERNS = {
    "clearance": re.compile(
        r"\b(?:(?:apparent(?:\s+(?:oral|plasma))?|blood|drug|plasma|renal|systemic|"
        r"total(?:\s+body)?)\s+)?clearance\b|\bCL(?:/F|r)?\b",
        re.IGNORECASE,
    ),
    "absolute_bioavailability": re.compile(
        r"\babsolute\s+(?:oral\s+)?bioavailability\b|\bAbsBio\b",
        re.IGNORECASE,
    ),
    "maximum_labeled_dose": re.compile(
        r"\b(?:maximum|maximal|max\.)\s+(?:(?:recommended|approved|single|"
        r"total|daily|adult|labeled|labelled|maintenance|treatment)\s+){0,5}"
        r"(?:dose|dosage)\b|\bnot\s+(?:to\s+)?exceed\b|"
        r"\bup\s+to\s+(?:a\s+)?maximum\b|\ba\s+maximum\s+of\b",
        re.IGNORECASE,
    ),
}
_ROUTES = (
    ("intravenous", re.compile(r"\b(?:intravenous(?:ly)?|IV)\b", re.IGNORECASE)),
    ("oral", re.compile(r"\b(?:oral(?:ly)?|PO)\b", re.IGNORECASE)),
    ("subcutaneous", re.compile(r"\bsubcutaneous(?:ly)?\b", re.IGNORECASE)),
    ("intramuscular", re.compile(r"\bintramuscular(?:ly)?\b", re.IGNORECASE)),
    ("inhalation", re.compile(r"\b(?:inhaled|inhalation)\b", re.IGNORECASE)),
    ("transdermal", re.compile(r"\btransdermal(?:ly)?\b", re.IGNORECASE)),
    ("intranasal", re.compile(r"\b(?:intranasal(?:ly)?|nasal)\b", re.IGNORECASE)),
)
_HUMAN = re.compile(
    r"\b(?:adult|adults|healthy|human|humans|male|males|female|females|patient|"
    r"patients|subject|subjects|volunteer|volunteers|participant|participants)\b",
    re.IGNORECASE,
)
_NON_HUMAN = re.compile(
    r"\b(?:animal|animals|dog|dogs|mice|mouse|monkey|monkeys|rabbit|rabbits|"
    r"rat|rats|rodent|rodents)\b",
    re.IGNORECASE,
)
_PEDIATRIC = re.compile(
    r"\b(?:pediatric|paediatric|child|children|adolescent|infant|neonate)\b",
    re.IGNORECASE,
)
_ADULT = re.compile(r"\b(?:adult|adults|geriatric|elderly)\b", re.IGNORECASE)
_RECOMMENDED_DOSE = re.compile(
    r"\b(?:recommended|recommendation|should\s+not\s+exceed|"
    r"must\s+not\s+exceed|do\s+not\s+exceed|not\s+(?:to\s+)?exceed)\b",
    re.IGNORECASE,
)
_COVARIATE_CLEARANCE = re.compile(
    r"\b(?:creatinine\s+clearance|CrCL|eGFR|GFR|glomerular\s+filtration)\b",
    re.IGNORECASE,
)
_STUDY_MAXIMUM = re.compile(
    r"\b(?:maximum\s+tolerated\s+dose|MTD|highest\s+(?:dose\s+)?studied|"
    r"dose[- ]escalation)\b",
    re.IGNORECASE,
)
_RELATIVE_BIOAVAILABILITY = re.compile(
    r"\brelative\s+(?:oral\s+)?bioavailability\b", re.IGNORECASE
)
_RANGE_OR_INEQUALITY = re.compile(
    r"(?:\d\s*(?:-|\u2013|\u2014|to)\s*\d)|(?:<=|>=|<|>)|(?:\+/-|\u00b1)",
    re.IGNORECASE,
)
_CLEARANCE_VALUE_IN_CLAUSE = re.compile(
    r"\b\d+(?:\.\d+)?\s*(?:mL|L)\s*/\s*(?:min|h|hr|hour)\b",
    re.IGNORECASE,
)
_PERCENT_VALUE_IN_CLAUSE = re.compile(
    r"\b\d+(?:\.\d+)?\s*(?:%|percent)(?!\w)", re.IGNORECASE
)
_DOSE_VALUE_IN_CLAUSE = re.compile(
    r"\b\d+(?:\.\d+)?\s*(?:mcg|ug|\u00b5g|mg|g)(?:\s*/\s*(?:kg|day))?\b",
    re.IGNORECASE,
)
_TRUNCATED_CLEARANCE_UNIT = re.compile(
    r"^\s*/\s*(?:kg|kilograms?|m\s*(?:2|\u00b2)|square\s+meters?)\b",
    re.IGNORECASE,
)


def _clean(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


_CLAUSE_BOUNDARY_RE = re.compile(r";|(?<!\d)[.!?](?!\d)|(?<=\d)[!?](?!\d)")


def _clause_bounds(text: str, start: int, end: int) -> tuple[int, int]:
    boundaries = list(_CLAUSE_BOUNDARY_RE.finditer(text))
    left = max(
        (match.end() for match in boundaries if match.end() <= start),
        default=0,
    )
    right = min(
        (match.end() for match in boundaries if match.start() >= end),
        default=len(text),
    )
    return left, right


def _nearest_match(
    pattern: re.Pattern[str], text: str, start: int, end: int
) -> re.Match[str] | None:
    matches = list(pattern.finditer(text))
    if not matches:
        return None
    return min(
        matches,
        key=lambda match: max(match.start() - end, start - match.end(), 0),
    )


def _drug_pattern(drug_id: str) -> re.Pattern[str] | None:
    tokens = re.findall(r"[a-z0-9]+", drug_id.casefold())
    if not tokens:
        return None
    common_salts = {
        "acetate",
        "calcium",
        "fumarate",
        "hydrochloride",
        "maleate",
        "mesylate",
        "phosphate",
        "potassium",
        "sodium",
        "succinate",
        "sulfate",
        "tartrate",
    }
    while len(tokens) > 1 and tokens[-1] in common_salts:
        tokens.pop()
    return re.compile(
        r"\b" + r"[\s-]+".join(map(re.escape, tokens)) + r"\b", re.IGNORECASE
    )


@lru_cache(maxsize=512)
def _load_payload(path: str) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("SPL cache payload is not an object")
    return value


def _records(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    value = payload.get("records")
    if isinstance(value, dict):
        return [value]
    return [record for record in value or [] if isinstance(record, dict)]


@lru_cache(maxsize=4096)
def _section_text(
    cache_file: str,
    source_record_id: str,
    spl_version: str,
    section: str,
    section_item_index: int,
) -> str:
    payload = _load_payload(cache_file)
    candidates = []
    for record in _records(payload):
        record_id = _clean(record.get("set_id") or record.get("id"))
        version = _clean(record.get("version"))
        if source_record_id and record_id != source_record_id:
            continue
        if spl_version and version and version != spl_version:
            continue
        candidates.append(record)
    if len(candidates) != 1:
        return ""
    items = [
        _clean(value) for value in _values(candidates[0].get(section)) if _clean(value)
    ]
    if section_item_index < 0 or section_item_index >= len(items):
        return ""
    return items[section_item_index]


def _endpoint_context(candidate_type: str, endpoint_text: str) -> str:
    value = endpoint_text.casefold()
    if candidate_type == "absolute_bioavailability":
        return "absolute_bioavailability"
    if candidate_type == "maximum_labeled_dose":
        if "single" in value:
            return "maximum_single_dose"
        if "daily" in value:
            return "maximum_daily_dose"
        if "not" in value and "exceed" in value:
            return "not_to_exceed_dose"
        return "maximum_labeled_dose_unspecified_interval"
    if "renal" in value or value.strip().casefold() == "clr":
        return "renal_clearance"
    if "apparent" in value and "plasma" in value:
        return "apparent_plasma_clearance"
    if "oral" in value or "apparent" in value or "cl/f" in value:
        return "apparent_oral_clearance"
    if "total" in value and "body" in value:
        return "total_body_clearance"
    for prefix in ("systemic", "plasma", "blood"):
        if prefix in value:
            return f"{prefix}_clearance"
    return "clearance_unspecified"


def _span(
    match: re.Match[str] | None, *, offset: int = 0
) -> tuple[str, int | str, int | str]:
    if match is None:
        return "", "", ""
    return match.group(0), offset + match.start(), offset + match.end()


def semantic_evidence_for_text(
    text: str,
    row: Mapping[str, Any],
) -> dict[str, Any]:
    reasons: list[str] = []
    try:
        value_start = int(row.get("match_start") or 0)
        value_end = int(row.get("match_end") or 0)
    except (TypeError, ValueError):
        value_start = value_end = 0
    if value_start < 0 or value_end <= value_start or value_end > len(text):
        reasons.append("value_span_invalid")
        value_start = max(0, min(value_start, len(text)))
        value_end = max(value_start, min(value_end, len(text)))
    clause_start, clause_end = _clause_bounds(text, value_start, value_end)
    clause = text[clause_start:clause_end]
    window_start = max(0, clause_start - 220)
    window_end = min(len(text), clause_end + 220)
    window = text[window_start:window_end]

    candidate_type = _clean(row.get("candidate_type"))
    endpoint_pattern = _ENDPOINT_PATTERNS.get(candidate_type)
    endpoint = (
        _nearest_match(
            endpoint_pattern,
            clause,
            value_start - clause_start,
            value_end - clause_start,
        )
        if endpoint_pattern
        else None
    )
    if endpoint is None:
        reasons.append("endpoint_not_in_value_clause")
    elif (
        max(
            endpoint.start() - (value_end - clause_start),
            (value_start - clause_start) - endpoint.end(),
            0,
        )
        > 48
    ):
        reasons.append("endpoint_too_far_from_value")

    raw_value = _clean(row.get("raw_value"))
    if _RANGE_OR_INEQUALITY.search(raw_value):
        reasons.append("range_inequality_or_variability_requires_manual_context")
    if candidate_type in {
        "clearance",
        "absolute_bioavailability",
    } and _RANGE_OR_INEQUALITY.search(clause):
        reasons.append(
            "range_or_variability_in_endpoint_clause_requires_manual_context"
        )
    if (
        candidate_type == "clearance"
        and len(_CLEARANCE_VALUE_IN_CLAUSE.findall(clause)) > 1
    ):
        reasons.append("multiple_clearance_values_in_clause")
    if (
        candidate_type == "absolute_bioavailability"
        and len(_PERCENT_VALUE_IN_CLAUSE.findall(clause)) > 1
    ):
        reasons.append("multiple_bioavailability_values_in_clause")
    if candidate_type in {"clearance", "absolute_bioavailability"} and re.search(
        r"\bvs\.?\b|\bversus\b", clause, re.IGNORECASE
    ):
        reasons.append("comparative_or_stratified_endpoint_requires_manual_context")
    if (
        candidate_type == "maximum_labeled_dose"
        and len(_DOSE_VALUE_IN_CLAUSE.findall(clause)) > 1
    ):
        reasons.append("multiple_dose_values_in_clause")
    trailing = text[value_end : min(len(text), value_end + 28)]
    if candidate_type == "clearance" and _TRUNCATED_CLEARANCE_UNIT.search(trailing):
        reasons.append("clearance_unit_truncated_after_numeric_span")
    if candidate_type == "clearance" and endpoint is not None:
        before = clause[max(0, endpoint.start() - 28) : endpoint.end()]
        between_start = min(endpoint.start(), value_start - clause_start)
        between_end = max(endpoint.end(), value_end - clause_start)
        binding_text = clause[
            max(0, between_start - 30) : min(len(clause), between_end + 30)
        ]
        if _COVARIATE_CLEARANCE.search(before) or _COVARIATE_CLEARANCE.search(
            binding_text
        ):
            reasons.append("subject_renal_function_not_drug_clearance")

    drug_pattern = _drug_pattern(_clean(row.get("drug_id")))
    analyte = drug_pattern.search(clause) if drug_pattern else None
    if analyte is None:
        reasons.append("parent_analyte_not_in_value_clause")
    if re.search(r"\bmetabolite\b", clause, re.IGNORECASE):
        reasons.append("metabolite_context_requires_separate_identity")

    same_clause_route_required = candidate_type in {
        "clearance",
        "absolute_bioavailability",
    }
    route_scope = clause if same_clause_route_required else window
    route_scope_start = clause_start if same_clause_route_required else window_start
    route_matches: list[tuple[str, re.Match[str]]] = []
    for name, pattern in _ROUTES:
        route_match = pattern.search(route_scope)
        if route_match is not None:
            route_matches.append((name, route_match))
    route = (
        min(
            route_matches,
            key=lambda item: max(
                item[1].start() - (value_end - route_scope_start),
                (value_start - route_scope_start) - item[1].end(),
                0,
            ),
        )
        if route_matches
        else None
    )
    human = _HUMAN.search(window)
    if _NON_HUMAN.search(window):
        reasons.append("non_human_context")
    elif candidate_type in {"clearance", "absolute_bioavailability"} and human is None:
        reasons.append("explicit_human_study_context_not_found")

    if candidate_type == "clearance" and route is None:
        reasons.append("administration_route_not_in_value_clause")
    semantic_endpoint = _endpoint_context(
        candidate_type, endpoint.group(0) if endpoint is not None else ""
    )
    if (
        candidate_type == "clearance"
        and semantic_endpoint
        in {
            "blood_clearance",
            "plasma_clearance",
            "systemic_clearance",
            "total_body_clearance",
        }
        and route is not None
        and route[0] != "intravenous"
    ):
        reasons.append("non_iv_clearance_not_explicitly_apparent")
    if candidate_type == "absolute_bioavailability":
        if _RELATIVE_BIOAVAILABILITY.search(clause):
            reasons.append("relative_bioavailability_not_absolute")
        routes = {name for name, _ in route_matches}
        if "intravenous" not in routes:
            reasons.append("intravenous_reference_not_in_value_clause")
        if not routes.difference({"intravenous"}):
            reasons.append("extravascular_test_route_not_in_value_clause")
    if candidate_type == "maximum_labeled_dose":
        label_names = [
            name.strip()
            for field in ("label_generic_names", "label_substance_names")
            for name in _clean(row.get(field)).split(";")
            if name.strip()
        ]
        if len({name.casefold() for name in label_names}) > 1:
            reasons.append(
                "combination_product_maximum_requires_manual_ingredient_binding"
            )
        if route is None:
            reasons.append("maximum_dose_route_not_in_source_context")
        if re.search(r"\binitial\s+dose\b", clause, re.IGNORECASE):
            reasons.append("initial_dose_cap_not_maximum_recommended_adult_dose")
        if not _ADULT.search(clause):
            reasons.append("explicit_adult_context_not_in_value_clause")
        if not _RECOMMENDED_DOSE.search(clause):
            reasons.append("explicit_recommended_dose_context_not_in_value_clause")
        if _PEDIATRIC.search(clause):
            reasons.append("pediatric_or_age_specific_maximum")
        if _STUDY_MAXIMUM.search(clause):
            reasons.append("highest_studied_or_mtd_not_recommended_maximum")
        if _clean(row.get("section")) != "dosage_and_administration":
            reasons.append("maximum_not_in_dosage_and_administration")

    endpoint_text, endpoint_start, endpoint_end = _span(endpoint, offset=clause_start)
    semantic_measurement_context = _endpoint_context(candidate_type, endpoint_text)
    analyte_text, analyte_start, analyte_end = _span(analyte, offset=clause_start)
    route_text = route[0] if route else ""
    route_start = route_scope_start + route[1].start() if route else ""
    route_end = route_scope_start + route[1].end() if route else ""
    route_evidence = [
        {
            "route": name,
            "text": match.group(0),
            "start": route_scope_start + match.start(),
            "end": route_scope_start + match.end(),
        }
        for name, match in sorted(route_matches, key=lambda item: item[1].start())
    ]
    human_text, human_start, human_end = _span(human, offset=window_start)
    return {
        "semantic_context_status": "verified" if not reasons else "not_verified",
        "semantic_measurement_context": semantic_measurement_context,
        "semantic_context_reasons": ";".join(dict.fromkeys(reasons)),
        "semantic_clause": clause,
        "semantic_clause_start": clause_start,
        "semantic_clause_end": clause_end,
        "endpoint_evidence_text": endpoint_text,
        "endpoint_evidence_start": endpoint_start,
        "endpoint_evidence_end": endpoint_end,
        "analyte_evidence_text": analyte_text,
        "analyte_evidence_start": analyte_start,
        "analyte_evidence_end": analyte_end,
        "route_evidence_text": route_text,
        "route_evidence_start": route_start,
        "route_evidence_end": route_end,
        "route_evidence_json": json.dumps(
            route_evidence,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ),
        "population_evidence_text": human_text,
        "population_evidence_start": human_start,
        "population_evidence_end": human_end,
        "adjudicated_analyte": _clean(row.get("drug_id")) if analyte else "",
    }


def enrich_spl_semantic_evidence(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    enriched: list[dict[str, Any]] = []
    for source in rows:
        row = dict(source)
        try:
            text = _section_text(
                _clean(row.get("cache_file")),
                _clean(row.get("source_record_id")),
                _clean(row.get("spl_version")),
                _clean(row.get("section")),
                int(row.get("section_item_index") or 0),
            )
        except (OSError, ValueError, json.JSONDecodeError):
            text = ""
        expected_hash = _clean(row.get("section_text_sha256"))
        verified = (
            bool(text)
            and hashlib.sha256(text.encode("utf-8")).hexdigest() == expected_hash
        )
        row["source_span_verified"] = _bool_text(verified)
        if verified:
            row.update(semantic_evidence_for_text(text, row))
        else:
            row.update(
                {
                    column: ""
                    for column in SEMANTIC_COLUMNS
                    if column != "source_span_verified"
                }
            )
            row["semantic_context_status"] = "not_verified"
            row["semantic_context_reasons"] = "source_section_or_hash_not_verified"
        enriched.append(row)
    return enriched


__all__ = [
    "SEMANTIC_COLUMNS",
    "enrich_spl_semantic_evidence",
    "semantic_evidence_for_text",
]
