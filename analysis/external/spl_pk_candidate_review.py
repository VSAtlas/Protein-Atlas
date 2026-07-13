from __future__ import annotations

from collections import Counter
from collections.abc import Collection, Iterator, Mapping
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Any

from analysis.external.openfda_pk import (
    _dose_context,
    _openfda_values,
    _values,
)


SOURCE_NAME = "DailyMed_openFDA_SPL_cache"
SCHEMA_VERSION = "atlas_spl_pk_candidate_review_v1"
DEFAULT_EXCERPT_CHARS = 480
MIN_EXCERPT_CHARS = 160
MAX_EXCERPT_CHARS = 800
BASE_EXCLUSION = "candidate_only_manual_review_required"

CANDIDATE_COLUMNS = (
    "candidate_id",
    "candidate_type",
    "candidate_status",
    "drug_id",
    "cache_file",
    "cache_file_sha256",
    "cache_search",
    "cache_status",
    "cache_last_updated",
    "source_name",
    "source_record_id",
    "source_set_id",
    "source_document_id",
    "spl_version",
    "effective_time",
    "source_url",
    "section",
    "section_item_index",
    "section_text_sha256",
    "match_start",
    "match_end",
    "source_excerpt",
    "raw_value",
    "raw_unit",
    "normalized_value",
    "normalized_unit",
    "normalization_status",
    "measurement_context",
    "analyte",
    "analyte_status",
    "label_generic_names",
    "label_substance_names",
    "route",
    "route_source",
    "formulation",
    "formulation_source",
    "population",
    "steady_state",
    "regimen",
    "confidence",
    "review_required",
    "exclusion_reason",
    "model_ready",
    "extraction_method",
)

_PK_SECTIONS = ("pharmacokinetics", "clinical_pharmacology")
_MAXIMUM_DOSE_SECTIONS = (
    "dosage_and_administration",
    "pharmacokinetics",
    "clinical_pharmacology",
)

_NUMBER = r"(?:\d+(?:,\d{3})*(?:\.\d+)?|\.\d+)"
_VALUE_EXPRESSION = (
    rf"(?P<value_expression>(?P<value>{_NUMBER})"
    rf"(?:\s*(?:(?P<range_connector>to|-|\u2013|\u2014)\s*"
    rf"(?P<range_value>{_NUMBER})|"
    rf"(?P<variability_connector>\+/-|\u00b1)\s*"
    rf"(?P<variability_value>{_NUMBER})))?)"
)

_DOSE_BASE_UNIT = (
    r"(?:micrograms?|milligrams?|nanograms?|grams?|mcg|ug|\u00b5g|ng|mg|g)"
)
_DOSE_DENOMINATOR = (
    r"(?:kg|kilograms?|m(?:2|\u00b2)|square\s+meters?|days?|hours?|hrs?|h|"
    r"minutes?|mins?|min|doses?)"
)
_DOSE_VALUE_RE = re.compile(
    rf"{_VALUE_EXPRESSION}\s*"
    rf"(?P<unit>{_DOSE_BASE_UNIT}(?:\s*(?:/|per)\s*"
    rf"{_DOSE_DENOMINATOR}){{0,2}})\b"
    r"(?!\s*(?:/|per)\s*(?:mL|L|dL)\b)"
    r"(?!\s*(?:[.*x\u00d7]\s*)?(?:h|hr|hours?)\s*/\s*(?:mL|L|dL)\b)"
    r"(?!\s*(?:/|per)\s*[A-Za-z])",
    re.IGNORECASE,
)

_CLEARANCE_VOLUME = r"(?:milliliters?|liters?|mL|L)"
_CLEARANCE_TIME = r"(?:seconds?|secs?|sec|minutes?|mins?|min|hours?|hrs?|hr|h|days?|day)"
_CLEARANCE_SCALE = r"(?:kg|kilograms?|m(?:2|\u00b2)|square\s+meters?)"
_CLEARANCE_UNIT = (
    rf"{_CLEARANCE_VOLUME}\s*(?:/|per)\s*(?:"
    rf"{_CLEARANCE_TIME}(?:\s*(?:/|per)\s*{_CLEARANCE_SCALE})?"
    rf"|{_CLEARANCE_SCALE}\s*(?:/|per)\s*{_CLEARANCE_TIME})"
    rf"(?:\s*(?:/|per)\s*1\.73\s*m(?:2|\u00b2))?"
)
_CLEARANCE_VALUE_RE = re.compile(
    rf"{_VALUE_EXPRESSION}\s*(?P<unit>{_CLEARANCE_UNIT})\b",
    re.IGNORECASE,
)
_BIOAVAILABILITY_VALUE_RE = re.compile(
    rf"{_VALUE_EXPRESSION}\s*(?P<unit>%|percent(?:age)?)?",
    re.IGNORECASE,
)

_CLEARANCE_KEYWORD_RE = re.compile(
    r"\b(?:(?:apparent|average|blood|body|drug|mean|oral|plasma|renal|"
    r"systemic|total)\s+){0,3}clearance\b|\bCL(?:/[A-Za-z]+)?\b",
    re.IGNORECASE,
)
_ABSOLUTE_BIOAVAILABILITY_RE = re.compile(
    r"\babsolute\s+(?:oral\s+)?bioavailability\b|\bAbsBio\b",
    re.IGNORECASE,
)
_DOSE_CUE_RE = re.compile(
    r"\b(?:administer(?:ed|ing|ation)?|dose(?:d|s)?|dosing|given|infus(?:ed|ion)|"
    r"inject(?:ed|ion)|received?)\b|\bfollowing\b|\bafter\b",
    re.IGNORECASE,
)
_STRONG_OBSERVED_DOSE_RE = re.compile(
    r"\b(?:administer(?:ed|ing|ation)?|dosed?|given|infus(?:ed|ion)|inject(?:ed|ion)|"
    r"received?)\b|\b(?:following|after)\b[^.;]{0,80}\b(?:dose|dosing|"
    r"administration|infusion|injection)\b|\b(?:subjects?|patients?|volunteers?)\b"
    r"[^.;]{0,100}\b(?:dose|dosed|received|administered)\b",
    re.IGNORECASE,
)
_MAXIMUM_DOSE_RE = re.compile(
    r"\b(?:maximum|maximal|max\.)\s+(?:(?:recommended|approved|single|total|daily|human|"
    r"labeled|labelled|tolerated|maintenance|treatment)\s+){0,5}(?:dose|dosage)\b|"
    r"\b(?:dose|dosage)\s+(?:(?:may|can|should|be|increased|titrated|to|a|the)\s+)"
    r"{0,8}(?:maximum|maximal)\b|\b(?:not\s+(?:to\s+)?exceed|up\s+to|"
    r"(?:to\s+)?a\s+maximum(?:\s+of)?)\b",
    re.IGNORECASE,
)

_METABOLITE_RE = re.compile(
    r"\b(?:active\s+|major\s+|primary\s+)?metabolites?\b",
    re.IGNORECASE,
)
_NON_HUMAN_RE = re.compile(
    r"\b(?:animal|animals|canine|dog|dogs|mice|mouse|monkey|monkeys|rabbit|"
    r"rabbits|rat|rats|rodent|rodents)\b",
    re.IGNORECASE,
)
_UNKNOWN_BIOAVAILABILITY_RE = re.compile(
    r"\babsolute\s+(?:oral\s+)?bioavailability\b[^.;]{0,45}\b"
    r"(?:not\s+known|unknown|has\s+not\s+been\s+determined)\b",
    re.IGNORECASE,
)

_ROUTE_PATTERNS = (
    ("intravenous", re.compile(r"\b(?:intravenous(?:ly)?|IV)\b", re.IGNORECASE)),
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
_FORMULATION_PATTERNS = (
    ("oral solution", re.compile(r"\boral\s+solution\b", re.IGNORECASE)),
    ("oral suspension", re.compile(r"\boral\s+suspension\b", re.IGNORECASE)),
    ("extended-release tablet", re.compile(r"\bextended[- ]release\s+tablets?\b", re.IGNORECASE)),
    ("extended-release capsule", re.compile(r"\bextended[- ]release\s+capsules?\b", re.IGNORECASE)),
    ("tablet", re.compile(r"\btablets?\b", re.IGNORECASE)),
    ("capsule", re.compile(r"\bcapsules?\b", re.IGNORECASE)),
    ("injection", re.compile(r"\binjections?\b", re.IGNORECASE)),
    ("infusion", re.compile(r"\binfusions?\b", re.IGNORECASE)),
    ("transdermal patch", re.compile(r"\b(?:transdermal\s+)?patch(?:es)?\b", re.IGNORECASE)),
    ("inhalation powder", re.compile(r"\binhalation\s+powder\b", re.IGNORECASE)),
    ("inhalation aerosol", re.compile(r"\binhalation\s+aerosol\b", re.IGNORECASE)),
    ("spray", re.compile(r"\bsprays?\b", re.IGNORECASE)),
    ("film", re.compile(r"\bfilms?\b", re.IGNORECASE)),
    ("implant", re.compile(r"\bimplants?\b", re.IGNORECASE)),
    ("cream", re.compile(r"\bcreams?\b", re.IGNORECASE)),
    ("gel", re.compile(r"\bgels?\b", re.IGNORECASE)),
    ("ointment", re.compile(r"\bointments?\b", re.IGNORECASE)),
)
_POPULATION_PATTERNS = (
    re.compile(
        r"\bhealthy(?:\s+(?:adult|male|female|young|elderly)){0,3}\s+"
        r"(?:subjects?|volunteers?|participants?)\b",
        re.IGNORECASE,
    ),
    re.compile(
        r"\b(?:adult|elderly|geriatric|pediatric)\s+(?:subjects?|patients?|"
        r"volunteers?|participants?)\b(?:\s+with\s+[^.;]{1,80})?",
        re.IGNORECASE,
    ),
    re.compile(r"\b(?:subjects?|patients?)\s+with\s+[^.;]{1,100}", re.IGNORECASE),
    re.compile(r"\bchildren(?:\s+(?:aged|ages?)\s+[^.;]{1,60})?", re.IGNORECASE),
    re.compile(r"\bpostmenopausal\s+women\b", re.IGNORECASE),
    re.compile(r"\b\d+\s+(?:subjects?|patients?|volunteers?|participants?)\b", re.IGNORECASE),
)


def _clean_text(value: Any) -> str:
    return " ".join(str(value or "").split())


def _unique(values: Collection[str]) -> list[str]:
    output: list[str] = []
    seen: set[str] = set()
    for value in values:
        clean = _clean_text(value)
        key = clean.casefold()
        if clean and key not in seen:
            output.append(clean)
            seen.add(key)
    return output


def _join(values: Collection[str]) -> str:
    return ";".join(_unique(values))


def _as_number(value: str) -> float:
    return float(value.replace(",", ""))


def _rounded(value: float) -> float:
    return float(f"{value:.12g}")


def _has_range_or_variability(match: re.Match[str]) -> bool:
    return bool(match.group("range_connector") or match.group("variability_connector"))


def _numeric_ambiguity(match: re.Match[str]) -> str:
    if match.group("range_connector"):
        return "numeric_range_requires_adjudication"
    if match.group("variability_connector"):
        return "summary_statistic_requires_adjudication"
    return ""


def _canonical_denominator(value: str) -> str | None:
    key = re.sub(r"\s+", "", value.casefold()).replace("\u00b2", "2")
    if key in {"kg", "kilogram", "kilograms"}:
        return "kg"
    if key in {"m2", "squaremeter", "squaremeters"}:
        return "m2"
    if key in {"day", "days"}:
        return "day"
    if key in {"h", "hr", "hrs", "hour", "hours"}:
        return "h"
    if key in {"min", "mins", "minute", "minutes"}:
        return "min"
    if key in {"dose", "doses"}:
        return "dose"
    return None


def _normalize_dose(
    match: re.Match[str],
) -> tuple[float | None, str, str, list[str]]:
    ambiguity = _numeric_ambiguity(match)
    if ambiguity:
        return None, "", "withheld_ambiguous_numeric_expression", [ambiguity]

    unit = match.group("unit")
    normalized = re.sub(r"\s*per\s*", "/", unit, flags=re.IGNORECASE)
    parts = [part.strip() for part in normalized.split("/")]
    base = parts[0].casefold().replace("\u00b5", "u")
    factors = {
        "g": 1000.0,
        "gram": 1000.0,
        "grams": 1000.0,
        "mg": 1.0,
        "milligram": 1.0,
        "milligrams": 1.0,
        "mcg": 0.001,
        "ug": 0.001,
        "microgram": 0.001,
        "micrograms": 0.001,
        "ng": 0.000001,
        "nanogram": 0.000001,
        "nanograms": 0.000001,
    }
    factor = factors.get(base)
    denominators = [_canonical_denominator(part) for part in parts[1:]]
    if factor is None or any(part is None for part in denominators):
        return None, "", "withheld_unsupported_unit", ["dose_unit_not_unambiguous"]
    value = _as_number(match.group("value"))
    if value <= 0:
        return None, "", "withheld_nonpositive_value", ["nonpositive_numeric_value"]
    suffix = "".join(f"/{part}" for part in denominators)
    return _rounded(value * factor), f"mg{suffix}", "normalized_unambiguous", []


def _canonical_time(value: str) -> tuple[str, float] | None:
    key = re.sub(r"\s+", "", value.casefold())
    if key in {"h", "hr", "hrs", "hour", "hours"}:
        return "h", 1.0
    if key in {"min", "mins", "minute", "minutes"}:
        return "h", 60.0
    if key in {"sec", "secs", "second", "seconds"}:
        return "h", 3600.0
    if key in {"day", "days"}:
        return "h", 1.0 / 24.0
    return None


def _normalize_clearance(
    match: re.Match[str],
) -> tuple[float | None, str, str, list[str]]:
    ambiguity = _numeric_ambiguity(match)
    if ambiguity:
        return None, "", "withheld_ambiguous_numeric_expression", [ambiguity]

    raw_unit = match.group("unit")
    if re.search(r"(?:/|per)\s*1\.73\s*m", raw_unit, re.IGNORECASE):
        return (
            None,
            "",
            "withheld_body_surface_area_reference",
            ["body_surface_area_reference_requires_adjudication"],
        )
    normalized = re.sub(r"\s*per\s*", "/", raw_unit, flags=re.IGNORECASE)
    parts = [part.strip() for part in normalized.split("/")]
    if len(parts) not in {2, 3}:
        return None, "", "withheld_unsupported_unit", ["clearance_unit_not_unambiguous"]

    volume_key = parts[0].casefold()
    if volume_key in {"l", "liter", "liters"}:
        volume_factor = 1.0
    elif volume_key in {"ml", "milliliter", "milliliters"}:
        volume_factor = 0.001
    else:
        return None, "", "withheld_unsupported_unit", ["clearance_unit_not_unambiguous"]

    time_part: tuple[str, float] | None = None
    scale = ""
    for part in parts[1:]:
        parsed_time = _canonical_time(part)
        parsed_scale = _canonical_denominator(part)
        if parsed_time is not None:
            if time_part is not None:
                return None, "", "withheld_unsupported_unit", ["multiple_time_denominators"]
            time_part = parsed_time
        elif parsed_scale in {"kg", "m2"}:
            if scale:
                return None, "", "withheld_unsupported_unit", ["multiple_scale_denominators"]
            scale = str(parsed_scale)
        else:
            return None, "", "withheld_unsupported_unit", ["clearance_unit_not_unambiguous"]
    if time_part is None:
        return None, "", "withheld_unsupported_unit", ["clearance_time_unit_missing"]

    value = _as_number(match.group("value"))
    if value <= 0:
        return None, "", "withheld_nonpositive_value", ["nonpositive_numeric_value"]
    _, time_factor = time_part
    suffix = f"/{scale}" if scale else ""
    return (
        _rounded(value * volume_factor * time_factor),
        f"L/h{suffix}",
        "normalized_unambiguous",
        [],
    )


def _normalize_bioavailability(
    match: re.Match[str],
) -> tuple[float | None, str, str, list[str]]:
    ambiguity = _numeric_ambiguity(match)
    if ambiguity:
        return None, "", "withheld_ambiguous_numeric_expression", [ambiguity]
    unit = _clean_text(match.group("unit"))
    if not unit:
        return None, "", "withheld_missing_unit", ["bioavailability_unit_missing"]
    value = _as_number(match.group("value"))
    if not 0 <= value <= 100:
        return None, "", "withheld_out_of_range", ["bioavailability_percent_out_of_range"]
    return _rounded(value), "%", "normalized_unambiguous", []


def _span_distance(first: re.Match[str], second: re.Match[str]) -> int:
    if first.end() < second.start():
        return second.start() - first.end()
    if second.end() < first.start():
        return first.start() - second.end()
    return 0


def _nearest_match(
    text: str,
    target: re.Match[str],
    pattern: re.Pattern[str],
    *,
    max_distance: int,
) -> re.Match[str] | None:
    start = max(0, target.start() - max_distance)
    end = min(len(text), target.end() + max_distance)
    matches = list(pattern.finditer(text, start, end))
    if not matches:
        return None
    nearest = min(matches, key=lambda item: _span_distance(item, target))
    return nearest if _span_distance(nearest, target) <= max_distance else None


def _short_excerpt(text: str, start: int, end: int, limit: int) -> str:
    limit = max(MIN_EXCERPT_CHARS, min(MAX_EXCERPT_CHARS, limit))
    anchor_length = min(end - start, limit)
    side = max(0, (limit - anchor_length) // 2)
    excerpt_start = max(0, start - side)
    excerpt_end = min(len(text), end + side)
    if excerpt_end - excerpt_start < limit:
        excerpt_start = max(0, excerpt_end - limit)
        excerpt_end = min(len(text), excerpt_start + limit)
    excerpt = text[excerpt_start:excerpt_end].strip()
    if excerpt_start > 0:
        excerpt = f"...{excerpt}"
    if excerpt_end < len(text):
        excerpt = f"{excerpt}..."
    return excerpt[:limit]


def _local_context(text: str, start: int, end: int, radius: int = 240) -> str:
    return text[max(0, start - radius) : min(len(text), end + radius)]


def _detect_analyte(
    context: str,
    record: dict[str, Any],
    drug_id: str,
) -> tuple[str, str, list[str], list[str], list[str]]:
    generic_names = _unique(_openfda_values(record, "generic_name"))
    substance_names = _unique(_openfda_values(record, "substance_name"))
    names = substance_names or generic_names
    if _METABOLITE_RE.search(context):
        return (
            "",
            "metabolite_context_ambiguous",
            generic_names,
            substance_names,
            ["metabolite_or_analyte_context_ambiguous"],
        )

    context_matches = [
        name
        for name in names
        if re.search(rf"\b{re.escape(name)}\b", context, re.IGNORECASE)
    ]
    if len(context_matches) == 1:
        return (
            context_matches[0],
            "context_label_match",
            generic_names,
            substance_names,
            [],
        )
    if len(names) == 1:
        return names[0], "single_label_analyte", generic_names, substance_names, []
    if len(names) > 1:
        return (
            "",
            "multiple_label_analytes",
            generic_names,
            substance_names,
            ["analyte_ambiguous"],
        )
    if drug_id:
        return drug_id, "cache_drug_id_fallback", generic_names, substance_names, []
    return "", "not_detected", generic_names, substance_names, ["analyte_not_detected"]


def _detect_route(context: str, record: dict[str, Any]) -> tuple[str, str]:
    detected = [name for name, pattern in _ROUTE_PATTERNS if pattern.search(context)]
    if detected:
        return _join(detected), "section_context"
    label_routes = [*_values(record.get("route")), *_openfda_values(record, "route")]
    return (_join(label_routes), "label_metadata") if label_routes else ("", "not_detected")


def _detect_formulation(context: str, record: dict[str, Any]) -> tuple[str, str]:
    detected = [name for name, pattern in _FORMULATION_PATTERNS if pattern.search(context)]
    if detected:
        return _join(detected), "section_context"
    label_forms = [
        *_values(record.get("dosage_form")),
        *_openfda_values(record, "dosage_form"),
    ]
    return (_join(label_forms), "label_metadata") if label_forms else ("", "not_detected")


def _detect_population(context: str) -> str:
    for pattern in _POPULATION_PATTERNS:
        match = pattern.search(context)
        if match:
            return _clean_text(match.group(0))[:120]
    animal = _NON_HUMAN_RE.search(context)
    return f"non-human: {animal.group(0).casefold()}" if animal else ""


def _regimen_context(context: str) -> tuple[str, str]:
    _, _, regimen, steady_state = _dose_context(context)
    if not regimen and re.search(r"\b(?:multiple|repeated)\s+doses?\b", context, re.IGNORECASE):
        regimen = "multiple dose"
    if re.search(r"\bsteady[ -]state\b", context, re.IGNORECASE):
        steady_state = "yes"
    elif re.search(r"\bsingle\s+dose\b", context, re.IGNORECASE):
        steady_state = "no"
    return regimen, steady_state or "unknown"


def _clearance_context(context: str) -> str:
    checks = (
        ("renal_clearance", r"\brenal\s+clearance\b|\bCLr\b"),
        ("apparent_oral_clearance", r"\bapparent\s+(?:oral\s+)?clearance\b|\bCL/F\b"),
        ("total_body_clearance", r"\btotal(?:\s+body)?\s+clearance\b"),
        ("systemic_clearance", r"\bsystemic\s+clearance\b"),
        ("plasma_clearance", r"\bplasma\s+clearance\b"),
        ("blood_clearance", r"\bblood\s+clearance\b"),
    )
    for label, pattern in checks:
        if re.search(pattern, context, re.IGNORECASE):
            return label
    return "clearance_unspecified"


def _maximum_dose_context(context: str) -> str:
    if re.search(r"\b(?:single|per)\s+(?:treatment\s+)?dose\b", context, re.IGNORECASE):
        return "maximum_single_dose"
    if re.search(r"\b(?:daily|per\s+day)\b", context, re.IGNORECASE):
        return "maximum_daily_dose"
    if re.search(r"\bnot\s+(?:to\s+)?exceed\b", context, re.IGNORECASE):
        return "not_to_exceed_dose"
    return "maximum_labeled_dose_unspecified_interval"


def _maximum_dose_directly_linked(
    text: str,
    value_match: re.Match[str],
    keyword_match: re.Match[str],
) -> bool:
    if _span_distance(value_match, keyword_match) > 60:
        return False
    cue = keyword_match.group(0).casefold()
    if value_match.end() <= keyword_match.start() and re.match(
        r"(?:not\s|up\s|(?:to\s+)?a\s+maximum)", cue
    ):
        return False
    if keyword_match.end() <= value_match.start():
        bridge = text[keyword_match.end() : value_match.start()]
    else:
        bridge = text[value_match.end() : keyword_match.start()]
    if re.search(_NUMBER, bridge):
        return False
    return not bool(
        re.search(
            r"[.;!?]|\b(?:AUC|C\s*max|concentration|exposure|T\s*max)\b",
            bridge,
            re.IGNORECASE,
        )
    )


def _bioavailability_directly_linked(
    text: str,
    value_match: re.Match[str],
    keyword_match: re.Match[str],
) -> bool:
    distance = _span_distance(value_match, keyword_match)
    if distance > 100:
        return False
    if keyword_match.end() <= value_match.start():
        bridge = text[keyword_match.end() : value_match.start()]
    else:
        bridge = text[value_match.end() : keyword_match.start()]
    if re.search(
        r"[.;!?%]|\b(?:AUC|clearance|C\s*max|half[- ]life|T\s*max)\b",
        bridge,
        re.IGNORECASE,
    ):
        return False
    if value_match.group("unit"):
        return len(bridge.split()) <= 12
    return bool(
        re.fullmatch(
            r"\s*(?:(?:is|was|of|=|:)\s*)?(?:approximately\s*|about\s*)?",
            bridge,
            re.IGNORECASE,
        )
        or (
            keyword_match.group(0).casefold() == "absbio"
            and len(bridge.split()) <= 3
        )
    )


def _candidate_row(
    *,
    candidate_type: str,
    measurement_context: str,
    payload: dict[str, Any],
    record: dict[str, Any],
    cache_file: str,
    cache_file_sha256: str,
    section: str,
    section_item_index: int,
    section_text: str,
    value_match: re.Match[str],
    keyword_match: re.Match[str],
    normalized_value: float | None,
    normalized_unit: str,
    normalization_status: str,
    normalization_reasons: Collection[str],
    extra_reasons: Collection[str] = (),
    excerpt_chars: int = DEFAULT_EXCERPT_CHARS,
) -> dict[str, Any]:
    anchor_start = min(value_match.start(), keyword_match.start())
    anchor_end = max(value_match.end(), keyword_match.end())
    context = _local_context(section_text, anchor_start, anchor_end)
    drug_id = _clean_text(payload.get("drug_id"))
    analyte, analyte_status, generic_names, substance_names, analyte_reasons = (
        _detect_analyte(context, record, drug_id)
    )
    route, route_source = _detect_route(context, record)
    formulation, formulation_source = _detect_formulation(context, record)
    population = _detect_population(context)
    regimen, steady_state = _regimen_context(context)
    reasons = [BASE_EXCLUSION, *normalization_reasons, *analyte_reasons, *extra_reasons]
    if _NON_HUMAN_RE.search(context):
        reasons.append("non_human_context")
    reasons = _unique(reasons)

    set_id = _clean_text(record.get("set_id"))
    document_id = _clean_text(record.get("id"))
    source_record_id = set_id or document_id
    spl_version = _clean_text(record.get("version"))
    excerpt = _short_excerpt(section_text, anchor_start, anchor_end, excerpt_chars)
    anchor_fingerprint = _clean_text(
        section_text[max(0, anchor_start - 100) : min(len(section_text), anchor_end + 100)]
    ).casefold()
    identity = "|".join(
        (
            drug_id.casefold(),
            source_record_id,
            spl_version,
            candidate_type,
            _clean_text(value_match.group("value_expression")),
            _clean_text(value_match.group("unit")),
            anchor_fingerprint,
        )
    )
    candidate_id = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:24]
    confidence = (
        "medium"
        if normalized_value is not None
        and not analyte_reasons
        and "non_human_context" not in reasons
        and not extra_reasons
        else "low"
    )
    return {
        "candidate_id": candidate_id,
        "candidate_type": candidate_type,
        "candidate_status": "manual_review_required",
        "drug_id": drug_id,
        "cache_file": cache_file,
        "cache_file_sha256": cache_file_sha256,
        "cache_search": _clean_text(payload.get("search")),
        "cache_status": _clean_text(payload.get("status")),
        "cache_last_updated": _clean_text(payload.get("last_updated")),
        "source_name": SOURCE_NAME,
        "source_record_id": source_record_id,
        "source_set_id": set_id,
        "source_document_id": document_id,
        "spl_version": spl_version,
        "effective_time": _clean_text(record.get("effective_time")),
        "source_url": (
            f"https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={set_id}"
            if set_id
            else ""
        ),
        "section": section,
        "section_item_index": section_item_index,
        "section_text_sha256": hashlib.sha256(section_text.encode("utf-8")).hexdigest(),
        "match_start": value_match.start(),
        "match_end": value_match.end(),
        "source_excerpt": excerpt,
        "raw_value": _clean_text(value_match.group("value_expression")),
        "raw_unit": _clean_text(value_match.group("unit")),
        "normalized_value": normalized_value,
        "normalized_unit": normalized_unit,
        "normalization_status": normalization_status,
        "measurement_context": measurement_context,
        "analyte": analyte,
        "analyte_status": analyte_status,
        "label_generic_names": _join(generic_names),
        "label_substance_names": _join(substance_names),
        "route": route,
        "route_source": route_source,
        "formulation": formulation,
        "formulation_source": formulation_source,
        "population": population,
        "steady_state": steady_state,
        "regimen": regimen,
        "confidence": confidence,
        "review_required": True,
        "exclusion_reason": ";".join(reasons),
        "model_ready": False,
        "extraction_method": SCHEMA_VERSION,
        "_dedupe_key": candidate_id,
    }


def _section_items(record: dict[str, Any], section: str) -> Iterator[tuple[int, str]]:
    for index, value in enumerate(_values(record.get(section))):
        text = _clean_text(value)
        if text:
            yield index, text


def _extract_clearance(
    *,
    text: str,
    section: str,
    section_item_index: int,
    payload: dict[str, Any],
    record: dict[str, Any],
    cache_file: str,
    cache_file_sha256: str,
    excerpt_chars: int,
) -> Iterator[dict[str, Any]]:
    for value_match in _CLEARANCE_VALUE_RE.finditer(text):
        keyword_match = _nearest_match(
            text,
            value_match,
            _CLEARANCE_KEYWORD_RE,
            max_distance=180,
        )
        if keyword_match is None:
            continue
        value, unit, status, reasons = _normalize_clearance(value_match)
        anchor_start = min(value_match.start(), keyword_match.start())
        anchor_end = max(value_match.end(), keyword_match.end())
        context = _local_context(text, anchor_start, anchor_end)
        yield _candidate_row(
            candidate_type="clearance",
            measurement_context=_clearance_context(context),
            payload=payload,
            record=record,
            cache_file=cache_file,
            cache_file_sha256=cache_file_sha256,
            section=section,
            section_item_index=section_item_index,
            section_text=text,
            value_match=value_match,
            keyword_match=keyword_match,
            normalized_value=value,
            normalized_unit=unit,
            normalization_status=status,
            normalization_reasons=reasons,
            excerpt_chars=excerpt_chars,
        )


def _extract_bioavailability(
    *,
    text: str,
    section: str,
    section_item_index: int,
    payload: dict[str, Any],
    record: dict[str, Any],
    cache_file: str,
    cache_file_sha256: str,
    excerpt_chars: int,
) -> Iterator[dict[str, Any]]:
    for value_match in _BIOAVAILABILITY_VALUE_RE.finditer(text):
        keyword_match = _nearest_match(
            text,
            value_match,
            _ABSOLUTE_BIOAVAILABILITY_RE,
            max_distance=120,
        )
        if keyword_match is None or not _bioavailability_directly_linked(
            text, value_match, keyword_match
        ):
            continue
        value, unit, status, reasons = _normalize_bioavailability(value_match)
        anchor_start = min(value_match.start(), keyword_match.start())
        anchor_end = max(value_match.end(), keyword_match.end())
        context = _local_context(text, anchor_start, anchor_end)
        extra_reasons = (
            ["absolute_bioavailability_stated_unknown"]
            if _UNKNOWN_BIOAVAILABILITY_RE.search(context)
            else []
        )
        if extra_reasons:
            value = None
            unit = ""
            status = "withheld_source_states_unknown"
        yield _candidate_row(
            candidate_type="absolute_bioavailability",
            measurement_context="absolute_bioavailability",
            payload=payload,
            record=record,
            cache_file=cache_file,
            cache_file_sha256=cache_file_sha256,
            section=section,
            section_item_index=section_item_index,
            section_text=text,
            value_match=value_match,
            keyword_match=keyword_match,
            normalized_value=value,
            normalized_unit=unit,
            normalization_status=status,
            normalization_reasons=reasons,
            extra_reasons=extra_reasons,
            excerpt_chars=excerpt_chars,
        )


def _extract_observed_doses(
    *,
    text: str,
    section: str,
    section_item_index: int,
    payload: dict[str, Any],
    record: dict[str, Any],
    cache_file: str,
    cache_file_sha256: str,
    excerpt_chars: int,
) -> Iterator[dict[str, Any]]:
    for value_match in _DOSE_VALUE_RE.finditer(text):
        keyword_match = _nearest_match(
            text,
            value_match,
            _DOSE_CUE_RE,
            max_distance=120,
        )
        if keyword_match is None:
            continue
        strong_context = _local_context(text, value_match.start(), value_match.end(), 140)
        extra_reasons = (
            []
            if _STRONG_OBSERVED_DOSE_RE.search(strong_context)
            else ["observed_administration_not_explicit"]
        )
        value, unit, status, reasons = _normalize_dose(value_match)
        yield _candidate_row(
            candidate_type="observed_pk_dose",
            measurement_context="observed_pk_study_dose",
            payload=payload,
            record=record,
            cache_file=cache_file,
            cache_file_sha256=cache_file_sha256,
            section=section,
            section_item_index=section_item_index,
            section_text=text,
            value_match=value_match,
            keyword_match=keyword_match,
            normalized_value=value,
            normalized_unit=unit,
            normalization_status=status,
            normalization_reasons=reasons,
            extra_reasons=extra_reasons,
            excerpt_chars=excerpt_chars,
        )


def _extract_maximum_doses(
    *,
    text: str,
    section: str,
    section_item_index: int,
    payload: dict[str, Any],
    record: dict[str, Any],
    cache_file: str,
    cache_file_sha256: str,
    excerpt_chars: int,
) -> Iterator[dict[str, Any]]:
    value_matches = list(_DOSE_VALUE_RE.finditer(text))
    for keyword_match in _MAXIMUM_DOSE_RE.finditer(text):
        keyword_text = keyword_match.group(0)
        if section != "dosage_and_administration" and not re.search(
            r"\b(?:recommended|approved|labeled|labelled)\b",
            keyword_text,
            re.IGNORECASE,
        ):
            continue
        nearby = [
            value_match
            for value_match in value_matches
            if _maximum_dose_directly_linked(text, value_match, keyword_match)
        ]
        if not nearby:
            continue
        value_match = min(
            nearby,
            key=lambda item: (
                _span_distance(item, keyword_match),
                abs(item.start() - keyword_match.end()),
            ),
        )
        value, unit, status, reasons = _normalize_dose(value_match)
        anchor_start = min(value_match.start(), keyword_match.start())
        anchor_end = max(value_match.end(), keyword_match.end())
        context = _local_context(text, anchor_start, anchor_end, 80)
        yield _candidate_row(
            candidate_type="maximum_labeled_dose",
            measurement_context=_maximum_dose_context(context),
            payload=payload,
            record=record,
            cache_file=cache_file,
            cache_file_sha256=cache_file_sha256,
            section=section,
            section_item_index=section_item_index,
            section_text=text,
            value_match=value_match,
            keyword_match=keyword_match,
            normalized_value=value,
            normalized_unit=unit,
            normalization_status=status,
            normalization_reasons=reasons,
            excerpt_chars=excerpt_chars,
        )


def extract_spl_record_candidates(
    record: Mapping[str, Any],
    *,
    payload: Mapping[str, Any] | None = None,
    cache_file: str | Path = "",
    cache_file_sha256: str = "",
    excerpt_chars: int = DEFAULT_EXCERPT_CHARS,
) -> list[dict[str, Any]]:
    """Extract review-only numeric PK candidates from one cached SPL record."""
    record_dict = dict(record)
    payload_dict = dict(payload or {})
    cache_path = str(cache_file)
    candidates: list[dict[str, Any]] = []
    seen: set[str] = set()

    for section in _PK_SECTIONS:
        for item_index, text in _section_items(record_dict, section):
            extractors = (
                _extract_clearance,
                _extract_bioavailability,
                _extract_observed_doses,
            )
            for extractor in extractors:
                for row in extractor(
                    text=text,
                    section=section,
                    section_item_index=item_index,
                    payload=payload_dict,
                    record=record_dict,
                    cache_file=cache_path,
                    cache_file_sha256=cache_file_sha256,
                    excerpt_chars=excerpt_chars,
                ):
                    key = str(row.pop("_dedupe_key"))
                    if key not in seen:
                        candidates.append(row)
                        seen.add(key)

    for section in _MAXIMUM_DOSE_SECTIONS:
        for item_index, text in _section_items(record_dict, section):
            for row in _extract_maximum_doses(
                text=text,
                section=section,
                section_item_index=item_index,
                payload=payload_dict,
                record=record_dict,
                cache_file=cache_path,
                cache_file_sha256=cache_file_sha256,
                excerpt_chars=excerpt_chars,
            ):
                key = str(row.pop("_dedupe_key"))
                if key not in seen:
                    candidates.append(row)
                    seen.add(key)
    return candidates


def _records(payload: dict[str, Any]) -> tuple[list[dict[str, Any]], bool]:
    raw_records = payload.get("records")
    if raw_records is None:
        return [], False
    if isinstance(raw_records, dict):
        return [raw_records], True
    if not isinstance(raw_records, list):
        return [], True
    return [record for record in raw_records if isinstance(record, dict)], False


def _write_candidates(path: Path, rows: Collection[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def review_spl_pk_candidates(
    *,
    cache_dir: str | Path,
    out_dir: str | Path,
    drug_ids: Collection[str] | None = None,
    max_cache_files: int = 0,
    excerpt_chars: int = DEFAULT_EXCERPT_CHARS,
) -> dict[str, Any]:
    """Write a provenance-first, manual-review-only SPL PK candidate queue."""
    if max_cache_files < 0:
        raise ValueError("max_cache_files must be zero or greater")
    excerpt_chars = max(MIN_EXCERPT_CHARS, min(MAX_EXCERPT_CHARS, excerpt_chars))
    cache = Path(cache_dir)
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    cache_files = sorted(cache.glob("*.json")) if cache.is_dir() else []
    requested_drugs = {str(value).strip().casefold() for value in drug_ids or [] if str(value).strip()}

    rows: list[dict[str, Any]] = []
    parse_failure_files: list[str] = []
    invalid_payload_files: list[str] = []
    cache_files_read = 0
    cache_files_selected = 0
    records_seen = 0
    record_containers_repaired = 0

    for cache_file in cache_files:
        if max_cache_files and cache_files_selected >= max_cache_files:
            break
        try:
            raw_payload = cache_file.read_bytes()
            payload_value = json.loads(raw_payload)
        except (OSError, json.JSONDecodeError, UnicodeDecodeError):
            parse_failure_files.append(str(cache_file))
            continue
        cache_files_read += 1
        if not isinstance(payload_value, dict):
            invalid_payload_files.append(str(cache_file))
            continue
        payload = payload_value
        drug_id = _clean_text(payload.get("drug_id")).casefold()
        if requested_drugs and drug_id not in requested_drugs:
            continue
        cache_files_selected += 1
        payload_records, repaired = _records(payload)
        record_containers_repaired += int(repaired)
        checksum = hashlib.sha256(raw_payload).hexdigest()
        for record in payload_records:
            records_seen += 1
            rows.extend(
                extract_spl_record_candidates(
                    record,
                    payload=payload,
                    cache_file=cache_file,
                    cache_file_sha256=checksum,
                    excerpt_chars=excerpt_chars,
                )
            )

    candidates_path = output / "spl_pk_candidate_review.csv"
    manifest_path = output / "spl_pk_candidate_review_manifest.json"
    _write_candidates(candidates_path, rows)

    type_counts = Counter(str(row["candidate_type"]) for row in rows)
    normalized_counts = Counter(
        str(row["candidate_type"])
        for row in rows
        if row.get("normalized_value") is not None
    )
    confidence_counts = Counter(str(row["confidence"]) for row in rows)
    additional_exclusion_counts: Counter[str] = Counter()
    for row in rows:
        for reason in str(row["exclusion_reason"]).split(";"):
            if reason and reason != BASE_EXCLUSION:
                additional_exclusion_counts[reason] += 1

    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_name": SOURCE_NAME,
        "cache_dir": str(cache),
        "cache_files_discovered": len(cache_files),
        "cache_files_read": cache_files_read,
        "cache_files_selected": cache_files_selected,
        "cache_parse_failures": len(parse_failure_files),
        "cache_parse_failure_files": parse_failure_files,
        "invalid_payloads": len(invalid_payload_files),
        "invalid_payload_files": invalid_payload_files,
        "record_containers_repaired": record_containers_repaired,
        "records_reviewed": records_seen,
        "records_with_candidates": len(
            {
                (row["cache_file_sha256"], row["source_record_id"], row["spl_version"])
                for row in rows
            }
        ),
        "drug_filter": sorted(requested_drugs),
        "max_cache_files": max_cache_files,
        "excerpt_character_limit": excerpt_chars,
        "candidate_output": str(candidates_path),
        "manifest_output": str(manifest_path),
        "candidate_rows": len(rows),
        "candidate_rows_by_type": dict(sorted(type_counts.items())),
        "normalized_rows_by_type": dict(sorted(normalized_counts.items())),
        "confidence_counts": dict(sorted(confidence_counts.items())),
        "review_required_rows": sum(bool(row["review_required"]) for row in rows),
        "model_ready_rows": sum(bool(row["model_ready"]) for row in rows),
        "rows_with_additional_exclusions": sum(
            str(row["exclusion_reason"]) != BASE_EXCLUSION for row in rows
        ),
        "additional_exclusion_counts": dict(sorted(additional_exclusion_counts.items())),
        "policy": {
            "candidate_use": "manual_review_queue_only",
            "automatic_truth_promotion": "prohibited",
            "observed_pk_dose": "extracted_only_from_PK_or_clinical_pharmacology_sections",
            "maximum_labeled_dose": "separate_candidate_type_requiring_explicit_maximum_language",
            "cmax": "not_extracted_or_linearly_rescaled",
            "normalization": "only_singular_values_with_unambiguous_supported_units",
            "confidence": "medium_is_still_review_required; no_candidate_is_high_confidence",
        },
        "provenance_note": (
            "match offsets refer to whitespace-normalized section text; source hashes and short "
            "excerpts are retained for SPL-version-specific adjudication"
        ),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest
