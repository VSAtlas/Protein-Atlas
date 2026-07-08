from __future__ import annotations

import re
from collections.abc import Iterable as IterableABC
from typing import Any, Dict, Iterable, Mapping, Sequence

from analysis.reporting.value_utils import normalize_side_effect_label


_LIST_SPLIT_RE = re.compile(r"[;\n|]+")

_SIDE_EFFECT_FIELDS = (
    "side_effects",
    "known_side_effect_examples",
    "side_effect_examples",
    "raw_adverse_event_examples",
    "adverse_event_examples",
    "direct_liability_examples",
    "adverse_reactions",
    "openfda_adverse_reactions",
    "fda_adverse_reactions",
    "sider_side_effects",
)
_SAFETY_BUCKET_FIELDS = ("safety_buckets", "drug_ae_buckets")
_DIRECT_SAFETY_FIELDS = ("direct_safety_buckets",)
_SECONDARY_SAFETY_FIELDS = ("secondary_safety_buckets",)
_SOURCE_FIELDS = (
    "known_side_effect_sources",
    "side_effect_sources",
    "adverse_event_sources",
    "sider_sources",
    "openfda_sources",
)
_SUMMARY_FIELDS = (
    "known_side_effect_summary",
    "side_effect_evidence_summary",
    "adverse_event_summary",
)


def _dedupe(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        clean = normalize_side_effect_label(value)
        key = clean.casefold()
        if not clean or key in seen:
            continue
        seen.add(key)
        out.append(clean)
    return out


def _list_values(source: Mapping[str, Any], fields: Sequence[str]) -> list[str]:
    values: list[str] = []
    for field in fields:
        raw = source.get(field)
        if raw is None:
            continue
        if isinstance(raw, str):
            values.extend(part for part in _LIST_SPLIT_RE.split(raw) if part.strip())
        elif isinstance(raw, IterableABC):
            values.extend(
                str(part) for part in raw if normalize_side_effect_label(part)
            )
        else:
            values.append(str(raw))
    return _dedupe(values)


def _first_text(source: Mapping[str, Any], fields: Sequence[str]) -> str:
    for field in fields:
        text = normalize_side_effect_label(source.get(field))
        if text:
            return text
    return ""


def _first_number(source: Mapping[str, Any], fields: Sequence[str]) -> float | None:
    for field in fields:
        raw = source.get(field)
        if raw in (None, ""):
            continue
        try:
            value = float(str(raw).strip())
        except (TypeError, ValueError):
            continue
        if value == value:
            return value
    return None


def collect_ligand_side_effect_fields(source: Mapping[str, Any]) -> Dict[str, Any]:
    """Collect optional ligand safety and exposure fields for report row metadata."""

    out: Dict[str, Any] = {}
    side_effects = _list_values(source, _SIDE_EFFECT_FIELDS)
    if side_effects:
        out["side_effects"] = side_effects

    safety_buckets = _list_values(source, _SAFETY_BUCKET_FIELDS)
    if safety_buckets:
        out["safety_buckets"] = safety_buckets
        out["primary_display_safety"] = normalize_side_effect_label(
            source.get("primary_display_safety")
        ) or safety_buckets[0]

    secondary_buckets = _list_values(source, _SECONDARY_SAFETY_FIELDS)
    if secondary_buckets:
        out["secondary_safety_buckets"] = secondary_buckets

    direct_buckets = _list_values(source, _DIRECT_SAFETY_FIELDS)
    if direct_buckets:
        out["direct_safety_buckets"] = direct_buckets

    sources = _list_values(source, _SOURCE_FIELDS)
    if sources:
        out["side_effect_sources"] = sources

    summary = _first_text(source, _SUMMARY_FIELDS)
    if summary:
        out["side_effect_evidence_summary"] = summary

    free_cmax_um = _first_number(
        source, ("free_cmax_um", "cmax_free_um", "unbound_cmax_um")
    )
    if free_cmax_um is None:
        total_cmax_um = _first_number(source, ("cmax_um", "total_cmax_um"))
        fraction_unbound = _first_number(
            source,
            (
                "fraction_unbound_plasma",
                "fu_plasma",
                "plasma_unbound_fraction",
                "unbound_fraction",
            ),
        )
        if total_cmax_um is not None and fraction_unbound is not None:
            free_cmax_um = total_cmax_um * fraction_unbound
            out["cmax_um"] = total_cmax_um
            out["fraction_unbound_plasma"] = fraction_unbound
    if free_cmax_um is not None:
        out["free_cmax_um"] = free_cmax_um

    exposure_source = _first_text(
        source, ("free_cmax_source", "cmax_source", "exposure_source")
    )
    if exposure_source:
        out["exposure_source"] = exposure_source

    return out


def ligand_side_effect_report_fields(source: Mapping[str, Any]) -> Dict[str, Any]:
    """Return normalized side-effect fields expected by the report client."""

    fields = collect_ligand_side_effect_fields(source)
    return {
        "side_effects": fields.get("side_effects", []),
        "safety_buckets": fields.get("safety_buckets", []),
        "primary_display_safety": fields.get("primary_display_safety", ""),
        "secondary_safety_buckets": fields.get("secondary_safety_buckets", []),
        "direct_safety_buckets": fields.get("direct_safety_buckets", []),
        "side_effect_sources": fields.get("side_effect_sources", []),
        "side_effect_evidence_summary": fields.get("side_effect_evidence_summary", ""),
        "free_cmax_um": fields.get("free_cmax_um", ""),
        "cmax_um": fields.get("cmax_um", ""),
        "fraction_unbound_plasma": fields.get("fraction_unbound_plasma", ""),
        "exposure_source": fields.get("exposure_source", ""),
    }


def target_side_effect_report_fields(source: Mapping[str, Any]) -> Dict[str, Any]:
    """Return target side-effect fields expected by the report client."""

    side_effects = _list_values(source, _SIDE_EFFECT_FIELDS)
    return {
        "side_effects": side_effects,
        "safety_evidence_terms": _list_values(source, ("safety_evidence_terms",)),
        "direct_liability_examples": _list_values(source, ("direct_liability_examples",)),
        "raw_adverse_event_examples": _list_values(source, ("raw_adverse_event_examples",)),
    }
