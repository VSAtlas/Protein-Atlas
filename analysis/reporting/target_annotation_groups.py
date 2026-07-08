from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Set, cast

import yaml  # type: ignore[import-untyped]

from analysis.reporting.target_safety_evidence import (
    HEURISTIC_CONFIDENCE,
    HEURISTIC_SOURCE,
    MANUAL_OVERRIDE_CONFIDENCE,
    MANUAL_OVERRIDE_SOURCE,
    SAFETY_DEFAULT,
    load_target_safety_profile,
)
from analysis.reporting.target_expression_cache import lookup_target_expression_entry
from analysis.reporting.value_utils import normalize_text

ADME_DEFAULT = "Non-ADME / Unassigned"


def _normalize_text(value: Any) -> str:
    return normalize_text(value)


def _normalize_upper(value: Any) -> str:
    return _normalize_text(value).upper()


def _normalize_lower(value: Any) -> str:
    return _normalize_text(value).lower()


@lru_cache(maxsize=1)
def _load_catalog(catalog_path: str) -> Dict[str, Any]:
    path = Path(catalog_path)
    if not path.exists():
        path = Path(__file__).resolve().parents[2] / "chemdb" / "target_organization_annotations.yaml"
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return cast(Dict[str, Any], payload)


@lru_cache(maxsize=4)
def _load_override_annotations(repo_root: Path) -> Dict[str, Any]:
    override_path = repo_root / "pathways" / "cache" / "target_annotation_overrides.yaml"
    if not override_path.exists():
        return {}
    payload = yaml.safe_load(override_path.read_text(encoding="utf-8")) or {}
    return cast(Dict[str, Any], payload)


def _catalog_section(catalog: Dict[str, Any], section_name: str) -> Dict[str, Any]:
    section = catalog.get(section_name)
    return cast(Dict[str, Any], section if isinstance(section, dict) else {})


def _section_order(section: Dict[str, Any], default_value: str) -> List[str]:
    order = [str(item).strip() for item in cast(Sequence[Any], section.get("order") or [])]
    order = [item for item in order if item]
    if default_value not in order:
        order.append(default_value)
    return order


def _rule_matches(
    rule: Dict[str, Any],
    uniprots: Set[str],
    gene_symbols: Set[str],
    search_text: str,
) -> bool:
    rule_uniprots = {
        _normalize_upper(value)
        for value in cast(Sequence[Any], rule.get("uniprot_ids") or [])
        if _normalize_text(value)
    }
    if rule_uniprots and (rule_uniprots & uniprots):
        return True
    rule_genes = {
        _normalize_upper(value)
        for value in cast(Sequence[Any], rule.get("gene_symbols") or [])
        if _normalize_text(value)
    }
    if rule_genes and (rule_genes & gene_symbols):
        return True
    keywords = [
        _normalize_lower(value)
        for value in cast(Sequence[Any], rule.get("target_name_keywords") or [])
        if _normalize_text(value)
    ]
    return any(keyword in search_text for keyword in keywords)


def _apply_exact_override(
    override_section: Dict[str, Any], uniprots: Set[str], gene_symbols: Set[str]
) -> Dict[str, Any]:
    uniprot_map = cast(Dict[str, Any], override_section.get("by_uniprot") or {})
    for accession in sorted(uniprots):
        match = uniprot_map.get(accession)
        if isinstance(match, dict):
            return cast(Dict[str, Any], match)
    gene_map = cast(Dict[str, Any], override_section.get("by_gene") or {})
    for symbol in sorted(gene_symbols):
        match = gene_map.get(symbol)
        if isinstance(match, dict):
            return cast(Dict[str, Any], match)
    return {}


def _first_present(source: Dict[str, Any], *keys: str) -> Any:
    for key in keys:
        value = source.get(key)
        if value not in (None, ""):
            return value
    return None


def _normalize_expression_profile(source: Dict[str, Any]) -> Dict[str, Any]:
    level = _normalize_text(
        _first_present(
            source, "target_tissue_expression", "tissue_expression_level", "level"
        )
    )
    score_raw = _first_present(
        source, "target_tissue_expression_score", "tissue_expression_score", "score"
    )
    score = ""
    if score_raw not in (None, ""):
        try:
            score_float = float(str(score_raw).strip())
        except (TypeError, ValueError):
            score_float = None
        if score_float is not None and score_float == score_float:
            score = f"{max(0.0, min(1.0, score_float)):.3f}"
    source_text = _normalize_text(
        source.get("target_tissue_expression_source")
        or source.get("source")
        or source.get("expression_source")
    )
    summary = _normalize_text(
        source.get("target_tissue_expression_summary")
        or source.get("summary")
        or source.get("expression_summary")
    )
    out: Dict[str, Any] = {}
    if level:
        out["target_tissue_expression"] = level
        out["target_tissue_expression_label"] = level
    if score:
        out["target_tissue_expression_score"] = score
    if source_text:
        out["target_tissue_expression_source"] = source_text
    if summary:
        out["target_tissue_expression_summary"] = summary
    return out


def _resolve_expression_profile(
    *,
    repo_root: Path,
    catalog: Dict[str, Any],
    overrides: Dict[str, Any],
    uniprots: Set[str],
    gene_symbols: Set[str],
    search_text: str,
) -> Dict[str, Any]:
    override = _apply_exact_override(
        _catalog_section(overrides, "expression"), uniprots, gene_symbols
    )
    if override:
        return _normalize_expression_profile(override)

    cached = lookup_target_expression_entry(
        repo_root, uniprots=uniprots, gene_symbols=gene_symbols
    )
    if cached:
        return _normalize_expression_profile(cached)

    expression_section = _catalog_section(catalog, "expression")
    for raw_rule in cast(Sequence[Any], expression_section.get("rules") or []):
        rule = cast(Dict[str, Any], raw_rule if isinstance(raw_rule, dict) else {})
        if _rule_matches(rule, uniprots, gene_symbols, search_text):
            return _normalize_expression_profile(rule)
    return {}


def _match_categories(
    section: Dict[str, Any],
    order: Sequence[str],
    uniprots: Set[str],
    gene_symbols: Set[str],
    search_text: str,
) -> List[str]:
    matched: List[str] = []
    rules = cast(Sequence[Any], section.get("rules") or [])
    for category in order:
        for raw_rule in rules:
            rule = cast(Dict[str, Any], raw_rule if isinstance(raw_rule, dict) else {})
            if _normalize_text(rule.get("category")) != category:
                continue
            if _rule_matches(rule, uniprots, gene_symbols, search_text):
                matched.append(category)
                break
    return matched


def resolve_target_annotations(
    repo_root: Path,
    *,
    target_name: str,
    uniprots: Iterable[str],
    gene_symbols: Iterable[str] | None = None,
) -> Dict[str, Any]:
    catalog = _load_catalog(
        str(repo_root / "chemdb" / "target_organization_annotations.yaml")
    )
    overrides = _load_override_annotations(repo_root)
    adme_section = _catalog_section(catalog, "adme")
    safety_section = _catalog_section(catalog, "safety")
    adme_order = _section_order(adme_section, ADME_DEFAULT)
    safety_order = _section_order(safety_section, SAFETY_DEFAULT)
    uniprot_set = {_normalize_upper(value) for value in uniprots if _normalize_text(value)}
    gene_symbol_set = {
        _normalize_upper(value) for value in (gene_symbols or []) if _normalize_text(value)
    }
    search_text = _normalize_lower(target_name)

    adme_override = _apply_exact_override(
        _catalog_section(overrides, "adme"), uniprot_set, gene_symbol_set
    )
    if _normalize_text(adme_override.get("category")):
        adme_category = _normalize_text(adme_override.get("category"))
    else:
        adme_matches = _match_categories(
            adme_section,
            adme_order,
            uniprot_set,
            gene_symbol_set,
            search_text,
        )
        adme_category = adme_matches[0] if adme_matches else ADME_DEFAULT

    safety_override = _apply_exact_override(
        _catalog_section(overrides, "safety"), uniprot_set, gene_symbol_set
    )
    safety_matches = _match_categories(
        safety_section,
        safety_order,
        uniprot_set,
        gene_symbol_set,
        search_text,
    )
    override_buckets = [
        _normalize_text(value)
        for value in cast(Sequence[Any], safety_override.get("buckets") or [])
        if _normalize_text(value)
    ]
    if override_buckets:
        safety_buckets = [
            bucket for bucket in safety_order if bucket in set(override_buckets)
        ] or [SAFETY_DEFAULT]
        primary_display_safety = safety_buckets[0] if safety_buckets else SAFETY_DEFAULT
        secondary_safety_buckets = safety_buckets[1:]
        direct_safety_buckets = []
        drug_ae_buckets = []
        safety_confidence = _normalize_text(safety_override.get("confidence")) or MANUAL_OVERRIDE_CONFIDENCE
        safety_sources = [
            _normalize_text(value)
            for value in cast(Sequence[Any], safety_override.get("sources") or [MANUAL_OVERRIDE_SOURCE])
            if _normalize_text(value)
        ] or [MANUAL_OVERRIDE_SOURCE]
        safety_bucket_scores = {
            bucket: "1.000"
            for bucket in safety_buckets
            if bucket != SAFETY_DEFAULT
        }
        safety_evidence_summary = _normalize_text(
            safety_override.get("summary")
        ) or "Manual override"
        direct_liability_examples = []
        raw_adverse_event_examples = []
        safety_evidence_terms = []
    else:
        cached_profile = load_target_safety_profile(
            repo_root,
            target_name=target_name,
            uniprots=uniprot_set,
            gene_symbols=gene_symbol_set,
        )
        if cached_profile:
            safety_buckets = [
                _normalize_text(value)
                for value in cast(Sequence[Any], cached_profile.get("safety_buckets") or [])
                if _normalize_text(value)
            ] or [SAFETY_DEFAULT]
            secondary_safety_buckets = [
                _normalize_text(value)
                for value in cast(
                    Sequence[Any], cached_profile.get("secondary_safety_buckets") or []
                )
                if _normalize_text(value)
            ]
            direct_safety_buckets = [
                _normalize_text(value)
                for value in cast(
                    Sequence[Any], cached_profile.get("direct_safety_buckets") or []
                )
                if _normalize_text(value)
            ]
            drug_ae_buckets = [
                _normalize_text(value)
                for value in cast(
                    Sequence[Any], cached_profile.get("drug_ae_buckets") or []
                )
                if _normalize_text(value)
            ]
            primary_display_safety = (
                _normalize_text(cached_profile.get("primary_display_safety"))
                or safety_buckets[0]
            )
            if safety_matches:
                preferred = safety_matches[0]
                if not direct_safety_buckets and preferred != SAFETY_DEFAULT:
                    primary_display_safety = preferred
                    if preferred not in safety_buckets:
                        safety_buckets = [preferred] + [
                            bucket for bucket in safety_buckets if bucket != preferred
                        ]
                elif preferred in direct_safety_buckets and primary_display_safety not in direct_safety_buckets:
                    primary_display_safety = preferred
            secondary_safety_buckets = [
                bucket
                for bucket in secondary_safety_buckets
                if bucket and bucket != primary_display_safety
            ]
            if primary_display_safety != SAFETY_DEFAULT:
                safety_buckets = [primary_display_safety] + [
                    bucket
                    for bucket in safety_buckets
                    if bucket
                    and bucket != primary_display_safety
                    and bucket != SAFETY_DEFAULT
                ]
            safety_confidence = (
                _normalize_text(cached_profile.get("safety_confidence")) or "unassigned"
            )
            safety_sources = [
                _normalize_text(value)
                for value in cast(Sequence[Any], cached_profile.get("safety_sources") or [])
                if _normalize_text(value)
            ]
            safety_bucket_scores = {
                _normalize_text(bucket): _normalize_text(score)
                for bucket, score in cast(
                    Dict[str, Any], cached_profile.get("safety_bucket_scores") or {}
                ).items()
                if _normalize_text(bucket) and _normalize_text(score)
            }
            safety_evidence_summary = _normalize_text(
                cached_profile.get("safety_evidence_summary")
            )
            safety_evidence_terms = [
                _normalize_text(value)
                for value in cast(
                    Sequence[Any], cached_profile.get("safety_evidence_terms") or []
                )
                if _normalize_text(value)
            ]
            direct_liability_examples = [
                _normalize_text(value)
                for value in cast(
                    Sequence[Any], cached_profile.get("direct_liability_examples") or []
                )
                if _normalize_text(value)
            ]
            raw_adverse_event_examples = [
                _normalize_text(value)
                for value in cast(
                    Sequence[Any], cached_profile.get("raw_adverse_event_examples") or []
                )
                if _normalize_text(value)
            ]
        else:
            safety_buckets = list(safety_matches)
            if not safety_buckets:
                safety_buckets = [SAFETY_DEFAULT]
            primary_display_safety = safety_buckets[0] if safety_buckets else SAFETY_DEFAULT
            secondary_safety_buckets = safety_buckets[1:]
            direct_safety_buckets = []
            drug_ae_buckets = []
            safety_confidence = HEURISTIC_CONFIDENCE
            safety_sources = [HEURISTIC_SOURCE]
            safety_bucket_scores = {
                bucket: "1.000" for bucket in safety_buckets if bucket != SAFETY_DEFAULT
            }
            safety_evidence_summary = "Local heuristic catalog"
            safety_evidence_terms = []
            direct_liability_examples = []
            raw_adverse_event_examples = []

    expression_profile = _resolve_expression_profile(
        repo_root=repo_root,
        catalog=catalog,
        overrides=overrides,
        uniprots=uniprot_set,
        gene_symbols=gene_symbol_set,
        search_text=search_text,
    )

    return {
        "adme_category": adme_category,
        "safety_buckets": safety_buckets,
        "primary_display_safety": primary_display_safety,
        "secondary_safety_buckets": secondary_safety_buckets,
        "direct_safety_buckets": direct_safety_buckets,
        "drug_ae_buckets": drug_ae_buckets,
        "safety_confidence": safety_confidence,
        "safety_sources": safety_sources,
        "safety_bucket_scores": safety_bucket_scores,
        "safety_evidence_summary": safety_evidence_summary,
        "safety_evidence_terms": safety_evidence_terms,
        "direct_liability_examples": direct_liability_examples,
        "raw_adverse_event_examples": raw_adverse_event_examples,
        **expression_profile,
    }


def target_annotation_block_summary(
    *,
    ordered_labels: Sequence[str],
    target_meta: Dict[str, Dict[str, Any]],
    field_name: str,
    order: Sequence[str],
    default_value: str,
) -> List[Dict[str, Any]]:
    counts: Dict[str, int] = {}
    for label in ordered_labels:
        value = _normalize_text(target_meta.get(label, {}).get(field_name)) or default_value
        counts[value] = counts.get(value, 0) + 1
    order_map = {name: idx for idx, name in enumerate(order)}
    return [
        {"key": name, "label": name, "count": count}
        for name, count in sorted(
            counts.items(),
            key=lambda item: (
                order_map.get(item[0], len(order_map)),
                item[0].lower(),
            ),
        )
    ]


def catalog_orders(repo_root: Path) -> Dict[str, List[str]]:
    catalog = _load_catalog(
        str(repo_root / "chemdb" / "target_organization_annotations.yaml")
    )
    return {
        "adme": _section_order(_catalog_section(catalog, "adme"), ADME_DEFAULT),
        "safety": _section_order(_catalog_section(catalog, "safety"), SAFETY_DEFAULT),
    }
