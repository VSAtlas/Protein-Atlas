from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Optional, Protocol, Sequence, Set, Tuple, cast

import requests  # type: ignore[import-untyped]
import yaml  # type: ignore[import-untyped]

from analysis.reporting.value_utils import normalize_text
from config.output_paths import run_output_dir
from protein_prep.pdb_records import local_uniprots_from_pdb

SAFETY_DEFAULT = "Unassigned"
HEURISTIC_CONFIDENCE = "low"
MANUAL_OVERRIDE_CONFIDENCE = "high"
TARGET_SAFETY_CACHE_VERSION = 1
OPEN_TARGETS_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

DIRECT_SAFETY_SOURCE = "Open Targets safety"
DRUG_AE_SOURCE = "Open Targets drug adverse events"
HEURISTIC_SOURCE = "Local heuristic catalog"
MANUAL_OVERRIDE_SOURCE = "Manual override"

DIRECT_TARGET_WEIGHT = 5.0
DRUG_AE_WEIGHT = 2.0
MAX_DRUG_BUCKET_WEIGHT_PER_DRUG = 3.5
DIRECT_EXAMPLE_LIMIT = 6
DRUG_EVENT_EXAMPLE_LIMIT = 8
SAFETY_EVIDENCE_TERM_LIMIT = 64

_LOG = logging.getLogger(__name__)
_UNIPROT_ACCESSION_RE = re.compile(
    r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b"
)


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _normalize_text(value: Any) -> str:
    return normalize_text(value)


def _normalize_upper(value: Any) -> str:
    return _normalize_text(value).upper()


def _normalize_lower(value: Any) -> str:
    return _normalize_text(value).lower()


def _normalize_number(value: Any) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except Exception:
        return None


@lru_cache(maxsize=4)
def _load_bucket_map(bucket_map_path: str) -> Dict[str, Any]:
    path = Path(bucket_map_path)
    if not path.exists():
        path = Path(__file__).resolve().parents[2] / "chemdb" / "target_safety_bucket_map.yaml"
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return cast(Dict[str, Any], payload)


def safety_bucket_order(repo_root: Path) -> List[str]:
    payload = _load_bucket_map(str(repo_root / "chemdb" / "target_safety_bucket_map.yaml"))
    order = [
        _normalize_text(item)
        for item in cast(Sequence[Any], payload.get("order") or [])
        if _normalize_text(item)
    ]
    if SAFETY_DEFAULT not in order:
        order.append(SAFETY_DEFAULT)
    return order


def _bucket_rules(repo_root: Path) -> Dict[str, Dict[str, List[str]]]:
    payload = _load_bucket_map(str(repo_root / "chemdb" / "target_safety_bucket_map.yaml"))
    raw_rules = cast(Dict[str, Any], payload.get("buckets") or {})
    rules: Dict[str, Dict[str, List[str]]] = {}
    for bucket_name, raw in raw_rules.items():
        bucket = _normalize_text(bucket_name)
        if not bucket:
            continue
        section = cast(Dict[str, Any], raw if isinstance(raw, dict) else {})
        rules[bucket] = {
            "strong_event_keywords": [
                _normalize_lower(item)
                for item in cast(Sequence[Any], section.get("strong_event_keywords") or [])
                if _normalize_text(item)
            ],
            "weak_event_keywords": [
                _normalize_lower(item)
                for item in cast(Sequence[Any], section.get("weak_event_keywords") or [])
                if _normalize_text(item)
            ],
            "exclude_event_keywords": [
                _normalize_lower(item)
                for item in cast(Sequence[Any], section.get("exclude_event_keywords") or [])
                if _normalize_text(item)
            ],
            "strong_tissue_keywords": [
                _normalize_lower(item)
                for item in cast(Sequence[Any], section.get("strong_tissue_keywords") or [])
                if _normalize_text(item)
            ],
            "weak_tissue_keywords": [
                _normalize_lower(item)
                for item in cast(Sequence[Any], section.get("weak_tissue_keywords") or [])
                if _normalize_text(item)
            ],
            "study_keywords": [
                _normalize_lower(item)
                for item in cast(Sequence[Any], section.get("study_keywords") or [])
                if _normalize_text(item)
            ],
        }
    return rules


def _keyword_hit_count(text: str, keywords: Sequence[str]) -> int:
    return sum(1 for keyword in keywords if keyword and keyword in text)


def classify_safety_bucket_scores(
    repo_root: Path,
    *,
    event_text: str = "",
    tissue_text: str = "",
    study_text: str = "",
    evidence_source: str = "direct",
) -> Dict[str, float]:
    rules = _bucket_rules(repo_root)
    order = safety_bucket_order(repo_root)
    event_lower = _normalize_lower(event_text)
    tissue_lower = _normalize_lower(tissue_text)
    study_lower = _normalize_lower(study_text)
    source = _normalize_lower(evidence_source) or "direct"
    scores: Dict[str, float] = {}
    for bucket in order:
        if bucket == SAFETY_DEFAULT:
            continue
        rule = rules.get(bucket, {})
        strong_event_hits = _keyword_hit_count(event_lower, rule.get("strong_event_keywords", []))
        weak_event_hits = _keyword_hit_count(event_lower, rule.get("weak_event_keywords", []))
        exclude_event_hits = _keyword_hit_count(
            event_lower, rule.get("exclude_event_keywords", [])
        )
        strong_tissue_hits = _keyword_hit_count(
            tissue_lower, rule.get("strong_tissue_keywords", [])
        )
        weak_tissue_hits = _keyword_hit_count(
            tissue_lower, rule.get("weak_tissue_keywords", [])
        )
        study_hits = _keyword_hit_count(study_lower, rule.get("study_keywords", []))
        score = 0.0
        if source == "drug_ae":
            if strong_event_hits:
                score += min(strong_event_hits, 2) * 2.5
            if weak_event_hits:
                score += min(weak_event_hits, 3) * 0.6
            if exclude_event_hits and strong_event_hits == 0:
                score = max(0.0, score - 1.0)
            if score < 2.5:
                continue
            if strong_event_hits == 0 and weak_event_hits < 2:
                continue
        else:
            if strong_event_hits:
                score += min(strong_event_hits, 2) * 3.0
            if weak_event_hits:
                score += min(weak_event_hits, 3) * 0.8
            if strong_tissue_hits:
                score += min(strong_tissue_hits, 2) * 2.0
            if weak_tissue_hits:
                score += min(weak_tissue_hits, 2) * 0.75
            if study_hits:
                score += min(study_hits, 2) * 0.75
            if exclude_event_hits and strong_event_hits == 0 and strong_tissue_hits == 0:
                score = max(0.0, score - 0.75)
            if score < 2.0:
                continue
        if score > 0.0:
            scores[bucket] = round(score, 6)
    return scores


def classify_safety_buckets(
    repo_root: Path,
    *,
    event_text: str = "",
    tissue_text: str = "",
    study_text: str = "",
    evidence_source: str = "direct",
) -> List[str]:
    order = safety_bucket_order(repo_root)
    scores = classify_safety_bucket_scores(
        repo_root,
        event_text=event_text,
        tissue_text=tissue_text,
        study_text=study_text,
        evidence_source=evidence_source,
    )
    return [bucket for bucket in order if float(scores.get(bucket, 0.0)) > 0.0]


def _confidence_rank(level: str) -> int:
    normalized = _normalize_lower(level)
    if normalized == "high":
        return 3
    if normalized == "medium":
        return 2
    if normalized == "low":
        return 1
    return 0


def _confidence_from_evidence(
    has_direct: bool, source_names: Set[str], manual_override: bool = False
) -> str:
    if manual_override:
        return MANUAL_OVERRIDE_CONFIDENCE
    if has_direct:
        return "high"
    if len(source_names) >= 2:
        return "medium"
    if len(source_names) == 1:
        return "low"
    return "unassigned"


def _stringify_sources(values: Iterable[str]) -> List[str]:
    return sorted({_normalize_text(value) for value in values if _normalize_text(value)})


def _bucket_scores_from_entries(
    order: Sequence[str], entries: Sequence[Dict[str, Any]]
) -> Dict[str, float]:
    scores: Dict[str, float] = {}
    for entry in entries:
        bucket = _normalize_text(entry.get("bucket"))
        if not bucket:
            continue
        weight = _normalize_number(entry.get("weight")) or 0.0
        scores[bucket] = scores.get(bucket, 0.0) + weight
    for bucket in order:
        scores.setdefault(bucket, 0.0)
    return scores


def _primary_bucket_from_scores(
    scores: Dict[str, float], order: Sequence[str], default_value: str
) -> str:
    populated = [(bucket, score) for bucket, score in scores.items() if score > 0]
    if not populated:
        return default_value
    order_map = {name: idx for idx, name in enumerate(order)}
    bucket, _score = min(
        populated,
        key=lambda item: (
            -item[1],
            order_map.get(item[0], len(order_map)),
            item[0].lower(),
        ),
    )
    return bucket


def _bucket_list_from_scores(
    scores: Dict[str, float], order: Sequence[str], default_value: str
) -> List[str]:
    buckets = [bucket for bucket in order if scores.get(bucket, 0.0) > 0]
    return buckets or [default_value]


def _secondary_bucket_list(
    scores: Dict[str, float], order: Sequence[str], primary: str, default_value: str
) -> List[str]:
    if primary == default_value:
        return []
    return [bucket for bucket in _bucket_list_from_scores(scores, order, default_value) if bucket != primary]


def _format_bucket_scores(scores: Dict[str, float], order: Sequence[str]) -> Dict[str, str]:
    return {
        bucket: f"{float(scores.get(bucket, 0.0)):.3f}"
        for bucket in order
        if float(scores.get(bucket, 0.0)) > 0.0
    }


def _merge_scores(
    order: Sequence[str], *score_maps: Dict[str, float]
) -> Dict[str, float]:
    merged: Dict[str, float] = {bucket: 0.0 for bucket in order}
    for score_map in score_maps:
        for bucket in order:
            merged[bucket] = float(merged.get(bucket, 0.0)) + float(score_map.get(bucket, 0.0))
    return merged


def _unique_example_events(
    entries: Sequence[Dict[str, Any]], source_name: str, limit: int
) -> List[str]:
    seen: Set[str] = set()
    examples: List[str] = []
    for entry in entries:
        if _normalize_lower(entry.get("source")) != _normalize_lower(source_name):
            continue
        event_name = _normalize_text(entry.get("event"))
        key = event_name.lower()
        if not event_name or key in seen:
            continue
        seen.add(key)
        examples.append(event_name)
        if len(examples) >= limit:
            break
    return examples


def _unique_evidence_terms(
    entries: Sequence[Dict[str, Any]], limit: int = SAFETY_EVIDENCE_TERM_LIMIT
) -> List[str]:
    seen: Set[str] = set()
    terms: List[str] = []
    for entry in entries:
        event_name = _normalize_text(entry.get("event"))
        key = event_name.lower()
        if not event_name or key in seen:
            continue
        seen.add(key)
        terms.append(event_name)
        if len(terms) >= limit:
            break
    return terms


def _build_evidence_summary(entries: Sequence[Dict[str, Any]]) -> str:
    counts: Dict[str, int] = {}
    for entry in entries:
        source = _normalize_text(entry.get("source")) or "Unknown source"
        counts[source] = counts.get(source, 0) + 1
    return "; ".join(
        f"{source}: {count}" for source, count in sorted(counts.items(), key=lambda item: item[0].lower())
    )


def _compose_safety_profile(
    order: Sequence[str],
    *,
    approved_symbol: str,
    approved_name: str,
    direct_scores: Dict[str, float],
    drug_scores: Dict[str, float],
    source_names: Set[str],
    evidence_entries: Sequence[Dict[str, Any]],
    confidence: str,
) -> Dict[str, Any]:
    direct_buckets = _bucket_list_from_scores(direct_scores, order, SAFETY_DEFAULT)
    drug_buckets = _bucket_list_from_scores(drug_scores, order, SAFETY_DEFAULT)
    total_scores = _merge_scores(order, direct_scores, drug_scores)
    primary_scores = direct_scores if any(score > 0.0 for score in direct_scores.values()) else total_scores
    primary = _primary_bucket_from_scores(primary_scores, order, SAFETY_DEFAULT)
    secondary = _secondary_bucket_list(total_scores, order, primary, SAFETY_DEFAULT)
    safety_buckets = [primary] + [bucket for bucket in secondary if bucket != primary]
    if primary == SAFETY_DEFAULT:
        safety_buckets = [SAFETY_DEFAULT]
    direct_examples = _unique_example_events(
        evidence_entries, DIRECT_SAFETY_SOURCE, DIRECT_EXAMPLE_LIMIT
    )
    drug_examples = _unique_example_events(
        evidence_entries, DRUG_AE_SOURCE, DRUG_EVENT_EXAMPLE_LIMIT
    )
    return {
        "approved_symbol": approved_symbol,
        "approved_name": approved_name,
        "safety_buckets": safety_buckets,
        "primary_display_safety": primary,
        "secondary_safety_buckets": secondary,
        "direct_safety_buckets": direct_buckets if direct_buckets != [SAFETY_DEFAULT] else [],
        "drug_ae_buckets": drug_buckets if drug_buckets != [SAFETY_DEFAULT] else [],
        "safety_confidence": confidence,
        "safety_sources": _stringify_sources(source_names),
        "safety_bucket_scores": _format_bucket_scores(total_scores, order),
        "direct_bucket_scores": _format_bucket_scores(direct_scores, order),
        "drug_ae_bucket_scores": _format_bucket_scores(drug_scores, order),
        "safety_evidence_summary": _build_evidence_summary(evidence_entries),
        "safety_evidence_terms": _unique_evidence_terms(evidence_entries),
        "direct_liability_examples": direct_examples,
        "raw_adverse_event_examples": drug_examples,
    }


def _normalize_cached_profile(
    repo_root: Path, profile: Dict[str, Any], *, fallback_name: str = ""
) -> Dict[str, Any]:
    order = safety_bucket_order(repo_root)
    raw_scores = cast(Dict[str, Any], profile.get("safety_bucket_scores") or {})
    raw_direct_scores = cast(Dict[str, Any], profile.get("direct_bucket_scores") or {})
    raw_drug_scores = cast(Dict[str, Any], profile.get("drug_ae_bucket_scores") or {})
    numeric_scores = {
        bucket: float(_normalize_number(raw_scores.get(bucket)) or 0.0) for bucket in order
    }
    direct_scores = {
        bucket: float(_normalize_number(raw_direct_scores.get(bucket)) or 0.0) for bucket in order
    }
    drug_scores = {
        bucket: float(_normalize_number(raw_drug_scores.get(bucket)) or 0.0) for bucket in order
    }
    safety_buckets = [
        _normalize_text(value)
        for value in cast(Sequence[Any], profile.get("safety_buckets") or [])
        if _normalize_text(value)
    ]
    if not safety_buckets:
        safety_buckets = _bucket_list_from_scores(numeric_scores, order, SAFETY_DEFAULT)
    primary = _normalize_text(profile.get("primary_display_safety")) or _primary_bucket_from_scores(
        numeric_scores, order, SAFETY_DEFAULT
    )
    evidence_entries = [
        cast(Dict[str, Any], item)
        for item in cast(Sequence[Any], profile.get("evidence") or [])
        if isinstance(item, dict)
    ]
    source_names = _stringify_sources(
        list(cast(Sequence[Any], profile.get("safety_sources") or []))
        + [str(entry.get("source") or "") for entry in evidence_entries]
    )
    confidence = _normalize_lower(profile.get("safety_confidence")) or _confidence_from_evidence(
        any(_normalize_lower(entry.get("source")) == _normalize_lower(DIRECT_SAFETY_SOURCE) for entry in evidence_entries),
        set(source_names),
    )
    if confidence not in {"high", "medium", "low", "unassigned"}:
        confidence = "unassigned"
    summary = _normalize_text(profile.get("safety_evidence_summary")) or _build_evidence_summary(
        evidence_entries
    )
    approved_symbol = _normalize_text(profile.get("approved_symbol"))
    approved_name = _normalize_text(profile.get("approved_name")) or fallback_name
    secondary = [
        _normalize_text(value)
        for value in cast(Sequence[Any], profile.get("secondary_safety_buckets") or [])
        if _normalize_text(value)
    ]
    if not secondary:
        secondary = _secondary_bucket_list(numeric_scores, order, primary, SAFETY_DEFAULT)
    direct_buckets = [
        _normalize_text(value)
        for value in cast(Sequence[Any], profile.get("direct_safety_buckets") or [])
        if _normalize_text(value)
    ]
    if not direct_buckets:
        direct_buckets = [
            bucket
            for bucket in _bucket_list_from_scores(direct_scores, order, SAFETY_DEFAULT)
            if bucket != SAFETY_DEFAULT
        ]
    drug_buckets = [
        _normalize_text(value)
        for value in cast(Sequence[Any], profile.get("drug_ae_buckets") or [])
        if _normalize_text(value)
    ]
    if not drug_buckets:
        drug_buckets = [
            bucket
            for bucket in _bucket_list_from_scores(drug_scores, order, SAFETY_DEFAULT)
            if bucket != SAFETY_DEFAULT
        ]
    return {
        "approved_symbol": approved_symbol,
        "approved_name": approved_name,
        "safety_buckets": safety_buckets or [SAFETY_DEFAULT],
        "primary_display_safety": primary or SAFETY_DEFAULT,
        "secondary_safety_buckets": secondary,
        "direct_safety_buckets": direct_buckets,
        "drug_ae_buckets": drug_buckets,
        "safety_confidence": confidence,
        "safety_sources": source_names,
        "safety_bucket_scores": _format_bucket_scores(numeric_scores, order),
        "direct_bucket_scores": _format_bucket_scores(direct_scores, order),
        "drug_ae_bucket_scores": _format_bucket_scores(drug_scores, order),
        "safety_evidence_summary": summary,
        "safety_evidence_terms": _unique_evidence_terms(
            cast(Sequence[Dict[str, Any]], profile.get("evidence") or [])
        ),
        "direct_liability_examples": [
            _normalize_text(value)
            for value in cast(Sequence[Any], profile.get("direct_liability_examples") or [])
            if _normalize_text(value)
        ],
        "raw_adverse_event_examples": [
            _normalize_text(value)
            for value in cast(Sequence[Any], profile.get("raw_adverse_event_examples") or [])
            if _normalize_text(value)
        ],
    }


@lru_cache(maxsize=4)
def _load_safety_cache(cache_path: str) -> Dict[str, Any]:
    path = Path(cache_path)
    if not path.exists():
        return {"version": TARGET_SAFETY_CACHE_VERSION, "entries": [], "by_uniprot": {}, "by_gene": {}}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return {"version": TARGET_SAFETY_CACHE_VERSION, "entries": [], "by_uniprot": {}, "by_gene": {}}
    return cast(Dict[str, Any], payload)


@lru_cache(maxsize=4)
def _load_local_target_adr_profiles(repo_root_text: str) -> Dict[str, Dict[str, Dict[str, Any]]]:
    repo_root = Path(repo_root_text)
    order = safety_bucket_order(repo_root)
    candidate_paths = [
        repo_root / "data" / "external" / "adrecs_target" / "target_adr_evidence.tsv",
        repo_root / "data" / "external" / "paper_sources" / "normalized_v2" / "combined_target_adr_evidence.tsv",
        repo_root / "data" / "external" / "opentargets" / "safety.tsv",
    ]
    raw_profiles: Dict[str, Dict[str, Any]] = {}
    for path in candidate_paths:
        if not path.exists():
            continue
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    label_state = _normalize_text(row.get("label_state"))
                    if label_state and label_state not in {"1", "1.0", "strict_positive"}:
                        continue
                    namespace = _normalize_upper(
                        row.get("target_id") or row.get("uniprot") or row.get("source_target_id")
                    )
                    gene_symbol = _normalize_upper(
                        row.get("gene_symbol") or row.get("_source_approved_symbol")
                    )
                    if not namespace and not gene_symbol:
                        continue
                    adr_term = _normalize_text(row.get("adr_term") or row.get("adr"))
                    if not adr_term:
                        continue
                    namespace = namespace or gene_symbol
                    source = _normalize_text(row.get("source_name") or row.get("source")) or path.stem
                    profile = raw_profiles.setdefault(
                        namespace,
                        {
                            "approved_symbol": gene_symbol,
                            "approved_name": _normalize_text(
                                row.get("target_name") or row.get("_source_approved_name")
                            ),
                            "source_names": set(),
                            "direct_scores": {bucket: 0.0 for bucket in order},
                            "drug_scores": {bucket: 0.0 for bucket in order},
                            "evidence": [],
                        },
                    )
                    if gene_symbol and not profile.get("approved_symbol"):
                        profile["approved_symbol"] = gene_symbol
                    cast(Set[str], profile["source_names"]).add(source)
                    scores = classify_safety_bucket_scores(
                        repo_root,
                        event_text=adr_term,
                        study_text=source,
                        evidence_source="direct",
                    )
                    if scores:
                        for bucket, score in scores.items():
                            direct_scores = cast(Dict[str, float], profile["direct_scores"])
                            direct_scores[bucket] = direct_scores.get(bucket, 0.0) + float(score)
                            cast(List[Dict[str, Any]], profile["evidence"]).append(
                                {
                                    "source": source,
                                    "event": adr_term,
                                    "bucket": bucket,
                                    "weight": float(score),
                                }
                            )
                    else:
                        cast(List[Dict[str, Any]], profile["evidence"]).append(
                            {
                                "source": source,
                                "event": adr_term,
                                "bucket": SAFETY_DEFAULT,
                                "weight": 0.0,
                            }
                        )
        except (OSError, csv.Error) as exc:
            _LOG.warning("local target ADR evidence load failed path=%s error=%s", path, exc)

    by_uniprot: Dict[str, Dict[str, Any]] = {}
    by_gene: Dict[str, Dict[str, Any]] = {}
    for namespace, raw in raw_profiles.items():
        profile = _compose_safety_profile(
            order,
            approved_symbol=_normalize_text(raw.get("approved_symbol")),
            approved_name=_normalize_text(raw.get("approved_name")) or namespace,
            direct_scores=cast(Dict[str, float], raw.get("direct_scores") or {}),
            drug_scores=cast(Dict[str, float], raw.get("drug_scores") or {}),
            source_names=cast(Set[str], raw.get("source_names") or set()),
            evidence_entries=cast(Sequence[Dict[str, Any]], raw.get("evidence") or []),
            confidence="medium",
        )
        profile["evidence"] = cast(Sequence[Dict[str, Any]], raw.get("evidence") or [])
        if _UNIPROT_ACCESSION_RE.fullmatch(namespace):
            by_uniprot[namespace] = profile
        else:
            by_gene[namespace] = profile
    return {"by_uniprot": by_uniprot, "by_gene": by_gene}


def clear_target_safety_cache() -> None:
    _load_safety_cache.cache_clear()
    _load_bucket_map.cache_clear()
    _load_local_target_adr_profiles.cache_clear()


def load_target_safety_profile(
    repo_root: Path,
    *,
    target_name: str,
    uniprots: Iterable[str],
    gene_symbols: Iterable[str] | None = None,
) -> Optional[Dict[str, Any]]:
    cache = _load_safety_cache(str(repo_root / "pathways" / "cache" / "target_safety_aggregated.json"))
    by_uniprot = cast(Dict[str, Any], cache.get("by_uniprot") or {})
    by_gene = cast(Dict[str, Any], cache.get("by_gene") or {})
    matches: List[Dict[str, Any]] = []
    for accession in sorted({_normalize_upper(value) for value in uniprots if _normalize_text(value)}):
        match = by_uniprot.get(accession)
        if isinstance(match, dict):
            matches.append(cast(Dict[str, Any], match))
    for symbol in sorted({_normalize_upper(value) for value in (gene_symbols or []) if _normalize_text(value)}):
        match = by_gene.get(symbol)
        if isinstance(match, dict):
            matches.append(cast(Dict[str, Any], match))
    if not matches:
        local_cache = _load_local_target_adr_profiles(str(repo_root.resolve()))
        local_by_uniprot = cast(Dict[str, Any], local_cache.get("by_uniprot") or {})
        local_by_gene = cast(Dict[str, Any], local_cache.get("by_gene") or {})
        for accession in sorted(
            {_normalize_upper(value) for value in uniprots if _normalize_text(value)}
        ):
            match = local_by_uniprot.get(accession)
            if isinstance(match, dict):
                matches.append(cast(Dict[str, Any], match))
        for symbol in sorted(
            {_normalize_upper(value) for value in (gene_symbols or []) if _normalize_text(value)}
        ):
            match = local_by_gene.get(symbol)
            if isinstance(match, dict):
                matches.append(cast(Dict[str, Any], match))
    if not matches:
        return None
    if len(matches) == 1:
        return _normalize_cached_profile(repo_root, matches[0], fallback_name=target_name)

    order = safety_bucket_order(repo_root)
    merged_entries: List[Dict[str, Any]] = []
    source_names: Set[str] = set()
    bucket_scores: Dict[str, float] = {bucket: 0.0 for bucket in order}
    direct_scores: Dict[str, float] = {bucket: 0.0 for bucket in order}
    drug_scores: Dict[str, float] = {bucket: 0.0 for bucket in order}
    approved_symbol = ""
    approved_name = target_name
    best_confidence = "unassigned"
    direct_examples: List[str] = []
    drug_examples: List[str] = []
    for match in matches:
        normalized = _normalize_cached_profile(repo_root, match, fallback_name=target_name)
        approved_symbol = approved_symbol or _normalize_text(normalized.get("approved_symbol"))
        approved_name = approved_name or _normalize_text(normalized.get("approved_name"))
        best_confidence = max(
            [best_confidence, _normalize_text(normalized.get("safety_confidence"))],
            key=_confidence_rank,
        )
        source_names.update(cast(List[str], normalized.get("safety_sources") or []))
        merged_entries.extend(
            cast(List[Dict[str, Any]], cast(Dict[str, Any], match).get("evidence") or [])
        )
        raw_scores = cast(Dict[str, Any], normalized.get("safety_bucket_scores") or {})
        raw_direct = cast(Dict[str, Any], normalized.get("direct_bucket_scores") or {})
        raw_drug = cast(Dict[str, Any], normalized.get("drug_ae_bucket_scores") or {})
        for bucket in order:
            bucket_scores[bucket] = bucket_scores.get(bucket, 0.0) + float(
                _normalize_number(raw_scores.get(bucket)) or 0.0
            )
            direct_scores[bucket] = direct_scores.get(bucket, 0.0) + float(
                _normalize_number(raw_direct.get(bucket)) or 0.0
            )
            drug_scores[bucket] = drug_scores.get(bucket, 0.0) + float(
                _normalize_number(raw_drug.get(bucket)) or 0.0
            )
        for example in cast(List[str], normalized.get("direct_liability_examples") or []):
            if example not in direct_examples:
                direct_examples.append(example)
        for example in cast(List[str], normalized.get("raw_adverse_event_examples") or []):
            if example not in drug_examples:
                drug_examples.append(example)

    composed = _compose_safety_profile(
        order,
        approved_symbol=approved_symbol,
        approved_name=approved_name or target_name,
        direct_scores=direct_scores,
        drug_scores=drug_scores,
        source_names=source_names,
        evidence_entries=merged_entries,
        confidence=best_confidence,
    )
    composed["direct_liability_examples"] = direct_examples[:DIRECT_EXAMPLE_LIMIT]
    composed["raw_adverse_event_examples"] = drug_examples[:DRUG_EVENT_EXAMPLE_LIMIT]
    return composed


@dataclass(frozen=True)
class TargetSpec:
    pdb_id: str
    target_name: str
    uniprots: Tuple[str, ...]
    gene_symbols: Tuple[str, ...] = ()


class OpenTargetsClient:
    def __init__(self, timeout: int = 20) -> None:
        self.timeout = timeout
        self.session = requests.Session()

    def _graphql(self, query: str, variables: Dict[str, Any]) -> Dict[str, Any]:
        response = self.session.post(
            OPEN_TARGETS_GRAPHQL_URL,
            json={"query": query, "variables": variables},
            timeout=self.timeout,
        )
        response.raise_for_status()
        payload = response.json()
        if isinstance(payload, dict) and payload.get("errors"):
            raise RuntimeError(str(payload.get("errors")))
        return cast(Dict[str, Any], payload.get("data") or {})

    def search_target(self, query_term: str) -> List[Dict[str, Any]]:
        query = (
            "query Search($q:String!){ "
            "search(queryString:$q, entityNames:[\"target\"], page:{index:0,size:5})"
            "{ hits { id name entity } } }"
        )
        data = self._graphql(query, {"q": query_term})
        search = cast(Dict[str, Any], data.get("search") or {})
        hits = cast(List[Any], search.get("hits") or [])
        return [cast(Dict[str, Any], hit) for hit in hits if isinstance(hit, dict)]

    def fetch_target(self, ensembl_id: str, known_drug_size: int = 40) -> Dict[str, Any]:
        query = (
            "query Target($id:String!, $size:Int!){ "
            "target(ensemblId:$id){ "
            "approvedSymbol approvedName "
            "proteinIds{ id source } "
            "safetyLiabilities{ event eventId datasource effects{ direction dosing } "
            "biosamples{ tissueLabel tissueId cellLabel cellFormat } studies{ name description type } } "
            "knownDrugs(size:$size){ count rows{ drugId prefName phase mechanismOfAction } } "
            "} }"
        )
        data = self._graphql(query, {"id": ensembl_id, "size": known_drug_size})
        return cast(Dict[str, Any], data.get("target") or {})

    def fetch_drug(self, chembl_id: str) -> Dict[str, Any]:
        query = (
            "query Drug($id:String!){ "
            "drug(chemblId:$id){ id name adverseEvents{ count rows{ name count logLR meddraCode } } } }"
        )
        data = self._graphql(query, {"id": chembl_id})
        return cast(Dict[str, Any], data.get("drug") or {})


class TargetSafetyClient(Protocol):
    def search_target(self, term: str) -> List[Dict[str, Any]]: ...

    def fetch_target(self, ensembl_id: str, known_drug_size: int = 40) -> Dict[str, Any]: ...

    def fetch_drug(self, chembl_id: str) -> Dict[str, Any]: ...


def _best_search_hit(
    client: TargetSafetyClient,
    *,
    target_name: str,
    uniprots: Sequence[str],
    gene_symbols: Sequence[str],
) -> Optional[Dict[str, Any]]:
    candidate_terms = list(uniprots) + list(gene_symbols)
    if _normalize_text(target_name):
        candidate_terms.append(target_name)
    seen: Set[str] = set()
    for term in candidate_terms:
        normalized_term = _normalize_text(term)
        if not normalized_term or normalized_term.lower() in seen:
            continue
        seen.add(normalized_term.lower())
        try:
            hits = client.search_target(normalized_term)
        except Exception as exc:
            _LOG.warning("[target-safety.search.failed] term=%s err=%s", normalized_term, exc)
            continue
        if hits:
            return hits[0]
    return None


def _string_terms(values: Iterable[Any]) -> str:
    return " | ".join(_normalize_text(value) for value in values if _normalize_text(value))


def _extract_direct_evidence(
    repo_root: Path, liabilities: Sequence[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    evidence: List[Dict[str, Any]] = []
    for liability in liabilities:
        event = _normalize_text(liability.get("event"))
        tissue_text = _string_terms(
            [cast(Dict[str, Any], item).get("tissueLabel") for item in cast(List[Any], liability.get("biosamples") or [])]
        )
        study_text = _string_terms(
            [cast(Dict[str, Any], item).get("type") for item in cast(List[Any], liability.get("studies") or [])]
        )
        bucket_scores = classify_safety_bucket_scores(
            repo_root,
            event_text=event,
            tissue_text=tissue_text,
            study_text=study_text,
            evidence_source="direct",
        )
        for bucket, match_score in bucket_scores.items():
            evidence.append(
                {
                    "bucket": bucket,
                    "source": DIRECT_SAFETY_SOURCE,
                    "weight": DIRECT_TARGET_WEIGHT + match_score,
                    "event": event,
                    "tissue": tissue_text,
                    "study": study_text,
                    "datasource": _normalize_text(liability.get("datasource")),
                }
            )
    return evidence


def _extract_drug_evidence(
    repo_root: Path,
    client: TargetSafetyClient,
    known_drugs: Sequence[Dict[str, Any]],
    *,
    max_approved_drugs: int,
) -> List[Dict[str, Any]]:
    approved: List[Dict[str, Any]] = []
    seen_ids: Set[str] = set()
    for drug in known_drugs:
        chembl_id = _normalize_text(drug.get("drugId"))
        if not chembl_id or chembl_id in seen_ids:
            continue
        phase = int(_normalize_number(drug.get("phase")) or 0)
        if phase < 4:
            continue
        seen_ids.add(chembl_id)
        approved.append(drug)
        if len(approved) >= max_approved_drugs:
            break

    evidence: List[Dict[str, Any]] = []
    for drug in approved:
        chembl_id = _normalize_text(drug.get("drugId"))
        try:
            payload = client.fetch_drug(chembl_id)
        except Exception as exc:
            _LOG.warning("[target-safety.drug.failed] drug=%s err=%s", chembl_id, exc)
            continue
        adverse = cast(Dict[str, Any], payload.get("adverseEvents") or {})
        rows = [
            cast(Dict[str, Any], row)
            for row in cast(List[Any], adverse.get("rows") or [])
            if isinstance(row, dict)
        ]
        rows = sorted(
            rows,
            key=lambda row: (
                -(float(_normalize_number(row.get("logLR")) or 0.0)),
                -(float(_normalize_number(row.get("count")) or 0.0)),
                _normalize_text(row.get("name")).lower(),
            ),
        )
        best_by_bucket: Dict[str, Dict[str, Any]] = {}
        for row in rows:
            event_name = _normalize_text(row.get("name"))
            if not event_name:
                continue
            bucket_scores = classify_safety_bucket_scores(
                repo_root,
                event_text=event_name,
                evidence_source="drug_ae",
            )
            if not bucket_scores:
                continue
            log_lr = float(_normalize_number(row.get("logLR")) or 0.0)
            count = float(_normalize_number(row.get("count")) or 0.0)
            base_weight = DRUG_AE_WEIGHT + max(0.0, min(log_lr, 3.0)) * 0.25
            if count >= 100:
                base_weight += 0.5
            for bucket, match_score in bucket_scores.items():
                weight = min(
                    MAX_DRUG_BUCKET_WEIGHT_PER_DRUG,
                    base_weight + min(match_score, 3.0) * 0.35,
                )
                candidate = {
                    "bucket": bucket,
                    "source": DRUG_AE_SOURCE,
                    "weight": weight,
                    "event": event_name,
                    "drug_id": chembl_id,
                    "drug_name": _normalize_text(payload.get("name"))
                    or _normalize_text(drug.get("prefName")),
                    "log_lr": f"{log_lr:.3f}",
                    "count": f"{count:.0f}",
                }
                current = best_by_bucket.get(bucket)
                if current is None or float(_normalize_number(current.get("weight")) or 0.0) < weight:
                    best_by_bucket[bucket] = candidate
        evidence.extend(best_by_bucket.values())
    return evidence


def build_target_safety_entry(
    repo_root: Path,
    target: TargetSpec,
    *,
    client: Optional[TargetSafetyClient] = None,
    include_drug_evidence: bool = True,
    max_approved_drugs: int = 12,
) -> Optional[Dict[str, Any]]:
    ot_client = client or OpenTargetsClient()
    hit = _best_search_hit(
        ot_client,
        target_name=target.target_name,
        uniprots=list(target.uniprots),
        gene_symbols=list(target.gene_symbols),
    )
    if not hit:
        return None
    ensembl_id = _normalize_text(hit.get("id"))
    if not ensembl_id:
        return None
    target_payload = ot_client.fetch_target(ensembl_id, known_drug_size=max(5, max_approved_drugs * 3))
    approved_symbol = _normalize_text(target_payload.get("approvedSymbol"))
    approved_name = _normalize_text(target_payload.get("approvedName")) or target.target_name
    protein_ids = [
        cast(Dict[str, Any], item)
        for item in cast(List[Any], target_payload.get("proteinIds") or [])
        if isinstance(item, dict)
    ]
    swissprot_ids = [
        _normalize_upper(item.get("id"))
        for item in protein_ids
        if _normalize_lower(item.get("source")) == "uniprot_swissprot" and _normalize_text(item.get("id"))
    ]
    gene_symbols = sorted(
        {_normalize_upper(value) for value in list(target.gene_symbols) + [approved_symbol] if _normalize_text(value)}
    )
    uniprots = sorted(
        {_normalize_upper(value) for value in list(target.uniprots) + swissprot_ids if _normalize_text(value)}
    )
    direct_evidence = _extract_direct_evidence(
        repo_root,
        [
            cast(Dict[str, Any], item)
            for item in cast(List[Any], target_payload.get("safetyLiabilities") or [])
            if isinstance(item, dict)
        ],
    )
    drug_evidence: List[Dict[str, Any]] = []
    if include_drug_evidence:
        drug_evidence = _extract_drug_evidence(
            repo_root,
            ot_client,
            [
                cast(Dict[str, Any], item)
                for item in cast(List[Any], cast(Dict[str, Any], target_payload.get("knownDrugs") or {}).get("rows") or [])
                if isinstance(item, dict)
            ],
            max_approved_drugs=max_approved_drugs,
        )
    order = safety_bucket_order(repo_root)
    evidence = direct_evidence + drug_evidence
    direct_scores = _bucket_scores_from_entries(order, direct_evidence)
    drug_scores = _bucket_scores_from_entries(order, drug_evidence)
    source_names = {str(entry.get("source") or "") for entry in evidence if _normalize_text(entry.get("source"))}
    confidence = _confidence_from_evidence(bool(direct_evidence), source_names)
    payload = _compose_safety_profile(
        order,
        approved_symbol=approved_symbol,
        approved_name=approved_name,
        direct_scores=direct_scores,
        drug_scores=drug_scores,
        source_names=source_names,
        evidence_entries=evidence,
        confidence=confidence,
    )
    payload.update(
        {
            "ensembl_id": ensembl_id,
            "uniprots": uniprots,
            "gene_symbols": gene_symbols,
            "evidence": evidence,
        }
    )
    return payload


def build_target_safety_cache(
    repo_root: Path,
    targets: Sequence[TargetSpec],
    *,
    include_drug_evidence: bool = True,
    max_approved_drugs: int = 12,
) -> Dict[str, Any]:
    client = OpenTargetsClient()
    entries: List[Dict[str, Any]] = []
    by_uniprot: Dict[str, Dict[str, Any]] = {}
    by_gene: Dict[str, Dict[str, Any]] = {}
    for target in targets:
        try:
            entry = build_target_safety_entry(
                repo_root,
                target,
                client=client,
                include_drug_evidence=include_drug_evidence,
                max_approved_drugs=max_approved_drugs,
            )
        except Exception as exc:
            _LOG.warning(
                "[target-safety.build.failed] pdb=%s target=%s err=%s",
                target.pdb_id,
                target.target_name,
                exc,
            )
            continue
        if not entry:
            continue
        entries.append(entry)
        for accession in cast(List[str], entry.get("uniprots") or []):
            by_uniprot[_normalize_upper(accession)] = entry
        for symbol in cast(List[str], entry.get("gene_symbols") or []):
            by_gene[_normalize_upper(symbol)] = entry

    return {
        "version": TARGET_SAFETY_CACHE_VERSION,
        "generated_at": _now_iso(),
        "source": "open_targets",
        "entries": entries,
        "by_uniprot": by_uniprot,
        "by_gene": by_gene,
    }


def write_target_safety_cache(repo_root: Path, payload: Dict[str, Any], out_path: Optional[Path] = None) -> Path:
    path = out_path or (repo_root / "pathways" / "cache" / "target_safety_aggregated.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    clear_target_safety_cache()
    return path


def extract_local_uniprots_from_pdb(repo_root: Path, pdb_id: str) -> List[str]:
    return local_uniprots_from_pdb(repo_root, pdb_id)


def load_run_targets(repo_root: Path, run_id: str) -> List[TargetSpec]:
    data_dir = run_output_dir(repo_root, "data", run_id)
    heatmap_path = data_dir / "heatmap_input.csv"
    if not heatmap_path.exists():
        raise FileNotFoundError(f"heatmap_input.csv not found for run {run_id}: {heatmap_path}")
    targets: Dict[str, TargetSpec] = {}
    with heatmap_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pdb_id = _normalize_upper(row.get("pdb_id"))
            if not pdb_id:
                continue
            target_name = _normalize_text(row.get("target_name"))
            current = targets.get(pdb_id)
            if current is None:
                uniprots = tuple(extract_local_uniprots_from_pdb(repo_root, pdb_id))
                targets[pdb_id] = TargetSpec(
                    pdb_id=pdb_id,
                    target_name=target_name,
                    uniprots=uniprots,
                    gene_symbols=(),
                )
                continue
            if target_name and not current.target_name:
                targets[pdb_id] = TargetSpec(
                    pdb_id=current.pdb_id,
                    target_name=target_name,
                    uniprots=current.uniprots,
                    gene_symbols=current.gene_symbols,
                )
    return sorted(targets.values(), key=lambda item: item.pdb_id.lower())
