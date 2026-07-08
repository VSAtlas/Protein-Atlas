from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple, cast

import yaml  # type: ignore[import-untyped]

from analysis.reporting.ligand_side_effect_cache import lookup_ligand_side_effect_entry
from analysis.reporting.side_effect_overlap_meta import collect_ligand_side_effect_fields
from analysis.reporting.value_utils import normalize_text

CHEMOTYPE_DEFAULT = "Other / Unassigned"
DRUG_EFFECT_DEFAULT = "Other / Unassigned"


def _normalize_text(value: Any) -> str:
    return normalize_text(value)


def _normalize_lower(value: Any) -> str:
    return _normalize_text(value).lower()


def _slugify(value: str) -> str:
    cleaned = re.sub(r"[^a-z0-9]+", "-", _normalize_lower(value))
    return cleaned.strip("-") or "unassigned"


@lru_cache(maxsize=1)
def _load_catalog(catalog_path: str) -> Dict[str, Any]:
    path = Path(catalog_path)
    if not path.exists():
        path = (
            Path(__file__).resolve().parents[2]
            / "chemdb"
            / "ligand_organization_annotations.yaml"
        )
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return cast(Dict[str, Any], payload)


@lru_cache(maxsize=4)
def _load_override_annotations(repo_root: Path) -> Dict[str, Any]:
    override_path = repo_root / "pathways" / "cache" / "ligand_annotation_overrides.yaml"
    if not override_path.exists():
        return {}
    payload = yaml.safe_load(override_path.read_text(encoding="utf-8")) or {}
    return cast(Dict[str, Any], payload)


def _catalog_section(catalog: Mapping[str, Any], section_name: str) -> Dict[str, Any]:
    section = catalog.get(section_name)
    return cast(Dict[str, Any], section if isinstance(section, dict) else {})


def _section_order(section: Mapping[str, Any], default_value: str) -> List[str]:
    order = [
        str(item).strip()
        for item in cast(Sequence[Any], section.get("order") or [])
        if str(item).strip()
    ]
    if default_value not in order:
        order.append(default_value)
    return order


def _rule_match_mode(
    rule: Mapping[str, Any],
    ligand_base: str,
    ligand_display: str,
    library: str,
    search_text: str,
) -> str:
    exact_names = {
        _normalize_lower(value)
        for value in cast(Sequence[Any], rule.get("exact_names") or [])
        if _normalize_text(value)
    }
    exact_bases = {
        _normalize_lower(value)
        for value in cast(Sequence[Any], rule.get("ligand_bases") or [])
        if _normalize_text(value)
    }
    library_keywords = {
        _normalize_lower(value)
        for value in cast(Sequence[Any], rule.get("libraries") or [])
        if _normalize_text(value)
    }
    if exact_names and (
        _normalize_lower(ligand_display) in exact_names
        or _normalize_lower(ligand_base) in exact_names
    ):
        return "exact"
    if exact_bases and _normalize_lower(ligand_base) in exact_bases:
        return "exact"
    if library_keywords and _normalize_lower(library) in library_keywords:
        return "exact"
    keywords = [
        _normalize_lower(value)
        for value in cast(Sequence[Any], rule.get("name_keywords") or [])
        if _normalize_text(value)
    ]
    if any(keyword in search_text for keyword in keywords):
        return "keyword"
    return ""


def _match_categories_with_quality(
    section: Mapping[str, Any],
    order: Sequence[str],
    ligand_base: str,
    ligand_display: str,
    library: str,
    search_text: str,
) -> List[Tuple[str, str]]:
    rules = cast(Sequence[Any], section.get("rules") or [])
    matched: List[Tuple[str, str]] = []
    for category in order:
        best_quality = ""
        for raw_rule in rules:
            rule = cast(Dict[str, Any], raw_rule if isinstance(raw_rule, dict) else {})
            if _normalize_text(rule.get("category")) != category:
                continue
            quality = _rule_match_mode(
                rule,
                ligand_base=ligand_base,
                ligand_display=ligand_display,
                library=library,
                search_text=search_text,
            )
            if quality == "exact":
                best_quality = "exact"
                break
            if quality == "keyword":
                best_quality = "keyword"
        if best_quality:
            matched.append((category, best_quality))
    return matched


def _match_motif_tags(
    section: Mapping[str, Any],
    ligand_base: str,
    ligand_display: str,
    library: str,
    search_text: str,
) -> List[str]:
    order = [
        _normalize_text(item)
        for item in cast(Sequence[Any], section.get("order") or [])
        if _normalize_text(item)
    ]
    rules = cast(Sequence[Any], section.get("rules") or [])
    matched: List[str] = []
    for tag in order:
        for raw_rule in rules:
            rule = cast(Dict[str, Any], raw_rule if isinstance(raw_rule, dict) else {})
            if _normalize_text(rule.get("tag")) != tag:
                continue
            if _rule_match_mode(
                rule,
                ligand_base=ligand_base,
                ligand_display=ligand_display,
                library=library,
                search_text=search_text,
            ):
                matched.append(tag)
                break
    return matched


def _normalize_memberships(
    values: Sequence[Any], order: Sequence[str], default_value: str
) -> List[str]:
    wanted = {_normalize_text(value) for value in values if _normalize_text(value)}
    ordered = [label for label in order if label in wanted]
    return ordered or [default_value]


def _merge_override(
    base: Dict[str, Any], override: Mapping[str, Any], *, chemotype_order: Sequence[str], effect_order: Sequence[str]
) -> Dict[str, Any]:
    if _normalize_text(override.get("chemotype_primary")):
        base["chemotype_primary"] = _normalize_text(override.get("chemotype_primary"))
    if _normalize_text(override.get("chemotype_series")):
        base["chemotype_series"] = _normalize_text(override.get("chemotype_series"))
    memberships = cast(Sequence[Any], override.get("drug_effect_memberships") or [])
    if memberships:
        base["drug_effect_memberships"] = _normalize_memberships(
            memberships, effect_order, DRUG_EFFECT_DEFAULT
        )
    if _normalize_text(override.get("drug_effect_primary")):
        base["drug_effect_primary"] = _normalize_text(override.get("drug_effect_primary"))
        if base["drug_effect_primary"] not in base["drug_effect_memberships"]:
            base["drug_effect_memberships"] = [
                base["drug_effect_primary"]
            ] + [
                item
                for item in cast(List[str], base["drug_effect_memberships"])
                if item != base["drug_effect_primary"]
            ]
    motif_tags = [
        _normalize_text(value)
        for value in cast(Sequence[Any], override.get("motif_tags") or [])
        if _normalize_text(value)
    ]
    if motif_tags:
        base["motif_tags"] = motif_tags
    if _normalize_text(override.get("confidence")):
        base["ligand_annotation_confidence"] = _normalize_text(
            override.get("confidence")
        )
    sources = [
        _normalize_text(value)
        for value in cast(Sequence[Any], override.get("sources") or [])
        if _normalize_text(value)
    ]
    if sources:
        base["ligand_annotation_sources"] = sources
    base.update(collect_ligand_side_effect_fields(override))
    base["chemotype_primary"] = (
        base["chemotype_primary"]
        if base["chemotype_primary"] in chemotype_order
        else CHEMOTYPE_DEFAULT
    )
    base["drug_effect_primary"] = (
        base["drug_effect_primary"]
        if base["drug_effect_primary"] in effect_order
        else DRUG_EFFECT_DEFAULT
    )
    return base


def resolve_ligand_annotations(
    repo_root: Path,
    *,
    ligand_name: str,
    ligand_base: str = "",
    library: str = "",
) -> Dict[str, Any]:
    catalog = _load_catalog(
        str(repo_root / "chemdb" / "ligand_organization_annotations.yaml")
    )
    cache_overrides = _load_override_annotations(repo_root)
    catalog_overrides = _catalog_section(catalog, "overrides")
    chemotype_section = _catalog_section(catalog, "chemotype")
    effect_section = _catalog_section(catalog, "drug_effect")
    motif_section = _catalog_section(catalog, "motif_tags")
    chemotype_order = _section_order(chemotype_section, CHEMOTYPE_DEFAULT)
    effect_order = _section_order(effect_section, DRUG_EFFECT_DEFAULT)
    search_text = " ".join(
        part
        for part in (
            _normalize_lower(ligand_name),
            _normalize_lower(ligand_base),
        )
        if part
    )

    chemotype_matches = _match_categories_with_quality(
        chemotype_section,
        chemotype_order,
        ligand_base,
        ligand_name,
        library,
        search_text,
    )
    effect_matches = _match_categories_with_quality(
        effect_section,
        effect_order,
        ligand_base,
        ligand_name,
        library,
        search_text,
    )
    motif_tags = _match_motif_tags(
        motif_section,
        ligand_base,
        ligand_name,
        library,
        search_text,
    )

    chemotype_primary = (
        chemotype_matches[0][0] if chemotype_matches else CHEMOTYPE_DEFAULT
    )
    drug_effect_memberships = (
        [item[0] for item in effect_matches] if effect_matches else [DRUG_EFFECT_DEFAULT]
    )
    drug_effect_primary = drug_effect_memberships[0]
    quality_hits = [quality for _category, quality in chemotype_matches + effect_matches]
    if "exact" in quality_hits:
        confidence = "medium"
    elif "keyword" in quality_hits or motif_tags:
        confidence = "low"
    else:
        confidence = "unassigned"
    sources = ["Ligand annotation catalog"] if quality_hits or motif_tags else []

    result: Dict[str, Any] = {
        "chemotype_primary": chemotype_primary,
        "chemotype_series": chemotype_primary,
        "motif_tags": motif_tags,
        "drug_effect_primary": drug_effect_primary,
        "drug_effect_memberships": drug_effect_memberships,
        "ligand_annotation_confidence": confidence,
        "ligand_annotation_sources": sources,
        "scaffold_key": _slugify(chemotype_primary),
        "pairability_key": _slugify(chemotype_primary),
        "normalized_smiles": "",
        "inchikey": "",
    }

    merged_overrides: Dict[str, Any] = {}
    for section in (catalog_overrides, cache_overrides):
        by_base = cast(Dict[str, Any], section.get("by_ligand_base") or {})
        by_name = cast(Dict[str, Any], section.get("by_ligand_name") or {})
        exact = by_base.get(_normalize_text(ligand_base)) or by_base.get(
            _normalize_lower(ligand_base)
        )
        if not isinstance(exact, dict):
            exact = by_name.get(_normalize_text(ligand_name)) or by_name.get(
                _normalize_lower(ligand_name)
            )
        if isinstance(exact, dict):
            merged_overrides.update(cast(Dict[str, Any], exact))
    if merged_overrides:
        result = _merge_override(
            result,
            merged_overrides,
            chemotype_order=chemotype_order,
            effect_order=effect_order,
        )

    result["motif_signature"] = "; ".join(cast(List[str], result["motif_tags"]))
    return result


def build_ligand_grouping_meta(
    repo_root: Path,
    rows: Iterable[Mapping[str, Any]],
    ligand_labels: Sequence[str],
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any]]:
    row_lookup: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        ligand_display = _normalize_text(row.get("ligand_display"))
        ligand_name = (
            _normalize_text(row.get("ligand_name"))
            or ligand_display
            or _normalize_text(row.get("ligand_base"))
        )
        if not ligand_name:
            continue
        entry = {
            "ligand_base": _normalize_text(row.get("ligand_base")),
            "ligand_display": ligand_display or ligand_name,
            "library": _normalize_text(row.get("library")),
            "side_effect_meta": collect_ligand_side_effect_fields(row),
        }
        for key in {ligand_name, ligand_display, _normalize_text(row.get("ligand_base"))}:
            if key and key not in row_lookup:
                row_lookup[key] = entry

    meta: Dict[str, Dict[str, Any]] = {}
    chemotype_counts: Dict[str, int] = {}
    effect_counts: Dict[str, int] = {}
    for ligand_label in ligand_labels:
        source = row_lookup.get(ligand_label, {})
        resolved = resolve_ligand_annotations(
            repo_root,
            ligand_name=_normalize_text(source.get("ligand_display")) or ligand_label,
            ligand_base=_normalize_text(source.get("ligand_base")),
            library=_normalize_text(source.get("library")),
        )
        resolved["ligand"] = ligand_label
        resolved["library"] = _normalize_text(source.get("library"))
        cached_side_effects = lookup_ligand_side_effect_entry(
            repo_root,
            ligand_label=ligand_label,
            ligand_base=_normalize_text(source.get("ligand_base")),
            ligand_display=_normalize_text(source.get("ligand_display")),
        )
        if cached_side_effects:
            resolved.update(collect_ligand_side_effect_fields(cached_side_effects))
        resolved.update(cast(Dict[str, Any], source.get("side_effect_meta") or {}))
        meta[ligand_label] = resolved
        chemotype_counts[resolved["chemotype_primary"]] = (
            chemotype_counts.get(resolved["chemotype_primary"], 0) + 1
        )
        effect_counts[resolved["drug_effect_primary"]] = (
            effect_counts.get(resolved["drug_effect_primary"], 0) + 1
        )

    catalog = _load_catalog(
        str(repo_root / "chemdb" / "ligand_organization_annotations.yaml")
    )
    chemotype_order = _section_order(_catalog_section(catalog, "chemotype"), CHEMOTYPE_DEFAULT)
    effect_order = _section_order(_catalog_section(catalog, "drug_effect"), DRUG_EFFECT_DEFAULT)
    summary = {
        "default_mode": "none",
        "chemotype_blocks": [
            {"key": label, "label": label, "count": chemotype_counts.get(label, 0)}
            for label in chemotype_order
            if chemotype_counts.get(label, 0) > 0 or label == CHEMOTYPE_DEFAULT
        ],
        "drug_effect_blocks": [
            {"key": label, "label": label, "count": effect_counts.get(label, 0)}
            for label in effect_order
            if effect_counts.get(label, 0) > 0 or label == DRUG_EFFECT_DEFAULT
        ],
    }
    return meta, summary
