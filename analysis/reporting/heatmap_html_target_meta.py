from __future__ import annotations

import logging
import re
from pathlib import Path
from types import ModuleType
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, cast

from analysis.reporting.target_annotation_groups import (
    ADME_DEFAULT,
    SAFETY_DEFAULT,
    catalog_orders,
    resolve_target_annotations,
    target_annotation_block_summary,
)
from analysis.reporting.value_utils import (
    normalize_side_effect_label,
    normalize_text,
    normalized_side_effect_list,
)
from analysis.reporting.heatmap_html_runtime import heatmap_call
from analysis.target_ids import parse_target_id, pdb_id_from_target_id
from protein_prep.pdb_records import local_uniprots_from_pdb

_LOG = logging.getLogger("heatmap-html")
_TARGET_TISSUE_EXPRESSION_KEYS = (
    "target_tissue_expression",
    "target_tissue_expression_label",
    "target_tissue_expression_score",
    "target_tissue_expression_source",
    "target_tissue_expression_summary",
)
_TARGET_FAMILY_ORDER = [
    "Kinase",
    "GPCR",
    "Ion channel",
    "Nuclear receptor",
    "Transporter",
    "Protease",
    "Enzyme",
    "Other",
    "Unassigned",
]
_UNIPROT_ACCESSION_RE = re.compile(
    r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])(?:-\d+)?\b"
)

def _infer_target_family(target_name: str) -> str:
    name = normalize_text(target_name)
    if not name:
        return "Unassigned"
    lowered = name.lower()
    if "channel" in lowered:
        return "Ion channel"
    try:
        from docking import druggability_runtime

        _tier, _triggers, protein_class = (
            druggability_runtime._infer_family_prior_from_metadata(  # type: ignore[attr-defined]
                uniprot_id=None,
                protein_name=name,
                protein_family=None,
                ec_numbers=[],
                pdb_header_text="",
                logger=_LOG,
            )
        )
    except Exception:
        protein_class = None
    class_token = normalize_text(protein_class).upper()
    if class_token == "KINASE":
        return "Kinase"
    if class_token == "GPCR":
        return "GPCR"
    if class_token == "NUCLEAR_RECEPTOR":
        return "Nuclear receptor"
    if class_token in {"ABC_TRANSPORTER", "TRANSPORTER"}:
        return "Transporter"
    if "PROTEASE" in class_token:
        return "Protease"
    if class_token in {
        "METALLOENZYME",
        "PHOSPHODIESTERASE",
        "BETA_LACTAMASE",
        "OTHER_ENZYME",
        "ENZYME",
        "CYSTEINE_PROTEASE",
    }:
        return "Enzyme"
    if "kinase" in lowered:
        return "Kinase"
    if "receptor" in lowered:
        return "Nuclear receptor" if "nuclear" in lowered else "Other"
    if any(
        keyword in lowered
        for keyword in (
            "enzyme",
            "dehydrogenase",
            "oxidase",
            "reductase",
            "hydrolase",
            "synthase",
            "isomerase",
            "transferase",
            "lyase",
            "carboxylase",
            "esterase",
        )
    ):
        return "Enzyme"
    return "Other"


def _import_pathway_resolver() -> ModuleType | None:
    """Return ``analysis.pathway_resolver`` if available."""
    try:
        from analysis import pathway_resolver

        return pathway_resolver
    except Exception as exc:
        _LOG.warning("[heatmap.target.pathways.import_failed] err=%s", exc)
        return None


def _fetch_target_uniprots_by_pdb(
    repo_root: Path, pdb_ids: Iterable[str]
) -> Dict[str, List[str]]:
    memberships: Dict[str, List[str]] = {}
    unique_pdb_ids = sorted(
        {normalize_text(pdb_id).upper() for pdb_id in pdb_ids if normalize_text(pdb_id)}
    )
    pathway_resolver = _import_pathway_resolver()
    if pathway_resolver is None:
        for pdb_id in unique_pdb_ids:
            memberships[pdb_id] = local_uniprots_from_pdb(repo_root, pdb_id)
        return memberships

    cache = pathway_resolver.Cache(
        cache_dir=repo_root / "pathways" / "cache",
        refresh=False,
        logger=_LOG,
    )
    http_client = pathway_resolver.HttpClient()
    for pdb_id in unique_pdb_ids:
        uniprot_ids: Set[str] = set(local_uniprots_from_pdb(repo_root, pdb_id))
        try:
            mapped = pathway_resolver.map_pdb_to_uniprots(pdb_id, cache, http_client)
            for uniprot in sorted({normalize_text(item).upper() for item in mapped}):
                if uniprot:
                    uniprot_ids.add(uniprot)
        except Exception as exc:
            _LOG.warning(
                "[heatmap.target.uniprots.failed] pdb=%s err=%s",
                pdb_id,
                exc,
            )
        memberships[pdb_id] = sorted(uniprot_ids, key=lambda value: value.lower())
    return memberships


def _fetch_pathway_memberships_by_pdb(
    repo_root: Path,
    pdb_ids: Iterable[str],
    *,
    uniprots_by_pdb: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, List[str]]:
    memberships: Dict[str, List[str]] = {}
    unique_pdb_ids = sorted(
        {normalize_text(pdb_id).upper() for pdb_id in pdb_ids if normalize_text(pdb_id)}
    )
    pathway_resolver = _import_pathway_resolver()
    if pathway_resolver is None:
        for pdb_id in unique_pdb_ids:
            memberships[pdb_id] = []
        return memberships

    cache = pathway_resolver.Cache(
        cache_dir=repo_root / "pathways" / "cache",
        refresh=False,
        logger=_LOG,
    )
    http_client = pathway_resolver.HttpClient()
    resolved_uniprots = uniprots_by_pdb or heatmap_call(
        "_fetch_target_uniprots_by_pdb", repo_root, unique_pdb_ids
    )
    for pdb_id in unique_pdb_ids:
        pathway_names: Set[str] = set()
        try:
            for uniprot in sorted(
                {normalize_text(item).upper() for item in resolved_uniprots.get(pdb_id, [])}
            ):
                if not uniprot:
                    continue
                entries = pathway_resolver.map_uniprot_to_reactome_pathways(
                    uniprot,
                    cache,
                    http_client,
                )
                for entry in entries:
                    name = normalize_text(entry.get("name") or entry.get("stId"))
                    if name:
                        pathway_names.add(name)
        except Exception as exc:
            _LOG.warning(
                "[heatmap.target.pathways.failed] pdb=%s err=%s",
                pdb_id,
                exc,
            )
        memberships[pdb_id] = sorted(pathway_names, key=lambda value: value.lower())
    return memberships


def _build_target_grouping_meta(
    repo_root: Path,
    rows: Iterable[Dict[str, Any]],
    target_labels: Iterable[str],
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any]]:
    ordered_labels = [label for label in target_labels if normalize_text(label)]
    label_set = set(ordered_labels)
    target_meta: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        target_label = normalize_text(row.get("target_label"))
        if target_label not in label_set:
            continue
        raw_target = normalize_text(row.get("target_id"))
        pdb_id = normalize_text(row.get("pdb_id")).upper() or pdb_id_from_target_id(
            raw_target
        )
        target_name = normalize_text(row.get("target_name"))
        current = target_meta.get(target_label)
        if current is None:
            target_meta[target_label] = {
                "target_raw": raw_target,
                "pdb": pdb_id,
                "target_name": target_name,
                "protein_family": _infer_target_family(target_name),
            }
            continue
        if target_name and not normalize_text(current.get("target_name")):
            current["target_name"] = target_name
            current["protein_family"] = _infer_target_family(target_name)
        if pdb_id and not normalize_text(current.get("pdb")):
            current["pdb"] = pdb_id
        if raw_target and not normalize_text(current.get("target_raw")):
            current["target_raw"] = raw_target

    orders = catalog_orders(repo_root)
    uniprots_by_pdb = heatmap_call(
        "_fetch_target_uniprots_by_pdb",
        repo_root,
        [str(meta.get("pdb") or "") for meta in target_meta.values()],
    )
    memberships_by_pdb = heatmap_call(
        "_fetch_pathway_memberships_by_pdb",
        repo_root,
        [str(meta.get("pdb") or "") for meta in target_meta.values()],
        uniprots_by_pdb=uniprots_by_pdb,
    )
    pathway_target_counts: Dict[str, int] = {}
    family_counts: Dict[str, int] = {}
    for label in ordered_labels:
        meta = target_meta.setdefault(label, {})
        family = normalize_text(meta.get("protein_family")) or "Unassigned"
        meta["protein_family"] = family
        family_counts[family] = family_counts.get(family, 0) + 1
        pdb_id = normalize_text(meta.get("pdb")).upper()
        uniprots = uniprots_by_pdb.get(pdb_id, [])
        meta["uniprot_accessions"] = list(uniprots)
        annotation_meta = resolve_target_annotations(
            repo_root,
            target_name=normalize_text(meta.get("target_name")),
            uniprots=uniprots,
            gene_symbols=[],
        )
        meta["adme_category"] = (
            normalize_text(annotation_meta.get("adme_category")) or ADME_DEFAULT
        )
        safety_buckets = normalized_side_effect_list(
            annotation_meta.get("safety_buckets") or []
        )
        if not safety_buckets:
            safety_buckets = [SAFETY_DEFAULT]
        meta["safety_buckets"] = safety_buckets
        pds = normalize_side_effect_label(
            annotation_meta.get("primary_display_safety")
        )
        meta["primary_display_safety"] = pds or safety_buckets[0]
        meta["secondary_safety_buckets"] = normalized_side_effect_list(
            annotation_meta.get("secondary_safety_buckets") or []
        )
        meta["direct_safety_buckets"] = normalized_side_effect_list(
            annotation_meta.get("direct_safety_buckets") or []
        )
        meta["drug_ae_buckets"] = normalized_side_effect_list(
            annotation_meta.get("drug_ae_buckets") or []
        )
        meta["safety_confidence"] = (
            normalize_text(annotation_meta.get("safety_confidence")) or "unassigned"
        )
        meta["safety_sources"] = [
            normalize_text(value)
            for value in cast(List[Any], annotation_meta.get("safety_sources") or [])
            if normalize_text(value)
        ]
        meta["safety_bucket_scores"] = cast(
            Dict[str, Any], annotation_meta.get("safety_bucket_scores") or {}
        )
        meta["safety_evidence_summary"] = normalize_side_effect_label(
            annotation_meta.get("safety_evidence_summary")
        )
        meta["safety_evidence_terms"] = normalized_side_effect_list(
            annotation_meta.get("safety_evidence_terms") or []
        )
        meta["direct_liability_examples"] = normalized_side_effect_list(
            annotation_meta.get("direct_liability_examples") or []
        )
        meta["raw_adverse_event_examples"] = normalized_side_effect_list(
            annotation_meta.get("raw_adverse_event_examples") or []
        )
        for expression_key in _TARGET_TISSUE_EXPRESSION_KEYS:
            expression_value = annotation_meta.get(expression_key)
            if expression_value not in (None, ""):
                meta[expression_key] = expression_value
        pathways = memberships_by_pdb.get(pdb_id, [])
        pathway_list = sorted(
            {normalize_text(name) for name in pathways if normalize_text(name)},
            key=lambda value: value.lower(),
        )
        meta["pathway_memberships"] = pathway_list
        for pathway_name in pathway_list:
            pathway_target_counts[pathway_name] = pathway_target_counts.get(pathway_name, 0) + 1

    pathway_block_counts: Dict[str, int] = {}
    for label in ordered_labels:
        meta = target_meta.setdefault(label, {})
        memberships = cast(List[str], meta.get("pathway_memberships") or [])
        if memberships:
            primary_display_pathway = min(
                memberships,
                key=lambda name: (
                    -int(pathway_target_counts.get(name, 0)),
                    name.lower(),
                ),
            )
        else:
            primary_display_pathway = "Unassigned"
        meta["primary_display_pathway"] = primary_display_pathway
        pathway_block_counts[primary_display_pathway] = (
            pathway_block_counts.get(primary_display_pathway, 0) + 1
        )

    family_order_map = {name: idx for idx, name in enumerate(_TARGET_FAMILY_ORDER)}
    family_blocks = [
        {"key": name, "label": name, "count": count}
        for name, count in sorted(
            family_counts.items(),
            key=lambda item: (
                family_order_map.get(item[0], len(family_order_map)),
                item[0].lower(),
            ),
        )
    ]
    pathway_blocks = [
        {"key": name, "label": name, "count": count}
        for name, count in sorted(
            pathway_block_counts.items(),
            key=lambda item: (
                item[0] == "Unassigned",
                -item[1],
                item[0].lower(),
            ),
        )
    ]
    adme_blocks = target_annotation_block_summary(
        ordered_labels=ordered_labels,
        target_meta=target_meta,
        field_name="adme_category",
        order=orders.get("adme", [ADME_DEFAULT]),
        default_value=ADME_DEFAULT,
    )
    safety_blocks = target_annotation_block_summary(
        ordered_labels=ordered_labels,
        target_meta=target_meta,
        field_name="primary_display_safety",
        order=orders.get("safety", [SAFETY_DEFAULT]),
        default_value=SAFETY_DEFAULT,
    )
    return target_meta, {
        "default_mode": "none",
        "family_blocks": family_blocks,
        "pathway_blocks": pathway_blocks,
        "adme_blocks": adme_blocks,
        "safety_blocks": safety_blocks,
    }


def _target_expression_report_fields(grouping_meta: Dict[str, Any]) -> Dict[str, str]:
    return {
        key: normalize_text(grouping_meta.get(key))
        for key in _TARGET_TISSUE_EXPRESSION_KEYS
    }


def _resolve_target_display_label(row: Dict[str, Any], target_id: str) -> str:
    pdb_id = normalize_text(row.get("pdb_id")).upper()
    if not pdb_id:
        pdb_id = pdb_id_from_target_id(target_id)
    target_name = normalize_text(row.get("target_name"))
    variant = normalize_text(row.get("variant"))
    ph_label = normalize_text(row.get("ph_label"))
    if (not variant or not ph_label) and "|" in target_id:
        tid = parse_target_id(target_id)
        if not variant:
            variant = normalize_text(tid.variant)
        if not ph_label:
            ph_label = normalize_text(tid.ph_label)
    variant_ph = " | ".join([part for part in (variant, ph_label) if part])
    if target_name and pdb_id:
        lines = [pdb_id, target_name]
        if variant_ph:
            lines.append(variant_ph)
        return "\n".join(lines)
    if pdb_id and variant_ph:
        return f"{pdb_id}\n{variant_ph}"
    if pdb_id:
        return pdb_id
    if target_name:
        return target_name
    return target_id


def _build_target_display_map(rows: List[Dict[str, Any]]) -> Dict[str, str]:
    labels: Dict[str, str] = {}
    has_name: Dict[str, bool] = {}
    for row in rows:
        target_id = normalize_text(row.get("target_id"))
        if not target_id:
            continue
        candidate = _resolve_target_display_label(row, target_id)
        candidate_has_name = bool(normalize_text(row.get("target_name")))
        if target_id not in labels:
            labels[target_id] = candidate
            has_name[target_id] = candidate_has_name
            continue
        # Prefer labels that include target_name metadata when available.
        if candidate_has_name and not has_name.get(target_id, False):
            labels[target_id] = candidate
            has_name[target_id] = True
    return labels
