from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from analysis.reporting.ligand_annotation_groups import resolve_ligand_annotations
from analysis.reporting.target_annotation_groups import resolve_target_annotations

from analysis._common import clean_text, parse_binary, parse_float, read_csv_rows, repo_root


_UNIPROT_RE = re.compile(
    r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])(?:-\d+)?\b"
)


@dataclass(frozen=True)
class MetadataResult:
    value: Any
    source: str = ""
    confidence: str = ""
    missing_reason: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "value": self.value,
            "source": self.source,
            "confidence": self.confidence,
            "missing_reason": self.missing_reason,
        }


class AnnotationIndex:
    def __init__(self, rows: list[dict[str, Any]] | None = None) -> None:
        self.rows = rows or []
        self.by_drug = self._index("drug_id", "ligand_id", "ligand_base", "ligand_display")
        self.by_target = self._index("target_id", "pdb_id", "gene_symbol", "uniprot")
        self.by_pair = self._pair_index()

    @classmethod
    def from_path(cls, path: Path | None) -> "AnnotationIndex":
        if path is None or not path.exists():
            return cls()
        rows: list[dict[str, Any]] = []
        paths = [path] if path.is_file() else sorted(p for p in path.iterdir() if p.suffix.lower() in {".csv", ".json"})
        for item in paths:
            if item.suffix.lower() == ".csv":
                rows.extend(read_csv_rows(item))
            elif item.suffix.lower() == ".json":
                payload = json.loads(item.read_text(encoding="utf-8"))
                if isinstance(payload, list):
                    rows.extend(row for row in payload if isinstance(row, dict))
                elif isinstance(payload, dict):
                    for key in ("rows", "annotations", "data"):
                        value = payload.get(key)
                        if isinstance(value, list):
                            rows.extend(row for row in value if isinstance(row, dict))
        return cls(rows)

    def _index(self, *keys: str) -> dict[str, dict[str, Any]]:
        out: dict[str, dict[str, Any]] = {}
        for row in self.rows:
            for key in keys:
                value = clean_text(row.get(key)).lower()
                if value and value not in out:
                    out[value] = row
        return out

    def _pair_index(self) -> dict[tuple[str, str], dict[str, Any]]:
        out: dict[tuple[str, str], dict[str, Any]] = {}
        for row in self.rows:
            drug = clean_text(row.get("drug_id") or row.get("ligand_id") or row.get("ligand_base")).lower()
            target = clean_text(row.get("target_id") or row.get("pdb_id")).lower()
            if drug and target:
                out[(drug, target)] = row
        return out

    def drug(self, drug_id: str) -> Mapping[str, Any]:
        return self.by_drug.get(clean_text(drug_id).lower(), {})

    def target(self, target_id: str, pdb_id: str = "") -> Mapping[str, Any]:
        return self.by_target.get(clean_text(target_id).lower(), {}) or self.by_target.get(clean_text(pdb_id).lower(), {})

    def pair(self, drug_id: str, target_id: str, pdb_id: str = "") -> Mapping[str, Any]:
        drug_key = clean_text(drug_id).lower()
        return self.by_pair.get((drug_key, clean_text(target_id).lower()), {}) or self.by_pair.get((drug_key, clean_text(pdb_id).lower()), {})


def _result(value: Any, source: str, confidence: str = "", missing_reason: str = "") -> dict[str, Any]:
    return MetadataResult(value=value, source=source, confidence=confidence, missing_reason=missing_reason).as_dict()


def get_drug_metadata(drug_id: str, annotations: AnnotationIndex | None = None) -> dict[str, Any]:
    row = (annotations or AnnotationIndex()).drug(drug_id)
    value = row if row else {"drug_id": drug_id}
    return _result(value, "annotation_table" if row else "input", "high" if row else "identity_only")


def get_target_metadata(target_id: str, pdb_id: str, annotations: AnnotationIndex | None = None) -> dict[str, Any]:
    ann = annotations or AnnotationIndex()
    row = dict(ann.target(target_id, pdb_id))
    target_text = " ".join(v for v in (target_id, pdb_id, clean_text(row.get("target_name"))) if v)
    uniprots = _UNIPROT_RE.findall(target_text)
    try:
        resolved = resolve_target_annotations(
            repo_root(),
            target_name=clean_text(row.get("target_name")) or target_id or pdb_id,
            uniprots=uniprots,
            gene_symbols=[clean_text(row.get("gene_symbol"))] if row.get("gene_symbol") else [],
        )
        row.update(resolved)
        source = "annotation_table+target_annotation_catalog" if ann.rows else "target_annotation_catalog"
        confidence = clean_text(resolved.get("safety_confidence")) or "heuristic"
    except Exception:
        source = "annotation_table" if row else ""
        confidence = "low" if row else ""
    return _result(row, source, confidence, "" if row else "not_found")


def get_pk_metadata(drug_id: str, annotations: AnnotationIndex | None = None) -> dict[str, Any]:
    row = (annotations or AnnotationIndex()).drug(drug_id)
    value = clean_text(row.get("free_cmax") or row.get("free_cmax_um") or row.get("cmax_free"))
    return _result(value, "annotation_table" if value else "", "medium" if value else "", "" if value else "missing")


def get_tissue_expression(target_id: str, annotations: AnnotationIndex | None = None) -> dict[str, Any]:
    row = (annotations or AnnotationIndex()).target(target_id)
    value = clean_text(
        row.get("tissue_expression")
        or row.get("target_tissue_expression_score")
        or row.get("target_tissue_expression")
    )
    return _result(value, "annotation_table" if value else "", "medium" if value else "", "" if value else "missing")


def get_target_adr_evidence(target_id: str, annotations: AnnotationIndex | None = None) -> dict[str, Any]:
    row = (annotations or AnnotationIndex()).target(target_id)
    explicit = parse_binary(row.get("target_adr_evidence"))
    if explicit is not None:
        return _result(explicit, "annotation_table", "high")
    safety = clean_text(row.get("primary_display_safety") or row.get("safety_buckets"))
    value = 1 if safety and safety.lower() not in {"unassigned", "none", "unknown"} else 0
    return _result(value, "target_annotation_catalog" if safety else "", "heuristic" if safety else "low")


def get_pathway_evidence(target_id: str, annotations: AnnotationIndex | None = None) -> dict[str, Any]:
    row = (annotations or AnnotationIndex()).target(target_id)
    explicit = parse_binary(row.get("pathway_evidence"))
    if explicit is not None:
        return _result(explicit, "annotation_table", "high")
    pathway = clean_text(row.get("pathway") or row.get("reactome_pathway") or row.get("aop_pathway"))
    return _result(1 if pathway else 0, "annotation_table" if pathway else "", "medium" if pathway else "low")


def get_ligand_chemotype(drug_id: str, annotations: AnnotationIndex | None = None) -> dict[str, Any]:
    row = (annotations or AnnotationIndex()).drug(drug_id)
    value = clean_text(row.get("ligand_chemotype") or row.get("chemotype_primary") or row.get("chemotype_series"))
    if value:
        return _result(value, "annotation_table", "high")
    try:
        resolved = resolve_ligand_annotations(repo_root(), ligand_name=drug_id, ligand_base=drug_id)
        value = clean_text(resolved.get("chemotype_primary"))
        return _result(value, "ligand_annotation_catalog", "heuristic", "" if value else "missing")
    except Exception:
        return _result("", "", "", "missing")


def get_literature_supported_label(
    drug_id: str, target_id: str, annotations: AnnotationIndex | None = None
) -> dict[str, Any]:
    row = (annotations or AnnotationIndex()).pair(drug_id, target_id)
    value = clean_text(row.get("literature_supported_label") or row.get("label") or row.get("supported"))
    return _result(value, "annotation_table" if value else "", "high" if value else "", "" if value else "missing")


def numeric_feature(value: Any) -> float | None:
    binary = parse_binary(value)
    if binary is not None:
        return float(binary)
    return parse_float(value)

