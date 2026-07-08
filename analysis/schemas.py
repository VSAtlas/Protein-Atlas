from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from analysis._common import clean_text, parse_float


REQUIRED_PAIR_COLUMNS = [
    "drug_id",
    "target_id",
    "pdb_id",
    "atlas_score",
    "mmgbsa_score",
    "free_cmax",
    "exposure_plausibility",
    "tissue_expression",
    "target_adr_evidence",
    "pathway_evidence",
    "structure_quality",
    "protein_class",
    "ligand_chemotype",
    "literature_supported_label",
]
INTERNAL_PREFIXES = ("_source_", "_missing_", "_confidence_", "_internal_")
LEAKY_FEATURES = {"drug_id", "target_id", "pdb_id", "literature_supported_label"}


def validate_pair_table(
    rows: Sequence[Mapping[str, Any]],
    *,
    require_atlas_score: bool = True,
) -> list[str]:
    errors: list[str] = []
    if not rows:
        return errors
    columns = set().union(*(row.keys() for row in rows))
    missing = [name for name in REQUIRED_PAIR_COLUMNS if name not in columns]
    if missing:
        errors.append(f"missing required columns: {', '.join(missing)}")
    for name in columns:
        if name in REQUIRED_PAIR_COLUMNS or name in {"exposure_adjusted_score", "priority_score"}:
            continue
        if name.startswith(INTERNAL_PREFIXES):
            continue
        errors.append(f"non-standard column lacks internal prefix: {name}")
    if require_atlas_score:
        bad = sum(1 for row in rows if parse_float(row.get("atlas_score")) is None)
        if bad == len(rows):
            errors.append("atlas_score has no numeric values")
    return errors


def pair_table_schema() -> dict[str, Any]:
    return {
        "required_columns": REQUIRED_PAIR_COLUMNS,
        "internal_prefixes": list(INTERNAL_PREFIXES),
        "score_directions": {"atlas_score": "higher", "mmgbsa_score": "lower"},
        "leakage_excluded_features": sorted(LEAKY_FEATURES),
    }


def summarize_pair_table(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    def distinct(name: str) -> int:
        return len({clean_text(row.get(name)) for row in rows if clean_text(row.get(name))})

    def missing_rate(name: str) -> float:
        if not rows:
            return 0.0
        return sum(1 for row in rows if not clean_text(row.get(name))) / len(rows)

    score_ranges: dict[str, dict[str, float | None]] = {}
    for name in ("atlas_score", "mmgbsa_score", "priority_score"):
        values = [value for row in rows if (value := parse_float(row.get(name))) is not None]
        score_ranges[name] = {
            "min": min(values) if values else None,
            "max": max(values) if values else None,
            "mean": (sum(values) / len(values)) if values else None,
        }
    return {
        "n_rows": len(rows),
        "n_drugs": distinct("drug_id"),
        "n_targets": distinct("target_id"),
        "n_pdbs": distinct("pdb_id"),
        "missingness": {name: missing_rate(name) for name in REQUIRED_PAIR_COLUMNS},
        "score_ranges": score_ranges,
    }


def write_pair_table_sidecars(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.with_suffix(".schema.json").write_text(json.dumps(pair_table_schema(), indent=2), encoding="utf-8")
    path.with_suffix(".summary.json").write_text(json.dumps(summarize_pair_table(rows), indent=2), encoding="utf-8")

