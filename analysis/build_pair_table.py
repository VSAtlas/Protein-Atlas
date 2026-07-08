from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Mapping

from analysis._common import (
    clean_text,
    first_value,
    format_float,
    parse_float,
    write_csv_rows,
)
from analysis.metadata import (
    AnnotationIndex,
    get_drug_metadata,
    get_ligand_chemotype,
    get_literature_supported_label,
    get_pathway_evidence,
    get_pk_metadata,
    get_target_adr_evidence,
    get_target_metadata,
    get_tissue_expression,
    numeric_feature,
)
from analysis.priority_score import DEFAULT_WEIGHTS, compute_priority_score
from analysis.schemas import REQUIRED_PAIR_COLUMNS, write_pair_table_sidecars


REQUIRED_COLUMNS = [
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
DERIVED_COLUMNS = ["exposure_adjusted_score", "priority_score", "priority_score_components_json", "priority_score_missing_components"]
LOWER_IS_BETTER_SCORE_ALIASES = ("selected_docking_score", "vina_score", "docking_score")


def _strip_ligand_suffix(value: str) -> str:
    text = clean_text(value)
    text = Path(text).name
    text = re.sub(r"\.(pdbqt|sdf|mol2|csv)$", "", text, flags=re.IGNORECASE)
    text = re.sub(r"\.sanitized$", "", text, flags=re.IGNORECASE)
    return text


def _target_id(row: Mapping[str, Any]) -> str:
    direct = first_value(row, ("target_id", "target", "uniprot", "gene_symbol"))
    if direct:
        return direct
    pdb = first_value(row, ("pdb_id", "pdb", "target_pdb"))
    variant = first_value(row, ("variant", "apo_holo"))
    ph_label = first_value(row, ("ph_label", "pH", "ph"))
    return "|".join(part for part in (pdb, variant, ph_label) if part)


def _atlas_score(row: Mapping[str, Any]) -> tuple[str, str]:
    value = first_value(
        row,
        (
            "atlas_score",
            "z_selected",
            "z_stage2",
            "z_stage1",
            "z_vs_decoys_blend",
            "z_vs_decoys_consensus",
            "z_vs_decoys_consensus_pre",
        ),
    )
    if value:
        return value, "higher_is_better_z_score"
    for name in LOWER_IS_BETTER_SCORE_ALIASES:
        raw = parse_float(row.get(name))
        if raw is not None:
            return format_float(-raw), f"converted_from_lower_is_better:{name}"
    return "", "missing"


def _exposure_plausibility(row: Mapping[str, Any], atlas_score: str, free_cmax: str) -> str:
    explicit = first_value(row, ("exposure_plausibility", "exposure_supported"))
    if explicit:
        return explicit
    cmax = parse_float(free_cmax)
    potency = parse_float(first_value(row, ("AC50", "ac50", "Ki", "ki", "ic50", "IC50")))
    if cmax is not None and cmax > 0 and potency is not None:
        margin = potency / cmax
        if margin <= 10:
            return "1"
        if margin <= 100:
            return "0.5"
        return "0"
    score = parse_float(atlas_score)
    if cmax is not None and cmax > 0 and score is not None and score >= 2.0:
        return "0.5"
    return ""


def _field_from_meta(meta: Mapping[str, Any], *names: str) -> str:
    value = meta.get("value")
    if isinstance(value, Mapping):
        return first_value(value, names)
    return clean_text(value)


def build_pair_feature_table(
    screening_results: Path,
    annotations: Path | None,
    out_path: Path | None = None,
    config: Mapping[str, Any] | None = None,
    *,
    include_debug_metadata: bool = False,
) -> Any:
    rows = []
    ann = AnnotationIndex.from_path(annotations)
    from analysis._common import read_csv_rows

    for raw in read_csv_rows(screening_results):
        drug_id = _strip_ligand_suffix(
            first_value(raw, ("drug_id", "ligand_id", "ligand_base", "ligand_display", "ligand", "ligand_file"))
        )
        target_id = _target_id(raw)
        pdb_id = first_value(raw, ("pdb_id", "pdb", "target_pdb"))
        atlas_score, atlas_source = _atlas_score(raw)
        drug_meta = get_drug_metadata(drug_id, ann)
        target_meta = get_target_metadata(target_id, pdb_id, ann)
        pk_meta = get_pk_metadata(drug_id, ann)
        free_cmax = first_value(raw, ("free_cmax", "free_cmax_um", "cmax_free")) or clean_text(pk_meta.get("value"))
        tissue_meta = get_tissue_expression(target_id, ann)
        adr_meta = get_target_adr_evidence(target_id, ann)
        pathway_meta = get_pathway_evidence(target_id, ann)
        chemotype_meta = get_ligand_chemotype(drug_id, ann)
        literature_meta = get_literature_supported_label(drug_id, target_id, ann)
        raw_target_value = target_meta.get("value")
        target_value: Mapping[str, Any] = raw_target_value if isinstance(raw_target_value, Mapping) else {}
        row: dict[str, Any] = {
            "drug_id": drug_id,
            "target_id": target_id,
            "pdb_id": pdb_id,
            "atlas_score": atlas_score,
            "mmgbsa_score": first_value(raw, ("mmgbsa_score", "MMGBSA", "mmgbsa", "delta_g_binding")),
            "free_cmax": free_cmax,
            "exposure_plausibility": _exposure_plausibility(raw, atlas_score, free_cmax),
            "tissue_expression": first_value(raw, ("tissue_expression", "target_tissue_expression_score", "target_tissue_expression"))
            or clean_text(tissue_meta.get("value")),
            "target_adr_evidence": first_value(raw, ("target_adr_evidence", "adr_evidence"))
            or clean_text(adr_meta.get("value")),
            "pathway_evidence": first_value(raw, ("pathway_evidence", "pathway_supported"))
            or clean_text(pathway_meta.get("value")),
            "structure_quality": first_value(raw, ("structure_quality", "pdb_quality_score", "resolution_quality")),
            "protein_class": first_value(raw, ("protein_class", "target_class"))
            or first_value(target_value, ("protein_class", "adme_category", "primary_display_safety")),
            "ligand_chemotype": first_value(raw, ("ligand_chemotype", "chemotype_primary"))
            or clean_text(chemotype_meta.get("value")),
            "literature_supported_label": first_value(raw, ("literature_supported_label", "supported_label"))
            or clean_text(literature_meta.get("value")),
        }
        exposure = numeric_feature(row["exposure_plausibility"])
        score = parse_float(row["atlas_score"])
        row["exposure_adjusted_score"] = format_float(score * exposure if score is not None and exposure is not None else None)
        if include_debug_metadata:
            row.update(
                {
                    "_atlas_score_source": atlas_source,
                    "_drug_metadata_source": clean_text(drug_meta.get("source")),
                    "_target_metadata_source": clean_text(target_meta.get("source")),
                    "_pk_metadata_source": clean_text(pk_meta.get("source")),
                }
            )
        rows.append(row)

    _add_priority_scores(rows, config or {})
    if out_path is not None:
        import pandas as pd

        out = Path(out_path)
        fieldnames = [*REQUIRED_PAIR_COLUMNS, *[name for name in DERIVED_COLUMNS if any(name in row for row in rows)]]
        debug = [
            name
            for name in ("_atlas_score_source", "_drug_metadata_source", "_target_metadata_source", "_pk_metadata_source")
            if any(name in row for row in rows)
        ]
        write_csv_rows(out, rows, [*fieldnames, *debug])
        write_pair_table_sidecars(out, rows)
        return pd.DataFrame(rows)
    return rows


def _add_priority_scores(rows: list[dict[str, Any]], config: Mapping[str, Any] | None = None) -> None:
    if not rows:
        return
    import pandas as pd

    priority_config = config.get("priority_score", {}) if isinstance(config, Mapping) else {}
    weights = priority_config.get("weights", DEFAULT_WEIGHTS) if isinstance(priority_config, Mapping) else DEFAULT_WEIGHTS
    scored = compute_priority_score(pd.DataFrame(rows), weights=weights)
    for idx, row in enumerate(rows):
        for name in ("priority_score", "priority_score_components_json", "priority_score_missing_components"):
            row[name] = format_float(scored.at[idx, name]) if name == "priority_score" else clean_text(scored.at[idx, name])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build Atlas target-ligand pair feature table.")
    parser.add_argument("--screening-results", required=True, type=Path)
    parser.add_argument("--annotations", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--json-out", type=Path, default=None)
    parser.add_argument("--include-debug-metadata", action="store_true")
    args = parser.parse_args(argv)

    config: dict[str, Any] = {}
    if args.config:
        from analysis.io import load_config

        config = load_config(args.config)
    rows = build_pair_feature_table(
        args.screening_results,
        args.annotations,
        config=config,
        include_debug_metadata=args.include_debug_metadata,
    )
    extra = [name for name in DERIVED_COLUMNS if any(clean_text(row.get(name)) for row in rows)]
    debug = [name for name in ("_atlas_score_source", "_drug_metadata_source", "_target_metadata_source", "_pk_metadata_source") if any(name in row for row in rows)]
    write_csv_rows(args.out, rows, [*REQUIRED_COLUMNS, *extra, *debug])
    write_pair_table_sidecars(args.out, rows)
    json_out = args.json_out or args.out.with_suffix(".json")
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
