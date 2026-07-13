from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.build_ml_dataset import add_standard_ml_columns, deduplicate_ml_rows
from analysis.ml.labels import binary_label_series
from analysis.ml.mechanism_label_projection import project_mechanism_labels_for_run_master
from analysis.ml.spd_label_enrichment import enrich_spd_labels_for_run_master
from analysis.ml.tissue_expression_projection import add_tissue_expression_context
from analysis.ml.source_benchmark_tables import ID_COLS, MODEL_READY_PROVENANCE_COLS

BINDING_LABEL_POLICY_VERSION = "spd_binding_censor_aware_v1"

SCORE_FEATURE_COLS = [
    "atlas_score",
    "consensus_score",
    "banana_score",
    "banana_binding_probability",
    "banana_score_normalized",
    "banana_atlas_blend_score",
    "binding_expert_score",
    "structure_quality",
    "protein_class",
    "target_family",
    "ligand_chemotype",
    "scaffold_key",
]
EXPOSURE_CONTEXT_COLS = [
    "free_cmax_um",
    "free_cmax_uM",
    "free_cmax",
    "cmax_um",
    "fraction_unbound_plasma",
    "exposure_margin",
    "spd_ac50_uM",
    "spd_exposure_relevant",
    "spd_exposure_weak",
    "spd_exposure_unlikely",
]
TISSUE_CONTEXT_COLS = [
    "adr_id",
    "adr_term",
    "adr_pt",
    "adr_soc",
    "meddra_pt",
    "meddra_soc",
    "adr_site_group",
    "adr_site_mapping_source",
    "adr_site_mapping_confidence",
    "adr_site_terms_matched",
    "adr_site_match_score",
    "target_adr_terms_aggregated",
    "target_adr_sources_aggregated",
    "target_adr_evidence_count",
    "site_id",
    "site_name",
    "uberon_id",
    "site_relevance_score",
    "site_relevance_score_by_adr",
    "site_safety_risk_score",
    "expression_presence_score",
    "site_specificity_score",
    "expression_concordance_score",
    "distribution_penalty",
    "sensitive_offsite_expression_score",
    "hpa_consensus_ntpm",
    "hpa_specificity_score",
    "gtex_median_tpm",
    "bgee_expression_score",
    "ot_expression_value",
    "ot_expression_zscore",
    "hpa_max_tissue",
    "gtex_max_tissue",
    "bgee_max_tissue",
    "ot_expression_max_tissue",
    "ot_safety_biosample_tissues",
    "ot_safety_liability_count",
]
MECHANISM_CONTEXT_COLS = [
    "mechanism_graph_score",
    "mechanism_path_count",
    "drug_adr_known",
    "target_adr_known",
    "target_pathway_adr_link",
    "triad_complete",
    "target_adr_evidence",
    "pathway_evidence",
    "mechanism_label_status",
    "negative_evidence_type",
    "negative_source",
    "negative_confidence",
    "mechanism_label_source",
    "mechanism_label_projection_key",
    "mechanism_label_projection_source",
    "mechanism_label_projection_source_rows",
    "projected_label_source",
    "projected_source_family",
    "projected_upstream_source",
    "projected_source_label_policy",
]
RAW_KEEP_COLS = [
    *ID_COLS,
    *SCORE_FEATURE_COLS,
    *EXPOSURE_CONTEXT_COLS,
    *TISSUE_CONTEXT_COLS,
    *MECHANISM_CONTEXT_COLS,
    "z_selected",
    "z_selected_source",
    "final_score",
    "SCORCH_score_used",
    "spd_label_status",
    "spd_assay_count",
    "spd_assay_ids",
    "spd_assay_name",
    "ml_binary_label",
    "ml_supervised_eligible",
    "ml_exclude_reason",
    "dedup_drug_key",
    "dedup_target_key",
]
MODEL_READY_COLS = [
    *ID_COLS,
    "source_objective",
    *MODEL_READY_PROVENANCE_COLS,
    *SCORE_FEATURE_COLS,
    *EXPOSURE_CONTEXT_COLS,
    *TISSUE_CONTEXT_COLS,
    *MECHANISM_CONTEXT_COLS,
]


def _numeric(series: pd.Series | None, index: pd.Index) -> pd.Series:
    if series is None:
        return pd.Series(pd.NA, index=index, dtype="Float64")
    return pd.to_numeric(series, errors="coerce").astype("Float64")


def _first_label(df: pd.DataFrame, columns: list[str]) -> pd.Series:
    label = pd.Series(pd.NA, index=df.index, dtype="Int64")
    for col in columns:
        if col not in df.columns:
            continue
        parsed = pd.to_numeric(df[col], errors="coerce")
        parsed = parsed.where(parsed.isin([-1, 0, 1]))
        label = label.where(label.notna(), parsed)
    return label.astype("Int64")


def _base_table(source: pd.DataFrame) -> pd.DataFrame:
    source = source.loc[:, ~source.columns.duplicated()].copy()
    out = source[[col for col in RAW_KEEP_COLS if col in source.columns]].copy()
    out = out.loc[:, ~out.columns.duplicated()].copy()
    if "drug_id" not in out.columns and "dedup_drug_key" in out.columns:
        out["drug_id"] = out["dedup_drug_key"]
    if "target_id" not in out.columns and "dedup_target_key" in out.columns:
        out["target_id"] = out["dedup_target_key"]
    if "atlas_score" not in out.columns and "z_selected" in out.columns:
        out["atlas_score"] = pd.to_numeric(out["z_selected"], errors="coerce")
    if "consensus_score" not in out.columns and "final_score" in out.columns:
        out["consensus_score"] = pd.to_numeric(out["final_score"], errors="coerce")
    if "free_cmax_um" not in out.columns:
        for col in ("free_cmax_uM", "free_cmax"):
            if col in out.columns:
                out["free_cmax_um"] = pd.to_numeric(out[col], errors="coerce")
                break
    return out


def _add_common_provenance(
    table: pd.DataFrame,
    *,
    objective: str,
    label_col: str,
    endpoint_type: str,
    policy: str,
) -> pd.DataFrame:
    out = table.copy()
    out["source_objective"] = objective
    out["label_source"] = objective
    out["assay_type"] = "secondary_pharmacology"
    out["endpoint_type"] = endpoint_type
    out["activity_type"] = endpoint_type
    out["source_family"] = "spd"
    out["upstream_source"] = "novartis_spd"
    out["source_label_policy"] = policy
    out["database_release_year"] = 2023
    out["database_release_year_source"] = "Sutherland et al. 2023 SPD publication year; coarse dataset provenance"
    out["availability_year"] = 2023
    out["availability_year_source"] = "SPD public release/publication year; not row-level assay timing"
    out["_expert_label_col"] = label_col
    return out


def _binding_label(source: pd.DataFrame, *, active_um: float, inactive_um: float) -> pd.Series:
    direct = _first_label(source, ["spd_binding_label", "spd_activity_label", "activity_label"])
    ac50 = _numeric(source.get("spd_ac50_uM"), source.index)
    relation = pd.Series("", index=source.index, dtype="string")
    for column in (
        "spd_activity_relation",
        "activity_relation",
        "standard_relation",
        "relation",
    ):
        if column not in source.columns:
            continue
        values = source[column].astype("string").str.strip()
        relation = relation.where(relation.ne(""), values)
    relation = relation.replace(
        {
            "≤": "<=",
            "≥": ">=",
            "==": "=",
            "eq": "=",
            "lt": "<",
            "gt": ">",
        }
    )

    derived = pd.Series(pd.NA, index=source.index, dtype="Int64")
    exact = relation.eq("=")
    upper_bound = relation.isin(["<", "<="])
    lower_bound = relation.isin([">", ">="])
    active = ((exact | upper_bound) & ac50.le(active_um)).fillna(False)
    inactive = ((exact | lower_bound) & ac50.ge(inactive_um)).fillna(False)
    derived = derived.mask(active, 1)
    derived = derived.mask(inactive, 0)
    interpretable_measurement = ac50.notna() & (exact | upper_bound | lower_bound)
    return direct.mask(interpretable_measurement, derived).astype("Int64")


def _exposure_label(source: pd.DataFrame) -> pd.Series:
    label = _first_label(source, ["spd_exposure_label", "ml_binary_label", "spd_exposure_relevant"])
    if label.notna().any():
        return label.mask(label.eq(-1), pd.NA).astype("Int64")
    margin = _numeric(source.get("exposure_margin"), source.index)
    derived = pd.Series(pd.NA, index=source.index, dtype="Int64")
    relevant = margin.le(10).fillna(False)
    not_relevant = margin.gt(10).fillna(False)
    derived = derived.mask(relevant, 1)
    derived = derived.mask(not_relevant, 0)
    return derived.astype("Int64")


def _tissue_label(source: pd.DataFrame) -> pd.Series:
    label = _first_label(
        source,
        ["tissue_site_label", "tissue_relevance_label", "site_relevance_label", "target_site_relevance_label"],
    )
    return label.mask(label.eq(-1), pd.NA).astype("Int64")


def _mechanism_label(source: pd.DataFrame) -> pd.Series:
    label = _first_label(
        source,
        ["mechanism_ml_label", "mechanism_label", "four_state_ml_label", "drug_target_adr_mechanism_label"],
    )
    return label.mask(label.eq(-1), pd.NA).astype("Int64")


def _write_raw_and_model_ready(
    table: pd.DataFrame,
    *,
    label_col: str,
    raw_path: Path,
    model_ready_path: Path,
) -> dict[str, Any]:
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(raw_path, index=False)
    label = binary_label_series(table[label_col]) if label_col in table.columns else pd.Series(pd.NA, index=table.index)
    summary: dict[str, Any] = {
        "raw_path": str(raw_path),
        "model_ready_path": None,
        "rows": int(len(table)),
        "label_col": label_col,
        "labelable_rows": int(label.notna().sum()),
        "positive_rows": int(label.eq(1).sum()),
        "negative_rows": int(label.eq(0).sum()),
        "unknown_rows": int(label.isna().sum()),
    }
    supervised = table[label.notna()].copy()
    if supervised.empty:
        summary["status"] = "raw_only_no_labelable_rows"
        return summary
    supervised[label_col] = label.loc[supervised.index].astype(int)
    supervised = add_standard_ml_columns(supervised)
    supervised = deduplicate_ml_rows(supervised, label_col)
    keep: list[str] = []
    seen: set[str] = set()
    for col in [*MODEL_READY_COLS, label_col]:
        if col in supervised.columns and col not in seen:
            keep.append(col)
            seen.add(col)
    model_ready = supervised[keep].copy()
    model_ready_path.parent.mkdir(parents=True, exist_ok=True)
    model_ready.to_csv(model_ready_path, index=False)
    summary.update(
        {
            "status": "model_ready_written",
            "model_ready_path": str(model_ready_path),
            "model_ready_rows": int(len(model_ready)),
        }
    )
    return summary


def build_spd_four_expert_tables(
    spd_table_path: str | Path,
    out_dir: str | Path,
    *,
    active_um: float = 1.0,
    inactive_um: float = 10.0,
    spd_panel_path: str | Path | None = None,
    target_map_path: str | Path | None = None,
    ligand_map_path: str | Path | None = None,
    target_metadata_path: str | Path | None = None,
    mechanism_label_path: str | Path | None = None,
    target_expression_path: str | Path | None = None,
    run_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Build objective-specific SPD tables for the four Atlas experts.

    The function is intentionally conservative: it writes raw objective tables even
    when a label is absent, but only writes model-ready tables for labelable rows.
    Missing tissue or mechanism labels remain unknown instead of being converted
    to negatives.
    """

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(spd_table_path, low_memory=False)
    source, enrichment_summary = enrich_spd_labels_for_run_master(
        source,
        spd_panel_path=spd_panel_path,
        target_map_path=target_map_path,
        ligand_map_path=ligand_map_path,
        target_metadata_path=target_metadata_path,
    )
    source, mechanism_projection_summary = project_mechanism_labels_for_run_master(
        source,
        mechanism_label_path=mechanism_label_path,
        run_dir=run_dir,
    )
    source, tissue_expression_summary = add_tissue_expression_context(
        source,
        target_expression_path=target_expression_path,
        run_dir=run_dir,
    )
    base = _base_table(source)
    tables: dict[str, pd.DataFrame] = {}

    binding = base.copy()
    binding["spd_binding_label"] = _binding_label(source, active_um=active_um, inactive_um=inactive_um)
    binding["spd_binding_label_policy_version"] = BINDING_LABEL_POLICY_VERSION
    binding = _add_common_provenance(
        binding,
        objective="spd_binding_activity",
        label_col="spd_binding_label",
        endpoint_type="AC50_potency",
        policy=(
            f"SPD censor-aware activity: exact or upper-bound AC50 <= {active_um:g} uM "
            f"positive; exact or lower-bound AC50 >= {inactive_um:g} uM negative; "
            "bounds that cross a threshold and the intermediate range remain unknown."
        ),
    )
    tables["binding"] = binding

    exposure = base.copy()
    exposure["spd_exposure_label"] = _exposure_label(source)
    exposure = _add_common_provenance(
        exposure,
        objective="spd_exposure_relevance",
        label_col="spd_exposure_label",
        endpoint_type="AC50/free_cmax_margin",
        policy="SPD exposure relevance: AC50/free-Cmax <= 10 positive; missing margin unknown. AC50/free-Cmax fields are label-definition fields, not predictors for nonleaky models.",
    )
    tables["exposure"] = exposure

    tissue = base.copy()
    tissue["tissue_site_label"] = _tissue_label(source)
    tissue = _add_common_provenance(
        tissue,
        objective="tissue_site_relevance",
        label_col="tissue_site_label",
        endpoint_type="target_expression_in_adr_site",
        policy="Tissue labels must come from explicit ADR-site or curated organ-panel labels. Absence of expression/ADR-site evidence remains unknown.",
    )
    tables["tissue"] = tissue

    mechanism = base.copy()
    mechanism["mechanism_ml_label"] = _mechanism_label(source)
    mechanism = _add_common_provenance(
        mechanism,
        objective="mechanism_adr_support",
        label_col="mechanism_ml_label",
        endpoint_type="drug_target_adr_mechanism_evidence",
        policy="Mechanism positives require curated target/drug/pathway ADR support. Negatives must be measured/reliable controls; missing mechanism evidence remains unknown.",
    )
    if "mechanism_label_source" in mechanism.columns:
        source_values = mechanism["mechanism_label_source"].fillna("").astype(str).str.strip()
        mechanism["label_source"] = mechanism["label_source"].where(source_values.eq(""), source_values)
    tables["mechanism"] = mechanism

    summaries: dict[str, Any] = {
        "source_table": str(spd_table_path),
        "out_dir": str(out),
        "purpose": "Prepare objective-specific SPD tables for Atlas binding, exposure, tissue, and mechanism experts.",
        "active_um": active_um,
        "inactive_um": inactive_um,
        "spd_auto_enrichment": enrichment_summary,
        "mechanism_label_projection": mechanism_projection_summary,
        "tissue_expression_projection": tissue_expression_summary,
        "experts": {},
    }
    for name, table in tables.items():
        label_col = str(table["_expert_label_col"].iloc[0])
        summaries["experts"][name] = _write_raw_and_model_ready(
            table.drop(columns=["_expert_label_col"]),
            label_col=label_col,
            raw_path=out / f"ml_spd_{name}_table.csv",
            model_ready_path=out / "model_ready" / f"ml_spd_{name}_model_ready.csv",
        )
    manifest_path = out / "spd_four_expert_tables_manifest.json"
    manifest_path.write_text(json.dumps(summaries, indent=2, sort_keys=True), encoding="utf-8")
    return summaries
