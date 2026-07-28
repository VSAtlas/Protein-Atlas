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
from analysis.spd_activity_policy import (
    label_spd_binding_observation,
    parse_spd_activity_interval,
)

BINDING_LABEL_POLICY_VERSION = "spd_binding_interval_aware_v2"

SCORE_FEATURE_COLS = [
    "atlas_score",
    "consensus_score",
    "final_score",
    "final_score_source",
    "ml_blend_mode",
    "ml_blend_scorch_weight_effective",
    "ml_blend_cnn_weight_effective",
    "z_vs_decoys_blend",
    "z_vs_decoys_consensus",
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
RECEPTOR_MAPPING_COLS = [
    "spd_receptor_mapping_policy_version",
    "spd_receptor_mapping_status",
    "spd_receptor_mapping_strict_eligible",
    "spd_receptor_mapping_exclusion_reason",
    "spd_receptor_mapping_contract_source",
    "spd_receptor_mapping_contract_sha256",
    "spd_receptor_mapping_changed",
    "spd_receptor_mapping_requires_target_derived_rebuild",
    "spd_receptor_mapping_legacy_target_id",
    "spd_receptor_mapping_legacy_target_gene",
    "spd_receptor_mapping_legacy_target_uniprot",
    "spd_receptor_mapping_revised_target_id",
    "spd_receptor_mapping_revised_target_gene",
    "spd_receptor_mapping_revised_target_uniprot",
    "spd_receptor_mapping_selected_chain",
    "spd_receptor_mapping_site_chain",
]
RAW_KEEP_COLS = [
    *ID_COLS,
    *RECEPTOR_MAPPING_COLS,
    *SCORE_FEATURE_COLS,
    *EXPOSURE_CONTEXT_COLS,
    *TISSUE_CONTEXT_COLS,
    *MECHANISM_CONTEXT_COLS,
    "z_selected",
    "z_selected_source",
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
    *RECEPTOR_MAPPING_COLS,
    "source_objective",
    *MODEL_READY_PROVENANCE_COLS,
    *SCORE_FEATURE_COLS,
    *EXPOSURE_CONTEXT_COLS,
    *TISSUE_CONTEXT_COLS,
    *MECHANISM_CONTEXT_COLS,
]

STRICT_ELIGIBILITY_COL = "spd_receptor_mapping_strict_eligible"
STRICT_REBUILD_COL = "spd_receptor_mapping_requires_target_derived_rebuild"
STRICT_TARGET_DERIVED_CONTEXT_COLS = [
    "target_family",
    *TISSUE_CONTEXT_COLS,
    *MECHANISM_CONTEXT_COLS,
    "tissue_site_label",
    "tissue_relevance_label",
    "site_relevance_label",
    "target_site_relevance_label",
    "mechanism_ml_label",
    "mechanism_label",
    "four_state_ml_label",
    "drug_target_adr_mechanism_label",
]


def _numeric(series: pd.Series | None, index: pd.Index) -> pd.Series:
    if series is None:
        return pd.Series(pd.NA, index=index, dtype="Float64")
    return pd.to_numeric(series, errors="coerce").astype("Float64")


def _bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.fillna(False).astype(bool)
    normalized = series.fillna("").astype(str).str.strip().str.lower()
    return normalized.isin({"1", "true", "yes"})


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
    if "final_score" in out.columns:
        out["final_score"] = pd.to_numeric(out["final_score"], errors="coerce")
    if "final_score_source" not in out.columns and "final_score" in out.columns:
        out["final_score_source"] = out["final_score"].notna().map(
            {True: "legacy_missing_provenance", False: ""}
        )
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
    derived = pd.Series(pd.NA, index=source.index, dtype="Int64")
    interpretable_measurement = pd.Series(False, index=source.index, dtype=bool)
    for position, (value, raw_relation) in enumerate(zip(ac50, relation)):
        interval = parse_spd_activity_interval(value, raw_relation, "uM")
        if not interval.is_valid:
            continue
        interpretable_measurement.iloc[position] = True
        result = label_spd_binding_observation(
            interval,
            active_um=active_um,
            inactive_um=inactive_um,
        )
        if result.numeric_label is not None:
            derived.iloc[position] = result.numeric_label
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
    enforce_strict_eligibility: bool = False,
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
    eligible = pd.Series(True, index=table.index, dtype=bool)
    if enforce_strict_eligibility:
        if STRICT_ELIGIBILITY_COL not in table.columns:
            raise ValueError(f"Strict output is missing {STRICT_ELIGIBILITY_COL}")
        eligible = _bool_series(table[STRICT_ELIGIBILITY_COL])
        summary["strict_eligible_rows"] = int(eligible.sum())
        summary["strict_ineligible_rows"] = int((~eligible).sum())
        summary["strict_ineligible_labelable_rows_excluded"] = int(
            (label.notna() & ~eligible).sum()
        )
    supervised = table[label.notna() & eligible].copy()
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
    receptor_mapping_mode: str = "legacy",
    receptor_mapping_contract_path: str | Path | None = None,
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
        receptor_mapping_mode=receptor_mapping_mode,
        receptor_mapping_contract_path=receptor_mapping_contract_path,
    )
    rebuild_mask = pd.Series(False, index=source.index, dtype=bool)
    if receptor_mapping_mode == "strict":
        if STRICT_REBUILD_COL not in source.columns:
            raise ValueError(f"Strict receptor mapping did not produce {STRICT_REBUILD_COL}")
        rebuild_mask = _bool_series(source[STRICT_REBUILD_COL])
        for col in STRICT_TARGET_DERIVED_CONTEXT_COLS:
            if col in source.columns:
                source.loc[rebuild_mask, col] = pd.NA
    if rebuild_mask.any():
        unchanged = source.loc[~rebuild_mask].copy()
        to_rebuild = source.loc[rebuild_mask].copy()
        to_rebuild["_strict_rebuild_row_id"] = to_rebuild.index
        rebuilt, mechanism_projection_summary = project_mechanism_labels_for_run_master(
            to_rebuild,
            mechanism_label_path=mechanism_label_path,
            run_dir=run_dir,
        )
        if len(rebuilt) != len(to_rebuild) or rebuilt["_strict_rebuild_row_id"].duplicated().any():
            raise ValueError("Strict mechanism reprojection multiplied source rows")
        rebuilt = rebuilt.set_index("_strict_rebuild_row_id")
        source = pd.concat([unchanged, rebuilt], axis=0, sort=False).sort_index(kind="stable")
        mechanism_projection_summary = {
            **mechanism_projection_summary,
            "strict_target_identity_rebuild": True,
            "strict_rebuilt_rows": int(rebuild_mask.sum()),
        }
    else:
        source, mechanism_projection_summary = project_mechanism_labels_for_run_master(
            source,
            mechanism_label_path=mechanism_label_path,
            run_dir=run_dir,
        )
    rows_before_tissue = len(source)
    source, tissue_expression_summary = add_tissue_expression_context(
        source,
        target_expression_path=target_expression_path,
        run_dir=run_dir,
    )
    if len(source) != rows_before_tissue:
        raise ValueError("Tissue context projection multiplied source rows")
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
            f"SPD interval-aware activity: exact or upper-bound AC50 <= {active_um:g} uM "
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
        policy=(
            "SPD exposure relevance: the complete AC50/free-Cmax margin interval "
            "must be <= 10 for positive or strictly > 10 for negative; intervals "
            "that can equal 10 and extend above it remain unknown. Exact margins "
            "in (10, 100] remain supervised weak negatives. AC50/free-Cmax fields "
            "are label-definition fields, not predictors for nonleaky models."
        ),
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
        "receptor_mapping_mode": receptor_mapping_mode,
        "receptor_mapping_contract_path": (
            str(receptor_mapping_contract_path) if receptor_mapping_contract_path else None
        ),
        "strict_target_derived_rebuild_rows": int(rebuild_mask.sum()),
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
            enforce_strict_eligibility=receptor_mapping_mode == "strict",
        )
    manifest_path = out / "spd_four_expert_tables_manifest.json"
    manifest_path.write_text(json.dumps(summaries, indent=2, sort_keys=True), encoding="utf-8")
    return summaries
