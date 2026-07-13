from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


SCORE_COLS = [
    "run_id",
    "pdb_id",
    "ligand_base",
    "ligand_display",
    "atlas_score",
    "consensus_score",
    "SCORCH_score_used",
    "final_score",
    "z_selected",
    "z_selected_source",
    "pose_valid_any",
    "pose_invalid_reason_top",
    "pocket_method",
    "fdr_score_field",
    "fdr_null_source",
    "fdr_n_decoys",
    "fdr_reliable",
    "fdr_p_empirical",
    "fdr_q_target",
    "scorch_fdr_n_decoys",
    "scorch_fdr_n_tested",
    "scorch_fdr_q_bh",
    "scorch_fdr_decoy_competition_q_plus1",
    "rescored_flag",
    "rescored_stage",
    "source_csv",
]

REFERENCE_VINA_COLS = [
    "reference_vina_status",
    "reference_vina_stage1_kcal_mol",
    "reference_vina_empirical_p",
    "reference_vina_q_bh_global",
    "reference_vina_q_bh_target",
    "reference_vina_empirical_tail_z",
    "reference_vina_top10pct",
    "reference_vina_config",
    "reference_vina_pose",
    "z_vs_compare_run_vina_stage1",
    "z_vs_compare_run_vina_stage1_source",
    "z_vs_compare_run_vina_stage1_run_id",
    "z_vs_compare_run_vina_stage1_run_sha256",
    "z_vs_compare_run_vina_stage1_null_sha256",
    "z_vs_compare_run_vina_stage1_decoy_n",
    "z_vs_compare_run_vina_stage1_decoy_mu",
    "z_vs_compare_run_vina_stage1_decoy_sigma",
    "z_vs_compare_run_vina_stage1_decoy_median",
    "z_vs_compare_run_vina_stage1_decoy_mad",
    "z_vs_compare_run_vina_stage1_decoy_unique_scores",
    "z_vs_compare_run_vina_stage1_score_space_sha256",
    "z_vs_compare_run_vina_stage1_receptor_sha256",
    "z_vs_compare_run_vina_stage1_ligand_pdbqt_sha256",
]

ADDON_SCORE_PROVENANCE_COLS = [
    "addon_run_consensus_rank",
    "addon_run_consensus_rank_source",
    "addon_run_target_row_count",
    "addon_run_score_comparability",
    "addon_run_scorch_score",
    *REFERENCE_VINA_COLS,
]

ADDON_PROVENANCE_COLS = [
    "_internal_spd_addons_score_pair_key",
    "_internal_spd_addons_run_mode",
    "_internal_spd_addons_library_original",
    "scenario_a_tier1_label",
    "scenario_a_tier1_label_status",
    "scenario_a_tier1_sources",
    "scenario_a_tier1_raw_sources",
    "scenario_a_tier1_min_activity_nM",
    "scenario_a_tier1_free_cmax_um",
    "scenario_a_tier1_external_exposure_margin",
    "scenario_a_tier1_external_exposure_relevant_candidate",
    "external_addon_source_drug_id",
    "external_addon_activity_relation",
    "external_addon_activity_nM",
    "external_addon_activity_uM",
    "external_addon_assay_id",
    "external_addon_activity_type",
    "external_addon_publication_year",
    "external_addon_ligand_domain",
    "external_addon_training_domain",
    "external_addon_selection_reason",
    "external_addon_source_priority",
    "external_addon_endpoint_priority",
    "external_addon_selection_tier",
    "external_addon_present_ligand_in_phase1",
    "external_addon_present_pair_in_phase1",
    "external_addon_compound_key",
    "external_addon_benchmark_only",
    "external_addon_production_truth_allowed",
    "external_addon_training_allowed",
    "external_addon_n_raw_assay_rows",
    "external_addon_n_direct_potency_rows",
    "external_addon_n_hts_rows",
    "external_addon_n_curated_positive_rows",
    "external_addon_label_policy",
    "pair_evidence_record_count",
    "pair_evidence_collapsed",
    "pair_evidence_sources_collapsed",
    "pair_evidence_assay_ids_collapsed",
    "pair_evidence_compound_keys_collapsed",
    *ADDON_SCORE_PROVENANCE_COLS,
]


def _binary_label(selected: pd.DataFrame) -> pd.Series:
    if "projected_label" in selected.columns:
        projected = pd.to_numeric(selected["projected_label"], errors="coerce")
        label = pd.Series(pd.NA, index=selected.index, dtype="Int64")
        label = label.mask(projected.eq(1).fillna(False), 1)
        label = label.mask(projected.eq(0).fillna(False), 0)
        if label.notna().any():
            return label.astype("Int64")

    activity = pd.to_numeric(selected.get("activity_uM"), errors="coerce")
    label = pd.Series(pd.NA, index=selected.index, dtype="Int64")
    label = label.mask(activity.le(1).fillna(False), 1)
    label = label.mask(activity.ge(10).fillna(False), 0)
    if "is_binding_positive_le1uM" in selected.columns:
        source_positive = selected["is_binding_positive_le1uM"].map({True: 1, "True": 1})
        label = label.where(label.notna(), source_positive)
    return label.astype("Int64")


def _ligand_key(frame: pd.DataFrame) -> pd.Series:
    ligand_base = frame.get("ligand_base", pd.Series("", index=frame.index)).fillna("").astype(str)
    if "ligand_file_stem" not in frame.columns:
        return ligand_base.str.casefold()
    file_stem = frame["ligand_file_stem"].fillna("").astype(str)
    return file_stem.where(file_stem.str.strip().ne(""), ligand_base).str.casefold()


def _row_key(frame: pd.DataFrame) -> pd.Series:
    return frame["pdb_id"].astype(str).str.upper() + "|" + _ligand_key(frame)


def _numeric_column(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(float("nan"), index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _first_nonempty(frame: pd.DataFrame, columns: tuple[str, ...]) -> pd.Series:
    """Return the first populated text value from each row."""

    result = pd.Series("", index=frame.index, dtype="object")
    for column in columns:
        if column not in frame.columns:
            continue
        candidate = frame[column].fillna("").astype(str).str.strip()
        result = result.where(result.str.strip().ne(""), candidate)
    return result.replace("", pd.NA)


def _collapse_duplicate_score_pairs(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, int]]:
    work = frame.copy()
    work["_score_pair_key"] = _row_key(work)
    work["_pair_original_order"] = range(len(work))
    duplicate_mask = work.duplicated("_score_pair_key", keep=False)
    duplicate_groups = int(work.loc[duplicate_mask, "_score_pair_key"].nunique())
    unique = work.loc[~duplicate_mask].copy()
    unique["pair_evidence_record_count"] = 1
    unique["pair_evidence_collapsed"] = False
    label_columns = [
        column
        for column in (
            "spd_binding_label",
            "spd_exposure_label",
            "combined_activity_ml_label",
            "external_four_state_ml_label",
        )
        if column in work.columns
    ]
    conflict_groups = 0
    rows: list[pd.Series] = []
    duplicates = work.loc[duplicate_mask]
    for _, group in duplicates.groupby(
        "_score_pair_key",
        sort=False,
        dropna=False,
    ):
        for column in label_columns:
            labels = pd.to_numeric(group[column], errors="coerce").dropna().unique()
            if len(labels) > 1:
                conflict_groups += 1
                break
        priority = pd.to_numeric(
            group.get(
                "external_addon_source_priority",
                pd.Series(float("nan"), index=group.index),
            ),
            errors="coerce",
        ).fillna(float("inf"))
        activity = pd.to_numeric(
            group.get(
                "external_addon_activity_nM",
                pd.Series(float("nan"), index=group.index),
            ),
            errors="coerce",
        ).fillna(float("inf"))
        representative = group.loc[
            pd.DataFrame(
                {"priority": priority, "activity": activity},
                index=group.index,
            )
            .sort_values(["priority", "activity"])
            .index[0]
        ].copy()
        representative["_pair_original_order"] = int(
            group["_pair_original_order"].min()
        )
        representative["pair_evidence_record_count"] = int(len(group))
        representative["pair_evidence_collapsed"] = True
        for destination, sources in {
            "pair_evidence_sources_collapsed": (
                "scenario_a_tier1_sources",
                "external_activity_source_family",
                "combined_activity_source_family",
            ),
            "pair_evidence_assay_ids_collapsed": (
                "external_addon_assay_id",
                "bindingdb_assay_ids",
            ),
            "pair_evidence_compound_keys_collapsed": (
                "external_addon_compound_key",
                "inchikey",
            ),
        }.items():
            values: set[str] = set()
            for source in sources:
                if source not in group.columns:
                    continue
                for value in group[source].dropna().astype(str):
                    values.update(
                        item.strip()
                        for item in value.split(";")
                        if item.strip()
                    )
            representative[destination] = ";".join(sorted(values))
        rows.append(representative)
    if conflict_groups:
        raise ValueError(
            f"{conflict_groups} duplicate score-pair group(s) have conflicting labels"
        )
    duplicate_representatives = pd.DataFrame(rows, columns=work.columns)
    collapsed = (
        pd.concat([unique, duplicate_representatives], ignore_index=True, sort=False)
        .sort_values("_pair_original_order")
        .drop(columns=["_score_pair_key", "_pair_original_order"])
        .reset_index(drop=True)
    )
    return collapsed, {
        "input_rows": int(len(work)),
        "output_rows": int(len(collapsed)),
        "duplicate_groups_collapsed": duplicate_groups,
        "rows_removed": int(len(work) - len(collapsed)),
        "conflict_groups": conflict_groups,
    }


def build_addon_rows(base: pd.DataFrame, master: pd.DataFrame, selected: pd.DataFrame, *, addon_name: str) -> pd.DataFrame:
    base = base.copy()
    for col in ADDON_PROVENANCE_COLS:
        if col not in base.columns:
            base[col] = pd.NA
    base_cols = list(base.columns)
    selected = selected.copy()
    master = master.copy()
    selected["_addon_pair_key"] = _row_key(selected)
    master["_addon_pair_key"] = _row_key(master)

    score_cols = [
        col for col in (*SCORE_COLS, *REFERENCE_VINA_COLS) if col in master.columns
    ]
    score = master[["_addon_pair_key", *score_cols]].drop_duplicates("_addon_pair_key", keep="first")
    joined = selected.merge(score, on="_addon_pair_key", how="left", suffixes=("", "_score"))

    rows = pd.DataFrame(index=joined.index, columns=base_cols)
    for col in rows.columns:
        if col in joined.columns:
            rows[col] = joined[col]

    rows["pdb_id"] = joined["pdb_id"].astype(str).str.upper()
    rows["ligand_base"] = _ligand_key(joined)
    rows["drug_id"] = joined.get("generic_name", joined.get("display_name", joined["ligand_base"]))
    rows["target_id"] = joined.get("target_uniprot", joined.get("target_gene", joined["pdb_id"]))
    rows["target_gene"] = joined.get("target_gene")
    rows["target_uniprot"] = joined.get("target_uniprot")
    rows["display_name"] = joined.get("display_name")
    rows["generic_name"] = joined.get("generic_name")
    rows["target_family"] = joined.get("target_family")
    rows["protein_class"] = joined.get("target_family")

    for source_col, dest_col in [
        ("candidate_chemotype", "ligand_chemotype"),
        ("candidate_scaffold_key", "scaffold_key"),
    ]:
        if source_col in joined.columns and dest_col in rows.columns:
            rows[dest_col] = joined[source_col]

    consensus_rank = _numeric_column(joined, "consensus_score")
    reference_z = _numeric_column(joined, "z_vs_compare_run_vina_stage1")
    rows["addon_run_consensus_rank"] = consensus_rank
    rows["addon_run_consensus_rank_source"] = consensus_rank.map(
        lambda value: (
            "within_addon_target_library_rank_percentile"
            if pd.notna(value)
            else pd.NA
        )
    )
    master_target_counts = (
        master["pdb_id"].fillna("").astype(str).str.upper().value_counts()
    )
    rows["addon_run_target_row_count"] = (
        joined["pdb_id"].fillna("").astype(str).str.upper().map(master_target_counts)
    )
    rows["addon_run_score_comparability"] = reference_z.map(
        lambda value: (
            "common_spd90_raw_vina_stage1_dud_score_space"
            if pd.notna(value)
            else "not_comparable_to_full_spd_dud_library"
        )
    )
    rows["addon_run_scorch_score"] = pd.to_numeric(
        joined.get("SCORCH_score_used"), errors="coerce"
    )

    for col in [
        "consensus_z_score",
        "consensus_z_score_source",
        "consensus_z_decoy_n",
        "consensus_z_decoy_unique",
        "consensus_z_decoy_zero_fraction",
        "consensus_z_decoy_mu",
        "consensus_z_decoy_sigma",
        "final_score",
        "z_selected",
        "z_selected_source",
        "pose_valid_any",
        "pose_invalid_reason_top",
        "pocket_method",
        "fdr_n_decoys",
        "fdr_reliable",
        "fdr_q_target",
        "scorch_fdr_n_decoys",
        "scorch_fdr_n_tested",
        "rescored_flag",
        "rescored_stage",
    ]:
        if col in joined.columns and col in rows.columns:
            rows[col] = joined[col]

    # Percentiles from a tiny add-on library are not interchangeable with the
    # full SPD/FDA/DUD score space. Preserve them above for audit, but do not
    # expose them as primary model features until a common-background rescore.
    for col in (
        "consensus_score",
        "consensus_score_raw",
        "SCORCH_score_used",
    ):
        if col in rows.columns:
            rows[col] = pd.NA

    if "atlas_score" in rows.columns:
        z_feature = pd.to_numeric(
            rows.get("consensus_z_score", pd.Series(index=rows.index)),
            errors="coerce",
        )
        selected = pd.to_numeric(
            rows.get("z_selected", pd.Series(index=rows.index)),
            errors="coerce",
        )
        selected_source = rows.get(
            "z_selected_source",
            pd.Series("", index=rows.index),
        ).fillna("").astype(str).str.lower()
        valid_selected = ~selected_source.str.contains(
            "fallback|missing_consensus_decoy_null",
            regex=True,
        )
        legacy_atlas = z_feature.where(
            z_feature.notna(), selected.where(valid_selected)
        )
        rows["atlas_score"] = reference_z.where(reference_z.notna(), legacy_atlas)
    if "final_score" in rows.columns:
        legacy_final = pd.to_numeric(rows["final_score"], errors="coerce")
        rows["final_score"] = reference_z.where(reference_z.notna(), legacy_final)
    if "final_score_source" in rows.columns:
        legacy_source = rows["final_score_source"]
        rows["final_score_source"] = legacy_source.where(
            reference_z.isna(), "comparison_run_vina_stage1"
        )

    label = _binary_label(joined)
    for col in ["scenario_a_tier1_label", "external_four_state_ml_label", "external_four_state_label", "combined_activity_label", "combined_activity_ml_label"]:
        if col in rows.columns:
            rows[col] = label

    # External direct assays have different semantics from the uniform SPD
    # AC50 panel. They may train the combined activity sensitivity model, but
    # they must never be represented as SPD assay or exposure truth.
    for col in (
        "spd_binding_label",
        "spd_binding_ml_label",
        "spd_activity_label",
        "spd_activity_ml_label",
        "spd_exposure_label",
        "spd_exposure_ml_label",
        "spd_exposure_relevant",
        "spd_exposure_weak",
        "spd_exposure_unlikely",
        "spd_ac50",
        "spd_ac50_nm",
        "spd_ac50_uM",
        "spd_exposure_margin",
        "spd_activity_relation",
        "spd_inchikey",
    ):
        if col in rows.columns:
            rows[col] = pd.NA

    if "scenario_a_tier1_label_status" in rows.columns:
        rows["scenario_a_tier1_label_status"] = label.map({1: "projected_positive", 0: "projected_negative"}).fillna("projected_unknown")
    if "external_four_state_label_status" in rows.columns:
        rows["external_four_state_label_status"] = rows.get("scenario_a_tier1_label_status")
    if "combined_activity_label_status" in rows.columns:
        rows["combined_activity_label_status"] = rows.get("scenario_a_tier1_label_status")
    if "scenario_a_tier1_sources" in rows.columns:
        rows["scenario_a_tier1_sources"] = joined.get("sources")
    if "scenario_a_tier1_raw_sources" in rows.columns:
        rows["scenario_a_tier1_raw_sources"] = joined.get("raw_sources")
    if "scenario_a_tier1_min_activity_nM" in rows.columns:
        rows["scenario_a_tier1_min_activity_nM"] = pd.to_numeric(joined.get("activity_uM"), errors="coerce") * 1000.0
    if "scenario_a_tier1_free_cmax_um" in rows.columns:
        rows["scenario_a_tier1_free_cmax_um"] = pd.to_numeric(joined.get("free_cmax_uM_for_margin"), errors="coerce")
    if "scenario_a_tier1_external_exposure_margin" in rows.columns:
        rows["scenario_a_tier1_external_exposure_margin"] = pd.to_numeric(joined.get("exposure_margin"), errors="coerce")
    if "scenario_a_tier1_external_exposure_relevant_candidate" in rows.columns:
        rows["scenario_a_tier1_external_exposure_relevant_candidate"] = joined.get("is_exposure_relevant_margin_le10")

    selected_provenance = {
        "external_addon_source_drug_id": "drug_id",
        "external_addon_activity_relation": "activity_relation",
        "external_addon_activity_nM": "activity_nM",
        "external_addon_activity_uM": "activity_uM",
        "external_addon_assay_id": "assay_id",
        "external_addon_activity_type": "activity_type",
        "external_addon_publication_year": "publication_year",
        "external_addon_ligand_domain": "ligand_domain",
        "external_addon_training_domain": "training_domain",
        "external_addon_selection_reason": "selection_reason",
        "external_addon_source_priority": "source_priority",
        "external_addon_endpoint_priority": "endpoint_priority",
        "external_addon_selection_tier": "selection_tier",
        "external_addon_present_ligand_in_phase1": "present_ligand_in_phase1",
        "external_addon_present_pair_in_phase1": "present_pair_in_phase1",
        "external_addon_compound_key": "compound_key",
        "external_addon_benchmark_only": "benchmark_only",
        "external_addon_production_truth_allowed": "production_truth_allowed",
        "external_addon_training_allowed": "training_allowed",
        "external_addon_n_raw_assay_rows": "n_raw_assay_rows",
        "external_addon_n_direct_potency_rows": "n_direct_potency_rows",
        "external_addon_n_hts_rows": "n_hts_rows",
        "external_addon_n_curated_positive_rows": "n_curated_positive_rows",
        "external_addon_label_policy": "addon_label_policy",
    }
    for dest_col, source_col in selected_provenance.items():
        if dest_col in rows.columns and source_col in joined.columns:
            rows[dest_col] = joined[source_col]

    provenance = {
        "source_objective": "external_binding_activity",
        "label_source": f"{addon_name}_external",
        "source_family": f"{addon_name}_external",
        "upstream_source": "BindingDB;ChEMBL;DrugCentral;ToxCast;SPD",
        "source_label_policy": "external_direct_binding_positive_or_strict_spd_exposure_candidate",
        "assay_type": "external_direct_binding_addon",
        "endpoint_type": "activity",
        "activity_type": "activity",
        "benchmark_only": False,
        "training_allowed": True,
    }
    for col, value in provenance.items():
        if col in rows.columns:
            current = rows[col]
            if current.notna().any():
                rows[col] = current.fillna(value)
            else:
                rows[col] = value

    row_source = _first_nonempty(
        joined,
        ("source_name", "source_family", "representative_source", "source", "sources"),
    )
    row_policy = _first_nonempty(
        joined,
        ("source_label_policy", "addon_label_policy", "measured_positive_policy"),
    )
    for column in ("source_family", "upstream_source"):
        if column in rows.columns:
            rows[column] = row_source.where(row_source.notna(), rows[column])
    if "label_source" in rows.columns:
        rows["label_source"] = row_source.map(
            lambda value: f"{value}_{addon_name}" if pd.notna(value) else pd.NA
        ).where(row_source.notna(), rows["label_source"])
    if "source_label_policy" in rows.columns:
        rows["source_label_policy"] = row_policy.where(
            row_policy.notna(), rows["source_label_policy"]
        )

    if "_internal_spd_addons_score_pair_key" in rows.columns:
        rows["_internal_spd_addons_score_pair_key"] = joined["_addon_pair_key"]
    if "_internal_spd_addons_run_mode" in rows.columns:
        rows["_internal_spd_addons_run_mode"] = master.get("run_mode", pd.Series(["dud"])).dropna().astype(str).iloc[0] if "run_mode" in master else "dud"
    if "_internal_spd_addons_library_original" in rows.columns:
        rows["_internal_spd_addons_library_original"] = "dud"

    if "structure_quality" in rows.columns and "structure_quality" in base.columns:
        quality = base[["pdb_id", "structure_quality"]].dropna().drop_duplicates("pdb_id")
        quality["pdb_id"] = quality["pdb_id"].astype(str).str.upper()
        rows = rows.drop(columns=["structure_quality"]).merge(quality, on="pdb_id", how="left")
        rows = rows.reindex(columns=base_cols)

    rows = rows.loc[:, base_cols]
    return rows


def _backfill_existing_reference_scores(
    base: pd.DataFrame,
    addon: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Fill missing score fields on existing pairs without touching labels/identity."""

    if addon.empty:
        return base, {"matched_rows": 0, "rows_updated": 0, "values_filled": 0}
    out = base.copy()
    out["_reference_backfill_key"] = _row_key(out)
    source = addon.copy()
    source["_reference_backfill_key"] = _row_key(source)
    if source["_reference_backfill_key"].duplicated().any():
        raise ValueError("reference-score backfill contains duplicate pair keys")
    lookup = source.set_index("_reference_backfill_key")
    matched = out["_reference_backfill_key"].isin(lookup.index)
    rows_updated = pd.Series(False, index=out.index)
    values_filled = 0
    fill_columns = [
        *REFERENCE_VINA_COLS,
        "addon_run_score_comparability",
        "final_score",
        "final_score_source",
        "atlas_score",
    ]
    for column in fill_columns:
        if column not in out or column not in lookup:
            continue
        candidate = out["_reference_backfill_key"].map(lookup[column])
        missing = out[column].isna()
        if out[column].dtype == object:
            missing |= out[column].fillna("").astype(str).str.strip().eq("")
        fill_mask = matched & missing & candidate.notna()
        if not fill_mask.any():
            continue
        out.loc[fill_mask, column] = candidate.loc[fill_mask]
        rows_updated |= fill_mask
        values_filled += int(fill_mask.sum())
    out = out.drop(columns=["_reference_backfill_key"])
    return out, {
        "matched_rows": int(matched.sum()),
        "rows_updated": int(rows_updated.sum()),
        "values_filled": values_filled,
    }


def merge_spd_addon_tables(
    *,
    base_table: Path,
    master_rows: Path,
    selected_pairs: Path,
    out: Path,
    addon_name: str = "spdaddon",
    addon_rows_out: Path | None = None,
    strict_selected_coverage: bool = False,
    backfill_existing_scores: bool = False,
) -> dict[str, object]:
    """Merge one scored add-on selection while preserving assay namespaces."""

    base = pd.read_csv(base_table, low_memory=False)
    for col in ADDON_PROVENANCE_COLS:
        if col not in base.columns:
            base[col] = pd.NA
    base, base_collapse = _collapse_duplicate_score_pairs(base)
    master = pd.read_csv(master_rows, low_memory=False)
    selected = pd.read_csv(selected_pairs, low_memory=False)
    selected_keys = set(_row_key(selected).astype(str))
    master_keys = set(_row_key(master).astype(str))
    missing_selected_keys = sorted(selected_keys - master_keys)
    if strict_selected_coverage and missing_selected_keys:
        preview = ", ".join(missing_selected_keys[:10])
        raise ValueError(
            f"{len(missing_selected_keys)} selected pair(s) have no scored master row; examples: {preview}"
        )
    addon = build_addon_rows(base, master, selected, addon_name=addon_name)
    scored = pd.Series(False, index=addon.index)
    for column in (
        "z_vs_compare_run_vina_stage1",
        "addon_run_consensus_rank",
        "addon_run_scorch_score",
    ):
        if column in addon:
            scored |= pd.to_numeric(addon[column], errors="coerce").notna()
    addon = addon.loc[scored].copy()

    base_keys = set(_row_key(base).astype(str)) if {"pdb_id", "ligand_base"}.issubset(base.columns) else set()
    addon_keys = _row_key(addon).astype(str)
    existing_addon = addon.loc[addon_keys.isin(base_keys)].copy()
    backfill = {"matched_rows": 0, "rows_updated": 0, "values_filled": 0}
    if backfill_existing_scores:
        base, backfill = _backfill_existing_reference_scores(base, existing_addon)
    addon = addon.loc[~addon_keys.isin(base_keys)].copy()

    merged = pd.concat([base, addon], ignore_index=True, sort=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, index=False)
    if addon_rows_out:
        addon_rows_out.parent.mkdir(parents=True, exist_ok=True)
        addon.to_csv(addon_rows_out, index=False)

    manifest: dict[str, object] = {
        "base_table": str(base_table),
        "master_rows": str(master_rows),
        "selected_pairs": str(selected_pairs),
        "out": str(out),
        "addon_rows_out": str(addon_rows_out) if addon_rows_out else None,
        "base_rows": int(len(base)),
        "base_pair_evidence_collapse": base_collapse,
        "selected_rows": int(len(selected)),
        "selected_unique_pairs": int(len(selected_keys)),
        "selected_pairs_without_score": int(len(missing_selected_keys)),
        "selected_pairs_without_score_examples": missing_selected_keys[:20],
        "scored_selected_rows": int(len(addon)),
        "merged_rows": int(len(merged)),
        "addon_positive_rows": int((pd.to_numeric(addon.get("combined_activity_ml_label"), errors="coerce") == 1).sum()),
        "addon_negative_rows": int((pd.to_numeric(addon.get("combined_activity_ml_label"), errors="coerce") == 0).sum()),
        "addon_common_reference_vina_rows": int(
            pd.to_numeric(
                addon.get("z_vs_compare_run_vina_stage1"), errors="coerce"
            ).notna().sum()
        ),
        "existing_pair_score_backfill": backfill,
    }
    out.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Merge a targeted SPD add-on docking run into a Phase 1 model-ready table.")
    parser.add_argument("--base-table", required=True, type=Path)
    parser.add_argument("--master-rows", required=True, type=Path)
    parser.add_argument("--selected-pairs", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--addon-name", default="spdaddon")
    parser.add_argument("--addon-rows-out", type=Path, default=None)
    parser.add_argument(
        "--strict-selected-coverage",
        action="store_true",
        help="Fail when a selected PDB-ligand pair has no corresponding scored master row.",
    )
    parser.add_argument(
        "--backfill-existing-scores",
        action="store_true",
        help=(
            "Fill only missing common-reference score fields for selected pairs "
            "already present in the base table; labels and identities are untouched."
        ),
    )
    args = parser.parse_args(argv)
    manifest = merge_spd_addon_tables(
        base_table=args.base_table,
        master_rows=args.master_rows,
        selected_pairs=args.selected_pairs,
        out=args.out,
        addon_name=args.addon_name,
        addon_rows_out=args.addon_rows_out,
        strict_selected_coverage=args.strict_selected_coverage,
        backfill_existing_scores=args.backfill_existing_scores,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
