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

ADDON_SCORE_PROVENANCE_COLS = [
    "addon_run_consensus_rank",
    "addon_run_consensus_rank_source",
    "addon_run_target_row_count",
    "addon_run_score_comparability",
    "addon_run_scorch_score",
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

    score_cols = [col for col in SCORE_COLS if col in master.columns]
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

    rows["addon_run_consensus_rank"] = pd.to_numeric(
        joined.get("consensus_score"), errors="coerce"
    )
    rows["addon_run_consensus_rank_source"] = (
        "within_addon_target_library_rank_percentile"
    )
    master_target_counts = (
        master["pdb_id"].fillna("").astype(str).str.upper().value_counts()
    )
    rows["addon_run_target_row_count"] = (
        joined["pdb_id"].fillna("").astype(str).str.upper().map(master_target_counts)
    )
    rows["addon_run_score_comparability"] = (
        "not_comparable_to_full_spd_dud_library"
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
        rows["atlas_score"] = z_feature.where(z_feature.notna(), selected.where(valid_selected))

    label = _binary_label(joined)
    for col in ["scenario_a_tier1_label", "external_four_state_ml_label", "external_four_state_label", "combined_activity_label", "combined_activity_ml_label"]:
        if col in rows.columns:
            rows[col] = label

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
        "source_objective": "spd_binding_activity",
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
    args = parser.parse_args(argv)

    base = pd.read_csv(args.base_table, low_memory=False)
    for col in ADDON_PROVENANCE_COLS:
        if col not in base.columns:
            base[col] = pd.NA
    base, base_collapse = _collapse_duplicate_score_pairs(base)
    master = pd.read_csv(args.master_rows, low_memory=False)
    selected = pd.read_csv(args.selected_pairs, low_memory=False)
    selected_keys = set(_row_key(selected).astype(str))
    master_keys = set(_row_key(master).astype(str))
    missing_selected_keys = sorted(selected_keys - master_keys)
    if args.strict_selected_coverage and missing_selected_keys:
        preview = ", ".join(missing_selected_keys[:10])
        raise ValueError(
            f"{len(missing_selected_keys)} selected pair(s) have no scored master row; examples: {preview}"
        )
    addon = build_addon_rows(base, master, selected, addon_name=args.addon_name)
    addon = addon.dropna(subset=["addon_run_consensus_rank"], how="all")

    base_keys = set(_row_key(base).astype(str)) if {"pdb_id", "ligand_base"}.issubset(base.columns) else set()
    addon_keys = _row_key(addon).astype(str)
    addon = addon.loc[~addon_keys.isin(base_keys)].copy()

    merged = pd.concat([base, addon], ignore_index=True, sort=False)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, index=False)
    if args.addon_rows_out:
        args.addon_rows_out.parent.mkdir(parents=True, exist_ok=True)
        addon.to_csv(args.addon_rows_out, index=False)

    manifest = {
        "base_table": str(args.base_table),
        "master_rows": str(args.master_rows),
        "selected_pairs": str(args.selected_pairs),
        "out": str(args.out),
        "addon_rows_out": str(args.addon_rows_out) if args.addon_rows_out else None,
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
    }
    args.out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
