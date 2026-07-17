from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.feature_metadata import enrich_ml_feature_metadata
from analysis.ml.tissue_site_labels import add_independent_tissue_site_labels


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _key(value: Any) -> str:
    return _clean(value).lower()


def _source_family_from_evidence(value: Any) -> str:
    text = _key(value)
    families = [
        label
        for token, label in (
            ("spd", "SPD"),
            ("toxcast", "ToxCast"),
            ("chembl", "ChEMBL"),
            ("bindingdb", "BindingDB"),
            ("iuphar", "IUPHAR"),
            ("guide to pharmacology", "IUPHAR"),
            ("ttd", "TTD"),
            ("pubchem", "PubChem"),
            ("papyrus", "Papyrus"),
            ("drugcentral", "DrugCentral"),
        )
        if token in text
    ]
    return ";".join(dict.fromkeys(families)) or "external_unspecified"


def _first_nonempty(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in cols:
        if col not in df.columns:
            continue
        values = df[col].map(_key)
        out = out.where(out.astype(str).str.len().gt(0), values)
    return out


def _make_merge_keys(df: pd.DataFrame, *, table_role: str) -> pd.DataFrame:
    out = df.copy()
    if table_role == "spd":
        drug_name = _first_nonempty(out, ["generic_name", "display_name", "mapped_drug_name", "drug_id"])
        target_gene = _first_nonempty(out, ["target_gene", "gene_symbol", "target_id"])
        target_uniprot = _first_nonempty(out, ["target_uniprot", "uniprot", "target_id"])
        smiles = _first_nonempty(out, ["canonical_smiles", "smiles", "ligand_smiles"])
        inchikey = _first_nonempty(out, ["inchikey", "standard_inchikey", "ligand_inchikey"])
        rdk = _first_nonempty(out, ["ligand_base", "rdk_id", "drug_id"])
    else:
        drug_name = _first_nonempty(out, ["spd_drug_name_id", "generic_name", "display_name", "drug_name", "drug_id"])
        target_gene = _first_nonempty(out, ["spd_target_gene_id", "target_gene", "gene_symbol"])
        target_uniprot = _first_nonempty(out, ["target_uniprot", "uniprot", "target_id"])
        smiles = _first_nonempty(out, ["canonical_smiles", "smiles", "ligand_smiles"])
        inchikey = _first_nonempty(out, ["inchikey", "standard_inchikey", "ligand_inchikey"])
        rdk = _first_nonempty(out, ["drug_id", "ligand_base", "rdk_id"])
    out["_join_drug_name"] = drug_name
    out["_join_target_gene"] = target_gene
    out["_join_target_uniprot"] = target_uniprot
    out["_join_smiles"] = smiles
    out["_join_inchikey"] = inchikey
    out["_join_rdk"] = rdk
    out["_key_name_gene"] = drug_name + "||" + target_gene
    out["_key_name_uniprot"] = drug_name + "||" + target_uniprot
    out["_key_smiles_gene"] = smiles + "||" + target_gene
    out["_key_smiles_uniprot"] = smiles + "||" + target_uniprot
    out["_key_inchikey_gene"] = inchikey + "||" + target_gene
    out["_key_inchikey_uniprot"] = inchikey + "||" + target_uniprot
    out["_key_rdk_gene"] = rdk + "||" + target_gene
    out["_key_rdk_uniprot"] = rdk + "||" + target_uniprot
    return out


def _valid_key(series: pd.Series) -> pd.Series:
    values = series.fillna("").astype(str)
    has_separator = values.str.contains("||", regex=False)
    parts = values.str.split("||", n=1, expand=True, regex=False)
    if parts.shape[1] < 2:
        return pd.Series(False, index=series.index)
    return has_separator & parts[0].str.len().gt(0) & parts[1].str.len().gt(0)


def _label_priority(row: pd.Series) -> tuple[int, str]:
    status = _key(row.get("four_state_label_status"))
    state = pd.to_numeric(pd.Series([row.get("four_state_ml_label")]), errors="coerce").iloc[0]
    if pd.notna(state) and int(state) == 1:
        return 0, "external_positive"
    if pd.notna(state) and int(state) == 0:
        return 1, "external_negative"
    if "conflict" in status or "excluded" in status:
        return 2, "external_excluded_or_conflicting"
    return 3, "external_unknown"


def _collapse_external(additive: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    keep_cols = [
        "four_state_ml_label",
        "four_state_label",
        "four_state_label_status",
        "evidence_sources",
        "parent_sources",
        "n_raw_evidence_rows",
        "n_production_evidence_rows",
        "n_benchmark_only_rows",
        "evidence_state_counts_json",
        "max_confidence",
        "source_family",
        "upstream_source",
        "drug_id",
        "target_id",
        "target_gene",
        "target_uniprot",
        "spd_drug_name_id",
        "spd_target_gene_id",
        "display_name",
        "generic_name",
        "inchikey",
        "smiles",
    ]
    work = additive[[col for col in keep_cols if col in additive.columns] + [col for col in additive.columns if col.startswith("_key_")]].copy()
    parsed = pd.to_numeric(work.get("four_state_ml_label"), errors="coerce")
    work["four_state_ml_label"] = parsed.where(parsed.isin([0, 1]))
    priorities = work.apply(_label_priority, axis=1, result_type="expand")
    work["_external_priority"] = priorities[0]
    work["external_label_basis"] = priorities[1]
    key_cols = [col for col in work.columns if col.startswith("_key_")]
    rows = []
    for key_col in key_cols:
        subset = work.loc[_valid_key(work[key_col]), [key_col, *[c for c in work.columns if not c.startswith("_key_")]]].copy()
        if subset.empty:
            continue
        subset = subset.sort_values([key_col, "_external_priority", "max_confidence"], ascending=[True, True, False], na_position="last")
        dedup = subset.drop_duplicates(key_col, keep="first").rename(columns={key_col: "_merge_key"})
        dedup["external_join_key_type"] = key_col.removeprefix("_key_")
        rows.append(dedup)
    collapsed = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if collapsed.empty:
        return collapsed, {"external_rows": int(len(additive)), "collapsed_rows": 0}
    # Prefix non-key source columns so they cannot collide with SPD labels/features.
    rename = {col: f"external_{col}" for col in collapsed.columns if col not in {"_merge_key", "external_join_key_type"}}
    collapsed = collapsed.rename(columns=rename)
    return collapsed, {
        "external_rows": int(len(additive)),
        "collapsed_key_rows": int(len(collapsed)),
        "join_key_types": sorted(collapsed["external_join_key_type"].dropna().unique().tolist()),
    }


def _join_external(spd: pd.DataFrame, external: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    out = spd.copy()
    if external.empty:
        out["external_four_state_ml_label"] = pd.NA
        return out, {"status": "no_external_rows"}
    key_priority = [
        "_key_rdk_uniprot",
        "_key_rdk_gene",
        "_key_inchikey_uniprot",
        "_key_inchikey_gene",
        "_key_smiles_uniprot",
        "_key_smiles_gene",
        "_key_name_uniprot",
        "_key_name_gene",
    ]
    ext_cols = [col for col in external.columns if col not in {"_merge_key", "external_join_key_type"}]
    explicit_addon = out.get(
        "external_addon_training_allowed",
        pd.Series(False, index=out.index),
    ).fillna(False).astype(str).str.lower().isin({"1", "true", "yes"})
    earlier_addon_label = pd.to_numeric(
        out.get("scenario_a_tier1_label", pd.Series(pd.NA, index=out.index)), errors="coerce"
    )
    addon_mask = explicit_addon | earlier_addon_label.isin([0, 1])
    for col in ext_cols:
        if col not in out.columns:
            out[col] = pd.NA
        else:
            out.loc[~addon_mask, col] = pd.NA
    for col in ("external_join_key_type", "external_join_key_value", "external_join_confidence"):
        if col not in out.columns:
            out[col] = pd.NA
        else:
            out.loc[~addon_mask, col] = pd.NA
    existing_label = pd.to_numeric(out.get("external_four_state_ml_label"), errors="coerce")
    filled = addon_mask & existing_label.isin([0, 1])
    preserved_addon_rows = int(filled.sum())
    join_counts: dict[str, int] = {}
    confidence = {
        "rdk_uniprot": "exact_identifier_target",
        "rdk_gene": "exact_identifier_target",
        "inchikey_uniprot": "exact_structure_target",
        "inchikey_gene": "exact_structure_target",
        "smiles_uniprot": "exact_structure_target",
        "smiles_gene": "exact_structure_target",
        "name_uniprot": "name_target_fallback",
        "name_gene": "name_target_fallback",
    }
    for key_col in key_priority:
        if key_col not in out.columns:
            continue
        subset = out.loc[~filled & _valid_key(out[key_col]), [key_col]].copy()
        if subset.empty:
            continue
        key_type = key_col.removeprefix("_key_")
        ext_subset = external[external["external_join_key_type"].eq(key_type)].copy()
        ext_label = pd.to_numeric(ext_subset.get("external_four_state_ml_label"), errors="coerce")
        ext_subset = ext_subset.loc[ext_label.isin([0, 1])].copy()
        if ext_subset.empty:
            continue
        merged = subset.reset_index(names="_row_index").merge(
            ext_subset.drop_duplicates("_merge_key", keep="first"),
            left_on=key_col,
            right_on="_merge_key",
            how="left",
        )
        matched = merged["_merge_key"].notna()
        if not matched.any():
            continue
        row_index = merged.loc[matched, "_row_index"]
        for col in ext_cols:
            if col not in out.columns:
                out[col] = pd.NA
            values = merged.loc[matched, col]
            out.loc[row_index.to_numpy(), col] = values.to_numpy()
        out.loc[row_index.to_numpy(), "external_join_key_type"] = key_type
        out.loc[row_index.to_numpy(), "external_join_key_value"] = merged.loc[matched, key_col].to_numpy()
        out.loc[row_index.to_numpy(), "external_join_confidence"] = confidence[key_type]
        filled.loc[row_index.to_numpy()] = True
        join_counts[key_type] = int(matched.sum())
    return out, {
        "matched_rows": int(filled.sum()),
        "joined_external_rows": int(sum(join_counts.values())),
        "preserved_addon_rows": preserved_addon_rows,
        "unmatched_rows": int((~filled).sum()),
        "join_counts": join_counts,
    }


def _combine_activity_labels(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    spd = pd.to_numeric(out.get("spd_binding_label"), errors="coerce")
    ext = pd.to_numeric(out.get("external_four_state_ml_label"), errors="coerce")
    external_source_text = _first_nonempty(
        out,
        [
            "external_evidence_sources",
            "scenario_a_tier1_sources",
            "external_parent_sources",
        ],
    )
    external_family = external_source_text.map(_source_family_from_evidence)
    combined = pd.Series(pd.NA, index=out.index, dtype="Float64")
    status = pd.Series("unknown_no_label", index=out.index, dtype="object")
    spd_only = spd.notna() & ext.isna()
    ext_only = spd.isna() & ext.notna()
    concordant = spd.notna() & ext.notna() & spd.eq(ext)
    conflict = spd.notna() & ext.notna() & ~spd.eq(ext)
    combined.loc[spd_only] = spd.loc[spd_only]
    status.loc[spd_only] = "spd_only"
    combined.loc[ext_only] = ext.loc[ext_only]
    status.loc[ext_only] = "external_only"
    combined.loc[concordant] = spd.loc[concordant]
    status.loc[concordant] = "spd_external_concordant"
    combined.loc[conflict] = -1.0
    status.loc[conflict] = "spd_external_conflict_excluded"
    out["combined_activity_label"] = combined
    out["combined_activity_ml_label"] = combined.mask(combined.eq(-1))
    out["combined_activity_label_status"] = status
    out["combined_activity_label_policy"] = (
        "SPD binding labels and external four-state measured activity labels are combined only when non-conflicting; "
        "comparable SPD/external disagreements become -1 and are excluded from combined_activity_ml_label."
    )
    combined_source_family = pd.Series("unknown", index=out.index, dtype="object")
    combined_label_source = pd.Series("", index=out.index, dtype="object")
    combined_source_family.loc[spd_only] = "SPD"
    combined_label_source.loc[spd_only] = "SPD"
    combined_source_family.loc[ext_only] = external_family.loc[ext_only]
    combined_label_source.loc[ext_only] = external_source_text.loc[ext_only]
    combined_source_family.loc[concordant] = "SPD;" + external_family.loc[concordant]
    combined_label_source.loc[concordant] = "SPD;" + external_source_text.loc[concordant]
    combined_source_family.loc[conflict] = "conflict;" + external_family.loc[conflict]
    combined_label_source.loc[conflict] = "conflict;SPD;" + external_source_text.loc[conflict]
    out["external_activity_source_family"] = external_family.where(ext.isin([0, 1]), "")
    out["combined_activity_source_family"] = combined_source_family
    out["combined_activity_label_source"] = combined_label_source
    return out


def build_spd_external_four_state_merged_table(
    spd_table_path: str | Path,
    four_state_table_path: str | Path,
    out_path: str | Path,
    *,
    run_dir: str | Path | None = None,
    repo_root: str | Path | None = None,
    add_tissue_labels: bool = True,
    refresh_metadata: bool = True,
) -> dict[str, Any]:
    spd_path = Path(spd_table_path)
    four_path = Path(four_state_table_path)
    out = Path(out_path)
    spd = pd.read_csv(spd_path, low_memory=False)
    additive = pd.read_csv(four_path, low_memory=False)
    spd = _make_merge_keys(spd, table_role="spd")
    additive = _make_merge_keys(additive, table_role="external")
    external, collapse_summary = _collapse_external(additive)
    merged, join_summary = _join_external(spd, external)
    merged = _combine_activity_labels(merged)
    tissue_summary: dict[str, Any] = {"status": "not_requested"}
    if add_tissue_labels:
        merged, tissue_summary = add_independent_tissue_site_labels(merged, run_dir=run_dir)
    metadata_summary: dict[str, Any] = {"status": "not_requested"}
    if refresh_metadata:
        merged, metadata_summary = enrich_ml_feature_metadata(
            merged,
            run_dir=Path(run_dir) if run_dir else None,
            repo_root=Path(repo_root).resolve() if repo_root else None,
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    # Preserve join audit keys for reproducibility, but put them at the end.
    merged.to_csv(out, index=False)
    spd_label = pd.to_numeric(merged.get("spd_binding_label"), errors="coerce")
    ext_label = pd.to_numeric(merged.get("external_four_state_ml_label"), errors="coerce")
    combined = pd.to_numeric(merged.get("combined_activity_ml_label"), errors="coerce")
    tissue = pd.to_numeric(merged.get("tissue_site_label"), errors="coerce") if "tissue_site_label" in merged else pd.Series(pd.NA, index=merged.index)
    summary = {
        "spd_table": str(spd_path),
        "four_state_table": str(four_path),
        "output": str(out),
        "rows": int(len(merged)),
        "columns": int(len(merged.columns)),
        "collapse_external": collapse_summary,
        "join_external": join_summary,
        "labels": {
            "spd_binding_label": {
                "labelable": int(spd_label.notna().sum()),
                "positive": int(spd_label.eq(1).sum()),
                "negative": int(spd_label.eq(0).sum()),
            },
            "external_four_state_ml_label": {
                "labelable": int(ext_label.notna().sum()),
                "positive": int(ext_label.eq(1).sum()),
                "negative": int(ext_label.eq(0).sum()),
            },
            "combined_activity_ml_label": {
                "labelable": int(combined.notna().sum()),
                "positive": int(combined.eq(1).sum()),
                "negative": int(combined.eq(0).sum()),
                "conflict_excluded": int(pd.to_numeric(merged.get("combined_activity_label"), errors="coerce").eq(-1).sum()),
            },
            "tissue_site_label": {
                "labelable": int(tissue.notna().sum()),
                "positive": int(tissue.eq(1).sum()),
                "negative": int(tissue.eq(0).sum()),
                "conflict_excluded": int(tissue.eq(-1).sum()),
            },
        },
        "combined_activity_status_counts": {str(k): int(v) for k, v in merged["combined_activity_label_status"].value_counts(dropna=False).items()},
        "tissue_label_projection": tissue_summary,
        "metadata_refresh": metadata_summary,
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return summary
