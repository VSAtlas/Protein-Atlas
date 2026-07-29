from __future__ import annotations

from decimal import Decimal
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.spd import SPD_LABEL_POLICY_VERSION
from analysis.ml.spd_receptor_mapping import apply_receptor_target_contract
from analysis.spd_activity_policy import (
    label_spd_exposure_observation,
    parse_spd_activity_interval,
)


DEFAULT_SPD_PANEL_CANDIDATES = [
    Path("data/spd_full_panel_latest/spd_full_panel_relaxed/spd_full_panel_assay_pairs.csv"),
    Path("data/spd_full_panel_latest/spd_full_panel_strict/spd_full_panel_assay_pairs.csv"),
    Path("data/spd_full_panel_extended_v3/spd_full_panel_relaxed/spd_full_panel_assay_pairs.csv"),
    Path("data/spd_full_panel_extended_v3/spd_full_panel_strict/spd_full_panel_assay_pairs.csv"),
    Path("data/spd_full_panel_extended/spd_full_panel_relaxed/spd_full_panel_assay_pairs.csv"),
    Path("data/spd_full_panel_extended/spd_full_panel_strict/spd_full_panel_assay_pairs.csv"),
]
DEFAULT_TARGET_MAP_CANDIDATES = [
    Path("data/spd_full_panel_latest/spd_full_panel_relaxed/spd_target_selected.csv"),
    Path("data/spd_full_panel_latest/spd_full_panel_strict/spd_target_selected.csv"),
    Path("data/spd_full_panel_extended_v3/spd_full_panel_relaxed/spd_target_selected.csv"),
    Path("data/spd_full_panel_extended_v3/spd_full_panel_strict/spd_target_selected.csv"),
]
DEFAULT_TARGET_METADATA_CANDIDATES = [
    Path("analysis/gene_list/spd_assays_all.csv"),
    Path("data/spd_prep_demo/spd_assay_target_metadata.csv"),
]
DEFAULT_LIGAND_MAP_CANDIDATES = [
    Path("chemdb/data/fda_mapping_from_pdbqt.csv"),
    Path("data/fda_mapping_from_pdbqt.csv"),
]

SPD_RECEPTOR_STRICT_ELIGIBLE_COL = "spd_receptor_mapping_strict_eligible"
SPD_RECEPTOR_REVISED_IDENTITY_COLS = {
    "target_id": "spd_receptor_mapping_revised_target_id",
    "target_gene": "spd_receptor_mapping_revised_target_gene",
    "target_uniprot": "spd_receptor_mapping_revised_target_uniprot",
}
SPD_LABEL_REFRESH_COLS = [
    "spd_drug_id",
    "spd_target_id",
    "spd_ac50_uM",
    "free_cmax_um",
    "cmax_um",
    "exposure_margin",
    "spd_exposure_label",
    "spd_exposure_relevant",
    "spd_exposure_weak",
    "spd_exposure_unlikely",
    "spd_label_status",
    "spd_missing_reason",
    "spd_label_policy_version",
    "spd_activity_relation",
    "spd_target_protein_class",
    "spd_drugcentral_struct_id",
    "spd_inchikey",
    "spd_assay_count",
    "spd_assay_ids",
    "spd_assay_name",
    "spd_source",
]
SPD_STRICT_DIRECT_LABEL_REFRESH_COLS = [
    "spd_binding_label",
    "spd_activity_label",
    "activity_label",
    "ml_binary_label",
]
SPD_STRICT_TARGET_CONTEXT_COLS = ["protein_class", "target_family"]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_existing(explicit: str | Path | None, candidates: list[Path]) -> Path | None:
    if explicit is not None:
        path = Path(explicit)
        return path if path.exists() else None
    roots = [Path.cwd(), _repo_root()]
    for root in roots:
        for rel in candidates:
            path = rel if rel.is_absolute() else root / rel
            if path.exists():
                return path
    return None


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _lower(value: Any) -> str:
    return _clean(value).lower()


def _upper(value: Any) -> str:
    return _clean(value).upper()


def _bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.fillna(False).astype(bool)
    normalized = series.fillna("").astype(str).str.strip().str.lower()
    return normalized.isin({"1", "true", "yes"})


def _first_nonempty_frame(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in cols:
        if col not in df.columns:
            continue
        values = df[col].map(_clean)
        out = out.where(out.astype(str).str.len().gt(0), values)
    return out


def _read_fda_ligand_map(path: Path) -> pd.DataFrame:
    wanted = {
        "path",
        "ligand_base",
        "generic_name",
        "display_name",
        "smiles",
        "inchikey",
        "rxnorm_rxcui",
        "drugcentral_id",
    }
    raw = pd.read_csv(path, usecols=lambda col: col in wanted, low_memory=False)
    if "ligand_base" not in raw.columns:
        if "path" not in raw.columns:
            raise ValueError(f"FDA ligand map {path} lacks path/ligand_base columns")
        raw["ligand_base"] = raw["path"].astype(str).str.extract(r"([^/\\]+)\.pdbqt$", expand=False)
    for col in wanted:
        if col not in raw.columns:
            raw[col] = pd.NA
    raw["_mapped_drug_id"] = _first_nonempty_frame(raw, ["generic_name", "display_name", "ligand_base"])
    raw["_mapped_drug_key"] = raw["_mapped_drug_id"].map(_lower)
    raw = raw.sort_values(["_mapped_drug_key", "ligand_base"], na_position="last")
    return raw[
        [
            "ligand_base",
            "_mapped_drug_id",
            "_mapped_drug_key",
            "generic_name",
            "display_name",
            "smiles",
            "inchikey",
            "rxnorm_rxcui",
            "drugcentral_id",
        ]
    ].drop_duplicates("ligand_base", keep="first")


def _read_target_metadata(path: Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=["gene", "protein_class"])
    raw = pd.read_csv(path, low_memory=False)
    gene_col = None
    for col in ["HumanEntrezGeneSymbol(representative)", "EntrezGeneSymbol", "gene", "target_gene"]:
        if col in raw.columns:
            gene_col = col
            break
    class_col = None
    for col in ["Protein/Target/ProteinClass", "protein_class", "target_class"]:
        if col in raw.columns:
            class_col = col
            break
    if gene_col is None or class_col is None:
        return pd.DataFrame(columns=["gene", "protein_class"])
    out = pd.DataFrame({"gene": raw[gene_col].map(_upper), "protein_class": raw[class_col].map(_clean)})
    out = out[out["gene"].astype(str).str.len().gt(0)]
    return out.drop_duplicates("gene", keep="first")


def _read_target_map(path: Path, metadata_path: Path | None) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    if "pdb_id" not in raw.columns or "gene" not in raw.columns:
        raise ValueError(f"SPD target map {path} must contain pdb_id and gene")
    out = raw[["pdb_id", "gene"]].copy()
    out["pdb_id"] = out["pdb_id"].map(_upper)
    out["_mapped_target_id"] = out["gene"].map(_upper)
    out["_mapped_target_key"] = out["_mapped_target_id"]
    meta = _read_target_metadata(metadata_path)
    if not meta.empty:
        out = out.merge(meta, left_on="_mapped_target_id", right_on="gene", how="left", suffixes=("", "_meta"))
        out = out.drop(columns=[col for col in ["gene_meta"] if col in out.columns])
    if "protein_class" not in out.columns:
        out["protein_class"] = pd.NA
    out = out.rename(columns={"gene": "target_gene"})
    return out[["pdb_id", "target_gene", "_mapped_target_id", "_mapped_target_key", "protein_class"]].drop_duplicates("pdb_id", keep="first")


def _first_nonempty_series(left: pd.Series, right: pd.Series) -> pd.Series:
    left_text = left.map(_clean) if left is not None else pd.Series("", index=right.index)
    right_text = right.map(_clean)
    return left_text.where(left_text.astype(str).str.len().gt(0), right_text)


def _read_spd_panel(path: Path) -> pd.DataFrame:
    wanted = {
        "drug_id",
        "target_id",
        "assay_id",
        "assay_name",
        "ac50_nM",
        "free_cmax_nM",
        "total_cmax_nM",
        "source",
        "assay_count",
        "source_assay_ids",
        "exposure_margin",
        "spd_assayed",
        "spd_label_status",
        "spd_missing_reason",
        "spd_exposure_relevant",
        "spd_exposure_weak",
        "spd_exposure_unlikely",
        "drug_match_key",
        "target_match_key",
        "spd_label_policy_version",
        "spd_activity_relation",
        "spd_target_protein_class",
        "spd_drugcentral_struct_id",
        "spd_inchikey",
    }
    raw = pd.read_csv(path, usecols=lambda col: col in wanted, low_memory=False)
    for col in wanted:
        if col not in raw.columns:
            raw[col] = pd.NA
    raw["drug_match_key"] = raw["drug_match_key"].where(raw["drug_match_key"].notna(), raw["drug_id"]).map(_lower)
    raw["target_match_key"] = raw["target_match_key"].where(raw["target_match_key"].notna(), raw["target_id"]).map(_upper)
    raw["_ac50_num"] = pd.to_numeric(raw["ac50_nM"], errors="coerce")
    raw = raw.sort_values(["drug_match_key", "target_match_key", "_ac50_num"], na_position="last")

    def join_unique(values: pd.Series) -> str:
        seen: list[str] = []
        for value in values:
            text = _clean(value)
            if text and text not in seen:
                seen.append(text)
        return ";".join(seen)

    rows = []
    for (drug_key, target_key), group in raw.groupby(["drug_match_key", "target_match_key"], dropna=False):
        first = group.iloc[0]
        ac50 = pd.to_numeric(first.get("ac50_nM"), errors="coerce")
        free_cmax = pd.to_numeric(first.get("free_cmax_nM"), errors="coerce")
        total_cmax = pd.to_numeric(first.get("total_cmax_nM"), errors="coerce")
        margin = pd.to_numeric(first.get("exposure_margin"), errors="coerce")
        if pd.isna(margin) and pd.notna(ac50) and pd.notna(free_cmax) and float(free_cmax) > 0:
            margin = float(ac50) / float(free_cmax)
        status = _lower(first.get("spd_label_status"))
        flag_text = _lower(first.get("spd_exposure_relevant"))
        exposure_label = pd.NA
        if status == "labeled_relevant":
            exposure_label = 1
        elif status in {"labeled_weak", "labeled_unlikely"}:
            exposure_label = 0
        elif status in {"labeled_not_relevant", "labeled_censored_not_relevant", "labeled_censored_unlikely"}:
            exposure_label = 0
        elif flag_text in {"true", "1", "1.0", "yes"}:
            exposure_label = 1
        elif flag_text in {"false", "0", "0.0", "no"} and status.startswith("labeled_"):
            exposure_label = 0
        elif not status.startswith("unknown_"):
            interval = parse_spd_activity_interval(
                first.get("ac50_nM"),
                first.get("spd_activity_relation"),
                "nM",
            )
            canonical_exposure = label_spd_exposure_observation(
                interval,
                (Decimal(str(free_cmax)) / Decimal("1000"))
                if pd.notna(free_cmax)
                else pd.NA,
            )
            if canonical_exposure.numeric_label is not None:
                exposure_label = canonical_exposure.numeric_label
        rows.append(
            {
                "_spd_drug_match_key": drug_key,
                "_spd_target_match_key": target_key,
                "spd_drug_id": first.get("drug_id"),
                "spd_target_id": first.get("target_id"),
                "spd_ac50_uM": (float(ac50) / 1000.0) if pd.notna(ac50) else pd.NA,
                "free_cmax_um": (float(free_cmax) / 1000.0) if pd.notna(free_cmax) else pd.NA,
                "cmax_um": (float(total_cmax) / 1000.0) if pd.notna(total_cmax) else pd.NA,
                "exposure_margin": margin,
                "spd_exposure_label": exposure_label,
                "spd_exposure_relevant": first.get("spd_exposure_relevant"),
                "spd_exposure_weak": first.get("spd_exposure_weak"),
                "spd_exposure_unlikely": first.get("spd_exposure_unlikely"),
                "spd_label_status": first.get("spd_label_status"),
                "spd_missing_reason": first.get("spd_missing_reason"),
                "spd_label_policy_version": first.get("spd_label_policy_version", SPD_LABEL_POLICY_VERSION),
                "spd_activity_relation": first.get("spd_activity_relation"),
                "spd_target_protein_class": first.get("spd_target_protein_class"),
                "spd_drugcentral_struct_id": first.get("spd_drugcentral_struct_id"),
                "spd_inchikey": first.get("spd_inchikey"),
                "spd_assay_count": int(pd.to_numeric(group.get("assay_count"), errors="coerce").fillna(1).sum()),
                "spd_assay_ids": join_unique(group.get("source_assay_ids", pd.Series(dtype=object))),
                "spd_assay_name": join_unique(group.get("assay_name", pd.Series(dtype=object))),
                "spd_source": join_unique(group.get("source", pd.Series(dtype=object))),
            }
        )
    return pd.DataFrame(rows)


def _has_current_spd_label_data(df: pd.DataFrame) -> bool:
    if "spd_label_policy_version" not in df.columns:
        return False
    version = df["spd_label_policy_version"].fillna("").astype(str)
    return bool(version.eq(SPD_LABEL_POLICY_VERSION).any())


def _has_any_spd_label_data(df: pd.DataFrame) -> bool:
    for col in ["spd_ac50_uM", "exposure_margin", "spd_exposure_label", "spd_binding_label", "spd_exposure_relevant"]:
        if col in df.columns and df[col].notna().any():
            return True
    return False


def enrich_spd_labels_for_run_master(
    source: pd.DataFrame,
    *,
    spd_panel_path: str | Path | None = None,
    target_map_path: str | Path | None = None,
    ligand_map_path: str | Path | None = None,
    target_metadata_path: str | Path | None = None,
    receptor_mapping_mode: str = "legacy",
    receptor_mapping_contract_path: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Join local SPD assay labels to an Atlas run master table when needed.

    The join is intentionally conservative: FDA rows without an assayed SPD
    drug-target match remain unknown, and DUD/decoy rows remain benchmark-only
    rather than production negatives.
    """

    contract_input = source.copy()
    if "_spd_contract_input_order" in contract_input.columns:
        raise ValueError("Reserved column _spd_contract_input_order already exists")
    contract_input["_spd_contract_input_order"] = range(len(contract_input))
    receptor_mapping_kwargs: dict[str, Any] = {"mode": receptor_mapping_mode}
    if receptor_mapping_contract_path is not None:
        receptor_mapping_kwargs["contract_path"] = receptor_mapping_contract_path
    mapped_source, receptor_mapping_summary = apply_receptor_target_contract(
        contract_input,
        **receptor_mapping_kwargs,
    )
    expected_order = list(range(len(source)))
    observed_order = mapped_source.get(
        "_spd_contract_input_order", pd.Series(dtype="int64")
    ).tolist()
    if len(mapped_source) != len(source) or observed_order != expected_order:
        raise ValueError(
            "SPD receptor mapping contract changed source row count or order "
            f"(input rows={len(source)}, output rows={len(mapped_source)})"
        )
    source = mapped_source.drop(columns=["_spd_contract_input_order"])
    strict_mode = receptor_mapping_mode == "strict"
    summary: dict[str, Any] = {
        "status": "skipped",
        "reason": "not_needed",
        "input_rows": int(len(source)),
        "receptor_mapping": receptor_mapping_summary,
    }
    has_current_spd_labels = _has_current_spd_label_data(source)
    has_any_spd_labels = _has_any_spd_label_data(source)
    if has_current_spd_labels and not strict_mode:
        return source, summary
    if "pdb_id" not in source.columns or not ({"ligand_base", "ligand", "ligand_file"} & set(source.columns)):
        requires_rebuild = _bool_series(
            source.get(
                "spd_receptor_mapping_requires_target_derived_rebuild",
                pd.Series(False, index=source.index),
            )
        )
        if strict_mode and requires_rebuild.any():
            raise ValueError(
                "Strict receptor mapping requires target-derived label rebuilding, "
                "but the source lacks PDB or ligand join columns"
            )
        summary.update({"reason": "missing_pdb_or_ligand_columns"})
        return source, summary

    panel_path = _resolve_existing(spd_panel_path, DEFAULT_SPD_PANEL_CANDIDATES)
    target_path = _resolve_existing(target_map_path, DEFAULT_TARGET_MAP_CANDIDATES)
    ligand_path = _resolve_existing(ligand_map_path, DEFAULT_LIGAND_MAP_CANDIDATES)
    metadata_path = _resolve_existing(target_metadata_path, DEFAULT_TARGET_METADATA_CANDIDATES)
    missing = [
        name
        for name, value in [("spd_panel", panel_path), ("target_map", target_path), ("ligand_map", ligand_path)]
        if value is None
    ]
    if missing:
        requires_rebuild = _bool_series(
            source.get(
                "spd_receptor_mapping_requires_target_derived_rebuild",
                pd.Series(False, index=source.index),
            )
        )
        if strict_mode and requires_rebuild.any():
            raise ValueError(
                "Strict receptor mapping requires target-derived label rebuilding, "
                f"but mapping inputs are unavailable: {', '.join(missing)}"
            )
        summary.update({"reason": "missing_default_mapping_files", "missing_files": missing})
        return source, summary
    assert panel_path is not None
    assert target_path is not None
    assert ligand_path is not None

    full_source = source.copy()
    full_source["_spd_source_row_order"] = range(len(full_source))
    strict_eligible = pd.Series(True, index=full_source.index, dtype=bool)
    rebuild = pd.Series(False, index=full_source.index, dtype=bool)
    if strict_mode:
        if SPD_RECEPTOR_STRICT_ELIGIBLE_COL not in full_source.columns:
            raise ValueError(
                "Strict SPD receptor mapping did not produce "
                f"{SPD_RECEPTOR_STRICT_ELIGIBLE_COL}"
            )
        strict_eligible = _bool_series(full_source[SPD_RECEPTOR_STRICT_ELIGIBLE_COL])
        rebuild = full_source.get(
            "spd_receptor_mapping_requires_target_derived_rebuild",
            pd.Series(False, index=full_source.index),
        )
        rebuild = _bool_series(rebuild)
    selective_strict_rebuild = (
        strict_mode and has_current_spd_labels and rebuild.any()
    )
    unchanged = (
        full_source.loc[~rebuild].copy()
        if selective_strict_rebuild
        else pd.DataFrame()
    )
    out = (
        full_source.loc[rebuild].copy()
        if selective_strict_rebuild
        else full_source.copy()
    )
    if strict_mode:
        for col in SPD_STRICT_TARGET_CONTEXT_COLS:
            if col in out.columns:
                out.loc[
                    _bool_series(
                        out[
                            "spd_receptor_mapping_requires_target_derived_rebuild"
                        ]
                    ),
                    col,
                ] = pd.NA
    if has_any_spd_labels:
        refresh_cols = list(SPD_LABEL_REFRESH_COLS)
        if strict_mode:
            refresh_cols.extend(SPD_STRICT_DIRECT_LABEL_REFRESH_COLS)
        out = out.drop(
            columns=[col for col in refresh_cols if col in out.columns],
            errors="ignore",
        )
    if "ligand_base" not in out.columns:
        for col in ["ligand_file", "ligand"]:
            if col in out.columns:
                out["ligand_base"] = out[col].map(_clean).str.replace(r"(_dud_stage[0-9]+|\.pdbqt)$", "", regex=True)
                break
    out["pdb_id"] = out["pdb_id"].map(_upper)
    out["ligand_base"] = out["ligand_base"].map(_clean)

    ligand_map = _read_fda_ligand_map(ligand_path)
    target_map = _read_target_map(target_path, metadata_path)
    spd_panel = _read_spd_panel(panel_path)
    if has_any_spd_labels:
        # Historical model-ready tables already carry columns emitted by these
        # mapping joins. Any policy-version refresh must remove that stale
        # projection first; otherwise pandas can create duplicate suffixed
        # columns (for example ``protein_class_target_map``) whose selection
        # yields a DataFrame instead of a Series.
        ligand_projection_cols = set(ligand_map.columns) - {"ligand_base"}
        target_projection_cols = set(target_map.columns) - {"pdb_id"}
        stale_mapping_cols = (
            ligand_projection_cols
            | target_projection_cols
            | {f"{col}_fda_map" for col in ligand_projection_cols}
            | {f"{col}_target_map" for col in target_projection_cols}
        )
        out = out.drop(
            columns=[col for col in stale_mapping_cols if col in out.columns],
            errors="ignore",
        )

    out = out.merge(
        ligand_map,
        on="ligand_base",
        how="left",
        suffixes=("", "_fda_map"),
        validate="many_to_one",
    )
    out = out.merge(
        target_map,
        on="pdb_id",
        how="left",
        suffixes=("", "_target_map"),
        validate="many_to_one",
    )

    existing_drug = out["drug_id"] if "drug_id" in out.columns else pd.Series("", index=out.index)
    out["drug_id"] = _first_nonempty_series(existing_drug, out["_mapped_drug_id"])
    existing_target = out["target_id"] if "target_id" in out.columns else pd.Series("", index=out.index)
    out["target_id"] = _first_nonempty_series(existing_target, out["_mapped_target_id"])
    if "target_gene" not in out.columns:
        out["target_gene"] = out["_mapped_target_id"]
    else:
        out["target_gene"] = _first_nonempty_series(out["target_gene"], out["_mapped_target_id"])
    for col in ["generic_name", "display_name", "smiles", "inchikey"]:
        mapped_col = col
        if mapped_col in out.columns:
            # Existing column may have come from the ligand map merge. Keep it.
            out[col] = out[mapped_col]
    if "protein_class_target_map" in out.columns:
        out["protein_class"] = _first_nonempty_series(out.get("protein_class", pd.Series("", index=out.index)), out["protein_class_target_map"])
    if strict_mode:
        for target_col, revised_col in SPD_RECEPTOR_REVISED_IDENTITY_COLS.items():
            if revised_col in out.columns:
                revised = out[revised_col]
                has_revised_identity = revised.notna() & revised.map(_clean).ne("")
                out.loc[has_revised_identity, target_col] = revised.loc[
                    has_revised_identity
                ]
        rebuild = _bool_series(
            out.get(
                "spd_receptor_mapping_requires_target_derived_rebuild",
                pd.Series(False, index=out.index),
            )
        )
        if rebuild.any() and "protein_class" in out.columns:
            revised_metadata = _read_target_metadata(metadata_path)
            class_by_gene = (
                revised_metadata.drop_duplicates("gene", keep="first")
                .set_index("gene")["protein_class"]
                if not revised_metadata.empty
                else pd.Series(dtype="object")
            )
            revised_gene = out["target_gene"].map(_upper)
            out.loc[rebuild, "protein_class"] = revised_gene.loc[rebuild].map(
                class_by_gene
            )

    out["_spd_drug_match_key"] = _first_nonempty_frame(out, ["drug_id", "generic_name", "display_name", "_mapped_drug_id"]).map(_lower)
    out["_spd_target_match_key"] = _first_nonempty_frame(out, ["target_gene", "target_id", "_mapped_target_id"]).map(_upper)
    enriched_input_rows = len(out)
    out = out.merge(
        spd_panel,
        on=["_spd_drug_match_key", "_spd_target_match_key"],
        how="left",
        validate="many_to_one",
    )
    if out["_spd_source_row_order"].duplicated().any() or len(out) != enriched_input_rows:
        raise ValueError(
            "SPD enrichment multiplied source rows; receptor and panel joins must be many-to-one"
        )
    out = out.sort_values("_spd_source_row_order", kind="stable")

    decoy = pd.to_numeric(out.get("is_decoy", pd.Series(0, index=out.index)), errors="coerce").fillna(0).eq(1)
    out["benchmark_only"] = decoy
    out["training_allowed"] = ~decoy
    out["spd_auto_enriched"] = True
    out["spd_auto_enrichment_source"] = str(panel_path)
    out["spd_target_map_source"] = str(target_path)
    out["spd_ligand_map_source"] = str(ligand_path)
    out["spd_join_available"] = out["spd_ac50_uM"].notna()
    out["spd_label_policy_version"] = out["spd_label_policy_version"].fillna(SPD_LABEL_POLICY_VERSION)

    if selective_strict_rebuild:
        out = pd.concat([unchanged, out], axis=0, ignore_index=False, sort=False)
    out = out.sort_values("_spd_source_row_order", kind="stable").reset_index(drop=True)
    if len(out) != len(source):
        raise ValueError(
            "Selective strict SPD enrichment changed source row count "
            f"from {len(source)} to {len(out)}"
        )
    out = out.drop(columns=["_spd_source_row_order"])
    out = out.loc[:, ~out.columns.duplicated()].copy()
    joined_mask = _bool_series(out["spd_join_available"])
    summary.update(
        {
            "status": "enriched",
            "reason": "joined_local_spd_panel",
            "replaced_existing_spd_labels": bool(has_any_spd_labels),
            "selective_strict_rebuild": bool(selective_strict_rebuild),
            "strict_target_derived_rebuilt_rows": int(rebuild.sum()),
            "spd_panel_path": str(panel_path),
            "target_map_path": str(target_path),
            "ligand_map_path": str(ligand_path),
            "target_metadata_path": str(metadata_path) if metadata_path else None,
            "rows": int(len(out)),
            "mapped_ligand_rows": int(out["_mapped_drug_key"].notna().sum()),
            "mapped_target_rows": int(out["_mapped_target_key"].notna().sum()),
            "spd_joined_rows": int(joined_mask.sum()),
            "unique_spd_joined_pairs": int(
                out.loc[
                    joined_mask,
                    ["_spd_drug_match_key", "_spd_target_match_key"],
                ]
                .drop_duplicates()
                .shape[0]
            ),
            "decoy_rows_kept_benchmark_only": int(decoy.sum()),
            "strict_model_eligible_rows": int(strict_eligible.sum()) if strict_mode else int(len(out)),
            "strict_model_ineligible_rows": int((~strict_eligible).sum()) if strict_mode else 0,
        }
    )
    return out, summary
