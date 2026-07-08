from __future__ import annotations

from pathlib import Path
from typing import Any, TypedDict

import pandas as pd

from analysis.reporting.value_utils import TRUTHY_PUBLICATION_STRINGS
from analysis.target_ids import build_target_id


class PublicationReadinessMetric(TypedDict):
    metric: str
    value: int


_TARGET_COMPONENT_COLUMNS = ("pdb_id", "variant", "ph_label")
_REQUIRED_PUBLICATION_COLUMNS = ("ligand_base",)


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _truthy(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(TRUTHY_PUBLICATION_STRINGS)
    )


def _numeric_present(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return pd.to_numeric(df[col], errors="coerce").notna()


def _series(df: pd.DataFrame, col: str, default: str = "") -> pd.Series:
    if col in df.columns:
        return df[col]
    return pd.Series(default, index=df.index)


def _target_id_series(df: pd.DataFrame) -> pd.Series:
    if "target_id" in df.columns:
        values = df["target_id"].fillna("").astype(str).str.strip()
        if values.ne("").any():
            return values
    return df.apply(
        lambda row: build_target_id(
            _text(row.get("pdb_id")),
            _text(row.get("variant")),
            _text(row.get("ph_label")),
        ),
        axis=1,
    )


def _merge_master_columns(heatmap: pd.DataFrame, master: pd.DataFrame) -> pd.DataFrame:
    heatmap = heatmap.copy()
    master = master.copy()
    heatmap["_publication_target_id"] = _target_id_series(heatmap)
    master["_publication_target_id"] = _target_id_series(master)
    key_cols = ["_publication_target_id", "ligand_base", "library"]
    if not set(key_cols).issubset(master.columns) or not set(key_cols).issubset(heatmap.columns):
        return heatmap
    master = master.drop_duplicates(key_cols, keep="first")
    extra_cols = [
        col
        for col in (
            "consensus_score",
            "consensus_score_pre",
            "final_score",
            "SCORCH_score_used",
            "scorch_composite",
            "free_cmax_um",
            "cmax_um",
            "fraction_unbound_plasma",
            "exposure_source",
            "spd_exposure_relevant",
            "spd_exposure_weak",
            "spd_exposure_unlikely",
            "spd_label_status",
            "ml_binary_label",
            "ml_supervised_eligible",
            "ml_exclude_reason",
        )
        if col in master.columns and col not in heatmap.columns
    ]
    if not extra_cols:
        return heatmap
    return heatmap.merge(master[key_cols + extra_cols], on=key_cols, how="left")


def _merge_optional_spd_columns(heatmap: pd.DataFrame, data_dir: Path) -> pd.DataFrame:
    candidates = [
        data_dir / "spd_atlas_benchmark_dedup_drug_target_ml_annotated.csv",
        data_dir / "spd_atlas_benchmark_dedup_drug_target.csv",
        data_dir / "spd_atlas_benchmark_ml_ready_supervised.csv",
        data_dir / "spd_atlas_benchmark.csv",
    ]
    spd_path = next((path for path in candidates if path.exists() and path.stat().st_size > 0), None)
    if spd_path is None:
        return heatmap
    spd = pd.read_csv(spd_path)
    keys = ["pdb_id", "variant", "ph_label", "ligand_base"]
    if not set(keys).issubset(heatmap.columns) or not set(keys).issubset(spd.columns):
        return heatmap
    spd_cols = [
        col
        for col in (
            "exposure_margin",
            "spd_exposure_relevant",
            "spd_exposure_weak",
            "spd_exposure_unlikely",
            "spd_label_status",
            "spd_ac50_uM",
            "spd_assay_count",
            "spd_assay_ids",
            "spd_assay_name",
            "ml_binary_label",
            "ml_supervised_eligible",
            "ml_exclude_reason",
        )
        if col in spd.columns and col not in heatmap.columns
    ]
    if not spd_cols:
        return heatmap
    spd = spd.drop_duplicates(keys, keep="first")
    return heatmap.merge(spd[keys + spd_cols], on=keys, how="left")


def build_publication_heatmap_input(
    heatmap_input_path: str | Path,
    master_rows_path: str | Path,
) -> pd.DataFrame:
    heatmap = pd.read_csv(heatmap_input_path)
    missing = [col for col in _REQUIRED_PUBLICATION_COLUMNS if col not in heatmap.columns]
    has_target_id = "target_id" in heatmap.columns
    has_target_components = all(col in heatmap.columns for col in _TARGET_COMPONENT_COLUMNS)
    if missing or (not has_target_id and not has_target_components):
        missing_msg = ", ".join(missing) if missing else "none"
        raise ValueError(
            "publication heatmap input missing required columns: "
            f"missing={missing_msg}; requires either target_id or "
            f"{', '.join(_TARGET_COMPONENT_COLUMNS)}"
        )
    master_path = Path(master_rows_path)
    if master_path.exists():
        heatmap = _merge_master_columns(heatmap, pd.read_csv(master_path))
        heatmap = _merge_optional_spd_columns(heatmap, master_path.parent)
    heatmap["atlas_score"] = pd.to_numeric(_series(heatmap, "z_selected"), errors="coerce")
    heatmap["atlas_score_source"] = _series(heatmap, "z_selected_source")
    lig_disp = _series(heatmap, "ligand_display")
    heatmap["drug_name"] = lig_disp.where(
        lig_disp.fillna("").astype(str).str.len() > 0,
        _series(heatmap, "ligand_base"),
    )
    heatmap["drug_ligand_id"] = _series(heatmap, "ligand_base")
    heatmap["full_library_bh_q"] = pd.to_numeric(_series(heatmap, "fdr_q_target"), errors="coerce")
    heatmap["scorch_conditional_q"] = pd.to_numeric(
        _series(heatmap, "scorch_fdr_decoy_competition_q_plus1"),
        errors="coerce",
    ).fillna(pd.to_numeric(_series(heatmap, "scorch_fdr_decoy_competition_q"), errors="coerce"))
    heatmap["has_full_library_fdr"] = _numeric_present(heatmap, "fdr_q_target").astype(int)
    heatmap["has_scorch_q"] = heatmap["scorch_conditional_q"].notna().astype(int)
    heatmap["has_free_cmax"] = _numeric_present(heatmap, "free_cmax_um").astype(int)
    heatmap["is_ml_supervised_eligible"] = (
        _truthy(heatmap["ml_supervised_eligible"])
        if "ml_supervised_eligible" in heatmap.columns
        else pd.Series(False, index=heatmap.index)
    ).astype(int)
    cols = [
        "drug_name",
        "drug_ligand_id",
        "target_id",
        "target_name",
        "pdb_id",
        "variant",
        "ph_label",
        "library",
        "rank",
        "pct_rank",
        "atlas_score",
        "atlas_score_source",
        "consensus_score",
        "final_score",
        "SCORCH_score_used",
        "full_library_bh_q",
        "fdr_null_source",
        "fdr_n_decoys",
        "scorch_conditional_q",
        "scorch_fdr_scope",
        "scorch_fdr_n_decoys",
        "scorch_fdr_n_tested",
        "free_cmax_um",
        "cmax_um",
        "fraction_unbound_plasma",
        "exposure_source",
        "spd_exposure_relevant",
        "spd_exposure_weak",
        "spd_exposure_unlikely",
        "spd_label_status",
        "ml_binary_label",
        "ml_supervised_eligible",
        "ml_exclude_reason",
        "has_full_library_fdr",
        "has_scorch_q",
        "has_free_cmax",
        "is_ml_supervised_eligible",
        "is_decoy",
        "is_control",
    ]
    return heatmap[[col for col in cols if col in heatmap.columns]].copy()


def build_heatmap_ml_readiness_summary(publication_df: pd.DataFrame) -> pd.DataFrame:
    is_decoy = _truthy(publication_df.get("is_decoy", pd.Series("", index=publication_df.index)))
    is_fda = ~is_decoy
    metrics: dict[str, int] = {
        "n_rows": len(publication_df),
        "n_fda_rows": int(is_fda.sum()),
        "n_decoy_rows": int(is_decoy.sum()),
        "n_targets": int(publication_df.get("target_id", pd.Series(dtype=object)).nunique()),
        "n_ligands": int(publication_df.get("drug_ligand_id", pd.Series(dtype=object)).nunique()),
        "n_atlas_score_present": int(_numeric_present(publication_df, "atlas_score").sum()),
        "n_full_library_fdr_present": int(_numeric_present(publication_df, "full_library_bh_q").sum()),
        "n_full_library_q10_hits": int((pd.to_numeric(publication_df.get("full_library_bh_q"), errors="coerce") <= 0.10).sum()),
        "n_scorch_q_present": int(_numeric_present(publication_df, "scorch_conditional_q").sum()),
        "n_scorch_q10_hits": int((pd.to_numeric(publication_df.get("scorch_conditional_q"), errors="coerce") <= 0.10).sum()),
        "n_free_cmax_present": int(_numeric_present(publication_df, "free_cmax_um").sum()),
        "n_spd_label_present": int(
            publication_df.get("spd_label_status", pd.Series("", index=publication_df.index)).fillna("").astype(str).str.startswith("labeled").sum()
        ),
        "n_ml_supervised_eligible": int(
            _truthy(publication_df.get("ml_supervised_eligible", pd.Series("", index=publication_df.index))).sum()
        ),
    }
    rows: list[PublicationReadinessMetric] = [
        {"metric": key, "value": value} for key, value in metrics.items()
    ]
    return pd.DataFrame(rows)


def write_publication_exports(
    heatmap_input_path: str | Path,
    master_rows_path: str | Path,
    out_dir: str | Path,
    *,
    suffix: str = "",
) -> tuple[Path, Path]:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    stem_suffix = f"_{suffix}" if suffix else ""
    publication_path = out_path / f"publication_heatmap_input{stem_suffix}.csv"
    summary_path = out_path / f"heatmap_ml_readiness_summary{stem_suffix}.csv"
    publication = build_publication_heatmap_input(heatmap_input_path, master_rows_path)
    publication.to_csv(publication_path, index=False)
    build_heatmap_ml_readiness_summary(publication).to_csv(summary_path, index=False)
    return publication_path, summary_path
