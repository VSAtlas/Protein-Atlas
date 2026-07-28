from __future__ import annotations

from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.mapping import apply_pair_mapping, normalize_columns, read_optional_mapping
from analysis.external.source_tables import read_source_table
from analysis.spd_activity_policy import (
    label_spd_exposure_observation,
    parse_spd_activity_interval,
)
from analysis.statistics import ranked_binary_metrics


STATISTICAL_FORMULAS = (
    {
        "title": "SPD exposure relevance",
        "formula": "exposure_margin = AC50 / free_Cmax; spd_exposure_relevant = 1 if exposure_margin <= margin_strong, otherwise 0 when the measured/censored bound proves margin > margin_strong",
        "notes": "The default strong margin is 10. Missing SPD assays, AC50, or free-Cmax remain unknown rather than negative. An open lower bound at 10 has feasible margins strictly above 10 and is negative; a closed lower bound at 10 remains unknown because equality is possible.",
        "source": "analysis.external.spd.aggregate_spd_assays",
    },
)


SPD_LABEL_POLICY_VERSION = "spd_interval_censor_aware_v3"


SPD_ALIASES = {
    "drug_id": ["drug_id", "compound_id", "drug", "compound", "molecule_id", "drugcentral_struct_id", "name"],
    "target_id": ["target_id", "target", "gene", "gene_symbol", "uniprot", "assay_target", "assay_group_name"],
    "assay_id": ["assay_id", "assay", "protocol_id"],
    "assay_name": ["assay_name", "assay_title", "assay_group_name", "name"],
    "ac50_nM": ["ac50_nM", "ac50_nm", "ac50", "ac50_um", "potency_nM", "summarized ic50"],
    "free_cmax_nM": ["free_cmax_nM", "free_cmax_nm", "free_cmax", "free_cmax_um", "combined_free_cmax_um"],
    "total_cmax_nM": ["total_cmax_nM", "total_cmax_nm", "total_cmax", "cmax_nm", "combined_cmax_um"],
    "source": ["source", "dataset"],
    "assay_type": ["assay_type", "format", "assay_format"],
    "activity_relation": ["activity_relation", "standard_relation", "relation", "ac50_relation", "summarized prefix"],
}


def _read_table(path: str | Path) -> pd.DataFrame:
    table_path = Path(path)
    if table_path.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        try:
            assays = pd.read_excel(table_path, sheet_name="S Data 1", header=2)
            try:
                targets = pd.read_excel(table_path, sheet_name="S Data 3", header=2)
            except Exception:
                targets = pd.DataFrame()
            if not targets.empty and "assay name" in targets.columns:
                target_cols = [
                    col
                    for col in [
                        "assay name",
                        "assay group ID",
                        "preferred assay ID",
                        "Format",
                        "Mode",
                        "Readout",
                        "Technology",
                        "EntrezGeneID",
                        "EntrezGeneSymbol",
                        "Protein/Target/Name",
                        "Protein/Target/ProteinClass",
                        "Protein/Target/Species",
                        "HumanEntrezGeneSymbol(representative)",
                    ]
                    if col in targets.columns
                ]
                assays = assays.merge(
                    targets[target_cols].drop_duplicates("assay name"),
                    left_on="assay_group_name",
                    right_on="assay name",
                    how="left",
                )
            return assays
        except Exception:
            return read_source_table(table_path)
    return read_source_table(table_path)


def normalize_spd_table(spd_path: str | Path, mapping_path: str | Path | None = None) -> pd.DataFrame:
    raw = _read_table(spd_path)
    df = normalize_columns(raw, SPD_ALIASES)
    df = apply_pair_mapping(df, read_optional_mapping(mapping_path))
    for col in ("ac50_nM", "free_cmax_nM", "total_cmax_nM"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if "summarized IC50" in df.columns:
        df["ac50_nM"] = pd.to_numeric(df["summarized IC50"], errors="coerce") * 1000.0
    if "combined_free_Cmax_uM" in df.columns:
        df["free_cmax_nM"] = pd.to_numeric(df["combined_free_Cmax_uM"], errors="coerce") * 1000.0
    if "combined_Cmax_uM" in df.columns:
        df["total_cmax_nM"] = pd.to_numeric(df["combined_Cmax_uM"], errors="coerce") * 1000.0
    if "HumanEntrezGeneSymbol(representative)" in df.columns:
        representative = df["HumanEntrezGeneSymbol(representative)"].fillna("").astype(str).str.strip()
        df["target_id"] = df["target_id"].where(representative.eq(""), representative)
    if "EntrezGeneSymbol" in df.columns:
        gene = df["EntrezGeneSymbol"].fillna("").astype(str).str.strip()
        current = df["target_id"].fillna("").astype(str).str.strip()
        df["target_id"] = df["target_id"].where(current.ne("") & current.ne("nan"), gene)
    if "Format" in df.columns:
        df["assay_type"] = df["assay_type"].where(df["assay_type"].notna(), df["Format"])
    if "Protein/Target/ProteinClass" in df.columns:
        df["spd_target_protein_class"] = df["Protein/Target/ProteinClass"]
    if "drugcentral_struct_id" in df.columns:
        df["spd_drugcentral_struct_id"] = df["drugcentral_struct_id"]
    if "inchi_key" in df.columns:
        df["spd_inchikey"] = df["inchi_key"]
    if "ac50_um" in df.columns and df["ac50_nM"].isna().any():
        df["ac50_nM"] = df["ac50_nM"].fillna(pd.to_numeric(df["ac50_um"], errors="coerce") * 1000.0)
    if "free_cmax_um" in df.columns and df["free_cmax_nM"].isna().any():
        df["free_cmax_nM"] = df["free_cmax_nM"].fillna(pd.to_numeric(df["free_cmax_um"], errors="coerce") * 1000.0)
    for col in ("ac50_nM", "free_cmax_nM", "total_cmax_nM"):
        source_um = col.replace("_nM", "_um")
        if source_um in df.columns:
            df[col] = df[col].where(df[col].notna(), pd.to_numeric(df[source_um], errors="coerce") * 1000.0)
    df["source"] = df["source"].fillna("SPD")
    return df


def _spd_label_status(
    row: pd.Series | dict[str, Any],
    margin_strong: float,
    margin_weak: float,
) -> str:
    ac50 = row.get("ac50_nM")
    free_cmax = row.get("free_cmax_nM")
    if pd.isna(ac50):
        return "unknown_missing_ac50"
    try:
        free_cmax_numeric = Decimal(str(free_cmax))
    except (InvalidOperation, TypeError, ValueError):
        return "unknown_missing_free_cmax"
    if not free_cmax_numeric.is_finite() or free_cmax_numeric <= 0:
        return "unknown_missing_free_cmax"
    free_cmax_um = free_cmax_numeric / Decimal("1000")

    interval = parse_spd_activity_interval(
        ac50,
        row.get("activity_relation"),
        "nM",
    )
    exposure = label_spd_exposure_observation(
        interval,
        free_cmax_um,
        threshold=margin_strong,
    )
    if not interval.is_valid:
        return f"unknown_invalid_activity_{interval.activity_parse_reason}"
    if exposure.numeric_label == 1:
        return "labeled_relevant"
    if exposure.numeric_label is None:
        if interval.normalized_relation in {">", ">="}:
            return "unknown_censored_ac50_gt_crosses_relevant_threshold"
        return "unknown_censored_activity_interval_crosses_relevant_threshold"

    weak_exposure = label_spd_exposure_observation(
        interval,
        free_cmax_um,
        threshold=margin_weak,
    )
    relation = interval.normalized_relation
    if relation in {">", ">="}:
        return (
            "labeled_censored_unlikely"
            if weak_exposure.numeric_label == 0
            else "labeled_censored_not_relevant"
        )
    return (
        "labeled_unlikely"
        if weak_exposure.numeric_label == 0
        else "labeled_weak"
    )


def aggregate_spd_assays(df: pd.DataFrame, margin_strong: float = 10.0, margin_weak: float = 100.0) -> pd.DataFrame:
    sortable = df.copy()
    sortable["_is_biochemical"] = sortable.get("assay_type", "").astype(str).str.lower().str.contains("biochem|target", regex=True, na=False)
    sortable = sortable.sort_values(["drug_id", "target_id", "_is_biochemical", "ac50_nM"], ascending=[True, True, False, True])
    rows: list[dict[str, object]] = []
    for (drug_id, target_id), group in sortable.groupby(["drug_id", "target_id"], dropna=False):
        best = group.iloc[0]
        free_cmax = best["free_cmax_nM"]
        ac50 = best["ac50_nM"]
        margin = ac50 / free_cmax if pd.notna(ac50) and pd.notna(free_cmax) and free_cmax > 0 else pd.NA
        label_row = best.copy()
        label_row["exposure_margin"] = margin
        status = _spd_label_status(label_row, margin_strong, margin_weak)
        labelable = status.startswith("labeled_")
        relevant = (
            True if status == "labeled_relevant" else (False if labelable else pd.NA)
        )
        weak = (
            status == "labeled_weak"
            if labelable and not status.startswith("labeled_censored_")
            else pd.NA
        )
        unlikely = (
            status in {"labeled_unlikely", "labeled_censored_unlikely"}
            if labelable and status != "labeled_censored_not_relevant"
            else pd.NA
        )
        rows.append(
            {
                "drug_id": drug_id,
                "target_id": target_id,
                "assay_id": best.get("assay_id", pd.NA),
                "assay_name": best.get("assay_name", pd.NA),
                "ac50_nM": ac50,
                "free_cmax_nM": free_cmax,
                "total_cmax_nM": best.get("total_cmax_nM", pd.NA),
                "source": best.get("source", "SPD"),
                "assay_count": len(group),
                "source_assay_ids": ";".join(sorted({str(v) for v in group["assay_id"].dropna()})),
                "exposure_margin": margin,
                "spd_assayed": True,
                "spd_label_status": status,
                "spd_missing_reason": "" if labelable else status,
                "spd_exposure_relevant": relevant,
                "spd_exposure_weak": weak,
                "spd_exposure_unlikely": unlikely,
                "spd_label_policy_version": SPD_LABEL_POLICY_VERSION,
                "spd_activity_relation": best.get("activity_relation", pd.NA),
                "spd_target_protein_class": best.get("spd_target_protein_class", pd.NA),
                "spd_drugcentral_struct_id": best.get("spd_drugcentral_struct_id", pd.NA),
                "spd_inchikey": best.get("spd_inchikey", pd.NA),
            }
        )
    return pd.DataFrame(rows)


def _first_nonempty(df: pd.DataFrame, candidates: tuple[str, ...]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in candidates:
        if col not in df.columns:
            continue
        values = df[col].fillna("").astype(str).str.strip()
        out = out.where(out.astype(str).str.len() > 0, values)
    return out


def _truthy_series(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"1", "1.0", "true", "yes", "y", "positive", "active"})


def _falsey_series(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"0", "0.0", "false", "no", "n", "negative", "inactive"})


def deduplicate_spd_benchmark(
    benchmark_table_path: str | Path,
    out_path: str | Path,
    *,
    ml_ready_out_path: str | Path | None = None,
) -> pd.DataFrame:
    """Collapse SPD benchmark rows to one row per drug-target pair for ML.

    The full benchmark may contain repeated ligand forms, PDB variants, and pH
    states. This keeps missing SPD evidence unknown, treats weak-margin rows as
    ambiguous for binary ML, and writes an optional supervised-eligible subset.
    """

    df = _read_table(benchmark_table_path)
    if df.empty:
        out = df.copy()
    else:
        work = df.copy()
        work["dedup_drug_key"] = _first_nonempty(
            work,
            (
                "drug_id",
                "drugcentral_id",
                "mapped_drug_name",
                "spd_drug_name",
                "ligand_display",
                "ligand_base",
            ),
        ).str.lower()
        work["dedup_target_key"] = _first_nonempty(
            work, ("target_id", "target_gene", "target_uniprot", "pdb_id")
        ).str.upper()
        work = work[(work["dedup_drug_key"] != "") & (work["dedup_target_key"] != "")]
        score_col = "z_selected" if "z_selected" in work.columns else "atlas_score"
        if score_col in work.columns:
            work["_dedup_score"] = pd.to_numeric(work[score_col], errors="coerce")
        else:
            work["_dedup_score"] = pd.NA
        work["_dedup_margin"] = pd.to_numeric(
            work.get("exposure_margin", pd.Series(pd.NA, index=work.index)),
            errors="coerce",
        )
        status = work.get("spd_label_status", pd.Series("", index=work.index)).astype(str)
        work["_dedup_label_rank"] = 4
        work["_dedup_relevant_true"] = _truthy_series(
            work.get("spd_exposure_relevant", pd.Series("", index=work.index))
        )
        work["_dedup_relevant_false"] = _falsey_series(
            work.get("spd_exposure_relevant", pd.Series("", index=work.index))
        )
        work["_dedup_weak_true"] = _truthy_series(
            work.get("spd_exposure_weak", pd.Series("", index=work.index))
        )
        work["_dedup_weak_false"] = _falsey_series(
            work.get("spd_exposure_weak", pd.Series("", index=work.index))
        )
        work["_dedup_unlikely_true"] = _truthy_series(
            work.get("spd_exposure_unlikely", pd.Series("", index=work.index))
        )
        work["_dedup_unlikely_false"] = _falsey_series(
            work.get("spd_exposure_unlikely", pd.Series("", index=work.index))
        )
        work.loc[work["_dedup_relevant_true"], "_dedup_label_rank"] = 0
        work.loc[
            (work["_dedup_label_rank"] > 0) & work["_dedup_unlikely_true"],
            "_dedup_label_rank",
        ] = 1
        work.loc[
            (work["_dedup_label_rank"] > 1) & work["_dedup_weak_true"],
            "_dedup_label_rank",
        ] = 2
        work.loc[(work["_dedup_label_rank"] > 2) & status.str.startswith("labeled"), "_dedup_label_rank"] = 3
        work = work.sort_values(
            ["dedup_drug_key", "dedup_target_key", "_dedup_label_rank", "_dedup_margin", "_dedup_score"],
            ascending=[True, True, True, True, False],
        )
        keys = ["dedup_drug_key", "dedup_target_key"]
        chosen = work.drop_duplicates(keys, keep="first").copy()
        agg_spec: dict[str, tuple[str, str]] = {
            "dedup_source_rows": ("dedup_target_key", "size"),
            "_relevant_true": ("_dedup_relevant_true", "max"),
            "_relevant_false": ("_dedup_relevant_false", "max"),
            "_weak_true": ("_dedup_weak_true", "max"),
            "_weak_false": ("_dedup_weak_false", "max"),
            "_unlikely_true": ("_dedup_unlikely_true", "max"),
            "_unlikely_false": ("_dedup_unlikely_false", "max"),
        }
        if "pdb_id" in work.columns:
            agg_spec["dedup_n_pdbs"] = ("pdb_id", "nunique")
        if "ligand_base" in work.columns:
            agg_spec["dedup_n_ligand_forms"] = ("ligand_base", "nunique")
        grouped = work.groupby(keys, dropna=False).agg(**agg_spec).reset_index()
        out = chosen.merge(grouped, on=keys, how="left")

        def label_from_flags(true_col: str, false_col: str) -> pd.Series:
            label = pd.Series("", index=out.index, dtype="object")
            label = label.mask(out[false_col].fillna(False), "0")
            label = label.mask(out[true_col].fillna(False), "1")
            return label

        out["spd_exposure_relevant"] = label_from_flags("_relevant_true", "_relevant_false")
        out["spd_exposure_weak"] = label_from_flags("_weak_true", "_weak_false")
        out["spd_exposure_unlikely"] = label_from_flags("_unlikely_true", "_unlikely_false")
        out["ml_binary_label"] = ""
        out["ml_supervised_eligible"] = "0"
        out["ml_exclude_reason"] = (
            out.get("spd_missing_reason", pd.Series("", index=out.index))
            .fillna("")
            .astype(str)
        )
        if "missing_reason" in out.columns:
            out["ml_exclude_reason"] = out["ml_exclude_reason"].where(
                out["ml_exclude_reason"].astype(str).str.len() > 0,
                out["missing_reason"].fillna("").astype(str),
            )
        out["ml_exclude_reason"] = out["ml_exclude_reason"].where(
            out["ml_exclude_reason"].astype(str).str.len() > 0,
            "unknown_no_spd_binary_label",
        )
        relevant_mask = out["spd_exposure_relevant"].eq("1")
        negative_mask = out["spd_exposure_relevant"].eq("0") & ~relevant_mask
        unlikely_mask = out["spd_exposure_unlikely"].eq("1") & ~relevant_mask
        weak_mask = out["spd_exposure_weak"].eq("1") & ~relevant_mask & ~unlikely_mask
        out.loc[relevant_mask, ["spd_label_status", "ml_binary_label", "ml_supervised_eligible", "ml_exclude_reason"]] = [
            "labeled_relevant",
            "1",
            "1",
            "",
        ]
        out.loc[negative_mask, ["spd_label_status", "ml_binary_label", "ml_supervised_eligible", "ml_exclude_reason"]] = [
            "labeled_not_relevant",
            "0",
            "1",
            "",
        ]
        out.loc[unlikely_mask, ["spd_label_status", "ml_binary_label", "ml_supervised_eligible", "ml_exclude_reason"]] = [
            "labeled_unlikely",
            "0",
            "1",
            "",
        ]
        out.loc[weak_mask, ["spd_label_status", "ml_binary_label", "ml_supervised_eligible", "ml_exclude_reason"]] = [
            "labeled_weak",
            "0",
            "1",
            "",
        ]
        out = out.drop(columns=[c for c in out.columns if c.startswith("_dedup_") or c in {"_relevant_true", "_relevant_false", "_weak_true", "_weak_false", "_unlikely_true", "_unlikely_false"}])
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    summary = pd.DataFrame(
        [
            {
                "n_dedup_rows": len(out),
                "n_supervised_eligible": int(out.get("ml_supervised_eligible", pd.Series(dtype=object)).astype(str).eq("1").sum()),
                "n_positive": int(out.get("ml_binary_label", pd.Series(dtype=object)).astype(str).eq("1").sum()),
                "n_negative": int(out.get("ml_binary_label", pd.Series(dtype=object)).astype(str).eq("0").sum()),
                "n_weak_ambiguous": int(out.get("ml_exclude_reason", pd.Series(dtype=object)).astype(str).eq("weak_margin_ambiguous").sum()),
            }
        ]
    )
    summary.to_csv(path.with_name(path.stem + ".summary.csv"), index=False)
    if ml_ready_out_path is not None:
        ml_ready = out[out.get("ml_supervised_eligible", pd.Series(dtype=object)).astype(str).eq("1")].copy()
        ml_path = Path(ml_ready_out_path)
        ml_path.parent.mkdir(parents=True, exist_ok=True)
        ml_ready.to_csv(ml_path, index=False)
    return out


def build_spd_benchmark(
    pair_table_path: str | Path,
    spd_path: str | Path,
    mapping_path: str | Path | None,
    out_path: str | Path,
    margin_strong: float = 10.0,
    margin_weak: float = 100.0,
) -> pd.DataFrame:
    pair = _read_table(pair_table_path)
    spd = aggregate_spd_assays(normalize_spd_table(spd_path, mapping_path), margin_strong, margin_weak)
    merged = pair.merge(spd, on=["drug_id", "target_id"], how="left", suffixes=("", "_spd"))
    no_match = merged["assay_count"].isna() if "assay_count" in merged else pd.Series(False, index=merged.index)
    merged.loc[no_match, "spd_assayed"] = False
    merged.loc[no_match, "spd_label_status"] = "unknown_no_spd_pair_match"
    merged.loc[no_match, "spd_missing_reason"] = "unknown_no_spd_pair_match"
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, index=False)
    deduplicate_spd_benchmark(
        out,
        out.with_name(out.stem + "_dedup_drug_target.csv"),
        ml_ready_out_path=out.with_name(out.stem + "_ml_ready.csv"),
    )
    summary = pd.DataFrame(
        [
            {
                "n_pair_rows": len(pair),
                "n_spd_pairs": len(spd),
                "n_joined_pairs": int(merged["assay_count"].notna().sum()) if "assay_count" in merged else 0,
                "n_labelable_pairs": int(merged["spd_exposure_relevant"].notna().sum()) if "spd_exposure_relevant" in merged else 0,
                "n_exposure_relevant": int(merged["spd_exposure_relevant"].eq(True).sum()) if "spd_exposure_relevant" in merged else 0,
                "label_status_counts": ";".join(
                    f"{key}:{value}" for key, value in merged.get("spd_label_status", pd.Series(dtype=object)).value_counts(dropna=False).items()
                ),
            }
        ]
    )
    summary.to_csv(out.with_name("spd_mapping_summary.csv"), index=False)
    return merged


def spd_metrics(df: pd.DataFrame, score_col: str = "atlas_score", label_col: str = "spd_exposure_relevant") -> dict[str, float]:
    scored = df[[score_col, label_col]].dropna()
    scores = pd.to_numeric(scored[score_col], errors="coerce")
    labels = scored[label_col].astype(bool).astype(int)
    keep = scores.notna()
    return ranked_binary_metrics(scores[keep].tolist(), labels[keep].tolist(), [0.01, 0.05, 0.10])


def write_spd_plots(df: pd.DataFrame, figures_dir: str | Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(figures_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if "exposure_margin" in df:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
        pd.to_numeric(df["exposure_margin"], errors="coerce").dropna().plot.hist(ax=ax, bins=40)
        ax.set_xlabel("AC50 / free Cmax")
        fig.tight_layout()
        fig.savefig(out_dir / "spd_exposure_margin_distribution.png")
        plt.close(fig)
    if {"atlas_score", "exposure_margin"}.issubset(df.columns):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=160)
        ax.scatter(pd.to_numeric(df["atlas_score"], errors="coerce"), pd.to_numeric(df["exposure_margin"], errors="coerce"), s=12, alpha=0.7)
        ax.set_xlabel("atlas_score")
        ax.set_ylabel("SPD exposure margin")
        fig.tight_layout()
        fig.savefig(out_dir / "atlas_score_vs_spd_margin.png")
        plt.close(fig)
