from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def _first_nonempty(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in cols:
        if col not in df.columns:
            continue
        vals = df[col].fillna("").astype(str).str.strip()
        out = out.where(out.astype(str).str.len() > 0, vals)
    return out


def _drug_key(df: pd.DataFrame) -> pd.Series:
    return _first_nonempty(
        df,
        ["generic_name", "mapped_drug_name", "display_name", "ligand_display", "drug_id"],
    ).str.lower()


def _top_mean(values: pd.Series, n: int = 5) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().sort_values(ascending=False)
    if numeric.empty:
        return float("nan")
    return float(numeric.head(n).mean())


def _aggregate_pair_features(pair_table: pd.DataFrame, mechanism_scores: pd.DataFrame) -> pd.DataFrame:
    pairs = pair_table.copy()
    pairs["_adr_panel_drug_key"] = _drug_key(pairs)
    if not mechanism_scores.empty and {"drug_id", "target_id"}.issubset(pairs.columns) and {
        "drug_id",
        "target_id",
    }.issubset(mechanism_scores.columns):
        score_cols = [col for col in mechanism_scores.columns if col not in pairs.columns and col not in {"drug_id", "target_id"}]
        pairs = pairs.merge(
            mechanism_scores[["drug_id", "target_id", *score_cols]],
            on=["drug_id", "target_id"],
            how="left",
        )
    numeric_cols = [
        col
        for col in [
            "atlas_score",
            "consensus_score",
            "mechanism_graph_score",
            "mechanism_path_count",
            "drug_adr_known",
            "target_adr_known",
            "target_pathway_adr_link",
            "drug_target_known",
            "triad_complete",
        ]
        if col in pairs.columns
    ]
    rows: list[dict[str, Any]] = []
    for drug_key, group in pairs.groupby("_adr_panel_drug_key", dropna=False):
        if not str(drug_key).strip():
            continue
        row: dict[str, Any] = {
            "_adr_panel_drug_key": drug_key,
            "n_targets_screened": int(group["target_id"].nunique()) if "target_id" in group.columns else int(len(group)),
            "n_pdbs_screened": int(group["pdb_id"].nunique()) if "pdb_id" in group.columns else int(len(group)),
        }
        for col in numeric_cols:
            vals = pd.to_numeric(group[col], errors="coerce")
            row[f"{col}_max"] = float(vals.max()) if vals.notna().any() else float("nan")
            row[f"{col}_mean"] = float(vals.mean()) if vals.notna().any() else float("nan")
            row[f"{col}_top5_mean"] = _top_mean(vals, 5)
        rows.append(row)
    return pd.DataFrame(rows)


def build_adr_mechanism_panel_table(
    pair_table_path: str | Path,
    four_state_labels_path: str | Path,
    out_path: str | Path,
    *,
    mechanism_scores_path: str | Path | None = None,
) -> pd.DataFrame:
    """Build a drug-ADR keyed table for clinical ADR control calibration.

    This is intentionally separate from drug-target mechanism training. OMOP/OHDSI
    and FAERS non-signals are drug-outcome controls, not biochemical target labels.
    """

    pairs = pd.read_csv(pair_table_path, low_memory=False)
    labels = pd.read_csv(four_state_labels_path, low_memory=False)
    labels = labels[labels["pair_type"].astype(str).eq("drug_adr")].copy()
    labels["four_state_ml_label"] = pd.to_numeric(labels["four_state_ml_label"], errors="coerce")
    scores = (
        pd.read_csv(mechanism_scores_path, low_memory=False)
        if mechanism_scores_path is not None and Path(mechanism_scores_path).exists()
        else pd.DataFrame()
    )
    features = _aggregate_pair_features(pairs, scores)
    labels["_adr_panel_drug_key"] = labels["drug_id"].fillna("").astype(str).str.strip().str.lower()
    out = labels.merge(features, on="_adr_panel_drug_key", how="left")
    out["adr_mechanism_panel_label"] = out["four_state_ml_label"]
    out["has_atlas_pair_features"] = out["n_targets_screened"].notna() if "n_targets_screened" in out.columns else False
    out["adr_mechanism_panel_ml_label"] = out["adr_mechanism_panel_label"].where(out["has_atlas_pair_features"])
    out["label_scope"] = "drug_adr"
    out["label_use_warning"] = (
        "Drug-ADR controls calibrate clinical outcome evidence. They are not strict target-mechanism labels."
    )
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    manifest = {
        "rows": int(len(out)),
        "labelable_rows": int(out["adr_mechanism_panel_label"].notna().sum()),
        "atlas_feature_labelable_rows": int(out["adr_mechanism_panel_ml_label"].notna().sum()),
        "positive_rows": int(out["adr_mechanism_panel_label"].eq(1).sum()),
        "negative_rows": int(out["adr_mechanism_panel_label"].eq(0).sum()),
        "missing_feature_rows": int(out["n_targets_screened"].isna().sum()) if "n_targets_screened" in out.columns else int(len(out)),
        "sources": {str(k): int(v) for k, v in out["evidence_sources"].value_counts(dropna=False).items()},
        "warning": (
            "Use this for ADR-level calibration/control and MoE monitoring. Do not interpret it as direct proof that a "
            "specific target mediates an ADR unless a target-ADR panel table is added."
        ),
    }
    path.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return out
