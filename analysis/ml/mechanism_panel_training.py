from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.mechanism_panel_audit import filter_mechanism_panel, panel_component_masks
from analysis.ml.mechanism_recovery import write_mechanism_topk_recovery


DEFAULT_QUARANTINED_FEATURES = [
    "target_family",
    "protein_class",
    "site_relevance_score",
    "site_relevance_score_by_adr",
    "mechanism_graph_score",
    "mechanism_path_count",
    "drug_adr_known",
    "target_adr_known",
    "target_pathway_adr_link",
    "triad_complete",
    "target_adr_evidence",
    "pathway_evidence",
    "site_safety_risk_score",
    "ot_safety_liability_count",
    "label_source",
    "source_family",
    "upstream_source",
    "adr_site_group",
    "adr_site_mapping_source",
    "adr_site_mapping_confidence",
    "adr_site_terms_matched",
    "adr_site_match_score",
    "mechanism_label_status",
    "mechanism_label_source",
    "mechanism_label_projection_source",
    "projected_label_source",
]


def _first_available_label(df: pd.DataFrame, preferred: str) -> str:
    candidates = [
        preferred,
        "mechanism_pu_label",
        "mechanism_ml_label_clean",
        "mechanism_ml_label",
        "four_state_ml_label",
        "mechanism_label",
    ]
    for col in candidates:
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any():
            return col
    raise ValueError(f"no usable mechanism label found among: {', '.join(candidates)}")


def _label_summary(df: pd.DataFrame, label_col: str) -> dict[str, Any]:
    label = pd.to_numeric(df.get(label_col, pd.Series(pd.NA, index=df.index)), errors="coerce")
    return {
        "rows": int(len(df)),
        "labelable_rows": int(label.notna().sum()),
        "positive_rows": int(label.eq(1).sum()),
        "negative_rows": int(label.eq(0).sum()),
        "unknown_rows": int(label.isna().sum()),
        "positive_rate_labelable": float(label.mean()) if label.notna().any() else None,
        "unique_drugs": int(df["drug_id"].dropna().astype(str).nunique()) if "drug_id" in df.columns else 0,
        "unique_targets": int(df["target_id"].dropna().astype(str).nunique()) if "target_id" in df.columns else 0,
        "unique_pdbs": int(df["pdb_id"].dropna().astype(str).nunique()) if "pdb_id" in df.columns else 0,
    }


def _source_balance(df: pd.DataFrame, label_col: str) -> pd.DataFrame:
    if "label_source" not in df.columns:
        return pd.DataFrame()
    label = pd.to_numeric(df[label_col], errors="coerce")
    rows: list[dict[str, Any]] = []
    for source, group in df.assign(_label=label).groupby("label_source", dropna=False):
        y = pd.to_numeric(group["_label"], errors="coerce")
        rows.append(
            {
                "label_source": "" if pd.isna(source) else str(source),
                "rows": int(len(group)),
                "labelable_rows": int(y.notna().sum()),
                "positive_rows": int(y.eq(1).sum()),
                "negative_rows": int(y.eq(0).sum()),
                "unknown_rows": int(y.isna().sum()),
                "positive_rate_labelable": float(y.mean()) if y.notna().any() else None,
            }
        )
    return pd.DataFrame(rows).sort_values(["labelable_rows", "positive_rows"], ascending=[False, False])




def _join_feature_source(df: pd.DataFrame, feature_source: str | Path | None) -> pd.DataFrame:
    if feature_source is None:
        return df
    source_path = Path(feature_source)
    if not source_path.exists():
        raise FileNotFoundError(source_path)
    source = pd.read_csv(source_path, low_memory=False)
    key_cols = [col for col in ["drug_id", "target_id", "pdb_id", "ligand_base"] if col in df.columns and col in source.columns]
    if not key_cols:
        raise ValueError(f"feature source {source_path} has no shared join keys")
    feature_cols = [
        col
        for col in [
            "banana_score",
            "banana_binding_probability",
            "banana_score_normalized",
            "banana_score_feature_source",
            "banana_binding_probability_feature_source",
        ]
        if col in source.columns
    ]
    if not feature_cols:
        return df
    slim = source[key_cols + feature_cols].drop_duplicates(subset=key_cols)
    out = df.merge(slim, on=key_cols, how="left", suffixes=("", "_feature_source_joined"))
    for col in feature_cols:
        joined = f"{col}_feature_source_joined"
        if joined not in out.columns:
            continue
        if col in out.columns:
            out[col] = out[col].where(out[col].notna(), out[joined])
            out = out.drop(columns=[joined])
        else:
            out = out.rename(columns={joined: col})
    out["mechanism_panel_feature_source_table"] = str(source_path)
    out["mechanism_panel_feature_source_join_keys"] = ";".join(key_cols)
    return out

def build_mechanism_panel_training_table(
    source_table: str | Path,
    out_table: str | Path,
    *,
    panel: str = "cardiac_qt",
    source_label_col: str = "mechanism_pu_label",
    label_col: str | None = None,
    topk_out: str | Path | None = None,
    quarantine_features: list[str] | None = None,
    feature_source_table: str | Path | None = None,
) -> dict[str, Any]:
    source_path = Path(source_table)
    df = pd.read_csv(source_path, low_memory=False)
    base_label_col = _first_available_label(df, source_label_col)
    site_mask, term_mask = panel_component_masks(df, panel)
    panel_df = filter_mechanism_panel(df, panel)
    if panel_df.empty:
        raise ValueError(f"panel {panel!r} has no matching rows in {source_path}")
    panel_df = _join_feature_source(panel_df, feature_source_table)
    label_name = label_col or f"{panel}_mechanism_label"
    panel_df["mechanism_panel"] = panel
    panel_df[label_name] = pd.to_numeric(panel_df[base_label_col], errors="coerce")
    panel_df["mechanism_panel_source_label_col"] = base_label_col
    panel_df["mechanism_panel_match_basis"] = "unknown"
    both_idx = panel_df.index[site_mask.loc[panel_df.index] & term_mask.loc[panel_df.index]]
    site_idx = panel_df.index[site_mask.loc[panel_df.index] & ~term_mask.loc[panel_df.index]]
    term_idx = panel_df.index[~site_mask.loc[panel_df.index] & term_mask.loc[panel_df.index]]
    panel_df.loc[both_idx, "mechanism_panel_match_basis"] = "site_and_term"
    panel_df.loc[site_idx, "mechanism_panel_match_basis"] = "site"
    panel_df.loc[term_idx, "mechanism_panel_match_basis"] = "term"
    out = Path(out_table)
    out.parent.mkdir(parents=True, exist_ok=True)
    panel_df.to_csv(out, index=False)

    topk_path = Path(topk_out) if topk_out is not None else out.with_name(out.stem + ".topk_recovery.csv")
    topk = write_mechanism_topk_recovery(
        out,
        topk_path,
        label_col=label_name,
        score_cols=[
            "ml_prediction_score",
            "binding_expert_score",
            "consensus_score",
            "banana_score_normalized",
            "mechanism_graph_score",
            "site_relevance_score",
        ],
        group_cols=["target_id", "pdb_id", "mechanism_panel_match_basis", "all"],
        ks=[5, 10, 20, 50],
    )
    source_balance = _source_balance(panel_df, label_name)
    source_balance_path = out.with_name(out.stem + ".source_balance.csv")
    source_balance.to_csv(source_balance_path, index=False)
    quarantine = sorted(set(quarantine_features or DEFAULT_QUARANTINED_FEATURES))
    quarantine_path = out.with_name(out.stem + ".quarantined_features.json")
    quarantine_path.write_text(json.dumps(quarantine, indent=2, sort_keys=True), encoding="utf-8")
    manifest = {
        "source_table": str(source_path),
        "out_table": str(out),
        "panel": panel,
        "label_col": label_name,
        "source_label_col": base_label_col,
        "feature_source_table": str(feature_source_table) if feature_source_table is not None else None,
        "summary": _label_summary(panel_df, label_name),
        "topk_recovery": str(topk_path),
        "topk_rows": int(len(topk)),
        "source_balance": str(source_balance_path),
        "quarantined_features": str(quarantine_path),
        "training_policy": (
            "Panel/source/ADR/evidence fields are preserved for audit/reporting, "
            "but clean models must exclude quarantined features unless the run is "
            "explicitly marked as evidence-augmented sensitivity."
        ),
    }
    manifest_path = out.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
