from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_WEIGHTS = {
    "binding_activity": 0.25,
    "exposure_relevance": 0.25,
    "tissue_site_relevance": 0.20,
    "mechanism_adr_support": 0.30,
}


def _read(path: str | Path | None) -> pd.DataFrame:
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    return pd.read_csv(path, low_memory=False)


def _join_keys(left: pd.DataFrame, right: pd.DataFrame) -> list[str]:
    keys = [key for key in ["drug_id", "target_id", "pdb_id"] if key in left.columns and key in right.columns]
    if len(keys) >= 2:
        return keys
    return [key for key in ["drug_id", "target_id"] if key in left.columns and key in right.columns]


def _first_nonempty(df: pd.DataFrame, columns: list[str]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in columns:
        if col not in df.columns:
            continue
        vals = df[col].fillna("").astype(str).str.strip()
        out = out.where(out.astype(str).str.len() > 0, vals)
    return out


def _with_moe_keys(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_moe_drug_key"] = _first_nonempty(
        out,
        ["generic_name", "mapped_drug_name", "display_name", "ligand_display", "drug_id"],
    ).str.lower()
    out["_moe_target_key"] = _first_nonempty(
        out,
        ["target_gene", "target_id", "target_uniprot"],
    ).str.upper()
    return out


def _deduplicate_for_keys(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    if not keys:
        return df
    score_cols = [
        col
        for col in [
            "ml_prediction_score",
            "source_calibrated_prediction_score",
            "raw_prediction_score",
            "binding_expert_score",
            "banana_atlas_blend_score",
            "banana_binding_probability",
            "mechanism_graph_score",
            "target_adr_evidence",
            "pathway_evidence",
        ]
        if col in df.columns
    ]
    if not score_cols:
        return df.drop_duplicates(keys, keep="first")
    work = df.copy()
    work["_moe_sort_score"] = pd.to_numeric(work[score_cols[0]], errors="coerce")
    work = work.sort_values([*keys, "_moe_sort_score"], ascending=[*([True] * len(keys)), False])
    return work.drop(columns=["_moe_sort_score"]).drop_duplicates(keys, keep="first")


def _normalize(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    lo = values.min(skipna=True)
    hi = values.max(skipna=True)
    if pd.isna(lo) or pd.isna(hi) or hi == lo:
        return pd.Series(pd.NA, index=series.index, dtype="Float64")
    out = (values - lo) / (hi - lo)
    if not higher_is_better:
        out = 1.0 - out
    return out.astype("Float64")


def _first_available(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series(pd.NA, index=df.index, dtype="Float64")
    for col in cols:
        if col not in df.columns:
            continue
        vals = pd.to_numeric(df[col], errors="coerce").astype("Float64")
        out = out.where(out.notna(), vals)
    return out


def _merge_expert(base: pd.DataFrame, expert: pd.DataFrame, prefix: str) -> pd.DataFrame:
    if expert.empty:
        base[f"{prefix}_expert_missing_reason"] = "missing_expert_table"
        return base
    base = _with_moe_keys(base)
    expert = _with_moe_keys(expert)
    keys = _join_keys(base, expert)
    if not keys:
        base[f"{prefix}_expert_missing_reason"] = "missing_join_keys"
        return base
    expert = _deduplicate_for_keys(expert, keys)
    rename = {col: f"{prefix}_{col}" for col in expert.columns if col not in keys}
    merged = base.merge(expert.rename(columns=rename), on=keys, how="left")
    marker_cols = [col for col in merged.columns if col.startswith(f"{prefix}_")]
    if marker_cols and not merged[marker_cols].notna().any(axis=1).any():
        alt_keys = ["_moe_drug_key", "_moe_target_key"]
        expert = _deduplicate_for_keys(expert, alt_keys)
        rename = {col: f"{prefix}_{col}" for col in expert.columns if col not in alt_keys}
        merged = base.merge(expert.rename(columns=rename), on=alt_keys, how="left")
        marker_cols = [col for col in merged.columns if col.startswith(f"{prefix}_")]
    if marker_cols:
        has_any = merged[marker_cols].notna().any(axis=1)
        merged[f"{prefix}_expert_missing_reason"] = ""
        merged.loc[~has_any, f"{prefix}_expert_missing_reason"] = "no_expert_row_match"
    return merged


def build_moe_adr_score(
    pair_table_path: str | Path,
    out_path: str | Path,
    *,
    bioactivity_predictions: str | Path | None = None,
    binding_expert: str | Path | None = None,
    exposure_predictions: str | Path | None = None,
    tissue_predictions: str | Path | None = None,
    mechanism_scores: str | Path | None = None,
    weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Combine target-specific experts into one ADR-implication evidence score.

    This is a mixture of available evidence, not a trained clinical ADR model.
    Missing experts are ignored and weights are rescaled over available scores.
    """

    base = _read(pair_table_path)
    if base.empty:
        raise ValueError(f"pair table is empty or missing: {pair_table_path}")
    keep = [
        col
        for col in [
            "drug_id",
            "target_id",
            "pdb_id",
            "generic_name",
            "display_name",
            "ligand_display",
            "target_gene",
            "target_uniprot",
            "atlas_score",
            "consensus_score",
        ]
        if col in base.columns
    ]
    out = base[keep].copy()
    out = _merge_expert(out, _read(bioactivity_predictions), "bioactivity")
    out = _merge_expert(out, _read(binding_expert), "binding")
    out = _merge_expert(out, _read(exposure_predictions), "exposure")
    out = _merge_expert(out, _read(tissue_predictions), "tissue")
    out = _merge_expert(out, _read(mechanism_scores), "mechanism")

    out["bioactivity_expert_score"] = _first_available(
        out,
        [
            "bioactivity_ml_prediction_score",
            "bioactivity_source_calibrated_prediction_score",
            "bioactivity_raw_prediction_score",
        ],
    )
    out["binding_activity_expert_score"] = _first_available(
        out,
        [
            "binding_binding_expert_score",
            "binding_banana_atlas_blend_score",
            "binding_banana_binding_probability",
            "bioactivity_ml_prediction_score",
            "bioactivity_source_calibrated_prediction_score",
            "bioactivity_raw_prediction_score",
        ],
    )
    out["exposure_relevance_expert_score"] = _first_available(
        out,
        [
            "exposure_ml_prediction_score",
            "exposure_source_calibrated_prediction_score",
            "exposure_raw_prediction_score",
        ],
    )
    if out["bioactivity_expert_score"].isna().all() and "atlas_score" in out.columns:
        out["bioactivity_expert_score"] = _normalize(out["atlas_score"])
        out["bioactivity_expert_missing_reason"] = out.get(
            "bioactivity_expert_missing_reason", "atlas_score_fallback"
        )
    if out["binding_activity_expert_score"].isna().all():
        out["binding_activity_expert_score"] = out["bioactivity_expert_score"]
    if out["binding_activity_expert_score"].isna().all() and "atlas_score" in out.columns:
        out["binding_activity_expert_score"] = _normalize(out["atlas_score"])
        out["binding_expert_missing_reason"] = out.get(
            "binding_expert_missing_reason", "atlas_score_fallback"
        )
    out["tissue_site_relevance_expert_score"] = _first_available(
        out,
        [
            "tissue_ml_prediction_score",
            "tissue_site_relevance_score",
            "tissue_site_safety_risk_score",
            "tissue_expression_presence_score",
            "tissue_target_tissue_expression_score",
        ],
    )
    if out["tissue_site_relevance_expert_score"].notna().any():
        out["tissue_site_relevance_expert_score"] = _normalize(out["tissue_site_relevance_expert_score"])
    out["mechanism_adr_support_expert_score"] = _first_available(
        out,
        [
            "mechanism_ml_prediction_score",
            "mechanism_mechanism_graph_score",
            "mechanism_target_adr_evidence",
            "mechanism_pathway_evidence",
        ],
    )
    if out["mechanism_adr_support_expert_score"].notna().any():
        out["mechanism_adr_support_expert_score"] = _normalize(out["mechanism_adr_support_expert_score"])

    weights = weights or DEFAULT_WEIGHTS
    expert_cols = {
        "binding_activity": "binding_activity_expert_score",
        "exposure_relevance": "exposure_relevance_expert_score",
        "tissue_site_relevance": "tissue_site_relevance_expert_score",
        "mechanism_adr_support": "mechanism_adr_support_expert_score",
    }
    numer = pd.Series(0.0, index=out.index)
    denom = pd.Series(0.0, index=out.index)
    availability: list[str] = []
    for name, col in expert_cols.items():
        w = float(weights.get(name, 0.0))
        vals = pd.to_numeric(out[col], errors="coerce")
        present = vals.notna()
        numer = numer.add(vals.fillna(0.0) * w)
        denom = denom.add(present.astype(float) * w)
        availability.append(name)
    out["moe_adr_implication_score"] = (numer / denom).where(denom.gt(0), pd.NA)
    missing_cols = list(expert_cols.values())
    out["moe_available_experts"] = out[missing_cols].notna().apply(
        lambda row: ";".join([name for name, present in zip(availability, row, strict=False) if present]),
        axis=1,
    )
    out["moe_missing_experts"] = out[missing_cols].notna().apply(
        lambda row: ";".join([name for name, present in zip(availability, row, strict=False) if not present]),
        axis=1,
    )
    out["moe_interpretation"] = (
        "Mixture-of-evidence score from separate binding/activity, exposure-relevance, "
        "tissue/site, and mechanism experts. Treat as prioritization evidence, not a clinical ADR probability."
    )

    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    manifest: dict[str, Any] = {
        "status": "written",
        "path": str(path),
        "rows": int(len(out)),
        "weights": weights,
        "n_with_bioactivity_expert": int(out["bioactivity_expert_score"].notna().sum()),
        "n_with_binding_activity_expert": int(out["binding_activity_expert_score"].notna().sum()),
        "n_with_exposure_expert": int(out["exposure_relevance_expert_score"].notna().sum()),
        "n_with_tissue_expert": int(out["tissue_site_relevance_expert_score"].notna().sum()),
        "n_with_mechanism_expert": int(out["mechanism_adr_support_expert_score"].notna().sum()),
        "n_with_any_moe_score": int(out["moe_adr_implication_score"].notna().sum()),
        "warning": (
            "A publishable trained evidence ensemble needs calibrated held-out predictions from all four expert tasks. "
            "This file is the conservative evidence-combination layer."
        ),
    }
    path.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return out
