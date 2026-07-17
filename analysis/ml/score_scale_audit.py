from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.materialize_final_score import (
    ACCEPTED_FINAL_SCORE_SOURCES,
    canonical_atlas_score_source,
)


MIXED_SCALE_SOURCE_TOKENS = (
    "consensus_score_fallback",
    "rank_percentile",
    "raw_consensus",
)
Z_DEPENDENT_FEATURES = {
    "atlas_score",
    "z_selected",
    "atlas_binding_prior",
    "consensus_z_score",
    "z_vs_decoys_consensus",
    "z_vs_compare_run_consensus",
    "z_vs_decoys_blend",
    "final_score",
}
SAFE_BINDING_PRIOR_SOURCES = {
    "banana_consensus_blend",
    "consensus_score_fallback",
}
SAFE_FINAL_SCORE_SOURCES = ACCEPTED_FINAL_SCORE_SOURCES


def _source_mask(series: pd.Series, tokens: Sequence[str]) -> pd.Series:
    text = series.fillna("").astype(str).str.strip().str.lower()
    mask = pd.Series(False, index=series.index)
    for token in tokens:
        mask |= text.str.contains(token, regex=False)
    return mask


def audit_score_scale_frame(
    frame: pd.DataFrame,
    *,
    feature_names: Sequence[str] | None = None,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    selected = sorted(set(feature_names or []))
    selected_z_features = sorted(set(selected) & Z_DEPENDENT_FEATURES)
    raw_prior_sensitive_selected = bool(
        set(selected) & {"atlas_score", "atlas_binding_prior"}
    )
    selected_binding_prior_features = sorted(
        set(selected) & {"binding_expert_score", "banana_atlas_blend_score"}
    )
    z_source = frame.get(
        "z_selected_source",
        pd.Series("", index=frame.index, dtype="object"),
    )
    final_source = frame.get(
        "final_score_source",
        pd.Series("", index=frame.index, dtype="object"),
    )
    fallback_mask = _source_mask(z_source, MIXED_SCALE_SOURCE_TOKENS)
    fallback_mask |= _source_mask(final_source, MIXED_SCALE_SOURCE_TOKENS)
    missing_null_mask = _source_mask(z_source, ("missing_consensus_decoy_null",))
    missing_null_mask |= _source_mask(final_source, ("missing_consensus_decoy_null",))

    prior_source = frame.get(
        "atlas_binding_prior_source",
        pd.Series("", index=frame.index, dtype="object"),
    )
    raw_prior_mask = _source_mask(prior_source, ("rank_percentile", "raw_consensus"))
    binding_source = (
        frame.get(
            "binding_expert_source",
            pd.Series("", index=frame.index, dtype="object"),
        )
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
    )
    unsafe_binding_prior = binding_source.ne("") & ~binding_source.isin(
        SAFE_BINDING_PRIOR_SOURCES
    )

    final_values = pd.to_numeric(
        frame.get("final_score", pd.Series(pd.NA, index=frame.index)),
        errors="coerce",
    )
    normalized_final_source = (
        final_source.fillna("").astype(str).str.strip().str.lower()
    )
    unsafe_final_score = final_values.notna() & ~normalized_final_source.isin(
        SAFE_FINAL_SCORE_SOURCES
    )
    atlas_values = pd.to_numeric(
        frame.get("atlas_score", pd.Series(pd.NA, index=frame.index)),
        errors="coerce",
    )
    atlas_source = frame.get(
        "atlas_score_source_for_ml",
        pd.Series("", index=frame.index, dtype="object"),
    )
    normalized_atlas_source = (
        atlas_source.fillna("").astype(str).str.strip().str.lower()
    )
    canonical_atlas_source = normalized_atlas_source.map(
        canonical_atlas_score_source
    )
    unsafe_atlas_score = atlas_values.notna() & canonical_atlas_source.eq("")

    z_values = pd.Series(float("nan"), index=frame.index, dtype="float64")
    for column in (
        "z_vs_decoys_consensus",
        "z_vs_compare_run_consensus",
        "consensus_z_score",
        "z_selected",
    ):
        if column in frame.columns:
            z_values = z_values.fillna(pd.to_numeric(frame[column], errors="coerce"))
    consensus_null_n = pd.Series(float("nan"), index=frame.index, dtype="float64")
    for column in (
        "consensus_z_decoy_n",
        "z_vs_compare_run_consensus_decoy_n",
        "consensus_n_decoys",
    ):
        if column in frame.columns:
            consensus_null_n = consensus_null_n.fillna(
                pd.to_numeric(frame[column], errors="coerce")
            )
    consensus_null_sigma = pd.Series(float("nan"), index=frame.index, dtype="float64")
    for column in (
        "consensus_z_decoy_sigma",
        "z_vs_compare_run_consensus_decoy_sigma",
        "consensus_sigma_decoy",
    ):
        if column in frame.columns:
            consensus_null_sigma = consensus_null_sigma.fillna(
                pd.to_numeric(frame[column], errors="coerce")
            )
    blend_null_n = pd.Series(float("nan"), index=frame.index, dtype="float64")
    for column in (
        "blend_n_decoys",
        "z_vs_compare_run_scorch_composite_decoy_n",
    ):
        if column in frame.columns:
            blend_null_n = blend_null_n.fillna(
                pd.to_numeric(frame[column], errors="coerce")
            )
    blend_null_sigma = pd.Series(float("nan"), index=frame.index, dtype="float64")
    for column in (
        "blend_sigma_decoy",
        "z_vs_compare_run_scorch_composite_decoy_sigma",
    ):
        if column in frame.columns:
            blend_null_sigma = blend_null_sigma.fillna(
                pd.to_numeric(frame[column], errors="coerce")
            )
    invalid_consensus_null = (
        consensus_null_n.lt(2)
        | consensus_null_n.isna()
        | consensus_null_sigma.le(0)
        | consensus_null_sigma.isna()
    )
    invalid_blend_null = (
        blend_null_n.lt(2)
        | blend_null_n.isna()
        | blend_null_sigma.le(0)
        | blend_null_sigma.isna()
    )
    final_consensus = final_values.notna() & normalized_final_source.eq(
        "z_vs_decoys_consensus"
    )
    final_blend = final_values.notna() & normalized_final_source.eq("z_vs_decoys_blend")
    missing_null_mask |= final_consensus & invalid_consensus_null
    missing_null_mask |= final_blend & invalid_blend_null
    consensus_feature_selected = bool(
        set(selected)
        & {
            "z_vs_decoys_consensus",
            "z_vs_compare_run_consensus",
            "consensus_z_score",
            "z_selected",
        }
    )
    if consensus_feature_selected:
        missing_null_mask |= z_values.notna() & invalid_consensus_null
    else:
        missing_null_mask |= (
            z_values.notna() & final_values.isna() & invalid_consensus_null
        )

    raw_consensus = pd.to_numeric(
        frame.get(
            "consensus_score_raw",
            frame.get("consensus_score", pd.Series(pd.NA, index=frame.index)),
        ),
        errors="coerce",
    )
    missing_z_with_raw = z_values.isna() & raw_consensus.notna()

    source_rows: list[dict[str, Any]] = []
    for source_column, source_series in (
        ("final_score_source", final_source),
        ("z_selected_source", z_source),
        ("atlas_score_source_for_ml", atlas_source),
    ):
        if source_column not in frame.columns:
            continue
        counts = source_series.fillna("missing").astype(str).value_counts(dropna=False)
        for value, count in counts.items():
            source_rows.append(
                {
                    "source_column": source_column,
                    "source_value": value,
                    "n_rows": int(count),
                    "fraction": float(count / len(frame)) if len(frame) else 0.0,
                }
            )
    if "atlas_binding_prior_source" in frame.columns:
        counts = prior_source.fillna("missing").astype(str).value_counts(dropna=False)
        for value, count in counts.items():
            source_rows.append(
                {
                    "source_column": "atlas_binding_prior_source",
                    "source_value": value,
                    "n_rows": int(count),
                    "fraction": float(count / len(frame)) if len(frame) else 0.0,
                }
            )

    null_quality_cols = [
        col
        for col in (
            "pdb_id",
            "blend_n_decoys",
            "blend_mu_decoy",
            "blend_sigma_decoy",
            "z_vs_compare_run_scorch_composite_decoy_n",
            "z_vs_compare_run_scorch_composite_decoy_mu",
            "z_vs_compare_run_scorch_composite_decoy_sigma",
            "z_vs_compare_run_consensus_decoy_n",
            "z_vs_compare_run_consensus_decoy_mu",
            "z_vs_compare_run_consensus_decoy_sigma",
            "z_vs_compare_run_consensus_null_sha256",
            "z_vs_compare_run_consensus_dud_sha256",
            "consensus_n_decoys",
            "consensus_mu_decoy",
            "consensus_sigma_decoy",
            "consensus_z_decoy_n",
            "consensus_z_decoy_unique",
            "consensus_z_decoy_zero_fraction",
            "consensus_z_decoy_mu",
            "consensus_z_decoy_sigma",
            "consensus_z_score_source",
        )
        if col in frame.columns
    ]
    null_quality = pd.DataFrame()
    if "pdb_id" in null_quality_cols:
        null_quality = (
            frame[null_quality_cols]
            .drop_duplicates("pdb_id", keep="first")
            .sort_values("pdb_id")
            .reset_index(drop=True)
        )

    flags: list[str] = []
    low_resolution_pdbs: list[str] = []
    if not null_quality.empty and {
        "consensus_z_decoy_n",
        "consensus_z_decoy_unique",
        "consensus_z_decoy_zero_fraction",
    }.issubset(null_quality.columns):
        n_decoys = pd.to_numeric(null_quality["consensus_z_decoy_n"], errors="coerce")
        n_unique = pd.to_numeric(
            null_quality["consensus_z_decoy_unique"],
            errors="coerce",
        )
        zero_fraction = pd.to_numeric(
            null_quality["consensus_z_decoy_zero_fraction"],
            errors="coerce",
        )
        low_resolution = n_unique.div(n_decoys).lt(0.5) | zero_fraction.gt(0.5)
        low_resolution_pdbs = (
            null_quality.loc[low_resolution, "pdb_id"].astype(str).tolist()
        )
        if low_resolution_pdbs:
            flags.append(
                "low-resolution decoy consensus nulls: "
                + ", ".join(low_resolution_pdbs)
            )
    n_fallback = int(fallback_mask.sum())
    n_raw_prior = int(raw_prior_mask.sum())
    n_missing_null = int(missing_null_mask.sum())
    n_missing_z_with_raw = int(missing_z_with_raw.sum())
    n_unsafe_final_score = int(unsafe_final_score.sum())
    n_unsafe_atlas_score = int(unsafe_atlas_score.sum())
    if n_fallback:
        flags.append(
            f"{n_fallback} rows place a raw consensus percentile on a z-score provenance path"
        )
    if n_raw_prior:
        flags.append(
            f"{n_raw_prior} rows use a raw consensus percentile as atlas_binding_prior"
        )
    if n_missing_null:
        flags.append(f"{n_missing_null} rows lack a same-run decoy null")
    if selected_z_features and n_fallback:
        flags.append(
            "selected z-dependent features contain raw-percentile fallback values: "
            + ", ".join(selected_z_features)
        )
    if selected_z_features and n_missing_z_with_raw:
        flags.append(
            f"{n_missing_z_with_raw} rows have raw consensus but no standardized consensus z-score"
        )
    if "final_score" in selected and n_unsafe_final_score:
        flags.append(
            f"{n_unsafe_final_score} finite final_score rows lack canonical z-vs-decoys provenance"
        )
    if "atlas_score" in selected and n_unsafe_atlas_score:
        flags.append(
            f"{n_unsafe_atlas_score} finite atlas_score rows lack accepted row provenance"
        )
    n_unsafe_binding_prior = int(unsafe_binding_prior.sum())
    if selected_binding_prior_features and n_unsafe_binding_prior:
        flags.append(
            "selected frozen binding-prior features use incompatible or legacy "
            "score provenance: " + ", ".join(selected_binding_prior_features)
        )

    z_scale_failure = bool(
        selected_z_features and (n_fallback or n_missing_z_with_raw or n_missing_null)
    )
    raw_prior_failure = bool(raw_prior_sensitive_selected and n_raw_prior)
    binding_prior_failure = bool(
        selected_binding_prior_features and n_unsafe_binding_prior
    )
    final_score_failure = bool("final_score" in selected and n_unsafe_final_score)
    atlas_score_failure = bool("atlas_score" in selected and n_unsafe_atlas_score)
    status = "passed"
    if (
        z_scale_failure
        or raw_prior_failure
        or binding_prior_failure
        or final_score_failure
        or atlas_score_failure
    ):
        status = "failed"
    elif flags:
        status = "warning"
    result = {
        "status": status,
        "publication_blocker": status == "failed",
        "n_rows": int(len(frame)),
        "selected_features": selected,
        "selected_z_dependent_features": selected_z_features,
        "selected_binding_prior_features": selected_binding_prior_features,
        "n_raw_consensus_fallback": n_fallback,
        "n_raw_atlas_prior": n_raw_prior,
        "n_unsafe_binding_prior": n_unsafe_binding_prior,
        "n_unsafe_final_score": n_unsafe_final_score,
        "n_unsafe_atlas_score": n_unsafe_atlas_score,
        "n_missing_consensus_decoy_null": n_missing_null,
        "n_missing_z_with_raw_consensus": n_missing_z_with_raw,
        "n_low_resolution_decoy_null_pdbs": len(low_resolution_pdbs),
        "low_resolution_decoy_null_pdbs": low_resolution_pdbs,
        "flags": flags,
        "policy": (
            "Raw consensus percentiles may be retained as consensus_score_raw for audit, "
            "but must not populate z_selected, atlas_score, atlas_binding_prior, or final_score."
        ),
    }
    return result, pd.DataFrame(source_rows), null_quality


def write_score_scale_audit(
    dataset: str | Path | pd.DataFrame,
    out_dir: str | Path,
    *,
    feature_names: Sequence[str] | None = None,
) -> dict[str, Any]:
    frame = (
        dataset.copy()
        if isinstance(dataset, pd.DataFrame)
        else pd.read_csv(dataset, low_memory=False)
    )
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    result, source_counts, null_quality = audit_score_scale_frame(
        frame,
        feature_names=feature_names,
    )
    source_counts.to_csv(out / "score_scale_source_counts.csv", index=False)
    null_quality.to_csv(out / "consensus_decoy_null_quality.csv", index=False)
    (out / "score_scale_integrity.json").write_text(
        json.dumps(result, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return {
        **result,
        "outputs": {
            "integrity": str(out / "score_scale_integrity.json"),
            "source_counts": str(out / "score_scale_source_counts.csv"),
            "decoy_null_quality": str(out / "consensus_decoy_null_quality.csv"),
        },
    }
