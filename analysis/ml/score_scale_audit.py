from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd

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
}
SAFE_BINDING_PRIOR_SOURCES = {
    "banana_consensus_blend",
    "consensus_score_fallback",
}


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
    selected_binding_prior_features = sorted(
        set(selected) & {"binding_expert_score", "banana_atlas_blend_score"}
    )
    source = frame.get(
        "z_selected_source",
        pd.Series("", index=frame.index, dtype="object"),
    )
    fallback_mask = _source_mask(source, MIXED_SCALE_SOURCE_TOKENS)
    missing_null_mask = _source_mask(source, ("missing_consensus_decoy_null",))

    prior_source = frame.get(
        "atlas_binding_prior_source",
        pd.Series("", index=frame.index, dtype="object"),
    )
    raw_prior_mask = _source_mask(prior_source, ("rank_percentile", "raw_consensus"))
    binding_source = frame.get(
        "binding_expert_source",
        pd.Series("", index=frame.index, dtype="object"),
    ).fillna("").astype(str).str.strip().str.lower()
    unsafe_binding_prior = (
        binding_source.ne("")
        & ~binding_source.isin(SAFE_BINDING_PRIOR_SOURCES)
    )

    z_values = pd.to_numeric(
        frame.get("z_selected", pd.Series(pd.NA, index=frame.index)),
        errors="coerce",
    )
    raw_consensus = pd.to_numeric(
        frame.get("consensus_score_raw", frame.get("consensus_score", pd.Series(pd.NA, index=frame.index))),
        errors="coerce",
    )
    missing_z_with_raw = z_values.isna() & raw_consensus.notna()

    source_rows: list[dict[str, Any]] = []
    if "z_selected_source" in frame.columns:
        counts = source.fillna("missing").astype(str).value_counts(dropna=False)
        for value, count in counts.items():
            source_rows.append(
                {
                    "source_column": "z_selected_source",
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
    if selected_z_features and (n_fallback or n_raw_prior):
        flags.append(
            "selected z-dependent features contain raw-percentile fallback values: "
            + ", ".join(selected_z_features)
        )
    n_unsafe_binding_prior = int(unsafe_binding_prior.sum())
    if selected_binding_prior_features and n_unsafe_binding_prior:
        flags.append(
            "selected frozen binding-prior features use incompatible or legacy "
            "score provenance: " + ", ".join(selected_binding_prior_features)
        )

    status = "passed"
    if (
        selected_z_features
        and (n_fallback or n_raw_prior)
        or selected_binding_prior_features
        and n_unsafe_binding_prior
    ):
        status = "failed"
    elif flags or n_missing_z_with_raw:
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
        "n_missing_consensus_decoy_null": n_missing_null,
        "n_missing_z_with_raw_consensus": n_missing_z_with_raw,
        "n_low_resolution_decoy_null_pdbs": len(low_resolution_pdbs),
        "low_resolution_decoy_null_pdbs": low_resolution_pdbs,
        "flags": flags,
        "policy": (
            "Raw consensus percentiles may be retained as consensus_score_raw for audit, "
            "but must not populate z_selected, atlas_score, or atlas_binding_prior."
        ),
    }
    return result, pd.DataFrame(source_rows), null_quality


def write_score_scale_audit(
    dataset: str | Path | pd.DataFrame,
    out_dir: str | Path,
    *,
    feature_names: Sequence[str] | None = None,
) -> dict[str, Any]:
    frame = dataset.copy() if isinstance(dataset, pd.DataFrame) else pd.read_csv(dataset, low_memory=False)
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
