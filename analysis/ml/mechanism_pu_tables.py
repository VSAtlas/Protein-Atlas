from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


POSITIVE_STATES = {"strict_positive", "projected_strict_positive", "strict_positive_mechanism_triad"}
NEGATIVE_HINTS = (
    "measured",
    "inactive",
    "reliable_negative",
    "negative_control",
    "omop_negative",
    "faers_nonsignal",
    "offsides_nonsignal",
)
EXCLUDED_HINTS = ("excluded", "ambiguous", "conflict", "conflicting", "benchmark_only")


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _lower(value: Any) -> str:
    return _clean(value).lower()


def _first_nonempty(row: pd.Series, columns: tuple[str, ...]) -> str:
    for col in columns:
        if col in row.index:
            text = _clean(row[col])
            if text:
                return text
    return ""


def _join_unique(values: pd.Series) -> str:
    seen: list[str] = []
    for value in values.dropna().astype(str):
        for token in value.split(";"):
            cleaned = _clean(token)
            if cleaned and cleaned not in seen:
                seen.append(cleaned)
    return ";".join(seen)


def _label_series(df: pd.DataFrame) -> pd.Series:
    out = pd.Series(pd.NA, index=df.index, dtype="Float64")
    for col in ("mechanism_ml_label_clean", "mechanism_ml_label", "four_state_ml_label", "mechanism_label"):
        if col not in df.columns:
            continue
        vals = pd.to_numeric(df[col], errors="coerce")
        vals = vals.where(vals.isin([-1, 0, 1]))
        out = out.where(out.notna(), vals)
    return out


def _state_for_row(row: pd.Series, label: float | None) -> str:
    status = _lower(row.get("mechanism_label_status", ""))
    negative_type = _lower(row.get("negative_evidence_type", ""))
    if label == 1:
        return "strict_positive"
    if label == 0:
        return "reliable_negative" if any(hint in f"{status} {negative_type}" for hint in NEGATIVE_HINTS) else "weak_negative"
    if label == -1 or any(hint in status for hint in EXCLUDED_HINTS):
        return "excluded_ambiguous_conflicting"
    return "unknown"


def _first_nonempty_frame(df: pd.DataFrame, columns: tuple[str, ...], *, upper: bool = False, lower: bool = False) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in columns:
        if col not in df.columns:
            continue
        values = df[col].map(_clean)
        out = out.where(out.astype(str).str.len().gt(0), values)
    if upper:
        out = out.str.upper()
    if lower:
        out = out.str.lower()
    return out


def _canonical_keys(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_pu_drug_key"] = _first_nonempty_frame(
        out,
        ("dedup_drug_key", "drug_id", "ligand_id", "drug_name", "compound_id"),
        lower=True,
    )
    out["_pu_target_key"] = _first_nonempty_frame(
        out,
        ("dedup_target_key", "target_gene", "gene_symbol", "target_id", "uniprot"),
        upper=True,
    )
    out["_pu_adr_key"] = _first_nonempty_frame(
        out,
        ("adr_id", "meddra_pt", "adr_pt", "adr_term", "mechanism_adr_key"),
        lower=True,
    )
    out["canonical_pair_key"] = (
        out["_pu_drug_key"].astype(str) + "|" + out["_pu_target_key"].astype(str) + "|" + out["_pu_adr_key"].astype(str)
    )
    return out


def _raw_states(df: pd.DataFrame, labels: pd.Series) -> pd.Series:
    status = df.get("mechanism_label_status", pd.Series("", index=df.index)).map(_lower)
    negative_type = df.get("negative_evidence_type", pd.Series("", index=df.index)).map(_lower)
    combined_negative = status + " " + negative_type
    states = pd.Series("unknown", index=df.index, dtype="object")
    positive = labels.eq(1).fillna(False).astype(bool)
    negative = labels.eq(0).fillna(False).astype(bool)
    excluded = labels.eq(-1).fillna(False).astype(bool) | status.map(
        lambda text: any(hint in text for hint in EXCLUDED_HINTS)
    )
    excluded = excluded.fillna(False).astype(bool)
    states = states.mask(positive, "strict_positive")
    reliable_negative = negative & combined_negative.map(lambda text: any(hint in text for hint in NEGATIVE_HINTS))
    reliable_negative = reliable_negative.fillna(False).astype(bool)
    states = states.mask(negative, "weak_negative")
    states = states.mask(reliable_negative, "reliable_negative")
    states = states.mask(excluded, "excluded_ambiguous_conflicting")
    return states


def _join_unique_by_group(df: pd.DataFrame, group_key: str, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(dtype="object")
    values = df[[group_key, col]].dropna(subset=[col]).copy()
    if values.empty:
        return pd.Series(dtype="object")
    values[col] = values[col].astype(str).str.split(";")
    values = values.explode(col)
    values[col] = values[col].map(_clean)
    values = values[values[col].astype(str).str.len().gt(0)]
    if values.empty:
        return pd.Series(dtype="object")
    values = values.drop_duplicates([group_key, col])
    return values.groupby(group_key, dropna=False)[col].agg(";".join)


def _collapse_usable(usable: pd.DataFrame) -> pd.DataFrame:
    priority = {
        "strict_positive": 0,
        "reliable_negative": 1,
        "weak_negative": 2,
        "unknown": 3,
        "excluded_ambiguous_conflicting": 4,
    }
    ranked = usable.assign(_collapse_priority=usable["_mechanism_pu_state_raw"].map(priority).fillna(9))
    representatives = ranked.sort_values(["canonical_pair_key", "_collapse_priority"]).drop_duplicates(
        "canonical_pair_key", keep="first"
    )
    label_raw = pd.to_numeric(usable["_mechanism_pu_label_raw"], errors="coerce")
    usable = usable.assign(
        _is_positive=label_raw.eq(1).fillna(False).astype(bool),
        _is_negative=label_raw.eq(0).fillna(False).astype(bool),
        _is_excluded=usable["_mechanism_pu_state_raw"].astype(str).eq("excluded_ambiguous_conflicting"),
        _is_reliable_negative=usable["_mechanism_pu_state_raw"].astype(str).eq("reliable_negative"),
        _is_weak_negative=usable["_mechanism_pu_state_raw"].astype(str).eq("weak_negative"),
    )
    flags = usable.groupby("canonical_pair_key", dropna=False).agg(
        mechanism_pu_source_rows=("canonical_pair_key", "size"),
        _has_positive=("_is_positive", "max"),
        _has_negative=("_is_negative", "max"),
        _has_excluded=("_is_excluded", "max"),
        _has_reliable_negative=("_is_reliable_negative", "max"),
        _has_weak_negative=("_is_weak_negative", "max"),
    )
    collapsed = representatives.drop(columns=["_collapse_priority"], errors="ignore").merge(
        flags,
        left_on="canonical_pair_key",
        right_index=True,
        how="left",
    )
    flag_cols = ["_has_positive", "_has_negative", "_has_excluded", "_has_reliable_negative", "_has_weak_negative"]
    for col in flag_cols:
        collapsed[col] = collapsed[col].fillna(False).astype(bool)
    conflict = collapsed["_has_positive"] & collapsed["_has_negative"]
    collapsed["mechanism_pu_state"] = "unknown"
    collapsed.loc[collapsed["_has_excluded"], "mechanism_pu_state"] = "excluded_ambiguous_conflicting"
    collapsed.loc[collapsed["_has_weak_negative"], "mechanism_pu_state"] = "weak_negative"
    collapsed.loc[collapsed["_has_reliable_negative"], "mechanism_pu_state"] = "reliable_negative"
    collapsed.loc[collapsed["_has_positive"], "mechanism_pu_state"] = "strict_positive"
    collapsed.loc[conflict, "mechanism_pu_state"] = "excluded_ambiguous_conflicting"
    collapsed["mechanism_pu_exclude_reason"] = ""
    source_excluded = collapsed["mechanism_pu_state"].eq("excluded_ambiguous_conflicting") & collapsed["_has_excluded"]
    collapsed.loc[source_excluded, "mechanism_pu_exclude_reason"] = "source_marked_excluded_or_ambiguous"
    collapsed.loc[conflict, "mechanism_pu_exclude_reason"] = "positive_negative_conflict"
    collapsed["mechanism_pu_label"] = pd.NA
    collapsed.loc[collapsed["mechanism_pu_state"].eq("strict_positive"), "mechanism_pu_label"] = 1
    collapsed.loc[collapsed["mechanism_pu_state"].isin(["reliable_negative", "weak_negative"]), "mechanism_pu_label"] = 0
    collapsed["mechanism_pu_training_role"] = collapsed["mechanism_pu_state"].map(
        {
            "strict_positive": "positive",
            "reliable_negative": "reliable_negative",
            "weak_negative": "weak_negative_sensitivity",
            "unknown": "unlabeled_background",
            "excluded_ambiguous_conflicting": "excluded",
        }
    )
    collapsed["_sample_weight"] = collapsed["mechanism_pu_state"].map(
        {
            "strict_positive": 1.0,
            "reliable_negative": 1.0,
            "weak_negative": 0.35,
            "unknown": 0.25,
            "excluded_ambiguous_conflicting": 0.0,
        }
    )
    for col in ("label_source", "source_family", "upstream_source", "mechanism_label_source", "negative_source"):
        joined = _join_unique_by_group(usable, "canonical_pair_key", col)
        if not joined.empty:
            collapsed = collapsed.merge(
                joined.rename(f"mechanism_pu_{col}_all"),
                left_on="canonical_pair_key",
                right_index=True,
                how="left",
            )
    return collapsed.drop(
        columns=[
            "_has_positive",
            "_has_negative",
            "_has_excluded",
            "_has_reliable_negative",
            "_has_weak_negative",
        ],
        errors="ignore",
    )

def build_mechanism_pu_table(
    mechanism_table_path: str | Path,
    out_path: str | Path,
    *,
    keep_weak_negatives_labeled: bool = False,
) -> pd.DataFrame:
    source = pd.read_csv(mechanism_table_path, low_memory=False)
    if source.empty:
        raise ValueError(f"mechanism table is empty: {mechanism_table_path}")
    work = _canonical_keys(source)
    labels = _label_series(work)
    work["_mechanism_pu_label_raw"] = labels
    work["_mechanism_pu_state_raw"] = _raw_states(work, labels)
    usable = work[work["canonical_pair_key"].astype(str).str.strip().ne("||")].copy()
    collapsed = _collapse_usable(usable).reset_index(drop=True)
    if not keep_weak_negatives_labeled:
        weak = collapsed["mechanism_pu_state"].eq("weak_negative")
        collapsed.loc[weak, "mechanism_pu_label"] = pd.NA
        collapsed.loc[weak, "mechanism_pu_training_role"] = "unlabeled_weak_negative_sensitivity"
    collapsed["mechanism_pu_label"] = pd.to_numeric(collapsed["mechanism_pu_label"], errors="coerce")
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    collapsed.to_csv(out, index=False)
    label = pd.to_numeric(collapsed["mechanism_pu_label"], errors="coerce")
    summary = {
        "source_table": str(mechanism_table_path),
        "out_path": str(out),
        "rows_in": int(len(source)),
        "rows_out_collapsed": int(len(collapsed)),
        "strict_positive_rows": int(collapsed["mechanism_pu_state"].eq("strict_positive").sum()),
        "reliable_negative_rows": int(collapsed["mechanism_pu_state"].eq("reliable_negative").sum()),
        "weak_negative_rows": int(collapsed["mechanism_pu_state"].eq("weak_negative").sum()),
        "unknown_rows": int(collapsed["mechanism_pu_state"].str.contains("unknown", na=False).sum()),
        "excluded_rows": int(collapsed["mechanism_pu_state"].str.contains("excluded", na=False).sum()),
        "positive_negative_conflict_rows": int(
            collapsed.get("mechanism_pu_exclude_reason", pd.Series("", index=collapsed.index))
            .eq("positive_negative_conflict")
            .sum()
        ),
        "source_marked_excluded_rows": int(
            collapsed.get("mechanism_pu_exclude_reason", pd.Series("", index=collapsed.index))
            .eq("source_marked_excluded_or_ambiguous")
            .sum()
        ),
        "labelable_rows": int(label.notna().sum()),
        "positive_rows": int(label.eq(1).sum()),
        "negative_rows": int(label.eq(0).sum()),
        "policy": "Strict positives, reliable negatives, unknowns, and excluded/conflicting rows are separated. Unknowns are available only as fold-local PU background, not global negatives.",
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return collapsed
