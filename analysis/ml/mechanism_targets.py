from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.source_tables import read_source_table
from analysis.graph.graph_scores import compute_pair_mechanism_scores
from analysis.ml.labels import binary_label_series, truthy_series


def build_mechanism_scores_from_edges(
    pair_table_path: str | Path,
    mechanism_edges_path: str | Path,
    out_path: str | Path,
) -> pd.DataFrame:
    """Build pair-level mechanism scores from a local mechanism edge table."""

    return compute_pair_mechanism_scores(pair_table_path, mechanism_edges_path, out_path)


def mechanism_label_status(df: pd.DataFrame) -> pd.Series:
    if "triad_complete" in df.columns:
        triad = pd.to_numeric(df["triad_complete"], errors="coerce").fillna(0).gt(0)
        return pd.Series("unknown_no_triad_evidence", index=df.index).mask(
            triad, "labeled_mechanism_triad"
        )
    if "mechanism_graph_score" in df.columns:
        score = pd.to_numeric(df["mechanism_graph_score"], errors="coerce")
        return pd.Series("unknown_no_mechanism_graph_score", index=df.index).mask(
            score.gt(0), "labeled_mechanism_evidence"
        )
    return pd.Series("unknown_no_mechanism_evidence_source", index=df.index)


def _base_positive_labels(df: pd.DataFrame) -> pd.Series:
    label = pd.Series(pd.NA, index=df.index, dtype="Int64")
    if "literature_supported_label" in df.columns:
        parsed = binary_label_series(df["literature_supported_label"])
        label = label.mask(parsed.eq(1), 1)
    if "triad_complete" in df.columns:
        label = label.mask(truthy_series(df["triad_complete"]), 1)
    elif "mechanism_graph_score" in df.columns:
        score = pd.to_numeric(df["mechanism_graph_score"], errors="coerce")
        label = label.mask(score.gt(0), 1)
    return label


def _merge_negative_scope(
    table: pd.DataFrame,
    evidence: pd.DataFrame,
    keys: list[str],
    suffix: str,
) -> pd.DataFrame:
    if not set(keys).issubset(table.columns) or not set(keys).issubset(evidence.columns):
        return table
    cols = [
        col
        for col in [
            *keys,
            "mechanism_label",
            "negative_evidence_type",
            "negative_source",
            "negative_confidence",
            "negative_label_status",
            "negative_selection_fold",
            "negative_selection_features",
        ]
        if col in evidence.columns
    ]
    neg = evidence[cols].dropna(subset=["mechanism_label"]).copy()
    if neg.empty:
        return table
    neg = neg.sort_values([*keys, "negative_confidence"], ascending=[*(True for _ in keys), False])
    neg = neg.drop_duplicates(keys, keep="first")
    return table.merge(neg, on=keys, how="left", suffixes=("", suffix))


def _merge_positive_scope(
    table: pd.DataFrame,
    evidence: pd.DataFrame,
    keys: list[str],
) -> pd.Series:
    if not set(keys).issubset(table.columns) or not set(keys).issubset(evidence.columns):
        return pd.Series(pd.NA, index=table.index, dtype="Int64")
    label = pd.to_numeric(evidence.get("mechanism_label", pd.Series(pd.NA, index=evidence.index)), errors="coerce")
    positives = evidence[label.eq(1)].dropna(subset=keys).copy()
    if positives.empty:
        return pd.Series(pd.NA, index=table.index, dtype="Int64")
    positives = positives.drop_duplicates(keys)
    marked = table[keys].merge(positives[keys].assign(_positive_mechanism_label=1), on=keys, how="left")
    return pd.to_numeric(marked["_positive_mechanism_label"], errors="coerce").astype("Int64")


def build_four_state_mechanism_labels(
    pair_table_path: str | Path,
    out_path: str | Path,
    *,
    mechanism_scores_path: str | Path | None = None,
    negative_evidence_path: str | Path | None = None,
    positive_evidence_path: str | Path | None = None,
) -> pd.DataFrame:
    """Build the four-state mechanism label table.

    Labels are 1 = strict positive triad, 0 = measured/reliable negative,
    NaN = unknown, and -1 = ambiguous/conflicting. Drug-ADR negatives only join
    when an ADR key is present in the pair/mechanism table.
    """

    base = pd.read_csv(pair_table_path, low_memory=False)
    if mechanism_scores_path is not None and Path(mechanism_scores_path).exists():
        scores = pd.read_csv(mechanism_scores_path, low_memory=False)
        if {"drug_id", "target_id"}.issubset(base.columns) and {"drug_id", "target_id"}.issubset(
            scores.columns
        ):
            score_cols = [
                col
                for col in scores.columns
                if col not in {"drug_id", "target_id"} and col not in base.columns
            ]
            base = base.merge(scores[["drug_id", "target_id", *score_cols]], on=["drug_id", "target_id"], how="left")
    labels = _base_positive_labels(base)
    if positive_evidence_path is not None and Path(positive_evidence_path).exists():
        positive = read_source_table(positive_evidence_path)
        drug_target_positive = positive[
            positive.get("evidence_scope", pd.Series("", index=positive.index)).astype(str).eq("drug_target")
        ]
        direct_positive = _merge_positive_scope(base, drug_target_positive, ["drug_id", "target_id"])
        labels = labels.mask(direct_positive.eq(1).fillna(False), 1)
    out = base.copy()
    if "mechanism_label" in out.columns:
        out = out.rename(columns={"mechanism_label": "_input_mechanism_label"})
    for col in (
        "negative_evidence_type",
        "negative_source",
        "negative_confidence",
        "negative_label_status",
        "negative_selection_fold",
        "negative_selection_features",
    ):
        if col not in out.columns:
            out[col] = pd.NA
    if negative_evidence_path is not None and Path(negative_evidence_path).exists():
        neg = pd.read_csv(negative_evidence_path, low_memory=False)
        neg_label = pd.to_numeric(neg.get("mechanism_label", pd.Series(pd.NA, index=neg.index)), errors="coerce")
        neg = neg[neg_label.isin([0, -1])].copy()
        neg["mechanism_label"] = neg_label.loc[neg.index].astype("Int64")
        if not neg.empty:
            drug_target = neg[neg.get("evidence_scope", "").astype(str).eq("drug_target")]
            out = _merge_negative_scope(out, drug_target, ["drug_id", "target_id"], "_drug_target")
            drug_adr = neg[neg.get("evidence_scope", "").astype(str).eq("drug_adr")]
            adr_key = "adr_id" if "adr_id" in out.columns else "adr_term" if "adr_term" in out.columns else None
            if adr_key and adr_key in drug_adr.columns:
                out = _merge_negative_scope(out, drug_adr, ["drug_id", adr_key], "_drug_adr")
    neg_label_cols = [
        col for col in out.columns if col == "mechanism_label" or col.startswith("mechanism_label_")
    ]
    neg_matrix = pd.DataFrame(index=out.index)
    for col in neg_label_cols:
        neg_matrix[col] = pd.to_numeric(out[col], errors="coerce")
    neg_values = pd.Series(pd.NA, index=out.index, dtype="Float64")
    if not neg_matrix.empty:
        any_ambiguous = neg_matrix.eq(-1).any(axis=1)
        any_negative = neg_matrix.eq(0).any(axis=1)
        neg_values = neg_values.mask(any_negative, 0)
        neg_values = neg_values.mask(any_ambiguous, -1)
    for stem in (
        "negative_evidence_type",
        "negative_source",
        "negative_confidence",
        "negative_label_status",
        "negative_selection_fold",
        "negative_selection_features",
    ):
        for col in [c for c in out.columns if c.startswith(f"{stem}_")]:
            fill_mask = out[stem].isna() & out[col].notna()
            if fill_mask.any():
                out.loc[fill_mask, stem] = out.loc[fill_mask, col]
    final = labels.copy()
    conflicts = labels.eq(1).fillna(False) & neg_values.isin([0, -1]).fillna(False)
    final = final.mask(neg_values.eq(0).fillna(False) & labels.isna(), 0)
    final = final.mask(neg_values.eq(-1).fillna(False), -1)
    final = final.mask(conflicts, -1)
    out["mechanism_label"] = final.astype("Int64")
    out["mechanism_ml_label"] = out["mechanism_label"].mask(out["mechanism_label"].eq(-1), pd.NA)
    status = pd.Series("unknown", index=out.index)
    status = status.mask(out["mechanism_label"].eq(1).fillna(False), "strict_positive_mechanism_triad")
    status = status.mask(out["mechanism_label"].eq(0).fillna(False), "measured_or_reliable_negative")
    status = status.mask(out["mechanism_label"].eq(-1).fillna(False), "excluded_ambiguous_or_conflicting")
    out["mechanism_label_status"] = status
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    return out
