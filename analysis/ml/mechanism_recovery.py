from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd

DEFAULT_SCORE_COLS = [
    "ml_prediction_score",
    "mechanism_graph_score",
    "site_relevance_score_by_adr",
    "site_relevance_score",
    "binding_expert_score",
    "consensus_score",
]


def _label_state(df: pd.DataFrame, label_col: str) -> pd.Series:
    label = pd.to_numeric(df.get(label_col, pd.Series(pd.NA, index=df.index)), errors="coerce")
    state = pd.Series("unknown", index=df.index, dtype="object")
    state = state.mask(label.eq(1), "positive")
    state = state.mask(label.eq(0), "negative")
    if "mechanism_pu_state" in df.columns:
        excluded = df["mechanism_pu_state"].fillna("").astype(str).str.contains("excluded|conflict|ambiguous", case=False)
        state = state.mask(excluded, "excluded")
    return state


def _ef(pos_top: int, total_top: int, pos_total: int, total_labelable: int) -> float | None:
    if total_top <= 0 or pos_total <= 0 or total_labelable <= 0:
        return None
    background = pos_total / total_labelable
    if background <= 0:
        return None
    return (pos_top / total_top) / background


def _topk_for_group(
    group: pd.DataFrame,
    *,
    label_col: str,
    score_col: str,
    group_name: str,
    group_value: str,
    ks: list[int],
) -> list[dict[str, Any]]:
    ranked = group.dropna(subset=[score_col]).copy()
    if ranked.empty:
        return []
    ranked["_label_state"] = _label_state(ranked, label_col)
    ranked = ranked[~ranked["_label_state"].eq("excluded")].sort_values(score_col, ascending=False)
    if ranked.empty:
        return []
    labelable = ranked[ranked["_label_state"].isin(["positive", "negative"])]
    pos_total = int(labelable["_label_state"].eq("positive").sum())
    rows: list[dict[str, Any]] = []
    scopes = [("all_rows", ranked), ("labelable_only", labelable)]
    for ranking_scope, scope_ranked in scopes:
        if scope_ranked.empty:
            continue
        for k in ks:
            top = scope_ranked.head(min(k, len(scope_ranked)))
            top_labelable = top[top["_label_state"].isin(["positive", "negative"])]
            pos_top = int(top["_label_state"].eq("positive").sum())
            neg_top = int(top["_label_state"].eq("negative").sum())
            unknown_top = int(top["_label_state"].eq("unknown").sum())
            denom_precision = pos_top + neg_top
            known_precision = float(pos_top / denom_precision) if denom_precision else math.nan
            positive_recall = float(pos_top / pos_total) if pos_total else math.nan
            known_f1 = (
                2 * known_precision * positive_recall / (known_precision + positive_recall)
                if pd.notna(known_precision) and pd.notna(positive_recall) and known_precision + positive_recall > 0
                else math.nan
            )
            rows.append(
                {
                    "group_name": group_name,
                    "group_value": group_value,
                    "score_col": score_col,
                    "ranking_scope": ranking_scope,
                    "k": int(k),
                    "n_ranked": int(len(ranked)),
                    "n_labelable": int(len(labelable)),
                    "n_positive_total": pos_total,
                    "n_top": int(len(top)),
                    "n_top_labelable": int(len(top_labelable)),
                    "positive_hits_at_k": pos_top,
                    "reliable_negative_hits_at_k": neg_top,
                    "unknown_at_k": unknown_top,
                    "known_precision_at_k": known_precision,
                    "positive_recall_at_k": positive_recall,
                    "known_f1_at_k": known_f1,
                    "unknown_fraction_at_k": float(unknown_top / len(top)) if len(top) else math.nan,
                    "enrichment_at_k": _ef(pos_top, max(1, len(top_labelable)), pos_total, len(labelable)),
                }
            )
    return rows

def mechanism_topk_recovery(
    df: pd.DataFrame,
    *,
    label_col: str = "mechanism_pu_label",
    score_cols: list[str] | None = None,
    group_cols: list[str] | None = None,
    ks: list[int] | None = None,
) -> pd.DataFrame:
    scores = [col for col in (score_cols or DEFAULT_SCORE_COLS) if col in df.columns]
    if not scores:
        return pd.DataFrame()
    groups = group_cols or ["target_id", "pdb_id", "adr_site_group", "all"]
    top_ks = ks or [5, 10, 20, 50]
    rows: list[dict[str, Any]] = []
    for score_col in scores:
        rows.extend(
            _topk_for_group(
                df,
                label_col=label_col,
                score_col=score_col,
                group_name="all",
                group_value="all",
                ks=top_ks,
            )
        )
        for group_col in groups:
            if group_col == "all" or group_col not in df.columns:
                continue
            for value, group in df.groupby(group_col, dropna=False):
                value_text = "" if pd.isna(value) else str(value)
                if not value_text:
                    continue
                rows.extend(
                    _topk_for_group(
                        group,
                        label_col=label_col,
                        score_col=score_col,
                        group_name=group_col,
                        group_value=value_text,
                        ks=top_ks,
                    )
                )
    return pd.DataFrame(rows)


def write_mechanism_topk_recovery(
    table_path: str | Path,
    out_path: str | Path,
    *,
    label_col: str = "mechanism_pu_label",
    score_cols: list[str] | None = None,
    group_cols: list[str] | None = None,
    ks: list[int] | None = None,
) -> pd.DataFrame:
    df = pd.read_csv(table_path, low_memory=False)
    result = mechanism_topk_recovery(
        df,
        label_col=label_col,
        score_cols=score_cols,
        group_cols=group_cols,
        ks=ks,
    )
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, index=False)
    manifest = {
        "source_table": str(table_path),
        "out_path": str(out),
        "label_col": label_col,
        "score_cols": score_cols or DEFAULT_SCORE_COLS,
        "group_cols": group_cols or ["target_id", "pdb_id", "adr_site_group", "all"],
        "ks": ks or [5, 10, 20, 50],
        "rows": int(len(result)),
        "policy": "Top-K recovery treats unknown rows as unknown, not false negatives. ranking_scope=all_rows reports operational triage with unknown fractions; ranking_scope=labelable_only reports recovery only among known positives/reliable negatives.",
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return result
