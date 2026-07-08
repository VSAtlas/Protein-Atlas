from __future__ import annotations

from collections.abc import Iterable
import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.audit_utils import grouped_label_balance, label_balance, source_label_balance
from analysis.ml.labels import binary_label_series


DEFAULT_GROUP_COLS = [
    "target_family",
    "protein_class",
    "target_id",
    "label_source",
    "source_family",
    "upstream_source",
    "adr_site_group",
    "mechanism_panel_match_basis",
]


def _available(columns: Iterable[str], df: pd.DataFrame) -> list[str]:
    return [col for col in columns if col in df.columns]


def _balance_with_flags(
    balance: pd.DataFrame,
    *,
    group_cols: list[str],
    min_positive: int,
    min_negative: int,
    extreme_rate_low: float,
    extreme_rate_high: float,
) -> pd.DataFrame:
    if balance.empty:
        return balance
    out = balance.copy()
    if "positive" not in out.columns and "pos" in out.columns:
        out["positive"] = out["pos"]
    if "negative" not in out.columns and "neg" in out.columns:
        out["negative"] = out["neg"]
    out["needs_positive_rows"] = (min_positive - pd.to_numeric(out["positive"], errors="coerce").fillna(0)).clip(lower=0).astype(int)
    out["needs_negative_rows"] = (min_negative - pd.to_numeric(out["negative"], errors="coerce").fillna(0)).clip(lower=0).astype(int)
    rate = pd.to_numeric(out["positive_rate"], errors="coerce")
    out["is_all_positive"] = rate.eq(1.0)
    out["is_all_negative"] = rate.eq(0.0)
    out["is_extreme_positive_rate"] = rate.ge(extreme_rate_high) | rate.le(extreme_rate_low)
    out["fails_minimum_fold_counts"] = out["needs_positive_rows"].gt(0) | out["needs_negative_rows"].gt(0)
    out["audit_flag"] = "ok"
    out.loc[out["is_all_positive"], "audit_flag"] = "all_positive"
    out.loc[out["is_all_negative"], "audit_flag"] = "all_negative"
    out.loc[out["fails_minimum_fold_counts"], "audit_flag"] = "insufficient_fold_class_counts"
    out["group_key"] = out[group_cols].astype(str).agg("|".join, axis=1) if group_cols else "all"
    return out


def _group_balance(
    data: pd.DataFrame,
    *,
    label_col: str,
    group_cols: list[str],
    min_positive: int,
    min_negative: int,
    extreme_rate_low: float,
    extreme_rate_high: float,
) -> pd.DataFrame:
    available = _available(group_cols, data)
    if not available:
        return pd.DataFrame()
    balance = grouped_label_balance(data, label_col, available)
    return _balance_with_flags(
        balance,
        group_cols=available,
        min_positive=min_positive,
        min_negative=min_negative,
        extreme_rate_low=extreme_rate_low,
        extreme_rate_high=extreme_rate_high,
    )


def _pairwise_group_balance(
    data: pd.DataFrame,
    *,
    label_col: str,
    primary_groups: list[str],
    context_cols: list[str],
    min_positive: int,
    min_negative: int,
    extreme_rate_low: float,
    extreme_rate_high: float,
) -> pd.DataFrame:
    rows = []
    for group_col in primary_groups:
        for context_col in context_cols:
            if group_col == context_col or group_col not in data.columns or context_col not in data.columns:
                continue
            balance = _group_balance(
                data,
                label_col=label_col,
                group_cols=[group_col, context_col],
                min_positive=min_positive,
                min_negative=min_negative,
                extreme_rate_low=extreme_rate_low,
                extreme_rate_high=extreme_rate_high,
            )
            if balance.empty:
                continue
            balance.insert(0, "primary_group_col", group_col)
            balance.insert(1, "context_col", context_col)
            rows.append(balance)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _acquisition_targets(group_balance: pd.DataFrame, source_balance: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if not group_balance.empty:
        for row in group_balance.itertuples(index=False):
            positive = int(getattr(row, "positive", 0))
            negative = int(getattr(row, "negative", 0))
            if positive == 0 or getattr(row, "needs_positive_rows", 0):
                rows.append(
                    {
                        "priority": "add_positives",
                        "group_col": getattr(row, "primary_group_col", None) or "group",
                        "group_key": getattr(row, "group_key", ""),
                        "current_positive": positive,
                        "current_negative": negative,
                        "needed_rows": int(getattr(row, "needs_positive_rows", 0)),
                        "reason": "group has too few positives for grouped validation",
                    }
                )
            if negative == 0 or getattr(row, "needs_negative_rows", 0):
                rows.append(
                    {
                        "priority": "add_negatives",
                        "group_col": getattr(row, "primary_group_col", None) or "group",
                        "group_key": getattr(row, "group_key", ""),
                        "current_positive": positive,
                        "current_negative": negative,
                        "needed_rows": int(getattr(row, "needs_negative_rows", 0)),
                        "reason": "group has too few negatives for grouped validation",
                    }
                )
    if not source_balance.empty:
        source_counts = source_balance.groupby("group_key", dropna=False)["label_source"].nunique() if "label_source" in source_balance else pd.Series(dtype=int)
        for group_key, n_sources in source_counts.items():
            if int(n_sources) < 2:
                rows.append(
                    {
                        "priority": "add_independent_source",
                        "group_col": "group_key",
                        "group_key": group_key,
                        "current_positive": None,
                        "current_negative": None,
                        "needed_rows": None,
                        "reason": "group has labels from fewer than two source lineages",
                    }
                )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.drop_duplicates().sort_values(["priority", "group_key"], kind="stable")


def _write_markdown_summary(
    out_path: Path,
    *,
    manifest: dict[str, Any],
    group_balance: pd.DataFrame,
    acquisition: pd.DataFrame,
) -> None:
    lines = [
        "# Label Balance Audit",
        "",
        f"- Dataset: `{manifest['dataset_path']}`",
        f"- Label: `{manifest['label_col']}`",
        f"- Labelable rows: `{manifest['label_balance'].get('n_labeled')}`",
        f"- Positives: `{manifest['label_balance'].get('positive')}`",
        f"- Negatives: `{manifest['label_balance'].get('negative')}`",
        f"- Unknown/missing: `{manifest['label_balance'].get('unknown_or_missing')}`",
        "",
        "## Findings",
    ]
    for finding in manifest["findings"]:
        lines.append(f"- {finding}")
    if not manifest["findings"]:
        lines.append("- No severe label-balance findings under the configured thresholds.")
    lines.extend(["", "## Highest Priority Acquisition Targets"])
    if acquisition.empty:
        lines.append("- No acquisition targets emitted.")
    else:
        for row in acquisition.head(25).itertuples(index=False):
            lines.append(f"- `{row.priority}` for `{row.group_key}`: {row.reason}")
    if not group_balance.empty:
        lines.extend(["", "## Most Imbalanced Groups"])
        cols = [col for col in ["group_key", "n", "positive", "negative", "positive_rate", "audit_flag"] if col in group_balance]
        preview = group_balance.sort_values(["audit_flag", "n"], ascending=[True, False]).head(20)
        lines.append(preview[cols].to_markdown(index=False))
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def audit_label_balance(
    dataset_path: str | Path,
    *,
    label_col: str,
    out_dir: str | Path,
    group_cols: list[str] | None = None,
    source_cols: list[str] | None = None,
    context_cols: list[str] | None = None,
    min_positive: int = 10,
    min_negative: int = 10,
    extreme_rate_low: float = 0.05,
    extreme_rate_high: float = 0.95,
) -> dict[str, Any]:
    df = pd.read_csv(dataset_path, low_memory=False)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    labels = binary_label_series(df[label_col])
    data = df.loc[labels.notna()].copy()
    data[label_col] = labels.loc[labels.notna()].astype(int)

    group_cols = group_cols or DEFAULT_GROUP_COLS
    source_cols = source_cols or ["label_source", "source_family", "upstream_source"]
    context_cols = context_cols or ["label_source", "source_family", "upstream_source", "assay_type", "endpoint_type"]
    available_groups = _available(group_cols, data)
    available_sources = _available(source_cols, data)

    grouped_frames = []
    for group_col in available_groups:
        grouped_frames.append(
            _group_balance(
                data,
                label_col=label_col,
                group_cols=[group_col],
                min_positive=min_positive,
                min_negative=min_negative,
                extreme_rate_low=extreme_rate_low,
                extreme_rate_high=extreme_rate_high,
            ).assign(primary_group_col=group_col)
        )
    group_balance = pd.concat(grouped_frames, ignore_index=True) if grouped_frames else pd.DataFrame()
    group_balance.to_csv(out / "label_balance_by_group.csv", index=False)

    source_frames = []
    for source_col in available_sources:
        source_frames.append(source_label_balance(data, label_col, source_col).rename(columns={source_col: "source_value"}).assign(source_col=source_col))
    source_summary = pd.concat(source_frames, ignore_index=True) if source_frames else pd.DataFrame()
    source_summary.to_csv(out / "label_balance_by_source.csv", index=False)

    pairwise = _pairwise_group_balance(
        data,
        label_col=label_col,
        primary_groups=[col for col in ["target_family", "protein_class", "target_id", "adr_site_group"] if col in available_groups],
        context_cols=_available(context_cols, data),
        min_positive=min_positive,
        min_negative=min_negative,
        extreme_rate_low=extreme_rate_low,
        extreme_rate_high=extreme_rate_high,
    )
    pairwise.to_csv(out / "label_balance_by_group_and_context.csv", index=False)
    acquisition = _acquisition_targets(group_balance, pairwise)
    acquisition.to_csv(out / "label_acquisition_targets.csv", index=False)

    findings = []
    if not group_balance.empty:
        n_all_positive = int(group_balance["is_all_positive"].sum())
        n_all_negative = int(group_balance["is_all_negative"].sum())
        n_failed = int(group_balance["fails_minimum_fold_counts"].sum())
        if n_all_positive:
            findings.append(f"{n_all_positive} groups are all-positive.")
        if n_all_negative:
            findings.append(f"{n_all_negative} groups are all-negative.")
        if n_failed:
            findings.append(f"{n_failed} groups fail min_positive={min_positive} or min_negative={min_negative}.")
    if not source_summary.empty and source_summary["positive_rate"].dropna().nunique() > 1:
        rate = pd.to_numeric(source_summary["positive_rate"], errors="coerce").dropna()
        if not rate.empty and rate.max() - rate.min() >= 0.50:
            findings.append("label prevalence differs by at least 50 percentage points across sources.")
    manifest = {
        "dataset_path": str(dataset_path),
        "label_col": label_col,
        "label_balance": label_balance(df, label_col),
        "group_cols": available_groups,
        "source_cols": available_sources,
        "context_cols": _available(context_cols, data),
        "min_positive": int(min_positive),
        "min_negative": int(min_negative),
        "extreme_rate_low": float(extreme_rate_low),
        "extreme_rate_high": float(extreme_rate_high),
        "findings": findings,
        "outputs": {
            "label_balance_by_group": str(out / "label_balance_by_group.csv"),
            "label_balance_by_source": str(out / "label_balance_by_source.csv"),
            "label_balance_by_group_and_context": str(out / "label_balance_by_group_and_context.csv"),
            "label_acquisition_targets": str(out / "label_acquisition_targets.csv"),
            "summary_md": str(out / "label_balance_audit.md"),
        },
    }
    (out / "label_balance_audit_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_markdown_summary(out / "label_balance_audit.md", manifest=manifest, group_balance=group_balance, acquisition=acquisition)
    return manifest
