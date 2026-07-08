from __future__ import annotations

import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.audit_utils import load_table, split_source_tokens
from analysis.ml.labels import binary_label_series


DEFAULT_LABELS = [
    "spd_binding_label",
    "spd_exposure_relevant",
    "combined_activity_ml_label",
    "external_four_state_ml_label",
    "mechanism_ml_label",
    "mechanism_ml_label_clean",
    "tissue_site_label",
]

DEFAULT_GROUP_COLS = [
    "target_family",
    "protein_class",
    "target_id",
    "pdb_id",
    "scaffold_key",
    "chemical_cluster",
    "ligand_chemotype",
    "drug_id",
    "label_source",
    "source_family",
    "upstream_source",
    "external_source_family",
    "external_upstream_source",
    "external_evidence_sources",
    "mechanism_label_source",
    "assay_type",
    "endpoint_type",
    "activity_type",
]

HOLDOUT_GROUP_CANDIDATES: dict[str, list[str]] = {
    "target_holdout": ["target_id", "pdb_id", "target_uniprot", "target_gene"],
    "target_family_holdout": ["target_family", "protein_class"],
    "scaffold_holdout": ["scaffold_key", "ligand_chemotype"],
    "chemical_cluster_holdout": ["chemical_cluster", "scaffold_key", "ligand_chemotype"],
    "source_holdout": [
        "external_source_family",
        "source_family",
        "label_source",
        "external_upstream_source",
        "upstream_source",
        "external_evidence_sources",
    ],
}

SOURCE_CALIBRATION_COLS = [
    "external_source_family",
    "source_family",
    "label_source",
    "external_upstream_source",
    "upstream_source",
    "external_evidence_sources",
]

SOURCE_LIKE_COLUMNS = {
    "external_evidence_sources",
    "external_parent_sources",
    "external_source_family",
    "external_upstream_source",
    "mechanism_label_source",
    "label_source",
    "source_family",
    "upstream_source",
}


def run_data_gap_report(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    labels: Sequence[str] | None = None,
    group_cols: Sequence[str] | None = None,
    min_group_labelable: int = 20,
    min_group_positives: int = 10,
    min_group_negatives: int = 10,
    min_holdout_test_positives: int = 10,
    min_holdout_test_negatives: int = 10,
    min_holdout_train_positives: int = 20,
    min_holdout_train_negatives: int = 20,
    min_calibration_positives: int = 20,
    min_calibration_negatives: int = 20,
    severe_positive_rate: float = 0.85,
    sparse_positive_rate: float = 0.05,
    max_unknown_fraction: float = 0.90,
) -> dict[str, Any]:
    """Write claim-oriented label/source/family/scaffold gap reports.

    The output is intentionally diagnostic. It should guide data collection,
    calibration, and split choice; it should not create predictive features.
    """

    dataset = Path(dataset_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = load_table(dataset)
    active_labels = [label for label in (labels or DEFAULT_LABELS) if label in df.columns]
    active_groups = [col for col in (group_cols or DEFAULT_GROUP_COLS) if col in df.columns]

    label_summary = _label_summary(df, active_labels)
    stratum_gaps = _stratum_gaps(
        df,
        active_labels,
        active_groups,
        min_group_labelable=min_group_labelable,
        min_group_positives=min_group_positives,
        min_group_negatives=min_group_negatives,
        severe_positive_rate=severe_positive_rate,
        sparse_positive_rate=sparse_positive_rate,
        max_unknown_fraction=max_unknown_fraction,
    )
    holdout_readiness = _holdout_readiness(
        df,
        active_labels,
        min_test_positives=min_holdout_test_positives,
        min_test_negatives=min_holdout_test_negatives,
        min_train_positives=min_holdout_train_positives,
        min_train_negatives=min_holdout_train_negatives,
    )
    source_calibration = _source_calibration_readiness(
        df,
        active_labels,
        min_positives=min_calibration_positives,
        min_negatives=min_calibration_negatives,
    )
    priorities = _data_collection_priorities(stratum_gaps, holdout_readiness, source_calibration)

    outputs = {
        "label_summary": out / "data_gap_label_summary.csv",
        "stratum_gaps": out / "ml_data_gap_report.csv",
        "holdout_readiness": out / "holdout_claim_readiness.csv",
        "source_calibration_readiness": out / "source_calibration_readiness.csv",
        "data_collection_priorities": out / "data_collection_priorities.csv",
    }
    label_summary.to_csv(outputs["label_summary"], index=False)
    stratum_gaps.to_csv(outputs["stratum_gaps"], index=False)
    holdout_readiness.to_csv(outputs["holdout_readiness"], index=False)
    source_calibration.to_csv(outputs["source_calibration_readiness"], index=False)
    priorities.to_csv(outputs["data_collection_priorities"], index=False)

    manifest: dict[str, Any] = {
        "dataset": str(dataset),
        "out_dir": str(out),
        "n_rows": int(len(df)),
        "labels": active_labels,
        "group_cols": active_groups,
        "thresholds": {
            "min_group_labelable": min_group_labelable,
            "min_group_positives": min_group_positives,
            "min_group_negatives": min_group_negatives,
            "min_holdout_test_positives": min_holdout_test_positives,
            "min_holdout_test_negatives": min_holdout_test_negatives,
            "min_holdout_train_positives": min_holdout_train_positives,
            "min_holdout_train_negatives": min_holdout_train_negatives,
            "min_calibration_positives": min_calibration_positives,
            "min_calibration_negatives": min_calibration_negatives,
            "severe_positive_rate": severe_positive_rate,
            "sparse_positive_rate": sparse_positive_rate,
            "max_unknown_fraction": max_unknown_fraction,
        },
        "policy": {
            "provenance_is_diagnostic_not_clean_feature": True,
            "blocked_holdout_means_no_claim_grade_generalization": True,
            "positive_unlabeled_labels_are_evidence_layers_unless_controls_exist": True,
        },
        "outputs": {key: str(path) for key, path in outputs.items()},
        "n_blocking_strata": int(stratum_gaps["status"].isin(["blocked", "diagnostic_only"]).sum())
        if "status" in stratum_gaps
        else 0,
        "n_blocked_holdouts": int(holdout_readiness["status"].eq("blocked").sum()) if "status" in holdout_readiness else 0,
    }
    (out / "data_gap_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    _write_markdown_summary(out / "data_gap_report.md", manifest, label_summary, priorities)
    return manifest


def _label_summary(df: pd.DataFrame, labels: Sequence[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for label in labels:
        series = pd.to_numeric(df[label], errors="coerce")
        binary = binary_label_series(df[label])
        positives = int(binary.eq(1).sum())
        negatives = int(binary.eq(0).sum())
        excluded = int(series.eq(-1).sum())
        unknown = int(series.isna().sum())
        labelable = positives + negatives
        status = "ok"
        if positives and not negatives and unknown:
            status = "positive_unlabeled_evidence_layer"
        elif positives == 0 and negatives == 0:
            status = "unlabeled"
        elif positives == 0 or negatives == 0:
            status = "single_class"
        rows.append(
            {
                "label": label,
                "n_rows": int(len(df)),
                "n_labelable": labelable,
                "n_positive": positives,
                "n_negative": negatives,
                "n_unknown": unknown,
                "n_excluded_or_conflicting": excluded,
                "positive_rate": positives / labelable if labelable else math.nan,
                "unknown_fraction": unknown / len(df) if len(df) else math.nan,
                "status": status,
                "recommended_use": _recommended_label_use(status),
            }
        )
    return pd.DataFrame(rows)


def _stratum_gaps(
    df: pd.DataFrame,
    labels: Sequence[str],
    group_cols: Sequence[str],
    *,
    min_group_labelable: int,
    min_group_positives: int,
    min_group_negatives: int,
    severe_positive_rate: float,
    sparse_positive_rate: float,
    max_unknown_fraction: float,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for label in labels:
        for group_col in group_cols:
            for group_value, group in _iter_groups(df, group_col):
                stats = _binary_stats(group[label], n_total=len(group))
                status, gap, action = _classify_group_gap(
                    stats,
                    min_group_labelable=min_group_labelable,
                    min_group_positives=min_group_positives,
                    min_group_negatives=min_group_negatives,
                    severe_positive_rate=severe_positive_rate,
                    sparse_positive_rate=sparse_positive_rate,
                    max_unknown_fraction=max_unknown_fraction,
                )
                rows.append(
                    {
                        "label": label,
                        "group_col": group_col,
                        "group_value": group_value,
                        **stats,
                        "status": status,
                        "gap_type": gap,
                        "recommended_action": action,
                    }
                )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return frame.sort_values(["status", "label", "group_col", "n_labelable"], ascending=[True, True, True, False])


def _holdout_readiness(
    df: pd.DataFrame,
    labels: Sequence[str],
    *,
    min_test_positives: int,
    min_test_negatives: int,
    min_train_positives: int,
    min_train_negatives: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for label in labels:
        total_stats = _binary_stats(df[label], n_total=len(df))
        for split_name, candidates in HOLDOUT_GROUP_CANDIDATES.items():
            group_col = _first_available(df, candidates)
            if group_col is None:
                rows.append(
                    {
                        "label": label,
                        "split": split_name,
                        "group_col": "",
                        "status": "blocked",
                        "reason": f"missing one of {candidates}",
                        **total_stats,
                    }
                )
                continue
            valid_groups = 0
            candidate_groups = 0
            blocker_examples: list[str] = []
            for group_value, group in _iter_groups(df, group_col):
                test_stats = _binary_stats(group[label], n_total=len(group))
                if test_stats["n_labelable"] == 0:
                    continue
                candidate_groups += 1
                train_stats = _subtract_stats(total_stats, test_stats)
                test_ok = (
                    test_stats["n_positive"] >= min_test_positives
                    and test_stats["n_negative"] >= min_test_negatives
                )
                train_ok = (
                    train_stats["n_positive"] >= min_train_positives
                    and train_stats["n_negative"] >= min_train_negatives
                )
                if test_ok and train_ok:
                    valid_groups += 1
                elif len(blocker_examples) < 5:
                    blocker_examples.append(
                        f"{group_value}: test_pos={test_stats['n_positive']} "
                        f"test_neg={test_stats['n_negative']} train_pos={train_stats['n_positive']} "
                        f"train_neg={train_stats['n_negative']}"
                    )
            if valid_groups >= 2:
                status = "ok"
                reason = "enough viable held-out groups for repeated/grouped evaluation"
            elif valid_groups == 1:
                status = "diagnostic_only"
                reason = "only one viable held-out group; use as case-study/diagnostic, not a broad claim"
            else:
                status = "blocked"
                reason = "no held-out group has enough positives and negatives"
            rows.append(
                {
                    "label": label,
                    "split": split_name,
                    "group_col": group_col,
                    "status": status,
                    "reason": reason,
                    "n_candidate_groups": candidate_groups,
                    "n_valid_groups": valid_groups,
                    "blocker_examples": "; ".join(blocker_examples),
                    **total_stats,
                }
            )
    return pd.DataFrame(rows)


def _source_calibration_readiness(
    df: pd.DataFrame,
    labels: Sequence[str],
    *,
    min_positives: int,
    min_negatives: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    source_cols = [col for col in SOURCE_CALIBRATION_COLS if col in df.columns]
    for label in labels:
        for source_col in source_cols:
            for source_value, group in _iter_groups(df, source_col):
                stats = _binary_stats(group[label], n_total=len(group))
                ready = stats["n_positive"] >= min_positives and stats["n_negative"] >= min_negatives
                status = "calibration_ready" if ready else "not_enough_two_class_rows"
                rows.append(
                    {
                        "label": label,
                        "source_col": source_col,
                        "source_value": source_value,
                        **stats,
                        "status": status,
                        "recommended_action": (
                            "fit source-specific sigmoid/isotonic calibrator on a held-out split"
                            if ready
                            else "report as subgroup only; do not fit a source-specific calibrator yet"
                        ),
                    }
                )
    return pd.DataFrame(rows)


def _data_collection_priorities(
    stratum_gaps: pd.DataFrame,
    holdout_readiness: pd.DataFrame,
    source_calibration: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if not stratum_gaps.empty:
        blockers = stratum_gaps.loc[stratum_gaps["status"].isin(["blocked", "diagnostic_only"])].copy()
        for _, row in blockers.head(100).iterrows():
            rows.append(
                {
                    "priority_type": "stratum_balance",
                    "label": row["label"],
                    "axis": row["group_col"],
                    "value": row["group_value"],
                    "status": row["status"],
                    "reason": row["gap_type"],
                    "recommended_action": row["recommended_action"],
                    "n_positive": row["n_positive"],
                    "n_negative": row["n_negative"],
                    "n_unknown": row["n_unknown"],
                }
            )
    if not holdout_readiness.empty:
        for _, row in holdout_readiness.loc[holdout_readiness["status"].ne("ok")].iterrows():
            rows.append(
                {
                    "priority_type": "holdout_claim",
                    "label": row["label"],
                    "axis": row["split"],
                    "value": row["group_col"],
                    "status": row["status"],
                    "reason": row["reason"],
                    "recommended_action": _holdout_action(str(row["split"]), str(row["group_col"])),
                    "n_positive": row["n_positive"],
                    "n_negative": row["n_negative"],
                    "n_unknown": row["n_unknown"],
                }
            )
    if not source_calibration.empty:
        no_cal = source_calibration.loc[source_calibration["status"].ne("calibration_ready")]
        for _, row in no_cal.head(50).iterrows():
            rows.append(
                {
                    "priority_type": "source_calibration",
                    "label": row["label"],
                    "axis": row["source_col"],
                    "value": row["source_value"],
                    "status": row["status"],
                    "reason": "source-specific calibration lacks two-class support",
                    "recommended_action": row["recommended_action"],
                    "n_positive": row["n_positive"],
                    "n_negative": row["n_negative"],
                    "n_unknown": row["n_unknown"],
                }
            )
    return pd.DataFrame(rows)


def _binary_stats(series: pd.Series, *, n_total: int) -> dict[str, Any]:
    labels = binary_label_series(series)
    numeric = pd.to_numeric(series, errors="coerce")
    positives = int(labels.eq(1).sum())
    negatives = int(labels.eq(0).sum())
    excluded = int(numeric.eq(-1).sum())
    unknown = int(n_total - positives - negatives - excluded)
    labelable = positives + negatives
    return {
        "n_rows": int(n_total),
        "n_labelable": labelable,
        "n_positive": positives,
        "n_negative": negatives,
        "n_unknown": unknown,
        "n_excluded_or_conflicting": excluded,
        "positive_rate": positives / labelable if labelable else math.nan,
        "unknown_fraction": unknown / n_total if n_total else math.nan,
    }


def _subtract_stats(total: dict[str, Any], part: dict[str, Any]) -> dict[str, Any]:
    n_rows = int(total["n_rows"]) - int(part["n_rows"])
    positives = int(total["n_positive"]) - int(part["n_positive"])
    negatives = int(total["n_negative"]) - int(part["n_negative"])
    unknown = int(total["n_unknown"]) - int(part["n_unknown"])
    excluded = int(total["n_excluded_or_conflicting"]) - int(part["n_excluded_or_conflicting"])
    labelable = positives + negatives
    return {
        "n_rows": n_rows,
        "n_labelable": labelable,
        "n_positive": positives,
        "n_negative": negatives,
        "n_unknown": unknown,
        "n_excluded_or_conflicting": excluded,
        "positive_rate": positives / labelable if labelable else math.nan,
        "unknown_fraction": unknown / n_rows if n_rows else math.nan,
    }


def _classify_group_gap(
    stats: dict[str, Any],
    *,
    min_group_labelable: int,
    min_group_positives: int,
    min_group_negatives: int,
    severe_positive_rate: float,
    sparse_positive_rate: float,
    max_unknown_fraction: float,
) -> tuple[str, str, str]:
    positives = int(stats["n_positive"])
    negatives = int(stats["n_negative"])
    labelable = int(stats["n_labelable"])
    unknown = int(stats["n_unknown"])
    unknown_fraction = float(stats["unknown_fraction"]) if not pd.isna(stats["unknown_fraction"]) else 1.0
    rate = float(stats["positive_rate"]) if labelable else math.nan
    if positives > 0 and negatives == 0 and unknown > 0:
        return (
            "diagnostic_only",
            "positive_unlabeled_no_controls",
            "keep as PU/evidence-ranking; add independent controls before supervised AUROC claims",
        )
    if labelable < min_group_labelable:
        return ("blocked", "too_few_labelable_rows", "add measured positives and negatives for this stratum")
    if positives < min_group_positives and negatives < min_group_negatives:
        return ("blocked", "too_few_two_class_examples", "add both positives and measured negatives")
    if positives < min_group_positives:
        return ("blocked", "needs_more_positives", "add measured positives in this family/source/scaffold")
    if negatives < min_group_negatives:
        return ("blocked", "needs_more_measured_negatives", "add measured inactive/control rows in this stratum")
    if not math.isnan(rate) and rate >= severe_positive_rate:
        return ("diagnostic_only", "positive_rich_stratum", "use source/family calibration or add measured negatives")
    if not math.isnan(rate) and rate <= sparse_positive_rate:
        return ("diagnostic_only", "negative_rich_stratum", "use stratified metrics or add positives")
    if unknown_fraction >= max_unknown_fraction:
        return ("diagnostic_only", "mostly_unknown", "use PU/background treatment and top-K recovery, not ordinary supervised metrics")
    return ("ok", "none", "no immediate data balance blocker")


def _iter_groups(df: pd.DataFrame, group_col: str) -> list[tuple[str, pd.DataFrame]]:
    if group_col not in df.columns:
        return []
    if group_col in SOURCE_LIKE_COLUMNS:
        exploded_rows: list[tuple[int, str]] = []
        for idx, value in df[group_col].items():
            tokens = split_source_tokens(value)
            if not tokens:
                tokens = ["<missing>"]
            for token in tokens:
                exploded_rows.append((idx, token))
        if not exploded_rows:
            return []
        exploded = pd.DataFrame(exploded_rows, columns=["_idx", "_group"])
        groups = []
        for value, sub in exploded.groupby("_group", dropna=False):
            groups.append((str(value), df.loc[sub["_idx"].unique()]))
        return groups
    out = []
    for value, group in df.groupby(df[group_col].fillna("<missing>").astype(str), dropna=False):
        out.append((str(value), group))
    return out


def _first_available(df: pd.DataFrame, candidates: Sequence[str]) -> str | None:
    for col in candidates:
        if col in df.columns and df[col].notna().any():
            return col
    return None


def _recommended_label_use(status: str) -> str:
    if status == "ok":
        return "supervised evaluation allowed if hard holdout readiness also passes"
    if status == "positive_unlabeled_evidence_layer":
        return "PU/top-K/evidence-ranking only until independent negatives or controls exist"
    if status == "single_class":
        return "not supervised-trainable as a two-class classifier"
    return "not labelable yet"


def _holdout_action(split: str, group_col: str) -> str:
    if split == "source_holdout":
        return "add labelable rows from independent sources with both positives and negatives"
    if split == "target_family_holdout":
        return "add positives and measured negatives across multiple target families"
    if split in {"chemical_cluster_holdout", "scaffold_holdout"}:
        return "add scaffold-diverse positives and measured negatives"
    if split == "target_holdout":
        return "add multiple targets with both classes so new-target tests have support"
    return f"add balanced label support for {group_col}"


def _write_markdown_summary(path: Path, manifest: dict[str, Any], label_summary: pd.DataFrame, priorities: pd.DataFrame) -> None:
    lines = [
        "# Atlas ML Data-Gap Report",
        "",
        f"Dataset: `{manifest['dataset']}`",
        f"Rows: {manifest['n_rows']}",
        "",
        "## Policy",
        "",
        "- Provenance columns explain source/family/scaffold failure modes; they are not clean predictive features.",
        "- Blocked holdout axes should block claim-grade generalization on that axis.",
        "- Positive-unlabeled labels should be evaluated with top-K recovery/enrichment unless independent controls exist.",
        "",
        "## Labels",
        "",
    ]
    if label_summary.empty:
        lines.append("No configured labels were present.")
    else:
        preview = label_summary[["label", "n_labelable", "n_positive", "n_negative", "n_unknown", "status"]].to_dict("records")
        lines.append("```json")
        lines.append(json.dumps(preview, indent=2, default=str))
        lines.append("```")
    lines.extend(["", "## Top Priorities", ""])
    if priorities.empty:
        lines.append("No blocking priorities were emitted under the configured thresholds.")
    else:
        for row in priorities.head(20).to_dict("records"):
            lines.append(
                f"- `{row['label']}` / `{row['axis']}` / `{row['value']}`: "
                f"{row['reason']} -> {row['recommended_action']}"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
