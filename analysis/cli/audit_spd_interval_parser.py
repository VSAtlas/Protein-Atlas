from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Iterable

import pandas as pd

from analysis.external.spd import (
    SPD_LABEL_POLICY_VERSION,
    _spd_label_status,
    aggregate_spd_assays,
    normalize_spd_table,
)
from analysis.ml.spd_four_expert_tables import (
    BINDING_LABEL_POLICY_VERSION,
    _binding_label,
    _exposure_label,
)
from analysis.ml.spd_panel import _extract_target_id, _first_nonempty
from analysis.spd_activity_policy import (
    SPD_ACTIVITY_PARSER_VERSION,
    interval_record,
    label_spd_binding_observation,
    label_spd_exposure_observation,
    parse_spd_activity_interval,
)


EXPECTED_GIT_SHA = "3554269ce088b709344f6a22dd9a9cd8958c6105"
EXPECTED_RELATIONS = {"=": 12_792, ">": 108_305}
EXPECTED_LEGACY_BINDING = {
    "positive": 2_484,
    "negative": 88_345,
    "unknown": 4_683,
}
EXPECTED_LEGACY_EXPOSURE = {
    "positive": 1_397,
    "negative": 37_638,
    "unknown": 56_477,
}

BEHAVIOR_COLUMNS = [
    "source_row_id",
    "source_excel_row",
    "drug_identity",
    "target_identity",
    "raw_relation",
    "raw_activity_value",
    "activity_unit",
    "assay_id",
    "assay_group",
    "assay_group_name",
    "assay_mode",
    "assay_type",
    "species",
    "source_representative",
    "production_selected",
    "materialized_binding_label",
    "legacy_binding_label",
    "canonical_binding_label",
    "materialized_exposure_label",
    "legacy_exposure_label",
    "canonical_exposure_label",
    "binding_changed",
    "exposure_changed",
    "legacy_binding_status",
    "canonical_binding_status",
    "legacy_exposure_status",
    "canonical_exposure_status",
    "canonical_binding_reason",
    "canonical_exposure_reason",
    "activity_value_uM",
    "activity_lower_bound_uM",
    "activity_upper_bound_uM",
    "lower_inclusive",
    "upper_inclusive",
    "normalized_relation",
    "activity_parse_status",
    "activity_parse_reason",
    "free_cmax_um",
    "canonical_margin_lower",
    "canonical_margin_upper",
]

RELATION_COLUMNS = [
    "input_stage",
    "relation_category",
    "raw_relation",
    "normalized_relation",
    "row_count",
    "nonmissing_activity_count",
    "examples",
]

BOUNDARY_COLUMNS = [
    "case_id",
    "activity_value",
    "activity_relation",
    "activity_unit",
    "free_cmax_um",
    "expected_parse_status",
    "expected_binding_label",
    "actual_binding_label",
    "expected_exposure_label",
    "actual_exposure_label",
    "passed",
]


def _git(root: Path, *args: str, binary: bool = False) -> str | bytes:
    value = subprocess.check_output(["git", *args], cwd=root)
    return value if binary else value.decode("utf-8", errors="replace")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _clean(value: Any) -> str:
    if value is None:
        return ""
    try:
        if bool(pd.isna(value)):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _key_text(value: Any) -> str:
    text = _clean(value)
    if text.endswith(".0"):
        try:
            return str(int(float(text)))
        except ValueError:
            return text
    return text


def _json_scalar(value: Any) -> Any:
    if value is None:
        return None
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _label_name(value: Any) -> str:
    if value is None:
        return "unknown"
    try:
        if bool(pd.isna(value)):
            return "unknown"
    except (TypeError, ValueError):
        pass
    text = str(value).strip().lower()
    if text in {"positive", "1", "1.0", "true"}:
        return "positive"
    if text in {"negative", "0", "0.0", "false"}:
        return "negative"
    if text in {"-1", "-1.0", "excluded", "conflict", "excluded_conflict"}:
        return "excluded_conflict"
    return "unknown"


def _binding_status(value: Any) -> str:
    return {
        "positive": "positive",
        "negative": "negative",
        "unknown": "unknown",
        "excluded_conflict": "excluded_conflict",
    }[_label_name(value)]


def _exposure_label_from_status(status: Any, margin: Any) -> int | None:
    text = _clean(status)
    if not text.startswith("labeled_"):
        return None
    numeric_margin = pd.to_numeric(pd.Series([margin]), errors="coerce").iloc[0]
    if pd.isna(numeric_margin):
        return None
    return 1 if float(numeric_margin) <= 10.0 else 0


def _bool_label(value: Any) -> int | None:
    text = _clean(value).lower()
    if text in {"true", "1", "1.0", "yes"}:
        return 1
    if text in {"false", "0", "0.0", "no"}:
        return 0
    return None


def _relation_category(value: Any) -> str:
    if value is None:
        return "null"
    try:
        if bool(pd.isna(value)):
            return "null"
    except (TypeError, ValueError):
        pass
    text = str(value)
    if text == "":
        return "blank"
    if text.strip() == "":
        return "whitespace"
    stripped = text.strip()
    if stripped in {"=", ">", ">=", "<", "<="}:
        return stripped
    return "unexpected"


def _serialize_bound(value: Any) -> Any:
    if value is None:
        return pd.NA
    try:
        if math.isinf(float(value)):
            return "inf"
    except (TypeError, ValueError):
        return value
    return value


def _relation_summary(
    frame: pd.DataFrame,
    *,
    stage: str,
    relation_col: str,
    activity_col: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    work = frame[[relation_col, activity_col]].copy()
    work["_category"] = work[relation_col].map(_relation_category)
    work["_activity_numeric"] = pd.to_numeric(work[activity_col], errors="coerce")
    for (category, raw_value), group in work.groupby(
        ["_category", relation_col], dropna=False, sort=False
    ):
        raw_display = "<NULL>" if category == "null" else repr(raw_value)
        normalized = (
            str(raw_value).strip() if category in {"=", ">", ">=", "<", "<="} else ""
        )
        examples = "; ".join(
            str(value)
            for value in frame.loc[group.index, activity_col].dropna().head(5).tolist()
        )
        rows.append(
            {
                "input_stage": stage,
                "relation_category": category,
                "raw_relation": raw_display,
                "normalized_relation": normalized,
                "row_count": int(len(group)),
                "nonmissing_activity_count": int(
                    group["_activity_numeric"].notna().sum()
                ),
                "examples": examples,
            }
        )
    present = Counter(work["_category"])
    for category in [
        "=",
        "blank",
        "null",
        ">",
        ">=",
        "<",
        "<=",
        "whitespace",
        "unexpected",
    ]:
        if present.get(category, 0):
            continue
        rows.append(
            {
                "input_stage": stage,
                "relation_category": category,
                "raw_relation": "",
                "normalized_relation": category
                if category in {"=", ">", ">=", "<", "<="}
                else "",
                "row_count": 0,
                "nonmissing_activity_count": 0,
                "examples": "",
            }
        )
    return rows


def _prepare_runtime_observations(workbook: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_excel(workbook, sheet_name="S Data 1", header=2)
    normalized = normalize_spd_table(workbook)
    if len(raw) != len(normalized):
        raise RuntimeError(
            f"normalization changed row count: workbook={len(raw)} normalized={len(normalized)}"
        )
    normalized = normalized.copy()
    normalized["_source_position"] = range(len(normalized))
    normalized["_source_excel_row"] = normalized["_source_position"] + 4
    normalized["_source_row_id"] = normalized.get(
        "RowId", pd.Series(normalized["_source_position"], index=normalized.index)
    ).map(_key_text)
    if (
        normalized["_source_row_id"].eq("").any()
        or normalized["_source_row_id"].duplicated().any()
    ):
        raise RuntimeError("workbook source RowId is missing or duplicated")
    normalized["drug_id"] = _first_nonempty(normalized, ["name", "drug_id"]).str.lower()
    normalized["target_id"] = _extract_target_id(normalized).astype(str).str.upper()
    normalized["activity_relation"] = normalized["activity_relation"].astype("object")
    normalized["_ac50_uM"] = (
        pd.to_numeric(normalized["ac50_nM"], errors="coerce") / 1000.0
    )
    normalized["_free_cmax_uM"] = (
        pd.to_numeric(normalized["free_cmax_nM"], errors="coerce") / 1000.0
    )
    normalized["_legacy_margin"] = (
        pd.to_numeric(normalized["ac50_nM"], errors="coerce")
        / pd.to_numeric(normalized["free_cmax_nM"], errors="coerce")
    ).where(pd.to_numeric(normalized["free_cmax_nM"], errors="coerce").gt(0))
    return raw, normalized


def _legacy_raw_results(frame: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series]:
    binding_source = pd.DataFrame(
        {
            "spd_ac50_uM": frame["_ac50_uM"],
            "spd_activity_relation": frame["activity_relation"],
        },
        index=frame.index,
    )
    binding = _binding_label(binding_source, active_um=1.0, inactive_um=10.0)
    status_frame = frame.copy()
    status_frame["exposure_margin"] = status_frame["_legacy_margin"]
    exposure_status = status_frame.apply(
        lambda row: _spd_label_status(row, 10.0, 100.0), axis=1
    )
    exposure = pd.Series(
        [
            _exposure_label_from_status(status, margin)
            for status, margin in zip(exposure_status, frame["_legacy_margin"])
        ],
        index=frame.index,
        dtype="Int64",
    )
    return binding, exposure, exposure_status


def _canonical_results(
    frame: pd.DataFrame,
) -> tuple[list[Any], list[Any], list[Any]]:
    intervals = []
    binding = []
    exposure = []
    for value, relation, free_cmax in zip(
        frame["_ac50_uM"],
        frame["activity_relation"],
        frame["_free_cmax_uM"],
    ):
        interval = parse_spd_activity_interval(value, relation, "uM")
        intervals.append(interval)
        binding.append(label_spd_binding_observation(interval))
        exposure.append(label_spd_exposure_observation(interval, free_cmax))
    return intervals, binding, exposure


def _selected_source_rows(frame: pd.DataFrame) -> pd.DataFrame:
    selected = frame.copy()
    selected["_is_biochemical"] = (
        selected.get("assay_type", "")
        .astype(str)
        .str.lower()
        .str.contains("biochem|target", regex=True, na=False)
    )
    selected = selected.sort_values(
        ["drug_id", "target_id", "_is_biochemical", "ac50_nM"],
        ascending=[True, True, False, True],
    )
    return selected.groupby(["drug_id", "target_id"], dropna=False, sort=False).head(1)


def _load_materialized_selection(path: Path) -> pd.DataFrame:
    usecols = [
        "code_state",
        "selection_stage",
        "drug_id",
        "target_id",
        "current_row_id",
        "current_binding_label",
        "current_exposure_label",
        "current_exposure_status",
    ]
    table = pd.read_csv(path, usecols=usecols, low_memory=False)
    table = table[
        table["code_state"].eq("committed_head")
        & table["selection_stage"].eq("raw_spd_pair_aggregation")
    ].copy()
    if table.duplicated(["drug_id", "target_id"]).any():
        raise RuntimeError(
            "materialized selection audit has duplicate committed pair keys"
        )
    table["drug_id"] = table["drug_id"].astype(str).str.lower()
    table["target_id"] = table["target_id"].astype(str).str.upper()
    table["current_row_id"] = table["current_row_id"].map(_key_text)
    return table


def _compare_pair_outputs(
    normalized: pd.DataFrame,
    materialized_pair_path: Path,
    materialized_selection_path: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    legacy_pairs = aggregate_spd_assays(
        normalized, margin_strong=10.0, margin_weak=100.0
    )
    materialized = pd.read_csv(materialized_pair_path, low_memory=False)
    selections = _load_materialized_selection(materialized_selection_path)
    for frame in (legacy_pairs, materialized):
        frame["drug_id"] = frame["drug_id"].astype(str).str.lower()
        frame["target_id"] = frame["target_id"].astype(str).str.upper()
    key = ["drug_id", "target_id"]
    if any(frame.duplicated(key).any() for frame in (legacy_pairs, materialized)):
        raise RuntimeError("pair comparison input contains duplicate drug-target keys")

    selected = _selected_source_rows(normalized)[
        [*key, "_source_row_id", "_ac50_uM", "_free_cmax_uM", "activity_relation"]
    ].copy()
    selected = selected.rename(columns={"activity_relation": "_selected_relation"})
    pair = (
        legacy_pairs.merge(
            materialized,
            on=key,
            how="outer",
            suffixes=("_legacy", "_materialized"),
            indicator=True,
            validate="one_to_one",
        )
        .merge(selections, on=key, how="left", validate="one_to_one")
        .merge(selected, on=key, how="left", validate="one_to_one")
    )
    pair["materialized_binding_label"] = pair["current_binding_label"].map(_label_name)
    binding_source = pd.DataFrame(
        {
            "spd_ac50_uM": pd.to_numeric(pair["ac50_nM_legacy"], errors="coerce")
            / 1000.0,
            "spd_activity_relation": pair["spd_activity_relation_legacy"],
        },
        index=pair.index,
    )
    pair["legacy_binding_label"] = _binding_label(
        binding_source, active_um=1.0, inactive_um=10.0
    ).map(_label_name)
    pair["legacy_exposure_label"] = [
        _label_name(_exposure_label_from_status(status, margin))
        for status, margin in zip(
            pair["spd_label_status_legacy"], pair["exposure_margin_legacy"]
        )
    ]
    pair["panel_materialized_exposure_label"] = (
        pair["spd_exposure_relevant_materialized"].map(_bool_label).map(_label_name)
    )
    pair["audit_materialized_exposure_label"] = pair["current_exposure_label"].map(
        _label_name
    )

    canonical_binding: list[str] = []
    canonical_exposure: list[str] = []
    canonical_status: list[str] = []
    for ac50, relation, free_cmax in zip(
        pair["ac50_nM_legacy"],
        pair["spd_activity_relation_legacy"],
        pair["free_cmax_nM_legacy"],
    ):
        interval = parse_spd_activity_interval(ac50, relation, "nM")
        binding = label_spd_binding_observation(interval)
        exposure = label_spd_exposure_observation(
            interval,
            None if pd.isna(free_cmax) else float(free_cmax) / 1000.0,
        )
        canonical_binding.append(_label_name(binding.numeric_label))
        canonical_exposure.append(_label_name(exposure.numeric_label))
        canonical_status.append(exposure.exposure_status)
    pair["canonical_binding_label"] = canonical_binding
    pair["canonical_exposure_label"] = canonical_exposure
    pair["canonical_exposure_status"] = canonical_status

    float_matches = {}
    for column in ("ac50_nM", "free_cmax_nM", "exposure_margin"):
        left = pd.to_numeric(pair[f"{column}_legacy"], errors="coerce")
        right = pd.to_numeric(pair[f"{column}_materialized"], errors="coerce")
        equal = (left.isna() & right.isna()) | (left - right).abs().le(
            1e-9
            * pd.concat([left.abs(), right.abs()], axis=1).max(axis=1).clip(lower=1.0)
        )
        float_matches[column] = int((~equal).sum())
    relation_match = pair["spd_activity_relation_legacy"].fillna("<NA>").astype(
        str
    ) == pair["spd_activity_relation_materialized"].fillna("<NA>").astype(str)
    status_match = pair["spd_label_status_legacy"].fillna("<NA>").astype(str) == pair[
        "spd_label_status_materialized"
    ].fillna("<NA>").astype(str)
    selected_row_match = (
        pair["_source_row_id"].map(_key_text).eq(pair["current_row_id"].map(_key_text))
    )

    summary = {
        "legacy_pair_rows": int(len(legacy_pairs)),
        "materialized_pair_rows": int(len(materialized)),
        "materialized_selection_rows": int(len(selections)),
        "pair_key_merge_counts": {
            str(k): int(v) for k, v in pair["_merge"].value_counts(dropna=False).items()
        },
        "materialized_vs_legacy_float_mismatches": float_matches,
        "materialized_vs_legacy_relation_mismatches": int((~relation_match).sum()),
        "materialized_vs_legacy_status_mismatches": int((~status_match).sum()),
        "materialized_audit_vs_legacy_status_mismatches": int(
            pair["current_exposure_status"]
            .fillna("<NA>")
            .astype(str)
            .ne(pair["spd_label_status_legacy"].fillna("<NA>").astype(str))
            .sum()
        ),
        "materialized_vs_replayed_selected_row_id_mismatches": int(
            (~selected_row_match).sum()
        ),
        "materialized_audit_vs_legacy_binding_mismatches": int(
            pair["materialized_binding_label"].ne(pair["legacy_binding_label"]).sum()
        ),
        "materialized_panel_vs_legacy_exposure_mismatches": int(
            pair["panel_materialized_exposure_label"]
            .ne(pair["legacy_exposure_label"])
            .sum()
        ),
        "materialized_audit_vs_panel_exposure_mismatches": int(
            pair["audit_materialized_exposure_label"]
            .ne(pair["panel_materialized_exposure_label"])
            .sum()
        ),
        "legacy_vs_canonical_binding_mismatches": int(
            pair["legacy_binding_label"].ne(pair["canonical_binding_label"]).sum()
        ),
        "legacy_vs_canonical_exposure_mismatches": int(
            pair["legacy_exposure_label"].ne(pair["canonical_exposure_label"]).sum()
        ),
    }
    for source in (
        "materialized_binding_label",
        "legacy_binding_label",
        "canonical_binding_label",
        "panel_materialized_exposure_label",
        "legacy_exposure_label",
        "canonical_exposure_label",
    ):
        summary[f"{source}_counts"] = {
            label: int(pair[source].eq(label).sum())
            for label in ("positive", "negative", "unknown", "excluded_conflict")
        }
    return pair, summary


def _compare_materialized_expert(path: Path) -> dict[str, Any]:
    wanted = [
        "drug_id",
        "target_id",
        "pdb_id",
        "spd_ac50_uM",
        "spd_activity_relation",
        "free_cmax_um",
        "spd_binding_label",
        "spd_exposure_label",
    ]
    header = pd.read_csv(path, nrows=0)
    missing = [column for column in wanted if column not in header.columns]
    if missing:
        raise RuntimeError(
            f"{path} is missing expert reconciliation columns: {missing}"
        )
    frame = pd.read_csv(path, usecols=wanted, low_memory=False)
    legacy_binding = _binding_label(frame, active_um=1.0, inactive_um=10.0)
    legacy_exposure = _exposure_label(frame)
    materialized_binding = pd.to_numeric(
        frame["spd_binding_label"], errors="coerce"
    ).astype("Int64")
    materialized_exposure = (
        pd.to_numeric(frame["spd_exposure_label"], errors="coerce")
        .replace(-1, pd.NA)
        .astype("Int64")
    )

    canonical_binding: list[int | None] = []
    canonical_exposure: list[int | None] = []
    changed_rows: list[dict[str, Any]] = []
    for row_number, row in frame.iterrows():
        interval = parse_spd_activity_interval(
            row["spd_ac50_uM"], row["spd_activity_relation"], "uM"
        )
        binding = label_spd_binding_observation(interval)
        exposure = label_spd_exposure_observation(interval, row["free_cmax_um"])
        canonical_binding.append(binding.numeric_label)
        canonical_exposure.append(exposure.numeric_label)
        old_binding = _label_name(materialized_binding.loc[row_number])
        old_exposure = _label_name(materialized_exposure.loc[row_number])
        new_binding = _label_name(binding.numeric_label)
        new_exposure = _label_name(exposure.numeric_label)
        if old_binding != new_binding or old_exposure != new_exposure:
            changed_rows.append(
                {
                    "source_row_number": int(row_number),
                    "drug_id": _json_scalar(row["drug_id"]),
                    "target_id": _json_scalar(row["target_id"]),
                    "pdb_id": _json_scalar(row["pdb_id"]),
                    "relation": _json_scalar(row["spd_activity_relation"]),
                    "ac50_uM": _json_scalar(row["spd_ac50_uM"]),
                    "free_cmax_um": _json_scalar(row["free_cmax_um"]),
                    "materialized_binding": old_binding,
                    "canonical_binding": new_binding,
                    "materialized_exposure": old_exposure,
                    "canonical_exposure": new_exposure,
                }
            )
    canonical_binding_series = pd.Series(canonical_binding, index=frame.index)
    canonical_exposure_series = pd.Series(canonical_exposure, index=frame.index)
    return {
        "path": str(path.resolve()),
        "sha256": _sha256(path),
        "rows": int(len(frame)),
        "row_key": ["drug_id", "pdb_id"],
        "unique_row_keys": int(frame[["drug_id", "pdb_id"]].drop_duplicates().shape[0]),
        "materialized_vs_legacy_binding_mismatches": int(
            materialized_binding.map(_label_name)
            .ne(legacy_binding.map(_label_name))
            .sum()
        ),
        "materialized_vs_legacy_exposure_mismatches": int(
            materialized_exposure.map(_label_name)
            .ne(legacy_exposure.map(_label_name))
            .sum()
        ),
        "materialized_vs_canonical_binding_mismatches": int(
            materialized_binding.map(_label_name)
            .ne(canonical_binding_series.map(_label_name))
            .sum()
        ),
        "materialized_vs_canonical_exposure_mismatches": int(
            materialized_exposure.map(_label_name)
            .ne(canonical_exposure_series.map(_label_name))
            .sum()
        ),
        "materialized_binding_counts": _counts(materialized_binding),
        "legacy_binding_counts": _counts(legacy_binding),
        "canonical_binding_counts": _counts(canonical_binding_series),
        "materialized_exposure_counts": _counts(materialized_exposure),
        "legacy_exposure_counts": _counts(legacy_exposure),
        "canonical_exposure_counts": _counts(canonical_exposure_series),
        "changed_rows": changed_rows,
    }


def _boundary_cases() -> pd.DataFrame:
    cases = [
        ("binding_exact_active", 1.0, "=", "uM", 1.0, "valid", 1, 1),
        ("binding_lower_unknown", 1.0, ">=", "uM", 1.0, "valid", None, None),
        ("binding_exact_inactive", 10.0, "=", "uM", 1.0, "valid", 0, 1),
        ("binding_open_inactive", 10.0, ">", "uM", 1.0, "valid", 0, 0),
        ("binding_closed_inactive", 10.0, ">=", "uM", 1.0, "valid", 0, None),
        ("exposure_upper_positive", 10.0, "<=", "uM", 1.0, "valid", None, 1),
        ("exposure_crossing_upper", 20.0, "<", "uM", 1.0, "valid", None, None),
        ("exposure_crossing_lower", 5.0, ">", "uM", 1.0, "valid", None, None),
        ("blank_relation", 1.0, "", "uM", 1.0, "invalid", None, None),
        ("missing_relation", 1.0, None, "uM", 1.0, "invalid", None, None),
        ("malformed_value", "bad", "=", "uM", 1.0, "invalid", None, None),
        ("unsupported_unit", 1.0, "=", "mg/L", 1.0, "invalid", None, None),
        ("negative_value", -1.0, "=", "uM", 1.0, "invalid", None, None),
        ("missing_cmax", 1.0, "=", "uM", None, "valid", 1, None),
        ("zero_cmax", 1.0, "=", "uM", 0.0, "valid", 1, None),
    ]
    rows = []
    for (
        case_id,
        value,
        relation,
        unit,
        free_cmax,
        expected_parse,
        expected_binding,
        expected_exposure,
    ) in cases:
        interval = parse_spd_activity_interval(value, relation, unit)
        binding = label_spd_binding_observation(interval)
        exposure = label_spd_exposure_observation(interval, free_cmax)
        passed = (
            interval.activity_parse_status == expected_parse
            and binding.numeric_label == expected_binding
            and exposure.numeric_label == expected_exposure
        )
        rows.append(
            {
                "case_id": case_id,
                "activity_value": value,
                "activity_relation": relation,
                "activity_unit": unit,
                "free_cmax_um": free_cmax,
                "expected_parse_status": expected_parse,
                "expected_binding_label": expected_binding,
                "actual_binding_label": binding.numeric_label,
                "expected_exposure_label": expected_exposure,
                "actual_exposure_label": exposure.numeric_label,
                "passed": passed,
            }
        )
    return pd.DataFrame(rows, columns=BOUNDARY_COLUMNS)


def _csv_profile(
    path: Path,
    *,
    stage: str,
    key_cols: Iterable[str],
    label_col: str | None = None,
    eligible_col: str | None = None,
) -> dict[str, Any]:
    header = pd.read_csv(path, nrows=0)
    wanted = [column for column in [*key_cols, label_col, eligible_col] if column]
    missing = [column for column in wanted if column not in header.columns]
    if missing:
        raise RuntimeError(f"{path} is missing lineage columns: {missing}")
    frame = pd.read_csv(path, usecols=wanted, low_memory=False)
    keys = list(key_cols)
    missing_key = pd.Series(False, index=frame.index)
    for column in keys:
        missing_key |= frame[column].isna() | frame[column].astype(str).str.strip().eq(
            ""
        )
    profile: dict[str, Any] = {
        "stage": stage,
        "path": str(path.resolve()),
        "sha256": _sha256(path),
        "rows": int(len(frame)),
        "row_key": keys,
        "missing_row_keys": int(missing_key.sum()),
        "unique_row_keys": int(frame[keys].drop_duplicates().shape[0]),
        "duplicate_row_key_rows": int(frame.duplicated(keys, keep=False).sum()),
    }
    if label_col:
        labels = pd.to_numeric(frame[label_col], errors="coerce")
        labelable = labels.isin([0, 1])
        profile["label_col"] = label_col
        profile["labelable_rows"] = int(labelable.sum())
        if eligible_col:
            eligible = (
                frame[eligible_col]
                .fillna(False)
                .astype(str)
                .str.strip()
                .str.lower()
                .isin({"1", "1.0", "true", "yes"})
            )
            profile["eligible_col"] = eligible_col
            profile["eligible_labelable_rows"] = int((eligible & labelable).sum())
    return profile


def _lineage_profiles(args: argparse.Namespace) -> list[dict[str, Any]]:
    specs = [
        (
            args.phase1_source,
            "frozen_phase1_source",
            ("drug_id", "pdb_id"),
            "spd_binding_label",
            None,
        ),
        (
            args.historical_master,
            "historical_run_master",
            ("pdb_id", "ligand_base"),
            None,
            None,
        ),
        (
            args.historical_binding_raw,
            "historical_binding_raw_expert",
            ("pdb_id", "ligand_base"),
            "spd_binding_label",
            None,
        ),
        (
            args.historical_binding_model_ready,
            "historical_binding_model_ready",
            ("drug_id", "target_id"),
            "spd_binding_label",
            None,
        ),
        (
            args.historical_exposure_raw,
            "historical_exposure_raw_expert",
            ("pdb_id", "ligand_base"),
            "spd_exposure_label",
            None,
        ),
        (
            args.historical_exposure_model_ready,
            "historical_exposure_model_ready",
            ("drug_id", "target_id"),
            "spd_exposure_label",
            None,
        ),
        (
            args.strict_binding_raw,
            "mapping_v9_strict_binding_raw",
            ("drug_id", "pdb_id"),
            "spd_binding_label",
            "spd_receptor_mapping_strict_eligible",
        ),
        (
            args.strict_binding_model_ready,
            "mapping_v9_strict_binding_model_ready",
            ("drug_id", "target_id"),
            "spd_binding_label",
            None,
        ),
        (
            args.strict_exposure_raw,
            "mapping_v9_strict_exposure_raw",
            ("drug_id", "pdb_id"),
            "spd_exposure_label",
            "spd_receptor_mapping_strict_eligible",
        ),
        (
            args.strict_exposure_model_ready,
            "mapping_v9_strict_exposure_model_ready",
            ("drug_id", "target_id"),
            "spd_exposure_label",
            None,
        ),
    ]
    profiles = []
    for path, stage, keys, label, eligible in specs:
        if path is None:
            continue
        profiles.append(
            _csv_profile(
                path,
                stage=stage,
                key_cols=keys,
                label_col=label,
                eligible_col=eligible,
            )
        )
    return profiles


def _counts(values: Iterable[Any]) -> dict[str, int]:
    labels = [_label_name(value) for value in values]
    return {
        label: int(labels.count(label))
        for label in ("positive", "negative", "unknown", "excluded_conflict")
    }


def _transition_records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    if frame.empty:
        return []
    columns = [
        "raw_relation",
        "canonical_margin_lower",
        "legacy_exposure_label",
        "canonical_exposure_label",
        "legacy_exposure_status",
        "canonical_exposure_reason",
        "drug_identity",
        "production_selected",
    ]
    records = []
    grouped = frame.groupby(columns, dropna=False, sort=True).size()
    for key, count in grouped.items():
        values = key if isinstance(key, tuple) else (key,)
        record = {column: _json_scalar(value) for column, value in zip(columns, values)}
        record["row_count"] = int(count)
        records.append(record)
    return records


def _write_summary(
    path: Path,
    *,
    git_sha: str,
    raw_rows: int,
    pair_summary: dict[str, Any],
    raw_binding_changes: int,
    raw_exposure_changes: int,
    relation_counts: dict[str, int],
    lineage: list[dict[str, Any]],
    expert_reconciliation: dict[str, Any] | None,
    changed_behavior: pd.DataFrame,
    boundary_passed: bool,
    gate_failures: list[str],
) -> None:
    lineage_lines = [
        "| Stage | Artifact | Rows | Row key | Labelable/eligible | Unique keys |",
        "|---|---|---:|---|---:|---:|",
    ]
    for item in lineage:
        labelable = item.get("eligible_labelable_rows", item.get("labelable_rows", "—"))
        lineage_lines.append(
            f"| {item['stage']} | `{item['path']}` | {item['rows']:,} | "
            f"`{' × '.join(item['row_key'])}` | {labelable} | {item['unique_row_keys']:,} |"
        )
    binding_counts = pair_summary["canonical_binding_label_counts"]
    exposure_counts = pair_summary["canonical_exposure_label_counts"]
    selected_changes = changed_behavior[
        changed_behavior["production_selected"].astype(bool)
    ]
    changed_drugs = (
        changed_behavior["drug_identity"].value_counts().sort_index().to_dict()
    )
    selected_changed_drugs = (
        selected_changes["drug_identity"].value_counts().sort_index().to_dict()
    )
    expert_line = ""
    if expert_reconciliation is not None:
        expert_line = (
            "\nMaterialized frozen-Phase1 reconciliation: "
            f"materialized-versus-committed binding mismatches="
            f"{expert_reconciliation['materialized_vs_legacy_binding_mismatches']:,}, "
            f"materialized-versus-canonical binding mismatches="
            f"{expert_reconciliation['materialized_vs_canonical_binding_mismatches']:,}, "
            f"materialized-versus-canonical exposure mismatches="
            f"{expert_reconciliation['materialized_vs_canonical_exposure_mismatches']:,}. "
            "These frozen-source differences are reported separately from the "
            "95,512-pair comparison; they indicate artifact/code drift, not a "
            "new parser-versus-committed binding difference."
        )
    gate = "PASS" if not gate_failures else "FAIL"
    text = f"""# SPD interval parser behavior comparison

- Git SHA: `{git_sha}`
- Parser: `{SPD_ACTIVITY_PARSER_VERSION}`
- Phase A gate: **{gate}**
- Raw observations compared: {raw_rows:,}
- Boundary cases: {"PASS" if boundary_passed else "FAIL"}

## Relation inventory

`=`: {relation_counts.get("=", 0):,}; `>`: {relation_counts.get(">", 0):,};
`<`: {relation_counts.get("<", 0):,}; `<=`: {relation_counts.get("<=", 0):,};
`>=`: {relation_counts.get(">=", 0):,}; blank/null/whitespace/unexpected:
{sum(relation_counts.get(key, 0) for key in ("blank", "null", "whitespace", "unexpected")):,}.

## Three-way reconciliation

Materialized production artifacts, committed legacy functions, and the canonical
parser were compared independently. Raw legacy-versus-canonical changes:
binding={raw_binding_changes:,}, exposure={raw_exposure_changes:,}.
Pair legacy-versus-canonical mismatches:
binding={pair_summary["legacy_vs_canonical_binding_mismatches"]:,},
exposure={pair_summary["legacy_vs_canonical_exposure_mismatches"]:,}.
{expert_line}

Canonical pair counts:

- Binding: positive={binding_counts["positive"]:,}, negative={binding_counts["negative"]:,}, unknown={binding_counts["unknown"]:,}
- Exposure: positive={exposure_counts["positive"]:,}, negative={exposure_counts["negative"]:,}, unknown={exposure_counts["unknown"]:,}

Gate failures: {", ".join(gate_failures) if gate_failures else "none"}.

Every current-data exposure mismatch is retained in
`parser_behavior_comparison.csv` with `exposure_changed=true`. There are
{len(changed_behavior):,} affected raw rows and {len(selected_changes):,}
production-selected pairs. Raw affected-drug counts: `{changed_drugs}`;
selected affected-drug counts: `{selected_changed_drugs}`. The affected
observations all have relation `>` and a lower exposure margin exactly 10:
legacy status is unknown, while the canonical open interval `(10,+∞)` is a
definite negative.

## Drug-PDB versus drug-target lineage

{chr(10).join(lineage_lines)}

The historical run master and raw expert tables contain 440,183 physical
`pdb_id × ligand_base` rows. The historical binding model-ready table contains
23,849 rows only **after** filtering and `drug_id × target_id` deduplication.
The frozen 23,849-row Phase 1 table was then reused as the input to the receptor
mapping comparison; it is locally pre-dedup for that second pass, but is not an
end-to-end pre-dedup artifact.

In the reviewed strict mapping run, the raw expert stage is
`strict/ml_spd_{{binding,exposure}}_table.csv`, keyed by `drug_id × pdb_id`.
The later `_write_raw_and_model_ready` stage invokes `deduplicate_ml_rows` on
`drug_id × target_id`. It removes 365 concordant binding and 156 concordant
exposure 9I52/9LLG DRD1 pairs. No deduplication behavior is changed by Prompt 2,
and the retained structure is not asserted to be biologically preferred.
"""
    path.write_text(text, encoding="utf-8")


def run_comparison(args: argparse.Namespace) -> dict[str, Any]:
    started = datetime.now(timezone.utc)
    root = args.repo_root.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    git_sha = str(_git(root, "rev-parse", "HEAD")).strip()
    if git_sha != args.expected_git_sha:
        raise RuntimeError(f"expected git SHA {args.expected_git_sha}, found {git_sha}")

    production_sources = [
        "analysis/external/spd.py",
        "analysis/ml/spd_four_expert_tables.py",
        "analysis/ml/spd_label_enrichment.py",
    ]
    source_hashes: dict[str, dict[str, Any]] = {}
    for relative in production_sources:
        working = (root / relative).read_bytes()
        committed = _git(root, "show", f"HEAD:{relative}", binary=True)
        if not isinstance(committed, bytes):
            raise AssertionError("binary git show unexpectedly returned text")
        source_hashes[relative] = {
            "working_sha256": _sha256_bytes(working),
            "committed_sha256": _sha256_bytes(committed),
            "matches_committed_head": working == committed,
        }
        if working != committed:
            raise RuntimeError(
                f"Phase A requires unmodified production derivation source: {relative}"
            )

    raw, normalized = _prepare_runtime_observations(args.spd_workbook)
    legacy_binding, legacy_exposure, legacy_exposure_status = _legacy_raw_results(
        normalized
    )
    intervals, canonical_binding, canonical_exposure = _canonical_results(normalized)
    pair, pair_summary = _compare_pair_outputs(
        normalized, args.pair_table, args.materialized_selection
    )
    expert_reconciliation = (
        _compare_materialized_expert(args.phase1_source)
        if args.phase1_source is not None
        else None
    )
    materialized_by_row = pair.set_index("current_row_id")[
        ["materialized_binding_label", "audit_materialized_exposure_label"]
    ].to_dict("index")
    selected_ids = set(materialized_by_row)

    rows = []
    for position, (
        index,
        source_row,
        old_binding,
        old_exposure,
        old_exposure_status,
        interval,
        new_binding,
        new_exposure,
    ) in enumerate(
        zip(
            normalized.index,
            normalized.to_dict("records"),
            legacy_binding,
            legacy_exposure,
            legacy_exposure_status,
            intervals,
            canonical_binding,
            canonical_exposure,
        )
    ):
        del index
        row_id = _key_text(source_row["_source_row_id"])
        materialized = materialized_by_row.get(row_id, {})
        interval_fields = interval_record(interval)
        rows.append(
            {
                "source_row_id": row_id,
                "source_excel_row": source_row["_source_excel_row"],
                "drug_identity": source_row["drug_id"],
                "target_identity": source_row["target_id"],
                "raw_relation": source_row["activity_relation"],
                "raw_activity_value": source_row["_ac50_uM"],
                "activity_unit": "uM",
                "assay_id": source_row.get("assay_id", pd.NA),
                "assay_group": source_row.get("assay_group", pd.NA),
                "assay_group_name": source_row.get("assay_group_name", pd.NA),
                "assay_mode": source_row.get("Mode", pd.NA),
                "assay_type": source_row.get("assay_type", pd.NA),
                "species": source_row.get("Protein/Target/Species", pd.NA),
                "source_representative": source_row.get(
                    "representative_result_drug_assay_group_pair", pd.NA
                ),
                "production_selected": row_id in selected_ids,
                "materialized_binding_label": materialized.get(
                    "materialized_binding_label", ""
                ),
                "legacy_binding_label": _label_name(old_binding),
                "canonical_binding_label": _label_name(new_binding.numeric_label),
                "materialized_exposure_label": materialized.get(
                    "audit_materialized_exposure_label", ""
                ),
                "legacy_exposure_label": _label_name(old_exposure),
                "canonical_exposure_label": _label_name(new_exposure.numeric_label),
                "binding_changed": _label_name(old_binding)
                != _label_name(new_binding.numeric_label),
                "exposure_changed": _label_name(old_exposure)
                != _label_name(new_exposure.numeric_label),
                "legacy_binding_status": _binding_status(old_binding),
                "canonical_binding_status": new_binding.binding_status,
                "legacy_exposure_status": old_exposure_status,
                "canonical_exposure_status": new_exposure.exposure_status,
                "canonical_binding_reason": new_binding.binding_reason,
                "canonical_exposure_reason": new_exposure.exposure_reason,
                "activity_value_uM": interval_fields["activity_value_uM"],
                "activity_lower_bound_uM": _serialize_bound(
                    interval_fields["activity_lower_bound_uM"]
                ),
                "activity_upper_bound_uM": _serialize_bound(
                    interval_fields["activity_upper_bound_uM"]
                ),
                "lower_inclusive": interval_fields["lower_inclusive"],
                "upper_inclusive": interval_fields["upper_inclusive"],
                "normalized_relation": interval_fields["normalized_relation"],
                "activity_parse_status": interval_fields["activity_parse_status"],
                "activity_parse_reason": interval_fields["activity_parse_reason"],
                "free_cmax_um": source_row["_free_cmax_uM"],
                "canonical_margin_lower": _serialize_bound(
                    new_exposure.exposure_margin_lower
                ),
                "canonical_margin_upper": _serialize_bound(
                    new_exposure.exposure_margin_upper
                ),
            }
        )
    behavior = pd.DataFrame(rows, columns=BEHAVIOR_COLUMNS)

    relation_rows = _relation_summary(
        raw,
        stage="raw_workbook_s_data_1",
        relation_col="summarized prefix",
        activity_col="summarized IC50",
    )
    relation_rows.extend(
        _relation_summary(
            normalized,
            stage="normalized_runtime",
            relation_col="activity_relation",
            activity_col="ac50_nM",
        )
    )
    relation_summary = pd.DataFrame(relation_rows, columns=RELATION_COLUMNS)
    boundary = _boundary_cases()
    lineage = _lineage_profiles(args)

    raw_binding_changes = int(behavior["binding_changed"].sum())
    raw_exposure_changes = int(behavior["exposure_changed"].sum())
    raw_categories = raw["summarized prefix"].map(_relation_category)
    relation_counts = {
        category: int(raw_categories.eq(category).sum())
        for category in [
            "=",
            ">",
            ">=",
            "<",
            "<=",
            "blank",
            "null",
            "whitespace",
            "unexpected",
        ]
    }
    gate_failures: list[str] = []
    if len(raw) != 121_097:
        gate_failures.append(f"raw_row_count={len(raw)}")
    if {key: relation_counts[key] for key in ("=", ">")} != EXPECTED_RELATIONS:
        gate_failures.append("raw_relation_counts")
    if any(
        relation_counts[key]
        for key in (">=", "<", "<=", "blank", "null", "whitespace", "unexpected")
    ):
        gate_failures.append("unexpected_current_relations")
    if raw_binding_changes:
        gate_failures.append(f"raw_binding_changes={raw_binding_changes}")
    if raw_exposure_changes:
        gate_failures.append(f"raw_exposure_changes={raw_exposure_changes}")
    for field in (
        "materialized_vs_legacy_relation_mismatches",
        "materialized_vs_legacy_status_mismatches",
        "materialized_vs_replayed_selected_row_id_mismatches",
        "materialized_audit_vs_legacy_binding_mismatches",
        "materialized_panel_vs_legacy_exposure_mismatches",
        "materialized_audit_vs_panel_exposure_mismatches",
        "materialized_audit_vs_legacy_status_mismatches",
        "legacy_vs_canonical_binding_mismatches",
        "legacy_vs_canonical_exposure_mismatches",
    ):
        if pair_summary[field]:
            gate_failures.append(f"{field}={pair_summary[field]}")
    if any(pair_summary["materialized_vs_legacy_float_mismatches"].values()):
        gate_failures.append("materialized_vs_legacy_float_mismatches")
    if pair_summary["canonical_binding_label_counts"] != {
        **EXPECTED_LEGACY_BINDING,
        "excluded_conflict": 0,
    }:
        gate_failures.append("binding_pair_counts")
    if pair_summary["canonical_exposure_label_counts"] != {
        **EXPECTED_LEGACY_EXPOSURE,
        "excluded_conflict": 0,
    }:
        gate_failures.append("canonical_exposure_counts_differ_from_locked_legacy")
    if expert_reconciliation is not None:
        for field in (
            "materialized_vs_legacy_binding_mismatches",
            "materialized_vs_legacy_exposure_mismatches",
            "materialized_vs_canonical_binding_mismatches",
            "materialized_vs_canonical_exposure_mismatches",
        ):
            if expert_reconciliation[field]:
                gate_failures.append(f"expert_{field}={expert_reconciliation[field]}")
    if not bool(boundary["passed"].all()):
        gate_failures.append("boundary_cases")

    behavior_path = out / "parser_behavior_comparison.csv"
    relation_path = out / "parser_relation_summary.csv"
    boundary_path = out / "parser_boundary_cases.csv"
    summary_path = out / "parser_summary.md"
    manifest_path = out / "parser_manifest.json"
    behavior.to_csv(behavior_path, index=False)
    relation_summary.to_csv(relation_path, index=False)
    boundary.to_csv(boundary_path, index=False)
    _write_summary(
        summary_path,
        git_sha=git_sha,
        raw_rows=len(raw),
        pair_summary=pair_summary,
        raw_binding_changes=raw_binding_changes,
        raw_exposure_changes=raw_exposure_changes,
        relation_counts=relation_counts,
        lineage=lineage,
        expert_reconciliation=expert_reconciliation,
        changed_behavior=behavior[behavior["exposure_changed"]].copy(),
        boundary_passed=bool(boundary["passed"].all()),
        gate_failures=gate_failures,
    )

    input_paths = [
        args.spd_workbook,
        args.pair_table,
        args.materialized_selection,
        *[
            path
            for path in (
                args.phase1_source,
                args.historical_master,
                args.historical_binding_raw,
                args.historical_binding_model_ready,
                args.historical_exposure_raw,
                args.historical_exposure_model_ready,
                args.strict_binding_raw,
                args.strict_binding_model_ready,
                args.strict_exposure_raw,
                args.strict_exposure_model_ready,
            )
            if path is not None
        ],
    ]
    implementation_paths = [
        root / "analysis/spd_activity_policy.py",
        root / "analysis/cli/audit_spd_interval_parser.py",
        root / "analysis/test_spd_activity_policy.py",
        root / "analysis/cli/test_audit_spd_interval_parser.py",
    ]
    implementation_digest = hashlib.sha256()
    implementation_hashes = []
    for path in implementation_paths:
        relative = path.relative_to(root).as_posix()
        content = path.read_bytes()
        implementation_digest.update(relative.encode("utf-8"))
        implementation_digest.update(b"\0")
        implementation_digest.update(content)
        implementation_digest.update(b"\0")
        implementation_hashes.append(
            {
                "path": relative,
                "sha256": _sha256_bytes(content),
                "size_bytes": len(content),
            }
        )
    tracked_diff = _git(root, "diff", "--binary", binary=True)
    if not isinstance(tracked_diff, bytes):
        raise AssertionError("binary git diff unexpectedly returned text")
    changed_behavior = behavior[behavior["exposure_changed"]].copy()
    selected_changed_behavior = changed_behavior[
        changed_behavior["production_selected"].astype(bool)
    ].copy()
    materialization_failure_fields = [
        field
        for field in (
            "materialized_vs_legacy_relation_mismatches",
            "materialized_vs_legacy_status_mismatches",
            "materialized_audit_vs_legacy_status_mismatches",
            "materialized_vs_replayed_selected_row_id_mismatches",
            "materialized_audit_vs_legacy_binding_mismatches",
            "materialized_panel_vs_legacy_exposure_mismatches",
            "materialized_audit_vs_panel_exposure_mismatches",
        )
        if pair_summary[field]
    ]
    if any(pair_summary["materialized_vs_legacy_float_mismatches"].values()):
        materialization_failure_fields.append("materialized_vs_legacy_float_mismatches")
    if expert_reconciliation is not None:
        for field in (
            "materialized_vs_legacy_binding_mismatches",
            "materialized_vs_legacy_exposure_mismatches",
        ):
            if expert_reconciliation[field]:
                materialization_failure_fields.append(f"expert_{field}")
    gate_sections = {
        "materialization_lineage": {
            "status": "passed" if not materialization_failure_fields else "failed",
            "failures": materialization_failure_fields,
        },
        "binding_legacy_compatibility": {
            "status": (
                "passed"
                if not raw_binding_changes
                and not pair_summary["legacy_vs_canonical_binding_mismatches"]
                else "failed"
            ),
            "raw_changes": raw_binding_changes,
            "pair_changes": pair_summary["legacy_vs_canonical_binding_mismatches"],
        },
        "exposure_legacy_compatibility": {
            "status": (
                "passed"
                if not raw_exposure_changes
                and not pair_summary["legacy_vs_canonical_exposure_mismatches"]
                else "failed"
            ),
            "raw_changes": raw_exposure_changes,
            "pair_changes": pair_summary["legacy_vs_canonical_exposure_mismatches"],
        },
        "canonical_exposure_semantic_delta": {
            "status": (
                "review_required"
                if pair_summary["legacy_vs_canonical_exposure_mismatches"]
                else "none"
            ),
            "decision": "not_selected",
        },
        "boundary_contract": {
            "status": "passed" if bool(boundary["passed"].all()) else "failed"
        },
        "production_integration": {
            "status": "blocked" if gate_failures else "eligible"
        },
    }
    artifact_paths = [behavior_path, relation_path, boundary_path, summary_path]
    manifest = {
        "schema_version": 1,
        "status": "passed" if not gate_failures else "failed",
        "phase": "A_standalone_parser_comparison",
        "git_sha": git_sha,
        "git_log": str(_git(root, "log", "-1", "--format=fuller")),
        "git_status_short_at_run": str(_git(root, "status", "--short")),
        "parser_version": SPD_ACTIVITY_PARSER_VERSION,
        "existing_label_policy_versions": {
            "binding": BINDING_LABEL_POLICY_VERSION,
            "exposure": SPD_LABEL_POLICY_VERSION,
        },
        "production_source_hashes": source_hashes,
        "implementation_sources": implementation_hashes,
        "implementation_bundle_sha256": implementation_digest.hexdigest(),
        "tracked_diff_sha256": _sha256_bytes(tracked_diff),
        "inputs": [
            {
                "path": str(path.resolve()),
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in input_paths
        ],
        "row_counts": {
            "raw_workbook_observations": int(len(raw)),
            "normalized_observations": int(len(normalized)),
            "materialized_pairs": int(pair_summary["materialized_pair_rows"]),
            "legacy_pairs": int(pair_summary["legacy_pair_rows"]),
        },
        "relation_counts": relation_counts,
        "raw_label_counts": {
            "legacy_binding": _counts(legacy_binding),
            "canonical_binding": _counts(
                result.numeric_label for result in canonical_binding
            ),
            "legacy_exposure": _counts(legacy_exposure),
            "canonical_exposure": _counts(
                result.numeric_label for result in canonical_exposure
            ),
        },
        "pair_reconciliation": pair_summary,
        "materialized_expert_reconciliation": expert_reconciliation,
        "changed_row_counts": {
            "raw_binding": raw_binding_changes,
            "raw_exposure": raw_exposure_changes,
            "pair_binding": pair_summary["legacy_vs_canonical_binding_mismatches"],
            "pair_exposure": pair_summary["legacy_vs_canonical_exposure_mismatches"],
        },
        "semantic_deltas": {
            "raw_exposure_transitions": _transition_records(changed_behavior),
            "selected_exposure_transitions": _transition_records(
                selected_changed_behavior
            ),
            "affected_source_row_ids": changed_behavior["source_row_id"]
            .astype(str)
            .tolist(),
            "affected_drugs": sorted(
                changed_behavior["drug_identity"].astype(str).unique().tolist()
            ),
            "affected_targets": sorted(
                changed_behavior["target_identity"].astype(str).unique().tolist()
            ),
            "selected_affected_drugs": sorted(
                selected_changed_behavior["drug_identity"].astype(str).unique().tolist()
            ),
            "selected_affected_targets": sorted(
                selected_changed_behavior["target_identity"]
                .astype(str)
                .unique()
                .tolist()
            ),
        },
        "expected_locked_pair_counts": {
            "binding": EXPECTED_LEGACY_BINDING,
            "exposure": EXPECTED_LEGACY_EXPOSURE,
            "basis": "materialized committed legacy production behavior",
        },
        "gate_sections": gate_sections,
        "gate_failures": gate_failures,
        "lineage_profiles": lineage,
        "commands": {
            "argv": [sys.executable, *sys.argv],
            "display_record": args.command_record or "",
        },
        "tests": args.test_record,
        "numeric_comparison": {
            "materialized_float_tolerance": "absolute delta <= 1e-9 * max(abs(left), abs(right), 1)",
            "boundary_classification": "canonical interval comparison; open lower bound at 10 is strictly above 10",
            "unbounded_csv_serialization": "inf",
            "missing_csv_serialization": "blank",
        },
        "environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "pandas": pd.__version__,
            "executable": sys.executable,
        },
        "artifacts": [
            {
                "path": str(path),
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in artifact_paths
        ],
        "timestamp_started_utc": started.isoformat(),
        "timestamp_completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False),
        encoding="utf-8",
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare materialized, committed-legacy, and canonical interval-aware "
            "SPD labels without regenerating production tables."
        )
    )
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--spd-workbook", required=True, type=Path)
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--materialized-selection", required=True, type=Path)
    parser.add_argument("--phase1-source", type=Path)
    parser.add_argument("--historical-master", type=Path)
    parser.add_argument("--historical-binding-raw", type=Path)
    parser.add_argument("--historical-binding-model-ready", type=Path)
    parser.add_argument("--historical-exposure-raw", type=Path)
    parser.add_argument("--historical-exposure-model-ready", type=Path)
    parser.add_argument("--strict-binding-raw", type=Path)
    parser.add_argument("--strict-binding-model-ready", type=Path)
    parser.add_argument("--strict-exposure-raw", type=Path)
    parser.add_argument("--strict-exposure-model-ready", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--expected-git-sha", default=EXPECTED_GIT_SHA)
    parser.add_argument("--command-record")
    parser.add_argument("--test-record", action="append", default=[])
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = run_comparison(args)
    print(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False))
    return 0 if manifest["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
