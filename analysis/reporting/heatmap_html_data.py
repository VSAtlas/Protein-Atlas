from __future__ import annotations

import csv
import logging
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, cast

from analysis.reporting.heatmap_html_scale import _is_finite
from analysis.reporting.value_utils import (
    as_float,
    is_truthy_matrix_cell,
    median_floats,
    normalize_text,
)
from analysis.target_ids import build_target_id, pdb_id_from_target_id

_LOG = logging.getLogger("heatmap-html")

def _sample_names(values: List[str], max_items: int = 5, max_len: int = 40) -> List[str]:
    sample: List[str] = []
    for value in values[:max_items]:
        text = str(value)
        if len(text) > max_len:
            text = text[: max_len - 3] + "..."
        sample.append(text)
    return sample


def _parse_fraction_config(value: Optional[str], default_value: float) -> float:
    if value is None:
        return default_value
    text = str(value).strip()
    if not text:
        return default_value
    try:
        if text.endswith("%"):
            parsed = float(text[:-1].strip()) / 100.0
        else:
            parsed = float(text)
    except ValueError:
        return default_value
    if not math.isfinite(parsed) or parsed < 0:
        return default_value
    return parsed


def _promiscuity_pose_valid(value: Any) -> bool:
    normalized = normalize_text(value)
    if normalized == "":
        return True
    return is_truthy_matrix_cell(value)


def _parse_float_or_none(value: Any) -> Optional[float]:
    return as_float(value)


def _normalized_entropy(weights: List[float]) -> float:
    k = len(weights)
    if k <= 1:
        return 0.0
    positive_weights = [max(float(weight), 0.0) for weight in weights]
    total = sum(positive_weights)
    if total <= 0:
        probs = [1.0 / k] * k
    else:
        probs = [weight / total for weight in positive_weights]
    entropy = 0.0
    for prob in probs:
        if prob > 0:
            entropy -= prob * math.log(prob)
    return entropy / math.log(k)


def _robust_z_scores(counts_by_name: Dict[str, int]) -> Dict[str, float]:
    if not counts_by_name:
        return {}
    values = [float(value) for value in counts_by_name.values()]
    median_value = median_floats(values)
    if median_value is None:
        return {name: 0.0 for name in counts_by_name}
    deviations = [abs(value - median_value) for value in values]
    mad = median_floats(deviations) or 0.0
    if mad > 0:
        return {
            name: 0.67449 * (float(count) - median_value) / mad
            for name, count in counts_by_name.items()
        }
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / len(values)
    std = math.sqrt(variance)
    if std > 0:
        return {
            name: (float(count) - mean_value) / std
            for name, count in counts_by_name.items()
        }
    return {name: 0.0 for name in counts_by_name}


def _promiscuity_hits(
    entries: Iterable[Tuple[str, Dict[str, Any]]],
    *,
    pct_threshold: float,
) -> List[Tuple[str, float]]:
    hits: List[Tuple[str, float]] = []
    for partner_name, cell in entries:
        pct_rank = as_float(cell.get("pct_rank"))
        if pct_rank is None or pct_rank > pct_threshold:
            continue
        if not _promiscuity_pose_valid(cell.get("pose_valid_any")):
            continue
        z_selected = as_float(cell.get("z_selected"))
        hits.append((partner_name, z_selected if z_selected is not None else 0.0))
    return hits


def _build_promiscuity_meta(
    agg: Dict[Tuple[str, str], Dict[str, Any]],
    *,
    pct_threshold: float,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, float]]]:
    ligand_entries: Dict[str, List[Tuple[str, Dict[str, Any]]]] = {}
    target_entries: Dict[str, List[Tuple[str, Dict[str, Any]]]] = {}
    all_ligands = sorted({ligand for ligand, _target in agg.keys()})
    all_targets = sorted({target for _ligand, target in agg.keys()})
    for (ligand, target), cell in agg.items():
        ligand_entries.setdefault(ligand, []).append((target, cell))
        target_entries.setdefault(target, []).append((ligand, cell))

    ligand_counts = {
        ligand: len(_promiscuity_hits(ligand_entries.get(ligand, []), pct_threshold=pct_threshold))
        for ligand in all_ligands
    }
    target_counts = {
        target: len(_promiscuity_hits(target_entries.get(target, []), pct_threshold=pct_threshold))
        for target in all_targets
    }
    ligand_z = _robust_z_scores(ligand_counts)
    target_z = _robust_z_scores(target_counts)

    total_targets = len(all_targets)
    total_ligands = len(all_ligands)
    ligand_meta: Dict[str, Dict[str, float]] = {}
    for ligand in all_ligands:
        hits = _promiscuity_hits(ligand_entries.get(ligand, []), pct_threshold=pct_threshold)
        weights = [weight for _partner, weight in hits]
        hit_count = len(hits)
        ligand_meta[ligand] = {
            "promiscuity_z": ligand_z.get(ligand, 0.0),
            "significant_hits": float(hit_count),
            "partner_total": float(total_targets),
            "hit_fraction": (float(hit_count) / total_targets) if total_targets else 0.0,
            "entropy": _normalized_entropy(weights),
            "threshold_pct_rank": pct_threshold,
        }

    target_meta: Dict[str, Dict[str, float]] = {}
    for target in all_targets:
        hits = _promiscuity_hits(target_entries.get(target, []), pct_threshold=pct_threshold)
        weights = [weight for _partner, weight in hits]
        hit_count = len(hits)
        target_meta[target] = {
            "promiscuity_z": target_z.get(target, 0.0),
            "significant_hits": float(hit_count),
            "partner_total": float(total_ligands),
            "hit_fraction": (float(hit_count) / total_ligands) if total_ligands else 0.0,
            "entropy": _normalized_entropy(weights),
            "threshold_pct_rank": pct_threshold,
        }
    return ligand_meta, target_meta


def _load_heatmap_rows_csv(
    input_csv: Path, include_decoys: bool, allowed_pdb_ids: Optional[Set[str]] = None
) -> List[Dict[str, Any]]:
    normalized_allowed_pdb_ids: Optional[Set[str]] = None
    if allowed_pdb_ids is not None:
        normalized_allowed_pdb_ids = {
            str(pdb_id).strip().upper() for pdb_id in allowed_pdb_ids if str(pdb_id).strip()
        }

    with input_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        selected_col = "z_selected" if "z_selected" in fieldnames else ""
        if not selected_col:
            raise ValueError("Missing required column: z_selected")
        has_target_id = "target_id" in fieldnames
        if not has_target_id:
            required_cols = {"pdb_id", "variant", "ph_label"}
            missing = required_cols - set(fieldnames)
            if missing:
                raise ValueError(
                    f"Missing required columns: {', '.join(sorted(missing))}"
                )
        has_ligand_display = "ligand_display" in fieldnames
        has_ligand_base = "ligand_base" in fieldnames
        if not has_ligand_display and not has_ligand_base:
            raise ValueError("Missing ligand_display/ligand_base columns in input CSV")

        rows: List[Dict[str, Any]] = []
        for row in reader:
            z_val = as_float(row.get(selected_col))
            if z_val is None:
                continue
            if not include_decoys and "is_decoy" in fieldnames:
                decoy_flag = is_truthy_matrix_cell(row.get("is_decoy"))
                if decoy_flag:
                    continue

            ligand_display = normalize_text(row.get("ligand_display")) if has_ligand_display else ""
            ligand_base = normalize_text(row.get("ligand_base")) if has_ligand_base else ""
            if not ligand_display or ligand_display.lower() == "nan":
                ligand_display = ligand_base
            ligand_name = ligand_display
            if not ligand_name:
                continue

            if has_target_id:
                target_id = normalize_text(row.get("target_id"))
                if not target_id and {"pdb_id", "variant", "ph_label"}.issubset(fieldnames):
                    target_id = build_target_id(
                        row.get("pdb_id"),
                        row.get("variant"),
                        row.get("ph_label"),
                    )
            else:
                target_id = build_target_id(
                    row.get("pdb_id"),
                    row.get("variant"),
                    row.get("ph_label"),
                )
            if not target_id:
                continue

            pdb_id = normalize_text(row.get("pdb_id")).upper() if "pdb_id" in fieldnames else ""
            if not pdb_id:
                pdb_id = pdb_id_from_target_id(target_id)
            if (
                normalized_allowed_pdb_ids is not None
                and pdb_id not in normalized_allowed_pdb_ids
            ):
                continue

            row["z_selected"] = z_val
            row["ligand_name"] = ligand_name
            row["target_id"] = target_id
            row["pdb_id"] = pdb_id
            rows.append(row)

    if not rows:
        raise ValueError("No usable z_selected values found in input CSV")
    return rows


def _load_heatmap_rows_parquet(
    input_path: Path, include_decoys: bool, allowed_pdb_ids: Optional[Set[str]] = None
) -> List[Dict[str, Any]]:
    normalized_allowed_pdb_ids: Optional[Set[str]] = None
    if allowed_pdb_ids is not None:
        normalized_allowed_pdb_ids = {
            str(pdb_id).strip().upper() for pdb_id in allowed_pdb_ids if str(pdb_id).strip()
        }

    try:
        import pyarrow.dataset as ds  # type: ignore[import-untyped]
    except Exception as exc:  # pragma: no cover - import error path
        raise RuntimeError(
            "pyarrow is required to read parquet heatmap inputs"
        ) from exc

    dataset = ds.dataset(str(input_path), format="parquet", partitioning="hive")
    schema_names = set(dataset.schema.names)
    selected_col = "z_selected" if "z_selected" in schema_names else ""
    if not selected_col:
        raise ValueError("Missing required column: z_selected")

    has_target_id = "target_id" in schema_names
    has_pdb_id = "pdb_id" in schema_names
    has_target_components = {"pdb_id", "variant", "ph_label"}.issubset(schema_names)
    if not has_target_id and not has_target_components:
        raise ValueError(
            "Missing required columns: target_id or pdb_id/variant/ph_label"
        )

    has_ligand_display = "ligand_display" in schema_names
    has_ligand_base = "ligand_base" in schema_names
    if not has_ligand_display and not has_ligand_base:
        raise ValueError("Missing ligand_display/ligand_base columns in input parquet")

    scan_columns: List[str] = [selected_col]
    for col in ("target_id", "pdb_id", "variant", "ph_label"):
        if col in schema_names and col not in scan_columns:
            scan_columns.append(col)
    if "target_name" in schema_names and "target_name" not in scan_columns:
        scan_columns.append("target_name")
    for col in (
        "ligand_display",
        "ligand_base",
        "is_decoy",
        "rank",
        "pct_rank",
        "pose_valid_any",
        "pose_invalid_reason_top",
        "library",
    ):
        if col in schema_names and col not in scan_columns:
            scan_columns.append(col)

    rows: List[Dict[str, Any]] = []
    scanner = dataset.scanner(columns=scan_columns, use_threads=True)
    for batch in scanner.to_batches():
        payload = batch.to_pydict()
        n_rows = batch.num_rows
        for idx in range(n_rows):
            z_val = as_float(payload[selected_col][idx])
            if z_val is None:
                continue

            if not include_decoys and "is_decoy" in payload:
                if is_truthy_matrix_cell(payload["is_decoy"][idx]):
                    continue

            ligand_display = (
                normalize_text(payload["ligand_display"][idx])
                if has_ligand_display
                else ""
            )
            ligand_base = (
                normalize_text(payload["ligand_base"][idx]) if has_ligand_base else ""
            )
            if not ligand_display or ligand_display.lower() == "nan":
                ligand_display = ligand_base
            ligand_name = ligand_display
            if not ligand_name:
                continue

            if has_target_id:
                target_id = normalize_text(payload["target_id"][idx])
                if not target_id and has_target_components:
                    target_id = build_target_id(
                        payload["pdb_id"][idx],
                        payload["variant"][idx],
                        payload["ph_label"][idx],
                    )
            else:
                target_id = build_target_id(
                    payload["pdb_id"][idx],
                    payload["variant"][idx],
                    payload["ph_label"][idx],
                )
            if not target_id:
                continue

            pdb_id = normalize_text(payload["pdb_id"][idx]).upper() if has_pdb_id else ""
            if not pdb_id:
                pdb_id = pdb_id_from_target_id(target_id)
            if (
                normalized_allowed_pdb_ids is not None
                and pdb_id not in normalized_allowed_pdb_ids
            ):
                continue

            row: Dict[str, Any] = {
                "z_selected": z_val,
                "ligand_name": ligand_name,
                "target_id": target_id,
                "pdb_id": pdb_id,
                "ligand_display": ligand_display,
                "ligand_base": ligand_base,
            }
            for col in (
                "rank",
                "pct_rank",
                "pose_valid_any",
                "pose_invalid_reason_top",
                "library",
                "target_name",
            ):
                if col in payload:
                    row[col] = payload[col][idx]
            rows.append(row)

    if not rows:
        raise ValueError("No usable z_selected values found in input parquet")
    return rows


def _load_heatmap_rows(
    input_path: Path,
    include_decoys: bool,
    allowed_pdb_ids: Optional[Set[str]] = None,
) -> List[Dict[str, Any]]:
    is_parquet_source = input_path.is_dir() or input_path.suffix.lower() == ".parquet"
    if is_parquet_source:
        return _load_heatmap_rows_parquet(
            input_path, include_decoys, allowed_pdb_ids=allowed_pdb_ids
        )
    return _load_heatmap_rows_csv(
        input_path, include_decoys, allowed_pdb_ids=allowed_pdb_ids
    )


def _aggregate_rows(
    rows: Iterable[Dict[str, Any]]
) -> Tuple[Dict[Tuple[str, str], Dict[str, Any]], Dict[str, float]]:
    agg: Dict[Tuple[str, str], Dict[str, Any]] = {}
    ligand_max: Dict[str, float] = {}
    for row in rows:
        ligand_name = row.get("ligand_name") or ""
        target_id = row.get("target_id") or ""
        target_label = row.get("target_label") or target_id
        z_val = row.get("z_selected")
        if not _is_finite(cast(Optional[float], z_val)):
            continue
        z_val = cast(float, z_val)
        key = (ligand_name, target_label)
        existing = agg.get(key)
        if existing is None or z_val > existing["z_selected"]:
            agg[key] = {
                "z_selected": z_val,
                "ligand_name": ligand_name,
                "target_id": target_label,
                "target_raw": normalize_text(target_id),
                "rank": normalize_text(row.get("rank")),
                "pct_rank": normalize_text(row.get("pct_rank")),
                "pose_valid_any": normalize_text(row.get("pose_valid_any")),
                "pose_invalid_reason_top": normalize_text(row.get("pose_invalid_reason_top")),
                "library": normalize_text(row.get("library")),
            }
        current_max = ligand_max.get(ligand_name)
        if current_max is None or z_val > current_max:
            ligand_max[ligand_name] = z_val
    return agg, ligand_max


def _build_matrix(
    agg: Dict[Tuple[str, str], Dict[str, Any]],
    ligand_order: List[str],
    target_order: List[str],
) -> List[List[Optional[float]]]:
    matrix: List[List[Optional[float]]] = []
    for ligand in ligand_order:
        row = []
        for target in target_order:
            entry = agg.get((ligand, target))
            row.append(entry["z_selected"] if entry else None)
        matrix.append(row)
    return matrix


def _filter_empty(
    matrix: List[List[Optional[float]]],
    row_labels: List[str],
    col_labels: List[str],
    *,
    drop_empty_cols: bool = True,
) -> Tuple[List[List[Optional[float]]], List[str], List[str]]:
    row_keep = [any(_is_finite(v) for v in row) for row in matrix]
    if drop_empty_cols:
        col_keep = []
        for col_idx in range(len(col_labels)):
            col_keep.append(
                any(
                    _is_finite(matrix[row_idx][col_idx])
                    for row_idx in range(len(row_labels))
                )
            )
    else:
        col_keep = [True for _ in col_labels]

    filtered_rows = [row for keep, row in zip(row_keep, matrix) if keep]
    filtered_row_labels = [label for keep, label in zip(row_keep, row_labels) if keep]
    if not filtered_row_labels:
        return [], [], []

    filtered_matrix: List[List[Optional[float]]] = []
    for row in filtered_rows:
        filtered_matrix.append([v for keep, v in zip(col_keep, row) if keep])
    filtered_col_labels = [label for keep, label in zip(col_keep, col_labels) if keep]
    return filtered_matrix, filtered_row_labels, filtered_col_labels
