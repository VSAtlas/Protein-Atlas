from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple


NUMERIC_PRIORITY = [
    "consensus_score",
    "consensus_score_pre",
    "z_vs_decoys_consensus",
    "z_vs_decoys_consensus_pre",
    "z_vs_decoys_blend",
    "t_vs_decoys_consensus",
    "t_vs_decoys_consensus_pre",
    "t_vs_decoys_blend",
    "blend_mu_decoy",
    "blend_sigma_decoy",
    "blend_n_decoys",
    "final_score",
    "scorch_composite",
    "SCORCH_score_used",
    "SCORCH_certainty_used",
    "rescored_flag",
    "consensus_rank",
    "final_rank",
    "p_energy",
    "p_vina",
    "p_gnina_energy",
    "p_ledock",
    "p_dock6",
    "p_cnn",
    "cnn_score_used",
    "cnn_affinity_used",
    "cnn_vs_used",
    "cnn_rescored_flag",
    "scorch_pct",
    "cnn_pct",
    "ml_blend_score",
    "n_engines_with_data",
    "selected_docking_score",
]
TEXT_PRIORITY = [
    "best_engine",
    "best_signal",
    "scorch_source_used",
    "run_mode",
    "selected_stage",
    "rescored_stage",
    "stage_match_flag",
    "stage_fallback_reason",
]


@dataclass(frozen=True)
class ScorchPick:
    pdb_id: str
    variant: str
    ph: str
    source: str
    ligand_id: str
    ligand_base: str
    stage_num: int
    scorch_score: Optional[float]
    scorch_certainty: Optional[float]
    run_mode: str
    selected_stage: str
    rescored_stage: str
    selected_docking_score: Optional[float]
    stage_match_flag: str
    stage_fallback_reason: str


@dataclass(frozen=True)
class CnnPick:
    ligand_base: str
    stage_num: int
    cnn_score: Optional[float]
    cnn_affinity: Optional[float]
    cnn_vs: Optional[float]


def find_consensus_csv(docked_combo_dir: Path) -> Optional[Path]:
    exact = docked_combo_dir / "consensus_docking_scores.csv"
    if exact.exists() and exact.stat().st_size > 0:
        return exact
    candidates = sorted(docked_combo_dir.glob("consensus*.csv")) + sorted(
        docked_combo_dir.glob("*_consensus_docking_scores.csv")
    )
    seen: set[Path] = set()
    for p in candidates:
        if p in seen:
            continue
        seen.add(p)
        name = p.name.lower()
        if p.exists() and p.stat().st_size > 0 and "dud_" not in name:
            return p
    return None


def find_dud_consensus_csv(docked_combo_dir: Path) -> Optional[Path]:
    exact = docked_combo_dir / "dud_consensus_docking_scores.csv"
    if exact.exists() and exact.stat().st_size > 0:
        return exact
    for p in sorted(docked_combo_dir.glob("dud_consensus*.csv")):
        if p.exists() and p.stat().st_size > 0:
            return p
    for p in sorted(docked_combo_dir.glob("*dud*consensus*.csv")):
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def _read_csv(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [dict(row) for row in reader]
        fields = list(reader.fieldnames or [])
    return rows, fields


def _read_consensus_rows_with_chunk_union(
    consensus_csv: Path,
    logger: logging.Logger,
    *,
    component: str,
    consensus_ligand_base: Callable[[str], str],
) -> Tuple[List[Dict[str, str]], List[str]]:
    rows, fields = _read_csv(consensus_csv)
    combo_dir = consensus_csv.parent
    primary_name = consensus_csv.name.lower()

    def _eligible_chunk_consensus(path: Path) -> bool:
        name = path.name.lower()
        if name == primary_name:
            return False
        if "dud_consensus" in name or "rerank" in name or "pretty" in name:
            return False
        return path.exists() and path.stat().st_size > 0

    chunk_files = [
        path
        for path in sorted(combo_dir.glob("*_consensus_docking_scores.csv"))
        if _eligible_chunk_consensus(path)
    ]
    if not chunk_files:
        return rows, fields

    merged: List[Dict[str, str]] = [dict(row) for row in rows]
    merged_fields = list(fields)
    index_by_key: Dict[Tuple[str, str, str, str, str], int] = {}

    def _row_key(row: Dict[str, Any]) -> Tuple[str, str, str, str, str]:
        run_id = str(row.get("run_id", "")).strip()
        pdb_id = str(row.get("pdb_id", "")).strip().upper()
        variant = str(row.get("variant", "")).strip().upper()
        ph_label = str(row.get("ph_label", row.get("ph", ""))).strip()
        lig_base = consensus_ligand_base(str(row.get("ligand", "")).strip())
        return (run_id, pdb_id, variant, ph_label, lig_base)

    for idx, row in enumerate(merged):
        key = _row_key(row)
        if key[-1]:
            index_by_key[key] = idx

    added_rows = 0
    replaced_rows = 0
    for chunk_file in chunk_files:
        try:
            chunk_rows, chunk_fields = _read_csv(chunk_file)
        except Exception as exc:
            logger.warning(
                "%s action=consensus status=degraded reason=chunk_read_failed path=%s error=%s",
                component,
                str(chunk_file),
                exc,
            )
            continue
        for name in chunk_fields:
            if name not in merged_fields:
                merged_fields.append(name)
        for row in chunk_rows:
            key = _row_key(row)
            if not key[-1]:
                continue
            existing_idx = index_by_key.get(key)
            if existing_idx is None:
                merged.append(dict(row))
                index_by_key[key] = len(merged) - 1
                added_rows += 1
                continue
            existing = merged[existing_idx]
            existing_score = str(existing.get("consensus_score", "")).strip()
            incoming_score = str(row.get("consensus_score", "")).strip()
            if not existing_score and incoming_score:
                merged[existing_idx] = {**existing, **dict(row)}
                replaced_rows += 1

    logger.info(
        "%s action=consensus status=ok reason=chunk_union primary=%s chunk_files=%d base_rows=%d added_rows=%d replaced_rows=%d merged_rows=%d",
        component,
        str(consensus_csv),
        len(chunk_files),
        len(rows),
        added_rows,
        replaced_rows,
        len(merged),
    )
    return merged, merged_fields


def _write_csv(
    path: Path,
    rows: List[Dict[str, Any]],
    fieldnames: List[str],
    preamble_lines: Optional[List[str]] = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        if preamble_lines:
            for line in preamble_lines:
                handle.write(line)
                if not line.endswith("\n"):
                    handle.write("\n")
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def _is_numeric_field(
    name: str,
    rows: List[Dict[str, Any]],
    *,
    numeric_priority: List[str],
) -> bool:
    lowered = name.lower()
    if lowered in {value.lower() for value in numeric_priority}:
        return True
    values = [row.get(name, "") for row in rows]
    non_empty = [value for value in values if str(value).strip() != ""]
    if not non_empty:
        return False
    success = 0
    for value in non_empty:
        try:
            float(str(value).strip())
            success += 1
        except Exception:
            continue
    return success >= 0.8 * len(non_empty)


def _ordered_fields(
    original_fields: List[str],
    *,
    numeric_priority: List[str],
    text_priority: List[str],
) -> List[str]:
    base_fields = [
        field
        for field in original_fields
        if field not in {"ligand_base", "scorch_stage_used"}
    ]
    extras = [
        "scorch_source_used",
        "SCORCH_score_used",
        "SCORCH_certainty_used",
        "scorch_composite",
        "selected_stage",
        "rescored_stage",
        "selected_docking_score",
        "stage_match_flag",
        "stage_fallback_reason",
        "cnn_score_used",
        "cnn_affinity_used",
        "cnn_vs_used",
        "cnn_rescored_flag",
        "scorch_pct",
        "cnn_pct",
        "ml_blend_score",
        "library",
        "run_mode",
        "consensus_score_pre",
        "z_vs_decoys_consensus",
        "z_vs_decoys_consensus_pre",
        "t_vs_decoys_consensus",
        "t_vs_decoys_consensus_pre",
        "consensus_mu_decoy",
        "consensus_sigma_decoy",
        "consensus_n_decoys",
        "z_vs_decoys_blend",
        "t_vs_decoys_blend",
        "blend_mu_decoy",
        "blend_sigma_decoy",
        "blend_n_decoys",
        "final_score",
        "rescored_flag",
        "consensus_rank",
        "final_rank",
    ]
    for extra in extras:
        if extra not in base_fields:
            base_fields.append(extra)

    dummy_rows: List[Dict[str, Any]] = []
    numeric_fields = []
    for field in base_fields:
        if _is_numeric_field(field, dummy_rows, numeric_priority=numeric_priority):
            numeric_fields.append(field)

    ordered: List[str] = []
    if "ligand" in base_fields:
        ordered.append("ligand")
    for field in numeric_priority:
        if field in base_fields and field not in ordered:
            ordered.append(field)
    for field in base_fields:
        if field in numeric_fields and field not in ordered:
            ordered.append(field)
    for field in text_priority:
        if field in base_fields and field not in ordered:
            ordered.append(field)
    for field in base_fields:
        if field not in ordered:
            ordered.append(field)
    return ordered


def _sort_rows(
    rows: List[Dict[str, Any]],
    *,
    as_float: Callable[[Any], Optional[float]],
) -> List[Dict[str, Any]]:
    def _rank_key(row: Dict[str, Any]) -> Tuple[int, float, float, float]:
        try:
            final_rank = int(str(row.get("final_rank", "")).strip())
            return (final_rank, 0.0, 0.0, 0.0)
        except Exception:
            comp = as_float(row.get("scorch_composite")) or -1e18
            scs = as_float(row.get("SCORCH_score_used")) or -1e18
            cs = as_float(row.get("consensus_score")) or -1e18
            return (int(1e9), -comp, -scs, -cs)

    return sorted(rows, key=_rank_key)
