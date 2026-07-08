from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

from analysis.reporting.value_utils import normalize_text
from config.output_paths import run_output_dir


_STALE_CACHE_DAYS = 30.0


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _parse_iso(value: str) -> datetime | None:
    text = normalize_text(value)
    if not text:
        return None
    try:
        return datetime.fromisoformat(text.replace("Z", "+00:00"))
    except Exception:
        return None


def _age_days(timestamp: datetime | None) -> float | None:
    if timestamp is None:
        return None
    now = datetime.now(timezone.utc)
    if timestamp.tzinfo is None:
        timestamp = timestamp.replace(tzinfo=timezone.utc)
    return max(0.0, (now - timestamp).total_seconds() / 86400.0)


def _safe_count_from_csv(path: Path, field_names: Sequence[str]) -> int | None:
    if not path.exists():
        return None
    seen: set[str] = set()
    try:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                for field_name in field_names:
                    value = normalize_text(row.get(field_name))
                    if value:
                        seen.add(value)
                        break
    except Exception:
        return None
    return len(seen)


def _cache_freshness(repo_root: Path) -> Dict[str, Any]:
    cache_path = repo_root / "pathways" / "cache" / "target_safety_aggregated.json"
    cache_payload = _load_json(cache_path)
    generated_at_text = normalize_text(cache_payload.get("generated_at"))
    generated_at = _parse_iso(generated_at_text)
    if generated_at is None and cache_path.exists():
        generated_at = datetime.fromtimestamp(cache_path.stat().st_mtime, tz=timezone.utc)
        generated_at_text = generated_at.isoformat()
    age_days = _age_days(generated_at)
    if not cache_path.exists():
        status = "missing"
    elif age_days is None:
        status = "unknown"
    elif age_days <= _STALE_CACHE_DAYS:
        status = "fresh"
    else:
        status = "stale"
    return {
        "path": str(cache_path),
        "generated_at": generated_at_text,
        "age_days": None if age_days is None else round(age_days, 3),
        "status": status,
        "max_age_days": _STALE_CACHE_DAYS,
        "entry_count": len(cache_payload.get("entries") or []),
    }


def _drift_summary(repo_root: Path) -> Dict[str, Any]:
    drift_path = repo_root / "pathways" / "cache" / "target_safety_drift_summary.json"
    drift_payload = _load_json(drift_path)
    return {
        "path": str(drift_path),
        "generated_at": normalize_text(drift_payload.get("generated_at")),
        "changed_count": int(drift_payload.get("changed_count") or 0),
        "changes": drift_payload.get("changes") or [],
    }


def build_heatmap_report_meta(
    *,
    repo_root: Path,
    run_id: str,
    asset_mode: str = "auto",
    rows: Iterable[Dict[str, Any]],
    full_row_labels: Sequence[str],
    full_col_labels: Sequence[str],
    display_top_k: int,
    target_group_meta: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    row_list = list(rows)
    generated_at = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    data_dir = run_output_dir(repo_root, "data", run_id)
    master_rows_path = data_dir / "master_rows.csv"
    full_run_ligands = _safe_count_from_csv(
        master_rows_path,
        ["ligand_base", "ligand_display", "ligand_name"],
    )
    full_run_targets = _safe_count_from_csv(
        master_rows_path,
        ["target_id", "target_label", "pdb_id"],
    )
    unresolved_targets: List[Dict[str, str]] = []
    for label in full_col_labels:
        meta = target_group_meta.get(label) or {}
        pdb_id = normalize_text(meta.get("pdb")).upper()
        if meta.get("uniprot_accessions"):
            continue
        unresolved_targets.append(
            {
                "pdb": pdb_id or normalize_text(meta.get("target_raw")).upper(),
                "target_name": normalize_text(meta.get("target_name")) or normalize_text(label),
            }
        )
    unresolved_targets = sorted(
        unresolved_targets,
        key=lambda item: (item["pdb"].lower(), item["target_name"].lower()),
    )
    coverage_summary = {
        "default_visible_ligands": min(max(1, int(display_top_k)), len(full_row_labels)),
        "interactive_ligands": len(full_row_labels),
        "interactive_targets": len(full_col_labels),
        "full_run_ligands": full_run_ligands or len(full_row_labels),
        "full_run_targets": full_run_targets or len(full_col_labels),
        "matrix_is_subset": bool(full_run_ligands and full_run_ligands > len(full_row_labels)),
        "row_source_rows": len(row_list),
    }
    return {
        "report_summary": {
            "run_id": run_id,
            "generated_at": generated_at,
            "asset_mode": normalize_text(asset_mode) or "unknown",
        },
        "coverage_summary": coverage_summary,
        "validation_summary": {
            "cache": _cache_freshness(repo_root),
            "missing_target_mappings": {
                "count": len(unresolved_targets),
                "targets": unresolved_targets[:16],
            },
            "safety_bucket_drift": _drift_summary(repo_root),
        },
    }
