from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from analysis.reporting.value_utils import (
    normalize_side_effect_label,
    normalize_text,
    normalized_side_effect_set,
)


def _entries_by_pdb(payload: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    entries = payload.get("entries") or []
    by_pdb: Dict[str, Dict[str, Any]] = {}
    if not isinstance(entries, list):
        return by_pdb
    for raw in entries:
        if not isinstance(raw, dict):
            continue
        pdb_id = normalize_text(raw.get("pdb_id")).upper()
        if pdb_id:
            by_pdb[pdb_id] = raw
    return by_pdb


def build_target_safety_drift_summary(
    previous_payload: Dict[str, Any] | None,
    current_payload: Dict[str, Any],
) -> Dict[str, Any]:
    previous = _entries_by_pdb(previous_payload or {})
    current = _entries_by_pdb(current_payload)
    all_pdbs = sorted(set(previous) | set(current))
    changes: List[Dict[str, Any]] = []
    for pdb_id in all_pdbs:
        old = previous.get(pdb_id) or {}
        new = current.get(pdb_id) or {}
        old_bucket = (
            normalize_side_effect_label(old.get("primary_display_safety")) or "Unassigned"
        )
        new_bucket = (
            normalize_side_effect_label(new.get("primary_display_safety")) or "Unassigned"
        )
        old_conf = normalize_text(old.get("safety_confidence")) or "unassigned"
        new_conf = normalize_text(new.get("safety_confidence")) or "unassigned"
        old_sources = sorted(normalized_side_effect_set(old.get("safety_sources") or []))
        new_sources = sorted(normalized_side_effect_set(new.get("safety_sources") or []))
        if (
            old_bucket == new_bucket
            and old_conf == new_conf
            and old_sources == new_sources
        ):
            continue
        changes.append(
            {
                "pdb_id": pdb_id,
                "target_name": normalize_text(new.get("target_name"))
                or normalize_text(old.get("target_name")),
                "old_primary_display_safety": old_bucket,
                "new_primary_display_safety": new_bucket,
                "old_safety_confidence": old_conf,
                "new_safety_confidence": new_conf,
                "old_safety_sources": old_sources,
                "new_safety_sources": new_sources,
            }
        )
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "previous_generated_at": normalize_text((previous_payload or {}).get("generated_at")),
        "current_generated_at": normalize_text(current_payload.get("generated_at")),
        "previous_entry_count": len(previous),
        "current_entry_count": len(current),
        "changed_count": len(changes),
        "changes": changes,
    }


def write_target_safety_drift_summary(
    repo_root: Path,
    payload: Dict[str, Any],
    out_path: Optional[Path] = None,
) -> Path:
    path = out_path or (repo_root / "pathways" / "cache" / "target_safety_drift_summary.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path
