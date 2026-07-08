from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path
from typing import Any, Mapping

from analysis.reporting.manifest_utils import load_run_manifest
from cli.status_dashboard_progress import lifecycle, overall_progress, stage_progress
from cli.status_dashboard_utils import (
    STAGE_ORDER,
    as_float,
    as_int,
    elapsed_seconds,
    mapping,
    median,
)
from config.output_paths import legacy_root, output_root


ETA_HISTORY_NAME = "status_eta_history.json"


def estimate_eta(
    status: Mapping[str, Any], history: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    progress = mapping(status.get("progress"))
    lifecycle_data = mapping(status.get("lifecycle"))
    timing = mapping(lifecycle_data.get("timing"))
    elapsed = elapsed_seconds(timing)
    scheduled = as_int(progress.get("scheduled")) or 0
    done = as_int(progress.get("done")) or 0
    remaining = max(0, scheduled - done)
    eta: dict[str, Any] = {
        "elapsed_sec": elapsed,
        "remaining_units": remaining,
        "eta_sec": None,
        "confidence": "low",
        "basis": "unknown",
    }
    history_eta = _eta_from_history(eta, remaining, history)
    if history_eta is not None:
        return history_eta
    distributed_eta = _eta_from_distributed(eta, mapping(status.get("distributed")), elapsed)
    if distributed_eta is not None:
        return distributed_eta
    return _eta_from_protein_completion(eta, scheduled, done, remaining, elapsed)


def _eta_from_history(
    eta: dict[str, Any],
    remaining: int,
    history: Mapping[str, Any] | None,
) -> dict[str, Any] | None:
    history = history if isinstance(history, Mapping) else {}
    per_protein = history.get("per_protein_wall_sec")
    if remaining <= 0:
        return None
    if not isinstance(per_protein, Mapping):
        return None
    median = per_protein.get("median")
    if not isinstance(median, (int, float)):
        return None
    sample_count = as_int(per_protein.get("sample_count")) or 0
    eta["eta_sec"] = round(float(remaining) * float(median), 1)
    eta["confidence"] = "high" if sample_count >= 10 else "medium"
    eta["basis"] = "history_per_protein"
    eta["history_samples"] = sample_count
    return eta


def _eta_from_distributed(
    eta: dict[str, Any],
    distributed: Mapping[str, Any],
    elapsed: float | None,
) -> dict[str, Any] | None:
    planned_chunks = as_int(distributed.get("planned_chunks")) or 0
    chunk_done = (as_int(distributed.get("completed_chunks")) or 0) + (
        as_int(distributed.get("failed_chunks")) or 0
    )
    if planned_chunks <= 0 or chunk_done <= 0 or not elapsed or elapsed <= 0:
        return None
    rate = float(chunk_done) / elapsed
    eta["remaining_units"] = max(0, planned_chunks - chunk_done)
    eta["eta_sec"] = round(float(eta["remaining_units"]) / rate, 1) if rate else None
    eta["confidence"] = "medium" if chunk_done >= 5 else "low"
    eta["basis"] = "distributed_chunks"
    return eta


def _eta_from_protein_completion(
    eta: dict[str, Any],
    scheduled: int,
    done: int,
    remaining: int,
    elapsed: float | None,
) -> dict[str, Any]:
    if scheduled > 0 and done > 0 and elapsed and elapsed > 0:
        rate = float(done) / elapsed
        eta["eta_sec"] = round(float(remaining) / rate, 1) if rate else None
        eta["confidence"] = "medium" if done >= 3 else "low"
        eta["basis"] = "protein_completion"
    return eta


def refresh_eta_history(root: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    seen_run_ids: set[str] = set()
    manifests_root = output_root(root, "manifests")
    for scan_root in (output_root(root, "manifests"), legacy_root(root, "manifests")):
        manifests_root = scan_root
        if not scan_root.exists():
            continue
        for run_dir in sorted(scan_root.iterdir()):
            if not run_dir.is_dir():
                continue
            if run_dir.name in seen_run_ids:
                continue
            manifest_path = run_dir / "run_manifest.yaml"
            if not manifest_path.exists():
                continue
            manifest, _manifest_path = load_run_manifest(root, run_dir.name)
            row = eta_history_row(run_dir.name, manifest)
            if row:
                rows.append(row)
                seen_run_ids.add(run_dir.name)
    per_protein_samples = [
        float(row["per_protein_wall_sec"])
        for row in rows
        if isinstance(row.get("per_protein_wall_sec"), (int, float))
    ]
    stage_samples: dict[str, list[float]] = {stage: [] for stage in STAGE_ORDER}
    for row in rows:
        stages = row.get("stages")
        if not isinstance(stages, Mapping):
            continue
        for stage in STAGE_ORDER:
            value = stages.get(stage)
            if isinstance(value, (int, float)):
                stage_samples[stage].append(float(value))
    history: dict[str, Any] = {
        "schema_version": 1,
        "updated_at": _dt.datetime.now(_dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "run_count": len(rows),
        "runs": rows,
        "per_protein_wall_sec": {
            "median": median(per_protein_samples),
            "sample_count": len(per_protein_samples),
        },
        "stage_wall_sec": {
            stage: {"median": median(values), "sample_count": len(values)}
            for stage, values in stage_samples.items()
        },
    }
    history_path = manifests_root / ETA_HISTORY_NAME
    try:
        manifests_root.mkdir(parents=True, exist_ok=True)
        history_path.write_text(json.dumps(history, indent=2, sort_keys=True), encoding="utf-8")
    except OSError:
        pass
    return history


def eta_history_row(run_id: str, manifest: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(manifest, Mapping):
        return None
    lifecycle_data = lifecycle(manifest)
    if str(lifecycle_data.get("status") or "").lower() != "completed":
        return None
    progress = overall_progress(manifest)
    scheduled = as_int(progress.get("scheduled")) or 0
    if scheduled <= 0:
        return None
    timing = lifecycle_data.get("timing")
    timing = timing if isinstance(timing, Mapping) else {}
    wall = as_float(timing.get("wall_time_sec"))
    if wall is None or wall <= 0:
        return None
    stages = stage_progress(manifest)
    return {
        "run_id": str(run_id),
        "scheduled": scheduled,
        "wall_time_sec": round(float(wall), 3),
        "per_protein_wall_sec": round(float(wall) / float(scheduled), 3),
        "stages": {
            stage: row.get("median_wall_time_sec")
            for stage, row in stages.items()
            if isinstance(row, Mapping) and isinstance(row.get("median_wall_time_sec"), (int, float))
        },
    }


def eta_history_summary(history: Mapping[str, Any]) -> dict[str, Any]:
    per_protein = history.get("per_protein_wall_sec")
    per_protein = per_protein if isinstance(per_protein, Mapping) else {}
    return {
        "run_count": as_int(history.get("run_count")) or 0,
        "per_protein_wall_sec_median": per_protein.get("median"),
        "per_protein_sample_count": per_protein.get("sample_count", 0),
    }
