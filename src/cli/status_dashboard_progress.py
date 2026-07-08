from __future__ import annotations

import datetime as _dt
from collections import Counter
from pathlib import Path
from typing import Any, Mapping

from cli.status_dashboard_utils import (
    STAGE_ORDER,
    as_float,
    as_int,
    chunk_ids_from_json,
    median,
    read_json,
)


def lifecycle(manifest: Mapping[str, Any] | None) -> dict[str, Any]:
    if not isinstance(manifest, Mapping):
        return {"status": "unknown", "timing": {}}
    timing = manifest.get("timing")
    return {
        "status": str(manifest.get("status") or "unknown"),
        "timing": timing if isinstance(timing, Mapping) else {},
    }


def overall_progress(manifest: Mapping[str, Any] | None) -> dict[str, Any]:
    proteins = manifest.get("proteins") if isinstance(manifest, Mapping) else None
    summary = manifest.get("summary") if isinstance(manifest, Mapping) else None
    summary = summary if isinstance(summary, Mapping) else {}
    counter: Counter[str] = Counter()
    failed_examples: list[str] = []
    running_examples: list[str] = []
    if isinstance(proteins, Mapping):
        for key, entry in proteins.items():
            if not isinstance(entry, Mapping):
                continue
            state = str(entry.get("status") or "unknown").lower()
            counter[state] += 1
            if state == "failed" and len(failed_examples) < 5:
                failed_examples.append(str(key))
            if state == "running" and len(running_examples) < 5:
                running_examples.append(str(key))
    scheduled = as_int(summary.get("total_proteins_scheduled"))
    if scheduled is None:
        scheduled = sum(counter.values())
    completed = as_int(summary.get("total_proteins_completed"))
    if completed is None:
        completed = counter.get("completed", 0)
    failed = as_int(summary.get("total_proteins_failed"))
    if failed is None:
        failed = counter.get("failed", 0)
    running = counter.get("running", 0)
    done = int(completed) + int(failed)
    pct = (100.0 * float(done) / float(scheduled)) if scheduled else 0.0
    return {
        "scheduled": int(scheduled or 0),
        "completed": int(completed or 0),
        "failed": int(failed or 0),
        "running": int(running),
        "done": int(done),
        "percent_done": round(pct, 1),
        "protein_status": dict(counter),
        "failed_examples": failed_examples,
        "running_examples": running_examples,
        "summary": dict(summary),
    }


def report_only_lifecycle(report_html: Path) -> dict[str, Any]:
    timing: dict[str, Any] = {}
    try:
        timing["updated_at"] = _dt.datetime.fromtimestamp(
            report_html.stat().st_mtime, _dt.timezone.utc
        ).replace(microsecond=0).isoformat()
    except OSError:
        pass
    return {
        "status": "completed",
        "timing": timing,
        "source": "report_html",
    }


def report_only_progress(report_html: Path) -> dict[str, Any]:
    return {
        "scheduled": 1,
        "completed": 1,
        "failed": 0,
        "running": 0,
        "done": 1,
        "percent_done": 100.0,
        "protein_status": {"completed": 1},
        "failed_examples": [],
        "running_examples": [],
        "summary": {
            "source": "report_html",
            "report_html": str(report_html),
        },
    }


def stage_progress(manifest: Mapping[str, Any] | None) -> dict[str, Any]:
    proteins = manifest.get("proteins") if isinstance(manifest, Mapping) else None
    out: dict[str, Any] = {}
    for stage in STAGE_ORDER:
        counter: Counter[str] = Counter()
        wall_samples: list[float] = []
        failed_examples: list[str] = []
        if isinstance(proteins, Mapping):
            for key, entry in proteins.items():
                if not isinstance(entry, Mapping):
                    continue
                stages = entry.get("stages")
                stage_entry = stages.get(stage) if isinstance(stages, Mapping) else None
                if not isinstance(stage_entry, Mapping):
                    state = "pending"
                else:
                    state = str(stage_entry.get("status") or "pending").lower()
                    timing = stage_entry.get("timing")
                    if isinstance(timing, Mapping):
                        wall = as_float(timing.get("wall_time_sec"))
                        if wall is not None:
                            wall_samples.append(float(wall))
                    if state == "failed" and len(failed_examples) < 5:
                        failed_examples.append(str(key))
                counter[state] += 1
        total = sum(counter.values())
        complete = counter.get("completed", 0) + counter.get("failed", 0)
        out[stage] = {
            "total": total,
            "completed": counter.get("completed", 0),
            "failed": counter.get("failed", 0),
            "running": counter.get("running", 0),
            "pending": counter.get("pending", 0),
            "percent_done": round((100.0 * complete / total) if total else 0.0, 1),
            "status_counts": dict(counter),
            "median_wall_time_sec": median(wall_samples),
            "failed_examples": failed_examples,
        }
    return out


def scorch_failure_summary(manifest: Mapping[str, Any] | None) -> dict[str, Any]:
    proteins = manifest.get("proteins") if isinstance(manifest, Mapping) else None
    errors: Counter[str] = Counter()
    returncodes: Counter[str] = Counter()
    examples: list[dict[str, Any]] = []
    failed = 0
    if not isinstance(proteins, Mapping):
        return {"failed": 0, "errors": {}, "returncodes": {}, "examples": []}
    for key, entry in proteins.items():
        if not isinstance(entry, Mapping):
            continue
        stages = entry.get("stages")
        post = stages.get("postprocessing") if isinstance(stages, Mapping) else None
        if not isinstance(post, Mapping):
            continue
        if str(post.get("status") or "").lower() != "failed":
            continue
        details = post.get("details")
        details = details if isinstance(details, Mapping) else {}
        error = str(post.get("error") or "postprocessing_failed")
        failed += 1
        errors[error] += 1
        returncode = details.get("scorch_returncode")
        if returncode is not None:
            returncodes[str(returncode)] += 1
        if len(examples) < 5:
            timing = post.get("timing")
            timing = timing if isinstance(timing, Mapping) else {}
            examples.append(
                {
                    "combo": str(key),
                    "error": error,
                    "returncode": returncode,
                    "mode": details.get("scorch_mode"),
                    "phase": details.get("scorch_phase"),
                    "device": details.get("scorch_device_effective"),
                    "wall_time_sec": timing.get("wall_time_sec"),
                }
            )
    return {
        "failed": failed,
        "errors": dict(errors),
        "returncodes": dict(returncodes),
        "examples": examples,
    }


def distributed_progress(dist_dir: Path) -> dict[str, Any]:
    out: dict[str, Any] = {
        "enabled": dist_dir.exists(),
        "planned_chunks": 0,
        "completed_chunks": 0,
        "failed_chunks": 0,
        "active_claims": 0,
        "percent_done": 0.0,
        "phase_markers": {},
        "combo_prep": {},
    }
    if not dist_dir.exists():
        return out
    plan_ids: set[str] = set()
    for path in sorted(dist_dir.glob("combo_chunks_*.json")):
        plan_ids.update(chunk_ids_from_json(path))
    result_dir = dist_dir / "chunk_results"
    result_status: Counter[str] = Counter()
    result_times: list[float] = []
    if result_dir.exists():
        for path in result_dir.glob("*.json"):
            payload = read_json(path)
            state = str(payload.get("status") or "completed").lower()
            result_status[state] += 1
            completed_at = as_float(payload.get("completed_at"))
            if completed_at is not None:
                result_times.append(float(completed_at))
    claims_dir = dist_dir / "chunk_claims"
    active_claims = sum(1 for _ in claims_dir.glob("*.json")) if claims_dir.exists() else 0
    phase_dir = dist_dir / "phase_markers"
    phases: dict[str, int] = {}
    if phase_dir.exists():
        for child in phase_dir.iterdir():
            if child.is_dir():
                phases[child.name] = sum(1 for _ in child.glob("*.json"))
    prep_dir = dist_dir / "combo_prep_state"
    prep_counts: Counter[str] = Counter()
    if prep_dir.exists():
        for path in prep_dir.glob("*.json"):
            payload = read_json(path)
            prep_counts[str(payload.get("status") or "unknown").lower()] += 1
    planned = len(plan_ids)
    completed = result_status.get("completed", 0) + result_status.get("ok", 0)
    failed = result_status.get("failed", 0)
    out.update(
        {
            "planned_chunks": planned,
            "completed_chunks": completed,
            "failed_chunks": failed,
            "active_claims": active_claims,
            "result_status_counts": dict(result_status),
            "percent_done": round((100.0 * (completed + failed) / planned) if planned else 0.0, 1),
            "phase_markers": phases,
            "combo_prep": dict(prep_counts),
            "first_chunk_completed_at": min(result_times) if result_times else None,
            "last_chunk_completed_at": max(result_times) if result_times else None,
        }
    )
    return out
