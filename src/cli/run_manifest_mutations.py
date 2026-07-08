"""Protein entry and stage mutation helpers for run manifests."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, MutableMapping, Optional

from cli.run_manifest_io import protein_key as _protein_key
from cli.run_manifest_support import (
    STAGE_KEYS,
    coerce_vec3 as _coerce_vec3,
    default_protein_entry as _default_protein_entry,
    default_stage_entry as _default_stage_entry,
    extract_error_from_fail_log as _extract_error_from_fail_log,
    utc_now_iso as _utc_now_iso,
)

def mutate_entry_druggability_and_engine_plan(
    entry: MutableMapping[str, Any],
    *,
    tier: Optional[str],
    use_gnina: bool,
    use_ledock: bool,
    use_dock6: bool,
) -> tuple[Optional[str], Dict[str, bool]]:
    stages = entry.setdefault("stages", {})

    pocket_stage = stages.get("pocket_detection")
    if not isinstance(pocket_stage, MutableMapping):
        pocket_stage = _default_stage_entry()
        stages["pocket_detection"] = pocket_stage

    pocket_details = pocket_stage.get("details")
    if not isinstance(pocket_details, MutableMapping):
        pocket_details = {}
    pocket_stage["details"] = pocket_details

    tier_norm = None
    if tier is not None:
        try:
            t = str(tier).strip().upper()
        except Exception:
            t = ""
        if t in {"A", "B", "C"}:
            tier_norm = t

    if tier_norm is not None:
        pocket_details["druggability_tier"] = tier_norm

    docking_stage = stages.get("docking")
    if not isinstance(docking_stage, MutableMapping):
        docking_stage = _default_stage_entry()
        stages["docking"] = docking_stage

    docking_details = docking_stage.get("details")
    if not isinstance(docking_details, MutableMapping):
        docking_details = {}
    docking_stage["details"] = docking_details

    engine_plan = docking_details.get("engine_plan")
    if not isinstance(engine_plan, MutableMapping):
        engine_plan = {}
    docking_details["engine_plan"] = engine_plan

    if tier_norm is not None:
        engine_plan["tier"] = tier_norm

    engines = {
        "vina": True,
        "gnina": bool(use_gnina),
        "ledock": bool(use_ledock),
        "dock6": bool(use_dock6),
    }
    engine_plan["engines"] = engines
    return tier_norm, engines

def ensure_protein(
    manifest: MutableMapping[str, Any],
    pdb_id: str,
    variant_label: Optional[str],
    ph_tag: Optional[str],
) -> Dict[str, Any]:
    proteins = manifest.setdefault("proteins", {})
    if not isinstance(proteins, MutableMapping):
        manifest["proteins"] = {}
        proteins = manifest["proteins"]

    key = _protein_key(pdb_id, variant_label, ph_tag)
    entry = proteins.get(key)
    if not isinstance(entry, MutableMapping):
        entry = _default_protein_entry()
        proteins[key] = entry
    else:
        # Ensure stages/timing keys exist
        entry.setdefault("pdb_id", str(pdb_id))
        entry.setdefault("library", None)
        entry.setdefault("variant", None)
        entry.setdefault("ph", ph_tag)
        entry.setdefault("status", "pending")
        entry.setdefault("error", None)
        entry.setdefault(
            "timing", {"started_at": None, "finished_at": None, "wall_time_sec": None}
        )
        stages = entry.setdefault("stages", {})
        for stage_key in STAGE_KEYS:
            stage_entry = stages.get(stage_key)
            if not isinstance(stage_entry, MutableMapping):
                stage_entry = _default_stage_entry()
                stages[stage_key] = stage_entry
            stage_entry.setdefault("status", "pending")
            stage_entry.setdefault("error", None)
            stage_entry.setdefault(
                "timing",
                {"started_at": None, "finished_at": None, "wall_time_sec": None},
            )
            stage_entry.setdefault("details", {})

    entry["pdb_id"] = str(pdb_id).upper()
    entry["variant"] = (variant_label or "legacy").strip().upper() or "LEGACY"
    entry["ph"] = ph_tag

    return entry  # type: ignore[return-value]
def mutate_entry_protein_start(
    entry: MutableMapping[str, Any],
    cfg: Mapping[str, Any],
    variant_label: str,
    library: Optional[str],
    ph_tag: Optional[str],
) -> None:
    entry["variant"] = (variant_label or "legacy").upper()
    if library is None:
        library = cfg.get("LIBRARY_SUBDIR_DEFAULT")
    entry["library"] = library
    entry["status"] = "running"
    entry.setdefault("ph", ph_tag)
    timing = entry.get("timing", {})
    if isinstance(timing, MutableMapping):
        if not timing.get("started_at"):
            timing["started_at"] = _utc_now_iso()
    else:
        entry["timing"] = {
            "started_at": _utc_now_iso(),
            "finished_at": None,
            "wall_time_sec": None,
        }

def mutate_entry_protein_success(
    entry: MutableMapping[str, Any],
    elapsed_sec: float,
) -> None:
    now = _utc_now_iso()
    entry["status"] = "completed"
    entry["error"] = None

    timing = entry.get("timing", {})
    if isinstance(timing, MutableMapping):
        timing.setdefault("started_at", now)
        timing["finished_at"] = now
        timing["wall_time_sec"] = round(float(elapsed_sec), 3)
    else:
        entry["timing"] = {
            "started_at": now,
            "finished_at": now,
            "wall_time_sec": round(float(elapsed_sec), 3),
        }

    stages = entry.get("stages", {})
    if isinstance(stages, MutableMapping):
        for stage_key in STAGE_KEYS:
            stage_entry = stages.get(stage_key)
            if not isinstance(stage_entry, MutableMapping):
                stages[stage_key] = _default_stage_entry()
                stage_entry = stages[stage_key]
            if str(stage_entry.get("status") or "").lower() != "failed":
                stage_entry["status"] = "completed"
                stage_timing = ensure_stage_timing(stage_entry)
                if not stage_timing.get("started_at"):
                    stage_timing["started_at"] = now
                if not stage_timing.get("finished_at"):
                    stage_timing["finished_at"] = now
            stage_entry.setdefault("error", None)

def mutate_entry_protein_failure(
    entry: MutableMapping[str, Any],
    variant_label: str,
    fail_log_path: Optional[Path],
) -> None:
    entry["variant"] = (variant_label or "legacy").upper()
    entry["status"] = "failed"
    entry["error"] = _extract_error_from_fail_log(fail_log_path)
    entry["fail_log"] = str(fail_log_path) if fail_log_path is not None else None
def mutate_entry_pocket_detection(
    entry: MutableMapping[str, Any],
    *,
    method: Optional[str],
    center: Optional[Any],
    box_size: Optional[Any],
) -> tuple[Optional[str], Optional[list[float]], Optional[list[float]]]:
    stages = entry.setdefault("stages", {})
    pocket_stage = stages.get("pocket_detection")
    if not isinstance(pocket_stage, MutableMapping):
        pocket_stage = _default_stage_entry()
        stages["pocket_detection"] = pocket_stage

    details = pocket_stage.get("details")
    if not isinstance(details, MutableMapping):
        details = {}
    pocket_stage["details"] = details

    existing_method = details.get("method")
    method_raw = (method or "").strip()
    if not method_raw and existing_method is not None:
        try:
            method_raw = str(existing_method).strip()
        except Exception:
            method_raw = ""
    method_str = method_raw or None
    center_vec = _coerce_vec3(center)
    box_vec = _coerce_vec3(box_size)

    if method_str is not None:
        details["method"] = method_str
    if center_vec is not None:
        details["center"] = [round(float(x), 3) for x in center_vec]
    if box_vec is not None:
        details["box_size"] = [round(float(x), 1) for x in box_vec]

    if method_str and center_vec and box_vec:
        pocket_stage.setdefault("status", "completed")
        pocket_stage.setdefault("error", None)
    return method_str, center_vec, box_vec
def ensure_stage_timing(
    stage_entry: MutableMapping[str, Any],
) -> MutableMapping[str, Any]:
    timing = stage_entry.get("timing")
    if not isinstance(timing, MutableMapping):
        timing = {"started_at": None, "finished_at": None, "wall_time_sec": None}
    else:
        timing.setdefault("started_at", None)
        timing.setdefault("finished_at", None)
        timing.setdefault("wall_time_sec", None)
    stage_entry["timing"] = timing
    return timing


def mutate_entry_stage_status(
    entry: MutableMapping[str, Any],
    *,
    stage_key: str,
    status: str,
    elapsed_sec: Optional[float],
    details: Optional[Mapping[str, Any]],
    error: Optional[str],
) -> bool:
    if stage_key not in STAGE_KEYS:
        return False
    now = _utc_now_iso()
    stages = entry.setdefault("stages", {})
    stage_entry = stages.get(stage_key)
    if not isinstance(stage_entry, MutableMapping):
        stage_entry = _default_stage_entry()
        stages[stage_key] = stage_entry

    timing = ensure_stage_timing(stage_entry)
    stage_details = stage_entry.get("details")
    if not isinstance(stage_details, MutableMapping):
        stage_details = {}
    stage_entry["details"] = stage_details

    status_token = str(status or "").strip().lower()
    if status_token == "running":
        stage_entry["status"] = "running"
        stage_entry["error"] = None
        if not timing.get("started_at"):
            timing["started_at"] = now
        if details:
            for key, value in details.items():
                stage_details[str(key)] = value
        return True

    if status_token in {"completed", "failed"}:
        if not timing.get("started_at"):
            timing["started_at"] = now
        timing["finished_at"] = now
        if elapsed_sec is not None:
            timing["wall_time_sec"] = round(float(elapsed_sec), 3)
        if details:
            for key, value in details.items():
                stage_details[str(key)] = value
        stage_entry["status"] = status_token
        if error is not None:
            stage_entry["error"] = str(error)
        elif status_token == "completed":
            stage_entry["error"] = None
        return True

    return False
def mutate_entry_docking_overall(
    entry: MutableMapping[str, Any],
    *,
    event: str,
    elapsed_sec: Optional[float],
    error: Optional[str],
) -> bool:
    now = _utc_now_iso()
    stages = entry.setdefault("stages", {})
    docking_stage = stages.get("docking")
    if not isinstance(docking_stage, MutableMapping):
        docking_stage = _default_stage_entry()
        stages["docking"] = docking_stage

    timing = ensure_stage_timing(docking_stage)

    if event == "start":
        docking_stage["status"] = "running"
        docking_stage["error"] = None
        if not timing.get("started_at"):
            timing["started_at"] = now
        return True

    if event in ("end", "fail"):
        if not timing.get("started_at"):
            timing["started_at"] = now
        timing["finished_at"] = now
        if elapsed_sec is not None:
            timing["wall_time_sec"] = round(float(elapsed_sec), 3)
        docking_stage["status"] = "failed" if event == "fail" else "completed"
        if error is not None:
            docking_stage["error"] = str(error)
        return True

    return False
def mutate_entry_docking_stage(
    entry: MutableMapping[str, Any],
    *,
    raw_name: str,
    subrun_label: str,
    base_stage_name: str,
    status: str,
    elapsed_sec: Optional[float],
    error: Optional[str],
) -> bool:
    now = _utc_now_iso()
    stages = entry.setdefault("stages", {})
    docking_stage = stages.get("docking")
    if not isinstance(docking_stage, MutableMapping):
        docking_stage = _default_stage_entry()
        stages["docking"] = docking_stage

    details = docking_stage.get("details")
    if not isinstance(details, MutableMapping):
        details = {}
    docking_stage["details"] = details

    per_stage = details.get("per_stage")
    if not isinstance(per_stage, MutableMapping):
        per_stage = {}
    details["per_stage"] = per_stage

    stage_entry = per_stage.get(raw_name)
    if not isinstance(stage_entry, MutableMapping):
        stage_entry = {
            "status": "pending",
            "error": None,
            "timing": {
                "started_at": None,
                "finished_at": None,
                "wall_time_sec": None,
            },
            "subrun": subrun_label,
            "stage_base_name": base_stage_name,
        }
    per_stage[raw_name] = stage_entry
    stage_entry["subrun"] = subrun_label
    stage_entry["stage_base_name"] = base_stage_name

    timing = stage_entry.get("timing")
    if not isinstance(timing, MutableMapping):
        timing = {"started_at": None, "finished_at": None, "wall_time_sec": None}
    else:
        timing.setdefault("started_at", None)
        timing.setdefault("finished_at", None)
        timing.setdefault("wall_time_sec", None)
    stage_entry["timing"] = timing

    by_subrun = details.get("by_subrun")
    if not isinstance(by_subrun, MutableMapping):
        by_subrun = {}
    details["by_subrun"] = by_subrun

    subrun_map = by_subrun.get(subrun_label)
    if not isinstance(subrun_map, MutableMapping):
        subrun_map = {}
    by_subrun[subrun_label] = subrun_map
    subrun_map[raw_name] = stage_entry

    if status == "running":
        if not timing.get("started_at"):
            timing["started_at"] = now
        stage_entry["status"] = "running"
        stage_entry["error"] = None
        return True

    if status in ("completed", "failed"):
        if not timing.get("started_at"):
            timing["started_at"] = now
        timing["finished_at"] = now
        if elapsed_sec is not None:
            timing["wall_time_sec"] = round(float(elapsed_sec), 3)
        stage_entry["status"] = status
        if error is not None:
            stage_entry["error"] = str(error)
        return True

    return False
def mutate_entry_postprocessing(
    entry: MutableMapping[str, Any],
    *,
    status: str,
    elapsed_sec: Optional[float],
    details: Optional[Mapping[str, Any]],
    error: Optional[str],
) -> bool:
    now = _utc_now_iso()
    stages = entry.setdefault("stages", {})
    post_stage = stages.get("postprocessing")
    if not isinstance(post_stage, MutableMapping):
        post_stage = _default_stage_entry()
        stages["postprocessing"] = post_stage

    timing = ensure_stage_timing(post_stage)
    stage_details = post_stage.get("details")
    if not isinstance(stage_details, MutableMapping):
        stage_details = {}
    post_stage["details"] = stage_details

    if status == "running":
        post_stage["status"] = "running"
        post_stage["error"] = None
        if not timing.get("started_at"):
            timing["started_at"] = now
        return True

    if status in ("completed", "failed"):
        if not timing.get("started_at"):
            timing["started_at"] = now
        timing["finished_at"] = now
        if elapsed_sec is not None:
            elapsed = round(float(elapsed_sec), 3)
            timing["wall_time_sec"] = elapsed
            stage_details["scorch_wall_time_sec"] = elapsed
        if details:
            for key, value in details.items():
                stage_details[str(key)] = value
        post_stage["status"] = status
        if error is not None:
            post_stage["error"] = str(error)
        elif status == "completed":
            post_stage["error"] = None
        return True

    return False
