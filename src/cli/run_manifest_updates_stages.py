"""Run-manifest stage and pocket update entrypoints."""

from __future__ import annotations

import json
import logging
from dataclasses import asdict
from typing import Any, Mapping, MutableMapping, Optional

from cli.run_manifest_distributed import distributed_mutate_protein_entry as _distributed_mutate_protein_entry
from cli.run_manifest_io import (
    get_manifest_paths,
    load_manifest as _load_manifest,
    manifest_lock_path as _manifest_lock_path,
    write_manifest as _write_manifest,
)
from cli.run_manifest_mutations import (
    ensure_protein as _ensure_protein,
    mutate_entry_docking_overall as _mutate_entry_docking_overall,
    mutate_entry_docking_stage as _mutate_entry_docking_stage,
    mutate_entry_pocket_detection as _mutate_entry_pocket_detection,
    mutate_entry_postprocessing as _mutate_entry_postprocessing,
    mutate_entry_stage_status as _mutate_entry_stage_status,
)
from cli.run_manifest_support import (
    PocketDetectionEvent,
    default_stage_entry as _default_stage_entry,
    distributed_enabled as _distributed_enabled,
    exclusive_file_lock as _exclusive_file_lock,
    refresh_summary as _refresh_summary,
)

def emit_pocket_detection_event(
    cfg: Mapping[str, Any], event: PocketDetectionEvent
) -> None:
    """
    Append a pocket-detection event to a JSONL log for later replay by the
    main process. Best-effort only; failures are logged as warnings.
    """
    ph_label = event.ph_tag if event.ph_tag is not None else "base"
    try:
        if _distributed_enabled(cfg):
            update_manifest_for_pocket_detection(
                cfg,
                run_id=event.run_id,
                pdb_id=event.pdb_id,
                variant_label=event.variant_label,
                ph_tag=event.ph_tag,
                method=event.method,
                center=event.center,
                box_size=event.box_size,
            )
            return

        if not event.run_id:
            logging.debug(
                "[run-manifest.pocket_detection.event.skip] reason=missing_run_id pdb=%s variant=%s ph=%s",
                event.pdb_id,
                event.variant_label,
                ph_label,
            )
            return

        _, manifest_path = get_manifest_paths(cfg, event.run_id)
        events_path = manifest_path.with_name("run_manifest_events.jsonl")
        payload = {"type": "pocket_detection", **asdict(event)}
        events_path.parent.mkdir(parents=True, exist_ok=True)
        with events_path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(payload, sort_keys=True) + "\n")

        logging.debug(
            "[run-manifest.pocket_detection.event] run_id=%s pdb=%s variant=%s ph=%s method=%s",
            event.run_id,
            event.pdb_id,
            event.variant_label,
            ph_label,
            event.method,
        )
    except Exception:
        logging.warning(
            "[run-manifest.pocket_detection.event.error] run_id=%s pdb=%s variant=%s ph=%s",
            event.run_id,
            event.pdb_id,
            event.variant_label,
            ph_label,
            exc_info=True,
        )


def apply_pocket_detection_events(cfg: Mapping[str, Any], run_id: str) -> None:
    """
    Replay pocket-detection events for a run and write them to the manifest in
    a single-threaded main-process context. Best-effort; errors are logged and
    ignored.
    """
    if not run_id:
        return
    if _distributed_enabled(cfg):
        return

    try:
        _, manifest_path = get_manifest_paths(cfg, run_id)
        events_path = manifest_path.with_name("run_manifest_events.jsonl")
        if not events_path.exists() or events_path.stat().st_size == 0:
            return

        try:
            with events_path.open("r", encoding="utf-8") as fh:
                lines = fh.readlines()
        except Exception:
            logging.warning(
                "[run-manifest.pocket_detection.events.read.error] run_id=%s path=%s",
                run_id,
                events_path,
                exc_info=True,
            )
            return

        events: list[dict[str, Any]] = []
        for raw in lines:
            line = raw.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except Exception:
                continue
            if payload.get("type") != "pocket_detection":
                continue
            if str(payload.get("run_id", "")) != str(run_id):
                continue
            events.append(payload)

        if not events:
            try:
                events_path.unlink()
            except Exception:
                pass
            return

        for payload in events:
            try:
                pdb_token = payload.get("pdb_id")
                if not pdb_token:
                    continue
                update_manifest_for_pocket_detection(
                    cfg,
                    run_id=str(payload.get("run_id") or run_id),
                    pdb_id=str(pdb_token),
                    variant_label=payload.get("variant_label"),
                    ph_tag=payload.get("ph_tag"),
                    method=payload.get("method"),
                    center=payload.get("center"),
                    box_size=payload.get("box_size"),
                )
            except Exception:
                logging.warning(
                    "[run-manifest.pocket_detection.events.apply.error] run_id=%s pdb=%s variant=%s ph=%s",
                    run_id,
                    payload.get("pdb_id"),
                    payload.get("variant_label"),
                    payload.get("ph_tag"),
                    exc_info=True,
                )

        try:
            events_path.unlink()
        except Exception:
            try:
                events_path.write_text("", encoding="utf-8")
            except Exception:
                pass
    except Exception:
        logging.warning(
            "[run-manifest.pocket_detection.events.error] run_id=%s",
            run_id,
            exc_info=True,
        )


def update_manifest_for_pocket_detection(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    *,
    ph_tag: Optional[str],
    method: Optional[str],
    center: Optional[Any],
    box_size: Optional[Any],
) -> None:
    """
    Record pocket-detection method + geometry under the pocket_detection stage
    for a given (pdb, variant, pH) context.

    Best-effort only: any error logs a WARNING and is otherwise ignored.
    """
    ph_label = ph_tag if ph_tag is not None else "base"
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.pocket_detection.skip] no run_id pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_label,
                ph_label,
            )
            return

        logging.debug(
            "[run-manifest.pocket_detection.request] run_id=%s pdb=%s variant=%s ph=%s method=%s center=%r box=%r",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            method,
            center,
            box_size,
        )

        if _distributed_mutate_protein_entry(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            lambda entry: _mutate_entry_pocket_detection(
                entry, method=method, center=center, box_size=box_size
            ),
        ):
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        with _exclusive_file_lock(_manifest_lock_path(manifest_path)):
            manifest = _load_manifest(manifest_path)
            if manifest is None:
                logging.warning(
                    "[run-manifest.pocket_detection.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s",
                    run_id,
                    manifest_path,
                    pdb_id,
                    variant_label,
                    ph_label,
                )
                return

            entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
            method_str, center_vec, box_vec = _mutate_entry_pocket_detection(
                entry, method=method, center=center, box_size=box_size
            )
            stages = entry.setdefault("stages", {})
            pocket_stage = stages.get("pocket_detection")
            if (
                isinstance(pocket_stage, MutableMapping)
                and method_str
                and center_vec
                and box_vec
            ):
                pocket_stage.setdefault("status", "completed")
                pocket_stage.setdefault("error", None)
                details = pocket_stage.get("details")
                if not isinstance(details, Mapping):
                    details = {}
                logging.info(
                    "[run-manifest.pocket_detection.ok] run_id=%s pdb=%s variant=%s ph=%s method=%s center=%s box=%s",
                    run_id,
                    pdb_id,
                    variant_label,
                    ph_label,
                    method_str,
                    details.get("center"),
                    details.get("box_size"),
                )
            else:
                logging.info(
                    "[run-manifest.pocket_detection.partial] run_id=%s pdb=%s variant=%s ph=%s method=%s center=%r box=%r",
                    run_id,
                    pdb_id,
                    variant_label,
                    ph_label,
                    method_str,
                    center,
                    box_size,
                )

            _refresh_summary(manifest)
            _write_manifest(manifest_path, manifest)

    except Exception:
        logging.warning(
            "[run-manifest.pocket_detection.error] run_id=%s pdb=%s variant=%s ph=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            exc_info=True,
        )

def update_manifest_for_stage(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    stage_key: str,
    status: str,
    *,
    ph_tag: Optional[str] = None,
    elapsed_sec: Optional[float] = None,
    details: Optional[Mapping[str, Any]] = None,
    error: Optional[str] = None,
) -> None:
    ph_label = ph_tag if ph_tag is not None else "base"
    stage_token = str(stage_key or "").strip()
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.stage.skip] no run_id pdb=%s variant=%s ph=%s stage=%s status=%s",
                pdb_id,
                variant_label,
                ph_label,
                stage_token,
                status,
            )
            return

        if _distributed_mutate_protein_entry(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            lambda entry: _mutate_entry_stage_status(
                entry,
                stage_key=stage_token,
                status=status,
                elapsed_sec=elapsed_sec,
                details=details,
                error=error,
            ),
        ):
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        with _exclusive_file_lock(_manifest_lock_path(manifest_path)):
            manifest = _load_manifest(manifest_path)
            if manifest is None:
                logging.warning(
                    "[run-manifest.stage.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s stage=%s",
                    run_id,
                    manifest_path,
                    pdb_id,
                    variant_label,
                    ph_label,
                    stage_token,
                )
                return
            entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
            applied = _mutate_entry_stage_status(
                entry,
                stage_key=stage_token,
                status=status,
                elapsed_sec=elapsed_sec,
                details=details,
                error=error,
            )
            if not applied:
                logging.debug(
                    "[run-manifest.stage.skip] unknown stage/status run_id=%s pdb=%s variant=%s ph=%s stage=%s status=%s",
                    run_id,
                    pdb_id,
                    variant_label,
                    ph_label,
                    stage_token,
                    status,
                )
                return
            _refresh_summary(manifest)
            _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.stage.error] run_id=%s pdb=%s variant=%s ph=%s stage=%s status=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            stage_token,
            status,
            exc_info=True,
        )


def update_manifest_for_prep_stage(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    status: str,
    *,
    ph_tag: Optional[str] = None,
    elapsed_sec: Optional[float] = None,
    details: Optional[Mapping[str, Any]] = None,
    error: Optional[str] = None,
) -> None:
    update_manifest_for_stage(
        cfg,
        run_id,
        pdb_id,
        variant_label,
        "prep",
        status,
        ph_tag=ph_tag,
        elapsed_sec=elapsed_sec,
        details=details,
        error=error,
    )
def update_manifest_for_docking_overall(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    *,
    ph_tag: Optional[str] = None,
    event: str = "start",
    elapsed_sec: Optional[float] = None,
    error: Optional[str] = None,
) -> None:
    """
    Record overall docking timing + status for a single (pdb, variant, pH).

    event == "start" : mark running + capture started_at if unset.
    event == "end"   : mark completed + capture finished_at/wall_time_sec.
    event == "fail"  : mark failed + capture finished_at/wall_time_sec/error.
    """
    ph_label = ph_tag if ph_tag is not None else "base"
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.docking-overall.skip] no run_id pdb=%s variant=%s ph=%s event=%s",
                pdb_id,
                variant_label,
                ph_label,
                event,
            )
            return

        if _distributed_mutate_protein_entry(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            lambda entry: _mutate_entry_docking_overall(
                entry,
                event=event,
                elapsed_sec=elapsed_sec,
                error=error,
            ),
        ):
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        with _exclusive_file_lock(_manifest_lock_path(manifest_path)):
            manifest = _load_manifest(manifest_path)
            if manifest is None:
                logging.warning(
                    "[run-manifest.docking-overall.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s",
                    run_id,
                    manifest_path,
                    pdb_id,
                    variant_label,
                    ph_label,
                )
                return

            entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
            stages = entry.setdefault("stages", {})
            docking_stage = stages.get("docking")
            if not isinstance(docking_stage, MutableMapping):
                docking_stage = _default_stage_entry()
                stages["docking"] = docking_stage

            applied = _mutate_entry_docking_overall(
                entry,
                event=event,
                elapsed_sec=elapsed_sec,
                error=error,
            )
            if not applied:
                logging.debug(
                    "[run-manifest.docking-overall.skip] unknown_event=%s run_id=%s pdb=%s variant=%s ph=%s",
                    event,
                    run_id,
                    pdb_id,
                    variant_label,
                    ph_label,
                )
                return

            _refresh_summary(manifest)
            _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.docking-overall.error] run_id=%s pdb=%s variant=%s ph=%s event=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            event,
            exc_info=True,
        )

def update_manifest_for_docking_stage(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    stage_name: str,
    status: str,
    *,
    ph_tag: Optional[str] = None,
    elapsed_sec: Optional[float] = None,
    error: Optional[str] = None,
) -> None:
    """
    Record status/timing for a single docking stage (e.g., "stage1", "stage2").
    """
    ph_label = ph_tag if ph_tag is not None else "base"
    stage_token = (stage_name or "stage").strip() or "stage"

    # --- normalize stage token into subrun + engine + base stage name ---
    lib_prefix = None
    core = stage_token
    if core.startswith("dud_"):
        lib_prefix = "dud"
        core = core[len("dud_") :] or core
    elif core.startswith("hmdb_"):
        lib_prefix = "hmdb"
        core = core[len("hmdb_") :] or core

    engine = None
    if core.startswith("vina_"):
        engine = "vina"
        core = core[len("vina_") :] or core
    elif core.startswith("gnina_"):
        engine = "gnina"
        core = core[len("gnina_") :] or core
    elif core.startswith("ledock_"):
        engine = "ledock"
        core = core[len("ledock_") :] or core

    # Handle engine-prefixed strings that still carry a library prefix (e.g., gnina_dud_stage1)
    if lib_prefix is None:
        if core.startswith("dud_"):
            lib_prefix = "dud"
            core = core[len("dud_") :] or core
        elif core.startswith("hmdb_"):
            lib_prefix = "hmdb"
            core = core[len("hmdb_") :] or core

    base_stage_name = core or stage_token
    engine = engine or "vina"
    subrun_label = lib_prefix if lib_prefix else "primary"

    if lib_prefix:
        raw_name = f"{lib_prefix}_{engine}_{base_stage_name}"
    else:
        raw_name = f"{engine}_{base_stage_name}"
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.docking-stage.skip] no run_id pdb=%s variant=%s ph=%s stage=%s",
                pdb_id,
                variant_label,
                ph_label,
                stage_token,
            )
            return

        if _distributed_mutate_protein_entry(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            lambda entry: _mutate_entry_docking_stage(
                entry,
                raw_name=raw_name,
                subrun_label=subrun_label,
                base_stage_name=base_stage_name,
                status=status,
                elapsed_sec=elapsed_sec,
                error=error,
            ),
        ):
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        with _exclusive_file_lock(_manifest_lock_path(manifest_path)):
            manifest = _load_manifest(manifest_path)
            if manifest is None:
                logging.warning(
                    "[run-manifest.docking-stage.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s stage=%s",
                    run_id,
                    manifest_path,
                    pdb_id,
                    variant_label,
                    ph_label,
                    stage_token,
                )
                return

            entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
            stages = entry.setdefault("stages", {})
            docking_stage = stages.get("docking")
            if not isinstance(docking_stage, MutableMapping):
                docking_stage = _default_stage_entry()
                stages["docking"] = docking_stage

            applied = _mutate_entry_docking_stage(
                entry,
                raw_name=raw_name,
                subrun_label=subrun_label,
                base_stage_name=base_stage_name,
                status=status,
                elapsed_sec=elapsed_sec,
                error=error,
            )
            if not applied:
                logging.debug(
                    "[run-manifest.docking-stage.skip] unknown_status=%s run_id=%s pdb=%s variant=%s ph=%s stage=%s",
                    status,
                    run_id,
                    pdb_id,
                    variant_label,
                    ph_label,
                    stage_token,
                )
                return

            _refresh_summary(manifest)
            _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.docking-stage.error] run_id=%s pdb=%s variant=%s ph=%s stage=%s status=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            stage_token,
            status,
            exc_info=True,
        )
def update_manifest_for_postprocessing(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    *,
    ph_tag: Optional[str] = None,
    status: str = "running",
    elapsed_sec: Optional[float] = None,
    details: Optional[Mapping[str, Any]] = None,
    error: Optional[str] = None,
) -> None:
    ph_label = ph_tag if ph_tag is not None else "base"
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.postprocessing.skip] no run_id pdb=%s variant=%s ph=%s status=%s",
                pdb_id,
                variant_label,
                ph_label,
                status,
            )
            return

        if _distributed_mutate_protein_entry(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            lambda entry: _mutate_entry_postprocessing(
                entry,
                status=status,
                elapsed_sec=elapsed_sec,
                details=details,
                error=error,
            ),
        ):
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        with _exclusive_file_lock(_manifest_lock_path(manifest_path)):
            manifest = _load_manifest(manifest_path)
            if manifest is None:
                logging.warning(
                    "[run-manifest.postprocessing.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s",
                    run_id,
                    manifest_path,
                    pdb_id,
                    variant_label,
                    ph_label,
                )
                return

            entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
            applied = _mutate_entry_postprocessing(
                entry,
                status=status,
                elapsed_sec=elapsed_sec,
                details=details,
                error=error,
            )
            if not applied:
                logging.debug(
                    "[run-manifest.postprocessing.skip] unknown_status=%s run_id=%s pdb=%s variant=%s ph=%s",
                    status,
                    run_id,
                    pdb_id,
                    variant_label,
                    ph_label,
                )
                return

            _refresh_summary(manifest)
            _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.postprocessing.error] run_id=%s pdb=%s variant=%s ph=%s status=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            status,
            exc_info=True,
        )
