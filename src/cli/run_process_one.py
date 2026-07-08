from __future__ import annotations

import logging
import os
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Optional

from cli.run_context import ConfigDict


_CHUNK_FAILURE_REASON_KEY = "_LAST_CHUNK_FAILURE_REASON"
_CHUNK_FAILURE_RETRYABLE_KEY = "_LAST_CHUNK_FAILURE_RETRYABLE"


def _chunk_failure_retryable(reason: str) -> bool:
    token = str(reason or "").strip().lower()
    if not token:
        return False
    if token.startswith("retryable_"):
        return True
    return (
        "missing_stage1_chunk_marker" in token
        and "missing_docking_summary" in token
    )


def _set_chunk_failure_metadata(
    cfg: ConfigDict,
    *,
    reason: str,
    retryable: bool,
) -> None:
    cfg[_CHUNK_FAILURE_REASON_KEY] = str(reason or "")
    cfg[_CHUNK_FAILURE_RETRYABLE_KEY] = bool(retryable)


def _as_int(value: Any) -> int:
    try:
        return int(value or 0)
    except Exception:
        return 0


def _as_bool(value: Any) -> bool:
    token = str(value or "").strip().lower()
    if token in {"1", "true", "yes", "on", "enabled"}:
        return True
    if token in {"0", "false", "no", "off", "disabled", "none", "null", ""}:
        return False
    return bool(value)


def _resume_manifest_entry(
    cfg: Mapping[str, Any],
    *,
    pdb_id: str,
    variant_label: str,
    ph_token: str,
) -> Mapping[str, Any] | None:
    proteins = cfg.get("_RESUME_MANIFEST_PROTEINS")
    if not isinstance(proteins, Mapping):
        return None
    variant = str(variant_label or "legacy").strip().upper() or "LEGACY"
    ph = str(ph_token or "base").strip() or "base"
    for key in (
        f"{str(pdb_id).upper()}|{variant}|{ph}",
        f"{str(pdb_id).upper()}|LEGACY|{ph}",
        f"{str(pdb_id).upper()}|BASE|{ph}",
    ):
        entry = proteins.get(key)
        if isinstance(entry, Mapping):
            return entry
    return None


def _stage_details(entry: Mapping[str, Any], *stage_names: str) -> Mapping[str, Any]:
    stages = entry.get("stages")
    if not isinstance(stages, Mapping):
        return {}
    for stage_name in stage_names:
        stage = stages.get(stage_name)
        if not isinstance(stage, Mapping):
            continue
        details = stage.get("details")
        if isinstance(details, Mapping):
            return details
    return {}


def _resume_docking_outputs_complete(entry: Mapping[str, Any]) -> bool:
    details = _stage_details(entry, "distributed_chunk_coverage", "docking")
    planned = _as_int(details.get("planned_chunks"))
    completed = _as_int(details.get("completed_chunks"))
    failed = _as_int(details.get("failed_chunks"))
    missing = _as_int(details.get("missing_result_chunks"))
    return (
        planned > 0
        and completed >= planned
        and failed == 0
        and missing == 0
        and _as_bool(details.get("has_consensus_csv"))
    )


def _resume_scorch_still_pending(entry: Mapping[str, Any]) -> bool:
    post = entry.get("stages", {})
    post_stage = post.get("postprocessing") if isinstance(post, Mapping) else None
    post_details = (
        post_stage.get("details")
        if isinstance(post_stage, Mapping)
        and isinstance(post_stage.get("details"), Mapping)
        else {}
    )
    post_status = (
        str(post_stage.get("status") or "").strip().lower()
        if isinstance(post_stage, Mapping)
        else ""
    )
    return post_status != "completed" or not _as_bool(
        post_details.get("has_scorch_reranked_csv")
    )


def _resume_scorch_required(cfg: Mapping[str, Any], entry: Mapping[str, Any]) -> bool:
    details = _stage_details(entry, "distributed_chunk_coverage", "docking")
    if _as_bool(details.get("require_scorch")):
        return True
    return _as_bool(cfg.get("USE_SCORCH")) or _as_bool(cfg.get("SCORCH_ENABLED"))


def _resume_allowed_count_hint(entry: Mapping[str, Any]) -> int | None:
    details = _stage_details(entry, "distributed_chunk_coverage", "docking")
    completed_ligands = _as_int(details.get("completed_ligands"))
    planned_ligands = _as_int(details.get("planned_ligands"))
    value = completed_ligands if completed_ligands > 0 else planned_ligands
    return value if value > 0 else None


def _should_resume_scorch_only(
    cfg: Mapping[str, Any],
    *,
    pdb_id: str,
    variant_label: str,
    ph_token: str,
) -> tuple[bool, Mapping[str, Any] | None]:
    entry = _resume_manifest_entry(
        cfg,
        pdb_id=pdb_id,
        variant_label=variant_label,
        ph_token=ph_token,
    )
    if entry is None:
        return False, None
    if not _resume_scorch_required(cfg, entry):
        return False, entry
    if not _resume_scorch_still_pending(entry):
        return False, entry
    return _resume_docking_outputs_complete(entry), entry


@dataclass
class ProcessOneContext:
    run_id: str
    label: str
    mode: str
    variant: Optional[str]
    is_resume: bool
    completed_combo_lookup: set[tuple[str, str, str]]
    dist_combo_chunk_mode: bool
    cfg_v: ConfigDict
    stages: Any
    params: Any
    tokens: list[str]
    global_start: float


@dataclass
class ProcessOneSharedState:
    failed_root: str
    failed_entries: list[tuple[str, str, str, str, str]]
    run_scope_completion_times: dict[tuple[str, str, str], float]
    retention_lock: Any
    pending_retention_combos: set[tuple[str, str, str]]
    pending_coverage_refresh: dict[tuple[str, str, str], dict[str, Any]]
    pending_retention_lock: Any
    scorch_queue_service: Any


@dataclass
class ProcessOneDeps:
    normalize_ph_tag_token: Callable[[Optional[str]], str]
    chunk_ligand_key: Callable[[str], str]
    process_one_protein: Callable[[ConfigDict, str, Any, Any], None]
    update_manifest_for_protein_start: Callable[..., Any]
    update_manifest_for_protein_success: Callable[..., Any]
    update_manifest_for_protein_failure: Callable[..., Any]
    verify_chunk_combo_outputs: Callable[..., tuple[bool, str]]
    resolve_combo_output_dir: Callable[..., Path]
    resolve_combo_docking_summary_csv: Callable[[Path, str], Path]
    load_scored_ligand_keys_from_summary: Callable[[Path], set[str]]
    update_combo_coverage_snapshot: Callable[..., Any]
    resolve_combo_post_consensus_csv: Callable[[Path, str], Path]
    load_post_scored_ligand_keys: Callable[[Path], set[str]]
    acquire_global_slot: Callable[..., Any]
    maybe_run_scorch_rescore_for_pdb: Callable[..., Any]
    maybe_run_artifact_retention_for_combo: Callable[..., Any]


def build_process_one_runner(
    *,
    context: ProcessOneContext,
    shared_state: ProcessOneSharedState,
    deps: ProcessOneDeps,
) -> Callable[..., bool]:
    def _process_one(
        pdb_file: str,
        cfg_for_pdb: ConfigDict,
        ph_tag: Optional[str],
        *,
        chunk_ligand_bases: Optional[list[str]] = None,
        chunk_id: Optional[str] = None,
        chunk_run_mode: Optional[str] = None,
        chunk_library_name: Optional[str] = None,
        chunk_library_root: Optional[str] = None,
    ) -> bool:
        pdb_id = os.path.splitext(os.path.basename(pdb_file))[0].upper()
        pdb_start = time.time()
        ph_token = deps.normalize_ph_tag_token(ph_tag)
        is_chunk_limited_work = chunk_ligand_bases is not None or bool(chunk_id)
        cfg_for_pdb.pop(_CHUNK_FAILURE_REASON_KEY, None)
        cfg_for_pdb.pop(_CHUNK_FAILURE_RETRYABLE_KEY, None)

        if context.is_resume:
            resume_key = (pdb_id, context.label, ph_token)
            if resume_key in context.completed_combo_lookup:
                logging.info(
                    "[resume.skip] pdb_id=%s variant=%s ph=%s already completed; skipping.",
                    pdb_id,
                    context.label,
                    ph_token,
                )
                return True

        try:
            library_name = str(chunk_library_name or "").strip() or None
            chunk_keys: set[str] = set()
            try:
                if not library_name and "dud" in context.tokens:
                    lib_map = cfg_for_pdb.get("_TEST_LIBRARY_CANONICAL", {}) or {}
                    library_name = lib_map.get(pdb_id)
                if not library_name:
                    library_name = cfg_for_pdb.get("LIBRARY_SUBDIR_DEFAULT")
            except Exception:
                if not library_name:
                    library_name = cfg_for_pdb.get("LIBRARY_SUBDIR_DEFAULT")

            scorch_only, resume_entry = _should_resume_scorch_only(
                cfg_for_pdb,
                pdb_id=pdb_id,
                variant_label=context.label,
                ph_token=ph_token,
            )
            if (
                context.is_resume
                and not is_chunk_limited_work
                and scorch_only
                and resume_entry is not None
            ):
                allowed_hint = _resume_allowed_count_hint(resume_entry)
                logging.info(
                    "[resume.scorch-only] pdb_id=%s variant=%s ph=%s reason=docking_chunks_complete allowed_count_hint=%s",
                    pdb_id,
                    context.label,
                    ph_token,
                    allowed_hint,
                )
                cfg_for_combo = cfg_for_pdb.copy()
                if ph_tag is not None:
                    cfg_for_combo["_PH_TAG_OVERRIDE"] = str(ph_tag)
                else:
                    cfg_for_combo.pop("_PH_TAG_OVERRIDE", None)
                queued_scorch = False
                if shared_state.scorch_queue_service is not None:
                    try:
                        queued_scorch = bool(
                            shared_state.scorch_queue_service.submit(
                                cfg=context.cfg_v,
                                pdb_id=pdb_id,
                                variant=None
                                if context.variant is None
                                else str(context.variant),
                                ph=ph_tag,
                                allowed_count_hint=allowed_hint,
                            )
                        )
                    except Exception:
                        queued_scorch = False
                        logging.warning(
                            "[resume.scorch-only] action=queue_submit_failed pdb_id=%s variant=%s ph=%s",
                            pdb_id,
                            context.label,
                            ph_token,
                            exc_info=True,
                        )
                if not queued_scorch:
                    try:
                        with deps.acquire_global_slot(
                            cfg_for_combo,
                            cores=1,
                            min_cores=1,
                            priority=2,
                            task_id=f"main:scorch-resume:pdb:{pdb_id}:{context.label}:{ph_token}",
                            task_type="scorch_score_chunk",
                        ):
                            deps.maybe_run_scorch_rescore_for_pdb(
                                context.cfg_v,
                                context.run_id,
                                pdb_id,
                                variant=None
                                if context.variant is None
                                else str(context.variant),
                                ph=ph_tag,
                                allowed_count_hint=allowed_hint,
                                verbose=False,
                            )
                    except Exception:
                        logging.warning(
                            "[resume.scorch-only] action=direct_hook_failed pdb_id=%s variant=%s ph=%s",
                            pdb_id,
                            context.label,
                            ph_token,
                            exc_info=True,
                        )
                if shared_state.scorch_queue_service is not None:
                    with shared_state.pending_retention_lock:
                        shared_state.pending_retention_combos.add(
                            (
                                str(pdb_id).upper(),
                                str(context.label or "").strip().upper() or "BASE",
                                deps.normalize_ph_tag_token(ph_tag),
                            )
                        )
                return True

            try:
                deps.update_manifest_for_protein_start(
                    cfg_for_pdb,
                    context.run_id,
                    pdb_id,
                    context.label,
                    library_name,
                    ph_tag=ph_tag,
                )
            except Exception:
                logging.warning(
                    "Failed to update run_manifest for start of %s (%s,%s)",
                    pdb_id,
                    context.label,
                    ph_token,
                    exc_info=True,
                )

            # Main per-combo work item (PDB + variant + pH).
            cfg_for_combo = cfg_for_pdb.copy()
            if ph_tag is not None:
                cfg_for_combo["_PH_TAG_OVERRIDE"] = str(ph_tag)
            else:
                cfg_for_combo.pop("_PH_TAG_OVERRIDE", None)
            if chunk_ligand_bases is not None:
                cfg_for_combo["_CHUNK_LIGAND_KEYS"] = list(chunk_ligand_bases)
                cfg_for_combo["_CHUNK_ID"] = str(chunk_id or "")
                if chunk_run_mode:
                    cfg_for_combo["_CHUNK_RUN_MODE"] = str(chunk_run_mode)
                else:
                    cfg_for_combo.pop("_CHUNK_RUN_MODE", None)
                if chunk_library_name:
                    cfg_for_combo["_CHUNK_LIBRARY_NAME"] = str(chunk_library_name)
                else:
                    cfg_for_combo.pop("_CHUNK_LIBRARY_NAME", None)
                if chunk_library_root:
                    cfg_for_combo["_CHUNK_LIBRARY_ROOT"] = str(chunk_library_root)
                else:
                    cfg_for_combo.pop("_CHUNK_LIBRARY_ROOT", None)
            else:
                cfg_for_combo.pop("_CHUNK_LIGAND_KEYS", None)
                cfg_for_combo.pop("_CHUNK_ID", None)
                cfg_for_combo.pop("_CHUNK_RUN_MODE", None)
                cfg_for_combo.pop("_CHUNK_LIBRARY_NAME", None)
                cfg_for_combo.pop("_CHUNK_LIBRARY_ROOT", None)
            deps.process_one_protein(cfg_for_combo, pdb_file, context.stages, context.params)
            try:
                softfail_count = int(cfg_for_combo.get("_SOFTFAIL_LIGAND_COUNT", 0) or 0)
            except Exception:
                softfail_count = 0
            if softfail_count > 0:
                prev_soft = int(cfg_for_pdb.get("_SOFTFAIL_LIGAND_COUNT", 0) or 0)
                cfg_for_pdb["_SOFTFAIL_LIGAND_COUNT"] = int(prev_soft + softfail_count)
                reasons_raw = cfg_for_combo.get("_SOFTFAIL_LIGAND_REASONS")
                if isinstance(reasons_raw, list):
                    reasons_existing = cfg_for_pdb.get("_SOFTFAIL_LIGAND_REASONS")
                    if not isinstance(reasons_existing, list):
                        reasons_existing = []
                    merged_reasons = [str(x) for x in reasons_existing if str(x).strip()]
                    for reason in reasons_raw:
                        reason_token = str(reason).strip()
                        if reason_token and reason_token not in merged_reasons:
                            merged_reasons.append(reason_token)
                    cfg_for_pdb["_SOFTFAIL_LIGAND_REASONS"] = merged_reasons
                logging.warning(
                    "[chunk.ligand-softfail] pdb=%s variant=%s ph=%s count=%d reasons=%s",
                    pdb_id,
                    context.label,
                    ph_token,
                    softfail_count,
                    ",".join(
                        str(x)
                        for x in cfg_for_pdb.get("_SOFTFAIL_LIGAND_REASONS", [])
                        if str(x).strip()
                    )
                    or "none",
                )
            # Preserve combo-prep reuse metadata discovered during owner prep so
            # waiting chunk workers can skip repeated control redock work.
            for prep_key in (
                "_COMBO_PREP_CENTER_BY_PH",
                "_COMBO_PREP_BOX_BY_PH",
                "_COMBO_PREP_SOURCE_BY_PH",
                "_COMBO_PREP_CONTROL_STEMS",
                "_COMBO_PREP_CONTROL_LOOKUP",
            ):
                prep_val = cfg_for_combo.get(prep_key)
                if prep_val:
                    cfg_for_pdb[prep_key] = prep_val
            if chunk_ligand_bases is not None:
                planned_chunk_keys = {
                    deps.chunk_ligand_key(str(x))
                    for x in chunk_ligand_bases
                    if str(x).strip()
                }
                effective_chunk_keys_raw = cfg_for_combo.get("_CHUNK_EFFECTIVE_LIGAND_KEYS")
                if isinstance(effective_chunk_keys_raw, (list, tuple, set)):
                    effective_chunk_keys = {
                        deps.chunk_ligand_key(str(x))
                        for x in effective_chunk_keys_raw
                        if str(x).strip()
                    }
                else:
                    effective_chunk_keys = set()
                chunk_keys = effective_chunk_keys if effective_chunk_keys else planned_chunk_keys
                if (
                    planned_chunk_keys
                    and effective_chunk_keys
                    and planned_chunk_keys != effective_chunk_keys
                ):
                    logging.info(
                        "[distributed.chunk-filter.effective] chunk_id=%s pdb=%s variant=%s ph=%s planned=%d effective=%d",
                        str(chunk_id or ""),
                        pdb_id,
                        context.label,
                        ph_token,
                        len(planned_chunk_keys),
                        len(effective_chunk_keys),
                    )
                verify_chunk_bases = sorted(chunk_keys) if chunk_keys else list(chunk_ligand_bases)
                chunk_ok, chunk_reason = deps.verify_chunk_combo_outputs(
                    cfg_for_combo,
                    run_id=str(context.run_id),
                    chunk_id=str(chunk_id or ""),
                    pdb_id=str(pdb_id),
                    variant_label=str(context.label),
                    ph_tag=ph_tag,
                    library_name=str(library_name or ""),
                    chunk_ligand_bases=verify_chunk_bases,
                    run_mode=str(chunk_run_mode or ""),
                )
                if not chunk_ok:
                    chunk_failure_reason = (
                        f"chunk_output_verification_failed:{chunk_reason}"
                    )
                    retryable_chunk_failure = _chunk_failure_retryable(chunk_reason)
                    _set_chunk_failure_metadata(
                        cfg_for_pdb,
                        reason=chunk_failure_reason,
                        retryable=retryable_chunk_failure,
                    )
                    if retryable_chunk_failure:
                        logging.warning(
                            "[distributed.chunk.verify.retryable] chunk_id=%s pdb=%s variant=%s ph=%s reason=%s",
                            str(chunk_id or ""),
                            pdb_id,
                            context.label,
                            ph_token,
                            chunk_failure_reason,
                        )
                        return False
                    raise RuntimeError(chunk_failure_reason)
                if chunk_reason != "ok":
                    logging.warning(
                        "[distributed.chunk.verify.degraded] chunk_id=%s pdb=%s variant=%s ph=%s reason=%s",
                        str(chunk_id or ""),
                        pdb_id,
                        context.label,
                        ph_token,
                        chunk_reason,
                    )
                combo_docked_dir = deps.resolve_combo_output_dir(
                    cfg_for_combo,
                    root_key="DOCKED_DIR",
                    default_name="docked",
                    run_id=str(context.run_id),
                    pdb_id=str(pdb_id),
                    variant_label=str(context.label),
                    ph_tag=ph_tag,
                )
                summary_csv = deps.resolve_combo_docking_summary_csv(
                    combo_docked_dir,
                    str(library_name or ""),
                )
                scored_keys = deps.load_scored_ligand_keys_from_summary(summary_csv)
                deps.update_combo_coverage_snapshot(
                    cfg_for_combo,
                    run_id=str(context.run_id),
                    pdb_id=str(pdb_id),
                    variant_label=str(context.label),
                    ph_tag=ph_tag,
                    library_name=str(library_name or ""),
                    expected_keys=chunk_keys,
                    docking_keys=(scored_keys & chunk_keys if chunk_keys else scored_keys),
                )
            # Queue SCORCH as mixed global work; fallback to legacy per-PDB direct path.
            # For distributed combo chunks, defer to one owner task per combo after
            # all chunks are complete to avoid duplicate tail submissions.
            defer_chunk_scorch = bool(
                context.dist_combo_chunk_mode
                and chunk_ligand_bases is not None
                and shared_state.scorch_queue_service is not None
            )
            if not defer_chunk_scorch:
                allowed_hint = (
                    max(1, len(chunk_ligand_bases))
                    if isinstance(chunk_ligand_bases, list) and chunk_ligand_bases
                    else None
                )
                queued_scorch = False
                if shared_state.scorch_queue_service is not None:
                    try:
                        queued_scorch = bool(
                            shared_state.scorch_queue_service.submit(
                                cfg=context.cfg_v,
                                pdb_id=pdb_id,
                                variant=None if context.variant is None else str(context.variant),
                                ph=ph_tag,
                                allowed_count_hint=allowed_hint,
                            )
                        )
                    except Exception:
                        queued_scorch = False
                        logging.warning(
                            "[scorch-rescore.queue] action=submit_failed pdb_id=%s variant=%s ph=%s",
                            pdb_id,
                            context.label,
                            ph_token,
                            exc_info=True,
                        )
                if not queued_scorch:
                    try:
                        with deps.acquire_global_slot(
                            cfg_for_combo,
                            cores=1,
                            min_cores=1,
                            priority=2,
                            task_id=f"main:scorch:pdb:{pdb_id}:{context.label}:{ph_token}",
                            task_type="scorch_score_chunk",
                        ):
                            deps.maybe_run_scorch_rescore_for_pdb(
                                context.cfg_v,
                                context.run_id,
                                pdb_id,
                                variant=None if context.variant is None else str(context.variant),
                                ph=ph_tag,
                                allowed_count_hint=allowed_hint,
                                verbose=False,
                            )
                            if chunk_ligand_bases is not None:
                                combo_post_dir = deps.resolve_combo_output_dir(
                                    cfg_for_combo,
                                    root_key="POST_DOCKED_DIR",
                                    default_name="post_docked",
                                    run_id=str(context.run_id),
                                    pdb_id=str(pdb_id),
                                    variant_label=str(context.label),
                                    ph_tag=ph_tag,
                                )
                                post_csv = deps.resolve_combo_post_consensus_csv(
                                    combo_post_dir,
                                    str(library_name or ""),
                                )
                                post_keys = deps.load_post_scored_ligand_keys(post_csv)
                                deps.update_combo_coverage_snapshot(
                                    cfg_for_combo,
                                    run_id=str(context.run_id),
                                    pdb_id=str(pdb_id),
                                    variant_label=str(context.label),
                                    ph_tag=ph_tag,
                                    library_name=str(library_name or ""),
                                    expected_keys=chunk_keys,
                                    post_keys=(post_keys & chunk_keys if chunk_keys else post_keys),
                                )
                    except Exception:
                        logging.warning(
                            "[scorch-rescore.hook] action=skip reason=unexpected_exception pdb_id=%s variant=%s ph=%s",
                            pdb_id,
                            context.label,
                            ph_token,
                            exc_info=True,
                        )

            pdb_elapsed = time.time() - pdb_start
            if is_chunk_limited_work:
                logging.debug(
                    "[distributed.chunk.manifest] action=skip_whole_target_success chunk_id=%s pdb=%s variant=%s ph=%s",
                    str(chunk_id or ""),
                    pdb_id,
                    context.label,
                    ph_token,
                )
            else:
                try:
                    deps.update_manifest_for_protein_success(
                        cfg_for_combo,
                        context.run_id,
                        pdb_id,
                        context.label,
                        pdb_elapsed,
                        ph_tag=ph_tag,
                    )
                except Exception:
                    logging.warning(
                        "Failed to update run_manifest for success of %s (%s,%s)",
                        pdb_id,
                        context.label,
                        ph_token,
                        exc_info=True,
                    )

            if shared_state.scorch_queue_service is not None:
                with shared_state.pending_retention_lock:
                    shared_state.pending_retention_combos.add(
                        (
                            str(pdb_id).upper(),
                            str(context.label or "").strip().upper() or "BASE",
                            deps.normalize_ph_tag_token(ph_tag),
                        )
                    )
                    if chunk_ligand_bases is not None:
                        combo_key = (
                            str(pdb_id).upper(),
                            str(context.label or "").strip().upper() or "BASE",
                            deps.normalize_ph_tag_token(ph_tag),
                        )
                        payload = shared_state.pending_coverage_refresh.setdefault(
                            combo_key,
                            {
                                "expected": set(),
                                "library": str(library_name or ""),
                            },
                        )
                        expected_payload = payload.get("expected")
                        if isinstance(expected_payload, set):
                            expected_payload.update(chunk_keys)
                        elif isinstance(expected_payload, list):
                            payload["expected"] = set(expected_payload) | set(chunk_keys)
                        else:
                            payload["expected"] = set(chunk_keys)
                        if not str(payload.get("library") or "").strip():
                            payload["library"] = str(library_name or "")
            else:
                try:
                    with deps.acquire_global_slot(
                        cfg_for_combo,
                        cores=1,
                        min_cores=1,
                        priority=0,
                        task_id=f"main:housekeeping:retention:{pdb_id}:{context.label}:{ph_token}",
                        task_type="housekeeping",
                    ):
                        deps.maybe_run_artifact_retention_for_combo(
                            cfg_for_combo,
                            context.run_id,
                            pdb_id,
                            str(context.label or "").strip().upper() or "BASE",
                            None if ph_token == "base" else str(ph_tag),
                            lock=shared_state.retention_lock,
                        )
                except Exception:
                    logging.warning(
                        "[artifact-retention.hook] action=skip reason=unexpected_exception pdb_id=%s variant=%s ph=%s",
                        pdb_id,
                        context.label,
                        ph_token,
                        exc_info=True,
                    )
            scope_metric_key = (
                str(pdb_id).upper(),
                str(context.label or "").strip().upper() or "BASE",
                deps.normalize_ph_tag_token(ph_tag),
            )
            if scope_metric_key not in shared_state.run_scope_completion_times:
                shared_state.run_scope_completion_times[scope_metric_key] = float(
                    max(0.0, time.time() - context.global_start)
                )
            return True

        except Exception as exc:
            if chunk_ligand_bases is not None:
                chunk_reason = str(
                    cfg_for_pdb.get(_CHUNK_FAILURE_REASON_KEY) or str(exc) or ""
                ).strip()
                retryable = bool(
                    cfg_for_pdb.get(_CHUNK_FAILURE_RETRYABLE_KEY)
                    or _chunk_failure_retryable(chunk_reason)
                )
                _set_chunk_failure_metadata(
                    cfg_for_pdb,
                    reason=chunk_reason,
                    retryable=retryable,
                )
                logging.warning(
                    "[distributed.chunk.failure] chunk_id=%s pdb=%s variant=%s ph=%s retryable=%s reason=%s",
                    str(chunk_id or ""),
                    pdb_id,
                    context.label,
                    ph_token,
                    str(bool(retryable)).lower(),
                    chunk_reason or str(exc),
                )
                return False

            # Per-PDB failure handling
            exc_type = type(exc).__name__
            exc_msg = str(exc)
            traceback_str = traceback.format_exc()

            fail_log_path = (
                Path(shared_state.failed_root)
                / f"{context.run_id}__{pdb_id}__{context.label}__{ph_token}.log"
            )
            tmp_fail_log_path = fail_log_path.with_suffix(fail_log_path.suffix + ".part")

            # Write a dedicated failure log for this PDB+variant
            with tmp_fail_log_path.open("w", encoding="utf-8") as fh:
                fh.write(
                    f"[FAILED PDB]\n"
                    f"  pdb_id        = {pdb_id}\n"
                    f"  variant_label = {context.label}\n"
                    f"  ph            = {ph_token}\n"
                    f"  mode          = {context.mode}\n"
                    f"  pdb_file      = {pdb_file} run_id={context.run_id}\n"
                    f"  exception     = {exc_type}: {exc_msg}\n\n"
                    f"[TRACEBACK]\n"
                    f"{traceback_str}\n"
                )
                traceback.print_exc(file=fh)
            os.replace(tmp_fail_log_path, fail_log_path)

            logging.error(
                "[apo-holo] FAILED pdb_id=%s variant_label=%s; see failure log at %s",
                pdb_id,
                context.label,
                fail_log_path,
            )

            shared_state.failed_entries.append(
                (pdb_id, context.label, str(fail_log_path), exc_type, exc_msg)
            )

            try:
                deps.update_manifest_for_protein_failure(
                    cfg_for_pdb,
                    context.run_id,
                    pdb_id,
                    context.label,
                    fail_log_path,
                    ph_tag=ph_tag,
                )
            except Exception:
                logging.warning(
                    "Failed to update run_manifest for failure of %s (%s,%s)",
                    pdb_id,
                    context.label,
                    ph_token,
                    exc_info=True,
                )
            return False

    return _process_one
