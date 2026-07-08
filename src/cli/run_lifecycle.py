from __future__ import annotations

import datetime
import logging
import os
import threading
import time
from typing import Any, Mapping, Optional, Sequence

from cli.cli_utils import _norm_pdb_id
from cli.distributed_context import (
    collect_failed_entries_from_markers,
    wait_for_all_markers,
    wait_for_all_phase_markers,
    write_completion_marker,
    write_json_atomic,
    write_phase_marker,
)
from cli.distributed_chunk_planner import (
    _chunk_ligand_key,
    _load_post_scored_ligand_keys,
    _reconcile_missing_post_scorch_outputs,
    _resolve_combo_output_dir,
    _resolve_combo_post_consensus_csv,
    _update_combo_coverage_snapshot,
)
from docking.global_scheduler import acquire_global_slot
from post_docking.mmgbsa.mmgbsa_pipeline import _maybe_run_mmgbsa_for_pdb
from cli.postrun_hooks_runtime import (
    _log_rescore_verification,
    _maybe_run_artifact_retention_for_combos,
    _maybe_run_dud_eval,
    _maybe_run_master_schema_export,
    _maybe_run_report_generation,
    _maybe_run_throughput_integrity,
    _write_run_efficiency_report,
)
from cli.run_manifest_runtime import finalize_run_manifest, load_run_manifest
from config.output_paths import runtime_root


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return float(default)
    try:
        value = float(raw)
    except Exception:
        return float(default)
    if value <= 0.0:
        return float(default)
    return float(value)


def _env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return bool(default)
    token = str(raw).strip().lower()
    if token in {"1", "true", "yes", "on"}:
        return True
    if token in {"0", "false", "no", "off"}:
        return False
    return bool(default)


def finalize_scheduler_and_retention(
    *,
    cfg: Mapping[str, Any],
    run_id: str,
    global_start: float,
    run_scope_completion_times: Mapping[tuple[str, str, str], float],
    run_chunk_rebalance_count: int,
    dist_ctx: Any,
    scorch_queue_service: Any,
    pending_coverage_refresh: Mapping[tuple[str, str, str], Mapping[str, Any]],
    pending_retention_combos: set[tuple[str, str, str]],
    retention_lock: threading.Lock,
) -> tuple[int, bool]:
    run_elapsed = float(max(0.0, time.time() - global_start))
    completion_samples = sorted(
        float(v) for v in run_scope_completion_times.values() if isinstance(v, (int, float))
    )
    first_complete = completion_samples[0] if completion_samples else None
    q25 = float(run_elapsed) * 0.25
    q50 = float(run_elapsed) * 0.50
    q75 = float(run_elapsed) * 0.75
    completion_at_25 = sum(1 for v in completion_samples if float(v) <= q25)
    completion_at_50 = sum(1 for v in completion_samples if float(v) <= q50)
    completion_at_75 = sum(1 for v in completion_samples if float(v) <= q75)
    scheduler_summary = {
        "run_id": str(run_id),
        "time_to_first_complete_protein_s": (
            None if first_complete is None else float(first_complete)
        ),
        "completed_proteins_at_25pct_wallclock": int(completion_at_25),
        "completed_proteins_at_50pct_wallclock": int(completion_at_50),
        "completed_proteins_at_75pct_wallclock": int(completion_at_75),
        "final_makespan_s": float(run_elapsed),
        "chunk_rebalance_count": int(max(0, run_chunk_rebalance_count)),
    }
    logging.info(
        "[execution.summary] run_id=%s first_complete_s=%s complete_25=%d complete_50=%d complete_75=%d makespan_s=%.2f chunk_rebalance_count=%d",
        str(run_id),
        (
            f"{float(first_complete):.2f}"
            if isinstance(first_complete, (int, float))
            else "none"
        ),
        int(completion_at_25),
        int(completion_at_50),
        int(completion_at_75),
        float(run_elapsed),
        int(max(0, run_chunk_rebalance_count)),
    )
    try:
        summary_path = (
            runtime_root(cfg, "MANIFESTS_DIR", "manifests")
            / str(run_id)
            / "scheduler_summary.json"
        )
        write_json_atomic(summary_path, scheduler_summary)
    except Exception:
        logging.warning(
            "[execution.summary] action=write_failed run_id=%s",
            run_id,
            exc_info=True,
        )

    scorch_queue_failed = 0
    retention_barrier_ok = True
    if scorch_queue_service is not None:
        try:
            scorch_queue_summary = scorch_queue_service.close_and_wait()
            scorch_queue_failed = int(scorch_queue_summary.get("failed", 0))
        except Exception:
            logging.warning(
                "[scorch-rescore.queue] action=close_failed run_id=%s",
                run_id,
                exc_info=True,
            )
        reconcile_allowed = True
        if dist_ctx.enabled:
            try:
                queue_barrier_timeout = _env_float(
                    "ATLAS_SCORCH_QUEUE_DRAIN_TIMEOUT_SEC",
                    min(
                        180.0,
                        float(
                            max(
                                30.0,
                                getattr(dist_ctx, "barrier_timeout_sec", 1800.0),
                            )
                        ),
                    ),
                )
                queue_barrier_poll = _env_float(
                    "ATLAS_SCORCH_QUEUE_DRAIN_POLL_SEC",
                    max(0.5, float(getattr(dist_ctx, "barrier_poll_sec", 5.0))),
                )
                write_phase_marker(
                    dist_ctx,
                    cfg,
                    phase="scorch_queue_closed",
                    payload={
                        "run_id": str(run_id),
                        "scorch_queue_failed": int(scorch_queue_failed),
                    },
                )
                wait_start = time.time()
                wait_deadline = wait_start + max(5.0, float(queue_barrier_timeout))
                queue_barrier_ok = False
                missing_task_ids: list[int] = []
                while time.time() < wait_deadline:
                    remaining = max(0.1, wait_deadline - time.time())
                    queue_barrier_ok, missing_task_ids = wait_for_all_phase_markers(
                        dist_ctx,
                        cfg,
                        phase="scorch_queue_closed",
                        timeout_sec=min(15.0, remaining),
                        poll_sec=max(0.5, float(queue_barrier_poll)),
                    )
                    if queue_barrier_ok:
                        break
                    logging.info(
                        "[scorch-rescore.queue-barrier] action=wait run_id=%s task_id=%d elapsed_s=%.1f missing_task_ids=%s",
                        run_id,
                        int(dist_ctx.task_id),
                        max(0.0, time.time() - wait_start),
                        ",".join(str(x) for x in missing_task_ids) or "none",
                    )
                if queue_barrier_ok:
                    logging.info(
                        "[scorch-rescore.queue-barrier] action=ok run_id=%s task_id=%d elapsed_s=%.1f",
                        run_id,
                        int(dist_ctx.task_id),
                        max(0.0, time.time() - wait_start),
                    )
                else:
                    reconcile_allowed = False
                    logging.warning(
                        "[scorch-rescore.queue-barrier] action=timeout_skip_reconcile run_id=%s task_id=%d elapsed_s=%.1f missing_task_ids=%s",
                        run_id,
                        int(dist_ctx.task_id),
                        max(0.0, time.time() - wait_start),
                        ",".join(str(x) for x in missing_task_ids) or "none",
                    )
            except Exception:
                reconcile_allowed = False
                logging.warning(
                    "[scorch-rescore.queue-barrier] action=error_skip_reconcile run_id=%s task_id=%d",
                    run_id,
                    int(getattr(dist_ctx, "task_id", 0)),
                    exc_info=True,
                )
        try:
            if reconcile_allowed:
                reconcile_summary = _reconcile_missing_post_scorch_outputs(
                    cfg,
                    run_id=str(run_id),
                    dist_ctx=dist_ctx,
                    coverage_refresh=pending_coverage_refresh,
                )
            else:
                reconcile_summary = {
                    "attempted": 0,
                    "succeeded": 0,
                    "failed": 0,
                    "skipped_owner": 0,
                }
        except Exception:
            reconcile_summary = {
                "attempted": 0,
                "succeeded": 0,
                "failed": 0,
                "skipped_owner": 0,
            }
            logging.warning(
                "[scorch-rescore.reconcile] action=error run_id=%s",
                run_id,
                exc_info=True,
            )
        scorch_queue_failed += int(reconcile_summary.get("failed", 0))
        if int(reconcile_summary.get("attempted", 0)) > 0:
            logging.info(
                "[scorch-rescore.reconcile] action=summary run_id=%s attempted=%d succeeded=%d failed=%d skipped_owner=%d",
                run_id,
                int(reconcile_summary.get("attempted", 0)),
                int(reconcile_summary.get("succeeded", 0)),
                int(reconcile_summary.get("failed", 0)),
                int(reconcile_summary.get("skipped_owner", 0)),
            )
        if dist_ctx.enabled:
            try:
                scorch_barrier_timeout = _env_float(
                    "ATLAS_SCORCH_BARRIER_TIMEOUT_SEC",
                    min(180.0, float(max(30.0, getattr(dist_ctx, "barrier_timeout_sec", 1800.0)))),
                )
                scorch_barrier_poll = _env_float(
                    "ATLAS_SCORCH_BARRIER_POLL_SEC",
                    max(0.5, float(getattr(dist_ctx, "barrier_poll_sec", 5.0))),
                )
                scorch_barrier_strict = _env_bool("ATLAS_SCORCH_BARRIER_STRICT", False)
                write_phase_marker(
                    dist_ctx,
                    cfg,
                    phase="scorch_drain",
                    payload={
                        "run_id": str(run_id),
                        "scorch_queue_failed": int(scorch_queue_failed),
                        "scorch_reconcile_attempted": int(
                            reconcile_summary.get("attempted", 0)
                        ),
                        "scorch_reconcile_failed": int(
                            reconcile_summary.get("failed", 0)
                        ),
                    },
                )
                wait_start = time.time()
                wait_deadline = wait_start + max(5.0, float(scorch_barrier_timeout))
                barrier_ok = False
                missing_task_ids: list[int] = []
                while time.time() < wait_deadline:
                    remaining = max(0.1, wait_deadline - time.time())
                    barrier_ok, missing_task_ids = wait_for_all_phase_markers(
                        dist_ctx,
                        cfg,
                        phase="scorch_drain",
                        timeout_sec=min(15.0, remaining),
                        poll_sec=max(0.5, float(scorch_barrier_poll)),
                    )
                    if barrier_ok:
                        break
                    logging.info(
                        "[scorch-rescore.barrier] action=wait run_id=%s task_id=%d elapsed_s=%.1f missing_task_ids=%s",
                        run_id,
                        int(dist_ctx.task_id),
                        max(0.0, time.time() - wait_start),
                        ",".join(str(x) for x in missing_task_ids) or "none",
                    )
                if not barrier_ok:
                    if scorch_barrier_strict:
                        retention_barrier_ok = False
                        logging.warning(
                            "[scorch-rescore.barrier] action=timeout_strict run_id=%s task_id=%d elapsed_s=%.1f missing_task_ids=%s",
                            run_id,
                            int(dist_ctx.task_id),
                            max(0.0, time.time() - wait_start),
                            ",".join(str(x) for x in missing_task_ids) or "none",
                        )
                    else:
                        logging.warning(
                            "[scorch-rescore.barrier] action=timeout_fail_open run_id=%s task_id=%d elapsed_s=%.1f missing_task_ids=%s",
                            run_id,
                            int(dist_ctx.task_id),
                            max(0.0, time.time() - wait_start),
                            ",".join(str(x) for x in missing_task_ids) or "none",
                        )
                else:
                    logging.info(
                        "[scorch-rescore.barrier] action=ok run_id=%s task_id=%d elapsed_s=%.1f",
                        run_id,
                        int(dist_ctx.task_id),
                        max(0.0, time.time() - wait_start),
                    )
            except Exception:
                retention_barrier_ok = False
                logging.warning(
                    "[scorch-rescore.barrier] action=error run_id=%s task_id=%d",
                    run_id,
                    int(getattr(dist_ctx, "task_id", 0)),
                    exc_info=True,
                )
        for coverage_key, payload in sorted(pending_coverage_refresh.items()):
            try:
                pdb_cov, variant_cov, ph_cov = coverage_key
                ph_tag_cov = None if str(ph_cov).strip().lower() == "base" else str(ph_cov)
                library_cov = str(payload.get("library") or "")
                expected_cov = payload.get("expected")
                if isinstance(expected_cov, set):
                    expected_keys_cov = set(expected_cov)
                elif isinstance(expected_cov, list):
                    expected_keys_cov = {
                        _chunk_ligand_key(str(x))
                        for x in expected_cov
                        if _chunk_ligand_key(str(x))
                    }
                else:
                    expected_keys_cov = set()
                variant_token = str(variant_cov or "").strip().upper()
                combo_post_dir = _resolve_combo_output_dir(
                    cfg,
                    root_key="POST_DOCKED_DIR",
                    default_name="post_docked",
                    run_id=str(run_id),
                    pdb_id=str(pdb_cov).upper(),
                    variant_label=str(variant_token or "BASE"),
                    ph_tag=ph_tag_cov,
                )
                post_csv = _resolve_combo_post_consensus_csv(combo_post_dir, library_cov)
                post_keys = _load_post_scored_ligand_keys(post_csv)
                _update_combo_coverage_snapshot(
                    cfg,
                    run_id=str(run_id),
                    pdb_id=str(pdb_cov).upper(),
                    variant_label=str(variant_token or "BASE"),
                    ph_tag=ph_tag_cov,
                    library_name=library_cov,
                    expected_keys=expected_keys_cov,
                    post_keys=(
                        post_keys & expected_keys_cov
                        if expected_keys_cov
                        else post_keys
                    ),
                )
            except Exception:
                logging.warning(
                    "[coverage-snapshot.refresh] run_id=%s combo=%s action=skip",
                    run_id,
                    coverage_key,
                    exc_info=True,
                )

    return int(scorch_queue_failed), bool(retention_barrier_ok)


def finalize_deferred_retention(
    *,
    cfg: Mapping[str, Any],
    run_id: str,
    retention_barrier_ok: bool,
    pending_retention_combos: set[tuple[str, str, str]],
    retention_lock: threading.Lock,
) -> None:
    retention_started = time.time()
    pending_combo_list: list[tuple[str, str, str]] = [
        (
            str(pending_pdb_id).upper(),
            str(pending_variant or "").strip().upper() or "LEGACY",
            "none" if str(pending_ph).strip().lower() == "base" else str(pending_ph),
        )
        for pending_pdb_id, pending_variant, pending_ph in sorted(pending_retention_combos)
    ]
    if pending_combo_list and not retention_barrier_ok:
        logging.warning(
            "[pipeline.phase] skip=retention_deferred run_id=%s reason=scorch_barrier_incomplete combos=%d",
            run_id,
            len(pending_combo_list),
        )
        return
    if not pending_combo_list:
        return

    logging.info(
        "[pipeline.phase] start=retention_deferred run_id=%s combos=%d",
        run_id,
        len(pending_combo_list),
    )
    try:
        with acquire_global_slot(
            cfg,
            cores=1,
            min_cores=1,
            priority=0,
            task_id=f"main:housekeeping:retention-batch:{run_id}",
            task_type="housekeeping",
        ):
            _maybe_run_artifact_retention_for_combos(
                cfg,
                run_id,
                pending_combo_list,
                lock=retention_lock,
            )
    except Exception:
        logging.warning(
            "[artifact-retention.hook] action=skip reason=unexpected_exception run_id=%s combos=%d",
            run_id,
            len(pending_combo_list),
            exc_info=True,
        )
    finally:
        logging.info(
            "[pipeline.phase] done=retention_deferred run_id=%s combos=%d elapsed_sec=%.2f",
            run_id,
            len(pending_combo_list),
            time.time() - retention_started,
        )


def finalize_post_run(
    *,
    cfg: Mapping[str, Any],
    run_id: str,
    global_start: float,
    failed_entries: list[tuple[str, str, str, str, str]],
    plan_only: bool,
    variants: Sequence[Optional[str]],
    pdb_files: Sequence[str],
    test_mode: str,
    dist_ctx: Any,
    dist_combo_chunk_mode: bool,
    distributed_assigned_pdb_ids: set[str],
    is_bench: bool,
    is_bench2: bool,
    is_bench_small: bool,
) -> bool:
    if not dist_ctx.enabled:
        try:
            finalize_run_manifest(
                cfg, run_id, start_time=global_start, failed_entries=failed_entries
            )
        except Exception:
            logging.warning("Failed to finalize run_manifest.yaml", exc_info=True)
    else:
        logging.info(
            "[run-manifest] finalize deferred to dependency finalizer job (task_id=%d)",
            dist_ctx.task_id,
        )

    if plan_only:
        return True

    try:
        mmgbsa_auto_run = str(cfg.get("MMGBSA_AUTO_RUN", "false")).strip().lower() in {
            "1",
            "true",
            "yes",
            "on",
        }
        if not mmgbsa_auto_run:
            logging.info(
                "[mmgbsa.pipeline] action=skip reason=auto_run_disabled hint=use_atlas_mmgbsa_for_report_hits"
            )
        failed_lookup = {(entry[0], entry[1]) for entry in failed_entries}
        if mmgbsa_auto_run:
            for variant in variants:
                label = "legacy" if variant is None else str(variant).lower()
                mmgbsa_variant_token: str | None = (
                    None if variant is None else str(variant).upper()
                )
                legacy_mode = variant is None
                for pdb_file in pdb_files:
                    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0].upper()
                    if (pdb_id, label) in failed_lookup:
                        continue
                    _maybe_run_mmgbsa_for_pdb(
                        cfg,
                        pdb_file,
                        pdb_id,
                        mmgbsa_variant_token,
                        run_id,
                        test_mode,
                        legacy_mode,
                    )
    except Exception:
        strict_mmgbsa = False
        try:
            strict_mmgbsa = bool(cfg.get("MMGBSA_STRICT", False))
        except Exception:
            strict_mmgbsa = False
        if strict_mmgbsa:
            raise
        logging.warning(
            "[mmgbsa.pipeline] action=skip reason=unexpected_exception", exc_info=True
        )

    if dist_ctx.enabled:
        assigned_ids: list[str] = []
        if dist_combo_chunk_mode and distributed_assigned_pdb_ids:
            assigned_ids = sorted({str(x).upper() for x in distributed_assigned_pdb_ids})
        else:
            for pdb_file in pdb_files:
                nid = _norm_pdb_id(pdb_file)
                if nid:
                    assigned_ids.append(nid)
        write_completion_marker(
            dist_ctx,
            cfg,
            assigned_pdb_ids=sorted(set(assigned_ids)),
            failed_entries=failed_entries,
        )
        if dist_ctx.is_leader:
            try:
                barrier_ok, missing_task_ids = wait_for_all_markers(dist_ctx, cfg)
                if not barrier_ok:
                    logging.warning(
                        "[distributed.finalize.fallback] barrier_timeout=true missing_task_ids=%s",
                        missing_task_ids,
                    )
                merged_failures = collect_failed_entries_from_markers(dist_ctx, cfg)
                if missing_task_ids:
                    merged_failures.append(
                        (
                            "DISTRIBUTED",
                            "distributed",
                            "-",
                            "MissingMarkers",
                            "missing_task_ids=" + ",".join(str(x) for x in missing_task_ids),
                        )
                    )
                manifest_started_epoch = float(global_start)
                try:
                    manifest_data = load_run_manifest(cfg, run_id) or {}
                    timing = (
                        manifest_data.get("timing", {})
                        if isinstance(manifest_data, dict)
                        else {}
                    )
                    started_iso = timing.get("started_at")
                    if started_iso:
                        started_dt = datetime.datetime.fromisoformat(
                            str(started_iso).replace("Z", "+00:00")
                        )
                        manifest_started_epoch = float(started_dt.timestamp())
                except Exception:
                    logging.warning(
                        "[distributed.finalize.fallback] run_id=%s action=manifest_started_at_parse_failed",
                        run_id,
                        exc_info=True,
                    )
                finalize_run_manifest(
                    cfg,
                    run_id,
                    start_time=manifest_started_epoch,
                    failed_entries=merged_failures,
                )
                _write_run_efficiency_report(cfg, run_id)
                try:
                    _maybe_run_dud_eval(cfg, run_id, list(pdb_files))
                except Exception:
                    logging.warning(
                        "[dud-eval.invoke] action=skip reason=unexpected_exception",
                        exc_info=True,
                    )
                try:
                    _maybe_run_master_schema_export(cfg, run_id)
                except Exception:
                    logging.warning(
                        "[master-export.invoke] action=skip reason=unexpected_exception",
                        exc_info=True,
                    )
                try:
                    _maybe_run_report_generation(cfg, run_id)
                except Exception:
                    logging.warning(
                        "[report.invoke] action=skip reason=unexpected_exception",
                        exc_info=True,
                    )

                throughput_integrity_strict = bool(is_bench or is_bench2 or is_bench_small)
                throughput_integrity_ok = True
                try:
                    throughput_integrity_ok = _maybe_run_throughput_integrity(
                        cfg, run_id, strict=throughput_integrity_strict
                    )
                except Exception:
                    throughput_integrity_ok = False
                    logging.warning(
                        "[throughput-integrity.invoke] action=skip reason=unexpected_exception",
                        exc_info=True,
                    )
                if throughput_integrity_strict and not throughput_integrity_ok:
                    logging.error(
                        "[throughput-integrity.enforce] run_id=%s strict=true status=failed",
                        run_id,
                    )
                    raise SystemExit(2)

                try:
                    _log_rescore_verification(run_id, cfg)
                except Exception:
                    logging.warning(
                        "[post-check.invoke] action=skip reason=unexpected_exception",
                        exc_info=True,
                    )
                logging.info(
                    "[distributed.finalize.fallback] run_id=%s status=done",
                    run_id,
                )
            except Exception:
                logging.warning(
                    "[distributed.finalize.fallback] run_id=%s action=failed",
                    run_id,
                    exc_info=True,
                )
        logging.info(
            "[distributed.worker.done] task_id=%d role=worker action=exit_after_marker",
            dist_ctx.task_id,
        )
        return True

    try:
        _maybe_run_dud_eval(cfg, run_id, list(pdb_files))
    except Exception:
        logging.warning(
            "[dud-eval.invoke] action=skip reason=unexpected_exception", exc_info=True
        )

    try:
        _maybe_run_master_schema_export(cfg, run_id)
    except Exception:
        logging.warning(
            "[master-export.invoke] action=skip reason=unexpected_exception",
            exc_info=True,
        )

    try:
        _maybe_run_report_generation(cfg, run_id)
    except Exception:
        logging.warning(
            "[report.invoke] action=skip reason=unexpected_exception", exc_info=True
        )
    try:
        _write_run_efficiency_report(cfg, run_id)
    except Exception:
        logging.warning(
            "[run-efficiency.invoke] action=skip reason=unexpected_exception",
            exc_info=True,
        )

    throughput_integrity_strict = bool(is_bench or is_bench2 or is_bench_small)
    throughput_integrity_ok = True
    try:
        throughput_integrity_ok = _maybe_run_throughput_integrity(
            cfg,
            run_id,
            strict=throughput_integrity_strict,
        )
    except Exception:
        throughput_integrity_ok = False
        logging.warning(
            "[throughput-integrity.invoke] action=skip reason=unexpected_exception",
            exc_info=True,
        )
    if throughput_integrity_strict and not throughput_integrity_ok:
        logging.error(
            "[throughput-integrity.enforce] run_id=%s strict=true status=failed",
            run_id,
        )
        raise SystemExit(2)

    try:
        _log_rescore_verification(run_id, cfg)
    except Exception:
        logging.warning(
            "[post-check.invoke] action=skip reason=unexpected_exception", exc_info=True
        )
    return False
