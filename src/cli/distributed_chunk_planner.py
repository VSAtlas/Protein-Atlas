from __future__ import annotations

import json
import logging
import os
import time
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from cli.distributed_context import write_json_atomic
from cli.planner_chunking import (
    _distributed_chunk_cpu_per_worker,
    _distributed_chunk_local_workers,
    _log_distributed_cpu_node,
    _resolve_runtime_cpu_settings,
    _split_ligands_into_chunks,
)
from cli.planner_combo_plan import (
    _distributed_combo_owner_task_id,
    _load_or_create_variant_chunk_plan,
)
from cli.planner_combo_work import (
    _build_combo_work_items,
    _combo_pdb_id_from_file,
    _combo_prep_key,
    _combo_receptor_ready,
    _local_combo_scope_key,
    _normalize_ph_tag_token,
    _pick_next_combo_index,
    _refresh_local_prep_ready_scopes,
    _scope_prep_key,
)
from cli.planner_manifest import (
    _chunk_ligand_key,
    _chunk_planner_manifest_mode,
    _chunk_planner_manifest_sample_size,
    _chunk_planner_scan_workers,
    _collect_combo_library_roots,
    _evaluate_manifest_preflight_roots,
    _manifest_prebuild_enabled,
    _rebuild_library_manifests,
    _resolve_combo_ligand_bases,
    _resolve_combo_ligand_streams,
    _validate_chunk_stream_resolution,
)
from cli.planner_validation import (
    _load_post_scored_ligand_keys,
    _load_scored_ligand_keys_from_summary,
    _resolve_combo_docking_summary_csv,
    _resolve_combo_output_dir,
    _resolve_combo_post_consensus_csv,
    _update_combo_coverage_snapshot,
    _verify_chunk_combo_outputs,
)
from cli.postrun_hooks_runtime import (
    _maybe_run_scorch_rescore_for_pdb,
    _resolve_scorch_execution_mode,
)
from config.output_paths import runtime_root
from docking.global_scheduler import global_scheduler_capacity
from post_docking.rescoring.scorch_coverage import (
    evaluate_scorch_coverage,
    remove_done_sentinel,
    resolve_top_fraction,
    write_coverage_summary,
)

__all__ = [
    "_build_combo_work_items",
    "_chunk_ligand_key",
    "_chunk_planner_manifest_mode",
    "_chunk_planner_manifest_sample_size",
    "_chunk_planner_scan_workers",
    "_combo_pdb_id_from_file",
    "_combo_prep_key",
    "_combo_receptor_ready",
    "_distributed_combo_owner_task_id",
    "_load_or_create_variant_chunk_plan",
    "_local_combo_scope_key",
    "_manifest_prebuild_enabled",
    "_normalize_ph_tag_token",
    "_pick_next_combo_index",
    "_refresh_local_prep_ready_scopes",
    "_resolve_combo_ligand_bases",
    "_resolve_combo_ligand_streams",
    "_validate_chunk_stream_resolution",
    "_resolve_combo_docking_summary_csv",
    "_resolve_combo_output_dir",
    "_resolve_combo_post_consensus_csv",
    "_scope_prep_key",
    "_load_scored_ligand_keys_from_summary",
    "_load_post_scored_ligand_keys",
    "_update_combo_coverage_snapshot",
    "_verify_chunk_combo_outputs",
    "_split_ligands_into_chunks",
    "_distributed_chunk_local_workers",
    "_distributed_chunk_cpu_per_worker",
    "_resolve_runtime_cpu_settings",
    "_log_distributed_cpu_node",
]

_CHUNK_TARGET_SIZE = 96
_CHUNK_MIN_SIZE = 8
_CHUNK_MAX_SIZE = 256
_CHUNK_OVERSHARD_FACTOR = 6
_CHUNK_LOCAL_WORKER_CPU_DIVISOR = 1
_CHUNK_LOCAL_WORKER_MAX = 512
_CHUNK_CLAIM_LEASE_SEC = 900.0
_CHUNK_CLAIM_HEARTBEAT_SEC = 60.0
_CHUNK_IDLE_WAIT_SEC = 0.5
_CHUNK_IDLE_TIMEOUT_SEC = 900.0
_CHUNK_STATE_RECONCILE_SEC = 10.0
_CHUNK_MAX_ATTEMPTS = 3
_COMBO_PREP_LEASE_SEC = 14400.0
_COMBO_PREP_WAIT_POLL_SEC = 2.0
_COMBO_PREP_WAIT_TIMEOUT_SEC = 14400.0
_MANIFEST_PREFLIGHT_WAIT_LOG_SEC = 15.0


def _resolve_global_scheduler_plan(
    cfg: Mapping[str, Any], *, pdb_count: int, cpu: int
) -> dict[str, Any]:
    scheduler_enabled = bool(cfg.get("ENABLE_GLOBAL_SCHEDULER", True))
    scheduler_policy = str(cfg.get("GLOBAL_ADMISSION_POLICY", "adaptive")).strip().lower()
    scheduler_cpus = global_scheduler_capacity(dict(cfg), cpu)
    min_parallel = max(1, int(cfg.get("GLOBAL_MIN_PROTEINS", 1) or 1))
    max_parallel_cfg = int(cfg.get("GLOBAL_MAX_PROTEINS", 0) or 0)

    if max_parallel_cfg > 0:
        max_parallel = max(
            min_parallel, min(max_parallel_cfg, int(pdb_count), scheduler_cpus)
        )
    else:
        max_parallel = max(min_parallel, min(int(pdb_count), scheduler_cpus))

    if scheduler_policy == "adaptive":
        initial_target = max(
            min_parallel, min(max_parallel, max(1, scheduler_cpus // 2))
        )
    else:
        initial_target = max_parallel

    return {
        "enabled": scheduler_enabled,
        "policy": scheduler_policy,
        "scheduler_cpus": scheduler_cpus,
        "min_parallel": min_parallel,
        "max_parallel": max_parallel,
        "initial_target": max(1, min(initial_target, max_parallel)),
    }


def _run_scoped_root_from_cfg(
    cfg: Mapping[str, Any],
    *,
    root_key: str,
    default_name: str,
    run_id: str,
) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    root = Path(str(cfg.get(root_key, repo_root / default_name) or repo_root / default_name)).resolve()
    token = str(run_id or "").strip()
    if token and root.name != token:
        root = root / token
    return root


def _scorch_coverage_complete_for_reconcile(
    cfg: Mapping[str, Any],
    *,
    run_id: str,
    pdb_id: str,
    variant_label: str,
    ph_tag: Optional[str],
    logger: logging.Logger,
) -> bool:
    docked_run_root = _run_scoped_root_from_cfg(
        cfg,
        root_key="DOCKED_DIR",
        default_name="docked",
        run_id=str(run_id),
    )
    post_run_root = _run_scoped_root_from_cfg(
        cfg,
        root_key="POST_DOCKED_DIR",
        default_name="post_docked",
        run_id=str(run_id),
    )
    summary = evaluate_scorch_coverage(
        docked_run_root=docked_run_root,
        post_run_root=post_run_root,
        pdb_id=pdb_id,
        variant=variant_label,
        ph=ph_tag or "",
        top_fraction=resolve_top_fraction(cfg),
        decoy_prefix=str(cfg.get("DECOY_PREFIX") or cfg.get("DUD_PREFIX") or "dud"),
    )
    try:
        write_coverage_summary(post_run_root, summary)
    except Exception:
        logger.debug(
            "[scorch-rescore.reconcile] action=coverage_summary status=failed run_id=%s pdb_id=%s",
            str(run_id),
            str(pdb_id),
            exc_info=True,
        )
    if summary.all_complete:
        return True
    try:
        remove_done_sentinel(post_run_root, summary)
    except Exception:
        logger.debug(
            "[scorch-rescore.reconcile] action=clear_done status=failed run_id=%s pdb_id=%s",
            str(run_id),
            str(pdb_id),
            exc_info=True,
        )
    logger.warning(
        "[scorch-rescore.reconcile] action=retry reason=scorch_coverage_incomplete run_id=%s pdb_id=%s variant=%s ph=%s details=%s",
        str(run_id),
        str(pdb_id),
        str(variant_label or "LEGACY"),
        str(ph_tag or "base"),
        ",".join(summary.reasons) or "unknown",
    )
    return False


def _reconcile_missing_post_scorch_outputs(
    cfg: Mapping[str, Any],
    *,
    run_id: str,
    dist_ctx: Any,
    coverage_refresh: Mapping[tuple[str, str, str], Mapping[str, Any]],
) -> dict[str, int]:
    logger = logging.getLogger("scorch-rescore-reconcile")
    if not coverage_refresh:
        return {"attempted": 0, "succeeded": 0, "failed": 0, "skipped_owner": 0}

    execution_mode = _resolve_scorch_execution_mode(cfg)
    attempted = 0
    succeeded = 0
    failed = 0
    skipped_owner = 0

    for coverage_key, payload in sorted(coverage_refresh.items()):
        pdb_cov, variant_cov, ph_cov = coverage_key
        pdb_id = str(pdb_cov).upper()
        variant_token = str(variant_cov or "").strip().upper() or "BASE"
        ph_norm = str(ph_cov).strip()
        ph_tag = None if ph_norm.lower() == "base" else ph_norm

        if bool(getattr(dist_ctx, "enabled", False)):
            owner_task = _distributed_combo_owner_task_id(
                dist_ctx,
                variant_label=str(variant_token),
                pdb_id=str(pdb_id),
                ph_tag=ph_tag,
            )
            if int(owner_task) != int(getattr(dist_ctx, "task_id", -1)):
                skipped_owner += 1
                continue

        library_name = str(payload.get("library") or "")
        combo_post_dir = _resolve_combo_output_dir(
            cfg,
            root_key="POST_DOCKED_DIR",
            default_name="post_docked",
            run_id=str(run_id),
            pdb_id=str(pdb_id),
            variant_label=str(variant_token),
            ph_tag=ph_tag,
        )
        post_csv = _resolve_combo_post_consensus_csv(combo_post_dir, library_name)
        done_sentinel = combo_post_dir / "scorch" / "_DONE"
        has_post = post_csv.exists() and post_csv.stat().st_size > 0
        has_done = done_sentinel.exists()
        coverage_ok = False
        if has_post and has_done:
            coverage_ok = _scorch_coverage_complete_for_reconcile(
                cfg,
                run_id=str(run_id),
                pdb_id=str(pdb_id),
                variant_label=str(variant_token),
                ph_tag=ph_tag,
                logger=logger,
            )
        if has_post and has_done and coverage_ok:
            continue

        expected_raw = payload.get("expected")
        expected_keys: set[str] = set()
        if isinstance(expected_raw, set):
            expected_keys = {str(x) for x in expected_raw if str(x).strip()}
        elif isinstance(expected_raw, list):
            expected_keys = {str(x) for x in expected_raw if str(x).strip()}
        allowed_hint = len(expected_keys) if expected_keys else None
        if has_post and not has_done:
            try:
                post_keys = _load_post_scored_ligand_keys(post_csv)
            except Exception:
                post_keys = set()
            if (
                post_keys
                and (not expected_keys or expected_keys.issubset(post_keys))
                and _scorch_coverage_complete_for_reconcile(
                    cfg,
                    run_id=str(run_id),
                    pdb_id=str(pdb_id),
                    variant_label=str(variant_token),
                    ph_tag=ph_tag,
                    logger=logger,
                )
            ):
                done_sentinel.parent.mkdir(parents=True, exist_ok=True)
                done_sentinel.write_text("ok\n", encoding="utf-8")
                succeeded += 1
                logger.info(
                    "[scorch-rescore.reconcile] action=mark_done run_id=%s pdb_id=%s variant=%s ph=%s post_csv=%s scored=%d expected=%d",
                    str(run_id),
                    str(pdb_id),
                    str(variant_token),
                    str(ph_tag or "base"),
                    str(post_csv),
                    len(post_keys),
                    len(expected_keys),
                )
                continue
        variant_arg = None if variant_token in {"BASE"} else str(variant_token)
        attempted += 1
        logger.warning(
            "[scorch-rescore.reconcile] action=retry run_id=%s pdb_id=%s variant=%s ph=%s has_post=%s has_done=%s expected=%d mode=%s",
            str(run_id),
            str(pdb_id),
            str(variant_token),
            str(ph_tag or "base"),
            bool(has_post),
            bool(has_done),
            len(expected_keys),
            str(execution_mode),
        )
        rc = _maybe_run_scorch_rescore_for_pdb(
            cfg,
            str(run_id),
            str(pdb_id),
            variant=variant_arg,
            ph=ph_tag,
            allowed_count_hint=allowed_hint,
            execution_mode=execution_mode,
            overwrite=False,
            verbose=False,
        )
        has_post_now = post_csv.exists() and post_csv.stat().st_size > 0
        has_done_now = done_sentinel.exists()
        coverage_ok_now = (
            has_post_now
            and has_done_now
            and _scorch_coverage_complete_for_reconcile(
                cfg,
                run_id=str(run_id),
                pdb_id=str(pdb_id),
                variant_label=str(variant_token),
                ph_tag=ph_tag,
                logger=logger,
            )
        )
        if rc is not None and int(rc) == 0 and coverage_ok_now:
            succeeded += 1
            logger.info(
                "[scorch-rescore.reconcile] action=ok run_id=%s pdb_id=%s variant=%s ph=%s",
                str(run_id),
                str(pdb_id),
                str(variant_token),
                str(ph_tag or "base"),
            )
            continue
        failed += 1
        logger.warning(
            "[scorch-rescore.reconcile] action=failed run_id=%s pdb_id=%s variant=%s ph=%s returncode=%s has_post=%s has_done=%s",
            str(run_id),
            str(pdb_id),
            str(variant_token),
            str(ph_tag or "base"),
            str(rc),
            bool(has_post_now),
            bool(has_done_now),
        )

    return {
        "attempted": int(attempted),
        "succeeded": int(succeeded),
        "failed": int(failed),
        "skipped_owner": int(skipped_owner),
    }


def _distributed_manifest_dir(cfg: Mapping[str, Any], *, run_id: str) -> Path:
    return runtime_root(cfg, "MANIFESTS_DIR", "manifests") / str(run_id) / "distributed"


def _ensure_distributed_manifest_preflight(
    cfg: Mapping[str, Any],
    *,
    run_id: str,
    dist_ctx: Any,
    variant_label: str,
    combo_items: Sequence[tuple[str, Optional[str]]],
    run_tokens: list[str],
) -> None:
    if not _manifest_prebuild_enabled(cfg, distributed_enabled=bool(dist_ctx.enabled)):
        return
    logger = logging.getLogger("distributed.chunk")
    dist_dir = _distributed_manifest_dir(cfg, run_id=run_id)
    token = str(variant_label or "legacy").strip().lower().replace(" ", "_")
    state_path = dist_dir / f"manifest_preflight_{token}.json"
    lock_path = state_path.with_suffix(state_path.suffix + ".lock")

    try:
        payload = json.loads(state_path.read_text(encoding="utf-8")) or {}
    except Exception:
        payload = {}
    if isinstance(payload, dict) and str(payload.get("status") or "").strip().lower() == "ready":
        return

    got_lock = False
    try:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        os.close(fd)
        got_lock = True
    except FileExistsError:
        got_lock = False

    if got_lock:
        started = time.perf_counter()
        try:
            roots = _collect_combo_library_roots(
                cfg,
                combo_items=combo_items,
                run_tokens=run_tokens,
                logger=logger,
            )
            broken_roots, preflight_summary = _evaluate_manifest_preflight_roots(
                cfg,
                roots=roots,
                logger=logger,
            )
            if broken_roots:
                logger.info(
                    "[distributed.chunk.manifest-preflight.rebuild] roots_broken=%d action=targeted_rebuild",
                    len(broken_roots),
                )
                rebuild_summary = _rebuild_library_manifests(
                    cfg, roots=broken_roots, logger=logger
                )
            else:
                rebuild_summary = {
                    "roots_total": 0,
                    "built": 0,
                    "rebuilt": 0,
                    "unchanged": 0,
                    "missing_root": 0,
                    "failed": 0,
                }
            elapsed = time.perf_counter() - started
            logger.info(
                "[distributed.chunk.manifest-preflight] variant=%s roots=%d healthy=%d broken=%d built=%d rebuilt=%d unchanged=%d missing_root=%d failed=%d elapsed_s=%.2f",
                str(variant_label).upper(),
                int(preflight_summary.get("roots_total", 0)),
                int(preflight_summary.get("healthy", 0)),
                int(preflight_summary.get("broken_total", 0)),
                int(rebuild_summary.get("built", 0)),
                int(rebuild_summary.get("rebuilt", 0)),
                int(rebuild_summary.get("unchanged", 0)),
                int(preflight_summary.get("missing_root", 0)),
                int(rebuild_summary.get("failed", 0)),
                elapsed,
            )
            merged_summary = dict(preflight_summary)
            merged_summary.update(
                {
                    "rebuilt_roots_total": int(rebuild_summary.get("roots_total", 0)),
                    "built": int(rebuild_summary.get("built", 0)),
                    "rebuilt": int(rebuild_summary.get("rebuilt", 0)),
                    "unchanged": int(rebuild_summary.get("unchanged", 0)),
                    "rebuild_missing_root": int(rebuild_summary.get("missing_root", 0)),
                    "failed": int(rebuild_summary.get("failed", 0)),
                }
            )
            write_json_atomic(
                state_path,
                {
                    "run_id": str(run_id),
                    "variant_label": str(variant_label),
                    "status": "ready",
                    "task_id": int(dist_ctx.task_id),
                    "task_count": int(dist_ctx.task_count),
                    "updated_at": float(time.time()),
                    "elapsed_s": float(elapsed),
                    "summary": merged_summary,
                },
            )
            return
        finally:
            try:
                lock_path.unlink()
            except Exception:
                pass

    deadline = time.time() + _COMBO_PREP_WAIT_TIMEOUT_SEC
    wait_started = time.perf_counter()
    next_wait_log_sec = float(_MANIFEST_PREFLIGHT_WAIT_LOG_SEC)
    while time.time() < deadline:
        try:
            payload = json.loads(state_path.read_text(encoding="utf-8")) or {}
        except Exception:
            payload = {}
        if isinstance(payload, dict) and str(payload.get("status") or "").strip().lower() == "ready":
            return
        waited = time.perf_counter() - wait_started
        if waited >= next_wait_log_sec:
            logger.info(
                "[distributed.chunk.manifest-preflight.wait] variant=%s task_id=%d elapsed_s=%.1f state=waiting_for_ready",
                str(variant_label).upper(),
                int(getattr(dist_ctx, "task_id", 0)),
                waited,
            )
            next_wait_log_sec += float(_MANIFEST_PREFLIGHT_WAIT_LOG_SEC)
        time.sleep(_COMBO_PREP_WAIT_POLL_SEC)

    logger.warning(
        "[distributed.chunk.manifest-preflight] variant=%s timeout_s=%.1f action=local_fallback",
        str(variant_label).upper(),
        float(_COMBO_PREP_WAIT_TIMEOUT_SEC),
    )
    roots = _collect_combo_library_roots(
        cfg,
        combo_items=combo_items,
        run_tokens=run_tokens,
        logger=logger,
    )
    broken_roots, _ = _evaluate_manifest_preflight_roots(
        cfg,
        roots=roots,
        logger=logger,
    )
    if broken_roots:
        _rebuild_library_manifests(cfg, roots=broken_roots, logger=logger)
