from __future__ import annotations

import logging
import os
import sys
import threading
from dataclasses import dataclass
from typing import Any, Callable, Optional

from tqdm import tqdm  # type: ignore[import-untyped]

from cli.run_context import ConfigDict


@dataclass
class VariantExecutionContext:
    run_id: str
    cfg: ConfigDict
    variants: list[Optional[str]]
    pdb_files: list[str]
    tokens: list[str]
    mode: str
    is_resume: bool
    dist_ctx: Any
    dist_combo_chunk_mode: bool
    dist_combo_hybrid_mode: bool
    execution_mode_decision: Any
    stages: Any
    params: Any
    global_start: float


@dataclass
class VariantExecutionSharedState:
    failed_root: str
    completed_combo_lookup: set[tuple[str, str, str]]
    failed_entries: list[tuple[str, str, str, str, str]]
    run_scope_completion_times: dict[tuple[str, str, str], float]
    recorded_terminal_chunk_ids: set[str]
    distributed_assigned_pdb_ids: set[str]
    run_chunk_rebalance_count: int
    retention_lock: threading.Lock
    pending_retention_combos: set[tuple[str, str, str]]
    pending_coverage_refresh: dict[tuple[str, str, str], dict[str, Any]]
    pending_retention_lock: threading.Lock
    scorch_queue_service: Any


@dataclass
class VariantExecutionDeps:
    to_bool: Callable[[Any], bool]
    build_combo_work_items: Callable[..., list[tuple[str, Optional[str]]]]
    resolve_global_scheduler_plan: Callable[..., dict[str, Any]]
    build_process_one_runner: Callable[..., Callable[[str, ConfigDict, Optional[str]], bool]]
    process_one_context_type: Any
    process_one_shared_state_type: Any
    process_one_deps_type: Any
    normalize_ph_tag_token: Callable[[Optional[str]], str]
    chunk_ligand_key: Callable[[str], str]
    process_one_protein: Callable[..., Any]
    update_manifest_for_protein_start: Callable[..., Any]
    update_manifest_for_protein_success: Callable[..., Any]
    update_manifest_for_protein_failure: Callable[..., Any]
    verify_chunk_combo_outputs: Callable[..., Any]
    resolve_combo_output_dir: Callable[..., Any]
    resolve_combo_docking_summary_csv: Callable[..., Any]
    load_scored_ligand_keys_from_summary: Callable[..., Any]
    update_combo_coverage_snapshot: Callable[..., Any]
    resolve_combo_post_consensus_csv: Callable[..., Any]
    load_post_scored_ligand_keys: Callable[..., Any]
    acquire_global_slot: Callable[..., Any]
    maybe_run_scorch_rescore_for_pdb: Callable[..., Any]
    maybe_run_artifact_retention_for_combo: Callable[..., Any]
    run_distributed_combo_variant: Callable[..., int]
    run_multi_pdb_single_ligand: Callable[..., None]
    run_global_scheduler: Callable[..., ConfigDict]
    run_serial_scheduler: Callable[..., None]
    apply_pocket_detection_events: Callable[..., Any]


@dataclass
class VariantExecutionResult:
    run_chunk_rebalance_count: int


def run_variant_execution(
    *,
    context: VariantExecutionContext,
    shared_state: VariantExecutionSharedState,
    deps: VariantExecutionDeps,
) -> VariantExecutionResult:
    for variant in context.variants:
        # Make variant visible to any module still reading env (legacy compatibility)
        os.environ["APO_HOLO_VARIANT"] = "" if variant is None else str(variant).upper()
        label = "legacy" if variant is None else str(variant).lower()
        logging.info(
            "[apo-holo] start_variant mode=%s variant_label=%s env_token=%s",
            context.mode,
            label,
            os.environ.get("APO_HOLO_VARIANT", ""),
        )
        variant_token = None if variant is None else str(variant).upper()

        with tqdm(
            total=0,
            desc=f"Processing Combos ({label})",
            unit="combo",
            position=0,
            dynamic_ncols=True,
            mininterval=0.2,
            leave=True,
            file=sys.stdout,
        ) as bar:
            cfg_v = context.cfg  # no per-variant mutation; variant is in env
            combo_items = deps.build_combo_work_items(
                cfg_v,
                context.pdb_files,
                variant_token=variant_token,
            )
            bar.reset(total=len(combo_items))
            logging.info(
                "[main.combo-plan] variant=%s proteins=%d combos=%d",
                label.upper(),
                len(context.pdb_files),
                len(combo_items),
            )

            single_ligand_mode = bool(context.cfg.get("_EFFECTIVE_SINGLE_LIGAND"))
            multi_pdb_single_ligand = single_ligand_mode and len(combo_items) > 1
            cpu = int(context.cfg.get("CPU", os.cpu_count() or 1))
            max_pdb_workers = max(1, min(cpu, len(combo_items))) if combo_items else 1
            scheduler_plan = deps.resolve_global_scheduler_plan(
                context.cfg,
                pdb_count=len(combo_items),
                cpu=cpu,
            )
            scheduler_enabled = bool(scheduler_plan["enabled"])
            scheduler_policy = str(scheduler_plan["policy"])
            scheduler_cpus = int(scheduler_plan["scheduler_cpus"])
            min_parallel_proteins = int(scheduler_plan["min_parallel"])
            max_parallel_proteins = int(scheduler_plan["max_parallel"])
            # Same PDB/variant pH members still share several prep-side work files
            # and ligand/control directories.  Keep them serialized by default to
            # avoid cross-pH file truncation races; advanced runs can opt into the
            # faster mode once their artifacts are fully pH-scoped.
            local_ph_hybrid_unlock = bool(
                deps.to_bool(context.cfg.get("LOCAL_PH_HYBRID_UNLOCK", False))
            )
            if not context.dist_combo_chunk_mode and len(combo_items) > 1:
                logging.info(
                    "[main.parallel.hybrid] mode=%s enabled=%s",
                    label.upper(),
                    str(local_ph_hybrid_unlock).lower(),
                )

            variant_small_task_cfg: ConfigDict = cfg_v
            process_one_context = deps.process_one_context_type(
                run_id=str(context.run_id),
                label=str(label),
                mode=str(context.mode),
                variant=None if variant is None else str(variant),
                is_resume=bool(context.is_resume),
                completed_combo_lookup=shared_state.completed_combo_lookup,
                dist_combo_chunk_mode=bool(context.dist_combo_chunk_mode),
                cfg_v=cfg_v,
                stages=context.stages,
                params=context.params,
                tokens=list(context.tokens),
                global_start=float(context.global_start),
            )
            process_one_shared_state = deps.process_one_shared_state_type(
                failed_root=str(shared_state.failed_root),
                failed_entries=shared_state.failed_entries,
                run_scope_completion_times=shared_state.run_scope_completion_times,
                retention_lock=shared_state.retention_lock,
                pending_retention_combos=shared_state.pending_retention_combos,
                pending_coverage_refresh=shared_state.pending_coverage_refresh,
                pending_retention_lock=shared_state.pending_retention_lock,
                scorch_queue_service=shared_state.scorch_queue_service,
            )
            process_one_deps = deps.process_one_deps_type(
                normalize_ph_tag_token=deps.normalize_ph_tag_token,
                chunk_ligand_key=deps.chunk_ligand_key,
                process_one_protein=deps.process_one_protein,
                update_manifest_for_protein_start=deps.update_manifest_for_protein_start,
                update_manifest_for_protein_success=deps.update_manifest_for_protein_success,
                update_manifest_for_protein_failure=deps.update_manifest_for_protein_failure,
                verify_chunk_combo_outputs=deps.verify_chunk_combo_outputs,
                resolve_combo_output_dir=deps.resolve_combo_output_dir,
                resolve_combo_docking_summary_csv=deps.resolve_combo_docking_summary_csv,
                load_scored_ligand_keys_from_summary=deps.load_scored_ligand_keys_from_summary,
                update_combo_coverage_snapshot=deps.update_combo_coverage_snapshot,
                resolve_combo_post_consensus_csv=deps.resolve_combo_post_consensus_csv,
                load_post_scored_ligand_keys=deps.load_post_scored_ligand_keys,
                acquire_global_slot=deps.acquire_global_slot,
                maybe_run_scorch_rescore_for_pdb=deps.maybe_run_scorch_rescore_for_pdb,
                maybe_run_artifact_retention_for_combo=deps.maybe_run_artifact_retention_for_combo,
            )
            process_one_runner = deps.build_process_one_runner(
                context=process_one_context,
                shared_state=process_one_shared_state,
                deps=process_one_deps,
            )

            if context.dist_combo_chunk_mode:
                shared_state.run_chunk_rebalance_count = deps.run_distributed_combo_variant(
                    cfg_v=cfg_v,
                    dist_ctx=context.dist_ctx,
                    run_id=context.run_id,
                    label=label,
                    combo_items=list(combo_items),
                    tokens=list(context.tokens),
                    bar=bar,
                    variant=None if variant is None else str(variant),
                    dist_combo_hybrid_mode=context.dist_combo_hybrid_mode,
                    execution_mode_decision=context.execution_mode_decision,
                    scorch_queue_service=shared_state.scorch_queue_service,
                    process_one=process_one_runner,
                    recorded_terminal_chunk_ids=shared_state.recorded_terminal_chunk_ids,
                    failed_entries=shared_state.failed_entries,
                    run_scope_completion_times=shared_state.run_scope_completion_times,
                    distributed_assigned_pdb_ids=shared_state.distributed_assigned_pdb_ids,
                    run_chunk_rebalance_count=shared_state.run_chunk_rebalance_count,
                )
            elif multi_pdb_single_ligand:
                deps.run_multi_pdb_single_ligand(
                    cfg_v=cfg_v,
                    combo_items=list(combo_items),
                    max_pdb_workers=max_pdb_workers,
                    variant_token=variant_token,
                    variant_label=label,
                    local_ph_hybrid_unlock=local_ph_hybrid_unlock,
                    process_one=process_one_runner,
                    bar=bar,
                )
            elif scheduler_enabled and len(combo_items) > 1:
                variant_small_task_cfg = deps.run_global_scheduler(
                    cfg_v=cfg_v,
                    combo_items=list(combo_items),
                    variant_token=variant_token,
                    variant_label=label,
                    scheduler_plan=scheduler_plan,
                    scheduler_policy=scheduler_policy,
                    scheduler_cpus=scheduler_cpus,
                    min_parallel_proteins=min_parallel_proteins,
                    max_parallel_proteins=max_parallel_proteins,
                    local_ph_hybrid_unlock=local_ph_hybrid_unlock,
                    process_one=process_one_runner,
                    bar=bar,
                )
            else:
                deps.run_serial_scheduler(
                    cfg_v=cfg_v,
                    combo_items=list(combo_items),
                    process_one=process_one_runner,
                    bar=bar,
                )

            try:
                with deps.acquire_global_slot(
                    variant_small_task_cfg,
                    cores=1,
                    min_cores=1,
                    priority=0,
                    task_id=f"main:housekeeping:pocket-events:{label}",
                    task_type="housekeeping",
                ):
                    deps.apply_pocket_detection_events(context.cfg, context.run_id)
            except Exception:
                logging.warning(
                    "[run-manifest.pocket_detection.events.apply] run_id=%s action=skip",
                    context.run_id,
                    exc_info=True,
                )

    return VariantExecutionResult(
        run_chunk_rebalance_count=int(shared_state.run_chunk_rebalance_count)
    )
