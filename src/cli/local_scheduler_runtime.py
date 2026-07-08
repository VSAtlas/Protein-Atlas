from __future__ import annotations

import logging
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from typing import Any, Callable, Optional

from cli.distributed_chunk_planner import (
    _combo_pdb_id_from_file,
    _local_combo_scope_key,
    _normalize_ph_tag_token,
    _pick_next_combo_index,
    _refresh_local_prep_ready_scopes,
)
from cli.run_context import ConfigDict
from docking.global_scheduler import GlobalAdmissionScheduler, sem_free_slots


def run_multi_pdb_single_ligand(
    *,
    cfg_v: ConfigDict,
    combo_items: list[tuple[str, Optional[str]]],
    max_pdb_workers: int,
    variant_token: Optional[str],
    variant_label: str,
    local_ph_hybrid_unlock: bool,
    process_one: Callable[[str, ConfigDict, Optional[str]], bool],
    bar: Any,
) -> None:
    logging.info(
        "[main.parallel.single_ligand] mode=%s n_items=%d max_workers=%d",
        variant_label.upper(),
        len(combo_items),
        max_pdb_workers,
    )
    local_prep_ready_floor = time.time()

    def _cfg_for_pdb() -> ConfigDict:
        cfg_local = cfg_v.copy()
        cfg_local["CPU"] = 1
        return cfg_local

    with ThreadPoolExecutor(max_workers=max_pdb_workers) as pool:
        pending_single = list(combo_items)
        running_single: dict[Any, tuple[str, Optional[str], str]] = {}
        inflight_single_scope_counts: dict[str, int] = {}
        prep_ready_single_scopes: set[str] = set()
        prep_scope_ready_after_single: dict[str, float] = {}

        def _refresh_single_ready_scopes() -> None:
            probe_items = list(pending_single)
            probe_items.extend(
                (pdb_file, ph_tag)
                for pdb_file, ph_tag, _scope_key in running_single.values()
            )
            _refresh_local_prep_ready_scopes(
                cfg_v,
                variant_token=variant_token,
                variant_label=variant_label,
                combos=probe_items,
                ready_scopes=prep_ready_single_scopes,
                scope_ready_after=prep_scope_ready_after_single,
                default_ready_after=local_prep_ready_floor,
            )

        def _refill_single_ligand() -> None:
            _refresh_single_ready_scopes()
            while pending_single and len(running_single) < max_pdb_workers:
                idx = _pick_next_combo_index(
                    pending_single,
                    inflight_scope_counts=inflight_single_scope_counts,
                    prep_ready_scopes=prep_ready_single_scopes,
                    allow_same_scope_when_ready=local_ph_hybrid_unlock,
                    variant_label=variant_label,
                )
                if idx is None:
                    return
                pdb_file, ph_tag = pending_single.pop(idx)
                pdb_id = _combo_pdb_id_from_file(pdb_file)
                scope_key = _local_combo_scope_key(
                    pdb_file=pdb_file,
                    variant_label=variant_label,
                )
                fut = pool.submit(process_one, pdb_file, _cfg_for_pdb(), ph_tag)
                running_single[fut] = (pdb_file, ph_tag, scope_key)
                prev_scope_inflight = int(
                    inflight_single_scope_counts.get(scope_key, 0) or 0
                )
                inflight_single_scope_counts[scope_key] = prev_scope_inflight + 1
                if scope_key not in prep_ready_single_scopes:
                    prep_scope_ready_after_single[scope_key] = time.time()
                logging.info(
                    "[main.parallel.single_ligand.schedule] action=submit pdb=%s ph=%s scope=%s scope_in_flight=%d prep_ready=%s total_in_flight=%d pending=%d",
                    pdb_id,
                    _normalize_ph_tag_token(ph_tag),
                    scope_key,
                    int(inflight_single_scope_counts.get(scope_key, 0)),
                    str(scope_key in prep_ready_single_scopes).lower(),
                    len(running_single),
                    len(pending_single),
                )

        _refill_single_ligand()
        while running_single:
            done, _ = wait(
                list(running_single.keys()),
                timeout=0.5,
                return_when=FIRST_COMPLETED,
            )
            if not done:
                _refill_single_ligand()
                continue
            for fut in done:
                pdb_file, ph_tag, scope_key = running_single.pop(fut)
                try:
                    fut.result()
                except Exception as exc:
                    logging.exception(
                        "[main.parallel.single_ligand.error] pdb=%s ph=%s error=%s",
                        pdb_file,
                        _normalize_ph_tag_token(ph_tag),
                        exc,
                    )
                finally:
                    prev = int(inflight_single_scope_counts.get(scope_key, 0) or 0)
                    if prev <= 1:
                        inflight_single_scope_counts.pop(scope_key, None)
                    else:
                        inflight_single_scope_counts[scope_key] = prev - 1
                    bar.update(1)
            _refill_single_ligand()


def run_global_scheduler(
    *,
    cfg_v: ConfigDict,
    combo_items: list[tuple[str, Optional[str]]],
    variant_token: Optional[str],
    variant_label: str,
    scheduler_plan: dict[str, Any],
    scheduler_policy: str,
    scheduler_cpus: int,
    min_parallel_proteins: int,
    max_parallel_proteins: int,
    local_ph_hybrid_unlock: bool,
    process_one: Callable[[str, ConfigDict, Optional[str]], bool],
    bar: Any,
) -> ConfigDict:
    global_admission = GlobalAdmissionScheduler(max(1, scheduler_cpus))
    target_parallel = int(scheduler_plan["initial_target"])
    logging.info(
        "[main.parallel.global] mode=%s policy=%s scheduler=central_backfill scheduler_cpus=%d max_parallel=%d target_initial=%d n_items=%d",
        variant_label.upper(),
        scheduler_policy,
        scheduler_cpus,
        max_parallel_proteins,
        target_parallel,
        len(combo_items),
    )
    local_prep_ready_floor = time.time()

    pending_global = list(combo_items)
    running_global: dict[Any, tuple[str, Optional[str], str]] = {}
    inflight_global_scope_counts: dict[str, int] = {}
    prep_ready_global_scopes: set[str] = set()
    prep_scope_ready_after_global: dict[str, float] = {}

    def _cfg_for_scheduler() -> ConfigDict:
        cfg_local = cfg_v.copy()
        cfg_local["GLOBAL_DOCKING_SCHEDULER"] = global_admission
        cfg_local["GLOBAL_DOCKING_SEM"] = global_admission
        # Keep compatibility with existing benchmark-style semaphore key.
        cfg_local["GLOBAL_LIGAND_SEM"] = global_admission
        cfg_local["CPU"] = max(int(cfg_local.get("CPU", 1) or 1), scheduler_cpus)
        return cfg_local

    with ThreadPoolExecutor(max_workers=max_parallel_proteins) as pool:

        def _refresh_global_ready_scopes() -> None:
            probe_items = list(pending_global)
            probe_items.extend(
                (pdb_file, ph_tag)
                for pdb_file, ph_tag, _scope_key in running_global.values()
            )
            _refresh_local_prep_ready_scopes(
                cfg_v,
                variant_token=variant_token,
                variant_label=variant_label,
                combos=probe_items,
                ready_scopes=prep_ready_global_scopes,
                scope_ready_after=prep_scope_ready_after_global,
                default_ready_after=local_prep_ready_floor,
            )

        def _refill() -> None:
            nonlocal target_parallel
            _refresh_global_ready_scopes()
            while pending_global and len(running_global) < target_parallel:
                idx = _pick_next_combo_index(
                    pending_global,
                    inflight_scope_counts=inflight_global_scope_counts,
                    prep_ready_scopes=prep_ready_global_scopes,
                    allow_same_scope_when_ready=local_ph_hybrid_unlock,
                    variant_label=variant_label,
                )
                if idx is None:
                    return
                pdb_file, ph_tag = pending_global.pop(idx)
                pdb_id = _combo_pdb_id_from_file(pdb_file)
                scope_key = _local_combo_scope_key(
                    pdb_file=pdb_file,
                    variant_label=variant_label,
                )
                fut = pool.submit(
                    process_one,
                    pdb_file,
                    _cfg_for_scheduler(),
                    ph_tag,
                )
                running_global[fut] = (pdb_file, ph_tag, scope_key)
                prev_scope_inflight = int(
                    inflight_global_scope_counts.get(scope_key, 0) or 0
                )
                inflight_global_scope_counts[scope_key] = prev_scope_inflight + 1
                if scope_key not in prep_ready_global_scopes:
                    prep_scope_ready_after_global[scope_key] = time.time()
                logging.info(
                    "[main.parallel.global.schedule] action=submit pdb=%s ph=%s scope=%s scope_in_flight=%d prep_ready=%s total_in_flight=%d pending=%d target=%d",
                    pdb_id,
                    _normalize_ph_tag_token(ph_tag),
                    scope_key,
                    int(inflight_global_scope_counts.get(scope_key, 0)),
                    str(scope_key in prep_ready_global_scopes).lower(),
                    len(running_global),
                    len(pending_global),
                    target_parallel,
                )

        variant_small_task_cfg = _cfg_for_scheduler()
        _refill()
        while running_global:
            done, _ = wait(
                list(running_global.keys()),
                timeout=0.5,
                return_when=FIRST_COMPLETED,
            )
            if not done:
                _refill()
                continue
            for fut in done:
                pdb_file, ph_tag, scope_key = running_global.pop(fut)
                try:
                    fut.result()
                except Exception as exc:
                    logging.exception(
                        "[main.parallel.global.error] pdb=%s ph=%s error=%s",
                        pdb_file,
                        _normalize_ph_tag_token(ph_tag),
                        exc,
                    )
                finally:
                    prev = int(inflight_global_scope_counts.get(scope_key, 0) or 0)
                    if prev <= 1:
                        inflight_global_scope_counts.pop(scope_key, None)
                    else:
                        inflight_global_scope_counts[scope_key] = prev - 1
                    bar.update(1)

            if scheduler_policy == "adaptive":
                free_slots = sem_free_slots(global_admission, scheduler_cpus)
                grow_threshold = max(1, scheduler_cpus // 3)
                if free_slots >= grow_threshold:
                    target_parallel = min(max_parallel_proteins, target_parallel + 1)
                elif free_slots == 0 and target_parallel > min_parallel_proteins:
                    target_parallel = max(min_parallel_proteins, target_parallel - 1)

            _refill()

    return variant_small_task_cfg


def run_serial_scheduler(
    *,
    cfg_v: ConfigDict,
    combo_items: list[tuple[str, Optional[str]]],
    process_one: Callable[[str, ConfigDict, Optional[str]], bool],
    bar: Any,
) -> None:
    for pdb_file, ph_tag in combo_items:
        try:
            process_one(pdb_file, cfg_v, ph_tag)
        finally:
            # Always advance the progress bar, even if this PDB failed
            bar.update(1)
