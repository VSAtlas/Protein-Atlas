from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Set, cast

from docking.global_scheduler import acquire_global_cores
from post_docking.rescoring.scorch_orchestration_types import ComboKey, OrchestrationDeps
from post_docking.rescoring.scorch_shards import (
    planned_completed_output_bases,
    planned_completed_output_csvs,
)
from post_docking.rescoring.scorch_types import StageSpec


def finalize_orchestration_combos(
    *,
    cfg: dict[str, object],
    combos: Set[ComboKey],
    combos_with_tasks: Set[ComboKey],
    combo_modes: Dict[ComboKey, Set[str]],
    combo_failed_local: Dict[ComboKey, bool],
    specs: List[StageSpec],
    run_root: Path,
    post_run_root: Path,
    overwrite: bool,
    decoy_prefix: str,
    component: str,
    logger: logging.Logger,
    deps: OrchestrationDeps,
    total_jobs: int,
    completed: int,
    failed_jobs: int,
    combos_attempted: int,
    skipped_missing_receptor: int,
    skipped_missing_consensus: int,
) -> Set[ComboKey]:
    local_completed_combos: Set[ComboKey] = set()
    prefix_value = decoy_prefix
    for combo in sorted(combos_with_tasks):
        modes = combo_modes.get(combo, {"fda"})
        fda_all = None
        dud_all = None
        plan_bound_modes: Set[str] = set()
        plan_incomplete_modes: Set[str] = set()
        for mode in sorted(modes):
            out_name = (
                "scorch_scores_all.csv"
                if mode == "fda"
                else f"{prefix_value}_scorch_scores_all.csv"
            )
            pdb_id, variant, ph = combo
            planned_csvs = planned_completed_output_csvs(
                cfg,
                run_root.name,
                combo,
                specs,
                mode,
                decoy_prefix,
            )
            planned_bases = None
            if planned_csvs is not None:
                plan_bound_modes.add(mode)
                if not planned_csvs:
                    plan_incomplete_modes.add(mode)
                else:
                    planned_bases = planned_completed_output_bases(
                        cfg,
                        run_root.name,
                        combo,
                        specs,
                        mode,
                        decoy_prefix,
                    )
            agg_task_id = f"scorch:aggregate:{pdb_id}:{variant}:{ph}:{mode}"
            deps.emit_task_event(
                "TASK_SUBMITTED",
                task_id=agg_task_id,
                task_type="consensus_rerank_annotate",
                want_cores=1,
                min_cores=1,
                est_duration_sec=1.0,
                pdb_id=pdb_id,
                variant=variant,
                ph=ph,
            )
            with acquire_global_cores(
                cfg,
                cores=1,
                min_cores=1,
                priority=1,
                est_duration_sec=1.0,
                task_id=agg_task_id,
                task_type="consensus_rerank_annotate",
            ) as granted:
                deps.emit_task_event(
                    "TASK_START",
                    task_id=agg_task_id,
                    task_type="consensus_rerank_annotate",
                    want_cores=1,
                    min_cores=1,
                    granted_cores=int(granted),
                    est_duration_sec=1.0,
                    pdb_id=pdb_id,
                    variant=variant,
                    ph=ph,
                )
                agg_path = deps.aggregate_combo(
                    post_run_root,
                    specs,
                    combo,
                    logger,
                    run_mode=mode,
                    output_name=out_name,
                    input_csvs=planned_csvs,
                    allowed_ligand_bases=planned_bases,
                )
                deps.emit_task_event(
                    "TASK_END",
                    task_id=agg_task_id,
                    task_type="consensus_rerank_annotate",
                    want_cores=1,
                    min_cores=1,
                    granted_cores=int(granted),
                    est_duration_sec=1.0,
                    status="ok" if agg_path is not None else "error",
                    pdb_id=pdb_id,
                    variant=variant,
                    ph=ph,
                )
            if mode == "fda":
                fda_all = agg_path
            elif mode == "dud":
                dud_all = agg_path

        if plan_incomplete_modes:
            logger.warning(
                "%s action=finalize status=incomplete reason=plan_bound_mode_missing_outputs pdb_id=%s variant=%s ph=%s modes=%s",
                component,
                combo[0],
                combo[1],
                combo[2],
                ",".join(sorted(plan_incomplete_modes)),
            )

        if fda_all is None and "fda" not in plan_bound_modes:
            candidate = post_run_root / combo[0] / combo[1] / combo[2] / "scorch_scores_all.csv"
            if candidate.exists():
                fda_all = candidate
        if dud_all is None and "dud" not in plan_bound_modes:
            candidate = (
                post_run_root
                / combo[0]
                / combo[1]
                / combo[2]
                / f"{prefix_value}_scorch_scores_all.csv"
            )
            if candidate.exists():
                dud_all = candidate

        chosen_all = fda_all or dud_all
        if chosen_all is None:
            for mode in sorted(modes):
                if mode in plan_bound_modes:
                    continue
                alt = (
                    "scorch_scores_all.csv"
                    if mode == "fda"
                    else f"{prefix_value}_scorch_scores_all.csv"
                )
                candidate = post_run_root / combo[0] / combo[1] / combo[2] / alt
                if candidate.exists():
                    chosen_all = candidate
                    break

        raw_complete = bool(
            chosen_all
            and not combo_failed_local.get(combo, False)
            and not plan_incomplete_modes
        )
        rerank_ok = False if raw_complete else True

        if (
            chosen_all
            and deps.rerank_consensus_with_scorch is not None
            and deps.find_consensus_csv is not None
            and not plan_incomplete_modes
        ):
            try:
                pdb_id, variant, ph = combo
                dock_combo_dir = run_root / pdb_id / variant / ph
                consensus_csv = deps.find_consensus_csv(dock_combo_dir)
                if not consensus_csv:
                    logger.warning(
                        "%s action=rerank status=skip reason=missing_consensus pdb_id=%s variant=%s ph=%s dock_dir=%s",
                        component,
                        pdb_id,
                        variant,
                        ph,
                        dock_combo_dir,
                    )
                    rerank_ok = False
                else:
                    out_csv = chosen_all.parent / "consensus_reranked_scorch.csv"
                    scorch_inputs = [p for p in (fda_all, dud_all) if p is not None]
                    if not scorch_inputs:
                        scorch_inputs = [chosen_all]
                    rerank_task_id = f"scorch:rerank:{pdb_id}:{variant}:{ph}:{prefix_value}"
                    deps.emit_task_event(
                        "TASK_SUBMITTED",
                        task_id=rerank_task_id,
                        task_type="consensus_rerank_annotate",
                        want_cores=1,
                        min_cores=1,
                        est_duration_sec=1.0,
                        pdb_id=pdb_id,
                        variant=variant,
                        ph=ph,
                    )
                    with acquire_global_cores(
                        cfg,
                        cores=1,
                        min_cores=1,
                        priority=1,
                        est_duration_sec=1.0,
                        task_id=rerank_task_id,
                        task_type="consensus_rerank_annotate",
                    ) as granted:
                        deps.emit_task_event(
                            "TASK_START",
                            task_id=rerank_task_id,
                            task_type="consensus_rerank_annotate",
                            want_cores=1,
                            min_cores=1,
                            granted_cores=int(granted),
                            est_duration_sec=1.0,
                            pdb_id=pdb_id,
                            variant=variant,
                            ph=ph,
                        )
                        reranked = bool(deps.rerank_consensus_with_scorch(
                            consensus_csv,
                            scorch_inputs,
                            out_csv,
                            logger,
                            overwrite=True,
                            decoy_prefix=prefix_value,
                        ))
                        expected_outputs = [out_csv]
                        if dud_all is not None:
                            expected_outputs.append(
                                out_csv.with_name(
                                    f"{prefix_value}_consensus_reranked_scorch.csv"
                                )
                            )
                            if str(prefix_value).lower() != "dud":
                                expected_outputs.append(
                                    out_csv.with_name(
                                        "dud_consensus_reranked_scorch.csv"
                                    )
                                )
                        missing_outputs = [
                            path
                            for path in expected_outputs
                            if not path.exists() or path.stat().st_size == 0
                        ]
                        rerank_ok = bool(reranked and not missing_outputs)
                        deps.emit_task_event(
                            "TASK_END",
                            task_id=rerank_task_id,
                            task_type="consensus_rerank_annotate",
                            want_cores=1,
                            min_cores=1,
                            granted_cores=int(granted),
                            est_duration_sec=1.0,
                            status="ok" if rerank_ok else "error",
                            pdb_id=pdb_id,
                            variant=variant,
                            ph=ph,
                        )
                        if not rerank_ok:
                            logger.warning(
                                "%s action=rerank status=failed combo=%s returned=%s missing_outputs=%s",
                                component,
                                combo,
                                reranked,
                                ",".join(str(path) for path in missing_outputs),
                            )
            except Exception as exc:
                rerank_ok = False
                logger.warning(
                    "%s action=rerank status=skip reason=exception combo=%s error=%s",
                    component,
                    combo,
                    exc,
                )
        elif chosen_all and raw_complete:
            logger.warning(
                "%s action=rerank status=skip reason=reranker_unavailable combo=%s",
                component,
                combo,
            )

        if fda_all and dud_all:
            annot_task_id = f"scorch:annotate:{combo[0]}:{combo[1]}:{combo[2]}"
            deps.emit_task_event(
                "TASK_SUBMITTED",
                task_id=annot_task_id,
                task_type="consensus_rerank_annotate",
                want_cores=1,
                min_cores=1,
                est_duration_sec=1.0,
                pdb_id=combo[0],
                variant=combo[1],
                ph=combo[2],
            )
            with acquire_global_cores(
                cfg,
                cores=1,
                min_cores=1,
                priority=1,
                est_duration_sec=1.0,
                task_id=annot_task_id,
                task_type="consensus_rerank_annotate",
            ) as granted:
                deps.emit_task_event(
                    "TASK_START",
                    task_id=annot_task_id,
                    task_type="consensus_rerank_annotate",
                    want_cores=1,
                    min_cores=1,
                    granted_cores=int(granted),
                    est_duration_sec=1.0,
                    pdb_id=combo[0],
                    variant=combo[1],
                    ph=combo[2],
                )
                deps.annotate_scorch_z_scores(cast(Path, fda_all), cast(Path, dud_all), logger)
                deps.emit_task_event(
                    "TASK_END",
                    task_id=annot_task_id,
                    task_type="consensus_rerank_annotate",
                    want_cores=1,
                    min_cores=1,
                    granted_cores=int(granted),
                    est_duration_sec=1.0,
                    status="ok",
                    pdb_id=combo[0],
                    variant=combo[1],
                    ph=combo[2],
                )

        if raw_complete and rerank_ok:
            local_completed_combos.add(combo)

    logger.info(
        "%s action=summary status=%s combos=%d combos_attempted=%d total_jobs=%d completed=%d failed=%d missing_receptor=%d missing_consensus=%d decoy_prefix=%s",
        component,
        "ok" if failed_jobs == 0 else "failed",
        len(combos),
        combos_attempted,
        total_jobs,
        completed,
        failed_jobs,
        skipped_missing_receptor,
        skipped_missing_consensus,
        prefix_value,
    )
    return local_completed_combos
