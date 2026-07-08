from __future__ import annotations

import logging
import os
from typing import Any, Callable, Mapping, Sequence

from cli.distributed_chunk_planner import _normalize_ph_tag_token, _resolve_combo_output_dir
from cli.distributed_chunk_runtime_state import chunk_terminal
from cli.hybrid_chunk_scheduler import chunk_scope_key
from docking.completion_markers import aggregate_combo_chunk_markers


ScopeKey = tuple[str, str, str]


def chunk_payload_scope(
    chunk_payload: Mapping[str, Any], *, variant_label: str
) -> ScopeKey:
    pdb_scope = (
        str(chunk_payload.get("pdb_id") or "").strip().upper()
        or os.path.splitext(os.path.basename(str(chunk_payload.get("pdb_file") or "")))[
            0
        ].upper()
    )
    ph_scope = _normalize_ph_tag_token(chunk_payload.get("ph_tag"))
    return (
        str(pdb_scope),
        str(variant_label or "").strip().upper() or "BASE",
        str(ph_scope),
    )


def build_scope_maps(
    *,
    chunk_plan: Sequence[Mapping[str, Any]],
    label: str,
    dist_ctx: Any,
    owner_task_resolver: Callable[..., int],
) -> tuple[
    dict[str, Mapping[str, Any]],
    dict[ScopeKey, int],
    dict[ScopeKey, int],
    dict[ScopeKey, int],
]:
    by_chunk_id = {
        str(item.get("chunk_id")): item
        for item in chunk_plan
        if str(item.get("chunk_id") or "").strip()
    }
    combo_owner_task_id: dict[ScopeKey, int] = {}
    combo_total_chunks: dict[ScopeKey, int] = {}
    combo_total_ligands: dict[ScopeKey, int] = {}

    for payload in by_chunk_id.values():
        scope = chunk_payload_scope(payload, variant_label=label)
        combo_total_chunks[scope] = int(combo_total_chunks.get(scope, 0)) + 1
        bases_raw = payload.get("ligand_bases") or []
        bases_n = sum(1 for x in bases_raw if str(x).strip())
        combo_total_ligands[scope] = int(combo_total_ligands.get(scope, 0)) + int(bases_n)

    try:
        task_ids = sorted({int(x) for x in dist_ctx.expected_task_ids()})
    except Exception:
        task_ids = []
    if len(task_ids) > 1 and combo_total_chunks:
        task_load = {int(task_id): 0 for task_id in task_ids}
        prep_group_scopes: dict[tuple[str, str], list[ScopeKey]] = {}
        prep_group_chunks: dict[tuple[str, str], int] = {}
        prep_group_ligands: dict[tuple[str, str], int] = {}
        for scope in sorted(combo_total_chunks):
            prep_group = (str(scope[0]), str(scope[1]))
            prep_group_scopes.setdefault(prep_group, []).append(scope)
            prep_group_chunks[prep_group] = int(prep_group_chunks.get(prep_group, 0)) + int(
                combo_total_chunks.get(scope, 0)
            )
            prep_group_ligands[prep_group] = int(prep_group_ligands.get(prep_group, 0)) + int(
                combo_total_ligands.get(scope, 0)
            )
        for prep_group in sorted(
            prep_group_scopes,
            key=lambda item: (
                -int(prep_group_chunks.get(item, 0)),
                -int(prep_group_ligands.get(item, 0)),
                item,
            ),
        ):
            owner_task = min(task_ids, key=lambda task_id: (task_load[int(task_id)], int(task_id)))
            for scope in prep_group_scopes.get(prep_group, []):
                combo_owner_task_id[scope] = int(owner_task)
            task_load[int(owner_task)] += int(prep_group_chunks.get(prep_group, 0))
    else:
        for scope in combo_total_chunks:
            combo_owner_task_id[scope] = owner_task_resolver(
                dist_ctx,
                variant_label=str(label),
                pdb_id=str(scope[0]),
                ph_tag=None if str(scope[2]).strip().lower() == "base" else str(scope[2]),
            )

    return by_chunk_id, combo_owner_task_id, combo_total_chunks, combo_total_ligands


def build_combo_weight_totals(
    *, by_chunk_id: Mapping[str, Mapping[str, Any]], label: str
) -> tuple[dict[ScopeKey, float], float]:
    combo_scope_total_weight: dict[ScopeKey, float] = {}
    for payload in by_chunk_id.values():
        scope = chunk_scope_key(
            payload,
            variant_label=(str(label or "").strip().upper() or "BASE"),
        )
        combo_scope_total_weight[scope] = combo_scope_total_weight.get(scope, 0.0) + float(
            payload.get("weight", 0.0) or 0.0
        )
    return combo_scope_total_weight, float(sum(combo_scope_total_weight.values()))


def collect_terminal_chunk_failures(
    *,
    by_chunk_id: Mapping[str, Mapping[str, Any]],
    label: str,
    dist_ctx: Any,
    max_attempts: int,
    read_chunk_result: Callable[..., Mapping[str, Any] | None],
    cfg_v: Mapping[str, Any],
    recorded_terminal_chunk_ids: set[str],
    failed_entries: list[tuple[str, str, str, str, str]],
) -> dict[ScopeKey, int]:
    combo_completed_chunks: dict[ScopeKey, int] = {}
    for chunk_id, chunk_payload in by_chunk_id.items():
        try:
            result_payload = read_chunk_result(dist_ctx, cfg_v, chunk_id=str(chunk_id))
        except Exception:
            logging.warning(
                "[distributed.chunk.result.read] action=skip chunk_id=%s reason=read_exception",
                str(chunk_id),
                exc_info=True,
            )
            continue
        if not isinstance(result_payload, dict):
            continue
        combo_scope = chunk_payload_scope(chunk_payload, variant_label=label)
        status = str(result_payload.get("status") or "").strip().lower()
        if status == "completed":
            combo_completed_chunks[combo_scope] = int(
                combo_completed_chunks.get(combo_scope, 0)
            ) + 1
        try:
            attempt = int(result_payload.get("attempt") or 1)
        except Exception:
            attempt = 1
        if status not in {"terminal_failed", "failed"}:
            continue
        if status == "failed" and int(attempt) < int(max_attempts):
            continue
        try:
            owner_task = int(result_payload.get("task_id") or -1)
        except Exception:
            owner_task = -1
        if owner_task != int(dist_ctx.task_id):
            continue
        if chunk_id in recorded_terminal_chunk_ids:
            continue
        recorded_terminal_chunk_ids.add(chunk_id)
        pdb_failed = (
            str(chunk_payload.get("pdb_id") or "").strip().upper()
            or os.path.splitext(os.path.basename(str(chunk_payload.get("pdb_file") or "")))[
                0
            ].upper()
        )
        ph_failed = str(chunk_payload.get("ph_tag") or "")
        failed_entries.append(
            (
                pdb_failed,
                label,
                "-",
                "ChunkTerminalFailed",
                f"chunk_id={chunk_id} ph={ph_failed} status={status} attempts={attempt}",
            )
        )
    return combo_completed_chunks


def collect_combo_chunk_progress(
    *,
    by_chunk_id: Mapping[str, Mapping[str, Any]],
    label: str,
    dist_ctx: Any,
    max_attempts: int,
    read_chunk_result: Callable[..., Mapping[str, Any] | None],
    cfg_v: Mapping[str, Any],
) -> dict[ScopeKey, tuple[int, int, int]]:
    combo_progress: dict[ScopeKey, list[int]] = {}
    for chunk_id, chunk_payload in by_chunk_id.items():
        try:
            result_payload = read_chunk_result(dist_ctx, cfg_v, chunk_id=str(chunk_id))
        except Exception:
            logging.warning(
                "[distributed.chunk.progress.read] action=skip chunk_id=%s reason=read_exception",
                str(chunk_id),
                exc_info=True,
            )
            continue
        if not isinstance(result_payload, dict):
            continue
        status = str(result_payload.get("status") or "").strip().lower()
        try:
            attempt = int(result_payload.get("attempt") or 1)
        except Exception:
            attempt = 1
        if not chunk_terminal(status, attempt, max_attempts=int(max_attempts)):
            continue
        combo_scope = chunk_payload_scope(chunk_payload, variant_label=label)
        counts = combo_progress.setdefault(combo_scope, [0, 0, 0])
        counts[0] += 1
        if status == "completed":
            counts[1] += 1
        else:
            counts[2] += 1
    return {
        scope: (int(counts[0]), int(counts[1]), int(counts[2]))
        for scope, counts in combo_progress.items()
    }


def aggregate_completion_markers_for_scope(
    *,
    combo_scope: ScopeKey,
    cfg_v: Mapping[str, Any],
    run_id: str,
    logger: logging.Logger | None = None,
) -> int:
    pdb_scope, variant_scope, ph_scope = combo_scope
    ph_value = None if str(ph_scope).strip().lower() == "base" else str(ph_scope)
    combo_docked_dir = _resolve_combo_output_dir(
        cfg_v,
        root_key="DOCKED_DIR",
        default_name="docked",
        run_id=str(run_id),
        pdb_id=str(pdb_scope).upper(),
        variant_label=str(variant_scope),
        ph_tag=ph_value,
    )
    return int(
        aggregate_combo_chunk_markers(
            combo_docked_dir,
            logger=logger,
        )
    )


def aggregate_completion_markers_for_owned_scopes(
    *,
    combo_owner_task_id: Mapping[ScopeKey, int],
    combo_total_chunks: Mapping[ScopeKey, int],
    dist_ctx: Any,
    cfg_v: Mapping[str, Any],
    run_id: str,
) -> int:
    owned_combo_scopes = [
        scope
        for scope, owner_task in sorted(combo_owner_task_id.items())
        if int(owner_task) == int(dist_ctx.task_id)
        and int(combo_total_chunks.get(scope, 0)) > 0
    ]
    completion_aggregate_count = 0
    completion_aggregate_logger = logging.getLogger("completion.aggregate")
    for combo_scope in owned_combo_scopes:
        try:
            completion_aggregate_count += int(
                aggregate_completion_markers_for_scope(
                    combo_scope=combo_scope,
                    cfg_v=cfg_v,
                    run_id=str(run_id),
                    logger=completion_aggregate_logger,
                )
            )
        except Exception:
            pdb_scope, variant_scope, ph_scope = combo_scope
            logging.warning(
                "[completion.aggregate] action=skip pdb=%s variant=%s ph=%s",
                str(pdb_scope),
                str(variant_scope),
                str(ph_scope),
                exc_info=True,
            )

    if completion_aggregate_count > 0:
        logging.info(
            "[completion.aggregate] run_id=%s task_id=%d combos_owned=%d files=%d",
            run_id,
            int(dist_ctx.task_id),
            len(owned_combo_scopes),
            int(completion_aggregate_count),
        )
    return int(completion_aggregate_count)
