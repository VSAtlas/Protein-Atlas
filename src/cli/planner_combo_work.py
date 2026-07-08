from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from cli.planner_hash import _hash_token
from config.normalize import _to_bool
from path_router.context_ph import select_ph_values_for_protonation
from path_router.path_router import load_ph_tags, make_paths, receptor_file


def _normalize_ph_tag_token(ph_tag: Optional[str]) -> str:
    token = str(ph_tag or "").strip()
    return token if token else "base"


def _format_ph_tag_from_value(value: float) -> str:
    rounded = max(0.0, min(14.0, round(float(value), 1)))
    return f"pH{rounded:.1f}".replace(".", "_")


def _discover_combo_ph_tags(
    cfg: Mapping[str, Any],
    *,
    pdb_file: str,
    pdb_id: str,
    variant_token: Optional[str],
) -> list[Optional[str]]:
    if not _to_bool(cfg.get("PH_ENSEMBLE", False)):
        return [None]

    tags = [
        str(tag).strip()
        for tag in (load_ph_tags(pdb_id, variant=variant_token) or [])
        if str(tag).strip()
    ]
    if not tags:
        pdb_path = Path(str(cfg.get("INPUT_DIR", "."))) / str(pdb_file)
        try:
            guessed_values = select_ph_values_for_protonation(str(pdb_path))
        except Exception:
            guessed_values = []
        tags = []
        for ph_value in guessed_values or []:
            try:
                tags.append(_format_ph_tag_from_value(float(ph_value)))
            except Exception:
                continue

    if not tags:
        return [None]

    seen: set[str] = set()
    ordered: list[Optional[str]] = []
    for raw in tags:
        token = str(raw).strip()
        if not token or token in seen:
            continue
        seen.add(token)
        ordered.append(token)
    return ordered or [None]


def _build_combo_work_items(
    cfg: Mapping[str, Any], pdb_files: list[str], *, variant_token: Optional[str]
) -> list[tuple[str, Optional[str]]]:
    items: list[tuple[str, Optional[str]]] = []
    for pdb_file in pdb_files:
        pdb_id = os.path.splitext(os.path.basename(str(pdb_file)))[0].upper()
        ph_tags = _discover_combo_ph_tags(
            cfg,
            pdb_file=str(pdb_file),
            pdb_id=pdb_id,
            variant_token=variant_token,
        )
        for ph_tag in ph_tags:
            items.append((str(pdb_file), ph_tag))
    return items


def _combo_pdb_id_from_file(pdb_file: str) -> str:
    return os.path.splitext(os.path.basename(str(pdb_file)))[0].upper()


def _combo_prep_key(
    *,
    run_id: str,
    variant_label: str,
    pdb_file: str,
    ph_tag: Optional[str],
) -> str:
    return _hash_token(
        "combo_prep",
        str(run_id),
        str(variant_label).upper(),
        _combo_pdb_id_from_file(pdb_file),
        _normalize_ph_tag_token(ph_tag),
    )[:24]


def _scope_prep_key(
    *,
    run_id: str,
    variant_label: str,
    pdb_file: str,
    ph_tag: Optional[str],
) -> str:
    return _hash_token(
        "scope_prep",
        str(run_id),
        str(variant_label).upper(),
        _combo_pdb_id_from_file(pdb_file),
        _normalize_ph_tag_token(ph_tag),
    )[:24]


def _local_combo_scope_key(*, pdb_file: str, variant_label: Optional[str]) -> str:
    variant_token = str(variant_label or "").strip().upper() or "BASE"
    return f"{_combo_pdb_id_from_file(pdb_file)}|{variant_token}"


def _pick_next_combo_index(
    pending: Sequence[tuple[str, Optional[str]]],
    *,
    inflight_pdb_ids: Optional[set[str]] = None,
    inflight_scope_counts: Optional[Mapping[str, int]] = None,
    prep_ready_scopes: Optional[set[str]] = None,
    allow_same_scope_when_ready: bool = False,
    variant_label: Optional[str] = None,
) -> Optional[int]:
    ready_scopes = prep_ready_scopes or set()
    for idx, (pdb_file, _ph_tag) in enumerate(pending):
        scope_key = _local_combo_scope_key(
            pdb_file=pdb_file,
            variant_label=variant_label,
        )
        if inflight_scope_counts is not None:
            inflight = int(inflight_scope_counts.get(scope_key, 0) or 0)
        elif inflight_pdb_ids is not None:
            inflight = (
                1 if _combo_pdb_id_from_file(pdb_file) in set(inflight_pdb_ids) else 0
            )
        else:
            inflight = 0
        if inflight <= 0:
            return idx
        if allow_same_scope_when_ready and scope_key in ready_scopes:
            return idx
        if inflight > 0:
            continue
    return None


def _combo_receptor_ready(
    cfg: Mapping[str, Any],
    *,
    pdb_file: str,
    variant_token: Optional[str],
    ph_tag: Optional[str],
    min_mtime: Optional[float] = None,
) -> bool:
    pdb_id = _combo_pdb_id_from_file(pdb_file)
    try:
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        ph_token = _normalize_ph_tag_token(ph_tag)
        ph_value = None if ph_token == "base" else str(ph_tag)
        receptor_path = receptor_file(
            paths.pdb_id,
            variant=None if variant_token is None else str(variant_token).upper(),
            ph_tag=ph_value,
            legacy=bool(cfg.get("_ROUTER_LEGACY", False)),
        )
        if not receptor_path.exists():
            return False
        stat = receptor_path.stat()
        if stat.st_size <= 0:
            return False
        if min_mtime is not None and float(stat.st_mtime) + 1e-6 < float(min_mtime):
            return False
        return True
    except Exception:
        return False


def _refresh_local_prep_ready_scopes(
    cfg: Mapping[str, Any],
    *,
    variant_token: Optional[str],
    variant_label: str,
    combos: Sequence[tuple[str, Optional[str]]],
    ready_scopes: set[str],
    scope_ready_after: Optional[Mapping[str, float]] = None,
    default_ready_after: Optional[float] = None,
) -> None:
    ready_after_lookup = scope_ready_after or {}
    for pdb_file, ph_tag in combos:
        scope_key = _local_combo_scope_key(
            pdb_file=pdb_file,
            variant_label=variant_label,
        )
        if scope_key in ready_scopes:
            continue
        min_mtime = ready_after_lookup.get(scope_key, default_ready_after)
        try:
            min_mtime = None if min_mtime is None else float(min_mtime)
        except Exception:
            min_mtime = default_ready_after
        if _combo_receptor_ready(
            cfg,
            pdb_file=pdb_file,
            variant_token=variant_token,
            ph_tag=ph_tag,
            min_mtime=min_mtime,
        ):
            ready_scopes.add(scope_key)
