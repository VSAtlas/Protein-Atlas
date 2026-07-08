"""State, cache, and active-site helpers for docking runtime."""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

from cli.distributed_context import write_json_atomic
from path_router.path_router import make_paths
from docking.active_site_detection import detect_active_site

_ACTIVE_SITE_CACHE_KEY = "_ACTIVE_SITE_CACHE"
_COMBO_PREP_CENTER_KEY = "_COMBO_PREP_CENTER_BY_PH"
_COMBO_PREP_BOX_KEY = "_COMBO_PREP_BOX_BY_PH"
_COMBO_PREP_SOURCE_KEY = "_COMBO_PREP_SOURCE_BY_PH"
_COMBO_PREP_CONTROL_STEMS_KEY = "_COMBO_PREP_CONTROL_STEMS"
_COMBO_PREP_CONTROL_LOOKUP_KEY = "_COMBO_PREP_CONTROL_LOOKUP"

def _decode_ph_key(raw: Any) -> Optional[str]:
    token = str(raw or "").strip()
    return None if token.lower() in {"", "none", "base"} else token

def _encode_ph_key(ph_key: Optional[str]) -> str:
    token = str(ph_key or "").strip()
    return token if token else "base"

def _decode_vec3_map(raw: Any) -> Dict[Optional[str], Tuple[float, float, float]]:
    if not isinstance(raw, dict):
        return {}
    out: Dict[Optional[str], Tuple[float, float, float]] = {}
    for raw_key, raw_vec in raw.items():
        if not isinstance(raw_vec, (list, tuple)) or len(raw_vec) != 3:
            continue
        try:
            out[_decode_ph_key(raw_key)] = (
                float(raw_vec[0]),
                float(raw_vec[1]),
                float(raw_vec[2]),
            )
        except Exception:
            continue
    return out

def _encode_vec3_map(
    raw: Dict[Optional[str], Tuple[float, float, float]],
) -> Dict[str, List[float]]:
    out: Dict[str, List[float]] = {}
    for ph_key, vec in (raw or {}).items():
        if not isinstance(vec, (list, tuple)) or len(vec) != 3:
            continue
        try:
            out[_encode_ph_key(ph_key)] = [float(vec[0]), float(vec[1]), float(vec[2])]
        except Exception:
            continue
    return out

def _decode_source_map(raw: Any) -> Dict[Optional[str], str]:
    if not isinstance(raw, dict):
        return {}
    out: Dict[Optional[str], str] = {}
    for raw_key, raw_val in raw.items():
        token = str(raw_val or "").strip()
        if not token:
            continue
        out[_decode_ph_key(raw_key)] = token
    return out

def _encode_source_map(raw: Dict[Optional[str], str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for ph_key, raw_val in (raw or {}).items():
        token = str(raw_val or "").strip()
        if not token:
            continue
        out[_encode_ph_key(ph_key)] = token
    return out

def _decode_control_lookup(raw: Any) -> Dict[str, Path]:
    if not isinstance(raw, dict):
        return {}
    out: Dict[str, Path] = {}
    for raw_key, raw_val in raw.items():
        key = str(raw_key or "").strip()
        path_token = str(raw_val or "").strip()
        if not key or not path_token:
            continue
        out[key] = Path(path_token)
    return out

def _signal_distributed_prep_ready(
    cfg: Dict,
    logger,
    *,
    pdb_id: str,
    variant_label: str,
    active_ph_label,
) -> None:
    signal = cfg.get("_DISTRIBUTED_PREP_READY_SIGNAL")
    if not isinstance(signal, Mapping):
        return
    state_path_raw = str(signal.get("state_path") or "").strip()
    if not state_path_raw:
        return
    center_map = cfg.get(_COMBO_PREP_CENTER_KEY)
    box_map = cfg.get(_COMBO_PREP_BOX_KEY)
    if not isinstance(center_map, dict) or not center_map:
        logger.info(
            "[prep.ready.signal] action=skip reason=missing_center_map pdb=%s variant=%s ph=%s",
            pdb_id,
            variant_label,
            _encode_ph_key(active_ph_label),
        )
        return
    if not isinstance(box_map, dict) or not box_map:
        logger.info(
            "[prep.ready.signal] action=skip reason=missing_box_map pdb=%s variant=%s ph=%s",
            pdb_id,
            variant_label,
            _encode_ph_key(active_ph_label),
        )
        return

    def _int_signal(key: str, default: int) -> int:
        try:
            return int(signal.get(key, default) or default)
        except Exception:
            return int(default)

    payload: dict[str, Any] = {
        "run_id": str(signal.get("run_id") or cfg.get("RUN_ID", "") or ""),
        "combo_key": str(signal.get("combo_key") or ""),
        "task_id": _int_signal("task_id", 0),
        "task_count": _int_signal("task_count", 1),
        "status": "ready",
        "updated_at": float(time.time()),
        "owner_task_id": _int_signal("owner_task_id", _int_signal("task_id", 0)),
        "pdb_id": str(signal.get("pdb_id") or pdb_id).upper(),
        "variant": str(signal.get("variant_label") or variant_label).upper(),
        "ph": str(signal.get("ph_tag") or _encode_ph_key(active_ph_label)),
        "note": "early_prep_ready",
        "center_by_ph": dict(center_map),
        "box_by_ph": dict(box_map),
    }
    source_map = cfg.get(_COMBO_PREP_SOURCE_KEY)
    if isinstance(source_map, dict) and source_map:
        payload["center_source_by_ph"] = dict(source_map)
    control_stems = cfg.get(_COMBO_PREP_CONTROL_STEMS_KEY)
    if isinstance(control_stems, list) and control_stems:
        payload["control_stems"] = [str(x) for x in control_stems if str(x).strip()]
    control_lookup = cfg.get(_COMBO_PREP_CONTROL_LOOKUP_KEY)
    if isinstance(control_lookup, dict) and control_lookup:
        payload["control_lookup"] = {
            str(k): str(v)
            for k, v in control_lookup.items()
            if str(k).strip() and str(v).strip()
        }

    try:
        write_json_atomic(Path(state_path_raw), payload)
        claim_path_raw = str(signal.get("claim_path") or "").strip()
        if claim_path_raw:
            try:
                Path(claim_path_raw).unlink()
            except FileNotFoundError:
                pass
        logger.info(
            "[prep.ready.signal] action=write pdb=%s variant=%s ph=%s centers=%d state=%s",
            payload["pdb_id"],
            payload["variant"],
            payload["ph"],
            len(center_map),
            state_path_raw,
        )
    except Exception:
        logger.warning(
            "[prep.ready.signal] action=skip reason=write_failed pdb=%s variant=%s ph=%s state=%s",
            pdb_id,
            variant_label,
            _encode_ph_key(active_ph_label),
            state_path_raw,
            exc_info=True,
        )

def _cache_active_site(
    cfg: Dict[str, Any],
    pdb_id: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    logger: Optional[logging.Logger] = None,
) -> None:
    """Record active-site center/box for reuse by downstream helpers (e.g., DOCK6 prep)."""
    try:
        cache = cfg.setdefault(_ACTIVE_SITE_CACHE_KEY, {})
        cache_key = (
            str(pdb_id).upper(),
            (variant_token or "HOLO"),
            (ph_label or "base"),
        )
        cache[cache_key] = (
            tuple(float(x) for x in center),
            tuple(float(x) for x in box_size),
        )
        if logger:
            logger.debug(
                "[active-site.cache.store] key=%s center=%s box=%s",
                cache_key,
                cache[cache_key][0],
                cache[cache_key][1],
            )
    except Exception:
        if logger:
            logger.debug("[active-site.cache.skip]")

def get_active_site_center_and_size(
    cfg: Dict[str, Any],
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    logger: logging.Logger,
) -> Optional[tuple[tuple[float, float, float], tuple[float, float, float]]]:
    """
    Returns (center, size) for the active site used by Vina/LeDock.

    center: (cx, cy, cz)
    size:   (sx, sy, sz) full box lengths (same values emitted to Vina configs)

    Prefers cached values populated during docking setup; falls back to running
    detect_active_site on the cleaned receptor.
    """
    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )
    variant_key = variant_token or "HOLO"
    ph_key = (str(ph_label).strip() or "") or "base"

    cache = cfg.get(_ACTIVE_SITE_CACHE_KEY)
    if isinstance(cache, dict):
        cache_key = (str(pdb_id).upper(), variant_key, ph_key)
        hit = cache.get(cache_key)
        if not hit and ph_key != "base":
            hit = cache.get((str(pdb_id).upper(), variant_key, "base"))
        if hit:
            logger.info(
                "[active-site.cache.hit] pdb=%s variant=%s ph=%s center=%s box=%s",
                pdb_id,
                variant_key,
                ph_key,
                hit[0],
                hit[1],
            )
            return hit  # type: ignore[return-value]

    try:
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        cleaned_pdb = paths.receptor_cleaned_pdb(variant_token)
    except Exception as exc:
        logger.warning(
            "[active-site.helper.error] pdb=%s variant=%s ph=%s reason=%s",
            pdb_id,
            variant_key,
            ph_key,
            exc,
        )
        return None

    if not cleaned_pdb or not Path(cleaned_pdb).exists():
        logger.warning(
            "[active-site.helper.skip] reason=missing_cleaned pdb=%s variant=%s ph=%s path=%s",
            pdb_id,
            variant_key,
            ph_key,
            cleaned_pdb,
        )
        return None

    center, box_size, source = detect_active_site(cleaned_pdb)
    if center and box_size:
        center_t = (float(center[0]), float(center[1]), float(center[2]))
        box_t = (float(box_size[0]), float(box_size[1]), float(box_size[2]))
        _cache_active_site(cfg, pdb_id, variant_token, ph_key, center_t, box_t, logger)
        logger.info(
            "[active-site.helper.detect] pdb=%s variant=%s ph=%s source=%s center=%s box=%s",
            pdb_id,
            variant_key,
            ph_key,
            source or "activesite",
            center_t,
            box_t,
        )
        return center_t, box_t

    logger.warning(
        "[active-site.helper.miss] pdb=%s variant=%s ph=%s source=%s",
        pdb_id,
        variant_key,
        ph_key,
        source if "source" in locals() else "unknown",
    )
    return None
