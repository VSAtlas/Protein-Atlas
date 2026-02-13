# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import yaml  # type: ignore[import-untyped]

LOGGER = logging.getLogger("manifest-utils")
CANONICAL_MANIFEST_NAME = "run_manifest.yaml"


def _as_float(val: Any) -> Optional[float]:
    try:
        f_val = float(val)
        if math.isfinite(f_val):
            return f_val
    except Exception:
        return None
    return None


def _clean_vector(vec: Any) -> Optional[list[float]]:
    if isinstance(vec, (list, tuple)) and len(vec) == 3:
        cleaned: list[float] = []
        for item in vec:
            f_val = _as_float(item)
            if f_val is None:
                return None
            cleaned.append(f_val)
        return cleaned
    return None


def _extract_details(entry: Dict[str, Any]) -> Dict[str, Optional[Any]]:
    details = entry.get("stages", {}).get("pocket_detection", {}).get("details", {})
    method = details.get("method") or details.get("pocket_method")
    center = _clean_vector(details.get("center"))
    if center is None:
        center = _clean_vector(
            [details.get("center_x"), details.get("center_y"), details.get("center_z")]
        )
    box = _clean_vector(details.get("box_size"))
    if box is None:
        box = _clean_vector(
            [details.get("box_x"), details.get("box_y"), details.get("box_z")]
        )

    if method is None and center is None and box is None:
        return {}

    return {
        "method": method,
        "center": center,
        "box": box,
    }


def load_run_manifest(
    repo_root: Path, run_id: str
) -> Tuple[Optional[Dict[str, Any]], Optional[Path]]:
    """
    Load the run manifest for a given run_id.

    Returns (manifest_dict, path_used) or (None, None) if not found/readable.
    """
    canonical = repo_root / "manifests" / run_id / CANONICAL_MANIFEST_NAME
    candidates = [
        canonical,
        repo_root / "post_docked" / run_id / CANONICAL_MANIFEST_NAME,
        repo_root / "docked" / run_id / CANONICAL_MANIFEST_NAME,
        repo_root / CANONICAL_MANIFEST_NAME,
    ]

    for path in candidates:
        if not path.exists() or path.stat().st_size <= 0:
            continue
        try:
            with path.open("r", encoding="utf-8") as fh:
                return yaml.safe_load(fh) or {}, path
        except Exception:
            LOGGER.debug(
                "[manifest.load] action=read_failed path=%s", path, exc_info=True
            )
    return None, None


def extract_pocket(
    manifest: Optional[Dict[str, Any]], pdb_id: str, variant: str, ph_label: str
) -> Dict[str, Optional[Any]]:
    """
    Extract pocket metadata for a (pdb_id, variant, ph_label) triple.

    Supports composite protein keys like "<pdb>|<variant>|<ph>" and falls back
    to matching entries with pdb_id/variant/ph fields if the composite key is missing.
    """
    result = {"method": None, "center": None, "box": None}
    if not manifest:
        return result

    proteins = manifest.get("proteins", {}) or {}
    if not isinstance(proteins, dict):
        return result

    pdb_norm = (pdb_id or "").upper()
    variant_norm = (variant or "").upper()
    ph_norm = ph_label or ""
    composite_key = f"{pdb_norm}|{variant_norm}|{ph_norm}"

    def _match_entry(entry: Dict[str, Any]) -> Dict[str, Optional[Any]]:
        info = _extract_details(entry)
        if info:
            return info

        variants = entry.get("variants", {})
        if isinstance(variants, dict) and variant_norm:
            var_entry = variants.get(variant_norm) or variants.get(variant)
            if isinstance(var_entry, dict):
                info = _extract_details(var_entry)
                if info:
                    return info
        return {}

    # 1) Exact composite key
    direct_entry = proteins.get(composite_key) or proteins.get(
        f"{pdb_id}|{variant}|{ph_label}"
    )
    if isinstance(direct_entry, dict):
        info = _match_entry(direct_entry)
        if info:
            return info

    # 2) Keys that decompose to the desired combo
    for prot_key, entry in proteins.items():
        if not isinstance(entry, dict):
            continue
        parts = str(prot_key).split("|")
        if len(parts) >= 3:
            key_pdb, key_variant, key_ph = parts[0].upper(), parts[1].upper(), parts[2]
            if (
                key_pdb == pdb_norm
                and key_variant == variant_norm
                and key_ph == ph_norm
            ):
                info = _match_entry(entry)
                if info:
                    return info

    # 3) Entries whose fields declare matching pdb/variant/ph
    for entry in proteins.values():
        if not isinstance(entry, dict):
            continue
        field_pdb = (entry.get("pdb_id") or "").upper()
        field_variant = (entry.get("variant") or "").upper()
        field_ph = entry.get("ph") or entry.get("ph_label") or ""
        if (
            field_pdb == pdb_norm
            and field_variant == variant_norm
            and field_ph == ph_norm
        ):
            info = _match_entry(entry)
            if info:
                return info

    return result
