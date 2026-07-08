# checkpoints.py
# Stage fingerprinting + checkpoint file helpers for docking runs.

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from path_router import docked_dir


def _checkpoint_path(
    cfg: Dict,
    pdb_id: str,
    stage_name: str,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
) -> Path:
    variant_token = (
        variant or os.environ.get("APO_HOLO_VARIANT", "") or ""
    ).strip().upper() or None
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    root = docked_dir(
        pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode
    )
    return root / f".ckpt_{stage_name}.json"


def checkpoint_should_skip(
    cfg: Dict,
    pdb_id: str,
    stage_name: str,
    fingerprint: Dict[str, Any],
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
    engine: Optional[str] = None,
) -> bool:
    p = _checkpoint_path(cfg, pdb_id, stage_name, ph_label=ph_label, variant=variant)
    if not p.exists():
        return False
    try:
        prev = json.loads(p.read_text())
    except Exception:
        return False
    if prev != fingerprint:
        return False

    if engine:
        legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
        variant_token = (
            variant or os.environ.get("APO_HOLO_VARIANT", "") or ""
        ).strip().upper() or None
        stage_dir = (
            docked_dir(
                pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode
            )
            / stage_name
        )
        marker = stage_dir / f"completion_{engine}.json"
        if not marker.exists():
            return False
        try:
            payload = json.loads(marker.read_text())
            if not payload or payload.get("success") is False:
                return False
            missing_after = payload.get("missing_count_after")
            if isinstance(missing_after, int) and missing_after > 0:
                return False
        except Exception:
            return False
    return True


def checkpoint_mark_done(
    cfg: Dict,
    pdb_id: str,
    stage_name: str,
    fingerprint: Dict[str, Any],
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
) -> None:
    p = _checkpoint_path(cfg, pdb_id, stage_name, ph_label=ph_label, variant=variant)
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        p.write_text(json.dumps(fingerprint, indent=2))
    except Exception:
        pass


def checkpoint_invalidate_from(
    cfg: Dict,
    pdb_id: str,
    stages: List[Dict],
    start_index: int,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
) -> None:
    for j in range(start_index, len(stages)):
        names = [stages[j]["name"]]
        try:
            base = stages[j]["name"]
            names.append(f"gnina_{base}")
        except Exception:
            pass
        for nm in names:
            try:
                _checkpoint_path(
                    cfg,
                    pdb_id,
                    nm,
                    ph_label=ph_label,
                    variant=variant,
                ).unlink(missing_ok=True)
            except Exception:
                pass
