from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, Optional, Tuple

from path_router.context_ph import select_ph_values_for_protonation
from docking.ph_ensemble_runtime import build_ph_ensemble
from path_router.path_router import Paths


def _parse_ph_override_values(label: str) -> list[float]:
    values: list[float] = []
    token = str(label or "").strip()
    if not token:
        return values
    core = token
    m = re.search(r"(?i)pH(.+)", token)
    if m:
        core = m.group(1)
    if "-" in core:
        core = core.split("-", 1)[0]
    for part in core.split("+"):
        part = part.strip().replace("_", ".")
        if not part:
            continue
        try:
            ph_value = float(part)
        except Exception:
            continue
        if 0.0 < ph_value < 15.0:
            values.append(ph_value)
    return values


def _phase5_ph_ensemble_global(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    variant_env: str,
    variant_token: Optional[str],
    legacy_mode: bool,
    cleaned_pdb: Optional[str],
    receptor_pdbqt: Optional[str],
    center: Optional[Tuple[float, float, float]],
    box_size: Optional[Tuple[float, float, float]],
) -> Tuple[Optional[str], Optional[str]]:
    active_ph_label = None
    if bool(cfg.get("PH_ENSEMBLE", False)):
        try:
            ph_override = str(cfg.get("_PH_TAG_OVERRIDE", "") or "").strip()
            if ph_override:
                raw_vals = _parse_ph_override_values(ph_override)
                logger.info(
                    "[ph.ctx.override] pdb=%s variant=%s ph_override=%s values=%s",
                    paths.pdb_id,
                    variant_token or "legacy",
                    ph_override,
                    ",".join(f"{v:.2f}" for v in (raw_vals or [])),
                )
            else:
                raw_vals = select_ph_values_for_protonation(str(paths.input_pdb_path))
                logger.info(
                    "[ph.ctx.raw] path=%s values=%s",
                    str(paths.input_pdb_path),
                    ",".join(f"{v:.2f}" for v in (raw_vals or [])),
                )

            ph_values = sorted(
                {max(3.0, min(10.5, round(float(x), 1))) for x in (raw_vals or [])}
            )
            if not ph_values:
                logger.warning("[ph.ctx.fallback] context list empty -> using [7.0]")
                ph_values = [7.0]

            logger.info(
                "[ph.list] n=%d values=%s",
                len(ph_values),
                ",".join(f"{v:.1f}" for v in ph_values),
            )

            ph_center = (0.0, 0.0, 0.0)
            ph_radius = 1_000_000.0
            if center is not None:
                ph_center = (
                    float(center[0]),
                    float(center[1]),
                    float(center[2]),
                )
                logger.info(
                    "[ph.center] source=precomputed center=%s radius=%s",
                    ph_center,
                    ph_radius,
                )
            else:
                logger.info(
                    "[ph.center] source=global_default center=%s radius=%s",
                    ph_center,
                    ph_radius,
                )

            if not cleaned_pdb:
                raise ValueError("missing_cleaned_receptor")
            manifest_path = build_ph_ensemble(
                pdb_id=paths.pdb_id,
                cleaned_receptor_pdb=str(Path(cleaned_pdb)),
                out_dir=str(Path(cfg["OUTPUT_DIR"])),
                center=ph_center,
                radius=ph_radius,
                ph_values=ph_values,
                member_index_start=0,
                variant=variant_token,
                legacy=legacy_mode,
            )
            logger.info("[ph_ensemble.manifest] path=%s", manifest_path)

        except Exception as e:
            logger.warning("[ph_ensemble.skip] error=%s", e)
    return receptor_pdbqt, active_ph_label
