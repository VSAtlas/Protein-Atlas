from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Tuple

from context_ph import select_ph_values_for_protonation
from ph_ensemble import build_ph_ensemble
from protein_functions import detect_active_site
from path_router import Paths


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
            raw_vals = select_ph_values_for_protonation(str(paths.input_pdb_path))
            logger.info("[ph.ctx.raw] path=%s values=%s", str(paths.input_pdb_path),
                        ",".join(f"{v:.2f}" for v in (raw_vals or [])))

            ph_values = sorted({max(3.0, min(10.5, round(float(x), 1))) for x in (raw_vals or [])})
            if not ph_values:
                logger.warning("[ph.ctx.fallback] context list empty -> using [7.0]")
                ph_values = [7.0]

            logger.info("[ph.list] n=%d values=%s", len(ph_values),
                        ",".join(f"{v:.1f}" for v in ph_values))

            ph_center = (0.0, 0.0, 0.0)
            ph_radius = 1_000_000.0
            if center is not None:
                ph_center = tuple(float(x) for x in center)
                logger.info(
                    "[ph.center] source=precomputed center=%s radius=%s",
                    ph_center,
                    ph_radius,
                )
            elif cleaned_pdb:
                c2, _b2, src = detect_active_site(cleaned_pdb)
                if c2:
                    ph_center = tuple(float(x) for x in c2)
                    logger.info(
                        "[ph.center] source=%s center=%s radius=%s",
                        src or "activesite",
                        ph_center,
                        ph_radius,
                    )
                else:
                    logger.warning(
                        "[ph.center] fallback=origin reason=detect_active_site_none radius=%s",
                        ph_radius,
                    )
            else:
                logger.warning(
                    "[ph.center] fallback=origin reason=missing_cleaned_pdb radius=%s",
                    ph_radius,
                )

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
