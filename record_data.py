# record_data.py
# Helpers for recording docking scores, ligand efficiency, and CSV outputs.

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from path_router import make_paths


def compute_ligand_efficiency(score: Optional[float], heavy_atoms: Optional[int]) -> Optional[float]:
    try:
        if score is None or heavy_atoms is None or int(heavy_atoms) <= 0:
            return None
        return float(-float(score) / int(heavy_atoms))
    except Exception:
        return None


def record_le(score_history: Dict[str, Dict[str, Dict]],
              stage_name: str,
              lig_path: str,
              score: Optional[float],
              heavy_atom_counts: Dict[str, int]) -> Optional[float]:
    ha = heavy_atom_counts.get(lig_path)
    le = compute_ligand_efficiency(score, ha)
    stage_map = score_history.setdefault(stage_name, {})
    rec = stage_map.setdefault(lig_path, {})
    rec["heavy_atoms"] = int(ha) if isinstance(ha, (int, float)) else None
    rec["le"] = le
    return le


def write_gnina_scores_csv(*args, **kwargs):
    from docking.docking_gnina import write_gnina_scores_csv as _impl
    return _impl(*args, **kwargs)


def _pose_path_for(csv_cfg: Dict, pdb_id: str, stage_name: str, lig_path: str,
                   ph_label: Optional[str] = None, variant: Optional[str] = None) -> str:
    """Build the expected pose path for a ligand at a given stage."""
    from pathlib import Path
    paths = make_paths(csv_cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_token = (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    stage_dir = paths.docked_stage_dir(variant_token, stage_name, ph_token)
    return str(stage_dir / f"{Path(lig_path).stem}_{stage_name}.pdbqt")


def write_scores_csv(cfg: Dict, pdb_id: str, score_history: Dict[str, Dict[str, Dict]],
                     ph_label: Optional[str] = None, variant: Optional[str] = None, *,
                     csv_prefix: str = "",
                     summary_basename: str = "docking_score_summary.csv",
                     long_basename: str = "docking_score_long.csv") -> str:
    from docking.docking_vina import write_scores_csv as _impl
    return _impl(
        cfg,
        pdb_id,
        score_history,
        ph_label=ph_label,
        variant=variant,
        csv_prefix=csv_prefix,
        summary_basename=summary_basename,
        long_basename=long_basename,
    )
