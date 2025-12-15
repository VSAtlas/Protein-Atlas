# record_data.py
# Helpers for recording docking scores, ligand efficiency, and CSV outputs.

from __future__ import annotations

import csv
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from path_router import make_paths
from input_and_export_functions import write_score_summary_to_csv
from pose_validation import compute_self_rmsd


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


def _pose_path_for(csv_cfg: Dict, pdb_id: str, stage_name: str, lig_path: str,
                   ph_label: Optional[str] = None, variant: Optional[str] = None) -> str:
    """Build the expected pose path for a ligand at a given stage."""
    from pathlib import Path
    # >>> DOCKED PATHS PATCH START
    paths = make_paths(csv_cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_token = (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    stage_dir = paths.docked_stage_dir(variant_token, stage_name, ph_token)
    return str(stage_dir / f"{Path(lig_path).stem}_{stage_name}.pdbqt")
    # >>> DOCKED PATHS PATCH END


def write_scores_csv(cfg: Dict, pdb_id: str, score_history: Dict[str, Dict[str, Dict]],
                     ph_label: Optional[str] = None, variant: Optional[str] = None, *,
                     csv_prefix: str = "",
                     summary_basename: str = "docking_score_summary.csv",
                     long_basename: str = "docking_score_long.csv") -> str:
    import csv, math

    # >>> DOCKED PATHS PATCH START
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_env = (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    dock_dir = paths.docked_variant_root(variant_token, ph_token)
    dock_dir.mkdir(parents=True, exist_ok=True)
    # >>> DOCKED PATHS PATCH END

    run_id_value = str(cfg.get("RUN_ID") or "")
    variant_value = variant_env
    include_variant = bool(variant_value)

    # --- Wide summary (unchanged shape) ---
    long_name = f"{csv_prefix}{long_basename}"
    summary_name = f"{csv_prefix}{summary_basename}"
    csv_out_wide = str(dock_dir / summary_name)
    flat = {}
    for stage_name, stage_map in score_history.items():
        flat[stage_name] = {}
        for lig, rec in stage_map.items():
            lig_key = os.path.basename(lig)
            s = rec.get("score", None)
            if rec.get("valid", False):
                flat[stage_name][lig_key] = s if s is not None else ""
            else:
                flat[stage_name][lig_key] = f"{s:.2f} (invalid)" if isinstance(s, (int, float)) else "(invalid)"
    write_score_summary_to_csv(
        flat,
        output_path=csv_out_wide,
        run_id=run_id_value,
        variant=variant_value if include_variant else None,
    )

    # --- Long format with self_rmsd added ---
    csv_out_long = str(dock_dir / long_name)
    with open(csv_out_long, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        header = ["run_id"]
        if include_variant:
            header.append("variant")
        header.extend(["stage", "ligand", "score", "valid", "reason", "heavy_atoms", "le", "self_rmsd", "pains_flag"])
        writer.writerow(header)

        for stage_name, stage_map in score_history.items():
            for lig, rec in stage_map.items():
                lig_key = os.path.basename(lig)
                score = rec.get("score", None)
                valid = bool(rec.get("valid", False))
                reason = rec.get("reason", "")

                ha = rec.get("heavy_atoms", None)
                le = rec.get("le", None)
                pains_hit = rec.get("pains_flag", False)

                # Compute self-RMSD from the saved pose for this stage (if present)
                pose_path = _pose_path_for(
                    cfg,
                    pdb_id,
                    stage_name,
                    lig,
                    ph_label=ph_token,
                    variant=variant_token,
                )
                if os.path.exists(pose_path):
                    try:
                        sr = compute_self_rmsd(pose_path)
                        sr_str = f"{sr:.2f}" if isinstance(sr, (int, float)) and math.isfinite(sr) else ""
                    except Exception:
                        sr_str = ""
                else:
                    sr_str = ""

                score_str = f"{score:.2f}" if isinstance(score, (int, float)) else ""
                ha_str = str(int(ha)) if isinstance(ha, (int, float)) else ""
                le_str = f"{le:.4f}" if isinstance(le, (int, float)) else ""
                reason_str = str(reason) if reason is not None else ""

                row = [run_id_value]
                if include_variant:
                    row.append(variant_value)
                row.extend([stage_name, lig_key, score_str, int(valid), reason_str, ha_str, le_str, sr_str, int(pains_hit)])
                writer.writerow(row)

    return csv_out_wide


def write_gnina_scores_csv(
    cfg: Dict,
    pdb_id: str,
    gnina_metrics_by_stage: Dict[str, Dict[str, Dict[str, Any]]],
    *,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
    csv_prefix: str = "",
) -> str:
    import csv, math

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_env = (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    dock_dir = paths.docked_variant_root(variant_token, ph_token)
    dock_dir.mkdir(parents=True, exist_ok=True)

    run_id_value = str(cfg.get("RUN_ID") or "")
    variant_value = variant_env
    include_variant = bool(variant_value)

    summary_name = "gnina_docking_score_summary.csv"
    long_name = f"{csv_prefix}gnina_docking_score_long.csv"
    csv_out_wide = str(dock_dir / summary_name)

    flat: Dict[str, Dict[str, Any]] = {}
    for stage_name, stage_map in gnina_metrics_by_stage.items():
        flat[stage_name] = {}
        for lig, rec in stage_map.items():
            lig_key = os.path.basename(lig)
            primary = rec.get("gnina_primary_score")
            if rec.get("valid"):
                flat[stage_name][lig_key] = f"{primary:.2f}" if isinstance(primary, (int, float)) else ""
            else:
                if isinstance(primary, (int, float)):
                    flat[stage_name][lig_key] = f"{primary:.2f} (invalid)"
                else:
                    flat[stage_name][lig_key] = "(invalid)"

    write_score_summary_to_csv(
        flat,
        output_path=csv_out_wide,
        run_id=run_id_value,
        variant=variant_value if include_variant else None,
    )

    csv_out_long = str(dock_dir / long_name)
    with open(csv_out_long, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        header = ["run_id"]
        if include_variant:
            header.append("variant")
        header.extend(
            [
                "stage",
                "ligand",
                "valid",
                "reason",
                "heavy_atoms",
                "le",
                "self_rmsd",
                "pains_flag",
                "gnina_minimized_affinity_kcal",
                "gnina_cnn_score",
                "gnina_cnn_affinity_pK",
                "gnina_primary_score",
            ]
        )
        writer.writerow(header)

        for stage_name, stage_map in gnina_metrics_by_stage.items():
            for lig, rec in stage_map.items():
                lig_key = os.path.basename(lig)
                valid = bool(rec.get("valid", False))
                reason = rec.get("reason", "")
                ha = rec.get("heavy_atoms", None)
                le = rec.get("le", None)
                sr = rec.get("self_rmsd", None)
                pains_hit = rec.get("pains_flag", False)
                minimized_affinity = rec.get("minimized_affinity_kcal")
                cnn_score = rec.get("cnn_score")
                cnn_affinity = rec.get("cnn_affinity_pK")
                primary = rec.get("gnina_primary_score")

                ha_str = str(int(ha)) if isinstance(ha, (int, float)) else ""
                le_str = f"{le:.4f}" if isinstance(le, (int, float)) else ""
                sr_str = f"{sr:.2f}" if isinstance(sr, (int, float)) and math.isfinite(sr) else ""
                reason_str = str(reason) if reason is not None else ""

                def _fmt(val: Any, places: int = 2) -> str:
                    return f"{val:.{places}f}" if isinstance(val, (int, float)) else ""

                row = [run_id_value]
                if include_variant:
                    row.append(variant_value)
                row.extend(
                    [
                        stage_name,
                        lig_key,
                        int(valid),
                        reason_str,
                        ha_str,
                        le_str,
                        sr_str,
                        int(bool(pains_hit)),
                        _fmt(minimized_affinity),
                        _fmt(cnn_score, places=3),
                        _fmt(cnn_affinity, places=3),
                        _fmt(primary, places=3),
                    ]
                )
                writer.writerow(row)

    return csv_out_wide
