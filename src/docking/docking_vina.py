from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional

from path_router.path_router import (
    make_paths,
    config_dir as router_config_dir,
    config_file as router_config_file,
    docked_dir as router_docked_dir,
    receptor_file as router_receptor_file,
)
from input_and_export_functions import write_score_summary_to_csv
from docking.pose_validation import compute_self_rmsd

# Engine-specific home for Vina config emission and CSV writing.


def emit_vina_config(
    cfg: Dict[str, Any],
    pdb_id: str,
    receptor_pdbqt: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    ligand_path: str,
    stage_name: str,
    stage_info: Dict[str, Any],
    cpu_per_job: int,
    logger: Optional[Any] = None,
    *,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    legacy: bool = False,
):
    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )
    ph_label = (str(ph_token).strip() or None) if ph_token is not None else None
    legacy_mode = bool(legacy)

    lig_base = Path(ligand_path).stem
    run_id = cfg["RUN_ID"]

    cfg_dir = router_config_dir(
        run_id,
        pdb_id,
        stage_name,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    cfg_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = router_config_file(
        run_id,
        pdb_id,
        stage_name,
        variant=variant_token,
        ph_tag=ph_label,
        name="vina.json",
        legacy=legacy_mode,
    )

    cfg_path = cfg_dir / f"{lig_base}_{stage_name}.txt"

    stage_root = router_docked_dir(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    out_dir = stage_root / stage_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{lig_base}_{stage_name}.pdbqt"

    expected_receptor = router_receptor_file(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    receptor_exists = expected_receptor.exists()
    receptor_for_config = str(expected_receptor)

    variant_display = variant_token or "None"
    ph_display = ph_label or "None"
    breadcrumb = (
        "[cfg.emit] run=%s pdb=%s stage=%s variant=%s ph=%s\n"
        "           cfg_dir=%s receptor=%s out_root=%s"
    )
    breadcrumb_args = (
        run_id,
        pdb_id,
        stage_name,
        variant_display,
        ph_display,
        str(cfg_dir),
        receptor_for_config,
        str(stage_root),
    )
    if logger:
        logger.info(breadcrumb, *breadcrumb_args)
    else:
        print(breadcrumb % breadcrumb_args)

    if not receptor_exists:
        msg = (
            f"[router.error] missing receptor for pdb={pdb_id} variant={variant_display} "
            f"ph={ph_display} -> {expected_receptor}"
        )
        if logger:
            logger.error(msg)
        else:
            print(msg)

    exhaustiveness = int(stage_info.get("exhaustiveness", 8))
    num_modes = int(stage_info.get("num_modes", 4))
    if cfg.get("FAST_MODE"):
        exhaustiveness = 1
        num_modes = 1

    lines = [
        f"receptor = {receptor_for_config}",
        f"ligand   = {ligand_path}",
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x   = {box_size[0]:.3f}",
        f"size_y   = {box_size[1]:.3f}",
        f"size_z   = {box_size[2]:.3f}",
        f"cpu      = {int(cpu_per_job)}",
        f"exhaustiveness = {exhaustiveness}",
        f"energy_range   = {int(stage_info.get('energy_range', 4))}",
        f"num_modes      = {num_modes}",
        f"verbosity      = {int(stage_info.get('verbosity', 0))}",
        f"out = {out_path}",
    ]

    if "seed" in stage_info:
        lines.append(f"seed = {int(stage_info['seed'])}")
    if logger:
        cx, cy, cz = center
        sx, sy, sz = box_size
        lig_name = os.path.basename(str(ligand_path))
        logger.info(
            "[vina.cfg] lig=%s center=(%.3f,%.3f,%.3f) size=(%.1f,%.1f,%.1f)",
            lig_name,
            cx,
            cy,
            cz,
            sx,
            sy,
            sz,
        )

    payload = ("\n".join(lines)).encode("utf-8")
    overwrite = cfg_path.exists()

    tmp = cfg_path.with_suffix(".part")
    with open(tmp, "wb") as f:
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, cfg_path)

    manifest_data: Dict[str, Any]
    entries_map: Dict[str, Dict[str, Any]]
    if manifest_path.exists():
        try:
            manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            manifest_data = {}
    else:
        manifest_data = {}

    entries = manifest_data.get("entries") if isinstance(manifest_data, dict) else None
    entries_map = {}
    if isinstance(entries, list):
        for item in entries:
            if isinstance(item, dict):
                lig = str(item.get("ligand", ""))
                if lig:
                    entries_map[lig] = item

    entry = {
        "ligand": lig_base,
        "config": str(cfg_path),
        "out": str(out_path),
        "receptor": receptor_for_config,
    }
    entries_map[lig_base] = entry

    manifest_data = {
        "run_id": run_id,
        "pdb_id": pdb_id,
        "stage": stage_name,
        "variant": variant_token,
        "ph": ph_label,
        "legacy": legacy_mode,
        "entries": [entries_map[k] for k in sorted(entries_map.keys())],
    }

    manifest_tmp = manifest_path.with_suffix(".part")
    try:
        manifest_tmp.parent.mkdir(parents=True, exist_ok=True)
        with open(manifest_tmp, "w", encoding="utf-8") as fh:
            json.dump(manifest_data, fh, indent=2, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(manifest_tmp, manifest_path)
    except FileNotFoundError:
        if logger:
            logger.warning(
                "[cfg.emit] manifest_tmp missing; skipping manifest update tmp=%s dest=%s stage=%s",
                str(manifest_tmp),
                str(manifest_path),
                stage_name,
            )
    except Exception:
        if logger:
            logger.exception(
                "[cfg.emit] manifest update failed; continuing without manifest stage=%s path=%s",
                stage_name,
                str(manifest_path),
            )

    emit_msg = (
        "[cfg.emit] run=%s pdb=%s variant=%s ph=%s stage=%s ligand=%s "
        "cfg_dir=%s docked_root=%s path=%s overwrite=%s bytes=%d"
    )
    emit_args = (
        run_id,
        pdb_id,
        variant_display,
        ph_display,
        stage_name,
        lig_base,
        str(cfg_dir),
        str(stage_root),
        str(cfg_path),
        str(overwrite).lower(),
        len(payload),
    )
    if logger:
        logger.info(emit_msg, *emit_args)
    else:
        print(emit_msg % emit_args)

    return str(cfg_path), str(out_path)


def _pose_path_for(
    csv_cfg: Dict,
    pdb_id: str,
    stage_name: str,
    lig_path: str,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
) -> str:
    """Build the expected pose path for a ligand at a given stage."""
    paths = make_paths(csv_cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_token = (
        variant or os.environ.get("APO_HOLO_VARIANT", "") or ""
    ).strip().upper() or None
    stage_dir = paths.docked_stage_dir(variant_token, stage_name, ph_token)
    return str(stage_dir / f"{Path(lig_path).stem}_{stage_name}.pdbqt")


def write_scores_csv(
    cfg: Dict,
    pdb_id: str,
    score_history: Dict[str, Dict[str, Dict]],
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
    *,
    csv_prefix: str = "",
    summary_basename: str = "docking_score_summary.csv",
    long_basename: str = "docking_score_long.csv",
) -> str:
    import csv as _csv
    import math as _math

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_env = (
        (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    )
    variant_token = variant_env or None
    dock_dir = paths.docked_variant_root(variant_token, ph_token)
    dock_dir.mkdir(parents=True, exist_ok=True)

    run_id_value = str(cfg.get("RUN_ID") or "")
    variant_value = variant_env
    include_variant = bool(variant_value)

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
                flat[stage_name][lig_key] = (
                    f"{s:.2f} (invalid)" if isinstance(s, (int, float)) else "(invalid)"
                )
    write_score_summary_to_csv(
        flat,
        output_path=csv_out_wide,
        run_id=run_id_value,
        variant=variant_value if include_variant else None,
    )

    csv_out_long = str(dock_dir / long_name)
    with open(csv_out_long, "w", newline="", encoding="utf-8") as f:
        writer = _csv.writer(f)
        header = ["run_id"]
        if include_variant:
            header.append("variant")
        header.extend(
            [
                "stage",
                "ligand",
                "score",
                "valid",
                "reason",
                "heavy_atoms",
                "le",
                "self_rmsd",
                "pains_flag",
            ]
        )
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
                        sr_str = (
                            f"{sr:.2f}"
                            if isinstance(sr, (int, float)) and _math.isfinite(sr)
                            else ""
                        )
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
                row.extend(
                    [
                        stage_name,
                        lig_key,
                        score_str,
                        int(valid),
                        reason_str,
                        ha_str,
                        le_str,
                        sr_str,
                        int(pains_hit),
                    ]
                )
                writer.writerow(row)

    return csv_out_wide
