# -*- coding: utf-8 -*-
import csv
import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from analysis.reporting.fda_name_map import resolve_ligand_display_name
from analysis.reporting.ligand_side_effect_cache import lookup_ligand_side_effect_entry
from analysis.reporting.side_effect_overlap_meta import ligand_side_effect_report_fields
from analysis.reporting.value_utils import as_float
from analysis.reporting.master_schema.constants import COMPONENT, _SCHEMA_CSV_ERRORS
from analysis.reporting.master_schema.decoy_helpers import (
    _canonical_ligand_base_with_prefix,
    _resolve_library_name,
    _row_is_decoy,
)
from analysis.reporting.master_schema.metadata_loaders import (
    _get_pocket_info,
    _load_control_bases,
)


def _collect_master_rows(
    *,
    consensus_files: List[Path],
    post_root: Path,
    repo_root: Path,
    run_id: str,
    decoy_prefix: str,
    manifest: Optional[Dict[str, Any]],
    processed_root: Path,
    pb_map: Dict[Tuple[str, str, str, str], Tuple[str, str]],
    dud_map: Dict[Tuple[str, str, str], Dict[str, str]],
    fda_index: Any,
    logger: logging.Logger,
) -> Tuple[List[Dict[str, Any]], Dict[str, Set[str]]]:
    control_caches: Dict[str, Set[str]] = {}
    master_rows: List[Dict[str, Any]] = []
    ligand_exposure_report_cache: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for csv_path in sorted(consensus_files):
        try:
            with csv_path.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
        except _SCHEMA_CSV_ERRORS as e:
            logger.warning(
                "%s action=read_error path=%s error=%s", COMPONENT, csv_path, e
            )
            continue

        if not rows:
            continue

        # Context from path or row
        # Path: post_docked/<runid>/<pdb>/<variant>/<ph>/...
        try:
            rel = csv_path.relative_to(post_root / run_id)
            pdb_id = rel.parts[0]
            variant = rel.parts[1]
            ph = rel.parts[2]
        except (ValueError, IndexError):
            # Fallback to first row
            pdb_id = rows[0].get("pdb_id", "")
            variant = rows[0].get("variant", "")
            ph = rows[0].get("ph_label", "") or rows[0].get("ph", "")
        pdb_id = str(pdb_id).upper()
        variant = str(variant).upper()
        ph = str(ph)

        if pdb_id not in control_caches:
            control_caches[pdb_id] = _load_control_bases(processed_root, pdb_id)
        controls = control_caches[pdb_id]

        pocket = _get_pocket_info(manifest, pdb_id, variant, ph)

        # Try finding DUD metrics
        metrics = dud_map.get((pdb_id, variant, ph)) or dud_map.get(
            (pdb_id, variant, ph.strip())
        )
        if not metrics and ph:
            metrics = dud_map.get((pdb_id, variant, ph.lower())) or dud_map.get(
                (pdb_id, variant, "")
            )
        if not metrics:
            metrics = {}
        ef1_val = str(metrics.get("ef1", "")).strip()
        roc_auc_val = str(metrics.get("roc_auc", "")).strip()
        roc_auc_adj_val = str(metrics.get("roc_auc_adj", "") or roc_auc_val).strip()
        if not ef1_val:
            ef1_val = "0"
        if not roc_auc_val:
            roc_auc_val = "0"
        if not roc_auc_adj_val:
            roc_auc_adj_val = roc_auc_val
        status_reason = str(metrics.get("dud_eval_status_reason", "") or "").strip()

        source_rel = str(csv_path.relative_to(repo_root))

        for row in rows:
            # Base identity
            lig_file = (
                row.get("ligand_file")
                or row.get("ligand", "")
                or row.get("Ligand_ID", "")
            )
            if not lig_file:
                continue
            lig_file_name = Path(lig_file).name
            base = _canonical_ligand_base_with_prefix(
                Path(lig_file_name).stem, decoy_prefix
            )
            ligand_display = base
            if fda_index is not None:
                ligand_display = resolve_ligand_display_name(
                    base, lig_file_name, fda_index
                )
            if not ligand_display:
                ligand_display = base

            # Decoy/Control
            is_control = base in controls
            is_decoy = False if is_control else _row_is_decoy(row, decoy_prefix)
            library_name = _resolve_library_name(
                row, source_rel, decoy_prefix, is_decoy=is_decoy
            )

            # Pose Validity
            pb_valid, pb_reason = pb_map.get((pdb_id, variant, ph, base), ("", ""))

            # Z-scores and selection
            z_stage2 = ""
            z_stage2_source = ""
            # prefer z_vs_decoys_blend, else z_vs_decoys_consensus
            z_blend = as_float(row.get("z_vs_decoys_blend"))
            z_cons = as_float(row.get("z_vs_decoys_consensus"))
            if z_blend is not None and math.isfinite(z_blend):
                z_stage2 = str(z_blend)
                z_stage2_source = "z_vs_decoys_blend"
            elif z_cons is not None and math.isfinite(z_cons):
                z_stage2 = str(z_cons)
                z_stage2_source = "z_vs_decoys_consensus"

            z_stage1_pre = row.get("z_vs_decoys_consensus_pre", "")
            z_stage1 = z_stage1_pre or row.get("z_vs_decoys_consensus", "")
            z_stage1_source = (
                "z_vs_decoys_consensus_pre"
                if z_stage1_pre
                else "z_vs_decoys_consensus"
            )

            z_selected = ""
            z_source = ""

            if z_stage2:
                z_selected = z_stage2
                z_source = z_stage2_source
            elif z_stage1:
                z_selected = z_stage1
                z_source = z_stage1_source
            else:
                raw_consensus = as_float(row.get("consensus_score"))
                if raw_consensus is not None and math.isfinite(raw_consensus):
                    z_source = "missing_consensus_decoy_null"

            if pb_valid == "0":
                z_source = f"{z_source}_pose_invalid"

            # Construct Master Row
            new_row = {
                "run_id": run_id,
                "pdb_id": pdb_id,
                "variant": variant,
                "ph_label": ph,
                "library": library_name,
                "run_mode": row.get("run_mode", ""),
                "ligand": row.get("ligand", ""),
                "ligand_file": lig_file_name,
                "ligand_base": base,
                "ligand_display": ligand_display,
                "is_decoy": 1 if is_decoy else 0,
                "is_control": 1 if is_control else 0,
                "pose_valid_any": pb_valid,
                "pose_invalid_reason_top": pb_reason,
                "pocket_method": pocket.get("pocket_method", ""),
                "center_x": pocket.get("center_x", ""),
                "center_y": pocket.get("center_y", ""),
                "center_z": pocket.get("center_z", ""),
                "box_x": pocket.get("box_x", ""),
                "box_y": pocket.get("box_y", ""),
                "box_z": pocket.get("box_z", ""),
                "ef1": ef1_val,
                "roc_auc": roc_auc_val,
                "roc_auc_adj": roc_auc_adj_val,
                "dud_eval_status_reason": status_reason,
                "z_stage1": z_stage1,
                "z_stage2": z_stage2,
                "z_selected": z_selected,
                "z_selected_source": z_source,
                "source_csv": source_rel,
            }

            # Add remaining columns from source
            for k, v in row.items():
                if k not in new_row:
                    new_row[k] = v

            if not is_decoy:
                _ex_key = (ligand_display, base)
                cached_exposure = ligand_exposure_report_cache.get(_ex_key)
                if cached_exposure is None:
                    cached_exposure = ligand_side_effect_report_fields(
                        lookup_ligand_side_effect_entry(
                            repo_root,
                            ligand_label=ligand_display,
                            ligand_base=base,
                            ligand_display=ligand_display,
                        )
                    )
                    ligand_exposure_report_cache[_ex_key] = cached_exposure
                for exposure_key in (
                    "free_cmax_um",
                    "cmax_um",
                    "fraction_unbound_plasma",
                    "exposure_source",
                ):
                    if cached_exposure.get(exposure_key) not in (None, ""):
                        new_row.setdefault(exposure_key, cached_exposure[exposure_key])

            master_rows.append(new_row)
    return master_rows, control_caches
