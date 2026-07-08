# -*- coding: utf-8 -*-
import csv
import logging
from pathlib import Path
from typing import Any, Dict, List, Tuple

from analysis.reporting.value_utils import as_float
from analysis.reporting.master_schema.constants import COMPONENT, _SCHEMA_WRITE_ERRORS


def _write_master_exports(
    *,
    master_rows: List[Dict[str, Any]],
    output_csv: Path,
    data_dir: Path,
    fdr_summary_rows: List[Dict[str, Any]],
    scorch_fdr_summary_rows: List[Dict[str, Any]],
    logger: logging.Logger,
) -> int:
    # Sort
    master_rows.sort(
        key=lambda r: (
            r["pdb_id"],
            r["variant"],
            r["ph_label"],
            r["library"],
            r["ligand_base"],
        )
    )

    # Write Master CSV
    base_cols = [
        "run_id",
        "pdb_id",
        "variant",
        "ph_label",
        "library",
        "run_mode",
        "ligand",
        "ligand_file",
        "ligand_base",
        "ligand_display",
        "is_decoy",
        "is_control",
        "pose_valid_any",
        "pose_invalid_reason_top",
        "pocket_method",
        "center_x",
        "center_y",
        "center_z",
        "box_x",
        "box_y",
        "box_z",
        "ef1",
        "roc_auc",
        "roc_auc_adj",
        "dud_eval_status_reason",
        "z_stage1",
        "z_stage2",
        "z_selected",
        "z_selected_source",
        "fdr_score_field",
        "fdr_null_source",
        "fdr_n_decoys",
        "fdr_unique_decoy_scores",
        "fdr_reliable",
        "fdr_p_empirical",
        "fdr_q_target",
        "fdr_hit_q05",
        "fdr_hit_q10",
        "fdr_decoy_competition_q",
        "fdr_decoy_competition_hit_q05",
        "fdr_decoy_competition_hit_q10",
        "scorch_fdr_scope",
        "scorch_fdr_n_decoys",
        "scorch_fdr_n_tested",
        "scorch_fdr_p_empirical",
        "scorch_fdr_q_bh",
        "scorch_fdr_decoy_competition_q",
        "scorch_fdr_decoy_competition_q_plus1",
        "scorch_fdr_decoy_competition_hit_q05",
        "scorch_fdr_decoy_competition_hit_q10",
        "scorch_fdr_decoy_competition_hit_q05_plus1",
        "scorch_fdr_decoy_competition_hit_q10_plus1",
        "source_csv",
    ]

    # Collect all dynamic columns
    all_keys = set().union(*(d.keys() for d in master_rows))
    extra_cols = sorted([k for k in all_keys if k not in base_cols])
    final_cols = base_cols + extra_cols

    try:
        with output_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=final_cols)
            writer.writeheader()
            writer.writerows(master_rows)
        logger.info(
            "%s action=write status=ok path=%s rows=%d",
            COMPONENT,
            output_csv,
            len(master_rows),
        )
    except _SCHEMA_WRITE_ERRORS as e:
        logger.error(
            "%s action=write status=failed path=%s error=%s", COMPONENT, output_csv, e
        )
        return 1

    # Optional: Target FDR Summary
    try:
        if fdr_summary_rows:
            fdr_summary_csv = data_dir / "target_fdr_summary.csv"
            with fdr_summary_csv.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(
                    f,
                    fieldnames=[
                        "pdb_id",
                        "variant",
                        "ph_label",
                        "fdr_score_field",
                        "fdr_null_source",
                        "fdr_n_decoys",
                        "fdr_unique_decoy_scores",
                        "fdr_reliable",
                        "n_tested",
                        "n_hits_q05",
                        "n_hits_q10",
                        "best_q",
                        "min_possible_bh_q",
                        "bh_resolvable_q05",
                        "bh_resolvable_q10",
                        "n_decoy_competition_hits_q05",
                        "n_decoy_competition_hits_q10",
                        "best_decoy_competition_q",
                    ],
                )
                writer.writeheader()
                writer.writerows(fdr_summary_rows)
            logger.info(
                "%s action=write status=ok path=%s rows=%d",
                COMPONENT,
                fdr_summary_csv,
                len(fdr_summary_rows),
            )
        if scorch_fdr_summary_rows:
            scorch_fdr_summary_csv = data_dir / "target_scorch_fdr_summary.csv"
            with scorch_fdr_summary_csv.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(
                    f,
                    fieldnames=[
                        "pdb_id",
                        "variant",
                        "ph_label",
                        "score_field",
                        "scope",
                        "n_decoys",
                        "n_tested",
                        "best_bh_q",
                        "min_possible_bh_q",
                        "best_decoy_competition_q",
                        "best_decoy_competition_q_plus1",
                        "n_decoy_competition_hits_q05",
                        "n_decoy_competition_hits_q10",
                        "n_decoy_competition_hits_q05_plus1",
                        "n_decoy_competition_hits_q10_plus1",
                    ],
                )
                writer.writeheader()
                writer.writerows(scorch_fdr_summary_rows)
            logger.info(
                "%s action=write status=ok path=%s rows=%d",
                COMPONENT,
                scorch_fdr_summary_csv,
                len(scorch_fdr_summary_rows),
            )
    except _SCHEMA_WRITE_ERRORS as e:
        logger.warning(
            "%s action=write_optional status=failed path=target_fdr_summary.csv error=%s",
            COMPONENT,
            e,
            exc_info=True,
        )

    # Optional: Ligand Top Targets
    try:
        ligand_targets: Dict[str, List[Dict[str, Any]]] = {}
        for r in master_rows:
            base = r["ligand_base"]
            ligand_targets.setdefault(base, []).append(r)

        top_targets_rows = []
        for base, targets in ligand_targets.items():
            # Sort by z_selected descending
            sorted_targets = sorted(
                targets,
                key=lambda x: as_float(x.get("z_selected")) or -1e18,
                reverse=True,
            )
            for i, target in enumerate(sorted_targets[:25], 1):
                top_targets_rows.append(
                    {
                        "ligand_base": base,
                        "rank": i,
                        "pdb_id": target["pdb_id"],
                        "variant": target["variant"],
                        "ph_label": target["ph_label"],
                        "z_selected": target["z_selected"],
                        "z_stage2": target["z_stage2"],
                        "pose_valid_any": target["pose_valid_any"],
                        "library": target["library"],
                    }
                )

        if top_targets_rows:
            top_targets_csv = data_dir / "ligand_top_targets.csv"
            with top_targets_csv.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(
                    f,
                    fieldnames=[
                        "ligand_base",
                        "rank",
                        "pdb_id",
                        "variant",
                        "ph_label",
                        "z_selected",
                        "z_stage2",
                        "pose_valid_any",
                        "library",
                    ],
                )
                writer.writeheader()
                writer.writerows(top_targets_rows)
            logger.info(
                "%s action=write status=ok path=%s rows=%d",
                COMPONENT,
                top_targets_csv,
                len(top_targets_rows),
            )
    except _SCHEMA_WRITE_ERRORS as e:
        logger.warning(
            "%s action=write_optional status=failed path=ligand_top_targets.csv error=%s",
            COMPONENT,
            e,
            exc_info=True,
        )

    # Optional: Target QC Summary
    try:
        target_groups: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {}
        for r in master_rows:
            key = (r["pdb_id"], r["variant"], r["ph_label"])
            target_groups.setdefault(key, []).append(r)

        qc_rows = []
        for (pdb, variant, ph), rows in target_groups.items():
            n_rows = len(rows)
            n_fda = sum(
                1
                for r in rows
                if r.get("library") == "FDA" or r.get("run_mode") == "fda"
            )
            n_decoy = sum(1 for r in rows if r.get("is_decoy") == 1)
            n_controls = sum(1 for r in rows if r.get("is_control") == 1)
            n_valid = sum(1 for r in rows if r.get("pose_valid_any") == "1")

            t_vals = [
                t_val
                for r in rows
                if (t_val := as_float(r.get("z_selected"))) is not None
            ]
            mean_t = sum(t_vals) / len(t_vals) if t_vals else 0.0
            max_t = max(t_vals) if t_vals else 0.0

            # Use metrics from first row for that target
            ef1 = rows[0].get("ef1", "")
            roc_auc_adj = rows[0].get("roc_auc_adj", "")

            qc_rows.append(
                {
                    "pdb_id": pdb,
                    "variant": variant,
                    "ph_label": ph,
                    "n_rows": n_rows,
                    "n_fda": n_fda,
                    "n_decoy": n_decoy,
                    "n_controls": n_controls,
                    "n_pose_valid": n_valid,
                    "pose_valid_rate": f"{n_valid / n_rows:.3f}"
                    if n_rows > 0
                    else "0.000",
                    "mean_z_selected": f"{mean_t:.3f}",
                    "max_z_selected": f"{max_t:.3f}",
                    "ef1": ef1,
                    "roc_auc_adj": roc_auc_adj,
                }
            )

        if qc_rows:
            qc_csv = data_dir / "target_qc_summary.csv"
            with qc_csv.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(
                    f,
                    fieldnames=[
                        "pdb_id",
                        "variant",
                        "ph_label",
                        "n_rows",
                        "n_fda",
                        "n_decoy",
                        "n_controls",
                        "n_pose_valid",
                        "pose_valid_rate",
                        "mean_z_selected",
                        "max_z_selected",
                        "ef1",
                        "roc_auc_adj",
                    ],
                )
                writer.writeheader()
                writer.writerows(qc_rows)
            logger.info(
                "%s action=write status=ok path=%s rows=%d",
                COMPONENT,
                qc_csv,
                len(qc_rows),
            )
    except _SCHEMA_WRITE_ERRORS as e:
        logger.warning(
            "%s action=write_optional status=failed path=target_qc_summary.csv error=%s",
            COMPONENT,
            e,
            exc_info=True,
        )

    return 0
