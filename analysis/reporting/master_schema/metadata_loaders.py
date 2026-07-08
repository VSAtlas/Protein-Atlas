# -*- coding: utf-8 -*-
import csv
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from analysis.reporting.manifest_utils import extract_pocket
from analysis.reporting.value_utils import as_float
from analysis.reporting.master_schema.constants import (
    COMPONENT,
    DECOY_PREFIX_DEFAULT,
    _SCHEMA_CSV_ERRORS,
)
from analysis.reporting.master_schema.decoy_helpers import _canonical_ligand_base_with_prefix
from analysis.reporting.master_schema.io_utils import _try_read_utf8_lines_logged

def _load_control_bases(processed_root: Path, pdb_id: str) -> Set[str]:
    roots = [
        processed_root / pdb_id / "ligands_raw",
        processed_root / f"{pdb_id}_NOLIG" / "ligands_raw",
    ]
    bases: Set[str] = set()
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            # Canonicalize control base name
            stem = path.stem.split("_stage")[0]
            stem = stem.split(".sanitized")[0]
            if stem:
                bases.add(stem)
    return bases


def _get_pocket_info(
    manifest: Optional[Dict[str, Any]], pdb_id: str, variant: str, ph: str
) -> Dict[str, str]:
    info = {
        "pocket_method": "",
        "center_x": "",
        "center_y": "",
        "center_z": "",
        "box_x": "",
        "box_y": "",
        "box_z": "",
    }
    pocket = extract_pocket(manifest, pdb_id, variant, ph)
    if (
        manifest
        and pocket.get("method") is None
        and pocket.get("center") is None
        and pocket.get("box") is None
    ):
        proteins = manifest.get("proteins", {}) or {}
        entry = proteins.get(pdb_id) or proteins.get(pdb_id.upper())
        if isinstance(entry, dict):
            details = (
                entry.get("stages", {}).get("pocket_detection", {}).get("details", {})
            )
            pocket = {
                "method": details.get("method") or details.get("pocket_method"),
                "center": details.get("center")
                or [
                    details.get("center_x"),
                    details.get("center_y"),
                    details.get("center_z"),
                ],
                "box": details.get("box")
                or details.get("box_size")
                or [details.get("box_x"), details.get("box_y"), details.get("box_z")],
            }
    method = pocket.get("method")
    center = pocket.get("center")
    box = pocket.get("box")

    if method is not None:
        info["pocket_method"] = str(method).strip()

    if isinstance(center, (list, tuple)) and len(center) == 3:
        info["center_x"] = str(center[0])
        info["center_y"] = str(center[1])
        info["center_z"] = str(center[2])

    if isinstance(box, (list, tuple)) and len(box) == 3:
        info["box_x"] = str(box[0])
        info["box_y"] = str(box[1])
        info["box_z"] = str(box[2])

    return info


def _load_posebusters_map(
    post_root: Path, run_id: str, decoy_prefix: str = DECOY_PREFIX_DEFAULT
) -> Dict[Tuple[str, str, str, str], Tuple[str, str]]:
    # Returns {(pdb, variant, ph, ligand_base): (pose_valid_any, reason_top)}
    # Note: Using ligand_base for join stability.

    pb_map: Dict[Tuple[str, str, str, str], Tuple[str, str]] = {}

    # Glob all posebusters_all_stages.csv under post_docked/<runid>
    run_dir = post_root / run_id
    if not run_dir.exists():
        return pb_map

    for csv_path in run_dir.rglob("posebusters_all_stages.csv"):
        # Infer combo from path
        try:
            rel = csv_path.relative_to(run_dir)
            parts = rel.parts
            if len(parts) < 3:
                continue
            pdb_id, variant, ph = parts[0], parts[1], parts[2]
        except ValueError:
            continue

        try:
            with csv_path.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
        except _SCHEMA_CSV_ERRORS:
            continue

        # Group by ligand_base to find 'any' validity and top reason
        # Fields: ligand_file, posebusters_pass (TRUE/FALSE), posebusters_reason

        ligand_groups: Dict[str, List[Dict[str, str]]] = {}
        for row in rows:
            lig_file = row.get("ligand_file", "")
            if not lig_file:
                continue
            base = _canonical_ligand_base_with_prefix(Path(lig_file).stem, decoy_prefix)
            if not base:
                continue
            ligand_groups.setdefault(base, []).append(row)

        for base, group in ligand_groups.items():
            valid_any = "0"
            reasons: Dict[str, int] = {}
            for r in group:
                valid = str(r.get("posebusters_pass", "")).strip().lower() in (
                    "true",
                    "1",
                    "yes",
                )
                if valid:
                    valid_any = "1"
                reason = str(r.get("posebusters_reason", "")).strip()
                if reason and reason != "None":
                    reasons[reason] = reasons.get(reason, 0) + 1

            top_reason = ""
            if reasons:
                top_reason = max(reasons.items(), key=lambda x: x[1])[0]

            pb_map[(pdb_id, variant, ph, base)] = (valid_any, top_reason)

    return pb_map


def _find_dud_eval_tsv(
    run_id: str, repo_root: Path, decoy_prefix: str = DECOY_PREFIX_DEFAULT
) -> Optional[Path]:
    analysis_eval = repo_root / "analysis" / f"{decoy_prefix}_eval" / run_id
    legacy_root = repo_root / f"{decoy_prefix}_eval" / run_id
    legacy_analysis = repo_root / "analysis" / "dud_eval" / run_id
    legacy_docked = repo_root / "docked" / f"{decoy_prefix}_eval" / run_id
    bases = []
    for base in (analysis_eval, legacy_root, legacy_analysis, legacy_docked):
        if base not in bases:
            bases.append(base)

    patterns = [
        f"consensus_reranked_scorch_summary_*{run_id}*.tsv",
        f"consensus_summary_*{run_id}*.tsv",
        f"consensus_reranked_scorch_summary_*{run_id}*_pretty.txt",
        f"consensus_summary_*{run_id}*_pretty.txt",
    ]
    for base in bases:
        if not base.exists():
            continue
        for pattern in patterns:
            matches = sorted(base.rglob(pattern))
            if matches:
                return matches[0]
    return None


def _parse_dud_eval_pretty(
    summary_path: Path, run_id: str, logger: logging.Logger
) -> Dict[Tuple[str, str, str], Dict[str, str]]:
    metrics_map: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    lines = _try_read_utf8_lines_logged(
        summary_path, logger, action="dud_eval_pretty"
    )
    if lines is None:
        return metrics_map

    for line in lines:
        parts = line.strip().split()
        if len(parts) < 6:
            continue
        if parts[0].lower().startswith("header"):
            continue
        variant = parts[0].strip().upper()
        ph = parts[1].strip()
        row_run = parts[2].strip()
        pdb_id = parts[3].strip().upper()
        if row_run and row_run != run_id:
            continue
        nums: List[float] = []
        rest: List[str] = []
        for tok in parts[4:]:
            try:
                nums.append(float(tok))
            except (TypeError, ValueError):
                rest.append(tok)
        if not nums:
            continue
        roc_auc_adj = nums[0] if len(nums) >= 1 else 0.0
        roc_auc = nums[1] if len(nums) >= 2 else roc_auc_adj
        ef1 = nums[3] if len(nums) >= 4 else 0.0
        status_reason = " ".join(rest).strip()
        metrics_map[(pdb_id, variant, ph)] = {
            "ef1": f"{ef1}",
            "roc_auc": f"{roc_auc}",
            "roc_auc_adj": f"{roc_auc_adj}",
            "dud_eval_status_reason": status_reason,
        }

    logger.info(
        "%s action=dud_eval_pretty status=loaded path=%s rows=%d",
        COMPONENT,
        summary_path,
        len(metrics_map),
    )
    return metrics_map


def _metric_with_default(
    row: Dict[str, Any], keys: Tuple[str, ...], default: Optional[float]
) -> str:
    for key in keys:
        if key in row:
            val = as_float(row.get(key))
            if val is not None:
                return f"{val}"
            raw = row.get(key)
            if raw is not None and str(raw).strip():
                return str(raw).strip()
    if default is None:
        return ""
    return f"{default}"


def _parse_dud_eval_summary(
    run_id: str,
    repo_root: Path,
    logger: logging.Logger,
    decoy_prefix: str = DECOY_PREFIX_DEFAULT,
) -> Dict[Tuple[str, str, str], Dict[str, str]]:
    # Returns {(pdb_id, variant, ph): {metrics...}}
    metrics_map: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    summary_path = _find_dud_eval_tsv(run_id, repo_root, decoy_prefix)
    if not summary_path:
        logger.info(
            "%s action=dud_eval_summary status=missing run_id=%s", COMPONENT, run_id
        )
        return metrics_map
    if summary_path.suffix.lower() == ".txt":
        return _parse_dud_eval_pretty(summary_path, run_id, logger)

    lines = _try_read_utf8_lines_logged(
        summary_path, logger, action="dud_eval_summary"
    )
    if lines is None:
        return metrics_map

    header_idx = None
    for idx, line in enumerate(lines):
        if "\t" not in line:
            continue
        lower = line.lower()
        if "variant" in lower and ("pdb" in lower or "target" in lower):
            header_idx = idx
            break

    if header_idx is None:
        logger.warning(
            "%s action=dud_eval_summary status=no_header path=%s",
            COMPONENT,
            summary_path,
        )
        return metrics_map

    reader = csv.DictReader(lines[header_idx:], delimiter="\t")
    for row in reader:
        if not row:
            continue
        row_run = (row.get("run_id") or row.get("run") or "").strip()
        if row_run and row_run != run_id:
            continue

        pdb_id = (row.get("pdb_id") or row.get("target_name") or "").strip().upper()
        variant = (row.get("variant") or row.get("Variant") or "").strip().upper()
        ph = (row.get("pH") or row.get("ph_label") or row.get("ph") or "").strip()
        if not pdb_id or not variant:
            continue

        ef1 = _metric_with_default(row, ("EF@1%", "EF1"), 0.0)
        roc_auc = _metric_with_default(row, ("ROC_AUC", "ROC-AUC", "roc_auc"), 0.0)
        roc_auc_adj = _metric_with_default(
            row, ("ROC_AUC_adj", "ROC-AUC_adj", "roc_auc_adj"), None
        )
        if not roc_auc_adj:
            roc_auc_adj = roc_auc
        status_reason = str(row.get("status_reason") or row.get("status") or "").strip()

        metrics_map[(pdb_id, variant, ph)] = {
            "ef1": ef1,
            "roc_auc": roc_auc,
            "roc_auc_adj": roc_auc_adj,
            "dud_eval_status_reason": status_reason,
        }

        if not ef1 or not roc_auc:
            logger.info(
                "%s action=dud_eval_summary status=fallback pdb=%s variant=%s ph=%s reason=missing_metrics",
                COMPONENT,
                pdb_id,
                variant,
                ph,
            )

    logger.info(
        "%s action=dud_eval_summary status=loaded path=%s rows=%d",
        COMPONENT,
        summary_path,
        len(metrics_map),
    )
    return metrics_map
