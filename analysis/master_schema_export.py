# -*- coding: utf-8 -*-
import argparse
import csv
import logging
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from analysis.manifest_utils import extract_pocket, load_run_manifest

# Try to import canonical_ligand_base from rescore_reranker
try:
    from rescore_reranker import canonical_ligand_base
except ImportError:
    # Fallback if import fails (should not happen if rescore_reranker is in path)
    def canonical_ligand_base(lig: str) -> str:
        s = str(lig or "").strip()
        s = s.replace(".sanitized", "")
        s = re.sub(r"(_dud_gnina_stage\d+)$", "", s)
        s = re.sub(r"(_gnina_dud_stage\d+)$", "", s)
        s = re.sub(r"(_dock6_dud_stage\d+)$", "", s)
        s = re.sub(r"(_dud_dock6_stage\d+)$", "", s)
        s = re.sub(r"(_dud_stage\d+)$", "", s)
        s = re.sub(r"(__dock6_dud_stage\d+)$", "", s)
        s = re.sub(r"(__ledock_stage\d+)$", "", s)
        s = re.sub(r"(__dock6_stage\d+)$", "", s)
        s = re.sub(r"(_gnina_stage\d+)$", "", s)
        s = re.sub(r"(_stage\d+)$", "", s)
        s = re.sub(r"\.(mol2|pdbqt)$", "", s, flags=re.IGNORECASE)
        s = s.replace("__", "_")
        s = re.sub(r"_+$", "", s)
        return s


COMPONENT = "[master-export]"
_DECOY_RE = re.compile(r"\bdecoys?_", re.IGNORECASE)
MIN_DECOYS_FOR_FDR = 200
MIN_UNIQUE_DECOY_SCORES = 10


def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("master-export")


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            return list(reader)
    except Exception:
        return []


def _is_decoy_file(filename: str) -> bool:
    return _DECOY_RE.search(filename or "") is not None


def _row_is_decoy(row: Dict[str, Any]) -> bool:
    library = str(row.get("library", "")).strip().upper()
    if library == "DECOY":
        return True
    if str(row.get("is_decoy", "")).strip() == "1":
        return True
    ligand_file = (
        row.get("ligand_file") or row.get("ligand") or row.get("Ligand_ID") or ""
    )
    try:
        lig_name = Path(str(ligand_file)).name
    except Exception:
        lig_name = str(ligand_file)
    return _is_decoy_file(lig_name)


def _choose_fdr_score_field(
    decoy_rows: List[Dict[str, Any]], min_decoys: int = MIN_DECOYS_FOR_FDR
) -> Optional[str]:
    counts: Dict[str, int] = {}
    for field in ("consensus_score", "consensus_score_pre"):
        vals = []
        for r in decoy_rows:
            v = _as_float(r.get(field))
            if v is None or not math.isfinite(v):
                continue
            vals.append(v)
        counts[field] = len(vals)
    preferred_field = None
    preferred_count = -1
    for field, count in counts.items():
        if count >= min_decoys and count > preferred_count:
            preferred_field = field
            preferred_count = count
    if preferred_field:
        return preferred_field

    if counts:
        fallback_field = max(counts.items(), key=lambda kv: kv[1])[0]
        if counts[fallback_field] > 0:
            return fallback_field

    return "consensus_score"


def _get_decoy_scores_for_combo(
    combo_dir: Path, fallback_rows: List[Dict[str, Any]], min_decoys: int = 1
) -> Tuple[Optional[str], List[float], bool, Dict[str, List[float]]]:
    dud_path = combo_dir / "dud_consensus_reranked_scorch.csv"
    fallback_used = False
    decoy_rows: List[Dict[str, Any]] = []
    if dud_path.exists():
        decoy_rows = _read_csv_rows(dud_path)
    if not decoy_rows:
        fallback_used = True
        decoy_rows = fallback_rows
    decoys_only = [r for r in decoy_rows if _row_is_decoy(r)]
    if not decoys_only:
        fallback_used = True
        decoys_only = [r for r in fallback_rows if _row_is_decoy(r)]
    scores_by_field: Dict[str, List[float]] = {}
    for field_name in ("consensus_score", "consensus_score_pre"):
        vals: List[float] = []
        for r in decoys_only:
            val = _as_float(r.get(field_name))
            if val is None or not math.isfinite(val):
                continue
            vals.append(val)
        scores_by_field[field_name] = vals

    field = _choose_fdr_score_field(decoys_only, min_decoys)
    scores = scores_by_field.get(field, []) if field else []
    return field, scores, fallback_used, scores_by_field


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
    post_root: Path, run_id: str
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
        except Exception:
            continue

        # Group by ligand_base to find 'any' validity and top reason
        # Fields: ligand_file, posebusters_pass (TRUE/FALSE), posebusters_reason

        ligand_groups: Dict[str, List[Dict[str, str]]] = {}
        for row in rows:
            lig_file = row.get("ligand_file", "")
            if not lig_file:
                continue
            base = canonical_ligand_base(Path(lig_file).stem)
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


def _find_dud_eval_tsv(run_id: str, repo_root: Path) -> Optional[Path]:
    base = repo_root / "analysis" / "dud_eval" / run_id
    if not base.exists():
        return None

    patterns = [
        f"consensus_reranked_scorch_summary_*{run_id}*.tsv",
        f"consensus_summary_*{run_id}*.tsv",
    ]
    for pattern in patterns:
        matches = sorted(base.rglob(pattern))
        if matches:
            return matches[0]
    return None


def _metric_with_default(
    row: Dict[str, Any], keys: Tuple[str, ...], default: Optional[float]
) -> str:
    for key in keys:
        if key in row:
            val = _as_float(row.get(key))
            if val is not None:
                return f"{val}"
            raw = row.get(key)
            if raw is not None and str(raw).strip():
                return str(raw).strip()
    if default is None:
        return ""
    return f"{default}"


def _parse_dud_eval_summary(
    run_id: str, repo_root: Path, logger: logging.Logger
) -> Dict[Tuple[str, str, str], Dict[str, str]]:
    # Returns {(pdb_id, variant, ph): {metrics...}}
    metrics_map: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    summary_path = _find_dud_eval_tsv(run_id, repo_root)
    if not summary_path:
        logger.info(
            "%s action=dud_eval_summary status=missing run_id=%s", COMPONENT, run_id
        )
        return metrics_map

    try:
        with summary_path.open("r", encoding="utf-8") as f:
            lines = f.readlines()
    except Exception as exc:
        logger.warning(
            "%s action=dud_eval_summary status=read_failed path=%s error=%s",
            COMPONENT,
            summary_path,
            exc,
        )
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


def _as_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        s = str(x).strip()
        if not s:
            return None
        return float(s)
    except Exception:
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Export master CSV for a docking run.")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    run_id = args.run_id

    post_root = repo_root / "post_docked"
    processed_root = repo_root / "processed_pdbs"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True, exist_ok=True)

    output_csv = data_dir / "master_rows.csv"
    if output_csv.exists() and not args.overwrite:
        logger.info("%s action=skip reason=exists path=%s", COMPONENT, output_csv)
        return 0

    # Discovery
    consensus_files = list((post_root / run_id).rglob("consensus_reranked_scorch.csv"))
    logger.info(
        "%s action=discover count=%d path_pattern=%s",
        COMPONENT,
        len(consensus_files),
        f"post_docked/{run_id}/**/consensus_reranked_scorch.csv",
    )

    if not consensus_files:
        logger.warning("%s action=exit reason=no_consensus_files", COMPONENT)
        return 0

    # Load metadata sources
    manifest, _manifest_path = load_run_manifest(repo_root, run_id)
    pb_map = _load_posebusters_map(post_root, run_id)
    dud_map = _parse_dud_eval_summary(run_id, repo_root, logger)
    control_caches: Dict[str, Set[str]] = {}

    master_rows: List[Dict[str, Any]] = []

    # Process files
    for csv_path in sorted(consensus_files):
        try:
            with csv_path.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
        except Exception as e:
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
            base = canonical_ligand_base(Path(lig_file_name).stem)

            # Decoy/Control
            is_decoy = _is_decoy_file(lig_file_name)
            is_control = base in controls

            # Pose Validity
            pb_valid, pb_reason = pb_map.get((pdb_id, variant, ph, base), ("", ""))

            # T-Scores and selection
            t_stage2 = ""
            # prefer t_vs_decoys_blend, else t_vs_decoys_consensus
            t_blend = _as_float(row.get("t_vs_decoys_blend"))
            t_cons = _as_float(row.get("t_vs_decoys_consensus"))
            if t_blend is not None and math.isfinite(t_blend):
                t_stage2 = str(t_blend)
            elif t_cons is not None and math.isfinite(t_cons):
                t_stage2 = str(t_cons)

            t_stage1 = row.get("t_vs_decoys_consensus_pre", "") or row.get(
                "t_vs_decoys_consensus", ""
            )
            if t_stage1 == t_stage2 and not row.get("t_vs_decoys_consensus_pre"):
                # If stage1 score is missing column-wise, sometimes it's implied same if not reranked?
                # Spec says: use t_vs_decoys_consensus_pre if present.
                pass

            t_selected = ""
            t_source = ""

            valid_bool = pb_valid == "1"

            if valid_bool and t_stage2:
                t_selected = t_stage2
                t_source = "stage2"
            elif t_stage1:
                t_selected = t_stage1
                t_source = "stage1"
            elif t_stage2:
                # Fallback if valid is false but stage2 exists? Spec says "if pose_valid_any==1 use stage2 else stage1"
                # So if invalid, use stage1.
                t_selected = t_stage2  # Wait, if invalid, we fallback to stage1.
                t_source = "stage2"  # Correction below

            # Correct logic:
            # if pose_valid_any==1, use t_stage2
            # else use t_stage1 (or blank)
            if valid_bool:
                t_selected = t_stage2
                t_source = "stage2"
            else:
                t_selected = t_stage1
                t_source = "stage1"

            if not t_selected:
                t_selected = ""
                t_source = ""

            # Construct Master Row
            new_row = {
                "run_id": run_id,
                "pdb_id": pdb_id,
                "variant": variant,
                "ph_label": ph,
                "library": row.get("library", ""),
                "run_mode": row.get("run_mode", ""),
                "ligand": row.get("ligand", ""),
                "ligand_file": lig_file_name,
                "ligand_base": base,
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
                "t_stage1": t_stage1,
                "t_stage2": t_stage2,
                "t_selected": t_selected,
                "t_selected_source": t_source,
                "source_csv": source_rel,
            }

            # Add remaining columns from source
            for k, v in row.items():
                if k not in new_row:
                    new_row[k] = v

            master_rows.append(new_row)

    if not master_rows:
        logger.warning("%s action=exit reason=no_rows_collected", COMPONENT)
        return 0

    # Initialize FDR fields
    for r in master_rows:
        r.setdefault("fdr_score_field", "")
        r.setdefault("fdr_n_decoys", "")
        r.setdefault("fdr_unique_decoy_scores", "")
        r.setdefault("fdr_reliable", "")
        r.setdefault("fdr_p_empirical", "")
        r.setdefault("fdr_q_target", "")
        r.setdefault("fdr_hit_q05", "")
        r.setdefault("fdr_hit_q10", "")

    combo_indices: Dict[Tuple[str, str, str], List[int]] = {}
    for idx, row in enumerate(master_rows):
        key = (row["pdb_id"], row["variant"], row["ph_label"])
        combo_indices.setdefault(key, []).append(idx)

    fdr_summary_rows: List[Dict[str, Any]] = []
    for (pdb_id, variant, ph), indices in combo_indices.items():
        combo_dir = post_root / run_id / pdb_id / variant / ph
        dud_path = combo_dir / "dud_consensus_reranked_scorch.csv"
        fallback_rows = [master_rows[i] for i in indices]
        (
            field,
            decoy_scores,
            used_fallback,
            _decoy_scores_by_field,
        ) = _get_decoy_scores_for_combo(combo_dir, fallback_rows, MIN_DECOYS_FOR_FDR)
        n_decoys = len(decoy_scores)
        unique_scores = len(set(decoy_scores))
        reliable = (
            n_decoys >= MIN_DECOYS_FOR_FDR and unique_scores >= MIN_UNIQUE_DECOY_SCORES
        )

        if used_fallback:
            logger.warning(
                "%s action=fdr_decoy_fallback pdb=%s variant=%s ph=%s path=%s reason=%s",
                COMPONENT,
                pdb_id,
                variant,
                ph,
                dud_path,
                "missing_or_empty" if not dud_path.exists() else "no_decoys_found",
            )

        field_for_rows = field or "consensus_score"

        non_decoy_indices = [
            idx for idx in indices if not _row_is_decoy(master_rows[idx])
        ]
        for idx in non_decoy_indices:
            master_rows[idx]["fdr_score_field"] = field_for_rows
            master_rows[idx]["fdr_n_decoys"] = str(n_decoys)
            master_rows[idx]["fdr_unique_decoy_scores"] = str(unique_scores)
            master_rows[idx]["fdr_reliable"] = "1" if reliable else "0"

        tested: List[Tuple[int, float]] = []
        for idx in non_decoy_indices:
            row = master_rows[idx]
            score_val = _as_float(row.get(field_for_rows))
            if score_val is None or not math.isfinite(score_val):
                continue
            p_val = (1 + sum(1 for d in decoy_scores if d >= score_val)) / (
                1 + n_decoys
            )
            tested.append((idx, p_val))
            row["fdr_p_empirical"] = f"{p_val:.6g}"

        hits05 = 0
        hits10 = 0
        best_q = 1.0
        tested_sorted: List[Tuple[int, float]] = []
        n_tested_rows = 0
        if tested:
            tested_sorted = sorted(tested, key=lambda x: x[1])
            m = len(tested_sorted)
            n_tested_rows = m
            raw_q: List[Tuple[int, float]] = []
            for rank, (idx, p_val) in enumerate(tested_sorted, start=1):
                raw_q.append((idx, p_val * m / rank))

            q_values: Dict[int, float] = {}
            prev = 1.0
            for idx, q_val in reversed(raw_q):
                adj = min(q_val, prev, 1.0)
                prev = adj
                q_values[idx] = adj

            for idx, q_val in q_values.items():
                master_rows[idx]["fdr_q_target"] = f"{q_val:.6g}"
                hit05 = q_val <= 0.05
                hit10 = q_val <= 0.10
                master_rows[idx]["fdr_hit_q05"] = "1" if hit05 else "0"
                master_rows[idx]["fdr_hit_q10"] = "1" if hit10 else "0"
                if hit05:
                    hits05 += 1
                if hit10:
                    hits10 += 1
                best_q = min(best_q, q_val)
        else:
            n_tested_rows = len(non_decoy_indices)
            for idx in non_decoy_indices:
                if not master_rows[idx].get("fdr_p_empirical"):
                    master_rows[idx]["fdr_p_empirical"] = "1"
                master_rows[idx]["fdr_q_target"] = (
                    master_rows[idx].get("fdr_q_target", "") or "1"
                )
                master_rows[idx]["fdr_hit_q05"] = (
                    master_rows[idx].get("fdr_hit_q05", "") or "0"
                )
                master_rows[idx]["fdr_hit_q10"] = (
                    master_rows[idx].get("fdr_hit_q10", "") or "0"
                )

        for idx in non_decoy_indices:
            if not master_rows[idx].get("fdr_p_empirical"):
                master_rows[idx]["fdr_p_empirical"] = "1"
            if not master_rows[idx].get("fdr_q_target"):
                master_rows[idx]["fdr_q_target"] = "1"
            if not master_rows[idx].get("fdr_hit_q05"):
                master_rows[idx]["fdr_hit_q05"] = "0"
            if not master_rows[idx].get("fdr_hit_q10"):
                master_rows[idx]["fdr_hit_q10"] = "0"

        fdr_summary_rows.append(
            {
                "pdb_id": pdb_id,
                "variant": variant,
                "ph_label": ph,
                "fdr_score_field": field_for_rows,
                "fdr_n_decoys": str(n_decoys),
                "fdr_unique_decoy_scores": str(unique_scores),
                "fdr_reliable": "1" if reliable else "0",
                "n_tested": str(n_tested_rows),
                "n_hits_q05": str(hits05),
                "n_hits_q10": str(hits10),
                "best_q": f"{best_q:.6g}" if math.isfinite(best_q) else "",
            }
        )

        logger.info(
            "%s action=fdr_compute pdb=%s variant=%s ph=%s field=%s n_decoys=%d n_tested=%d n_hits_q05=%d",
            COMPONENT,
            pdb_id,
            variant,
            ph,
            field_for_rows,
            n_decoys,
            n_tested_rows,
            hits05,
        )

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
        "t_stage1",
        "t_stage2",
        "t_selected",
        "t_selected_source",
        "fdr_score_field",
        "fdr_n_decoys",
        "fdr_unique_decoy_scores",
        "fdr_reliable",
        "fdr_p_empirical",
        "fdr_q_target",
        "fdr_hit_q05",
        "fdr_hit_q10",
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
    except Exception as e:
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
                        "fdr_n_decoys",
                        "fdr_unique_decoy_scores",
                        "fdr_reliable",
                        "n_tested",
                        "n_hits_q05",
                        "n_hits_q10",
                        "best_q",
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
    except Exception as e:
        logger.warning(
            "%s action=write_optional status=failed path=target_fdr_summary.csv error=%s",
            COMPONENT,
            e,
        )

    # Optional: Ligand Top Targets
    try:
        ligand_targets: Dict[str, List[Dict[str, Any]]] = {}
        for r in master_rows:
            base = r["ligand_base"]
            ligand_targets.setdefault(base, []).append(r)

        top_targets_rows = []
        for base, targets in ligand_targets.items():
            # Sort by t_selected descending
            sorted_targets = sorted(
                targets,
                key=lambda x: _as_float(x.get("t_selected")) or -1e18,
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
                        "t_selected": target["t_selected"],
                        "t_stage2": target["t_stage2"],
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
                        "t_selected",
                        "t_stage2",
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
    except Exception as e:
        logger.warning(
            "%s action=write_optional status=failed path=ligand_top_targets.csv error=%s",
            COMPONENT,
            e,
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
                _as_float(r.get("t_selected"))
                for r in rows
                if _as_float(r.get("t_selected")) is not None
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
                    "mean_t_selected": f"{mean_t:.3f}",
                    "max_t_selected": f"{max_t:.3f}",
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
                        "mean_t_selected",
                        "max_t_selected",
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
    except Exception as e:
        logger.warning(
            "%s action=write_optional status=failed path=target_qc_summary.csv error=%s",
            COMPONENT,
            e,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
