# -*- coding: utf-8 -*-
import argparse
import csv
import logging
import math
import re
import sys
import yaml
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

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


def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("master-export")


def _is_decoy_file(filename: str) -> bool:
    return _DECOY_RE.search(filename or "") is not None


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


def _find_manifest(run_id: str, repo_root: Path) -> Optional[Dict[str, Any]]:
    candidates = [
        repo_root / "post_docked" / run_id / "run_manifest.yaml",
        repo_root / "docked" / run_id / "run_manifest.yaml",
        repo_root / "run_manifest.yaml",
    ]
    for path in candidates:
        if path.exists() and path.stat().st_size > 0:
            try:
                with path.open("r", encoding="utf-8") as f:
                    return yaml.safe_load(f)
            except Exception:
                pass
    return None


def _get_pocket_info(manifest: Optional[Dict[str, Any]], pdb_id: str, variant: str, ph: str) -> Dict[str, str]:
    info = {
        "pocket_method": "",
        "center_x": "",
        "center_y": "",
        "center_z": "",
        "box_x": "",
        "box_y": "",
        "box_z": "",
    }
    if not manifest:
        return info
    
    # Try direct pdb access or variant access
    # Structure varies: proteins -> <pdb> -> stages -> pocket_detection -> details
    # Or proteins -> <pdb> -> variants -> <variant> -> ...
    
    # Flatten variant if needed or check typical paths
    proteins = manifest.get("proteins", {})
    prot_data = proteins.get(pdb_id)
    if not prot_data:
        return info

    # Check for pocket details at protein level (common if no variants or shared)
    stages = prot_data.get("stages", {})
    pocket = stages.get("pocket_detection", {})
    details = pocket.get("details", {})
    
    # If empty, check variant specific? (Less common in current pipeline but possible)
    if not details:
        variants = prot_data.get("variants", {})
        var_data = variants.get(variant)
        if var_data:
            details = var_data.get("stages", {}).get("pocket_detection", {}).get("details", {})

    if details:
        info["pocket_method"] = str(details.get("pocket_method", "")).strip()
        for k in ["center_x", "center_y", "center_z", "box_x", "box_y", "box_z"]:
            val = details.get(k)
            if val is not None:
                info[k] = str(val)
    
    return info


def _load_posebusters_map(post_root: Path, run_id: str) -> Dict[Tuple[str, str, str, str], Tuple[str, str]]:
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
                valid = str(r.get("posebusters_pass", "")).strip().lower() in ("true", "1", "yes")
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


def _parse_dud_eval_summary(run_id: str, repo_root: Path) -> Dict[Tuple[str, str, str], Dict[str, str]]:
    # Returns {(variant, ph, pdb_id): {metrics...}}
    # Note: Keying by (variant, ph, pdb_id) to match rows.
    
    summary_path = repo_root / "analysis" / "dud_eval" / run_id / "post_docked" / f"consensus_reranked_scorch_summary_{run_id}_pretty.txt"
    metrics_map: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    
    if not summary_path.exists():
        return metrics_map
        
    with summary_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        
    # Skip header
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 8:
            continue
            
        # Parse from left
        variant = parts[0]
        ph = parts[1]
        rid = parts[2]
        if rid != run_id:
            continue
            
        # Parse from right
        # Order: ... BEDROC_alpha_20 EF@1% EF@2% EF@5% EF@10% Reason
        # Before BEDROC usually: ROC_AUC_adj ROC_AUC (or just one of them)
        
        status_reason = parts[-1]
        try:
            ef10 = parts[-2]
            ef5 = parts[-3]
            ef2 = parts[-4]
            ef1 = parts[-5]
            # bedroc = parts[-6] # Not storing bedroc for now per spec, or could store it
            
            # Look for ROC_AUC
            # Tokens before -6 are potential ROC/Target
            remaining = parts[3:-6]
            
            roc_auc = ""
            roc_auc_adj = ""
            target_guess = ""
            
            # Heuristic: traverse remaining from right to find numbers
            numeric_vals = []
            non_numeric_vals = []
            
            for token in reversed(remaining):
                try:
                    float(token)
                    numeric_vals.append(token)
                except ValueError:
                    non_numeric_vals.append(token)
            
            # Expect roc_auc_adj then roc_auc (reversed order in numeric_vals)
            if len(numeric_vals) >= 1:
                roc_auc_adj = numeric_vals[0]
            if len(numeric_vals) >= 2:
                roc_auc = numeric_vals[1]
                
            # Target guess is whatever is left in non-numeric, joined
            if non_numeric_vals:
                target_guess = " ".join(reversed(non_numeric_vals))
            else:
                # If everything was numeric? unlikely for target name
                pass
                
        except IndexError:
            continue
            
        # Store
        # If target_guess matches a PDB format, use it.
        # Otherwise, relying on variant/ph uniqueness might be risky if multiple targets share run_id/variant/ph?
        # Typically run_id is unique per study, but multiple targets can be in one run.
        # If target_guess is empty, we might skip.
        
        combo_key = (variant, ph, target_guess) # Ideally target_guess is PDB ID
        
        data = {
            "ef1": ef1,
            "ef2": ef2,
            "ef5": ef5,
            "ef10": ef10,
            "roc_auc": roc_auc,
            "roc_auc_adj": roc_auc_adj,
            "dud_eval_status_reason": status_reason
        }
        
        metrics_map[combo_key] = data
        
    return metrics_map


def _as_float(x: Any) -> Optional[float]:
    try:
        if x is None: return None
        s = str(x).strip()
        if not s: return None
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
    logger.info("%s action=discover count=%d path_pattern=%s", COMPONENT, len(consensus_files), f"post_docked/{run_id}/**/consensus_reranked_scorch.csv")
    
    if not consensus_files:
        logger.warning("%s action=exit reason=no_consensus_files", COMPONENT)
        return 0
        
    # Load metadata sources
    manifest = _find_manifest(run_id, repo_root)
    pb_map = _load_posebusters_map(post_root, run_id)
    dud_map = _parse_dud_eval_summary(run_id, repo_root)
    control_caches: Dict[str, Set[str]] = {}
    
    master_rows: List[Dict[str, Any]] = []
    
    # Process files
    for csv_path in sorted(consensus_files):
        try:
            with csv_path.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
        except Exception as e:
            logger.warning("%s action=read_error path=%s error=%s", COMPONENT, csv_path, e)
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

        if pdb_id not in control_caches:
            control_caches[pdb_id] = _load_control_bases(processed_root, pdb_id)
        controls = control_caches[pdb_id]
        
        pocket = _get_pocket_info(manifest, pdb_id, variant, ph)
        
        # Try finding DUD metrics
        # First try exact PDB match
        metrics = dud_map.get((variant, ph, pdb_id))
        if not metrics:
            # Fallback: if only one target in DUD map for this variant/ph? 
            # Risk of collision. For now, strict match.
            metrics = {}

        source_rel = str(csv_path.relative_to(repo_root))
        
        for row in rows:
            # Base identity
            lig_file = row.get("ligand_file") or row.get("ligand", "") or row.get("Ligand_ID", "")
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
                
            t_stage1 = row.get("t_vs_decoys_consensus_pre", "") or row.get("t_vs_decoys_consensus", "")
            if t_stage1 == t_stage2 and not row.get("t_vs_decoys_consensus_pre"):
                # If stage1 score is missing column-wise, sometimes it's implied same if not reranked? 
                # Spec says: use t_vs_decoys_consensus_pre if present.
                pass
            
            t_selected = ""
            t_source = ""
            
            valid_bool = (pb_valid == "1")
            
            if valid_bool and t_stage2:
                t_selected = t_stage2
                t_source = "stage2"
            elif t_stage1:
                t_selected = t_stage1
                t_source = "stage1"
            elif t_stage2:
                # Fallback if valid is false but stage2 exists? Spec says "if pose_valid_any==1 use stage2 else stage1"
                # So if invalid, use stage1.
                t_selected = t_stage2 # Wait, if invalid, we fallback to stage1.
                t_source = "stage2" # Correction below
            
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
                "ef1": metrics.get("ef1", ""),
                "roc_auc": metrics.get("roc_auc", ""),
                "roc_auc_adj": metrics.get("roc_auc_adj", ""),
                "dud_eval_status_reason": metrics.get("dud_eval_status_reason", ""),
                "t_stage1": t_stage1,
                "t_stage2": t_stage2,
                "t_selected": t_selected,
                "t_selected_source": t_source,
                "source_csv": source_rel
            }
            
            # Add remaining columns from source
            for k, v in row.items():
                if k not in new_row:
                    new_row[k] = v
            
            master_rows.append(new_row)
            
    if not master_rows:
        logger.warning("%s action=exit reason=no_rows_collected", COMPONENT)
        return 0
        
    # Sort
    master_rows.sort(key=lambda r: (r["pdb_id"], r["variant"], r["ph_label"], r["library"], r["ligand_base"]))
    
    # Write Master CSV
    base_cols = [
        "run_id", "pdb_id", "variant", "ph_label", "library", "run_mode",
        "ligand", "ligand_file", "ligand_base", "is_decoy", "is_control",
        "pose_valid_any", "pose_invalid_reason_top",
        "pocket_method", "center_x", "center_y", "center_z", "box_x", "box_y", "box_z",
        "ef1", "roc_auc", "roc_auc_adj", "dud_eval_status_reason",
        "t_stage1", "t_stage2", "t_selected", "t_selected_source",
        "source_csv"
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
        logger.info("%s action=write status=ok path=%s rows=%d", COMPONENT, output_csv, len(master_rows))
    except Exception as e:
        logger.error("%s action=write status=failed path=%s error=%s", COMPONENT, output_csv, e)
        return 1

    # Optional: Ligand Top Targets
    try:
        ligand_targets: Dict[str, List[Dict[str, Any]]] = {}
        for r in master_rows:
            base = r["ligand_base"]
            ligand_targets.setdefault(base, []).append(r)
        
        top_targets_rows = []
        for base, targets in ligand_targets.items():
            # Sort by t_selected descending
            sorted_targets = sorted(targets, key=lambda x: _as_float(x.get("t_selected")) or -1e18, reverse=True)
            for i, target in enumerate(sorted_targets[:25], 1):
                top_targets_rows.append({
                    "ligand_base": base,
                    "rank": i,
                    "pdb_id": target["pdb_id"],
                    "variant": target["variant"],
                    "ph_label": target["ph_label"],
                    "t_selected": target["t_selected"],
                    "t_stage2": target["t_stage2"],
                    "pose_valid_any": target["pose_valid_any"],
                    "library": target["library"]
                })
        
        if top_targets_rows:
            top_targets_csv = data_dir / "ligand_top_targets.csv"
            with top_targets_csv.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["ligand_base", "rank", "pdb_id", "variant", "ph_label", "t_selected", "t_stage2", "pose_valid_any", "library"])
                writer.writeheader()
                writer.writerows(top_targets_rows)
            logger.info("%s action=write status=ok path=%s rows=%d", COMPONENT, top_targets_csv, len(top_targets_rows))
    except Exception as e:
        logger.warning("%s action=write_optional status=failed path=ligand_top_targets.csv error=%s", COMPONENT, e)

    # Optional: Target QC Summary
    try:
        target_groups: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {}
        for r in master_rows:
            key = (r["pdb_id"], r["variant"], r["ph_label"])
            target_groups.setdefault(key, []).append(r)
            
        qc_rows = []
        for (pdb, variant, ph), rows in target_groups.items():
            n_rows = len(rows)
            n_fda = sum(1 for r in rows if r.get("library") == "FDA" or r.get("run_mode") == "fda")
            n_decoy = sum(1 for r in rows if r.get("is_decoy") == 1)
            n_controls = sum(1 for r in rows if r.get("is_control") == 1)
            n_valid = sum(1 for r in rows if r.get("pose_valid_any") == "1")
            
            t_vals = [_as_float(r.get("t_selected")) for r in rows if _as_float(r.get("t_selected")) is not None]
            mean_t = sum(t_vals) / len(t_vals) if t_vals else 0.0
            max_t = max(t_vals) if t_vals else 0.0
            
            # Use metrics from first row for that target
            ef1 = rows[0].get("ef1", "")
            roc_auc_adj = rows[0].get("roc_auc_adj", "")
            
            qc_rows.append({
                "pdb_id": pdb,
                "variant": variant,
                "ph_label": ph,
                "n_rows": n_rows,
                "n_fda": n_fda,
                "n_decoy": n_decoy,
                "n_controls": n_controls,
                "n_pose_valid": n_valid,
                "pose_valid_rate": f"{n_valid/n_rows:.3f}" if n_rows > 0 else "0.000",
                "mean_t_selected": f"{mean_t:.3f}",
                "max_t_selected": f"{max_t:.3f}",
                "ef1": ef1,
                "roc_auc_adj": roc_auc_adj
            })
            
        if qc_rows:
            qc_csv = data_dir / "target_qc_summary.csv"
            with qc_csv.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["pdb_id", "variant", "ph_label", "n_rows", "n_fda", "n_decoy", "n_controls", "n_pose_valid", "pose_valid_rate", "mean_t_selected", "max_t_selected", "ef1", "roc_auc_adj"])
                writer.writeheader()
                writer.writerows(qc_rows)
            logger.info("%s action=write status=ok path=%s rows=%d", COMPONENT, qc_csv, len(qc_rows))
    except Exception as e:
        logger.warning("%s action=write_optional status=failed path=target_qc_summary.csv error=%s", COMPONENT, e)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
