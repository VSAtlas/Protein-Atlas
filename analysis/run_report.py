# -*- coding: utf-8 -*-
import argparse
import csv
import datetime
import logging
import math
import sys
import yaml
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

COMPONENT = "[run-report]"

def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("run-report")

def _as_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        s = str(x).strip()
        if not s:
            return None
        v = float(s)
        if math.isnan(v):
            return None
        return v
    except Exception:
        return None

def _as_int(x: Any) -> int:
    try:
        return int(float(x))
    except Exception:
        return 0

def _pick_repeated_value(rows: List[Dict[str, Any]], col: str, is_int: bool = False) -> Optional[Any]:
    for r in rows:
        val = r.get(col)
        if val is None:
            continue
        s = str(val).strip()
        if not s or s.lower() == "nan":
            continue
        try:
            f = float(s)
            if math.isfinite(f):
                return int(f) if is_int else f
        except Exception:
            continue
    return None

def _load_scorch_stats(repo_root: Path, run_id: str, pdb_id: str, variant: str, ph: str) -> Dict[str, Optional[Any]]:
    stats = {"mu": None, "sigma": None, "n": None}
    combo_dir = repo_root / "post_docked" / run_id / pdb_id / variant / ph
    
    candidates = [
        combo_dir / "scorch_scores_all.csv",
        combo_dir / "dud_scorch_scores_all.csv"
    ]
    
    for path in candidates:
        if path.exists() and path.stat().st_size > 0:
            try:
                with path.open("r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    # Just need first row with data
                    for row in reader:
                        mu = _as_float(row.get("scorch_mu_decoy"))
                        sigma = _as_float(row.get("scorch_sigma_decoy"))
                        n = _as_int(row.get("scorch_n_decoys"))
                        if mu is not None or sigma is not None or n > 0:
                            stats["mu"] = mu
                            stats["sigma"] = sigma
                            stats["n"] = n if n > 0 else None
                            return stats
            except Exception:
                continue
    return stats

def _as_bool(x: Any) -> bool:
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(x)
    s = str(x).strip().lower()
    return s in ("1", "true", "yes", "on")


def _compute_fdr_stats(group_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    non_decoys = [r for r in group_rows if not _as_bool(r.get("is_decoy"))]
    score_field = None
    for r in non_decoys:
        val = (r.get("fdr_score_field") or "").strip()
        if val:
            score_field = val
            break
    if not score_field:
        return {}

    n_decoys = None
    for r in non_decoys:
        dec = r.get("fdr_n_decoys")
        if dec:
            try:
                n_decoys = int(float(str(dec)))
                break
            except Exception:
                continue

    n_tested = 0
    n_hits_q05 = 0
    n_hits_q10 = 0
    best_q = None

    for r in non_decoys:
        q_val = _as_float(r.get("fdr_q_target"))
        if q_val is None:
            continue
        n_tested += 1
        best_q = q_val if best_q is None else min(best_q, q_val)
        if q_val <= 0.05:
            n_hits_q05 += 1
        if q_val <= 0.10:
            n_hits_q10 += 1

    if n_tested == 0:
        return {}

    return {
        "score_field": score_field,
        "n_decoys": n_decoys,
        "n_tested": n_tested,
        "n_hits_q05": n_hits_q05,
        "n_hits_q10": n_hits_q10,
        "best_q": best_q
    }

def _find_manifest(run_id: str, repo_root: Path) -> Optional[Path]:
    candidates = [
        repo_root / "post_docked" / run_id / "run_manifest.yaml",
        repo_root / "docked" / run_id / "run_manifest.yaml",
        repo_root / "run_manifest.yaml",
    ]
    for path in candidates:
        if path.exists() and path.stat().st_size > 0:
            return path
    return None

def _load_manifest_pocket_map(manifest_path: Optional[Path]) -> Dict[Tuple[str, str, str], Dict[str, Any]]:
    # Returns (pdb_id, variant, ph) -> {method, center, box}
    pocket_map = {}
    if not manifest_path:
        return pocket_map
    
    try:
        with manifest_path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
    except Exception:
        return pocket_map

    proteins = data.get("proteins", {})
    if not isinstance(proteins, dict):
        return pocket_map

    for pdb_key, prot_data in proteins.items():
        if not isinstance(prot_data, dict):
            continue
        
        # Check top-level stages
        stages = prot_data.get("stages", {})
        details = stages.get("pocket_detection", {}).get("details", {})
        
        # Default info if found at top level
        default_info = {}
        if details:
            default_info = {
                "method": details.get("pocket_method"),
                "center": [details.get("center_x"), details.get("center_y"), details.get("center_z")],
                "box": [details.get("box_x"), details.get("box_y"), details.get("box_z")]
            }
            # Clean up Nones in lists
            if any(x is None for x in default_info["center"]):
                default_info["center"] = None
            if any(x is None for x in default_info["box"]):
                default_info["box"] = None

        # Iterate variants? If manifest structure supports it.
        # Often manifest is structured: proteins -> <PDB> -> variants -> <VARIANT> or similar?
        # The prompt says: match keys like "<PDB>|<VARIANT>|base" or "<PDB>|<VARIANT>|pH7_0"
        # Since we don't have perfect manifest iteration logic in spec, we will assume
        # the manifest matches the run structure if possible.
        #
        # However, for this implementation, we will store (pdb_id, *, *) -> default_info
        # and if specific variants exist, we'd use them.
        # Let's support simple lookups by PDB first.
        
        # We'll store normalized PDB key.
        pdb_norm = str(pdb_key).upper()
        pocket_map[(pdb_norm, "*", "*")] = default_info

        # If there are specific variants in manifest (not standard in current `run_manifest.yaml` structure shown previously, 
        # which tracks processing events linearly), we might miss them. 
        # We'll stick to the "default_info" found under stages->pocket_detection.

    return pocket_map

def _resolve_pocket(pocket_map: Dict, pdb: str, variant: str, ph: str) -> Dict[str, Any]:
    # Try exact match (not implemented fully above due to manifest ambiguity)
    # Fallback to (pdb, *, *)
    info = pocket_map.get((pdb.upper(), "*", "*"))
    if info:
        return info
    return {"method": None, "center": None, "box": None}

def build_report(run_id: str, repo_root: Path, top_n: int = 5) -> Dict[str, Any]:
    master_csv = repo_root / "data" / run_id / "master_rows.csv"
    if not master_csv.exists():
        raise FileNotFoundError(f"Master CSV not found: {master_csv}")

    manifest_path = _find_manifest(run_id, repo_root)
    pocket_map = _load_manifest_pocket_map(manifest_path)

    rows = []
    with master_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)

    # Group by target combo
    # Key: "<pdb_id>|<variant>|<ph_label>"
    target_groups: Dict[str, List[Dict[str, Any]]] = {}
    ligand_hits: Dict[str, List[Dict[str, Any]]] = {} # ligand_base -> list of hit info

    unique_ligands = set()
    
    for r in rows:
        pdb = r.get("pdb_id", "").strip()
        variant = r.get("variant", "").strip()
        ph = r.get("ph_label", "").strip()
        key = f"{pdb}|{variant}|{ph}"
        target_groups.setdefault(key, []).append(r)
        
        lig_base = r.get("ligand_base", "").strip()
        if lig_base:
            unique_ligands.add(lig_base)

    summary = {
        "n_target_combos": len(target_groups),
        "n_rows": len(rows),
        "n_unique_ligands": len(unique_ligands)
    }

    targets_out = {}
    
    sorted_keys = sorted(target_groups.keys())
    
    for key in sorted_keys:
        group_rows = target_groups[key]
        pdb, variant, ph = key.split("|")
        
        # QC Stats
        n_rows = len(group_rows)
        n_decoys = sum(1 for r in group_rows if _as_bool(r.get("is_decoy")))
        n_controls = sum(1 for r in group_rows if _as_bool(r.get("is_control")))
        
        # Extract metrics from first row (assuming they are repeated/consistent per target)
        first = group_rows[0]
        ef1 = _as_float(first.get("ef1"))
        roc_auc = _as_float(first.get("roc_auc"))
        roc_auc_adj = _as_float(first.get("roc_auc_adj"))
        dud_reason = first.get("dud_eval_status_reason") or None
        
        # Pocket info
        # Prefer row data if present (it was joined in master export)
        # master export columns: pocket_method, center_x, etc.
        p_method = first.get("pocket_method") or None
        if p_method:
            cx = _as_float(first.get("center_x"))
            cy = _as_float(first.get("center_y"))
            cz = _as_float(first.get("center_z"))
            bx = _as_float(first.get("box_x"))
            by = _as_float(first.get("box_y"))
            bz = _as_float(first.get("box_z"))
            center = [cx, cy, cz] if (cx is not None and cy is not None and cz is not None) else None
            box = [bx, by, bz] if (bx is not None and by is not None and bz is not None) else None
        else:
            # Fallback to manifest map (though master export should have handled this)
            p_info = _resolve_pocket(pocket_map, pdb, variant, ph)
            p_method = p_info.get("method")
            center = p_info.get("center")
            box = p_info.get("box")

        # Top N Ligands
        # Sort by t_selected desc, ligand_base asc
        def sort_key(x):
            t = _as_float(x.get("t_selected"))
            base = x.get("ligand_base", "")
            # t desc (None is lowest), base asc
            return (t if t is not None else -float('inf'), base) # wait, we want t desc.
            # Python sorts tuples element-wise.
            # We want t descending -> use -t? or reverse?
            # But we want base ascending as tie breaker.
            # If we reverse the whole thing: t ascending (bad), base descending (bad).
            # Let's use key properly.
            
        # To support (Desc, Asc):
        # key = (-t, base) works if t is float.
        # But we need to handle None for t.
        def rank_key(x):
            t = _as_float(x.get("t_selected"))
            if t is None:
                t_val = -float('inf')
            else:
                t_val = t
            return (-t_val, x.get("ligand_base", ""))

        sorted_rows = sorted(group_rows, key=rank_key)
        top_rows = sorted_rows[:top_n]
        
        top_list = []
        for r in top_rows:
            lig_base = r.get("ligand_base", "")
            t_sel = _as_float(r.get("t_selected"))
            t_src = r.get("t_selected_source") or None
            
            entry = {
                "ligand_base": lig_base,
                "t_selected": t_sel,
                "t_selected_source": t_src,
                "is_control": _as_bool(r.get("is_control")),
                "is_decoy": _as_bool(r.get("is_decoy"))
            }
            top_list.append(entry)
            
            # Register for ligand section
            if lig_base:
                ligand_hits.setdefault(lig_base, []).append({
                    "target": key,
                    "t_selected": t_sel,
                    "t_selected_source": t_src
                })

        fdr_stats = _compute_fdr_stats(group_rows)

        targets_out[key] = {
            "qc": {
                "ef1": ef1,
                "roc_auc": roc_auc,
                "roc_auc_adj": roc_auc_adj,
                "dud_eval_status_reason": dud_reason,
                "n_rows": n_rows,
                "n_decoys": n_decoys,
                "n_controls": n_controls,
                "decoy_stats": {
                    "consensus": {
                        "mu": _pick_repeated_value(group_rows, "consensus_mu_decoy"),
                        "sigma": _pick_repeated_value(group_rows, "consensus_sigma_decoy"),
                        "n": _pick_repeated_value(group_rows, "consensus_n_decoys", is_int=True)
                    },
                    "blend": {
                        "mu": _pick_repeated_value(group_rows, "blend_mu_decoy"),
                        "sigma": _pick_repeated_value(group_rows, "blend_sigma_decoy"),
                        "n": _pick_repeated_value(group_rows, "blend_n_decoys", is_int=True)
                    },
                    "scorch": _load_scorch_stats(repo_root, run_id, pdb, variant, ph)
                },
                "fdr": fdr_stats if fdr_stats else None
            },
            "pocket": {
                "method": p_method,
                "center": center,
                "box": box
            },
            "top5_ligands": top_list
        }

    # Ligands section
    ligands_out = {}
    for base in sorted(ligand_hits.keys()):
        hits = ligand_hits[base]
        # Sort hits by target key for determinism
        hits.sort(key=lambda x: x["target"])
        ligands_out[base] = {
            "targets_in_top5": hits
        }

    report = {
        "run_id": run_id,
        "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        "sources": {
            "master_rows_csv": str(master_csv.relative_to(repo_root)),
            "manifest_yaml": str(manifest_path.relative_to(repo_root)) if manifest_path else None
        },
        "summary": summary,
        "targets": targets_out,
        "ligands": ligands_out
    }
    
    return report

def write_yaml(report: Dict[str, Any], out_path: Path) -> None:
    # Use a compact representation for lists of numbers (center/box) if possible?
    # PyYAML default dump is okay.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        yaml.dump(report, f, sort_keys=False, default_flow_style=False)

def main() -> int:
    parser = argparse.ArgumentParser(description="Generate run report YAML")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    run_id = args.run_id
    
    out_path = repo_root / "data" / run_id / "report.yaml"
    if out_path.exists() and not args.overwrite:
        logger.info("%s action=skip reason=exists path=%s", COMPONENT, out_path)
        return 0

    try:
        report = build_report(run_id, repo_root)
        write_yaml(report, out_path)
        logger.info("%s action=write status=ok path=%s", COMPONENT, out_path)
    except Exception as e:
        logger.error("%s action=fail error=%s", COMPONENT, e)
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
