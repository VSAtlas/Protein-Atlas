# -*- coding: utf-8 -*-
import argparse
import csv
import datetime
import logging
import math
import sys
import yaml
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from analysis.manifest_utils import extract_pocket, load_run_manifest
except ModuleNotFoundError:
    REPO_ROOT = Path(__file__).resolve().parents[1]
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    SRC_ROOT = REPO_ROOT / "src"
    if str(SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(SRC_ROOT))
    from analysis.manifest_utils import extract_pocket, load_run_manifest

COMPONENT = "[run-report]"
DECOY_PREFIX_KEY = "DECOY_PREFIX"
DUD_PREFIX_KEY = "DUD_PREFIX"
DECOY_PREFIX_DEFAULT = "dud"
MIN_DECOYS_FOR_FDR = 200
MIN_UNIQUE_DECOY_SCORES = 10


def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("run-report")


def _strip_quotes(value: str) -> str:
    stripped = str(value or "").strip()
    if (
        len(stripped) >= 2
        and stripped[0] == stripped[-1]
        and stripped[0] in ("'", '"')
    ):
        return stripped[1:-1].strip()
    return stripped


def _read_decoy_prefix_from_file(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    decoy_value = None
    dud_value = None
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                key, raw_value = stripped.split("=", 1)
                key = key.strip()
                value = _strip_quotes(raw_value)
                if not value:
                    continue
                if key == DECOY_PREFIX_KEY:
                    decoy_value = value
                elif key == DUD_PREFIX_KEY:
                    dud_value = value
    except Exception:
        return None
    return decoy_value or dud_value


def _resolve_decoy_prefix(
    repo_root: Path, run_id: str, cli_value: Optional[str] = None
) -> str:
    if cli_value is not None:
        candidate = _strip_quotes(str(cli_value))
        if candidate:
            return candidate

    candidates = [
        repo_root / "data" / run_id / "config.txt",
        repo_root / "data" / run_id / "config_snapshot.txt",
        repo_root / "manifests" / run_id / "config.txt",
        repo_root / "config.txt",
    ]
    for path in candidates:
        value = _read_decoy_prefix_from_file(path)
        if value:
            return value
    return DECOY_PREFIX_DEFAULT


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


def _pick_repeated_value(
    rows: List[Dict[str, Any]], col: str, is_int: bool = False
) -> Optional[Any]:
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


def _load_scorch_stats(
    repo_root: Path,
    run_id: str,
    pdb_id: str,
    variant: str,
    ph: str,
    decoy_prefix: str = DECOY_PREFIX_DEFAULT,
) -> Dict[str, Optional[Any]]:
    stats = {"mu": None, "sigma": None, "n": None}
    combo_dir = repo_root / "post_docked" / run_id / pdb_id / variant / ph

    candidates = [
        combo_dir / f"{decoy_prefix}_scorch_scores_all.csv",
        combo_dir / "dud_scorch_scores_all.csv",
        combo_dir / "scorch_scores_all.csv",
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
    if not non_decoys:
        return {}

    score_field = None
    for r in non_decoys:
        val = (r.get("fdr_score_field") or "").strip()
        if val:
            score_field = val
            break

    n_decoys = None
    unique_decoys = None
    for r in non_decoys:
        dec = _as_float(r.get("fdr_n_decoys"))
        if dec is not None:
            n_decoys = int(dec)
            break
    for r in non_decoys:
        uniq = _as_float(r.get("fdr_unique_decoy_scores"))
        if uniq is not None:
            unique_decoys = int(uniq)
            break

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

    if best_q is None:
        best_q = 1.0

    reliable = (n_decoys is not None and n_decoys >= MIN_DECOYS_FOR_FDR) and (
        unique_decoys is not None and unique_decoys >= MIN_UNIQUE_DECOY_SCORES
    )
    return {
        "score_field": score_field,
        "n_decoys": n_decoys,
        "n_tested": n_tested,
        "n_hits_q05": n_hits_q05,
        "n_hits_q10": n_hits_q10,
        "best_q": best_q,
        "unique_decoy_scores": unique_decoys,
        "reliable": reliable,
    }


def _vector_from_row(row: Dict[str, Any], prefix: str) -> Optional[List[float]]:
    vals = [
        _as_float(row.get(f"{prefix}_x")),
        _as_float(row.get(f"{prefix}_y")),
        _as_float(row.get(f"{prefix}_z")),
    ]
    if all(v is not None and math.isfinite(v) for v in vals):
        return [float(v) for v in vals]  # type: ignore
    return None


def _resolve_manifest_pocket(
    manifest: Optional[Dict[str, Any]], pdb: str, variant: str, ph: str
) -> Dict[str, Any]:
    pocket = extract_pocket(manifest, pdb, variant, ph)
    center = pocket.get("center")
    box = pocket.get("box")

    normalized_center = None
    normalized_box = None
    if isinstance(center, (list, tuple)) and len(center) == 3:
        center_vals = [_as_float(v) for v in center]
        if all(v is not None for v in center_vals):
            normalized_center = [float(v) for v in center_vals]  # type: ignore

    if isinstance(box, (list, tuple)) and len(box) == 3:
        box_vals = [_as_float(v) for v in box]
        if all(v is not None for v in box_vals):
            normalized_box = [float(v) for v in box_vals]  # type: ignore

    return {
        "method": pocket.get("method"),
        "center": normalized_center,
        "box": normalized_box,
    }


def build_report(
    run_id: str,
    repo_root: Path,
    top_n: int = 5,
    decoy_prefix: str = DECOY_PREFIX_DEFAULT,
) -> Dict[str, Any]:
    master_csv = repo_root / "data" / run_id / "master_rows.csv"
    if not master_csv.exists():
        raise FileNotFoundError(f"Master CSV not found: {master_csv}")

    canonical_manifest_rel = Path("manifests") / run_id / "run_manifest.yaml"
    manifest_data, _manifest_path = load_run_manifest(repo_root, run_id)

    rows = []
    with master_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)

    # Group by target combo
    # Key: "<pdb_id>|<variant>|<ph_label>"
    target_groups: Dict[str, List[Dict[str, Any]]] = {}
    ligand_hits: Dict[str, List[Dict[str, Any]]] = {}  # ligand_base -> list of hit info

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
        "n_unique_ligands": len(unique_ligands),
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
        roc_auc_adj = _as_float(first.get("roc_auc_adj")) or roc_auc
        if ef1 is None:
            ef1 = 0.0
        if roc_auc is None:
            roc_auc = 0.0
        if roc_auc_adj is None:
            roc_auc_adj = roc_auc
        dud_reason = first.get("dud_eval_status_reason") or None

        # Pocket info
        # Prefer row data if present (it was joined in master export)
        # master export columns: pocket_method, center_x, etc.
        p_method = (first.get("pocket_method") or "").strip() or None
        center = _vector_from_row(first, "center")
        box = _vector_from_row(first, "box")

        manifest_pocket = _resolve_manifest_pocket(
            manifest_data, pdb.upper(), variant.upper(), ph
        )
        if not p_method:
            p_method = manifest_pocket.get("method")
        if center is None:
            center = manifest_pocket.get("center")
        if box is None:
            box = manifest_pocket.get("box")

        # Top N Ligands
        # Sort by t_selected desc, ligand_base asc
        def rank_key(x):
            t = _as_float(x.get("t_selected"))
            if t is None:
                t_val = -float("inf")
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
                "is_decoy": _as_bool(r.get("is_decoy")),
            }
            top_list.append(entry)

            # Register for ligand section
            if lig_base:
                ligand_hits.setdefault(lig_base, []).append(
                    {"target": key, "t_selected": t_sel, "t_selected_source": t_src}
                )

        fdr_stats = _compute_fdr_stats(group_rows)
        fdr_value: Optional[Any] = None
        if fdr_stats:
            best_q = _as_float(fdr_stats.get("best_q"))
            if best_q is None:
                best_q = 1.0
            reliable = bool(fdr_stats.get("reliable"))
            fdr_value = float(best_q) if reliable else f"SMALL {best_q:.6g}"

        targets_out[key] = {
            "qc": {
                "ef1": ef1,
                "roc_auc": roc_auc,
                "roc_auc_adj": roc_auc_adj,
                "dud_eval_status_reason": dud_reason,
                "n_rows": n_rows,
                "n_decoys": n_decoys,
                "n_controls": n_controls,
                "fdr": fdr_value,
                "decoy_stats": {
                    "consensus": {
                        "mu": _pick_repeated_value(group_rows, "consensus_mu_decoy"),
                        "sigma": _pick_repeated_value(
                            group_rows, "consensus_sigma_decoy"
                        ),
                        "n": _pick_repeated_value(
                            group_rows, "consensus_n_decoys", is_int=True
                        ),
                    },
                    "blend": {
                        "mu": _pick_repeated_value(group_rows, "blend_mu_decoy"),
                        "sigma": _pick_repeated_value(group_rows, "blend_sigma_decoy"),
                        "n": _pick_repeated_value(
                            group_rows, "blend_n_decoys", is_int=True
                        ),
                    },
                    "scorch": _load_scorch_stats(
                        repo_root,
                        run_id,
                        pdb,
                        variant,
                        ph,
                        decoy_prefix=decoy_prefix,
                    ),
                },
                "fdr_stats": fdr_stats if fdr_stats else None,
            },
            "pocket": {"method": p_method, "center": center, "box": box},
            "top5_ligands": top_list,
        }

    # Ligands section
    ligands_out = {}
    for base in sorted(ligand_hits.keys()):
        hits = ligand_hits[base]
        # Sort hits by target key for determinism
        hits.sort(key=lambda x: x["target"])
        ligands_out[base] = {"targets_in_top5": hits}

    report = {
        "run_id": run_id,
        "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        "sources": {
            "master_rows_csv": str(master_csv.relative_to(repo_root)),
            "manifest_yaml": str(canonical_manifest_rel),
        },
        "summary": summary,
        "targets": targets_out,
        "ligands": ligands_out,
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
    parser.add_argument("--decoy-prefix", default=None)
    args = parser.parse_args()

    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    run_id = args.run_id
    decoy_prefix = _resolve_decoy_prefix(repo_root, run_id, args.decoy_prefix)
    logger.info("%s action=preflight decoy_prefix=%s", COMPONENT, decoy_prefix)

    out_path = repo_root / "data" / run_id / "report.yaml"
    if out_path.exists() and not args.overwrite:
        logger.info("%s action=skip reason=exists path=%s", COMPONENT, out_path)
        return 0

    try:
        report = build_report(run_id, repo_root, decoy_prefix=decoy_prefix)
        write_yaml(report, out_path)
        logger.info("%s action=write status=ok path=%s", COMPONENT, out_path)
    except Exception as e:
        logger.error("%s action=fail error=%s", COMPONENT, e)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
