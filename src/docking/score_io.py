from __future__ import annotations

import csv
import math
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

# Robust AutoDock Vina REMARK scanning (spacing/case); shared by run pipeline and parsers.
VINA_RESULT_SCORE_RE = re.compile(
    r"REMARK\s+VINA\s+RESULT[:\s]+(-?\d+(?:\.\d+)?)", re.IGNORECASE
)
VINA_MODEL_START_RE = re.compile(r"^\s*MODEL\b", re.IGNORECASE)
LIGAND_STAGE_FILENAME_RE = re.compile(r"^(.*?)(_stage\d+)?\.pdbqt$", re.IGNORECASE)


def _fold_vina_remark_into_best(line: str, best: Optional[float]) -> Optional[float]:
    """If ``line`` matches a Vina REMARK score, return the lower of ``best`` and that score."""
    m = VINA_RESULT_SCORE_RE.search(line)
    if not m:
        return best
    try:
        e = float(m.group(1))
    except ValueError:
        return best
    if best is None or e < best:
        return e
    return best


def extract_best_vina_score_from_lines(lines: Iterable[str]) -> Optional[float]:
    """Minimum Vina affinity from in-memory PDBQT ``lines`` (one MODEL block)."""
    best: Optional[float] = None
    for ln in lines:
        best = _fold_vina_remark_into_best(ln, best)
    return best


def extract_best_vina_score(pdbqt_path: str | os.PathLike[str]) -> tuple[Optional[float], int]:
    """
    Stream PDBQT for the lowest Vina affinity.

    Returns ``(best_score, n_models)`` for logging; ``n_models`` assumes 1 when no
    ``MODEL`` records are present (single logical pose file).
    """
    p = Path(pdbqt_path)
    if not p.exists() or p.stat().st_size < 32:
        return None, 0
    best: Optional[float] = None
    n_models = 0
    saw_model_tag = False
    try:
        with open(p, "r", encoding="utf-8", errors="ignore") as fh:
            buf: list[str] = []
            for ln in fh:
                if VINA_MODEL_START_RE.match(ln):
                    if buf:
                        n_models += 1
                    buf = [ln]
                    saw_model_tag = True
                else:
                    buf.append(ln)
                best = _fold_vina_remark_into_best(ln, best)
            if saw_model_tag:
                n_models += 1
    except OSError:
        return None, n_models
    if n_models == 0:
        n_models = 1
    return best, n_models

def write_score_summary_to_csv(
    score_history, output_path="docking_score_summary.csv", run_id=None, variant=None
):
    ligand_scores = defaultdict(dict)

    for stage, stage_scores in score_history.items():
        for path, score in stage_scores.items():
            ligand_filename = os.path.basename(path)
            match = LIGAND_STAGE_FILENAME_RE.match(ligand_filename)
            if match:
                ligand_core = f"{match.group(1)}.pdbqt"
                ligand_scores[ligand_core][stage] = score

    all_ligands = sorted(ligand_scores.keys())
    all_stages = sorted(score_history.keys())

    run_id_value = "" if run_id is None else str(run_id)
    variant_value = (variant or "").strip()
    include_variant = bool(variant is not None and variant_value)

    header = ["run_id"]
    if include_variant:
        header.append("variant")
    header.append("Ligand")
    header.extend(all_stages)
    rows = []
    for ligand in all_ligands:
        row = [run_id_value]
        if include_variant:
            row.append(variant_value)
        row.append(ligand)
        for stage in all_stages:
            score = ligand_scores[ligand].get(stage, "")
            if isinstance(score, (float, int)):
                row.append(f"{score:.2f}")
            else:
                row.append(str(score) if score else "")
        rows.append(row)

    with open(output_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerows(rows)

    print(f"\nScore summary written to: {output_path}")


def extract_best_score(docked_pdbqt_path: str | os.PathLike[str]) -> Optional[float]:
    best, _n = extract_best_vina_score(docked_pdbqt_path)
    return best


def extract_gnina_scores(docked_pdbqt_path: str) -> Dict[str, Optional[float]]:
    """
    Parse GNINA REMARK lines from a docked PDBQT.

    Returns keys:
      - minimized_affinity_kcal
      - cnn_score
      - cnn_affinity_pK
    Values are float or None if not present/parseable.
    """
    metrics: Dict[str, Optional[float]] = {
        "minimized_affinity_kcal": None,
        "cnn_score": None,
        "cnn_affinity_pK": None,
    }
    try:
        with open(docked_pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.startswith("REMARK"):
                    continue
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                token = parts[1].lower()
                try:
                    val = float(parts[-1])
                except Exception:
                    val = None
                if token == "minimizedaffinity":
                    metrics["minimized_affinity_kcal"] = val
                elif token == "cnnscore":
                    metrics["cnn_score"] = val
                elif token == "cnnaffinity":
                    metrics["cnn_affinity_pK"] = val
    except Exception:
        pass
    return metrics


# -------------------------
# Vina config writer, uses shim(lazy yeah)
# -------------------------
def generate_config(
    output_dir,
    pdb_id,
    receptor_pdbqt,
    center,
    box_size,
    ligand_path,
    stage,
    stage_info,
    cpu_per_job,
    docked_dir=None,
):
    # Legacy shim: derive a minimal cfg for older callers
    from docking.docking_vina import emit_vina_config  # type: ignore[import-not-found]

    cfg = {
        "OVERALL_DIR": output_dir,
        "CONFIGS_DIR": os.path.join(output_dir, "configs"),
        "DOCKED_DIR": docked_dir or os.path.join(output_dir, "docked"),
        "RUN_ID": "legacy",
        "RESET_CONFIGS": False,
    }
    cfg["CONFIG_RUN_DIR"] = os.path.join(cfg["CONFIGS_DIR"], cfg["RUN_ID"])
    Path(cfg["CONFIG_RUN_DIR"]).mkdir(parents=True, exist_ok=True)
    if not str(ligand_path).lower().endswith(".pdbqt"):
        raise ValueError(f"Ligand is not a .pdbqt file: {ligand_path}")
    return emit_vina_config(
        cfg,
        pdb_id,
        receptor_pdbqt,
        center,
        box_size,
        ligand_path,
        stage,
        stage_info,
        cpu_per_job,
        logger=None,
        legacy=True,
    )


# -------------------------
# Score bookkeeping
# -------------------------
def record_score(score_history, stage_name, ligand, score: Any, valid, reason=None):
    def coerce_score(x):
        if isinstance(x, (int, float)):
            return float(x)
        if isinstance(x, str) and x.lower().endswith(".pdbqt"):
            try:
                return extract_best_score(x)
            except Exception:
                return None
        try:
            return float(x)
        except Exception:
            return None

    s = coerce_score(score)
    score_history[stage_name][ligand] = {
        "score": s,
        "valid": bool(valid),
        "reason": reason,
    }


def score_key(item):
    """item = (ligand, rec). Sort by numeric score; invalid or None go to bottom."""
    _, rec = item
    s = rec.get("score")
    if s is None:
        return math.inf
    try:
        s = float(s)
    except (TypeError, ValueError):
        return math.inf
    return s if rec.get("valid", False) else s + 1e-6

