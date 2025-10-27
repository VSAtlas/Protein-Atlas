#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DUD/DUD-E style evaluator for Atlas docking outputs (filename-labeled actives/decoys).

Assumptions:
- Input: docked/<PDB_ID>/docking_score_long.csv (per-POSE rows)
- Lower score = better (Vina-style)
- Ligand filename contains token 'active' or 'decoy' (case-insensitive), e.g.
    ABL1_active_0123.pdbqt  or  ABL1_decoy_0456.pdbqt
- We derive labels from the filename; no chemdb/actives files are used.

Outputs:
- atlas/analysis/dud_eval/out/<PDB_ID>/
    - metrics.tsv
    - ef_curve.png, roc.png, pr.png, score_hist.png
- atlas/analysis/dud_eval/out/summary.tsv (per-target table)
- atlas/analysis/dud_eval/out/summary_macro.tsv (macro average over targets)

Metrics:
- EF@{1,2,5,10}%
- ROC-AUC, PR-AUC
- logAUC (lambda=1e-3), adjusted logAUC = logAUC - 0.14462
- BEDROC(alpha=20.0)

Notes:
- We "collapse" poses by taking the minimum score per ligand (best pose).
- We auto-detect score and ligand filename columns via heuristics; CLI flags override.
"""

import argparse
import math
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc

SCORE_CANDIDATES = [
    "score","docking_score","vina_score","affinity","pose_score",
    "dockscore","energy","gnina_score","cnnscore","cnn_affinity"
]
LIGFILE_CANDIDATES = [
    "ligand_file","ligand_path","ligand","ligand_name","ligand","pose_file",
    "pose_path","output_ligand","output_ligand_path","file","filepath","filename",
    "pdbqt_path","pdbqt","out_path"
]

def guess_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    lower = {c.lower(): c for c in df.columns}
    for key in candidates:
        if key in lower:
            return lower[key]
    for frag in candidates:
        for c in df.columns:
            if frag in c.lower():
                return c
    return None

def guess_score_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    col = guess_col(df, SCORE_CANDIDATES)
    if col:
        return col
    # last resort: first numeric column that smells like score
    for c in df.columns:
        if df[c].dtype.kind in "fi" and any(k in c.lower() for k in ("score","affin","dock","energy")):
            return c
    raise ValueError("Could not find score column. Use --score-col.")

def guess_ligfile_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    col = guess_col(df, LIGFILE_CANDIDATES)
    if col:
        return col
    # fallback: first object column
    for c in df.columns:
        if df[c].dtype == object:
            return c
    raise ValueError("Could not find ligand filename/path column. Use --lig-col.")

# --- filename parsing: derive ligand_id root + is_active from filename ---

TOKEN_RE = re.compile(r'(?<![A-Za-z0-9])(active|decoy)(?![A-Za-z0-9])', re.IGNORECASE)
POSE_TAIL_RE = re.compile(r'(?:_pose\d+|_mode\d+|_conf\d+|_rank\d+|_cluster\d+|_p\d+)$', re.IGNORECASE)
EXT_RE = re.compile(r'\.(pdbqt|sdf|mol2)(?:\.gz)?$', re.IGNORECASE)

def parse_name_and_label(path_str: str) -> Tuple[str, Optional[int]]:
    """Return (ligand_root_id, is_active) from a path or filename.
       is_active: 1 for active, 0 for decoy, None if not found."""
    base = os.path.basename(str(path_str))
    stem = EXT_RE.sub("", base)               # drop .pdbqt/.sdf/.mol2(.gz)
    stem = POSE_TAIL_RE.sub("", stem)         # drop trailing _pose1 etc.
    # Keep the full stem as the ligand identifier (unique per ligand)
    lig_id = stem
    m = TOKEN_RE.search(base)
    if not m:
        return lig_id, None
    tok = m.group(1).lower()
    return lig_id, 1 if tok == "active" else 0

# --- metrics ---

def ef_at_fractions(y_true: np.ndarray, scores_low_is_better: np.ndarray,
                    fractions=(0.01,0.02,0.05,0.10)) -> Dict[str,float]:
    N = len(y_true)
    n_act = int(y_true.sum())
    out = {}
    if N == 0 or n_act == 0:
        return {f"EF@{int(fr*100)}%": float("nan") for fr in fractions}
    order = np.argsort(scores_low_is_better)    # lowest score first
    y_sorted = y_true[order]
    cum_act = np.cumsum(y_sorted)
    for fr in fractions:
        k = max(1, int(round(fr * N)))
        found = int(cum_act[k-1])
        hit_rate = found / k
        base_rate = n_act / N
        out[f"EF@{int(fr*100)}%"] = float(hit_rate / base_rate) if base_rate > 0 else float("nan")
    return out

def pr_auc(y_true: np.ndarray, y_score_high_is_better: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    precision, recall, _ = precision_recall_curve(y_true, y_score_high_is_better)
    return float(auc(recall, precision)), precision, recall

def log_auc_from_roc(fpr: np.ndarray, tpr: np.ndarray, lam: float=1e-3) -> float:
    order = np.argsort(fpr)
    fpr = fpr[order]; tpr = tpr[order]
    mask = fpr >= lam
    if not np.any(mask):
        return 0.0
    # insert point at lam if needed
    if not np.isclose(fpr[mask][0], lam):
        i = np.searchsorted(fpr, lam)
        x0,x1 = fpr[i-1], fpr[i]; y0,y1 = tpr[i-1], tpr[i]
        ylam = y0 + (y1 - y0) * (lam - x0) / (x1 - x0)
        fpr = np.insert(fpr, i, lam)
        tpr = np.insert(tpr, i, ylam)
        mask = fpr >= lam
    xf = fpr[mask]; yf = tpr[mask]
    logx = np.log10(xf)
    num = np.sum((logx[1:] - logx[:-1]) * (yf[1:] + yf[:-1]) / 2.0)
    denom = math.log10(1.0/lam)
    return float(num/denom) if denom > 0 else float("nan")

def bedroc(y_true: np.ndarray, y_score_high_is_better: np.ndarray, alpha: float=20.0) -> float:
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score_high_is_better, dtype=float)
    N = len(y_true); n = int(y_true.sum())
    if N == 0 or n == 0 or n == N:
        return float("nan")
    order = np.argsort(-y_score)
    ranks = np.nonzero(y_true[order]==1)[0]  # 0-based ranks among sorted list
    s = float(np.sum(np.exp(-alpha * ranks / N)))
    ra = n / N
    # constants per Truchon & Bayly (2007), matching common implementations
    k1 = (ra * (1.0 - math.exp(-alpha))) / (math.exp(alpha / N) - 1.0)
    k2 = (ra * math.sinh(alpha/2.0)) / (math.cosh(alpha/2.0) - math.cosh(alpha/2.0 - alpha*ra))
    return float((s / k1) * k2)

# --- evaluation ---

def evaluate_target(pdb_id: str,
                    csv_path: Path,
                    out_dir: Path,
                    lig_col_cli: Optional[str],
                    score_col_cli: Optional[str],
                    bedroc_alpha: float,
                    logauc_lambda: float) -> Optional[pd.Series]:
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"[WARN] {pdb_id}: failed to read {csv_path}: {e}")
        return None

    lig_col = guess_ligfile_col(df, lig_col_cli)
    score_col = guess_score_col(df, score_col_cli)
    # Coerce scores to numeric; drop NaN/±inf early to avoid NaNs in metrics
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    before_nf = len(df)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna(subset=[score_col])
    if len(df) < before_nf:
        print(f"[INFO] {pdb_id}: dropped {before_nf - len(df)} rows with non-finite scores")

    # derive ligand_id + label from filename
    parsed = df[lig_col].astype(str).apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])

    # drop rows where label couldn't be inferred
    before = len(df)
    df = df.dropna(subset=["is_active"])
    if len(df) < before:
        print(f"[INFO] {pdb_id}: dropped {before - len(df)} rows without 'active/decoy' token in filename")

    # collapse to best score per ligand (min score)
    best = df.groupby("lig_id", as_index=False).agg(
        best_score=(score_col, "min"),
        is_active=("is_active", "max"),   # any active -> active
    )
    # Some ligands may still have all-NaN scores → best_score=NaN; drop them
    before_best = len(best)
    best = best.dropna(subset=["best_score"])
    if len(best) < before_best:
        print(f"[INFO] {pdb_id}: dropped {before_best - len(best)} ligands with NaN best_score")

    y_true = best["is_active"].astype(int).to_numpy()
    y_low = best["best_score"].to_numpy()
    y_high = -y_low

    N = len(best); n_act = int(y_true.sum())
    if N == 0 or n_act == 0 or n_act == N:
        print(f"[WARN] {pdb_id}: degenerate set (N={N}, actives={n_act}). Skipping.")
        return None

    # metrics
    ef = ef_at_fractions(y_true, y_low, fractions=(0.01,0.02,0.05,0.10))
    rocAUC = float(roc_auc_score(y_true, y_high))
    fpr, tpr, _ = roc_curve(y_true, y_high)
    prAUC, precision, recall = pr_auc(y_true, y_high)
    lAUC = log_auc_from_roc(fpr, tpr, lam=logauc_lambda)
    lAUC_adj = lAUC - 0.14462
    bed = bedroc(y_true, y_high, alpha=bedroc_alpha)

    # plots
    out_dir.mkdir(parents=True, exist_ok=True)

    # EF curve
    order = np.argsort(y_low)
    y_sorted = y_true[order]
    cum_pos = np.cumsum(y_sorted)
    k = np.arange(1, N+1)
    ef_curve = (N / max(n_act,1)) * (cum_pos / k)
    frac = k / N
    plt.figure()
    plt.plot(frac, ef_curve)
    plt.xlabel("Fraction screened")
    plt.ylabel("Enrichment factor (EF)")
    plt.title(f"{pdb_id} - EF curve")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "ef_curve.png", dpi=200)
    plt.close()

    # ROC
    plt.figure()
    plt.plot(fpr, tpr, label=f"AUC={rocAUC:.3f}")
    plt.plot([0,1],[0,1],"k--",lw=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title(f"{pdb_id} - ROC")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "roc.png", dpi=200)
    plt.close()

    # PR
    base = n_act / N
    plt.figure()
    plt.plot(recall, precision, label=f"PR-AUC={prAUC:.3f}")
    plt.hlines(base, 0, 1, colors="k", linestyles="--", linewidth=1, label=f"Baseline={base:.3f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"{pdb_id} - Precision–Recall")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "pr.png", dpi=200)
    plt.close()

    # Score distributions
    plt.figure()
    plt.hist(best.loc[best["is_active"]==1, "best_score"], bins=40, alpha=0.6, label="Actives")
    plt.hist(best.loc[best["is_active"]==0, "best_score"], bins=40, alpha=0.6, label="Decoys")
    plt.xlabel("Best docking score (lower = better)")
    plt.ylabel("Count")
    plt.title(f"{pdb_id} - Score distributions")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "score_hist.png", dpi=200)
    plt.close()

    row = {
        "pdb_id": pdb_id,
        "N": N,
        "n_actives": n_act,
        "actives_fraction": base,
        "ROC_AUC": rocAUC,
        "PR_AUC": prAUC,
        "logAUC": lAUC,
        "logAUC_adj": lAUC_adj,
        f"BEDROC_alpha_{bedroc_alpha:g}": bed,
        **ef,
    }
    pd.DataFrame([row]).to_csv(out_dir / "metrics.tsv", sep="\t", index=False)
    return pd.Series(row)

def main():
    ap = argparse.ArgumentParser(description="Atlas VS benchmark evaluator (filename-labeled actives/decoys).")
    ap.add_argument("--docked-root", type=str, default="docked",
                    help="Root folder (or a single target folder) to scan for docking_score_long.csv.")
    ap.add_argument("--out-dir", type=str, default="atlas/analysis/dud_eval/out",
                    help="Output root directory.")
    ap.add_argument("--lig-col", type=str, default=None,
                    help="Column containing ligand filename/path (auto-detected if omitted).")
    ap.add_argument("--score-col", type=str, default=None,
                    help="Score column (lower is better). Auto-detected if omitted.")
    ap.add_argument("--bedroc-alpha", type=float, default=20.0)
    ap.add_argument("--logauc-lambda", type=float, default=1e-3)
    args = ap.parse_args()

    docked_root = Path(args.docked_root)
    out_root = Path(args.out_dir)

    # discover targets
    targets: List[Tuple[str, Path]] = []
    if (docked_root / "docking_score_long.csv").exists():
        targets.append((docked_root.name, docked_root / "docking_score_long.csv"))
    else:
        for sub in sorted(docked_root.iterdir()):
            if not sub.is_dir():
                continue
            csvp = sub / "docking_score_long.csv"
            if csvp.exists():
                targets.append((sub.name, csvp))

    if not targets:
        print(f"[ERR] No docking_score_long.csv under {docked_root}")
        raise SystemExit(2)

    rows = []
    for pdb_id, csvp in targets:
        series = evaluate_target(
            pdb_id=pdb_id,
            csv_path=csvp,
            out_dir=out_root / pdb_id,
            lig_col_cli=args.lig_col,
            score_col_cli=args.score_col,
            bedroc_alpha=args.bedroc_alpha,
            logauc_lambda=args.logauc_lambda,
        )
        if series is not None:
            rows.append(series)

    if not rows:
        print("[ERR] Nothing evaluated.")
        raise SystemExit(3)

    df = pd.DataFrame(rows).sort_values("pdb_id")
    out_root.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_root / "summary.tsv", sep="\t", index=False)

    metric_cols = [c for c in df.columns if c not in ("pdb_id","N","n_actives","actives_fraction")]
    macro = df[metric_cols].mean(numeric_only=True).to_dict()
    pd.DataFrame([{"pdb_id":"macro_avg", **{k: macro[k] for k in metric_cols}}]) \
      .to_csv(out_root / "summary_macro.tsv", sep="\t", index=False)

    print(f"[OK] Wrote per-target results to: {out_root}")
    print(f"[OK] Summary: {out_root / 'summary.tsv'}")

if __name__ == "__main__":
    main()
