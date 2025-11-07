#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DUD/DUD-E style evaluator for Atlas docking outputs (filename-labeled actives/decoys).

Assumptions:
- Input: docked/<PDB_ID>/docking_score_long.csv (per-POSE rows)
- Lower score = better (Vina-style)
- Ligand filename contains token 'active/actives' or 'decoy/decoys' (case-insensitive)
    ABL1_active_0123.pdbqt  or  ABL1_decoy_0456.pdbqt
- We derive labels from the filename; no chemdb/actives files are used.

Outputs:
- <OVERALL_DIR>/analysis/out/<PDB_ID>/
    - metrics.tsv
    - ef_curve.png, roc.png, pr.png, score_hist.png
- <OVERALL_DIR>/analysis/out/summary.tsv (per-target table)
- <OVERALL_DIR>/analysis/out/summary_macro.tsv (macro average over targets)

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
import json
import math
import os
import re
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc

# >>> PATHS IMPORT START
from pathlib import Path
# --- ensure repo root is importable when running from analysis/ ---
import sys
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from path_router import make_paths, expand_variants


def _dedup(seq: Iterable[str]) -> List[str]:
    """Preserve order while removing duplicates."""
    seen: Set[str] = set()
    out: List[str] = []
    for item in seq:
        if item and item not in seen:
            seen.add(item)
            out.append(item)
    return out


def extract_compnd_molecules(pdb_lines: List[str]) -> List[str]:
    """Pull COMPND.MOLECULE entries from a PDB header."""
    molecules: List[str] = []
    pending: Optional[str] = None
    last_key: Optional[str] = None

    def _flush() -> None:
        nonlocal pending
        if pending:
            molecules.append(pending.strip().rstrip(","))
            pending = None

    for line in pdb_lines:
        if not line.startswith("COMPND"):
            continue
        payload = line[10:].strip()
        if not payload:
            continue
        segments = [seg.strip() for seg in payload.split(";")]
        for seg in segments:
            if not seg:
                continue
            if ":" in seg:
                key, val = seg.split(":", 1)
                key = key.strip().upper()
                val = val.strip()
                if key != last_key and last_key == "MOLECULE":
                    _flush()
                last_key = key
                if key == "MOLECULE":
                    pending = f"{pending} {val}".strip() if pending else val
                elif pending and key != "MOLECULE":
                    _flush()
            else:
                if last_key == "MOLECULE":
                    pending = f"{pending} {seg}".strip() if pending else seg
    if last_key != "MOLECULE":
        _flush()
    _flush()
    return _dedup(molecules)


def extract_uniprot_from_dbref(pdb_lines: List[str]) -> Tuple[List[str], List[str]]:
    """Return (entry_names, accessions) discovered from DBREF UNP rows."""
    entry_names: List[str] = []
    accessions: List[str] = []
    for line in pdb_lines:
        if not line.startswith("DBREF"):
            continue
        tokens = line.split()
        if len(tokens) < 7:
            continue
        db = tokens[5].upper()
        if db not in {"UNP", "UNIPROT"}:
            continue
        accession = tokens[6].strip()
        entry = tokens[7].strip() if len(tokens) >= 8 else ""
        if accession:
            accessions.append(accession)
        if entry:
            entry_names.append(entry)
    return _dedup(entry_names), _dedup(accessions)


def choose_target_name(compnd_mols: List[str],
                       uniprot_entries: List[str],
                       uniprot_accessions: List[str],
                       prefer: str = "auto") -> str:
    """Pick the best target label respecting preference order."""
    prefer_key = (prefer or "auto").lower()
    compnd_choice = ", ".join(compnd_mols) if compnd_mols else ""
    uniprot_entry = uniprot_entries[0] if uniprot_entries else ""
    uniprot_acc = uniprot_accessions[0] if uniprot_accessions else ""

    if prefer_key == "compnd":
        return compnd_choice or uniprot_entry or uniprot_acc
    if prefer_key == "uniprot":
        return uniprot_entry or uniprot_acc or compnd_choice
    # auto
    return compnd_choice or uniprot_entry or uniprot_acc


def _read_pdb_header_lines(pdb_path: Path) -> List[str]:
    lines: List[str] = []
    try:
        with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")):
                    break
                lines.append(line.rstrip("\n"))
    except Exception:
        return []
    return lines


def infer_library_name(target_id: str,
                       ligand_basenames: Set[str],
                       *,
                       prepped_root_override: Optional[Path],
                       cfg: Dict) -> str:
    """Infer library folder name via manifest or filename overlap."""
    candidate_roots: List[Path] = []
    if prepped_root_override:
        candidate_roots.append(prepped_root_override)
    if cfg:
        try:
            paths = make_paths(cfg, base_id=target_id, pdb_file=f"{target_id}.pdb")
            candidate_roots.append(paths.prepped_root)
        except Exception:
            pass

    seen: Set[Path] = set()
    roots: List[Path] = []
    for root in candidate_roots:
        if root and root not in seen:
            seen.add(root)
            roots.append(root)

    def _manifest_label(manifest_path: Path) -> Optional[str]:
        try:
            with manifest_path.open("r", encoding="utf-8") as handle:
                data = json.load(handle)
        except Exception:
            return None
        for key in ("library_name", "name", "library", "label"):
            val = data.get(key)
            if isinstance(val, str) and val.strip():
                return val.strip()
        return None

    def _collect_basenames(folder: Path) -> Set[str]:
        allowed_suffixes = {".pdbqt", ".sdf", ".mol2"}
        names: Set[str] = set()
        for entry in folder.iterdir():
            if entry.is_file() and entry.suffix.lower() in allowed_suffixes:
                names.add(entry.name)
        return names

    for root in roots:
        target_dir = root / target_id
        if not target_dir.exists():
            continue
        best_name: Optional[str] = None
        best_score = -1
        tied = False
        for library_dir in sorted(target_dir.iterdir()):
            if not library_dir.is_dir():
                continue
            manifest = library_dir / "manifest.json"
            if manifest.exists():
                label = _manifest_label(manifest)
                if label:
                    return label
            basenames = _collect_basenames(library_dir)
            if not basenames:
                continue
            score = len(ligand_basenames & basenames) if ligand_basenames else 0
            if score > best_score:
                best_score = score
                best_name = library_dir.name
                tied = False
            elif score == best_score and score >= 0:
                tied = True
        if best_name and best_score > 0 and not tied:
            return best_name
        if best_name and tied:
            print(f"[WARN] {target_id}: multiple libraries tied for filename overlap; leaving library_name blank")
    return ""


def derive_target_name(target_id: str,
                       *,
                       prefer: str,
                       pdb_root_override: Optional[Path],
                       cfg: Dict) -> str:
    candidates: List[Path] = []
    seen: Set[Path] = set()

    def _append_candidate(path: Optional[Path]) -> None:
        if path and path.exists() and path not in seen:
            seen.add(path)
            candidates.append(path)

    def _collect_from_dir(base: Path) -> None:
        if not base.exists():
            return
        for pattern in ("*.pdb", "*/*.pdb"):
            for found in sorted(base.glob(pattern)):
                _append_candidate(found)

    if pdb_root_override:
        base = pdb_root_override / target_id
        _collect_from_dir(base)

    if cfg:
        try:
            paths = make_paths(cfg, base_id=target_id, pdb_file=f"{target_id}.pdb")
        except Exception:
            paths = None
        if paths is not None:
            _collect_from_dir(paths.root_pdb_dir)
            _append_candidate(paths.input_pdb_path)

    if not candidates:
        return ""

    header_lines = _read_pdb_header_lines(candidates[0])
    if not header_lines:
        return ""
    compnd = extract_compnd_molecules(header_lines)
    entries, accessions = extract_uniprot_from_dbref(header_lines)
    return choose_target_name(compnd, entries, accessions, prefer=prefer)
# >>> PATHS IMPORT END


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
TOKEN_RE = re.compile(r'(?<![A-Za-z0-9])(active|decoy)s?(?![A-Za-z0-9])', re.IGNORECASE)
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
                    logauc_lambda: float,
                    run_id: Optional[str]) -> Optional[Tuple[pd.Series, Set[str]]]:
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"[WARN] {pdb_id}: failed to read {csv_path}: {e}")
        return None

    if run_id:
        total_rows = len(df)
        if "run_id" in df.columns:
            mask = df["run_id"].astype(str) == str(run_id)
            df = df.loc[mask].copy()
            print(f"[eval] Using run_id={run_id} (rows kept: {len(df)}/{total_rows})")
        else:
            print("[eval] --run-id specified but CSV has no 'run_id'; proceeding unfiltered.")

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
    lig_paths = df[lig_col].astype(str)
    parsed = lig_paths.apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])
    used_ligand_basenames: Set[str] = {os.path.basename(p) for p in lig_paths.tolist() if p}

    # drop rows where label couldn't be inferred
    before = len(df)
    df = df.dropna(subset=["is_active"])
    if len(df) < before:
        print(f"[INFO] {pdb_id}: dropped {before - len(df)} rows without 'active(s)/decoy(s)' token in filename")

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

    variant_label: Optional[str] = None
    if "variant" in df.columns:
        variant_series = df["variant"].dropna()
        if not variant_series.empty:
            normalized = variant_series.astype(str).str.strip()
            normalized = normalized[normalized != ""]
            if not normalized.empty:
                upper = normalized.str.upper()
                has_apo = bool((upper == "APO").any())
                has_holo = bool((upper == "HOLO").any())
                if has_apo and not has_holo:
                    variant_label = "APO"
                elif has_holo and not has_apo:
                    variant_label = "HOLO"
                elif has_apo or has_holo:
                    variant_label = "MIXED"

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
    if variant_label:
        row["variant"] = variant_label
    pd.DataFrame([row]).to_csv(out_dir / "metrics.tsv", sep="\t", index=False)
    return pd.Series(row), used_ligand_basenames

def _load_default_cfg() -> Dict:
    try:
        from input_and_export_functions import load_config, validate_config
    except Exception:
        return {}

    try:
        cfg = load_config("config.txt") or {}
        if cfg:
            validate_config(cfg)
        return cfg
    except Exception:
        return {}

def main():
    ap = argparse.ArgumentParser(description="Atlas VS benchmark evaluator (filename-labeled actives/decoys).")
    ap.add_argument("--docked-root", type=str, default="docked",
                    help="Root folder (or a single target folder) to scan for docking_score_long.csv.")
    ap.add_argument("--out-dir", type=str, default="atlas/analysis/out",
                    help="Output root directory.")
    ap.add_argument("--lig-col", type=str, default=None,
                    help="Column containing ligand filename/path (auto-detected if omitted).")
    ap.add_argument("--score-col", type=str, default=None,
                    help="Score column (lower is better). Auto-detected if omitted.")
    ap.add_argument("--run-id", type=str, default=None,
                    help="Filter docking_score_long.csv rows to a specific run identifier.")
    ap.add_argument("--bedroc-alpha", type=float, default=20.0)
    ap.add_argument("--logauc-lambda", type=float, default=1e-3)
    ap.add_argument("--target-name-from-pdb", action="store_true",
                    help="If set, add a target_name column derived from PDB headers.")
    ap.add_argument("--target-name-prefer", type=str, default="auto",
                    choices=("auto", "compnd", "uniprot"),
                    help="Preference order when selecting target_name (default: auto).")
    ap.add_argument("--pdb-root", type=str, default=None,
                    help="Optional override root for processed PDB folders (processed_pdbs/<target>/).")
    ap.add_argument("--report-library", action="store_true",
                    help="If set, attempt to infer library_name from prepped_ligands.")
    ap.add_argument("--prepped-root", type=str, default=None,
                    help="Optional override root for prepped ligands (prepped_ligands/<target>/).")
    args = ap.parse_args()

    cfg = _load_default_cfg()
    docked_root = Path(args.docked_root)
    out_root = Path(args.out_dir)
    pdb_root_override = Path(args.pdb_root) if args.pdb_root else None
    prepped_root_override = Path(args.prepped_root) if args.prepped_root else None

    analysis_root = out_root
    if cfg and "OVERALL_DIR" in cfg:
        # >>> ANALYSIS OUT ROOT PATCH START
        analysis_root = Path(cfg["OVERALL_DIR"]) / "analysis" / "out"
        # >>> ANALYSIS OUT ROOT PATCH END
    analysis_root.mkdir(parents=True, exist_ok=True)

    def _resolve_docking_csv(pdb_id: str, fallback_dir: Path) -> Optional[Path]:
        if not cfg:
            candidate = fallback_dir / "docking_score_long.csv"
            return candidate if candidate.exists() else None

        try:
            # >>> PATHS INIT START
            mode    = str(cfg.get("APO_HOLO_MODE", "")).strip()
            variants = expand_variants(mode)  # returns [None] | ["APO","HOLO"]
            paths   = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
            ph_token = cfg.get("PH_TOKEN") or None
            # >>> PATHS INIT END
        except Exception:
            candidate = fallback_dir / "docking_score_long.csv"
            return candidate if candidate.exists() else None

        # >>> DOCKED INPUT PATHS PATCH START
        docked_root_cfg = paths.docked_pdb_root()
        summary_csv = paths.docking_score_summary_csv()
        long_csv    = paths.docking_score_long_csv()
        # >>> DOCKED INPUT PATHS PATCH END

        if not docked_root_cfg.exists():
            return None

        if long_csv.exists():
            return long_csv

        for variant in variants:
            variant_root = paths.docked_variant_root(variant)
            candidate = variant_root / "docking_score_long.csv"
            if candidate.exists():
                return candidate

        candidate = fallback_dir / "docking_score_long.csv"
        if candidate.exists():
            return candidate

        # Touch summary path to exercise router (no fallback to summary file for eval)
        summary_csv.exists()

        return None

    # discover targets
    targets: List[Tuple[str, Path]] = []

    if (docked_root / "docking_score_long.csv").exists():
        csv_path = _resolve_docking_csv(docked_root.name, docked_root)
        if csv_path is not None:
            targets.append((docked_root.name, csv_path))
    elif docked_root.is_dir():
        for sub in sorted(docked_root.iterdir()):
            if not sub.is_dir():
                continue
            csv_path = _resolve_docking_csv(sub.name, sub)
            if csv_path is not None:
                targets.append((sub.name, csv_path))

    if not targets:
        print(f"[ERR] No docking_score_long.csv under {docked_root}")
        raise SystemExit(2)

    rows: List[pd.Series] = []
    ligand_basenames_by_target: Dict[str, Set[str]] = {}
    for pdb_id, csvp in targets:
        evaluated = evaluate_target(
            pdb_id=pdb_id,
            csv_path=csvp,
            out_dir=analysis_root / pdb_id,
            lig_col_cli=args.lig_col,
            score_col_cli=args.score_col,
            bedroc_alpha=args.bedroc_alpha,
            logauc_lambda=args.logauc_lambda,
            run_id=args.run_id,
        )
        if evaluated is not None:
            series, ligand_basenames = evaluated
            rows.append(series)
            ligand_basenames_by_target[pdb_id] = ligand_basenames

    if not rows:
        print("[ERR] Nothing evaluated.")
        raise SystemExit(3)

    df = pd.DataFrame(rows).sort_values("pdb_id")

    if args.target_name_from_pdb:
        names: Dict[str, str] = {}
        for pdb_id in df["pdb_id"].tolist():
            # Prefer COMPND.MOLECULE -> UniProt entry -> accession, honoring CLI preference.
            names[pdb_id] = derive_target_name(
                pdb_id,
                prefer=args.target_name_prefer,
                pdb_root_override=pdb_root_override,
                cfg=cfg,
            )
        df["target_name"] = df["pdb_id"].map(names).fillna("")

    if args.report_library:
        libs: Dict[str, str] = {}
        for pdb_id in df["pdb_id"].tolist():
            lig_basenames = ligand_basenames_by_target.get(pdb_id, set())
            # Match ligands against prepped_ligands/<target>/<library>/ (or manifest.json).
            libs[pdb_id] = infer_library_name(
                pdb_id,
                lig_basenames,
                prepped_root_override=prepped_root_override,
                cfg=cfg,
            )
        df["library_name"] = df["pdb_id"].map(libs).fillna("")

    analysis_root.mkdir(parents=True, exist_ok=True)

    summary_name = "summary.tsv"
    if args.run_id:
        run_suffix = str(args.run_id).replace(" ", "")
        summary_name = f"summary{run_suffix}.tsv"
    summary_path = analysis_root / summary_name

    variant_col_present = "variant" in df.columns
    variant_upper = df["variant"].astype(str).str.upper() if variant_col_present else None
    apo_mask = (variant_upper == "APO") if variant_col_present else None
    holo_mask = (variant_upper == "HOLO") if variant_col_present else None
    has_sections = bool(variant_col_present and ((apo_mask is not None and apo_mask.any()) or (holo_mask is not None and holo_mask.any())))

    with open(summary_path, "w", newline="") as fh:
        if has_sections:
            apo_rows = df.loc[apo_mask].sort_values("pdb_id") if apo_mask is not None else pd.DataFrame()
            holo_rows = df.loc[holo_mask].sort_values("pdb_id") if holo_mask is not None else pd.DataFrame()
            n_apo = len(apo_rows)
            n_holo = len(holo_rows)

            if n_apo:
                fh.write("Apo\n")
                apo_rows.drop(columns=["variant"], errors="ignore").to_csv(fh, sep="\t", index=False)
            if n_holo:
                if n_apo:
                    fh.write("\n")
                fh.write("Holo\n")
                holo_rows.drop(columns=["variant"], errors="ignore").to_csv(fh, sep="\t", index=False)

            leftover_mask = ~(apo_mask | holo_mask) if (apo_mask is not None and holo_mask is not None) else pd.Series(False, index=df.index)
            leftover_rows = df.loc[leftover_mask].sort_values("pdb_id") if not df.empty else pd.DataFrame()
            if not leftover_rows.empty:
                if n_apo or n_holo:
                    fh.write("\n")
                leftover_rows.drop(columns=["variant"], errors="ignore").to_csv(fh, sep="\t", index=False)

            print(f"[eval] sections: Apo={n_apo} Holo={n_holo}")
        else:
            df.drop(columns=["variant"], errors="ignore").to_csv(fh, sep="\t", index=False)

    excluded_cols = {"pdb_id", "N", "n_actives", "actives_fraction", "variant", "target_name", "library_name"}
    metric_cols = [c for c in df.columns if c not in excluded_cols]
    macro = df[metric_cols].mean(numeric_only=True).to_dict()
    pd.DataFrame([{"pdb_id":"macro_avg", **{k: macro[k] for k in metric_cols}}]) \
      .to_csv(analysis_root / "summary_macro.tsv", sep="\t", index=False)

    print(f"[OK] Wrote per-target results to: {analysis_root}")
    print(f"[eval] summary_out={summary_path}")

if __name__ == "__main__":
    main()
