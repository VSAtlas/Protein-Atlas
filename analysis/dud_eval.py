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
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc

LOG_LEVELS = {"DEBUG": 10, "INFO": 20, "WARN": 30, "ERROR": 40}
LEVEL_ABBREV = {"DEBUG": "DBG", "INFO": "INF", "WARN": "WRN", "ERROR": "ERR"}
_current_log_level = LOG_LEVELS["INFO"]
_BACKEND_LOGGED = False

FULL_RUN_MIN_LIGANDS = int(os.environ.get("FULL_RUN_MIN_LIGANDS", 10))


def set_log_level(level_name: Optional[str]) -> None:
    global _current_log_level
    key = (level_name or "INFO").strip().upper()
    if key not in LOG_LEVELS:
        key = "INFO"
    _current_log_level = LOG_LEVELS[key]


def dbg(level: str, tag: str, message: str) -> None:
    upper = (level or "INFO").strip().upper()
    if upper not in LOG_LEVELS:
        upper = "INFO"
    if LOG_LEVELS[upper] < _current_log_level:
        return
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    label = LEVEL_ABBREV.get(upper, upper[:3])
    print(f"[{stamp}][{label}][{tag}] {message}")


CONTROL_CENTERS_RE = re.compile(
    r"\[control-centers\]\s+"
    r"n=(\d+)\s+"
    r"max\?=([0-9]+(?:\.[0-9]+)?)\s*A?\b\s+"
    r"policy=(.+)$",
    re.IGNORECASE,
)

CONTROL_REDOCK_RE = re.compile(
    r"\[control-redock\]\s+"
    r"lig=([^\s]+)\s+"
    r"rmsd=([0-9]+(?:\.[0-9]+)?)\s*A?\b\s+"
    r"score=([0-9.\-]+)",
    re.IGNORECASE,
)

_CONTROL_PATTERNS_LOGGED = False

# >>> PATHS IMPORT START
from pathlib import Path
# --- ensure repo root is importable when running from analysis/ ---
import sys
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from path_router import make_paths, expand_variants



def _fmt_counts_and_round(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with big counts comma-formatted and floats rounded to 3 dp (for display only)."""
    out = df.copy()
    # Comma-format counts if present
    for c in ("N", "n_actives"):
        if c in out.columns:
            out[c] = out[c].map(lambda v: f"{int(v):,}" if pd.notna(v) else "")
    # Round floats to 3 decimals (display only)
    for c in out.select_dtypes(include=["float", "float64"]).columns:
        if c not in ("actives_fraction",):  # optional: keep full precision here if you prefer
            out[c] = out[c].map(lambda v: f"{v:.3f}" if pd.notna(v) else "")
    return out

def _df_to_pretty_text(df: pd.DataFrame, title: str | None = None) -> str:
    """
    Produce an aligned, readable table using pandas' to_string.
    Assumes df is already ordered and formatted for display.
    """
    # Choose widths implicitly; pandas to_string aligns numbers right by default
    body = df.to_string(index=False)
    return (title + "\n" + body) if title else body

def _write_pretty_summary(sections: list[tuple[str | None, pd.DataFrame]], out_path: Path) -> None:
    """
    sections: list of (title, df) pairs. title can be None for single-table case.
    Writes a single text file with optional section headers separated by blank lines.
    """
    parts = []
    for title, df in sections:
        if df is None or df.empty:
            continue
        df_disp = _fmt_counts_and_round(df)
        parts.append(_df_to_pretty_text(df_disp, title=title))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n\n".join(parts) + ("\n" if parts else ""))


def _write_pretty_table_noformat(df: pd.DataFrame, out_path: Path, title: str | None = None) -> None:
    """
    Write a readable, aligned text table without changing any values.
    This preserves rows, columns, order, and raw numeric formatting.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        if title:
            f.write(f"{title}\n")
        f.write(df.to_string(index=False))
        f.write("\n")





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
    """
    Infer library folder name when libraries live at:
        <prepped_root>/<library>/*files*
    No per-PDB subdirectory is expected.

    Preference:
      1) (currently no json but might make one in the future) manifest.json at library root with "library_name"/"name"/"library"/"label"
      2) filename-overlap between docked basenames and files in <library> (non-recursive; optional shallow)
    """

    # --- gather candidate roots from CLI override + config
    candidate_roots: List[Path] = []
    if prepped_root_override:
        candidate_roots.append(prepped_root_override)
    if cfg:
        try:
            paths = make_paths(cfg, base_id=target_id, pdb_file=f"{target_id}.pdb")
            candidate_roots.append(paths.prepped_root)
        except Exception:
            pass

    roots: List[Path] = []
    seen: Set[Path] = set()
    for r in candidate_roots:
        if r and r not in seen:
            seen.add(r)
            roots.append(r)

    # DEBUG: show roots + docked basename sample
    if roots:
        dbg("DEBUG", "library", f"pdb={target_id} search_roots={';'.join(str(r) for r in roots)}")
    else:
        dbg("WARN", "library", f"pdb={target_id} no candidate roots for library inference")

    if ligand_basenames:
        dbg("DEBUG", "library",
            f"pdb={target_id} docked_basenames_n={len(ligand_basenames)} "
            f"sample={list(sorted(ligand_basenames))[:3]}")
    else:
        dbg("DEBUG", "library", f"pdb={target_id} docked_basenames_n=0")

    def _manifest_label(manifest_path: Path) -> Optional[str]:
        try:
            with manifest_path.open("r", encoding="utf-8") as fh:
                data = json.load(fh)
        except Exception:
            return None
        for key in ("library_name", "name", "library", "label"):
            val = data.get(key)
            if isinstance(val, str) and val.strip():
                return val.strip()
        return None

    def _collect_basenames_top(folder: Path) -> Set[str]:
        """Non-recursive: only files directly under <library>."""
        allowed = {".pdbqt", ".sdf", ".mol2"}
        out: Set[str] = set()
        try:
            for e in folder.iterdir():
                if e.is_file() and e.suffix.lower() in allowed:
                    out.add(e.name)
        except FileNotFoundError:
            pass
        return out

    # (Optional) shallow peek—uncomment if you want actives/decoys support without full recursion
    def _collect_basenames_shallow(folder: Path) -> Set[str]:
        allowed = {".pdbqt", ".sdf", ".mol2"}
        out: Set[str] = set()
        # top level
        out |= _collect_basenames_top(folder)
        # shallow subdirs commonly used
        for sub in ("actives", "decoys"):
            subdir = folder / sub
            if subdir.is_dir():
                try:
                    for e in subdir.iterdir():
                        if e.is_file() and e.suffix.lower() in allowed:
                            out.add(e.name)
                except FileNotFoundError:
                    pass
        return out

    best_name: Optional[str] = None
    best_score = -1
    have_tie = False

    for root in roots:
        dbg("DEBUG", "library", f"pdb={target_id} scanning_root={root}")
        if not root.exists():
            dbg("WARN", "library", f"pdb={target_id} root_missing={root}")
            continue

        # iterate each <library> folder at the root
        for library_dir in sorted(root.iterdir()):
            if not library_dir.is_dir():
                continue

            dbg("DEBUG", "library", f"pdb={target_id} candidate_library={library_dir.name}")

            # 1) manifest preference
            manifest = library_dir / "manifest.json"
            if manifest.exists():
                label = _manifest_label(manifest)
                dbg("DEBUG", "library", f"pdb={target_id} manifest_found={manifest} label={label or '<none>'}")
                if label:
                    dbg("INFO", "library",
                        f"pdb={target_id} picked='{label}' method='manifest' root={library_dir}")
                    return label

            # 2) filename-overlap (non-recursive; flip to _collect_basenames_shallow if needed)
            basenames = _collect_basenames_top(library_dir)
            dbg("DEBUG", "library",
                f"pdb={target_id} library={library_dir.name} files_seen={len(basenames)} "
                f"sample={list(sorted(basenames))[:3] if basenames else []}")

            if not basenames:
                # Uncomment next two lines to consider shallow subdirs if top-level empty
                # basenames = _collect_basenames_shallow(library_dir)
                # dbg("DEBUG", "library", f"pdb={target_id} library={library_dir.name} shallow_files_seen={len(basenames)}")
                if not basenames:
                    continue

            score = len(ligand_basenames & basenames) if ligand_basenames else 0
            dbg("DEBUG", "library", f"pdb={target_id} library={library_dir.name} overlap_score={score}")
            if score > best_score:
                best_score = score
                best_name = library_dir.name
                have_tie = False
            elif score == best_score and score >= 0:
                have_tie = True

    if best_name and best_score > 0 and not have_tie:
        dbg("INFO", "library",
            f"pdb={target_id} picked='{best_name}' method='filename_overlap' score={best_score}")
        return best_name
    if best_name and have_tie:
        dbg("WARN", "library", f"pdb={target_id} filename_overlap ties score={best_score}")
    dbg("WARN", "library", f"pdb={target_id} no library match found in roots")
    return ""


@dataclass
class TargetEvaluation:
    metrics: Optional[pd.Series]
    ligand_basenames: Set[str]
    has_run_id_column: bool
    run_ids: Set[str]


def derive_target_name(target_id: str,
                       *,
                       prefer: str,
                       pdb_root_override: Optional[Path],
                       cfg: Dict) -> str:
    pdb_path = Path("input_pdbs") / f"{target_id}.pdb"
    exists = pdb_path.exists()
    print(f"[dbg.target.source] pdb={target_id} path=input_pdbs/{target_id}.pdb exists={str(exists)}")
    if not exists:
        return ""

    header_lines = _read_pdb_header_lines(pdb_path)
    compnd = extract_compnd_molecules(header_lines)
    entries, accessions = extract_uniprot_from_dbref(header_lines)
    compnd_label = compnd[0] if compnd else ""
    uniprot_label = accessions[0] if accessions else (entries[0] if entries else "")
    print(f"[dbg.target.extract] pdb={target_id} compnd='{compnd_label}' uniprot='{uniprot_label}' prefer={prefer}")
    choice = choose_target_name(compnd, entries, accessions, prefer=prefer)
    source = ""
    compnd_choice = ", ".join(compnd) if compnd else ""
    if choice:
        if compnd_choice and choice == compnd_choice:
            source = "PDB:COMPND"
        elif entries and choice == entries[0]:
            source = "PDB:DBREF_ENTRY"
        elif accessions and choice == accessions[0]:
            source = "PDB:DBREF_ACCESSION"
    source_label = source if source else "none"
    dbg("DEBUG", "target", f"pdb={target_id} name='{choice}' source='{source_label}'")
    return choice
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

# >>> RUN-SELECTION START
def select_default_run_id(targets: List[Tuple[str, Path]],
                          csv_paths_by_target: Dict[str, Path],
                          lig_col_cli: Optional[str],
                          score_col_cli: Optional[str]) -> Optional[str]:
    run_full_map: Dict[str, Set[str]] = {}
    run_mtimes: Dict[str, float] = {}

    digit_re = re.compile(r"^\d{8,}$")
    iso_re = re.compile(r"^\d{4}-\d{2}-\d{2}([Tt _].*)?$")

    def _parse_isoish(token: str) -> Optional[float]:
        text = token.strip()
        if not text:
            return None
        cleaned = text.rstrip("Z").rstrip("z")
        candidates = [cleaned]
        if "T" in cleaned:
            candidates.append(cleaned.replace("T", " "))
        if "_" in cleaned:
            candidates.append(cleaned.replace("_", " "))
        if cleaned.endswith("T"):
            candidates.append(cleaned.rstrip("T"))
        for cand in candidates:
            cand = cand.strip()
            if not cand:
                continue
            try:
                dt = datetime.fromisoformat(cand)
            except Exception:
                continue
            try:
                return dt.replace(tzinfo=dt.tzinfo or timezone.utc).timestamp()
            except Exception:
                try:
                    return dt.timestamp()
                except Exception:
                    continue
        return None

    for pdb_id, _ in targets:
        csv_path = csv_paths_by_target.get(pdb_id)
        if not csv_path or not csv_path.exists():
            continue
        try:
            df = pd.read_csv(csv_path)
        except Exception as exc:
            dbg("WARN", "run", f"pdb={pdb_id} auto_select_read_err={exc}")
            continue

        if "run_id" not in df.columns:
            continue

        try:
            lig_col = guess_ligfile_col(df, lig_col_cli)
        except Exception as exc:
            dbg("WARN", "run", f"pdb={pdb_id} auto_select_ligcol_err={exc}")
            continue

        lig_series = df[lig_col].astype(str)
        parsed = lig_series.apply(parse_name_and_label)
        lig_ids = parsed.apply(lambda t: t[0])

        run_series = df["run_id"].dropna()
        if run_series.empty:
            continue
        run_ids = run_series.index
        df_subset = df.loc[run_ids].copy()
        df_subset["__run_id"] = run_series.astype(str).str.strip()
        df_subset["__lig_id"] = lig_ids.loc[run_ids]
        df_subset = df_subset[df_subset["__run_id"] != ""]
        if df_subset.empty:
            continue

        lig_counts = df_subset.groupby("__run_id")["__lig_id"].nunique()
        for run_label, count in lig_counts.items():
            if count >= FULL_RUN_MIN_LIGANDS:
                run_full_map.setdefault(run_label, set()).add(pdb_id)
            try:
                stat = csv_path.stat()
            except Exception:
                continue
            run_mtimes[run_label] = max(run_mtimes.get(run_label, 0.0), getattr(stat, "st_mtime", 0.0))

    if not run_full_map:
        return None

    scored = {run_id: len(pdbs) for run_id, pdbs in run_full_map.items()}
    best_score = max(scored.values())
    candidates = [run_id for run_id, score in scored.items() if score == best_score]
    if not candidates:
        return None

    parsed_values: Dict[str, float] = {}
    all_parseable = True
    for run_id in candidates:
        token = run_id.strip()
        parsed_val: Optional[float] = None
        if digit_re.match(token):
            try:
                parsed_val = float(int(token))
            except Exception:
                parsed_val = None
        elif iso_re.match(token):
            parsed_val = _parse_isoish(token)
        if parsed_val is None:
            all_parseable = False
            break
        parsed_values[run_id] = parsed_val

    if all_parseable and parsed_values:
        return max(parsed_values, key=lambda k: (parsed_values[k], k))

    return max(candidates, key=lambda k: (run_mtimes.get(k, 0.0), k))
# >>> RUN-SELECTION END

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
                    run_id: Optional[str]) -> Optional[TargetEvaluation]:
    try:
        df = pd.read_csv(csv_path)
        dbg("INFO", "csv", f"pdb={pdb_id} path={csv_path} rows={len(df)}")
    except Exception as e:
        dbg("ERROR", "csv", f"pdb={pdb_id} path={csv_path} err={e}")
        return None

    has_run_id_col = "run_id" in df.columns
    if run_id:
        total_rows = len(df)
        if has_run_id_col:
            mask = df["run_id"].astype(str) == str(run_id)
            df = df.loc[mask].copy()
            dbg("INFO", "filter", f"pdb={pdb_id} run_id={run_id} kept={len(df)}/{total_rows}")
        else:
            dbg("WARN", "filter", f"pdb={pdb_id} run_id={run_id} column_missing proceeding_unfiltered")

    dbg("DEBUG", "schema", f"pdb={pdb_id} cols={list(df.columns)} has_variant={'variant' in df.columns} has_run_id={'run_id' in df.columns}")

    lig_col = guess_ligfile_col(df, lig_col_cli)
    score_col = guess_score_col(df, score_col_cli)
    dbg("DEBUG", "schema", f"pdb={pdb_id} ligand_col={lig_col} source={'CLI' if lig_col_cli else 'auto'} score_col={score_col} source={'CLI' if score_col_cli else 'auto'}")
    # Coerce scores to numeric; drop NaN/±inf early to avoid NaNs in metrics
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    before_nf = len(df)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna(subset=[score_col])
    drop_nonfinite = before_nf - len(df)

    raw_lig = df[lig_col]
    missing_names = int(raw_lig.isna().sum())
    empty_names = int(raw_lig.astype(str).str.strip().eq("").sum())
    # derive ligand_id + label from filename
    lig_paths = raw_lig.astype(str)
    parsed = lig_paths.apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])
    used_ligand_basenames: Set[str] = {os.path.basename(p) for p in lig_paths.tolist() if p}

    run_ids_present: Set[str] = set()
    if has_run_id_col and "run_id" in df.columns:
        try:
            run_ids_present = set(df["run_id"].dropna().astype(str).unique().tolist())
        except Exception:
            run_ids_present = set()

    # drop rows where label couldn't be inferred
    before = len(df)
    df = df.dropna(subset=["is_active"])
    dropped_token = before - len(df)

    # collapse to best score per ligand (min score)
    best = df.groupby("lig_id", as_index=False).agg(
        best_score=(score_col, "min"),
        is_active=("is_active", "max"),   # any active -> active
    )
    # Some ligands may still have all-NaN scores → best_score=NaN; drop them
    before_best = len(best)
    best = best.dropna(subset=["best_score"])
    drop_nan_best = before_best - len(best)

    y_true = best["is_active"].astype(int).to_numpy()
    y_low = best["best_score"].to_numpy()
    y_high = -y_low

    N = len(best); n_act = int(y_true.sum())
    missing_name_detected = missing_names + empty_names
    dbg("DEBUG", "screen", f"pdb={pdb_id} missing_name_detected={missing_name_detected} kept_rows={len(df)} ligands={N}")
    if N == 0 or n_act == 0 or n_act == N:
        dbg("WARN", "metrics", f"pdb={pdb_id} degenerate_set N={N} actives={n_act}")
        return TargetEvaluation(metrics=None,
                                ligand_basenames=used_ligand_basenames,
                                has_run_id_column=has_run_id_col,
                                run_ids=run_ids_present)

    # metrics
    ef = ef_at_fractions(y_true, y_low, fractions=(0.01,0.02,0.05,0.10))
    dbg("DEBUG", "metrics", f"pdb={pdb_id} start N={N} n_actives={n_act}")
    try:
        rocAUC = float(roc_auc_score(y_true, y_high))
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} metric=ROC_AUC err={exc}")
        rocAUC = float("nan")
    try:
        fpr, tpr, _ = roc_curve(y_true, y_high)
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} metric=ROC_curve err={exc}")
        fpr = np.array([0.0, 1.0])
        tpr = np.array([0.0, 1.0])
    try:
        prAUC, precision, recall = pr_auc(y_true, y_high)
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} metric=PR_AUC err={exc}")
        prAUC = float("nan")
        precision = np.array([0.0, 1.0])
        recall = np.array([0.0, 1.0])
    try:
        lAUC = log_auc_from_roc(fpr, tpr, lam=logauc_lambda)
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} metric=logAUC err={exc}")
        lAUC = float("nan")
    lAUC_adj = lAUC - 0.14462
    try:
        bed = bedroc(y_true, y_high, alpha=bedroc_alpha)
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} metric=BEDROC err={exc}")
        bed = float("nan")

    # plots
    out_dir.mkdir(parents=True, exist_ok=True)
    global _BACKEND_LOGGED
    if not _BACKEND_LOGGED:
        try:
            backend = plt.get_backend()
        except Exception as exc:
            dbg("WARN", "mpl", f"pdb={pdb_id} backend_detect_failed err={exc}")
        else:
            dbg("DEBUG", "mpl", f"backend={backend}")
        _BACKEND_LOGGED = True

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
    dbg("INFO", "screen", f"pdb={pdb_id} drop_nonfinite={drop_nonfinite} missing_name_detected={missing_name_detected} drop_missing_token={dropped_token} drop_nan_best_score={drop_nan_best} kept_ligands={N} actives={n_act} active_fraction={base:.3f} token_regex='active|decoy'")
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
                unique_vals = sorted(set(upper.tolist()))
                counts = {val: int((upper == val).sum()) for val in unique_vals}
                dbg("DEBUG", "variant", f"pdb={pdb_id} unique={unique_vals} counts={counts}")

    row = {
        "run_id": str(run_id) if run_id else "(none)",
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
    dbg("DEBUG", "metrics", f"pdb={pdb_id} ROC_AUC={rocAUC:.3f} PR_AUC={prAUC:.3f} BEDROC={bed:.3f}")
    return TargetEvaluation(metrics=pd.Series(row),
                            ligand_basenames=used_ligand_basenames,
                            has_run_id_column=has_run_id_col,
                            run_ids=run_ids_present)

def _load_default_cfg() -> Dict:
    cfg_path = Path("config.txt")
    try:
        from input_and_export_functions import load_config, validate_config
    except Exception:
        dbg("WARN", "config", "input_and_export_functions unavailable; skipping config.txt")
        return {}

    try:
        if not cfg_path.exists():
            dbg("DEBUG", "config", "config.txt not found; using CLI defaults")
            return {}
        cfg = load_config("config.txt") or {}
        if cfg:
            validate_config(cfg)
            dbg("DEBUG", "config", f"config.txt loaded keys={sorted(cfg.keys())}")
        else:
            dbg("DEBUG", "config", "config.txt empty; using CLI defaults")
        return cfg
    except Exception as exc:
        dbg("WARN", "config", f"failed to load config.txt err={exc}")
        return {}


def _resolve_run_label(pdb_id: str,
                       meta: Optional[TargetEvaluation],
                       cli_run_id: Optional[str],
                       active_run_id: Optional[str]) -> str:
    if cli_run_id:
        if meta and meta.has_run_id_column:
            return str(cli_run_id)
        return "(none)"
    if meta and meta.has_run_id_column:
        if meta.run_ids:
            ordered = sorted(str(x) for x in meta.run_ids if x)
            if not ordered:
                return ""
            return ordered[0] if len(ordered) == 1 else ",".join(ordered)
        return ""
    if active_run_id:
        return str(active_run_id)
    return "(none)"


def _candidate_protein_logs(pdb_id: str,
                            csv_path: Path,
                            docked_root: Path,
                            log_root_override: Optional[Path],
                            cfg: Dict) -> Tuple[Optional[Path], List[Path]]:
    candidates: List[Path] = []
    seen: Set[str] = set()

    def _add(path: Path) -> None:
        try:
            key = str(path.resolve())
        except Exception:
            key = str(path)
        if key not in seen:
            seen.add(key)
            candidates.append(path)

    if log_root_override:
        _add(log_root_override / pdb_id / "protein.log")
    else:
        _add(Path("docked") / pdb_id / "protein.log")

    if docked_root.is_dir():
        if docked_root.name.upper() == pdb_id.upper():
            _add(docked_root / "protein.log")
        _add(docked_root / pdb_id / "protein.log")
    else:
        _add(docked_root / pdb_id / "protein.log")

    csv_parent = csv_path.parent
    _add(csv_parent / "protein.log")
    if csv_parent.name.upper() != pdb_id.upper():
        _add(csv_parent.parent / "protein.log")

    if cfg:
        try:
            paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        except Exception as exc:
            dbg("WARN", "control", f"pdb={pdb_id} router_err={exc}")
        else:
            _add(paths.docked_root / paths.pdb_id / "protein.log")

    selected: Optional[Path] = None
    for cand in candidates:
        if cand.exists():
            selected = cand
            break
    return selected, candidates


def _build_control_records(pdb_id: str,
                           run_id_label: str,
                           target_name: str,
                           library_name: str,
                           log_path: Optional[Path]) -> Tuple[List[Dict], Dict, Dict[str, int]]:
    long_rows: List[Dict] = []
    controls_found = 0
    max_spread: Optional[float] = None
    policy = ""
    centers_detected = False

    if not log_path or not log_path.exists():
        summary_row = {
            "pdb_id": pdb_id,
            "run_id": run_id_label,
            "target_name": target_name,
            "library_name": library_name,
            "controls_found": 0,
            "control_id": "",
            "rmsd_to_crystal_A": None,
            "redock_best_energy_kcal_mol": None,
            "chosen": 0,
            "max_spread_A": None,
            "policy": "",
        }
        meta = {"centers_found": 0, "redock_lines": 0, "report_rows": 0}
        return long_rows, summary_row, meta

    try:
        with log_path.open("r", encoding="utf-8", errors="ignore") as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line:
                    continue
                m_centers = CONTROL_CENTERS_RE.search(line)
                if m_centers:
                    centers_detected = True
                    try:
                        controls_found = int(m_centers.group(1))
                    except Exception:
                        controls_found = 0
                    try:
                        max_spread = float(m_centers.group(2))
                    except Exception:
                        max_spread = None
                    policy = m_centers.group(3)
                    continue
                m_control = CONTROL_REDOCK_RE.search(line)
                if not m_control:
                    continue
                control_id = m_control.group(1)
                try:
                    rmsd = float(m_control.group(2))
                except Exception:
                    rmsd = float("nan")
                try:
                    score = float(m_control.group(3))
                except Exception:
                    score = float("nan")
                long_rows.append({
                    "pdb_id": pdb_id,
                    "run_id": run_id_label,
                    "target_name": target_name,
                    "library_name": library_name,
                    "controls_found": 0,  # fill later
                    "control_id": control_id,
                    "rmsd_to_crystal_A": rmsd,
                    "redock_best_energy_kcal_mol": score,
                    "chosen": 0,
                    "max_spread_A": None,
                    "policy": "",
                })
    except Exception as exc:
        dbg("WARN", "control", f"pdb={pdb_id} log_read_err={exc}")
        long_rows.clear()
        controls_found = 0
        max_spread = None
        policy = ""

    if controls_found == 0 and long_rows:
        controls_found = len(long_rows)

    chosen_index: Optional[int] = None
    if long_rows:
        # Fill common fields first
        for row in long_rows:
            row["controls_found"] = controls_found
            row["max_spread_A"] = max_spread
            row["policy"] = policy

        def _rmsd_key(val: Optional[float]) -> float:
            if val is None:
                return float("inf")
            try:
                return float(val) if not math.isnan(float(val)) else float("inf")
            except Exception:
                return float("inf")

        ordered = sorted(
            enumerate(long_rows),
            key=lambda pair: (
                _rmsd_key(pair[1]["rmsd_to_crystal_A"]),
                str(pair[1]["control_id"]),
            ),
        )
        if ordered:
            chosen_index = ordered[0][0]
            long_rows[chosen_index]["chosen"] = 1
            chosen_row = long_rows[chosen_index]
            dbg(
                "DEBUG",
                "control",
                (
                    f"pdb={pdb_id} chosen_control={chosen_row['control_id']} "
                    f"best_rmsd_A={chosen_row['rmsd_to_crystal_A']} "
                    f"best_score_kcal={chosen_row['redock_best_energy_kcal_mol']}"
                ),
            )

    summary_row = {
        "pdb_id": pdb_id,
        "run_id": run_id_label,
        "target_name": target_name,
        "library_name": library_name,
        "controls_found": controls_found,
        "control_id": "",
        "rmsd_to_crystal_A": None,
        "redock_best_energy_kcal_mol": None,
        "chosen": 0,
        "max_spread_A": max_spread,
        "policy": policy,
    }

    if chosen_index is not None:
        chosen = long_rows[chosen_index]
        summary_row.update({
            "control_id": chosen["control_id"],
            "rmsd_to_crystal_A": chosen["rmsd_to_crystal_A"],
            "redock_best_energy_kcal_mol": chosen["redock_best_energy_kcal_mol"],
            "chosen": 1,
        })

    meta = {
        "centers_found": 1 if centers_detected else 0,
        "redock_lines": len(long_rows),
        "report_rows": len(long_rows),
    }

    return long_rows, summary_row, meta

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
    ap.add_argument("--log-level", type=str, default="INFO",
                    choices=("DEBUG", "INFO", "WARN", "ERROR"),
                    help="Logging verbosity (default: INFO).")
    ap.add_argument("--target-name-from-pdb", action="store_true", default=True,
                    help="If set, add a target_name column derived from PDB headers.")
    ap.add_argument("--target-name-prefer", type=str, default="auto", 
                    choices=("auto", "compnd", "uniprot"),
                    help="Preference order when selecting target_name (default: auto).")
    ap.add_argument("--pdb-root", type=str, default=None,
                    help="Optional override root for processed PDB folders (processed_pdbs/<target>/).")
    ap.add_argument("--report-library", action="store_true", default=True,
                    help="If set, attempt to infer library_name from prepped_ligands.")
    ap.add_argument("--prepped-root", type=str, default=None,
                    help="Optional override root for prepped ligands (prepped_ligands/<target>/).")
    ap.add_argument("--emit-control-report", action="store_true", default=True,
                    help="If set, parse protein.log control entries and emit control_redock TSVs.")
    ap.add_argument("--log-root", type=str, default=None,
                    help="Optional override root containing <PDB>/protein.log (default: docked/).")
    try:
        ap.add_argument("--pretty-summary", action=argparse.BooleanOptionalAction, default=True,
                        help="Also write a human-readable aligned text summary (default: on).")
    except Exception:
        ap.add_argument("--pretty-summary", dest="pretty_summary", action="store_true",
                        help="Write a human-readable aligned text summary.")
        ap.add_argument("--no-pretty-summary", dest="pretty_summary", action="store_false")
        ap.set_defaults(pretty_summary=True)

    args = ap.parse_args()

    set_log_level(args.log_level)
    dbg("DEBUG", "args", f"log_level={args.log_level} target_name_from_pdb={'ON' if args.target_name_from_pdb else 'OFF'} report_library={'ON' if args.report_library else 'OFF'} control_report={'ON' if args.emit_control_report else 'OFF'} run_id={args.run_id or 'none'}")

    cfg = _load_default_cfg()
    docked_root = Path(args.docked_root)
    out_root = Path(args.out_dir)
    pdb_root_override = Path(args.pdb_root) if args.pdb_root else None
    prepped_root_override = Path(args.prepped_root) if args.prepped_root else None
    log_root_override = Path(args.log_root) if getattr(args, "log_root", None) else None

    analysis_root = out_root
    if cfg and "OVERALL_DIR" in cfg:
        # >>> ANALYSIS OUT ROOT PATCH START
        analysis_root = Path(cfg["OVERALL_DIR"]) / "analysis" / "out"
        # >>> ANALYSIS OUT ROOT PATCH END
        dbg("DEBUG", "config", f"OVERALL_DIR override applied analysis_root={analysis_root}")
    analysis_root.mkdir(parents=True, exist_ok=True)
    dbg("DEBUG", "paths", f"docked_root={docked_root} out_root={out_root} analysis_root={analysis_root} pdb_root={pdb_root_override or 'none'} prepped_root={prepped_root_override or 'none'}")

    def _resolve_docking_csv(pdb_id: str, fallback_dir: Path) -> Optional[Path]:
        candidates_tried: List[str] = []
        picked: Optional[Path] = None

        def _record(path: Path) -> bool:
            nonlocal picked
            exists = path.exists()
            state = "hit" if exists else "miss"
            candidates_tried.append(f"{path} ({state})")
            if exists and picked is None:
                picked = path
                return True
            return False

        if not cfg:
            candidate = fallback_dir / "docking_score_long.csv"
            _record(candidate)
        else:
            try:
                # >>> PATHS INIT START
                mode    = str(cfg.get("APO_HOLO_MODE", "")).strip()
                variants = expand_variants(mode)  # returns [None] | ["APO","HOLO"]
                paths   = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
                ph_token = cfg.get("PH_TOKEN") or None
                # >>> PATHS INIT END
            except Exception as exc:
                dbg("WARN", "resolve", f"pdb={pdb_id} router_err={exc}")
                candidate = fallback_dir / "docking_score_long.csv"
                _record(candidate)
            else:
                # >>> DOCKED INPUT PATHS PATCH START
                docked_root_cfg = paths.docked_pdb_root()
                summary_csv = paths.docking_score_summary_csv()
                long_csv    = paths.docking_score_long_csv()
                # >>> DOCKED INPUT PATHS PATCH END

                if docked_root_cfg.exists():
                    if _record(long_csv):
                        pass
                    else:
                        for variant in variants:
                            variant_root = paths.docked_variant_root(variant)
                            candidate = variant_root / "docking_score_long.csv"
                            if _record(candidate):
                                break
                    if picked is None:
                        candidate = fallback_dir / "docking_score_long.csv"
                        _record(candidate)
                    # Touch summary path to exercise router (no fallback to summary file for eval)
                    summary_csv.exists()
                else:
                    dbg("WARN", "resolve", f"pdb={pdb_id} docked_root_missing root={docked_root_cfg}")
                    candidate = fallback_dir / "docking_score_long.csv"
                    _record(candidate)

        if picked is not None:
            dbg("DEBUG", "resolve", f"pdb={pdb_id} tried={len(candidates_tried)} candidates={candidates_tried}")
            dbg("INFO", "resolve", f"pdb={pdb_id} picked={picked}")
            return picked

        dbg("WARN", "resolve", f"pdb={pdb_id} no_csv_found tried={candidates_tried or ['<none>']} search_root={fallback_dir}")
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
        dbg("ERROR", "discover", f"no docking_score_long.csv under {docked_root}")
        raise SystemExit(2)

    dbg("INFO", "discover", f"targets={len(targets)} root={docked_root}")

    rows: List[pd.Series] = []
    ligand_basenames_by_target: Dict[str, Set[str]] = {}
    target_eval_results: Dict[str, TargetEvaluation] = {}
    csv_paths_by_target: Dict[str, Path] = {}
    for pdb_id, csvp in targets:
        csv_paths_by_target[pdb_id] = csvp

    # >>> ACTIVE-RUN PICK START
    active_run_id = args.run_id
    if not active_run_id:
        try:
            active_run_id = select_default_run_id(targets, csv_paths_by_target, args.lig_col, args.score_col)
        except Exception as _exc:
            active_run_id = None
            dbg("WARN", "run", f"auto_select_failed err={_exc}")
    dbg("INFO", "run", f"active={active_run_id or '(none)'} source={'CLI' if args.run_id else 'auto'}")
    # >>> ACTIVE-RUN PICK END

    for pdb_id, csvp in targets:
        evaluated = evaluate_target(
            pdb_id=pdb_id,
            csv_path=csvp,
            out_dir=analysis_root / pdb_id,
            lig_col_cli=args.lig_col,
            score_col_cli=args.score_col,
            bedroc_alpha=args.bedroc_alpha,
            logauc_lambda=args.logauc_lambda,
            run_id=active_run_id,
        )
        if evaluated is not None:
            target_eval_results[pdb_id] = evaluated
            ligand_basenames_by_target[pdb_id] = evaluated.ligand_basenames
            if evaluated.metrics is not None:
                rows.append(evaluated.metrics)

    if not rows:
        dbg("ERROR", "metrics", "no targets produced evaluable rows")
        raise SystemExit(3)

    df_all = pd.DataFrame(rows).sort_values("pdb_id")

    control_targets = [pdb for pdb, _ in targets if pdb in target_eval_results]

    names: Dict[str, str] = {p: "" for p in control_targets}
    if args.target_name_from_pdb:
        for pdb_id in control_targets:
            names[pdb_id] = derive_target_name(
                pdb_id,
                prefer=args.target_name_prefer,
                pdb_root_override=pdb_root_override,
                cfg=cfg,
            )
        if not df_all.empty:
            df_all["target_name"] = df_all["pdb_id"].map(names).fillna("")
    else:
        dbg("DEBUG", "target", "target_name_from_pdb=OFF")

    libs: Dict[str, str] = {p: "" for p in control_targets}
    if args.report_library:
        for pdb_id in control_targets:
            lig_basenames = ligand_basenames_by_target.get(pdb_id, set())
            libs[pdb_id] = infer_library_name(
                pdb_id,
                lig_basenames,
                prepped_root_override=prepped_root_override,
                cfg=cfg,
            )
        if not df_all.empty:
            df_all["library_name"] = df_all["pdb_id"].map(libs).fillna("")
    else:
        dbg("DEBUG", "library", "report_library=OFF")

    df = df_all.copy()
    if not df.empty and "target_name" not in df.columns:
        df["target_name"] = df["pdb_id"].map(names).fillna("")
    if not df.empty and "library_name" not in df.columns:
        df["library_name"] = df["pdb_id"].map(libs).fillna("")

    analysis_root.mkdir(parents=True, exist_ok=True)

    # --- Exclude rows where the inferred library is FDA or the PDB-specific library ---
    # Normalize library names and PDB IDs for robust comparison
    lib_norm = df["library_name"].astype(str).str.strip()
    pdb_norm = df["pdb_id"].astype(str).str.strip()
    # FDA library (accept both "fda_library" and "fda")
    mask_fda = lib_norm.str.lower().isin({"fda_library", "fda"})
    # Per-PDB library: library name equals the PDB code (case-insensitive)
    mask_pdb = lib_norm.str.upper() == pdb_norm.str.upper()
    exclude_mask = mask_fda | mask_pdb
    # Emit a visible exclusions report with reasons
    if exclude_mask.any():
        excluded = df.loc[exclude_mask, ["pdb_id", "library_name"]].copy()
        # Reason per row
        excluded["reason"] = np.where(
            mask_fda.loc[exclude_mask],
            "library name = fda library",
            "library name = pdb",
        )
        excl_path = analysis_root / "excluded.tsv"
        excluded.to_csv(excl_path, sep="\t", index=False)
        dbg("INFO", "exclude", f"excluded={len(excluded)} out={excl_path}")
        # Also print each line for quick visibility
        for _, r in excluded.iterrows():
            print(f"[exclude] pdb={r['pdb_id']} library={r['library_name']} reason={r['reason']}")
    else:
        dbg("DEBUG", "exclude", "excluded=0")

    # Keep only the non-excluded rows for downstream metrics/summary
    df = df.loc[~exclude_mask].copy()

    control_long_records: List[Dict] = []
    control_summary_records: List[Dict] = []
    scanned_ok: List[str] = []
    missing_logs: List[str] = []
    parsed_counts: Dict[str, Dict[str, int]] = {}
    if getattr(args, "emit_control_report", False):
        global _CONTROL_PATTERNS_LOGGED
        if not _CONTROL_PATTERNS_LOGGED:
            dbg(
                "DEBUG",
                "control",
                f"searching_pattern.control_centers=\"{CONTROL_CENTERS_RE.pattern}\"",
            )
            dbg(
                "DEBUG",
                "control",
                f"searching_pattern.control_redock=\"{CONTROL_REDOCK_RE.pattern}\"",
            )
            _CONTROL_PATTERNS_LOGGED = True
        for pdb_id in control_targets:
            eval_meta = target_eval_results.get(pdb_id)
            run_label = _resolve_run_label(pdb_id, eval_meta, args.run_id, active_run_id)
            target_label = names.get(pdb_id, "")
            library_label = libs.get(pdb_id, "")
            csv_path = csv_paths_by_target.get(pdb_id)
            if csv_path is None:
                continue
            parsed_counts[pdb_id] = {"centers": 0, "redock_lines": 0, "report_rows": 0}
            selected_log, candidates = _candidate_protein_logs(
                pdb_id,
                csv_path,
                docked_root,
                log_root_override,
                cfg,
            )
            dbg("DEBUG", "control", f"pdb={pdb_id} log_candidates={[str(c) for c in candidates]}")
            if selected_log:
                dbg("INFO", "control", f"pdb={pdb_id} protein_log={selected_log}")
                scanned_ok.append(pdb_id)
            else:
                missing_logs.append(pdb_id)
                missing_hint = str(candidates[0]) if candidates else str(Path("docked") / pdb_id / "protein.log")
                dbg(
                    "WARN",
                    "control",
                    f"pdb={pdb_id} protein_log_missing={missing_hint} action=skip_control_parse",
                )

            long_rows, summary_row, meta = _build_control_records(
                pdb_id,
                run_label,
                target_label,
                library_label,
                selected_log,
            )
            if pdb_id not in parsed_counts:
                parsed_counts[pdb_id] = {"centers": 0, "redock_lines": 0, "report_rows": 0}
            parsed_counts[pdb_id]["centers"] = meta.get("centers_found", 0)
            parsed_counts[pdb_id]["redock_lines"] = meta.get("redock_lines", 0)
            parsed_counts[pdb_id]["report_rows"] = meta.get("report_rows", 0)
            dbg(
                "INFO",
                "control",
                (
                    f"pdb={pdb_id} centers_found={parsed_counts[pdb_id]['centers']} "
                    f"redock_lines={parsed_counts[pdb_id]['redock_lines']} "
                    f"report_rows={parsed_counts[pdb_id]['report_rows']}"
                ),
            )
            if long_rows:
                control_long_records.extend(long_rows)
            control_summary_records.append(summary_row)
        for pdb_id in control_targets:
            if pdb_id not in parsed_counts:
                parsed_counts[pdb_id] = {"centers": 0, "redock_lines": 0, "report_rows": 0}
                if pdb_id not in missing_logs and pdb_id not in scanned_ok:
                    missing_logs.append(pdb_id)

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

    preferred_cols = ("run_id", "target_name", "library_name", "pdb_id")
    excluded_cols = {"run_id", "target_name", "library_name", "pdb_id", "N", "n_actives", "actives_fraction", "variant"}

    with open(summary_path, "w", newline="") as fh:
        if has_sections:
            _df = df.drop(columns=["variant"], errors="ignore")
            cols = [c for c in preferred_cols if c in _df.columns] \
                   + [c for c in _df.columns if c not in preferred_cols]
            apo_rows = df.loc[apo_mask].sort_values("pdb_id") if apo_mask is not None else pd.DataFrame()
            holo_rows = df.loc[holo_mask].sort_values("pdb_id") if holo_mask is not None else pd.DataFrame()
            n_apo = len(apo_rows)
            n_holo = len(holo_rows)

            sections = []
            leftover_rows = pd.DataFrame()
            if n_apo:
                apo_out = apo_rows.drop(columns=["variant"], errors="ignore")
                a_cols = [c for c in preferred_cols if c in apo_out.columns] \
                         + [c for c in apo_out.columns if c not in preferred_cols]
                fh.write("Apo\n");
                apo_out[a_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("Apo", apo_out[a_cols]))

            if n_holo:
                if n_apo: fh.write("\n")
                holo_out = holo_rows.drop(columns=["variant"], errors="ignore")
                h_cols = [c for c in preferred_cols if c in holo_out.columns] \
                         + [c for c in holo_out.columns if c not in preferred_cols]
                fh.write("Holo\n");
                holo_out[h_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("Holo", holo_out[h_cols]))

            leftover_mask = ~(apo_mask | holo_mask) if (apo_mask is not None and holo_mask is not None) else pd.Series(False, index=df.index)
            leftover_rows = df.loc[leftover_mask].sort_values("pdb_id") if not df.empty else pd.DataFrame()
            if not leftover_rows.empty:
                if n_apo or n_holo:
                    fh.write("\n")
                lo_base = leftover_rows.drop(columns=["variant"], errors="ignore")
                l_cols = [c for c in preferred_cols if c in lo_base.columns] \
                         + [c for c in lo_base.columns if c not in preferred_cols]
                lo_base[l_cols].to_csv(fh, sep="\t", index=False)
                sections.append(("Unlabeled", lo_base[l_cols]))
            dbg("INFO", "summary", f"sections=Apo:{n_apo} Holo:{n_holo} other={len(leftover_rows)}")

            # Pretty summary file (same run-id naming as TSV)
            if args.pretty_summary and sections:
                pretty_name = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary(sections, pretty_name)
                print(f"[eval] pretty_summary_out={pretty_name}")
        else:
            _df = df.drop(columns=["variant"], errors="ignore")
            cols = [c for c in preferred_cols if c in _df.columns] \
                   + [c for c in _df.columns if c not in preferred_cols]
            _df[cols].to_csv(fh, sep="\t", index=False)
            # Pretty summary alongside TSV
            if args.pretty_summary:  # if you added the flag; otherwise remove the 'if'
                pretty_name = summary_path.with_name(summary_path.stem + "_pretty.txt")
            _write_pretty_summary([(None, _df[cols])], pretty_name)
            print(f"[eval] pretty_summary_out={pretty_name}")
    metric_cols = [c for c in df.columns if c not in excluded_cols]
    macro = df[metric_cols].mean(numeric_only=True).to_dict()
    macro_df = pd.DataFrame([{"pdb_id": "macro_avg", **{k: macro[k] for k in metric_cols}}])
    macro_path = analysis_root / "summary_macro.tsv"
    macro_df.to_csv(macro_path, sep="\t", index=False)

    if getattr(args, "emit_control_report", False):
        long_columns = [
            "pdb_id",
            "run_id",
            "target_name",
            "library_name",
            "controls_found",
            "control_id",
            "rmsd_to_crystal_A",
            "redock_best_energy_kcal_mol",
            "chosen",
            "max_spread_A",
            "policy",
        ]
        report_df = pd.DataFrame(control_long_records, columns=long_columns)
        control_report_path = analysis_root / "control_redock_report.tsv"
        report_df.to_csv(control_report_path, sep="\t", index=False)

        control_summary_df = pd.DataFrame(control_summary_records, columns=long_columns)
        control_summary_path = analysis_root / "control_redock_summary.tsv"
        control_summary_df.to_csv(control_summary_path, sep="\t", index=False)

        dbg("INFO", "control", f"report_rows={len(report_df)} out={control_report_path}")
        dbg("INFO", "control", f"summary_rows={len(control_summary_df)} out={control_summary_path}")
        # Pretty control summary (no schema/value changes)
        if getattr(args, "pretty_summary", True):
            control_pretty_path = control_summary_path.with_name(control_summary_path.stem + "_pretty.txt")
            _write_pretty_table_noformat(control_summary_df, control_pretty_path, title="Control Redock Summary")
            print(f"[control] pretty_summary_out={control_pretty_path}")

    total_targets = len(control_targets)
    total_scanned = len(scanned_ok)
    total_missing = len(missing_logs)
    scanned_label = ",".join(scanned_ok) if scanned_ok else "(none)"
    missing_label = ",".join(missing_logs) if missing_logs else "(none)"
    total_redock = sum(stats.get("redock_lines", 0) for stats in parsed_counts.values())
    total_centers = sum(stats.get("centers", 0) for stats in parsed_counts.values())
    total_reports = sum(stats.get("report_rows", 0) for stats in parsed_counts.values())

    dbg("INFO", "summary", f"n_targets_total={total_targets}")
    dbg(
        "INFO",
        "summary",
        f"n_logs_scanned={total_scanned} n_logs_missing={total_missing}",
    )
    dbg("INFO", "summary", f"scanned_ok={scanned_label}")
    dbg("INFO", "summary", f"missing_logs={missing_label}")
    dbg(
        "INFO",
        "summary",
        (
            "totals: "
            f"redock_lines={total_redock} "
            f"control_centers_found={total_centers} "
            f"report_rows={total_reports}"
        ),
    )

    scan_summary_path = analysis_root / "control_redock_scan_summary.txt"
    per_target_rows: List[Dict[str, int | str]] = []
    for pdb_id in control_targets:
        stats = parsed_counts.get(pdb_id, {"centers": 0, "redock_lines": 0, "report_rows": 0})
        per_target_rows.append(
            {
                "pdb_id": pdb_id,
                "centers_found": stats.get("centers", 0),
                "redock_lines": stats.get("redock_lines", 0),
                "report_rows": stats.get("report_rows", 0),
            }
        )
    with open(scan_summary_path, "w", encoding="utf-8") as fh:
        fh.write("Control Log Scan Summary\n")
        fh.write(f"n_targets_total: {total_targets}\n")
        fh.write(f"n_logs_scanned: {total_scanned}\n")
        fh.write(f"n_logs_missing: {total_missing}\n")
        fh.write(f"scanned_ok: {scanned_label}\n")
        fh.write(f"missing_logs: {missing_label}\n")
        fh.write(
            (
                "totals: "
                f"redock_lines={total_redock} "
                f"control_centers_found={total_centers} "
                f"report_rows={total_reports}\n"
            )
        )
        fh.write("\nPer-target counts\n")
        if per_target_rows:
            counts_df = pd.DataFrame(per_target_rows, columns=[
                "pdb_id",
                "centers_found",
                "redock_lines",
                "report_rows",
            ])
            fh.write(counts_df.to_string(index=False))
            fh.write("\n")
        else:
            fh.write("(none)\n")

    dbg("INFO", "summary", f"scan_summary_pretty={scan_summary_path}")

    print(f"[dbg.summary] header={list(_df[cols].columns)}")
    dbg("INFO", "summary", f"out={summary_path} columns={metric_cols}")
    dbg("DEBUG", "summary", f"targets_written={len(df)} macro_path={macro_path}")

if __name__ == "__main__":
    main()
