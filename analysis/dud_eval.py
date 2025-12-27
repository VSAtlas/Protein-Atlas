#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DUD/DUD-E style evaluator for Atlas docking outputs (filename-labeled actives/decoys).

Assumptions:
- Input: docked/<PDB_ID>/docking_score_long.csv (per-POSE rows)
- Lower score = better (Vina-style)
- Ligand filename contains token 'active/actives' or 'decoy/decoys' (case-insensitive)
    ABL1_active_0123.pdbqt  or  ABL1_decoy_0456.pdbqt
- derive labels from the filename

Outputs:
- <OVERALL_DIR>/<out_dir>/<RUN_ID>/<PDB_ID>/
    - metrics.tsv
    - ef_curve.png, roc.png, pr.png, score_hist.png
- <OVERALL_DIR>/<out_dir>/<RUN_ID>/summary.tsv (per-target table)
- <OVERALL_DIR>/<out_dir>/<RUN_ID>/summary_macro.tsv (macro average over targets)

Metrics:
- EF@{1,2,5,10}%
- ROC-AUC, PR-AUC
- logAUC (lambda=1e-3), adjusted logAUC = logAUC - 0.14462
- BEDROC(alpha=20.0)

Notes:
- We "collapse" poses and microstates by taking the minimum score per normalized ligand id.
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
from library_index import LibraryIndex
from run_manifest import get_manifest_paths, _load_manifest  # type: ignore

_LIB_INDEX_CACHE: Dict[Tuple[str, str], LibraryIndex] = {}


def _get_library_index(root: Path, manifest_filename: str) -> LibraryIndex:
    """
    Return a cached LibraryIndex for a given prepped-ligands root + manifest filename.
    This avoids re-reading the same manifest for every target.
    """
    key = (str(root.resolve()), manifest_filename)
    idx = _LIB_INDEX_CACHE.get(key)
    if idx is None:
        idx = LibraryIndex(manifest_filename=manifest_filename)
        idx.load([root])
        _LIB_INDEX_CACHE[key] = idx
    return idx


def _is_probable_run_id_dirname(name: str) -> bool:
    """
    Heuristic to detect run-id style directory names (e.g., 20251205_225826, 2025-12-18T...).
    Avoid matching classic 4-char PDB IDs.
    """
    if not name or len(name) <= 6:
        return False
    upper = name.upper()
    if len(upper) == 4 and upper.isalnum():
        return False
    if "_" in name:
        return True
    if re.fullmatch(r"\d{8,}", name):
        return True
    if "T" in name and re.match(r"\d{4}-?\d{2}-?\d{2}", name):
        return True
    return False


def _resolve_scan_roots(docked_root: Path, run_id: Optional[str]) -> List[Path]:
    """
    Determine which roots to scan for docking outputs.

    Priority:
      - If run_id is provided and docked_root/run_id exists, scan that first,
        then fall back to docked_root for legacy layouts.
      - If run_id provided but missing, fall back to docked_root.
      - If no run_id: scan run-like subdirs (mtime desc) if present, then docked_root.
      - If docked_root is not a dir, just return [docked_root].
    """
    if not docked_root.is_dir():
        return [docked_root]

    if run_id:
        candidate = docked_root / run_id
        if candidate.is_dir():
            roots: List[Path] = [candidate]
            if docked_root != candidate:
                roots.append(docked_root)
            return roots
        return [docked_root]

    run_subdirs = [
        p for p in docked_root.iterdir()
        if p.is_dir() and _is_probable_run_id_dirname(p.name)
    ]
    if run_subdirs:
        run_subdirs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return run_subdirs + [docked_root]

    return [docked_root]


@dataclass
class TargetSpec:
    target_key: str
    pdb_id: str
    variant: Optional[str]
    ph_tag: Optional[str]
    csv_path: Optional[Path]
    source: str  # "manifest" | "scan" | "legacy"


def make_target_key(pdb_id: str,
                    variant: Optional[str],
                    ph_tag: Optional[str]) -> str:
    """
    Build a path-safe identifier for a target.
    Backward compatibility: if both variant and pH are blank, return the pdb_id unchanged.
    """
    variant_norm = (variant or "").strip()
    ph_norm = (ph_tag or "").strip()
    if not variant_norm and not ph_norm:
        return pdb_id
    variant_norm = variant_norm.upper().replace(" ", "_")
    ph_norm = ph_norm.replace(" ", "_")
    for sep in (os.sep, os.altsep):
        if sep:
            variant_norm = variant_norm.replace(sep, "_")
            ph_norm = ph_norm.replace(sep, "_")
    return f"{pdb_id}__{variant_norm}__{ph_norm}"


def _resolve_reranked_scorch_path(
    spec: TargetSpec,
    docked_root: Path,
    post_docked_root: Path,
    run_id: Optional[str],
    run_dir: Optional[Path],
) -> Optional[Path]:
    """
    Best-effort resolution of consensus_reranked_scorch.csv for a target.
    Prefers new layout (post_docked/<run_id>/...) but falls back to legacy.
    """
    candidates: List[Path] = []
    seen: Set[str] = set()

    if spec.csv_path:
        try:
            rel = spec.csv_path.parent.relative_to(docked_root)
            cand = post_docked_root / rel / RERANKED_SCORCH_BASENAME
            if str(cand) not in seen:
                candidates.append(cand)
                seen.add(str(cand))
        except Exception:
            pass

    if run_dir and spec.csv_path and _is_under(spec.csv_path, run_dir):
        try:
            rel_new = spec.csv_path.parent.relative_to(run_dir)
            base = post_docked_root / (run_id or "")
            cand = base / rel_new / RERANKED_SCORCH_BASENAME
            if str(cand) not in seen:
                candidates.append(cand)
                seen.add(str(cand))
        except Exception:
            pass

    if run_id:
        base = post_docked_root / run_id
        combos = [
            base / spec.pdb_id / (spec.variant or "") / (spec.ph_tag or "") / RERANKED_SCORCH_BASENAME,
            base / spec.pdb_id / (spec.variant or "") / RERANKED_SCORCH_BASENAME,
            base / spec.pdb_id / RERANKED_SCORCH_BASENAME,
        ]
        for cand in combos:
            if str(cand) not in seen:
                candidates.append(cand)
                seen.add(str(cand))

    legacy_base = post_docked_root
    legacy_candidates = [
        legacy_base / spec.pdb_id / (spec.variant or "") / (spec.ph_tag or "") / RERANKED_SCORCH_BASENAME,
        legacy_base / spec.pdb_id / (spec.variant or "") / RERANKED_SCORCH_BASENAME,
        legacy_base / spec.pdb_id / RERANKED_SCORCH_BASENAME,
    ]
    for cand in legacy_candidates:
        if str(cand) not in seen:
            candidates.append(cand)
            seen.add(str(cand))

    for cand in candidates:
        if cand.exists() and cand.is_file():
            return cand
    return None



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


def _normalize_pdb_ids(tokens: Optional[Iterable[str]]) -> List[str]:
    """
    Normalize tokens to uppercase 4-character PDB IDs.
    Strips extensions and non-alphanumerics; drops too-short tokens.
    """
    if tokens is None:
        return []
    normalized: List[str] = []
    for raw in tokens:
        if raw is None:
            continue
        text = str(raw).strip()
        if not text:
            continue
        text = re.sub(r"\.pdb(?:\.gz)?$", "", text, flags=re.IGNORECASE)
        alnum = re.sub(r"[^A-Za-z0-9]", "", text).upper()
        if len(alnum) < 4:
            dbg("WARN", "pdb-filter", f"token={text} reason=too_short")
            continue
        normalized.append(alnum[:4])
    return _dedup(normalized)


def _compute_analysis_root(args: argparse.Namespace, cfg: Optional[dict], run_id_override: Optional[str] = None) -> Path:
    """
    Resolve the final analysis root honoring OVERALL_DIR, --out-dir, and --run-id.
    Relative out_dir values are anchored under OVERALL_DIR when available.
    """
    out_root = Path(args.out_dir)
    analysis_root = out_root
    if cfg and "OVERALL_DIR" in cfg:
        base = Path(cfg["OVERALL_DIR"])
        analysis_root = out_root if out_root.is_absolute() else base / out_root
    run_id_val = run_id_override if run_id_override is not None else getattr(args, "run_id", None)
    if run_id_val:
        analysis_root = analysis_root / str(run_id_val)
    analysis_root.mkdir(parents=True, exist_ok=True)
    return analysis_root

def _format_run_label(run_id: Optional[str]) -> str:
    """
    Sanitize run_id for filenames. Falls back to 'none' when run_id is missing.
    """
    if not run_id:
        return "none"
    label = str(run_id).strip()
    for sep in (os.sep, os.altsep):
        if sep:
            label = label.replace(sep, "_")
    return label.replace(" ", "") or "none"

def _normalize_run_id_token(val: str) -> str:
    """Normalize run_id tokens for comparison."""
    if val is None:
        return ""
    text = str(val).strip().strip("'\"")
    return re.sub(r"[^0-9A-Za-z]+", "", text).lower()

def _is_under(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except Exception:
        return False

def filter_df_by_run_id(df: pd.DataFrame,
                        cli_run_id: Optional[str],
                        *,
                        pdb_id: str,
                        csv_path: Path,
                        strict: bool,
                        layout_label: str = "legacy") -> tuple[pd.DataFrame, str]:
    """
    Filter a dataframe by run_id with robust normalization.

    Returns (filtered_df, status_reason).
    status_reason is "ok" on match, or a mismatch/fallback label otherwise.
    """
    if not cli_run_id:
        return df, "ok"

    run_col = None
    for cand in RUN_ID_COL_CANDIDATES:
        if cand in df.columns:
            run_col = cand
            break
    if run_col is None:
        dbg("WARN", "run.filter", f"pdb={pdb_id} path={csv_path} reason=no_run_id_col")
        return df, "ok_no_run_id_col"

    target_norm = _normalize_run_id_token(cli_run_id)
    series_raw = df[run_col].astype(str).str.strip().str.strip("'\"")
    series_norm = series_raw.str.replace(r"[^0-9A-Za-z]+", "", regex=True).str.lower()
    mask = series_norm == target_norm
    kept = int(mask.sum())
    total = len(df)
    if kept > 0:
        dbg("INFO", "run.filter",
            f"pdb={pdb_id} path={csv_path} run_id={cli_run_id} kept={kept}/{total} col={run_col} layout={layout_label}")
        return df.loc[mask].copy(), "ok"

    value_counts = series_norm.value_counts(dropna=False)
    sample = value_counts.head(5).to_dict()
    dbg("WARN", "run.filter",
        f"pdb={pdb_id} path={csv_path} run_id={cli_run_id} matched=0/{total} layout={layout_label} sample={sample}")

    if strict:
        dbg("ERROR", "run.filter",
            f"pdb={pdb_id} path={csv_path} run_id={cli_run_id} matched=0 action=abort layout={layout_label}")
        return df.loc[mask].copy(), "run_id_mismatch_new_layout"

    if len(value_counts) == 1:
        dbg("WARN", "run.filter",
            f"pdb={pdb_id} layout=legacy_fallback action=use_singleton run_id_val={value_counts.index[0]} rows={total}")
        return df.copy(), "run_id_mismatch_legacy_fallback"

    top_value = value_counts.index[0]
    top_rows = int(value_counts.iloc[0])
    dbg("WARN", "run.filter",
        f"pdb={pdb_id} layout=legacy_fallback action=use_top run_id_val={top_value} rows={top_rows}/{total}")
    top_mask = series_norm == top_value
    return df.loc[top_mask].copy(), "run_id_mismatch_legacy_fallback"


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


def _infer_library_name_via_manifest(
    target_id: str,
    ligand_basenames: Set[str],
    roots: List[Path],
    cfg: Dict,
) -> str:
    """
    Fast path: infer library name using a pre-built manifest under prepped_ligands.

    Strategy:
      - For each existing root that has a manifest file (default: _manifest.json),
        load a LibraryIndex for that root.
      - For a bounded sample of ligand basenames, look up each basename in the
        manifest via LibraryIndex.lookup_filename.
      - Each hit yields a relative path like "abl1/ABL1_active_10.pdbqt";
        treat the first path component ("abl1") as the library name and count votes.
      - If there is a unique best library with >0 votes, return it.
      - Otherwise, return "" so the caller can fall back to the legacy scan.
    """
    if not ligand_basenames:
        return ""

    manifest_filename = cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json")
    if not manifest_filename:
        manifest_filename = "_manifest.json"
    manifest_filename = str(manifest_filename)

    sample_names = sorted(ligand_basenames)
    max_sample = 200
    if len(sample_names) > max_sample:
        sample_names = sample_names[:max_sample]

    overall_counts: Dict[str, int] = {}

    for root in roots:
        if not root or not root.exists():
            continue
        manifest_path = root / manifest_filename
        if not manifest_path.exists():
            dbg("DEBUG", "library",
                f"pdb={target_id} manifest_missing={manifest_path}")
            continue

        try:
            idx = _get_library_index(root, manifest_filename)
        except Exception as exc:
            dbg("WARN", "library",
                f"pdb={target_id} manifest_load_failed root={root} err={exc}")
            continue

        for name in sample_names:
            base = os.path.basename(name)
            if not base:
                continue
            hit = idx.lookup_filename(base, roots=[root])
            if not hit:
                continue
            try:
                rel = hit.relative_to(root)
            except Exception:
                rel = hit
            parts = rel.as_posix().split("/")
            if not parts:
                continue
            lib_name = parts[0].strip()
            if not lib_name:
                continue
            overall_counts[lib_name] = overall_counts.get(lib_name, 0) + 1

    if not overall_counts:
        return ""

    best_lib, best_count = max(overall_counts.items(), key=lambda kv: (kv[1], kv[0]))
    n_best = sum(1 for c in overall_counts.values() if c == best_count)
    if best_count <= 0 or n_best != 1:
        dbg("WARN", "library",
            f"pdb={target_id} manifest_votes ambiguous best_count={best_count} "
            f"candidates={len(overall_counts)}")
        return ""

    dbg("INFO", "library",
        f"pdb={target_id} picked='{best_lib}' method='manifest_entries' "
        f"matches={best_count}")
    return best_lib


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
      1) Global prepped-ligands manifest via LibraryIndex (fast path).
      2) manifest.json at library root with "library_name"/"name"/"library"/"label"
      3) filename-overlap between docked basenames and files in <library> (non-recursive; optional shallow)
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

    # 1) Fast path: try prepped-ligands manifest via LibraryIndex (if present).
    manifest_choice = _infer_library_name_via_manifest(
        target_id,
        ligand_basenames,
        roots,
        cfg,
    )
    if manifest_choice:
        return manifest_choice

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
    status_reason: str = "ok"


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
VALID_COL_CANDIDATES = [
    "valid",
    "pose_valid",
    "is_valid",
    "passed_pose_validation",
    "pose_is_valid",
    "valid_pose",
    "pose_ok",
    "ok",
    "passed",
    "pass",
]
VALID_TRUE_STRINGS = {"true", "t", "yes", "y", "pass", "passed", "ok", "valid"}
# We now support both DUD-specific and legacy score CSV names.
# Prefer the DUD-prefixed name when both exist.
CSV_BASENAMES = (
    "dud_docking_score_long.csv",  # new DUD runs (preferred)
    "docking_score_long.csv",      # legacy name (fallback)
)
CONSENSUS_CSV_BASENAME = "consensus_docking_scores.csv"
RERANKED_SCORCH_BASENAME = "consensus_reranked_scorch.csv"
CONSENSUS_SCORE_CANDIDATES = ["consensus_score", "score", "consensus", "final_score"]
RUN_ID_COL_CANDIDATES = ["run_id", "runid", "run"]

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

def guess_consensus_score_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    lower = {c.lower(): c for c in df.columns}
    for key in ("consensus_score", "score"):
        if key in lower:
            return lower[key]
    for cand in CONSENSUS_SCORE_CANDIDATES:
        if cand in lower:
            return lower[cand]
    col = guess_col(df, CONSENSUS_SCORE_CANDIDATES)
    if col:
        return col
    return guess_score_col(df, override=None)


def guess_reranked_scorch_score_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    for cand in ("scorch_composite", "SCORCH_score_used", "final_rank", "consensus_score"):
        if cand in df.columns:
            return cand
    lower = {c.lower(): c for c in df.columns}
    for cand in ("scorch_composite", "scorch_score_used", "final_rank", "consensus_score"):
        if cand in lower:
            return lower[cand]
    return guess_consensus_score_col(df, None)


def read_reranked_scorch_csv(path: Path) -> pd.DataFrame:
    """
    Read a reranked SCORCH CSV, tolerating optional metadata preamble lines.

    The file may start with comment-style metadata (e.g., '# run_id=...').
    We first attempt a normal read; if required columns are missing, re-read
    after skipping metadata lines.
    """
    def _detect_preamble(lines: List[str]) -> int:
        count = 0
        meta_prefixes = ("run_id", "pdb_id", "variant", "ph", "ph_label")
        for line in lines:
            stripped = line.strip()
            if not stripped:
                count += 1
                continue
            if stripped.startswith("#"):
                count += 1
                continue
            low = stripped.lower()
            if any(low.startswith(f"{p}=") or low.startswith(f"{p}:") for p in meta_prefixes):
                count += 1
                continue
            break
        return count

    try:
        return pd.read_csv(path)
    except Exception:
        pass

    try:
        with path.open("r", encoding="utf-8") as fh:
            head: List[str] = []
            for _ in range(5):
                try:
                    head.append(next(fh))
                except StopIteration:
                    break
    except Exception:
        return pd.read_csv(path)

    skip = _detect_preamble(head)
    if skip <= 0:
        return pd.read_csv(path)
    return pd.read_csv(path, skiprows=skip)

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

def guess_consensus_ligfile_col(df: pd.DataFrame, override: Optional[str]) -> str:
    if override and override in df.columns:
        return override
    lower = {c.lower(): c for c in df.columns}
    for key in ("ligand", "ligand_file"):
        if key in lower:
            return lower[key]
    return guess_ligfile_col(df, override=None)

def resolve_valid_col(df: pd.DataFrame, override: Optional[str], enabled: bool) -> Optional[str]:
    if not enabled:
        return None
    if override:
        if override in df.columns:
            return override
        raise ValueError(
            f"Valid-only requested but --valid-col '{override}' is missing. "
            "Pass --valid-col with an existing column name."
        )
    for name in VALID_COL_CANDIDATES:
        if name in df.columns:
            return name
    raise ValueError("Valid-only requested but no validity column found. Use --valid-col.")

def parse_valid_mask(series: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        numeric = pd.to_numeric(series, errors="coerce")
        return numeric.notna() & (numeric != 0)
    text = series.fillna("").astype(str).str.strip()
    numeric = pd.to_numeric(text, errors="coerce")
    numeric_mask = numeric.notna() & (numeric != 0)
    token_mask = text.str.lower().isin(VALID_TRUE_STRINGS)
    return numeric_mask | token_mask

def _filter_consensus_no_data(df: pd.DataFrame, score_col: str) -> tuple[pd.DataFrame, int]:
    """
    Drop rows that clearly have no engine data.

    Priority:
      1) n_engines_with_data > 0
      2) Boolean-like validity column (valid/ok/has_data)
      3) All engine columns zero (excluding consensus_score + obvious metadata)
    """
    start = len(df)
    if start == 0:
        return df, 0

    if "n_engines_with_data" in df.columns:
        counts = pd.to_numeric(df["n_engines_with_data"], errors="coerce").fillna(0)
        mask = counts > 0
        filtered = df.loc[mask].copy()
        dropped = start - len(filtered)
        dbg("INFO", "consensus.filter", f"method=n_engines_with_data kept={len(filtered)}/{start}")
        return filtered, dropped

    bool_candidates = [c for c in df.columns if c.lower() in {"valid", "ok", "has_data", "usable", "is_valid", "available"}]
    for col in bool_candidates:
        try:
            mask = parse_valid_mask(df[col])
        except Exception:
            continue
        filtered = df.loc[mask].copy()
        dropped = start - len(filtered)
        dbg("INFO", "consensus.filter", f"method=bool_col col={col} kept={len(filtered)}/{start}")
        return filtered, dropped

    numeric_cols = list(df.select_dtypes(include=["number", "bool"]).columns)
    exclude_tokens = {"rank", "percent", "perc", "stage", "pose", "order"}
    engine_cols: List[str] = []
    for col in numeric_cols:
        low = col.lower()
        if low == score_col.lower():
            continue
        if low in {"n_engines_with_data", "run_id", "is_active"}:
            continue
        if any(tok in low for tok in exclude_tokens):
            continue
        engine_cols.append(col)

    if engine_cols:
        engines = df[engine_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        consensus_scores = pd.to_numeric(df[score_col], errors="coerce").fillna(0)
        mask_zero = (engines == 0).all(axis=1) & (consensus_scores == 0)
        filtered = df.loc[~mask_zero].copy()
        dropped = int(mask_zero.sum())
        if dropped > 0:
            dbg("INFO", "consensus.filter", f"method=engine_zero engine_cols={engine_cols} dropped={dropped} kept={len(filtered)}/{start}")
        else:
            dbg("DEBUG", "consensus.filter", f"method=engine_zero engine_cols={engine_cols} dropped=0 kept={len(filtered)}/{start}")
        return filtered, dropped

    dbg("WARN", "consensus.filter", "method=engine_zero reason=no_numeric_candidates")
    return df, 0

# >>> RUN-SELECTION START
def select_default_run_id(targets: List[TargetSpec],
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

    for target in targets:
        csv_path = csv_paths_by_target.get(target.target_key)
        if not csv_path or not csv_path.exists():
            continue
        try:
            df = pd.read_csv(csv_path)
        except Exception as exc:
            dbg("WARN", "run", f"pdb={target.pdb_id} auto_select_read_err={exc}")
            continue

        if "run_id" not in df.columns:
            continue

        try:
            lig_col = guess_ligfile_col(df, lig_col_cli)
        except Exception as exc:
            dbg("WARN", "run", f"pdb={target.pdb_id} auto_select_ligcol_err={exc}")
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
                run_full_map.setdefault(run_label, set()).add(target.target_key)
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
MICROSTATE_TAIL_RE = re.compile(r'(?:__ms|_ms)_[A-Za-z0-9]+$', re.IGNORECASE)
EXT_RE = re.compile(r'\.(pdbqt|sdf|mol2)(?:\.gz)?$', re.IGNORECASE)

def parse_name_and_label(path_str: str) -> Tuple[str, Optional[int]]:
    """Return (ligand_root_id, is_active) from a path or filename.
       is_active: 1 for active, 0 for decoy, None if not found."""
    base = os.path.basename(str(path_str))
    stem = EXT_RE.sub("", base)               # drop .pdbqt/.sdf/.mol2(.gz)
    stem = POSE_TAIL_RE.sub("", stem)         # drop trailing _pose1 etc.
    stem = MICROSTATE_TAIL_RE.sub("", stem)   # drop trailing microstate token
    lig_id = stem                              # normalized ligand identifier (pose+microstate collapsed)
    m = TOKEN_RE.search(base)
    if not m:
        return lig_id, None
    tok = m.group(1).lower()
    return lig_id, 1 if tok == "active" else 0


def compute_decoy_stats_from_long_csv(
    csv_path: Path | str,
    lig_col_override: Optional[str] = None,
    score_col_override: Optional[str] = None,
) -> tuple[float, float, int]:
    """
    Compute (mean, std, n) of best decoy scores from a docking_score_long.csv-style file.

    - csv_path: path to the CSV file (per-POSE rows).
    - lig_col_override / score_col_override: optional explicit column names.

    Decoy/active labels are derived from the ligand filename using parse_name_and_label,
    and best scores are defined as the minimum score per ligand.

    Returns (mean_decoy_score, std_decoy_score, n_decoys).
    If no decoys are found or the file is missing, returns (nan, nan, 0).
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return float("nan"), float("nan"), 0

    df = pd.read_csv(csv_path)

    lig_col = guess_ligfile_col(df, lig_col_override)
    score_col = guess_score_col(df, score_col_override)

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    lig_series = df[lig_col].astype(str)
    parsed = lig_series.apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])

    best = df.groupby("lig_id", as_index=False).agg(
        best_score=(score_col, "min"),
        is_active=("is_active", "max"),
    )

    decoys = best[best["is_active"] == 0]
    n_decoys = int(len(decoys))
    if n_decoys == 0:
        return float("nan"), float("nan"), 0

    mu = float(decoys["best_score"].mean())
    sigma = float(decoys["best_score"].std(ddof=1))
    return mu, sigma, n_decoys

def _make_placeholder_row(pdb_id: str,
                          variant: Optional[str],
                          ph_tag: Optional[str],
                          run_id: Optional[str],
                          bedroc_alpha: float,
                          status_reason: str) -> pd.Series:
    base_row = {
        "run_id": str(run_id) if run_id else "(none)",
        "pdb_id": pdb_id,
        "variant": variant or "",
        "pH": ph_tag or "",
        "N": 0,
        "n_actives": 0,
        "actives_fraction": float("nan"),
        "ROC_AUC": float("nan"),
        "PR_AUC": float("nan"),
        "logAUC": float("nan"),
        "logAUC_adj": float("nan"),
        f"BEDROC_alpha_{bedroc_alpha:g}": float("nan"),
        "EF@1%": float("nan"),
        "EF@2%": float("nan"),
        "EF@5%": float("nan"),
        "EF@10%": float("nan"),
        "status_reason": status_reason,
    }
    return pd.Series(base_row)

# --- aggregation helpers ---

def _collapse_best_scores(df: pd.DataFrame, score_col: str, *, best_is_min: bool = True) -> tuple[pd.DataFrame, int]:
    """
    Collapse to one row per lig_id using either min or max score.

    Returns (best_df, drop_nan_best) where drop_nan_best counts ligands removed
    due to NaN best_score.
    """
    agg_fn = "min" if best_is_min else "max"
    best = df.groupby("lig_id", as_index=False).agg(
        best_score=(score_col, agg_fn),
        is_active=("is_active", "max"),   # any active -> active
    )
    before_best = len(best)
    best = best.dropna(subset=["best_score"])
    drop_nan_best = before_best - len(best)
    return best, drop_nan_best

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

def _compute_metrics_and_plots(
    *,
    pdb_id: str,
    best: pd.DataFrame,
    out_dir: Path,
    bedroc_alpha: float,
    logauc_lambda: float,
    score_high_is_better: bool,
    missing_name_detected: int,
    drop_nonfinite: int,
    dropped_token: int,
    drop_nan_best: int,
    token_regex: str,
    variant_label: Optional[str],
    run_id: Optional[str],
    ligand_basenames: Set[str],
    has_run_id_col: bool,
    run_ids_present: Set[str],
    hist_xlabel: str,
    title_suffix: str = "",
    mode_label: str = "docking",
    status_reason: str = "ok",
) -> TargetEvaluation:
    y_true = best["is_active"].astype(int).to_numpy()
    if score_high_is_better:
        y_high = best["best_score"].to_numpy()
        y_low = -y_high
    else:
        y_low = best["best_score"].to_numpy()
        y_high = -y_low

    N = len(best); n_act = int(y_true.sum())
    dbg("DEBUG", "screen",
        f"pdb={pdb_id} mode={mode_label} missing_name_detected={missing_name_detected} kept_rows={len(best)} ligands={N}")
    if N == 0 or n_act == 0 or n_act == N:
        dbg("WARN", "metrics", f"pdb={pdb_id} mode={mode_label} degenerate_set N={N} actives={n_act}")
        reason = status_reason if status_reason != "ok" else "degenerate_set"
        return TargetEvaluation(metrics=None,
                                ligand_basenames=ligand_basenames,
                                has_run_id_column=has_run_id_col,
                                run_ids=run_ids_present,
                                status_reason=reason)

    ef = ef_at_fractions(y_true, y_low, fractions=(0.01,0.02,0.05,0.10))
    dbg("DEBUG", "metrics", f"pdb={pdb_id} mode={mode_label} start N={N} n_actives={n_act}")
    try:
        rocAUC = float(roc_auc_score(y_true, y_high))
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} mode={mode_label} metric=ROC_AUC err={exc}")
        rocAUC = float("nan")
    try:
        fpr, tpr, _ = roc_curve(y_true, y_high)
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} mode={mode_label} metric=ROC_curve err={exc}")
        fpr = np.array([0.0, 1.0])
        tpr = np.array([0.0, 1.0])
    try:
        prAUC, precision, recall = pr_auc(y_true, y_high)
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} mode={mode_label} metric=PR_AUC err={exc}")
        prAUC = float("nan")
        precision = np.array([0.0, 1.0])
        recall = np.array([0.0, 1.0])
    try:
        lAUC = log_auc_from_roc(fpr, tpr, lam=logauc_lambda)
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} mode={mode_label} metric=logAUC err={exc}")
        lAUC = float("nan")
    lAUC_adj = lAUC - 0.14462
    try:
        bed = bedroc(y_true, y_high, alpha=bedroc_alpha)
    except Exception as exc:
        dbg("WARN", "metrics", f"pdb={pdb_id} mode={mode_label} metric=BEDROC err={exc}")
        bed = float("nan")

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
    plt.title(f"{pdb_id} - EF curve{title_suffix}")
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
    plt.title(f"{pdb_id} - ROC{title_suffix}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "roc.png", dpi=200)
    plt.close()

    # PR
    base = n_act / N
    dbg("INFO", "screen",
        f"pdb={pdb_id} mode={mode_label} drop_nonfinite={drop_nonfinite} missing_name_detected={missing_name_detected} drop_missing_token={dropped_token} drop_nan_best_score={drop_nan_best} kept_ligands={N} actives={n_act} active_fraction={base:.3f} token_regex='{token_regex}'")
    plt.figure()
    plt.plot(recall, precision, label=f"PR-AUC={prAUC:.3f}")
    plt.hlines(base, 0, 1, colors="k", linestyles="--", linewidth=1, label=f"Baseline={base:.3f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"{pdb_id} - Precision–Recall{title_suffix}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "pr.png", dpi=200)
    plt.close()

    # Score distributions (always plotted as low-is-better values)
    mask_active = best["is_active"] == 1
    mask_decoy = best["is_active"] == 0
    plt.figure()
    plt.hist(y_low[mask_active.to_numpy()], bins=40, alpha=1, label="Actives", zorder=2)
    plt.hist(y_low[mask_decoy.to_numpy()], bins=40, alpha=0.5, label="Decoys", zorder=1)
    plt.xlabel(hist_xlabel)
    plt.ylabel("Count")
    plt.title(f"{pdb_id} - Score distributions{title_suffix}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "score_hist.png", dpi=200)
    plt.close()

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
    row["status_reason"] = status_reason if status_reason else "ok"
    if variant_label:
        row["variant"] = variant_label
    pd.DataFrame([row]).to_csv(out_dir / "metrics.tsv", sep="\t", index=False)
    dbg("DEBUG", "metrics",
        f"pdb={pdb_id} mode={mode_label} ROC_AUC={rocAUC:.3f} PR_AUC={prAUC:.3f} BEDROC={bed:.3f}")
    return TargetEvaluation(metrics=pd.Series(row),
                            ligand_basenames=ligand_basenames,
                            has_run_id_column=has_run_id_col,
                            run_ids=run_ids_present,
                            status_reason=row["status_reason"])

# --- evaluation ---

def evaluate_target(pdb_id: str,
                    csv_path: Path,
                    out_dir: Path,
                    lig_col_cli: Optional[str],
                    score_col_cli: Optional[str],
                    bedroc_alpha: float,
                    logauc_lambda: float,
                    run_id: Optional[str],
                    valid_only: bool = False,
                    valid_col_cli: Optional[str] = None,
                    new_layout_active: bool = False,
                    layout_label: str = "legacy") -> Optional[TargetEvaluation]:
    status_reason = "ok"
    try:
        df = pd.read_csv(csv_path)
        dbg("INFO", "csv", f"pdb={pdb_id} path={csv_path} rows={len(df)}")
    except Exception as e:
        dbg("ERROR", "csv", f"pdb={pdb_id} path={csv_path} err={e}")
        return TargetEvaluation(metrics=None,
                                ligand_basenames=set(),
                                has_run_id_column=False,
                                run_ids=set(),
                                status_reason="read_error")

    has_run_id_col = "run_id" in df.columns
    df, filter_reason = filter_df_by_run_id(df, run_id,
                                            pdb_id=pdb_id,
                                            csv_path=csv_path,
                                            strict=new_layout_active,
                                            layout_label=layout_label)
    status_reason = filter_reason if filter_reason != "ok" else "ok"
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(metrics=None,
                                ligand_basenames=set(),
                                has_run_id_column=has_run_id_col,
                                run_ids=set(),
                                status_reason=final_reason)

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
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(metrics=None,
                                ligand_basenames=set(),
                                has_run_id_column=has_run_id_col,
                                run_ids=set(),
                                status_reason=final_reason)

    # Apply pose validity filtering before best-pose aggregation.
    if valid_only:
        valid_col = resolve_valid_col(df, valid_col_cli, enabled=True)
        before_valid = len(df)
        valid_mask = parse_valid_mask(df[valid_col])
        df = df.loc[valid_mask].copy()
        dropped_valid = before_valid - len(df)
        dbg("INFO", "filter", f"pdb={pdb_id} valid_only=ON valid_col={valid_col} rows_kept={len(df)} rows_dropped={dropped_valid}")

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
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(metrics=None,
                                ligand_basenames=used_ligand_basenames,
                                has_run_id_column=has_run_id_col,
                                run_ids=run_ids_present,
                                status_reason=final_reason)

    best, drop_nan_best = _collapse_best_scores(df, score_col, best_is_min=True)
    missing_name_detected = missing_names + empty_names

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

    return _compute_metrics_and_plots(
        pdb_id=pdb_id,
        best=best,
        out_dir=out_dir,
        bedroc_alpha=bedroc_alpha,
        logauc_lambda=logauc_lambda,
        score_high_is_better=False,
        missing_name_detected=missing_name_detected,
        drop_nonfinite=drop_nonfinite,
        dropped_token=dropped_token,
        drop_nan_best=drop_nan_best,
        token_regex="active|decoy",
        variant_label=variant_label,
        run_id=run_id,
        ligand_basenames=used_ligand_basenames,
        has_run_id_col=has_run_id_col,
        run_ids_present=run_ids_present,
        hist_xlabel="Best docking score (lower = better)",
        title_suffix="",
        mode_label="docking",
        status_reason=status_reason,
    )

def evaluate_target_consensus(
    pdb_id: str,
    csv_path: Path,
    out_dir: Path,
    lig_col_cli: Optional[str],
    score_col_cli: Optional[str],
    bedroc_alpha: float,
    logauc_lambda: float,
    run_id: Optional[str],
    variant_hint: Optional[str] = None,
    ph_tag: Optional[str] = None,  # kept for symmetry/future use
    new_layout_active: bool = False,
    layout_label: str = "legacy",
) -> Optional[TargetEvaluation]:
    status_reason = "ok"
    try:
        df = pd.read_csv(csv_path)
        dbg("INFO", "consensus.csv", f"pdb={pdb_id} path={csv_path} rows={len(df)}")
    except Exception as e:
        dbg("ERROR", "consensus.csv", f"pdb={pdb_id} path={csv_path} err={e}")
        return TargetEvaluation(metrics=None,
                                ligand_basenames=set(),
                                has_run_id_column=False,
                                run_ids=set(),
                                status_reason="read_error")

    has_run_id_col = "run_id" in df.columns
    df, filter_reason = filter_df_by_run_id(df, run_id,
                                            pdb_id=pdb_id,
                                            csv_path=csv_path,
                                            strict=new_layout_active,
                                            layout_label=layout_label)
    status_reason = filter_reason if filter_reason != "ok" else "ok"
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(metrics=None,
                                ligand_basenames=set(),
                                has_run_id_column=has_run_id_col,
                                run_ids=set(),
                                status_reason=final_reason)

    lig_col = guess_consensus_ligfile_col(df, lig_col_cli)
    score_col = guess_consensus_score_col(df, score_col_cli)
    dbg("DEBUG", "consensus.schema",
        f"pdb={pdb_id} ligand_col={lig_col} score_col={score_col} cols={list(df.columns)}")

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df, dropped_no_data = _filter_consensus_no_data(df, score_col)

    before_nf = len(df)
    df = df.dropna(subset=[score_col])
    drop_nonfinite = before_nf - len(df)
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(metrics=None,
                                ligand_basenames=set(),
                                has_run_id_column=has_run_id_col,
                                run_ids=set(),
                                status_reason=final_reason)

    raw_lig = df[lig_col]
    missing_names = int(raw_lig.isna().sum())
    empty_names = int(raw_lig.astype(str).str.strip().eq("").sum())
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

    before = len(df)
    df = df.dropna(subset=["is_active"])
    dropped_token = before - len(df)
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(metrics=None,
                                ligand_basenames=used_ligand_basenames,
                                has_run_id_column=has_run_id_col,
                                run_ids=run_ids_present,
                                status_reason=final_reason)

    best, drop_nan_best = _collapse_best_scores(df, score_col, best_is_min=False)
    missing_name_detected = missing_names + empty_names

    variant_label: Optional[str] = variant_hint.upper() if variant_hint else None
    if variant_label is None and "variant" in df.columns:
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
                dbg("DEBUG", "consensus.variant", f"pdb={pdb_id} unique={unique_vals} counts={counts}")

    dbg("DEBUG", "consensus.filter",
        f"pdb={pdb_id} drop_no_data={dropped_no_data} drop_nonfinite={drop_nonfinite} drop_missing_token={dropped_token} drop_nan_best={drop_nan_best}")

    return _compute_metrics_and_plots(
        pdb_id=pdb_id,
        best=best,
        out_dir=out_dir,
        bedroc_alpha=bedroc_alpha,
        logauc_lambda=logauc_lambda,
        score_high_is_better=True,
        missing_name_detected=missing_name_detected,
        drop_nonfinite=drop_nonfinite,
        dropped_token=dropped_token,
        drop_nan_best=drop_nan_best,
        token_regex="active|decoy",
        variant_label=variant_label,
        run_id=run_id,
        ligand_basenames=used_ligand_basenames,
        has_run_id_col=has_run_id_col,
        run_ids_present=run_ids_present,
        hist_xlabel="Consensus score (higher = better; negated for plot)",
        title_suffix=" (consensus)",
        mode_label="consensus",
        status_reason=status_reason,
    )


def evaluate_target_post_docked_reranked_scorch(
    pdb_id: str,
    csv_path: Path,
    out_dir: Path,
    lig_col_cli: Optional[str],
    score_col_cli: Optional[str],
    bedroc_alpha: float,
    logauc_lambda: float,
    run_id: Optional[str],
    variant_hint: Optional[str] = None,
    ph_tag: Optional[str] = None,
    new_layout_active: bool = False,
    layout_label: str = "post_docked",
) -> Optional[TargetEvaluation]:
    status_reason = "ok"
    try:
        df = read_reranked_scorch_csv(csv_path)
        dbg("INFO", "reranked.csv", f"pdb={pdb_id} path={csv_path} rows={len(df)}")
    except Exception as e:
        dbg("ERROR", "reranked.csv", f"pdb={pdb_id} path={csv_path} err={e}")
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=False,
            run_ids=set(),
            status_reason="read_error",
        )

    has_run_id_col = "run_id" in df.columns
    df, filter_reason = filter_df_by_run_id(
        df,
        run_id,
        pdb_id=pdb_id,
        csv_path=csv_path,
        strict=new_layout_active,
        layout_label=layout_label,
    )
    status_reason = filter_reason if filter_reason != "ok" else "ok"
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=has_run_id_col,
            run_ids=set(),
            status_reason=final_reason,
        )

    lig_col = guess_consensus_ligfile_col(df, lig_col_cli)
    score_col = guess_reranked_scorch_score_col(df, score_col_cli)
    dbg(
        "DEBUG",
        "reranked.schema",
        f"pdb={pdb_id} ligand_col={lig_col} score_col={score_col} cols={list(df.columns)}",
    )

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    score_eval_col = "_score_eval"
    if score_col.lower() == "final_rank":
        df[score_eval_col] = -pd.to_numeric(df[score_col], errors="coerce")
    else:
        df[score_eval_col] = df[score_col]

    before_nf = len(df)
    df = df.dropna(subset=[score_eval_col])
    drop_nonfinite = before_nf - len(df)
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=has_run_id_col,
            run_ids=set(),
            status_reason=final_reason,
        )

    raw_lig = df[lig_col]
    missing_names = int(raw_lig.isna().sum())
    empty_names = int(raw_lig.astype(str).str.strip().eq("").sum())
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

    before = len(df)
    df = df.dropna(subset=["is_active"])
    dropped_token = before - len(df)
    if df.empty:
        final_reason = status_reason if status_reason != "ok" else "empty_after_filter"
        return TargetEvaluation(
            metrics=None,
            ligand_basenames=used_ligand_basenames,
            has_run_id_column=has_run_id_col,
            run_ids=run_ids_present,
            status_reason=final_reason,
        )

    best, drop_nan_best = _collapse_best_scores(df, score_eval_col, best_is_min=False)
    missing_name_detected = missing_names + empty_names

    variant_label: Optional[str] = variant_hint.upper() if variant_hint else None
    if variant_label is None and "variant" in df.columns:
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
                dbg("DEBUG", "reranked.variant", f"pdb={pdb_id} unique={unique_vals} counts={counts}")

    return _compute_metrics_and_plots(
        pdb_id=pdb_id,
        best=best,
        out_dir=out_dir,
        bedroc_alpha=bedroc_alpha,
        logauc_lambda=logauc_lambda,
        score_high_is_better=True,
        missing_name_detected=missing_name_detected,
        drop_nonfinite=drop_nonfinite,
        dropped_token=dropped_token,
        drop_nan_best=drop_nan_best,
        token_regex="active|decoy",
        variant_label=variant_label,
        run_id=run_id,
        ligand_basenames=used_ligand_basenames,
        has_run_id_col=has_run_id_col,
        run_ids_present=run_ids_present,
        hist_xlabel="Reranked SCORCH score (higher = better)",
        title_suffix=" (reranked_scorch)",
        mode_label="post_docked_scorch",
        status_reason=status_reason,
    )

def emit_consensus_summary(df: pd.DataFrame, analysis_root: Path, run_label: str, pretty_enabled: bool) -> None:
    if df is None or df.empty:
        dbg("INFO", "consensus.summary", "no consensus rows to write")
        return

    preferred_cols = ("variant", "pH", "run_id", "target_name", "library_name", "pdb_id")
    summary_path = analysis_root / f"consensus_summary_{run_label}.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    variant_col_present = "variant" in df.columns
    variant_upper = df["variant"].astype(str).str.upper() if variant_col_present else None
    apo_mask = (variant_upper == "APO") if variant_col_present else None
    holo_mask = (variant_upper == "HOLO") if variant_col_present else None
    has_sections = bool(variant_col_present and ((apo_mask is not None and apo_mask.any()) or (holo_mask is not None and holo_mask.any())))

    with open(summary_path, "w", newline="") as fh:
        sections: List[tuple[str | None, pd.DataFrame]] = []
        if has_sections:
            _df = df.copy()
            cols = [c for c in preferred_cols if c in _df.columns] \
                   + [c for c in _df.columns if c not in preferred_cols]
            apo_rows = df.loc[apo_mask].sort_values("pdb_id") if apo_mask is not None else pd.DataFrame()
            holo_rows = df.loc[holo_mask].sort_values("pdb_id") if holo_mask is not None else pd.DataFrame()
            n_apo = len(apo_rows)
            n_holo = len(holo_rows)

            if n_apo:
                apo_out = apo_rows
                a_cols = [c for c in preferred_cols if c in apo_out.columns] \
                         + [c for c in apo_out.columns if c not in preferred_cols]
                fh.write("Apo\n")
                apo_out[a_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("CONSENSUS - Apo", apo_out[a_cols]))

            if n_holo:
                if n_apo:
                    fh.write("\n")
                holo_out = holo_rows
                h_cols = [c for c in preferred_cols if c in holo_out.columns] \
                         + [c for c in holo_out.columns if c not in preferred_cols]
                fh.write("Holo\n")
                holo_out[h_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("CONSENSUS - Holo", holo_out[h_cols]))

            leftover_mask = ~(apo_mask | holo_mask) if (apo_mask is not None and holo_mask is not None) else pd.Series(False, index=df.index)
            leftover_rows = df.loc[leftover_mask].sort_values("pdb_id") if not df.empty else pd.DataFrame()
            if not leftover_rows.empty:
                if n_apo or n_holo:
                    fh.write("\n")
                lo_base = leftover_rows
                l_cols = [c for c in preferred_cols if c in lo_base.columns] \
                         + [c for c in lo_base.columns if c not in preferred_cols]
                lo_base[l_cols].to_csv(fh, sep="\t", index=False)
                sections.append(("CONSENSUS - Unlabeled", lo_base[l_cols]))

            dbg("INFO", "consensus.summary", f"sections=Apo:{n_apo} Holo:{n_holo} other={len(leftover_rows)} out={summary_path}")
            if pretty_enabled and sections:
                pretty_path = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary(sections, pretty_path)
                print(f"[consensus] pretty_summary_out={pretty_path}")
        else:
            _df = df.copy()
            cols = [c for c in preferred_cols if c in _df.columns] \
                   + [c for c in _df.columns if c not in preferred_cols]
            _df[cols].to_csv(fh, sep="\t", index=False)
            dbg("INFO", "consensus.summary", f"out={summary_path} rows={len(_df)}")
            if pretty_enabled:
                pretty_path = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary([("CONSENSUS", _df[cols])], pretty_path)
                print(f"[consensus] pretty_summary_out={pretty_path}")


def emit_reranked_scorch_summary(df: pd.DataFrame, analysis_root: Path, run_label: str, pretty_enabled: bool) -> None:
    if df is None or df.empty:
        dbg("INFO", "reranked.summary", "no reranked rows to write")
        return

    preferred_cols = ("variant", "pH", "run_id", "target_name", "library_name", "pdb_id")
    summary_root = analysis_root / "post_docked"
    summary_path = summary_root / f"consensus_reranked_scorch_summary_{run_label}.tsv"
    summary_root.mkdir(parents=True, exist_ok=True)

    variant_col_present = "variant" in df.columns
    variant_upper = df["variant"].astype(str).str.upper() if variant_col_present else None
    apo_mask = (variant_upper == "APO") if variant_col_present else None
    holo_mask = (variant_upper == "HOLO") if variant_col_present else None
    has_sections = bool(
        variant_col_present
        and ((apo_mask is not None and apo_mask.any()) or (holo_mask is not None and holo_mask.any()))
    )

    with open(summary_path, "w", newline="") as fh:
        sections: List[tuple[str | None, pd.DataFrame]] = []
        if has_sections:
            apo_rows = df.loc[apo_mask].sort_values("pdb_id") if apo_mask is not None else pd.DataFrame()
            holo_rows = df.loc[holo_mask].sort_values("pdb_id") if holo_mask is not None else pd.DataFrame()
            n_apo = len(apo_rows)
            n_holo = len(holo_rows)

            if n_apo:
                apo_out = apo_rows
                a_cols = [c for c in preferred_cols if c in apo_out.columns] + [c for c in apo_out.columns if c not in preferred_cols]
                fh.write("Apo\n")
                apo_out[a_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("POST_DOCKED_SCORCH - Apo", apo_out[a_cols]))

            if n_holo:
                if n_apo:
                    fh.write("\n")
                holo_out = holo_rows
                h_cols = [c for c in preferred_cols if c in holo_out.columns] + [c for c in holo_out.columns if c not in preferred_cols]
                fh.write("Holo\n")
                holo_out[h_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("POST_DOCKED_SCORCH - Holo", holo_out[h_cols]))

            leftover_mask = ~(apo_mask | holo_mask) if (apo_mask is not None and holo_mask is not None) else pd.Series(False, index=df.index)
            leftover_rows = df.loc[leftover_mask].sort_values("pdb_id") if not df.empty else pd.DataFrame()
            if not leftover_rows.empty:
                if n_apo or n_holo:
                    fh.write("\n")
                lo_base = leftover_rows
                l_cols = [c for c in preferred_cols if c in lo_base.columns] + [c for c in lo_base.columns if c not in preferred_cols]
                lo_base[l_cols].to_csv(fh, sep="\t", index=False)
                sections.append(("POST_DOCKED_SCORCH - Unlabeled", lo_base[l_cols]))

            dbg("INFO", "reranked.summary", f"sections=Apo:{n_apo} Holo:{n_holo} other={len(leftover_rows)} out={summary_path}")
            if pretty_enabled and sections:
                pretty_path = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary(sections, pretty_path)
                print(f"[post_docked_scorch] pretty_summary_out={pretty_path}")
        else:
            _df = df.copy()
            cols = [c for c in preferred_cols if c in _df.columns] + [c for c in _df.columns if c not in preferred_cols]
            _df[cols].to_csv(fh, sep="\t", index=False)
            dbg("INFO", "reranked.summary", f"out={summary_path} rows={len(_df)}")
            if pretty_enabled:
                pretty_path = summary_path.with_name(summary_path.stem + "_pretty.txt")
                _write_pretty_summary([("POST_DOCKED_SCORCH", _df[cols])], pretty_path)
                print(f"[post_docked_scorch] pretty_summary_out={pretty_path}")

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
    # Primary: mirror the metrics row's run_id, if it exists.
    if meta is not None and meta.metrics is not None:
        if "run_id" in meta.metrics.index:
            return str(meta.metrics["run_id"])

    # Fallback: previous behaviour, but never join multiple IDs.
    if cli_run_id:
        if meta and meta.has_run_id_column:
            return str(cli_run_id)
        return "(none)"

    if meta and meta.has_run_id_column:
        if meta.run_ids:
            ordered = sorted(str(x) for x in meta.run_ids if x)
            if not ordered:
                return ""
            if active_run_id and active_run_id in ordered:
                return active_run_id
            return ordered[0]
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
        # Prefer the provided docked_root, fall back to legacy docked/ root last.
        _add(docked_root / pdb_id / "protein.log")
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


def _manifest_library_for_target(
    manifest: dict,
    pdb_id: str,
    variant_label: Optional[str],
    ph_tag: Optional[str],
) -> str:
    """
    Pull the library label for a given target directly from the run manifest.
    Falls back to a pdb_id-only search if the exact variant/pH key is missing.
    """
    if not manifest:
        return ""

    proteins = manifest.get("proteins") or {}
    if not isinstance(proteins, dict):
        return ""

    pdb_norm = str(pdb_id).upper()
    variant_norm = (variant_label or "legacy").strip().upper() or "LEGACY"
    ph_norm = (ph_tag or "base").strip()
    ph_token = ph_norm if ph_norm else "base"
    key = f"{pdb_norm}|{variant_norm}|{ph_token}"

    entry = proteins.get(key)
    if isinstance(entry, dict):
        lib = entry.get("library")
        if isinstance(lib, str) and lib.strip():
            return lib.strip()

    for _k, e in proteins.items():
        if not isinstance(e, dict):
            continue
        if str(e.get("pdb_id", "")).upper() != pdb_norm:
            continue
        lib = e.get("library")
        if isinstance(lib, str) and lib.strip():
            return lib.strip()

    return ""


def _manifest_protein_entries(manifest: dict) -> List[Dict[str, str]]:
    """
    Return the manifest proteins entries in order, without collapsing by pdb_id.
    Each entry is a dict with pdb_id, variant, and ph (strings, may be empty).
    """
    entries: List[Dict[str, str]] = []
    if not manifest:
        return entries
    proteins = manifest.get("proteins") or {}
    if not isinstance(proteins, dict):
        return entries
    for _key, entry in proteins.items():
        if not isinstance(entry, dict):
            continue
        pdb_id = str(entry.get("pdb_id", "")).strip()
        if not pdb_id:
            continue
        variant = str(entry.get("variant", "") or "").strip()
        ph = str(entry.get("ph", "") or "").strip()
        entries.append({"pdb_id": pdb_id, "variant": variant, "ph": ph})
    return entries


def _manifest_proteins_by_pdb(manifest: dict) -> Dict[str, Dict[str, str]]:
    """
    Return a mapping:
        pdb_id -> {"variant": <variant or \"\">, "ph": <ph or \"\">}
    using the manifest's proteins section.

    If multiple manifest entries share the same pdb_id, prefer the first.
    We treat the manifest's 'variant' and 'ph' fields as canonical labels
    and do NOT try to parse or split composite pH tags like 'pH7_2+7_7-dup4'.
    """
    out: Dict[str, Dict[str, str]] = {}
    if not manifest:
        return out
    proteins = manifest.get("proteins") or {}
    if not isinstance(proteins, dict):
        return out
    for _key, entry in proteins.items():
        if not isinstance(entry, dict):
            continue
        pdb_id = str(entry.get("pdb_id", "")).strip()
        if not pdb_id:
            continue
        if pdb_id in out:
            continue
        variant = str(entry.get("variant", "") or "").strip()
        ph = str(entry.get("ph", "") or "").strip()
        out[pdb_id] = {"variant": variant, "ph": ph}
    return out


def _build_control_records(
    pdb_id: str,
    run_id_label: str,
    target_name: str,
    library_name: str,
    log_path: Optional[Path],
    variant: Optional[str] = None,
    ph_tag: Optional[str] = None,
) -> Tuple[List[Dict], Dict, Dict[str, int]]:
    long_rows: List[Dict] = []
    controls_found = 0
    max_spread: Optional[float] = None
    policy = ""
    centers_detected = False

    if not log_path or not log_path.exists():
        summary_row = {
            "variant": variant,
            "pH": ph_tag,
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
                    "variant": variant,
                    "pH": ph_tag,
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
        "variant": variant,
        "pH": ph_tag,
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
    ap.add_argument("--out-dir", type=str, default="analysis/dud_eval",
                    help="Output root directory (run_id is appended automatically when provided).")
    ap.add_argument("--post-docked-root", type=str, default="post_docked",
                    help="Root folder for post-docked outputs (consensus reranked with SCORCH).")
    ap.add_argument("--lig-col", type=str, default=None,
                    help="Column containing ligand filename/path (auto-detected if omitted).")
    ap.add_argument("--score-col", type=str, default=None,
                    help="Score column (lower is better). Auto-detected if omitted.")
    ap.add_argument("--run-id", type=str, default=None,
                    help="Filter docking_score_long.csv rows to a specific run identifier.")
    ap.add_argument("--pdb-id", action="append", default=None,
                    help="Restrict evaluation to specific PDB IDs (repeatable or comma-separated).")
    ap.add_argument("-valid", "--valid", dest="valid_only", action="store_true", default=False,
                    help="Only evaluate poses marked valid in the docking_score_long.csv file.")
    ap.add_argument("--valid-col", type=str, default=None,
                    help="Validity column name to use with --valid (auto-detected if omitted).")
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
    pdb_id_filter = _normalize_pdb_ids(args.pdb_id)
    active_run_id: Optional[str] = args.run_id
    csv_basenames: Tuple[str, ...] = CSV_BASENAMES

    cfg = _load_default_cfg()
    docked_root = Path(args.docked_root)
    out_root = Path(args.out_dir)
    post_docked_root = Path(args.post_docked_root)
    pdb_root_override = Path(args.pdb_root) if args.pdb_root else None
    prepped_root_override = Path(args.prepped_root) if args.prepped_root else None
    log_root_override = Path(args.log_root) if getattr(args, "log_root", None) else None

    analysis_root = _compute_analysis_root(args, cfg)
    dbg("DEBUG", "paths", f"docked_root={docked_root} out_root={out_root} analysis_root={analysis_root} pdb_root={pdb_root_override or 'none'} prepped_root={prepped_root_override or 'none'}")
    scan_roots = _resolve_scan_roots(docked_root, args.run_id)
    dbg("INFO", "paths", f"scan_roots={[str(p) for p in scan_roots]}")
    manifest_docked_root_primary = scan_roots[0] if scan_roots else docked_root
    legacy_root = docked_root

    if args.run_id:
        os.environ["ATLAS_RUN_ID"] = str(args.run_id)
        if cfg is not None:
            cfg["RUN_ID"] = str(args.run_id)

    def _infer_variant_ph_from_csv_path(pdb_id: str, csv_path: Path) -> tuple[Optional[str], Optional[str]]:
        """
        Infer (variant, pH_tag) from the chosen docking_score_long.csv path.

        We treat 'docked/<PDB>/' from path_router as the root and look at
        the extra directory segments between that root and the CSV's parent.
        Any 'APO'/'HOLO' segment becomes the variant, any segment starting
        with 'pH'/'ph' becomes the pH tag.
        """
        try:
            paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
            base = paths.docked_pdb_root()  # docked/<PDB>/
        except Exception:
            base = docked_root / pdb_id

        try:
            rel = csv_path.parent.relative_to(base)
            parts = [p for p in rel.parts if p]
        except Exception:
            parts = []

        variant: Optional[str] = None
        ph_tag: Optional[str] = None
        for part in parts:
            up = part.upper()
            if variant is None and up in {"APO", "HOLO"}:
                variant = up
                continue
            low = part.lower()
            if ph_tag is None and low.startswith("ph"):
                ph_tag = part
        return variant, ph_tag

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

        def _record_basenames(base_path: Path) -> None:
            """
            Try all CSV_BASENAMES at base_path, e.g. base_path / name,
            in order. The first existing file wins via _record().
            """
            nonlocal picked
            for basename in csv_basenames:
                if picked is not None:
                    break
                candidate = base_path / basename
                _record(candidate)

        def _scan_tree(root: Path) -> None:
            """
            Look for CSV_BASENAMES starting at root.

            Search order:
            1) root/<basename> for each basename in CSV_BASENAMES.
            2) root/<variant>/<basename> (variant is any subdir name).
            3) root/<variant>/<pH>/<basename> (pH is any sub-subdir name).
            We stop as soon as _record() finds a hit.
            """
            nonlocal picked

            # 1) root/<basename>
            if picked is None:
                _record_basenames(root)

            # 2) root/<variant>/<basename> and 3) root/<variant>/<pH>/<basename>
            if picked is not None:
                return

            if not root.is_dir():
                return

            for variant_dir in sorted(root.iterdir()):
                if picked is not None:
                    break
                if not variant_dir.is_dir():
                    continue

                # 2) docked/<PDB>/<VARIANT>/<basename>
                _record_basenames(variant_dir)

                if picked is not None:
                    break

                # 3) docked/<PDB>/<VARIANT>/<pH>/<basename>
                for ph_dir in sorted(variant_dir.iterdir()):
                    if picked is not None:
                        break
                    if not ph_dir.is_dir():
                        continue
                    _record_basenames(ph_dir)

        if not cfg:
            _record_basenames(fallback_dir)
            if picked is None:
                _scan_tree(fallback_dir)
        else:
            try:
                mode    = str(cfg.get("APO_HOLO_MODE", "")).strip()
                variants = expand_variants(mode)  # returns [None] | ["APO","HOLO"]
                paths   = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
                ph_token = cfg.get("PH_TOKEN") or None
            except Exception as exc:
                dbg("WARN", "resolve", f"pdb={pdb_id} router_err={exc}")
                _record_basenames(fallback_dir)
                if picked is None:
                    _scan_tree(fallback_dir)
            else:
                docked_root_cfg = paths.docked_pdb_root()
                summary_csv = paths.docking_score_summary_csv()
                long_csv    = paths.docking_score_long_csv()

                if docked_root_cfg.exists():
                    # 1) Try path_router's long_csv location, but with both basenames
                    for basename in csv_basenames:
                        if picked is not None:
                            break
                        candidate = long_csv.with_name(basename)
                        _record(candidate)

                    # 2) If still nothing, try variant roots (and their pH subdirs)
                    if picked is None:
                        for variant in variants:
                            variant_root = paths.docked_variant_root(variant)
                            _record_basenames(variant_root)
                            if picked is not None:
                                break
                            if variant_root.is_dir():
                                for ph_dir in sorted(variant_root.iterdir()):
                                    if picked is not None:
                                        break
                                    if ph_dir.is_dir():
                                        _record_basenames(ph_dir)

                    # 3) Final fallback: CLI fallback_dir (includes nested scan)
                    if picked is None:
                        _scan_tree(fallback_dir)

                    # Touch summary path to exercise router (no fallback to summary file for eval)
                    try:
                        summary_csv.exists()
                    except Exception:
                        pass
                else:
                    dbg("WARN", "resolve", f"pdb={pdb_id} docked_root_missing root={docked_root_cfg}")
                    _record_basenames(fallback_dir)

        if picked is not None:
            dbg("DEBUG", "resolve", f"pdb={pdb_id} tried={len(candidates_tried)} candidates={candidates_tried}")
            dbg("INFO", "resolve", f"pdb={pdb_id} picked={picked}")

            return picked

        dbg("WARN", "resolve", f"pdb={pdb_id} no_csv_found tried={candidates_tried or ['<none>']} search_root={fallback_dir}")
        return None

    def _resolve_docking_csv_for_manifest(
        pdb_id: str,
        variant_hint: Optional[str],
        ph_hint: Optional[str],
        docked_root: Path,
        cfg: Dict,
        fallback_root: Optional[Path] = None,
    ) -> Optional[Path]:
        """
        Variant/pH-aware wrapper around _resolve_docking_csv.

        Strategy:
          1) If cfg and variant_hint/ph_hint are provided, try the exact
             docked/<PDB>/<VARIANT>/<PH>/<CSV_BASENAME> locations first.
          2) If that fails, fall back to _resolve_docking_csv(pdb_id, docked_root / pdb_id).
          3) If still missing and fallback_root provided, try legacy fallback_root/pdb_id.
        """
        try:
            paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
            base_root = paths.docked_pdb_root()
        except Exception:
            paths = None
            base_root = docked_root / pdb_id

        if variant_hint or ph_hint:
            var = (variant_hint or "").strip()
            ph = (ph_hint or "").strip()
            candidate_roots: List[Path] = []
            if var:
                try:
                    variant_root = paths.docked_variant_root(var) if paths else base_root / var
                except Exception:
                    variant_root = base_root / var
                if ph:
                    candidate_roots.append(variant_root / ph)
                else:
                    candidate_roots.append(variant_root)
            elif ph:
                candidate_roots.append(base_root / ph)

            for root in candidate_roots:
                for basename in csv_basenames:
                    candidate = root / basename
                    if candidate.exists():
                        dbg("INFO", "resolve", f"pdb={pdb_id} picked={candidate} via=manifest_hint")
                        return candidate

        picked = _resolve_docking_csv(pdb_id, docked_root / pdb_id)
        if picked is None and fallback_root and fallback_root != docked_root:
            picked = _resolve_docking_csv(pdb_id, fallback_root / pdb_id)
        return picked

    def _discover_docking_csvs(
        pdb_id: str,
        docked_root: Path,
        cfg: Optional[dict],
        csv_basenames: Tuple[str, ...],
    ) -> List[TargetSpec]:
        """
        Return all docking CSVs for a pdb_id across variants/pH subfolders.
        """
        try:
            paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb") if cfg else None
        except Exception:
            paths = None

        if paths:
            base_root = paths.docked_pdb_root()
        else:
            base_root = docked_root / pdb_id
            if docked_root.name.upper() == pdb_id.upper():
                base_root = docked_root

        seen_paths: Set[str] = set()
        ranked: List[Tuple[TargetSpec, int]] = []
        basename_rank = {name: idx for idx, name in enumerate(csv_basenames)}

        def _add(candidate: Path, variant: Optional[str], ph_tag: Optional[str], source: str) -> None:
            if not candidate.exists():
                return
            try:
                key = str(candidate.resolve())
            except Exception:
                key = str(candidate)
            if key in seen_paths:
                return
            seen_paths.add(key)
            variant_norm = variant.upper() if variant else None
            rank = basename_rank.get(candidate.name, len(csv_basenames))
            ranked.append((
                TargetSpec(
                    target_key=make_target_key(pdb_id, variant_norm, ph_tag),
                    pdb_id=pdb_id,
                    variant=variant_norm,
                    ph_tag=ph_tag,
                    csv_path=candidate,
                    source=source,
                ),
                rank,
            ))

        # Legacy root: docked/<PDB>/<basename>
        for basename in csv_basenames:
            _add(base_root / basename, None, None, "legacy")

        if base_root.is_dir():
            for ph_dir in sorted(base_root.iterdir()):
                if not ph_dir.is_dir():
                    continue
                ph_label = ph_dir.name
                for basename in csv_basenames:
                    _add(ph_dir / basename, None, ph_label, "scan")

        for variant in ("APO", "HOLO"):
            variant_dir = base_root / variant
            for basename in csv_basenames:
                _add(variant_dir / basename, variant, None, "scan")
            if variant_dir.is_dir():
                for ph_dir in sorted(variant_dir.iterdir()):
                    if not ph_dir.is_dir():
                        continue
                    ph_label = ph_dir.name
                    for basename in csv_basenames:
                        _add(ph_dir / basename, variant, ph_label, "scan")

        ranked.sort(key=lambda pair: (
            pair[0].pdb_id,
            pair[0].variant or "",
            pair[0].ph_tag or "",
            pair[1],
            str(pair[0].csv_path),
        ))
        specs = [spec for spec, _rank in ranked]
        pairs = [f"{spec.variant or 'legacy'}|{spec.ph_tag or 'base'}" for spec in specs]
        variant_label = ";".join(pairs) if pairs else "<none>"
        dbg("INFO", "discover", f"pdb={pdb_id} targets={len(specs)} variants={variant_label}")
        return specs

    manifest_data: dict | None = None
    manifest_entries: List[Dict[str, str]] = []
    manifest_proteins: Dict[str, Dict[str, str]] = {}
    manifest_driven = False
    effective_manifest_root = docked_root
    if cfg and active_run_id:
        if len(scan_roots) == 1 and scan_roots[0].is_dir() and _is_probable_run_id_dirname(scan_roots[0].name):
            effective_manifest_root = scan_roots[0]
        try:
            _manifest_dir, manifest_path = get_manifest_paths(cfg, str(active_run_id))
            manifest_data = _load_manifest(manifest_path) or {}
            manifest_entries = _manifest_protein_entries(manifest_data or {})
            manifest_proteins = _manifest_proteins_by_pdb(manifest_data or {})
            manifest_driven = bool(active_run_id and manifest_entries)
            test_mode_enable = ""
            cmd_block = manifest_data.get("command") if isinstance(manifest_data, dict) else {}
            if isinstance(cmd_block, dict):
                test_mode_enable = str(cmd_block.get("TEST_MODE_ENABLE", "")).strip().lower()
            if test_mode_enable == "dud":
                csv_basenames = ("docking_score_long.csv", "dud_docking_score_long.csv")
                dbg("INFO", "manifest", f"run_id={active_run_id} test_mode=dud csv_priority={';'.join(csv_basenames)}")
            dbg("INFO", "manifest", f"run_id={active_run_id} manifest={manifest_path} proteins={len(manifest_entries)}")
        except Exception as exc:
            dbg("WARN", "manifest", f"run_id={active_run_id} load_err={exc}")
            manifest_data = {}
            manifest_entries = []
            manifest_proteins = {}
            manifest_driven = False
    else:
        dbg("DEBUG", "manifest", "pre-discover: no cfg or active_run_id; manifest lookup disabled")

    targets: List[TargetSpec] = []

    if manifest_driven:
        for entry in manifest_entries:
            pdb_id = str(entry.get("pdb_id", "")).strip()
            if not pdb_id:
                continue
            variant_hint = str(entry.get("variant", "") or "").strip()
            ph_hint = str(entry.get("ph", "") or "").strip()
            variant_upper = variant_hint.upper() if variant_hint else None
            csv_path = _resolve_docking_csv_for_manifest(
                pdb_id=pdb_id,
                variant_hint=variant_hint or None,
                ph_hint=ph_hint or None,
                docked_root=manifest_docked_root_primary,
                cfg=cfg,
                fallback_root=legacy_root,
            )
            target_key = make_target_key(pdb_id, variant_upper, ph_hint or None)
            targets.append(TargetSpec(
                target_key=target_key,
                pdb_id=pdb_id,
                variant=variant_upper,
                ph_tag=ph_hint or None,
                csv_path=csv_path,
                source="manifest",
            ))
            if csv_path is None:
                dbg("WARN", "discover", f"pdb={pdb_id} variant={variant_upper or 'base'} ph={ph_hint or 'base'} reason=no_docked_csv_for_manifest_entry")
    elif docked_root.is_dir():
        # Iterate over candidate scan roots (run-scoped and legacy).
        for root in scan_roots:
            if not root.is_dir():
                # Treat as single-target fallback
                pdb_id = root.name
                targets.extend(_discover_docking_csvs(pdb_id, root, cfg, csv_basenames))
                continue

            any_added = False
            for sub in sorted(root.iterdir()):
                if not sub.is_dir():
                    continue
                pdb_id = sub.name

                pdb_path = Path("input_pdbs") / f"{pdb_id}.pdb"
                if not pdb_path.exists():
                    dbg("WARN", "discover.skip", f"pdb={pdb_id} reason=no_input_pdb path={pdb_path}")
                    continue

                targets.extend(_discover_docking_csvs(pdb_id, root, cfg, csv_basenames))
                any_added = True

            if not any_added:
                pdb_id = root.name
                targets.extend(_discover_docking_csvs(pdb_id, root, cfg, csv_basenames))
    else:
        pdb_id = docked_root.name
        targets.extend(_discover_docking_csvs(pdb_id, docked_root, cfg, csv_basenames))

    unique_targets: List[TargetSpec] = []
    seen_keys: Set[str] = set()
    for spec in targets:
        if spec.target_key in seen_keys:
            dbg("WARN", "discover", f"target_key={spec.target_key} csv={spec.csv_path} source={spec.source} action=skip_duplicate")
            continue
        seen_keys.add(spec.target_key)
        unique_targets.append(spec)
    targets = unique_targets

    if pdb_id_filter:
        filter_set = set(pdb_id_filter)
        before = len(targets)
        targets = [t for t in targets if t.pdb_id.upper() in filter_set]
        dbg("INFO", "dud-eval", f"filter pdb_id={sorted(filter_set)} kept={len(targets)}/{before}")

    if not targets:
        if pdb_id_filter:
            dbg("ERROR", "discover", f"no targets after pdb_id filter pdb_id={pdb_id_filter}")
        else:
            dbg("ERROR", "discover", f"no docking_score_long.csv / dud_docking_score_long.csv under {docked_root}")
        raise SystemExit(2)

    dbg("INFO", "discover", f"targets={len(targets)} root={docked_root}")

    rows: List[pd.Series] = []
    consensus_rows: List[pd.Series] = []
    reranked_rows: List[pd.Series] = []
    target_eval_results: Dict[str, TargetEvaluation] = {}
    csv_paths_by_target: Dict[str, Path] = {}
    for spec in targets:
        if spec.csv_path:
            csv_paths_by_target[spec.target_key] = spec.csv_path
    targets_with_csv = [t for t in targets if t.csv_path]

    if not active_run_id:
        try:
            active_run_id = select_default_run_id(targets_with_csv, csv_paths_by_target, args.lig_col, args.score_col)
        except Exception as _exc:
            active_run_id = None
            dbg("WARN", "run", f"auto_select_failed err={_exc}")
    dbg("INFO", "run", f"active={active_run_id or '(none)'} source={'CLI' if args.run_id else 'auto'}")
    run_dir = Path(args.docked_root) / str(active_run_id) if active_run_id else None
    if active_run_id:
        analysis_root = _compute_analysis_root(args, cfg, run_id_override=active_run_id)
        dbg("INFO", "paths", f"analysis_root_updated run_id={active_run_id} path={analysis_root}")
    new_layout_root_exists = bool(run_dir and run_dir.exists())
    dbg("INFO", "run.layout", f"run_id={active_run_id or '(none)'} new_layout_root_exists={new_layout_root_exists}")
    run_label_for_files = _format_run_label(active_run_id)

    if manifest_data is None and cfg and active_run_id:
        try:
            _manifest_dir, manifest_path = get_manifest_paths(cfg, str(active_run_id))
            manifest_data = _load_manifest(manifest_path) or {}
            manifest_entries = _manifest_protein_entries(manifest_data or {})
            manifest_proteins = _manifest_proteins_by_pdb(manifest_data or {})
            manifest_driven = manifest_driven or bool(active_run_id and manifest_entries)
            test_mode_enable = ""
            cmd_block = manifest_data.get("command") if isinstance(manifest_data, dict) else {}
            if isinstance(cmd_block, dict):
                test_mode_enable = str(cmd_block.get("TEST_MODE_ENABLE", "")).strip().lower()
            if test_mode_enable == "dud":
                csv_basenames = ("docking_score_long.csv", "dud_docking_score_long.csv")
                dbg("INFO", "manifest", f"run_id={active_run_id} test_mode=dud csv_priority={';'.join(csv_basenames)}")
            dbg("INFO", "manifest", f"run_id={active_run_id} manifest={manifest_path} proteins={len(manifest_entries)}")
        except Exception as exc:
            dbg("WARN", "manifest", f"run_id={active_run_id} load_err={exc}")
            manifest_data = {}
            manifest_entries = []
            manifest_proteins = {}
            manifest_driven = bool(active_run_id and manifest_proteins)
    elif manifest_data is None:
        dbg("DEBUG", "manifest", "no cfg or active_run_id; manifest lookup disabled")

    if not manifest_entries:
        manifest_entries = _manifest_protein_entries(manifest_data or {})
    if not manifest_proteins:
        manifest_proteins = _manifest_proteins_by_pdb(manifest_data or {})
    manifest_driven = manifest_driven or bool(active_run_id and manifest_entries)

    for spec in targets:
        evaluated: Optional[TargetEvaluation] = None
        status_reason = "missing_csv"
        out_dir = analysis_root / spec.target_key
        if spec.csv_path:
            evaluated = evaluate_target(
                pdb_id=spec.pdb_id,
                csv_path=spec.csv_path,
                out_dir=out_dir,
                lig_col_cli=args.lig_col,
                score_col_cli=args.score_col,
                bedroc_alpha=args.bedroc_alpha,
                logauc_lambda=args.logauc_lambda,
                run_id=active_run_id,
                valid_only=args.valid_only,
                valid_col_cli=args.valid_col,
                new_layout_active=bool(run_dir and spec.csv_path and _is_under(spec.csv_path, run_dir)),
                layout_label="new_run" if (run_dir and spec.csv_path and _is_under(spec.csv_path, run_dir)) else "legacy",
            )
            if evaluated is not None:
                status_reason = evaluated.status_reason or "ok"
            else:
                status_reason = "read_error"
        consensus_candidate: Optional[Path] = None
        if spec.csv_path:
            candidate = spec.csv_path.parent / CONSENSUS_CSV_BASENAME
            if candidate.exists():
                consensus_candidate = candidate
        if consensus_candidate is not None:
            cons_eval = evaluate_target_consensus(
                pdb_id=spec.pdb_id,
                csv_path=consensus_candidate,
                out_dir=analysis_root / "consensus" / spec.target_key,
                lig_col_cli=args.lig_col,
                score_col_cli=args.score_col,
                bedroc_alpha=args.bedroc_alpha,
                logauc_lambda=args.logauc_lambda,
                run_id=active_run_id,
                variant_hint=spec.variant,
                ph_tag=spec.ph_tag,
                new_layout_active=bool(run_dir and consensus_candidate and _is_under(consensus_candidate, run_dir)),
                layout_label="new_run" if (run_dir and consensus_candidate and _is_under(consensus_candidate, run_dir)) else "legacy",
            )
            if cons_eval is not None and cons_eval.metrics is not None:
                cons_metrics = cons_eval.metrics.copy()
                for col in ("variant", "pH"):
                    if col in cons_metrics.index:
                        cons_metrics.drop(index=col, inplace=True)
                cons_metrics["pdb_id"] = spec.pdb_id
                cons_metrics["variant"] = spec.variant or ""
                cons_metrics["pH"] = spec.ph_tag or ""
                if "status_reason" not in cons_metrics:
                    cons_metrics["status_reason"] = cons_eval.status_reason or "ok"
                consensus_rows.append(cons_metrics)
            elif cons_eval is not None:
                placeholder_cons = _make_placeholder_row(
                    spec.pdb_id,
                    spec.variant,
                    spec.ph_tag,
                    active_run_id,
                    args.bedroc_alpha,
                    cons_eval.status_reason or "empty_after_filter",
                )
                consensus_rows.append(placeholder_cons)
            else:
                dbg("WARN", "consensus.eval", f"pdb={spec.pdb_id} path={consensus_candidate} reason=evaluation_failed")
        reranked_candidate: Optional[Path] = _resolve_reranked_scorch_path(
            spec, docked_root, post_docked_root, active_run_id, run_dir
        )
        if reranked_candidate is not None:
            post_run_root = post_docked_root / str(active_run_id) if active_run_id else None
            new_layout_reranked = bool(post_run_root and _is_under(reranked_candidate, post_run_root))
            reranked_eval = evaluate_target_post_docked_reranked_scorch(
                pdb_id=spec.pdb_id,
                csv_path=reranked_candidate,
                out_dir=analysis_root / "post_docked" / "consensus_reranked_scorch" / spec.target_key,
                lig_col_cli=args.lig_col,
                score_col_cli=args.score_col,
                bedroc_alpha=args.bedroc_alpha,
                logauc_lambda=args.logauc_lambda,
                run_id=active_run_id,
                variant_hint=spec.variant,
                ph_tag=spec.ph_tag,
                new_layout_active=new_layout_reranked,
                layout_label="post_docked_new" if new_layout_reranked else "post_docked_legacy",
            )
            if reranked_eval is not None and reranked_eval.metrics is not None:
                rer_metrics = reranked_eval.metrics.copy()
                for col in ("variant", "pH"):
                    if col in rer_metrics.index:
                        rer_metrics.drop(index=col, inplace=True)
                rer_metrics["pdb_id"] = spec.pdb_id
                rer_metrics["variant"] = spec.variant or ""
                rer_metrics["pH"] = spec.ph_tag or ""
                rer_metrics["run_id"] = active_run_id or ""
                if "status_reason" not in rer_metrics:
                    rer_metrics["status_reason"] = reranked_eval.status_reason or "ok"
                reranked_rows.append(rer_metrics)
            elif reranked_eval is not None:
                placeholder_reranked = _make_placeholder_row(
                    spec.pdb_id,
                    spec.variant,
                    spec.ph_tag,
                    active_run_id,
                    args.bedroc_alpha,
                    reranked_eval.status_reason or "empty_after_filter",
                )
                reranked_rows.append(placeholder_reranked)
        if evaluated is not None:
            target_eval_results[spec.target_key] = evaluated
            if evaluated.metrics is not None:
                metrics = evaluated.metrics.copy()
                for col in ("variant", "pH"):
                    if col in metrics.index:
                        metrics.drop(index=col, inplace=True)
                metrics["pdb_id"] = spec.pdb_id
                metrics["variant"] = spec.variant or ""
                metrics["pH"] = spec.ph_tag or ""
                if "status_reason" not in metrics:
                    metrics["status_reason"] = status_reason
                rows.append(metrics)
            else:
                placeholder = _make_placeholder_row(
                    spec.pdb_id,
                    spec.variant,
                    spec.ph_tag,
                    active_run_id,
                    args.bedroc_alpha,
                    status_reason,
                )
                rows.append(placeholder)
        else:
            placeholder = _make_placeholder_row(
                spec.pdb_id,
                spec.variant,
                spec.ph_tag,
                active_run_id,
                args.bedroc_alpha,
                status_reason,
            )
            rows.append(placeholder)

    if not rows:
        dbg("ERROR", "metrics", "no targets produced evaluable rows")
        raise SystemExit(3)

    df_all = pd.DataFrame(rows)
    df_consensus_all = pd.DataFrame(consensus_rows) if consensus_rows else pd.DataFrame()
    df_reranked_all = pd.DataFrame(reranked_rows) if reranked_rows else pd.DataFrame()
    for col in ("variant", "pH"):
        if col not in df_all.columns:
            df_all[col] = ""
    has_variant = "variant" in df_all.columns and df_all["variant"].astype(str).str.strip().ne("").any()
    has_ph = "pH" in df_all.columns and df_all["pH"].astype(str).str.strip().ne("").any()
    if not has_variant and "variant" in df_all.columns:
        df_all = df_all.drop(columns=["variant"])
    if not has_ph and "pH" in df_all.columns:
        df_all = df_all.drop(columns=["pH"])
    sort_cols = [c for c in ("pdb_id", "variant", "pH") if c in df_all.columns]
    df_all = df_all.sort_values(sort_cols or ["pdb_id"])

    if not df_consensus_all.empty:
        for col in ("variant", "pH"):
            if col not in df_consensus_all.columns:
                df_consensus_all[col] = ""
        cons_has_variant = "variant" in df_consensus_all.columns and df_consensus_all["variant"].astype(str).str.strip().ne("").any()
        cons_has_ph = "pH" in df_consensus_all.columns and df_consensus_all["pH"].astype(str).str.strip().ne("").any()
        if not cons_has_variant and "variant" in df_consensus_all.columns:
            df_consensus_all = df_consensus_all.drop(columns=["variant"])
        if not cons_has_ph and "pH" in df_consensus_all.columns:
            df_consensus_all = df_consensus_all.drop(columns=["pH"])
        sort_cols_cons = [c for c in ("pdb_id", "variant", "pH") if c in df_consensus_all.columns]
        df_consensus_all = df_consensus_all.sort_values(sort_cols_cons or ["pdb_id"])

    if manifest_driven:
        expected_keys = [
            make_target_key(
                str(entry.get("pdb_id", "")).strip(),
                (entry.get("variant") or "").strip() or None,
                (entry.get("ph") or "").strip() or None,
            )
            for entry in manifest_entries
            if str(entry.get("pdb_id", "")).strip()
        ]
        observed_keys: List[str] = []
        for _, row in df_all.iterrows():
            variant_val = ""
            ph_val = ""
            if "variant" in df_all.columns:
                raw_var = row.get("variant", "")
                variant_val = "" if pd.isna(raw_var) else str(raw_var)
            if "pH" in df_all.columns:
                raw_ph = row.get("pH", "")
                ph_val = "" if pd.isna(raw_ph) else str(raw_ph)
            observed_keys.append(make_target_key(str(row.get("pdb_id", "")), variant_val, ph_val))
        missing = sorted(set(expected_keys) - set(observed_keys))
        extra = sorted(set(observed_keys) - set(expected_keys))
        dbg("INFO", "manifest.coverage", f"expected={len(expected_keys)} observed={len(observed_keys)} missing={len(missing)} extra={len(extra)}")
        if missing:
            suffix = "..." if len(missing) > 5 else ""
            dbg("WARN", "manifest.coverage", f"missing_keys={missing[:5]}{suffix}")

    annotate_targets = list(df_all["pdb_id"].astype(str).unique())
    if not df_consensus_all.empty:
        for pdb_id in df_consensus_all["pdb_id"].astype(str).unique():
            if pdb_id not in annotate_targets:
                annotate_targets.append(pdb_id)

    names: Dict[str, str] = {p: "" for p in annotate_targets}
    if args.target_name_from_pdb:
        for pdb_id in annotate_targets:
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

    if args.report_library:
        if manifest_data:
            df_all["library_name"] = df_all.apply(
                lambda r: _manifest_library_for_target(
                    manifest=manifest_data,
                    pdb_id=str(r.get("pdb_id", "")),
                    variant_label=(
                        "" if "variant" not in df_all.columns else ("" if pd.isna(r.get("variant", "")) else str(r.get("variant", "")))
                    ),
                    ph_tag=(
                        "" if "pH" not in df_all.columns else ("" if pd.isna(r.get("pH", "")) else str(r.get("pH", "")))
                    ),
                ),
                axis=1,
            )
            df_all["library_name"] = df_all["library_name"].fillna("")
        else:
            df_all["library_name"] = ""
    else:
        dbg("DEBUG", "library", "report_library=OFF")
        if "library_name" not in df_all.columns:
            df_all["library_name"] = ""

    if not df_all.empty and "target_name" not in df_all.columns:
        df_all["target_name"] = df_all["pdb_id"].map(names).fillna("")
    if "library_name" not in df_all.columns:
        df_all["library_name"] = ""

    if not df_consensus_all.empty:
        if args.report_library:
            if manifest_data:
                df_consensus_all["library_name"] = df_consensus_all.apply(
                    lambda r: _manifest_library_for_target(
                        manifest=manifest_data,
                        pdb_id=str(r.get("pdb_id", "")),
                        variant_label=(
                            "" if "variant" not in df_consensus_all.columns else ("" if pd.isna(r.get("variant", "")) else str(r.get("variant", "")))
                        ),
                        ph_tag=(
                            "" if "pH" not in df_consensus_all.columns else ("" if pd.isna(r.get("pH", "")) else str(r.get("pH", "")))
                        ),
                    ),
                    axis=1,
                )
                df_consensus_all["library_name"] = df_consensus_all["library_name"].fillna("")
            else:
                df_consensus_all["library_name"] = df_consensus_all.get("library_name", pd.Series("", index=df_consensus_all.index))
        else:
            if "library_name" not in df_consensus_all.columns:
                df_consensus_all["library_name"] = ""
        if "target_name" not in df_consensus_all.columns:
            df_consensus_all["target_name"] = df_consensus_all["pdb_id"].map(names).fillna("")

    if not df_reranked_all.empty:
        for col in ("variant", "pH"):
            if col not in df_reranked_all.columns:
                df_reranked_all[col] = ""
        reranked_has_variant = (
            "variant" in df_reranked_all.columns and df_reranked_all["variant"].astype(str).str.strip().ne("").any()
        )
        reranked_has_ph = "pH" in df_reranked_all.columns and df_reranked_all["pH"].astype(str).str.strip().ne("").any()
        if not reranked_has_variant and "variant" in df_reranked_all.columns:
            df_reranked_all = df_reranked_all.drop(columns=["variant"])
        if not reranked_has_ph and "pH" in df_reranked_all.columns:
            df_reranked_all = df_reranked_all.drop(columns=["pH"])
        sort_cols_reranked = [c for c in ("pdb_id", "variant", "pH") if c in df_reranked_all.columns]
        df_reranked_all = df_reranked_all.sort_values(sort_cols_reranked or ["pdb_id"])

        if args.report_library:
            if manifest_data:
                df_reranked_all["library_name"] = df_reranked_all.apply(
                    lambda r: _manifest_library_for_target(
                        manifest=manifest_data,
                        pdb_id=str(r.get("pdb_id", "")),
                        variant_label=(
                            "" if "variant" not in df_reranked_all.columns else ("" if pd.isna(r.get("variant", "")) else str(r.get("variant", "")))
                        ),
                        ph_tag=(
                            "" if "pH" not in df_reranked_all.columns else ("" if pd.isna(r.get("pH", "")) else str(r.get("pH", "")))
                        ),
                    ),
                    axis=1,
                )
                df_reranked_all["library_name"] = df_reranked_all["library_name"].fillna("")
            else:
                df_reranked_all["library_name"] = df_reranked_all.get("library_name", pd.Series("", index=df_reranked_all.index))
        else:
            if "library_name" not in df_reranked_all.columns:
                df_reranked_all["library_name"] = ""
        if "target_name" not in df_reranked_all.columns:
            df_reranked_all["target_name"] = df_reranked_all["pdb_id"].map(names).fillna("")

    control_targets = [t for t in targets if t.target_key in target_eval_results]

    df = df_all.copy()
    df_consensus = df_consensus_all.copy()

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
    if not df_consensus.empty:
        cons_lib_norm = df_consensus["library_name"].astype(str).str.strip()
        cons_pdb_norm = df_consensus["pdb_id"].astype(str).str.strip()
        cons_mask_fda = cons_lib_norm.str.lower().isin({"fda_library", "fda"})
        cons_mask_pdb = cons_lib_norm.str.upper() == cons_pdb_norm.str.upper()
        cons_exclude_mask = cons_mask_fda | cons_mask_pdb
        dropped_consensus = int(cons_exclude_mask.sum())
        if dropped_consensus:
            dbg("INFO", "consensus.exclude", f"excluded={dropped_consensus} reason=library_filter")
        df_consensus = df_consensus.loc[~cons_exclude_mask].copy()
        sort_cols_cons = [c for c in ("pdb_id", "variant", "pH") if c in df_consensus.columns]
        df_consensus = df_consensus.sort_values(sort_cols_cons or ["pdb_id"])
    if not df_reranked_all.empty:
        rer_lib_norm = df_reranked_all["library_name"].astype(str).str.strip()
        rer_pdb_norm = df_reranked_all["pdb_id"].astype(str).str.strip()
        rer_mask_fda = rer_lib_norm.str.lower().isin({"fda_library", "fda"})
        rer_mask_pdb = rer_lib_norm.str.upper() == rer_pdb_norm.str.upper()
        rer_exclude_mask = rer_mask_fda | rer_mask_pdb
        if rer_exclude_mask.any():
            dbg("INFO", "reranked.exclude", f"excluded={int(rer_exclude_mask.sum())} reason=library_filter")
        df_reranked_all = df_reranked_all.loc[~rer_exclude_mask].copy()
        sort_cols_rer = [c for c in ("pdb_id", "variant", "pH") if c in df_reranked_all.columns]
        df_reranked_all = df_reranked_all.sort_values(sort_cols_rer or ["pdb_id"])

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
        for spec in control_targets:
            eval_meta = target_eval_results.get(spec.target_key)
            run_label = _resolve_run_label(spec.pdb_id, eval_meta, args.run_id, active_run_id)
            target_label = names.get(spec.pdb_id, "")
            library_label = ""
            if args.report_library and manifest_data:
                library_label = _manifest_library_for_target(
                    manifest=manifest_data,
                    pdb_id=spec.pdb_id,
                    variant_label=spec.variant,
                    ph_tag=spec.ph_tag,
                )
            if not library_label and not df.empty and "library_name" in df.columns:
                mask = df["pdb_id"].astype(str) == str(spec.pdb_id)
                if "variant" in df.columns:
                    mask &= df["variant"].astype(str) == (spec.variant or "")
                if "pH" in df.columns:
                    mask &= df["pH"].astype(str) == (spec.ph_tag or "")
                if mask.any():
                    library_label = str(df.loc[mask, "library_name"].iloc[0])

            csv_path = csv_paths_by_target.get(spec.target_key)
            if csv_path is None:
                continue
            parsed_counts[spec.target_key] = {"centers": 0, "redock_lines": 0, "report_rows": 0}
            selected_log, candidates = _candidate_protein_logs(
                spec.pdb_id,
                csv_path,
                docked_root,
                log_root_override,
                cfg,
            )
            dbg(
                "DEBUG",
                "control",
                f"target={spec.target_key} log_candidates={[str(c) for c in candidates]}",
            )
            if selected_log:
                dbg("INFO", "control", f"target={spec.target_key} protein_log={selected_log}")
                scanned_ok.append(spec.target_key)
            else:
                missing_logs.append(spec.target_key)
                missing_hint = str(candidates[0]) if candidates else str(Path("docked") / spec.pdb_id / "protein.log")
                dbg(
                    "WARN",
                    "control",
                    f"target={spec.target_key} protein_log_missing={missing_hint} action=skip_control_parse",
                )

            long_rows, summary_row, meta = _build_control_records(
                spec.pdb_id,
                run_label,
                target_label,
                library_label,
                selected_log,
                spec.variant,
                spec.ph_tag,
            )
            if spec.target_key not in parsed_counts:
                parsed_counts[spec.target_key] = {"centers": 0, "redock_lines": 0, "report_rows": 0}
            parsed_counts[spec.target_key]["centers"] = meta.get("centers_found", 0)
            parsed_counts[spec.target_key]["redock_lines"] = meta.get("redock_lines", 0)
            parsed_counts[spec.target_key]["report_rows"] = meta.get("report_rows", 0)
            dbg(
                "INFO",
                "control",
                (
                    f"target={spec.target_key} centers_found={parsed_counts[spec.target_key]['centers']} "
                    f"redock_lines={parsed_counts[spec.target_key]['redock_lines']} "
                    f"report_rows={parsed_counts[spec.target_key]['report_rows']}"
                ),
            )
            if long_rows:
                control_long_records.extend(long_rows)
            control_summary_records.append(summary_row)
        for spec in control_targets:
            if spec.target_key not in parsed_counts:
                parsed_counts[spec.target_key] = {"centers": 0, "redock_lines": 0, "report_rows": 0}
                if spec.target_key not in missing_logs and spec.target_key not in scanned_ok:
                    missing_logs.append(spec.target_key)

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

    preferred_cols = ("variant", "pH", "run_id", "target_name", "library_name", "pdb_id")
    excluded_cols = {"run_id", "target_name", "library_name", "pdb_id", "N", "n_actives", "actives_fraction", "variant", "pH", "status_reason"}

    with open(summary_path, "w", newline="") as fh:
        if has_sections:
            _df = df.copy()
            cols = [c for c in preferred_cols if c in _df.columns] \
                   + [c for c in _df.columns if c not in preferred_cols]
            apo_rows = df.loc[apo_mask].sort_values("pdb_id") if apo_mask is not None else pd.DataFrame()
            holo_rows = df.loc[holo_mask].sort_values("pdb_id") if holo_mask is not None else pd.DataFrame()
            n_apo = len(apo_rows)
            n_holo = len(holo_rows)

            sections = []
            leftover_rows = pd.DataFrame()
            if n_apo:
                apo_out = apo_rows
                a_cols = [c for c in preferred_cols if c in apo_out.columns] \
                         + [c for c in apo_out.columns if c not in preferred_cols]
                fh.write("Apo\n");
                apo_out[a_cols].to_csv(fh, sep="\t", index=False, float_format="%.3f")
                sections.append(("Apo", apo_out[a_cols]))

            if n_holo:
                if n_apo: fh.write("\n")
                holo_out = holo_rows
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
                lo_base = leftover_rows
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
            _df = df.copy()
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
    emit_consensus_summary(df_consensus, analysis_root, run_label_for_files, getattr(args, "pretty_summary", True))
    emit_reranked_scorch_summary(df_reranked_all, analysis_root, run_label_for_files, getattr(args, "pretty_summary", True))

    if getattr(args, "emit_control_report", False):
        long_columns = [
            "variant",
            "pH",
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
        for col in ("variant", "pH"):
            if col in control_summary_df.columns:
                s = control_summary_df[col]
                if s.isna().all() or (s.astype(str).str.strip() == "").all():
                    control_summary_df.drop(columns=[col], inplace=True)
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
    for spec in control_targets:
        stats = parsed_counts.get(spec.target_key, {"centers": 0, "redock_lines": 0, "report_rows": 0})
        per_target_rows.append(
            {
                "target_key": spec.target_key,
                "pdb_id": spec.pdb_id,
                "variant": spec.variant or "",
                "pH": spec.ph_tag or "",
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
                "target_key",
                "pdb_id",
                "variant",
                "pH",
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
