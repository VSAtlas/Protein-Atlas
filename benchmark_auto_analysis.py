# --- snip header / your docstring stays the same ---
"""
Benchmark Auto-Analysis — IDE-friendly, self-contained script (v0.3)

What this does
- Crawls: DOCKED/<PDB>/bench_pocketX_single/ (e.g., bench_pocket1_single, bench_pocket2_single)
- Finds controls like: NIL_A601_bench_pocket2_single.best.pdb
- Finds RDK library poses like: rdk_0003792_bench_pocket2_single.pdbqt
- Parses Vina scores from REMARK lines, computes pose centroids, and an approximate RMSD.
- Pairs each CONTROL to the best RDK match (prefer identity when mapping name is available, else centroid-nearest).
- Scores each pair up to **4 points**:
    1) |Δ score| ≤ 2.0 kcal/mol
    2) Centroid distance ≤ 1.0 Å
    3) Pose–pose RMSD ≤ 3.0 Å
    4) Identity match (RDKit ligand resolves to the *same drug* as the control via mapping or alias)
- Writes two CSVs under DOCKED/_analysis/ by default:
    • benchmark_analysis_details.csv  (one row per control pairing; now includes pass_category)
    • benchmark_analysis_summary.csv  (one row per PDB; now includes passed/near/fail flags)

Notes
- This file is **standalone**: no project-internal imports are required. Only dependency is NumPy.
- Identity match (#4) improves when you provide your FDA mapping CSV. If omitted, #4 usually = 0 (fallback to nearest).
- RMSD uses a robust, order-agnostic approximation (nearest-neighbor pairing + Kabsch alignment) for heavy atoms only.
- “Pass / Near / Fail” in summary is computed from **best_total_points** per PDB:
    • Pass  = best_total_points == max_points (4 with identity, 3 without)
    • Near  = 1 < best_total_points < max_points
    • Fail  = best_total_points == 0
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import base64
import html as _html
from dataclasses import dataclass
from pathlib import Path
import re, csv, json, math, logging, itertools, statistics
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Iterable, Set, Sequence, Mapping, List, Sequence
import numpy as np
from chemdb.chem_alias_db import alias_list_for_het as _alias_list_for_het



# -----------------------------
# Portable defaults: prefer env, then relative to this file
from pathlib import Path
_SCRIPT_ROOT = Path(__file__).resolve().parent

DEFAULT_DOCKED_ROOT  = os.environ.get("DOCKED_DIR") \
    or os.environ.get("OUTPUT_DIR") \
    or str(_SCRIPT_ROOT / "docked")

DEFAULT_MAPPING_CSV  = os.environ.get("MAPPING_CSV") \
    or str(_SCRIPT_ROOT / "fda_mapping_from_pdbqt.csv")
DEFAULT_ONLY_PDB = None  # e.g., "5MO4"
DEFAULT_SCORE_TOL = 2   # kcal/mol
DEFAULT_CENTER_TOL = 1.0  # Å
DEFAULT_RMSD_TOL = 3.0    # Å
DEFAULT_OUTDIR = None     # None -> <DOCKED>\_analysis
DEFAULT_INCLUDE_IDENTITY = True  # toggle 4th criterion
# Prefer explicit analysis cutoff if provided; else use the run's stamped epoch
DEFAULT_SINCE_EPOCH = os.environ.get("BENCH_ANALYSIS_SINCE_EPOCH") or os.environ.get("BENCH_RUN_START_EPOCH")

try:
    from input_and_export_functions import load_config, validate_config
    _cfg = load_config("config.txt") or {}
    try:
        validate_config(_cfg)
    except Exception:
        pass
except Exception:
    _cfg = {}

DEFAULT_DOCKED_ROOT = _cfg.get("DOCKED_DIR") or DEFAULT_DOCKED_ROOT
DEFAULT_MAPPING_CSV = _cfg.get("FDA_MAPPING_CSV") or DEFAULT_MAPPING_CSV
DEFAULT_ONLY_PDB     = _cfg.get("ANALYSIS_ONLY_PDB") or DEFAULT_ONLY_PDB
DEFAULT_SCORE_TOL    = float(_cfg.get("SCORE_TOL_KCAL", DEFAULT_SCORE_TOL))
DEFAULT_CENTER_TOL   = float(_cfg.get("CENTER_TOL_A",   DEFAULT_CENTER_TOL))
DEFAULT_RMSD_TOL     = float(_cfg.get("RMSD_TOL_A",     DEFAULT_RMSD_TOL))
DEFAULT_OUTDIR       = _cfg.get("BENCH_ANALYSIS_OUTDIR") or DEFAULT_OUTDIR
VARIANT = (os.environ.get("APO_HOLO_MODE") or str(_cfg.get("APO_HOLO_MODE") or "holo")).strip().lower()

# -----------------------------
# Filename patterns
# -----------------------------
SCORE_PAT = re.compile(r"REMARK\s+VINA\s+RESULT[:\s]+(-?\d+\.\d+)")
CONTROL_PAT = re.compile(r"^(?P<het>[A-Za-z0-9]{3})_\w\d+_bench_pocket\d+_single\.best\.(?:pdb|pdbqt)$", re.I)
RDK_PAT = re.compile(r"^(rdk_\d{6,8})_bench_pocket\d+_single\.(?:pdbqt|pdb)$", re.I)


# Expected RDK names by PDB (case-insensitive match will be applied on use)
EXPECTED_RDK_BY_PDB: Dict[str, List[str]] = {
  "1T46": ["imatinib", "gleevec"],
  "1M17": ["erlotinib"],
  "2E2B": ["bafetinib"],
  "2HYY": ["imatinib", "gleevec"],
  "2GQG": ["dasatinib"],
  "3CS9": ["nilotinib", "tasigna"],
  "3ERT": ["4-hydroxytamoxifen", "tamoxifen"],
  "3OG7": ["vemurafenib"],
  "3QX3": ["etoposide"],
  "4RT7": ["quizartinib"],
  "4XUF": ["quizartinib"],
  "6O0L": ["venetoclax"],
  "6JQR": ["gilteritinib"],
  "6WTN": ["ruxolitinib"],
  "4XV2": ["dabrafenib"],
  "3OXZ": ["ponatinib"],
  "5I96": ["enasidenib"],
  "4U5J": ["ruxolitinib"],
  "3LXK": ["tofacitinib", "xeljanz"],
  "3ZOS": ["ponatinib", "Ponatinib"],
  "3WZD": ["lenvatinib"],
  "2WGJ": ["crizotinib"],
  "2XP2": ["crizotinib"],
  "3DZY": ["Rosiglitazone"],
  "3WZE": ["sorafenib"],
  "3ZBF": ["Crizotinib"],
  "4AG8": ["axitinib"],
  "4ASD": ["sorafenib"],
  "5L7I": ["vismodegib", "erivedge"],
  "6O0K": ["venetoclax"],
  "5MO4": ["asciminib", "nilotinib"],
  "6U4J": ["olutasidenib", "FT-2102"],
}



def _is_expected_name(pdb_id: str, rdk_name: str) -> bool:
    exp = [_norm_text(x) for x in EXPECTED_RDK_BY_PDB.get((pdb_id or '').upper(), [])]
    nm = _norm_text(rdk_name or '')
    if not exp or not nm:
        return True  # fail-open
    return any(en and (en in nm or nm in en) for en in exp)

# name normalization for robust identity checks
from typing import Optional as _Optional
def _norm_text(s: _Optional[str]) -> str:
    return re.sub(r'[^a-z0-9]+', '', (s or '').lower())

# Debug flag (env or CLI --debug)
DEBUG = bool(int(os.environ.get("BENCH_DEBUG", "0")))

def _pdbqt_element_histogram(path: Path) -> Dict[str, int]:
    hist: Dict[str, int] = {}
    try:
        with open(path, "r", errors="ignore") as f:
            for ln in f:
                if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                    continue
                el = (ln[76:78].strip() if len(ln) >= 78 else "").upper()
                if not el:
                    # fallback: atom name heuristic
                    an = ln[12:16].strip().upper()
                    el = (an[0] if an else "")
                if el:
                    hist[el] = hist.get(el, 0) + 1
    except Exception:
        pass
    return hist

def _hist_equal(a: Dict[str,int], b: Dict[str,int]) -> bool:
    # ignore H in strict heavy atom match; compare H optionally
    ah = {k:v for k,v in a.items() if k != "H"}
    bh = {k:v for k,v in b.items() if k != "H"}
    return ah == bh
def alignment_metrics(ctrl: np.ndarray, rdk: np.ndarray, nn_cap: int = 128):
    """
    Returns (rmsd, raw_centroid, aligned_centroid, n_ctrl, n_rdk).
    Aligns RDK -> CTRL using nearest-neighbor pairing + Kabsch.
    """
    n_ctrl = int(ctrl.shape[0])
    n_rdk  = int(rdk.shape[0])

    rmsd = None

    if n_ctrl >= 1 and n_rdk >= 1:
        c_ctr = ctrl.mean(axis=0)
        r_ctr = rdk.mean(axis=0)

    if n_ctrl < 3 or n_rdk < 3:
        return rmsd, n_ctrl, n_rdk

    # build NN pairs (A=ctrl reference, B=rdk moved to A)
    A = ctrl
    B = rdk
    if A.shape[0] > nn_cap:
        idx = np.linspace(0, A.shape[0] - 1, nn_cap, dtype=int)
        A = A[idx]

    used = set()
    pairs_a, pairs_b = [], []
    for a in A:
        best_j, best_d2 = None, 1e99
        for j, b in enumerate(B):
            if j in used:
                continue
            d2 = ((a - b) ** 2).sum()
            if d2 < best_d2:
                best_d2, best_j = d2, j
        if best_j is not None:
            used.add(best_j)
            pairs_a.append(a)
            pairs_b.append(B[best_j])

    if len(pairs_a) < 3:
        return rmsd, n_ctrl, n_rdk

    P = np.vstack(pairs_a)  # ctrl
    Q = np.vstack(pairs_b)  # rdk

    # R that maps Q -> P (note the reversed order):
    R, _ = _kabsch(Q, P)
    Qc = Q - Q.mean(axis=0)
    Pc = P - P.mean(axis=0)
    # RMSD after alignment
    diff2 = ((Qc @ R) - Pc) ** 2
    rmsd = float(np.sqrt(diff2.sum(axis=1).mean()))

    # Apply transform to ALL rdk atoms for centroid-after-alignment
    rdk_aligned = (rdk - Q.mean(axis=0)) @ R + P.mean(axis=0)
    return rmsd, n_ctrl, n_rdk


# -----------------------------
# Mapping CSV (optional, no pandas needed)
# -----------------------------
@dataclass
class MapRow:
    path: str
    display_name: str = ""
    generic_name: str = ""
    brand_names: str = ""
    pubchem_name: str = ""
    pubchem_record_title: str = ""
    pubchem_iupac_name: str = ""
    pubchem_synonyms: str = ""
    rxnorm_generic_name: str = ""
    rxnorm_brand_names: str = ""
    drugcentral_generic_name: str = ""
    drugcentral_brand_names: str = ""
    remark_name: str = ""
    sdf_title: str = ""
    inchikey: str = ""

    def best_name(self) -> str:
        for f in (
            self.display_name,
            self.generic_name,
            self.remark_name,
            self.sdf_title,
            self.brand_names,
            self.pubchem_name,
            self.pubchem_record_title,
        ):
            if f and str(f).strip():
                return str(f).strip()
        return ""


class MappingIndex:
    def __init__(self, csv_path: Optional[Path]):
        self.rows: List[MapRow] = []
        if not csv_path or not Path(csv_path).is_file():
            return
        try:
            with open(csv_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
                reader = csv.DictReader(f)
                for r in reader:
                    self.rows.append(
                        MapRow(
                            path=str(r.get("path", "")),
                            display_name=str(r.get("display_name", "")),
                            generic_name=str(r.get("generic_name", "")),
                            brand_names=str(r.get("brand_names", "")),
                            pubchem_name=str(r.get("pubchem_name", "")),
                            pubchem_record_title=str(r.get("pubchem_record_title", "")),
                            pubchem_iupac_name=str(r.get("pubchem_iupac_name", "")),
                            pubchem_synonyms=str(r.get("pubchem_synonyms", "")),
                            rxnorm_generic_name=str(r.get("rxnorm_generic_name", "")),
                            rxnorm_brand_names=str(r.get("rxnorm_brand_names", "")),
                            drugcentral_generic_name=str(r.get("drugcentral_generic_name", "")),
                            drugcentral_brand_names=str(r.get("drugcentral_brand_names", "")),
                            remark_name=str(r.get("remark_name", "")),
                            sdf_title=str(r.get("sdf_title", "")),
                            inchikey=str(r.get("inchikey", "")),
                        )
                    )
        except Exception as e:
            print(f"[analysis] WARNING: failed to parse mapping CSV: {e}")

    def resolve_name_for_rdk(self, rdk_id: str) -> Optional[str]:
        rid = (rdk_id or "").strip().lower()
        if not rid:
            return None
        best: Optional[str] = None
        for row in self.rows:
            p = str(row.path or "").lower()
            if rid in p:
                n = row.best_name()
                if n:
                    best = n
                    break
        return best


# -----------------------------
# PDB/PDBQT parsing helpers
# -----------------------------

def parse_vina_score(path: Path) -> Optional[float]:
    try:
        with open(path, "r", errors="ignore") as f:
            for line in f:
                m = SCORE_PAT.search(line)
                if m:
                    return float(m.group(1))
    except Exception:
        pass
    return None


def parse_coords(path: Path) -> np.ndarray:
    """Return Nx3 coords for heavy atoms from PDB/PDBQT (order-agnostic)."""
    coords: List[Tuple[float, float, float]] = []
    try:
        with open(path, "r", errors="ignore") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                atom_name = line[12:16].strip() if len(line) >= 16 else ""
                if atom_name.upper().startswith("H"):
                    continue
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except Exception:
                    parts = line.split()
                    flts = []
                    for p in parts:
                        try:
                            flts.append(float(p))
                        except Exception:
                            pass
                    if len(flts) >= 3:
                        x, y, z = flts[-3:]
                    else:
                        continue
                coords.append((x, y, z))
    except Exception:
        pass
    if not coords:
        return np.zeros((0, 3), dtype=float)
    return np.asarray(coords, dtype=float)


# --- Prepped ligand path resolver (module-scope) ------------------------------

def _collapse_core_stem(stem: str) -> str:
    # strip stage suffixes like "_bench_pocket1_single" and optional ".best"
    s = re.sub(r"_bench_pocket\d+_single(?:\.best)?$", "", stem, flags=re.I)
    # strip trailing ".sanitized" tokens if present
    s = re.sub(r"(?:\.sanitized)+$", "", s, flags=re.I)
    return s

def _prepped_search_dirs_for(pdb_id: str) -> list[Path]:
    """
    Likely directories that contain the *input* prepped ligands for a given PDB.
    Uses config/env fallbacks consistent with the rest of this script.
    """
    base_overall = Path(DEFAULT_DOCKED_ROOT).parent
    processed_root = Path(_cfg.get("OUTPUT_DIR") or (base_overall / "processed_pdbs"))
    lib_root = Path(_cfg.get("OUTPUT_LIGANDS_DIR") or (base_overall / "prepped_ligands"))

    dirs: list[Path] = []
    # per-PDB prepped ligands
    dirs.append(processed_root / pdb_id / "prepped_ligands")
    # library roots (if user keeps FDA libraries here)
    for sub in ("", "fda_library", "fda_library2", "fda_new"):
        d = lib_root if not sub else (lib_root / sub)
        if d.is_dir():
            dirs.append(d)

    # de-dup, keep order
    out, seen = [], set()
    for d in dirs:
        try:
            key = str(d.resolve())
        except Exception:
            key = str(d)
        if key not in seen:
            seen.add(key)
            out.append(d)
    return out

def guess_prepped_from_pose(pose_path: Path, pdb_id: str) -> Optional[Path]:
    """
    Best-effort guess for the *input* prepped ligand (.pdbqt) that produced a pose file.
    Works for both controls and RDKs by stem matching and searching known roots.
    """
    core = _collapse_core_stem(pose_path.stem)

    # 1) try exact file name in likely dirs
    for d in _prepped_search_dirs_for(pdb_id):
        q = d / f"{core}.pdbqt"
        if q.is_file():
            return q
        q2 = d / f"{core}.sanitized.pdbqt"
        if q2.is_file():
            return q2

    # 2) looser glob as fallback
    for d in _prepped_search_dirs_for(pdb_id):
        hits = sorted(d.glob(f"*{core}*.pdbqt"))
        if hits:
            return hits[0]
    return None

def extract_raw_text_block(path: Path, max_chars: int = 20000, only_header: bool = False) -> str:
    """
    Return a 'raw text' snippet from a PDB/PDBQT pose.

    If only_header=True:
      - Keep header-ish lines (MODEL, REMARK, TORSDOF, ENDMDL; plus ROOT marker for PDBQT).
      - Skip ATOM/HETATM coordinate blocks to stay compact.

    If only_header=False (default for our HTML/CSV fields now):
      - Include the *entire* file contents (including ATOM/HETATM/BRANCH blocks).
      - Still enforce a max_chars cap to keep reports reasonable.

    The result is newline-preserved; caller can wrap in <pre> for HTML.
    """
    try:
        ext = (path.suffix or "").lower()
        is_pdbqt = (ext == ".pdbqt")

        if not only_header:
            # Full text path: include everything
            with open(path, "r", errors="ignore") as f:
                raw = f.read()
            raw = raw.strip()
            if len(raw) > max_chars:
                raw = raw[:max_chars].rstrip() + " …"
            return raw
        
        # Header-only path (kept for possible re-use)
        keep_prefixes = ("MODEL", "REMARK", "TORSDOF", "ENDMDL")
        out_lines = []
        with open(path, "r", errors="ignore") as f:
            for ln in f:
                if any(ln.startswith(pfx) for pfx in keep_prefixes):
                    out_lines.append(ln.rstrip("\r\n"))
                    continue
                if is_pdbqt and ln.startswith("ROOT"):
                    out_lines.append("ROOT")
                    continue
                # header-only mode deliberately skips ATOM/HETATM
        raw = "\n".join(out_lines).strip()
        if len(raw) > max_chars:
            raw = raw[:max_chars].rstrip() + " …"
        return raw
    except Exception:
        return ""



def centroid(pts: np.ndarray) -> Optional[np.ndarray]:
    if pts is None or pts.size == 0:
        return None
    return pts.mean(axis=0)


# -----------------------------
# Approximate RMSD (robust, order-agnostic)
# -----------------------------

def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, float]:
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    P_rot = Pc @ R
    diff2 = ((P_rot - Qc) ** 2).sum(axis=1)
    rmsd = float(np.sqrt(diff2.mean())) if len(diff2) else float("inf")
    return R, rmsd


def approximate_rmsd(coords_a: np.ndarray, coords_b: np.ndarray, nn_cap: int = 128) -> Optional[float]:
    if coords_a.shape[0] < 3 or coords_b.shape[0] < 3:
        return None
    # Use smaller set as reference; cap to nn_cap for speed
    A, B = (coords_a, coords_b) if coords_a.shape[0] <= coords_b.shape[0] else (coords_b, coords_a)
    if A.shape[0] > nn_cap:
        idx = np.linspace(0, A.shape[0] - 1, nn_cap, dtype=int)
        A = A[idx]
    used = set()
    pairs_a: List[np.ndarray] = []
    pairs_b: List[np.ndarray] = []
    for a in A:
        best_j = None
        best_d2 = 1e99
        for j, b in enumerate(B):
            if j in used:
                continue
            d2 = ((a - b) ** 2).sum()
            if d2 < best_d2:
                best_d2 = d2
                best_j = j
        if best_j is None:
            continue
        used.add(best_j)
        pairs_a.append(a)
        pairs_b.append(B[best_j])
    if len(pairs_a) < 3:
        return None
    P = np.vstack(pairs_a)
    Q = np.vstack(pairs_b)
    _, rmsd = _kabsch(P, Q)
    return rmsd


# -----------------------------
# Identity helpers
# -----------------------------

def _rdk_id_from_stem(stem: str) -> Optional[str]:
    m = re.search(r"(rdk_\d{6,8})", stem.lower())
    return m.group(1) if m else None


def resolve_rdk_name_from_mapping(rdk_id: str, mapping: MappingIndex) -> Optional[str]:
    if not mapping:
        return None
    try:
        return mapping.resolve_name_for_rdk(rdk_id)
    except Exception:
        return None


def expected_names_for_het(het: str) -> List[str]:
    """
    Return the synonym list for a 3–5 letter HET code, backed by chemdb/aliases.yaml.
    Tokens are already normalized to lowercase by chemdb.chem_alias_db.
    """
    return _alias_list_for_het(het or "")



# -----------------------------
# Data classes
# -----------------------------
@dataclass
class Pose:
    path: Path
    score: Optional[float]
    coords: np.ndarray
    rdk_name: Optional[str] = None  # resolved drug name for RDKs


# -----------------------------
# Screenshot helpers (PyMOL renders)
# -----------------------------

def find_pymol_screenshots(pdb_id: str, pocket_dir: Path) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """
    Looks for three PNGs in the given pocket directory:
      <PDB>_cleaned__CONTROL+RDKclosest_side.png
      <PDB>_cleaned__CONTROL+RDKclosest_front.png
      <PDB>_cleaned__CONTROL+RDKclosest_top.png
    Returns (side, front, top) Paths if found, else None entries.
    """
    base = f"{pdb_id}_cleaned__CONTROL+RDKclosest"
    views = ["side", "front", "top"]
    out: List[Optional[Path]] = []
    for v in views:
        p = pocket_dir / f"{base}_{v}.png"
        if p.is_file():
            out.append(p)
            continue
        # Fallback: looser glob in case of minor naming differences
        cands = sorted(pocket_dir.glob(f"*{pdb_id}*CONTROL+RDKclosest*{v}*.png"))
        out.append(cands[0] if cands else None)
    return out[0], out[1], out[2]
def find_pair_only_screenshots(pdb_id: str, pocket_dir: Path, rdk_id: Optional[str]) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """
    Looks for three PNGs named like:
      <PDB>_cleaned__CONTROL+<RDKID>_PAIR_side.png
      <PDB>_cleaned__CONTROL+<RDKID>_PAIR_front.png
      <PDB>_cleaned__CONTROL+<RDKID>_PAIR_top.png
    Falls back to a loose glob if names vary.
    """
    rid = (rdk_id or "RDKclosest")
    base = f"{pdb_id}_cleaned__CONTROL+{rid}_PAIR"
    views = ["side", "front", "top"]
    out: List[Optional[Path]] = []
    for v in views:
        p = pocket_dir / f"{base}_{v}.png"
        if p.is_file():
            out.append(p)
            continue
        # Fallback: allow slight naming differences, but keep rid in the filename if possible
        cands = [q for q in pocket_dir.glob(f"*{pdb_id}*CONTROL+*PAIR*{v}*.png") if rid.lower() in q.name.lower()]
        if not cands:
            cands = list(pocket_dir.glob(f"*{pdb_id}*CONTROL+*PAIR*{v}*.png"))
        out.append(sorted(cands)[0] if cands else None)
    return out[0], out[1], out[2]


@dataclass
class PairEval:
    pdb_id: str
    pocket: str
    control_id: str
    control_het: str
    control_file: Path
    rdk_file: Path
    rdk_id: str

    # fields with defaults must come after all required fields
    control_prepped: Optional[Path] = None
    rdk_prepped: Optional[Path] = None
    rdk_name: Optional[str] = None
    control_score: Optional[float] = None
    rdk_score: Optional[float] = None

    delta_kcal: Optional[float] = None
    rmsd: Optional[float] = None
    flag_score: int = 0
    flag_center: int = 0
    flag_rmsd: int = 0
    flag_identity: int = 0
    total_points: int = 0
    n_atoms_ctrl: int = 0
    n_atoms_rdk: int = 0
    pick_reason: str = ""

    png_side: Optional[Path] = None
    png_front: Optional[Path] = None
    png_top: Optional[Path] = None

    # interpretability/readability
    control_display_name: Optional[str] = None
    flags_sum: int = 0
    confidence_label: str = ""
    confident_match: int = 0
    good_control: int = 0
    png_pair_side: Optional[Path] = None
    png_pair_front: Optional[Path] = None
    png_pair_top: Optional[Path] = None
    expected_rdks: str = ""
    matched_success: int = 0
    rdk_raw_text: str = ""
    control_raw_text: str = ""
    identity_strict: int = 0
    ligprep_errors: str = ""
    proteinprep_errors: str = ""
    is_expected: int = 0


# -----------------------------
# Crawling & pairing
# -----------------------------

def find_pocket_dirs(pdb_dir: Path) -> List[Path]:
    out: List[Path] = []
    if not pdb_dir.is_dir():
        return out
    for child in pdb_dir.iterdir():
        if child.is_dir() and child.name.startswith("bench_pocket") and child.name.endswith("_single"):
            out.append(child)
    return sorted(out)

def load_controls_and_rdks(
    pocket_dir: Path,
    mapping: MappingIndex,
    since_epoch: Optional[float] = None
) -> Tuple[List[Tuple[str, Pose]], List[Pose]]:
    controls: List[Tuple[str, Pose]] = []
    rdks: List[Pose] = []
    cutoff = float(since_epoch) if since_epoch is not None else None
    for p in pocket_dir.iterdir():
        if not p.is_file():
            continue
        if cutoff is not None:
            try:
                if float(p.stat().st_mtime) < cutoff:
                    continue
            except Exception:
                pass
        name = p.name
        if CONTROL_PAT.match(name):
            control_id = name.split("_bench_")[0]  # NIL_A601
            coords = parse_coords(p)
            sc = parse_vina_score(p)
            controls.append((control_id, Pose(path=p, score=sc, coords=coords)))
        elif RDK_PAT.match(name):
            coords = parse_coords(p)
            sc = parse_vina_score(p)
            rdk_id = _rdk_id_from_stem(Path(name).stem)
            rname = resolve_rdk_name_from_mapping(rdk_id or "", mapping)
            rdks.append(Pose(path=p, score=sc, coords=coords, rdk_name=rname))
    return controls, rdks


#should stay unused
def _closest_by_centroid(ctrl: Pose, rdks: Sequence[Pose]) -> Optional[Pose]:
    if ctrl.centroid is None:
        return None
    best = None
    best_d = 1e99
    for r in rdks:
        if r.centroid is None:
            continue
        d = float(np.linalg.norm(ctrl.centroid - r.centroid))
        if d < best_d:
            best_d = d
            best = r
    return best



def _match_by_identity(ctrl_het: str, rdks: Sequence[Pose]) -> List[Pose]:
    expected_norm = [_norm_text(s) for s in expected_names_for_het(ctrl_het)]
    if not expected_norm:
        return []
    hits: List[Pose] = []
    for r in rdks:
        nm_norm = _norm_text(r.rdk_name)
        if nm_norm and any(en and (en in nm_norm or nm_norm in en) for en in expected_norm):
            hits.append(r)
    return hits


def evaluate_pairs(
    pdb_id: str,
    pocket_name: str,
    pocket_dir: Path,
    controls: Sequence[Tuple[str, Pose]],
    rdks: Sequence[Pose],
    score_tol: float,
    center_tol: float,
    rmsd_tol: float,
    include_identity: bool,
) -> List[PairEval]:
    # Locate the PyMOL screenshots once per pocket
    png_side, png_front, png_top = find_pymol_screenshots(pdb_id, pocket_dir)

    rows: List[PairEval] = []

    # lazy load sidecar status TSVs from sibling dirs (optional/minimal wiring)
    def _lookup_errors(p: Path) -> Tuple[str, str]:
        # try to find <pdb>/ligands_raw/ligand_prep_status.tsv and <pdb>/_NOLIG/receptor/*.log
        # Minimal: read a TSV that has 'ligprep_errors' and 'proteinprep_errors' keyed by basename
        return "", ""  # no-op if not available


    for control_id, ctrl in controls:
        m = re.match(r"^(?P<het>[A-Za-z0-9]{3})_(?P<chain>\w)(?P<res>\d+)$", control_id)
        het = (m.group("het") if m else control_id.split("_")[0]).upper()

        # Prefer identity match; fallback to score-nearest if none
        picked = None
        identity_flag = 0
        candidates = _match_by_identity(het, rdks)
        if candidates:
            identity_flag = 1
            if ctrl.score is not None:
                picked = min(candidates, key=lambda r: abs((r.score or 9e9) - (ctrl.score or 0.0)))
            else:
                picked = min(candidates, key=lambda r: (r.score or 9e9))
        else:
            # score-based fallback (no identity): keep rows flowing without centroid
            if not rdks:
                continue  # nothing to pair with
            if ctrl.score is not None:
                picked = min(rdks, key=lambda r: abs((r.score or 9e9) - (ctrl.score or 0.0)))
                pick_reason = "score-nearest"
            else:
                picked = min(rdks, key=lambda r: (r.score or 9e9))
                pick_reason = "best-score"

        if picked is None:
            continue

        pick_reason = "identity"
        identity_strict = 0
        try:
            ctrl_hist = _pdbqt_element_histogram(ctrl.path)
            rdk_hist = _pdbqt_element_histogram(picked.path)
            if _hist_equal(ctrl_hist, rdk_hist):
                identity_strict = 1
        except Exception:
            identity_strict = 0
        # alias for consistency (fixes NameError in confidence calc)
        flag_identity = identity_flag

        # Score delta
        # Score delta
        delta_kcal: Optional[float] = None
        flag_score = 0
        if ctrl.score is not None and picked.score is not None:
            delta_kcal = abs(picked.score - ctrl.score)
            flag_score = int(delta_kcal <= score_tol)

        # Geometry metrics (RMSD only; no centroid flag)
        rmsd_val, n_ctrl, n_rdk = alignment_metrics(ctrl.coords, picked.coords)
        flag_rmsd = int((rmsd_val is not None) and (rmsd_val <= rmsd_tol))

        # No centroid flag in this version
        flag_center = 0

        if DEBUG and rmsd_val is not None and flag_rmsd:
            print(f"[debug] {pdb_id}/{pocket_name}/{control_id}: RMSD ok ({rmsd_val:.2f} Å) "
                  f"Using aligned for scoring. pick={picked.path.name} via {pick_reason}")

        rdk_id = _rdk_id_from_stem(picked.path.stem) or ""
        # Expected RDK names (by PDB)
        expected_list = EXPECTED_RDK_BY_PDB.get(pdb_id.upper(), [])
        expected_rdks = ", ".join(expected_list)

        # HTML filter flag: does picked.rdk_name match the allow-list for this PDB?
        exp_norm = [_norm_text(x) for x in EXPECTED_RDK_BY_PDB.get(pdb_id.upper(), [])]
        nm_norm = _norm_text(picked.rdk_name)
        is_expected_bool = 0
        if not exp_norm or not nm_norm:
            # fail-open: if no allow-list for this PDB or name missing, include row
            is_expected_bool = 1
        else:
            for en in exp_norm:
                if en and (en in nm_norm or nm_norm in en):
                    is_expected_bool = 1
                    break

        # interpretability / labels (resolve control name first)
        ctrl_guess_list = expected_names_for_het(het)
        control_display_name = (ctrl_guess_list[0] if ctrl_guess_list else "")

        # Raw text straight from pose files INCLUDING ATOM/HETATM (capped for size)
        rdk_raw_text = extract_raw_text_block(picked.path, only_header=False)
        control_raw_text = extract_raw_text_block(ctrl.path, only_header=False)


        # Success criterion (ignore centroid): identity & RMSD & score
        matched_success = int(bool(flag_identity and flag_rmsd and flag_score))

        pair_side, pair_front, pair_top = find_pair_only_screenshots(pdb_id, pocket_dir, rdk_id or "RDKclosest")

        # Confidence rule use only score & RMSD (+ identity for "confident")
        _two = (flag_rmsd + flag_score)
        confident = int(bool(flag_identity and flag_rmsd and flag_score))
        confidence_label = "confident" if confident else ("plausible" if _two >= 2 else "weak")


        good_control = int((ctrl.score is not None) and (int(ctrl.coords.shape[0]) >= 10) and bool(control_display_name))
        flags_sum = _two + (identity_flag if include_identity else 0)
        if DEBUG and (ctrl.score is not None) and (picked.score is not None):
            logging.debug(
                "[Δkcal-check] %s/%s ctrl=%.3f rdk=%.3f Δ=%.3f tol=%.3f -> score✔=%d",
                pdb_id, control_id, ctrl.score, picked.score, delta_kcal, score_tol, flag_score
            )

        total_points = flag_score + flag_center + flag_rmsd + (identity_flag if include_identity else 0)
        lig_e, prot_e = _lookup_errors(picked.path)
        # resolve the actual prepped input ligands (best-effort)
        _control_prepped = guess_prepped_from_pose(ctrl.path, pdb_id)
        _rdk_prepped     = guess_prepped_from_pose(picked.path, pdb_id)

        rows.append(
            PairEval(
                pdb_id=pdb_id,
                pocket=pocket_name,
                control_id=control_id,
                control_het=het,
                control_file=ctrl.path,
                rdk_file=picked.path,
                control_prepped=_control_prepped,
                rdk_prepped=_rdk_prepped,
                rdk_id=rdk_id,
                rdk_name=picked.rdk_name,
                control_score=ctrl.score,
                rdk_score=picked.score,
                delta_kcal=delta_kcal,
                rmsd=rmsd_val,
                flag_score=flag_score,
                flag_center=flag_center,
                flag_rmsd=flag_rmsd,
                flag_identity=identity_flag,
                total_points=total_points,
                n_atoms_ctrl=int(ctrl.coords.shape[0]),
                n_atoms_rdk=int(picked.coords.shape[0]),

                pick_reason=pick_reason,
                png_side=png_side,
                png_front=png_front,
                png_top=png_top,
                control_display_name=control_display_name,
                flags_sum=flags_sum,
                confidence_label=confidence_label,
                confident_match=confident,
                good_control=good_control,
                png_pair_side=pair_side,
                png_pair_front=pair_front,
                png_pair_top=pair_top,
                expected_rdks=expected_rdks,
                matched_success=matched_success,
                rdk_raw_text=rdk_raw_text,
                control_raw_text=control_raw_text,
                identity_strict=identity_strict,
                ligprep_errors=lig_e, proteinprep_errors=prot_e,
                is_expected=int(is_expected_bool),

            )
        )
    return rows


# -----------------------------
# CSV writing
# -----------------------------

def _fmt(x: Optional[float]) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return ""
    try:
        return f"{float(x):.3f}"
    except Exception:
        return ""


def _pass_category(total_points: int, include_identity: bool) -> str:
    max_points = 4 if include_identity else 3
    if total_points == max_points:
        return "pass"
    if 1 < total_points < max_points:
        return "near"
    if total_points == 0:
        return "fail"
    return ""  # exactly 1 -> no category requested


def write_details_csv(out_dir: Path, rows: Sequence[PairEval], include_identity: bool) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "benchmark_analysis_details.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "variant", "pdb_id", "pocket", "control_id", "control_het",
            "control_file", "control_prepped", "rdk_file", "rdk_prepped", "rdk_id",
            "control_drug_guess", "rdk_suspected_fda_name",
            "control_score_kcal", "rdk_score_kcal", "delta_kcal",
            "rmsd_ctrl_to_rdk_A", "rmsd_rdk_to_ctrl_A",
            "score_within_tol", "rmsd_within_tol", "identity_match",
            "n_atoms_ctrl", "n_atoms_rdk", "pick_reason",
            "png_side", "png_front", "png_top",
            "pair_png_side", "pair_png_front", "pair_png_top",
            "flags_sum", "confident_match", "confidence_label", "good_control",
            "total_points", "pass_category",
            "ligprep_errors", "proteinprep_errors", "identity_strict",
            "expected_rdks", "matched_success", "rdk_raw_text", "control_raw_text",
        ])
        for r in rows:
            w.writerow([
                VARIANT,
                r.pdb_id, r.pocket, r.control_id, r.control_het,
                str(r.control_file), (str(r.control_prepped) if r.control_prepped else ""),
                str(r.rdk_file), (str(r.rdk_prepped) if r.rdk_prepped else ""),
                r.rdk_id,
                (r.control_display_name or ""), (r.rdk_name or ""),
                _fmt(r.control_score), _fmt(r.rdk_score), _fmt(r.delta_kcal),
                _fmt(r.rmsd), _fmt(r.rmsd),
                r.flag_score,  r.flag_rmsd, r.flag_identity,
                r.n_atoms_ctrl, r.n_atoms_rdk, r.pick_reason,
                (str(r.png_side) if r.png_side else ""),
                (str(r.png_front) if r.png_front else ""),
                (str(r.png_top) if r.png_top else ""),
                (str(r.png_pair_side) if r.png_pair_side else ""),
                (str(r.png_pair_front) if r.png_pair_front else ""),
                (str(r.png_pair_top) if r.png_pair_top else ""),
                r.flags_sum, r.confident_match, r.confidence_label, r.good_control,
                r.total_points, _pass_category(r.total_points, include_identity),
                (r.ligprep_errors or ""), (r.proteinprep_errors or ""), r.identity_strict,
                r.expected_rdks, r.matched_success, r.rdk_raw_text, r.control_raw_text,

            ])
    return path


# =========================
# Summary writers (drop-in)
# =========================

def write_summary_csv(out_dir: Path, rows: Sequence[PairEval], include_identity: bool) -> Path:
    """
    Per-PDB rollup with flags and counts, including 'good_control' and 'confident_match'.
    Writes: <out_dir>/benchmark_analysis_summary.csv
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    max_points = 4 if include_identity else 3
    path = out_dir / "benchmark_analysis_summary.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "pdb_id",
            "best_total_points",
            "num_controls",
            "num_pockets",
            "any_score_within_tol",
            "any_rmsd_within_tol",
            "any_identity_match",
            "num_good_controls",
            "num_confident_matches",
            "passed_4_or_max",
            "near_pass_2to(max-1)",
            "fail_0",
        ])
        by_pdb: Dict[str, List[PairEval]] = {}
        for r in rows:
            by_pdb.setdefault(r.pdb_id, []).append(r)

        for pdb_id in sorted(by_pdb.keys()):
            items = by_pdb[pdb_id]
            best_total = max((it.total_points for it in items), default=0)
            num_controls = len({it.control_id for it in items})
            num_pockets = len({it.pocket for it in items})
            any_score = int(any(it.flag_score for it in items))
            any_center = int(any(it.flag_center for it in items))
            any_rmsd = int(any(it.flag_rmsd for it in items))
            any_ident = int(any(it.flag_identity for it in items))
            num_good_ctrl = sum(int(it.good_control) for it in items)
            num_confident = sum(int(it.confident_match) for it in items)
            passed = int(best_total == max_points)
            near = int(1 < best_total < max_points)
            fail0 = int(best_total == 0)

            w.writerow([
                pdb_id, best_total, num_controls, num_pockets,
                any_score, any_center, any_rmsd, any_ident,
                num_good_ctrl, num_confident,
                passed, near, fail0,
            ])
    return path


# -----------------------------
# Visual reports (images embedded)
# -----------------------------

def _html_escape(s: Optional[str]) -> str:
    return _html.escape("" if s is None else str(s))


def _img_to_data_uri_or_link(p: Optional[Path], max_width_px: int = 280) -> str:
    if not p or not Path(p).is_file():
        return ""
    try:
        data = Path(p).read_bytes()
        b64 = base64.b64encode(data).decode("ascii")
        return f'<img src="data:image/png;base64,{b64}" style="max-width:{max_width_px}px; height:auto; border:1px solid #ddd;" />'
    except Exception:
        try:
            return f'<img src="{Path(p).as_uri()}" style="max-width:{max_width_px}px; height:auto; border:1px solid #ddd;" />'
        except Exception:
            return _html_escape(str(p))

def _mtime_ok(p: Optional[Path], since_epoch: Optional[float]) -> bool:
    if not p or not Path(p).is_file():
        return False
    if since_epoch is None:
        return True
    try:
        return float(Path(p).stat().st_mtime) >= float(since_epoch)
    except Exception:
        return False

def write_details_html(
    out_dir: Path,
    rows: Sequence[PairEval],
    include_identity: bool,
    *,
    html_success_only: Optional[bool] = None,
    since_epoch: Optional[float] = None,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "benchmark_analysis_details_with_images.html"

    def _badge(label: str) -> str:
        cls = "bad"
        if label == "confident":
            cls = "good"
        elif label == "plausible":
            cls = "ok"
        return f"<span class='badge {cls}'>{_html_escape(label)}</span>"

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("""
<!doctype html>
<html>
<head>
<meta charset="utf-8" />
<title>Benchmark Analysis — Details (with images)</title>
<style>
  body{font-family:Segoe UI,Arial,sans-serif;margin:24px}
  table{border-collapse:collapse;width:100%}
  th,td{border:1px solid #e5e5e5;padding:6px 8px;vertical-align:top;font-size:13px}
  th{background:#fafafa;position:sticky;top:0;z-index:1}
  /* Raw text cells: make tiny and clamp height so images decide row size */
  pre.raw{
    white-space:pre-wrap;
    max-width:420px;
    overflow-wrap:anywhere;
    margin:0;
    font-size:9px;
    line-height:1.05;
    max-height:96px;
    overflow:auto;
  }
  td.rawcol{ padding:2px 4px; }
  .num{text-align:right;white-space:nowrap}
  .imgcell{white-space:nowrap}
  .muted{color:#777}
  .badge{display:inline-block;padding:2px 8px;border-radius:12px;border:1px solid #ddd;font-size:12px}
  .good{background:#e7f6ec;border-color:#c9e8d0;color:#137333}
  .ok{background:#fff8e1;border-color:#f7e0a3;color:#8a6d3b}
  .bad{background:#fdecea;border-color:#f5c6cb;color:#a61b1b}
</style>
</head>
<body>
<h2>Benchmark Analysis — Details (with images)</h2>
<p class="muted">Columns mirror <code>benchmark_analysis_details.csv</code>, plus raw pose text and six image columns.</p>
<table>
<thead><tr>
  <th>variant</th><th>pdb_id</th><th>pocket</th><th>control_id</th><th>rdk_id</th>
  <th>rdk_prepped</th><th>control_prepped</th>
  <th>control_drug</th><th>rdk_suspected_fda</th>
  <th class="num">ctrl_score</th><th class="num">rdk_score</th><th class="num">Δkcal</th>
  <th class="num">RMSD_A</th>
  <th>expected_rdks</th><th>matched_success</th><th>rdk_raw_text</th><th>control_raw_text</th>
  <th>pick</th><th>confidence</th><th>good control</th>
  <th class="num">score✔</th><th class="num">rmsd✔</th><th class="num">ident✔</th>
  <th class="imgcell">side</th><th class="imgcell">front</th><th class="imgcell">top</th>
  <th class="imgcell">pair_side</th><th class="imgcell">pair_front</th><th class="imgcell">pair_top</th>
</tr></thead>
<tbody>
""")

        # Do NOT filter rows: render everything provided in `rows`.
        for r in rows:
            # Gate images by mtime ≥ since_epoch (rows are never filtered)
            def _maybe(p: Optional[Path]) -> Optional[Path]:
                return p if _mtime_ok(p, since_epoch) else None
            side = _maybe(r.png_side)
            front = _maybe(r.png_front)
            top = _maybe(r.png_top)
            pside = _maybe(r.png_pair_side)
            pfront = _maybe(r.png_pair_front)
            ptop = _maybe(r.png_pair_top)

            fh.write("<tr>")
            fh.write(f"<td>{_html_escape(VARIANT.upper())}</td>")
            fh.write(f"<td>{_html_escape(r.pdb_id)}</td>")
            fh.write(f"<td>{_html_escape(r.pocket)}</td>")
            fh.write(f"<td>{_html_escape(r.control_id)}</td>")
            fh.write(f"<td>{_html_escape(r.rdk_id)}</td>")
            fh.write(f"<td>{_html_escape(str(r.rdk_prepped or ''))}</td>")
            fh.write(f"<td>{_html_escape(str(r.control_prepped or ''))}</td>")
            fh.write(f"<td>{_html_escape(r.control_display_name or '')}</td>")

            _nm = r.rdk_name or ''
            def _norm(s: str) -> str:
                return re.sub(r'[^a-z0-9]+', '', (s or '').lower())
            _exp = [_norm(x) for x in EXPECTED_RDK_BY_PDB.get((r.pdb_id or '').upper(), [])]
            _nm_norm = _norm(_nm)
            _match = (not _exp) or (not _nm_norm) or any(en and (en in _nm_norm or _nm_norm in en) for en in _exp)
            note = '' if _match else " <span class='muted'>(not-in-expected)</span>"
            fh.write(f"<td>{_html_escape(_nm)}{note}</td>")

            fh.write(f"<td class='num'>{_html_escape(_fmt(r.control_score))}</td>")
            fh.write(f"<td class='num'>{_html_escape(_fmt(r.rdk_score))}</td>")
            fh.write(f"<td class='num'>{_html_escape(_fmt(r.delta_kcal))}</td>")
            fh.write(f"<td class='num'>{_html_escape(_fmt(r.rmsd))}</td>")
            fh.write(f"<td>{_html_escape(r.expected_rdks)}</td>")
            fh.write(f"<td>{'True' if r.matched_success else 'False'}</td>")

            fh.write(f"<td class='rawcol'><pre class='raw'>{_html_escape(r.rdk_raw_text)}</pre></td>")
            fh.write(f"<td class='rawcol'><pre class='raw'>{_html_escape(r.control_raw_text)}</pre></td>")

            fh.write(f"<td>{_html_escape(r.pick_reason)}</td>")
            fh.write(f"<td>{_badge(r.confidence_label)}</td>")
            fh.write(f"<td class='num'>{'✔' if r.good_control else '—'}</td>")
            fh.write(f"<td class='num'>{r.flag_score}</td>")
            fh.write(f"<td class='num'>{r.flag_rmsd}</td>")
            fh.write(f"<td class='num'>{r.flag_identity}</td>")

            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(side)}</td>")
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(front)}</td>")
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(top)}</td>")
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(pside)}</td>")
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(pfront)}</td>")
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(ptop)}</td>")
            fh.write("</tr>")

        fh.write("""
</tbody>
</table>
</body>
</html>
""")
    return path


def write_details_xlsx(out_dir: Path, rows: Sequence[PairEval], include_identity: bool) -> Optional[Path]:
    """Try to write an .xlsx with embedded images. Uses xlsxwriter if available, else openpyxl. Returns path or None."""
    # First try xlsxwriter (no Pillow dependency required)
    try:
        import xlsxwriter  # type: ignore
        xlsx_path = out_dir / "benchmark_analysis_details.xlsx"
        wb = xlsxwriter.Workbook(str(xlsx_path))
        ws = wb.add_worksheet("details")
        # Column headers
        headers = [
            "pdb_id","pocket","control_id","control_het","rdk_id","rdk_suspected_fda","control_drug_guess",
            "ctrl_score","rdk_score","delta_kcal","rmsd_A","pick_reason",
            "confidence","good_control","score✔","center✔","rmsd✔","ident✔",
            "image_side","image_front","image_top","pair_side","pair_front","pair_top",
        ]
        for c,h in enumerate(headers):
            ws.write(0, c, h)
        # Set widths
        ws.set_column(0, 1, 10)
        ws.set_column(2, 3, 16)
        ws.set_column(4, 6, 20)
        ws.set_column(7, 12, 12)
        ws.set_column(13, 18, 11)
        ws.set_column(19, 24, 32)  # cover all 6 image columns
        # Data rows
        row = 1
        for r in rows:
            ws.write_row(row, 0, [
                r.pdb_id, r.pocket, r.control_id, r.control_het, r.rdk_id, r.rdk_name or "", r.control_display_name or "",
                _fmt(r.control_score), _fmt(r.rdk_score), _fmt(r.delta_kcal), _fmt(r.rmsd), r.pick_reason,
                r.confidence_label, r.good_control, r.flag_score, r.flag_center, r.flag_rmsd, r.flag_identity
            ])
            # Insert overlay images (scaled down)
            for offset, p in enumerate([r.png_side, r.png_front, r.png_top]):
                if p and Path(p).is_file():
                    try:
                        ws.insert_image(row, 19 + offset, str(p), {"x_scale":0.35, "y_scale":0.35})
                    except Exception:
                        pass
            # Insert pair-only images (scaled down)
            for offset, p in enumerate([r.png_pair_side, r.png_pair_front, r.png_pair_top]):
                if p and Path(p).is_file():
                    try:
                        ws.insert_image(row, 22 + offset, str(p), {"x_scale":0.35, "y_scale":0.35})
                    except Exception:
                        pass
            row += 1
        wb.close()
        return xlsx_path
    except Exception:
        pass

    # Fallback: openpyxl (requires Pillow for PNG sizing)
    try:
        from openpyxl import Workbook  # type: ignore
        from openpyxl.drawing.image import Image as XLImage  # type: ignore
        from openpyxl.utils import get_column_letter  # type: ignore
        xlsx_path = out_dir / "benchmark_analysis_details.xlsx"
        wb = Workbook()
        ws = wb.active
        ws.title = "details"
        headers = [
            "pdb_id","pocket","control_id","control_het","rdk_id","rdk_suspected_fda","control_drug_guess",
            "ctrl_score","rdk_score","delta_kcal","rmsd_A","pick_reason",
            "confidence","good_control","score✔","center✔","rmsd✔","ident✔",
            "image_side","image_front","image_top","pair_side","pair_front","pair_top",
        ]
        ws.append(headers)
        for r in rows:
            ws.append([
                r.pdb_id, r.pocket, r.control_id, r.control_het, r.rdk_id, r.rdk_name or "", r.control_display_name or "",
                _fmt(r.control_score), _fmt(r.rdk_score), _fmt(r.delta_kcal), _fmt(r.rmsd), r.pick_reason,
                r.confidence_label, r.good_control, r.flag_score, r.flag_center, r.flag_rmsd, r.flag_identity,
                "", "", "", "", "", ""  # placeholders for 6 images
            ])
            row_idx = ws.max_row
            # overlay
            for i, p in enumerate([r.png_side, r.png_front, r.png_top], start=20):
                try:
                    if p and Path(p).is_file():
                        img = XLImage(str(p))
                        ws.add_image(img, f"{get_column_letter(i)}{row_idx}")
                except Exception:
                    continue
            # pair-only
            for i, p in enumerate([r.png_pair_side, r.png_pair_front, r.png_pair_top], start=23):
                try:
                    if p and Path(p).is_file():
                        img = XLImage(str(p))
                        ws.add_image(img, f"{get_column_letter(i)}{row_idx}")
                except Exception:
                    continue
        wb.save(str(xlsx_path))
        return xlsx_path
    except Exception:
        return None



def write_details_visual_report(
    out_dir: Path,
    rows: Sequence[PairEval],
    include_identity: bool,
    *,
    html_success_only: Optional[bool] = None,
    since_epoch: Optional[float] = None
) -> Path:
    """Create a visual report attempting XLSX (with images) first, else HTML with inline images."""
    xlsx = write_details_xlsx(out_dir, rows, include_identity)
    if xlsx is not None:
        return xlsx
    return write_details_html(out_dir, rows, include_identity,
                              html_success_only=html_success_only, since_epoch=since_epoch)



# -----------------------------
# Main
# -----------------------------

def find_pdb_dirs(docked_root: Path, only_pdb: Optional[str]) -> List[Path]:
    pdb_dirs: List[Path] = []
    for child in docked_root.iterdir():
        if not child.is_dir():
            continue
        if only_pdb and child.name.upper() != only_pdb.upper():
            continue
        pdb_dirs.append(child)
    pdb_dirs.sort()
    return pdb_dirs

def run_analysis(
    docked_root: Path,
    mapping_csv: Optional[Path],
    only_pdb: Optional[str],
    score_tol: float,
    center_tol: float,
    rmsd_tol: float,
    include_identity: bool,
    out_dir: Optional[Path] = None,
    since_epoch: Optional[float] = None,
    html_success_only: Optional[bool] = None,
) -> Tuple[Optional[Path], Optional[Path]]:
    if not docked_root.is_dir():
        print(f"[analysis] DOCKED root not found: {docked_root}")
        return None, None

    mapping = MappingIndex(mapping_csv if mapping_csv and Path(mapping_csv).is_file() else None)

    all_rows: List[PairEval] = []
    for pdb_dir in find_pdb_dirs(docked_root, only_pdb):
        pdb_id = pdb_dir.name.upper()
        for pocket_dir in find_pocket_dirs(pdb_dir):
            # Skip pockets explicitly marked as having no RDKs this run
            sentinel = pocket_dir / ".rdk_skipped"
            if sentinel.exists():
                try:
                    has_rdk = any(RDK_PAT.match(ch.name) for ch in pocket_dir.iterdir() if ch.is_file())
                except Exception:
                    has_rdk = False
                if not has_rdk:
                    if DEBUG:
                        print(f"[debug] Skipping {pdb_id}/{pocket_dir.name} due to .rdk_skipped sentinel")
                    continue
            controls, rdks = load_controls_and_rdks(pocket_dir, mapping, since_epoch=since_epoch)
            if not controls or not rdks:
                if DEBUG:
                    print(f"[debug] Skipping {pdb_id}/{pocket_dir.name} (controls={len(controls)} rdks={len(rdks)})")
                continue
            rows = evaluate_pairs(
                pdb_id,
                pocket_dir.name,
                pocket_dir,
                controls,
                rdks,
                score_tol=score_tol,
                center_tol=center_tol,
                rmsd_tol=rmsd_tol,
                include_identity=include_identity,
            )
            all_rows.extend(rows)

    if not all_rows:
        print("[analysis] No matches found. Are the bench_pocketX_single folders populated?")
        return None, None

    out_root = out_dir if out_dir else (docked_root / "_analysis")
    details_path = write_details_csv(out_root, all_rows, include_identity)
    summary_path = write_summary_csv(out_root, all_rows, include_identity)
    visual_path = write_details_visual_report(out_root, all_rows, include_identity,
                                              html_success_only=html_success_only,
                                              since_epoch=since_epoch)

    max_points = 4 if include_identity else 3
    by_pdb_best: Dict[str, int] = {}
    for r in all_rows:
        by_pdb_best[r.pdb_id] = max(by_pdb_best.get(r.pdb_id, 0), r.total_points)

    passed_pdbs = sorted([p for p, v in by_pdb_best.items() if v == max_points])
    near_pdbs   = sorted([p for p, v in by_pdb_best.items() if 1 < v < max_points])
    fail0_pdbs  = sorted([p for p, v in by_pdb_best.items() if v == 0])

    global_best = max(by_pdb_best.values()) if by_pdb_best else 0
    best_pdbs   = sorted([p for p, v in by_pdb_best.items() if v == global_best])

    print(f"[analysis] Wrote details CSV: {details_path}")
    print(f"[analysis] Wrote summary CSV: {summary_path}")
    print(f"[analysis] Wrote visual details: {visual_path}")
    print(f"[analysis] Passed (=={max_points}): {len(passed_pdbs)} -> {', '.join(passed_pdbs) if passed_pdbs else '-'}")
    print(f"[analysis] Near (2..{max_points-1}): {len(near_pdbs)} -> {', '.join(near_pdbs) if near_pdbs else '-'}")
    print(f"[analysis] Fail (==0): {len(fail0_pdbs)} -> {', '.join(fail0_pdbs) if fail0_pdbs else '-'}")
    print(f"[analysis] PDBs with BEST score = {global_best}: {', '.join(best_pdbs) if best_pdbs else '-'}")

    return details_path, summary_path


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Auto-interpret benchmark results across proteins/pockets.")
    p.add_argument("--docked-root", default=DEFAULT_DOCKED_ROOT)
    p.add_argument("--mapping", default=DEFAULT_MAPPING_CSV)
    p.add_argument("--only-pdb", default=DEFAULT_ONLY_PDB)
    p.add_argument("--score-tol", type=float, default=DEFAULT_SCORE_TOL,
                   help="Score tolerance (kcal/mol) used for 'score✔'.")
    p.add_argument("--center-tol", type=float, default=DEFAULT_CENTER_TOL)
    p.add_argument("--rmsd-tol", type=float, default=DEFAULT_RMSD_TOL)
    p.add_argument("--no-identity", action="store_true", help="Exclude identity-match from scoring (max 3 points)")
    p.add_argument("--out", default=DEFAULT_OUTDIR)
    p.add_argument("--debug", action="store_true", help="Verbose per-pair diagnostics (or set BENCH_DEBUG=1)")
    p.add_argument("--analysis-score-tol", type=float, dest="score_tol",
                   help="Alias for --score-tol (kcal/mol).")
    p.add_argument("--html-expected-only", dest="html_expected_only", action="store_true", default=True,
                   help="Show only expected FDA RDKs in the HTML report (env ANALYSIS_HTML_EXPECTED_ONLY=1).")
    p.add_argument("--no-html-expected-only", dest="html_expected_only", action="store_false",
                   help="Disable expected-only filtering in HTML (env ANALYSIS_HTML_EXPECTED_ONLY=0).")
    p.add_argument("--since-epoch", type=float,
                   default=(float(DEFAULT_SINCE_EPOCH) if DEFAULT_SINCE_EPOCH else None),
                   help="Only ingest files with mtime >= this UNIX epoch (env BENCH_ANALYSIS_SINCE_EPOCH or BENCH_RUN_START_EPOCH).")
    p.add_argument("--html-success-only", dest="html_success_only", action="store_true", default=True,
                   help="HTML shows only success rows (score✔ & RMSD✔ & identity✔). Env ANALYSIS_HTML_SUCCESS_ONLY=1.")
    p.add_argument("--no-html-success-only", dest="html_success_only", action="store_false",
                   help="Disable success-only filtering for HTML (env ANALYSIS_HTML_SUCCESS_ONLY=0).")

    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_argparser().parse_args(argv)
    docked_root = Path(args.docked_root)
    mapping_csv = Path(args.mapping) if args.mapping else None
    out = Path(args.out) if args.out else None

    # honor CLI --debug
    global DEBUG
    DEBUG = DEBUG or bool(args.debug)
    _default_expected_only = str(_cfg.get("ANALYSIS_HTML_EXPECTED_ONLY", "true")).lower() in ("1", "true", "yes")

    run_analysis(
        docked_root=docked_root,
        mapping_csv=mapping_csv,
        only_pdb=(args.only_pdb or None),
        score_tol=float(args.score_tol),
        center_tol=float(args.center_tol),
        rmsd_tol=float(args.rmsd_tol),
        include_identity=not bool(args.no_identity),
        out_dir=out,
        since_epoch=(float(args.since_epoch) if args.since_epoch is not None else None),
        html_success_only=bool(args.html_success_only),
    )


if __name__ == "__main__":
    # Allow running directly in an IDE without supplying args
    try:
        main(None)
    except SystemExit:
        run_analysis(
            docked_root=Path(DEFAULT_DOCKED_ROOT),
            mapping_csv=Path(DEFAULT_MAPPING_CSV) if DEFAULT_MAPPING_CSV else None,
            only_pdb=DEFAULT_ONLY_PDB,
            score_tol=DEFAULT_SCORE_TOL,
            center_tol=DEFAULT_CENTER_TOL,
            rmsd_tol=DEFAULT_RMSD_TOL,
            include_identity=DEFAULT_INCLUDE_IDENTITY,
            out_dir=Path(DEFAULT_OUTDIR) if DEFAULT_OUTDIR else None,
            since_epoch=(float(DEFAULT_SINCE_EPOCH) if DEFAULT_SINCE_EPOCH else None),
            html_success_only=True,
        )
