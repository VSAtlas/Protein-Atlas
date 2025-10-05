# benchmark_mode.py
"""
Benchmark mode: single ultra-stage docking of likely co-crystal FDA ligands, per pocket.

This version:
- Resolves INPUT/OUTPUT/PREPPED from YAML config or CLI; no Windows paths.
- Resolves FDA mapping CSV from CLI -> YAML (FDA_MAPPING_CSV) -> alongside library -> OVERALL_DIR -> CWD
  -> fixed BRCF fallback: /stor/home/mpg2352/atlas/code/protein_automation/fda_mapping_from_pdbqt.csv
- Sources CHEMCOMP_ALIAS / EXCLUDE_* / PER_PDB_HINTS / HARD_FDA_CONTROL_BY_PDB from chemdb.chem_alias_db
- Keeps optional alias augmentation from prior benchmark details CSVs
- Uses  logic from main.py for docking/validation/switching
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import threading
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# Rendering helpers
from capture_pose import (
    _render_native_on_original_pdb,
    _render_three_views_with_pymol,
    _safe_open_csv_for_write,
    pick_control_and_nearest_rdk,
)

from concurrent.futures import ThreadPoolExecutor, as_completed

# Config
from input_and_export_functions import load_inputs, validate_config

# Library resolution helpers
from metabolite_resolver import (
    ensure_parent_drugs_for_controls,
    load_library_index,
    resolve_corresponding_name_for_rdk,
    resolve_corresponding_name_from_text,
)

# Fallback active-site detector
from protein_functions import detect_active_site

# Unified chem/alias config (YAML-backed)
from chemdb.chem_alias_db import (
    CHEMCOMP_ALIAS,
    EXCLUDE_HET_IDS,
    EXCLUDE_HET_NAME_KEYWORDS,
    HARD_FDA_CONTROL_BY_PDB,
    PER_PDB_HINTS,
)

# Limit simultaneous PyMOL renders (1 by default; override via env PYMOL_PARALLEL)
_RENDER_LOCK = threading.Semaphore(int(os.environ.get("PYMOL_PARALLEL", "1")))

# Import core pipeline pieces from main.py (reuses logic verbatim)
from main import (  # noqa: E402
    BudgetGuard,
    CenterSelector,
    GlobalCenterGuard,
    build_control_lookup,
    checkpoint_invalidate_from,
    checkpoint_mark_done,
    checkpoint_should_skip,
    detect_pocket,
    extract_ligands_to_nolig,
    final_pose_validation_and_screenshots,
    get_recenter_params,
    make_paths,
    make_protein_logger,
    prepare_receptor,
    record_le,
    record_score,
    run_one_stage,
    _count_heavy_atoms_from_pdbqt,
    _fingerprint_stage,
    early_recenter_decision,
    RetryManager,
)
from prep_ligands import prep_ligands_from_pdb



# ----- Deferred render queue (global) -----
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List
import threading
from contextlib import contextmanager

@contextmanager
def _acquire(sem: threading.Semaphore):
    sem.acquire()
    try:
        yield
    finally:
        sem.release()




# ---- Re-entrant patch for main.ThreadPoolExecutor ----
import threading
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor as _RealTPE

_GLOBAL_LIGAND_SEM = None  # unchanged
class _GuardedTPE(_RealTPE):
    def submit(self, fn, *args, **kwargs):
        sem = _GLOBAL_LIGAND_SEM
        if sem is None:
            return super().submit(fn, *args, **kwargs)
        def _wrapped(*a, **k):
            with _acquire(sem):
                return fn(*a, **k)
        return super().submit(_wrapped, *args, **kwargs)

# re-entrant state
_TP_PATCH_LOCK   = threading.Lock()
_TP_PATCH_COUNT  = 0
_TP_PREV_EXEC    = None

@contextmanager
def _patch_main_threadpool(sem):
    """
    Temporarily replace main.ThreadPoolExecutor with a guarded one.
    Re-entrant and thread-safe: multiple overlapping uses are OK.
    """
    import main as _main_mod
    global _GLOBAL_LIGAND_SEM, _TP_PATCH_COUNT, _TP_PREV_EXEC
    with _TP_PATCH_LOCK:
        _GLOBAL_LIGAND_SEM = sem
        if _TP_PATCH_COUNT == 0:
            _TP_PREV_EXEC = getattr(_main_mod, "ThreadPoolExecutor", None)
            setattr(_main_mod, "ThreadPoolExecutor", _GuardedTPE)
        _TP_PATCH_COUNT += 1
    try:
        yield
    finally:
        with _TP_PATCH_LOCK:
            _TP_PATCH_COUNT -= 1
            if _TP_PATCH_COUNT == 0:
                setattr(_main_mod, "ThreadPoolExecutor", _TP_PREV_EXEC)
                _GLOBAL_LIGAND_SEM = None


# ---- Re-entrant deferral for capture_pose PyMOL functions ----
import capture_pose as _cap
from concurrent.futures import ThreadPoolExecutor as _RealTPE2
from contextlib import contextmanager

_DEFERRED_PYMOL_CALLS = []
_ORIG_RENDER_THREE = getattr(_cap, "_render_three_views_with_pymol")
_ORIG_RENDER_NATIVE = getattr(_cap, "_render_native_on_original_pdb")

_PYMOL_DEFER_LOCK  = threading.Lock()
_PYMOL_DEFER_COUNT = 0

@contextmanager
def _defer_pymol_capture_calls(enable: bool = True):
    """
    Queue PyMOL renders instead of running them immediately.
    Re-entrant: only restore originals when the outermost context exits.
    """
    global _PYMOL_DEFER_COUNT
    if not enable:
        yield
        return
    with _PYMOL_DEFER_LOCK:
        if _PYMOL_DEFER_COUNT == 0:
            def _enqueue_three(*args, **kwargs):
                _DEFERRED_PYMOL_CALLS.append(("three", args, kwargs))
            def _enqueue_native(*args, **kwargs):
                _DEFERRED_PYMOL_CALLS.append(("native", args, kwargs))
            _cap._render_three_views_with_pymol = _enqueue_three
            _cap._render_native_on_original_pdb = _enqueue_native
        _PYMOL_DEFER_COUNT += 1
    try:
        yield
    finally:
        with _PYMOL_DEFER_LOCK:
            _PYMOL_DEFER_COUNT -= 1
            if _PYMOL_DEFER_COUNT == 0:
                _cap._render_three_views_with_pymol = _ORIG_RENDER_THREE
                _cap._render_native_on_original_pdb = _ORIG_RENDER_NATIVE

def run_deferred_captures(max_workers: int):
    """Replay queued mid-run captures in parallel."""
    if not _DEFERRED_PYMOL_CALLS:
        print("[render] no mid-run captures to replay.")
        return
    print(f"[render] replaying deferred mid-run captures: n={len(_DEFERRED_PYMOL_CALLS)} | workers={max_workers}")
    def _do(task):
        kind, args, kwargs = task
        if kind == "three":
            return _ORIG_RENDER_THREE(*args, **kwargs)
        elif kind == "native":
            return _ORIG_RENDER_NATIVE(*args, **kwargs)
    with _RealTPE2(max_workers=max_workers) as pool:
        list(pool.map(_do, _DEFERRED_PYMOL_CALLS))
    _DEFERRED_PYMOL_CALLS.clear()
    print("[render] deferred mid-run captures complete.")



_DEFERRED_RENDER_TASKS: List[Dict[str, str]] = []

def _enqueue_render_task(
    *,
    pdb_id: str,
    stage_name: str,
    root_project: str,
    cleaned_pdb_path: str,
    original_pdb_path: str,
    ctrl_pose_path: str | None,
    rdk_pose_path: str | None,
    exclude_resns: List[str],
) -> None:
    _DEFERRED_RENDER_TASKS.append(
        {
            "pdb_id": pdb_id,
            "stage_name": stage_name,
            "root_project": root_project,
            "cleaned_pdb_path": cleaned_pdb_path,
            "original_pdb_path": original_pdb_path,
            "ctrl_pose_path": ctrl_pose_path or "",
            "rdk_pose_path": rdk_pose_path or "",
            "exclude_resns": list(exclude_resns),
        }
    )

def _execute_render_task(task: Dict[str, str]) -> None:
    pdb_id = task["pdb_id"]
    stage_name = task["stage_name"]
    stage_dir_target = Path(task["root_project"]) / "docked" / pdb_id / stage_name
    stage_dir_target.mkdir(parents=True, exist_ok=True)

    # Always render native on original PDB
    _render_native_on_original_pdb(
        original_pdb=task["original_pdb_path"],
        outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-native__NATIVE"),
        exclude_resns=sorted(task["exclude_resns"]),
    )

    ctrl = task["ctrl_pose_path"] or ""
    rdk  = task["rdk_pose_path"] or ""
    cleaned = task["cleaned_pdb_path"]
    orig    = task["original_pdb_path"]

    if ctrl:
        _render_three_views_with_pymol(
            receptor_path=cleaned,
            ligand_paths_and_colors=[(ctrl, "control", "green")],
            outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL"),
        )

    if rdk:
        _render_three_views_with_pymol(
            receptor_path=cleaned,
            ligand_paths_and_colors=[(rdk, "rdk_closest", "magenta")],
            outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__RDKclosest"),
        )
        _render_three_views_with_pymol(
            receptor_path=orig,
            ligand_paths_and_colors=[(rdk, "rdk_closest", "magenta")],
            outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-rdk__RDKclosest"),
        )
        if ctrl:
            _render_three_views_with_pymol(
                receptor_path=cleaned,
                ligand_paths_and_colors=[(ctrl, "control", "green"), (rdk, "rdk_closest", "magenta")],
                outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL+RDKclosest"),
                label_top_n_res=5,
                label_cutoff=5.0,
            )

def run_deferred_renders(max_workers: int) -> None:
    if not _DEFERRED_RENDER_TASKS:
        print("[render] nothing to render (queue empty).")
        return
    print(f"[render] starting deferred renders: {len(_DEFERRED_RENDER_TASKS)} tasks | workers={max_workers}")
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futs = [pool.submit(_execute_render_task, t) for t in _DEFERRED_RENDER_TASKS]
        for _ in as_completed(futs):
            pass
    print("[render] all deferred renders complete.")



# =============================
# Small text/path helper utils
# =============================

def _norm(s: Optional[str]) -> str:
    """Lowercase, trim, collapse spaces, and strip punctuation-like separators."""
    if s is None:
        return ""
    s = str(s).strip().lower()
    s = s.replace("\u00A0", " ")
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"[-_/\\,;:\|\[\]\(\)\{\}\.\+\*'\"]+", " ", s)
    return s

def _tokenize(s: str) -> List[str]:
    """Tokenize to [a-z0-9+] words used by search/ranking."""
    if not s:
        return []
    toks = re.split(r"[^a-z0-9\+]+", s.lower())
    return [t for t in toks if t]

def _extract_rdk_id(text: str) -> Optional[str]:
    """Extract rdk_XXXXXX id from a filename-like string."""
    if not text:
        return None
    m = re.search(r"(rdk_\d{6,8})", Path(text).stem.lower())
    return m.group(1) if m else None

def _norm_text(s: Optional[str]) -> str:
    """Letters+digits only, for approximate comparisons."""
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())

def _looks_like_brand_or_generic(name: str) -> bool:
    n = (name or "").strip()
    if len(n) < 3 or len(n) > 64:
        return False
    if not re.search(r"[a-zA-Z]", n):
        return False
    noisy = len(re.findall(r"[0-9\[\]\(\)\/\.\-\,;:]", n))
    return (noisy / max(1, len(n))) < 0.40

def _clean_alias(name: str) -> Optional[str]:
    """Light cleanup and guardrails for text harvested from CSVs."""
    if not name:
        return None
    s = name.strip()
    if s.startswith("?"):
        s = s.lstrip("?\uFF1F").strip()
    s = s.replace("\u2019", "'").replace("\u00AE", "").replace("\u2122", "")
    s = re.sub(r"\s+", " ", s)
    return s.lower() if _looks_like_brand_or_generic(s) else None

def _expand_globs(paths: Sequence[str]) -> List[Path]:
    out: List[Path] = []
    for p in paths:
        if any(ch in p for ch in "*?[]"):
            out.extend(Path(x) for x in glob.glob(p))
        else:
            out.append(Path(p))
    return out

def _dedupe_str(seq: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for s in seq:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out

def _path_is_within(child: Path, parent: Path) -> bool:
    """Return True if `child` is within `parent` (no strict existence requirement)."""
    try:
        child = child.resolve(strict=False)
        parent = parent.resolve(strict=False)
        child.relative_to(parent)
        return True
    except Exception:
        return False

# =============================
# PDB parsing / hint extraction
# =============================

def parse_pdb_het_hints(pdb_path: Path) -> Tuple[List[str], List[str]]:
    """
    Parse HET / HETNAM / HETATM records from a PDB to extract plausible ligand HET codes
    and human-readable names, applying EXCLUDE_* filters.
    """
    pdb_path = Path(pdb_path)
    het_ids: set[str] = set()
    hetnam_map: Dict[str, List[Tuple[int, str]]] = defaultdict(list)

    if not pdb_path.is_file():
        return [], []

    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for raw in fh:
                rec = raw[:6].strip().upper()

                if rec == "HET":
                    het = (raw[7:10].strip() if len(raw) >= 10 else "").upper()
                    if not het:
                        toks = raw.split()
                        if len(toks) >= 2:
                            het = toks[1].upper()
                    if het:
                        het_ids.add(het)

                elif rec == "HETNAM":
                    toks = raw.split()
                    cont = 1
                    het, name = "", ""
                    if len(toks) >= 3 and toks[1].isdigit():
                        cont = int(toks[1]); het = toks[2].upper(); name = " ".join(toks[3:])
                    elif len(toks) >= 2:
                        het = toks[1].upper(); name = " ".join(toks[2:])
                    if het and name:
                        hetnam_map[het].append((cont, name.strip()))

                elif rec == "HETATM":
                    resn = raw[17:20].strip().upper()
                    if resn:
                        het_ids.add(resn)
    except Exception:
        return [], []

    STANDARD_RES = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET",
        "PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL"
    }

    het_ids = {
        h for h in het_ids
        if h not in EXCLUDE_HET_IDS and h not in STANDARD_RES and 2 <= len(h) <= 5
    }

    het_names: List[str] = []
    for het in sorted(het_ids):
        if het in hetnam_map:
            parts = [t for _, t in sorted(hetnam_map[het], key=lambda x: x[0])]
            nm = " ".join(parts).strip()
            if not nm:
                continue
            nm_up = nm.upper()
            if any(kw in nm_up for kw in EXCLUDE_HET_NAME_KEYWORDS):
                continue
            cleaned = _clean_alias(nm)
            het_names.append(cleaned if cleaned else nm)

    # dedupe preserving order
    seen = set()
    dedup_names: List[str] = []
    for n in het_names:
        nn = _norm(n)
        if nn and nn not in seen:
            seen.add(nn)
            dedup_names.append(n)

    return sorted(het_ids), dedup_names

# =============================
# Alias augmentation (details CSV)
# =============================

def load_aliases_from_details_csv(paths: Sequence[str]) -> Dict[str, List[str]]:
    """
    Harvest (HET -> alias) from benchmark_analysis_details.csv files.
    Returns aliases lowercased and deduped per HET code.
    """
    exp: Dict[str, set] = {}
    for csv_path in _expand_globs(paths):
        if not csv_path.is_file():
            continue
        try:
            with open(csv_path, "r", encoding="utf-8", errors="ignore", newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    het = (row.get("control_het") or "").strip().upper()
                    nm = _clean_alias(row.get("rdk_name") or "")
                    if het and nm:
                        exp.setdefault(het, set()).add(nm)
        except Exception as e:
            print(f"[alias] WARN: failed to parse {csv_path}: {e}")
    return {k: sorted(v) for k, v in exp.items()}

def merge_aliases_into_chemcomp(base: Dict[str, List[str]], extra: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """Merge extra aliases into an existing chemcomp map. Output values are lowercased & deduped."""
    out: Dict[str, List[str]] = {k: sorted(set(vv.lower() for vv in vs)) for k, vs in base.items()}
    for het, names in (extra or {}).items():
        merged = set(out.get(het, []))
        for nm in names:
            nm_norm = (nm or "").strip().lower()
            if nm_norm:
                merged.add(nm_norm)
        out[het] = sorted(merged)
    return out

# =============================
# Mapping index (FDA mapping CSV)
# =============================

@dataclass
class MappingRow:
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

    def all_name_fields(self) -> List[Tuple[str, str]]:
        return [
            ("display_name", self.display_name),
            ("generic_name", self.generic_name),
            ("brand_names", self.brand_names),
            ("pubchem_name", self.pubchem_name),
            ("pubchem_record_title", self.pubchem_record_title),
            ("pubchem_iupac_name", self.pubchem_iupac_name),
            ("pubchem_synonyms", self.pubchem_synonyms),
            ("rxnorm_generic_name", self.rxnorm_generic_name),
            ("rxnorm_brand_names", self.rxnorm_brand_names),
            ("drugcentral_generic_name", self.drugcentral_generic_name),
            ("drugcentral_brand_names", self.drugcentral_brand_names),
            ("remark_name", self.remark_name),
            ("sdf_title", self.sdf_title),
        ]

class MappingIndex:
    def __init__(self, csv_path: Path):
        import pandas as pd

        self.csv_path = Path(csv_path)
        if not self.csv_path.exists():
            raise FileNotFoundError(f"Mapping CSV not found: {self.csv_path}")
        self.df = pd.read_csv(self.csv_path)
        if "path" not in self.df.columns:
            raise ValueError("Mapping CSV must include a 'path' column to .pdbqt files")

        self.rows: List[MappingRow] = []
        for _, r in self.df.iterrows():
            self.rows.append(
                MappingRow(
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

    def rows_by_rdk_id(self, rdk_id: str) -> List[MappingRow]:
        rid = (rdk_id or "").strip().lower()
        if not rid:
            return []
        out = []
        for row in self.rows:
            stem = Path(row.path).stem.lower()
            if rid in stem:
                out.append(row)
        return out

    def search(
        self,
        hints: Sequence[str],
        inchikey: Optional[str] = None,
        max_results: int = 6,
    ) -> List[Tuple[MappingRow, int, str]]:
        hints_norm = [_norm(h) for h in hints if h]
        hint_tokens = set(t for h in hints_norm for t in _tokenize(h))
        raw_het_codes = {h.strip().upper() for h in hints if h and 2 <= len(h.strip()) <= 5}
        out: List[Tuple[MappingRow, int, str]] = []

        for row in self.rows:
            if not row.path or not Path(row.path).exists():
                continue

            best = 0
            why = ""
            p = Path(row.path)
            stem_upper = p.stem.upper()
            base_upper = p.name.upper()
            parents_upper = " ".join([pp.name.upper() for pp in p.parents])

            # Priority 1: exact InChIKey
            if inchikey and row.inchikey and row.inchikey.strip().upper() == inchikey.strip().upper():
                best, why = 100, "inchikey_exact"

            # Priority 2: HET codes visible in path tokens
            if best < 100 and raw_het_codes:
                for code in raw_het_codes:
                    if re.search(rf"\b{re.escape(code)}\b", stem_upper) or re.search(
                        rf"\b{re.escape(code)}\b", parents_upper
                    ):
                        sc = 93
                        rs = f"path_token:{code}"
                    elif code in base_upper:
                        sc = 88
                        rs = f"path_substr:{code}"
                    else:
                        sc = 0
                        rs = ""
                    if sc > best:
                        best, why = sc, rs

            # Priority 3: name-field matches
            if best < 100 and hints_norm:
                for field, value in row.all_name_fields():
                    v = _norm(value)
                    if not v:
                        continue
                    if v in hints_norm:
                        sc = 95
                        rs = f"{field}_exact"
                    elif any(h in v for h in hints_norm if len(h) >= 3):
                        sc = 85
                        rs = f"{field}_substr"
                    else:
                        vtok = set(_tokenize(v))
                        overlap = len(vtok & hint_tokens)
                        sc = 70 if overlap >= 2 else (65 if overlap == 1 else 0)
                        rs = f"{field}_tokens:{overlap}"
                    if sc > best:
                        best, why = sc, rs

            if best > 0:
                out.append((row, best, why))

        out.sort(key=lambda t: t[1], reverse=True)
        return out[:max_results]

# =============================
# Selection / orchestration
# =============================

def _resolve_name_for_path_or_text(p_or_text: str, fda_index) -> str:
    """Prefer RDK-based resolution; fallback to filename text matching."""
    rdk = _extract_rdk_id(p_or_text or "")
    if rdk:
        name = resolve_corresponding_name_for_rdk(rdk, fda_index) or ""
        if name:
            return name
    stem = Path(p_or_text).stem if os.path.exists(p_or_text) else str(p_or_text)
    name = resolve_corresponding_name_from_text(stem, fda_index) or ""
    return name

def _remap_to_prepped(p: Path, prepped_dir: Path) -> Path:
    """If p isn't within prepped_dir, try to map it to a file that *is* there (by basename or rdk_* id)."""
    if p.is_file() and p.resolve().parent == prepped_dir.resolve():
        return p
    # 1) Same basename inside prepped
    q = prepped_dir / p.name
    if q.is_file():
        return q
    # 2) rdk_* id match inside prepped
    rid = _extract_rdk_id(p.name)
    if rid:
        hits = list(prepped_dir.glob(f"*{rid}*.pdbqt"))
        if hits:
            return hits[0]
    return p  # unchanged if we couldn't remap

def select_candidates_for_protein(
    mapping: MappingIndex,
    hints: Sequence[str],
    prepped_dir: Path,
    extra_parent_ids: Sequence[str],
    max_candidates: int,
) -> List[Tuple[MappingRow, int, str]]:
    candidates = mapping.search(hints, inchikey=None, max_results=max(6, max_candidates * 4))

    parent_rows: List[Tuple[MappingRow, int, str]] = []
    prepped_dir_resolved = prepped_dir.resolve()
    cand_rows: List[Tuple[MappingRow, int, str]] = []

    _d_exist_fail = 0
    _d_scope_fail = 0

    for row, sc, why in (candidates + parent_rows):
        p = Path(row.path)

        # NEW: remap to the current prepped_dir if needed
        p = _remap_to_prepped(p, prepped_dir_resolved)

        # Keep only files that exist *and* live inside prepped_dir
        if not p.is_file():
            _d_exist_fail += 1
            continue
        if p.resolve().parent != prepped_dir_resolved and not _path_is_within(p, prepped_dir_resolved):
            _d_scope_fail += 1
            continue

        cand_rows.append((MappingRow(path=str(p),  # keep the remapped, absolute path
                                     display_name=row.display_name,
                                     generic_name=row.generic_name,
                                     brand_names=row.brand_names,
                                     pubchem_name=row.pubchem_name,
                                     pubchem_record_title=row.pubchem_record_title,
                                     pubchem_iupac_name=row.pubchem_iupac_name,
                                     pubchem_synonyms=row.pubchem_synonyms,
                                     rxnorm_generic_name=row.rxnorm_generic_name,
                                     rxnorm_brand_names=row.rxnorm_brand_names,
                                     drugcentral_generic_name=row.drugcentral_generic_name,
                                     drugcentral_brand_names=row.drugcentral_brand_names,
                                     remark_name=row.remark_name,
                                     sdf_title=row.sdf_title,
                                     inchikey=row.inchikey),
                           sc, why))

    if _d_exist_fail or _d_scope_fail:
        print(f"[DEBUG] select_candidates: dropped not_exist={_d_exist_fail}, out_of_scope={_d_scope_fail}, kept={len(cand_rows)}")

    for pid in extra_parent_ids or []:
        rid = _extract_rdk_id(pid or "")
        if rid:
            parent_rows.extend((r, 99, "parent_from_metabolite") for r in mapping.rows_by_rdk_id(rid))
        else:
            for (r, sc, why) in mapping.search([pid], inchikey=None, max_results=2):
                parent_rows.append((r, max(sc, 92), f"{why}|parent_from_metabolite"))

    # best per path
    best_by_path: Dict[str, Tuple[MappingRow, int, str]] = {}
    for row, sc, why in (candidates + parent_rows):
        p_abs = str(Path(row.path).resolve())
        if (p_abs not in best_by_path) or (sc > best_by_path[p_abs][1]):
            best_by_path[p_abs] = (row, sc, why)

    prepped_dir_resolved = prepped_dir.resolve()
    cand_rows: List[Tuple[MappingRow, int, str]] = []
    for row, sc, why in best_by_path.values():
        p = Path(row.path)
        if p.is_file() and _path_is_within(p, prepped_dir_resolved):
            cand_rows.append((row, sc, why))

    cand_rows.sort(key=lambda t: t[1], reverse=True)
    return cand_rows[:max_candidates]

def split_controls_and_whitelist(
    whitelist_paths: Sequence[str],
    prepped_control_pdbqts: Sequence[str],
    control_stems_lower: set,
    heavy_atom_counts: Dict[str, int],
    cfg: Dict,
    logger,
) -> Tuple[List[str], List[str]]:
    """Divide ligands into control vs non-control (with caps and HA minimum)."""
    ctrl_blacklist = {t.strip().upper() for t in str(cfg.get("CONTROL_BLACKLIST", "")).split(",") if t.strip()}
    min_ha = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))

    def _is_control_path(p: str) -> bool:
        stem = Path(p).stem.split("_stage")[0]
        if stem.upper() in ctrl_blacklist:
            return False
        if stem.lower() not in control_stems_lower:
            return False
        ha = heavy_atom_counts.get(p)
        return (ha is None) or (ha >= min_ha)

    ctrls = [p for p in prepped_control_pdbqts if _is_control_path(p)]
    non_ctrls = [p for p in whitelist_paths if not _is_control_path(p)]
    ctrls = _dedupe_str(ctrls)
    non_ctrls = _dedupe_str(non_ctrls)

    max_ctrls = int(cfg.get("BENCH_MAX_CONTROLS", 8))
    if len(ctrls) > max_ctrls:
        if logger:
            logger.warning(f"[SANITY] controls={len(ctrls)} looks high; capping to {max_ctrls}. Offenders will be skipped.")
            for p in ctrls[max_ctrls:]:
                logger.warning(f"[SANITY-OFFENDER] {p}")
        ctrls = ctrls[:max_ctrls]
    return ctrls, non_ctrls

def _merge_per_pdb_hints(pdb_id: str, hints: List[str]) -> List[str]:
    """Union (deduped) of existing hints with per-PDB and hard-control hints."""
    out = list(hints or [])
    add: List[str] = []
    pid = pdb_id.upper()
    if pid in PER_PDB_HINTS:
        add.extend(PER_PDB_HINTS[pid])
    if pid in HARD_FDA_CONTROL_BY_PDB:
        add.extend(HARD_FDA_CONTROL_BY_PDB[pid])
    seen = set()
    merged: List[str] = []
    for h in (out + add):
        hn = _norm(h)
        if hn and hn not in seen:
            seen.add(hn)
            merged.append(h)
    return merged

def _promote_forced_controls_to_ctrls(
    pdb_id: str,
    mapping: MappingIndex,
    prepped_dir: Path,
    names: Sequence[str],
    logger=None,
) -> Tuple[List[str], List[str]]:
    """Promote explicitly requested controls (from HARD_FDA_CONTROL_BY_PDB) to the control pool."""
    if not names:
        return [], []
    rows = mapping.search(names, inchikey=None, max_results=max(8, len(names) * 3))
    prepped_dir_resolved = prepped_dir.resolve()
    promoted_paths: List[str] = []
    promoted_stems_lower: List[str] = []
    for row, sc, why in rows:
        p = Path(row.path)
        if p.is_file() and p.resolve().parent == prepped_dir_resolved:
            promoted_paths.append(str(p.resolve()))
            stem0 = p.stem.split("_stage")[0].lower()
            promoted_stems_lower.append(stem0)
            if logger:
                logger.info(f"[forced-control] {pdb_id}: promoting '{stem0}' as control (reason={why}, score={sc})")
    promoted_paths = _dedupe_str(promoted_paths)
    promoted_stems_lower = _dedupe_str(promoted_stems_lower)
    if logger and not promoted_paths:
        logger.warning(
            f"[forced-control] {pdb_id}: requested {list(names)} but none were found in prepped library '{prepped_dir}'."
        )
    return promoted_paths, promoted_stems_lower

def collect_prepped_controls_for_protein(paths, control_stems: Sequence[str], logger=None) -> Tuple[List[str], set]:
    """Find control pdbqts colocated in the protein’s prepped ligands directory."""
    control_stems_lower = {s.lower() for s in control_stems if s}
    prepped_control_pdbqts: List[str] = []
    prepped_dir_resolved = paths.prepped_ligands_dir.resolve()
    if paths.prepped_ligands_dir.exists():
        for p in paths.prepped_ligands_dir.glob("*.pdbqt"):
            stem0 = p.stem.split("_stage")[0].lower()
            if stem0 in control_stems_lower and p.parent.resolve() == prepped_dir_resolved:
                prepped_control_pdbqts.append(str(p.resolve()))
    prepped_control_pdbqts = sorted(set(prepped_control_pdbqts))
    if logger:
        leaks = [pp for pp in prepped_control_pdbqts if Path(pp).resolve().parent != prepped_dir_resolved]
        if leaks:
            logger.warning(f"[CONTROL-LEAK] Found {len(leaks)} controls outside {prepped_dir_resolved}")
            for pp in leaks[:10]:
                logger.warning(f"  leak -> {pp}")
    return prepped_control_pdbqts, control_stems_lower

# =============================
# Benchmark driver
# =============================

def run_benchmark_for_protein(
    cfg: Dict,
    mapping: MappingIndex,
    pdb_file: str,
    prepped_dir: Path,
    out_root: Path,
    exhaustiveness: int,
    num_modes: int,
    max_candidates: int,
    manual_hints: Optional[List[str]] = None,
    fda_index=None,
) -> None:
    base_id = os.path.splitext(pdb_file)[0]
    pdb_id = base_id.replace("_cleaned", "").upper()

    # Logger + ASCII filter for mixed terminals
    logger = make_protein_logger(cfg["DOCKED_DIR"], pdb_id, cfg)
    import logging, sys

    def _sanitize_msg(s: str) -> str:
        return (
            s.replace("≤", "<=")
             .replace("≥", ">=")
             .replace("Å", " Angstrom")
             .replace("µ", "u")
             .replace("°", " deg")
        )

    class _AsciiFilter(logging.Filter):
        def filter(self, record):
            if isinstance(record.msg, str):
                record.msg = _sanitize_msg(record.msg)
            return True

    for h in logger.handlers:
        h.addFilter(_AsciiFilter())
    for stream_name in ("stdout", "stderr"):
        try:
            getattr(sys, stream_name).reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    paths = make_paths(cfg, base_id, pdb_file)
    # --- DEBUG visibility ---
    try:
        logger.info("[DEBUG] prepped_dir=%s", str(paths.prepped_ligands_dir.resolve()))
        # Mapping stats
        try:
            n_rows = len(mapping.rows)
        except Exception:
            n_rows = -1
        logger.info("[DEBUG] mapping_rows=%d", n_rows)
        if n_rows > 0:
            # sample a few paths and whether they exist & are within prepped_dir
            from pathlib import Path as _P
            _pd = paths.prepped_ligands_dir.resolve()
            sample = mapping.rows[:12]  # first dozen
            bad_exist = 0
            bad_scope = 0
            for r in sample:
                p = _P(r.path)
                if not p.exists():
                    bad_exist += 1
                else:
                    try:
                        if p.resolve().parent != _pd and _pd not in p.resolve().parents:
                            bad_scope += 1
                    except Exception:
                        bad_scope += 1
            logger.info("[DEBUG] mapping sample: not_exist=%d, out_of_scope=%d (first 12)", bad_exist, bad_scope)
    except Exception as _e:
        logger.warning("[DEBUG] mapping-precheck failed: %s", _e)

    # 1) Extract controls & build nolig
    _, control_stems = extract_ligands_to_nolig(paths, logger)
    # Always prep the extracted crystal ligands with the same sanitizer/ADT flags (-A hydrogens)
    prep_ligands_from_pdb(
        ligand_output_dir=paths.ligand_output_dir,
        ligands_mol2_dir=paths.ligands_mol2_dir,
        prepped_ligands_dir=paths.prepped_ligands_dir,
    )

    # 1b) metabolites -> parents
    extra_parent_ids: List[str] = []
    try:
        extra_parent_ids = ensure_parent_drugs_for_controls(
            pdb_code=pdb_id,
            ligands_raw_dir=str(paths.ligand_output_dir),
            fda_index=fda_index,
            max_additions=3,
        ) or []
        if extra_parent_ids:
            logger.info(f"[metabolite→parent] Will also consider parents: {', '.join(extra_parent_ids)}")
    except Exception as e:
        logger.warning(f"[metabolite→parent] resolver failed: {e}")

    # 2) Prepare receptor
    cleaned_pdb, receptor_pdbqt = prepare_receptor(cfg, paths, logger)
    if not cleaned_pdb or not receptor_pdbqt:
        logger.warning("[benchmark] receptor prep failed; skipping protein")
        return

    # 3) Detect pocket (main logic)
    center, detected_box, src = detect_pocket(cleaned_pdb, paths.ligand_output_dir, logger)
    if center:
        pockets = [("pocket1", center)]
        pocket_box_size = tuple(min(28.0, float(s)) for s in (detected_box or (24.0, 24.0, 24.0)))
    else:
        logger.error("[benchmark] no pocket could be detected; skipping protein")
        return
    if not pockets:
        c, _, _ = detect_active_site(cleaned_pdb)
        if c:
            pockets = [("pocket1", c)]
        else:
            logger.error("[benchmark] no pocket could be detected; skipping protein")
            return

    # 4) Build hints
    het_ids, het_names = parse_pdb_het_hints(Path(paths.pdb_path))
    hints: List[str] = []
    hints.extend(het_names)
    hints.extend(het_ids)
    for het in het_ids:
        hints.extend(CHEMCOMP_ALIAS.get(het.upper(), []))
    if manual_hints:
        hints.extend(manual_hints)
    hints.extend(extra_parent_ids)

    # dedupe normalized
    seen = set()
    final_hints: List[str] = []
    for h in hints:
        hn = _norm(h)
        if hn and hn not in seen:
            final_hints.append(h)
            seen.add(hn)
    final_hints = _merge_per_pdb_hints(pdb_id, final_hints)

    logger.info("[DEBUG] hint_count=%d | examples=%s", len(final_hints), ", ".join(map(str, final_hints[:8])))
    # 5) Candidate selection
    cand_rows = select_candidates_for_protein(
        mapping=mapping,
        hints=final_hints,
        prepped_dir=prepped_dir,
        extra_parent_ids=extra_parent_ids,
        max_candidates=max_candidates,
    )
    cand_rows.sort(key=lambda t: t[1], reverse=True)

    # Audit
    out_dir = Path(out_root) / pdb_id
    fh, cand_csv_path = _safe_open_csv_for_write(out_dir / f"benchmark_candidates_{pdb_id}.csv")
    with fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "rank",
                "score",
                "reason",
                "path",
                "display_name",
                "generic_name",
                "brand_names",
                "inchikey",
                "corresponding_name",
            ]
        )
        for i, (row, sc, why) in enumerate(cand_rows, start=1):
            corr = _resolve_name_for_path_or_text(row.path, fda_index) or row.generic_name or row.display_name
            w.writerow([i, sc, why, row.path, row.display_name, row.generic_name, row.brand_names, row.inchikey, corr])
    if str(cand_csv_path) != str(out_dir / f"benchmark_candidates_{pdb_id}.csv"):
        logger.info(f"[benchmark] candidates CSV was locked; wrote to {cand_csv_path.name}")

    whitelist = [row.path for (row, _, _) in cand_rows]
    if not whitelist:
        logger.warning("[benchmark] no candidate ligands matched mapping; skipping protein")
        return

    # Heavy atom counts
    heavy_atom_counts: Dict[str, int] = {}
    for p in whitelist:
        try:
            heavy_atom_counts[p] = int(_count_heavy_atoms_from_pdbqt(Path(p)))
        except Exception:
            heavy_atom_counts[p] = 0

    # Control lookup
    control_lookup = build_control_lookup(paths)

    # Controls present in THIS protein’s prepped dir
    prepped_control_pdbqts, control_stems_lower = collect_prepped_controls_for_protein(
        paths=paths, control_stems=control_stems, logger=logger
    )
    if pdb_id.upper() in HARD_FDA_CONTROL_BY_PDB:
        forced_names = HARD_FDA_CONTROL_BY_PDB[pdb_id.upper()]
        forced_paths, forced_stems_lower = _promote_forced_controls_to_ctrls(
            pdb_id=pdb_id,
            mapping=mapping,
            prepped_dir=prepped_dir,
            names=forced_names,
            logger=logger,
        )
        if forced_paths:
            prepped_control_pdbqts = _dedupe_str(prepped_control_pdbqts + forced_paths)
            control_stems_lower = set(control_stems_lower) | set(forced_stems_lower)

    for p in prepped_control_pdbqts:
        try:
            heavy_atom_counts[p] = int(_count_heavy_atoms_from_pdbqt(Path(p)))
        except Exception:
            heavy_atom_counts[p] = heavy_atom_counts.get(p, 0)

    # -------- single-stage per pocket --------
    pocket_strength_rows: List[List[str]] = []

    for pocket_name, center in pockets[:2]:
        stage = {
            "name": f"bench_{pocket_name}_single",
            "exhaustiveness": int(exhaustiveness),
            "num_modes": int(num_modes),
            "energy_range": int(cfg.get("ENERGY_RANGE", 9)),
        }

        params = get_recenter_params(cfg)
        guard = GlobalCenterGuard(max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2)))
        selector = CenterSelector(
            cfg,
            logger,
            control_stems=set(control_stems),
            heavy_atom_counts=heavy_atom_counts,
            initial_center=center,
        )
        box_size: Tuple[float, float, float] = tuple(min(28.0, float(s)) for s in (24.0, 24.0, 24.0))
        # prefer detect_pocket size
        box_size = tuple(min(28.0, float(s)) for s in (box_size if not isinstance(center, tuple) else box_size))
        if 'detected_box' in locals():
            box_size = tuple(min(28.0, float(s)) for s in (detected_box or box_size))

        # Ligand pool
        ctrls, non_ctrls = split_controls_and_whitelist(
            whitelist_paths=whitelist,
            prepped_control_pdbqts=prepped_control_pdbqts,
            control_stems_lower=control_stems_lower,
            heavy_atom_counts=heavy_atom_counts,
            cfg=cfg,
            logger=logger,
        )

        max_total = int(cfg.get("BENCH_MAX_TOTAL_LIGANDS", 5))
        if max_total > 0:
            allowed_non_ctrls = max(0, max_total - len(ctrls))
            if allowed_non_ctrls < len(non_ctrls):
                if logger:
                    logger.info(
                        f"[CAP] {pdb_id}: limiting non-controls from {len(non_ctrls)} to {allowed_non_ctrls} "
                        f"(controls={len(ctrls)}, cap={max_total})"
                    )
                non_ctrls = non_ctrls[:allowed_non_ctrls]

        ligands = _dedupe_str(ctrls + non_ctrls)
        stage1_original = ligands[:]
        logger.info(f"[DEBUG] pools: ctrls={len(ctrls)}, whitelist={len(non_ctrls)}, total={len(ligands)}")

        # Per-ligand budget guards (persist across restarts/waves)
        per_ligand_seconds = float(cfg.get("BENCH_MAX_SECONDS", 9000.0))
        guards_all: Dict[str, BudgetGuard] = {lig: BudgetGuard(per_ligand_seconds) for lig in ligands}

        score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
        validated_ligands_last: List[str] = []
        recenter_attempts = 0

        # Optional checkpoint skip
        if bool(cfg.get("CHECKPOINT_ENABLE", False)):
            fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
            if checkpoint_should_skip(cfg, pdb_id, stage["name"], fp):
                logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                continue

        restarting = True
        while restarting:
            restarting = False
            guard.reset_stage()

            if not ligands:
                logger.warning(f"No ligands to dock at {stage['name']}; skipping pocket.")
                break

            logger.info(f"[benchmark] {pdb_id} | {pocket_name} | starting {stage['name']} with {len(ligands)} ligands")

            if ctrls and non_ctrls:
                logger.info(f"Single-stage two-wave: {len(ctrls)} controls first, then {len(non_ctrls)} whitelist.")

                # Wave A — controls
                with _patch_main_threadpool(cfg.get("GLOBAL_LIGAND_SEM")):
                    s1, v1, d1, rd1, inv1 = run_one_stage(
                        cfg=cfg,
                        pdb_id=pdb_id,
                        receptor_pdbqt=receptor_pdbqt,
                        center=center,
                        box_size=box_size,
                        stage=stage,
                        ligands=ctrls,
                        logger=logger,
                        retry_mgr=RetryManager(),
                        control_lookup=control_lookup,
                        budget_guards={lig: guards_all[lig] for lig in ctrls},
                    )

                # CenterSelector after controls
                try:
                    dec = selector.consider_switch(stage["name"], s1, v1, rd1, receptor_pdbqt, center, guard)
                    if dec.promoted and dec.new_center is not None:
                        old = center
                        center = dec.new_center
                        guard.mark_switch()
                        if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                            checkpoint_invalidate_from(cfg, pdb_id, [stage], start_index=0)
                        logger.info(f"[CENTER] switched (controls wave): {old} -> {center} ({dec.reason}) [global switch]")
                except Exception as e:
                    logger.warning(f"CenterSelector (controls) failed: {e}")

                # Optional control lock
                lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", -6.0))
                lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
                lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))
                qualified_controls = []
                for lig in v1:
                    stem = Path(lig).stem.split("_stage")[0].lower()
                    if stem not in control_stems_lower:
                        continue
                    sc = s1.get(lig)
                    if sc is None or not isinstance(sc, (int, float)):
                        continue
                    if sc > lock_score_max:
                        continue
                    pose_path = rd1.get(lig)
                    if not pose_path:
                        continue
                    c = CenterSelector._pdbqt_centroid(pose_path)
                    if c is None:
                        continue
                    import numpy as np
                    if np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                        qualified_controls.append(lig)
                if len(qualified_controls) >= lock_min_hits and not guard.locked:
                    guard.lock()
                    logger.info(
                        f"[CONTROL-LOCK] n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Angstrom"
                    )

                # Wave B — whitelist
                with _patch_main_threadpool(cfg.get("GLOBAL_LIGAND_SEM")):
                    s2, v2, d2, rd2, inv2 = run_one_stage(
                        cfg=cfg,
                        pdb_id=pdb_id,
                        receptor_pdbqt=receptor_pdbqt,
                        center=center,
                        box_size=box_size,
                        stage=stage,
                        ligands=non_ctrls,
                        logger=logger,
                        retry_mgr=RetryManager(),
                        control_lookup=control_lookup,
                        budget_guards={lig: guards_all[lig] for lig in non_ctrls},
                    )

                # Merge
                scores, validated, distances = ({**s1, **s2}, v1 + v2, d1 + d2)
                raw_docked = {**rd1, **rd2}
                invalids = {**inv1, **inv2}
            else:
                # Single wave
                with _patch_main_threadpool(cfg.get("GLOBAL_LIGAND_SEM")):
                    scores, validated, distances, raw_docked, invalids = run_one_stage(
                        cfg=cfg,
                        pdb_id=pdb_id,
                        receptor_pdbqt=receptor_pdbqt,
                        center=center,
                        box_size=box_size,
                        stage=stage,
                        ligands=ligands,
                        logger=logger,
                        retry_mgr=RetryManager(),
                        control_lookup=control_lookup,
                        budget_guards={lig: guards_all[lig] for lig in ligands},
                    )

            validated_ligands_last = validated

            # Record
            for lig, sc in scores.items():
                record_score(score_history, stage["name"], lig, sc, True)
                record_le(score_history, stage["name"], lig, sc, heavy_atom_counts)
            for lig, (sc, reason) in invalids.items():
                record_score(score_history, stage["name"], lig, sc, False, reason=str(reason))
                record_le(score_history, stage["name"], lig, sc, heavy_atom_counts)

            # CenterSelector after full run
            try:
                decision = selector.consider_switch(
                    stage["name"], scores, validated, raw_docked, receptor_pdbqt, center, guard
                )
                if decision.promoted and decision.new_center is not None:
                    old = center
                    center = decision.new_center
                    guard.mark_switch()
                    if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                        checkpoint_invalidate_from(cfg, pdb_id, [stage], start_index=0)
                    logger.info(f"[CENTER] {pdb_id} {pocket_name}: {old} -> {center} ({decision.reason}) [global switch]")
            except Exception as e:
                logger.warning(f"CenterSelector failed gracefully: {e}")

            # Early recenter (single-stage flavor)
            restart, center, box_size, redo_ligands, recenter_attempts = early_recenter_decision(
                0,
                scores,
                distances,
                box_size,
                center,
                stage1_original,
                recenter_attempts,
                params,
                cfg,
                pdb_id,
                receptor_pdbqt,
                logger,
                raw_docked,
                guard,
                control_anchor_hit=False,
            )
            if restart:
                ligands = redo_ligands
                if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                    checkpoint_invalidate_from(cfg, pdb_id, [stage], start_index=0)
                restarting = True
                continue

            # Adaptive shrink
            try:
                if cfg.get("ADAPTIVE_SHRINK_ENABLE", True) and validated:
                    import numpy as np
                    med = float(
                        np.median([d for d in distances if isinstance(d, (int, float))])
                    ) if distances else None
                    if (med is not None) and (med < float(cfg.get("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0))):
                        dec = float(cfg.get("ADAPTIVE_SHRINK_DEC", 4.0))
                        min_box = float(cfg.get("ADAPTIVE_SHRINK_MIN_BOX", 14.0))
                        new_box = tuple(max(min_box, s - dec) for s in box_size)
                        if new_box != box_size:
                            logger.info(f"Adaptive shrink: median {med:.2f} Angstrom -> box {box_size} -> {new_box}")
                            box_size = new_box
            except Exception as _e:
                logger.warning(f"Adaptive shrink skipped: {_e}")

            # Checkpoint done
            if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                try:
                    fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                    checkpoint_mark_done(cfg, pdb_id, stage["name"], fp)
                except Exception:
                    pass

            # ---- outputs ----
            out_dir.mkdir(parents=True, exist_ok=True)
            long_csv = out_dir / f"docking_score_long_{pocket_name}.csv"
            results = score_history.get(stage["name"], {})
            ordered = sorted(
                results.items(),
                key=lambda kv: (
                    kv[1].get("score", 1e9) if isinstance(kv[1].get("score"), (int, float)) else 1e9,
                    os.path.basename(kv[0]),
                ),
            )
            with open(long_csv, "w", newline="", encoding="utf-8") as f:
                f.write(f"# pocket={pocket_name}, center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f})\n")
                w = csv.writer(f)
                w.writerow(["stage", "ligand", "library_id", "corresponding_name", "score", "valid", "reason"])
                for lig, rec in ordered:
                    sc = rec.get("score")
                    sc_str = f"{sc:.2f}" if isinstance(sc, (int, float)) else ""
                    lib_id = _extract_rdk_id(lig) or ""
                    corr = _resolve_name_for_path_or_text(lig, fda_index)
                    w.writerow(
                        [
                            stage["name"],
                            os.path.basename(lig),
                            lib_id,
                            corr,
                            sc_str,
                            int(bool(rec.get("valid", False))),
                            str(rec.get("reason", "")),
                        ]
                    )

            if validated_ligands_last:
                final_pose_validation_and_screenshots(
                    cfg,
                    pdb_id,
                    [stage],
                    receptor_pdbqt,
                    center,
                    validated_ligands_last,
                    score_history,
                    cleaned_pdb,
                    cfg.get("DOCKING_MODE", "benchmark"),
                    logger,
                )

            # Optional PyMOL renders (when DOCKING_MODE == "benchmark")
            if str(cfg.get("DOCKING_MODE", "")).lower() == "benchmark":
                stage_name = stage["name"]
                root_project = Path(cfg.get("OVERALL_DIR", str(Path(cfg["OUTPUT_DIR"]).parent)))

                results_for_stage = score_history.get(stage_name, {})
                best_ctrl_lig, nearest_rdk_lig = pick_control_and_nearest_rdk(
                    results_for_stage, raw_docked, control_stems_lower
                )

                cleaned_pdb_path = str(cleaned_pdb)
                original_pdb_path = str(Path(paths.pdb_path))
                ctrl_pose_path = raw_docked.get(best_ctrl_lig) if best_ctrl_lig else None
                rdk_pose_path = raw_docked.get(nearest_rdk_lig) if nearest_rdk_lig else None

                if bool(cfg.get("DEFER_PYMOL", True)):
                    _enqueue_render_task(
                        pdb_id=pdb_id,
                        stage_name=stage_name,
                        root_project=str(root_project),
                        cleaned_pdb_path=cleaned_pdb_path,
                        original_pdb_path=original_pdb_path,
                        ctrl_pose_path=ctrl_pose_path,
                        rdk_pose_path=rdk_pose_path,
                        exclude_resns=sorted(list(EXCLUDE_HET_IDS)),
                    )
                else:
                    with _RENDER_LOCK:
                        stage_dir_target = root_project / "docked" / pdb_id / stage_name
                        stage_dir_target.mkdir(parents=True, exist_ok=True)

                        _render_native_on_original_pdb(
                            original_pdb=original_pdb_path,
                            outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-native__NATIVE"),
                            exclude_resns=sorted(list(EXCLUDE_HET_IDS)),
                        )
                        if ctrl_pose_path and Path(ctrl_pose_path).is_file():
                            _render_three_views_with_pymol(
                                receptor_path=cleaned_pdb_path,
                                ligand_paths_and_colors=[(ctrl_pose_path, "control", "green")],
                                outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL"),
                            )
                        if rdk_pose_path and Path(rdk_pose_path).is_file():
                            _render_three_views_with_pymol(
                                receptor_path=cleaned_pdb_path,
                                ligand_paths_and_colors=[(rdk_pose_path, "rdk_closest", "magenta")],
                                outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__RDKclosest"),
                            )
                            _render_three_views_with_pymol(
                                receptor_path=original_pdb_path,
                                ligand_paths_and_colors=[(rdk_pose_path, "rdk_closest", "magenta")],
                                outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-rdk__RDKclosest"),
                            )
                            if ctrl_pose_path and Path(ctrl_pose_path).is_file():
                                _render_three_views_with_pymol(
                                    receptor_path=cleaned_pdb_path,
                                    ligand_paths_and_colors=[(ctrl_pose_path, "control", "green"),
                                                             (rdk_pose_path, "rdk_closest", "magenta")],
                                    outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL+RDKclosest"),
                                    label_top_n_res=5,
                                    label_cutoff=5.0,
                                )

            scores_only = [rec.get("score") for rec in results.values() if isinstance(rec.get("score"), (int, float))]
            best = min(scores_only) if scores_only else None
            n_valid = sum(1 for rec in results.values() if rec.get("valid"))
            n_tested = len(results)
            pocket_strength_rows.append(
                [
                    pocket_name,
                    f"{center[0]:.3f}",
                    f"{center[1]:.3f}",
                    f"{center[2]:.3f}",
                    (f"{best:.2f}" if isinstance(best, (int, float)) else ""),
                    n_valid,
                    n_tested,
                ]
            )

    if pocket_strength_rows:
        summary_csv = out_dir / "benchmark_pocket_strength.csv"
        with open(summary_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["pocket", "center_x", "center_y", "center_z", "best_score", "n_valid", "n_tested"])
            for row in pocket_strength_rows:
                w.writerow(row)

# =============================
# Analysis wrapper
# =============================

def _run_auto_analysis(
    docked_root: Path,
    mapping_csv: Path,
    only_pdb: Optional[str],
    score_tol: float,
    center_tol: float,
    rmsd_tol: float,
    include_identity: bool,
    out_dir: Optional[Path] = None,
) -> Tuple[Optional[Path], Optional[Path]]:
    try:
        import benchmark_auto_analysis as ana
    except Exception as e:
        print(f"[analysis] ERROR: could not import benchmark_auto_analysis.py: {e}")
        return None, None
    if not hasattr(ana, "run_analysis"):
        print("[analysis] ERROR: benchmark_auto_analysis.py is missing run_analysis().")
        return None, None
    try:
        return ana.run_analysis(
            docked_root=docked_root,
            mapping_csv=mapping_csv,
            only_pdb=only_pdb,
            score_tol=score_tol,
            center_tol=center_tol,
            rmsd_tol=rmsd_tol,
            include_identity=include_identity,
            out_dir=out_dir,
        )
    except Exception as e:
        print(f"[analysis] ERROR running run_analysis(): {e}")
        return None, None

# =============================
# CLI
# =============================

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Benchmark mode: single ultra-stage docking of likely co-crystal FDA ligands, per pocket."
    )
    # Defaults are None; we resolve from YAML or fallbacks in main()
    p.add_argument("--input-dir", default=None)
    p.add_argument("--out-root", default=None)
    p.add_argument("--prepped", default=None)
    p.add_argument("--mapping", default=None)
    p.add_argument("--max-candidates", type=int, default=999)
    p.add_argument("--exhaustiveness", type=int, default=1)
    p.add_argument("--num-modes", type=int, default=1)
    p.add_argument("--hints", help="Optional manual comma-separated hints (e.g., 'imatinib,STI571')")

    # Outer (proteins) parallelism
    p.add_argument(
        "--jobs",
        type=int,
        default=max(1, (os.cpu_count() or 4) // 2),
        help="Number of proteins to process in parallel",
    )

    # New CPU/parallel shaping flags
    p.add_argument(
        "--total-cpus",
        type=int,
        default=None,
        help="Override: total logical CPUs available to this run (auto-detected if omitted).",
    )
    p.add_argument(
        "--ligand-workers",
        type=int,
        default=None,
        help="Max ligands per protein to dock concurrently. If omitted, auto = floor(total_cpus / jobs).",
    )
    p.add_argument(
        "--blas-threads",
        type=int,
        default=1,
        help="Set OMP/BLAS thread env (OMP/MKL/OPENBLAS/etc) to avoid oversubscription (default: 1).",
    )

    p.add_argument(
        "--alias-from-details",
        nargs="*",
        default=[],
        help=("Path(s) or glob(s) to benchmark_analysis_details.csv to harvest aliases"),
    )

    # Budget CLI flag
    p.add_argument(
        "--max-seconds",
        type=float,
        default=9000.0,
        help="Per-ligand wall-clock budget in seconds for retries/validation (default: 9000.0).",
    )

    # Auto-Analysis options
    p.add_argument(
        "--skip-analysis",
        dest="run_analysis",
        action="store_false",
        help="Skip the post-run Auto-Analysis.",
    )
    p.set_defaults(run_analysis=True)
    p.add_argument("--analysis-only-pdb", default=None, help="Optional: restrict analysis to one PDB (e.g., 5MO4).")
    p.add_argument("--analysis-score-tol", type=float, default=1.0, help="Score tolerance (kcal/mol) for a point.")
    p.add_argument("--analysis-center-tol", type=float, default=1.0, help="Centroid distance tolerance (Å).")
    p.add_argument("--analysis-rmsd-tol", type=float, default=3.0, help="RMSD tolerance (Å).")
    p.add_argument("--analysis-no-identity", action="store_true", help="Exclude identity-match from scoring.")
    p.add_argument("--analysis-out", default=None, help="Optional output dir for analysis CSVs.")
    p.add_argument("--only-pdb", default=None, help="Comma-separated PDB IDs to dock, e.g. '6GQO,4XUF'.")
    return p









# =============================
# CPU detection helpers
# =============================

def _env_int(*names: str) -> Optional[int]:
    """Return the first present env var (int) among names, else None."""
    for n in names:
        v = os.environ.get(n)
        if v:
            try:
                return int(v)
            except Exception:
                pass
    return None

def _detect_available_cpus() -> int:
    """
    Conservative CPU detector for clusters/desktops.
    Priority:
      1) SLURM_CPUS_PER_TASK
      2) SLURM_CPUS_ON_NODE
      3) NSLOTS / PBS_NP
      4) OMP_NUM_THREADS
      5) os.cpu_count()
    """
    for v in (
        _env_int("SLURM_CPUS_PER_TASK"),
        _env_int("SLURM_CPUS_ON_NODE"),
        _env_int("NSLOTS", "PBS_NP"),
        _env_int("OMP_NUM_THREADS"),
        os.cpu_count() or 1,
    ):
        if v and v > 0:
            return int(v)
    return 1

def _set_thread_env(n: int) -> None:
    """Pin common threaded libs to n to avoid hidden oversubscription."""
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_MAX_THREADS",
        "BLIS_NUM_THREADS",
    ):
        os.environ.setdefault(key, str(n))




# =============================
# Main (Linux/BRCF friendly)
# =============================

def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_argparser().parse_args(argv)

    # Load YAML config first so we can use BRCF paths by default
    cfg = load_inputs()
    validate_config(cfg)

    # Helper: first existing path in a list (returns Path or None)
    def _first_existing(*cands) -> Optional[Path]:
        for c in cands:
            if not c:
                continue
            p = Path(str(c))
            if p.exists():
                return p
        return None

    # 1) Resolve INPUT / OUTPUT / PREPPED
    input_dir = Path(args.input_dir or cfg.get("INPUT_DIR", ""))
    out_root  = Path(args.out_root  or cfg.get("DOCKED_DIR") or cfg.get("OUTPUT_DIR", ""))
    prepped   = Path(args.prepped   or cfg.get("OUTPUT_LIGANDS_DIR", ""))

    # 2) Resolve mapping CSV in priority order
    brcf_fixed_mapping = Path("/stor/home/mpg2352/atlas/code/protein_automation/fda_mapping_from_pdbqt.csv")
    mapping_path = (
        Path(args.mapping) if args.mapping else
        _first_existing(
            cfg.get("FDA_MAPPING_CSV"),
            prepped / "fda_mapping_from_pdbqt.csv",
            Path(cfg.get("OVERALL_DIR", Path(out_root).parent if out_root else Path.cwd())) / "fda_mapping_from_pdbqt.csv",
            Path.cwd() / "fda_mapping_from_pdbqt.csv",
            brcf_fixed_mapping,
        )
    )

    # 3) Sanity checks with helpful messages
    if not input_dir.exists():
        print(f"[benchmark] INPUT_DIR does not exist: {input_dir}")
        return

    if not prepped.exists():
        print(f"[benchmark] PREPPED library does not exist: {prepped}")
        return

    if not out_root:
        out_root = Path(cfg.get("OUTPUT_DIR", "")) if cfg.get("OUTPUT_DIR") else (Path.cwd() / "benchmarks")
    out_root.mkdir(parents=True, exist_ok=True)

    if not mapping_path or not Path(mapping_path).exists():
        print(
            "[benchmark] Mapping CSV not found.\n"
            "  Tried (in order):\n"
            f"    --mapping={args.mapping}\n"
            f"    cfg['FDA_MAPPING_CSV']={cfg.get('FDA_MAPPING_CSV')}\n"
            f"    {prepped / 'fda_mapping_from_pdbqt.csv'}\n"
            f"    {Path(cfg.get('OVERALL_DIR', Path(out_root).parent)) / 'fda_mapping_from_pdbqt.csv'}\n"
            f"    {Path.cwd() / 'fda_mapping_from_pdbqt.csv'}\n"
            f"    {brcf_fixed_mapping}\n"
            "  Fix by passing --mapping /path/to/fda_mapping_from_pdbqt.csv or setting FDA_MAPPING_CSV in your YAML."
        )
        return

    # 4) Build MappingIndex now that we know the file exists
    mapping = MappingIndex(Path(mapping_path))

    # 5) Freeze resolved paths back into cfg for downstream calls
    cfg = dict(cfg)
    cfg["INPUT_DIR"] = str(input_dir)
    cfg["DOCKED_DIR"] = str(out_root)
    cfg["OUTPUT_DIR"] = str(out_root)
    cfg["OUTPUT_LIGANDS_DIR"] = str(prepped)
    cfg["DOCKING_MODE"] = "benchmark"
    cfg.setdefault("OVERALL_DIR", str(Path(out_root).parent))
    cfg["BENCH_MAX_SECONDS"] = float(args.max_seconds)

    # Optionally harvest aliases from prior results
    alias_sources = args.alias_from_details or []
    if alias_sources:
        try:
            extra_aliases = load_aliases_from_details_csv(alias_sources)
            if extra_aliases:
                global CHEMCOMP_ALIAS
                CHEMCOMP_ALIAS = merge_aliases_into_chemcomp(CHEMCOMP_ALIAS, extra_aliases)
                # keep values normalized+deduped
                for _k, _vals in list(CHEMCOMP_ALIAS.items()):
                    CHEMCOMP_ALIAS[_k] = sorted(set(v.lower() for v in _vals))
                total_added = sum(len(v) for v in extra_aliases.values())
                print(f"[alias] Expanded CHEMCOMP_ALIAS with {total_added} aliases across {len(extra_aliases)} HET codes.")
        except Exception as e:
            print(f"[alias] WARN: could not augment aliases: {e}")

    # enumerate PDBs
    pdb_files = [f for f in os.listdir(cfg["INPUT_DIR"]) if f.lower().endswith(".pdb") and "_nolig" not in f.lower()]
    print(f"[benchmark] proteins queued: {len(pdb_files)} from {cfg['INPUT_DIR']}")

    # ONLY-PDB filter (optional)
    if args.only_pdb:
        wanted = {s.strip().upper() for s in str(args.only_pdb).split(",") if s.strip()}
        before = len(pdb_files)
        pdb_files = [f for f in pdb_files if os.path.splitext(f)[0].upper() in wanted]
        print(f"[benchmark] --only-pdb={sorted(wanted)} -> matched {len(pdb_files)} of {before} files:")
        for f in pdb_files:
            print(f"  - {f}")

    fda_index = load_library_index(Path(mapping_path))
    manual_hints = [h.strip() for h in (args.hints or "").split(",") if h.strip()] or None

    # ---- Dynamic CPU-saturating scheduler (outer proteins × inner ligands) ----
    from collections import deque

    total_cpus = int(args.total_cpus) if args.total_cpus else _detect_available_cpus()
    _set_thread_env(int(args.blas_threads))

    # We’ll defer PyMOL to the very end by default so docking never waits on renders.
    cfg["DEFER_PYMOL"] = True
    # <<< global token pool used by ALL ligand tasks (live rebalancing of cpu usage) >>>
    global_ligand_sem = threading.BoundedSemaphore(total_cpus)
    cfg["GLOBAL_LIGAND_SEM"] = global_ligand_sem
    # Quick, cheap prescan: estimate per-protein candidate counts to weight scheduling.
    # (We intentionally skip metabolite→parent expansion here for speed.)
    prescan: List[Tuple[str, int]] = []
    for pdb in pdb_files:
        pdb_path = Path(cfg["INPUT_DIR"]) / pdb
        try:
            het_ids, het_names = parse_pdb_het_hints(pdb_path)
        except Exception:
            het_ids, het_names = [], []
        hints = list(het_names) + list(het_ids)
        for het in het_ids:
            hints.extend(CHEMCOMP_ALIAS.get(het.upper(), []))
        if manual_hints:
            hints.extend(manual_hints)
        try:
            est = len(select_candidates_for_protein(
                mapping=mapping,
                hints=hints,
                prepped_dir=prepped,
                extra_parent_ids=[],            # keep fast
                max_candidates=int(args.max_candidates),
            ))
        except Exception:
            est = 8
        prescan.append((pdb, max(1, est)))

    # Big jobs first helps keep CPUs busy as small jobs finish.
    prescan.sort(key=lambda t: t[1], reverse=True)
    q = deque(prescan)

    # Outer cap: how many proteins can run at once
    max_protein_jobs = max(1, min(int(args.jobs) if args.jobs else total_cpus, len(pdb_files), total_cpus))

    print(f"[parallel] total_cpus={total_cpus} | max_proteins={max_protein_jobs} | queued={len(prescan)}")

    def _work(pdb_file: str, inner_for_this: int):
        cfg_local = dict(cfg)
        cfg_local["GLOBAL_LIGAND_SEM"] = cfg["GLOBAL_LIGAND_SEM"]
        # Set a high per-protein cap so the semaphore (not the pool size) is the real limiter.
        cfg_local["MAX_PARALLEL_JOBS"] = max(int(cfg_local.get("MAX_PARALLEL_JOBS") or 0), total_cpus)
        with _defer_pymol_capture_calls(True):
            run_benchmark_for_protein(
                cfg=cfg_local,
                mapping=mapping,
                pdb_file=pdb_file,
                prepped_dir=prepped,
                out_root=out_root,
                exhaustiveness=int(args.exhaustiveness),
                num_modes=int(args.num_modes),
                max_candidates=int(args.max_candidates),
                manual_hints=manual_hints,
                fda_index=fda_index,
            )
        return pdb_file

    import traceback

    # Submit loop: we keep track of "CPU slots" allocated to each running protein.
    running: Dict[object, Tuple[str, int]] = {}
    sum_inner = 0

    def _alloc_plan(remaining_slots: int, remaining_jobs: int) -> int:
        """Greedy: give each new protein at least floor, first few get +1 to absorb remainder."""
        base = max(1, remaining_slots // max(1, remaining_jobs))
        extra = remaining_slots - base * max(1, remaining_jobs)
        return base + (1 if extra > 0 else 0)

    with ThreadPoolExecutor(max_workers=max_protein_jobs) as pool:
        # Seed as many as we can to fill CPU budget
        while q and len(running) < max_protein_jobs and sum_inner < total_cpus:
            remaining_slots = total_cpus - sum_inner
            remaining_jobs = max_protein_jobs - len(running)
            inner = _alloc_plan(remaining_slots, remaining_jobs)
            pdb, est = q.popleft()
            fut = pool.submit(_work, pdb, inner)
            running[fut] = (pdb, inner)
            sum_inner += inner

        done = failed = 0
        total = len(pdb_files)

        while running:
            fut = next(as_completed(list(running.keys())))
            pdb, inner = running.pop(fut)
            sum_inner -= inner
            try:
                fut.result()
                done += 1
                print(f"[benchmark] ✔ finished {pdb} (used_inner={inner}) ({done}/{total})")
            except Exception as e:
                failed += 1
                tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
                print(f"[benchmark] ✖ error on {pdb} (used_inner={inner}): {e} ({done + failed}/{total})\n{tb}")

            # Refill immediately, redistributing freed CPU to the next proteins
            while q and len(running) < max_protein_jobs and sum_inner < total_cpus:
                remaining_slots = total_cpus - sum_inner
                remaining_jobs  = max_protein_jobs - len(running)
                inner = _alloc_plan(remaining_slots, remaining_jobs)
                next_pdb, _est = q.popleft()
                fut2 = pool.submit(_work, next_pdb, inner)
                running[fut2] = (next_pdb, inner)
                sum_inner += inner

    # ----- Post-run: render everything in parallel (no docking blocked) -----
    # Default workers: env PYMOL_PARALLEL (if set) else 4 as a safe HPC default
    try:
        render_workers = int(os.environ.get("PYMOL_PARALLEL", ""))
        if render_workers <= 0:
            render_workers = 4
    except Exception:
        render_workers = 4
    run_deferred_captures(render_workers)
    run_deferred_renders(render_workers)
    # optional post-run analysis
    if args.run_analysis:
        docked_root = Path(cfg.get("OVERALL_DIR", str(Path(out_root).parent))) / "docked"
        if not docked_root.is_dir():
            print(f"[analysis] Docked root not found at {docked_root} — skipping.")
        else:
            details, summary = _run_auto_analysis(
                docked_root=docked_root,
                mapping_csv=Path(mapping_path),
                only_pdb=(args.analysis_only_pdb or None),
                score_tol=float(args.analysis_score_tol),
                center_tol=float(args.analysis_center_tol),
                rmsd_tol=float(args.analysis_rmsd_tol),
                include_identity=not bool(args.analysis_no_identity),
                out_dir=(Path(args.analysis_out) if args.analysis_out else None),
            )
            if details and summary:
                print(f"[analysis] ✅ Details: {details}")
                print(f"[analysis] ✅ Summary: {summary}")
            else:
                print("[analysis] ❌ Analysis did not produce outputs.")

if __name__ == "__main__":
    main()
