# -*- coding: utf-8 -*-
"""
Automated protein preparation pipeline (Windows/WSL-friendly)

This module prepares a receptor structure from a raw PDB by performing:
  1) Alternate-conformation filtering (altLoc) with deterministic selection
  2) Removal of nonstandard residues (policy-aware keep/drop of cofactors/metals/waters)
  3) Element-column fixes to ensure PDB compliance (via activesite/YAML)
  4) Optional loop/residue completion via MODELLER
  5) Sanity checks for incomplete residues
  6) Phenix cleaning passes (auto-detect Linux Phenix, Windows Phenix, or skip—WSL-safe)
  7) Hydrogen sanity cleanup (CONECT- and geometry-based)
  8) Chain validation (ensure at least one CA-containing chain)
  9) Protonation via Reduce with automatic fallbacks (temp retry/OpenBabel)
 10) Final polish via phenix.pdbtools when available; otherwise degrade gracefully
 11) Receptor PDBQT preparation (Meeko preferred; ADT fallback)

Design notes
------------
• Single source of truth for element handling: **activesite** helpers (YAML-driven). No local element inference.
• WSL-aware external calls.
• Logging instead of prints; short, numbered steps for easier debugging.
• Never touch PDBQT with any PDB element-fixing logic.

Prerequisites (recommended)
---------------------------
• Phenix (Linux: phenix.pdbtools; Windows: phenix.pdbtools.bat, phenix.python.bat)
• Reduce (reduce.exe or reduce) — often bundled with Phenix
• MODELLER (licensed) — optional but recommended for loop filling
• Open Babel (obabel) — fallback for hydrogen addition
• Biopython — for PDB parsing
• MGLTools (for prepare_receptor4.py) and/or Meeko (mk_prepare_receptor)
"""
from __future__ import annotations

import json
import platform
from collections import Counter, defaultdict
from math import sqrt
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union
import os, sys, shutil, subprocess, logging

# -------- Project-local imports --------
from installation import load_config
from logger_setup import setup_logger
# Single source of truth (YAML-backed) — do not re-implement locally
from activesite import (
    fix_element_columns_in_file,
    scan_helium_counts,
    assert_no_helium_in_hydrogen_names,
    get_atom_rules,
    fix_pdb_elements,
    scan_helium_counts_with_hits,
    format_pdb_atom_debug,
    rules_version,
)
ALIASES = get_atom_rules()
# Normalize to a dict so existing RULES.get(...) calls work
RULES = ALIASES.__dict__ if hasattr(ALIASES, "__dict__") else dict(ALIASES)

# Toggle: default ON; set HELIUM_POSTWRITE_VERBOSE=0 to suppress the info line
_HE_POSTWRITE_VERBOSE = (os.environ.get("HELIUM_POSTWRITE_VERBOSE", "1") != "0")

def load_aliases(): #little shim, to fix later 
    return get_atom_rules()

def _flatten_semicolons(items):
    out = []
    for item in items or []:
        # items might be lines like "A; B; C"
        parts = [p.strip().upper() for p in str(item).split(";") if p.strip()]
        out.extend(parts)
    return out

# Sets from YAML
_RETAIN = set(_flatten_semicolons(RULES.get("retain_in_receptor_resnames", [])))

# Optional: recognize “waters” using the retain list (makes policy consistent)
_WATER_NAMES = {w for w in _RETAIN if w in {"HOH","WAT","DOD","H2O","TIP","TIP3","SOL"}}

# For diagnostics only (not for logic): which entries look like elemental ions
# We infer “metals/halides” from retain entries that are 1–2 char and in element tables.
# --- Canonical element tokens (upper-cased 1�2 letter symbols) ---
def _canonize_two_letter(xs) -> set[str]:
    out = set()
    for line in (xs or []):
        for tok in str(line).split(";"):
            t = tok.strip()
            if t:
                out.add(t.upper())
    return out

_ES      = RULES.get("element_sets", {}) or {}
_ONE     = {str(s).upper() for s in (_ES.get("one_letter_elements") or [])}
_TWORAW  = _ES.get("two_letter_elements", []) or []
_TWO     = _canonize_two_letter(_TWORAW)
# All allowed element tokens
_ELEM_CANON = _ONE | _TWO

def _is_element_token(sym):
    s = str(sym).strip().upper()
    return (len(s) in (1,2)) and (s in _ELEM_CANON)

# Treat as “ion-like” if it looks like an element token and is in the retain set
def _is_retained_ion(resname: str) -> bool:
    r = (resname or "").upper()
    return r in _RETAIN and _is_element_token(r)


# =============================
# Configuration helpers
# =============================
config = load_config()  # load once at import

# Case-insensitive lookup with alias support
# minimal _cfg using your existing installation.load_config()
import os
from installation import load_config

_CONFIG = load_config()

def _cfg(key: str, default: str = "", legacy_key: str | None = None) -> str:
    """env > config (key) > config (legacy_key) > default"""
    v = os.environ.get(key)
    if v not in (None, ""):
        return v
    if key in _CONFIG and str(_CONFIG.get(key)) != "":
        return str(_CONFIG.get(key))
    if legacy_key and legacy_key in _CONFIG and str(_CONFIG.get(legacy_key)) != "":
        return str(_CONFIG.get(legacy_key))
    return default

PHENIX_DIR         = _cfg("PHENIX_DIR", "", "phenix_dir")
PHENIX_LIB_PATH    = _cfg("PHENIX_LIB_PATH", "", "phenix_lib_path")
PHENIX_CLEAN_SCRIPT= _cfg("PHENIX_CLEAN_SCRIPT", "", "phenix_clean_script")
INPUT_DIR          = _cfg("INPUT_DIR", ".")

MGLTOOLS_PYTHON         = _cfg("MGLTOOLS_PYTHON", "")
PREPARE_RECEPTOR_SCRIPT = _cfg("PREPARE_RECEPTOR_SCRIPT", "")
USE_MEEKO = str(_cfg("USE_MEEKO", "")).lower() in ("1", "true", "yes")

# prefer explicit env var; else fall back to obabel on PATH
OPENBABEL_PATH = os.environ.get("OPENBABEL_PATH") or shutil.which("obabel") or ""

# Windows .bat wrappers (if using Windows Phenix)
PDBTOOLS_BAT     = str(Path(PHENIX_DIR) / "phenix.pdbtools.bat") if PHENIX_DIR else ""
MOLPROBITY_BAT   = str(Path(PHENIX_DIR) / "phenix.molprobity.bat") if PHENIX_DIR else ""
PHENIX_PYTHON_BAT= str(Path(PHENIX_DIR) / "phenix.python.bat")   if PHENIX_DIR else ""

def _pick_reduce_exe() -> str:
    """
    Choose the actual 'reduce' binary robustly.
    Precedence:
      1) $REDUCE_EXE (env) or config value
      2) User's local build: /stor/home/mpg2352/atlas/tools/reduce/reduce_src/reduce
      3) Phenix conda_base/bin/reduce (next to PHENIX_DIR)
      4) reduce on PATH
    Avoids 'reduce.python' wrapper entirely.
    """
    # 1) Explicit env/config
    explicit = (
        os.environ.get("REDUCE_EXE")
        or os.environ.get("REDUCE_BIN")
        or _cfg("REDUCE_EXE", "", "reduce_exe")
    )
    if explicit and Path(explicit).exists():
        return explicit

    # 2) Known local build (user-specific)
    user_local = Path("/stor/home/mpg2352/atlas/tools/reduce/reduce_src/reduce")
    if user_local.exists() and os.access(str(user_local), os.X_OK):
        return str(user_local)

    # 3) Phenix conda-base candidate near PHENIX_DIR
    if PHENIX_DIR:
        phenix_root = Path(PHENIX_DIR).resolve().parent  # .../phenix-1.21.2-5419
        cb = phenix_root / "conda_base" / "bin" / "reduce"
        if cb.exists() and os.access(str(cb), os.X_OK):
            return str(cb)

    # 4) PATH fallback
    which = shutil.which("reduce") or shutil.which("reduce.exe")
    if which:
        return which

    # Last resort (let subprocess resolve; will likely fail clearly)
    return "reduce"

REDUCE_EXE = _pick_reduce_exe()
logging.info("Using Reduce at: %s", REDUCE_EXE)



def _post_write_element_guard(step: str, pdb_path: str | Path) -> None:
    """
    Normalize element columns after each receptor write step and log helium status.
    """
    p = Path(pdb_path)
    try:
        fixed = fix_element_columns_in_file(p, dst_path=p, rewrite_atoms=False)
    except Exception as e:
        logging.warning(
            "[helium] stage=post_write step=%s file=%s He->H=? residual_He=? note=elemfix_error:%s",
            step, p, e
        )
        return

    txt = p.read_text(encoding="utf-8", errors="ignore")
    he_count = scan_helium_counts(txt)
    logging.info(
        "[helium] stage=post_write step=%s file=%s He->H=%s residual_He=%d",
        step, p, fixed if isinstance(fixed, int) else -1, he_count
    )


def _meeko_preflight_or_fail(pdb_input: str | Path, work_dir: str | Path) -> Path:
    """
    Ensure the exact PDB Meeko will read contains no 'He'.
    Saves copy to work/meeko_input_pre_sanitize.pdb for forensics.
    """
    src = Path(pdb_input).resolve()
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    try:
        fix_element_columns_in_file(src, dst_path=src, rewrite_atoms=False)
    except Exception as e:
        logging.warning("[helium] stage=preflight file=%s note=elemfix_exception:%s", src, e)

    txt = src.read_text(encoding="utf-8", errors="ignore")
    he_count, first_hits = scan_helium_counts_with_hits(txt, max_hits=5)

    # Save what Meeko will actually see
    prefile = work / "meeko_input_pre_sanitize.pdb"
    prefile.write_text(txt, encoding="utf-8")

    if he_count > 0:
        logging.error(
            "[helium] stage=preflight file=%s He_count=%d rules=%s",
            src, he_count, rules_version(),
        )
        for hit in first_hits:
            logging.error("[helium] offender %s", hit)
        raise RuntimeError("helium_preflight_failed")

    return src




def _helium_postwrite_counter(step_name: str, pdb_path: str | Path) -> None:
    """
    Run the element fixer on the freshly written receptor PDB and emit:
      [helium] stage=post_write step=<name> file=<path> He->H=<n> residual_He=<n>
    Fail fast if residual_He > 0.
    """
    p = Path(pdb_path)
    try:
        before = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        before = ""

    # Count BEFORE
    he_before = scan_helium_counts(before.splitlines()) if before else 0

    # Fix elements in-place on the exact file we pass forward
    fix_element_columns_in_file(p, p, rewrite_atoms=True)

    try:
        after = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        after = ""

    # Count AFTER and compute delta
    he_after = scan_helium_counts(after.splitlines()) if after else 0
    delta = max(0, he_before - he_after)

    if _HE_POSTWRITE_VERBOSE:
        logging.info(
            "[helium] stage=post_write step=%s file=%s He->H=%d residual_He=%d",
            step_name, str(p), delta, he_after
        )

    if he_after > 0:
        # Fail fast so we can see which step leaked helium
        raise RuntimeError("helium_residual_post_write")




try:
    from rdkit import Chem  # type: ignore
    _HAS_RDKIT = True
except Exception:
    _HAS_RDKIT = False

# =============================
# OS / WSL helpers
# =============================

def _is_wsl() -> bool:
    try:
        return "microsoft" in platform.release().lower() or "WSL_INTEROP" in os.environ
    except Exception:
        return False

def _bin_on_path(name: str) -> bool:
    return shutil.which(name) is not None

def _win_path(p: Union[str, Path]) -> str:
    s = str(p)
    if not _is_wsl():
        return s
    try:
        out = subprocess.check_output(["wslpath", "-w", s], text=True).strip()
        return out
    except Exception:
        return s

def _powershell(cmd: str) -> subprocess.CompletedProcess:
    return subprocess.run([
        "powershell.exe", "-NoProfile", "-NonInteractive", "-Command", cmd
    ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)





def _as_path(p) -> Path:
    return p if isinstance(p, Path) else Path(p)

def _cfg_env_or_default(key: str, default: Optional[str] = None) -> Optional[str]:
    """Prefer env var, then config.txt (sibling of this file), else default."""
    v = os.environ.get(key)
    if v:
        return v
    try:
        root = Path(__file__).resolve().parent
        cfg = root / "config.txt"
        if cfg.is_file():
            for line in cfg.read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, val = line.split("=", 1)
                if k.strip() == key:
                    return val.strip()
    except Exception:
        pass
    return default

def _canon_base(output_root: Path, pdb_id: str) -> Path:
    """Canonical per-protein base: processed_pdbs/<PDB>"""
    output_root = _as_path(output_root)
    return (output_root / pdb_id.upper()).resolve()

def _merge_dir(src: Path, dst: Path) -> None:
    """Merge src directory into dst; remove src after moving."""
    src, dst = src.resolve(), dst.resolve()
    if not src.exists():
        return
    dst.mkdir(parents=True, exist_ok=True)
    for root, dirs, files in os.walk(src):
        r = Path(root)
        rel = r.relative_to(src)
        (dst / rel).mkdir(parents=True, exist_ok=True)
        for d in dirs:
            (dst / rel / d).mkdir(parents=True, exist_ok=True)
        for f in files:
            s = r / f
            t = (dst / rel / f)
            if t.exists():
                try:
                    if s.stat().st_size == t.stat().st_size:
                        continue
                except Exception:
                    pass
            shutil.move(str(s), str(t))
    try:
        shutil.rmtree(src)
    except Exception:
        pass

def fold_legacy_layout(pdb_id: str, output_root) -> None:
    """
    Extended: migrate uppercase legacy dirs into canonical tree:
      <PDB>_NOLIG            -> processed_pdbs/<PDB>/nolig
      <PDB>_CLEANED_LIGANDS  -> processed_pdbs/<PDB>/ligands_raw
      <PDB>_nolig(.pdb)      -> processed_pdbs/<PDB>/nolig/<PDB>_nolig_phenix_clean.pdb
    """
    try:
        root = _as_path(output_root).resolve()
        pdb_idU = pdb_id.upper()
        base = _canon_base(root, pdb_idU)
        (base / "nolig").mkdir(parents=True, exist_ok=True)
        (base / "ligands_raw").mkdir(parents=True, exist_ok=True)

        legacy_dirs = [
            (root / f"{pdb_id}_nolig",                   base / "nolig"),
            (root / f"{pdb_id.lower()}_nolig",           base / "nolig"),
            (root / f"{pdb_idU}_NOLIG",                  base / "nolig"),
            (root / f"{pdb_id}_cleaned_ligands",         base / "ligands_raw"),
            (root / f"{pdb_id.lower()}_cleaned_ligands", base / "ligands_raw"),
            (root / f"{pdb_idU}_CLEANED_LIGANDS",        base / "ligands_raw"),
        ]
        for src, dst in legacy_dirs:
            if src.exists():
                logging.info("Migrating legacy directory %s -> %s", src, dst)
                _merge_dir(src, dst)

        candidates_files = [
            (root / f"{pdb_id}_nolig.pdb",         base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb"),
            (root / f"{pdb_id.lower()}_nolig.pdb", base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb"),
        ]
        for src, dst in candidates_files:
            if src.exists() and not dst.exists():
                logging.info("Moving legacy file %s -> %s", src, dst)
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src), str(dst))
    except Exception as e:
        logging.warning("fold_legacy_layout (extended) failed for %s: %s", pdb_id, e)


# =============================
# Ligand pristine reference (optional; for RMSD downstream)
# =============================

def _write_pristine_reference(pdb_lig_path: Path) -> None:
    """Write pristine ligand copy next to ligands_raw (SDF + light atom map JSON)."""
    ref_dir = pdb_lig_path.parent.parent / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    ref_sdf = ref_dir / (pdb_lig_path.stem + ".sdf")
    ref_json = ref_dir / (pdb_lig_path.stem + ".map.json")

    if _HAS_RDKIT:
        m = Chem.MolFromPDBFile(str(pdb_lig_path), sanitize=False, removeHs=False)
        if m is not None:
            w = Chem.SDWriter(str(ref_sdf)); w.write(m); w.close()
            amap = []
            with open(pdb_lig_path, "r", encoding="utf-8", errors="ignore") as f:
                for ln in f:
                    if ln.startswith(("ATOM","HETATM")):
                        amap.append({
                            "name": ln[12:16].strip(),
                            "element": (ln[76:78].strip() or ln[12:16].strip()[:1].upper())
                        })
            json.dump({"atoms": amap}, open(ref_json, "w"), indent=2)
            return

    # Fallback minimal SDF stub
    with open(ref_sdf, "w", encoding="utf-8") as out:
        out.write(f"{pdb_lig_path.stem}\n  -Pristine-\n\n")
        out.write("$$$$\n")

# =============================
# Canonical per-protein directory layout
# =============================

def canon_paths(pdb_id: str, output_root: Union[str, Path]) -> Dict[str, Path]:
    root = Path(output_root).resolve()
    base = root / pdb_id.upper()
    return {
        "protein_root": base,
        "raw":          base / "raw",
        "work":         base / "work",
        "ligands_raw":  base / "ligands_raw",
        "nolig":        base / "nolig",
        "receptor":     base / "receptor",
    }

# =============================
# Ligand extraction (single source lives here)
# =============================

def extract_ligands_from_filtered(filtered_pdb: Union[str, Path], out_dir: Union[str, Path]) -> List[Path]:
    """Extract non-water HETATM residues into individual PDBs (RES_CHAINRESI.pdb),
    and run text-level element repair on each (YAML rules)."""
    outd = Path(out_dir)
    outd.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    cur_key = None
    bucket: List[str] = []

    def flush():
        nonlocal bucket, cur_key, written
        if not bucket or cur_key is None:
            return
        resname, chain, resseq, icode = cur_key
        name = f"{resname}_{chain}{int(resseq)}"
        outp = outd / f"{name}.pdb"
        with open(outp, "w") as w:
            w.write(f"REMARK Extracted {name}\n")
            for ln in bucket:
                w.write(ln)
            w.write("TER\nEND\n")
        # YAML-backed element repair (PDB only)
        try:
            fix_element_columns_in_file(outp, outp)
        except Exception as e:
            logging.warning("Element-fix skipped for %s: %s", outp.name, e)
        written.append(outp)
        bucket = []
        cur_key = None
        try:
            _write_pristine_reference(outp)
        except Exception as e:
            logging.warning("Could not write pristine reference for %s: %s", outp.name, e)

    with open(filtered_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith("HETATM"):
                continue
            resname = ln[17:20].strip().upper()
            if resname in _WATER_NAMES:
                continue  # never extract waters
            if resname in _RETAIN:
                continue  # don't extract cofactors/metals/ions you keep with protein
            chain = ln[21]
            resseq = ln[22:26].strip() or "0"
            icode = ln[26]
            key = (resname, chain, resseq, icode)
            if cur_key is None:
                cur_key = key
            if key != cur_key:
                flush()
                cur_key = key
            bucket.append(ln)
    flush()
    logging.info("Extracted %d ligand residues to %s", len(written), out_dir)
    return written



def element_fix_all_in_dir(dir_path: Union[str, Path], rewrite_atoms: bool = False) -> int:
    """
    Run the text-level element column fixer on every *.pdb under dir_path.
    Returns the number of files rewritten.
    """
    d = Path(dir_path)
    if not d.exists():
        return 0
    n = 0
    for p in d.rglob("*.pdb"):
        try:
            fix_element_columns_in_file(p, dst_path=p, rewrite_atoms=rewrite_atoms)
            n += 1
        except Exception as e:
            logging.warning("element_fix_all_in_dir skipped %s: %s", p, e)
    logging.info("Element column sweep fixed %d PDB files under %s", n, d)
    return n

def expose_ligand_intermediates_for_debug(src_dir: Union[str, Path],
                                          link_dir: Union[str, Path]) -> None:
    """
    Create/refresh a symlink 'link_dir' -> 'src_dir' for easy browsing.
    On Windows, tries directory junction fallback if symlink fails.
    """
    src = Path(src_dir).resolve()
    dst = Path(link_dir)

    try:
        if dst.is_symlink() or dst.exists():
            try:
                if dst.is_symlink():
                    dst.unlink()
                else:
                    shutil.rmtree(dst)
            except Exception:
                pass

        if os.name == "nt":
            try:
                os.symlink(src, dst, target_is_directory=True)
            except OSError:
                cmd = ['cmd', '/c', 'mklink', '/J', str(dst), str(src)]
                subprocess.run(cmd, check=True)
        else:
            os.symlink(src, dst, target_is_directory=True)
        logging.info("Exposed ligand intermediates: %s -> %s", dst, src)
    except Exception as e:
        logging.warning("Could not create debug link %s -> %s: %s", dst, src, e)


# =============================
# Element & Hydrogen utilities (PDB only)
# =============================

def hydrogenation_status(pdb_path: Union[str, Path]) -> Tuple[str, int, int, float]:
    h = heavy = 0
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            el = line[76:78].strip().upper()
            if el == "H":
                h += 1
            else:
                heavy += 1
    ratio = h / max(heavy, 1)
    if h == 0:
        return "NO_H", h, heavy, ratio
    if ratio < 0.20:
        return "SUSPECT_LOW_H", h, heavy, ratio
    return "HAS_H", h, heavy, ratio


def file_contains_hydrogens(pdb_path: Union[str, Path]) -> bool:
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                el = line[76:78].strip().upper()
                if el in {"H", "D", "T"}:
                    return True
    except Exception as e:
        logging.warning("Could not read %s to check for H atoms: %s", pdb_path, e)
    return False



def strip_monoatomic_ions_inplace(pdb_path: Union[str, Path],
                                  keep_resnames: Optional[Set[str]] = None) -> int:
    """
    Remove single-atom HET residues that look like elemental ions (Na, Cl, Zn, ...),
    UNLESS resname is in keep_resnames (e.g., YAML retain list).
    Returns number of residues removed.
    """
    keep_resnames = set(keep_resnames or [])
    by_res = {}
    lines = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                resname = ln[17:20].strip().upper()
                key = (ln[21], ln[22:26], ln[26], resname)
                by_res.setdefault(key, []).append(ln)
            else:
                lines.append(ln)

    removed = 0
    for key, atms in by_res.items():
        chain, resi, icode, resname = key
        if resname in keep_resnames:
            lines.extend(atms); continue
        if len(atms) == 1 and _is_element_token(resname):
            removed += 1
            continue
        lines.extend(atms)

    with open(pdb_path, "w", encoding="utf-8") as w:
        w.writelines(lines)

    if removed:
        logging.info("Stripped %d monoatomic ions from %s", removed, pdb_path)
    return removed

def detect_catalytic_metals(pdb_path: Union[str, Path]) -> Set[str]:
    """
    Return the set of retained ion-like resnames present in the PDB (e.g., ZN, MG).
    Uses YAML retain list + element tokens. For diagnostics only.
    """
    RETAIN = set()
    for item in RULES.get("retain_in_receptor_resnames", []):
        for t in str(item).split(";"):
            tok = t.strip().upper()
            if tok:
                RETAIN.add(tok)

    ionic = set()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith("HETATM"):
                continue
            resname = ln[17:20].strip().upper()
            if (resname in RETAIN) and _is_element_token(resname):
                ionic.add(resname)

    if ionic:
        logging.info("Catalytic/retained ions present: %s", sorted(ionic))
    return ionic

def compare_ion_presence_between_pdb_and_pdbqt(clean_pdb: Union[str, Path],
                                                receptor_pdbqt: Union[str, Path]) -> None:
    """
    Log a warning if an ion present in the cleaned PDB is missing in receptor PDBQT.
    Uses YAML retain list + element tokens as "ions" definition.
    """
    def ions_in_pdb(p: Union[str, Path]) -> Set[Tuple[str, str]]:
        ions = set()
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith("HETATM"):
                    continue
                resname = ln[17:20].strip().upper()
                chain = ln[21]
                resi  = ln[22:26].strip()
                if _is_element_token(resname):
                    ions.add((resname, f"{chain}:{resi}"))
        return ions

    def ions_in_pdbqt(p: Union[str, Path]) -> Set[Tuple[str, str]]:
        ions = set()
        if not Path(p).exists():
            return ions
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                resname = ln[17:20].strip().upper()
                chain = ln[21]
                resi  = ln[22:26].strip()
                if _is_element_token(resname):
                    ions.add((resname, f"{chain}:{resi}"))
        return ions

    pdb_ions   = ions_in_pdb(clean_pdb)
    pdbqt_ions = ions_in_pdbqt(receptor_pdbqt)
    missing = pdb_ions - pdbqt_ions
    if missing:
        logging.warning("Ions present in PDB but not in receptor PDBQT: %s", sorted(missing))

def log_possible_metal_mislabels(pdb_path: Union[str, Path]) -> None:
    """
    Heuristic: flag HET residues whose resname is not an element token, but
    whose element column shows an element token and the residue has very few atoms.
    """
    suspects = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        by_res = {}
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                key = (ln[21], ln[22:26], ln[26], ln[17:20].strip().upper())
                by_res.setdefault(key, []).append(ln)

    for (chain, resi, icode, resname), atms in by_res.items():
        if len(atms) > 4:
            continue
        if _is_element_token(resname):
            continue
        elem_tokens = {ln[76:78].strip().upper() for ln in atms if len(ln) >= 78}
        if any(_is_element_token(e) for e in elem_tokens):
            suspects.append((resname, f"{chain}:{resi}", sorted(elem_tokens)))

    if suspects:
        logging.warning("Possible metal mislabels (check resname vs element cols): %s", suspects)


# =============================
# Structural cleaning utilities
# =============================

def filter_altlocs(pdb_input_path: Union[str, Path], pdb_output_path: Union[str, Path]) -> None:
    """Filter alternate locations (altLoc) deterministically and log decisions."""
    lines: List[str] = []
    atoms = {}
    removed_count = 0

    with open(pdb_input_path, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                atom_name = line[12:16]
                altLoc = line[16]
                chainID = line[21]
                resSeq = line[22:26].strip()
                iCode = line[26]
                key = (chainID, resSeq, iCode, atom_name.strip())
                atoms.setdefault(key, {})[altLoc] = line
            else:
                lines.append(line)

    filtered_atoms: List[str] = []
    for key, altloc_dict in atoms.items():
        altLocs = list(altloc_dict.keys())
        if len(altLocs) > 1:
            logging.info("AltLocs found for %s: %s", key, altLocs)
        if ' ' in altloc_dict:
            selected = ' '
        elif 'A' in altloc_dict:
            selected = 'A'
        else:
            selected = sorted(altloc_dict.keys())[0]
        if len(altloc_dict) > 1:
            removed = [alt for alt in altLocs if alt != selected]
            removed_count += len(removed)
            logging.info("Keeping altLoc '%s' for atom %s, removed %s", selected, key, removed)
        filtered_atoms.append(altloc_dict[selected])

    def _int_safe(s: str) -> int:
        s = s.strip()
        return int(s) if s and s.lstrip("-").isdigit() else 0

    filtered_atoms.sort(key=lambda l: (l[21], _int_safe(l[22:26]), l[26], l[12:16].strip()))

    with open(pdb_output_path, 'w', encoding="utf-8") as f:
        for line in lines:
            f.write(line)
        for atom_line in filtered_atoms:
            f.write(atom_line)
    _post_write_element_guard("altloc_filter", pdb_output_path)
    logging.info("Filtered altLocs in %s → %s (removed %d alternates)", pdb_input_path, pdb_output_path, removed_count)


def build_missing_loops(input_pdb: Union[str, Path], output_dir: Union[str, Path]) -> str:
    """Fill missing loops/residues using MODELLER; return output PDB path."""
    try:
        from modeller import environ, log
        from modeller.scripts import complete_pdb
    except Exception as e:
        logging.warning("MODELLER not available; skipping loop completion: %s", e)
        return str(input_pdb)

    output_pdb = os.path.join(str(output_dir), "modeller_filled.pdb")
    log.none()
    logging.info("Running MODELLER to complete missing parts of %s", input_pdb)

    try:
        env = environ()
        env.io.hetatm = True
        env.io.water = True
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        # MODELLER loop building
        mdl = complete_pdb(env, str(input_pdb))
        mdl.write(file=output_pdb)
        _post_write_element_guard("MODELLER", output_pdb)

        if os.path.exists(output_pdb):
            logging.info("MODELLER filled PDB saved to %s", output_pdb)
            return output_pdb
        logging.warning("MODELLER did not produce expected output: %s", output_pdb)
        return str(input_pdb)
    except Exception as e:
        logging.warning("MODELLER failed on %s: %s", input_pdb, e)
        return str(input_pdb)


def filter_invalid_chains(pdb_path: Union[str, Path], output_path: Union[str, Path]) -> None:
    """Remove entire chains lacking backbone atoms (CA/N/C/O)."""
    chains: Dict[str, List[str]] = defaultdict(list)
    valid_chains: Set[str] = set()

    with open(pdb_path, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                chain_id = line[21]
                atom_name = line[12:16].strip()
                chains[chain_id].append(line)
                if atom_name in {"CA", "N", "C", "O"}:
                    valid_chains.add(chain_id)
            else:
                chains["HEADER"].append(line)

    with open(output_path, 'w', encoding="utf-8") as f:
        for chain_id in chains:
            if chain_id == "HEADER" or chain_id in valid_chains:
                f.writelines(chains[chain_id])
            else:
                logging.warning("Skipping invalid chain '%s' (no CA atoms)", chain_id)

# --- Hydrogen cleanup (geometry + CONECT) ---

def _parse_xyz(line: str):
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except Exception:
        return None


def conect_coverage(pdb_path: Union[str, Path]) -> float:
    atom_ids, conect_ids = set(), set()
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_ids.add(line[6:11].strip())
            elif line.startswith("CONECT"):
                parts = line.split()
                conect_ids.update(parts[1:])
    return len(conect_ids & atom_ids) / max(len(atom_ids), 1)


def remove_unbonded_atoms(pdb_path: Union[str, Path]) -> None:
    bonded_atoms: Set[str] = set()
    all_atoms: List[str] = []
    with open(pdb_path, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("CONECT"):
                parts = line.split()
                for atom_serial in parts[1:]:
                    bonded_atoms.add(atom_serial)
            elif line.startswith(("ATOM", "HETATM")):
                all_atoms.append(line)

    filtered: List[str] = []
    for line in all_atoms:
        atom_serial = line[6:11].strip()
        if atom_serial in bonded_atoms or line[76:78].strip() != 'H':
            filtered.append(line)
        else:
            logging.info("Removed unbonded hydrogen: %s", line.strip())

    with open(pdb_path, 'w', encoding="utf-8") as f:
        f.writelines(filtered)


def remove_implausible_hydrogens_by_distance(pdb_path: Union[str, Path]) -> None:
    atoms: List[Tuple[str, str, Optional[Tuple[float, float, float]]]] = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    atoms.append((line, "UNK", (None, None, None)))
                    continue
                el = line[76:78].strip() or line[12:16].strip()[:1]
                atoms.append((line, el.upper(), (x, y, z)))

    kept: List[str] = []
    heavy_coords = [a[2] for a in atoms if a[1] != "H" and a[2][0] is not None]

    def near_heavy(coord: Optional[Tuple[float, float, float]]) -> bool:
        if coord[0] is None:
            return True
        x, y, z = coord
        for X, Y, Z in heavy_coords:
            dx = x - X; dy = y - Y; dz = z - Z
            if (dx*dx + dy*dy + dz*dz) <= (1.35 * 1.35):
                return True
        return False

    for line, el, coord in atoms:
        if el != "H" or near_heavy(coord):
            kept.append(line)
        else:
            logging.info("Removed implausible H: %s", line.strip())

    with open(pdb_path, "w", encoding="utf-8") as out:
        out.writelines(kept)


def clean_hydrogens(pdb_path: Union[str, Path], use_conect_if_reliable: bool = True, conect_min_cov: float = 0.6) -> None:
    cov = conect_coverage(pdb_path) if use_conect_if_reliable else 0.0
    if cov >= conect_min_cov:
        remove_unbonded_atoms(pdb_path)
    remove_implausible_hydrogens_by_distance(pdb_path)
    _post_write_element_guard("hydrogen_cleanup", pdb_path)


# =============================
# Cofactors / waters policy
# =============================

def _is_metal(resname: str) -> bool:
    """Treat as 'metal/ion' iff YAML retain list contains this resname AND it looks like an element token."""
    rn = (resname or "").upper()
    return (rn in _RETAIN) and _is_element_token(rn)

def _cofactor_policy_keep(resname: str) -> bool:
    """
    YAML-only policy:
      • Keep anything listed under 'retain_in_receptor_resnames'
      • Keep standard residues (handled elsewhere)
      • Water handling is done in strip_nonstandard_residues(), so return False here for waters.
      • Everything else → remove
    """
    rn = (resname or "").upper()
    if rn in _RETAIN:
        return True
    if rn in _WATER_NAMES:
        return False
    return False


def strip_nonstandard_residues(input_pdb: Union[str, Path], output_pdb: Union[str, Path]) -> Tuple[int, str]:
    """
    Remove nonstandard residues while keeping what the YAML says to keep
    (retain_in_receptor_resnames) and standard amino acids.
    Water handling still honors your numeric policy (radius/B-factor) but
    the water names themselves come from the YAML.
    """
    standard_residues = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
        "HID","HIE","HIP","SEC","PYL","MSE",
    }

    removed: Set[str] = set()
    kept_lines: List[str] = []

    # Estimate pocket center from any YAML-retained cofactors/metals present
    cofm_xyz: List[Tuple[float, float, float]] = []
    with open(input_pdb, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip().upper()
            if resname in _RETAIN:
                xyz = _parse_xyz(line)
                if xyz: cofm_xyz.append(xyz)

    pocket_center: Optional[Tuple[float, float, float]] = None
    if cofm_xyz:
        from statistics import fmean
        xs, ys, zs = zip(*cofm_xyz)
        pocket_center = (fmean(xs), fmean(ys), fmean(zs))

    water_policy = (config.get("WATER_POLICY", "site_only") or "site_only").lower()
    water_radius = float(config.get("WATER_SITE_RADIUS_ANG", 6.0))
    water_bmax   = float(config.get("WATER_MAX_BFACTOR", 60.0))

    with open(input_pdb, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM  "):
                resname = line[17:20].strip().upper()
                if resname in standard_residues:
                    kept_lines.append(line)
                else:
                    removed.add(resname)
                continue

            if line.startswith("HETATM"):
                resname = line[17:20].strip().upper()

                # Water names come from YAML (no hardcoded list here)
                if resname in _WATER_NAMES:
                    if water_policy == "keep_all":
                        kept_lines.append(line); continue
                    if water_policy == "remove_all":
                        removed.add(resname);      continue
                    # site_only / auto: keep if near pocket center and not too mobile
                    if pocket_center is not None:
                        xyz = _parse_xyz(line)
                        keep = False
                        if xyz:
                            dx = xyz[0]-pocket_center[0]; dy = xyz[1]-pocket_center[1]; dz = xyz[2]-pocket_center[2]
                            dist2 = dx*dx + dy*dy + dz*dz
                            if dist2 <= water_radius*water_radius:
                                try:
                                    b = float(line[60:66]); keep = (b <= water_bmax)
                                except Exception:
                                    keep = True
                        if keep:
                            kept_lines.append(line)
                        else:
                            removed.add(resname)
                    else:
                        # no pocket estimate: default to remove unless keep_all
                        removed.add(resname)
                    continue

                # Retain anything listed in YAML retain block (cofactors, metals, ions, etc.)
                if resname in _RETAIN:
                    kept_lines.append(line)
                else:
                    removed.add(resname)
                continue

            # non-coordinate records pass through
            kept_lines.append(line)

    with open(output_pdb, 'w', encoding="utf-8") as f:
        f.writelines(kept_lines)
    _post_write_element_guard("strip_nonstandard", output_pdb)

    logging.info("Removed nonstandard residues (YAML-driven): %s", sorted(removed))
    if pocket_center:
        logging.info("Estimated pocket center: (%.2f, %.2f, %.2f)", *pocket_center)
    return len(removed), str(output_pdb)


def quick_element_histogram(pdb_path: Union[str, Path]) -> None:
    cnt = Counter()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM","HETATM")):
                el = ln[76:78].strip().upper()
                cnt[el or ""] += 1
    logging.info("[Elem histogram %s] %s", os.path.basename(str(pdb_path)), dict(sorted(cnt.items())))

def assert_no_metal_in_peptidic(pdb_path: Union[str, Path]) -> None:
    # Build a peptide-like name set from YAML
    pep_name_lines = RULES["element_sets"].get("peptide_like_names", [])
    peptidey = set()
    for line in pep_name_lines:
        for tok in str(line).split(";"):
            t = tok.strip().upper()
            if t:
                peptidey.add(t)

    # Element tokens considered "ionic" from YAML context (retain + element list)
    ionic_tokens = {r for r in _RETAIN if _is_element_token(r)}
    # If you want to exclude halides from the warning, do:
    # halides = set(RULES.get("element_sets", {}).get("halide_resnames", []))
    # ionic_tokens = ionic_tokens - halides

    res_atoms = defaultdict(list)
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM  ","HETATM")):
                key = (ln[21], ln[22:26], ln[26], ln[17:20].strip().upper())
                res_atoms[key].append(ln)

    offenders = []
    for key, lines in res_atoms.items():
        names = [ln[12:16].strip().upper() for ln in lines]
        if not names:
            continue
        hits = sum((n in peptidey) or (n[:2] in {"OE","NE","OD","ND","SD"}) for n in names)
        pep_like = hits >= max(4, 0.6*len(names))
        if not pep_like:
            continue
        # if any atom's element is an ionic token (from YAML), flag
        if any((ln[76:78].strip().upper() in ionic_tokens) for ln in lines):
            offenders.append(key)

    if offenders:
        logging.warning("Peptide-like residues contain ionic elements (check labeling): %s", offenders)


# ============================
# End-to-end Cleaning Pipeline
# =============================
def clean_pdb(pdb_file: Union[str, Path], output_root: Union[str, Path]) -> Optional[str]:
    """Run the full cleaning pipeline and return path to final cleaned PDB (receptor)."""
    output_root = str(output_root)
    Path(output_root).mkdir(parents=True, exist_ok=True)

    pdb_id = os.path.splitext(os.path.basename(str(pdb_file)))[0].upper()
    paths = canon_paths(pdb_id, output_root)
    for d in ["protein_root", "raw", "work", "ligands_raw", "nolig", "receptor"]:
        paths[d].mkdir(parents=True, exist_ok=True)

    # (1) Working copy → raw/
    working_pdb = paths["raw"] / f"{pdb_id}_working.pdb"
    shutil.copyfile(str(pdb_file), working_pdb)
    _helium_postwrite_counter("copy_working", working_pdb)

    # (2) AltLoc filtering → raw/filtered.pdb
    filtered_pdb = paths["raw"] / f"{pdb_id}_filtered.pdb"
    filter_altlocs(working_pdb, filtered_pdb)

    # (2a) EARLY text-level element fix (YAML-driven), before any heavy tools
    try:
        fix_element_columns_in_file(filtered_pdb, filtered_pdb, rewrite_atoms=True)
        _helium_postwrite_counter("elemfix_filtered", filtered_pdb)
        logging.info("Early text-level element fix applied to %s", filtered_pdb)
    except Exception as e:
        logging.warning("Early text-level element fix skipped for %s: %s", filtered_pdb, e)

    # (3) Extract ligands now (controls live here), with YAML element repair per-file
    _ = extract_ligands_from_filtered(filtered_pdb, paths["ligands_raw"])

    # (4) Strip nonstandard from protein (policy aware) → work/stripped.pdb
    stripped_pdb = paths["work"] / f"{pdb_id}_stripped.pdb"
    removed_count, _out = strip_nonstandard_residues(filtered_pdb, stripped_pdb)
    logging.info("Removed %d nonstandard residue lines.", removed_count)

    # (5) Element fix → MODELLER → element fix again (PDB only)
    elemfix_pdb = paths["work"] / f"{pdb_id}_elemfix.pdb"
    fix_pdb_elements(stripped_pdb, elemfix_pdb)
    _helium_postwrite_counter("elemfix_before_modeller", elemfix_pdb)
    loop_fixed_pdb = build_missing_loops(elemfix_pdb, paths["work"])
    fix_pdb_elements(loop_fixed_pdb, loop_fixed_pdb)
    _helium_postwrite_counter("elemfix_after_modeller", loop_fixed_pdb)

    # (6) Optional external Phenix polish (non-fatal if missing)
    receptor_pdb = paths["receptor"] / f"{pdb_id}_cleaned.pdb"
    if run_phenix_pdbtools(loop_fixed_pdb, receptor_pdb, remove_waters=True):
        pass  # already written by phenix.pdbtools
    else:
        shutil.copyfile(loop_fixed_pdb, receptor_pdb)
    _helium_postwrite_counter("phenix_or_copy_receptor", receptor_pdb)

    # (7) Hydrogen cleanup & chain validation
    debulked_pdb = paths["work"] / f"{pdb_id}_debulked.pdb"
    shutil.copyfile(receptor_pdb, debulked_pdb)
    clean_hydrogens(debulked_pdb, use_conect_if_reliable=True, conect_min_cov=0.6)

    chain_validated_pdb = paths["work"] / f"{pdb_id}_validated.pdb"
    filter_invalid_chains(debulked_pdb, chain_validated_pdb)
    _helium_postwrite_counter("chain_validate", chain_validated_pdb)

    try:
        _txt_before = Path(chain_validated_pdb).read_text(encoding="utf-8", errors="ignore")
    except Exception:
        _txt_before = ""

    try:
        # Rewrite PDB element columns (77–78) using the unified rules (also for ATOM when rewrite_atoms=True)
        fix_element_columns_in_file(chain_validated_pdb, chain_validated_pdb, rewrite_atoms=True)

        # Secondary invariant: if any H-named atom still carries He, fix and summarize
        _txt_after = Path(chain_validated_pdb).read_text(encoding="utf-8", errors="ignore")
        _fixed_text, _nname = assert_no_helium_in_hydrogen_names(_txt_after)
        if _nname > 0:
            Path(chain_validated_pdb).write_text(_fixed_text, encoding="utf-8")

        # Grep-friendly one-liner with He→H delta
        _before = scan_helium_counts(_txt_before)
        _after = scan_helium_counts(Path(chain_validated_pdb).read_text(encoding="utf-8", errors="ignore"))
        _delta = max(0, _before - _after)
        logging.info(f"[elem-fix] file={Path(chain_validated_pdb).name} stage=preflight He->H={_delta}")
    except Exception as _e:
        logging.warning(f"[elements] receptor preflight failed for {Path(chain_validated_pdb).name}: {_e}")
        
    # (8) Protonation (Reduce when safe; else Open Babel fallback)
    def _present_resnames(pdb_path: Path) -> set[str]:
        res = set()
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if ln.startswith(("ATOM  ", "HETATM")):
                    res.add(ln[17:20].strip().upper())
        return res

    present_resnames = _present_resnames(chain_validated_pdb)
    # Cofactors where Reduce tends to misinterpret atom names as elements (phosphorus variants etc.)
    nucleotide_like = {
        "NAP","NAD","NADP","NMN","FAD","FMN",
        "ATP","ADP","AMP","GTP","GDP","CTP","UTP","TTP"
    }
    use_reduce = not any(r in nucleotide_like for r in present_resnames)

    reduced_pdb = paths["work"] / f"{pdb_id}_reduced.pdb"
    assign_protonation_states(
        chain_validated_pdb,
        reduced_pdb,
        reduce_exe=REDUCE_EXE if use_reduce else None,  # skip Reduce for nucleotide cofactors
    )
    _helium_postwrite_counter("reduce_or_fallback", reduced_pdb)

    # (9) Final element fix and sanity on the protonated file
    fix_pdb_elements(reduced_pdb)
    _helium_postwrite_counter("elemfix_after_reduce", reduced_pdb)

    quick_element_histogram(reduced_pdb)

    assert_no_metal_in_peptidic(reduced_pdb)

    # (10) Move to receptor and re-fix (post-step edits)
    shutil.copyfile(reduced_pdb, receptor_pdb)
    _helium_postwrite_counter("promote_receptor_copy", receptor_pdb)

    fix_pdb_elements(receptor_pdb)
    _helium_postwrite_counter("elemfix_final_receptor", receptor_pdb)

    quick_element_histogram(receptor_pdb)
    assert file_contains_hydrogens(receptor_pdb), f"[FATAL] Cleaned file lost hydrogens: {receptor_pdb}"

    logging.info("Cleaned receptor: %s", receptor_pdb)
    return str(receptor_pdb)


# --- Histidine helpers used by Meeko wrapper ---
from collections import defaultdict

def _classify_and_rename_histidines(pdb_in: Union[str, Path], pdb_out: Union[str, Path]) -> None:
    """
    Inspect each HIS residue's side-chain hydrogens and rename:
      - HD1 only  -> HID
      - HE2 only  -> HIE
      - both      -> HIP
      - neither   -> keep HIS (caller may rewrite later)
    """
    pdb_in, pdb_out = str(pdb_in), str(pdb_out)

    residues = defaultdict(list)
    other_lines = []
    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                chain = ln[21]
                resseq = ln[22:26]
                icode  = ln[26]
                resn   = ln[17:20]
                key = (chain, resseq, icode, resn)
                residues[key].append(ln)
            else:
                other_lines.append(ln)

    out_lines = []
    for key, atoms in residues.items():
        chain, resseq, icode, resn = key
        resname = resn.strip()
        if resname != "HIS":
            out_lines.extend(atoms)
            continue

        names = {ln[12:16].strip().upper() for ln in atoms}
        has_hd1 = "HD1" in names  # proton on ND1
        has_he2 = "HE2" in names  # proton on NE2

        if has_hd1 and has_he2:
            new = "HIP"
        elif has_hd1 and not has_he2:
            new = "HID"
        elif has_he2 and not has_hd1:
            new = "HIE"
        else:
            new = None

        if new:
            for ln in atoms:
                if ln.startswith(("ATOM  ", "HETATM")):
                    out_lines.append(ln[:17] + f"{new:>3}" + ln[20:])
                else:
                    out_lines.append(ln)
        else:
            out_lines.extend(atoms)

    with open(pdb_out, "w", encoding="utf-8") as w:
        for ln in other_lines:
            w.write(ln)
        for ln in out_lines:
            w.write(ln)


def _rewrite_his_default(pdb_in: Union[str, Path], pdb_out: Union[str, Path], default: str = "HIE") -> None:
    pdb_in, pdb_out = str(pdb_in), str(pdb_out)
    default = default.upper()
    assert default in {"HIE","HID","HIP"}
    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as f, \
         open(pdb_out, "w", encoding="utf-8") as w:
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")) and ln[17:20] == "HIS":
                w.write(ln[:17] + default + ln[20:])
            else:
                w.write(ln)

# =============================
# External Tools (Meeko/ADT, Phenix, OpenBabel, Reduce)
# =============================
def run_prepare_receptor(input_pdb: Union[str, Path], output_pdbqt: Union[str, Path], cfg: dict) -> bool:
    """
    Robust Meeko/ADT wrapper with:
      â€¢ Early check for HIS tautomer ambiguity on the modern (--read_pdb) call
      â€¢ Auto -n mapping for all truly ambiguous HIS (no HD1/HE2)
      â€¢ Optional -a (allow_bad_res) retry on template-mismatch
      â€¢ Heavy-handed HIS?<default> fallback
      â€¢ ADT fallback
    """
    import re
    from tempfile import NamedTemporaryFile

    input_pdb = str(input_pdb)
    output_pdbqt = str(output_pdbqt)
    # --- Helium preflight on the EXACT file we’re about to feed into Meeko pipeline
    # Work dir is the canonical processed tree root/<PDB>/work (output_pdbqt = .../receptor/<PDB>.pdbqt)
    _work_dir = Path(output_pdbqt).resolve().parent.parent / "work"
    try:
        # This both fixes element columns last-mile and FAILS FAST with offenders if He>0.
        input_pdb = str(_meeko_preflight_or_fail(input_pdb, _work_dir))
    except RuntimeError as e:
        if str(e) == "helium_preflight_failed":
            logging.error("Abort: helium detected pre-Meeko; see [helium] offenders above.")
            raise
        raise

    his_default = str(cfg.get("HIS_DEFAULT", "HIE")).upper()
    if his_default not in {"HIE", "HID", "HIP"}:
        his_default = "HIE"

    # improved truthiness parsing
    allow_bad_res = str(cfg.get("MEEKO_ALLOW_BAD_RES", "true")).lower() in ("1", "true", "yes")
    default_altloc = (cfg.get("MEEKO_DEFAULT_ALTLOC") or "").strip()  # e.g. "A" or ""

    def _meeko_cmd() -> list[str]:
        # 1) PATH shims (with/without .py)
        for name in ("mk_prepare_receptor", "mk_prepare_receptor.py"):
            exe = shutil.which(name)
            if exe:
                return [exe]
        # 2) Scripts directory of the current interpreter
        try:
            import sysconfig, os
            scripts = sysconfig.get_path("scripts")
            cand = os.path.join(scripts, "mk_prepare_receptor.py")
            if os.path.exists(cand):
                return [sys.executable, cand]  # run via python to be safe
        except Exception:
            pass
        raise FileNotFoundError(
            "Meeko CLI not found. Try adding the Python 'scripts' dir to PATH or call mk_prepare_receptor.py directly."
        )
    def _run(cmd: List[str]) -> subprocess.CompletedProcess:
        logging.info("Meeko: %s", " ".join(cmd))
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        (logging.info if r.returncode == 0 else logging.warning)(
            "rc=%s\nSTDOUT:\n%s\nSTDERR:\n%s", r.returncode, (r.stdout or ""), (r.stderr or "")
        )
        return r

    def _build_his_override_mapping(pdb_path: str) -> str:
        """Return a single -n mapping like 'A:27,A:94=HIE' for HIS with no HD1/HE2."""
        from collections import defaultdict
        by_res = defaultdict(set)
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                if ln[17:20] != "HIS":
                    continue
                chain = ln[21]
                resi = ln[22:26].strip()
                aname = ln[12:16].strip().upper()
                by_res[(chain, resi)].add(aname)
        targets = []
        for (chain, resi), names in by_res.items():
            if ("HD1" not in names) and ("HE2" not in names):
                try:
                    rnum = int(resi)
                except Exception:
                    rnum = int(resi.strip() or "0")
                targets.append(f"{chain}:{rnum}")
        return (",".join(sorted(targets)) + f"={his_default}") if targets else ""

    def _modern_meeko(pdb_path: str, extra_flags: Optional[List[str]] = None) -> subprocess.CompletedProcess:
        try:
            cmd = _meeko_cmd() + ["--read_pdb", pdb_path, "-p", output_pdbqt]
        except FileNotFoundError as e:
            from types import SimpleNamespace
            return SimpleNamespace(returncode=127, stdout="", stderr=str(e))
        flags = extra_flags or []
        return _run(cmd + flags)

    def _legacy_meeko(pdb_path: str, extra_flags: Optional[List[str]] = None) -> bool:
        try:
            base = _meeko_cmd()
        except FileNotFoundError:
            return False
        flags = extra_flags or []
        tried = False
        r = _run(base + ["-i", pdb_path, "-p", output_pdbqt] + flags); tried = True
        if r.returncode == 0 and Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size > 0:
            return True
        r = _run(base + ["-r", pdb_path, "-o", output_pdbqt] + flags)
        return r.returncode == 0 and Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size > 0

    # --- Step 0: classify HIS by existing hydrogens (best-case, no coord change)
    with NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp1:
        tmp1_path = tmp1.name
    try:
        _classify_and_rename_histidines(input_pdb, tmp1_path)  # leaves non-diagnostic HIS as HIS
    except Exception as e:
        logging.warning("HIS classify/rename step failed (continuing with original): %s", e)
        tmp1_path = input_pdb

    # Preflight again on the HIS-classified file we’re actually handing to Meeko now
    try:
        _meeko_preflight_or_fail(tmp1_path, Path(output_pdbqt).resolve().parent.parent / "work")
    except RuntimeError as e:
        if str(e) == "helium_preflight_failed":
            logging.error("Abort: helium detected in HIS-classified PDB pre-Meeko; see [helium] offenders above.")
            raise
        raise

    # --- Step 1: Modern Meeko first
    r0 = _modern_meeko(tmp1_path)
    if r0.returncode == 0 and Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size > 0:
        logging.info("Receptor prepared (Meeko).")
        if tmp1_path != input_pdb:
            try: os.remove(tmp1_path)
            except Exception: pass
        return True

    se0 = (r0.stderr or "")
    # (1a) HIS tie? -> build a -n mapping and retry once
    his_tie = ("tied for fewest missing H" in se0) and ("HIE" in se0 and "HID" in se0)
    if his_tie:
        mapping = _build_his_override_mapping(tmp1_path)
        if not mapping:
            m = re.search(r"residue_key='([A-Za-z]):(\d+)'", se0)
            if m:
                mapping = f"{m.group(1)}:{int(m.group(2))}={his_default}"
        if mapping:
            logging.warning("Meeko histidine ambiguity -> retry with -n %s", mapping)
            r1 = _modern_meeko(tmp1_path, extra_flags=["-n", mapping])
            if r1.returncode == 0 and Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size > 0:
                logging.info("Receptor prepared after tautomer override (Meeko).")
                if tmp1_path != input_pdb:
                    try: os.remove(tmp1_path)
                    except Exception: pass
                return True

    # (1b) Template mismatch? -> retry with -a (and default altloc if provided)
    templ_fail = ("No template matched for residue_key" in se0) or ("Template matched failed" in se0) \
                 or ("Template matching failed" in se0)
    if allow_bad_res and (templ_fail or "allow_bad_res" in se0 or "recommendations" in se0.lower()):
        extra = ["-a"]
        if default_altloc:
            extra += ["--default_altloc", default_altloc]
        logging.warning("Template mismatch -> retry with %s", " ".join(extra))
        r2 = _modern_meeko(tmp1_path, extra_flags=extra)
        if r2.returncode == 0 and Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size > 0:
            logging.info("Receptor prepared after -a (allow_bad_res).")
            if tmp1_path != input_pdb:
                try: os.remove(tmp1_path)
                except Exception: pass
            return True

    # --- Step 2: Heavy-handed fallback: rewrite residual HIS -> default and retry modern Meeko (+-a)
    logging.warning("Retrying Meeko after rewriting residual HIS -> %s", his_default)
    with NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp2:
        tmp2_path = tmp2.name
    try:
        _rewrite_his_default(input_pdb, tmp2_path, default=his_default)
        extra = ["-a"] if allow_bad_res else []
        if default_altloc:
            extra += ["--default_altloc", default_altloc]
        r3 = _modern_meeko(tmp2_path, extra_flags=extra)
        if r3.returncode == 0 and Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size > 0:
            logging.info("Receptor prepared after HIS rewrite (Meeko).")
            return True
    finally:
        try: os.remove(tmp2_path)
        except Exception: pass

    # --- Step 3: Legacy variants
    extra = ["-a"] if allow_bad_res else []
    if default_altloc:
        extra += ["--default_altloc", default_altloc]
    if _legacy_meeko(tmp1_path, extra_flags=extra):
        logging.info("Receptor prepared (Meeko legacy).")
        if tmp1_path != input_pdb:
            try: os.remove(tmp1_path)
            except Exception: pass
        return True

    if tmp1_path != input_pdb:
        try: os.remove(tmp1_path)
        except Exception: pass

    # --- Step 4: ADT fallback
    mgltools_python = cfg.get("MGLTOOLS_PYTHON")
    prepare_script   = cfg.get("PREPARE_RECEPTOR_SCRIPT")
    if not mgltools_python or not prepare_script or not os.path.exists(prepare_script):
        logging.error("ADT receptor prep unavailable (MGLTOOLS_PYTHON or PREPARE_RECEPTOR_SCRIPT missing/not found)")
        return False

    cmd = [
        mgltools_python, prepare_script,
        "-r", input_pdb, "-o", output_pdbqt,
        "-A", "none", "-U", "nphs_lps_nonstdres"
    ]
    logging.info("Running prepare_receptor4.py: %s", " ".join(cmd))
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        logging.error("prepare_receptor4 failed:\n%s", r.stderr or r.stdout or "")
        return False
    logging.info("Receptor prepared (ADT).")
    return True

def _phenix_detect() -> tuple[str, list[str] | None, dict]:
    """
    Return (mode, cmd_list, env) for phenix.pdbtools detection.
    mode ∈ {'linux-path','linux-explicit','wsl-bat','unavailable'}
    """
    env = dict(os.environ)
    if PHENIX_LIB_PATH:
        env["PYTHONPATH"] = PHENIX_LIB_PATH

    # 1) On PATH (native Linux installs)
    if shutil.which("phenix.pdbtools"):
        return ("linux-path", ["phenix.pdbtools"], env)

    # 2) Explicit Linux path via PHENIX_DIR (your case)
    if PHENIX_DIR:
        explicit = Path(PHENIX_DIR) / "phenix.pdbtools"
        if explicit.exists() and os.access(str(explicit), os.X_OK):
            return ("linux-explicit", [str(explicit)], env)

    # 3) Windows .bat (usable from WSL via PowerShell)
    if PDBTOOLS_BAT and Path(PDBTOOLS_BAT).exists():
        return ("wsl-bat", [PDBTOOLS_BAT], env)

    return ("unavailable", None, env)


def _run_and_log(cmd: list[str], *, env: dict, check: bool = False) -> subprocess.CompletedProcess:
    logging.info("[phenix] exec: %s", " ".join(cmd))
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if cp.stdout:
        logging.debug("[phenix][stdout]\n%s", cp.stdout)
    if cp.stderr:
        logging.debug("[phenix][stderr]\n%s", cp.stderr)
    if check and cp.returncode != 0:
        raise subprocess.CalledProcessError(cp.returncode, cmd, cp.stdout, cp.stderr)
    return cp


def run_phenix_pdbtools(input_pdb: Union[str, Path],
                        output_pdb: Union[str, Path],
                        remove_waters: bool = True) -> bool:
    """
    Run Phenix pdbtools if available (Linux, or Windows via PowerShell under WSL).
    Return True if successful. Always log branch as: [phenix] mode=...
    """
    inp  = os.path.abspath(str(input_pdb)).replace("\\", "/")
    outp = os.path.abspath(str(output_pdb)).replace("\\", "/")

    mode, base_cmd, env = _phenix_detect()
    logging.info("[phenix] mode=%s", mode)

    if mode == "unavailable" or not base_cmd:
        logging.warning("Phenix not available; skipping pdbtools polish")
        return False

    cmd = base_cmd + [inp, f"output.file_name={outp}"]
    if remove_waters:
        cmd.append('remove="resname HOH"')

    try:
        cp = _run_and_log(cmd, env=env, check=True)
        if os.path.exists(outp) and os.path.getsize(outp) > 0:
            logging.info("phenix.pdbtools wrote %s", outp)
            return True
        logging.warning("phenix.pdbtools finished but did not create expected output: %s", outp)
        return False
    except Exception as e:
        logging.error("phenix.pdbtools failed: %s", e)
        # Diagnostic: phenix.python import iotbx (same env)
        py_cmd = None
        if PHENIX_DIR and (Path(PHENIX_DIR) / "phenix.python").exists():
            py_cmd = [str(Path(PHENIX_DIR) / "phenix.python")]
        elif shutil.which("phenix.python"):
            py_cmd = ["phenix.python"]
        elif PHENIX_PYTHON_BAT and Path(PHENIX_PYTHON_BAT).exists():
            py_cmd = [PHENIX_PYTHON_BAT]
        if py_cmd:
            try:
                _run_and_log(py_cmd + ["-c", "from iotbx import pdb; print('iotbx_ok')"], env=env, check=False)
            except Exception:
                pass
        return False


def run_windows_phenix_clean_script(loop_fixed_pdb: Union[str, Path], nolig_dir: Union[str, Path]) -> int:
    """Call your PHENIX_CLEAN_SCRIPT using Windows Phenix Python from WSL. Returns returncode."""
    if not PHENIX_PYTHON_BAT or not Path(PHENIX_PYTHON_BAT).exists() or not PHENIX_CLEAN_SCRIPT:
        return 127
    win_script = _win_path(PHENIX_CLEAN_SCRIPT)
    win_in = _win_path(loop_fixed_pdb)
    win_outdir = _win_path(nolig_dir)
    ps = (
        "$ErrorActionPreference='Stop';"
        f"if (!(Test-Path '{win_outdir}')) {{ New-Item -ItemType Directory -Force -Path '{win_outdir}' | Out-Null }};"
        f"Push-Location '{win_outdir}';"
        f"& '{_win_path(PHENIX_PYTHON_BAT)}' '{win_script}' '{win_in}' '{win_outdir}';"
        "Pop-Location;"
    )
    r = _powershell(ps)
    return r.returncode

def run_openbabel_add_h(input_pdb: Union[str, Path], output_pdb: Union[str, Path]) -> None:
    obabel = OPENBABEL_PATH or shutil.which("obabel")
    if not obabel or not shutil.which(Path(obabel).name):
        raise RuntimeError("Open Babel not found. Install or set OPENBABEL_PATH.")
    cmd = [obabel, str(input_pdb), "-O", str(output_pdb), "-h"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Open Babel failed:\n{r.stderr}")
    logging.info("Open Babel hydrogenation succeeded")
from typing import Union, Optional
from pathlib import Path
import shutil, subprocess, logging, os


def assign_protonation_states(input_pdb: Union[str, Path],
                              output_pdb: Union[str, Path],
                              reduce_exe: Optional[str] = None) -> str:
    import uuid

    input_pdb  = os.path.abspath(str(input_pdb)).replace("\\", "/")
    output_pdb = os.path.abspath(str(output_pdb)).replace("\\", "/")

    has_h = hydrogenation_status(input_pdb)[0] != "NO_H"
    reduce_flags = ["-quiet"] if has_h else ["-BUILD", "-quiet"]
    exe = reduce_exe or REDUCE_EXE
    exe_dir = os.path.dirname(exe) or None

    def run_reduce(in_pdb: str, stage_name: str) -> str:
        env = dict(os.environ)
        # Allow config-provided dict even if user didn't export it
        het = env.get("REDUCE_HET_DICT") or _cfg("REDUCE_HET_DICT", "")
        if het:
            env["REDUCE_HET_DICT"] = het

        with open(output_pdb, "w", encoding="utf-8") as out:
            cp = subprocess.run([exe] + reduce_flags + [in_pdb],
                                stdout=out, stderr=subprocess.PIPE, text=True,
                                cwd=exe_dir, env=env)
        if cp.returncode != 0:
            raise RuntimeError(f"{stage_name} reduce failed: {cp.stderr.strip()}")
        return output_pdb

    status_before, h0, hv0, r0 = hydrogenation_status(input_pdb)
    logging.info("[H-Scan before] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_before, h0, hv0, r0)

    try:
        run_reduce(input_pdb, "Reduce#1")
    except Exception as e:
        logging.warning("Reduce#1 failed: %s", e)
        try:
            # Repair element columns (YAML-driven fixer under the hood)
            fix_pdb_elements(input_pdb)
        except Exception:
            pass
        try:
            tmp_in = os.path.splitext(input_pdb)[0] + f"_retry_{uuid.uuid4().hex}.pdb"
            shutil.copy(input_pdb, tmp_in)
            try:
                run_reduce(tmp_in, "Reduce#2(temp)")
            finally:
                try:
                    os.remove(tmp_in)
                except Exception:
                    pass
        except Exception as e2:
            logging.error("Reduce temp retry failed: %s", e2)
            try:
                if has_h:
                    shutil.copy(input_pdb, output_pdb)
                    logging.warning("Reduce failed; keeping existing hydrogens (no rebuild).")
                else:
                    run_openbabel_add_h(input_pdb, output_pdb)
                    logging.info("Open Babel used to add hydrogens.")
            except Exception as babel_error:
                logging.error("OpenBabel fallback failed: %s", babel_error)
                shutil.copy(input_pdb, output_pdb)
                logging.warning("Hydrogenation skipped; copied input to output.")

    # everything below stays INSIDE the function
    status_after, h1, hv1, r1 = hydrogenation_status(output_pdb)
    logging.info("[H-Scan after ] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_after, h1, hv1, r1)

    if status_after == "NO_H":
        logging.warning("Reduce produced no hydrogens; trying Open Babel fallback.")
        try:
            tmp_babel = output_pdb + ".babel.pdb"
            run_openbabel_add_h(input_pdb, tmp_babel)
            shutil.move(tmp_babel, output_pdb)
            status_after, h1, hv1, r1 = hydrogenation_status(output_pdb)
            logging.info("[H-Scan babel] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_after, h1, hv1, r1)
        except Exception as e:
            logging.error("OpenBabel fallback failed after Reduce: %s", e)

        if hydrogenation_status(output_pdb)[0] == "NO_H":
            if str(os.environ.get("ALLOW_NO_HYDROGENS", "")).lower() in ("1","true","yes"):
                logging.warning("Continuing with NO_H due to ALLOW_NO_HYDROGENS env override.")
                return output_pdb
            raise RuntimeError("Protonation produced no hydrogens.")

    return output_pdb



def run_molprobity_validate(pdb_path: Union[str, Path], work_dir: Optional[Union[str, Path]] = None) -> int:
    """
    Try running MolProbity validation via 'phenix.molprobity' if available.
    Returns the process return code (0 means success). Non-fatal if missing.
    """
    exe = shutil.which("phenix.molprobity")
    if not exe:
        logging.info("MolProbity not available; skipping validation.")
        return 127
    work = Path(work_dir or Path(pdb_path).with_suffix(".molprobity")).resolve()
    work.mkdir(parents=True, exist_ok=True)
    cmd = [exe, str(Path(pdb_path).resolve())]
    logging.info("Running MolProbity: %s", " ".join(cmd))
    r = subprocess.run(cmd, cwd=str(work), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        logging.warning("MolProbity failed (rc=%s)\nSTDERR:\n%s", r.returncode, (r.stderr or ""))
    else:
        logging.info("MolProbity finished; results under %s", work)
    return r.returncode
# =============================
# Module Entrypoint
# =============================

def main(pdb_filename: str, output_dir: Union[str, Path] = r"./processed_pdbs"):
    """High-level wrapper: clean a PDB and prepare the receptor PDBQT."""
    try:
        # Resolve input path (abs path takes precedence)
        if os.path.isabs(pdb_filename) and os.path.isfile(pdb_filename):
            pdb_path = pdb_filename
        else:
            pdb_path = os.path.join(_cfg("INPUT_DIR", "."), pdb_filename)

        if not os.path.isfile(pdb_path):
            logging.error("ERROR: File does not exist: %s", pdb_path)
            return None

        # Run cleaning → returns the final cleaned receptor PDB path
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        cleaned_pdb = clean_pdb(pdb_path, output_dir)
        if not cleaned_pdb:
            logging.error("ERROR: Cleaning failed for %s", pdb_filename)
            return None

        # Prepare receptor PDBQT next to the cleaned tree
        pdb_id = Path(pdb_filename).stem.upper()
        paths = canon_paths(pdb_id, output_dir)
        output_pdbqt = str((paths["receptor"] / f"{pdb_id}.pdbqt").resolve())
        fix_element_columns_in_file(cleaned_pdb, cleaned_pdb, rewrite_atoms=True)

        if not run_prepare_receptor(cleaned_pdb, output_pdbqt, config):
            logging.error("ERROR: Failed to prepare receptor PDBQT for %s", pdb_id)
            return None

        logging.info("Prepared receptor PDBQT: %s", output_pdbqt)
        return cleaned_pdb, output_pdbqt

    except Exception as e:
        logging.exception("[FATAL] automate_protein_prep.main() failed: %s", e)
        return None

# Example usage
# if __name__ == "__main__":
#     setup_logger()
#     test_pdb = "1a3n.pdb"
#     result = main(test_pdb)
#     if result:
#         print("Final cleaned PDB:", result[0])
