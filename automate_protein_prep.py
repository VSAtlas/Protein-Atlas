"""
Automated protein preparation pipeline (Windows/WSL-friendly)

This module prepares a receptor structure from a raw PDB by performing:
  1) Alternate-conformation filtering (altLoc) with deterministic selection
  2) Removal of nonstandard residues (policy-aware keep/drop of cofactors/metals/waters)
  3) Element-column fixes to ensure PDB compliance
  4) Optional loop/residue completion via MODELLER
  5) Sanity checks for incomplete residues
  6) Phenix cleaning passes (auto-detect Linux Phenix, Windows Phenix, or skip—WSL-safe)
  7) Hydrogen sanity cleanup (CONECT- and geometry-based)
  8) Chain validation (ensure at least one CA-containing chain)
  9) Protonation via Reduce with automatic fallbacks (temp retry/OpenBabel)
 10) Final polish via phenix.pdbtools when available; otherwise degrade gracefully
 11) Receptor PDBQT preparation (Meeko if available/selected; else ADT)

Design notes
------------
• WSL-aware: Phenix calls are auto-selected: Linux binaries, Windows .bat via PowerShell + wslpath, or skipped.
• Config keys are case-insensitive (PHENIX_LIB_PATH vs phenix_lib_path, etc.).
• Logging instead of prints; short, numbered steps for readability and easier debugging.
• No functionality intentionally removed; behavior is clarified and guarded.

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
import logging
import os
import platform
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from math import sqrt
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

# -------- Project-local imports (assumed available in your repo) --------
from installation import load_config
from activesite import fix_pdb_elements  # element column rewriter (trusted external impl)
from logger_setup import setup_logger  # assumed to be called by the entrypoint

# =============================
# Configuration helpers
# =============================
config = load_config()  # load once at import time (matches original behavior)

# Case-insensitive lookup with alias support
def _cfg(key: str, default=None, *aliases: str):
    if key in config:
        return config[key]
    for a in aliases:
        if a in config:
            return config[a]
    # case-insensitive fallbacks
    lk = key.lower()
    uk = key.upper()
    if lk in config:
        return config[lk]
    if uk in config:
        return config[uk]
    for a in aliases:
        if a.lower() in config:
            return config[a.lower()]
        if a.upper() in config:
            return config[a.upper()]
    return default

# --- Cofactor / metal / water policy defaults (override in your config) ---
config.setdefault("COFACTOR_POLICY", "auto")  # ["auto","keep","remove"]
config.setdefault(
    "COFACTOR_KEEP_LIST",
    "HEM,FAD,NAD,NADH,NADP,NADPH,FMN,PLP,TPP,BIO,SAH,SAM,COA,H4B,B6P,ATP,ADP",
)
config.setdefault(
    "COFACTOR_REMOVE_LIST",
    "GOL,EDO,PG4,MPD,ACT,TRS,SO4,PO4,PEG,CL,BR,NA,K,CA",
)
config.setdefault("KEEP_METALS", True)
config.setdefault("METAL_LIST", "ZN,MG,FE,MN,CU,NI,CO,MO,CD,NA,K,CA")
config.setdefault("WATER_POLICY", "site_only")  # ["auto","keep_all","remove_all","site_only"]
config.setdefault("WATER_SITE_RADIUS_ANG", 6.0)
config.setdefault("WATER_MAX_BFACTOR", 60.0)

# Fast lookup sets
_COFACTOR_KEEP = {s.strip().upper() for s in config.get("COFACTOR_KEEP_LIST", "").split(",") if s.strip()}
_COFACTOR_DROP = {s.strip().upper() for s in config.get("COFACTOR_REMOVE_LIST", "").split(",") if s.strip()}
_METALS = {s.strip().upper() for s in config.get("METAL_LIST", "").split(",") if s.strip()}

# Tool paths (robust to case and missing keys)
PHENIX_DIR = _cfg("PHENIX_DIR", "", "phenix_dir")
PHENIX_LIB_PATH = _cfg("PHENIX_LIB_PATH", "", "phenix_lib_path")
PHENIX_CLEAN_SCRIPT = _cfg("PHENIX_CLEAN_SCRIPT", "", "phenix_clean_script")
INPUT_DIR = _cfg("INPUT_DIR", ".")
OPENBABEL_PATH = _cfg("OPENBABEL_PATH", shutil.which("obabel") or "")
MGLTOOLS_PYTHON = _cfg("MGLTOOLS_PYTHON", "")
PREPARE_RECEPTOR_SCRIPT = _cfg("PREPARE_RECEPTOR_SCRIPT", "")
USE_MEEKO = str(_cfg("USE_MEEKO", "")).lower() in ("1", "true", "yes")

# Make Phenix Python importable if available
if PHENIX_LIB_PATH:
    sys.path.insert(0, str(PHENIX_LIB_PATH))
if PHENIX_DIR:
    sys.path.insert(0, str(PHENIX_DIR))

# Windows .bat wrappers (if using Windows Phenix)
PDBTOOLS_BAT = str(Path(PHENIX_DIR) / "phenix.pdbtools.bat") if PHENIX_DIR else ""
MOLPROBITY_BAT = str(Path(PHENIX_DIR) / "phenix.molprobity.bat") if PHENIX_DIR else ""
PHENIX_PYTHON_BAT = str(Path(PHENIX_DIR) / "phenix.python.bat") if PHENIX_DIR else ""

# Prefer explicit Reduce from config; else try PHENIX_DIR; else PATH
REDUCE_EXE = _cfg("REDUCE_EXE", "", "reduce_exe")
if not REDUCE_EXE:
    if PHENIX_DIR and Path(PHENIX_DIR, "reduce.exe").exists():
        REDUCE_EXE = str(Path(PHENIX_DIR, "reduce.exe"))
    else:
        REDUCE_EXE = shutil.which("reduce") or shutil.which("reduce.exe") or "reduce"

# Two-letter elements commonly encountered in PDBs (for element inference)
_TWO_LETTER = {
    "ZN","FE","MG","MN","CL","NA","CA","CU","CO","BR","SI","PT","PD","NI","AL",
    "AG","AU","IR","SR","BA","BE","LI","RB","CS","MO","SE","TE","TI","CR","VD",
    "HG","PB","SN","SB","CD","GA","GE","ZR","Y","NB","W","RE","OS","RH","RU","I"
}

# RDKit optional
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

def _win_path(p: str | Path) -> str:
    """Convert WSL path to Windows path (no-op outside WSL)."""
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

# =============================
# Ligand pristine reference (for RMSD downstream)
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

def canon_paths(pdb_id: str, output_root: str | Path) -> dict[str, Path]:
    root = Path(output_root).resolve()
    base = root / pdb_id.upper()
    return {
        "protein_root": base,                      # processed_pdbs/1IEP
        "raw":          base / "raw",              # working copy, altloc-filtered
        "work":         base / "work",             # intermediates
        "ligands_raw":  base / "ligands_raw",      # extracted ligands
        "nolig":        base / "nolig",            # ligand-stripped protein + phenix outputs
        "receptor":     base / "receptor",         # final cleaned receptor & PDBQT
    }

def fold_legacy_layout(pdb_id: str, output_root: str | Path) -> None:
    """Best-effort migration of legacy sibling dirs into the canonical tree (non-fatal)."""
    root = Path(output_root).resolve()
    base = root / pdb_id.upper()
    base.mkdir(parents=True, exist_ok=True)

    candidates = [
        (root / f"{pdb_id}_nolig", base / "nolig"),
        (root / f"{pdb_id.lower()}_nolig", base / "nolig"),
        (root / f"{pdb_id}_cleaned_ligands", base / "ligands_raw"),
        (root / f"{pdb_id.lower()}_cleaned_ligands", base / "ligands_raw"),
        (root / f"{pdb_id}_nolig.pdb", base / "nolig" / f"{pdb_id}_nolig_phenix_clean.pdb"),
        (root / f"{pdb_id.lower()}_nolig.pdb", base / "nolig" / f"{pdb_id}_nolig_phenix_clean.pdb"),
    ]
    for src, dst in candidates:
        try:
            if src.is_dir():
                dst.parent.mkdir(parents=True, exist_ok=True)
                for p in src.rglob("*"):
                    rel = p.relative_to(src)
                    target = dst / rel
                    if p.is_dir():
                        target.mkdir(parents=True, exist_ok=True)
                    else:
                        if not target.exists():
                            target.parent.mkdir(parents=True, exist_ok=True)
                            shutil.move(str(p), str(target))
                shutil.rmtree(src, ignore_errors=True)
            elif src.is_file():
                dst.parent.mkdir(parents=True, exist_ok=True)
                if not dst.exists():
                    shutil.move(str(src), str(dst))
        except Exception:
            pass

# =============================
# Ligand extraction
# =============================

def extract_ligands_from_filtered(filtered_pdb: str | Path, out_dir: str | Path) -> list[Path]:
    """Extract non-water HETATM residues into individual PDBs (RES_CHAINRESI.pdb)."""
    outd = Path(out_dir)
    outd.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    cur_key = None
    bucket: list[str] = []

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
            if resname == "HOH":
                continue
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
    logging.info("Extracted %d ligand residues to %s", len(written), outd)
    return written

# =============================
# Element & Hydrogen utilities
# =============================

def looks_like_hydrogen_name(line: str) -> bool:
    name = line[12:16]
    return name[0] == " " and name[1].upper() in {"H", "D", "T"} or name.lstrip()[:1].upper() in {"H","D","T"}

def infer_element(line: str) -> str:
    el = line[76:78].strip().upper()
    if el:
        return "H" if el in {"D", "T"} else el
    name = line[12:16]
    c1 = name[1] if name[0] == " " else name[0]
    c2 = (name[2] if name[0] == " " else name[1]).strip()
    cand2 = (c1 + c2).upper()
    return cand2 if cand2 in _TWO_LETTER else c1.upper()

def count_atoms_by_element(pdb_path: str | Path) -> Tuple[int, int]:
    h = heavy = 0
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            el = infer_element(line)
            if el == "H":
                # If element field blank and name doesn't look like H, treat as heavy
                if line[76:78].strip() == "" and not looks_like_hydrogen_name(line):
                    heavy += 1
                else:
                    h += 1
            else:
                heavy += 1
    return h, heavy

def hydrogenation_status(pdb_path: str | Path) -> Tuple[str, int, int, float]:
    h, heavy = count_atoms_by_element(pdb_path)
    ratio = h / max(heavy, 1)
    if h == 0:
        return "NO_H", h, heavy, ratio
    if ratio < 0.20:
        return "SUSPECT_LOW_H", h, heavy, ratio
    return "HAS_H", h, heavy, ratio

def file_contains_hydrogens(pdb_path: str | Path) -> bool:
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                el = line[76:78].strip().upper()
                if el in {"H", "D", "T"}:
                    return True
                # also detect from atom name if element field is blank/misplaced
                name = line[12:16]
                if (name[0] == " " and name[1].upper() in {"H","D","T"}) or name.lstrip()[:1].upper() in {"H","D","T"}:
                    return True
    except Exception as e:
        logging.warning("Could not read %s to check for H atoms: %s", pdb_path, e)
    return False

def file_was_reduced(pdb_path: str | Path) -> bool:
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if "Reduce" in line and "protonation" in line.lower():
                    return True
    except Exception as e:
        logging.warning("Could not read %s to check for Reduce header: %s", pdb_path, e)
    return False

def fix_element_columns_in_file(src: str | Path, dst: str | Path) -> None:
    METALS = {"ZN","FE","MG","MN","CU","NI","CO","NA","K","CA","CL","BR","SR","BA","CD","HG"}
    PEPTIDEY = {"N","CA","C","O","OXT","CB","CG","CD","CE","CZ","SG",
                "ND","NE","OD","OE","SD","NZ","OH","CH","CZ1","CZ2","CD1","CD2","CE1","CE2","CE3"}

    def derive_from_name(aname: str) -> str:
        an = aname.strip().upper()
        if an[:2] in {"OE","NE","OD","ND","SD"}:
            return an[0]
        if an.startswith("OXT"):
            return "O"
        if an and an[0].isalpha():
            return an[0]
        for ch in an:
            if ch.isalpha():
                return ch
        return "C"

    residues = defaultdict(list)  # key -> list of (idx,line)
    lines = open(src, "r", encoding="utf-8", errors="ignore").read().splitlines(True)

    for i, line in enumerate(lines):
        if line.startswith(("ATOM  ","HETATM")):
            chain = line[21]
            resseq = line[22:26]
            icode = line[26]
            resname = line[17:20].strip().upper()
            key = (chain, resseq, icode, resname)
            residues[key].append((i, line))

    def is_peptidic_like(atom_lines):
        names = [(ln[12:16].strip().upper()) for _, ln in atom_lines]
        if not names:
            return False
        hits = sum((n in PEPTIDEY) or (n[:2] in {"OE","NE","OD","ND","SD"}) for n in names)
        return hits >= max(4, 0.6*len(names))

    def is_ion_like(resname, atom_lines):
        return (resname in METALS) and (len(atom_lines) <= 2)

    for key, atom_lines in residues.items():
        chain, resseq, icode, resname = key
        pep_like = is_peptidic_like(atom_lines)
        ion_like = is_ion_like(resname, atom_lines)
        for idx, old in atom_lines:
            if not old.startswith(("ATOM  ","HETATM")):
                continue
            aname = old[12:16]
            if pep_like:
                el = derive_from_name(aname)
            elif ion_like:
                el = resname[:2] if len(resname) >= 2 else derive_from_name(aname)
            else:
                el = derive_from_name(aname)
            lines[idx] = old[:76] + f"{el:>2}" + old[78:]

    with open(dst, "w", encoding="utf-8") as out:
        out.writelines(lines)

# =============================
# Structural cleaning utilities
# =============================

def filter_altlocs(pdb_input_path: str | Path, pdb_output_path: str | Path) -> None:
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

    logging.info("Filtered altLocs in %s → %s (removed %d alternates)", pdb_input_path, pdb_output_path, removed_count)

def build_missing_loops(input_pdb: str | Path, output_dir: str | Path) -> str:
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
        mdl = complete_pdb(env, str(input_pdb))
        mdl.write(file=output_pdb)
        if os.path.exists(output_pdb):
            logging.info("MODELLER filled PDB saved to %s", output_pdb)
            return output_pdb
        logging.warning("MODELLER did not produce expected output: %s", output_pdb)
        return str(input_pdb)
    except Exception as e:
        logging.warning("MODELLER failed on %s: %s", input_pdb, e)
        return str(input_pdb)

def find_invalid_atoms(pdb_path: str | Path) -> List[Tuple[str, str, int]]:
    """Flag residues that contain CA but have ≤2 atoms (likely malformed)."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", str(pdb_path))

    flagged: List[Tuple[str, str, int]] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                atoms = list(residue.get_atoms())
                atom_names = [a.get_name() for a in atoms]
                if 'CA' in atom_names and len(atoms) <= 2:
                    flagged.append((chain.id, residue.get_resname(), residue.id[1]))
    return flagged

def filter_invalid_chains(pdb_path: str | Path, output_path: str | Path) -> None:
    """Remove entire chains lacking backbone atoms (CA/N/C/O)."""
    chains: dict[str, List[str]] = defaultdict(list)
    valid_chains: set[str] = set()

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

def remove_implausible_hydrogens_by_distance(pdb_path: str | Path) -> None:
    """Remove H atoms >1.35Å from any heavy atom (conservative geometry filter)."""
    atoms: List[Tuple[str, str, Tuple[float | None, float | None, float | None]]] = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    atoms.append((line, "UNK", (None, None, None)))
                    continue
                el = line[76:78].strip()
                if not el:
                    name = line[12:16]
                    el = (name[1] if name[0] == " " else name[0])
                atoms.append((line, el.upper(), (x, y, z)))

    kept: List[str] = []
    heavy_coords = [a[2] for a in atoms if a[1] != "H" and a[2][0] is not None]

    def near_heavy(coord: Tuple[float | None, float | None, float | None]) -> bool:
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

def conect_coverage(pdb_path: str | Path) -> float:
    atom_ids, conect_ids = set(), set()
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_ids.add(line[6:11].strip())
            elif line.startswith("CONECT"):
                parts = line.split()
                conect_ids.update(parts[1:])
    return len(conect_ids & atom_ids) / max(len(atom_ids), 1)

def remove_unbonded_atoms(pdb_path: str | Path) -> None:
    bonded_atoms: set[str] = set()
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

def clean_hydrogens(pdb_path: str | Path, use_conect_if_reliable: bool = True, conect_min_cov: float = 0.6) -> None:
    cov = conect_coverage(pdb_path) if use_conect_if_reliable else 0.0
    if cov >= conect_min_cov:
        remove_unbonded_atoms(pdb_path)
    remove_implausible_hydrogens_by_distance(pdb_path)

def has_valid_chain(pdb_path: str | Path) -> bool:
    with open(pdb_path, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                return True
    return False

def _classify_and_rename_histidines(pdb_in: str | Path, pdb_out: str | Path) -> None:
    """
    Inspect each HIS residue's side-chain hydrogens and rename:
      - HD1 only  -> HID
      - HE2 only  -> HIE
      - both      -> HIP
      - neither   -> (leave as HIS; caller may choose to rewrite to a default)
    Works on plain PDB text; no dependencies.
    """
    pdb_in, pdb_out = str(pdb_in), str(pdb_out)

    # Collect atom lines by residue
    from collections import defaultdict
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
            # No clue from hydrogens (likely pre-Reduce). Keep HIS; caller may override.
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
        # Keep any header/TER/etc lines in relative order (pessimistic but safe)
        for ln in other_lines:
            w.write(ln)
        for ln in out_lines:
            w.write(ln)

def _rewrite_his_default(pdb_in: str | Path, pdb_out: str | Path, default: str = "HIE") -> None:
    """Rewrite any residual HIS → <default> (HIE/HID/HIP) without inspecting hydrogens."""
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
def run_prepare_receptor(input_pdb: str | Path, output_pdbqt: str | Path, cfg: dict) -> bool:
    """
    Robust Meeko/ADT wrapper with:
      • Early check for HIS tautomer ambiguity on the modern (--read_pdb) call
      • Auto -n mapping for all truly ambiguous HIS (no HD1/HE2)
      • Optional -a (allow_bad_res) retry on template-mismatch
      • Heavy-handed HIS→<default> fallback
      • ADT fallback
    """
    import re
    from tempfile import NamedTemporaryFile

    input_pdb = str(input_pdb)
    output_pdbqt = str(output_pdbqt)

    his_default = str(cfg.get("HIS_DEFAULT", "HIE")).upper()
    if his_default not in {"HIE", "HID", "HIP"}:
        his_default = "HIE"

    # improved truthiness parsing
    allow_bad_res = str(cfg.get("MEEKO_ALLOW_BAD_RES", "true")).lower() in ("1", "true", "yes")
    default_altloc = (cfg.get("MEEKO_DEFAULT_ALTLOC") or "").strip()  # e.g. "A" or ""

    def _meeko_cmd() -> list[str]:
        # Prefer a “bare” executable on PATH
        exe = shutil.which("mk_prepare_receptor") or shutil.which("mk_prepare_receptor.py")
        if exe:
            return [exe]
        # Fallback to calling the script with this Python
        abs_py = "/home/michael/miniconda3/envs/docking-env/bin/mk_prepare_receptor.py"
        if os.path.exists(abs_py):
            return [sys.executable, abs_py]
        raise FileNotFoundError("Meeko not found: mk_prepare_receptor(.py) not on PATH and no known absolute path.")
    
    def _run(cmd: list[str]) -> subprocess.CompletedProcess:
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

    def _modern_meeko(pdb_path: str, extra_flags: list[str] | None = None) -> subprocess.CompletedProcess:
        try:
            cmd = _meeko_cmd() + ["--read_pdb", pdb_path, "-p", output_pdbqt]
        except FileNotFoundError as e:
            class _R:  # little stub so the calling code keeps working
                returncode = 127; stdout = ""; stderr = str(e)
            return _R()
        flags = extra_flags or []
        return _run(cmd + flags)

    def _legacy_meeko(pdb_path: str, extra_flags: list[str] | None = None) -> bool:
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

def run_phenix_pdbtools(input_pdb: str | Path, output_pdb: str | Path, remove_waters: bool = True) -> bool:
    """Run Phenix pdbtools if available (Linux, or Windows via PowerShell under WSL). Return True if successful."""
    inp = os.path.abspath(str(input_pdb)).replace("\\", "/")
    outp = os.path.abspath(str(output_pdb)).replace("\\", "/")
    remove_arg = ' remove="resname HOH"' if remove_waters else ""

    # 1) Native Linux phenix.pdbtools
    if _bin_on_path("phenix.pdbtools"):
        cmd = ["phenix.pdbtools", inp, f"output.file_name={outp}"]
        if remove_waters:
            cmd.append('remove="resname HOH"')
        logging.info("Running (Linux) phenix.pdbtools: %s", " ".join(cmd))
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.returncode == 0:
            logging.info("phenix.pdbtools completed successfully")
            return True
        logging.warning("phenix.pdbtools failed (rc=%s):\n%s", r.returncode, r.stderr)
        return False

    # 2) Windows Phenix .bat (usable from WSL via PowerShell)
    if _is_wsl() and PHENIX_DIR and Path(PDBTOOLS_BAT).exists():
        win_in = _win_path(inp)
        win_out = _win_path(outp)
        out_dir = _win_path(Path(outp).parent)
        ps = (
            "$ErrorActionPreference='Stop';"
            f"if (!(Test-Path '{out_dir}')) {{ New-Item -ItemType Directory -Force -Path '{out_dir}' | Out-Null }};"
            f"Push-Location '{out_dir}';"
            f"& '{_win_path(PDBTOOLS_BAT)}' '{win_in}' output.file_name='{win_out}'{remove_arg};"
            "Pop-Location;"
        )
        logging.info("Running (WSL→Windows) phenix.pdbtools.bat via PowerShell")
        r = _powershell(ps)
        if r.returncode == 0:
            logging.info("phenix.pdbtools (Windows) completed successfully")
            return True
        logging.warning("phenix.pdbtools (Windows) failed (rc=%s):\n%s", r.returncode, r.stderr)
        return False

    logging.warning("Phenix not available; skipping pdbtools polish")
    return False

def run_windows_phenix_clean_script(loop_fixed_pdb: str | Path, nolig_dir: str | Path) -> int:
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

def run_openbabel_add_h(input_pdb: str | Path, output_pdb: str | Path) -> None:
    obabel = OPENBABEL_PATH or shutil.which("obabel")
    if not obabel or not shutil.which(Path(obabel).name):
        raise RuntimeError("Open Babel not found. Install or set OPENBABEL_PATH.")
    cmd = [obabel, str(input_pdb), "-O", str(output_pdb), "-h"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Open Babel failed:\n{r.stderr}")
    logging.info("Open Babel hydrogenation succeeded")

def assign_protonation_states(input_pdb: str | Path, output_pdb: str | Path, reduce_exe: str | None = None) -> str:
    import uuid
    input_pdb = os.path.abspath(str(input_pdb)).replace("\\", "/")
    output_pdb = os.path.abspath(str(output_pdb)).replace("\\", "/")

    has_h = hydrogenation_status(input_pdb)[0] != "NO_H"

    reduce_flags = ["-quiet"] if has_h else ["-BUILD", "-quiet"]
    exe = reduce_exe or REDUCE_EXE
    exe_dir = os.path.dirname(exe) or None

    def run_reduce(in_pdb: str, stage_name: str) -> str:
        with open(output_pdb, "w", encoding="utf-8") as out:
            cp = subprocess.run([exe] + reduce_flags + [in_pdb], stdout=out, stderr=subprocess.PIPE, text=True, cwd=exe_dir)
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
            fix_pdb_elements(input_pdb)  # repair element columns; helps Reduce
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

    # ⬇⬇⬇ everything below stays INSIDE the function ⬇⬇⬇
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
            # Optional: allow pipeline to continue if you want (unset to keep strict)
            if str(config.get("ALLOW_NO_HYDROGENS", "")).lower() in ("1","true","yes"):
                logging.warning("Continuing with NO_H due to ALLOW_NO_HYDROGENS config.")
                return output_pdb
            raise RuntimeError("Protonation produced no hydrogens.")

    # make sure we return the path on success
    return output_pdb


# =============================
# Residue & Policy utilities
# =============================

def _is_metal(resname: str) -> bool:
    return (resname or "").upper() in _METALS

def _cofactor_policy_keep(resname: str) -> bool:
    rn = (resname or "").upper()
    pol = (config.get("COFACTOR_POLICY", "auto") or "auto").lower()
    if rn in _COFACTOR_KEEP:
        return True
    if rn in _COFACTOR_DROP:
        return False
    if _is_metal(rn):
        return bool(config.get("KEEP_METALS", True))
    if pol == "keep":
        return rn not in _COFACTOR_DROP
    if pol == "remove":
        return rn in _COFACTOR_KEEP or (_is_metal(rn) and config.get("KEEP_METALS", True))
    # auto: keep likely cofactors unless blacklisted
    return rn not in _COFACTOR_DROP

def _parse_xyz(line: str):
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except Exception:
        return None

def _should_keep_water(line: str, pocket_center: Tuple[float, float, float] | None) -> bool:
    pol = (config.get("WATER_POLICY", "site_only") or "site_only").lower()
    if pol == "remove_all":
        return False
    if pol == "keep_all":
        return True
    if pocket_center is None:
        return pol == "auto"  # without a center, 'auto' defaults to dropping waters
    xyz = _parse_xyz(line)
    if not xyz:
        return False
    dx, dy, dz = xyz[0]-pocket_center[0], xyz[1]-pocket_center[1], xyz[2]-pocket_center[2]
    if sqrt(dx*dx + dy*dy + dz*dz) > float(config.get("WATER_SITE_RADIUS_ANG", 6.0)):
        return False
    try:
        b = float(line[60:66])
        if b > float(config.get("WATER_MAX_BFACTOR", 60.0)):
            return False
    except Exception:
        pass
    return True

def strip_nonstandard_residues(input_pdb: str | Path, output_pdb: str | Path) -> Tuple[int, str]:
    """Remove nonstandard residues while keeping cofactors/metals/waters per policy.

    Returns (#distinct residue names removed, output path).
    """
    standard_residues = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
        "HID","HIE","HIP","SEC","PYL","MSE",
    }

    removed: set[str] = set()
    kept_lines: List[str] = []

    # Estimate pocket center from kept cofactors/metals
    cofm_xyz: List[Tuple[float, float, float]] = []
    with open(input_pdb, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip().upper()
            if line.startswith("HETATM") and (_cofactor_policy_keep(resname) or _is_metal(resname)):
                xyz = _parse_xyz(line)
                if xyz:
                    cofm_xyz.append(xyz)

    pocket_center: Tuple[float, float, float] | None = None
    if cofm_xyz:
        xs, ys, zs = zip(*cofm_xyz)
        from statistics import fmean
        pocket_center = (fmean(xs), fmean(ys), fmean(zs))

    with open(input_pdb, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM  "):
                resname = line[17:20].strip().upper()
                if resname in standard_residues:
                    kept_lines.append(line)
                else:
                    removed.add(resname)
            elif line.startswith("HETATM"):
                resname = line[17:20].strip().upper()
                if resname == "HOH":
                    if _should_keep_water(line, pocket_center):
                        kept_lines.append(line)
                    else:
                        removed.add(resname)
                elif _is_metal(resname):
                    if bool(config.get("KEEP_METALS", True)):
                        kept_lines.append(line)
                    else:
                        removed.add(resname)
                elif _cofactor_policy_keep(resname):
                    kept_lines.append(line)
                else:
                    removed.add(resname)
            else:
                kept_lines.append(line)

    with open(output_pdb, 'w', encoding="utf-8") as f:
        f.writelines(kept_lines)

    logging.info("Removed nonstandard residues: %s", sorted(removed))
    if pocket_center:
        logging.info("Estimated pocket center: (%.2f, %.2f, %.2f)", *pocket_center)
    return len(removed), str(output_pdb)

def log_metal_mislabels(pdb_path: str | Path):
    metal_names = {"NA","K","CA","MG","MN","FE","CO","NI","CU","ZN","CL","BR"}
    mis = {}
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith(("ATOM","HETATM")):
                continue
            resname = ln[17:20].strip().upper()
            elem = ln[76:78].strip().upper()
            if resname in metal_names and elem not in {resname[:2], "H"}:
                mis.setdefault(resname, 0)
                mis[resname] += 1
    if mis:
        logging.warning("Possible mislabels in metal residues: %s", mis)

def quick_element_histogram(pdb_path: str | Path) -> None:
    cnt = Counter()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM","HETATM")):
                el = ln[76:78].strip().upper()
                cnt[el or ""] += 1
    logging.info("[Elem histogram %s] %s", os.path.basename(str(pdb_path)), dict(sorted(cnt.items())))

def assert_no_metal_in_peptidic(pdb_path: str | Path) -> None:
    METALS = {"NA","K","CA","MG","MN","FE","CO","NI","CU","ZN","CL","BR"}
    PEPTIDEY = {"N","CA","C","O","OXT","CB","CG","CD","CE","CZ","SG",
                "ND","NE","OD","OE","SD","NZ","OH","CH","CZ1","CZ2","CD1","CD2","CE1","CE2","CE3"}

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
        hits = sum((n in PEPTIDEY) or (n[:2] in {"OE","NE","OD","ND","SD"}) for n in names)
        pep_like = hits >= max(4, 0.6*len(names))
        if not pep_like:
            continue
        if any(ln[76:78].strip().upper() in METALS for ln in lines):
            offenders.append(key)

    if offenders:
        logging.warning("Peptide-like residues contain metal elements (check upstream labeling): %s", offenders)

# =============================
# End-to-end Cleaning Pipeline
# =============================

def clean_pdb(pdb_file: str | Path, output_root: str | Path) -> str | None:
    """Run the full cleaning pipeline and return path to final cleaned PDB (receptor)."""
    output_root = str(output_root)
    # Ensure base directory exists instead of pre-checking writability
    Path(output_root).mkdir(parents=True, exist_ok=True)

    pdb_id = os.path.splitext(os.path.basename(str(pdb_file)))[0].upper()
    fold_legacy_layout(pdb_id, output_root)
    paths = canon_paths(pdb_id, output_root)
    for d in ["protein_root", "raw", "work", "ligands_raw", "nolig", "receptor"]:
        paths[d].mkdir(parents=True, exist_ok=True)

    # (1) Working copy → raw/
    working_pdb = paths["raw"] / f"{pdb_id}_working.pdb"
    shutil.copyfile(str(pdb_file), working_pdb)

    # (2) AltLoc filtering → raw/filtered.pdb
    filtered_pdb = paths["raw"] / f"{pdb_id}_filtered.pdb"
    filter_altlocs(working_pdb, filtered_pdb)

    # (3) Extract ligands now (controls live here)
    _ = extract_ligands_from_filtered(filtered_pdb, paths["ligands_raw"])

    # (4) Strip nonstandard from protein (policy aware) → work/stripped.pdb
    stripped_pdb = paths["work"] / f"{pdb_id}_stripped.pdb"
    removed_count, _out = strip_nonstandard_residues(filtered_pdb, stripped_pdb)
    logging.info("Removed %d nonstandard residue lines.", removed_count)

    # (5) Element fix → MODELLER → element fix again
    elemfix_pdb = paths["work"] / f"{pdb_id}_elemfix.pdb"
    fix_pdb_elements(stripped_pdb, elemfix_pdb)
    log_metal_mislabels(elemfix_pdb)
    loop_fixed_pdb = build_missing_loops(elemfix_pdb, paths["work"])
    fix_pdb_elements(loop_fixed_pdb, loop_fixed_pdb)
    log_metal_mislabels(loop_fixed_pdb)

    # (6) Phenix clean (kept inside nolig/) — robust across Linux/Windows/WSL; non-fatal on failure
    phenix_out_pdb = paths["nolig"] / f"{pdb_id}_nolig_phenix_clean.pdb"

    ran_clean = False
    # Prefer system python3 for your phenix_clean.py (it’s Python 3 code)
    if PHENIX_CLEAN_SCRIPT:
        if shutil.which("python3"):
            r = subprocess.run(
                ["python3", PHENIX_CLEAN_SCRIPT, str(loop_fixed_pdb), str(paths["nolig"])],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            if r.returncode != 0:
                logging.warning("Phenix clean script failed (python3) rc=%s:\n%s", r.returncode, r.stderr)
            else:
                ran_clean = True
        elif _bin_on_path("phenix.python"):
            r = subprocess.run(
                ["phenix.python", PHENIX_CLEAN_SCRIPT, str(loop_fixed_pdb), str(paths["nolig"])],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            if r.returncode != 0:
                logging.warning("Phenix clean script failed (phenix.python) rc=%s:\n%s", r.returncode, r.stderr)
            else:
                ran_clean = True
        elif _is_wsl():
            rc = run_windows_phenix_clean_script(loop_fixed_pdb, paths["nolig"])
            if rc != 0:
                logging.warning("Phenix clean script failed (Windows via PowerShell) rc=%s", rc)
            else:
                ran_clean = True
        else:
            logging.info("No Phenix interpreter available for PHENIX_CLEAN_SCRIPT; skipping script stage")
    else:
        logging.info("No PHENIX_CLEAN_SCRIPT configured; skipping script stage")

    src_for_step = phenix_out_pdb if ran_clean and phenix_out_pdb.is_file() else loop_fixed_pdb

    # (7) Hydrogen cleanup & chain validation
    debulked_pdb = paths["work"] / f"{pdb_id}_debulked.pdb"
    shutil.copyfile(src_for_step, debulked_pdb)
    clean_hydrogens(debulked_pdb, use_conect_if_reliable=True, conect_min_cov=0.6)

    chain_validated_pdb = paths["work"] / f"{pdb_id}_validated.pdb"
    filter_invalid_chains(debulked_pdb, chain_validated_pdb)

    # (8) Reduce (with fallbacks)
    reduced_pdb = paths["work"] / f"{pdb_id}_reduced.pdb"
    assign_protonation_states(chain_validated_pdb, reduced_pdb)
    if not file_contains_hydrogens(reduced_pdb):
        logging.error("[FATAL] Reduced file missing hydrogens: %s", reduced_pdb)
        return None

    # (9) Final element fix & optional Phenix polish → receptor/
    fix_pdb_elements(reduced_pdb)
    quick_element_histogram(reduced_pdb)

    receptor_pdb = paths["receptor"] / f"{pdb_id}_cleaned.pdb"
    remove_waters_flag = (str(config.get("WATER_POLICY", "site_only")).lower() == "remove_all")

    polished = run_phenix_pdbtools(reduced_pdb, receptor_pdb, remove_waters=remove_waters_flag)
    if not polished:
        # Degrade gracefully: copy reduced file to receptor and continue
        shutil.copyfile(reduced_pdb, receptor_pdb)
        logging.warning("phenix.pdbtools unavailable/failed — using Reduce output as final receptor")

    fix_pdb_elements(receptor_pdb)
    quick_element_histogram(receptor_pdb)
    assert_no_metal_in_peptidic(receptor_pdb)
    assert file_contains_hydrogens(receptor_pdb), f"[FATAL] Cleaned file lost hydrogens: {receptor_pdb}"

    # (10) MolProbity (non-blocking)
    try:
        if _bin_on_path("phenix.molprobity"):
            subprocess.run(["phenix.molprobity", str(receptor_pdb)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        elif _is_wsl() and PHENIX_DIR and Path(MOLPROBITY_BAT).exists():
            _powershell(f"& '{_win_path(MOLPROBITY_BAT)}' '{_win_path(receptor_pdb)}'")
    except Exception as e:
        logging.warning("MolProbity stage failed (non-fatal): %s", e)

    logging.info("Cleaned receptor: %s", receptor_pdb)
    return str(receptor_pdb)

# =============================
# Module Entrypoint
# =============================

def main(pdb_filename: str, output_dir: str | Path = r"./processed_pdbs") -> tuple[str, str] | None:
    """High-level wrapper: clean a PDB and prepare the receptor PDBQT."""
    try:
        # Resolve input path (abs path takes precedence)
        if os.path.isabs(pdb_filename) and os.path.isfile(pdb_filename):
            pdb_path = pdb_filename
        else:
            pdb_path = os.path.join(_cfg("INPUT_DIR", "."), pdb_filename)
        pdb_id = os.path.splitext(os.path.basename(pdb_filename))[0].upper()

        if not os.path.isfile(pdb_path):
            logging.error("ERROR: File does not exist: %s", pdb_path)
            return None

        os.makedirs(output_dir, exist_ok=True)
        cleaned_pdb = clean_pdb(pdb_path, output_dir)
        if not cleaned_pdb:
            logging.error("ERROR: Cleaning failed for %s", pdb_filename)
            return None

        # Prepare receptor PDBQT
        paths = canon_paths(pdb_id, output_dir)
        output_pdbqt = str((paths["receptor"] / f"{pdb_id}.pdbqt").resolve())
        if not run_prepare_receptor(cleaned_pdb, output_pdbqt, config):
            logging.error("ERROR: Failed to prepare receptor PDBQT for %s", pdb_id)
            return None

        # Warn if metals expected but missing in the PDBQT
        if config.get("KEEP_METALS", True):
            try:
                txt = Path(output_pdbqt).read_text(encoding="utf-8", errors="ignore").upper()
                missing = [m for m in _METALS if (" " + m + " ") not in txt and (" " + m + "\n") not in txt]
                if missing:
                    logging.warning("Metals expected but not present in receptor PDBQT: %s", missing)
            except Exception:
                pass

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
