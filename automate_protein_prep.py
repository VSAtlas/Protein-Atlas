"""
Automated protein preparation pipeline (Windows-friendly)

This module prepares a receptor structure from a raw PDB by performing:
  1) Alternate-conformation filtering (altLoc) with deterministic selection
  2) Removal of nonstandard residues (while retaining common cofactors/metals)
  3) Element-column fixes to ensure PDB compliance
  4) Optional loop/residue completion via MODELLER
  5) Sanity checks for incomplete residues
  6) Phenix cleaning pass (robust on Windows)
  7) Hydrogen sanity cleanup (CONECT- and geometry-based)
  8) Chain validation (ensure at least one CA-containing chain)
  9) Protonation via Reduce with automatic fallbacks (temp retry/OpenBabel)
 10) Final polishing via phenix.pdbtools and MolProbity report
 11) AutoDockTools receptor preparation (PDBQT)

Design notes
------------
• Function-level docstrings explain the *why* behind each step.
• All path handling is normalized for Windows subprocess calls.
• Logging is used instead of prints (except where a tool’s console output is useful).
• No functionality is intentionally removed; behavior is clarified and guarded.

Prerequisites
-------------
• Phenix (phenix.python.bat, phenix.pdbtools.bat, phenix.molprobity.bat)
• Reduce (reduce.exe) — usually inside the Phenix bin directory
• MODELLER (licensed) — optional but recommended for loop filling
• Open Babel (obabel) — fallback for hydrogen addition
• Biopython — for PDB parsing
• MGLTools — for prepare_receptor4.py
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

# Project-local imports
from installation import load_config
from activesite import fix_pdb_elements  # element column rewriter (trusted external impl)
from activesite import ensure_model_records  # imported for parity with original file (unused here)
from logger_setup import setup_logger  # assumed to be called by the entrypoint

# =============================
# Configuration & Tool Paths
# =============================

# Load config once at import time (original behavior)
config = load_config()

# --- Cofactor / metal / water policy defaults (safe; change in your config to enable) ---  # << NEW
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

# Build fast-lookup sets once  # << NEW
_COFACTOR_KEEP = {s.strip().upper() for s in config.get("COFACTOR_KEEP_LIST", "").split(",") if s.strip()}
_COFACTOR_DROP = {s.strip().upper() for s in config.get("COFACTOR_REMOVE_LIST", "").split(",") if s.strip()}
_METALS        = {s.strip().upper() for s in config.get("METAL_LIST", "").split(",") if s.strip()}

# Required paths from config
PHENIX_DIR: str = config["PHENIX_DIR"]
PHENIX_LIB_PATH: str = config["phenix_lib_path"]
PHENIX_CLEAN_SCRIPT: str = config["PHENIX_CLEAN_SCRIPT"]
INPUT_DIR: str = config["INPUT_DIR"]

# Make Phenix' Python available for imports (mirrors original behavior)
sys.path.insert(0, PHENIX_LIB_PATH)
sys.path.insert(0, PHENIX_DIR)

PHENIX_DIR_PATH = Path(PHENIX_DIR)

# Prefer explicit REDUCE_EXE from config; otherwise use default inside PHENIX_DIR
REDUCE_EXE_PATH = Path(config.get("REDUCE_EXE", PHENIX_DIR_PATH / "reduce.exe"))

# Batch wrappers (Windows)
PDBTOOLS_BAT = PHENIX_DIR_PATH / "phenix.pdbtools.bat"
MOLPROBITY_BAT = PHENIX_DIR_PATH / "phenix.molprobity.bat"
PHENIX_PYTHON_BAT = PHENIX_DIR_PATH / "phenix.python.bat"

# Normalize for subprocess on Windows (prefer forward slashes)
REDUCE_EXE = str(REDUCE_EXE_PATH).replace("\\", "/")
PDBTOOLS_BAT = str(PDBTOOLS_BAT).replace("\\", "/")
MOLPROBITY_BAT = str(MOLPROBITY_BAT).replace("\\", "/")
PHENIX_PYTHON_BAT = str(PHENIX_PYTHON_BAT).replace("\\", "/")

logging.debug("[DEBUG] PHENIX_DIR = %s", PHENIX_DIR)
logging.debug("[DEBUG] REDUCE_EXE = %s", REDUCE_EXE)
logging.debug("[DEBUG] PDBTOOLS_BAT = %s (exists=%s)", PDBTOOLS_BAT, os.path.exists(PDBTOOLS_BAT))

# Two-letter elements commonly encountered in PDBs (for element inference)
_TWO_LETTER = {
    "ZN","FE","MG","MN","CL","NA","CA","CU","CO","BR","SI","PT","PD","NI","AL",
    "AG","AU","IR","SR","BA","BE","LI","RB","CS","MO","SE","TE","TI","CR","VD",
    "HG","PB","SN","SB","CD","GA","GE","ZR","Y","NB","W","RE","OS","RH","RU","I"
}

# =============================
# Element & Hydrogen Utilities
# =============================

def looks_like_hydrogen_name(line: str) -> bool:
    """Return True if the *atom name field* looks like a hydrogen."""
    name = line[12:16]
    return name[0] == " " and name[1].upper() == "H"

def infer_element(line: str) -> str:
    """Infer the element symbol for an ATOM/HETATM line."""
    el = line[76:78].strip().upper()
    if el:
        return "H" if el in {"D", "T"} else el
    name = line[12:16]
    c1 = name[1] if name[0] == " " else name[0]
    c2 = (name[2] if name[0] == " " else name[1]).strip()
    cand2 = (c1 + c2).upper()
    return cand2 if cand2 in _TWO_LETTER else c1.upper()

def count_atoms_by_element(pdb_path: str | Path) -> Tuple[int, int]:
    """Count hydrogens vs. heavy atoms in a PDB file."""
    h = heavy = 0
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            el = infer_element(line)
            if el == "H":
                if line[76:78].strip() == "" and not looks_like_hydrogen_name(line):
                    heavy += 1
                else:
                    h += 1
            else:
                heavy += 1
    return h, heavy

def hydrogenation_status(pdb_path: str | Path) -> Tuple[str, int, int, float]:
    """Categorize hydrogenation level: NO_H, SUSPECT_LOW_H (<20%), or HAS_H."""
    h, heavy = count_atoms_by_element(pdb_path)
    ratio = h / max(heavy, 1)
    if h == 0:
        return "NO_H", h, heavy, ratio
    if ratio < 0.20:
        return "SUSPECT_LOW_H", h, heavy, ratio
    return "HAS_H", h, heavy, ratio

def file_contains_hydrogens(pdb_path: str | Path) -> bool:
    """Fast check for presence of H/D/T in the element column."""
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    el = line[76:78].strip().upper()
                    if el in {"H", "D", "T"}:
                        return True
    except Exception as e:
        logging.warning("Could not read %s to check for H atoms: %s", pdb_path, e)
    return False

def file_was_reduced(pdb_path: str | Path) -> bool:
    """Heuristic: detect a Reduce header indicating prior protonation."""
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if "Reduce" in line and "protonation" in line.lower():
                    return True
    except Exception as e:
        logging.warning("Could not read %s to check for Reduce header: %s", pdb_path, e)
    return False

def fix_element_columns_in_file(src: str | Path, dst: str | Path) -> None:
    """Fill empty element columns using atom-name heuristics."""
    out_lines: List[str] = []
    with open(src) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                el = line[76:78].strip()
                if not el:
                    name = line[12:16]
                    c = name[1] if name[0] == " " else name[0]
                    c2 = (name[2] if name[0] == " " else name[1]).strip()
                    cand = (c + c2).upper()
                    el = cand if cand in _TWO_LETTER else c.upper()
                line = line[:76] + f"{el:>2}" + line[78:]
            out_lines.append(line)
    with open(dst, "w") as out:
        out.writelines(out_lines)

# =============================
# Structural Cleaning Utilities
# =============================

def filter_altlocs(pdb_input_path: str | Path, pdb_output_path: str | Path) -> None:
    """Filter alternate locations (altLoc) deterministically and log decisions."""
    lines: List[str] = []
    atoms = {}
    removed_count = 0

    with open(pdb_input_path, 'r') as f:
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

    with open(pdb_output_path, 'w') as f:
        for line in lines:
            f.write(line)
        for atom_line in filtered_atoms:
            f.write(atom_line)

    logging.info("Filtered altLocs in %s → %s (removed %d alternates)", pdb_input_path, pdb_output_path, removed_count)

def build_missing_loops(input_pdb: str | Path, output_dir: str | Path) -> str:
    """Fill missing loops/residues using MODELLER; return output PDB path."""
    from modeller import environ, log
    from modeller.scripts import complete_pdb

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
    from collections import defaultdict

    chains: dict[str, List[str]] = defaultdict(list)
    valid_chains: set[str] = set()

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                chain_id = line[21]
                atom_name = line[12:16].strip()
                chains[chain_id].append(line)
                if atom_name in {"CA", "N", "C", "O"}:
                    valid_chains.add(chain_id)
            else:
                chains["HEADER"].append(line)

    with open(output_path, 'w') as f:
        for chain_id in chains:
            if chain_id == "HEADER" or chain_id in valid_chains:
                f.writelines(chains[chain_id])
            else:
                logging.warning("Skipping invalid chain '%s' (no CA atoms)", chain_id)

def remove_implausible_hydrogens_by_distance(pdb_path: str | Path) -> None:
    """Remove H atoms >1.35Å from any heavy atom (very conservative geometry filter)."""
    atoms: List[Tuple[str, str, Tuple[float | None, float | None, float | None]]] = []
    with open(pdb_path) as f:
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

    with open(pdb_path, "w") as out:
        out.writelines(kept)

def conect_coverage(pdb_path: str | Path) -> float:
    """Return fraction of atoms that appear in any CONECT record."""
    atom_ids, conect_ids = set(), set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_ids.add(line[6:11].strip())
            elif line.startswith("CONECT"):
                parts = line.split()
                conect_ids.update(parts[1:])
    return len(conect_ids & atom_ids) / max(len(atom_ids), 1)

def remove_unbonded_atoms(pdb_path: str | Path) -> None:
    """Drop hydrogens that never appear in any CONECT bond list."""
    bonded_atoms: set[str] = set()
    all_atoms: List[str] = []
    with open(pdb_path, 'r') as f:
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

    with open(pdb_path, 'w') as f:
        f.writelines(filtered)

def clean_hydrogens(pdb_path: str | Path, use_conect_if_reliable: bool = True, conect_min_cov: float = 0.6) -> None:
    """Two-pass H cleanup."""
    cov = conect_coverage(pdb_path) if use_conect_if_reliable else 0.0
    if cov >= conect_min_cov:
        remove_unbonded_atoms(pdb_path)
    remove_implausible_hydrogens_by_distance(pdb_path)

def has_valid_chain(pdb_path: str | Path) -> bool:
    """Return True if at least one chain has a CA atom (very loose validity)."""
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                return True
    return False

# =============================
# External Tools (ADT, Phenix, OpenBabel, Reduce)
# =============================

def run_prepare_receptor(input_pdb: str | Path, output_pdbqt: str | Path, config: dict) -> bool:
    """Run MGLTools' prepare_receptor4.py to generate a receptor PDBQT."""
    mgltools_python = config["MGLTOOLS_PYTHON"]
    prepare_script = config["PREPARE_RECEPTOR_SCRIPT"]

    cmd = [
        mgltools_python,
        prepare_script,
        "-r", str(input_pdb),
        "-o", str(output_pdbqt),
        "-A", "none",
        "-U", "nphs_lps_nonstdres",
    ]
    logging.info("Running prepare_receptor4.py: %s", " ".join(cmd))

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        logging.error("prepare_receptor4 failed:\n%s", result.stderr)
        return False

    logging.info("Receptor prepared: %s", output_pdbqt)
    return True

def run_phenix_pdbtools(input_pdb: str | Path, output_pdb: str | Path, remove_waters: bool = True) -> bool:
    """Run phenix.pdbtools to polish a PDB after Reduce."""
    input_pdb = os.path.abspath(str(input_pdb)).replace("\\", "/")
    output_pdb = os.path.abspath(str(output_pdb)).replace("\\", "/")
    remove_arg = ' remove="resname HOH"' if remove_waters else ""
    cmd = f'"{PDBTOOLS_BAT}" "{input_pdb}" output.file_name="{output_pdb}"{remove_arg}'

    logging.info("Running: %s", cmd)
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        logging.warning("phenix.pdbtools failed (rc=%s):\n%s", result.returncode, result.stderr)
        return False

    logging.info("phenix.pdbtools completed successfully")
    return True

def run_openbabel_add_h(input_pdb: str | Path, output_pdb: str | Path) -> None:
    """Fallback hydrogenation via Open Babel (`obabel -h`)."""
    obabel = shutil.which("obabel") or config.get("OPENBABEL_PATH")
    if not obabel or not os.path.exists(obabel):
        raise RuntimeError("Open Babel not found. Install or set OPENBABEL_PATH.")

    cmd = [obabel, str(input_pdb), "-O", str(output_pdb), "-h"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Open Babel failed:\n{result.stderr}")
    logging.info("Open Babel hydrogenation succeeded")

def assign_protonation_states(input_pdb: str | Path, output_pdb: str | Path, reduce_exe: str | None = None) -> str:
    """Add/optimize hydrogens with Reduce; retry/fallback automatically."""
    import uuid

    input_pdb = os.path.abspath(str(input_pdb)).replace("\\", "/")
    output_pdb = os.path.abspath(str(output_pdb)).replace("\\", "/")

    already_reduced = file_was_reduced(input_pdb)
    has_h = hydrogenation_status(input_pdb)[0] != "NO_H"

    reduce_flags = ["-quiet"] if has_h else ["-BUILD", "-quiet"]
    exe = reduce_exe or REDUCE_EXE
    exe_dir = os.path.dirname(exe) or None

    def run_reduce(in_pdb: str, stage_name: str) -> str:
        with open(output_pdb, "w") as out:
            cp = subprocess.run([exe] + reduce_flags + [in_pdb], stdout=out, stderr=subprocess.PIPE, text=True, cwd=exe_dir)
        if cp.returncode != 0:
            raise RuntimeError(f"{stage_name} reduce failed: {cp.stderr.strip()}")
        return output_pdb

    status_before, h0, hv0, r0 = hydrogenation_status(input_pdb)
    logging.info("[H-Scan before] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_before, h0, hv0, r0)

    # Attempt 1: direct Reduce
    try:
        run_reduce(input_pdb, "Reduce#1")
    except Exception as e:
        logging.warning("Reduce#1 failed: %s", e)
        fix_pdb_elements(input_pdb)  # repair element columns; helps Reduce
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

    status_after, h1, hv1, r1 = hydrogenation_status(output_pdb)
    logging.info("[H-Scan after ] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_after, h1, hv1, r1)
    if status_after == "NO_H":
        raise RuntimeError("Protonation produced no hydrogens.")

    return output_pdb

# =============================
# Residue & Policy Utilities
# =============================

def guess_element_from_atom_name(name: str) -> str:
    """Very small helper for element inference from atom *name*."""
    name = name.strip()
    if len(name) == 4:
        name = name[1:] if name[0].isdigit() or name[0] == ' ' else name
    name = name.upper()
    two_letter = {"ZN", "FE", "MG", "MN", "CL", "NA", "CA", "CU", "CO"}
    if name[:2] in two_letter:
        return name[:2]
    if name[:1] in {"C", "H", "O", "N", "S", "P"}:
        return name[0]
    return "C"

def strip_incomplete_residues(pdb_path: str | Path, output_path: str | Path, min_atoms: int = 3) -> None:
    """Write a copy of the structure lacking residues with < min_atoms atoms."""
    from Bio.PDB import PDBParser, PDBIO, Select

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", str(pdb_path))

    class GoodResidueSelect(Select):
        def accept_residue(self, residue):
            return len(list(residue.get_atoms())) >= min_atoms

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), GoodResidueSelect())

# ---- Cofactor / Metal / Water helpers ----  # << NEW
def _is_metal(resname: str) -> bool:
    return (resname or "").upper() in _METALS

def _cofactor_policy_keep(resname: str) -> bool:
    """Return True if this HETATM residue should be kept in the receptor, under policy."""
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
    # auto: keep likely cofactors unless blacklisted; simple name-based rule
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
    # site_only / auto: keep HOH within radius of a tentative pocket center
    if pocket_center is None:
        return pol == "auto"  # without a center, 'auto' defaults to dropping waters
    xyz = _parse_xyz(line)
    if not xyz:
        return False
    from math import sqrt
    dx, dy, dz = xyz[0]-pocket_center[0], xyz[1]-pocket_center[1], xyz[2]-pocket_center[2]
    if sqrt(dx*dx + dy*dy + dz*dz) > float(config.get("WATER_SITE_RADIUS_ANG", 6.0)):
        return False
    # Optional B-factor screen; ignore if blank
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

    # First pass: gather candidate cofactor/metal coordinates to estimate a pocket center
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
        import statistics
        xs, ys, zs = zip(*cofm_xyz)
        pocket_center = (statistics.fmean(xs), statistics.fmean(ys), statistics.fmean(zs))

    # Second pass: write filtered file according to policy
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
        logging.info("Estimated pocket center from kept cofactors/metals: (%.2f, %.2f, %.2f)",
                     pocket_center[0], pocket_center[1], pocket_center[2])
    return len(removed), str(output_pdb)

# =============================
# End-to-end Cleaning Pipeline
# =============================

def clean_pdb(pdb_file: str | Path, output_root: str | Path) -> str | None:
    """Run the full cleaning pipeline and return path to final cleaned PDB."""
    output_root = str(output_root)
    if not os.access(output_root, os.W_OK):
        raise PermissionError(f"Cannot write to output directory: {output_root}")

    pdb_id = os.path.splitext(os.path.basename(str(pdb_file)))[0]
    if pdb_id.endswith("_cleaned"):
        pdb_id = pdb_id[:-8]

    output_dir = os.path.join(output_root, pdb_id)
    os.makedirs(output_dir, exist_ok=True)

    # 1) Working copy & altLoc filtering
    working_pdb = os.path.join(output_dir, f"{pdb_id}_working.pdb")
    shutil.copyfile(str(pdb_file), working_pdb)

    filtered_pdb = os.path.join(output_dir, f"{pdb_id}_filtered.pdb")
    filter_altlocs(working_pdb, filtered_pdb)

    # 1.5) Remove nonstandard residues (policy-aware keep of cofactors/metals/waters)
    stripped_pdb = os.path.join(output_dir, f"{pdb_id}_stripped.pdb")
    removed_count, _ = strip_nonstandard_residues(filtered_pdb, stripped_pdb)
    logging.info("Removed %d nonstandard residue lines.", removed_count)

    # 2) Element fix → MODELLER loop fill → element fix again (defensive)
    fixed_elements_pdb = os.path.join(output_dir, f"{pdb_id}_elemfix.pdb")
    fix_element_columns_in_file(stripped_pdb, fixed_elements_pdb)

    loop_fixed_pdb = build_missing_loops(fixed_elements_pdb, output_dir)
    fix_element_columns_in_file(loop_fixed_pdb, loop_fixed_pdb)

    # 3) Incomplete residue visibility
    invalid_residues = find_invalid_atoms(loop_fixed_pdb)
    if invalid_residues:
        for chain_id, resname, resnum in invalid_residues:
            logging.warning("Incomplete residue: %s %s%d has ≤2 atoms incl. CA", resname, chain_id, resnum)
    else:
        logging.info("No incomplete CA residues found.")

    # 4) Phenix clean pass (Windows-friendly)
    phenix_input = os.path.abspath(loop_fixed_pdb).replace("\\", "/")
    phenix_out_dir = os.path.abspath(os.path.join(output_root, f"{pdb_id}_nolig")).replace("\\", "/")
    os.makedirs(phenix_out_dir, exist_ok=True)
    phenix_out_pdb = os.path.join(phenix_out_dir, f"{pdb_id}_nolig_phenix_clean.pdb").replace("\\", "/")

    r = subprocess.run(
        ["cmd.exe", "/c", PHENIX_PYTHON_BAT, PHENIX_CLEAN_SCRIPT, phenix_input, phenix_out_dir],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if r.returncode != 0:
        logging.warning("Phenix clean script failed (rc=%s). stderr:\n%s", r.returncode, r.stderr)

    src_for_step4 = phenix_out_pdb if os.path.isfile(phenix_out_pdb) else loop_fixed_pdb

    # Handle transient file locks on Windows
    if src_for_step4 == phenix_out_pdb:
        import time
        for _ in range(20):  # ~4 seconds total
            try:
                with open(src_for_step4, "rb"):
                    break
            except PermissionError:
                time.sleep(0.2)
        else:
            logging.warning("Phenix output stays locked; falling back to loop_fixed_pdb.")
            src_for_step4 = loop_fixed_pdb

    # 5) Hydrogen cleanup & chain validation
    debulked_pdb = os.path.join(output_dir, f"{pdb_id}_debulked.pdb")
    shutil.copyfile(src_for_step4, debulked_pdb)

    clean_hydrogens(debulked_pdb, use_conect_if_reliable=True, conect_min_cov=0.6)

    chain_validated_pdb = os.path.join(output_dir, f"{pdb_id}_validated.pdb")
    filter_invalid_chains(debulked_pdb, chain_validated_pdb)

    # 6) Protonation & verification
    reduced_pdb = os.path.join(output_dir, f"{pdb_id}_reduced.pdb")
    assign_protonation_states(chain_validated_pdb, reduced_pdb)
    if not file_contains_hydrogens(reduced_pdb):
        logging.error("[FATAL] Reduced file missing hydrogens: %s", reduced_pdb)
        return None

    # Final element fix (+ external fix for anything pdbtools changed)
    fix_pdb_elements(reduced_pdb)

    cleaned_pdb = os.path.abspath(os.path.join(output_dir, f"{pdb_id}_cleaned.pdb")).replace("\\", "/")
    # Respect WATER_POLICY here so Phenix doesn't strip waters you decided to keep  # << UPDATED
    remove_waters_flag = (config.get("WATER_POLICY", "site_only").lower() == "remove_all")
    if not run_phenix_pdbtools(reduced_pdb, cleaned_pdb, remove_waters=remove_waters_flag):
        return None

    # One more element sanity pass in case pdbtools touched columns
    fix_pdb_elements(cleaned_pdb)
    assert file_contains_hydrogens(cleaned_pdb), f"[FATAL] Cleaned file lost hydrogens: {cleaned_pdb}"

    # 7) MolProbity report (non-blocking)
    molprobity_log = os.path.join(output_dir, f"{pdb_id}_molprobity.log")
    with open(molprobity_log, "w") as out:
        subprocess.run([MOLPROBITY_BAT, cleaned_pdb], stdout=out)

    logging.info("Cleaned: %s", cleaned_pdb)
    logging.info("Input PDB: %s", pdb_file)
    logging.info("Cleaned PDB will be written to: %s", cleaned_pdb)

    # 8) Final parse sanity check
    from Bio.PDB import PDBParser
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("validate", cleaned_pdb)
        atoms = list(structure.get_atoms())
        assert len(atoms) > 0, f"[FATAL] Final cleaned PDB has no atoms: {cleaned_pdb}"
    except Exception as e:
        logging.error("[FATAL] Could not parse cleaned PDB: %s", e)

    return cleaned_pdb

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
            pdb_path = os.path.join(INPUT_DIR, pdb_filename)
        pdb_id = os.path.splitext(os.path.basename(pdb_filename))[0].upper()

        if not os.path.isfile(pdb_path):
            logging.error("ERROR: File does not exist: %s", pdb_path)
            return None

        os.makedirs(output_dir, exist_ok=True)
        cleaned_pdb = clean_pdb(pdb_path, output_dir)
        if not cleaned_pdb:
            logging.error("ERROR: Cleaning failed for %s", pdb_filename)
            return None

        # Prepare receptor PDBQT with ADT
        output_pdbqt = os.path.join(output_dir, pdb_id, f"{pdb_id}.pdbqt")
        if not run_prepare_receptor(cleaned_pdb, output_pdbqt, config):
            logging.error("ERROR: Failed to prepare receptor PDBQT for %s", pdb_id)
            return None

        # Optional sanity: warn if metals were expected but missing in the PDBQT  # << NEW
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

# Example usage (leave commented for library mode)
# if __name__ == "__main__":
#     setup_logger()
#     test_pdb = "1a3n.pdb"
#     result = main(test_pdb)
#     if result:
#         print("Final cleaned PDB:", result[0])
