#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand preparation utilities (WSL/Windows/Linux friendly)

Highlights
----------
- Meeko vs ADT selection (config: USE_MEEKO, PREPARE_LIGAND_SCRIPT, MGLTOOLS_PYTHON)
- Optional pH pre-protonation via Open Babel BEFORE prep (config: LIGAND_PH), never passed to Meeko
- Meeko charge model control (config: MEEKO_CHARGE_MODEL = gasteiger|espaloma|zero)
- Detailed, reproducible logging (command lines, rc/stdout/stderr)
- Robust fallbacks & validation (quarantine on invalid outputs)
- RDKit ETKDG 3D for SDF batch path (optional)
- Windows/WSL path helpers for Open Babel & MGLTools

Add to config.txt (examples)
----------------------------
USE_MEEKO = true
PREPARE_LIGAND_SCRIPT = /path/to/mk_prepare_ligand.py
MGLTOOLS_PYTHON = /path/to/python
OPENBABEL_PATH = obabel

# optional:
LIGAND_PH = 7.4
MEEKO_CHARGE_MODEL = gasteiger
MEEKO_BAD_CHARGE_OK = false

Paths (optional, batch mode):
LIGAND_EXTRACTED_DIR = ./inputs/sdf
LIGANDS_MOL2_DIR     = ./work/mol2
OUTPUT_LIGANDS_DIR   = ./out/pdbqt
"""

from __future__ import annotations

import os
import sys
import ctypes
import shutil
import math
import logging
import subprocess
from datetime import datetime
from collections import defaultdict
from ctypes import wintypes, create_unicode_buffer
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any

# ---- Project-local imports ----
from activesite import fix_pdb_elements  # keep your existing helper

# ---- RDKit imports ----
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize as _std
    _HAS_STD = True
except Exception:
    _std = None
    _HAS_STD = False

# =========================
# Constants & configuration
# =========================

RUN_TAG = datetime.now().strftime("%Y%m%d_%H%M%S")
MALFORMED_LOG = Path(f"malformed_ligands_{RUN_TAG}.txt")
QUARANTINE_DIRNAME = "quarantine"

# tuning
RESUME_SKIP        = False
USE_RDKIT_FOR_3D   = True
OBABEL_THREADS     = 16
OBABEL_TIMEOUT_S   = 900
CHUNK_SIZE         = 200

ALLOWED_ELEMENTS = {
    "H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I",
    "B", "Si", "Se", "Zn", "Mg", "Ca", "Mn", "Fe", "K", "Na"
}
MAX_HEAVY_ATOMS       = 1200
MIN_ATOMS_FOR_DOCKING = 5
MIN_PARENT_HEAVY      = 8

STANDARD_AMINO_ACIDS = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "HID","HIE","HIP","SEC","PYL","MSE"
}

EXCLUDE_CRYSTAL_ADDITIVES = {
    "HOH","CIT","TAR","SO4","PO4","CA","NA","K","MG","MN","ZN",
    "GOL","EDO","PEG","MPD","TRS","MES","HEPES","ACET","ACT","FMT",
    "MAL","DMS","IPA","CLU","NAG","BOG",
    "TOS","BES","PTS","OTF","TRF","TFA","BF4","PF6","CL","BR","I",
}

MONOATOMIC_IONS = {
    "Cl","Br","I","F","Na","K","Ca","Mg","Zn","Mn","Fe","Cu","Co","Ni",
    "Al","Ag","Au","Pt","Li","Ba","Sr","Cs","Rb"
}

_REM = SaltRemover()

# =========================
# Config & small utilities
# =========================

def read_config(path: str = "config.txt") -> Dict[str, str]:
    cfg: Dict[str, str] = {}
    if not Path(path).exists():
        return cfg
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                cfg[k.strip()] = v.strip()
    return cfg

def _is_true(x: Any) -> bool:
    return str(x).strip().lower() in {"1", "true", "yes", "on"}

def get_short_path_name(long_name: str) -> str:
    """Return Windows short path if on Windows, else unchanged."""
    if sys.platform != "win32":
        return long_name
    _GetShortPathNameW = ctypes.windll.kernel32.GetShortPathNameW
    _GetShortPathNameW.argtypes = [wintypes.LPCWSTR, wintypes.LPWSTR, wintypes.DWORD]
    _GetShortPathNameW.restype = wintypes.DWORD
    buf_size = 260
    buf = create_unicode_buffer(buf_size)
    ret = _GetShortPathNameW(long_name, buf, buf_size)
    return long_name if ret == 0 or ret > buf_size else buf.value

def resolve_obabel_exe(obabel_cfg: str) -> str:
    """
    Accepts either:
      - the literal binary name on PATH (e.g., "obabel"),
      - a full path to the binary (Linux or Windows),
      - or a directory containing the binary.
    Returns a usable string to execute.
    """
    p = Path(obabel_cfg)
    if p.exists() and p.is_file():
        return str(p)
    if p.exists() and p.is_dir():
        for name in ("obabel", "obabel.exe"):
            cand = p / name
            if cand.exists():
                return str(cand)
    return obabel_cfg  # assume on PATH

def _maybe_set_babel_datadir(obabel_exe: str) -> None:
    """If the OpenBabel data/ directory is adjacent to the binary, export BABEL_DATADIR."""
    if os.environ.get("BABEL_DATADIR"):
        return
    try:
        obdir = Path(obabel_exe).resolve().parent
        data = obdir / "data"
        if data.exists():
            os.environ["BABEL_DATADIR"] = str(data)
    except Exception:
        pass

# =========================
# Logging helpers
# =========================

def _log_malformed(path: Path, reason: str) -> None:
    try:
        MALFORMED_LOG.parent.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    with open(MALFORMED_LOG, "a", encoding="utf-8") as fh:
        fh.write(f"{path}\t{reason}\n")

# =========================
# Chemistry helpers
# =========================

def _quick_filters(mol: Chem.Mol) -> Tuple[bool, str]:
    heavy = mol.GetNumHeavyAtoms()
    if heavy == 0:
        return False, "no_heavy_atoms"
    if heavy > MAX_HEAVY_ATOMS:
        return False, f"too_large({heavy})"
    if mol.GetNumAtoms() < MIN_ATOMS_FOR_DOCKING:
        return False, f"too_few_atoms({mol.GetNumAtoms()})"
    for a in mol.GetAtoms():
        if a.GetSymbol() not in ALLOWED_ELEMENTS:
            return False, f"disallowed_element:{a.GetSymbol()}"
    return True, ""

def _standardize_then_sanitize(mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], str]:
    if not _HAS_STD:
        return None, "std_module_missing"
    try:
        md = _std.MetalDisconnector()
        fr = _std.FragmentRemover()
        lf = _std.LargestFragmentChooser(preferOrganic=True)
        uc = _std.Uncharger()
        m  = uc.uncharge(lf.choose(fr.RemoveFragments(md.Disconnect(mol))))
        Chem.SanitizeMol(m)
        ok, why = _quick_filters(m)
        return (m, "") if ok else (None, why)
    except Exception as e:
        return None, f"std_resanitize_fail:{e}"

def _element_from_adt(adt: str) -> str:
    """Map ADT atom type tail token to element."""
    t = (adt or "").strip()
    if not t:
        return "C"
    u = t.upper()
    MAP = {
        "C":"C","A":"C", "N":"N","NA":"N", "O":"O","OA":"O", "S":"S","SA":"S",
        "H":"H","HD":"H", "F":"F","CL":"Cl","BR":"Br","I":"I", "P":"P",
        "B":"B","SI":"Si","SE":"Se", "ZN":"Zn","MG":"Mg","CA":"Ca","MN":"Mn","FE":"Fe",
        "K":"K","NA+":"Na","NA_":"Na","NA ":"Na",
    }
    if u in MAP: return MAP[u]
    if u == "CL": return "Cl"
    if u == "BR": return "Br"
    return u[0]

# ----- counterions / buffers -----

_COUNTERION_SMARTS = {
    "mesylate":    Chem.MolFromSmarts("[CH3]-S(=O)(=O)[O-]"),
    "tosylate":    Chem.MolFromSmarts("c1cccc(c1)S(=O)(=O)[O-]"),
    "triflate":    Chem.MolFromSmarts("C(F)(F)F-S(=O)(=O)[O-]"),
    "sulfonate":   Chem.MolFromSmarts("S(=O)(=O)[O-]"),
    "phosphate":   Chem.MolFromSmarts("P(=O)([O-])([O-])[O-]"),
    "sulfate":     Chem.MolFromSmarts("S(=O)(=O)([O-])[O-]"),
    "formate":     Chem.MolFromSmarts("[#6](=O)[O-]"),
    "acetate":     Chem.MolFromSmarts("CC(=O)[O-]"),
    "lactate":     Chem.MolFromSmarts("CC(O)C(=O)[O-]"),
    "citrate_like": Chem.MolFromSmarts("[CX4](-[CH2]-C(=O)[O-])(-[CH2]-C(=O)[O-])(-C(=O)[O-])O"),
    "tartrate_like": Chem.MolFromSmarts("OC([CH](O)C(=O)[O-])C(=O)[O-]"),
}

def _looks_like_buffer_salt(m: Chem.Mol) -> bool:
    from rdkit.Chem import rdMolDescriptors as rdmd
    hac = m.GetNumHeavyAtoms()
    rings = rdmd.CalcNumRings(m)
    arom  = rdmd.CalcNumAromaticRings(m)
    o     = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'O')
    n     = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'N')
    return hac <= 15 and o >= 6 and n == 0 and rings == 0 and arom == 0

def _matches_counterion(m: Chem.Mol) -> Optional[str]:
    try:
        hac = m.GetNumHeavyAtoms()
        if hac == 0:
            return None
        for name in ("mesylate", "tosylate", "triflate"):
            patt = _COUNTERION_SMARTS[name]
            if patt and m.HasSubstructMatch(patt):
                return name
        for name in ("citrate_like", "tartrate_like"):
            patt = _COUNTERION_SMARTS[name]
            if patt and m.HasSubstructMatch(patt):
                return name
        for name in ("phosphate", "sulfate", "sulfonate"):
            patt = _COUNTERION_SMARTS[name]
            if patt and m.HasSubstructMatch(patt) and hac <= 14:
                return name
        if hac <= 4 and _COUNTERION_SMARTS["formate"] and m.HasSubstructMatch(_COUNTERION_SMARTS["formate"]):
            return "formate"
        if hac <= 5 and _COUNTERION_SMARTS["acetate"] and m.HasSubstructMatch(_COUNTERION_SMARTS["acetate"]):
            return "acetate"
        if hac <= 6 and _COUNTERION_SMARTS["lactate"] and m.HasSubstructMatch(_COUNTERION_SMARTS["lactate"]):
            return "lactate"
        s = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "S")
        o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")
        c = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")
        if hac <= 8 and s == 1 and o >= 3 and c <= 2:
            return "small_sulfonate_like"
        return None
    except Exception:
        return None

_CARBOXYLATE = Chem.MolFromSmarts("[CX3](=O)[O-]")
_CARBOXYLIC  = Chem.MolFromSmarts("[CX3](=O)O")

def _is_polyacidic_buffer_like(m: Chem.Mol) -> bool:
    try:
        from rdkit.Chem import rdMolDescriptors as rdmd
        hac = m.GetNumHeavyAtoms()
        if hac == 0: return False
        rings = rdmd.CalcNumRings(m)
        arom  = rdmd.CalcNumAromaticRings(m)
        o     = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")
        n     = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N")
        o_ratio = (o / float(hac)) if hac else 0.0
        na = 0
        if _CARBOXYLATE: na += len(m.GetSubstructMatches(_CARBOXYLATE))
        if _CARBOXYLIC:  na += len(m.GetSubstructMatches(_CARBOXYLIC))
        has_sulfate   = bool(_COUNTERION_SMARTS["sulfate"]   and m.HasSubstructMatch(_COUNTERION_SMARTS["sulfate"]))
        has_phosphate = bool(_COUNTERION_SMARTS["phosphate"] and m.HasSubstructMatch(_COUNTERION_SMARTS["phosphate"]))
        return (rings == 0 and arom == 0 and o_ratio >= 0.35 and (na >= 3 or has_sulfate or has_phosphate) and n <= 1)
    except Exception:
        return False

def _buffer_like_by_counts_from_mol(m: Chem.Mol) -> bool:
    hac = m.GetNumHeavyAtoms()
    if hac == 0: return False
    o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'O')
    n = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'N')
    return (hac >= 10 and (o / float(hac)) >= 0.40 and n <= 1)

def _element_counts_from_pdb(pdb_path: Path) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fin:
            for ln in fin:
                if not ln.startswith(("ATOM", "HETATM")):
                    continue
                parts = ln.split()
                if not parts: continue
                elem = _element_from_adt(parts[-1])
                counts[elem] += 1
    except Exception:
        pass
    return counts

def _buffer_like_by_counts_from_pdbfile(pdb_path: Path) -> bool:
    c = _element_counts_from_pdb(pdb_path)
    hac = sum(v for k, v in c.items() if k != "H")
    o   = c.get("O", 0)
    n   = c.get("N", 0)
    return hac >= 10 and hac > 0 and (o / float(hac)) >= 0.40 and n <= 1

def _polyacidic_by_counts_from_pdbfile(pdb_path: Path) -> bool:
    c = _element_counts_from_pdb(pdb_path)
    hac = sum(v for k, v in c.items() if k != "H")
    o   = c.get("O", 0)
    n   = c.get("N", 0)
    return hac >= 10 and hac > 0 and (o / float(hac)) >= 0.35 and n <= 1

def _to_parent_mol(m: Chem.Mol) -> Optional[Chem.Mol]:
    try:
        if _HAS_STD and _std is not None:
            p = _std.ChargeParent(_std.LargestFragmentChooser(preferOrganic=True).choose(_std.Cleanup(m)))
            if _looks_like_buffer_salt(p) or _matches_counterion(p) or _is_polyacidic_buffer_like(p) or _buffer_like_by_counts_from_mol(p):
                p = _std.ChargeParent(max(Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True), key=lambda x: x.GetNumHeavyAtoms()))
        else:
            p = _REM.StripMol(m, dontRemoveEverything=True)
            p = max(Chem.GetMolFrags(p, asMols=True, sanitizeFrags=True), key=lambda x: x.GetNumHeavyAtoms())

        Chem.SanitizeMol(p)
        try:
            Chem.SetAromaticity(p, Chem.AromaticityModel.AROMATICITY_RDKIT)
        except Exception:
            pass

        if not any(a.GetSymbol() == "C" for a in p.GetAtoms()): return None
        if _is_polyacidic_buffer_like(p) or _buffer_like_by_counts_from_mol(p): return None
        if p.GetNumHeavyAtoms() < MIN_PARENT_HEAVY: return None
        if _looks_like_buffer_salt(p) or _matches_counterion(p): return None
        return p
    except Exception:
        return None

# =========================
# Pre-protonation (Open Babel) & CLI builders
# =========================

def preprotonate_with_obabel(in_path: Path, out_path: Path, obabel_exe: str, ph: float) -> Path:
    """
    Pre-protonate the ligand with Open Babel at the requested pH, writing out_path.
    Returns out_path if successful; otherwise returns in_path unchanged.
    """
    try:
        fmt_in  = in_path.suffix.lower().lstrip(".") or "sdf"
        fmt_out = out_path.suffix.lower().lstrip(".") or fmt_in
        cmd = [obabel_exe, f"-i{fmt_in}", get_short_path_name(str(in_path.resolve())),
                          f"-o{fmt_out}", "-O", get_short_path_name(str(out_path.resolve())),
                          "-p", str(ph)]
        logging.info("Pre-protonating via Open Babel: %s", " ".join(cmd))
        cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=300)
        if cp.returncode == 0 and out_path.exists() and out_path.stat().st_size > 0:
            logging.info("Open Babel pH pre-protonation ok")
            return out_path
        logging.warning("Open Babel pH pre-protonation failed (rc=%s):\nSTDOUT:\n%s\nSTDERR:\n%s",
                        cp.returncode, cp.stdout, cp.stderr)
    except subprocess.TimeoutExpired:
        logging.warning("Open Babel pH pre-protonation timed out")
    except Exception as e:
        logging.warning("Open Babel pH pre-protonation error: %s", e)
    return in_path

def build_ligprep_cmd(cfg: dict, in_path: Path, out_pdbqt: Path) -> List[str]:
    """
    Choose the ligand-prep CLI based on config (Meeko vs ADT).
    - Meeko (mk_prepare_ligand.py): supports --charge_model {gasteiger,espaloma,zero}
    - ADT   (prepare_ligand4.py):   classic flags preserved
    """
    use_meeko = _is_true(cfg.get("USE_MEEKO", False))
    py  = cfg.get("MGLTOOLS_PYTHON") or sys.executable
    tool = cfg.get("PREPARE_LIGAND_SCRIPT") or ""
    if not tool:
        raise RuntimeError("Missing PREPARE_LIGAND_SCRIPT in config.txt")

    py_short   = get_short_path_name(py)
    tool_short = get_short_path_name(tool)

    if use_meeko:
        # NOTE: DO NOT pass any pH flag to Meeko; do pH work before this call if desired.
        charge_model = (cfg.get("MEEKO_CHARGE_MODEL") or "gasteiger").strip().lower()
        if charge_model not in {"gasteiger", "espaloma", "zero"}:
            charge_model = "gasteiger"
        args = [py_short, tool_short, "-i", str(in_path), "-o", str(out_pdbqt), "--charge_model", charge_model]
        if _is_true(cfg.get("MEEKO_BAD_CHARGE_OK", False)):
            args.append("--bad_charge_ok")
        return args

    # ADT path
    return [py_short, tool_short, "-l", str(in_path), "-o", str(out_pdbqt), "-U", "nphs_lps", "-A", "checkhydrogens"]

# =========================
# OBabel helpers
# =========================

def _write_obabel_friendly_sdf(mol: Chem.Mol, out_path: Path) -> bool:
    """Write SDF robustly (try kekulized and aromatic models)."""
    try:
        m = Chem.Mol(mol); Chem.SanitizeMol(m); Chem.Kekulize(m, clearAromaticFlags=True)
        w = Chem.SDWriter(str(out_path)); w.write(m); w.close()
        return True
    except Exception:
        try:
            m2 = Chem.Mol(mol); Chem.SanitizeMol(m2)
            Chem.SetAromaticity(m2, Chem.AromaticityModel.AROMATICITY_RDKIT)
            w = Chem.SDWriter(str(out_path))
            try: w.SetKekulize(False)
            except Exception: pass
            w.write(m2); w.close()
            return True
        except Exception:
            return False

def _run_obabel(cmd: List[str], timeout_sec: int) -> bool:
    print("Running OpenBabel:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, timeout=timeout_sec)
        return True
    except subprocess.TimeoutExpired:
        print("OpenBabel timed out.")
        return False
    except subprocess.CalledProcessError as e:
        print("OpenBabel failed with code:", e.returncode)
        return False

def _attempt_obabel_series(base_cmd: List[str], timeout_sec: int, threads_list: List[int], use_fast_first: bool = True) -> bool:
    attempts: List[List[str]] = []
    if use_fast_first:
        for j in threads_list:
            attempts.append(base_cmd + ["--fast", "-j", str(j)])
    for j in threads_list:
        cmd = [c for c in base_cmd if c != "--fast"] + ["-j", str(j)]
        attempts.append(cmd)

    for i, cmd in enumerate(attempts, 1):
        print(f"[OBabel attempt {i}/{len(attempts)}]")
        if _run_obabel(cmd, timeout_sec):
            return True
        print("Retrying with a more conservative setting…")
    return False

# =========================
# SDF splitting & conversion
# =========================

def count_sdf_records(sdf_path: Path) -> int:
    n = 0
    with open(sdf_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith("$$$$"):
                n += 1
    return n

def split_sdf_into_chunks(sdf_path: Path, out_dir: Path, chunk_size: int = 1000) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    chunks: List[Path] = []
    idx, mols_in_chunk = 0, 0
    mol_buf: List[str] = []

    def flush_chunk():
        nonlocal idx, mol_buf, mols_in_chunk
        if not mol_buf: return
        idx += 1
        out_path = out_dir / f"{sdf_path.stem}_chunk{idx:04d}.sdf"
        with open(out_path, "w", encoding="utf-8") as out:
            out.write("".join(mol_buf))
        chunks.append(out_path)
        mol_buf, mols_in_chunk = [], 0

    with open(sdf_path, "r", errors="ignore") as fh:
        cur: List[str] = []
        for line in fh:
            cur.append(line)
            if line.startswith("$$$$"):
                mol_buf.extend(cur); cur = []
                mols_in_chunk += 1
                if mols_in_chunk >= chunk_size:
                    flush_chunk()
        if cur: mol_buf.extend(cur)
        flush_chunk()

    return chunks

def convert_sdf_to_mol2_split_parallel(
    sdf_path: Path, mol2_output_dir: Path, obabel_exe: str, threads: int = 32,
    timeout_sec: int = 36000, chunk_size: int = 1000
) -> List[Path]:
    mol2_output_dir.mkdir(parents=True, exist_ok=True)
    chunk_dir = mol2_output_dir / "_sdf_chunks"
    if chunk_dir.exists(): shutil.rmtree(chunk_dir)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    print(f"SDF has ~{count_sdf_records(sdf_path)} molecules")
    chunk_paths = split_sdf_into_chunks(sdf_path, chunk_dir, chunk_size=chunk_size)
    print(f"Split into {len(chunk_paths)} chunk(s) of up to {chunk_size} molecules")

    all_out: List[Path] = []
    t1 = threads if threads >= 16 else 16
    t2 = max(16, t1 // 2)
    threads_list = []
    for t in (t1, t2, 16):
        if t not in threads_list:
            threads_list.append(t)

    for ci, chunk in enumerate(chunk_paths, 1):
        prefix = mol2_output_dir / f"mol2_chunk{ci:04d}_"
        base_cmd = [obabel_exe, "-isdf", str(chunk), "--gen3d", "-omol2", "-m", "-O", str(prefix) + ".mol2"]
        ok = _attempt_obabel_series(base_cmd, timeout_sec=timeout_sec, threads_list=threads_list, use_fast_first=True)
        if not ok:
            print(f"Chunk {ci} failed entirely; moving on.")
            continue
        out_files = sorted(mol2_output_dir.glob(f"mol2_chunk{ci:04d}_*.mol2"))
        print(f"Chunk {ci}: wrote {len(out_files)} mol2 files")
        all_out.extend(out_files)

    print(f"Total MOL2 files: {len(all_out)}")
    return all_out

# =========================
# RDKit ETKDG 3D generation
# =========================

def rdkit_embed_sdf_to_mol2(
    sdf_in: Path, mol2_out_dir: Path, obabel_exe: str, max_workers: int = 8
) -> List[Path]:
    from rdkit.Chem import AllChem
    suppl = Chem.SDMolSupplier(str(sdf_in), removeHs=False, sanitize=False)
    mols = [(i, m) for i, m in enumerate(suppl) if m is not None]
    print(f"RDKit: loaded {len(mols)} molecules from {sdf_in.name}")

    sdf_tmp_dir = mol2_out_dir / "_rdkit_embedded_sdf"
    sdf_tmp_dir.mkdir(parents=True, exist_ok=True)

    def _embed_one(i_m):
        i, m = i_m
        try:
            Chem.SanitizeMol(m)
        except Exception:
            return None
        orig_heavy = m.GetNumHeavyAtoms()
        p = _to_parent_mol(m)
        if p is None:
            _log_malformed(Path(f"rdk_{i:07d}"), "no_parent_or_too_small_after_desalting")
            return None
        if p is not m and p.GetNumHeavyAtoms() != orig_heavy:
            logging.info(f"[parent-pick] rdk_{i:07d}: {orig_heavy}→{p.GetNumHeavyAtoms()} heavy atoms")
        m = p

        try:
            params = AllChem.ETKDGv3(); params.randomSeed = 0xC0FFEE
            ok = AllChem.EmbedMolecule(m, params)
            if ok != 0: return None
            try: AllChem.UFFOptimizeMolecule(m, maxIters=200)
            except Exception: pass
            out = sdf_tmp_dir / f"rdk_{i:07d}.sdf"
            if _write_obabel_friendly_sdf(m, out):
                return out
            _log_malformed(out, "write_obabel_friendly_sDF_failed")
            return None
        except Exception:
            return None

    paths: List[Path] = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        for res in ex.map(_embed_one, mols):
            if res: paths.append(res)

    print(f"RDKit: embedded {len(paths)} molecules; converting to MOL2 …")
    out_files: List[Path] = []
    for pth in paths:
        out = mol2_out_dir / (pth.stem + ".mol2")
        if RESUME_SKIP and out.exists() and out.stat().st_size > 100:
            out_files.append(out); continue
        if _run_obabel([obabel_exe, "-isdf", str(pth), "-omol2", "-O", str(out)], timeout_sec=120):
            out_files.append(out)
    print(f"RDKit: wrote {len(out_files)} MOL2 files")
    return out_files

def convert_sdf_to_mol2_split_parallel_via_rdkit(
    sdf_path: Path, mol2_output_dir: Path, obabel_exe: str, max_workers: int
) -> List[Path]:
    embedded = rdkit_embed_sdf_to_mol2(sdf_in=sdf_path, mol2_out_dir=mol2_output_dir, obabel_exe=obabel_exe, max_workers=max_workers)
    return embedded or []

# =========================
# Validation & simple detectors
# =========================

def _looks_like_monoatomic_ion_pdbqt(lines: List[str]) -> bool:
    # very small helper for is_valid_ligand; relies on last token mapping
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    if len(atom_lines) <= 2:
        elems = [_element_from_adt(ln.split()[-1]) for ln in atom_lines]
        if len(set(e for e in elems if e != "H")) == 1 and any(e in MONOATOMIC_IONS for e in elems):
            return True
    return False

def _is_probable_water_from_pdbqt_lines(lines: List[str]) -> bool:
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    if not atom_lines: return False
    if any(" HOH " in ln for ln in atom_lines): return True
    heavy, has_O = 0, False
    for ln in atom_lines:
        parts = ln.split()
        if not parts: continue
        elem = _element_from_adt(parts[-1])
        if elem == "O": has_O = True
        if elem != "H": heavy += 1
    return has_O and heavy <= 1

def is_valid_ligand(path: Path, log_dir: Path) -> bool:
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        atom_lines   = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
        torsion_lines= [ln for ln in lines if ln.startswith("TORSDOF")]
        if not atom_lines: raise ValueError("No ATOM or HETATM lines found.")
        if not torsion_lines: raise ValueError("No torsion info (TORSDOF) found.")
        if _looks_like_monoatomic_ion_pdbqt(lines): raise ValueError("monoatomic_counterion")
        if _is_probable_water_from_pdbqt_lines(lines): raise ValueError("looks_like_water_or_hydroxide")
        heavy = sum(1 for ln in atom_lines if _element_from_adt(ln.split()[-1]) != "H")
        if heavy < 3: raise ValueError(f"too_few_heavy_atoms_in_pdbqt({heavy})")
        return True
    except Exception as e:
        malformed_log = Path(log_dir) / MALFORMED_LOG.name
        with open(malformed_log, "a", encoding="utf-8") as f:
            f.write(f"{path.name} - PDBQT validation failed: {e}\n")
        return False

def _count_aromatic_atoms_in_mol2(mol2_path: Path) -> int:
    m = Chem.MolFromMol2File(str(mol2_path), sanitize=True, removeHs=False)
    return -1 if m is None else sum(int(a.GetIsAromatic()) for a in m.GetAtoms())

def _count_aromatic_ad_types_in_pdbqt(pdbqt_path: Path) -> int:
    n = 0
    try:
        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if ln.startswith(("ATOM", "HETATM")):
                    t = ln.split()[-1].upper()
                    if t in ("A", "NA"): n += 1
    except Exception:
        return -1
    return n

def _pdbqt_from_mol2_via_obabel(mol2_file: Path, out_pdbqt: Path, obabel_exe: str) -> bool:
    tmp = out_pdbqt.with_suffix(".obabel_tmp.pdbqt")
    cmd = [obabel_exe, "-imol2", str(mol2_file), "-opdbqt", "-O", str(tmp)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
    except subprocess.CalledProcessError:
        return False
    if tmp.exists() and tmp.stat().st_size > 100:
        try: tmp.replace(out_pdbqt); return True
        except Exception: return False
    return False

# =========================
# Preparation primitive
# =========================

def _prepare_one(mol2_file: Path, pdbqt_path: Path, cfg: dict, obabel_exe: Optional[str] = None) -> Tuple[str, str]:
    """
    Prepare a single ligand (MOL2 as input for ADT; Meeko can handle PDB/SDF/MOL2).
    - Optional pre-protonation via Open Babel if LIGAND_PH is set.
    - Clear logging of command and outputs.
    - ADT aromaticity-rescue preserved (tries alt flags or obabel conversion as backup).
    """
    use_meeko = _is_true(cfg.get("USE_MEEKO", False))
    try:
        if pdbqt_path.exists():
            pdbqt_path.unlink()
    except Exception:
        pass

    # Baseline aromaticity (ADT only; Meeko not needed)
    src_arom = _count_aromatic_atoms_in_mol2(mol2_file) if (not use_meeko and mol2_file.suffix.lower()==".mol2") else -1

    # Optional pre-protonation BEFORE Meeko/ADT (never pass pH to Meeko)
    prep_input = mol2_file
    ligand_ph = cfg.get("LIGAND_PH")
    if ligand_ph and obabel_exe:
        try:
            ph_val = float(ligand_ph)
            preprot = mol2_file.with_suffix(".pp.mol2") if mol2_file.suffix.lower() == ".mol2" else mol2_file.with_suffix(".pp" + mol2_file.suffix)
            prep_input = preprotonate_with_obabel(mol2_file, preprot, obabel_exe, ph_val)
        except Exception as e:
            logging.warning("Ignoring LIGAND_PH (%r) due to parse/error: %s", ligand_ph, e)

    # Build command and run
    cmd = build_ligprep_cmd(cfg, prep_input, pdbqt_path)
    logging.info("Ligand prep cmd: %s", " ".join(map(str, cmd)))
    try:
        cp = subprocess.run([str(x) for x in cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            cwd=str(mol2_file.parent), timeout=900)
        logging.log(logging.INFO if cp.returncode==0 else logging.WARNING,
                    "rc=%s\nSTDOUT:\n%s\nSTDERR:\n%s", cp.returncode, cp.stdout or "", cp.stderr or "")
    except subprocess.TimeoutExpired:
        return (mol2_file.name, "timeout")
    except subprocess.CalledProcessError as e:
        logging.warning("Ligand prep raised CalledProcessError (rc=%s)", e.returncode)

    # ADT-only aromaticity rescue path if initial prep succeeded but aromaticity looks off
    if not use_meeko and pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
        adt_arom = _count_aromatic_ad_types_in_pdbqt(pdbqt_path)
        if src_arom >= 0 and adt_arom >= 0 and adt_arom < src_arom:
            logging.warning(f"[arom-mismatch] {mol2_file.name}: MOL2_arom={src_arom} > PDBQT_arom={adt_arom} (trying rescue)")
            alt_pdbqt = pdbqt_path.with_suffix(".alt.pdbqt")
            alt_cmd = build_ligprep_cmd(cfg, prep_input, alt_pdbqt)
            # tweak: try without removing lone pairs if using ADT
            if "-U" in alt_cmd:
                idx = alt_cmd.index("-U"); alt_cmd[idx+1] = "nphs"
            try:
                cp2 = subprocess.run([str(x) for x in alt_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=600)
                logging.log(logging.INFO if cp2.returncode==0 else logging.WARNING,
                            "ALT rc=%s\nSTDOUT:\n%s\nSTDERR:\n%s", cp2.returncode, cp2.stdout or "", cp2.stderr or "")
                alt_arom = _count_aromatic_ad_types_in_pdbqt(alt_pdbqt)
            except Exception:
                alt_arom = -1

            ob_arom = -1
            if obabel_exe:
                ob_tmp = pdbqt_path.with_suffix(".ob.pdbqt")
                if _pdbqt_from_mol2_via_obabel(mol2_file, ob_tmp, obabel_exe):
                    ob_arom = _count_aromatic_ad_types_in_pdbqt(ob_tmp)

            candidates = [(adt_arom, pdbqt_path)]
            if alt_arom >= 0: candidates.append((alt_arom, alt_pdbqt))
            if ob_arom  >= 0: candidates.append((ob_arom,  pdbqt_path.with_suffix(".ob.pdbqt")))
            best_arom, best_path = max(candidates, key=lambda t: t[0])

            if best_path != pdbqt_path and Path(best_path).exists():
                try: pdbqt_path.unlink(missing_ok=True)
                except Exception: pass
                Path(best_path).replace(pdbqt_path)

            # cleanup
            for tmp in (alt_pdbqt, pdbqt_path.with_suffix(".ob.pdbqt")):
                try:
                    if tmp.exists(): tmp.unlink()
                except Exception: pass

    # Final validation / quarantine
    if not is_valid_ligand(pdbqt_path, log_dir=pdbqt_path.parent):
        quarantine = pdbqt_path.parent / QUARANTINE_DIRNAME
        quarantine.mkdir(exist_ok=True)
        try:
            if pdbqt_path.exists():
                pdbqt_path.replace(quarantine / pdbqt_path.name)
        except Exception:
            pass
        return (mol2_file.name, "postcheck_fail")

    return (mol2_file.name, "ok")

# =========================
# Crystal ligand path (PDB → [optional MOL2] → PDBQT)
# =========================

def prep_ligands_from_pdb(ligand_output_dir: Path, ligands_mol2_dir: Path, prepped_ligands_dir: Path):
    """
    Prepare ligands extracted as individual PDBs from a PDB structure (crystal-safe path).
    Uses Meeko directly on sanitized PDB if USE_MEEKO, else PDB -> MOL2 (OpenBabel) -> ADT.
    """
    logging.info("Starting ligand preparation from PDB files (crystal-safe)")

    cfg = read_config()
    if not cfg.get("PREPARE_LIGAND_SCRIPT"):
        raise RuntimeError("Missing PREPARE_LIGAND_SCRIPT in config.txt")
    if not cfg.get("OPENBABEL_PATH"):
        raise RuntimeError("Missing OPENBABEL_PATH in config.txt")

    obabel_exe = resolve_obabel_exe(cfg["OPENBABEL_PATH"])
    _maybe_set_babel_datadir(obabel_exe)
    use_meeko = _is_true(cfg.get("USE_MEEKO", False))

    prepped_ligands_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = list(ligand_output_dir.glob("*.pdb"))
    logging.info(f"Found {len(pdb_files)} PDB ligand file(s)")

    def _sanitize_pdb(in_pdb: Path, out_pdb: Path) -> bool:
        wrote_any = False
        with open(in_pdb, "r", encoding="utf-8", errors="ignore") as fin, open(out_pdb, "w", encoding="utf-8") as fout:
            for ln in fin:
                if not ln.startswith(("ATOM", "HETATM")):
                    fout.write(ln); continue
                altloc = ln[16].strip() if len(ln) > 16 else ""
                if altloc and altloc.upper() not in ("", "A"):
                    continue
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                    if any(math.isnan(v) for v in (x, y, z)) or max(abs(x), abs(y), abs(z)) > 1e6:
                        continue
                except Exception:
                    continue
                fout.write(ln); wrote_any = True
        if not wrote_any:
            logging.warning(f"[control-prep] {in_pdb.name}: no safe ATOM/HETATM lines kept.")
        return wrote_any and out_pdb.exists() and out_pdb.stat().st_size > 0

    for pdb_file in pdb_files:
        logging.info(f"Processing: {pdb_file.name}")

        residue_name = pdb_file.stem.split("_")[0].upper()
        if residue_name in STANDARD_AMINO_ACIDS:
            logging.info(f"Skipping standard amino acid residue: {pdb_file.name}")
            continue
        if residue_name in EXCLUDE_CRYSTAL_ADDITIVES:
            _log_malformed(pdb_file, f"excluded_crystal_additive:{residue_name}")
            logging.info(f"Skipping crystallization additive: {pdb_file.name}")
            continue

        with open(pdb_file, "r", encoding="utf-8", errors="ignore") as f:
            atom_lines = [ln for ln in f if ln.startswith(("HETATM", "ATOM"))]
        if len(atom_lines) < MIN_ATOMS_FOR_DOCKING:
            _log_malformed(pdb_file, f"tiny_ligand_fewer_than_{MIN_ATOMS_FOR_DOCKING}_atoms")
            logging.info(f"Skipping tiny ligand: {pdb_file.name}")
            continue

        try:
            fix_pdb_elements(str(pdb_file))
        except Exception:
            pass

        sanitized = pdb_file.with_suffix(".sanitized.pdb")
        if not _sanitize_pdb(pdb_file, sanitized):
            _log_malformed(pdb_file, "sanitize_kept_no_atoms")
            logging.warning(f"Sanitization produced no atoms: {pdb_file.name}")
            sanitized.unlink(missing_ok=True)
            continue

        # quick counter-ion/buffer checks
        flag_buffer, reason = False, None
        try:
            m_chk = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
            if m_chk is not None:
                if _buffer_like_by_counts_from_mol(m_chk):
                    flag_buffer, reason = True, "buffer_like_by_counts"
                if not flag_buffer:
                    try: Chem.SanitizeMol(m_chk)
                    except Exception: pass
                    reason = _matches_counterion(m_chk) or (_looks_like_buffer_salt(m_chk) and "buffer_like")
                    if reason: flag_buffer = True
                    if not flag_buffer and residue_name in {"UNL", "LIG"} and _is_polyacidic_buffer_like(m_chk):
                        flag_buffer, reason = True, "polyacidic_buffer_like"
        except Exception:
            if _buffer_like_by_counts_from_pdbfile(sanitized) or (residue_name in {"UNL", "LIG"} and _polyacidic_by_counts_from_pdbfile(sanitized)):
                flag_buffer, reason = True, "counts_only_polyacidic"

        if flag_buffer:
            _log_malformed(pdb_file, f"counterion_or_buffer:{reason or 'unknown'}")
            logging.info(f"Skipping likely counter-ion/buffer ({reason or 'unknown'}): {pdb_file.name}")
            sanitized.unlink(missing_ok=True)
            continue

        # save reference SDF (optional)
        try:
            m = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
            ref_dir = sanitized.parent.parent / "reference"; ref_dir.mkdir(parents=True, exist_ok=True)
            ref_sdf = ref_dir / (pdb_file.stem + ".sdf")
            w = Chem.SDWriter(str(ref_sdf))
            try: w.SetKekulize(False)
            except Exception: pass
            if m is not None: w.write(m)
            w.close()
        except Exception:
            pass

        pdbqt_path = prepped_ligands_dir / f"{pdb_file.stem}.pdbqt"
        if RESUME_SKIP and pdbqt_path.exists() and pdbqt_path.stat().st_size > 100 and is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir):
            logging.info(f"[resume] Valid PDBQT already exists, skipping: {pdbqt_path.name}")
            sanitized.unlink(missing_ok=True)
            continue

        # Decide prep input (Meeko can take PDB; ADT prefers MOL2)
        prep_input = sanitized
        tmp_mol2 = sanitized.with_suffix(".mol2")
        if not use_meeko:
            ob_cmd = [obabel_exe, "-ipdb", get_short_path_name(str(sanitized.resolve())), "-omol2",
                      "-O", get_short_path_name(str(tmp_mol2.resolve()))]
            logging.info("Converting sanitized PDB -> MOL2 (no gen3d): " + " ".join(map(str, ob_cmd)))
            try:
                res_ob = subprocess.run(ob_cmd, check=True, capture_output=True, text=True, timeout=300)
                logging.info(res_ob.stdout)
                prep_input = tmp_mol2
            except subprocess.CalledProcessError as e:
                logging.warning(f"OBabel PDB->MOL2 failed for {sanitized.name}:\n{e.stderr}")
                prep_input = sanitized  # if USE_MEEKO toggled, still okay; otherwise fallback later

        # Optional pre-protonation BEFORE prep
        ligand_ph = cfg.get("LIGAND_PH")
        if ligand_ph and obabel_exe:
            try:
                ph_val = float(ligand_ph)
                pp = prep_input.with_suffix(".pp" + prep_input.suffix)
                prep_input = preprotonate_with_obabel(prep_input, pp, obabel_exe, ph_val)
            except Exception as e:
                logging.warning("Ignoring LIGAND_PH for crystal ligand due to error: %s", e)

        # Run prep
        mgl_ok = False
        try:
            cmd = build_ligprep_cmd(cfg, prep_input, pdbqt_path)
            logging.info("Preparing ligand: " + " ".join(map(str, cmd)))
            result = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=str(sanitized.parent), timeout=900)
            logging.info(result.stdout)
            mgl_ok = True
        except subprocess.TimeoutExpired:
            logging.error(f"Timeout preparing {prep_input.name}")
        except subprocess.CalledProcessError as e:
            logging.warning(f"Ligand prep failed for {prep_input.name} (will try fallback):\n{e.stderr}")

        needs_fallback = (
            not mgl_ok
            or not pdbqt_path.exists()
            or pdbqt_path.stat().st_size < 100
            or not is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir)
        )

        if needs_fallback:
            logging.info(f"Primary prep failed/invalid for {sanitized.name}; trying RDKit→SDF→OBabel→MOL2→prep.")
            tmp_sdf = sanitized.with_suffix(".tmp.sdf")
            try:
                mol = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
                if mol is None: raise RuntimeError("RDKit failed to read sanitized PDB")

                try: Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
                except Exception: pass
                try: mol = Chem.AddHs(mol, addCoords=True)
                except Exception: pass

                if not _write_obabel_friendly_sdf(mol, tmp_sdf):
                    raise RuntimeError("Failed to write OBabel-friendly SDF in fallback")

                ob_cmd2 = [obabel_exe, "-isdf", get_short_path_name(str(tmp_sdf.resolve())),
                           "-omol2", "-O", get_short_path_name(str(tmp_mol2.resolve()))]
                logging.info("Converting SDF -> MOL2 (no gen3d): " + " ".join(map(str, ob_cmd2)))
                try:
                    res_ob2 = subprocess.run(ob_cmd2, check=True, capture_output=True, text=True, timeout=300)
                    logging.info(res_ob2.stdout)
                except subprocess.CalledProcessError as e:
                    logging.error(f"OBabel SDF->MOL2 fallback failed for {sanitized.name}:\n{e.stderr}")

                # prep again using the builder
                if tmp_mol2.exists() and tmp_mol2.stat().st_size > 100:
                    prepare_cmd2 = build_ligprep_cmd(cfg, tmp_mol2, pdbqt_path)
                    logging.info("Preparing ligand (fallback MOL2): " + " ".join(map(str, prepare_cmd2)))
                    try:
                        res2 = subprocess.run(prepare_cmd2, check=True, capture_output=True, text=True, timeout=900)
                        logging.info(res2.stdout)
                    except subprocess.CalledProcessError as e:
                        logging.error(f"Fallback ligand prep failed for {sanitized.name}:\n{e.stderr}")
                    except subprocess.TimeoutExpired:
                        logging.error(f"Fallback timeout for {sanitized.name}")
            finally:
                try: tmp_sdf.unlink(missing_ok=True)
                except Exception: pass

        try: sanitized.unlink(missing_ok=True)
        except Exception: pass
        try:
            if pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
                tmp_mol2.unlink(missing_ok=True)
        except Exception:
            pass

        if not pdbqt_path.exists() or pdbqt_path.stat().st_size < 100 or not is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir):
            _log_malformed(pdbqt_path, "pdbqt_postcheck_fail_or_small")
            quarantine = prepped_ligands_dir / QUARANTINE_DIRNAME
            quarantine.mkdir(exist_ok=True)
            try:
                if pdbqt_path.exists():
                    pdbqt_path.replace(quarantine / pdbqt_path.name)
            except Exception:
                pass
            continue

        logging.info(f"Created PDBQT: {pdbqt_path.name}")

# =========================
# Batch SDF → MOL2 → PDBQT
# =========================

def prep_ligands():
    """Main batch mode: iterate SDF files, produce MOL2 with RDKit/OBabel, then prep to PDBQT in parallel."""
    print("Starting ligand preparation")

    cfg = read_config()
    if not cfg.get("PREPARE_LIGAND_SCRIPT"):
        raise RuntimeError("Missing PREPARE_LIGAND_SCRIPT in config.txt")
    if not cfg.get("OPENBABEL_PATH"):
        raise RuntimeError("Missing OPENBABEL_PATH in config.txt")

    ligand_extracted_dir = Path(cfg["LIGAND_EXTRACTED_DIR"]).resolve()
    ligands_mol2_dir     = Path(cfg["LIGANDS_MOL2_DIR"]).resolve()
    output_ligands_dir   = Path(cfg["OUTPUT_LIGANDS_DIR"]).resolve()
    output_ligands_dir.mkdir(parents=True, exist_ok=True)

    obabel_exe = resolve_obabel_exe(cfg["OPENBABEL_PATH"])
    _maybe_set_babel_datadir(obabel_exe)

    # sanity on paths
    for label, p in [
        ("OPENBABEL_PATH",  obabel_exe),
        ("LIGAND_EXTRACTED_DIR", ligand_extracted_dir),
        ("LIGANDS_MOL2_DIR",     ligands_mol2_dir),
        ("OUTPUT_LIGANDS_DIR",   output_ligands_dir),
    ]:
        if not str(p).strip():
            raise RuntimeError(f"Config value missing/empty: {label}")

    sdf_files = list(ligand_extracted_dir.glob("*.sdf"))
    print(f"Found {len(sdf_files)} SDF file(s)")
    if not sdf_files:
        return

    for sdf_file in sdf_files:
        print(f"\n=== Processing SDF: {sdf_file.name} ===")
        sdf_abs = sdf_file.resolve()
        max_workers = max(1, min(8, os.cpu_count() or 8))

        if USE_RDKIT_FOR_3D:
            print("Using RDKit ETKDG for 3D with parent-picking; OpenBabel only for format conversion …")
            mol2_files = rdkit_embed_sdf_to_mol2(sdf_abs, ligands_mol2_dir, obabel_exe=obabel_exe, max_workers=max_workers)
        else:
            print("Using OpenBabel --gen3d; pre-cleaning SDF to parent-only …")
            cleaned_sdf = ligands_mol2_dir / (sdf_abs.stem + "_parents.sdf")
            # implement _write_parent_only_sdf if you want this path; omitted for brevity
            # n_kept = _write_parent_only_sdf(sdf_abs, cleaned_sdf)
            # if n_kept == 0: continue
            mol2_files = convert_sdf_to_mol2_split_parallel(sdf_abs, ligands_mol2_dir, obabel_exe,
                                                            threads=OBABEL_THREADS, timeout_sec=OBABEL_TIMEOUT_S, chunk_size=CHUNK_SIZE)

        if not mol2_files:
            print("No MOL2 files produced; skipping this SDF.")
            continue

        print(f"Preparing {len(mol2_files)} MOL2 files (parallel)…")
        ok_count = fail_count = 0
        futures = []
        with ThreadPoolExecutor(max_workers=max(1, min(8, os.cpu_count() or 8))) as ex:
            for mol2_file in mol2_files:
                pdbqt_path = output_ligands_dir / f"{mol2_file.stem}.pdbqt"
                if RESUME_SKIP and pdbqt_path.exists() and pdbqt_path.stat().st_size > 100 and is_valid_ligand(pdbqt_path, log_dir=output_ligands_dir):
                    print(f"[resume] Skipping already-valid PDBQT: {pdbqt_path.name}")
                    continue
                futures.append(ex.submit(_prepare_one, mol2_file, pdbqt_path, cfg, obabel_exe))

            total = len(futures)
            for i, fut in enumerate(as_completed(futures), 1):
                name, status = fut.result()
                if status == "ok": ok_count += 1
                else: fail_count += 1
                if i % 100 == 0 or status != "ok":
                    print(f"[{i}/{total}] {name}: {status}")

        print(f"Done: {ok_count} ok, {fail_count} failed/quarantined for {sdf_file.name}")

# =========================
# Standalone execution
# =========================

if __name__ == "__main__":
    prep_ligands()
