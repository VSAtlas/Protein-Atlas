import sys
import ctypes
import os
from ctypes import wintypes, create_unicode_buffer
from pathlib import Path
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Tuple, Dict, Any

from activesite import fix_pdb_elements

# --- RDKit / Standardization imports ---
from rdkit import Chem
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
    _HAS_STD = True
except Exception:
    _HAS_STD = False

# Ensure Open Babel can find its data (adjust if needed)
os.environ.setdefault("BABEL_DATADIR", "E:/OpenBabel-3.1.1/data")

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "SEC", "PYL", "MSE"
}

# --- Salvage / logging config ---
MALFORMED_LOG = Path("malformed_ligands.txt")
QUARANTINE_DIRNAME = "quarantine"

ALLOWED_ELEMENTS = {"H","C","N","O","F","P","S","Cl","Br","I","B","Si","Se","Zn","Mg","Ca","Mn","Fe","K","Na"}
MAX_HEAVY_ATOMS = 1200
MIN_ATOMS_FOR_DOCKING = 5  # replaces inline "5" and logs a reason
# ----------- tuning switches -----------
USE_RDKIT_FOR_3D = True          # <-- turn RDKit ETKDG on by default
OBABEL_THREADS    = 16           # do not exceed this; stable on Windows
OBABEL_TIMEOUT_S  = 900          # per OBabel call
CHUNK_SIZE        = 200          # smaller micro-chunks improve parallelism
# --------------------------------------


def _log_malformed(path: Path, reason: str):
    try:
        MALFORMED_LOG.parent.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    with open(MALFORMED_LOG, "a", encoding="utf-8") as fh:
        fh.write(f"{path}\t{reason}\n")


def _quick_filters(mol):
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


def _standardize_then_sanitize(mol):
    if not _HAS_STD:
        return None, "std_module_missing"
    try:
        md = rdMolStandardize.MetalDisconnector()
        fr = rdMolStandardize.FragmentRemover()
        lf = rdMolStandardize.LargestFragmentChooser()
        uc = rdMolStandardize.Uncharger()

        m = md.Disconnect(mol)
        m = fr.RemoveFragments(m)
        m = lf.choose(m)
        m = uc.uncharge(m)
        Chem.SanitizeMol(m)
        ok, why = _quick_filters(m)
        if not ok:
            return None, why
        return m, ""
    except Exception as e:
        return None, f"std_resanitize_fail:{e}"


def _reserialize_mol_via_obabel(mol, obabel_exe_short: str, target_mol2: Path):
    """Write RDKit mol to SDF, then make a fresh MOL2 with Open Babel."""
    tmp_sdf = target_mol2.with_suffix(".std.sdf")
    try:
        w = Chem.SDWriter(str(tmp_sdf))
        w.write(mol)
        w.close()
    except Exception as e:
        _log_malformed(target_mol2, f"rdkit_sdf_write_fail:{e}")
        return None

    fresh_mol2 = target_mol2.with_suffix(".std.mol2")
    cmd = [obabel_exe_short, "-isdf", str(tmp_sdf), "--gen3d", "-omol2", "-O", str(fresh_mol2)]
    try:
        # Capture here is fine; this is a rare salvage path
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        _log_malformed(target_mol2, f"obabel_reserialize_fail:{getattr(e, 'stderr', '')[:200]}")
        try:
            tmp_sdf.unlink(missing_ok=True)
        except Exception:
            pass
        return None

    try:
        tmp_sdf.unlink(missing_ok=True)
    except Exception:
        pass

    if not fresh_mol2.exists() or fresh_mol2.stat().st_size < 100:
        _log_malformed(fresh_mol2, "fresh_mol2_empty_or_small")
        return None
    return fresh_mol2


def read_config(path="config.txt"):
    config: Dict[str, str] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                config[key.strip()] = value.strip()
    return config


def get_short_path_name(long_name):
    if sys.platform != 'win32':
        return long_name

    _GetShortPathNameW = ctypes.windll.kernel32.GetShortPathNameW
    _GetShortPathNameW.argtypes = [wintypes.LPCWSTR, wintypes.LPWSTR, wintypes.DWORD]
    _GetShortPathNameW.restype = wintypes.DWORD

    output_buf_size = 260
    output_buf = create_unicode_buffer(output_buf_size)
    ret = _GetShortPathNameW(long_name, output_buf, output_buf_size)

    if ret == 0 or ret > output_buf_size:
        return long_name
    return output_buf.value


# ---------- Robust OBabel execution ----------

def _run_obabel(cmd: List[str], timeout_sec: int) -> bool:
    """
    Run obabel with live console output (no capture) so you see progress bars.
    Returns True on success, False on non-zero exit or timeout.
    """
    print("Running Open Babel:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, timeout=timeout_sec)
        return True
    except subprocess.TimeoutExpired:
        print("Open Babel timed out.")
        return False
    except subprocess.CalledProcessError as e:
        print("Open Babel failed with code:", e.returncode)
        return False


def _attempt_obabel_series(base_cmd: List[str], timeout_sec: int, threads_list: List[int], use_fast_first: bool = True) -> bool:
    """
    Try a sequence of OBabel invocations:
      1) with --fast (if requested) at given thread counts
      2) without --fast at given thread counts
    Returns True on first success.
    """
    attempts: List[List[str]] = []
    # with --fast
    if use_fast_first:
        for j in threads_list:
            attempts.append(base_cmd + ["--fast", "-j", str(j)])
    # without --fast
    for j in threads_list:
        # ensure --fast not present
        cmd = [c for c in base_cmd if c != "--fast"]
        cmd += ["-j", str(j)]
        attempts.append(cmd)

    total = len(attempts)
    for i, cmd in enumerate(attempts, 1):
        print(f"[OBabel attempt {i}/{total}]")
        if _run_obabel(cmd, timeout_sec):
            return True
        print("Retrying with a more conservative setting...")
    return False


# ---------- SDF chunking ----------

def count_sdf_records(sdf_path: Path) -> int:
    """Counts '$$$$' records to estimate molecules in SDF."""
    n = 0
    with open(sdf_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith("$$$$"):
                n += 1
    return n


def split_sdf_into_chunks(sdf_path: Path, out_dir: Path, chunk_size: int = 1000) -> List[Path]:
    """
    Splits a large SDF into many smaller SDFs with at most `chunk_size` molecules each.
    Returns the list of chunk paths.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    chunks: List[Path] = []
    idx = 0
    mol_buf: List[str] = []
    mols_in_chunk = 0

    def flush_chunk():
        nonlocal idx, mol_buf, mols_in_chunk
        if not mol_buf:
            return
        idx += 1
        out_path = out_dir / f"{sdf_path.stem}_chunk{idx:04d}.sdf"
        with open(out_path, "w", encoding="utf-8") as out:
            out.write("".join(mol_buf))
        chunks.append(out_path)
        mol_buf = []
        mols_in_chunk = 0

    with open(sdf_path, "r", errors="ignore") as fh:
        cur: List[str] = []
        for line in fh:
            cur.append(line)
            if line.startswith("$$$$"):
                mol_buf.extend(cur)
                cur = []
                mols_in_chunk += 1
                if mols_in_chunk >= chunk_size:
                    flush_chunk()
        # trailing (malformed) last block without $$$$?
        if cur:
            mol_buf.extend(cur)
        flush_chunk()
    return chunks


# ---------- SDF → MOL2 with chunking + thread fallback (not below 16) ----------

def convert_sdf_to_mol2_split_parallel(
    sdf_path: Path,
    mol2_output_dir: Path,
    obabel_exe: str,
    threads: int = 32,
    timeout_sec: int = 36000,
    chunk_size: int = 1000
) -> List[Path]:
    """
    Robust SDF → per-molecule MOL2 using OBabel:
      • Splits huge SDF into smaller chunks (default 1000 mols)
      • Tries threads [threads, max(16, threads//2), 16] (never below 16)
      • Falls back from --fast to conservative if needed
    """
    mol2_output_dir.mkdir(parents=True, exist_ok=True)
    n_mols = count_sdf_records(sdf_path)
    print(f"SDF has ~{n_mols} molecules")

    # Chunk the SDF to avoid crashing the whole run on one bad molecule
    chunk_dir = mol2_output_dir / "_sdf_chunks"
    chunk_paths = split_sdf_into_chunks(sdf_path, chunk_dir, chunk_size=chunk_size)
    print(f"Split into {len(chunk_paths)} chunk(s) of up to {chunk_size} molecules")

    all_out: List[Path] = []

    # Thread attempts: never drop below 16
    t1 = threads if threads >= 16 else 16
    t2 = max(16, t1 // 2)
    threads_list = [t1, t2, 16]
    # Deduplicate while preserving order
    seen = set()
    threads_list = [t for t in threads_list if not (t in seen or seen.add(t))]

    for ci, chunk in enumerate(chunk_paths, 1):
        prefix = mol2_output_dir / f"mol2_chunk{ci:04d}_"
        base_cmd = [
            obabel_exe, "-isdf", str(chunk),
            "--gen3d", "-omol2",
            "-m", "-O", str(prefix) + ".mol2"
        ]
        ok = _attempt_obabel_series(
            base_cmd,
            timeout_sec=timeout_sec,
            threads_list=threads_list,
            use_fast_first=True
        )
        if not ok:
            print(f"Chunk {ci} failed entirely; moving on.")
            continue

        # Collect outputs for this chunk
        out_files = sorted(mol2_output_dir.glob(f"mol2_chunk{ci:04d}_*.mol2"))
        print(f"Chunk {ci}: wrote {len(out_files)} mol2 files")
        all_out.extend(out_files)

    print(f"Total MOL2 files: {len(all_out)}")
    return all_out
def _obabel_convert_chunk(chunk_sdf: Path, out_prefix: Path, obabel_exe: str) -> int:
    """
    Convert one SDF chunk to many MOL2s (no --gen3d).
    Returns number of MOL2 files written.
    """
    cmd = [
        obabel_exe, "-isdf", str(chunk_sdf),
        "-omol2", "-m",
        "-O", str(out_prefix) + ".mol2",
        "-j", str(OBABEL_THREADS)
    ]
    print("Running Open Babel:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, timeout=OBABEL_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        print(f"Open Babel timed out on {chunk_sdf.name}")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Open Babel failed ({e.returncode}) on {chunk_sdf.name}")
        return 0

    written = len(list(out_prefix.parent.glob(out_prefix.name + "_*.mol2")))
    return written


def convert_sdf_to_mol2_split_parallel_via_rdkit(
    sdf_path: Path, mol2_output_dir: Path, obabel_exe: str, max_workers: int
) -> List[Path]:
    """
    RDKit ETKDG -> SDF micro-chunks -> parallel OBabel SDF->MOL2 (no gen3d).
    """
    # 1) RDKit embed once to a temp dir of per-mol SDFs
    embedded = rdkit_embed_sdf_to_mol2(  # writes mol2 via obabel at the end
        sdf_in=sdf_path, mol2_out_dir=mol2_output_dir, obabel_exe=obabel_exe, max_workers=max_workers
    )
    if embedded:
        return embedded
    return []

# ---------- RDKit ETKDG fallback (no OBabel gen3d) ----------

def rdkit_embed_sdf_to_mol2(
    sdf_in: Path, mol2_out_dir: Path, obabel_exe: str, max_workers: int = 8
) -> List[Path]:
    """
    ETKDG embeds with RDKit (parallel), writes per-molecule SDF, then converts SDF→MOL2 without --gen3d.
    Much more stable than OBabel --gen3d for gnarly molecules.
    """
    from rdkit import Chem
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
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = 0xC0FFEE
            ok = AllChem.EmbedMolecule(m, params)
            if ok != 0:
                return None
            try:
                AllChem.UFFOptimizeMolecule(m, maxIters=200)
            except Exception:
                pass
            out = sdf_tmp_dir / f"rdk_{i:07d}.sdf"
            w = Chem.SDWriter(str(out))
            w.write(m)
            w.close()
            return out
        except Exception:
            return None

    paths: List[Path] = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        for res in ex.map(_embed_one, mols):
            if res:
                paths.append(res)

    print(f"RDKit: embedded {len(paths)} molecules; converting to MOL2 …")
    out_files: List[Path] = []
    for p in paths:
        out = mol2_out_dir / (p.stem + ".mol2")
        cmd = [obabel_exe, "-isdf", str(p), "-omol2", "-O", str(out)]
        if _run_obabel(cmd, timeout_sec=120):
            out_files.append(out)
    print(f"RDKit: wrote {len(out_files)} MOL2 files")
    return out_files


# ---------- MGLTools preparation ----------

def is_valid_ligand(path: Path, log_dir: Path) -> bool:
    """
    Minimal check that ensures the .pdbqt is not malformed.
    Returns True if basic structure seems intact.
    """
    try:
        with open(path, 'r') as f:
            lines = f.readlines()

        atom_lines = [line for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]
        torsion_lines = [line for line in lines if line.startswith("TORSDOF")]

        if not atom_lines:
            raise ValueError("No ATOM or HETATM lines found.")
        if not torsion_lines:
            raise ValueError("No torsion info (TORSDOF) found.")

        return True

    except Exception as e:
        malformed_log = Path(log_dir) / "malformed_ligands.txt"
        with open(malformed_log, "a", encoding="utf-8") as f:
            f.write(f"{path.name} - PDBQT validation failed: {e}\n")
        return False


def _prepare_one(mgltools_python_short: str, prepare_script_short: str, mol2_file: Path, pdbqt_path: Path) -> Tuple[str, str]:
    mol2_short = mol2_file.name  # since cwd will be set to mol2_file.parent
    pdbqt_short = get_short_path_name(str(pdbqt_path.resolve()))
    cmd = [
        mgltools_python_short,
        prepare_script_short,
        "-l", mol2_short,
        "-o", pdbqt_short,
        "-A", "hydrogens"
    ]
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,   # avoid ADT chatter flooding the console
            text=True,
            cwd=str(mol2_file.parent),
            timeout=600            # 10 minutes per ligand
        )
        # Post-check
        if not is_valid_ligand(pdbqt_path, log_dir=pdbqt_path.parent):
            quarantine = pdbqt_path.parent / QUARANTINE_DIRNAME
            quarantine.mkdir(exist_ok=True)
            try:
                pdbqt_path.replace(quarantine / pdbqt_path.name)
            except Exception:
                pass
            return (mol2_file.name, "postcheck_fail")
        return (mol2_file.name, "ok")
    except subprocess.TimeoutExpired:
        return (mol2_file.name, "timeout")
    except subprocess.CalledProcessError as e:
        return (mol2_file.name, f"prepare_fail:{e.returncode}")


def load_mol2_lenient(path, logger):
    mol = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
    if mol is None:
        _log_malformed(Path(path), "rdkit_read_fail")
        return None
    try:
        Chem.SanitizeMol(mol)
        ok, why = _quick_filters(mol)
        if not ok:
            _log_malformed(Path(path), why)
            return None
        return mol
    except Exception as e:
        std_mol, why = _standardize_then_sanitize(mol)
        if std_mol is None:
            _log_malformed(Path(path), f"sanitize_fail:{e}|{why}")
            if logger:
                logger.warning(f"RDKit failed to sanitize: {path} ({e})")
            return None
        return std_mol


# ---------- PDB ligand path (unchanged, but kept for completeness) ----------
def prep_ligands_from_pdb(ligand_output_dir: Path, ligands_mol2_dir: Path, prepped_ligands_dir: Path):
    """
    Prepare extracted crystal ligands directly from their PDB coordinates.
    - No OBabel --gen3d on crystals (avoids NaN explosions)
    - MGLTools prepare_ligand4 on sanitized PDB
    - Fallback: RDKit → write temp PDB → MGLTools (still no gen3d)
    """
    logging.info("Starting ligand preparation from PDB files (crystal-safe path)")

    config = read_config()
    mgltools_python = config.get("MGLTOOLS_PYTHON")
    mgltools_path = config.get("MGLTOOLS_PATH")
    obabel_exe_cfg = config.get("OPENBABEL_PATH")

    if not mgltools_python or not mgltools_path or not obabel_exe_cfg:
        raise RuntimeError("Missing paths in config.txt: MGLTOOLS_PYTHON, MGLTOOLS_PATH, OPENBABEL_PATH")

    # Accept either folder or direct exe for OpenBabel
    obabel_exe = obabel_exe_cfg if obabel_exe_cfg.lower().endswith(".exe") else str(Path(obabel_exe_cfg) / "obabel.exe")

    mgltools_python_short = get_short_path_name(mgltools_python)
    prepare_script = Path(mgltools_path) / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
    if not prepare_script.exists():
        raise FileNotFoundError(f"prepare_ligand4.py not found at {prepare_script}")
    prepare_script_short = get_short_path_name(str(prepare_script.resolve()))
    obabel_exe_short = get_short_path_name(obabel_exe)

    pdb_files = list(ligand_output_dir.glob("*.pdb"))
    logging.info(f"Found {len(pdb_files)} PDB ligand file(s)")

    prepped_ligands_dir.mkdir(parents=True, exist_ok=True)

    def _sanitize_pdb(in_pdb: Path, out_pdb: Path) -> bool:
        wrote_any = False
        with open(in_pdb, "r", encoding="utf-8", errors="ignore") as fin, \
             open(out_pdb, "w", encoding="utf-8") as fout:
            for ln in fin:
                if not ln.startswith(("ATOM", "HETATM")):
                    fout.write(ln); continue
                # keep only primary altloc (blank or 'A')
                altloc = ln[16].strip() if len(ln) > 16 else ""
                if altloc and altloc.upper() not in ("", "A"):
                    continue
                # filter NaNs / absurd coords
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                    if any([(x != x), (y != y), (z != z)]) or max(abs(x), abs(y), abs(z)) > 1e6:
                        continue
                except Exception:
                    continue
                fout.write(ln); wrote_any = True
        if not wrote_any:
            logging.warning(f"[control-prep] {in_pdb.name}: no safe ATOM/HETATM lines kept.")
        return wrote_any and out_pdb.exists() and out_pdb.stat().st_size > 0

    for pdb_file in pdb_files:
        logging.info(f"Processing: {pdb_file.name}")

        # Skip standard amino acids (peptides etc.)
        residue_name = pdb_file.stem.split("_")[0].upper()
        if residue_name in STANDARD_AMINO_ACIDS:
            logging.info(f"Skipping standard amino acid residue: {pdb_file.name}")
            continue

        # Quick atom count filter
        with open(pdb_file, "r", encoding="utf-8", errors="ignore") as f:
            atom_lines = [line for line in f if line.startswith(("HETATM", "ATOM"))]
        if len(atom_lines) < MIN_ATOMS_FOR_DOCKING:
            _log_malformed(pdb_file, f"tiny_ligand_fewer_than_{MIN_ATOMS_FOR_DOCKING}_atoms")
            logging.info(f"Skipping tiny ligand: {pdb_file.name}")
            continue

        # Fix element columns before any tool reads it
        try:
            fix_pdb_elements(str(pdb_file))
        except Exception:
            pass

        sanitized = pdb_file.with_suffix(".sanitized.pdb")
        # Ensure a pristine reference exists (for older runs / direct starts)
        ref_dir = sanitized.parent.parent / "reference"
        ref_dir.mkdir(parents=True, exist_ok=True)
        ref_sdf = ref_dir / (pdb_file.stem + ".sdf")
        if not ref_sdf.exists():
            try:
                # prefer RDKit; falls back to a stub if unavailable
                from rdkit import Chem
                m = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
                if m is not None:
                    w = Chem.SDWriter(str(ref_sdf));
                    w.write(m);
                    w.close()
            except Exception:
                # last-ditch: create an empty placeholder so downstream code has a path
                with open(ref_sdf, "w") as out:
                    out.write(pdb_file.stem + "\n$$$$\n")

        ok = _sanitize_pdb(pdb_file, sanitized)
        if not ok:
            _log_malformed(pdb_file, "sanitize_kept_no_atoms")
            logging.warning(f"Sanitization produced no atoms: {pdb_file.name}")
            sanitized.unlink(missing_ok=True)
            continue
        # write pristine reference right here (post-altLoc filter)
        ref_dir = sanitized.parent.parent / "reference"
        ref_dir.mkdir(parents=True, exist_ok=True)
        ref_sdf = ref_dir / (pdb_file.stem + ".sdf")
        if not ref_sdf.exists():
            try:
                m = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
                if m:
                    w = Chem.SDWriter(str(ref_sdf));
                    w.write(m);
                    w.close()
            except Exception:
                pass
        # Output PDBQT
        pdbqt_path = prepped_ligands_dir / f"{pdb_file.stem}.pdbqt"

        # 1) Primary path: MGLTools directly on sanitized PDB (no gen3d)
        prepare_cmd = [
            mgltools_python_short,
            prepare_script_short,
            "-l", get_short_path_name(str(sanitized.resolve())),
            "-o", get_short_path_name(str(pdbqt_path.resolve())),
            "-U", "nphs_lps",          # mild neutralization
            "-A", "checkhydrogens"
        ]
        logging.info(f"Preparing ligand (crystal coords): {' '.join(prepare_cmd)}")
        mgl_ok = False
        try:
            result = subprocess.run(
                prepare_cmd,
                check=True,
                capture_output=True,
                text=True,
                cwd=str(sanitized.parent),
                timeout=600
            )
            logging.info(result.stdout)
            mgl_ok = True
        except subprocess.TimeoutExpired:
            logging.error(f"Timeout preparing {sanitized.name}")
        except subprocess.CalledProcessError as e:
            logging.warning(f"MGLTools failed for {sanitized.name} (will try fallback):\n{e.stderr}")

        # 2) Fallback: only if primary MGLTools step failed or produced a bad/small PDBQT
        needs_fallback = (
            not mgl_ok
            or not pdbqt_path.exists()
            or pdbqt_path.stat().st_size < 100
            or not is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir)
        )

        if needs_fallback:
            logging.info(f"Primary prep failed/invalid for {sanitized.name}; trying RDKit→PDB→MGLTools (no gen3d).")
            tmp_pdb = sanitized.with_suffix(".tmp.pdb")
            try:
                mol = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
                if mol is None:
                    raise RuntimeError("RDKit failed to read sanitized PDB")

                # Keep sanitization light; we only need a reasonable PDB to pass to MGLTools
                try:
                    # avoid heavy sanitization; we just ensure the molecule object is usable
                    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
                except Exception:
                    pass

                # Add Hs with coords if possible (helps MGLTools hydrogens/charges)
                try:
                    mol = Chem.AddHs(mol, addCoords=True)
                except Exception:
                    pass

                Chem.MolToPDBFile(mol, str(tmp_pdb))

                prepare_cmd2 = [
                    mgltools_python_short,
                    prepare_script_short,
                    "-l", get_short_path_name(str(tmp_pdb.resolve())),
                    "-o", get_short_path_name(str(pdbqt_path.resolve())),
                    "-U", "nphs_lps",
                    "-A", "checkhydrogens",
                ]
                logging.info("Fallback with MGLTools on RDKit-cleaned PDB: " + " ".join(prepare_cmd2))

                # capture output to avoid console spam, but keep it for logs on failure
                try:
                    res2 = subprocess.run(
                        prepare_cmd2,
                        check=True,
                        capture_output=True,
                        text=True,
                        timeout=600
                    )
                    logging.info(res2.stdout)
                except subprocess.CalledProcessError as e:
                    logging.error(f"Fallback MGLTools failed for {sanitized.name}:\n{e.stderr}")
                except subprocess.TimeoutExpired:
                    logging.error(f"Fallback timeout for {sanitized.name}")

            finally:
                # Always try to clean up temp PDB
                try:
                    tmp_pdb.unlink(missing_ok=True)
                except Exception:
                    pass


        # Final validation & quarantine on failure
        try:
            sanitized.unlink(missing_ok=True)
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



# ---------- Main SDF→MOL2→PDBQT pipeline ----------

def prep_ligands_with_mgltools():
    print("Starting ligand preparation")

    cfg = read_config()  # read_config keeps original case; config uses UPPERCASE

    # --- Directories ---
    ligand_extracted_dir = Path(cfg["LIGAND_EXTRACTED_DIR"]).resolve()
    ligands_mol2_dir     = Path(cfg["LIGANDS_MOL2_DIR"]).resolve()
    output_ligands_dir   = Path(cfg["OUTPUT_LIGANDS_DIR"]).resolve()
    output_ligands_dir.mkdir(parents=True, exist_ok=True)

    # --- Tool paths ---
    mgltools_python = cfg["MGLTOOLS_PYTHON"]
    mgltools_path = cfg["MGLTOOLS_PATH"]
    obabel_cfg = cfg["OPENBABEL_PATH"]

    # Normalize to the actual .exe if a folder was provided
    obabel_exe = obabel_cfg if obabel_cfg.lower().endswith(".exe") else str(Path(obabel_cfg) / "obabel.exe")
    obabel_exe_short = get_short_path_name(obabel_exe)

    # Optional: set BABEL_DATADIR from Open Babel install if not set
    if not os.environ.get("BABEL_DATADIR"):
        obabel_dir = Path(obabel_exe).resolve().parent
        data_dir = obabel_dir / "data"
        if data_dir.exists():
            os.environ["BABEL_DATADIR"] = str(data_dir)

    # Sanity checks with clear messages
    for label, p in [
        ("MGLTOOLS_PYTHON", mgltools_python),
        ("MGLTOOLS_PATH",   mgltools_path),
        ("OPENBABEL_PATH",  obabel_exe),
        ("LIGAND_EXTRACTED_DIR", ligand_extracted_dir),
        ("LIGANDS_MOL2_DIR",     ligands_mol2_dir),
        ("OUTPUT_LIGANDS_DIR",   output_ligands_dir),
    ]:
        if not str(p).strip():
            raise RuntimeError(f"Config value missing/empty: {label}")

    mgltools_python_short = get_short_path_name(mgltools_python)
    prepare_script = Path(mgltools_path) / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
    if not prepare_script.exists():
        raise FileNotFoundError(f"prepare_ligand4.py not found at {prepare_script}")
    prepare_script_short = get_short_path_name(str(prepare_script.resolve()))
    obabel_exe_short = get_short_path_name(obabel_exe)

    sdf_files = list(ligand_extracted_dir.glob("*.sdf"))
    print(f"Found {len(sdf_files)} SDF file(s)")
    if not sdf_files:
        return

    for sdf_file in sdf_files:
        print(f"\n=== Processing SDF: {sdf_file.name} ===")
        sdf_abs = sdf_file.resolve()

        # Step 1: Convert SDF to per-molecule MOL2s with OBabel (chunked, visible progress)
        max_workers = max(1, min(8, os.cpu_count() or 8))

        if USE_RDKIT_FOR_3D:
            print("Using RDKit ETKDG for 3D, OBabel only for format conversion …")
            mol2_files = rdkit_embed_sdf_to_mol2(
                sdf_abs, ligands_mol2_dir, obabel_exe=obabel_exe_short, max_workers=max_workers
            )
        else:
            print("Using OBabel --gen3d (may be less stable on Windows) …")
            mol2_files = convert_sdf_to_mol2_split_parallel(
                sdf_abs, ligands_mol2_dir, obabel_exe_short,
                threads=OBABEL_THREADS, timeout_sec=OBABEL_TIMEOUT_S, chunk_size=CHUNK_SIZE
            )

        if not mol2_files:
            print("No MOL2 files produced; skipping this SDF.")
            continue

        # If OBabel gen3d path yields nothing, try RDKit ETKDG fallback
        if not mol2_files:
            print("OBabel --gen3d failed across chunks; trying RDKit ETKDG fallback …")
            mol2_files = rdkit_embed_sdf_to_mol2(
                sdf_abs, ligands_mol2_dir, obabel_exe=obabel_exe_short,
                max_workers=min(8, os.cpu_count() or 8)
            )

        if not mol2_files:
            print("No MOL2 files produced; skipping this SDF.")
            continue

        # Step 2: Prepare each MOL2 with MGLTools in parallel
        print(f"Preparing {len(mol2_files)} MOL2 files with MGLTools (parallel)…")
        max_workers = max(1, min(8, os.cpu_count() or 8))
        futures = []
        ok_count = 0
        fail_count = 0
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            for mol2_file in mol2_files:
                pdbqt_path = output_ligands_dir / f"{mol2_file.stem}.pdbqt"
                futures.append(ex.submit(
                    _prepare_one, mgltools_python_short, prepare_script_short, mol2_file, pdbqt_path
                ))

            total = len(futures)
            for i, fut in enumerate(as_completed(futures), 1):
                name, status = fut.result()
                if status == "ok":
                    ok_count += 1
                else:
                    fail_count += 1
                if i % 100 == 0 or status != "ok":
                    print(f"[{i}/{total}] {name}: {status}")

        print(f"Done: {ok_count} ok, {fail_count} failed/quarantined for {sdf_file.name}")


if __name__ == "__main__":
    prep_ligands_with_mgltools()
