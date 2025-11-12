
import os, re, json, subprocess, shutil, math, tempfile, hashlib
from pathlib import Path
from typing import Dict, Tuple, Optional, Set
import logging
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio import BiopythonWarning
from activesite import fix_pdb_elements

#this filters out the occupancy messages, not  really relevant to us. occupancy from my understanding
# is not used by reduce, vina or anything else here 
warnings.filterwarnings("ignore", category=BiopythonWarning, module="Bio.PDB.PDBIO")

try:
    # Prefer the newer loader/validator if present
    from input_and_export_functions import load_config as _load_cfg_new, validate_config as _validate_cfg
    _CFG = _load_cfg_new("config.txt")
    try:
        _validate_cfg(_CFG)
    except Exception:
        pass
except Exception:
    # Legacy fallback
    try:
        from installation import load_config as _load_cfg_legacy
        _CFG = _load_cfg_legacy()
    except Exception:
        _CFG = {}

def _cfg(key: str, default: str = "", legacy_key: str | None = None) -> str:
    v = os.environ.get(key)
    if v not in (None, ""):
        return v
    if key in _CFG and str(_CFG.get(key)) != "":
        return str(_CFG.get(key))
    if legacy_key and legacy_key in _CFG and str(_CFG.get(legacy_key)) != "":
        return str(_CFG.get(legacy_key))
    return default

def _pick_reduce_exe() -> str:
    # 1) explicit config/env
    explicit = (os.environ.get("REDUCE_EXE") or os.environ.get("REDUCE_BIN") or
                _cfg("REDUCE_EXE", "", "reduce_exe"))
    if explicit and Path(explicit).exists():
        return explicit
    # 2) repo-local build
    here = Path(__file__).resolve().parent
    # Also try the repository root's tools/ directory (../../tools/reduce/...)
    repo_root = here.parent.parent
    cand_root = repo_root / "tools" / "reduce" / "reduce_src" / "reduce"
    if cand_root.exists() and os.access(str(cand_root), os.X_OK):
        return str(cand_root)
    cand = here / "tools" / "reduce" / "reduce_src" / "reduce"
    if cand.exists() and os.access(str(cand), os.X_OK):
        return str(cand)
    # 3) PATH
    which = shutil.which("reduce") or shutil.which("reduce.exe")
    return which or "reduce"

def _het_dict_path() -> str | None:
    hd = os.environ.get("REDUCE_HET_DICT") or _cfg("REDUCE_HET_DICT", "", "reduce_het_dict")
    if hd and Path(hd).is_file():
        return hd
    # repo default if bundled
    here = Path(__file__).resolve().parent
    repo_root = here.parent.parent
    default_root = repo_root / "tools" / "reduce" / "reduce_wwPDB_het_dict.txt"
    if default_root.is_file():
        return str(default_root)
    default = here / "tools" / "reduce" / "reduce_wwPDB_het_dict.txt"
    return str(default) if default.is_file() else None

REDUCE_EXE = _pick_reduce_exe()
logging.info("propka_wire: Using Reduce at %s", REDUCE_EXE)

import hashlib
def _count_atoms_pdb(path: Path) -> int:
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            return sum(1 for ln in fh if ln.startswith(("ATOM  ", "HETATM")))
    except Exception:
        return 0

def _sha1_of_file(p: Path) -> str:
    h = hashlib.sha1()
    try:
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return ""

def _peek_lines(p: Path, n: int = 12) -> list[str]:
    try:
        return p.read_text(errors="ignore").splitlines()[:n]
    except Exception:
        return []

# --------------------
# Small exec helpers
# --------------------
def _which(name: str) -> Optional[str]:
    p = shutil.which(name)
    return p if p else None

def _has_exe(name: str) -> Optional[str]:
    return _which(name)

def _ensure_exec(env_var: str, fallback_names: list[str]) -> Optional[str]:
    cand = os.environ.get(env_var, "")
    if cand and Path(cand).exists():
        return cand
    for nm in fallback_names:
        p = _which(nm)
        if p:
            return p
    return None

# --------------------
# PDB2PQR / PROPKA
# --------------------
def _strip_pqr_to_pdb(pqr_path: Path, pdb_out: Path) -> None:
    """Convert PQR to PDB by dropping charge/radius columns, keeping coords & Hs."""
    lines = []
    for line in Path(pqr_path).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            # Keep columns 1..54 (xyz), plus element in cols 77..78 if present
            core = line[:54]
            element = line[76:78] if len(line) >= 78 else ""
            # Re-pad to a minimal valid PDB ATOM/HETATM line
            new = f"{core:54s}{'':6s}{'':6s}{element:>2s}\n"
            lines.append(new)
        else:
            lines.append(line + ("\n" if not line.endswith("\n") else ""))
    pdb_out.write_text("".join(lines))

def _copy_if_exists(src: Path, dst: Path) -> Optional[Path]:
    if src and src.exists():
        shutil.copy2(src, dst)
        return dst
    return None
def _pick_propka_exe() -> Optional[str]:
    """
    Choose PROPKA CLI with precedence:
      1) config.txt (key: PROPKA_EXE, via installation.load_config)
      2) env var PROPKA_EXE
      3) PATH: propka31 / propka30 / propka
    """
    # 1) config.txt
    cfg_path = Path(__file__).resolve().parent / "config.txt"
    try:
        from installation import load_config as _load_cfg
        cfg = _load_cfg() if cfg_path.exists() else {}
    except Exception:
        cfg = {}
    if isinstance(cfg, dict):
        cand = str(cfg.get("PROPKA_EXE", "")).strip()
        if cand and Path(cand).exists():
            return cand

    # 2) env
    env_cand = os.environ.get("PROPKA_EXE", "").strip()
    if env_cand and Path(env_cand).exists():
        return env_cand

    # 3) PATH fallbacks
    return (_has_exe("propka31") or _has_exe("propka30") or
            _has_exe("propka3")  or _has_exe("propka"))




def pdb2pqr_protonate(
    pdb_in: str,
    target_ph: float,
    out_dir: str | Path,
    ff: str = "amber",
    keep_waters: bool = True
) -> Tuple[Optional[str], Optional[str]]:
    """
    Run PDB2PQR (with PROPKA) at the requested pH and return (protonated_pdb, propka_log).
    - Writes a temporary .pqr then converts back to .pdb (coords + Hs, charges dropped).
    - If pdb2pqr not found, returns (None, None).
    """
    pdb2pqr = _ensure_exec("PDB2PQR_EXE", ["pdb2pqr"])
    if not pdb2pqr:
        logger.warning("pdb2pqr not found on PATH; skipping pre-protonation.")
        return None, None

    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    base = Path(pdb_in).stem
    tag = str(target_ph).replace(".", "_")
    pqr_out = out_dir / f"{base}.p{tag}.pqr"
    pdb_out = out_dir / f"{base}.p{tag}.pdb"
    pk_log  = out_dir / f"{base}.propka_pka.txt"
    # PDB2PQR 3.x expects uppercase FF names; default to AMBER if unknown
    ALLOWED_FF = {"AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON"}
    ff_norm = (ff or "AMBER").upper()
    if ff_norm not in ALLOWED_FF:
        ff_norm = "AMBER"
    cmd = [
        pdb2pqr,
        f"--ff={ff_norm}",
        f"--with-ph={target_ph:.2f}",
        "--keep-chain",
        "--titration-state-method=propka",
        # note: default is to keep waters; only add --drop-water if requested
    ]
    variant = os.environ.get("VARIANT", None)
    if variant == "HOLO":
        keep_hetero = True
    else:
        keep_hetero = False
    if keep_hetero:
        cmd.append("--keep-hetero")
    if not keep_waters:
        cmd.append("--drop-water")
    cmd.extend([str(pdb_in), str(pqr_out)])
    logger.info(f"[pdb2pqr] keep_hetero={keep_hetero} in={str(pdb_in)} out={str(pdb_out)}")

    try:
        res = subprocess.run(cmd, check=True, text=True, capture_output=True, cwd=str(out_dir))
        _strip_pqr_to_pdb(pqr_out, pdb_out)
        fix_pdb_elements(str(pdb_out))
        logger.info(f"[pdb2pqr.fix] applied elemfix to restore element fields in {str(pdb_out)}")
        # Copy PROPKA table if it was emitted near the PQR (cwd was set to out_dir)
        pka_candidate = next((p for p in Path(out_dir).glob("*.propka*")), None)
        _copy_if_exists(pka_candidate, pk_log)
        logger.info("[pdb2pqr] ph=%.2f out=%s pkas=%s", float(target_ph), str(pdb_out), str(pk_log.exists()))
        return str(pdb_out), (str(pk_log) if pk_log.exists() else None)

    except subprocess.CalledProcessError as e:
        # Write full stdout/stderr so we can see why pdb2pqr failed
        fail_txt = out_dir / f"{base}.pdb2pqr.failed.txt"
        try:
            with open(fail_txt, "w", encoding="utf-8") as fh:
                if e.stdout:
                    fh.write(e.stdout)
                    fh.write("\n")
                fh.write("--- STDERR ---\n")
                if e.stderr:
                    fh.write(e.stderr)
        except Exception:
            pass
        logger.error("pdb2pqr failed at pH %.2f rc=%s; see %s", target_ph, str(e.returncode), str(fail_txt))
        return None, None
    except Exception as e:
        logger.error("pdb2pqr failed at pH %.2f: %s", target_ph, e)
        return None, None

# --------------------
# Local titration by PROPKA states
# --------------------
def _collect_residues_within(pdb_in: str, center_xyz: Tuple[float, float, float], radius_A: float) -> Set[Tuple[str,int,str]]:
    """Collect (chain, resi, resname) for residues with any atom within radius_A."""
    cx, cy, cz = center_xyz
    out: Set[Tuple[str,int,str]] = set()
    for line in Path(pdb_in).read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            resn = line[17:20].strip().upper()
            chain = line[21].strip() or "?"
            resi = int(line[22:26])
        except Exception:
            continue
        if (x-cx)**2 + (y-cy)**2 + (z-cz)**2 <= radius_A**2:
            out.add((chain, resi, resn))
    return out
# --- BEGIN: PROPKA-only prestate generator (no hydrogens here) ---

import os, re, shutil, subprocess, logging
from pathlib import Path
from typing import Dict, Tuple, Optional

_log = logging.getLogger("propka_wire")
if not _log.handlers:
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("[propka_wire] %(message)s"))
    _log.addHandler(_h)
_log.setLevel(logging.INFO)
logger = _log #lazy little  shim 
# Residue rename rules
ACIDS   = {"ASP", "GLU", "CYS", "TYR"}
BASES   = {"HIS", "LYS", "ARG"}
RENAMES = {
    ("ASP", True): "ASH", ("ASP", False): "ASP",
    ("GLU", True): "GLH", ("GLU", False): "GLU",
    ("CYS", False): "CYM",  # deprotonated sulfur; protonated stays CYS
    ("LYS", False): "LYN",  # deprotonated lysine; protonated stays LYS
    # HIS needs tautomer; default neutral = HIE; protonated = HIP; forced alternative neutral = HID (rarely)
}

# Accept blank/lowercase chain and optional insertion code next to resseq.
_PKA_RE = re.compile(
    r"^\s*(\d+)\s+([A-Za-z]{3})\s+(\S?)\s+(\d+)([A-Za-z]?)\s+([-\d\.]+)"
)
# A more permissive matcher for newer PROPKA tables (3.5.x can vary spacing/columns)
_PKA_RE2 = re.compile(
    r"^\s*\d+\s+([A-Za-z]{3})\s+(\S?)\s+(\d+)[A-Za-z]?\s+([-+]?\d+(?:\.\d+)?)\s*$"
)
def _run_propka(pdb_in: Path, ph: float, out_dir: Path) -> Optional[Path]:
    """
    Run propka with an explicit pH and return a per-pH .pka:
      <basename>.pH{ph_tag}.pka  (e.g., 6LYZ_cleaned.pH7_00.pka)

    Note: PROPKA 3.5.x sometimes omits an explicit “pH=X” banner in the .pka header.
    We therefore no longer reject based on header content; we accept any non-empty .pka.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    ph_tag = f"{ph:.2f}".replace(".", "_")
    pka_file = out_dir / f"{pdb_in.stem}.pH{ph_tag}.pka"
    dbg_file = out_dir / (pdb_in.stem + f".pH{ph_tag}.propka.debug.txt")

    exe = _pick_propka_exe()
    if not exe:
        _log.warning("[propka.exec] no propka exe found; skipping")
        return None

    cmd = [exe, str(pdb_in), "--pH", f"{ph:.2f}"]

    tried = []
    try:
        # Run in out_dir so propka emits its default-named outputs there
        cp = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True)
        tried.append((cmd, cp.returncode))

        # PROPKA typically writes <basename>.pka — capture that, then copy to our per-pH name
        default_pka = out_dir / f"{pdb_in.stem}.pka"
        if default_pka.exists() and default_pka.stat().st_size > 0:
            shutil.copy2(default_pka, pka_file)

        # Write a detailed debug trace (first 200 lines for format triage)
        with open(dbg_file, "w", encoding="utf-8") as fh:
            fh.write(f"[propka.exec] cwd={out_dir}\n")
            for tcmd, rc in tried:
                fh.write(f"[propka.try] rc={rc} cmd={' '.join(tcmd)}\n")
            cand = pka_file if pka_file.exists() else (out_dir / f"{pdb_in.stem}.pka")
            if cand.exists():
                fh.write(f"[propka.out] exists=1 size={cand.stat().st_size} sha1={_sha1_of_file(cand)}\n")
                try:
                    lines = cand.read_text(errors="ignore").splitlines()
                    for ln in lines[:200]:
                        fh.write(f"[propka.pka.peek] {ln}\n")
                except Exception as e:
                    fh.write(f"[propka.peek.error] {e}\n")
            else:
                fh.write("[propka.out] exists=0\n")

    except FileNotFoundError:
        _log.warning("[propka.exec] exe missing; tried=%s", exe)
    except Exception as e:
        _log.error("[propka.exec] ph=%.2f failed: %s", ph, e)

    if pka_file.exists() and pka_file.stat().st_size > 0:
        _log.info("[propka.exec] ph=%.2f wrote_pka=1 sha1=%s size=%d cwd=%s",
                  ph, _sha1_of_file(pka_file), pka_file.stat().st_size, str(out_dir))
        # Consistency: filename pH tag vs runtime pH
        ph_from_name = parse_ph_from_name(pka_file.name)
        if ph_from_name is not None and abs(ph_from_name - float(ph)) > 1e-2:
            _log.warning("[propka.exec] pH tag mismatch: name=%.2f runtime=%.2f file=%s",
                         ph_from_name, float(ph), pka_file.name)
        # Quick fingerprint (rows via tolerant parser)
        try:
            rows = len(_parse_pka_table(pka_file))
            _log.info("[propka.fingerprint] ph=%.2f file=%s size=%d sha1=%s rows=%d",
                      float(ph), pka_file.name, pka_file.stat().st_size, _sha1_of_file(pka_file), rows)
        except Exception:
            pass
        return pka_file

    _log.warning("[propka.exec] ph=%.2f wrote_pka=0 cwd=%s", ph, str(out_dir))
    return None

def parse_ph_from_name(name: str) -> Optional[float]:
    """
    Extract pH from filenames like ..._pH8_5.pka, ..._pH7.00.pka, ..._pH9.pka (case-insensitive).
    Returns float pH if found, else None.
    """
    s = (name or "").strip()
    m = re.search(r'(?i)(?:^|[._-])pH\s*([0-9]+)(?:[_\. ]([0-9]+))?', s)
    if not m:
        return None
    major = m.group(1)
    minor = m.group(2) or ""
    try:
        if minor != "":
            return float(f"{major}.{minor}")
        return float(major)
    except Exception:
        return None


def _parse_pka_table(pka_path: Path) -> Dict[Tuple[str,int,str], float]:
    """
    Robust PROPKA parser for v3.1–3.5.1 tables.

    Returns a dict keyed by (RESN, RESID, CHAIN) -> pKa.
    Also stores chain-agnostic fallback keys (CHAIN="").
    """
    table: Dict[Tuple[str,int,str], float] = {}
    if not pka_path or not pka_path.exists():
        _log.info("[propka.pka.stats] parsed_rows=0 unique_keys=0 (no file)")
        return table

    raw = pka_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    # --- Try legacy strictly-formatted pattern first (for older propka dumps) ---
    legacy_rows = 0
    for line in raw:
        m = _PKA_RE.match(line)
        if not m:
            continue
        _idx, resn, chain, resi, icode, pka = m.groups()
        try:
            resn_u  = resn.upper()
            if not re.fullmatch(r"[A-Z]{3}", resn_u):
                continue  # skip termini like C-, N+
            chain_u = (chain or "").upper() or "_"
            resi_i  = int(resi)   # insertion code is ignored
            val     = float(pka)  # tolerate sentinel 99.99 etc.
        except Exception:
            continue
        table[(resn_u, resi_i, chain_u)] = val
        table[(resn_u, resi_i, "")]      = val
        legacy_rows += 1

    if legacy_rows > 0:
        _log.info("[propka.pka.stats] parsed_rows=%d unique_keys=%d (legacy pattern)",
                  legacy_rows, len(table))
        items = list(table.items())[:10]
        sample = "; ".join([f"{k[2]}:{k[1]}:{k[0]}={v:.2f}" for k, v in items])
        _log.info("[propka.pka.sample] %s", sample)
        return table

    # --- SUMMARY block parser (tolerant across 3.x variants) ---
    in_summary = False
    parsed_rows = 0
    for line in raw:
        U = line.upper()
        if not in_summary:
            if "SUMMARY OF THIS PREDICTION" in U:
                in_summary = True
            continue

        # Stop at the first major separator/footer after we've parsed rows
        if re.search(r"-{5,}", line) and parsed_rows > 0:
            break
        if ("FREE ENERGY" in U or "PROTEIN CHARGE" in U or "REFERENCES" in U) and parsed_rows > 0:
            break

        s = line.strip()
        if not s or set(s) <= {"-"}:
            continue
        if s.upper().startswith(("GROUP", "RESIDUE")):
            continue

        toks = s.split()
        if len(toks) < 3:
            continue

        # Token 0 must be a 3-letter residue name (skip terminal groups like C-, N+)
        resn_tok = toks[0].upper()
        if not re.fullmatch(r"[A-Z]{3}", resn_tok):
            continue

        # Token 1 is residue index, possibly with insertion code (e.g., 35A)
        m_resi = re.match(r"(\d+)", toks[1])
        if not m_resi:
            continue
        resi_i = int(m_resi.group(1))

        # Token 2 is either a chain ID (single letter) or the first numeric pKa
        chain_u = "_"
        pka_val: Optional[float] = None
        idx_from = 2
        if len(toks[2]) == 1 and toks[2].isalpha():
            chain_u = toks[2].upper()
            idx_from = 3

        # Find first numeric token from idx_from onward as the pKa
        for tk in toks[idx_from:]:
            if re.match(r"^[-+]?\d+(?:\.\d+)?$", tk):
                try:
                    pka_val = float(tk)
                except Exception:
                    pka_val = None
                break

        if pka_val is None:
            continue

        key1 = (resn_tok, resi_i, chain_u)
        key2 = (resn_tok, resi_i, "")
        table[key1] = pka_val
        table[key2] = pka_val
        parsed_rows += 1

    _log.info("[propka.pka.stats] parsed_rows=%d unique_keys=%d", parsed_rows, len(table))
    if parsed_rows:
        items = list(table.items())[:10]
        sample = "; ".join([f"{k[2]}:{k[1]}:{k[0]}={v:.2f}" for k, v in items])
        _log.info("[propka.pka.sample] %s", sample)
    else:
        _log.warning("[propka.pka.sample] no rows parsed; verify .pka format")

    return table



def _within_sphere(x: float, y: float, z: float, cx: float, cy: float, cz: float, r2: float) -> bool:
    dx, dy, dz = x - cx, y - cy, z - cz
    return (dx*dx + dy*dy + dz*dz) <= r2

def _choose_his_name(delta: float) -> str:
    """
    HIS policy with gentle hysteresis:
      delta = pKa - pH
      if delta >= +0.5 -> HIP (+1)
      else -> neutral HIE (dock-friendly default)
    """
    return "HIP" if delta >= 0.5 else "HIE"


def _protonated_is_true_for(resname: str, pka_minus_ph: float) -> Optional[bool]:
    """
    Gentle hysteresis to avoid flip-flopping:
      ACIDS (ASP,GLU):   neutral (ASH/GLH) if (pKa - pH) >= +0.5 -> protonated=True
                         else deprotonated (ASP/GLU) -> protonated=False
      CYS:               deprotonate to CYM only when (pH - pKa) >= 1.0
      TYR:               keep TYR; only consider deprot when (pH - pKa) >= 1.0 (we keep name 'TYR')
      BASES (LYS):       deprotonate to LYN only when (pH - pKa) >= 1.0
      ARG:               practically always protonated for docking; leave as-is (True)
    """
    r = resname.upper()
    d = pka_minus_ph  # = pKa - pH

    if r in {"ASP", "GLU"}:
        # protonated (neutralized name) when delta >= +0.5
        return True if d >= 0.5 else False

    if r == "CYS":
        # deprotonate only when pH - pKa >= 1.0  =>  delta <= -1.0
        return False if d <= -1.0 else True  # True==protonated (stay CYS), False==deprot (-> CYM)

    if r == "TYR":
        # keep protonated unless strong evidence to deprotonate
        return True if d > -1.0 else False

    if r == "LYS":
        # deprotonate (-> LYN) only when pH - pKa >= 1.0
        return False if d <= -1.0 else True

    if r == "ARG":
        return True  # keep protonated for docking

    return None  # HIS handled separately


def _rename_line(line: str, new3: str) -> str:
    # PDB residue name columns 18-20 (1-indexed), i.e., [17:20] 0-indexed
    return line[:17] + f"{new3:>3}" + line[20:]

def _extract_xyz(line: str) -> Tuple[float,float,float]:
    # X at cols 31-38, Y at 39-46, Z at 47-54 (1-indexed). Use robust slicing.
    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
    return x,y,z

def apply_propka_states(
    cleaned_receptor_pdb: str,
    center: Tuple[float,float,float],
    radius: float,
    ph: float,
    out_dir: str,
    tag: str
) -> Tuple[str, str, int]:
    """
    Produce PROPKA-informed, pocket-localized renames only (or whole-protein in GLOBAL mode).
    Inputs:
      cleaned_receptor_pdb: path to base receptor (no hydrogens required)
      center: (cx,cy,cz), sphere center in Å
      radius: sphere radius in Å; if radius >= 1e6 or math.isinf(radius) => GLOBAL mode
      ph: target environmental pH
      out_dir: directory for outputs
      tag: filename tag (e.g., '6LYZ_pH7_0')
    Outputs:
      (prestate_pdb_path, pka_path, n_renamed)
    Side effects:
      - Writes <tag>.pka (if propka ran)
      - Writes <tag>.prestate.pdb (never empty; guard copies input when needed)
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    pdb_in  = Path(cleaned_receptor_pdb)
    prestate_pdb = out_dir / f"{tag}.prestate.pdb"
    pka_out      = out_dir / f"{tag}.pka"

    # GLOBAL titration sentinel (preserves call sites; pass a huge radius to enable)
    global_mode = (radius >= 1e6) or math.isinf(radius)
    cx, cy, cz = center
    r2 = float("inf") if global_mode else (radius * radius)

    def in_scope(x: float, y: float, z: float) -> bool:
        # Whole protein in GLOBAL mode; otherwise pocket sphere
        if global_mode:
            return True
        dx, dy, dz = x - cx, y - cy, z - cz
        return (dx*dx + dy*dy + dz*dz) <= r2

    # 1) Run PROPKA and persist a copy of the pKa table as <tag>.pka
    pka_file = _run_propka(pdb_in, ph, out_dir)
    if pka_file and pka_file.exists():
        try:
            shutil.copy2(pka_file, pka_out)
        except Exception:
            if pka_file.resolve() != pka_out.resolve():
                shutil.copyfile(pka_file, pka_out)
    else:
        # No pKa table; proceed with rename heuristics = none (n_renamed may be 0)
        pka_out.touch(exist_ok=True)

    pka_table = _parse_pka_table(pka_out) if pka_out.exists() else {}
    # --- Diagnostics: fingerprint and quick peek of the PROPKA table for this pH ---
    if pka_out.exists() and pka_out.stat().st_size > 0:
        pka_sha1 = _sha1_of_file(pka_out)
        _log.info("[propka.fingerprint] ph=%.2f file=%s size=%d sha1=%s rows=%d",
                  ph, pka_out.name, pka_out.stat().st_size, pka_sha1, len(pka_table))
        for ln in _peek_lines(pka_out, 8):
            _log.info("[propka.table.peek] %s", ln)
    else:
        _log.info("[propka.fingerprint] ph=%.2f file=%s size=0 sha1= rows=%d (no table)",
                  ph, pka_out.name, len(pka_table))

    # --- Diagnostics: enumerate residues in-scope and rename choices (cap to avoid spam) ---
    _diag_shown = 0
    _diag_cap = 40  # don’t flood logs; raise if you want more
    try:
        with open(cleaned_receptor_pdb, "r", errors="ignore") as _fh:
            seen_keys = set()
            for _ln in _fh:
                if not _ln.startswith(("ATOM  ", "HETATM")):
                    continue
                try:
                    x = float(_ln[30:38]); y = float(_ln[38:46]); z = float(_ln[46:54])
                except Exception:
                    continue
                if not in_scope(x, y, z):
                    continue
                resn = _ln[17:20].strip().upper()
                chain = (_ln[21].strip() or " ")
                resi = int(_ln[22:26])
                k = (resn, resi, chain)
                if k in seen_keys:
                    continue
                seen_keys.add(k)

                pka = (pka_table.get((resn, resi, chain))
                       or pka_table.get((resn, resi, ""))
                       or pka_table.get((resn, resi, "?"))
                       or None)

                if resn == "HIS":
                    new3 = _choose_his_name((pka - ph) if pka is not None else -999.0) if pka is not None else None
                elif resn in ACIDS or resn in BASES:
                    if pka is not None:
                        prot = _protonated_is_true_for(resn, pka - ph)
                        new3 = RENAMES.get((resn, True)) if prot is True else RENAMES.get(
                            (resn, False)) if prot is False else None
                    else:
                        new3 = None
                else:
                    new3 = None

                if _diag_shown < _diag_cap:
                    if pka is None:
                        _log.info("[propka.choice] ph=%.2f %s%d:%s pKa=? delta=? -> %s",
                                  ph, chain, resi, resn, str(new3 or resn))
                    else:
                        _log.info("[propka.choice] ph=%.2f %s%d:%s pKa=%.2f delta=%.2f -> %s",
                                  ph, chain, resi, resn, pka, (pka - ph), str(new3 or resn))
                    _diag_shown += 1
    except Exception:
        pass

    # 2) Renaming pass (in-scope only)
    n_renamed = 0
    wrote_atoms = 0
    with open(pdb_in, "r", errors="ignore") as fh_in, open(prestate_pdb, "w") as fh_out:
        for line in fh_in:
            rec = line[:6]
            if rec == "ATOM  " or rec == "HETATM":
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except Exception:
                    fh_out.write(line)
                    continue

                if in_scope(x, y, z):
                    resn  = line[17:20].strip().upper()
                    chain = line[21].strip() or " "
                    resi  = int(line[22:26])
                    new3  = None

                    if resn in ACIDS or resn in BASES or resn == "HIS":
                        pka = (pka_table.get((resn, resi, chain))
                               or pka_table.get((resn, resi, ""))
                               or pka_table.get((resn, resi, "?"))
                               or None)
                        if resn == "HIS":
                            if pka is not None:
                                new3 = _choose_his_name(pka - ph)
                        else:
                            if pka is not None:
                                prot = _protonated_is_true_for(resn, pka - ph)
                                if prot is True:
                                    new3 = RENAMES.get((resn, True), resn)
                                elif prot is False:
                                    new3 = RENAMES.get((resn, False), resn)

                    if new3 and new3 != resn:
                        line = _rename_line(line, new3)
                        n_renamed += 1

                fh_out.write(line)
                wrote_atoms += 1
            else:
                fh_out.write(line)

    # 3) Guard: never leave an empty prestate
    if wrote_atoms == 0:
        shutil.copy2(pdb_in, prestate_pdb)
        _log.warning("[prestate.guard] wrote_atoms=0 \u2192 copied input\u2192output: %s", prestate_pdb)

    # 4) Compact rename summary (same logic; scope via in_scope)
    try:
        def _res_map(path: Path):
            m = {}
            with open(path, "r", errors="ignore") as f:
                for ln in f:
                    if ln.startswith(("ATOM  ", "HETATM")):
                        try:
                            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                        except Exception:
                            continue
                        if in_scope(x, y, z):
                            chain = (ln[21].strip() or " ")
                            resi = int(ln[22:26])
                            resn = ln[17:20].strip().upper()
                            m.setdefault((chain, resi), resn)
            return m

        before = _res_map(pdb_in)
        after  = _res_map(prestate_pdb)
        ren_list = []
        for k in sorted(set(before) | set(after)):
            b = before.get(k); a = after.get(k)
            if b and a and a != b:
                ren_list.append((k[0], k[1], b, a))
        if ren_list:
            _log.info("[propka.renames] ph=%.2f within_r=%s n=%d %s",
                      ph, ("ALL" if global_mode else f"{radius:.1f}"),
                      len(ren_list),
                      " ".join([f"{c}{i}:{b}->{a}" for (c, i, b, a) in ren_list]))
    except Exception:
        pass

    try:
        matched = sum(1 for (k, v) in pka_table.items() if k[2] != "")
        _log.info("[propka.match] pKa_keys=%d (with_chain) + %d (chainless)",
                  matched, sum(1 for (k, v) in pka_table.items() if k[2] == ""))
    except Exception:
        pass

    # 5) State summary + fingerprint (explicit scope marker)
    _log.info("[propka.states] scope=%s center=(%.3f,%.3f,%.3f) r=%s ph=%.2f n_renamed=%d",
              ("GLOBAL" if global_mode else "SPHERE"), cx, cy, cz,
              ("ALL" if global_mode else f"{radius:.2f}"), ph, n_renamed)
    try:
        pre_sha1 = _sha1_of_file(prestate_pdb)
        n_atoms = _count_atoms_pdb(prestate_pdb)
        _log.info("[prestate.sha1] ph=%.2f file=%s atoms=%d sha1=%s size=%d",
                  ph, prestate_pdb.name, n_atoms, pre_sha1, prestate_pdb.stat().st_size)
    except Exception:
        pass

    # 6) Fixed-order residue-name bins; in GLOBAL mode these are whole-protein counts
    try:
        counts = {
            "ASP": 0, "ASH": 0, "GLU": 0, "GLH": 0,
            "HIS": 0, "HID": 0, "HIE": 0, "HIP": 0,
            "LYS": 0, "LYN": 0, "CYS": 0, "CYM": 0, "TYR": 0
        }
        with open(prestate_pdb, "r", errors="ignore") as fh:
            seen = set()
            for ln in fh:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                except Exception:
                    continue
                if not in_scope(x, y, z):
                    continue
                chain = (ln[21].strip() or " ")
                resi = int(ln[22:26])
                resn = ln[17:20].strip().upper()
                key = (chain, resi)
                if key in seen:
                    continue
                seen.add(key)
                if resn in counts:
                    counts[resn] += 1

        _log.info(
            "[propka.states.bin] ph=%.2f ASP=%d ASH=%d GLU=%d GLH=%d HIS=%d HID=%d HIE=%d HIP=%d CYS=%d CYM=%d LYS=%d LYN=%d TYR=%d",
            ph,
            counts["ASP"], counts["ASH"], counts["GLU"], counts["GLH"],
            counts["HIS"], counts["HID"], counts["HIE"], counts["HIP"],
            counts["CYS"], counts["CYM"], counts["LYS"], counts["LYN"], counts["TYR"]
        )
    except Exception:
        pass

    return (str(prestate_pdb), str(pka_out), n_renamed)





# --------------------
# Reduce
# --------------------
def reduce_add_hydrogens(pdb_in: str, pdb_out: str, mode: str | None = None) -> None:
    """
    Add/optimize hydrogens with Reduce, with strict empty-output detection and fallbacks.
    mode=None         -> auto: if input has H, do flip-only; else full build
    mode="flip_only"  -> always use -FLIP -Quiet
    mode="full_build" -> always use -BUILD -Quiet
    """
    import subprocess, shutil, os as os
    # ---- Guard: input must exist before we try Reduce/OpenBabel ----
    pdb_in_path = Path(pdb_in)
    if not pdb_in_path.exists():
        logger.error("[reduce] input_missing in=%s; aborting hydrogenation", pdb_in)
        raise FileNotFoundError(pdb_in)
    def _count_atoms(path: str) -> int:
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                return sum(1 for ln in fh if ln.startswith(("ATOM  ", "HETATM")))
        except Exception:
            return 0

    def _has_h(path: str) -> bool:
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if ln.startswith(("ATOM  ", "HETATM")) and ln[12:16].strip().startswith("H"):
                        return True
        except Exception:
            pass
        return False

    # Decide flags
    if mode == "flip_only":
        flags = ["-FLIP", "-Quiet"]; assume_h = True
    elif mode == "full_build":
        flags = ["-BUILD", "-Quiet"]; assume_h = False
    else:
        assume_h = _has_h(pdb_in)
        flags = ["-FLIP", "-Quiet"] if assume_h else ["-BUILD", "-Quiet"]

    out_dir = Path(pdb_out).resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    het_dict = _het_dict_path()
    if het_dict:
        env["REDUCE_HET_DICT"] = het_dict

    exe = REDUCE_EXE or (shutil.which("reduce") or shutil.which("reduce.exe") or "reduce")

    def _run(stage: str, use_flags: list[str]) -> tuple[int, int, str]:
        """return (rc, wrote_atoms, stderr)"""
        cmd = [exe] + use_flags + [pdb_in]
        stderr_path = out_dir / f"{stage}.stderr.txt"
        with open(pdb_out, "w", encoding="utf-8") as out:
            cp = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, env=env)
        try:
            stderr_path.write_text(cp.stderr or "", encoding="utf-8")
        except Exception:
            pass
        wrote = _count_atoms(pdb_out)
        logger.info("[reduce] stage=%s rc=%s flags=%s wrote_atoms=%d", stage, cp.returncode, " ".join(use_flags), wrote)
        return cp.returncode, wrote, cp.stderr or ""

    # 1) Primary attempt
    rc, wrote, _ = _run("Reduce#1", flags)
    if wrote == 0:
        # 2) Alternate flag attempt (flip <-> build)
        alt = ["-BUILD","-Quiet"] if "-FLIP" in flags else ["-FLIP","-Quiet"]
        rc2, wrote2, _ = _run("Reduce#retry", alt)
        if wrote2 == 0:
            # 3) OpenBabel fallback
            ob = shutil.which("obabel") or shutil.which("obabel.exe")
            if ob:
                cmd = [ob, "-i", "pdb", pdb_in, "-o", "pdb", "-O", pdb_out, "-h"]
                cp3 = subprocess.run(cmd, text=True, capture_output=True)
                wrote3 = _count_atoms(pdb_out)
                logger.warning("[fallback] openbabel in=%s out=%s rc=%s wrote_atoms=%d", pdb_in, pdb_out, cp3.returncode, wrote3)
                (out_dir / "OpenBabel.stderr.txt").write_text(cp3.stderr or "", encoding="utf-8")
                if wrote3 == 0:
                    shutil.copy2(pdb_in, pdb_out)
                    logger.error("[reduce] all attempts failed; copied input→output (wrote_atoms=0)")
            else:
                shutil.copy2(pdb_in, pdb_out)
                logger.error("[reduce] no OpenBabel found; copied input→output (wrote_atoms=0)")
    # Success/acceptance (even if rc==1) is purely “has atoms”
    return




