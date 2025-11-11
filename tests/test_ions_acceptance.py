# tests/test_ions_acceptance.py

import os
import re
import sys
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import urllib.request
import socket

import pytest

# --- Repo path setup (adjust if needed) ---------------------------------------
ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = ROOT / "code" / "protein_automation"
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

# Try to import the module; tests degrade gracefully if only CLI exists.
try:
    import automate_protein_prep  # type: ignore
except Exception:
    automate_protein_prep = None  # will use CLI fallback


# --- Minimal PDB metal parser (pure stdlib) -----------------------------------

# Conservative list of metal elements commonly observed in PDBs
METAL_ELEMENTS = {
    "LI","NA","MG","AL","K","CA","V","CR","MN","FE","CO","NI","CU","ZN","GA",
    "RB","SR","Y","ZR","NB","MO","RU","RH","PD","AG","CD","IN","SN","SB",
    "CS","BA","LA","CE","PR","ND","SM","EU","GD","TB","DY","HO","ER","TM",
    "YB","LU","HF","TA","W","RE","OS","IR","PT","AU","HG","TL","PB","BI"
}

_METAL_RE = re.compile(r"^(?:HETATM|ATOM)\s")  # some structures mislabel ions as ATOM

CACHE_DIR = ROOT / "tests" / "data"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
CACHE_PDB = CACHE_DIR / "1BN1.pdb"
RCSB_URL = "https://files.rcsb.org/download/1BN1.pdb"

def _try_download_1bn1(dest: Path, timeout: float = 20.0) -> bool:
    try:
        with urllib.request.urlopen(RCSB_URL, timeout=timeout) as r:
            data = r.read()
        if not data or len(data) < 5000:  # sanity check
            return False
        dest.write_bytes(data)
        return True
    except (OSError, socket.timeout):
        return False
    except Exception:
        return False

def _element_from_pdb_line(line: str) -> Optional[str]:
    """
    Extract element:
      (1) PDB fixed-width cols 77-78 (1-based),
      (2) fallback: first alphabetic chars of atom name (cols 13-16).
    """
    if len(line) >= 78:
        elem = line[76:78].strip().upper()
        if elem:
            return elem
    if len(line) >= 16:
        name = line[12:16].strip().upper()
        # Take 1-2 leading letters (handles 'ZN', 'MG', 'FE', 'CL', 'NA', 'K', etc.)
        m = re.match(r"[A-Z]{1,2}", name)
        if m:
            return m.group(0)
    return None

def _resname_from_pdb_line(line: str) -> Optional[str]:
    if len(line) >= 20:
        return line[17:20].strip().upper()
    return None

def _res_id_from_pdb_line(line: str) -> str:
    chain = line[21].strip() if len(line) > 21 else ""
    resseq = line[22:26].strip() if len(line) >= 26 else ""
    icode = line[26].strip() if len(line) >= 27 else ""
    return f"{chain or '-'}:{resseq or '?'}{icode or ''}"

def parse_pdb_metals(pdb_path: Path) -> Tuple[int, Dict[str, int], List[Tuple[str, str, str]]]:
    """
    Returns:
      total_count,
      counts_by_resname_or_elem (keyed by RESNAME if present else ELEMENT),
      list of (resname, element, resid) for each detected metal ion.
    """
    counts: Dict[str, int] = {}
    items: List[Tuple[str, str, str]] = []
    if not pdb_path.exists():
        return 0, counts, items

    with pdb_path.open("r", errors="ignore") as fh:
        for line in fh:
            if not _METAL_RE.match(line):
                continue
            elem = _element_from_pdb_line(line) or ""
            resn = _resname_from_pdb_line(line) or ""
            elem = elem.upper()
            resn = resn.upper()
            # Normalize: drop trailing digits in resname (e.g., ZN2 -> ZN)
            resn_norm = re.sub(r"\d+$", "", resn)
            is_metal = (elem in METAL_ELEMENTS) or (resn_norm in METAL_ELEMENTS)
            if not is_metal:
                continue
            key = resn_norm or elem or "UNK"
            counts[key] = counts.get(key, 0) + 1
            items.append((resn_norm or resn, elem or "?", _res_id_from_pdb_line(line)))

    total = sum(counts.values())
    return total, counts, items


# --- Test runner hook (Codex can align this to your actual entrypoint) --------

def _run_prep_python_entrypoint(pdb_path: Path, outdir: Path, variant: str) -> None:
    """
    Preferred path: call a Python function if available. Codex can adapt this
    to whatever entrypoint you expose (e.g., prepare_receptor or run()).
    """
    if automate_protein_prep is None:
        raise RuntimeError("automate_protein_prep import failed")

    # Try a few conventional entrypoints; Codex can collapse to your canonical one.
    if hasattr(automate_protein_prep, "prepare_receptor"):
        # COMMON SIGNATURES (Codex: align to your code)
        # prepare_receptor(pdb_path: str, outdir: str, variant: str)
        try:
            automate_protein_prep.prepare_receptor(str(pdb_path), str(outdir), variant)  # type: ignore[attr-defined]
            return
        except TypeError:
            pass
    if hasattr(automate_protein_prep, "run"):
        automate_protein_prep.run(str(pdb_path), str(outdir), variant)  # type: ignore[attr-defined]
        return
    if hasattr(automate_protein_prep, "main"):
        # Some mains accept argv; others read env (APO_HOLO_VARIANT)
        os.environ["APO_HOLO_VARIANT"] = variant
        try:
            automate_protein_prep.main([str(pdb_path), "--out", str(outdir)])  # type: ignore[attr-defined]
            return
        except TypeError:
            automate_protein_prep.main()  # type: ignore[attr-defined]
            return

    raise RuntimeError("No suitable Python entrypoint found; use CLI fallback.")


def _run_prep_cli(pdb_path: Path, outdir: Path, variant: str) -> None:
    """
    Fallback: invoke your CLI. Codex: adjust this command to match your CLI.
    Examples you've used:
      - `python main.py --1bn1 --holo --fast` (uses embedded target shortcuts)
      - or flags that accept an explicit path + variant.
    """
    env = os.environ.copy()
    env["APO_HOLO_VARIANT"] = variant

    # Codex: EDIT the command below to your real CLI.
    # Option A (module): python -m protein_automation.main --pdb <file> --variant <HOLO|APO> --fast
    cmd = [
        sys.executable,
        str(CODE_DIR / "main.py"),
        "--pdb", str(pdb_path),
        "--variant", variant.lower(),
        "--fast",
    ]

    # If your main.py doesn't support these flags, Codex should map to your actual flags.
    subprocess.run(cmd, check=True, env=env, cwd=str(ROOT))


def run_prep_variant(pdb_path: Path, variant: str) -> Path:
    """
    Runs the prep for one variant into a temporary directory and returns the
    path to the cleaned receptor PDB (glob for '*_cleaned.pdb').
    """
    tmp = Path(tempfile.mkdtemp(prefix=f"ions_{variant.lower()}_"))
    outdir = tmp / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        try:
            _run_prep_python_entrypoint(pdb_path, outdir, variant)
        except Exception:
            _run_prep_cli(pdb_path, outdir, variant)

        # Find a cleaned receptor PDB
        cleaned = list(outdir.rglob("*_cleaned.pdb"))
        if not cleaned:
            # Some pipelines place it under processed_pdbs/<PDB>/receptor/
            cleaned = list(tmp.rglob("*_cleaned.pdb"))
        assert cleaned, f"No cleaned receptor PDB found in {outdir} (or {tmp})"
        return cleaned[0]
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise


# --- Fixtures & tests ---------------------------------------------------------

def _get_test_pdb() -> Path:
    env_path = os.environ.get("ATLAS_TEST_PDB", "").strip()
    if env_path:
        return Path(env_path).resolve()
    if CACHE_PDB.exists():
        return CACHE_PDB.resolve()
    if _try_download_1bn1(CACHE_PDB):
        return CACHE_PDB.resolve()
    pytest.skip(
        f"ATLAS_TEST_PDB not set and could not download {RCSB_URL}. "
        "Set ATLAS_TEST_PDB or vendor tests/data/1BN1.pdb."
    )


class TestIonRetention:
    @pytest.mark.ions_acceptance
    def test_input_has_metals_or_skip(self):
        pdb_path = _get_test_pdb()
        total_in, _, _ = parse_pdb_metals(pdb_path)
        if total_in == 0:
            pytest.skip(f"Input PDB has no metals: {pdb_path}")

    @pytest.mark.ions_acceptance
    def test_holo_retains_metals_vs_apo_reduces(self):
        pdb_path = _get_test_pdb()
        total_in, _, _ = parse_pdb_metals(pdb_path)
        if total_in == 0:
            pytest.skip(f"Input PDB has no metals: {pdb_path}")
        holo_clean = run_prep_variant(pdb_path, "HOLO")
        holo_total, *_ = parse_pdb_metals(holo_clean)
        assert holo_total >= 1, f"HOLO should retain metals; got 0 in {holo_clean}"
        apo_clean = run_prep_variant(pdb_path, "APO")
        holo_total, *_ = parse_pdb_metals(holo_clean)
        apo_total, *_ = parse_pdb_metals(apo_clean)
        assert apo_total < holo_total, (
            f"APO should strip metals vs HOLO (apo={apo_total}, holo={holo_total})"
        )
