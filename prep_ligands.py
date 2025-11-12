import sys
import ctypes
import os
from ctypes import wintypes, create_unicode_buffer
from pathlib import Path
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Tuple, Dict, Set, Any, Union
import re
from datetime import datetime
import shutil
from collections import defaultdict
import argparse
from input_and_export_functions import load_config, validate_config
# >>> PATHS IMPORT START
from path_router import make_paths
# >>> PATHS IMPORT END

from activesite import (
    fix_pdb_elements,
    fix_element_columns_in_file,  # unify element rewrite for columns 77–78
    _fix_ligand_element_columns_in_memory,
    derive_element,
    scan_helium_counts,
    rules_version,
    assert_no_helium_in_hydrogen_names,
    assert_no_helium_in_pdbqt,
    get_atom_rules,
)
# --- RDKit / Standardization imports ---
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover

try:
    from rdkit.Chem.MolStandardize import rdMolStandardize as _std  # unified handle

    _HAS_STD = True
except Exception:
    _std = None
    _HAS_STD = False
    
    
LIGPREP_PH = 7.4
KEEP_NONPOLAR_H = 1
OBABEL_TIMEOUT_S = 900
STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "SEC", "PYL", "MSE"
}

# Exclude common crystallization additives/buffers/metals (+ porphyrins / modified residues)
try:
    # Prefer chemdb if available
    from chemdb.chem_alias_db import EXCLUDE_HET_IDS as EXCLUDE_CRYSTAL_ADDITIVES
except Exception:
    EXCLUDE_CRYSTAL_ADDITIVES = {
        "HOH","CIT","TAR","SO4","PO4","CA","NA","K","MG","MN","ZN","GOL","EDO","PEG","MPD","TRS",
        "MES","HEPES","ACET","ACT","FMT","MAL","DMS","IPA","CLU","NAG","BOG","TOS","BES","PTS","OTF",
        "TRF","TFA","BF4","PF6","CL","BR","I","HEM","HEC","HEA","HEB","HEO","HEG","HEF","HEH","PTR",
        "TPO","SEP"
    }


# --- Salvage / logging config ---
RUN_TAG = datetime.now().strftime("%Y%m%d_%H%M%S")
MALFORMED_LOG = Path("malformed_ligands.txt")
QUARANTINE_DIRNAME = "quarantine"

# --- Resume mode toggle ---
RESUME_SKIP = False


MAX_HEAVY_ATOMS = 1200
MIN_ATOMS_FOR_DOCKING = 5
MIN_PARENT_HEAVY = 8

# ----------- tuning switches -----------
USE_RDKIT_FOR_3D = True
OBABEL_THREADS = 50
OBABEL_TIMEOUT_S = 900
CHUNK_SIZE = 200


# --------------------------------------

# --- salt remover (fallback path also uses this) ---
_REM = SaltRemover()
# =========================
# Alias/rules centralization
# =========================
from typing import Optional

# Effective, alias-derived sets/maps populated at import:
RULES = None
ALLOWED_ELEMENTS: Set[str] = set()
MONOATOMIC_IONS: Set[str] = set()
# AD4_TYPES already exists later in the file; we will reconcile it here.

def _alias_token_to_element_symbol(tok: str) -> Optional[str]:
    """
    Convert an alias token from YAML (usually uppercase like 'CL','NA','ZN')
    into an RDKit-style element symbol ('Cl','Na','Zn'). Falls back to
    first letter uppercase for non-standard tokens.
    """
    if not tok:
        return None
    u = tok.strip().upper()
    if len(u) == 1:
        return u  # C, N, O, F, H, B, P, I
    # two-letter common case
    if len(u) == 2 and u.isalpha():
        return u[0] + u[1].lower()
    # Handle a few known multi-letter symbols that appear as 2 chars in YAML already
    return u[0] + (u[1:].lower() if len(u) > 1 else "")

def allowed_elements_from_aliases(rules) -> Set[str]:
    """
    Build RDKit-style element symbols from YAML-driven one- and two-letter sets.
    """
    out: Set[str] = set()
    for e in (rules.one_letter or []):
        if e: out.add(e.strip().upper())
    for e2 in (rules.two_letter or []):
        sym = _alias_token_to_element_symbol(e2)
        if sym: out.add(sym)
    return out

def free_ion_elements_from_aliases(rules) -> Set[str]:
    """
    Translate meeko.drop_free_ions and halide resnames into element symbols for
    'one-atom free ion' filtering during ligand prep.
    """
    out: Set[str] = set()
    for tok in (getattr(rules, "meeko_drop_free_ions", []) or []):
        sym = _alias_token_to_element_symbol(tok)
        if sym: out.add(sym)
    for hal in (getattr(rules, "halide_resnames", []) or []):
        sym = _alias_token_to_element_symbol(hal)
        if sym: out.add(sym)
    return out

def _init_alias_rules_cache() -> None:
    """
    Resolve RULES from activesite.get_atom_rules() and derive:
      - ALLOWED_ELEMENTS (RDKit symbols, canonical)
      - MONOATOMIC_IONS (elements considered free ions for 1-atom ligands)
      - AD4_TYPES (prefer YAML-driven; fallback to existing set)
    Also emit quiet DEBUG logs when LIGPREP_DEBUG=1.
    """
    global RULES, ALLOWED_ELEMENTS, MONOATOMIC_IONS, AD4_TYPES
    try:
        RULES = get_atom_rules()
    except Exception as e:
        logging.warning("[ligprep] get_atom_rules() failed; using built-ins only (%s)", e)
        RULES = None

    if RULES is not None:
        try:
            # Centralized allow-list
            ALLOWED_ELEMENTS = allowed_elements_from_aliases(RULES)
            # Derive free-ion set for monatomic ligand guard
            derived_free_ions = free_ion_elements_from_aliases(RULES)
            # Keep behavior stable: union with historical set if present
            try:
                _legacy = MONOATOMIC_IONS if isinstance(MONOATOMIC_IONS, set) else set()
            except NameError:
                _legacy = set()
            MONOATOMIC_IONS = (derived_free_ions or set()) | _legacy

            # Prefer canonical AD4 types from YAML if available
            if hasattr(RULES, "ad4_types") and RULES.ad4_types:
                AD4_TYPES = set(RULES.ad4_types)

            if os.environ.get("LIGPREP_DEBUG", "") == "1":
                logging.debug("[ligprep] aliases version=%s", rules_version())
                logging.debug("[ligprep] ALLOWED_ELEMENTS=%s", ", ".join(sorted(ALLOWED_ELEMENTS)))
                if AD4_TYPES:
                    logging.debug("[ligprep] AD4_TYPES(sample)=%s", ", ".join(sorted(list(AD4_TYPES))[:12]))
                if MONOATOMIC_IONS:
                    logging.debug("[ligprep] MONOATOMIC_IONS=%s", ", ".join(sorted(MONOATOMIC_IONS)))
        except Exception as e:
            logging.warning("[ligprep] alias init error; falling back to local lists (%s)", e)

# Initialize once at import
_init_alias_rules_cache()

# =========================
# Utility / logging helpers
# =========================

def _looks_like_monoatomic_ion_pdbqt(lines: List[str]) -> bool:
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    if len(atom_lines) != 1:
        return False
    parts = atom_lines[0].split()
    elem = _element_from_adt(parts[-1]) if parts else ""
    return elem in MONOATOMIC_IONS


# =========================
# Test-mode / ONLY-set helpers
# =========================
_ONLY_TOKEN_RE = re.compile(r"[0-9]{1,7}")

def _normalize_only_token(tok: str) -> Optional[str]:
    """
    Accepts things like 'rdk_0004931', 'rdk_4931', '0004931', '4931'.
    Returns canonical 'rdk_0004931' or None if it can't be parsed.
    """
    if not tok:
        return None
    t = tok.strip().lower().replace(",", " ")
    if not t:
        return None
    # pull first 1–7 digit run
    m = _ONLY_TOKEN_RE.search(t)
    if not m:
        return None
    n = m.group(0)
    try:
        i = int(n)
    except Exception:
        return None
    if i < 0 or i > 9999999:
        return None
    return f"rdk_{i:07d}"

def _collect_only_from_env_and_cli(cli_only: Optional[List[str]] = None) -> Set[str]:
    """
    Merge LIGPREP_ONLY (env), LIGPREP_ONLY_FILE (env path), and --only (CLI list).
    Normalize all tokens to 'rdk_0000000'. Empty/invalid tokens are ignored.
    """
    out: Set[str] = set()

    # env var, allow comma/space separated
    env_only = os.environ.get("LIGPREP_ONLY", "")
    if env_only:
        for raw in re.split(r"[,\s]+", env_only.strip()):
            norm = _normalize_only_token(raw)
            if norm:
                out.add(norm)

    # env file (one per line; also allow commas/spaces within lines)
    env_file = os.environ.get("LIGPREP_ONLY_FILE", "")
    if env_file:
        try:
            with open(env_file, "r", encoding="utf-8") as fh:
                for line in fh:
                    for raw in re.split(r"[,\s]+", line.strip()):
                        norm = _normalize_only_token(raw)
                        if norm:
                            out.add(norm)
        except Exception as e:
            logging.warning("[test-mode] could not read LIGPREP_ONLY_FILE=%s: %s", env_file, e)

    # CLI list (repeatable/variadic)
    if cli_only:
        for raw in cli_only:
            # each item may itself contain commas/spaces
            for tok in re.split(r"[,\s]+", raw.strip()):
                norm = _normalize_only_token(tok)
                if norm:
                    out.add(norm)

    return out

MALFORMED_DIR: Path | None = None
def _log_malformed(p: Path, reason: str, log_dir: Path | None = None) -> None:
    """
    Append a one-line reason for a malformed ligand into the per-protein log
    under prepped_ligands/<PDB>/, never the CWD.
    """
    try:
        line = f"{p.name}\t{reason}\n"
        base = Path(log_dir) if log_dir else (Path(MALFORMED_DIR) if MALFORMED_DIR else Path.cwd())
        dest = base / MALFORMED_LOG.name  # preserve original filename/timestamp pattern
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "a", encoding="utf-8") as out:
            out.write(line)
    except Exception:
        pass


def standardize_mol_with_activesite(mol: Chem.Mol) -> Chem.Mol:
    """
    Apply canonical ligand standardization before filtering:
      - RDKit metal disconnection, salt removal, largest fragment, uncharge
      - (light) sanitize to normalize valences where possible
      - Do NOT change protonation (keep existing behavior)
    Uses the same internal passes as _standardize_then_sanitize when available.
    """
    if mol is None:
        return mol
    m = mol
    # Prefer the existing standardization pipeline if available
    try:
        std, why = _standardize_then_sanitize(m)
        if std is not None:
            m = std
    except Exception:
        # fall through; keep original m
        pass

    # Guardrail: strip clearly invalid/dummy atoms if any slipped through
    try:
        bad = [a for a in m.GetAtoms() if (a.GetAtomicNum() == 0 or a.GetSymbol() == "*")]
        if bad:
            emsg = f"dummy_atoms({len(bad)})"
            m.SetProp("_ligprep_filter_hint", emsg)
    except Exception:
        pass

    return m
def _quick_filters(mol: Chem.Mol) -> Tuple[bool, str]:
    """
    Fast ligand gate: run *after* canonical standardization to avoid alias/case drift.
    Preserves legacy thresholds/messages; adds clearer diagnostics for element rejects.
    """
    # 1) Standardize first (aliases/salts) so checks run on a normalized molecule
    try:
        mol = standardize_mol_with_activesite(mol) or mol
    except Exception:
        # keep going; downstream checks are defensive
        pass

    # 2) Size/atom-count thresholds (unchanged)
    heavy = mol.GetNumHeavyAtoms()
    if heavy == 0:
        return False, "no_heavy_atoms"
    if heavy > MAX_HEAVY_ATOMS:
        return False, f"too_large({heavy})"
    total_atoms = mol.GetNumAtoms()
    if total_atoms < MIN_ATOMS_FOR_DOCKING:
        return False, f"too_few_atoms({total_atoms})"

    # 3) Element allow-list using aliases-derived set
    offending: Set[str] = set()
    try:
        for a in mol.GetAtoms():
            sym = a.GetSymbol()
            if sym not in ALLOWED_ELEMENTS:
                offending.add(sym)
    except Exception:
        # if RDKit access fails for any atom, be conservative and report failure
        return False, "element_scan_error"

    if offending:
        # Provide a hint for common alias/case drift
        sym = sorted(offending)[0]
        hint = ""
        if sym.upper() != sym and sym.capitalize() in ALLOWED_ELEMENTS:
            hint = f" (did_you_mean:{sym.capitalize()})"
        elif sym.upper() in (getattr(RULES, "two_letter", set()) or set()):
            hint = f" (did_you_mean:{sym[0].upper()}{sym[1:].lower()})"
        return False, f"disallowed_element:{sym}{hint}"

    # 4) Monatomic free-ion guard via aliases (canonical set; unions legacy)
    try:
        if mol.GetNumAtoms() == 1:
            s = mol.GetAtomWithIdx(0).GetSymbol()
            if s in MONOATOMIC_IONS:
                return False, f"free_monoatomic_ion:{s}"
    except Exception:
        pass

    return True, ""



def _standardize_then_sanitize(mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], str]:
    if not _HAS_STD:
        return None, "std_module_missing"
    try:
        md = _std.MetalDisconnector()
        fr = _std.FragmentRemover()
        lf = _std.LargestFragmentChooser(preferOrganic=True)
        uc = _std.Uncharger()

        m = md.Disconnect(mol)
        m = fr.RemoveFragments(m)
        m = lf.choose(m)
        m = uc.uncharge(m)

        try:
            m = Chem.AddHs(m, addCoords=True)
        except Exception as _e_hs:
            logging.warning("[ligprep] AddHs pre-sanitize failed: %s", _e_hs)

        Chem.SanitizeMol(m)

        ok, why = _quick_filters(m)
        if not ok:
            return None, why
        return m, ""
    except Exception as e:
        return None, f"std_resanitize_fail:{e}"


def _looks_like_buffer_salt(m: Chem.Mol) -> bool:
    from rdkit.Chem import rdMolDescriptors as rdmd
    hac = m.GetNumHeavyAtoms()
    o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'O')
    n = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'N')
    rings = rdmd.CalcNumRings(m)
    arom = rdmd.CalcNumAromaticRings(m)
    return (hac <= 15 and o >= 6 and n == 0 and rings == 0 and arom == 0)


# --- Counter-ion SMARTS (compiled once) ---
_COUNTERION_SMARTS = {
    "mesylate": Chem.MolFromSmarts("[CH3]-S(=O)(=O)[O-]"),
    "tosylate": Chem.MolFromSmarts("c1cccc(c1)S(=O)(=O)[O-]"),
    "triflate": Chem.MolFromSmarts("C(F)(F)F-S(=O)(=O)[O-]"),
    "sulfonate": Chem.MolFromSmarts("S(=O)(=O)[O-]"),
    "phosphate": Chem.MolFromSmarts("P(=O)([O-])([O-])[O-]"),
    "sulfate": Chem.MolFromSmarts("S(=O)(=O)([O-])[O-]"),
    "formate": Chem.MolFromSmarts("[#6](=O)[O-]"),
    "acetate": Chem.MolFromSmarts("CC(=O)[O-]"),
    "lactate": Chem.MolFromSmarts("CC(O)C(=O)[O-]"),
}
_COUNTERION_SMARTS.update({
    "citrate_like": Chem.MolFromSmarts("[CX4](-[CH2]-C(=O)[O-])(-[CH2]-C(=O)[O-])(-C(=O)[O-])O"),
    "tartrate_like": Chem.MolFromSmarts("IC([CH](O)C(=O)[O-])C(=O)[O-]".replace("I", "O")),
})


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
            if patt and m.HasSubstructMatch(patt):
                if hac <= 14:
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


# =========================
# Poly-acidic & buffer-like detection helpers
# =========================

_CARBOXYLATE = Chem.MolFromSmarts("[CX3](=O)[O-]")
_CARBOXYLIC = Chem.MolFromSmarts("[CX3](=O)O")


def _is_polyacidic_buffer_like(m: Chem.Mol) -> bool:
    try:
        from rdkit.Chem import rdMolDescriptors as rdmd
        hac = m.GetNumHeavyAtoms()
        if hac == 0:
            return False
        rings = rdmd.CalcNumRings(m)
        arom = rdmd.CalcNumAromaticRings(m)
        o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")
        n = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N")
        o_ratio = (o / float(hac)) if hac else 0.0

        na = 0
        if _CARBOXYLATE: na += len(m.GetSubstructMatches(_CARBOXYLATE))
        if _CARBOXYLIC:  na += len(m.GetSubstructMatches(_CARBOXYLIC))
        has_sulfate = bool(_COUNTERION_SMARTS["sulfate"] and m.HasSubstructMatch(_COUNTERION_SMARTS["sulfate"]))
        has_phosphate = bool(_COUNTERION_SMARTS["phosphate"] and m.HasSubstructMatch(_COUNTERION_SMARTS["phosphate"]))

        if rings == 0 and arom == 0 and o_ratio >= 0.35 and (na >= 3 or has_sulfate or has_phosphate):
            if n <= 1:
                return True
        return False
    except Exception:
        return False


def _element_from_adt(adt: str) -> str:
    t = (adt or "").strip()
    if not t:
        return "C"
    u = t.upper()
    MAP = {
        "C": "C", "A": "C",
        "N": "N", "NA": "N",
        "O": "O", "OA": "O",
        "S": "S", "SA": "S",
        "H": "H", "HD": "H",
        "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
        "P": "P",
        "B": "B", "SI": "Si", "SE": "Se",
        "ZN": "Zn", "MG": "Mg", "CA": "Ca", "MN": "Mn", "FE": "Fe",
        "K": "K", "NA+": "Na", "NA_": "Na", "NA ": "Na"
    }
    if u in MAP:
        return MAP[u]
    if u == "CL":
        return "Cl"
    if u == "BR":
        return "Br"
    return u[0]


def _get_ad4_types_from_aliases() -> set[str]:
    # Canonical AD4 atom type symbols (what ADT/Meeko/OBabel actually write)
    # Includes aromatic and hetero variants and standard polar hydrogens.
    return {
        # organics
        "C", "A",
        "N", "NA",
        "O", "OA",
        "S", "SA",
        "H", "HD",
        "P", "B",

        # halogens
        "F", "CL", "BR", "I",

        # metalloids / chalcogens
        "SI", "SE",

        # metals (commonly seen in PDBQTs)
        "ZN", "MG", "CA", "MN", "FE",
        "K", "NA", "CU", "CO", "NI", "AL", "AG", "AU", "PT",
        "LI", "BA", "SR", "CS", "RB",
    }


# Canonical AD4 type set (module-scope to avoid NameError anywhere)
try:
    AD4_TYPES: set[str] = _get_ad4_types_from_aliases()
except Exception:
    AD4_TYPES = {
        "C","A","N","NA","O","OA","S","SA","H","HD","F","CL","BR","I","P","B","SI","SE",
        "ZN","MG","CA","MN","FE","K","NA","CU","CO","NI","AL","AG","AU","PT","LI","BA","SR","CS","RB"
    }


# --- PDBQT metrics + invariants ---------------------------------------------
def _parse_pdbqt_metrics(pdbqt_path: Path) -> dict:
    """
    Lightweight scan of a PDBQT to extract invariants.
    - Only require polar H (HD) when chemically expected (NA/OA present).
    """
    try:
        n_atoms = 0
        typed = 0
        charges_present = 0
        torsdof = -1
        n_HD = 0
        n_NA = 0
        n_OA = 0
        bad_samples = []  # collect up to 10 non-AD4 tokens (token, line)

        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if ln.startswith("TORSDOF"):
                    try:
                        torsdof = int(ln.split()[-1])
                    except Exception:
                        torsdof = -1
                    continue
                if not ln.startswith(("ATOM", "HETATM")):
                    continue
                parts = ln.split()
                if len(parts) < 3:
                    continue

                n_atoms += 1

                # --- AD4 token: normalize common aliases and stray punctuation
                adt = parts[-1].strip().upper()
                # Strip trailing non-alphanumerics (rare stray chars)
                while adt and not adt[-1].isalnum():
                    adt = adt[:-1]
                # Normalize a couple of noisy aliases we’ve seen in the wild
                if adt in {"NA+", "NA_", "NA"}:
                    adt = "NA"

                if adt in AD4_TYPES:
                    typed += 1
                else:
                    if len(bad_samples) < 10:
                        bad_samples.append((adt, ln.strip()))

                if adt == "HD":
                    n_HD += 1
                elif adt == "NA":
                    n_NA += 1
                elif adt == "OA":
                    n_OA += 1

                # per-atom partial charge should be the penultimate token
                try:
                    _ = float(parts[-2])
                    charges_present += 1
                except Exception:
                    pass

        has_ad4_types = (n_atoms > 0 and typed == n_atoms)
        has_charges = (n_atoms > 0 and charges_present == n_atoms)
        polar_expected = (n_NA + n_OA) > 0
        has_polar_H = (n_HD > 0) or (not polar_expected)

        # Debug/audit when types are missing
        if (not has_ad4_types) and bad_samples and os.environ.get("LIGPREP_DEBUG", "").strip().lower() in {"1","true","yes","y"}:
            sample_str = "; ".join([f"{tkn}:{ln}" for tkn, ln in bad_samples[:10]])
            logging.warning("[ad4.audit] %s bad_types sample=%s", pdbqt_path.stem, sample_str)

        return {
            "has_polar_H": has_polar_H,
            "polarH_expected": polar_expected,
            "has_ad4_types": has_ad4_types,
            "has_charges": has_charges,
            "torsdof": torsdof,
            # expose counts when debugging downstream
            "n_atoms": n_atoms, "n_HD": n_HD, "n_NA": n_NA, "n_OA": n_OA, "typed": typed, "charges_present": charges_present,
        }
    except Exception as e:
        # Try to include a tiny context without risking new exceptions
        ctx = ""
        try:
            lines = []
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if ln.startswith(("ATOM","HETATM")):
                        lines.append(ln.strip())
                        if len(lines) >= 3:
                            break
            if lines:
                ctx = " | ctx=" + " | ".join(lines)
        except Exception:
            pass
        logging.warning("[ligprep] metrics parse failed for %s: %s%s", pdbqt_path.name, e, ctx)
        return {
            "has_polar_H": False,
            "polarH_expected": False,
            "has_ad4_types": False,
            "has_charges": False,
            "torsdof": -1,
        }



def validate_pdbqt_invariants(pdbqt_path: Path) -> tuple[bool, dict, list[str]]:
    m = _parse_pdbqt_metrics(pdbqt_path)
    failures: list[str] = []

    # Stop failing on missing polar H; donors are optional at accept stage.
    # if m.get("polarH_expected", False) and not m.get("has_polar_H", False):
    #     failures.append("missing_polar_H")

    if not m.get("has_ad4_types", False):
        failures.append("missing_AD4_types")
    if not m.get("has_charges", False):
        failures.append("missing_charges")

    td = m.get("torsdof", -1)
    # Allow rigid ligands (TORSDOF == 0). Fail only if missing or negative.
    if not (isinstance(td, int) and td >= 0):
        failures.append("torsdof_missing_or_negative")

    if os.environ.get("LIGPREP_DEBUG", "").strip() in {"1","true","yes","y"}:
        logging.info(
            "[invariants] lig=%s has_HD=%s donors_seen=%s has_polar_H=%s has_ad4_types=%s has_charges=%s "
            "torsdof=%s fail=%s",
            pdbqt_path.stem, m.get("n_HD","NA"), m.get("n_NA","NA"), m.get("has_polar_H"),
            m.get("has_ad4_types"), m.get("has_charges"), m.get("torsdof"), failures
        )
    return (len(failures) == 0), m, failures





def _adt_retype_and_normalize(
        mgltools_python_short: str,
        prepare_script_short: str,
        mol2_for_mgl: Path,
        out_pdbqt: Path,
        *,
        torsion_rule_label: str = "adt_normalized"
) -> tuple[bool, str]:
    """
    Re-run prepare_ligand4 to enforce AD4 types and ADT torsion ROOT/BRANCH, then replace out_pdbqt on success.
    """
    try:
        tmp = out_pdbqt.with_suffix(".adt_norm.pdbqt")
        cmd = [
            mgltools_python_short, prepare_script_short,
            "-l", get_short_path_name(str(mol2_for_mgl.resolve())),
            "-o", get_short_path_name(str(tmp.resolve())),
            "-U", "nphs_lps",
            "-A", "hydrogens",
        ]
        res = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=str(mol2_for_mgl.parent), timeout=600)
        if res.stderr:
            logging.warning("[mgltools retype stderr] %s", res.stderr.strip()[:300])

        if tmp.exists() and tmp.stat().st_size > 100:
            try:
                if out_pdbqt.exists():
                    out_pdbqt.unlink()
            except Exception:
                pass
            tmp.replace(out_pdbqt)

            m = _parse_pdbqt_metrics(out_pdbqt)
            typed, n_atoms = m.get("typed", 0), m.get("n_atoms", 0)
            logging.info("[retype.ok] replaced=True n_atoms=%d typed_ratio=%d/%d ad4_ok=%s torsdof=%s",
                         n_atoms, typed, n_atoms, m.get("has_ad4_types"), m.get("torsdof"))
            return True, torsion_rule_label

    except subprocess.CalledProcessError as e:
        logging.warning("[torsion-normalize] ADT retype failed: %s", (e.stderr or "")[-180:])
    except Exception as e:
        logging.warning("[torsion-normalize] exception: %s", e)

    return False, torsion_rule_label


def _reprep_via_meeko_path(mol2_path: Path, out_pdbqt: Path) -> bool:
    """
    Use Meeko from PATH (no config): prefer `mk_prepare_ligand.py`; fallback to `python -m meeko`.
    """
    exe = shutil.which("mk_prepare_ligand.py")
    if exe:
        cmd = [exe, "-i", str(mol2_path), "-o", str(out_pdbqt)]
    else:
        # last-resort: try the active interpreter with module form
        cmd = [sys.executable, "-m", "meeko", "-i", str(mol2_path), "-o", str(out_pdbqt)]
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
        if res.stderr:
            logging.warning("[meeko stderr] %s", res.stderr.strip()[:300])
        return out_pdbqt.exists() and out_pdbqt.stat().st_size > 100
    except Exception as e:
        logging.warning("[meeko] failed: %s", e)
        return False


def _element_counts_from_pdb(pdb_path: Path) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fin:
            for ln in fin:
                if not ln.startswith(("ATOM", "HETATM")):
                    continue
                parts = ln.split()
                if not parts:
                    continue
                adt = parts[-1]
                elem = _element_from_adt(adt)
                counts[elem] += 1
    except Exception:
        pass
    return counts


import logging
from collections import Counter


def _elem_hist(m):
    try:
        return dict(Counter(a.GetSymbol() for a in m.GetAtoms()))
    except Exception:
        return {}


def _formal_charge(m):
    try:
        return sum(a.GetFormalCharge() for a in m.GetAtoms())
    except Exception:
        return 0


def _has_explicit_Hs(m):
    try:
        return any(a.GetAtomicNum() == 1 for a in m.GetAtoms())
    except Exception:
        return False


def _buffer_like_by_counts_from_mol(m: Chem.Mol) -> bool:
    hac = m.GetNumHeavyAtoms()
    if hac == 0:
        return False
    o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'O')
    n = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'N')
    return (hac >= 10 and (o / float(hac)) >= 0.40 and n <= 1)


def _buffer_like_by_counts_from_pdbfile(pdb_path: Path) -> bool:
    counts = _element_counts_from_pdb(pdb_path)
    hac = sum(v for k, v in counts.items() if k != "H")
    o = counts.get("O", 0)
    n = counts.get("N", 0)
    return (hac >= 10 and hac > 0 and (o / float(hac)) >= 0.40 and n <= 1)


def _polyacidic_by_counts_from_pdbfile(pdb_path: Path) -> bool:
    counts = _element_counts_from_pdb(pdb_path)
    hac = sum(v for k, v in counts.items() if k != "H")
    o = counts.get("O", 0)
    n = counts.get("N", 0)
    return (hac >= 10 and hac > 0 and (o / float(hac)) >= 0.35 and n <= 1)


def _write_aromatic_sdf(mol: Chem.Mol, out_path: Path) -> bool:
    try:
        m = Chem.Mol(mol)
        try:
            Chem.SanitizeMol(m)
        except Exception:
            # If sanitize still fails, attempt minimal aromatic flag set and write anyway
            pass
        try:
            Chem.SetAromaticity(m, Chem.AromaticityModel.AROMATICITY_RDKIT)
        except Exception:
            pass

        w = Chem.SDWriter(str(out_path))
        try:
            w.SetKekulize(False)  # keep aromatic bonds marked if possible
        except Exception:
            pass
        w.write(m)
        w.close()
        return True
    except Exception as e:
        logging.warning(f"[ligprep] write aromatic SDF failed for {out_path.name}: {e}")
        return False


from rdkit import Chem
from rdkit.Chem import rdchem


def _re_aromatize_mol2_in_place(mol2_path: Path, obabel_exe_short: str) -> bool:
    """
    Try to recover aromaticity deterministically by: RDKit (no sanitize) -> SDF write with aromatic flags
    -> OBabel SDF->MOL2 round-trip. Emit a compact audit once.
    """
    m = None
    try:
        m = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
    except Exception as e:
        logging.warning("[ligprep] RDKit failed to read MOL2 (sanitize=False) %s: %s", mol2_path.name, e)

    if m is None:
        logging.info("[ligprep] re_arom skip (RDKit load failed) mol2=%s", mol2_path.name)
        return False

    # Count aromatic atoms BEFORE rescue
    try:
        m.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    try:
        Chem.GetSymmSSSR(m)
    except Exception:
        pass
    try:
        Chem.SetAromaticity(m, Chem.AromaticityModel.AROMATICITY_RDKIT)
    except Exception:
        pass
    try:
        before_arom = sum(int(a.GetIsAromatic()) for a in m.GetAtoms())
    except Exception:
        before_arom = -1

    tmp_sdf = mol2_path.with_suffix(".arom.sdf")
    if not _write_aromatic_sdf(m, tmp_sdf):
        logging.info("[ligprep] re_arom SDF write failed for %s", mol2_path.name)
        return False

    tmp_mol2 = mol2_path.with_suffix(".arom.mol2")
    ok_ob = _run_obabel([obabel_exe_short, "-isdf", str(tmp_sdf), "-omol2", "-O", str(tmp_mol2)], timeout_sec=600)
    if not ok_ob:
        logging.info("[ligprep] re_arom OBabel round-trip failed for %s", mol2_path.name)
        return False

    if tmp_mol2.exists() and tmp_mol2.stat().st_size > 100:
        try:
            mol2_path.unlink(missing_ok=True)
            tmp_mol2.replace(mol2_path)
            # post-replace: recount aromatics
            after_arom = _count_aromatic_atoms_in_mol2(mol2_path)
            used = bool(after_arom >= 0 and before_arom >= 0 and after_arom != before_arom)
            logging.info("[arom-rescue] file=%s before=%d after=%d used=%s",
                         mol2_path.name, before_arom, after_arom, used)
            return True
        except Exception as e:
            logging.info("[ligprep] re_arom replace failed for %s: %s", mol2_path.name, e)
            return False

    logging.info("[ligprep] re_arom produced empty MOL2 for %s", mol2_path.name)
    return False



_SANITIZED_RE = re.compile(r'(?:\.sanitized)+(?=\.)')


def collapse_sanitized_once(p: Path) -> Path:
    return p.with_name(_SANITIZED_RE.sub('.sanitized', p.name))


def add_hydrogens_mol2(in_path: Path, out_path: Path, obabel_exe: str) -> tuple[bool, str]:
    """
    Ensure explicit H before MGLTools. Returns (ok, stderr_text).
    Uses obabel -h to add hydrogens *without* changing atom order more than necessary.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [obabel_exe, "-imol2", str(in_path), "-omol2", "-O", str(out_path), "-h"]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True, (res.stderr or "").strip()
    except subprocess.CalledProcessError as e:
        return False, (e.stderr or "").strip()
def _log_std_diff(log_dir: Path, lig_name: str, branch: str, old_smiles: str, new_smiles: str) -> None:
    """Append SMILES diffs caused by standardization so we can spot chemistry changes."""
    try:
        path = Path(log_dir) / "standardization_diffs.tsv"
        if not path.exists():
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("ligand\tbranch\told_smiles\tnew_smiles\n")
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(f"{lig_name}\t{branch}\t{old_smiles}\t{new_smiles}\n")
    except Exception:
        pass


def quick_pdbqt_validate(path: Path) -> tuple[bool, int, str]:
    """
    Lightweight validity: file exists, non-trivial size, has ATOM/HETATM, no blatantly invalid ADT types.
    Returns (ok, atom_count, reason_if_fail)
    """
    if not path.exists() or path.stat().st_size < 100:
        return False, 0, "empty_or_missing"
    atoms = 0
    bad = 0
    with open(path, "r", errors="ignore") as fh:
        for ln in fh:
            if ln.startswith("ATOM") or ln.startswith("HETATM"):
                atoms += 1
                # simple ADT type sanity: last token should exist
                if len(ln.split()) < 2:  # too sparse
                    bad += 1
    if atoms == 0:
        return False, 0, "no_atoms"
    if bad > 0:
        return False, atoms, "malformed_records"
    return True, atoms, ""



# =========================
# Parent chooser
# =========================

def _to_parent_mol(m: Chem.Mol) -> Optional[Chem.Mol]:
    try:
        if _HAS_STD and _std is not None:
            p = _std.Cleanup(m)
            chooser = _std.LargestFragmentChooser(preferOrganic=True)
            p = chooser.choose(p)
            p = _std.ChargeParent(p)
            if _looks_like_buffer_salt(p) or _matches_counterion(p) or _is_polyacidic_buffer_like(
                    p) or _buffer_like_by_counts_from_mol(p):
                frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
                p = max(frags, key=lambda x: x.GetNumHeavyAtoms())
                p = _std.ChargeParent(p)
        else:
            p = _REM.StripMol(m, dontRemoveEverything=True)
            frags = Chem.GetMolFrags(p, asMols=True, sanitizeFrags=True)
            p = max(frags, key=lambda x: x.GetNumHeavyAtoms())

        Chem.SanitizeMol(p)
        try:
            Chem.SetAromaticity(p, Chem.AromaticityModel.AROMATICITY_RDKIT)
        except Exception:
            pass

        if not any(a.GetSymbol() == "C" for a in p.GetAtoms()):
            return None
        if _is_polyacidic_buffer_like(p) or _buffer_like_by_counts_from_mol(p):
            return None
        if p.GetNumHeavyAtoms() < MIN_PARENT_HEAVY:
            return None
        if _looks_like_buffer_salt(p) or _matches_counterion(p):
            return None
        return p
    except Exception:
        return None


# Intial SCAM Filter (post-cleaning)
from rdkit.Chem import Crippen

# Optional pKa/logD dependency: try to import, else fall back to cLogP
try:
    from pkasolver import pkasolver as _pka

    _HAS_PKASOLVER = True
except Exception:
    _pka = None
    _HAS_PKASOLVER = False

# SCAM substructure alerts
SCAM_SMARTS: Dict[str, str] = {
    # electrophiles / reactive
    "epoxide": "[OX2r3]",
    "aziridine": "[NX3r3]",
    "alkyl_halide": "[CX4;H0,H1,H2][Cl,Br,I,F]",
    "michael_acceptor": "[C,c]=[C,c]-[C,S](=O)[O,N,S] | [C,c]=[C,c]-C(=O)[O,N,S]",
    "acrylamide": "C=CC(=O)N",
    "isothiocyanate": "N=C=S",
    "sulfonyl_fluoride": "S(=O)(=O)F",

    # redox / interference
    "p_quinone": "O=C1C=CC(=O)C=C1",
    "o_quinone": "O=c1ccc(=O)[cH][cH]1",
    "phenothiazine_like": "n2c1ccccn1Sc3ccccc23",

    # chelators / aggregators
    "catechol": "c1cc(O)c(O)cc1",
    "hydroxamate": "C(=O)N[OH]",
    "8_hydroxyquinoline": "Oc1cccc2ncccc12",
    "tannin_polyphenol": "c(O)c(O)c(O)",

    # nucleophiles / potentially reactive
    "hydrazine": "NN",
    "hydroxylamine": "N[OH]",
    "thiol": "[SH]",
    "dithiol": "SCCS",

    # rhodanine / known hitters
    "rhodanine": "O=C1NC(=S)SC1",
    "barbiturate_like": "O=C1NC(=O)NC(=O)1",
}
SCAM_QUERIES = {name: Chem.MolFromSmarts(s) for name, s in SCAM_SMARTS.items()}


def scam_flags(mol: Chem.Mol) -> List[str]:
    flags = []
    for name, patt in SCAM_QUERIES.items():
        if patt is not None and mol.HasSubstructMatch(patt):
            flags.append(f"hard: {name}")
    return flags


def predict_logD(mol: Chem.Mol, ph: float = 7.4) -> float:
    """
    Prefer pkasolver logD if available, otherwise gracefully fall back to cLogP.
    """
    try:
        if _HAS_PKASOLVER and _pka is not None:
            smiles = Chem.MolToSmiles(mol)
            return float(_pka.calculate_logd(smiles, ph=ph))
    except Exception as e:
        print(f"[WARN] logD prediction failed for {Chem.MolToSmiles(mol)}: {e}")
    return float(Crippen.MolLogP(mol))


def annotate_ligand_with_scam(lig_path: str, ligand_record: Dict) -> Dict:
    """
    Annotates the ligand record with SCAM filter results (cLogP, logD7.4, SCAM_Flags, SMILES).
    """
    try:
        # load molecule from file
        if lig_path.endswith(".sdf") or lig_path.endswith(".mol"):
            mol = Chem.MolFromMolFile(lig_path, sanitize=True)
        elif lig_path.endswith(".mol2"):
            mol = Chem.MolFromMol2File(lig_path, sanitize=True)
        else:
            raise ValueError(f"Unsupported ligand format: {lig_path}")

        if mol is None:
            ligand_record["SCAM_Flags"] = "InvalidMol"
            return ligand_record

        # get ligand descriptors
        smiles = Chem.MolToSmiles(mol)
        clogp = Crippen.MolLogP(mol)
        logd = predict_logD(mol, ph=7.4)

        # check for typical SCAM substructures
        flags = scam_flags(mol)

        # soft flags based on property cutoffs
        if clogp is not None and clogp > 3.5:
            flags.append(f"soft:high_logP({clogp:.2f})")

        if logd is not None:
            if logd < 0:
                flags.append(f"soft:low_logD({logd:.2f})")
            elif logd > 3.5:
                flags.append(f"soft:high_logD({logd:.2f})")
            if logd > 5:
                flags.append(f"hard:very_high_logD({logd:.2f})")

        flag_str = ";".join(flags) if flags else "None"

        ligand_record.update({
            "SMILES": smiles,
            "cLogP": clogp,
            "logD7.4": logd,
            "SCAM_Flags": flag_str
        })

    except Exception as e:
        ligand_record["SCAM_Flags"] = f"Error:{e}"

    return ligand_record


def _append_prep_status(
        status_log_path: Path,
        ligand_name: str,
        status: str,
        reason: str = "",
        relpath: str = "",
        *,
        stage: str = "",
        failure_code: str = "",
        failure_detail: str = "",
        fixes_count: int = 0,
        rules_ver: str = "",
        writer_final: str = "",
        rescue_used: str = "",
        torsion_root_rule: str = "",
        polarH: str = "",
        ad4_types_ok: str = "",
        charges_ok: str = "",
        torsdof: str = "",
):
    """
    Backward-compatible TSV: old columns + appended provenance/metrics.
    """
    header = (
        "ligand\tstatus\treason\tpdbqt_rel\t"
        "stage\tfailure_code\tfailure_detail\tfixes_count\trules_version\t"
        "writer_final\trescue_used\ttorsion_root_rule\tpolarH\tad4_types_ok\tcharges_ok\ttorsdof\n"
    )
    if not status_log_path.exists():
        status_log_path.write_text(header, encoding="utf-8")
    with status_log_path.open("a", encoding="utf-8") as fh:
        row = [
            ligand_name, status, reason, relpath,
            stage, failure_code, failure_detail, str(fixes_count), rules_ver,
            writer_final, rescue_used, torsion_root_rule, polarH, ad4_types_ok, charges_ok, str(torsdof or "")
        ]
        fh.write("\t".join(row) + "\n")


# --- logging helper for element-fix one-liners ---
def _log_elem_fix_summary(file_path: Path, stage: str, before_text: str = None, after_text: str = None):
    """
    Emit a single compact line showing how many He->H (ADT/element) conversions were performed
    by comparing counts before/after. If before_text is None, read file.
    """
    try:
        if before_text is None:
            before_text = Path(file_path).read_text(encoding="utf-8", errors="ignore")
        n_before = scan_helium_counts(before_text.splitlines())
        if after_text is None:
            after_text = Path(file_path).read_text(encoding="utf-8", errors="ignore")
        n_after = scan_helium_counts(after_text.splitlines())
        # If we *re-wrote* columns (He->H), the *remaining* He count should drop. Report delta.
        delta = max(0, n_before - n_after)
        logging.info(f"[elem-fix] file={file_path.name} stage={stage} He->H={delta}")
    except Exception as e:
        logging.warning(f"[elem-fix] summary failed for {file_path}: {e}")


def _helium_postwrite_guard(pdbqt_path: Path, ligand_name: str, prepped_ligands_dir: Path) -> Tuple[bool, int, str]:
    try:
        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
        new_lines, fixes, q_reason = assert_no_helium_in_pdbqt(lines, ligand_name)
        logging.info("[helium] ligand=%s fixes=%d quarantine=%s rules=%s",
                     ligand_name, fixes, str(bool(q_reason)), rules_version())

        if fixes > 0:
            with open(pdbqt_path, "w", encoding="utf-8") as out:
                out.writelines(new_lines)

        ok = (q_reason == "")
        print(f"[helium] ligand={ligand_name} fixes={fixes} quarantine={not ok} rules={rules_version()}")
        return ok, fixes, q_reason
    except Exception as e:
        print(
            f"[helium] ligand={ligand_name} fixes=0 quarantine=True rules={rules_version()} note=postwrite_exception:{e}")
        return False, 0, f"postwrite_exception:{e}"


# =========================
# OBabel-friendly SDF writer
# =========================

def _write_obabel_friendly_sdf(mol: Chem.Mol, out_path: Path) -> bool:
    try:
        m = Chem.Mol(mol)
        Chem.SanitizeMol(m)
        from rdkit.Chem import AllChem
        try:
            # Helps fix some O valence issues due to charge misassignment
            Chem.Kekulize(m, clearAromaticFlags=True)
        except Exception:
            pass
        try:
            rdMolOps.AssignFormalCharges(m)
        except Exception:
            pass
        w = Chem.SDWriter(str(out_path))
        w.write(m)
        w.close()
        return True
    except Exception:
        try:
            m2 = Chem.Mol(mol)
            Chem.SanitizeMol(m2)
            Chem.SetAromaticity(m2, Chem.AromaticityModel.AROMATICITY_RDKIT)
            w = Chem.SDWriter(str(out_path))
            try:
                w.SetKekulize(False)
            except Exception:
                pass
            w.write(m2)
            w.close()
            return True
        except Exception:
            return False


# =========================
# OBabel helpers
# =========================

def _reserialize_mol_via_obabel(mol, obabel_exe_short: str, target_mol2: Path):
    tmp_sdf = target_mol2.with_suffix(".std.sdf")
    if not _write_obabel_friendly_sdf(mol, tmp_sdf):
        _log_malformed(target_mol2, "rdkit_sdf_write_fail:kekulize_or_aromaticity", log_dir=MALFORMED_DIR)
        return None

    fresh_mol2 = target_mol2.with_suffix(".std.mol2")
    cmd = [obabel_exe_short, "-isdf", str(tmp_sdf), "--gen3d", "-omol2", "-O", str(fresh_mol2)]
    try:
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


def read_config(path="config.txt") -> Dict[str, str]:
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


# =========================
# Open Babel runners
# =========================

def _run_obabel(cmd: List[str], timeout_sec: int) -> bool:
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


def _attempt_obabel_series(base_cmd: List[str], timeout_sec: int, threads_list: List[int],
                           use_fast_first: bool = True) -> bool:
    attempts: List[List[str]] = []
    if use_fast_first:
        for j in threads_list:
            attempts.append(base_cmd + ["--fast", "-j", str(j)])
    for j in threads_list:
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


# =========================
# SDF chunking
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
        if cur:
            mol_buf.extend(cur)
        flush_chunk()
    return chunks


# =========================
# SDF ? MOL2 conversions
# =========================

def convert_sdf_to_mol2_split_parallel(
        sdf_path: Path,
        mol2_output_dir: Path,
        obabel_exe: str,
        threads: int = 32,
        timeout_sec: int = 36000,
        chunk_size: int = 1000
) -> List[Path]:
    mol2_output_dir.mkdir(parents=True, exist_ok=True)

    chunk_dir = mol2_output_dir / "_sdf_chunks"
    if chunk_dir.exists():
        shutil.rmtree(chunk_dir)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    n_mols = count_sdf_records(sdf_path)
    print(f"SDF has ~{n_mols} molecules")

    chunk_paths = split_sdf_into_chunks(sdf_path, chunk_dir, chunk_size=chunk_size)
    print(f"Split into {len(chunk_paths)} chunk(s) of up to {chunk_size} molecules")

    all_out: List[Path] = []

    t1 = threads if threads >= 16 else 16
    t2 = max(16, t1 // 2)
    threads_list = []
    for t in [t1, t2, 16]:
        if t not in threads_list:
            threads_list.append(t)

    for ci, chunk in enumerate(chunk_paths, 1):
        prefix = mol2_output_dir / f"mol2_chunk{ci:04d}_"
        base_cmd = [
            obabel_exe, "-isdf", str(chunk),
            "--gen3d", "-omol2",
            "-m", "-O", str(prefix) + ".mol2"
        ]
        print(f"[ligprep] pre-MOL2-write: chunk={ci} sdf={chunk.name} -> prefix={prefix.name}")
        ok = _attempt_obabel_series(
            base_cmd, timeout_sec=timeout_sec, threads_list=threads_list, use_fast_first=True
        )
        if not ok:
            print(f"Chunk {ci} failed entirely; moving on.")
            continue

        out_files = sorted(mol2_output_dir.glob(f"mol2_chunk{ci:04d}_*.mol2"))
        print(f"Chunk {ci}: wrote {len(out_files)} mol2 files")
        for _mol2 in out_files[:4]:  # cap spam
            try:
                n_mol2_atoms = 0
                with open(_mol2, "r", errors="ignore") as fh:
                    in_atoms = False
                    for ln in fh:
                        s = ln.strip()
                        if s.startswith("@<TRIPOS>ATOM"):
                            in_atoms = True
                            continue
                        if s.startswith("@<TRIPOS>") and in_atoms:
                            break
                        if in_atoms and s and s[0].isdigit():
                            n_mol2_atoms += 1
                print(f"[ligprep] MOL2_written≈{n_mol2_atoms} file={_mol2.name}")
            except Exception as e:
                print(f"[ligprep] MOL2_count_error {_mol2.name}: {e}")
        all_out.extend(out_files)

    print(f"Total MOL2 files: {len(all_out)}")
    return all_out


def _obabel_convert_chunk(chunk_sdf: Path, out_prefix: Path, obabel_exe: str) -> int:
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


# =========================
# RDKit ETKDG path (with parent picking)
# =========================
def rdkit_embed_sdf_to_mol2(
        sdf_in: Path, mol2_out_dir: Path, obabel_exe: str, max_workers: int = 8
) -> List[Path]:
    # Ensure both imports exist in this scope
    from rdkit import Chem
    from rdkit.Chem import AllChem

    suppl = Chem.SDMolSupplier(str(sdf_in), removeHs=False, sanitize=False)
    mols = [(i, m) for i, m in enumerate(suppl) if m is not None]
    print(f"RDKit: loaded {len(mols)} molecules from {sdf_in.name}")

    sdf_tmp_dir = mol2_out_dir / "_rdkit_embedded_sdf"
    sdf_tmp_dir.mkdir(parents=True, exist_ok=True)

    def _embed_one(i_m):
        i, m = i_m
        from rdkit.Chem import rdmolops as rdMolOps

        logging.info("[bulkSDF] pre-standardize ligand=rdk_%07d atoms=%d charge=%+d elements=%s has_conf=%s",
                     i,
                     (m.GetNumAtoms() if m else -1),
                     _formal_charge(m),
                     _elem_hist(m),
                     bool(m and m.GetNumConformers() > 0))

        try:
            Chem.SanitizeMol(m)
            logging.info("[bulkSDF] post-sanitize ligand=rdk_%07d ok=True charge=%+d", i, _formal_charge(m))
        except Exception as e:
            print(f"[ligprep] SanitizeMol hard fail: {e} (trying staged flags)")
            try:
                rdMolOps.SanitizeMol(m, sanitizeOps=rdMolOps.SanitizeFlags.SANITIZE_FINDRADICALS |
                                                    rdMolOps.SanitizeFlags.SANITIZE_KEKULIZE |
                                                    rdMolOps.SanitizeFlags.SANITIZE_SETAROMATICITY |
                                                    rdMolOps.SanitizeFlags.SANITIZE_ADJUSTHS)
                print("[ligprep] staged-sanitize succeeded")
                logging.info("[bulkSDF] post-sanitize ligand=rdk_%07d ok=True(staged) charge=%+d", i, _formal_charge(m))
            except Exception as e2:
                _dbg_atom_stats(m, "sanitize-fail-snapshot")
                logging.info("[bulkSDF] post-sanitize ligand=rdk_%07d ok=False err=%s", i, str(e2)[:180])
                raise  # keep failing fast, but with context

        orig_heavy = m.GetNumHeavyAtoms()
        p = _to_parent_mol(m)
        if p is None:
            _log_malformed(Path(f"rdk_{i:07d}"), "no_parent_or_too_small_after_desalting")
            return None
        if p is not m and p.GetNumHeavyAtoms() != orig_heavy:
            logging.info(f"[parent-pick] rdk_{i:07d}: {orig_heavy}?{p.GetNumHeavyAtoms()} heavy atoms")

        m = p

        logging.info("[bulkSDF] post-standardize ligand=rdk_%07d atoms=%d charge=%+d elements=%s",
                     i, m.GetNumAtoms(), _formal_charge(m), _elem_hist(m))

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
            if _write_obabel_friendly_sdf(m, out):
                logging.info("[paths] sdf_out=%s", str(out.resolve()))

                return out

            else:
                _log_malformed(out, "write_obabel_friendly_sdf_failed")
                return None
        except Exception:
            return None

    paths: List[Path] = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        for res in ex.map(_embed_one, mols):
            if res:
                paths.append(res)

    print(f"RDKit: embedded {len(paths)} molecules; converting to MOL2")
    out_files: List[Path] = []
    for pth in paths:
        out = mol2_out_dir / (pth.stem + ".mol2")
        print(f"[ligprep] pre-MOL2-write: rdkit_embed sdf={pth.name} -> mol2={out.name}")
        logging.info("[bulkSDF] write.mol2 ligand=%s via=obabel out=%s", pth.stem, str(out.resolve()))

        if out.exists() and out.stat().st_size > 100:
            out_files.append(out)
            try:
                size = out.stat().st_size
                logging.info("[bulkSDF] write.mol2.done ligand=%s rc=0(size_only) size=%d", pth.stem, size)
            except Exception:
                pass
            continue

        cmd = [obabel_exe, "-isdf", str(pth), "-omol2", "-O", str(out)]
        if _run_obabel(cmd, timeout_sec=120):
            out_files.append(out)
            try:
                size = out.stat().st_size if out.exists() else 0
                logging.info("[bulkSDF] write.mol2.done ligand=%s rc=%s size=%d",
                             pth.stem, "0" if size > 0 else "unknown", size)
            except Exception:
                pass

            try:
                n_mol2_atoms = 0
                with open(out, "r", errors="ignore") as fh:
                    in_atoms = False
                    for ln in fh:
                        s = ln.strip()
                        if s.startswith("@<TRIPOS>ATOM"):
                            in_atoms = True
                            continue
                        if s.startswith("@<TRIPOS>") and in_atoms:
                            break
                        if in_atoms and s and s[0].isdigit():
                            n_mol2_atoms += 1
                print(f"[ligprep] MOL2_written≈{n_mol2_atoms} file={out.name}")
            except Exception as e:
                print(f"[ligprep] MOL2_count_error {out.name}: {e}")

    print(f"RDKit: wrote {len(out_files)} MOL2 files")
    return out_files


# =========================
# Optional: pre-clean entire SDF to parents only
# =========================

def _write_parent_only_sdf(src_sdf: Path, dst_sdf: Path) -> int:
    suppl = Chem.SDMolSupplier(str(src_sdf), removeHs=False, sanitize=False)
    w = Chem.SDWriter(str(dst_sdf))
    try:
        w.SetKekulize(False)
    except Exception:
        pass
    kept = 0
    for m in suppl:
        if m is None:
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            continue
        p = _to_parent_mol(m)
        if p is None:
            continue
        w.write(p)
        kept += 1
    w.close()
    return kept


# =========================
# MGLTools: prepare_ligand4 & validation
# =========================

def _is_probable_water_from_pdbqt_lines(lines: List[str]) -> bool:
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    if not atom_lines:
        return False
    if any(" HOH " in ln for ln in atom_lines):
        return True
    heavy = 0
    has_O = False
    for ln in atom_lines:
        parts = ln.split()
        if not parts:
            continue
        elem = _element_from_adt(parts[-1])
        if elem == 'O':
            has_O = True
        if elem != 'H':
            heavy += 1
    return has_O and heavy <= 1


def is_valid_ligand(path: Path, log_dir: Path) -> bool:
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        atom_lines = [line for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]
        torsion_lines = [line for line in lines if line.startswith("TORSDOF")]

        if not atom_lines:
            raise ValueError("No ATOM or HETATM lines found.")
        if not torsion_lines:
            raise ValueError("No torsion info (TORSDOF) found.")
        if _looks_like_monoatomic_ion_pdbqt(lines):
            raise ValueError("monoatomic_counterion")
        if _is_probable_water_from_pdbqt_lines(lines):
            raise ValueError("looks_like_water_or_hydroxide")
        heavy = 0
        for ln in atom_lines:
            parts = ln.split()
            if not parts:
                continue
            if _element_from_adt(parts[-1]) != 'H':
                heavy += 1
        if heavy < 3:
            raise ValueError(f"too_few_heavy_atoms_in_pdbqt({heavy})")

        return True
    except Exception as e:
        malformed_log = Path(log_dir) / MALFORMED_LOG.name
        with open(malformed_log, "a", encoding="utf-8") as f:
            f.write(f"{path.name} - PDBQT validation failed: {e}\n")
        return False

def _audit_protonation_metrics(tag: str,
                               mol_or_path: Union["Chem.Mol", str, Path],
                               context: str = "pre") -> Dict[str, Any]:
    """
    Returns {'ok': bool, 'h_count': int, 'formal_charge': int, 'has_partial_charges': bool,
             'aromatic_atoms': int, 'ad4_types_ok': Optional[bool], 'notes': List[str]}
    Works on RDKit Mol or path to PDB/PDBQT/MOL2/SDF.
    """
    notes: List[str] = []
    h_count = 0
    formal_charge = 0
    has_partial = False
    aromatic_atoms = 0
    ad4_ok: Optional[bool] = None

    mol = None
    try:
        from rdkit import Chem
        if hasattr(mol_or_path, "GetNumAtoms"):
            mol = mol_or_path
        else:
            p = str(mol_or_path)
            ext = os.path.splitext(p)[1].lower()
            if ext in {".sdf"}:
                suppl = Chem.SDMolSupplier(p, sanitize=False, removeHs=False)
                mol = next((m for m in suppl if m), None)
            elif ext in {".mol2"}:
                mol = Chem.MolFromMol2File(p, sanitize=False, removeHs=False)
            elif ext in {".pdb"}:
                mol = Chem.MolFromPDBFile(p, sanitize=False, removeHs=False)
        if mol:
            try:
                Chem.SanitizeMol(mol, catchErrors=True)
            except Exception:
                notes.append("sanitize_warn")
            try:
                formal_charge = Chem.GetFormalCharge(mol)
            except Exception:
                notes.append("formal_charge_calc_warn")
            try:
                h_count = sum(a.GetNumExplicitHs() + a.GetTotalNumHs() for a in mol.GetAtoms())
            except Exception:
                notes.append("h_count_warn")
            try:
                aromatic_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
            except Exception:
                notes.append("aroma_warn")
    except Exception:
        notes.append("rdkit_missing")

    # Fallback for PDBQT: detect AD4 types + partial charges heuristically
    if not mol and isinstance(mol_or_path, (str, Path)):
        try:
            has_q = False
            ad4_seen = False
            h_count_fallback = 0
            with open(mol_or_path, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                        continue
                    an = ln[12:16].strip().upper()
                    el = ln[76:78].strip().upper()
                    if an.startswith("H") or el == "H":
                        h_count_fallback += 1
                    toks = ln.split()
                    if len(toks) >= 9:
                        try:
                            float(toks[8]); has_q = True
                        except Exception:
                            pass
                    if "A" in ln[-6:]:
                        ad4_seen = True
            h_count = h_count or h_count_fallback
            has_partial = has_q
            ad4_ok = ad4_seen
        except Exception:
            notes.append("pdbqt_parse_warn")

    return {
        "ok": True,
        "h_count": int(h_count),
        "formal_charge": int(formal_charge),
        "has_partial_charges": bool(has_partial),
        "aromatic_atoms": int(aromatic_atoms),
        "ad4_types_ok": ad4_ok,
        "notes": notes,
        "tag": tag,
        "context": context,
    }
# === AROMATICITY AUDIT & RESCUE ===

def _count_aromatic_atoms_in_mol2(mol2_path: Path) -> int:
    """
    Count aromatic atoms from MOL2 without full sanitize to avoid RDKit precondition violations
    on under-hydrogenated inputs.
    """
    try:
        m = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
        if m is None:
            return -1
        # Guard: ensure property cache and ring info exist, then set aromaticity
        try:
            m.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        try:
            Chem.GetSymmSSSR(m)
        except Exception:
            pass
        try:
            Chem.SetAromaticity(m, Chem.AromaticityModel.AROMATICITY_RDKIT)
        except Exception:
            pass
        return sum(int(a.GetIsAromatic()) for a in m.GetAtoms())
    except Exception:
        return -1



def _count_aromatic_ad_types_in_pdbqt(pdbqt_path: Path) -> int:
    n = 0
    try:
        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                    continue
                t = ln.split()[-1].upper()
                if t in ("A", "NA"):
                    n += 1
    except Exception:
        return -1
    return n

def _pdbqt_from_mol2_via_obabel(mol2_file: Path, out_pdbqt: Path, obabel_exe_short: str) -> bool:
    """
    OBabel path with explicit polar hydrogens and Gasteiger partial charges.
    """
    try:
        tmp = out_pdbqt.with_suffix(".obabel_tmp.pdbqt")
        cmd = [
            obabel_exe_short,
            "-imol2", str(mol2_file),
            "-opdbqt",
            "-O", str(tmp),
            "-h",
            "--partialcharge", "gasteiger",
        ]
        ok = _run_obabel(cmd, timeout_sec=600)
        if ok and tmp.exists():
            try:
                if out_pdbqt.exists():
                    out_pdbqt.unlink()
            except Exception:
                pass
            tmp.replace(out_pdbqt)
            return True
    except Exception as e:
        logging.error("[obabel] MOL2->PDBQT failed: %s", e)
    return False

def _pdbqt_has_H(pdbqt_path: Path) -> bool:
    try:
        with open(pdbqt_path, "r", errors="ignore") as f:
            for line in f:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[76:78].strip() == "H":
                    return True
    except Exception:
        pass
    return False



def _count_explicit_H_in_mol2(path: Path) -> int:
    try:
        from rdkit import Chem
        m = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
        if not m:
            return -1
        return sum(1 for a in m.GetAtoms() if a.GetSymbol() == "H")
    except Exception:
        return -1

def _prepare_one(
        mgltools_python_short: str,
        prepare_script_short: str,
        mol2_file: Path,
        pdbqt_path: Path,
        obabel_exe_short: Optional[str] = None,
        *,
        status_log_dir: Path
) -> Tuple[str, str]:
    # use the H-enriched input

    lig_id = mol2_file.stem
    writer_final = "mgltools"
    rescue_used = False
    torsion_rule_label = "adt_default"
    logging.info("[debug] _prepare_one lig=%s mol2=%s out=%s obabel=%s",
                 lig_id, mol2_file.name, pdbqt_path.name, bool(obabel_exe_short))

    try:
        if obabel_exe_short:
            try:
                _pre_arom = _count_aromatic_atoms_in_mol2(mol2_file)
                _m_pre = Chem.MolFromMol2File(str(mol2_file), sanitize=False, removeHs=False)
                _preH = sum(1 for a in (_m_pre.GetAtoms() if _m_pre else []) if a.GetSymbol() == "H")
                logging.info("[ligprep] re_arom pre  arom=%d H=%d file=%s", _pre_arom, _preH, mol2_file.name)
            except Exception:
                logging.info("[ligprep] re_arom pre  arom=? H=? file=%s", mol2_file.name)

            pre_H = _count_explicit_H_in_mol2(mol2_file)
            _re_aromatize_mol2_in_place(mol2_file, obabel_exe_short)
            post_H = _count_explicit_H_in_mol2(mol2_file)
            logging.info(f"[re-arom] lig={lig_id} H_pre={pre_H} H_post={post_H} file={mol2_file.name}")

            if pre_H > 0 and post_H == 0:
                logging.warning(f"[re-arom] LOST all explicit H! Will re-add via obabel -h.")
                add_hydrogens_mol2(mol2_in, mol2_in, obabel_exe_short)  # in-place re-add
                post2_H = _count_explicit_H_in_mol2(mol2_in)
                logging.info(f"[re-arom] after re-add: H={post2_H}")

            try:
                _post_arom = _count_aromatic_atoms_in_mol2(mol2_file)
                _m_post = Chem.MolFromMol2File(str(mol2_file), sanitize=False, removeHs=False)
                _postH = sum(1 for a in (_m_post.GetAtoms() if _m_post else []) if a.GetSymbol() == "H")
                logging.info("[ligprep] re_arom post arom=%d H=%d file=%s", _post_arom, _postH, mol2_file.name)
            except Exception:
                logging.info("[ligprep] re_arom post arom=? H=? file=%s", mol2_file.name)
    except Exception as e:
        logging.info("[ligprep] re_arom skipped for %s: %s", mol2_file.name, e)

    try:
        if pdbqt_path.exists():
            pdbqt_path.unlink()
    except Exception:
        pass
    # Ensure explicit hydrogens: mol2_in -> mol2_h
    mol2_in = Path(mol2_file)
    mol2_h = mol2_in.with_name(mol2_in.stem + ".withH.mol2")

    preH0 = _count_explicit_H_in_mol2(mol2_in)
    okH, hstderr = add_hydrogens_mol2(mol2_in, mol2_h, obabel_exe_short)
    postH0 = _count_explicit_H_in_mol2(mol2_h) if mol2_h.exists() else -1

    logging.info("[ligprep] addHs stage=primary in=%s ok=%s preH=%s postH=%s",
                 mol2_in.name, okH, preH0, postH0)
    if (hstderr or "").strip():
        logging.warning("[ligprep] addHs stderr (primary) %s", hstderr.splitlines()[-1][:200])

    # Treat 'success but zero H' as failure so we fall back
    if not okH or postH0 <= 0 or (preH0 >= 0 and postH0 <= preH0):
        logging.warning("[ligprep] obabel_h_charge: H not added (pre=%s post=%s); forcing RDKit AddHs fallback: %s",
                        preH0, postH0, mol2_in.name)
        mol2_for_mgl = mol2_in
    else:
        mol2_for_mgl = mol2_h

    # If still H-free, attempt an RDKit-based AddHs roundtrip
    if mol2_for_mgl is mol2_in:
        logging.warning("[ligprep] OBabel -h failed (or no H added); attempting RDKit AddHs fallback: %s", mol2_in.name)
        try:
            _m = Chem.MolFromMol2File(str(mol2_in), sanitize=True, removeHs=False)
            if _m is not None:
                _mh = Chem.AddHs(_m)
                _tmp_sdf = mol2_in.with_suffix(".rdkH.sdf")
                _tmp_mol2 = mol2_in.with_suffix(".rdkH.mol2")
                with Chem.SDWriter(str(_tmp_sdf)) as w:
                    w.write(_mh)
                ok2 = _run_obabel([obabel_exe_short, "-isdf", str(_tmp_sdf), "-omol2", "-O", str(_tmp_mol2)],
                                  timeout_sec=600)
                if ok2 and _tmp_mol2.exists():
                    mol2_for_mgl = _tmp_mol2
                    logging.info("[ligprep] RDKit AddHs fallback produced MOL2: %s", _tmp_mol2.name)
        except Exception as e:
            logging.warning("[ligprep] RDKit AddHs fallback failed: %s", e)

    # Final ADT input sanity (now also logs H count explicitly)
    _m_chk = Chem.MolFromMol2File(str(mol2_for_mgl), sanitize=False, removeHs=False)
    curH = sum(1 for a in (_m_chk.GetAtoms() if _m_chk else []) if a.GetSymbol() == "H")
    if curH > 0:
        logging.info("[ligprep] ADT input Hs=%d path=%s", curH, mol2_for_mgl)
    else:
        logging.warning("[ligprep] ADT will add hydrogens internally (input has 0 explicit H): %s", mol2_for_mgl)

    try:
        _chk = Chem.MolFromMol2File(str(mol2_for_mgl), sanitize=False, removeHs=False)
        logging.info("[bulkSDF] post-AddHs ligand=%s atoms=%d explicit_H=%s",
                     lig_id,
                     (_chk.GetNumAtoms() if _chk else -1),
                     (any(a.GetAtomicNum() == 1 for a in _chk.GetAtoms()) if _chk else False))



    except Exception:
        pass

    def _dbg_atom_stats(mol2_path: Path, tag: str):
        """
        RDKit-safe debug of atoms/rings/charges on raw MOL2 without raising
        precondition violations. Adds greppable tag [rdkit-guarded].
        """
        try:
            from rdkit import Chem
            m = None
            if mol2_path.is_file():
                m = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
            if m:
                # Guard: property cache + ring info population
                try:
                    m.UpdatePropertyCache(strict=False)  # [rdkit-guarded]
                except Exception:
                    pass
                try:
                    Chem.GetSymmSSSR(m)  # [rdkit-guarded]
                except Exception:
                    pass

                n_atoms = m.GetNumAtoms()
                try:
                    ri = m.GetRingInfo()
                    n_rings = int(ri.NumRings()) if ri is not None else 0
                except Exception:
                    n_rings = 0

                try:
                    charge_sum = int(sum(a.GetFormalCharge() for a in m.GetAtoms()))
                except Exception:
                    charge_sum = 0

                try:
                    has_explicit_H = any(a.GetSymbol() == "H" for a in m.GetAtoms())
                except Exception:
                    has_explicit_H = False

                print(f"[rdkit-guarded] {tag}: atoms={n_atoms} rings={n_rings} "
                      f"charge_sum={charge_sum} has_explicit_H={has_explicit_H}")
            else:
                print(f"[ligprep] {tag}: RDKit failed to parse {mol2_path.name}")
        except Exception as e:
            print(f"[ligprep] {tag}: dbg_error {e}")

    _dbg_atom_stats(mol2_for_mgl, "pre-ADT/prepare_ligand4")
    # Compact pre-ADT banner (greppable)
    try:
        _m_pre = Chem.MolFromMol2File(str(mol2_for_mgl), sanitize=False, removeHs=False)
        _Hc = sum(1 for a in (_m_pre.GetAtoms() if _m_pre else []) if a.GetSymbol() == "H")
        # Optional: count aromatic atoms if you have a helper; else report -1
        try:
            _arom_atoms = _count_aromatic_atoms_in_mol2(mol2_for_mgl)
        except Exception:
            _arom_atoms = -1
        print(
            f"[pre-adt] lig={lig_id} in_mol2={mol2_for_mgl.name} "
            f"atoms={_m_pre.GetNumAtoms() if _m_pre else -1} H={_Hc} arom_atoms={_arom_atoms}"
        )
    except Exception:
        pass

    # aromatic baseline from source MOL2 (keep original for comparison)
    src_arom = _count_aromatic_atoms_in_mol2(mol2_file)
    # Final pre-ADT hydrogen sanity; add via OBabel if input is H-free
    H_in = _count_explicit_H_in_mol2(mol2_for_mgl)
    logging.info(f"[ADT] input={mol2_for_mgl.name} H={H_in}")

    if H_in == 0:
        logging.warning("[ADT] No explicit H in input MOL2; forcing obabel -h first.")
        add_hydrogens_mol2(mol2_for_mgl, mol2_for_mgl, obabel_exe_short)

    # use the H-enriched file for MGLTools (cwd is mol2_file.parent)
    mol2_short = mol2_for_mgl.name
    pdbqt_short = get_short_path_name(str(pdbqt_path.resolve()))

    # Environment toggle: KEEP_NONPOLAR_H=1 keeps nphs; default True for extracted runs
    _keep_nphs = str(os.environ.get("KEEP_NONPOLAR_H", "1")).lower() not in {"0", "false", "no"}
    _u_flag = ["-U", "lps"] if _keep_nphs else ["-U", "nphs_lps"]
    # Guard: ensure the ADT input MOL2 actually has explicit H
    H_in = _count_explicit_H_in_mol2(mol2_for_mgl)
    logging.info(f"[ADT] input={Path(mol2_for_mgl).name} H={H_in}")
    if H_in == 0 and obabel_exe_short:
        logging.warning("[ADT] No explicit H in input MOL2; forcing obabel -h first.")
        # in-place enrich (safe because we immediately read it)
        add_hydrogens_mol2(mol2_for_mgl, mol2_for_mgl, obabel_exe_short)

    cmd = [
        mgltools_python_short, prepare_script_short,
        "-l", get_short_path_name(str(mol2_for_mgl.resolve())),
        "-o", pdbqt_short,  
        *_u_flag,
        "-A", "hydrogens",
    ]

    print(f"[ligprep] Calling prepare_ligand4 on {mol2_file.name} -> {Path(pdbqt_path).name}")
    logging.info("[paths] pdbqt_out=%s (from mol2=%s)", str(pdbqt_path.resolve()), str(mol2_file.resolve()))

    try:
        _ = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            cwd=str(mol2_file.parent),
            timeout=600
        )
        try:
            size = pdbqt_path.stat().st_size if pdbqt_path.exists() else 0
            logging.info("[bulkSDF] write.pdbqt.done ligand=%s rc=0 size=%d via=mgltools", lig_id, size)
        except Exception:
            pass

        try:
            ok, n_atoms, reason = quick_pdbqt_validate(pdbqt_path)
            print(f"[ligprep] post-ADT PDBQT atoms={n_atoms} ok={ok} reason={reason}")
            print(f"[post-adt] lig={lig_id} out_atoms={n_atoms}")
            print(f"[fallback-obabel] lig={lig_id} out_atoms={n_atoms} stderr_last=''")

            # Initialize candidate list for deterministic writer selection
            candidates: list[tuple[str, Path, dict, int]] = []  # (label, path, metrics, arom_count)
            try:
                m0 = _parse_pdbqt_metrics(pdbqt_path)
                a0 = _count_aromatic_ad_types_in_pdbqt(pdbqt_path)
                candidates.append(("mgltools_primary", pdbqt_path, m0, a0))
            except Exception:
                pass

        except Exception:
            ok, n_atoms, reason = False, 0, "validate_exception"

        # Record atom counts for fallback comparisons (needed before the ADT→OBabel decision)
        _in_atoms = -1
        try:
            _in_atoms = _m_chk.GetNumAtoms() if _m_chk else -1  # or a MOL2 atom counter if _m_chk isn't set
        except Exception:
            pass
        _out_atoms = n_atoms
        # If ADT wrote noticeably fewer atoms than the MOL2 input, rescue via OBabel MOL2->PDBQT
        try:
            def _count_mol2_atoms(mol2_path: Path) -> int:
                n, in_atoms = 0, False
                with open(mol2_path, "r", encoding="utf-8", errors="ignore") as fh:
                    for ln in fh:
                        if ln.startswith("@<TRIPOS>ATOM"): in_atoms = True; continue
                        if ln.startswith("@<TRIPOS>BOND"): break
                        if in_atoms: n += 1
                return n

            _in_atoms = _m_chk.GetNumAtoms() if _m_chk else _count_mol2_atoms(mol2_for_mgl)
        except Exception:
            _in_atoms = _count_mol2_atoms(mol2_for_mgl)
        print(
            f"[ligprep] writer={writer_final} lig={lig_id} in_atoms={_in_atoms} out_atoms={n_atoms} ok={ok} reason={reason}")

        if _out_atoms > 0 and _in_atoms > 0 and _out_atoms < int(0.9 * _in_atoms):
            # One light retry: keep hydrogens; do not strip lone pairs
            try:
                alt_pdbqt = pdbqt_path.with_suffix(".loss_retry.pdbqt")
                alt_cmd = [
                    mgltools_python_short, prepare_script_short,
                    "-l", mol2_short,
                    "-o", get_short_path_name(str(alt_pdbqt.resolve())),
                    "-U", "lps",
                    "-A", "hydrogens",
                ]
                _ = subprocess.run(
                    alt_cmd, check=True, capture_output=True, text=True,
                    cwd=str(mol2_file.parent), timeout=600
                )
                ok_alt, n_alt, _ = quick_pdbqt_validate(alt_pdbqt)
                if ok_alt and n_alt >= int(0.9 * _in_atoms):
                    # Promote retry result to main path
                    try:
                        if pdbqt_path.exists():
                            pdbqt_path.unlink()
                    except Exception:
                        pass
                    alt_pdbqt.replace(pdbqt_path)
                    writer_final = "mgltools"
                    print(f"[post-adt] lig={lig_id} out_atoms={n_alt} (loss-retry)")
                    _out_atoms = n_alt
                else:
                    try:
                        alt_pdbqt.unlink(missing_ok=True)
                    except Exception:
                        pass
            except Exception:
                pass

            logging.warning("[ligprep] ADT wrote fewer atoms (%d -> %d); invoking OBabel MOL2->PDBQT fallback",
                            _in_atoms, _out_atoms)
            if obabel_exe_short:
                ok_ob = _pdbqt_from_mol2_via_obabel(Path(proto_mol2), pdbqt_path, obabel_exe_short)
                logging.info("[ligprep] ADT altpath (via OBabel) OK=%s", ok_ob)
                writer_final = "obabel"
                rescue_used = True

                # refresh validation counters for downstream decisions
                try:
                    # refresh validation counters for downstream decisions (post-OBabel rescue)
                    ok, n_atoms, reason = quick_pdbqt_validate(pdbqt_path)
                    print(f"[ligprep] post-OBabel rescue PDBQT atoms={n_atoms} ok={ok} reason={reason}")
                    # Force AD4 typing/normalization after OBabel rescue
                    print(f"[retype] lig={lig_id} reason=post-obabel path={pdbqt_path.name}")
                    ok_norm, torsion_rule_label = _adt_retype_and_normalize(
                        mgltools_python_short, prepare_script_short, mol2_for_mgl, pdbqt_path,
                        torsion_rule_label="adt_normalized_after_obabel"
                    )
                    # Probe invariants again to surface AD4 typing status
                    _, __metrics, inv_probe_fail = validate_pdbqt_invariants(pdbqt_path)
                    try:
                        # If your validator returns (ok, metrics, fails), grab metrics from the second element
                        inv_ok_probe, metrics_probe, _fails_dummy = validate_pdbqt_invariants(pdbqt_path)
                        print(f"[retype.result] lig={lig_id} ad4={metrics_probe.get('has_ad4_types')} "
                              f"torsdof={metrics_probe.get('torsdof')}")
                        # Capture OBabel+ADT-normalized as another candidate
                        try:
                            m1 = _parse_pdbqt_metrics(pdbqt_path)
                            a1 = _count_aromatic_ad_types_in_pdbqt(pdbqt_path)
                            candidates.append(("obabel+adt_norm", pdbqt_path, m1, a1))
                        except Exception:
                            pass

                    except Exception:
                        print(f"[retype.result] lig={lig_id} ad4=? torsdof=?")
                    if ok_norm:
                        writer_final = "mgltools"

                    _in_atoms = -1
                    try:
                        _in_atoms = _m_chk.GetNumAtoms() if _m_chk else -1
                    except Exception:
                        pass
                    _out_atoms = n_atoms
                except Exception:
                    pass

        adt_arom = _count_aromatic_ad_types_in_pdbqt(pdbqt_path)

        if src_arom >= 0 and adt_arom >= 0 and adt_arom < src_arom:
            logging.warning(
                f"[arom-mismatch] {mol2_file.name}: MOL2_arom={src_arom} > PDBQT_arom={adt_arom} (trying rescue)")

            # 1) Try prepare_ligand4 again WITHOUT removing lone pairs: "-U nphs"
            alt_pdbqt = pdbqt_path.with_suffix(".alt.pdbqt")
            alt_cmd = [
                mgltools_python_short, prepare_script_short,
                "-l", mol2_short,
                "-o", get_short_path_name(str(alt_pdbqt.resolve())),
                "-U", "nphs",
                "-A", "hydrogens",  # <-- fixed
            ]
            try:
                _ = subprocess.run(
                    alt_cmd,
                    check=True,
                    capture_output=True,
                    text=True,
                    cwd=str(mol2_file.parent),
                    timeout=600
                )
                alt_arom = _count_aromatic_ad_types_in_pdbqt(alt_pdbqt)
            except Exception:
                alt_arom = -1

            # 2) Try OBabel -> PDBQT as a final rescue
            ob_arom = -1
            if obabel_exe_short:
                ob_tmp = pdbqt_path.with_suffix(".ob.pdbqt")
                if _pdbqt_from_mol2_via_obabel(mol2_file, ob_tmp, obabel_exe_short):
                    ob_arom = _count_aromatic_ad_types_in_pdbqt(ob_tmp)

            # Choose the best among {original, alt, obabel} by aromatic count
            candidates = [(adt_arom, pdbqt_path)]
            if alt_arom >= 0:
                candidates.append((alt_arom, alt_pdbqt))
            if ob_arom >= 0:
                candidates.append((ob_arom, ob_tmp))
            best_arom, best_path = max(candidates, key=lambda t: t[0])

            # If best is not the current path, replace
            if best_path != pdbqt_path and best_path.exists():
                try:
                    pdbqt_path.unlink(missing_ok=True)
                except Exception:
                    pass
                best_path.replace(pdbqt_path)

            # cleanup temps
            for tmp in [alt_pdbqt, pdbqt_path.with_suffix(".ob.pdbqt")]:
                try:
                    if tmp.exists():
                        tmp.unlink()
                except Exception:
                    pass

        # Ensure the final PDBQT has explicit H; if not, re-write via OBabel with -h
        if pdbqt_path.exists() and obabel_exe_short and not _pdbqt_has_H(pdbqt_path):
            logging.warning(f"[post] {pdbqt_path.name} has no H; re-writing via OBabel with -h")
            _pdbqt_from_mol2_via_obabel(mol2_file, pdbqt_path, obabel_exe_short)
            logging.info(f"[post] re-write complete; has_H={_pdbqt_has_H(pdbqt_path)}")

        try:
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
                lines = fh.readlines()
            new_lines, fixes, q_reason = assert_no_helium_in_pdbqt(lines, lig_id)
            if q_reason:
                quarantine = pdbqt_path.parent / QUARANTINE_DIRNAME
                quarantine.mkdir(exist_ok=True)
                qpath = quarantine / pdbqt_path.name
                logging.info("[ligprep] quarantine reason=%s ligand=%s", q_reason, lig_id)
                try:
                    if pdbqt_path.exists():
                        pdbqt_path.replace(qpath)
                except Exception:
                    pass
                _append_prep_status(
                    status_log_path=status_log_dir / "ligand_prep_status.tsv",
                    ligand_name=lig_id,
                    status="quarantined",
                    reason="adt_helium_inconsistent",
                    relpath=str(qpath.name),
                    stage="validate_pdbqt",
                    failure_code="adt_helium_inconsistent",
                    failure_detail=q_reason,
                    fixes_count=0,
                    rules_ver=rules_version(),
                )
                logging.warning(
                    f"[elem-validate] file={pdbqt_path.name} stage=postwrite quarantine ligand={lig_id} code=adt_helium_inconsistent")
                return (mol2_file.name, "postcheck_fail")
            if fixes > 0:
                with open(pdbqt_path, "w", encoding="utf-8") as out:
                    out.writelines(new_lines)
                _append_prep_status(
                    status_log_path=status_log_dir / "ligand_prep_status.tsv",
                    ligand_name=lig_id,
                    status="ok",
                    reason="",
                    relpath=str(pdbqt_path.name),
                    stage="write_pdbqt",
                    failure_code="",
                    failure_detail="",
                    fixes_count=fixes,
                    rules_ver=rules_version(),
                )
                logging.info(
                    f"[elem-validate] file={pdbqt_path.name} stage=postwrite ok ligand={lig_id} fixes={fixes}")
            else:
                _append_prep_status(
                    status_log_path=status_log_dir / "ligand_prep_status.tsv",
                    ligand_name=lig_id,
                    status="ok",
                    reason="",
                    relpath=str(pdbqt_path.name),
                    stage="write_pdbqt",
                    failure_code="",
                    failure_detail="",
                    fixes_count=0,
                    rules_ver=rules_version(),
                )
        except Exception as e:
            logging.warning(f"[elem-validate] failed to post-check PDBQT for {lig_id}: {e}")

        # --- lightweight post-write PDBQT validity gate (primary) ---
        okV, n_atoms, whyV = quick_pdbqt_validate(pdbqt_path)
        logging.info("[ligprep] validate_pdbqt result=%s atoms=%d ligand=%s",
                     "ok" if okV else "fail", n_atoms, mol2_file.stem)
        if not okV:
            # quarantine + TSV with failure_code
            quarantine = pdbqt_path.parent / QUARANTINE_DIRNAME
            quarantine.mkdir(exist_ok=True)
            qpath = quarantine / pdbqt_path.name
            logging.info("[ligprep] quarantine reason=%s ligand=%s", whyV, mol2_file.stem)
            try:
                if pdbqt_path.exists():
                    pdbqt_path.replace(qpath)
            except Exception:
                pass
            _append_prep_status(
                status_log_path=status_log_dir / "ligand_prep_status.tsv",
                ligand_name=mol2_file.stem,
                status="quarantined",
                reason="invalid_pdbqt_postwrite",
                relpath=str(qpath.name if qpath.exists() else ""),
                stage="validate_pdbqt",
                failure_code="invalid_pdbqt_postwrite",
                failure_detail=whyV,
                fixes_count=0,
                rules_ver=rules_version(),
            )
            return (mol2_file.name, "postcheck_fail")

        # === Invariant validation + normalization ===
        inv_ok, metrics, inv_fail = validate_pdbqt_invariants(pdbqt_path)
        print(f"[invariants.summary] lig={lig_id} typed={metrics.get('typed')}/{metrics.get('n_atoms')} "
              f"charges_ok={metrics.get('has_charges')} ad4_ok={metrics.get('has_ad4_types')} "
              f"torsdof={metrics.get('torsdof')} writer={writer_final}")
        # Deterministic winner selection across accumulated candidates
        try:
            cur_m = _parse_pdbqt_metrics(pdbqt_path)
            cur_a = _count_aromatic_ad_types_in_pdbqt(pdbqt_path)
            have = {(p.resolve() if isinstance(p, Path) else p) for _, p, _, _ in candidates}
            if pdbqt_path.resolve() not in have:
                candidates.append((writer_final, pdbqt_path, cur_m, cur_a))

            def _score(entry):
                _label, _path, _m, _a = entry
                typed = int(_m.get("typed", 0));
                nat = int(_m.get("n_atoms", 0))
                typed_ratio = (typed / nat) if nat > 0 else 0.0
                return (typed_ratio, nat, _a)

            if candidates:
                best = max(candidates, key=_score)
                chosen_lbl, chosen_path, chosen_m, chosen_a = best
                if chosen_path.resolve() != pdbqt_path.resolve() and chosen_path.exists():
                    try:
                        tmp_best = pdbqt_path.with_suffix(".winner.pdbqt")
                        chosen_path.replace(tmp_best)
                        if pdbqt_path.exists():
                            pdbqt_path.unlink(missing_ok=True)
                        tmp_best.replace(pdbqt_path)
                    except Exception:
                        pass
                logging.info("[writer-select] chosen=%s candidates=%s",
                             chosen_lbl,
                             ",".join(f"{lbl}:{m.get('typed', 0)}/{m.get('n_atoms', 0)}@arom{a}"
                                      for (lbl, _, m, a) in candidates[:5]))
        except Exception:
            pass

        if not inv_ok:
            logging.warning("[invariants] %s failed: %s", pdbqt_path.name, ",".join(inv_fail))

            # If OBabel wrote it or AD4 types are missing, enforce ADT typing + torsion normalization
            if (writer_final == "obabel") or ("missing_AD4_types" in inv_fail):
                ok_norm, torsion_rule_label = _adt_retype_and_normalize(
                    mgltools_python_short, prepare_script_short, mol2_for_mgl, pdbqt_path,
                    torsion_rule_label="adt_normalized"
                )
                if ok_norm:
                    writer_final = "mgltools"
                    inv_ok, metrics, inv_fail = validate_pdbqt_invariants(pdbqt_path)

        if not inv_ok:
            # Try Meeko from PATH, then normalize back to ADT-style torsions
            meeko_ok = _reprep_via_meeko_path(mol2_for_mgl, pdbqt_path)
            if meeko_ok:
                writer_final = "meeko"
                rescue_used = True
                ok_norm, torsion_rule_label = _adt_retype_and_normalize(
                    mgltools_python_short, prepare_script_short, mol2_for_mgl, pdbqt_path,
                    torsion_rule_label="adt_normalized_after_meeko"
                )
                if ok_norm:
                    writer_final = "mgltools"
                inv_ok, metrics, inv_fail = validate_pdbqt_invariants(pdbqt_path)

        if not inv_ok:
            # Fail-closed: quarantine with precise failure codes
            quarantine = pdbqt_path.parent / QUARANTINE_DIRNAME
            quarantine.mkdir(exist_ok=True)
            qpath = quarantine / pdbqt_path.name
            try:
                if pdbqt_path.exists():
                    pdbqt_path.replace(qpath)
            except Exception:
                pass
            _append_prep_status(
                status_log_path=status_log_dir / "ligand_prep_status.tsv",
                ligand_name=lig_id,
                status="quarantined",
                reason="FAIL_CLOSED",
                relpath=str(qpath.name),
                stage="validate_invariants",
                failure_code="FAIL_CLOSED",
                failure_detail=",".join(inv_fail),
                fixes_count=0,
                rules_ver=rules_version(),
                writer_final=writer_final,
                rescue_used=str(rescue_used),
                torsion_root_rule=torsion_rule_label,
                polarH=str(bool(metrics.get("has_polar_H", False))),
                ad4_types_ok=str(bool(metrics.get("has_ad4_types", False))),
                charges_ok=str(bool(metrics.get("has_charges", False))),
                torsdof=str(metrics.get("torsdof", -1)),
            )
            return (mol2_file.name, "postcheck_fail")

        # If invariants OK, append enriched TSV row now
        try:
            _append_prep_status(
                status_log_path=status_log_dir / "ligand_prep_status.tsv",
                ligand_name=lig_id,
                status="ok",
                reason="",
                relpath=str(pdbqt_path.name),
                stage="write_pdbqt",
                failure_code="",
                failure_detail="",
                fixes_count=0,
                rules_ver=rules_version(),
                writer_final=writer_final,
                rescue_used=str(rescue_used),
                torsion_root_rule=torsion_rule_label,
                polarH=str(bool(metrics.get("has_polar_H", False))),
                ad4_types_ok=str(bool(metrics.get("has_ad4_types", False))),
                charges_ok=str(bool(metrics.get("has_charges", False))),
                torsdof=str(metrics.get("torsdof", -1)),
            )
        except Exception as e:
            logging.warning("[invariants] TSV append failed for %s: %s", lig_id, e)


        # Validation / quarantine
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
        try:
            print(f"[retry-timeout] lig={lig_id} threads={os.environ.get('OBABEL_THREADS', '')} timeout=1200")
            _ = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                cwd=str(mol2_file.parent),
                timeout=1200
            )
            okV, n_atoms, _ = quick_pdbqt_validate(pdbqt_path)
            if okV and n_atoms > 0:
                inv_ok, metrics, _ = validate_pdbqt_invariants(pdbqt_path)
                print(f"[invariants.summary] lig={lig_id} typed={metrics.get('typed')}/{metrics.get('n_atoms')} "
                      f"charges_ok={metrics.get('has_charges')} ad4_ok={metrics.get('has_ad4_types')} "
                      f"torsdof={metrics.get('torsdof')} writer={writer_final}")
                return (mol2_file.name, "ok" if inv_ok else "postcheck_fail")
        except Exception:
            pass
        return (mol2_file.name, "timeout_after_retry")

    except subprocess.CalledProcessError as e:
        return (mol2_file.name, f"prepare_fail:{e.returncode}")


def load_mol2_lenient(path, logger):
    mol = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
    if mol is None:
        _log_malformed(Path(path), "rdkit_read_fail")
        return None
    try:
        # --- DEBUG: pre-sanitize fingerprint ---
        def _dbg_atom_stats(mol2_path: Path, tag: str):
            try:
                from rdkit import Chem
                m = None
                if mol2_path.is_file():
                    m = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
                if m:
                    try:
                        m.UpdatePropertyCache(strict=False)
                    except Exception:
                        pass
                    # Safely derive ring info without spamming errors on unsanitized mols
                    try:
                        _ = Chem.GetSymmSSSR(m)  # populate ring info if possible
                        ri = m.GetRingInfo()
                        n_rings = (ri.NumRings() if ri is not None else 0)
                    except Exception:
                        n_rings = -1
                    try:
                        n_atoms = m.GetNumAtoms()
                    except Exception:
                        n_atoms = -1
                    try:
                        charge = int(sum(a.GetFormalCharge() for a in (m.GetAtoms() if m else [])))
                    except Exception:
                        charge = 0
                    try:
                        has_H = any(a.GetSymbol() == "H" for a in (m.GetAtoms() if m else []))
                    except Exception:
                        has_H = False
                    print(
                        f"[ligprep] {tag}: atoms={n_atoms} rings={n_rings} charge_sum={charge} has_explicit_H={has_H}")
                else:
                    print(f"[ligprep] {tag}: RDKit failed to parse {mol2_path.name}")
            except Exception as e:
                # compress noisy stack traces into a single line
                print(f"[ligprep] {tag}: dbg_skip ({e.__class__.__name__})")

                # Make valence/ring queries safe
                try:
                    m.UpdatePropertyCache(strict=False)
                except Exception:
                    pass

                atoms = list(m.GetAtoms())
                elems = Counter(a.GetSymbol() for a in atoms)

                # Oxygen valence histogram (guarded)
                o_vals = Counter()
                for a in atoms:
                    if a.GetSymbol() == "O":
                        try:
                            v = int(a.GetExplicitValence())
                        except Exception:
                            try:
                                v = int(a.GetTotalValence())
                            except Exception:
                                v = -1
                        o_vals[v] += 1

                # Ring count (guarded)
                try:
                    ri = m.GetRingInfo()
                    ring_count = ri.NumRings() if ri is not None else 0
                except Exception as e:
                    ring_count = -1
                    print(f"[ligprep][{tag}] dbg: GetRingInfo failed for {src_name}: {e}")

                chg_sum = 0
                try:
                    chg_sum = sum(int(a.GetFormalCharge()) for a in atoms)
                except Exception:
                    pass

                has_exp_h = any(a.GetNumExplicitHs() > 0 for a in atoms)

                print(
                    f"[ligprep][{tag}] file={src_name} atoms={len(atoms)} elems={dict(elems)} "
                    f"O_val_hist={dict(o_vals)} rings={ring_count} charge_sum={chg_sum} explicitH={int(has_exp_h)}"
                )

            except Exception as e:
                print(f"[ligprep][{tag}] DEBUG_STATS_ERROR: {e}")

        _dbg_atom_stats(mol, "pre-sanitize")
        Chem.SanitizeMol(mol)
        _dbg_atom_stats(mol, "post-sanitize")
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
                # On sanitize/valence error: dump quick forensics
                try:
                    smi = Chem.MolToSmiles(mol, isomericSmiles=True) if mol is not None else "None"
                except Exception:
                    smi = "MolToSmiles_failed"
                try:
                    from collections import Counter
                    elem_hist = dict(Counter(a.GetSymbol() for a in mol.GetAtoms())) if mol is not None else {}
                    fcharge = sum(int(a.GetFormalCharge()) for a in (mol.GetAtoms() if mol else []))
                except Exception:
                    elem_hist, fcharge = {}, 0
                logging.warning("[ligprep] sanitize_failed elem=%s formal_charge=%d smiles=%s", elem_hist, fcharge, smi)
            return None
        return std_mol


# =========================
# PDB-control ligand path (crystal-safe)
# =========================

def _resolve_obabel_exe(obabel_cfg: str) -> str:
    """Linux-safe: if cfg is a file use it; if it's a dir use dir/obabel; else pass through."""
    p = Path(obabel_cfg)
    if p.is_file():
        return str(p)
    if p.is_dir():
        cand = p / "obabel"
        return str(cand) if cand.exists() else str(p)
    return obabel_cfg


def _resolve_prepare_ligand4(mgltools_path: str, cfg: Dict[str, str]) -> Path:
    """Support PREPARE_LIGAND4 override; probe Windows-style and Linux-style locations."""
    override = cfg.get("PREPARE_LIGAND4")
    if override:
        pp = Path(override)
        if pp.exists():
            return pp
    mp = Path(mgltools_path)
    win_probe = mp / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
    lin_probe = mp / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
    print("using prepare_ligand4.py at ", lin_probe)
    if win_probe.exists():
        return win_probe
    return lin_probe


def prep_ligands_from_pdb(ligand_output_dir: Path, ligands_mol2_dir: Path, prepped_ligands_dir: Path):
    """
    Crystal-safe path to prepare ligands that were extracted from PDBs (processed_pdbs/*/ligands_raw/*.pdb).
    Hardened to mirror bulk path resilience:
      - layered protonation (RDKit AddHs -> OBabel -h -> ADT/Meeko during typing)
      - unified writer via _prepare_one(...) with AD4 normalization + invariant checks
      - crystal-safety scrubs + helium guard
      - concise one-line summary per ligand
      - test-mode subset selection via EXTRACT_ONLY/EXTRACT_TEST (and mirrored CLI)
    """
    logging.info("Starting ligand preparation from PDB files (crystal-safe, MOL2-first path)")

    config = read_config()
    mgltools_python = config.get("MGLTOOLS_PYTHON")
    mgltools_path = config.get("MGLTOOLS_PATH")
    obabel_exe_cfg = config.get("OPENBABEL_PATH")

    if not mgltools_python or not mgltools_path or not obabel_exe_cfg:
        raise RuntimeError("Missing paths in config.txt: MGLTOOLS_PYTHON, MGLTOOLS_PATH, OPENBABEL_PATH")

    obabel_exe = obabel_exe_cfg
    obabel_exe_short = get_short_path_name(obabel_exe)

    mgltools_python_short = get_short_path_name(mgltools_python)
    prepare_script = _resolve_prepare_ligand4(mgltools_path, config)
    if not prepare_script.exists():
        raise FileNotFoundError(f"prepare_ligand4.py not found at {prepare_script} (set PREPARE_LIGAND4 in config.txt)")
    prepare_script_short = get_short_path_name(str(prepare_script.resolve()))

    # --- discovery root diagnostics ---
    root = ligand_output_dir.resolve()
    print(f"[ligprep-extracted] discovery_root={root} exists={root.exists()} is_dir={root.is_dir()}")

    # --- collect EXTRACT_ONLY tokens (basename / filename / full path) ---
    raw_only = os.environ.get("EXTRACT_ONLY", "").strip()
    requested_only: Set[str] = set()
    if raw_only:
        for tok in raw_only.replace(",", " ").split():
            t = tok.strip()
            if t:
                requested_only.add(t)

    # If EXTRACT_ONLY contains any *existing* .pdb paths, use them directly
    direct_paths: List[Path] = []
    for tok in list(requested_only):
        p = Path(tok)
        if p.suffix.lower() == ".pdb" and p.exists():
            direct_paths.append(p.resolve())

    if direct_paths:
        pdb_files = sorted(set(direct_paths))
        print(f"[ligprep-extracted] direct-path mode count={len(pdb_files)} keep={','.join(p.stem for p in pdb_files)}")
        missing = set()  # nothing to match; we used the paths verbatim
    else:
        # Fall back to recursive discovery under the configured root
        pdb_files = sorted(root.rglob("*.pdb"))
        logging.info(f"Found {len(pdb_files)} PDB ligand file(s) under {root}")

        # If subset requested by name/stem, filter against discovered set
        missing: Set[str] = set()
        if requested_only:
            by_stem: Dict[str, List[Path]] = {}
            by_name: Dict[str, List[Path]] = {}
            by_abs: Dict[str, Path] = {}
            for p in pdb_files:
                by_stem.setdefault(p.stem, []).append(p)
                by_name.setdefault(p.name, []).append(p)
                by_abs[str(p.resolve())] = p

            keep: List[Path] = []
            for r in requested_only:
                if r in by_abs:
                    keep.append(by_abs[r]);
                    continue
                if r in by_name and by_name[r]:
                    keep.extend(by_name[r]);
                    continue
                stem = Path(r).stem if Path(r).suffix else r
                # ^
                if stem in by_stem and by_stem[stem]:
                    keep.extend(by_stem[stem]);
                    continue
                missing.add(r)

            pdb_files = sorted(set(keep))

        print(
            f"[ligprep-extracted] test-mode enabled count={len(pdb_files)} missing={len(missing)} keep={','.join(Path(p).stem for p in pdb_files)}")
        if missing:
            print(f"[ligprep-extracted] warn: requested_not_found={','.join(sorted(missing))}")

    # Verbose toggle like bulk
    if str(os.environ.get("EXTRACT_TEST", "0")).lower() not in {"0", "false", "no"}:
        os.environ["LIGPREP_DEBUG"] = "1"

    prepped_ligands_dir.mkdir(parents=True, exist_ok=True)
    status_log = prepped_ligands_dir / "ligand_prep_status.tsv"

    # Intermediates for easier debugging + symlink back to source dir (per-protein; no cfg)
    try:
        inter_dirname = "intermediates"
        ref_dirname = "reference"
        quar_dirname = "quarantine"
        inter_dir = prepped_ligands_dir / inter_dirname
        ref_dir = prepped_ligands_dir / ref_dirname
        quar_dir = prepped_ligands_dir / quar_dirname
        inter_dir.mkdir(parents=True, exist_ok=True)
        ref_dir.mkdir(parents=True, exist_ok=True)
        quar_dir.mkdir(parents=True, exist_ok=True)

        link = prepped_ligands_dir / "intermediates_src"
        if link.exists() or link.is_symlink():
            try:
                if link.is_dir() and not link.is_symlink():
                    shutil.rmtree(link)
                else:
                    link.unlink()
            except Exception:
                pass
        link.symlink_to(ligand_output_dir.resolve(), target_is_directory=True)
        logging.info("[intermediates] symlinked source -> intermediates_src")
    except Exception as e:
        logging.debug("intermediates link skipped: %s", e)

    for pdb_file in pdb_files:
        logging.info(f"Processing: {pdb_file.name}")

        lig_id = pdb_file.stem  # e.g., RXT_A1204
        tmp_pdb = str(pdb_file)  # updated to sanitized below
        residue_name = lig_id.split("_")[0].upper()

        # Skip standard amino acids and known crystallization additives
        if residue_name in STANDARD_AMINO_ACIDS:
            logging.info(f"Skipping standard amino acid residue: {pdb_file.name}")
            continue
        if residue_name in EXCLUDE_CRYSTAL_ADDITIVES:
            _log_malformed(pdb_file, f"excluded_crystal_additive:{residue_name}", log_dir=prepped_ligands_dir)
            logging.info(f"Skipping crystallization additive: {pdb_file.name}")
            continue

        # Quick size gate
        with open(pdb_file, "r", encoding="utf-8", errors="ignore") as f:
            atom_lines = [line for line in f if line.startswith(("HETATM", "ATOM"))]
        if len(atom_lines) < MIN_ATOMS_FOR_DOCKING:
            _log_malformed(pdb_file, f"tiny_ligand_fewer_than_{MIN_ATOMS_FOR_DOCKING}_atoms", log_dir=prepped_ligands_dir)
            logging.info(f"Skipping tiny ligand: {pdb_file.name}")
            continue

        # Sanitize PDB: keep only the ligand's HETATM triplet; drop all ATOM lines (protein)
        sanitized = pdb_file.with_suffix(".sanitized.pdb")

        def _sanitize_pdb(in_pdb: Path, out_pdb: Path) -> bool:
            target = None  # (resn, chain, resi)
            wrote_any = False

            def _triplet(ln: str):
                resn = (ln[17:20] if len(ln) >= 20 else "").strip()
                chain = (ln[21] if len(ln) >= 22 else " ").strip()
                resi = (ln[22:26] if len(ln) >= 26 else "").strip()
                return resn, chain, resi

            with open(in_pdb, "r", encoding="utf-8", errors="ignore") as fin:
                lines = fin.readlines()

            for ln in lines:
                if ln.startswith("HETATM"):
                    target = _triplet(ln); break

            if target:
                logging.info(f"[sanitize] {in_pdb.name}: keeping ligand {target[0]} chain {target[1]} resi {target[2]}")
            else:
                logging.warning(f"[sanitize] {in_pdb.name}: no HETATM found; nothing to keep")

            with open(out_pdb, "w", encoding="utf-8") as fout:
                for ln in lines:
                    if not ln.startswith(("ATOM", "HETATM")):
                        fout.write(ln); continue
                    if ln.startswith("ATOM"):
                        continue
                    if target is None:
                        continue
                    resn, chain, resi = _triplet(ln)
                    if (resn, chain, resi) != target:
                        continue
                    altloc = ln[16] if len(ln) > 16 else " "
                    keep_alt = {"", " ", "A"} | set("0123456789")
                    a = altloc.strip()
                    if a and (a not in keep_alt and a.upper() not in keep_alt):
                        continue
                    try:
                        x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                        if any([(x != x), (y != y), (z != z)]) or max(abs(x), abs(y), abs(z)) > 1e6:
                            continue
                    except Exception:
                        continue
                    fout.write(ln); wrote_any = True

            if not wrote_any:
                logging.warning(f"[control-prep] {in_pdb.name}: no HETATM kept after residue filtering.")
            return wrote_any and out_pdb.exists() and out_pdb.stat().st_size > 0

        ok_san = _sanitize_pdb(pdb_file, sanitized)
        if not ok_san:
            _log_malformed(pdb_file, "sanitize_kept_no_atoms")
            logging.warning(f"Sanitization produced no atoms: {pdb_file.name}")
            continue
        tmp_pdb = str(sanitized)

        # Element-field fixes (column-accurate) and helium guard preflight
        try:
            pre_txt = sanitized.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            pre_txt = None
        try:
            fix_element_columns_in_file(str(sanitized), rewrite_atoms=False)
            fixed_txt, nname = assert_no_helium_in_hydrogen_names(
                sanitized.read_text(encoding="utf-8", errors="ignore"))
            if nname > 0:
                sanitized.write_text(fixed_txt, encoding="utf-8")
            _log_elem_fix_summary(sanitized, stage="preflight", before_text=pre_txt,
                                  after_text=fixed_txt if nname > 0 else None)
        except Exception as e:
            logging.warning(f"[elements] preflight failed for {sanitized.name}: {e}")
        try:
            fix_pdb_elements(str(sanitized))
            logging.info(f"[elements] fixed element fields: {sanitized.name}")
        except Exception as e:
            logging.warning(f"[elements] could not fix elements for {sanitized.name}: {e}")

        # Counter-ion/buffer heuristics (quick bail-outs on obvious salts/buffers)
        flag_buffer = False
        reason = None
        try:
            m_chk = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
            if m_chk is not None:
                if _buffer_like_by_counts_from_mol(m_chk):
                    flag_buffer, reason = True, "buffer_like_by_counts"
                if not flag_buffer:
                    try:
                        Chem.SanitizeMol(m_chk)
                    except Exception:
                        pass
                    reason = _matches_counterion(m_chk) or (_looks_like_buffer_salt(m_chk) and "buffer_like")
                    if reason:
                        flag_buffer = True
                    if not flag_buffer and residue_name in {"UNL", "LIG"} and _is_polyacidic_buffer_like(m_chk):
                        flag_buffer, reason = True, "polyacidic_buffer_like"
        except Exception:
            if _buffer_like_by_counts_from_pdbfile(sanitized) or (
                residue_name in {"UNL", "LIG"} and _polyacidic_by_counts_from_pdbfile(sanitized)
            ):
                flag_buffer, reason = True, "counts_only_polyacidic"

        if flag_buffer:
            _log_malformed(pdb_file, f"counterion_or_buffer:{reason or 'unknown'}")
            logging.info(f"Skipping likely counter-ion/buffer ({reason or 'unknown'}): {pdb_file.name}")
            continue

        # Additional UNL O-rich ringless fragment guard
        try:
            if 'm_chk' in locals() and m_chk is not None and residue_name in {"UNL", "LIG"}:
                from rdkit.Chem import rdMolDescriptors as rdmd
                rings = rdmd.CalcNumRings(m_chk)
                arom = rdmd.CalcNumAromaticRings(m_chk)
                o = sum(1 for a in m_chk.GetAtoms() if a.GetSymbol() == "O")
                hac = m_chk.GetNumHeavyAtoms()
                if rings == 0 and arom == 0 and hac >= 12 and (o / float(hac)) >= 0.40:
                    _log_malformed(pdb_file, "UNL_O_rich_ringless", log_dir=prepped_ligands_dir)
                    logging.info(f"Skipping UNL O-rich ringless fragment: {pdb_file.name}")
                    continue
        except Exception:
            pass

        # Write a pristine SDF snapshot for reference/debug
        ref_dir = sanitized.parent.parent / "reference"
        ref_dir.mkdir(parents=True, exist_ok=True)
        ref_sdf = ref_dir / (pdb_file.stem + ".sdf")
        try:
            mref = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
            if mref is not None:
                w = Chem.SDWriter(str(ref_sdf))
                try: w.SetKekulize(False)
                except Exception: pass
                w.write(mref); w.close()
        except Exception:
            with open(ref_sdf, "w") as out:
                out.write(pdb_file.stem + "\n$$$$\n")

        # Resume check
        pdbqt_path = prepped_ligands_dir / f"{pdb_file.stem}.pdbqt"
        if pdbqt_path.exists() and pdbqt_path.stat().st_size > 100 and is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir):
            logging.info(f"[resume] Valid PDBQT already exists, skipping: {pdbqt_path.name}")
            continue

        # Primary PDB->MOL2 (no H) as a baseline file; later protonation layers may produce protoA/protoB
        tmp_mol2 = sanitized.with_suffix(".mol2")
        try:
            ob_cmd = [
                obabel_exe_short,
                "-ipdb", get_short_path_name(str(sanitized.resolve())),
                "-omol2", "-O", get_short_path_name(str(tmp_mol2.resolve()))
            ]
            logging.info("Converting sanitized PDB -> MOL2 (no gen3d): " + " ".join(map(str, ob_cmd)))
            res_ob = subprocess.run(ob_cmd, check=True, capture_output=True, text=True, timeout=300)
            if res_ob.stderr:
                logging.warning("[obabel primary stderr] %s", res_ob.stderr.strip())
        except subprocess.CalledProcessError as e:
            logging.warning(f"OBabel PDB->MOL2 failed for {sanitized.name}:\n{e.stderr}")

        if tmp_mol2.exists() and tmp_mol2.stat().st_size > 100:
            # Normalize aromaticity on MOL2 so downstream typing is stable
            try:
                _re_aromatize_mol2_in_place(tmp_mol2, obabel_exe_short)
            except Exception as e:
                logging.error("[ligprep] re_arom crash for %s: %s", tmp_mol2.name, e)

            # --- PROTONATION: layered strategy on the extracted fragment ---
            pre_metrics = _audit_protonation_metrics(tag=lig_id, mol_or_path=tmp_pdb, context="pre")
            h_strategy_label: List[str] = []
            fallbacks_used: List[str] = []

            # Strategy A: RDKit AddHs on sanitized PDB -> MOL2
            proto_mol2: Optional[str] = None
            try:
                mH = Chem.MolFromPDBFile(tmp_pdb, sanitize=False, removeHs=False)
                if mH is None:
                    raise ValueError("rdkit_from_pdb_failed")
                mH = Chem.AddHs(mH, addCoords=True)
                try:
                    Chem.Kekulize(mH, clearAromaticFlags=True)
                except Exception:
                    pass
                protoA = sanitized.with_suffix(".protoA.mol2")
                Chem.MolToMol2File(mH, str(protoA))
                if protoA.exists() and protoA.stat().st_size > 100:
                    proto_mol2 = str(protoA)
                    h_strategy_label.append("RDKit-AddHs")
            except Exception:
                fallbacks_used.append("RDKit-AddHs-fail")

            # Try Strategy B (OBabel -p -h) unconditionally
            protoB = sanitized.with_suffix(".protoB.mol2")
            ph = os.environ.get("LIGPREP_PH", "7.4")
            okB = _run_obabel([
                obabel_exe_short, "-ipdb", str(sanitized),
                "-omol2", "-O", str(protoB),
                "-p", str(ph), "--partialcharge", "gasteiger"
            ], timeout_sec=600)

            # Decide best candidate by explicit-H count
            cands = []
            for p in (protoA, protoB):
                if p.exists() and p.stat().st_size > 100:
                    cands.append((p, _count_explicit_H_in_mol2(p)))
            if cands:
                best, bestH = max(cands, key=lambda t: t[1])
                proto_mol2 = str(best)
                logging.info(f"[choose-mol2] picked={best.name} H={bestH} "
                             f"others={[(x[0].name, x[1]) for x in cands]}")
            else:
                # fallback to the baseline tmp_mol2 below
                proto_mol2 = None

            # Strategy C: no extra H; rely on ADT/Meeko typing additions
            if proto_mol2 is None:
                proto_mol2 = str(tmp_mol2)
                h_strategy_label.append("no-extra-H (ADT adds)")
                fallbacks_used.append("ADT/Meeko-rescue-possible")
                
            # --- Choose protoA/protoB by explicit H count if either exists ---
            protoA_path = sanitized.with_suffix(".protoA.mol2")
            protoB_path = sanitized.with_suffix(".protoB.mol2")

            # Collect existing candidates
            cands: List[Path] = [p for p in (protoA_path, protoB_path) if p and p.exists()]

            proto_mol2: str = ""  # string path we pass into RDKit/OBabel helpers
            mol2_for_mgl_path: Path = None  # Path object used by downstream helpers
            bestH: int = -1

            if cands:
                # Pick the proto with the highest explicit-H count
                counts = [(p, _count_explicit_H_in_mol2(p)) for p in cands]
                mol2_for_mgl_path, bestH = max(counts, key=lambda t: t[1])
                proto_mol2 = str(mol2_for_mgl_path)
                logging.info(
                    "[choose-mol2] picked=%s H=%s others=%s",
                    mol2_for_mgl_path.name, bestH,
                    [(p.name, h) for p, h in counts]
                )
            else:
                # Fallback: if neither proto exists, use the baseline tmp_mol2 produced earlier
                # (created as sanitized.with_suffix('.mol2') in the primary PDB→MOL2 step)
                fallback = tmp_mol2 if 'tmp_mol2' in locals() and tmp_mol2.exists() else sanitized.with_suffix(".mol2")
                mol2_for_mgl_path = Path(fallback)
                proto_mol2 = str(mol2_for_mgl_path)
                logging.warning(
                    "[choose-mol2] no protoA/protoB found; falling back to %s (exists=%s size=%s)",
                    mol2_for_mgl_path.name, mol2_for_mgl_path.exists(),
                    (mol2_for_mgl_path.stat().st_size if mol2_for_mgl_path.exists() else -1),
                )

            # Ensure the ADT/Meeko helpers get the expected variable name
            mol2_for_mgl: Path = Path(proto_mol2)

            # Keep aromaticity consistent on proto_mol2 as well
            try:
                _re_aromatize_mol2_in_place(mol2_for_mgl, obabel_exe_short)
            except Exception as e:
                logging.error("[ligprep] re_arom crash for %s: %s", mol2_for_mgl.name, e)

            # --- Enforce active-site standardization on extracted ligands (parity with bulk) ---
            try:
                # Parse the chosen proto_mol2 safely
                m_raw = Chem.MolFromMol2File(str(proto_mol2), sanitize=False, removeHs=False)
                old_smiles = ""
                try:
                    m_old = Chem.MolFromMol2File(str(proto_mol2), sanitize=True, removeHs=False)
                    if m_old:
                        old_smiles = Chem.MolToSmiles(m_old, isomericSmiles=True)
                except Exception:
                    pass

                m_std = standardize_mol_with_activesite(m_raw)
                if m_std is not None:
                    # Audit: if chemistry changed, log it
                    new_smiles = ""
                    try:
                        new_smiles = Chem.MolToSmiles(m_std, isomericSmiles=True)
                    except Exception:
                        pass
                    if old_smiles and new_smiles and old_smiles != new_smiles:
                        logging.warning(
                            "[std:audit][extracted] %s: SMILES changed %s -> %s",
                            Path(proto_mol2).name, old_smiles, new_smiles
                        )
                        _log_std_diff(prepped_ligands_dir, lig_id, "extracted", old_smiles, new_smiles)
            except Exception as e:
                logging.warning("[std:extracted] standardization skipped due to error: %s", e)

            # ... (writer selection remains unchanged) ...

            name, status = _prepare_one(
                mgltools_python_short, prepare_script_short,
                Path(proto_mol2), pdbqt_path,
                obabel_exe_short=obabel_exe_short,
                status_log_dir=prepped_ligands_dir
            )
            ok_write = (status == "ok")
            chosen_writer = "auto"
            reason = status


            # Post-write guard (extracted path): ensure explicit H present; OBabel re-write if not
            if pdbqt_path.exists() and obabel_exe_short and not _pdbqt_has_H(pdbqt_path):
                logging.warning("[post-extracted] %s has no H; OBabel -h re-write", pdbqt_path.name)
                _pdbqt_from_mol2_via_obabel(Path(proto_mol2), pdbqt_path, obabel_exe_short)

                logging.info("[post-extracted] re-write complete; has_H=%s", _pdbqt_has_H(pdbqt_path))

            # Post-write guards: helium + invariants (waters/ions/torsdof/charges/types)
            _he_ok, _he_fixes, _he_q = _helium_postwrite_guard(pdbqt_path, lig_id, prepped_ligands_dir)
            valid_ok = is_valid_ligand(pdbqt_path, prepped_ligands_dir)

            post_metrics = _audit_protonation_metrics(tag=lig_id, mol_or_path=pdbqt_path, context="post")

            # Concise one-line summary for logs
            print("[ligprep-extracted] "
                  f"lig={lig_id} input_fmt=PDB "
                  f"H_strategy={'+'.join(h_strategy_label) or 'none'} "
                  f"fallbacks={fallbacks_used or '[]'} "
                  f"writer={chosen_writer} ok={ok_write and valid_ok and not bool(_he_q)} reason={reason} "
                  f"final: H={post_metrics['h_count']} formal_charge={post_metrics['formal_charge']} "
                  f"has_partial_charges={post_metrics['has_partial_charges']} "
                  f"ad4_types_ok={post_metrics['ad4_types_ok']}")

            # TSV status/provenance row
            try:
                _append_prep_status(
                    status_log,
                    ligand_name=lig_id,
                    status=("OK" if (ok_write and valid_ok and not bool(_he_q)) else "FAIL"),
                    reason=reason or "",
                    relpath=str(pdbqt_path.name),
                    stage="write_pdbqt",
                    failure_code=("helium_rules" if _he_q else ""),
                    failure_detail=str(_he_q or ""),
                    fixes_count=int(_he_fixes or 0),
                    rules_ver=rules_version(),
                    writer_final=str(chosen_writer),
                    rescue_used=(",".join(fallbacks_used) if fallbacks_used else ""),
                    torsion_root_rule="",
                    polarH="",
                    ad4_types_ok=("yes" if post_metrics.get("ad4_types_ok") else "no"),
                    charges_ok=("yes" if post_metrics.get("has_partial_charges") else "no"),
                    torsdof="",
                )
            except Exception as _e:
                logging.warning("[status-log] append failed for %s: %s", lig_id, _e)

        # --- extra postcheck: ensure no *protein residues* slipped into ligand PDBQT
        try:
            has_protein_res = False
            if pdbqt_path.exists():
                with open(pdbqt_path, "r", errors="ignore") as fh:
                    for ln in fh:
                        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                            continue
                        resn = (ln[17:20] if len(ln) >= 20 else "").strip().upper()
                        if resn in STANDARD_AMINO_ACIDS:
                            has_protein_res = True
                            break
            if has_protein_res:
                _log_malformed(pdb_file, "protein_residue_in_pdbqt", log_dir=prepped_ligands_dir)
                quarantine = prepped_ligands_dir / QUARANTINE_DIRNAME
                quarantine.mkdir(exist_ok=True)
                try:
                    if pdbqt_path.exists():
                        pdbqt_path.replace(quarantine / pdbqt_path.name)
                except Exception:
                    pass
                try:
                    _append_prep_status(
                        status_log,
                        sanitized.name,
                        "FAIL",
                        "protein_residue_in_pdbqt",
                        str((quarantine / pdbqt_path.name).relative_to(prepped_ligands_dir))
                        if (quarantine / pdbqt_path.name).exists() else ""
                    )
                except Exception:
                    pass
                continue
        except Exception:
            logging.warning(f"[postcheck] unable to scan for protein ATOM in {pdbqt_path.name}")

        # Clean temp MOL2 if we succeeded in writing a nontrivial PDBQT
        try:
            if pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
                tmp_mol2.unlink(missing_ok=True)
        except Exception:
            pass

        # Final validation + quarantine (use correct log_dir)
        if (not pdbqt_path.exists()) or (pdbqt_path.stat().st_size < 100) or (not is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir)):
            _log_malformed(pdb_file, "pdbqt_postcheck_fail_or_small", log_dir=prepped_ligands_dir)
            quarantine = prepped_ligands_dir / QUARANTINE_DIRNAME
            quarantine.mkdir(exist_ok=True)
            try:
                if pdbqt_path.exists():
                    pdbqt_path.replace(quarantine / pdbqt_path.name)
            except Exception:
                pass
            _append_prep_status(
                status_log, sanitized.name, "FAIL", "postcheck_fail_or_small",
                str((quarantine / pdbqt_path.name).relative_to(prepped_ligands_dir))
                if (quarantine / pdbqt_path.name).exists() else ""
            )
            continue

        logging.info(f"Created PDBQT: {pdbqt_path.name}")
        _append_prep_status(status_log, sanitized.name, "OK", "", str(pdbqt_path.relative_to(prepped_ligands_dir)))



# =========================
# Main SDF ? MOL2 ? PDBQT pipeline
# =========================

def _valid_pdbqt(path: Path, log_dir: Path) -> bool:
    return path.exists() and path.stat().st_size > 100 and is_valid_ligand(path, log_dir=log_dir)


def prep_ligands_with_mgltools(*, force: bool = False, only: Optional[Set[str]] = None):
    # also honor an env var as a fallback (useful in batch/HPC)
    if not force:
        env_force = os.environ.get("LIGPREP_FORCE", "").strip().lower()
        force = env_force in {"1", "true", "yes", "y"}

    print("Starting ligand preparation")

    cfg = load_config("config.txt")
    validate_config(cfg)
    pdb_token = (
        os.environ.get("PDB_ID")
        or cfg.get("PDB_ID")
        or cfg.get("TARGET_PDB")
        or cfg.get("PDB")
        or cfg.get("INPUT_PDB")
        or cfg.get("PDB_FILE")
        or ""
    )
    pdb_id = Path(str(pdb_token)).stem.upper() if pdb_token else "LIGPREP"
    # >>> PATHS INIT START
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    # >>> PATHS INIT END
    # ---------- resolve numeric/env knobs you already had ----------
    MIN_TORS_DOF = _env_int("MIN_TORS_DOF", cfg.get("MIN_TORS_DOF", 0))
    MIN_PARENT_HEAVY = _env_int("MIN_PARENT_HEAVY", cfg.get("MIN_PARENT_HEAVY", 8))

    # ---------- NEW: resolve overrideable effective paths ----------
    # read ENV once here; CLI (if used) will populate these env vars in __main__
    in_sdf_env = (os.environ.get("LIGPREP_IN_SDF", "") or "").strip() or None
    in_sdf_dir_env = (os.environ.get("LIGPREP_IN_SDF_DIR", "") or "").strip() or None
    in_pdb_dir_env = (os.environ.get("LIGPREP_IN_PDB_DIR", "") or "").strip() or None
    mol2_dir_env = (os.environ.get("LIGPREP_MOL2_DIR", "") or "").strip() or None
    out_dir_env = (os.environ.get("LIGPREP_OUT_DIR", "") or "").strip() or None
    status_log_env = (os.environ.get("LIGPREP_STATUS_LOG", "") or "").strip() or None

    # default to config when not overridden
    # >>> LIGAND PATHS PATCH START
    ligands_raw_dir  = paths.ligand_output_dir
    prepped_lig_dir  = paths.prepped_ligands_dir
    ligands_mol2_dir = paths.ligands_mol2_dir
    # >>> LIGAND PATHS PATCH END

    ligand_extracted_dir = Path(in_sdf_dir_env).resolve() if in_sdf_dir_env else ligands_raw_dir
    ligands_mol2_dir = Path(mol2_dir_env).resolve() if mol2_dir_env else ligands_mol2_dir
    output_ligands_dir = Path(out_dir_env).resolve() if out_dir_env else prepped_lig_dir
    prepped_ligands_dir = output_ligands_dir
    # direct all malformed logs for this protein to its prepped_ligands dir
    global MALFORMED_DIR
    MALFORMED_DIR = prepped_ligands_dir
    # create working/output dirs
    output_ligands_dir.mkdir(parents=True, exist_ok=True)
    ligands_mol2_dir.mkdir(parents=True, exist_ok=True)

    # emit a compact audit banner (single line)
    print(
        "[paths.effective]"
        f" in_sdf={in_sdf_env or 'None'}"
        f" in_sdf_dir={ligand_extracted_dir}"
        f" in_pdb_dir={in_pdb_dir_env or 'None'}"
        f" mol2_dir={ligands_mol2_dir}"
        f" out_pdbqt_dir={output_ligands_dir}"
        f" status_log={(status_log_env or cfg.get('LIGAND_STATUS_LOG_BASENAME', 'ligand_prep_status.tsv'))}"
    )

    # env > config precedence for toggles/ints
    def _env_bool(name, default):
        v = os.environ.get(name)
        return (str(v).strip().lower() in {"1", "true", "yes", "on"}) if v is not None else bool(default)

    def _env_int(name, default):
        v = os.environ.get(name)
        return int(v) if v not in (None, "") else int(default)

    def _env_float(name, default):
        v = os.environ.get(name)
        return float(v) if v not in (None, "") else float(default)

    USE_RDKIT_FOR_3D = _env_bool("USE_RDKIT_FOR_3D", cfg.get("USE_RDKIT_FOR_3D", True))
    OBABEL_THREADS = _env_int("OBABEL_THREADS", cfg.get("OBABEL_THREADS", 50))
    OBABEL_TIMEOUT_S = _env_int("OBABEL_TIMEOUT_S", cfg.get("OBABEL_TIMEOUT_S", 900))
    CHUNK_SIZE = _env_int("LIGPREP_CHUNK_SIZE", cfg.get("LIGPREP_CHUNK_SIZE", 200))
    LIGPREP_PH = _env_float("LIGPREP_PH", cfg.get("LIGPREP_PH", 7.4))
    KEEP_NONPOLAR_H = _env_int("KEEP_NONPOLAR_H", cfg.get("KEEP_NONPOLAR_H", 1))
    MAX_HEAVY_ATOMS = _env_int("MAX_HEAVY_ATOMS", cfg.get("MAX_HEAVY_ATOMS", 1200))
    MIN_ATOMS_FOR_DOCKING = _env_int("MIN_ATOMS_FOR_DOCKING", cfg.get("MIN_ATOMS_FOR_DOCKING", 5))
    MIN_PARENT_HEAVY = _env_int("MIN_PARENT_HEAVY", cfg.get("MIN_PARENT_HEAVY", 8))
    output_ligands_dir.mkdir(parents=True, exist_ok=True)
    # If someone pointed MOL2s at 'prepped_ligands', redirect to 'ligands_mol2' (compat warning).
    if "prepped_ligands" in str(ligands_mol2_dir):
        suggested = Path(str(ligands_mol2_dir).replace("prepped_ligands", "ligands_mol2")).resolve()
        logging.warning("[compat] LIGANDS_MOL2_DIR points at prepped_ligands; redirecting to %s", suggested)
        ligands_mol2_dir = suggested

    ligands_mol2_dir.mkdir(parents=True, exist_ok=True)
    # Banner: show exactly where things will land for this run
    print(f"[ligprep] sdf_dir={ligands_mol2_dir / '_rdkit_embedded_sdf'} "
          f"mol2_dir={ligands_mol2_dir} "
          f"pdbqt_dir={output_ligands_dir}")

    mgltools_python = cfg["MGLTOOLS_PYTHON"]
    mgltools_path = cfg["MGLTOOLS_PATH"]
    obabel_cfg = cfg["OPENBABEL_PATH"]

    # Linux-safe obabel resolution
    obabel_exe = obabel_cfg
    obabel_exe_short = get_short_path_name(obabel_exe)

    if not os.environ.get("BABEL_DATADIR"):
        obabel_dir = Path(obabel_exe).resolve().parent
        data_dir = obabel_dir / "data"
        if data_dir.exists():
            os.environ["BABEL_DATADIR"] = str(data_dir)

    for label, p in [
        ("MGLTOOLS_PYTHON", mgltools_python),
        ("MGLTOOLS_PATH", mgltools_path),
        ("OPENBABEL_PATH", obabel_exe),
        ("LIGAND_EXTRACTED_DIR", ligand_extracted_dir),
        ("LIGANDS_MOL2_DIR", ligands_mol2_dir),
        ("OUTPUT_LIGANDS_DIR", output_ligands_dir),
    ]:
        if not str(p).strip():
            raise RuntimeError(f"Config value missing/empty: {label}")

    mgltools_python_short = get_short_path_name(mgltools_python)
    prepare_script = _resolve_prepare_ligand4(mgltools_path, cfg)
    if not prepare_script.exists():
        raise FileNotFoundError(f"prepare_ligand4.py not found at {prepare_script} (set PREPARE_LIGAND4 in config.txt)")
    prepare_script_short = get_short_path_name(str(prepare_script.resolve()))

    # Per-target status log (SDF pipeline)
    status_log = (
        Path(status_log_env).resolve()
        if (status_log_env and Path(status_log_env).suffix)
        else (output_ligands_dir / (status_log_env or cfg.get("LIGAND_STATUS_LOG_BASENAME", "ligand_prep_status.tsv")))
    )

    # =========================
    # TRUE TEST-MODE: prefer per-ligand SDFs and filter by ONLY
    # =========================
    rdkit_unit_sdf_dir = ligands_mol2_dir / "_rdkit_embedded_sdf"
    unit_sdfs = sorted(rdkit_unit_sdf_dir.glob("*.sdf")) if rdkit_unit_sdf_dir.is_dir() else []

    if unit_sdfs:
        print(f"[test-mode] per-ligand SDFs detected dir={rdkit_unit_sdf_dir}")
        # Build the candidate SDF list
        selected_sdfs: List[Path]
        if only:
            # only is already normalized to canonical rdk_{7d} by _collect_only_from_env_and_cli
            requested_ids = sorted(only)
            # map requested ids to files (only those that exist)
            selected_sdfs = [rdkit_unit_sdf_dir / f"{rid}.sdf" for rid in requested_ids if
                             (rdkit_unit_sdf_dir / f"{rid}.sdf").exists()]
            missing = [rid for rid in requested_ids if not (rdkit_unit_sdf_dir / f"{rid}.sdf").exists()]

            # diagnostics (bounded)
            print(f"[test-mode] enabled only_count={len(requested_ids)} sample={','.join(requested_ids[:10])}")
            print(f"[test-mode] selected_sdf_count={len(selected_sdfs)} missing={len(missing)}"
                  + (f" first_missing={','.join(missing[:10])}" if missing else ""))

            if len(selected_sdfs) == 0:
                print("[test-mode] No requested per-ligand SDFs found. Nothing to do; exiting cleanly.")
                return
        else:
            # ONLY not set: prefer all per-ligand SDFs instead of bulk
            selected_sdfs = unit_sdfs
        # Optional cleanup of existing PDBQTs for the selected subset (before scheduling)
        clean_env = os.environ.get("LIGPREP_CLEAN", "").strip().lower() in {"1", "true", "yes", "y"}
        if clean_env and selected_sdfs:
            cleaned = []
            for sdf_in in selected_sdfs:
                out_pdbqt = output_ligands_dir / (sdf_in.stem + ".pdbqt")
                if out_pdbqt.exists():
                    try:
                        out_pdbqt.unlink()
                        cleaned.append(out_pdbqt.name)
                    except Exception:
                        pass
            if cleaned:
                logging.info("[clean] removed_existing=%d sample=%s", len(cleaned), ",".join(cleaned[:5]))
                print(f"[clean] removed={len(cleaned)} sample={','.join(cleaned[:5])}")

        # Convert each selected per-ligand SDF → MOL2 (idempotent unless --force)
        mol2_files: List[Path] = []
        for sdf_in in selected_sdfs:
            out_mol2 = ligands_mol2_dir / (sdf_in.stem + ".mol2")
            if out_mol2.exists() and out_mol2.stat().st_size > 100 and not force:
                mol2_files.append(out_mol2)
                continue

            cmd = [obabel_exe_short, "-isdf", get_short_path_name(str(sdf_in.resolve())), "-omol2",
                   "-O", get_short_path_name(str(out_mol2.resolve()))]
            print(f"[ligprep] pre-MOL2-write: per-ligand sdf={sdf_in.name} -> mol2={out_mol2.name}")
            ok = _run_obabel(cmd, timeout_sec=600)
            if ok and out_mol2.exists() and out_mol2.stat().st_size > 100:
                mol2_files.append(out_mol2)

        if not mol2_files:
            print("[test-mode] No MOL2 files produced from per-ligand SDFs; exiting.")
            return

        # If ONLY set, show concise banner of the subset that will be scheduled
        if only:
            sample = ",".join(sorted({Path(p).stem for p in mol2_files})[:10])
            print(f"[test-mode] scheduling_only={len(mol2_files)} sample={sample}")

        # ===== Schedule only the selected MOL2s for PDBQT (same logic as bulk path) =====
        print(f"Preparing {len(mol2_files)} MOL2 files with MGLTools (parallel) into {output_ligands_dir}")
        max_workers = max(1, int(os.environ.get("CPU", "8")))
        futures = []
        ok_count = 0
        fail_count = 0
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            for mol2_file in mol2_files:
                norm_mol2 = collapse_sanitized_once(Path(mol2_file))
                if norm_mol2 != mol2_file:
                    try:
                        norm_mol2.write_bytes(Path(mol2_file).read_bytes())
                        mol2_file = norm_mol2
                        logging.info("[tidy] normalized double-sanitized -> %s", mol2_file.name)
                    except Exception as e:
                        logging.warning("[tidy] unable to normalize %s: %s", mol2_file, e)

                pdbqt_path = output_ligands_dir / f"{Path(mol2_file).stem}.pdbqt"

                resume_skip = (
                        not force
                        and pdbqt_path.exists()
                        and pdbqt_path.stat().st_size > 100
                        and is_valid_ligand(pdbqt_path, log_dir=output_ligands_dir)
                )
                if resume_skip:
                    logging.info("[resume] Valid PDBQT exists, skipping: %s", pdbqt_path.name)
                    resume_skips = (resume_skips + 1) if 'resume_skips' in locals() else 1
                    if 'resume_examples' not in locals(): resume_examples = []
                    if len(resume_examples) < 10: resume_examples.append(pdbqt_path.name)
                    continue

                if force and pdbqt_path.exists():
                    logging.info("[force] Overwriting existing PDBQT: %s", pdbqt_path.name)
                # --- Enforce active-site standardization on bulk MOL2s (explicit, before _prepare_one) ---
                try:
                    m_raw = Chem.MolFromMol2File(str(mol2_file), sanitize=False, removeHs=False)
                    old_smiles = ""
                    try:
                        m_old = Chem.MolFromMol2File(str(mol2_file), sanitize=True, removeHs=False)
                        if m_old:
                            old_smiles = Chem.MolToSmiles(m_old, isomericSmiles=True)
                    except Exception:
                        pass

                    m_std = standardize_mol_with_activesite(m_raw)
                    if m_std is not None:
                        new_smiles = ""
                        try:
                            new_smiles = Chem.MolToSmiles(m_std, isomericSmiles=True)
                        except Exception:
                            pass
                        if old_smiles and new_smiles and old_smiles != new_smiles:
                            logging.warning("[std:audit][bulk] %s: SMILES changed %s -> %s",
                                            Path(mol2_file).name, old_smiles, new_smiles)
                            _log_std_diff(output_ligands_dir, Path(mol2_file).stem, "bulk", old_smiles, new_smiles)

                        std_path = Path(mol2_file).with_suffix(".std.mol2")
                        Chem.MolToMol2File(m_std, str(std_path))
                        mol2_file = std_path  # pass standardized path downstream
                except Exception as e:
                    logging.warning("[std][bulk] skip for %s: %s", Path(mol2_file).name, e)
                # --- end standardization parity block ---

                futures.append(ex.submit(
                    _prepare_one,
                    mgltools_python_short, prepare_script_short,
                    Path(mol2_file), pdbqt_path,
                    obabel_exe_short=obabel_exe_short,
                    status_log_dir=output_ligands_dir
                ))

            total = len(futures)
            for i, fut in enumerate(as_completed(futures), 1):
                name, status = fut.result()
                lig_stem = Path(name).stem
                if status == "ok":
                    ok_count += 1
                    try:
                        _append_prep_status(status_log, lig_stem, "OK", "", f"{lig_stem}.pdbqt")
                    except Exception:
                        pass
                else:
                    fail_count += 1
                    try:
                        _append_prep_status(status_log, lig_stem, "FAIL", status, f"{lig_stem}.pdbqt")
                    except Exception:
                        pass
                if i % 100 == 0 or status != "ok":
                    print(f"[{i}/{total}] {name}: {status}")
            try:
                print(
                    f"[ligprep] scheduled={len(futures)} resume_skips={resume_skips if 'resume_skips' in locals() else 0}"
                    f" force={force} examples_skipped={(resume_examples if 'resume_examples' in locals() else [])[:5]}")
            except Exception:
                pass

        # Done with per-ligand “test-mode” path; avoid touching bulk SDFs.
        return
    # =========================
    # END test-mode, fall back to legacy bulk-SDF path below
    # =========================

    # : crystal-safe dispatch (takes precedence over bulk scan when in_sdf not set)
    if in_pdb_dir_env and not in_sdf_env:
        return prep_ligands_from_pdb(Path(in_pdb_dir_env).resolve(), ligands_mol2_dir, output_ligands_dir)

    #  single-file override; else scan directory as before
    sdf_files = [Path(in_sdf_env).resolve()] if in_sdf_env else list(ligand_extracted_dir.glob("*.sdf"))
    print(f"Found {len(sdf_files)} SDF file(s)")

    if not sdf_files:
        return

    for sdf_file in sdf_files:
        print(f"=== Processing SDF: {sdf_file.name} ===")

        # immediately normalize SDF names to avoid '.sanitized.sanitized'
        norm_sdf = collapse_sanitized_once(Path(sdf_file))
        if norm_sdf != Path(sdf_file):
            try:
                norm_sdf.write_bytes(Path(sdf_file).read_bytes())
                logging.info("[tidy] normalized double-sanitized SDF -> %s", norm_sdf.name)
                sdf_file = norm_sdf
            except Exception as e:
                logging.warning("[tidy] unable to normalize SDF %s: %s", sdf_file, e)
        sdf_abs = Path(sdf_file).resolve()

        sdf_abs = sdf_file.resolve()

        max_workers = max(1, int(os.environ.get("CPU", "8")))

        if USE_RDKIT_FOR_3D:
            print("Using RDKit ETKDG for 3D with parent-picking; OBabel only for format conversion ")
            mol2_files = rdkit_embed_sdf_to_mol2(
                sdf_abs, ligands_mol2_dir, obabel_exe=obabel_exe_short, max_workers=max_workers
            )
        else:
            print("Using OBabel --gen3d; pre-cleaning SDF to parent-only ")
            cleaned_sdf = ligands_mol2_dir / (sdf_abs.stem + "_parents.sdf")
            n_kept = _write_parent_only_sdf(sdf_abs, cleaned_sdf)
            print(f"Parent-only SDF kept {n_kept} records")
            if n_kept == 0:
                print("No parent molecules survived desalting; skipping.")
                continue
            mol2_files = convert_sdf_to_mol2_split_parallel(
                cleaned_sdf, ligands_mol2_dir, obabel_exe_short,
                threads=OBABEL_THREADS, timeout_sec=OBABEL_TIMEOUT_S, chunk_size=CHUNK_SIZE
            )

        if not mol2_files:
            for _root in (ligands_mol2_dir, output_ligands_dir):
                probe = sorted(_root.glob("*.mol2"))
                if probe:
                    logging.warning("[compat] Found pre-existing MOL2s under %s (DEPRECATED layout); continuing.",
                                    _root)
                    mol2_files = probe
                    break

        if not mol2_files:
            print("No MOL2 files produced; skipping this SDF.")
            continue


        # ----- TEST MODE FILTERING (ONLY-set) -----
        if only:
            stems = {Path(p).stem for p in mol2_files}
            requested = set(only)
            # Filter down
            mol2_files = [p for p in mol2_files if Path(p).stem in requested]

            # Diagnostics banner (single, concise line)
            sample = sorted(list(requested))[:10]
            remain = len(mol2_files)
            logging.info("[test-mode] enabled count=%d remain=%d sample=%s", len(requested), remain, ",".join(sample))
            print(f"[test-mode] enabled count={len(requested)} remain={remain} sample={','.join(sample)}")

            # Warn once about missing
            missing = sorted(list(requested - stems))
            if missing:
                logging.warning("[test-mode] requested ligands not found among MOL2s: %s",
                                ",".join(missing[:20]) + ("..." if len(missing) > 20 else ""))

            # Early exit if nothing remains
            if not mol2_files:
                print("[test-mode] No requested ligands were found. Nothing to do; exiting cleanly.")
                return



        print(f"Preparing {len(mol2_files)} MOL2 files with MGLTools (parallel) into {output_ligands_dir}")
        max_workers = max(1, int(os.environ.get("CPU", "8")))
        futures = []
        ok_count = 0
        fail_count = 0
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            for mol2_file in mol2_files:
                norm_mol2 = collapse_sanitized_once(Path(mol2_file))
                if norm_mol2 != mol2_file:
                    try:
                        norm_mol2.write_bytes(Path(mol2_file).read_bytes())
                        mol2_file = norm_mol2
                        logging.info("[tidy] normalized double-sanitized -> %s", mol2_file.name)
                    except Exception as e:
                        logging.warning("[tidy] unable to normalize %s: %s", mol2_file, e)

                # ← now back at the for-loop level
                pdbqt_path = output_ligands_dir / f"{mol2_file.stem}.pdbqt"

                resume_skip = (
                        not force
                        and pdbqt_path.exists()
                        and pdbqt_path.stat().st_size > 100
                        and is_valid_ligand(pdbqt_path, log_dir=output_ligands_dir)
                )
                if resume_skip:
                    logging.info("[resume] Valid PDBQT exists, skipping: %s", pdbqt_path.name)
                    resume_skips = (resume_skips + 1) if 'resume_skips' in locals() else 1
                    if 'resume_examples' not in locals(): resume_examples = []
                    if len(resume_examples) < 10: resume_examples.append(pdbqt_path.name)
                    continue

                if force and pdbqt_path.exists():
                    logging.info("[force] Overwriting existing PDBQT: %s", pdbqt_path.name)

                futures.append(ex.submit(
                    _prepare_one,
                    mgltools_python_short, prepare_script_short,
                    mol2_file, pdbqt_path,
                    obabel_exe_short=obabel_exe_short,
                    status_log_dir=output_ligands_dir
                ))

            total = len(futures)
            for i, fut in enumerate(as_completed(futures), 1):
                name, status = fut.result()
                lig_stem = Path(name).stem
                if status == "ok":
                    ok_count += 1
                    try:
                        _append_prep_status(status_log, lig_stem, "OK", "", f"{lig_stem}.pdbqt")
                    except Exception:
                        pass
                else:
                    fail_count += 1
                    try:
                        _append_prep_status(status_log, lig_stem, "FAIL", status, f"{lig_stem}.pdbqt")
                    except Exception:
                        pass
                if i % 100 == 0 or status != "ok":
                    print(f"[{i}/{total}] {name}: {status}")
            try:
                print(
                    f"[ligprep] scheduled={len(futures)} resume_skips={resume_skips if 'resume_skips' in locals() else 0}"
                    f" force={force} examples_skipped={(resume_examples if 'resume_examples' in locals() else [])[:5]}")
            except Exception:
                pass

        def _as_path(p) -> Path:
            return p if isinstance(p, Path) else Path(p)

        def _cfg_env_or_default(key: str, default: Optional[str] = None) -> Optional[str]:
            """Lightweight config reader that prefers env, then config.txt next to this file, else default."""
            v = os.environ.get(key)
            if v:
                return v
            # try config.txt next to this file (your project already uses this pattern)
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
            """Canonical per-protein base dir: processed_pdbs/<PDB>"""
            output_root = _as_path(output_root)
            return (output_root / pdb_id.upper()).resolve()

        def _merge_dir(src: Path, dst: Path) -> None:
            """Merge src directory into dst (mkdirs as needed); removes src if emptied."""
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
                        # prefer keeping existing canonical artifacts; only overwrite if target is missing
                        try:
                            # If same file, skip; else overwrite (safe in our case)
                            if s.stat().st_size == t.stat().st_size:
                                continue
                        except Exception:
                            pass
                    shutil.move(str(s), str(t))
            # try to remove empty src tree
            try:
                shutil.rmtree(src)
            except Exception:
                pass

        def fold_legacy_layout(pdb_id: str, output_root) -> None:
            """
            Override with extended handling to pull uppercase legacy dirs into the canonical tree:
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

                # legacy directories (both lower and UPPER)
                legacy_dirs = [
                    (root / f"{pdb_id}_nolig", base / "nolig"),
                    (root / f"{pdb_id.lower()}_nolig", base / "nolig"),
                    (root / f"{pdb_idU}_NOLIG", base / "nolig"),
                    (root / f"{pdb_id}_cleaned_ligands", base / "ligands_raw"),
                    (root / f"{pdb_id.lower()}_cleaned_ligands", base / "ligands_raw"),
                    (root / f"{pdb_idU}_CLEANED_LIGANDS", base / "ligands_raw"),
                ]
                for src, dst in legacy_dirs:
                    if src.exists():
                        logging.info("Migrating legacy directory %s -> %s", src, dst)
                        _merge_dir(src, dst)

                # legacy files
                candidates_files = [
                    (root / f"{pdb_id}_nolig.pdb", base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb"),
                    (root / f"{pdb_id.lower()}_nolig.pdb", base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb"),
                ]
                for src, dst in candidates_files:
                    if src.exists() and not dst.exists():
                        logging.info("Moving legacy file %s -> %s", src, dst)
                        dst.parent.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(src), str(dst))
            except Exception as e:
                logging.warning("fold_legacy_layout (extended) failed for %s: %s", pdb_id, e)

        def _expose_ligand_intermediates_for_debug(pdb_id: str, ligands_raw: Path) -> None:
            """
            If PREPPED_LIGANDS_ROOT is configured, create/update:
              <PREPPED_LIGANDS_ROOT>/<PDB>/intermediates -> <processed_pdbs>/<PDB>/ligands_raw
            so intermediates (sanitized PDB, MOL2) are visible next to final PDBQTs.
            """
            try:
                prepped_root = _cfg_env_or_default("PREPPED_LIGANDS_ROOT", "")
                if not prepped_root:
                    return
                dst_dir = _as_path(prepped_root) / pdb_id.upper()
                dst_dir.mkdir(parents=True, exist_ok=True)
                link = dst_dir / "intermediates"
                if link.is_symlink() or link.exists():
                    try:
                        if link.is_dir() and not link.is_symlink():
                            shutil.rmtree(link)
                        else:
                            link.unlink()
                    except Exception:
                        pass
                link.symlink_to(ligands_raw.resolve(), target_is_directory=True)
                logging.info("Debug symlink: %s -> %s", link, ligands_raw)
            except Exception as e:
                logging.warning("Could not create debug symlink for %s: %s", pdb_id, e)

        def _element_fix_all_in_dir(ligands_raw: Path) -> None:
            """
            Run your element-column fixer on every PDB in ligands_raw if the function exists.
            This reduces 'Unknown atom name F29 -> C' noise and improves template matches downstream.
            """
            try:
                fixer = globals().get("fix_element_columns_in_file", None)
                if fixer is None:
                    return
                for p in ligands_raw.glob("*.pdb"):
                    try:
                        fixer(p, p)
                    except Exception as e:
                        logging.warning("Element-fix skipped for %s: %s", p.name, e)
            except Exception as e:
                logging.warning("Bulk element-fix failed in %s: %s", ligands_raw, e)

        # Wrap/extend clean_pdb to: (1) fold legacy for this PDB, (2) expose intermediates, (3) element-fix newly extracted ligands.
        if "clean_pdb" in globals():
            _orig_clean_pdb = clean_pdb  # type: ignore[misc]

            def clean_pdb(pdb_file, output_root):
                pdb_file = _as_path(pdb_file)
                pdb_id = pdb_file.stem.upper()
                # 1) Make sure any legacy dirs for this PDB are migrated before we proceed
                fold_legacy_layout(pdb_id, output_root)
                # 2) Run the original pipeline
                out = _orig_clean_pdb(pdb_file, output_root)
                fold_legacy_layout(pdb_id, output_root)

                # 3) Expose intermediates via a symlink next to prepped_ligands/<PDB>
                try:
                    # Use the same canonical path helper your code already uses, if present
                    if "canon_paths" in globals():
                        paths = globals()["canon_paths"](pdb_id, output_root)  # type: ignore[index]
                        lig_raw = paths.get("ligands_raw", None)
                        if lig_raw:
                            _element_fix_all_in_dir(lig_raw)
                            _expose_ligand_intermediates_for_debug(pdb_id, lig_raw)
                except Exception as e:
                    logging.warning("post-clean_pdb expose failed for %s: %s", pdb_id, e)
                return out

if __name__ == "__main__":
    import argparse, sys

    parser = argparse.ArgumentParser(description="Bulk ligand preparation (SDF→MOL2→PDBQT)")
    parser.add_argument("--force", action="store_true",
                        help="Do not skip existing PDBQTs; overwrite outputs even if a valid PDBQT exists.")
    parser.add_argument("--migrate-legacy", metavar="PROCESSED_ROOT",
                        help="Scan processed_pdbs and migrate *_NOLIG/*_CLEANED_LIGANDS into canonical layout.")
    parser.add_argument("--expose-intermediates", nargs=2, metavar=("PROCESSED_ROOT", "PREPPED_LIGANDS_ROOT"),
                        help="Create/refresh prepped_ligands/<PDB>/intermediates symlinks for all PDBs.")
    parser.add_argument("--only", nargs="+",
                        help="Limit run to specific ligands (accepts rdk_0004931, rdk_4931, 0004931, 4931; comma/space OK)")
    parser.add_argument("--extracted", action="store_true",
                        help="Run ligand prep for extracted-from-PDB ligands (ligands_raw/*.pdb).")
    parser.add_argument("--only-extracted", nargs="+", default=None,
                        help="Limit extracted prep to specific ligand basenames (e.g., RXT_A1204). Mirrors EXTRACT_ONLY.")
    parser.add_argument("--test-extracted", action="store_true",
                        help="Verbose test mode for extracted ligands (mirrors EXTRACT_TEST=1).")
    # --- optional I/O override flags (non-breaking) ---
    parser.add_argument("--in-sdf", help="Path to a single SDF file (process only this file)")
    parser.add_argument("--in-sdf-dir", help="Directory of bulk SDFs (replaces default extracted-SDF root)")
    parser.add_argument("--in-pdb-dir", help="Directory of extracted ligand PDBs (crystal-safe path)")
    parser.add_argument("--mol2-dir", help="Directory for MOL2 intermediates")
    parser.add_argument("--out-pdbqt-dir", help="Destination directory for final PDBQTs")
    parser.add_argument("--status-log", help="Basename or full path for the TSV status log")

    args = parser.parse_args()


    # Promote CLI overrides to ENV so the core pipeline (and your config precedence)
    # can pick them up without changing existing logic. Precedence stays: CLI > ENV > config.
    def _set_env(k, v):
        if v is not None and str(v).strip() != "":
            os.environ[k] = str(v)


    _set_env("LIGPREP_IN_SDF", args.in_sdf)
    _set_env("LIGPREP_IN_SDF_DIR", args.in_sdf_dir)
    _set_env("LIGPREP_IN_PDB_DIR", args.in_pdb_dir)
    _set_env("LIGPREP_MOL2_DIR", args.mol2_dir)
    _set_env("LIGPREP_OUT_DIR", args.out_pdbqt_dir)
    _set_env("LIGPREP_STATUS_LOG", args.status_log)

    # Compact one-line audit of the effective overrides seen at startup
    print(
        "[paths.effective]",
        f"in_sdf={os.environ.get('LIGPREP_IN_SDF', 'None')}",
        f"in_sdf_dir={os.environ.get('LIGPREP_IN_SDF_DIR', 'None')}",
        f"in_pdb_dir={os.environ.get('LIGPREP_IN_PDB_DIR', 'None')}",
        f"mol2_dir={os.environ.get('LIGPREP_MOL2_DIR', 'None')}",
        f"out_pdbqt_dir={os.environ.get('LIGPREP_OUT_DIR', 'None')}",
        f"status_log={os.environ.get('LIGPREP_STATUS_LOG', 'None')}",
    )

    # Build ONLY set from env + file + CLI
    only_set = _collect_only_from_env_and_cli(args.only if hasattr(args, "only") else None)
    # --- extracted path entrypoint ---
    if getattr(args, "extracted", False):
        # surface ONLY/test via env for the function to read
        if getattr(args, "only_extracted", None):
            os.environ["EXTRACT_ONLY"] = " ".join(args.only_extracted)
        if getattr(args, "test_extracted", False):
            os.environ["EXTRACT_TEST"] = "1"

        cfg = read_config()
        ligand_extracted_dir = Path(cfg["LIGAND_EXTRACTED_DIR"]).resolve()
        ligands_mol2_dir = Path(cfg["LIGANDS_MOL2_DIR"]).resolve()
        output_ligands_dir = Path(cfg["OUTPUT_LIGANDS_DIR"]).resolve()

        # prep_ligands_from_pdb expects (ligand_output_dir, ligands_mol2_dir, prepped_ligands_dir)
        prep_ligands_from_pdb(ligand_extracted_dir, ligands_mol2_dir, output_ligands_dir)
        sys.exit(0)

    # Handle utility modes first, then exit
    if args.migrate_legacy:
        root = Path(args.migrate_legacy).resolve()
        for d in sorted(root.iterdir()):
            if not d.is_dir():
                continue
            m = re.match(r"^([A-Za-z0-9]{4})(?:_.+)?$", d.name)
            if not m:
                continue
            pdb_id = m.group(1).upper()
            try:
                fold_legacy_layout(pdb_id, root)  # defined above
            except Exception as e:
                logging.warning("migrate-legacy skip %s: %s", pdb_id, e)
        sys.exit(0)

    if args.expose_intermediates:
        proc_root = Path(args.expose_intermediates[0]).resolve()
        prepped_root = Path(args.expose_intermediates[1]).resolve()
        os.environ["PREPPED_LIGANDS_ROOT"] = str(prepped_root)
        for pdb_dir in sorted(proc_root.iterdir()):
            if not pdb_dir.is_dir():
                continue
            pdb_id = pdb_dir.name.split("_")[0].upper()
            lig_raw = pdb_dir / "ligands_raw"
            if lig_raw.is_dir():
                _element_fix_all_in_dir(lig_raw)
                _expose_ligand_intermediates_for_debug(pdb_id, lig_raw)
        print("Exposed intermediates under:", prepped_root)
        sys.exit(0)

    # Normal run: single call, preserve ONLY
    prep_ligands_with_mgltools(force=args.force, only=(only_set if only_set else None))


