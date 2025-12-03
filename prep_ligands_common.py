"""Common ligand prep helpers shared between bulk and crystal workflows."""
import ctypes
import logging
import os
import shutil
import subprocess
import sys
from ctypes import create_unicode_buffer, wintypes
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from collections import defaultdict

from activesite import assert_no_helium_in_pdbqt, get_atom_rules, rules_version, scan_helium_counts
from rdkit import Chem
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize as _std  # unified handle

    _HAS_STD = True
except Exception:
    _std = None
    _HAS_STD = False

from prep_ligands_bulk_sdf import _run_obabel

logger = logging.getLogger(__name__)

# Shared constants and global state
STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "SEC", "PYL", "MSE"
}

try:
    from chemdb.chem_alias_db import EXCLUDE_HET_IDS as EXCLUDE_CRYSTAL_ADDITIVES
except Exception:
    EXCLUDE_CRYSTAL_ADDITIVES = {
        "HOH","CIT","TAR","SO4","PO4","CA","NA","K","MG","MN","ZN","GOL","EDO","PEG","MPD","TRS",
        "MES","HEPES","ACET","ACT","FMT","MAL","DMS","IPA","CLU","NAG","BOG","TOS","BES","PTS","OTF",
        "TRF","TFA","BF4","PF6","CL","BR","I","HEM","HEC","HEA","HEB","HEO","HEG","HEF","HEH","PTR",
        "TPO","SEP"
    }

MAX_HEAVY_ATOMS = 1200
MIN_ATOMS_FOR_DOCKING = 3
QUARANTINE_DIRNAME = "quarantine"

MALFORMED_LOG = Path("malformed_ligands.txt")
MALFORMED_DIR: Path | None = None

RULES = None
ALLOWED_ELEMENTS: Set[str] = set()
MONOATOMIC_IONS: Set[str] = set()


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

_init_alias_rules_cache()

# --- Counter-ion SMARTS (shared between filters) ---
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

_CARBOXYLATE = Chem.MolFromSmarts("[CX3](=O)[O-]")
_CARBOXYLIC = Chem.MolFromSmarts("[CX3](=O)O")


def _looks_like_monoatomic_ion_pdbqt(lines: List[str]) -> bool:
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    if len(atom_lines) != 1:
        return False
    parts = atom_lines[0].split()
    elem = _element_from_adt(parts[-1]) if parts else ""
    return elem in MONOATOMIC_IONS


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


try:
    AD4_TYPES: set[str] = _get_ad4_types_from_aliases()
except Exception:
    AD4_TYPES = {
        "C","A","N","NA","O","OA","S","SA","H","HD","F","CL","BR","I","P","B","SI","SE",
        "ZN","MG","CA","MN","FE","K","NA","CU","CO","NI","AL","AG","AU","PT","LI","BA","SR","CS","RB"
    }

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


def add_hydrogens_mol2(
    in_path: Path,
    out_path: Path,
    obabel_exe: str,
    ph: Optional[float] = None,
) -> tuple[bool, str]:
    """
    Ensure explicit H before MGLTools. Returns (ok, stderr_text).
    Uses obabel -h to add hydrogens *without* changing atom order more than necessary.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [obabel_exe, "-imol2", str(in_path), "-omol2", "-O", str(out_path)]
    if ph is not None:
        cmd.extend(["-p", str(ph)])
    cmd.append("-h")
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


def read_config(path="config.txt") -> Dict[str, str]:
    config: Dict[str, str] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                config[key.strip()] = value.strip()
    return config


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
        status_log_dir: Path,
        ph: Optional[float] = None,
        copy_targets: Optional[List[Path]] = None,
) -> Tuple[str, str]:
    # use the H-enriched input

    lig_id = mol2_file.stem
    writer_final = "mgltools"
    rescue_used = False
    torsion_rule_label = "adt_default"

    def _relpath_for_status(target: Path) -> str:
        try:
            return str(target.relative_to(status_log_dir))
        except Exception:
            return target.name

    copy_targets = copy_targets or []

    logger.debug(
        "prep_ligands._prepare_one: ligand=%s ph=%.2f mol2=%s pdbqt=%s copy_targets=%d",
        lig_id,
        ph if ph is not None else float("nan"),
        mol2_file,
        pdbqt_path,
        len(copy_targets),
    )

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
                add_hydrogens_mol2(mol2_in, mol2_in, obabel_exe_short, ph=ph)  # in-place re-add
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
    okH, hstderr = add_hydrogens_mol2(mol2_in, mol2_h, obabel_exe_short, ph=ph)
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
        add_hydrogens_mol2(mol2_for_mgl, mol2_for_mgl, obabel_exe_short, ph=ph)

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
        add_hydrogens_mol2(mol2_for_mgl, mol2_for_mgl, obabel_exe_short, ph=ph)

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
                ok_ob = _pdbqt_from_mol2_via_obabel(mol2_for_mgl, pdbqt_path, obabel_exe_short)
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
                    relpath=_relpath_for_status(qpath),
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
                    relpath=_relpath_for_status(pdbqt_path),
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
                    relpath=_relpath_for_status(pdbqt_path),
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
                relpath=_relpath_for_status(qpath if qpath.exists() else pdbqt_path),
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
                relpath=_relpath_for_status(qpath),
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
                relpath=_relpath_for_status(pdbqt_path),
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
        for extra_target in copy_targets:
            try:
                if extra_target.resolve() == pdbqt_path.resolve():
                    continue
                extra_target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(pdbqt_path, extra_target)
            except Exception as e:
                logging.warning("[microstate] pdbqt_copy_failed src=%s dest=%s err=%s", pdbqt_path, extra_target, e)
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


__all__ = [
    "EXCLUDE_CRYSTAL_ADDITIVES",
    "MIN_ATOMS_FOR_DOCKING",
    "QUARANTINE_DIRNAME",
    "STANDARD_AMINO_ACIDS",
    "_append_prep_status",
    "_audit_protonation_metrics",
    "_buffer_like_by_counts_from_mol",
    "_buffer_like_by_counts_from_pdbfile",
    "_count_aromatic_atoms_in_mol2",
    "_count_explicit_H_in_mol2",
    "_helium_postwrite_guard",
    "_is_polyacidic_buffer_like",
    "_log_elem_fix_summary",
    "_log_malformed",
    "_log_std_diff",
    "_matches_counterion",
    "_pdbqt_from_mol2_via_obabel",
    "_pdbqt_has_H",
    "_polyacidic_by_counts_from_pdbfile",
    "_prepare_one",
    "_re_aromatize_mol2_in_place",
    "_resolve_obabel_exe",
    "_resolve_prepare_ligand4",
    "_run_obabel",
    "_looks_like_buffer_salt",
    "_write_aromatic_sdf",
    "get_short_path_name",
    "is_valid_ligand",
    "read_config",
    "standardize_mol_with_activesite",
]
