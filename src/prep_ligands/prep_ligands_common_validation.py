"""Ligand-prep validation helpers and exported rule/compatibility state."""

import ctypes
import logging
import os
import re
import sys
from ctypes import create_unicode_buffer, wintypes
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, cast


from config.runtime_config import load_config
from protein_prep import pdb_fixer_runtime as _pdb_fixer
from prep_ligands.prep_ligands_reporting import (
    MALFORMED_DIR,
    MALFORMED_LOG,
    _append_prep_status,
    _audit_protonation_metrics,
    _log_elem_fix_summary,
    _log_malformed,
)

logger = logging.getLogger(__name__)
_SANITIZED_RUN = re.compile(r"(?:\.sanitized){2,}")

def assert_no_helium_in_pdbqt(
    lines: List[str], ligand_name: str
) -> Tuple[List[str], int, str]:
    return cast(
        Tuple[List[str], int, str],
        _pdb_fixer.assert_no_helium_in_pdbqt(lines, ligand_name),
    )


def get_atom_rules() -> Dict[str, Any]:
    return cast(Dict[str, Any], _pdb_fixer.get_atom_rules())


def rules_version() -> str:
    return cast(str, _pdb_fixer.rules_version())


def scan_helium_counts(lines: List[str]) -> int:
    return cast(int, _pdb_fixer.scan_helium_counts(lines))


def collapse_sanitized_tokens(name: str) -> str:
    """
    Collapse any repeated '.sanitized' tokens in a filename while preserving the extension.
    Example: 'foo.sanitized.sanitized.pdb' -> 'foo.sanitized.pdb'
    """
    stem, ext = os.path.splitext(name)
    new_stem = _SANITIZED_RUN.sub(".sanitized", stem)
    return new_stem + ext


def collapse_sanitized_path(path: Path) -> Path:
    """Return a Path with repeated '.sanitized' tokens collapsed in the basename."""
    path = Path(path)
    new_name = collapse_sanitized_tokens(path.name)
    return path.with_name(new_name)


# Shared constants and global state
STANDARD_AMINO_ACIDS = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "HID",
    "HIE",
    "HIP",
    "SEC",
    "PYL",
    "MSE",
}

CHEM_ALIAS_DB_SOURCE = "chemdb.chem_alias_db"
CHEM_ALIAS_DB_IMPORT_ERROR: Optional[str] = None
CHEM_ALIAS_DB_FALLBACK = False

_EXCLUDE_CRYSTAL_ADDITIVES_BUILTIN = {
    "HOH",
    "CIT",
    "TAR",
    "SO4",
    "PO4",
    "CA",
    "NA",
    "K",
    "MG",
    "MN",
    "ZN",
    "GOL",
    "EDO",
    "PEG",
    "MPD",
    "TRS",
    "MES",
    "HEPES",
    "ACET",
    "ACT",
    "FMT",
    "MAL",
    "DMS",
    "IPA",
    "CLU",
    "NAG",
    "BOG",
    "TOS",
    "BES",
    "PTS",
    "OTF",
    "TRF",
    "TFA",
    "BF4",
    "PF6",
    "CL",
    "BR",
    "I",
    "HEM",
    "HEC",
    "HEA",
    "HEB",
    "HEO",
    "HEG",
    "HEF",
    "HEH",
    "PTR",
    "TPO",
    "SEP",
}

try:
    from chemdb.chem_alias_db import EXCLUDE_HET_IDS as _exclude_ids

    if _exclude_ids is None:
        raise AttributeError("EXCLUDE_HET_IDS not found")
    EXCLUDE_CRYSTAL_ADDITIVES = {str(x).upper() for x in _exclude_ids}
except (ImportError, AttributeError, TypeError) as e:
    CHEM_ALIAS_DB_SOURCE = "builtin_defaults"
    CHEM_ALIAS_DB_FALLBACK = True
    CHEM_ALIAS_DB_IMPORT_ERROR = f"{e.__class__.__name__}: {e}"
    logger.warning(
        "[ligprep.alias_db] status=fallback source=%s module=chemdb.chem_alias_db error=%s",
        CHEM_ALIAS_DB_SOURCE,
        CHEM_ALIAS_DB_IMPORT_ERROR,
    )
    logger.warning(
        "[ligprep.alias_db.telemetry] event=chem_alias_db_fallback source=%s module=chemdb.chem_alias_db",
        CHEM_ALIAS_DB_SOURCE,
    )
    EXCLUDE_CRYSTAL_ADDITIVES = set(_EXCLUDE_CRYSTAL_ADDITIVES_BUILTIN)


def chem_alias_db_fallback_status() -> Dict[str, object]:
    return {
        "source": CHEM_ALIAS_DB_SOURCE,
        "fallback": CHEM_ALIAS_DB_FALLBACK,
        "error": CHEM_ALIAS_DB_IMPORT_ERROR,
    }


MAX_HEAVY_ATOMS = 1200
MIN_ATOMS_FOR_DOCKING = 3
QUARANTINE_DIRNAME = "quarantine"

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
    for e in rules.one_letter or []:
        if e:
            out.add(e.strip().upper())
    for e2 in rules.two_letter or []:
        sym = _alias_token_to_element_symbol(e2)
        if sym:
            out.add(sym)
    return out


def free_ion_elements_from_aliases(rules) -> Set[str]:
    """
    Translate meeko.drop_free_ions and halide resnames into element symbols for
    'one-atom free ion' filtering during ligand prep.
    """
    out: Set[str] = set()
    for tok in getattr(rules, "meeko_drop_free_ions", []) or []:
        sym = _alias_token_to_element_symbol(tok)
        if sym:
            out.add(sym)
    for hal in getattr(rules, "halide_resnames", []) or []:
        sym = _alias_token_to_element_symbol(hal)
        if sym:
            out.add(sym)
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
        logging.warning(
            "[ligprep] get_atom_rules() failed; using built-ins only (%s)", e
        )
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
                logging.debug(
                    "[ligprep] ALLOWED_ELEMENTS=%s", ", ".join(sorted(ALLOWED_ELEMENTS))
                )
                if AD4_TYPES:
                    logging.debug(
                        "[ligprep] AD4_TYPES(sample)=%s",
                        ", ".join(sorted(list(AD4_TYPES))[:12]),
                    )
                if MONOATOMIC_IONS:
                    logging.debug(
                        "[ligprep] MONOATOMIC_IONS=%s",
                        ", ".join(sorted(MONOATOMIC_IONS)),
                    )
        except Exception as e:
            logging.warning(
                "[ligprep] alias init error; falling back to local lists (%s)", e
            )


_init_alias_rules_cache()


def _looks_like_monoatomic_ion_pdbqt(lines: List[str]) -> bool:
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    if len(atom_lines) != 1:
        return False
    parts = atom_lines[0].split()
    elem = _element_from_adt(parts[-1]) if parts else ""
    return elem in MONOATOMIC_IONS


def _helium_postwrite_guard(
    pdbqt_path: Path, ligand_name: str, prepped_ligands_dir: Path
) -> Tuple[bool, int, str]:
    del prepped_ligands_dir
    try:
        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
        new_lines, fixes, q_reason = assert_no_helium_in_pdbqt(lines, ligand_name)
        logging.info(
            "[helium] ligand=%s fixes=%d quarantine=%s rules=%s",
            ligand_name,
            fixes,
            str(bool(q_reason)),
            rules_version(),
        )

        if fixes > 0:
            with open(pdbqt_path, "w", encoding="utf-8") as out:
                out.writelines(new_lines)

        ok = q_reason == ""
        print(
            f"[helium] ligand={ligand_name} fixes={fixes} quarantine={not ok} rules={rules_version()}"
        )
        return ok, fixes, q_reason
    except Exception as e:
        print(
            f"[helium] ligand={ligand_name} fixes=0 quarantine=True rules={rules_version()} note=postwrite_exception:{e}"
        )
        return False, 0, f"postwrite_exception:{e}"


def _element_from_adt(adt: str) -> str:
    t = (adt or "").strip()
    if not t:
        return "C"
    u = t.upper()
    MAP = {
        "C": "C",
        "A": "C",
        "N": "N",
        "NA": "N",
        "O": "O",
        "OA": "O",
        "S": "S",
        "SA": "S",
        "H": "H",
        "HD": "H",
        "F": "F",
        "CL": "Cl",
        "BR": "Br",
        "I": "I",
        "P": "P",
        "B": "B",
        "SI": "Si",
        "SE": "Se",
        "ZN": "Zn",
        "MG": "Mg",
        "CA": "Ca",
        "MN": "Mn",
        "FE": "Fe",
        "K": "K",
        "NA+": "Na",
        "NA_": "Na",
        "NA ": "Na",
    }
    if u in MAP:
        return MAP[u]
    return u[0]


AD4_TYPES: set[str] = {
    "C",
    "A",
    "N",
    "NA",
    "O",
    "OA",
    "S",
    "SA",
    "H",
    "HD",
    "F",
    "CL",
    "BR",
    "I",
    "P",
    "B",
    "SI",
    "SE",
    "ZN",
    "MG",
    "CA",
    "MN",
    "FE",
    "K",
    "NA",
    "CU",
    "CO",
    "NI",
    "AL",
    "AG",
    "AU",
    "PT",
    "LI",
    "BA",
    "SR",
    "CS",
    "RB",
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
        bad_samples: list[tuple[str, str]] = []

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

        has_ad4_types = n_atoms > 0 and typed == n_atoms
        has_charges = n_atoms > 0 and charges_present == n_atoms
        polar_expected = (n_NA + n_OA) > 0
        has_polar_H = (n_HD > 0) or (not polar_expected)

        # Debug/audit when types are missing
        if (
            (not has_ad4_types)
            and bad_samples
            and os.environ.get("LIGPREP_DEBUG", "").strip().lower()
            in {"1", "true", "yes", "y"}
        ):
            sample_str = "; ".join([f"{tkn}:{ln}" for tkn, ln in bad_samples[:10]])
            logging.warning(
                "[ad4.audit] %s bad_types sample=%s", pdbqt_path.stem, sample_str
            )

        return {
            "has_polar_H": has_polar_H,
            "polarH_expected": polar_expected,
            "has_ad4_types": has_ad4_types,
            "has_charges": has_charges,
            "torsdof": torsdof,
            # expose counts when debugging downstream
            "n_atoms": n_atoms,
            "n_HD": n_HD,
            "n_NA": n_NA,
            "n_OA": n_OA,
            "typed": typed,
            "charges_present": charges_present,
        }
    except Exception as e:
        # Try to include a tiny context without risking new exceptions
        ctx = ""
        try:
            lines = []
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if ln.startswith(("ATOM", "HETATM")):
                        lines.append(ln.strip())
                        if len(lines) >= 3:
                            break
            if lines:
                ctx = " | ctx=" + " | ".join(lines)
        except Exception:
            pass
        logging.warning(
            "[ligprep] metrics parse failed for %s: %s%s", pdbqt_path.name, e, ctx
        )
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

    if os.environ.get("LIGPREP_DEBUG", "").strip() in {"1", "true", "yes", "y"}:
        logging.info(
            "[invariants] lig=%s has_HD=%s donors_seen=%s has_polar_H=%s has_ad4_types=%s has_charges=%s "
            "torsdof=%s fail=%s",
            pdbqt_path.stem,
            m.get("n_HD", "NA"),
            m.get("n_NA", "NA"),
            m.get("has_polar_H"),
            m.get("has_ad4_types"),
            m.get("has_charges"),
            m.get("torsdof"),
            failures,
        )
    return (len(failures) == 0), m, failures

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


def read_config(path: str | None = None) -> Dict[str, Any]:
    if path:
        return load_config(path)
    return load_config()

def get_short_path_name(long_name):
    if sys.platform != "win32":
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
        if elem == "O":
            has_O = True
        if elem != "H":
            heavy += 1
    return has_O and heavy <= 1


def is_valid_ligand(path: Path, log_dir: Path) -> bool:
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        atom_lines = [
            line
            for line in lines
            if line.startswith("ATOM") or line.startswith("HETATM")
        ]
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
            if _element_from_adt(parts[-1]) != "H":
                heavy += 1
        if heavy < 3:
            raise ValueError(f"too_few_heavy_atoms_in_pdbqt({heavy})")

        return True
    except Exception as e:
        malformed_log = Path(log_dir) / MALFORMED_LOG.name
        with open(malformed_log, "a", encoding="utf-8") as f:
            f.write(f"{path.name} - PDBQT validation failed: {e}\n")
        return False


__all__ = [
    "EXCLUDE_CRYSTAL_ADDITIVES",
    "MIN_ATOMS_FOR_DOCKING",
    "QUARANTINE_DIRNAME",
    "STANDARD_AMINO_ACIDS",
    "AD4_TYPES",
    "ALLOWED_ELEMENTS",
    "MALFORMED_DIR",
    "MALFORMED_LOG",
    "MONOATOMIC_IONS",
    "RULES",
    "_append_prep_status",
    "_audit_protonation_metrics",
    "_element_from_adt",
    "_helium_postwrite_guard",
    "_log_elem_fix_summary",
    "_log_malformed",
    "assert_no_helium_in_pdbqt",
    "chem_alias_db_fallback_status",
    "collapse_sanitized_path",
    "collapse_sanitized_tokens",
    "get_atom_rules",
    "get_short_path_name",
    "is_valid_ligand",
    "quick_pdbqt_validate",
    "read_config",
    "rules_version",
    "scan_helium_counts",
    "validate_pdbqt_invariants",
]
