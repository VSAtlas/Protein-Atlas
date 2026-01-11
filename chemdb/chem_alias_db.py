from __future__ import annotations
import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Set, Optional

try:
    import yaml  # PyYAML
except Exception:  # pragma: no cover
    yaml = None

# ---------------- path resolution ----------------

_ENV = os.environ.get("CHEM_ALIASES_YAML") or os.environ.get("CHEMDB_ALIASES")


def _find_yaml() -> Optional[Path]:
    # 1) explicit env var
    if _ENV:
        p = Path(_ENV)
        if p.is_file():
            return p
    # 2) CWD
    p = Path("aliases.yaml")
    if p.is_file():
        return p
    # 3) alongside this module
    p = Path(__file__).with_name("aliases.yaml")
    return p if p.is_file() else None


# ---------------- normalization ------------------

_norm_rx = re.compile(r"[^0-9a-z]+")


def _norm(s: str) -> str:
    if not s:
        return ""
    return _norm_rx.sub("", s.strip().lower())


# ---------------- loaders -----------------------


def _load_yaml_blocks() -> Dict[str, object]:
    """
    Returns a dict with keys we care about; empty/defaults if YAML missing.
    Expected keys in YAML:
      - alias_map: {HET or name: [synonyms...]}
      - canonical: {raw_name_or_code: canonical_name}
      - exclude_het_ids: [HET,...]
      - exclude_het_name_keywords: [strings...]
      - per_pdb_hints: {PDBID: {...}}
      - hard_fda_control_by_pdb: {PDBID: "name"}
    """
    p = _find_yaml()
    if not p or not p.is_file() or yaml is None:
        if yaml is None:
            logging.warning("[chem-alias] PyYAML not available; using built-ins only")
        else:
            logging.warning("[chem-alias] aliases.yaml not found; using built-ins only")
        return {}

    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception as e:
        logging.warning("[chem-alias] failed to read %s: %s", p, e)
        return {}


def _build_chemcomp_alias(
    alias_map: Dict[str, List[str]] | None, canonical: Dict[str, str] | None
) -> Dict[str, List[str]]:
    """
    Build a mapping like {'STI': ['gleevec','imatinib','sti571',...], ...}
    Keys are kept UPPER for HET codes if they look like 2-4 letters; values are lowercase.
    Also folds canonical remaps into the same synonym lists.
    """
    out: Dict[str, Set[str]] = {}

    alias_map = alias_map or {}
    for raw_key, syns in alias_map.items():
        key = raw_key.strip()
        key_u = key.upper()
        bucket = out.setdefault(key_u, set())
        if isinstance(syns, (list, tuple, set)):
            for s in syns:
                n = _norm(str(s))
                if n:
                    bucket.add(n)
        elif isinstance(syns, str):
            n = _norm(syns)
            if n:
                bucket.add(n)
        # include a normalized self-name if it looks like a name (not only a code)
        nself = _norm(key)
        if nself:
            bucket.add(nself)

    canonical = canonical or {}
    for raw, can in canonical.items():
        r_u = str(raw).strip().upper()
        c_norm = _norm(str(can))
        if not r_u or not c_norm:
            continue
        bucket = out.setdefault(r_u, set())
        bucket.add(c_norm)

    # finalize as lists (sorted for repeatability)
    return {k: sorted(v) for k, v in out.items() if v}


# ---------------- built-in fallback --------------

_BUILTIN_CHEMCOMP_ALIAS: Dict[str, List[str]] = {
    # just enough to keep the benchmark sane when YAML is missing
    "STI": ["imatinib", "gleevec", "sti571", "sti"],
    "NIL": ["nilotinib"],
    "OHT": ["4hydroxytamoxifen", "oht", "hydroxytamoxifen"],
    "BAX": ["bazedoxifene", "bax"],
    "032": ["crizotinib", "pf02341066", "032"],
    "LEV": ["levonorgestrel", "lev"],
    "VGH": ["dasatinib", "sprycel", "vgh"],
    "MI1": ["mi1"],  # placeholder
    "0LI": ["lapatinib", "0li"],  # example from your table
    "P30": ["palbociclib", "p30"],
    "PTR": ["pertuzumab", "ptr"],  # placeholder
    # extend here as needed
}
_BUILTIN_EXCLUDE_HET_IDS: Set[str] = {
    "HOH",
    "WAT",
    "DOD",
    "NA",
    "K",
    "CL",
    "BR",
    "I",
    "CA",
    "MG",
    "MN",
    "ZN",
    "FE",
    "CU",
    "CO",
    "NI",
    "MO",
}
_BUILTIN_EXCLUDE_HET_NAME_KEYWORDS: Set[str] = {
    # uppercase substrings checked against names
    "WATER",
    "SOLVENT",
    "BUFFER",
    "GLYCEROL",
    "TRIS",
    "HEPES",
    "MES",
    "PEG",
    "ETHYLENE",
    "DMSO",
}
_BUILTIN_PER_PDB_HINTS: Dict[str, dict] = {}
_BUILTIN_HARD_FDA_BY_PDB: Dict[str, str] = {}

# ---------------- module globals (API) ----------


def _init():
    data = _load_yaml_blocks()
    alias_map = data.get("alias_map") or {}
    canonical = data.get("canonical") or data.get("canon") or {}
    chem_alias = _build_chemcomp_alias(alias_map, canonical)

    exclude_ids = data.get("exclude_het_ids") or []
    exclude_kw = data.get("exclude_het_name_keywords") or []
    per_pdb = data.get("per_pdb_hints") or {}
    hard = data.get("hard_fda_control_by_pdb") or {}

    # Exported names (keep EXACT identifiers expected by external code)
    globals().update(
        {
            "CHEMCOMP_ALIAS": chem_alias or _BUILTIN_CHEMCOMP_ALIAS,
            "EXCLUDE_HET_IDS": set(map(str.upper, exclude_ids))
            or _BUILTIN_EXCLUDE_HET_IDS,
            "EXCLUDE_HET_NAME_KEYWORDS": set(map(str.upper, exclude_kw))
            or _BUILTIN_EXCLUDE_HET_NAME_KEYWORDS,
            "PER_PDB_HINTS": per_pdb or _BUILTIN_PER_PDB_HINTS,
            "HARD_FDA_CONTROL_BY_PDB": hard or _BUILTIN_HARD_FDA_BY_PDB,
        }
    )


_init()

# ---------------- convenience helpers -----------


def resolve_name_tokens(raw: str) -> List[str]:
    """
    Given any text (HET code, brand, generic), return possible normalized tokens
    that can be matched against CHEMCOMP_ALIAS values.
    """
    n = _norm(raw)
    if not n:
        return []
    return [n]


def alias_list_for_het(het_code: str) -> List[str]:
    """Return the synonym list for a three/four-letter HET code (lowercased tokens)."""
    return CHEMCOMP_ALIAS.get(str(het_code).upper(), [])


# Keep a minimal CLI for sanity checks
if __name__ == "__main__":
    import sys

    if len(sys.argv) == 1:
        print(
            "keys:",
            ", ".join(
                [
                    "CHEMCOMP_ALIAS",
                    "EXCLUDE_HET_IDS",
                    "EXCLUDE_HET_NAME_KEYWORDS",
                    "PER_PDB_HINTS",
                    "HARD_FDA_CONTROL_BY_PDB",
                ]
            ),
        )
        print("CHEMCOMP_ALIAS entries:", len(CHEMCOMP_ALIAS))
        print("example STI ->", alias_list_for_het("STI"))
        sys.exit(0)
    for arg in sys.argv[1:]:
        print(arg, "=>", alias_list_for_het(arg))
