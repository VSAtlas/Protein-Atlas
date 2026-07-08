# metabolite_resolver.py
from __future__ import annotations
import csv
import logging
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Mapping, Sequence

from config.runtime_config import load_config, validate_config

Chem: Any
Descriptors: Any
MurckoScaffold: Any
AllChem: Any
rdMolDescriptors: Any
try:
    from rdkit import Chem as _Chem
    from rdkit.Chem import rdMolDescriptors as _Descriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold as _MurckoScaffold
    from rdkit.Chem import AllChem as _AllChem, rdMolDescriptors as _rdMolDescriptors
except Exception:  # RDKit optional; name-based still works
    Chem = MurckoScaffold = AllChem = rdMolDescriptors = Descriptors = None
else:
    Chem = _Chem
    Descriptors = _Descriptors
    MurckoScaffold = _MurckoScaffold
    AllChem = _AllChem
    rdMolDescriptors = _rdMolDescriptors

log = logging.getLogger("metabolite_resolver")

_cfg = load_config()
validate_config(_cfg)

# optional chemdb data: prefer external alias/rule tables if available
if TYPE_CHECKING:
    _EXT_PARENT_ALIASES: Optional[Mapping[str, str]] = None
    _EXT_NAME_RULES: Optional[Sequence[Tuple[Any, str]]] = None
else:
    try:
        import chemdb.chem_alias_db as _chem_alias_db
    except Exception:
        _chem_alias_db = None
    _EXT_PARENT_ALIASES = (
        getattr(_chem_alias_db, "PARENT_ALIASES", None) if _chem_alias_db else None
    )
    _EXT_NAME_RULES = (
        getattr(_chem_alias_db, "NAME_RULES", None) if _chem_alias_db else None
    )

# RDKit fingerprint params (keep current defaults)
_METABO_FP_RADIUS = int(_cfg.get("METABO_FINGERPRINT_RADIUS", 2))
_METABO_FP_BITS = int(_cfg.get("METABO_FINGERPRINT_BITS", 2048))
_USE_SCAFFOLD = str(_cfg.get("METABO_USE_SCAFFOLD_MATCH", "true")).lower() in (
    "1",
    "true",
    "yes",
    "on",
)

# --------------------------
# Normalization helpers
# --------------------------
_WS = re.compile(r"\s+")
_PUNC = re.compile(r"[-_/\\,;:\|\[\]\(\)\{\}\.\+\*'\"]+")
_RDK_RE = re.compile(r"rdk[_-]?(\d+)", re.IGNORECASE)
_BAD_SYNONYM_PREFIXES = (
    "CHEMBL",
    "CID",
    "UNII",
    "SCHEMBL",
    "CHEBI",
    "ZINC",
    "PUBCHEM",
    "NCI",
    "NSC",
)
_NAME_FIELDS_PRIMARY = [
    "rxnorm_generic_name",
    "drugcentral_generic_name",
    "generic_name",
    "display_name",
    "pubchem_record_title",
    "pubchem_name",
    "pubchem_iupac_name",
    "remark_name",
    "sdf_title",
]
_NAME_FIELDS_FALLBACK = [
    "fda_name",
    "name",
    "drug_name",
    "preferred_name",
    "international_nonproprietary_name",
    "inn",
    "brand_name",
    "brand_names",
    "trade_name",
    "label_name",
]
_SYNONYM_FIELDS_EXTRA = [
    "synonyms",
    "alias",
    "alts",
    "pubchem_synonyms",
    "display_name",
    "brand_names",
    "rxnorm_brand_names",
    "drugcentral_brand_names",
    "generic_name",
    "rxnorm_generic_name",
    "drugcentral_generic_name",
    "pubchem_name",
    "pubchem_record_title",
    "pubchem_iupac_name",
    "remark_name",
    "sdf_title",
    "fda_name",
    "name",
    "drug_name",
    "preferred_name",
    "international_nonproprietary_name",
    "inn",
    "brand_name",
    "trade_name",
    "label_name",
]
_NAME_FALLBACK_EXCLUDE = ("filename", "file_name", "pathname", "path_name")


def _norm(s: Optional[str]) -> str:
    if not s:
        return ""
    s = str(s).strip().lower().replace("\u00a0", " ")
    s = _WS.sub(" ", s)
    s = _PUNC.sub(" ", s)
    return s


def _tokenize(s: str) -> List[str]:
    return [t for t in re.split(r"[^a-z0-9\+]+", (s or "").lower()) if t]


def _extract_rdk_id(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    match = _RDK_RE.search(str(value))
    if not match:
        return None
    return f"rdk_{match.group(1).zfill(7)}"


def _row_value(row: Dict[str, Any], key: str) -> str:
    value = row.get(key)
    if value is None:
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _pick_first(row: Dict[str, Any], *keys: str) -> str:
    for key in keys:
        value = _row_value(row, key)
        if value:
            return value
    return ""


def _collect_values(row: Dict[str, Any], *keys: str) -> List[str]:
    values = []
    for key in keys:
        value = _row_value(row, key)
        if value:
            values.append(value)
    return values


def _digit_ratio(text: str) -> float:
    if not text:
        return 1.0
    digits = sum(ch.isdigit() for ch in text)
    return digits / len(text)


def _is_registry_like(text: str) -> bool:
    if not text:
        return False
    upper = text.strip().upper()
    for prefix in _BAD_SYNONYM_PREFIXES:
        if upper.startswith(prefix):
            return True
    if re.match(r"^\d{2,7}-\d{2}-\d$", text.strip()):
        return True
    return False


def _is_iupac_like(text: str) -> bool:
    value = str(text or "").strip()
    if not value:
        return False
    length = len(value)
    if length < 25:
        return False
    digits = sum(ch.isdigit() for ch in value)
    punct = sum(1 for ch in value if not ch.isalnum() and not ch.isspace())
    if length >= 50:
        return True
    if length >= 35 and (punct / length) >= 0.12:
        return True
    if length >= 35 and (digits / length) >= 0.2:
        return True
    if length >= 30 and punct >= 6 and (digits / length) >= 0.1:
        return True
    return False


def _choose_pubchem_synonym(raw: str) -> str:
    if not raw:
        return ""
    for part in raw.split(";"):
        candidate = part.strip()
        if not candidate:
            continue
        if len(candidate) > 32:
            continue
        if not re.search(r"[A-Za-z]", candidate):
            continue
        if _digit_ratio(candidate) > 0.4:
            continue
        if _is_registry_like(candidate):
            continue
        if _is_iupac_like(candidate):
            continue
        return candidate
    return ""


def _canonical_fieldname(name: str) -> str:
    return str(name or "").strip().lower()


def _infer_name_fields(fieldnames: List[str]) -> List[str]:
    field_map = {_canonical_fieldname(name): name for name in fieldnames if name}
    ordered: List[str] = []

    for name in _NAME_FIELDS_PRIMARY + _NAME_FIELDS_FALLBACK:
        key = _canonical_fieldname(name)
        actual = field_map.get(key)
        if actual and actual not in ordered:
            ordered.append(actual)

    for name in fieldnames:
        key = _canonical_fieldname(name)
        if not key or name in ordered:
            continue
        if "name" in key and not any(bad in key for bad in _NAME_FALLBACK_EXCLUDE):
            ordered.append(name)

    return ordered


def _infer_synonym_fields(fieldnames: List[str], name_fields: List[str]) -> List[str]:
    field_map = {_canonical_fieldname(name): name for name in fieldnames if name}
    ordered: List[str] = []

    for name in _SYNONYM_FIELDS_EXTRA:
        key = _canonical_fieldname(name)
        actual = field_map.get(key)
        if actual and actual not in ordered:
            ordered.append(actual)

    for name in name_fields:
        if name not in ordered:
            ordered.append(name)

    return ordered


def _select_preferred_name(
    row: Dict[str, Any], rdk_id: str, name_fields: List[str]
) -> str:
    name = _pick_first(row, *name_fields)
    if not name:
        return rdk_id or "unknown"

    if _is_iupac_like(name):
        for field in name_fields:
            candidate = _row_value(row, field)
            if candidate and not _is_iupac_like(candidate):
                return candidate
        synonym = _choose_pubchem_synonym(_pick_first(row, "pubchem_synonyms"))
        if synonym:
            return synonym

    return name


# --------------------------
# Library structures
# --------------------------
@dataclass
class DrugRec:
    rdk_id: str  # e.g. rdk_0000229
    name: str  # canonical preferred name
    synonyms: List[str]  # any synonyms
    smiles: Optional[str]  # SMILES if available
    inchikey: Optional[str]  # InChIKey if available


class LibraryIndex:
    def __init__(self):
        self.id_to_rec: Dict[str, DrugRec] = {}
        self.name_index: Dict[str, List[str]] = {}  # norm(name/syn) -> [rdk_id]
        self.scaffold_index: Dict[str, List[str]] = {}  # murcko SMILES -> [rdk_id]
        self.firstblock_index: Dict[
            str, List[str]
        ] = {}  # InChIKey first block -> [rdk_id]

    def add(self, rec: DrugRec):
        self.id_to_rec[rec.rdk_id] = rec
        for key in {_norm(rec.name), *(_norm(s) for s in rec.synonyms)}:
            if key:
                self.name_index.setdefault(key, []).append(rec.rdk_id)
        if rec.inchikey:
            fb = rec.inchikey.split("-")[0]
            self.firstblock_index.setdefault(fb, []).append(rec.rdk_id)
        if Chem and rec.smiles:
            try:
                m = Chem.MolFromSmiles(rec.smiles)
                if m:
                    scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
                    if scaf:
                        self.scaffold_index.setdefault(scaf, []).append(rec.rdk_id)
            except Exception:
                pass


def load_library_index(mapping_csv: str) -> LibraryIndex:
    """
    Flexible reader for your FDA mapping CSV (column names are best-effort).
    Expected columns (any casing): rdk_id, drug/name, synonyms, smiles, inchikey.
    'synonyms' can be pipe- or comma-separated.
    """
    idx = LibraryIndex()
    with open(mapping_csv, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        fieldnames = [name for name in (r.fieldnames or []) if name]
        name_fields = _infer_name_fields(fieldnames)
        synonym_fields = _infer_synonym_fields(fieldnames, name_fields)
        for row in r:
            scheme = _pick_first(row, "scheme")
            rdk_id = _extract_rdk_id(
                _pick_first(row, "rdk_id", "rdkid", "ligand_id")
            )
            if not rdk_id:
                rdk_id = _extract_rdk_id(
                    _pick_first(row, "path", "sdf_title", "remark_name")
                )
            if not rdk_id and str(scheme).strip().lower() == "rdk":
                file_num = _pick_first(row, "file_num")
                if file_num:
                    try:
                        rdk_id = f"rdk_{int(float(str(file_num))):07d}"
                    except ValueError:
                        rdk_id = None
            if not rdk_id:
                raw_id = _pick_first(row, "id")
                rdk_id = _extract_rdk_id(raw_id) or (raw_id or "")

            name = _select_preferred_name(row, rdk_id, name_fields)
            syn_fields = []
            syn_fields.extend(_collect_values(row, *synonym_fields))
            path_val = _pick_first(row, "path")
            if path_val:
                base = os.path.basename(str(path_val))
                stem = os.path.splitext(base)[0]
                if stem:
                    syn_fields.append(stem)
            syns_raw = "|".join(str(s) for s in syn_fields if s)
            smiles = _pick_first(row, "smiles", "smile", "smiles_rdkit")
            inchikey = _pick_first(row, "inchikey", "inchi_key", "ikey")
            synonyms = [s.strip() for s in re.split(r"[|,;]", syns_raw) if s.strip()]
            if rdk_id:
                idx.add(DrugRec(rdk_id, str(name), synonyms, smiles, inchikey))
    log.info(
        "Loaded FDA library: %d entries (names=%d, scaffolds=%d, inchikeyFB=%d)",
        len(idx.id_to_rec),
        len(idx.name_index),
        len(idx.scaffold_index),
        len(idx.firstblock_index),
    )
    return idx


# --------------------------
# Metabolite name heuristics
# --------------------------
# Given a metabolite string, propose plausible parent names.
# Simple and safe rules first; you can extend this list freely.
# prefer external rules/aliases if present; otherwise keep in-file defaults
_NAME_RULES = _EXT_NAME_RULES or [
    (re.compile(r"\b4[- ]?hydroxy(tamoxifen)\b", re.I), r"\1"),
    (re.compile(r"\b17[- ]?des?acetyl(\w+)\b", re.I), r"\1"),
    (re.compile(r"\bN[- ]?des?methyl(\w+)\b", re.I), r"\1"),
    (re.compile(r"\bdes?methyl[- ]?(\w+)\b", re.I), r"\1"),
    (re.compile(r"\bO[- ]?des?methyl(\w+)\b", re.I), r"\1"),
    (re.compile(r"\bN[- ]?oxide\b", re.I), r""),
    (re.compile(r"\b[- ]?glucuronide\b", re.I), r""),
    (re.compile(r"\b[- ]?sulfate\b", re.I), r""),
    (re.compile(r"\b[- ]?phosphate\b", re.I), r""),
    (re.compile(r"\bcarboxylic acid\b", re.I), r""),
]

_KNOWN_PARENT_ALIASES = _EXT_PARENT_ALIASES or {
    "oht": "tamoxifen",
    "4-hydroxytamoxifen": "tamoxifen",
    "bax": "baxitinib",
    "lev": "levofloxacin",
    "axi": "axitinib",
    "nil": "nilotinib",
    "p31": "prostaglandin e1",
    "pg6": "prostaglandin g2",
    "pge": "prostaglandin e2",
}


def _propose_parent_names(metabolite_like: str) -> List[str]:
    base = metabolite_like or ""
    s = base
    out = set()
    s_norm = _norm(s)
    toks = _tokenize(s)
    # direct alias by 3-letter code etc.
    if s_norm in _KNOWN_PARENT_ALIASES:
        out.add(_KNOWN_PARENT_ALIASES[s_norm])
    for t in toks:
        if t in _KNOWN_PARENT_ALIASES:
            out.add(_KNOWN_PARENT_ALIASES[t])
    # regex rules
    s2 = s
    for pat, rep in _NAME_RULES:
        s2 = pat.sub(rep, s2)
    s2 = _norm(s2)
    if s2 and s2 != s_norm:
        out.add(s2)
    # also strip common suffixes
    s3 = re.sub(
        r"\b(metabolite|acid|sulfate|glucuronide|phosphate|oxide)\b", "", s, flags=re.I
    )
    s3 = _norm(s3)
    if s3 and s3 != s_norm:
        out.add(s3)
    return [o for o in out if o]


# --------------------------
# RDKit helpers (optional)
# --------------------------
def _largest_fragment(m: Chem.Mol) -> Chem.Mol:
    frags = list(Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False))
    frags.sort(key=lambda x: x.GetNumAtoms(), reverse=True)
    return frags[0]


def _murcko(m: Chem.Mol) -> Optional[str]:
    try:
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
        return scaf or None
    except Exception:
        return None


def _fp(m: Chem.Mol):
    return AllChem.GetMorganFingerprintAsBitVect(
        m, radius=_METABO_FP_RADIUS, nBits=_METABO_FP_BITS
    )


def _tanimoto(a, b) -> float:
    return rdMolDescriptors.TanimotoSimilarity(a, b)


# --------------------------
# Core API
# --------------------------
def ensure_parent_drugs_for_controls(
    pdb_code: str,
    ligands_raw_dir: str,
    fda_index: LibraryIndex,
    max_additions: int = int(_cfg.get("METABO_MAX_PARENT_ADDITIONS", 2)),
) -> List[str]:
    """
    Look at extracted co-crystal ligands for this PDB.
    If any looks like a metabolite, try to locate the parent drug in the FDA library.
    Returns a (deduped) list of rdk_ids to add to your whitelist.
    """
    candidates: List[Tuple[str, float]] = []  # (rdk_id, score)
    raw_files = []
    if os.path.isdir(ligands_raw_dir):
        raw_files = [
            os.path.join(ligands_raw_dir, f)
            for f in os.listdir(ligands_raw_dir)
            if f.lower().endswith((".pdb", ".sdf", ".mol"))
        ]

    # pass 1: name-based guesses from filenames
    for path in raw_files:
        base = os.path.basename(path)
        stem = os.path.splitext(base)[0]
        # typical pattern: OHT_A600.sanitized.pdb — get left token(s)
        name_guess = _norm(re.split(r"[_\.]", stem)[0])
        for proposed in _propose_parent_names(name_guess):
            key = _norm(proposed)
            if key in fda_index.name_index:
                for rdk_id in fda_index.name_index[key]:
                    candidates.append((rdk_id, 0.70))  # heuristic score

    # pass 2: structure-based scaffold match (RDKit only, if library has SMILES)
    if (
        _USE_SCAFFOLD
        and Chem
        and any(rec.smiles for rec in fda_index.id_to_rec.values())
    ):
        for path in raw_files:
            m = None
            try:
                if path.lower().endswith(".sdf"):
                    suppl = Chem.SDMolSupplier(path, removeHs=False)
                    m = next((x for x in suppl if x), None)
                elif path.lower().endswith(".mol"):
                    m = Chem.MolFromMolFile(path, removeHs=False)
                else:
                    m = Chem.MolFromPDBFile(path, removeHs=False)
            except Exception:
                m = None
            if not m:
                continue
            try:
                m = _largest_fragment(m)
                Chem.SanitizeMol(m)
            except Exception:
                pass
            scaf = _murcko(m)
            if not scaf:
                continue
            lib_ids = fda_index.scaffold_index.get(scaf, [])
            if not lib_ids:
                continue
            try:
                qfp = _fp(m)
            except Exception:
                qfp = None
            for rid in lib_ids:
                rec = fda_index.id_to_rec[rid]
                if not rec.smiles:
                    candidates.append((rid, 0.72))
                    continue
                try:
                    mol_lib = Chem.MolFromSmiles(rec.smiles)
                    if not mol_lib:
                        continue
                    score = 0.75
                    if qfp:
                        lfp = _fp(mol_lib)
                        score = _tanimoto(qfp, lfp)
                    candidates.append((rid, float(score)))
                except Exception:
                    pass

    # rank, dedupe, cap
    picked = []
    seen = set()
    for rid, sc in sorted(candidates, key=lambda x: x[1], reverse=True):
        if rid not in seen:
            seen.add(rid)
            picked.append(rid)
        if len(picked) >= max_additions:
            break
    log.info(
        "[Metabolite→Parent] %s: adding %d parent(s): %s", pdb_code, len(picked), picked
    )
    return picked


def resolve_corresponding_name_for_rdk(rdk_id: str, fda_index: LibraryIndex) -> str:
    """Return a nice, canonical drug name for any rdk_* id (falls back to the id)."""
    rec = fda_index.id_to_rec.get(rdk_id)
    return rec.name if rec else rdk_id


def resolve_corresponding_name_from_text(
    text: str, fda_index: LibraryIndex
) -> Optional[str]:
    """Map a free-text ligand name/code to a library drug name, if possible."""
    key = _norm(text)
    if key in fda_index.name_index:
        rid = fda_index.name_index[key][0]
        return fda_index.id_to_rec[rid].name
    # try token subset
    toks = _tokenize(text)
    for t in toks:
        k = _norm(t)
        if k in fda_index.name_index:
            rid = fda_index.name_index[k][0]
            return fda_index.id_to_rec[rid].name
    # last resort: known aliases
    if key in _KNOWN_PARENT_ALIASES:
        alias = _KNOWN_PARENT_ALIASES[key]
        alias_k = _norm(alias)
        if alias_k in fda_index.name_index:
            rid = fda_index.name_index[alias_k][0]
            return fda_index.id_to_rec[rid].name
        return alias
    return None
