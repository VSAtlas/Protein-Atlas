import re
import json
import time
import logging
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any

import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.SaltRemover import SaltRemover

# ---- RDKit standardize import (robust across builds) ----
try:
    # Your env (2025.03.5) works here:
    from rdkit.Chem.MolStandardize import rdMolStandardize as _std  # type: ignore
except Exception:
    try:
        # Some builds expose it here:
        from rdkit.Chem import rdMolStandardize as _std  # type: ignore
    except Exception:
        _std = None  # we'll gracefully fall back in neutralize_smiles()
from input_and_export_functions import load_inputs, validate_config
cfg = load_inputs(); validate_config(cfg)

# Logs & CSVs (paths)
from pathlib import Path
logs_root = Path(cfg.get("LOGS_DIR", cfg.get("DOCKED_DIR", ".")))
logs_root.mkdir(parents=True, exist_ok=True)
LOG_PATH = logs_root / cfg.get("IDENTIFY_LOG_BASENAME", "fda_mapping.log")

# Where to write mapping CSVs — default to DOCKED_DIR (keeps behavior close to outputs)
out_root = Path(cfg.get("DOCKED_DIR", "."))
base_csv = out_root / cfg.get("IDENTIFY_MAPPING_CSV", "fda_mapping_from_pdbqt.csv")
bad_csv  = out_root / cfg.get("IDENTIFY_UNMAPPED_CSV", "fda_mapping_unmapped_or_no_smiles.csv")
# Try to import alias tables from chemdb; fallback to in-code defaults below
try:
    from chemdb.chem_alias_db import SALT_WORDS as _SALT, HYDRATE_WORDS as _HYDR, FORM_WORDS as _FORM, GENERIC_SUFFIX_BONUS as _GENBONUS
    SALT_WORDS = set(_SALT) if _SALT else set()
    HYDRATE_WORDS = set(_HYDR) if _HYDR else set()
    FORM_WORDS = set(_FORM) if _FORM else set()
    GENERIC_SUFFIX_BONUS = tuple(_GENBONUS) if _GENBONUS else ()
except Exception:
    # fallbacks below remain as-is
    pass

# ========= User config =========
SDF_PATH = Path(r"E:\PythonProject\protein_automation\extracted_ligands\merged_dedup.sdf")
SEARCH_ROOT = Path(r"E:\PythonProject\protein_automation\prepped_ligands")
ONE_BASED = True            # set False if your fda_000 matches SDF index 0

# Lookups (toggle per source)
LOOKUP_PUBCHEM = True
LOOKUP_RXNORM = True
LOOKUP_DRUGCENTRAL = True
LOOKUP_OPENFDA = True

# Caches
CACHE_PATH_PUBCHEM = Path("pubchem_name_cache.json")
CACHE_PATH_RXNORM = Path("rxnorm_cache.json")  # retained (used implicitly for rate limiting if extended)
CACHE_PATH_DRUGCENTRAL = Path("drugcentral_cache.json")  # legacy (not used with TSV path; kept for compatibility)
CACHE_PATH_OPENFDA = Path("openfda_cache.json")

RATE_LIMIT_SEC = 0.20       # be gentle to remote APIs
MAX_LOOKUPS = None          # e.g., 200 for trial; None = all
LOG_PATH = Path("fda_mapping.log")

# DrugCentral local TSV (download their SMILES/InChI file and point here)
DRUGCENTRAL_TSV = Path(r"E:\PythonProject\protein_automation\ref\drugcentral_structures.tsv")

# Optional: copy these SDF props if present (we’ll mine them for hints)
SDF_FIELDS = [
    "DrugBank_ID", "DRUGBANK_ID", "DRUGBANK", "DRUG_NAME", "DRUGNAME",
    "PUBCHEM_CID", "PUBCHEM", "UNII", "CAS", "ZINC_ID", "ZINC", "ID", "NAME"
]

# ---- Generic/brand naming helpers ----
SALT_WORDS = {
    "hydrochloride","hydrobromide","chloride","bromide","mesylate","maleate","succinate",
    "fumarate","tartrate","citrate","acetate","phosphate","nitrate","sulfate","sulphate",
    "tosylate","besylate","pamoate","bitartrate","lactate","oxalate","carbonate","bicarbonate",
    "sodium","potassium","calcium","magnesium","zinc","aluminum","arginate","arginine"
}
HYDRATE_WORDS = {"anhydrous","hydrate","hemihydrate","monohydrate","dihydrate","trihydrate"}

# dosage/form words to strip from display names
FORM_WORDS = {
    "tablet","tablets","capsule","capsules","solution","suspension","injection","injectable",
    "oral","ophthalmic","topical","cream","ointment","gel","xr","er","sr","dr","extended",
    "prolonged","delayed","release","chewable","lozenges","patch","transdermal","nasal",
    "spray","aerosol","inhalation","powder","granules","syrup","elixir","drops","kit","for"
}

GENERIC_SUFFIX_BONUS = (
    "ine","ide","ol","one","ate","vir","mab","nib","pril","sartan","azole","caine","mycin","cillin","dipine","vastatin"
)

# ========= requests session with retries =========
def make_session() -> requests.Session:
    sess = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=3)
    sess.mount("http://", adapter)
    sess.mount("https://", adapter)
    return sess

REQ = make_session()

def _norm_ws(s: str) -> str:
    return " ".join(s.split())

def _strip_marks(s: str) -> str:
    # drop ™, ® and similar marks for scoring
    return s.replace("™","").replace("®","").replace("©","").strip()

def _strip_form_words(name: str) -> str:
    tok = name.replace("-", " ").replace("/", " ").replace(",", " ").split()
    keep = [t for t in tok if t.lower() not in FORM_WORDS]
    return " ".join(keep)

def _strip_salt_hydrate_tail(name: str) -> str:
    # remove trailing salt/hydrate tokens
    tokens = (
        name.replace(",", " ")
            .replace("-", " ")
            .replace("/", " ")
            .lower().split()
    )
    while tokens and (tokens[-1] in SALT_WORDS or tokens[-1] in HYDRATE_WORDS):
        tokens.pop()
    return " ".join(tokens)

def _clean_display_name(s: Optional[str]) -> Optional[str]:
    if not s:
        return s
    s = _norm_ws(_strip_marks(s))
    s = _strip_salt_hydrate_tail(s)
    s = _strip_form_words(s)
    return _norm_ws(s)

def neutralize_smiles(smi: Optional[str]) -> Optional[str]:
    """Return neutral, largest-fragment SMILES if possible (with graceful fallback)."""
    if not isinstance(smi, str) or not smi.strip():
        return None
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None

        if _std is not None:
            mol = _std.LargestFragmentChooser().choose(mol)
            mol = _std.Uncharger().uncharge(mol)
        else:
            # Fallback: remove salts and keep largest fragment
            mol = SaltRemover().StripMol(mol, dontRemoveEverything=True)
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if frags:
                mol = max(frags, key=lambda m: m.GetNumAtoms())

        Chem.SanitizeMol(mol, catchErrors=True)
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None

def safe_smiles(mol) -> Optional[str]:
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True)
    except Exception:
        try:
            return Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            return None

def safe_inchikey(mol) -> Optional[str]:
    if mol is None:
        return None
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None

# ========= Logging setup =========
LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.FileHandler(LOG_PATH, encoding="utf-8"),
        # console: quiet if requested
        (logging.StreamHandler())
    ]
)
# Adjust console level dynamically
for h in logging.getLogger().handlers:
    if isinstance(h, logging.StreamHandler):
        h.setLevel(logging.WARNING if bool(cfg.get("QUIET_CONSOLE", False)) else logging.INFO)

log = logging.getLogger("fda-map")

# ========= SDF & file scanning =========
def load_sdf_rows(sdf_path: Path, one_based: bool) -> pd.DataFrame:
    log.info(f"Parsing SDF: {sdf_path}")
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    rows = []
    for i, mol in enumerate(suppl):
        idx = i + 1 if one_based else i
        if mol is None:
            rows.append({"sdf_index": idx, "sdf_title": None, "smiles": None, "inchikey": None})
            continue

        title = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        row: Dict[str, Any] = {
            "sdf_index": idx,
            "sdf_title": title,
            "smiles": safe_smiles(mol),
            "inchikey": safe_inchikey(mol),
        }

        props = set(mol.GetPropNames())
        # store any helpful fields in lowercase
        for key in SDF_FIELDS:
            if key in props:
                row[key.lower()] = mol.GetProp(key)

        # extra: attempt to compute neutralized SMILES to aid dedup & lookup
        if row.get("smiles"):
            row["smiles_neutral"] = neutralize_smiles(row["smiles"])
        else:
            row["smiles_neutral"] = None

        rows.append(row)

    df = pd.DataFrame(rows)
    log.info(f"SDF records parsed: {len(df)} "
             f"(mol=None: {int(df['sdf_title'].isna().sum())})")
    return df

def _merge_coverage(files_df: pd.DataFrame, sdf_df: pd.DataFrame, one_based: bool) -> float:
    """
    Compute how many rows match when merging file_num (adjust for indexing) to sdf_index.
    Returns coverage fraction in [0,1].
    """
    if files_df.empty or sdf_df.empty:
        return 0.0
    # Adjust indices based on hypothesis
    adj = files_df.copy()
    adj["match_index"] = adj["file_num"] + (1 if one_based else 0)
    merged = adj.merge(sdf_df[["sdf_index"]], left_on="match_index", right_on="sdf_index", how="left")
    hits = int(merged["sdf_index"].notna().sum())
    total = len(adj)
    return (hits / total) if total else 0.0

def infer_indexing(files_df: pd.DataFrame, sdf_df: pd.DataFrame, default_one_based: bool) -> bool:
    """
    Try both zero-based and one-based; pick the one with higher coverage.
    If a tie, fall back to the provided default_one_based.
    """
    # files_df must contain column 'file_num'
    tmp = files_df[["file_num"]].drop_duplicates().copy()
    cov_zero = _merge_coverage(tmp, sdf_df, one_based=False)
    cov_one  = _merge_coverage(tmp, sdf_df, one_based=True)
    log.info(f"Indexing coverage test → zero-based: {cov_zero:.3f}, one-based: {cov_one:.3f}")

    if cov_zero > cov_one:
        log.info("Auto-selected ZERO-based indexing (rdk_* style).")
        return False
    if cov_one > cov_zero:
        log.info("Auto-selected ONE-based indexing (fda_* style).")
        return True

    log.info(f"Coverage tie; using configured ONE_BASED={default_one_based}.")
    return default_one_based

def sanity_check_index_merge(mapping: pd.DataFrame) -> None:
    total = len(mapping)
    missing = int(mapping["sdf_title"].isna().sum())
    pct = 0 if total == 0 else 100.0 * (total - missing) / total
    log.info(f"Merge coverage: {total - missing}/{total} ({pct:.1f}%)")
    if pct < 80:
        log.warning("Low coverage. Check that your merged_dedup.sdf ordering aligns with prepped ligand numbering.")

# ---------- cache utils ----------
def load_cache(path: Path) -> dict:
    if path.exists():
        try:
            return json.load(open(path, "r", encoding="utf-8"))
        except Exception:
            log.warning(f"Cache exists but failed to load ({path}); starting fresh.")
    return {}

def save_cache(cache: dict, path: Path) -> None:
    try:
        json.dump(cache, open(path, "w", encoding="utf-8"))
    except Exception as e:
        log.error(f"Failed to save cache to {path}: {e}")

def _calc_formula_heavy_unpack(args):  # new
    """
    Unpack (path, smiles, ...) and forward to _calc_formula_heavy.
    Exists only to avoid lambdas so ProcessPool can pickle it on Windows.
    """
    return _calc_formula_heavy(*args)

# ===================== PDBQT scanning & identity ========================
def scan_prepped_files(root: Path) -> pd.DataFrame:
    """
    Find prepped ligand files under `root`. Supports:
      • rdk_000123.pdbqt
      • rdk_000123_stage5.pdbqt
      • fda_1158.pdbqt
      • fda_1158_stage3.pdbqt

    Returns a dataframe with:
      scheme: 'rdk' or 'fda'
      file_num: int
      stage: optional int (None if not present)
      path: full file path
    """
    log.info(f"Scanning for prepped ligands under: {root}")

    rx_rdk = re.compile(r"(?i)\brdk_(\d+)(?:_stage(\d+))?\.pdbqt$")
    rx_fda = re.compile(r"(?i)\bfda_(\d+)(?:_stage(\d+))?\.pdbqt$")
    recs: List[Dict[str, Any]] = []

    total = 0
    for p in root.rglob("*.pdbqt"):
        if not p.is_file():
            continue
        name = p.name
        m = rx_rdk.search(name)
        if m:
            total += 1
            recs.append({
                "scheme": "rdk",
                "file_num": int(m.group(1)),
                "stage": int(m.group(2)) if m.group(2) is not None else None,
                "path": str(p)
            })
            continue
        m = rx_fda.search(name)
        if m:
            total += 1
            recs.append({
                "scheme": "fda",
                "file_num": int(m.group(1)),
                "stage": int(m.group(2)) if m.group(2) is not None else None,
                "path": str(p)
            })

    df = pd.DataFrame(recs)
    if df.empty:
        log.warning("No prepped ligand files found. Check SEARCH_ROOT and filename patterns.")
        return df

    # 🔧 make stage a numeric nullable int to avoid object-dtype fillna warnings later
    df["stage"] = pd.to_numeric(df["stage"], errors="coerce").astype("Int16")

    log.info(
        f"Found {total} .pdbqt files; unique (scheme,file_num) pairs: "
        f"{df.drop_duplicates(['scheme', 'file_num']).shape[0]}"
    )
    return df

# Atom-type to element mapping for AutoDock/Vina PDBQT
_ELEMENT_MAP = {
    "A": "C",   # aromatic carbon
    "C": "C",
    "HD": "H", "H": "H",
    "N": "N", "NA": "N",
    "O": "O", "OA": "O",
    "S": "S", "SA": "S",
    "P": "P",
    "F": "F",
    "CL": "Cl", "Cl": "Cl", "cl": "Cl",
    "BR": "Br", "Br": "Br", "br": "Br",
    "I": "I",
    "B": "B", "SI": "Si", "Si": "Si", "se": "Se", "Se": "Se"
}

def _hill_formula_from_counts(counts: Dict[str, int]) -> str:
    """Return Hill-system formula string from element counts dict."""
    norm = {}
    for k, v in counts.items():
        if not v:
            continue
        kk = k[0].upper() + (k[1:].lower() if len(k) > 1 else "")
        norm[kk] = norm.get(kk, 0) + int(v)

    c = norm.pop("C", 0)
    h = norm.pop("H", 0)
    parts = []
    if c:
        parts.append(f"C{c if c>1 else ''}")
        if h:
            parts.append(f"H{h if h>1 else ''}")
    else:
        if h:
            parts.append(f"H{h if h>1 else ''}")

    for el in sorted(norm.keys()):
        n = norm[el]
        parts.append(f"{el}{n if n>1 else ''}")
    return "".join(parts) if parts else ""

def _parse_pdbqt_element_counts(lines: List[str]) -> Dict[str, int]:
    """Count elements from PDBQT ATOM/HETATM lines using the last token (AutoDock type)."""
    counts: Dict[str, int] = {}
    for ln in lines:
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        toks = ln.split()
        if not toks:
            continue
        ad_type = toks[-1]  # AutoDock type is typically the last token
        el = _ELEMENT_MAP.get(ad_type, None)
        if el is None:
            # Fallback: try atom-name field token
            atom_name = toks[2] if len(toks) > 2 else ""
            prefix = "".join([ch for ch in atom_name if ch.isalpha()])[:2]
            el = _ELEMENT_MAP.get(prefix, None) or (prefix.capitalize() if prefix else None)
        if el:
            counts[el] = counts.get(el, 0) + 1
    return counts

_rx_remark_smiles = re.compile(r"(?i)^REMARK\s+(?:SMILES|SMI)\s*[:=]\s*(\S+)")
_rx_remark_inchikey = re.compile(r"(?i)^REMARK.*?(INCHIKEY|InChIKey)\s*[:=]\s*([A-Z0-9\-]+)")
_rx_remark_name = re.compile(r"(?i)^REMARK\s+(?:NAME|Title)\s*[:=]\s*(.+)$")

def parse_pdbqt_identity(p: Path) -> Dict[str, Any]:
    """Extract identity hints from a PDBQT file.

    Returns:
      {
        'path': str,
        'remark_smiles': Optional[str],
        'remark_inchikey': Optional[str],
        'remark_name': Optional[str],
        'pdbqt_formula': Optional[str],
        'pdbqt_heavy_atoms': Optional[int],
      }
    """
    out = {
        "path": str(p), "remark_smiles": None, "remark_inchikey": None, "remark_name": None,
        "pdbqt_formula": None, "pdbqt_heavy_atoms": None
    }
    try:
        txt = p.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return out

    # REMARK fields
    for ln in txt:
        m = _rx_remark_smiles.search(ln)
        if m and not out["remark_smiles"]:
            out["remark_smiles"] = m.group(1).strip()

        m2 = _rx_remark_inchikey.search(ln)
        if m2 and not out["remark_inchikey"]:
            out["remark_inchikey"] = m2.group(2).strip()

        m3 = _rx_remark_name.search(ln)
        if m3 and not out["remark_name"]:
            out["remark_name"] = m3.group(1).strip()

    # Composition fingerprint
    counts = _parse_pdbqt_element_counts(txt)
    if counts:
        out["pdbqt_formula"] = _hill_formula_from_counts(counts)
        # "heavy" excludes hydrogens
        heavy = sum(v for k, v in counts.items() if k.upper() != "H")
        out["pdbqt_heavy_atoms"] = heavy

    return out

def select_representative_pdbqt(files_df: pd.DataFrame) -> pd.DataFrame:
    """
    One row per (scheme, file_num), preferring the highest stage if present.
    Adds parsed PDBQT identity columns.
    """
    if files_df.empty:
        return files_df.assign(remark_smiles=None, remark_inchikey=None, remark_name=None,
                               pdbqt_formula=None, pdbqt_heavy_atoms=None)

    # Pick best path per group
    def _pick(group: pd.DataFrame) -> pd.Series:
        g = group.copy()
        # Prefer higher stage; if none, just first
        g["stage_rank"] = pd.to_numeric(g["stage"], errors="coerce").fillna(-1).astype("int32")
        g = g.sort_values(["stage_rank"], ascending=[False])
        return g.iloc[0]


    cols = [c for c in files_df.columns if c not in ("scheme", "file_num")]
    reps = (
        files_df
        .groupby(["scheme", "file_num"], group_keys=False)[cols]
        .apply(_pick)
        .reset_index(drop=False)  # scheme/file_num come back from the index
    )
    # Parse identity for each chosen path
    id_rows = []
    for _, r in reps.iterrows():
        id_rows.append(parse_pdbqt_identity(Path(r["path"])))
    id_df = pd.DataFrame(id_rows)

    reps = reps.merge(id_df, on="path", how="left")
    return reps

# SDF composition features
try:
    from rdkit.Chem import rdMolDescriptors as rdMD
except Exception:
    rdMD = None

from concurrent.futures import ProcessPoolExecutor

def _calc_formula_heavy(smi: Optional[str], smi_neutral: Optional[str]) -> Tuple[Optional[str], Optional[int]]:
    try:
        s = smi_neutral or smi
        if not s: return None, None
        mol = Chem.MolFromSmiles(s)
        if not mol: return None, None
        try:
            from rdkit.Chem import rdMolDescriptors as rdMD
            f = rdMD.CalcMolFormula(mol)
        except Exception:
            # manual fallback
            counts = {}
            for a in mol.GetAtoms():
                el = a.GetSymbol()
                counts[el] = counts.get(el, 0) + 1
            f = _hill_formula_from_counts(counts) or None
        heavy = sum(1 for a in mol.GetAtoms() if a.GetSymbol() != "H")
        return f, int(heavy)
    except Exception:
        return None, None

import os
def sdf_add_formula_features(df_sdf: pd.DataFrame) -> pd.DataFrame:
    inputs = list(zip(df_sdf.get("smiles"), df_sdf.get("smiles_neutral")))
    # Processes ≈ CPU cores is fine
    max_workers = int(cfg.get("IDENTIFY_MAX_WORKERS", 0)) or (os.cpu_count() or 2)
    with ProcessPoolExecutor(max_workers=max_workers) as ex:

        results = list(ex.map(_calc_formula_heavy_unpack, inputs))
    formulas, heavies = zip(*results) if results else ([], [])
    out = df_sdf.copy()
    out["sdf_formula"] = list(formulas)
    out["sdf_heavy_atoms"] = list(heavies)
    return out


def match_using_pdbqt_features(reps_df: pd.DataFrame, sdf_df: pd.DataFrame, default_one_based: bool) -> pd.DataFrame:
    """
    For each representative PDBQT:
      1) Try InChIKey exact match to SDF inchikey
      2) Try SMILES → InChIKey match
      3) Try composition match: (formula, heavy atoms)
      4) Fallback: index alignment (infer one-based vs zero-based)
    Returns a merged mapping with a 'match_method' column.
    """
    if reps_df.empty:
        raise RuntimeError("No representative PDBQT files to match.")

    # Decide index fallback globally
    inferred_one_based = infer_indexing(
        files_df=reps_df[["file_num"]].copy(),
        sdf_df=sdf_df[["sdf_index"]],
        default_one_based=default_one_based
    )

    # Fast lookups
    by_inchikey = {str(ik): i for i, ik in zip(sdf_df["sdf_index"], sdf_df["inchikey"].fillna(""))
                   if isinstance(ik, str) and ik}
    ik14 = {str(ik)[:14]: i for i, ik in zip(sdf_df["sdf_index"], sdf_df["inchikey"].fillna(""))
            if isinstance(ik, str) and ik}

    comp_index: Dict[Tuple[str, int], List[int]] = {}
    for i, f, h in zip(sdf_df["sdf_index"], sdf_df["sdf_formula"], sdf_df["sdf_heavy_atoms"]):
        if not isinstance(f, str) or pd.isna(f) or pd.isna(h):
            continue
        comp_index.setdefault((f, int(h)), []).append(int(i))

    matched_rows = []
    for _, row in reps_df.iterrows():
        method = None
        target_sdf_index = None

        # 1) InChIKey from REMARK
        ik = row.get("remark_inchikey")
        if isinstance(ik, str) and ik in by_inchikey:
            target_sdf_index = by_inchikey[ik]
            method = "inchikey"
        else:
            # 2) SMILES from REMARK (compute IK and match)
            smi = row.get("remark_smiles")
            if isinstance(smi, str) and smi.strip():
                try:
                    mol = Chem.MolFromSmiles(smi)
                    if mol:
                        ik2 = inchi.MolToInchiKey(mol)
                        if ik2 in by_inchikey:
                            target_sdf_index = by_inchikey[ik2]
                            method = "smiles→inchikey"
                        elif ik2[:14] in ik14:
                            target_sdf_index = ik14[ik2[:14]]
                            method = "smiles→IK14"
                except Exception:
                    pass

        # 3) Composition match
        if target_sdf_index is None:
            pf = row.get("pdbqt_formula")
            ph = row.get("pdbqt_heavy_atoms")
            if isinstance(pf, str) and isinstance(ph, (int, float)):
                candidates = comp_index.get((pf, int(ph)), [])
                if len(candidates) == 1:
                    target_sdf_index = candidates[0]
                    method = "composition"
                elif len(candidates) > 1:
                    idx0 = int(row["file_num"] + (1 if inferred_one_based else 0))
                    if idx0 in candidates:
                        target_sdf_index = idx0
                        method = "composition+index"
                    else:
                        target_sdf_index = min(candidates)
                        method = f"composition(ambig:{len(candidates)})"

        # 4) Fallback: index-only merge
        if target_sdf_index is None:
            target_sdf_index = int(row["file_num"] + (1 if inferred_one_based else 0))
            method = "index-fallback"

        matched_rows.append({
            **row.to_dict(),
            "sdf_index": target_sdf_index,
            "match_method": method
        })

    matched_df = pd.DataFrame(matched_rows)
    mapping = matched_df.merge(sdf_df, on="sdf_index", how="left")
    return mapping

# ============================ Enrichment: PubChem ============================

def parent_cid_from_cid(cid: Optional[int]) -> Optional[int]:
    """Ask PubChem for the parent CID of a given CID (PUG REST)."""
    if cid is None:
        return None
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{int(cid)}/cids/JSON?cids_type=parent"
        r = REQ.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()
        info = data.get("InformationList", {}).get("Information", [])
        if info and "CID" in info[0] and info[0]["CID"]:
            return int(info[0]["CID"][0])
    except Exception as e:
        log.debug(f"Parent CID lookup failed for {cid}: {e}")
    return None

def _is_generic_like(s: str) -> bool:
    if not s: return False
    t = s.strip()
    if any(ch.isdigit() for ch in t): return False
    if any(c in t for c in [",",";","(",")","[","]","{","}","/","\\"]): return False
    if t.count(" ") > 2: return False  # allow one or two words
    return t == t.lower()  # prefer lowercase generics

def _is_brand_like(s: str) -> bool:
    if not s: return False
    t = _strip_marks(s).strip()
    if any(ch.isdigit() for ch in t): return False
    if any(c in t for c in [",",";","(",")","[","]","{","}","/","\\"]): return False
    parts = t.split()
    if len(parts) > 2: return False
    if t == t.lower(): return False  # brands usually TitleCase / ALLCAPS
    return True

def _score_generic_candidate(s: str) -> float:
    # lower score = better
    base = len(s)
    if " " in s: base += 3
    if "-" in s: base += 1
    if any(s.endswith(suf) for suf in GENERIC_SUFFIX_BONUS):
        base -= 2
    return base

def choose_best_generic(synonyms: Optional[List[str]], record_title: Optional[str] = None) -> Optional[str]:
    syns = list(dict.fromkeys([_norm_ws(_strip_marks(x)) for x in (synonyms or []) if isinstance(x, str)]))
    candidates: List[str] = []
    for s in syns:
        cleaned = _strip_salt_hydrate_tail(s)
        if _is_generic_like(cleaned):
            candidates.append(cleaned)
    if record_title:
        rt = _strip_salt_hydrate_tail(_norm_ws(_strip_marks(record_title)))
        if _is_generic_like(rt):
            candidates.append(rt)
    if not candidates:
        return None
    candidates = sorted(set(candidates), key=lambda x: (_score_generic_candidate(x), len(x), x))
    return candidates[0]

def collect_brand_names(synonyms: Optional[List[str]], limit: int = 3) -> List[str]:
    syns = list(dict.fromkeys([_norm_ws(_strip_marks(x)) for x in (synonyms or []) if isinstance(x, str)]))
    brands = [s for s in syns if _is_brand_like(s)]
    seen = set()
    out: List[str] = []
    for b in brands:
        key = b.lower()
        if key not in seen:
            out.append(b)
            seen.add(key)
        if len(out) >= limit:
            break
    return out

def fetch_synonyms_and_title(cid: Optional[int]) -> Tuple[List[str], Optional[str]]:
    if cid is None:
        return [], None
    syns: List[str] = []
    title: Optional[str] = None
    try:
        u = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{int(cid)}/synonyms/JSON"
        r = REQ.get(u, timeout=10)
        r.raise_for_status()
        data = r.json()
        info = data.get("InformationList", {}).get("Information", [])
        if info and "Synonym" in info[0]:
            syns = info[0]["Synonym"]
        for s in syns:
            ss = s.strip()
            if any(c in ss for c in [",",";","(",")","[","]","{","}","/","\\"]):
                continue
            if any(ch.isdigit() for ch in ss):
                continue
            if ss.lower() == ss and 1 <= ss.count(" ") <= 2:
                title = ss
                break
        if not title and syns:
            title = syns[0]
    except Exception:
        pass
    return syns, title

def fetch_pubchem_view_by_cid(cid: Optional[int]) -> Dict[str, Any]:
    if cid is None:
        return {}
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?response_type=display"
    out = {"record_title": None, "iupac_name": None, "synonyms": [], "sources": {}}
    try:
        r = REQ.get(url, timeout=10)
        r.raise_for_status()
        j = r.json()
        def walk(sec):
            if not isinstance(sec, dict): return
            for inf in sec.get("Information", []):
                nm = inf.get("Name") or ""
                val = inf.get("Value", {})
                if "IUPAC Name" in nm:
                    vals = val.get("StringWithMarkup", [])
                    if vals and not out["iupac_name"]:
                        out["iupac_name"] = vals[0].get("String")
                if "Record Title" in nm:
                    vals = val.get("StringWithMarkup", [])
                    if vals and not out["record_title"]:
                        out["record_title"] = vals[0].get("String")
                if "Synonyms" in nm:
                    out["synonyms"].extend(val.get("String", []) or [])
                if any(k in nm for k in ["DrugBank", "UNII", "RXCUI", "RxCUI"]):
                    arr = []
                    if "String" in val: arr += val["String"]
                    if "StringWithMarkup" in val:
                        arr += [x.get("String") for x in val["StringWithMarkup"] if isinstance(x, dict)]
                    key = nm.split()[0].upper()
                    out["sources"].setdefault(key, [])
                    out["sources"][key].extend([x for x in arr if isinstance(x, str)])
            for s in sec.get("Section", []):
                walk(s)
        for top in j.get("Record", {}).get("Section", []):
            walk(top)
        out["synonyms"] = list(dict.fromkeys([_norm_ws(_strip_marks(s)) for s in out["synonyms"] if isinstance(s, str)]))
        # normalize UNIIs
        if "UNII" in out["sources"]:
            out["sources"]["UNII"] = [re.sub(r"[^A-Z0-9]", "", u.upper()) for u in out["sources"]["UNII"] if isinstance(u, str)]
        return out
    except Exception:
        return out

def merge_name_buckets(base_syns: List[str], view_syns: List[str]) -> List[str]:
    syns = base_syns + view_syns
    cleaned = []
    for s in syns:
        if not isinstance(s, str):
            continue
        if sum(1 for c in s if c in ",;()[]{}\\/") > 3:
            continue
        cleaned.append(s)
    return list(dict.fromkeys(cleaned))

def pick_best_display_name(candidates: List[str], record_title: Optional[str], iupac: Optional[str]) -> Optional[str]:
    pool = []
    for c in candidates:
        cc = _clean_display_name(c or "")
        if cc: pool.append(cc)
    if record_title:
        rt = _clean_display_name(record_title)
        if rt: pool.append(rt)
    pool = list(dict.fromkeys(pool))
    generics = [p for p in pool if _is_generic_like(p)]
    if generics:
        generics = sorted(set(generics), key=lambda x: (_score_generic_candidate(x), len(x), x))
        return generics[0]
    if iupac:
        return _norm_ws(iupac)
    return pool[0] if pool else None

def guess_generic_from_title(title: Optional[str]) -> Optional[str]:
    if not title:
        return None
    t = _strip_marks(_norm_ws(title))
    if re.match(r"^(ZINC|CHEMBL|DB|CID|UNII|CAS)[_:\-\s0-9A-Z]+$", t, re.I):
        return None
    t2 = _strip_salt_hydrate_tail(_strip_form_words(t))
    if _is_generic_like(t2):
        return t2
    return None

def pubchem_lookup(smiles: Optional[str], inchikey_val: Optional[str], cache: dict) -> Dict[str, Any]:
    import pubchempy as pcp
    def cache_get(k):
        v = cache.get(k)
        return v if isinstance(v, dict) else None
    def cache_put(k, entry):
        cache[k] = entry
    def extract_entry(cmpds):
        if not cmpds: return None
        c = cmpds[0]
        syns = getattr(c, "synonyms", None) or []
        title = getattr(c, "title", None) or None
        name  = getattr(c, "iupac_name", None) or (syns[0] if syns else None)
        cid   = getattr(c, "cid", None)
        view = fetch_pubchem_view_by_cid(cid)
        view_syns = view.get("synonyms", []) or []
        view_title = view.get("record_title")
        view_iupac = view.get("iupac_name")
        syns = merge_name_buckets(syns, view_syns)
        if not title and view_title:
            title = view_title
        if (not name) and view_iupac:
            name = view_iupac
        parent = parent_cid_from_cid(cid)
        if parent and parent != cid:
            p = pcp.get_compounds(parent, "cid")
            if p:
                c2 = p[0]
                syns2 = getattr(c2, "synonyms", None) or []
                title2 = getattr(c2, "title", None)
                view2 = fetch_pubchem_view_by_cid(parent)
                syns2 = merge_name_buckets(syns2, view2.get("synonyms", []) or [])
                title2 = title2 or view2.get("record_title")
                iupac2 = view2.get("iupac_name")
                merged = merge_name_buckets(syns, syns2)
                generic = choose_best_generic(merged, title2 or title)
                brands = collect_brand_names(merged)
                return {
                    "name": iupac2 or name,
                    "cid": parent,
                    "title": title2 or title,
                    "synonyms": merged,
                    "sources": view2.get("sources", {}),
                    "generic": generic,
                    "brands": brands,
                    "iupac_name": iupac2 or name,
                }
        generic = choose_best_generic(syns, title)
        brands = collect_brand_names(syns)
        return {
            "name": view_iupac or name,
            "cid": cid,
            "title": title,
            "synonyms": syns,
            "sources": view.get("sources", {}),
            "generic": generic,
            "brands": brands,
            "iupac_name": view_iupac or name,
        }
    def attempt(query, namespace):
        delay = RATE_LIMIT_SEC
        for attempt_no in range(5):
            try:
                cs = pcp.get_compounds(query, namespace)
                return extract_entry(cs)
            except Exception as e:
                if attempt_no < 2:
                    log.warning(f"PubChem query failed ({namespace}='{query}'): {e}")
                else:
                    log.warning(f"PubChem retry {attempt_no+1} for {namespace}")
                time.sleep(delay)
                delay = min(delay * 2, 5.0)
        return None
    key_smiles = smiles if isinstance(smiles, str) else None
    key_inchi  = f"IK:{inchikey_val}" if isinstance(inchikey_val, str) else None
    for k in (key_smiles, key_inchi):
        if k:
            c = cache_get(k)
            if c:
                if (not c.get("synonyms")) or (c.get("title") is None) or (c.get("iupac_name") is None):
                    view = fetch_pubchem_view_by_cid(c.get("cid"))
                    merged_syns = merge_name_buckets(c.get("synonyms", []), view.get("synonyms", []) or [])
                    c["synonyms"] = merged_syns
                    c["title"] = c.get("title") or view.get("record_title")
                    c["iupac_name"] = c.get("iupac_name") or view.get("iupac_name")
                    c["generic"] = choose_best_generic(c.get("synonyms"), c.get("title"))
                    c["brands"]  = collect_brand_names(c.get("synonyms"))
                return c
    entry = attempt(key_smiles, "smiles") if key_smiles else None
    if entry:
        cache_put(key_smiles, entry)
        neut = neutralize_smiles(key_smiles)
        if neut: cache_put(neut, entry)
        if entry.get("cid"): cache_put(f"CID:{entry['cid']}", entry)
        return entry
    entry = attempt(inchikey_val, "inchikey") if inchikey_val else None
    if entry:
        cache_put(key_inchi, entry)
        if entry.get("cid"): cache_put(f"CID:{entry['cid']}", entry)
        return entry
    neut = neutralize_smiles(key_smiles)
    if neut and neut != key_smiles:
        entry = attempt(neut, "smiles")
        if entry:
            cache_put(neut, entry)
            if entry.get("cid"): cache_put(f"CID:{entry['cid']}", entry)
            return entry
    return {"name": None, "cid": None, "title": None, "synonyms": [], "generic": None, "brands": [], "iupac_name": None}

from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

def enrich_with_pubchem(mapping: pd.DataFrame) -> pd.DataFrame:
    cache = load_cache(CACHE_PATH_PUBCHEM)
    log.info(f"Loaded PubChem cache entries: {len(cache)}")
    # Normalize/repair old cache entries
    changed = False
    for k, v in list(cache.items()):
        if isinstance(v, dict) and v.get("cid") and (not v.get("synonyms") or v.get("title") is None or v.get("iupac_name") is None):
            syns, ttl = fetch_synonyms_and_title(v["cid"])
            view = fetch_pubchem_view_by_cid(v.get("cid"))
            merged_syns = merge_name_buckets(syns, view.get("synonyms", []) or [])
            if merged_syns: v["synonyms"] = merged_syns
            if ttl and not v.get("title"): v["title"] = ttl
            if not v.get("title") and view.get("record_title"): v["title"] = view.get("record_title")
            if not v.get("iupac_name") and view.get("iupac_name"): v["iupac_name"] = view.get("iupac_name")
            v["generic"] = choose_best_generic(v.get("synonyms"), v.get("title"))
            v["brands"]  = collect_brand_names(v.get("synonyms"))
            changed = True
        elif isinstance(v, dict) and ("generic" not in v or "brands" not in v):
            v["generic"] = choose_best_generic(v.get("synonyms"), v.get("title"))
            v["brands"]  = collect_brand_names(v.get("synonyms"))
            changed = True
    if changed:
        save_cache(cache, CACHE_PATH_PUBCHEM)

    # Deduplicate queries (prefer SMILES; fallback to neutral → InChIKey)
    queries = []
    for _, row in mapping.iterrows():
        s  = row.get("smiles", None)
        sn = row.get("smiles_neutral", None)
        ik = row.get("inchikey", None)
        if isinstance(s, str) and s not in cache:
            queries.append(("smiles", s, s, ik)); continue
        if isinstance(sn, str) and sn not in cache:
            queries.append(("smiles", sn, sn, ik)); continue
        if (not isinstance(s, str)) and isinstance(ik, str) and f"IK:{ik}" not in cache:
            queries.append(("inchikey", f"IK:{ik}", s, ik))

    if MAX_LOOKUPS is not None:
        queries = queries[:MAX_LOOKUPS]
    log.info(f"PubChem lookups to perform (deduped): {len(queries)}")

    # Worker returns (cache_key, entry, optional neutral_smiles, cid)
    def _worker(item):
        _, cache_key, s, ik = item
        # polite spacing per worker; tiny sleep helps smooth bursts
        time.sleep(RATE_LIMIT_SEC)
        entry = pubchem_lookup(s, ik, {})  # use a local dict; merge later
        neut = neutralize_smiles(s) if isinstance(s, str) else None
        cid  = entry.get("cid")
        return cache_key, entry, neut, cid

    # Small pool for I/O-bound API calls
    results_local: Dict[str, dict] = {}
    extra_keys: List[Tuple[str, dict]] = []
    http_workers = int(cfg.get("IDENTIFY_HTTP_WORKERS", 4))
    with ThreadPoolExecutor(max_workers=http_workers) as ex:
        for cache_key, entry, neut, cid in ex.map(_worker, queries):
            results_local[cache_key] = entry
            if isinstance(neut, str):
                extra_keys.append((neut, entry))
            if cid:
                extra_keys.append((f"CID:{cid}", entry))

    # Merge once at the end (thread-safe, single thread)
    for k, v in results_local.items():
        cache[k] = v
    for k, v in extra_keys:
        cache[k] = v

    save_cache(cache, CACHE_PATH_PUBCHEM)
    log.info("Saved PubChem cache (parallel phase complete).")

    for i, (_, cache_key, s, ik) in enumerate(queries, start=1):
        time.sleep(RATE_LIMIT_SEC)
        entry = pubchem_lookup(s, ik, cache)
        cache[cache_key] = entry
        if isinstance(s, str):
            neut = neutralize_smiles(s)
            if neut: cache[neut] = entry
        if entry.get("cid"):
            cache[f"CID:{entry['cid']}"] = entry
        if i % 50 == 0:
            save_cache(cache, CACHE_PATH_PUBCHEM)
            log.info(f"PubChem progress: {i}/{len(queries)} cached…")

    save_cache(cache, CACHE_PATH_PUBCHEM)
    log.info("Saved PubChem cache.")

    # Attach resolved fields (do NOT overwrite existing SDF)
    def get_entry(row):
        s = row.get("smiles", None)
        sn = row.get("smiles_neutral", None)
        ik = row.get("inchikey", None)
        cid_val = row.get("pubchem_cid", None)
        if isinstance(s, str) and s in cache: return cache[s]
        if isinstance(sn, str) and sn in cache: return cache[sn]
        if isinstance(ik, str) and f"IK:{ik}" in cache: return cache[f"IK:{ik}"]
        try:
            cid_int = int(cid_val) if pd.notna(cid_val) else None
        except Exception:
            cid_int = None
        if cid_int and f"CID:{cid_int}" in cache: return cache[f"CID:{cid_int}"]
        return {}

    mapping["pubchem_name"] = mapping.apply(lambda r: get_entry(r).get("name"), axis=1)
    mapping["pubchem_cid_resolved"] = mapping.apply(lambda r: get_entry(r).get("cid"), axis=1)
    mapping["pubchem_record_title"] = mapping.apply(lambda r: get_entry(r).get("title"), axis=1)
    mapping["generic_name"] = mapping.apply(lambda r: get_entry(r).get("generic"), axis=1)
    mapping["brand_names"] = mapping.apply(lambda r: "; ".join(get_entry(r).get("brands", [])), axis=1)
    mapping["pubchem_iupac_name"] = mapping.apply(lambda r: get_entry(r).get("iupac_name"), axis=1)
    mapping["pubchem_synonyms"] = mapping.apply(
        lambda r: "; ".join((get_entry(r).get("synonyms") or [])) or None, axis=1
    )
    mapping["pubchem_unii_list"] = mapping.apply(
        lambda r: "; ".join((get_entry(r).get("sources", {}) or {}).get("UNII", [])) or None, axis=1
    )

    # normalize CID as int for readability
    def _as_int(x):
        try: return int(x)
        except Exception: return None
    mapping["pubchem_cid_resolved"] = mapping["pubchem_cid_resolved"].apply(_as_int)

    # Rescue generic if still missing
    def rescue_generic(row):
        if pd.notna(row.get("generic_name")) and row.get("generic_name"):
            return row.get("generic_name")
        guess = guess_generic_from_title(row.get("sdf_title"))
        if guess:
            return guess
        for k in ("drug_name","drugname","name"):
            v = row.get(k)
            if isinstance(v, str):
                g = guess_generic_from_title(v)
                if g:
                    return g
        return None
    mapping["generic_name"] = mapping.apply(rescue_generic, axis=1)

    # Friendly display name
    def pick_display(row):
        return (_clean_display_name(row.get("generic_name"))
                or _clean_display_name(row.get("pubchem_record_title"))
                or _clean_display_name(row.get("pubchem_iupac_name"))
                or _clean_display_name(row.get("pubchem_name"))
                or row.get("sdf_title"))
    mapping["display_name"] = mapping.apply(pick_display, axis=1)

    hits = int(mapping["pubchem_cid_resolved"].notna().sum())
    log.info(f"PubChem hits: {hits}/{len(mapping)}")
    return mapping

# ============================ Enrichment: RxNorm =============================

RXN_BASE = "https://rxnav.nlm.nih.gov/REST"

def rxn_get(url: str, params: Dict[str, Any] = None) -> Optional[dict]:
    try:
        r = REQ.get(url, params=params or {}, timeout=20)
        r.raise_for_status()
        return r.json()
    except Exception:
        return None

def rxnorm_from_unii(unii: Optional[str]) -> Optional[str]:
    """Use findRxcuiById with idtype=UNII."""
    if not unii:
        return None
    j = rxn_get(f"{RXN_BASE}/findRxcuiById.json", {"idtype": "UNII", "id": unii, "allsrc": 1})
    try:
        ids = j.get("idGroup", {}).get("rxnormId", [])
        return ids[0] if ids else None
    except Exception:
        return None

def rxnorm_from_name(name: Optional[str]) -> Optional[str]:
    """Try normalized match, then approximate, then getDrugs path."""
    if not name:
        return None
    # 1) Normalized exact-ish
    j1 = rxn_get(f"{RXN_BASE}/rxcui.json", {"name": name, "allsrc": 1})
    try:
        ids = j1.get("idGroup", {}).get("rxnormId", [])
        if ids:
            return ids[0]
    except Exception:
        pass
    # 2) Approximate
    j2 = rxn_get(f"{RXN_BASE}/approximateTerm.json", {"term": name, "maxEntries": 1})
    try:
        cands = j2.get("approximateGroup", {}).get("candidate", []) or []
        if cands:
            return cands[0].get("rxcui")
    except Exception:
        pass
    # 3) getDrugs by name
    j3 = rxn_get(f"{RXN_BASE}/drugs.json", {"name": name})
    try:
        drugs = j3.get("drugGroup", {}).get("conceptGroup", []) or []
        for grp in drugs:
            for c in grp.get("conceptProperties", []) or []:
                rx = c.get("rxcui")
                if rx:
                    return rx
    except Exception:
        pass
    return None

def rxnorm_names_for_rxcui(rxcui: str) -> Dict[str, Any]:
    out = {"rxcui": rxcui, "generic": None, "brands": []}
    j = rxn_get(f"{RXN_BASE}/rxcui/{rxcui}/allrelated.json")
    if not j:
        return out
    concept = j.get("allRelatedGroup", {}).get("conceptGroup", []) or []
    generics, brands = [], []
    for cg in concept:
        tty = cg.get("tty")
        for c in cg.get("conceptProperties", []) or []:
            nm = c.get("name")
            if not isinstance(nm, str):
                continue
            if tty in ("IN", "PIN"):
                generics.append(nm)
            if tty in ("BN", "SBD"):
                brands.append(nm)
    if generics:
        generics_sorted = sorted(generics, key=lambda x: (x != x.lower(), len(x)))
        out["generic"] = generics_sorted[0]
    seen = set(); out_brands = []
    for b in brands:
        lb = b.lower()
        if lb in seen:
            continue
        seen.add(lb); out_brands.append(b)
        if len(out_brands) >= 5:
            break
    out["brands"] = out_brands
    return out

def _rx_candidates(obj: Dict[str, Any]) -> List[str]:
    """Build clean 1–2 word candidates for RxNorm/OpenFDA (trimmed)."""
    raw: List[Optional[str]] = []
    raw += [obj.get("generic") if "generic" in obj else obj.get("generic_name")]
    raw += [obj.get("pubchem_name")]
    raw += [obj.get("display_name")]
    syns = obj.get("pubchem_synonyms")
    if isinstance(syns, str) and syns.strip():
        raw += syns.split(";")[:10]
    cands: List[str] = []
    for nm in raw:
        if not isinstance(nm, str):
            continue
        nm2 = _clean_display_name(nm)
        if not nm2:
            continue
        if nm2.count(" ") > 2:
            continue
        cands.append(nm2)
    seen = set(); out=[]
    for c in cands:
        lc=c.lower()
        if lc in seen:
            continue
        seen.add(lc); out.append(c)
    return out

def enrich_with_rxnorm(df: pd.DataFrame) -> pd.DataFrame:
    rx_rxcui, rx_generic, rx_brands = [], [], []

    def first_unii(row) -> Optional[str]:
        pu = row.get("pubchem_unii_list")
        if isinstance(pu, str) and pu.strip():
            return pu.split(";")[0].strip()
        return None

    for _, row in df.iterrows():
        rxcui_val = None

        # 1) UNII → RxCUI
        rxcui_val = rxnorm_from_unii(first_unii(row))

        # 2) fallback: names
        if not rxcui_val:
            nm_obj = {
                "generic_name": row.get("generic_name"),
                "pubchem_name": row.get("pubchem_name"),
                "display_name": row.get("display_name"),
                "pubchem_synonyms": row.get("pubchem_synonyms"),
            }
            for nm in _rx_candidates(nm_obj):
                rxcui_val = rxnorm_from_name(nm)
                if rxcui_val:
                    break

        generic = None; brands = []
        if rxcui_val:
            info = rxnorm_names_for_rxcui(rxcui_val)
            generic = info.get("generic")
            brands = info.get("brands") or []

        rx_rxcui.append(rxcui_val)
        rx_generic.append(generic)
        rx_brands.append("; ".join(brands) if brands else None)

        time.sleep(RATE_LIMIT_SEC)

    df["rxnorm_rxcui"] = rx_rxcui
    df["rxnorm_generic_name"] = rx_generic
    df["rxnorm_brand_names"] = rx_brands

    hits = int(pd.Series(rx_rxcui).notna().sum())
    log.info(f"RxNorm hits: {hits}/{len(df)}")
    return df

# =========================== Enrichment: DrugCentral =========================

def _load_drugcentral_index(tsv_path: Path):
    if not tsv_path.exists():
        log.warning(f"DrugCentral TSV not found at {tsv_path}; skipping DrugCentral enrichment.")
        return None, None
    df = pd.read_csv(tsv_path, sep="\t|,", engine="python")
    # normalize columns
    cols = {c.lower(): c for c in df.columns}
    ik_col = cols.get("inchikey") or cols.get("std_inchi_key") or cols.get("inchi_key")
    nm_col = cols.get("name") or cols.get("inn") or cols.get("generic_name")
    id_col = cols.get("id") or cols.get("drugcentral_id")
    if not ik_col:
        log.warning("DrugCentral TSV missing InChIKey column; skipping.")
        return None, None
    df = df.dropna(subset=[ik_col]).copy()
    df["__ik"] = df[ik_col].astype(str).str.strip()
    df["__ik14"] = df["__ik"].str[:14]
    return df, {"ik": "__ik", "ik14": "__ik14", "name": nm_col, "id": id_col}

_DRUGCENTRAL_DF, _DC_COLS = _load_drugcentral_index(DRUGCENTRAL_TSV)

def enrich_with_drugcentral(df: pd.DataFrame) -> pd.DataFrame:
    if _DRUGCENTRAL_DF is None:
        df["drugcentral_id"] = None
        df["drugcentral_generic_name"] = None
        df["drugcentral_brand_names"] = None
        log.info(f"DrugCentral hits: 0/{len(df)}")
        return df

    hits = 0
    dc_ids, dc_generic = [], []
    for _, row in df.iterrows():
        ik = row.get("inchikey")
        rec = None
        if isinstance(ik, str) and ik:
            # exact InChIKey
            m = _DRUGCENTRAL_DF[_DRUGCENTRAL_DF[_DC_COLS["ik"]] == ik]
            if m.empty:
                # fallback: connectivity block (salt/tautomer tolerant)
                m = _DRUGCENTRAL_DF[_DRUGCENTRAL_DF[_DC_COLS["ik14"]] == ik[:14]]
            if not m.empty:
                rec = m.iloc[0]

        if rec is not None:
            dc_ids.append(rec.get(_DC_COLS["id"]))
            dc_generic.append(rec.get(_DC_COLS["name"]))
            hits += 1
        else:
            # optional name-based assist if InChIKey was missing
            nm = _clean_display_name(row.get("generic_name") or row.get("pubchem_name") or row.get("display_name"))
            if nm and _DC_COLS["name"] in _DRUGCENTRAL_DF.columns:
                mm = _DRUGCENTRAL_DF[_DRUGCENTRAL_DF[_DC_COLS["name"]].astype(str).str.lower() == nm.lower()]
                if not mm.empty:
                    rec = mm.iloc[0]
                    dc_ids.append(rec.get(_DC_COLS["id"]))
                    dc_generic.append(rec.get(_DC_COLS["name"]))
                    hits += 1
                else:
                    dc_ids.append(None); dc_generic.append(None)
            else:
                dc_ids.append(None); dc_generic.append(None)

    df["drugcentral_id"] = dc_ids
    df["drugcentral_generic_name"] = dc_generic
    df["drugcentral_brand_names"] = None  # TSV typically lacks brand names

    log.info(f"DrugCentral hits: {hits}/{len(df)}")
    return df

# ============================= Enrichment: OpenFDA ===========================

OPENFDA_BASE = "https://api.fda.gov/drug/label.json"

def openfda_query(search: str) -> Optional[List[dict]]:
    params = {"search": search, "limit": 5}
    try:
        r = REQ.get(OPENFDA_BASE, params=params, timeout=20)
        if r.status_code == 404:
            return []  # no results
        r.raise_for_status()
        j = r.json()
        return j.get("results", []) or []
    except Exception:
        return None

def _extract_openfda_names(items: List[dict]) -> Dict[str, Any]:
    gen = []; brand = []; unii = []; rxcui = []
    for it in items:
        of = it.get("openfda", {}) if isinstance(it, dict) else {}
        if isinstance(of.get("generic_name"), list): gen += of["generic_name"]
        if isinstance(of.get("brand_name"), list): brand += of["brand_name"]
        if isinstance(of.get("unii"), list): unii += of["unii"]
        if isinstance(of.get("rxcui"), list): rxcui += of["rxcui"]
    # dedupe & compact
    def _dedup(arr, cap=10):
        out=[]; seen=set()
        for a in arr:
            if not isinstance(a, str): continue
            key=a.strip().lower()
            if not key or key in seen: continue
            seen.add(key); out.append(a.strip())
            if len(out)>=cap: break
        return out
    return {
        "generic": _dedup(gen, cap=3),
        "brands": _dedup(brand, cap=10),
        "unii": _dedup(unii, cap=5),
        "rxcui": _dedup(rxcui, cap=5),
    }

def _rx_candidates_for_openfda(obj: Dict[str, Any]) -> List[str]:
    return _rx_candidates(obj)

def enrich_with_openfda(df: pd.DataFrame) -> pd.DataFrame:
    cache = load_cache(CACHE_PATH_OPENFDA)
    results_generic, results_brand, results_unii, results_rxcui = [], [], [], []
    hits = 0

    for _, row in df.iterrows():
        time.sleep(RATE_LIMIT_SEC)
        # Build queries in priority order
        queries = []

        # 1) UNII from PubChem PUG-View
        unii_pubchem = None
        pu = row.get("pubchem_unii_list")
        if isinstance(pu, str) and pu.strip():
            unii_pubchem = pu.split(";")[0].strip()
        if unii_pubchem:
            queries.append(f'openfda.unii:"{unii_pubchem}"')

        # 2) RXCUI if we already got one
        rxcui_val = row.get("rxnorm_rxcui")
        if isinstance(rxcui_val, str) and rxcui_val.strip():
            queries.append(f'openfda.rxcui:"{rxcui_val.strip()}"')

        # 3) Names (generic, pubchem_name, display) + top PubChem synonyms
        nm_obj = {
            "generic_name": row.get("generic_name"),
            "pubchem_name": row.get("pubchem_name"),
            "display_name": row.get("display_name"),
            "pubchem_synonyms": row.get("pubchem_synonyms"),
        }
        for nm in _rx_candidates_for_openfda(nm_obj):
            queries.append(f'openfda.generic_name:"{nm}"')
            queries.append(f'openfda.brand_name:"{nm}"')

        # Execute until first good hit
        agg = {"generic": [], "brands": [], "unii": [], "rxcui": []}
        got = False
        for q in queries[:30]:  # cap number of remote calls per row
            if q in cache:
                items = cache[q]
            else:
                items = openfda_query(q)
                cache[q] = items
            if items is None:
                continue
            names = _extract_openfda_names(items)
            # accumulate
            agg["generic"] += names["generic"]
            agg["brands"]  += names["brands"]
            agg["unii"]    += names["unii"]
            agg["rxcui"]   += names["rxcui"]
            # consider it's a "hit" if we got at least a generic or brand
            if names["generic"] or names["brands"]:
                got = True
                break

        # finalize row
        gen = list(dict.fromkeys(agg["generic"]))[:1]  # pick best one
        brs = list(dict.fromkeys(agg["brands"]))[:10]
        uniis = list(dict.fromkeys(agg["unii"]))[:3]
        rxs = list(dict.fromkeys(agg["rxcui"]))[:3]

        results_generic.append(gen[0] if gen else None)
        results_brand.append("; ".join(brs) if brs else None)
        results_unii.append("; ".join(uniis) if uniis else None)
        results_rxcui.append("; ".join(rxs) if rxs else None)
        if got: hits += 1

    save_cache(cache, CACHE_PATH_OPENFDA)

    df["openfda_generic_name"] = results_generic
    df["openfda_brand_names"] = results_brand
    df["openfda_unii"] = results_unii
    df["openfda_rxcuis"] = results_rxcui

    log.info(f"OpenFDA rows with names: {hits}/{len(df)}")
    return df

# =============================== main ========================================

def main():
    # 1) Build SDF index → metadata table
    df_sdf = load_sdf_rows(SDF_PATH, ONE_BASED)

    # 2) Scan prepped ligands (supports rdk_* and fda_* and *_stageX)
    files_df = scan_prepped_files(SEARCH_ROOT)
    if files_df.empty:
        raise RuntimeError("No prepped ligand files found; cannot proceed.")

    # 3) Select one representative .pdbqt per ligand and parse identity
    reps_df = select_representative_pdbqt(files_df)

    # Add SDF composition features to enable composition matching
    df_sdf = sdf_add_formula_features(df_sdf)

    # 4) Match representative PDBQTs to SDF rows using identity hints
    mapping = match_using_pdbqt_features(reps_df, df_sdf, ONE_BASED)

    # 5) Sanity + write base mapping immediately
    sanity_check_index_merge(mapping)
    base_csv = Path("fda_mapping_from_pdbqt.csv")
    mapping.to_csv(base_csv, index=False)
    log.info(f"Wrote base mapping (with identity & match_method): {base_csv.resolve()}")

    # 6) PubChem enrichment
    if LOOKUP_PUBCHEM:
        log.info("Starting PubChem enrichment…")
        mapping = enrich_with_pubchem(mapping)
        mapping.to_csv(base_csv, index=False)
        log.info("Updated mapping with PubChem names/CIDs and generic/brand fields.")

    # 7) RxNorm enrichment
    if LOOKUP_RXNORM:
        log.info("Starting RxNorm enrichment…")
        mapping = enrich_with_rxnorm(mapping)
        mapping.to_csv(base_csv, index=False)
        log.info("Added RxNorm RXCUI, generic and brand names (source-specific).")

    # 8) DrugCentral enrichment
    if LOOKUP_DRUGCENTRAL:
        log.info("Starting DrugCentral enrichment…")
        mapping = enrich_with_drugcentral(mapping)
        mapping.to_csv(base_csv, index=False)
        log.info("Added DrugCentral IDs and generic names (source-specific).")

    # 9) OpenFDA enrichment
    if LOOKUP_OPENFDA:
        log.info("Starting OpenFDA enrichment…")
        mapping = enrich_with_openfda(mapping)
        mapping.to_csv(base_csv, index=False)
        log.info("Added OpenFDA generic/brand/UNII/RXCUI (source-specific).")

    # 10) Write problem rows for quick inspection (missing SDF title or SMILES)
    unmapped = mapping[mapping["sdf_title"].isna() | mapping["smiles"].isna()]
    if len(unmapped) > 0:
        bad_csv = Path("fda_mapping_unmapped_or_no_smiles.csv")
        unmapped.to_csv(bad_csv, index=False)
        log.warning(f"Wrote {len(unmapped)} problematic rows to: {bad_csv.resolve()}")

    # 11) Final preview
    log.info("Preview:")
    with pd.option_context("display.max_columns", None, "display.width", 220):
        log.info("\n" + mapping.head(10).to_string(index=False))

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log.exception(f"Fatal error: {e}")
