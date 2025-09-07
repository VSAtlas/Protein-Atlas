

from __future__ import annotations

import argparse
import csv
import glob
import logging
import os
import re
import sys
import threading
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path 
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
# -------- Optional heavy deps are imported lazily where needed --------

# ---- limit simultaneous PyMOL renders (override via env PYMOL_PARALLEL) ----
_RENDER_LOCK = threading.Semaphore(int(os.environ.get("PYMOL_PARALLEL", "1")))

# ---- import pipeline pieces from your main so we reuse logic verbatim ----
from main import (  # noqa: E402
    run_one_stage,
    make_protein_logger,
    make_paths,
    extract_ligands_to_nolig,
    robust_prepare_controls,
    prepare_receptor,
    _count_heavy_atoms_from_pdbqt,
    CenterSelector,
    GlobalCenterGuard,
    get_recenter_params,
    early_recenter_decision,
    fallback_recentering_if_empty,  # imported but kept for parity
    record_score,
    record_le,
    final_pose_validation_and_screenshots,
    RetryManager,
    build_control_lookup,
    _fingerprint_stage,
    checkpoint_should_skip,
    checkpoint_mark_done,
    checkpoint_invalidate_from,
    detect_pocket,
)

from capture_pose import (  # noqa: E402
    pick_control_and_nearest_rdk,
    _render_native_on_original_pdb,
    _render_three_views_with_pymol,
    _safe_open_csv_for_write,
)

from metabolite_resolver import (  # noqa: E402
    load_library_index,
    ensure_parent_drugs_for_controls,
    resolve_corresponding_name_for_rdk,
    resolve_corresponding_name_from_text,
)

from input_and_export_functions import load_inputs, validate_config  # noqa: E402
from protein_functions import detect_active_site  # noqa: E402


# =============================================================================
# Path helpers (WSL-friendly)
# =============================================================================

_WIN_DRIVE_RX = re.compile(r"^[A-Za-z]:[\\/]")
def _is_windowsish_path(s: str) -> bool:
    return bool(_WIN_DRIVE_RX.match(s or ""))

def _wslify(s: str) -> str:
    """
    Convert 'E:\\path\\to\\thing' -> '/mnt/e/path/to/thing' (best effort).
    Only applied if running on POSIX and the path doesn't exist.
    """
    if not s:
        return s
    s2 = os.path.expandvars(os.path.expanduser(s))
    if os.name == "posix" and _is_windowsish_path(s2):
        drive = s2[0].lower()
        tail = s2[2:].replace("\\", "/")
        guess = f"/mnt/{drive}/{tail.lstrip('/')}"
        return guess
    return s2

def _resolve_existing_path(p: str) -> Path:
    """
    Expand ~ and env vars; if missing on POSIX and Windowsy, try WSL-ify.
    Return a Path (existing or not), but logs intent.
    """
    p1 = os.path.expandvars(os.path.expanduser(p or ""))
    p2 = p1
    if os.name == "posix" and not os.path.exists(p1) and _is_windowsish_path(p1):
        p2 = _wslify(p1)
    return Path(p2)

def _echo_path(tag: str, p: Path):
    exists = "✅" if p.exists() else "❌"
    kind = "dir" if p.is_dir() else ("file" if p.is_file() else "path")
    print(f"[paths] {tag}: {p}  ({kind}) {exists}")


# =============================================================================
# Defaults (Windows-friendly paths; override on CLI)
# =============================================================================

# Defaults remain Windows-style; we auto-convert if you're on Linux/WSL.
DEFAULT_INPUT_DIR = r"E:\PythonProject\protein_automation\input_pdbs"
DEFAULT_OUT_ROOT = r"E:\PythonProject\protein_automation\benchmarks"
DEFAULT_PREPPED_DIR = r"E:\PythonProject\protein_automation\prepped_ligands"
DEFAULT_MAPPING_CSV = r"E:\PythonProject\protein_automation\fda_mapping_from_pdbqt.csv"


# =============================================================================
# Normalization utilities
# =============================================================================

def _norm(s: Optional[str]) -> str:
    """Normalize a free-text string for fuzzy matching."""
    if s is None:
        return ""
    s = str(s).strip().lower()
    s = s.replace("\u00A0", " ")
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"[-_/\\,;:\|\[\]\(\)\{\}\.\+\*'\"]+", " ", s)
    return s

def _tokenize(s: str) -> List[str]:
    """Tokenize into alnum chunks (plus '+') for token-overlap scoring."""
    if not s:
        return []
    toks = re.split(r"[^a-z0-9\+]+", s.lower())
    return [t for t in toks if t]

def _norm_text(s: Optional[str]) -> str:
    """Aggressive normalization (alphanum only)."""
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())

def _dedupe_str(seq: Sequence[str]) -> List[str]:
    """Stable dedupe for lists of strings (keeps first occurrence)."""
    seen: set = set()
    out: List[str] = []
    for s in seq:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out

def _path_is_within(child: Path, parent: Path) -> bool:
    """True if 'child' is inside 'parent' (Windows-safe)."""
    try:
        child = child.resolve(strict=False)
        parent = parent.resolve(strict=False)
        child.relative_to(parent)
        return True
    except Exception:
        return False

def _expand_globs(paths: Sequence[str]) -> List[Path]:
    """Expand any glob patterns into Paths."""
    out: List[Path] = []
    for p in paths:
        if any(ch in p for ch in "*?[]"):
            out.extend(Path(x) for x in glob.glob(p))
        else:
            out.append(Path(p))
    return out

def _extract_rdk_id(text: str) -> Optional[str]:
    """Extract an 'rdk_######' id from a filename/path if present."""
    if not text:
        return None
    m = re.search(r"(rdk_\d{6,8})", Path(text).stem.lower())
    return m.group(1) if m else None


# =============================================================================
# Exclusions (expanded) — HET filters and name keyword screens
# =============================================================================

EXCLUDE_HET_IDS: set = {
    # waters
    "HOH", "WAT", "DOD",

    # monoatomic / simple ions
    "F", "CL", "BR", "I",
    "LI", "NA", "K", "RB", "CS", "FR",
    "MG", "CA", "SR", "BA", "RA",
    "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
    "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "HGA", "HG", "PT", "AU",
    "AL", "GA", "IN", "TL", "PB", "BI", "SN", "SB", "AS", "SE",
    "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",

    # salts / crystallization additives
    "SO4", "SUL", "PO4", "HPO", "DPO", "CO3", "NO3", "SCN", "CYN", "BCT", "CAC", "IOD", "BR", "NCO", "CLO",

    # buffers / pH agents / precipitants / reducing agents
    "TRS", "TES", "BES", "HEP", "HEZ", "MES", "MOP", "PIP", "ADA", "CHES", "CAPS", "BTP", "POP", "MOPS", "PIPES",
    "DTT", "BME", "TCEP", "IMD", "IPR", "IPA", "ACN", "ACT", "ACE", "FMT", "CIT", "TAR", "TLA", "GLY",
    # standard residues often used as additives
    "HIS", "ARG", "LYS",

    # cryo / solvents / polyols
    "GOL", "EDO", "PGO", "MPD", "DMS", "DMF", "EOH", "ETOH", "MEO", "PEO", "BU3", "TBA", "2MO", "1PE",
    "PE8", "PEU", "PG4", "PG5", "PG6", "PGE", "P33",

    # PEG fragments & glymes
    "PEG", "1PG", "2PG", "3PG", "TEG", "P4G", "P6G", "P8G", "P10", "P20",

    # sugars & glycans
    "NAG", "NDG", "BMA", "MAN", "GLC", "GAL", "FUC", "SIA", "SLB", "BGC", "BOG", "TRE", "SUC", "MAL", "LMT", "XYP", "ARA", "RIB", "FRU",

    # lipids / fatty acids / sterols / detergents
    "CLR", "CHL", "OLA", "OLE", "STE", "PLM", "MYR", "PAL", "LDA", "LDAO", "DOD", "C8E", "C10E", "C12", "C14",

    # nucleotides / cofactors (not benchmarking small-molecule controls)
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "CTP", "CDP", "CMP", "UTP", "UDP", "UMP", "TTP", "TDP", "TMP",
    "NAD", "NAP", "NADH", "NADP", "FAD", "FMN", "PLP", "TPP", "COA", "ACP", "SAM", "SAH", "THF", "FOL",

    # porphyrins / quinones / retinal
    "HEM", "HEC", "HEB", "HEA", "UQ1", "UQ2", "UQ5", "MEN", "RET", "BCR",

    # detergents
    "CHO", "CHX", "CHD", "TRI", "CHAP", "CHAPSO", "DDM", "DM", "UDM", "MNG", "LMG", "NG",
}
EXCLUDE_HET_IDS |= {"2PE"}  # often phenylethanol/PEG fragment in co-crystals

EXCLUDE_HET_NAME_KEYWORDS: set = {
    "SULFATE", "PHOSPHATE", "CHLORIDE", "BROMIDE", "IODIDE", "NITRATE", "CARBONATE", "BICARBONATE",
    "THIOCYANATE", "CACODYLATE", "ACETATE", "FORMATE", "CITRATE", "TARTRATE", "IMIDAZOLE", "AMMONIUM",
    "ETHYLENE GLYCOL", "DIETHYLENE GLYCOL", "TRIETHYLENE GLYCOL", "GLYCEROL", "MPD", "ISOPROPANOL",
    "ETHANOL", "METHANOL", "DMSO", "DMF", "PEG", "POLYETHYLENE GLYCOL", "GLYME",
    "GLUCOSE", "GALACTOSE", "MANNOSE", "FRUCTOSE", "MALTOSE", "TREHALOSE", "SUCROSE", "N-ACETYLGLUCOSAMINE",
    "DODECYL", "MALTOSIDE", "NEOPENTYL GLYCOL", "DIGITONIN", "TWEEN", "TRITON",
    "OCTYLGLUCOSIDE", "LAURYLDIMETHYLAMINE-OXIDE", "LDAO", "CHOLESTEROL", "OLEATE", "PALMITATE", "STEARATE",
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "NAD", "NADP", "FAD", "FMN", "COENZYME A",
    "S-ADENOSYLMETHIONINE", "PYRIDOXAL PHOSPHATE", "THIAMINE PYROPHOSPHATE",
    "HEME", "UBIQUINONE", "RETINAL",
    "SCHEMBL", "CARBOXYLATE", "CARBOXAMIDE", "PROPANOATE", "BENZAMIDE",
    "GUANIDIN", "IMINIUM", "PIPERIDINIUM", "PYRIDINIUM", "AZONIA", "OXAN",
}


# =============================================================================
# Per-PDB manual hints and hard-coded FDA controls
# =============================================================================

PER_PDB_HINTS: Dict[str, List[str]] = {
    "1M17": ["erlotinib", "gefitinib", "afatinib", "osimertinib"],
    "2E2B": ["zidovudine", "azt", "acyclovir", "ganciclovir"],
    "2GQG": ["sunitinib", "sorafenib", "pazopanib"],
    "2ITN": ["sunitinib", "pazopanib", "regorafenib"],
    "2RGC": ["sotorasib", "adagrasib"],
    "3OG7": ["dasatinib", "ponatinib", "sorafenib"],
    "3QX3": ["etoposide", "doxorubicin", "topotecan", "irinotecan"],
    "4RT7": ["quizartinib", "gilteritinib", "midostaurin"],
    "4XUF": ["quizartinib", "gilteritinib", "midostaurin"],
    "5C7X": ["dinoprostone", "misoprostol"],
    "6ADQ": ["atovaquone"],
    "6GQO": ["axitinib", "sunitinib", "pazopanib", "cabozantinib", "regorafenib"],
    "6O0L": ["venetoclax", "navitoclax"],
}

HARD_FDA_CONTROL_BY_PDB: Dict[str, List[str]] = {
    "1M17": ["erlotinib", "gefitinib", "afatinib", "osimertinib"],
    "2E2B": ["zidovudine", "acyclovir"],
    "2GQG": ["sunitinib", "sorafenib"],
    "2ITN": ["sunitinib", "pazopanib"],
    "2RGC": ["sotorasib", "adagrasib"],
    "3OG7": ["dasatinib", "ponatinib"],
    "3QX3": ["etoposide", "doxorubicin"],
    "4RT7": ["quizartinib", "gilteritinib"],
    "4XUF": ["quizartinib", "gilteritinib"],
    "5C7X": ["dinoprostone"],
    "6ADQ": ["atovaquone"],
    "6GQO": ["axitinib", "sunitinib"],
    "6O0L": ["venetoclax"],
}


# =============================================================================
# ChemComp alias map (expanded)
# =============================================================================

CHEMCOMP_ALIAS: Dict[str, List[str]] = {
    # ABL / KIT / VEGFR TKIs
    "STI": ["imatinib", "gleevec", "sti571"],
    "NIL": ["nilotinib", "tasigna"],
    "ABL": ["asciminib", "abl001", "scemblix"],
    "DAS": ["dasatinib", "sprycel"],
    "BOS": ["bosutinib", "bosulif"],
    "PON": ["ponatinib", "iclusig"],
    "AXI": ["axitinib", "inlyta"],
    "SFB": ["sorafenib", "nexavar", "bay 43-9006", "bay439006"],
    "SU1": ["sunitinib", "sutent", "su11248"],
    "PAZ": ["pazopanib", "votrient"],
    "CAB": ["cabozantinib", "cabometyx", "cometriq"],
    "REG": ["regorafenib", "stivarga"],
    "VAN": ["vandetanib", "zd6474", "caprelsa"],

    # EGFR / ERBB2
    "ERL": ["erlotinib", "tarceva"],
    "GEF": ["gefitinib", "iressa"],
    "AFN": ["afatinib", "gilotrif"],
    "OSM": ["osimertinib", "tagrisso", "azd9291"],
    "LAP": ["lapatinib", "tykerb"],

    # ALK/ROS1/MET/RET
    "CRZ": ["crizotinib", "xalkori"],
    "CER": ["ceritinib", "zykadia"],
    "ALE": ["alectinib", "alecenza"],
    "LOR": ["lorlatinib", "lorbrena", "lorviqua"],
    "BRG": ["brigatinib", "ap26113"],
    "ENT": ["entrectinib", "rxdx-101"],
    "CAP": ["capmatinib", "tabrecta", "capecitabine", "xeloda"],
    "SELr": ["selpercatinib", "rxdx-105", "ret inhibitor"],
    "PRT": ["pralsetinib", "blud-667", "gavripranib", "gprc"],

    # RAS/RAF/MEK
    "VEM": ["vemurafenib", "zelboraf"],
    "DAB": ["dabrafenib", "tafinlar"],
    "ENC": ["encorafenib", "braftovi"],
    "COB": ["cobimetinib", "cotellic"],
    "BIN": ["binimetinib", "mektovi"],
    "TRM": ["trametinib", "mekinist"],
    "SEL": ["selumetinib", "koselugo"],

    # JAK
    "RUX": ["ruxolitinib", "jakafi"],
    "TOF": ["tofacitinib", "xeljanz"],
    "BAR": ["baricitinib", "olumiant"],
    "UPA": ["upadacitinib", "rinvoq"],
    "FED": ["fedratinib", "inoma", "indra", "indra-280"],

    # PI3K/mTOR
    "IDA": ["idelalisib", "zydelig"],
    "DUV": ["duvelisib", "copiktra"],
    "COP": ["copanlisib", "aliqopa"],
    "API": ["alpelisib", "piqray", "byl719"],
    "EVR": ["everolimus", "afinitor"],
    "TMS": ["temsirolimus", "torisel"],
    "RAP": ["rapamycin", "sirolimus"],

    # CDK4/6
    "P31": ["palbociclib", "pd-0332991", "ibrance"],
    "RIB": ["ribociclib", "lee011", "kiskali"],
    "ABE": ["abemaciclib", "ly2835219", "verzenio"],

    # BCL2 / apoptosis
    "ABT": ["venetoclax", "abt-199", "venclexta"],
    "NAV": ["navitoclax", "abt-263"],

    # BTK
    "IBR": ["ibrutinib", "imbruvica"],
    "ACB": ["acalabrutinib", "calquence"],
    "ZAN": ["zanubrutinib", "brukinsa"],

    # FLT3
    "QUI": ["quizartinib", "ac220", "vantictumab"],
    "GIL": ["gilteritinib", "asp2215", "xospata"],
    "CRE": ["crenolanib", "cp-868596"],
    "MID": ["midostaurin", "pkc412", "rydapt"],
    "LST": ["lestaurtinib", "cep-701"],

    # IDH
    "ENA": ["enasidenib", "ag-221", "idhifa"],
    "IVO": ["ivosidenib", "ag-120", "tibsovo"],

    # Hedgehog
    "GLB": ["glasdegib", "pf-04449913", "daurismo"],
    "VIS": ["vismodegib", "erivedge", "gdc-0449"],
    "SON": ["sonidegib", "odenzo", "lde225"],

    # HMAs / cytotoxics
    "AZA": ["azacitidine", "vidaza"],
    "DAC": ["decitabine", "dacogen"],
    "ATO": ["arsenic trioxide", "trisenox"],
    "DNR": ["daunorubicin"],
    "IDR": ["idarubicin"],
    "DXR": ["doxorubicin", "adriamycin"],
    "ETO": ["etoposide", "vp-16"],
    "TPT": ["topotecan"],
    "IRI": ["irinotecan", "cpt-11"],

    # PARP inhibitors
    "OLP": ["olaparib", "lynparza"],
    "NIR": ["niraparib", "zejula"],
    "RUC": ["rucaparib", "rubraca"],
    "TLZ": ["talazoparib", "talzenna"],
    "VLP": ["veliparib", "abt-888"],

    # HDAC inhibitors
    "48D": ["vorinostat", "saha", "zolinza"],
    "PNB": ["panobinostat", "farydak"],
    "BEL": ["belinostat", "beleodaq"],
    "ROM": ["romidepsin", "istodax"],

    # HRT / ER / AR axis
    "TAM": ["tamoxifen"],
    "OHT": ["4-hydroxytamoxifen", "hydroxytamoxifen", "endoxifen", "tamoxifen"],
    "BAX": ["bazedoxifene", "conbriza", "duavive"],
    "FUL": ["fulvestrant", "faslodex"],
    "E2": ["estradiol", "17beta-estradiol", "estrogen"],
    "EST": ["estradiol", "estrogen"],
    "E1": ["estrone"],
    "DHT": ["dihydrotestosterone", "androstanolone"],
    "TES": ["testosterone"],
    "PRG": ["progesterone"],
    "LET": ["letrozole", "femara"],
    "ANA": ["anastrozole", "arimidex"],
    "EXE": ["exemestane", "aromasin"],
    "BIC": ["bicalutamide", "casodex"],
    "ENZ": ["enzalutamide", "xtandi"],
    "APA": ["apalutamide", "erleada"],
    "DAR": ["darolutamide", "nubeqa"],

    # antimetabolites
    "MTX": ["methotrexate"],
    "5FU": ["5-fluorouracil", "fluorouracil"],
    "GEM": ["gemcitabine", "gemzar"],
    "FLUa": ["fludarabine", "f-ara-a"],
    "CLD": ["cladribine", "2-cda"],

    # proteasome
    "BOR": ["bortezomib", "velcade"],
    "CFZ": ["carfilzomib", "kyprolis"],
    "IXA": ["ixazomib", "ninlaro"],

    # HIV antivirals
    "RTV": ["ritonavir"],
    "LPV": ["lopinavir"],
    "ATV": ["atazanavir"],
    "EFV": ["efavirenz"],

    # others
    "MET": ["metformin"],
    "DXN": ["dexamethasone"],
    "CPT": ["camptothecin", "topotecan", "irinotecan"],

    # endocrine/other (alt key)
    "EVE": ["everolimus", "afinitor"],
}
CHEMCOMP_ALIAS.update({
    "0LI": ["ponatinib", "iclusig"],  # 3ZOS
    "1BQ": ["ensartinib"],            # 4I4E
    "69Q": ["enasidenib"],            # 5I96
    "B49": ["xenazine"],              # 3G0E
    "BRL": ["orbenin"],               # 3DZY
    "C6F": ["lampren"],               # 6JQR
    "CXS": ["lampren"],               # 6JQR
    "FLC": ["pempidine"],             # 5TQH
    "LQQ": ["palbociclib isethiolate", "palbociclib"],  # 5L2I
    "MI1": ["xeljanz", "tofacitinib"],                  # 3LXK
    "P06": ["valtrex"],                                  # 4XV2
    "RXT": ["jakafi", "ruxolitinib"],                    # 4U5J, 6WTN
    "VGH": ["crizotinib"],                               # 2WGJ, 2XP2, 3ZBF
})
CHEMCOMP_ALIAS.update({
    "LEV": ["lenvatinib", "lenvima"],
    "ANT": ["antimycin a"],
    "STG": ["stigmatellin"],
    "ATO": ["atovaquone"],
    "F0K": ["bi-2852"],  # 6GJ8
})
CHEMCOMP_ALIAS.update({
    "SOT": ["sotorasib", "lumakras"],
    "ADA": ["adagrasib", "krazati"],
    "ACV": ["acyclovir", "zovirax"],
    "ZDV": ["zidovudine", "azt"],
})
for _k, _vals in list(CHEMCOMP_ALIAS.items()):
    CHEMCOMP_ALIAS[_k] = sorted(set(v.lower() for v in _vals))


# =============================================================================
# Alias harvesting from benchmark_analysis_details.csv
# =============================================================================

def _looks_like_brand_or_generic(name: str) -> bool:
    n = (name or "").strip()
    if len(n) < 3 or len(n) > 64:
        return False
    if not re.search(r"[a-zA-Z]", n):
        return False
    noisy = len(re.findall(r"[0-9\[\]\(\)\/\.\-\,;:]", n))
    return (noisy / max(1, len(n))) < 0.40

def _clean_alias(name: str) -> Optional[str]:
    if not name:
        return None
    s = name.strip()
    if s.startswith("?"):
        s = s.lstrip("?\uFF1F").strip()
    s = s.replace("\u2019", "'").replace("\u00AE", "").replace("\u2122", "")
    s = re.sub(r"\s+", " ", s)
    s_low = s.lower()
    if _looks_like_brand_or_generic(s):
        return s_low
    return None

def load_aliases_from_details_csv(paths: Sequence[str]) -> Dict[str, List[str]]:
    exp: Dict[str, set] = {}
    for csv_path in _expand_globs(paths):
        if not csv_path.is_file():
            continue
        try:
            with open(csv_path, "r", encoding="utf-8", errors="ignore", newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    het = (row.get("control_het") or "").strip().upper()
                    nm = _clean_alias(row.get("rdk_name") or "")
                    if not het or not nm:
                        continue
                    exp.setdefault(het, set()).add(nm)
        except Exception as e:
            print(f"[alias] WARN: failed to parse {csv_path}: {e}")
    return {k: sorted(v) for k, v in exp.items()}

def merge_aliases_into_chemcomp(
    base: Dict[str, List[str]],
    extra: Dict[str, List[str]],
) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {k: sorted(set(vv.lower() for vv in vs)) for k, vs in base.items()}
    for het, names in (extra or {}).items():
        merged = set(out.get(het, []))
        for nm in names:
            nm_norm = (nm or "").strip().lower()
            if nm_norm:
                merged.add(nm_norm)
        out[het] = sorted(merged)
    return out


# =============================================================================
# PDB hint parsing
# =============================================================================

def parse_pdb_het_hints(pdb_path: Path) -> Tuple[List[str], List[str]]:
    pdb_path = Path(pdb_path)
    het_ids: set = set()
    hetnam_map: Dict[str, List[Tuple[int, str]]] = defaultdict(list)

    if not pdb_path.is_file():
        return [], []

    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for raw in fh:
                rec = raw[:6].strip().upper()

                if rec == "HET":
                    het = (raw[7:10].strip() if len(raw) >= 10 else "").upper()
                    if not het:
                        toks = raw.split()
                        if len(toks) >= 2:
                            het = toks[1].upper()
                    if het:
                        het_ids.add(het)

                elif rec == "HETNAM":
                    toks = raw.split()
                    cont = 1
                    het, name = "", ""
                    if len(toks) >= 3 and toks[1].isdigit():
                        cont = int(toks[1])
                        het = toks[2].upper()
                        name = " ".join(toks[3:])
                    elif len(toks) >= 2:
                        het = toks[1].upper()
                        name = " ".join(toks[2:])
                    if het and name:
                        hetnam_map[het].append((cont, name.strip()))

                elif rec == "HETATM":
                    resn = raw[17:20].strip().upper()  # residue name columns 18–20
                    if resn:
                        het_ids.add(resn)
    except Exception:
        return [], []

    STANDARD_RES = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
        "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "PYL"
    }

    het_ids = {
        h for h in het_ids
        if h not in EXCLUDE_HET_IDS and h not in STANDARD_RES and 2 <= len(h) <= 5
    }

    het_names: List[str] = []
    for het in sorted(het_ids):
        if het in hetnam_map:
            parts = [t for _, t in sorted(hetnam_map[het], key=lambda x: x[0])]
            nm = " ".join(parts).strip()
            if not nm:
                continue
            nm_up = nm.upper()
            if any(kw in nm_up for kw in EXCLUDE_HET_NAME_KEYWORDS):
                continue
            cleaned = _clean_alias(nm)
            het_names.append(cleaned if cleaned else nm)

    seen: set = set()
    dedup_names: List[str] = []
    for n in het_names:
        nn = _norm(n)
        if nn and nn not in seen:
            seen.add(nn)
            dedup_names.append(n)

    return sorted(het_ids), dedup_names


# =============================================================================
# Mapping index (reads your fda_mapping_from_pdbqt.csv)
# =============================================================================

@dataclass
class MappingRow:
    path: str
    display_name: str = ""
    generic_name: str = ""
    brand_names: str = ""
    pubchem_name: str = ""
    pubchem_record_title: str = ""
    pubchem_iupac_name: str = ""
    pubchem_synonyms: str = ""
    rxnorm_generic_name: str = ""
    rxnorm_brand_names: str = ""
    drugcentral_generic_name: str = ""
    drugcentral_brand_names: str = ""
    remark_name: str = ""
    sdf_title: str = ""
    inchikey: str = ""

    def all_name_fields(self) -> List[Tuple[str, str]]:
        return [
            ("display_name", self.display_name),
            ("generic_name", self.generic_name),
            ("brand_names", self.brand_names),
            ("pubchem_name", self.pubchem_name),
            ("pubchem_record_title", self.pubchem_record_title),
            ("pubchem_iupac_name", self.pubchem_iupac_name),
            ("pubchem_synonyms", self.pubchem_synonyms),
            ("rxnorm_generic_name", self.rxnorm_generic_name),
            ("rxnorm_brand_names", self.rxnorm_brand_names),
            ("drugcentral_generic_name", self.drugcentral_generic_name),
            ("drugcentral_brand_names", self.drugcentral_brand_names),
            ("remark_name", self.remark_name),
            ("sdf_title", self.sdf_title),
        ]

class MappingIndex:
    """In-memory index over the mapping CSV for efficient candidate selection."""
    def __init__(self, csv_path: Path):
        import pandas as pd
        self.csv_path = Path(csv_path)

        # Helpful error if Windows path on Linux
        if not self.csv_path.exists():
            guess = _resolve_existing_path(str(self.csv_path))
            if guess != self.csv_path and guess.exists():
                print(f"[mapping] Adjusted mapping path -> {guess}")
                self.csv_path = guess

        if not self.csv_path.exists():
            raise FileNotFoundError(
                f"Mapping CSV not found: {self.csv_path}\n"
                "Tip: pass --mapping with a Linux/WSL path (e.g. /mnt/e/...)\n"
                "or keep Windows path and this script will try to convert."
            )

        try:
            self.df = pd.read_csv(self.csv_path)
        except Exception as e:
            raise RuntimeError(f"Failed to read mapping CSV '{self.csv_path}': {e}")

        if "path" not in self.df.columns:
            raise ValueError(
                f"Mapping CSV must include a 'path' column to .pdbqt files. Columns present: {list(self.df.columns)}"
            )

        self.rows: List[MappingRow] = []
        for _, r in self.df.iterrows():
            self.rows.append(MappingRow(
                path=str(r.get("path", "")),
                display_name=str(r.get("display_name", "")),
                generic_name=str(r.get("generic_name", "")),
                brand_names=str(r.get("brand_names", "")),
                pubchem_name=str(r.get("pubchem_name", "")),
                pubchem_record_title=str(r.get("pubchem_record_title", "")),
                pubchem_iupac_name=str(r.get("pubchem_iupac_name", "")),
                pubchem_synonyms=str(r.get("pubchem_synonyms", "")),
                rxnorm_generic_name=str(r.get("rxnorm_generic_name", "")),
                rxnorm_brand_names=str(r.get("rxnorm_brand_names", "")),
                drugcentral_generic_name=str(r.get("drugcentral_generic_name", "")),
                drugcentral_brand_names=str(r.get("drugcentral_brand_names", "")),
                remark_name=str(r.get("remark_name", "")),
                sdf_title=str(r.get("sdf_title", "")),
                inchikey=str(r.get("inchikey", "")),
            ))

    def rows_by_rdk_id(self, rdk_id: str) -> List[MappingRow]:
        rid = (rdk_id or "").strip().lower()
        if not rid:
            return []
        out: List[MappingRow] = []
        for row in self.rows:
            stem = Path(row.path).stem.lower()
            if rid in stem:
                out.append(row)
        return out

    def search(
        self,
        hints: Sequence[str],
        inchikey: Optional[str] = None,
        max_results: int = 6,
    ) -> List[Tuple[MappingRow, int, str]]:
        hints_norm = [_norm(h) for h in hints if h]
        hint_tokens = set(t for h in hints_norm for t in _tokenize(h))
        raw_het_codes = {h.strip().upper() for h in hints if h and 2 <= len(h.strip()) <= 5}

        out: List[Tuple[MappingRow, int, str]] = []
        for row in self.rows:
            # Normalize each mapping path similarly (Windows->WSL) then existence check
            p = _resolve_existing_path(row.path)
            if not str(p):
                continue
            if not p.exists():
                # One more try: if path was relative to CSV folder
                p_rel = (self.csv_path.parent / row.path).resolve()
                if p_rel.exists():
                    p = p_rel
                else:
                    continue

            best = 0
            why = ""
            stem_upper = p.stem.upper()
            base_upper = p.name.upper()
            parents_upper = " ".join([pp.name.upper() for pp in p.parents])

            # 0) InChIKey exact
            if inchikey and row.inchikey and row.inchikey.strip().upper() == inchikey.strip().upper():
                best, why = 100, "inchikey_exact"

            # 1) Boost if filename/parents encode a 3–5 letter HET code
            if best < 100 and raw_het_codes:
                for code in raw_het_codes:
                    if re.search(rf"\b{re.escape(code)}\b", stem_upper) or re.search(rf"\b{re.escape(code)}\b", parents_upper):
                        sc = 93; rs = f"path_token:{code}"
                    elif code in base_upper:
                        sc = 88; rs = f"path_substr:{code}"
                    else:
                        sc = 0; rs = ""
                    if sc > best:
                        best, why = sc, rs

            # 2) Name-field scoring
            if best < 100 and hints_norm:
                for field, value in row.all_name_fields():
                    v = _norm(value)
                    if not v:
                        continue
                    if v in hints_norm:
                        sc = 95; rs = f"{field}_exact"
                    elif any(h in v for h in hints_norm if len(h) >= 3):
                        sc = 85; rs = f"{field}_substr"
                    else:
                        vtok = set(_tokenize(v))
                        overlap = len(vtok & hint_tokens)
                        sc = 70 if overlap >= 2 else (65 if overlap == 1 else 0)
                        rs = f"{field}_tokens:{overlap}"
                    if sc > best:
                        best, why = sc, rs

            if best > 0:
                # Store the resolved absolute path back into row.path for downstream use
                row.path = str(Path(p).resolve())
                out.append((row, best, why))

        out.sort(key=lambda t: t[1], reverse=True)
        return out[:max_results]


# =============================================================================
# Candidate & control selection helpers
# =============================================================================

def _resolve_name_for_path_or_text(p_or_text: str, fda_index) -> str:
    rdk = _extract_rdk_id(p_or_text or "")
    if rdk:
        name = resolve_corresponding_name_for_rdk(rdk, fda_index) or ""
        if name:
            return name
    stem = Path(p_or_text).stem if os.path.exists(p_or_text) else str(p_or_text)
    name = resolve_corresponding_name_from_text(stem, fda_index) or ""
    return name

def select_candidates_for_protein(
    mapping: MappingIndex,
    hints: Sequence[str],
    prepped_dir: Path,
    extra_parent_ids: Sequence[str],
    max_candidates: int,
) -> List[Tuple[MappingRow, int, str]]:
    candidates = mapping.search(hints, inchikey=None, max_results=max(6, max_candidates * 4))

    parent_rows: List[Tuple[MappingRow, int, str]] = []
    for pid in extra_parent_ids or []:
        rid = _extract_rdk_id(pid or "")
        if rid:
            for r in mapping.rows_by_rdk_id(rid):
                parent_rows.append((r, 99, "parent_from_metabolite"))
        else:
            for (r, sc, why) in mapping.search([pid], inchikey=None, max_results=2):
                parent_rows.append((r, max(sc, 92), f"{why}|parent_from_metabolite"))

    best_by_path: Dict[str, Tuple[MappingRow, int, str]] = {}
    for row, sc, why in (candidates + parent_rows):
        p_abs = str(Path(row.path).resolve())
        if (p_abs not in best_by_path) or (sc > best_by_path[p_abs][1]):
            best_by_path[p_abs] = (row, sc, why)

    prepped_dir_resolved = prepped_dir.resolve()
    cand_rows: List[Tuple[MappingRow, int, str]] = []
    for row, sc, why in best_by_path.values():
        p = Path(row.path)
        if p.is_file() and _path_is_within(p, prepped_dir_resolved):
            cand_rows.append((row, sc, why))

    cand_rows.sort(key=lambda t: t[1], reverse=True)
    return cand_rows[:max_candidates]

def split_controls_and_whitelist(
    whitelist_paths: Sequence[str],
    prepped_control_pdbqts: Sequence[str],
    control_stems_lower: set,
    heavy_atom_counts: Dict[str, int],
    cfg: Dict,
    logger,
) -> Tuple[List[str], List[str]]:
    ctrl_blacklist = {t.strip().upper() for t in str(cfg.get("CONTROL_BLACKLIST", "")).split(",") if t.strip()}
    min_ha = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))

    def _is_control_path(p: str) -> bool:
        stem = Path(p).stem.split("_stage")[0]
        if stem.upper() in ctrl_blacklist:
            return False
        if stem.lower() not in control_stems_lower:
            return False
        ha = heavy_atom_counts.get(p)
        return (ha is None) or (ha >= min_ha)

    ctrls = [p for p in prepped_control_pdbqts if _is_control_path(p)]
    non_ctrls = [p for p in whitelist_paths if not _is_control_path(p)]

    ctrls = _dedupe_str(ctrls)
    non_ctrls = _dedupe_str(non_ctrls)

    max_ctrls = int(cfg.get("BENCH_MAX_CONTROLS", 8))
    if len(ctrls) > max_ctrls:
        if logger:
            logger.warning(f"[SANITY] controls={len(ctrls)} looks high; capping to {max_ctrls}. Extra will be skipped.")
            for p in ctrls[max_ctrls:]:
                logger.warning(f"[SANITY-OFFENDER] {p}")
        ctrls = ctrls[:max_ctrls]

    return ctrls, non_ctrls

def _merge_per_pdb_hints(pdb_id: str, hints: List[str]) -> List[str]:
    out = list(hints or [])
    add: List[str] = []
    u = pdb_id.upper()
    if u in PER_PDB_HINTS:
        add.extend(PER_PDB_HINTS[u])
    if u in HARD_FDA_CONTROL_BY_PDB:
        add.extend(HARD_FDA_CONTROL_BY_PDB[u])

    seen: set = set()
    merged: List[str] = []
    for h in (out + add):
        hn = _norm(h)
        if hn and hn not in seen:
            seen.add(hn)
            merged.append(h)
    return merged

def _promote_forced_controls_to_ctrls(
    pdb_id: str,
    mapping: MappingIndex,
    prepped_dir: Path,
    names: Sequence[str],
    logger=None,
) -> Tuple[List[str], List[str]]:
    if not names:
        return [], []
    rows = mapping.search(names, inchikey=None, max_results=max(8, len(names) * 3))
    prepped_dir_resolved = prepped_dir.resolve()
    promoted_paths: List[str] = []
    promoted_stems_lower: List[str] = []

    for row, sc, why in rows:
        p = Path(row.path)
        if p.is_file() and p.resolve().parent == prepped_dir_resolved:
            promoted_paths.append(str(p.resolve()))
            stem0 = p.stem.split("_stage")[0].lower()
            promoted_stems_lower.append(stem0)
            if logger:
                logger.info(f"[forced-control] {pdb_id}: promoting '{stem0}' as control (reason={why}, score={sc})")

    promoted_paths = _dedupe_str(promoted_paths)
    promoted_stems_lower = _dedupe_str(promoted_stems_lower)
    if logger and not promoted_paths:
        logger.warning(f"[forced-control] {pdb_id}: requested {list(names)} but none found in '{prepped_dir}'.")
    return promoted_paths, promoted_stems_lower

def collect_prepped_controls_for_protein(
    paths,
    control_stems: Sequence[str],
    logger=None,
) -> Tuple[List[str], set]:
    control_stems_lower = {s.lower() for s in control_stems if s}
    prepped_control_pdbqts: List[str] = []
    prepped_dir_resolved = paths.prepped_ligands_dir.resolve()

    if paths.prepped_ligands_dir.exists():
        for p in paths.prepped_ligands_dir.glob("*.pdbqt"):
            stem0 = p.stem.split("_stage")[0].lower()
            if stem0 in control_stems_lower and p.parent.resolve() == prepped_dir_resolved:
                prepped_control_pdbqts.append(str(p.resolve()))

    prepped_control_pdbqts = sorted(set(prepped_control_pdbqts))

    if logger:
        leaks = [pp for pp in prepped_control_pdbqts
                 if Path(pp).resolve().parent != prepped_dir_resolved]
        if leaks:
            logger.warning(f"[CONTROL-LEAK] Found {len(leaks)} controls outside {prepped_dir_resolved}")
            for pp in leaks[:10]:
                logger.warning(f"  leak -> {pp}")

    return prepped_control_pdbqts, control_stems_lower


# =============================================================================
# Benchmark driver — single ultra-stage per pocket
# =============================================================================

def run_benchmark_for_protein(
    cfg: Dict,
    mapping: MappingIndex,
    pdb_file: str,
    prepped_dir: Path,
    out_root: Path,
    exhaustiveness: int,
    num_modes: int,
    max_candidates: int,
    manual_hints: Optional[List[str]] = None,
    fda_index=None,
) -> None:
    import numpy as np  # used a few times; import once here
    import logging

    t0 = time.time()
    base_id = os.path.splitext(pdb_file)[0]
    pdb_id = base_id.replace("_cleaned", "").upper()
    logger = make_protein_logger(cfg["DOCKED_DIR"], pdb_id, cfg)

    # Keep Windows consoles readable if they don't like unicode
    def _sanitize_msg(s: str) -> str:
        return (s.replace("≤", "<=").replace("≥", ">=").replace("Å", " Angstrom")
                 .replace("µ", "u").replace("°", " deg"))

    class _AsciiFilter(logging.Filter):
        def filter(self, record):
            if isinstance(record.msg, str):
                record.msg = _sanitize_msg(record.msg)
            return True

    class _ConsoleFilter(logging.Filter):
        """
        Allow: progress & docking notes.
        Drop: Reduce/phenix/OpenBabel spam, histograms, H-Scan, chain breaks, 'H atom too close', etc.
        In --silent, also drop non-allowlisted WARNING (ERRORs still pass).
        """
        ALLOW_PREFIX = (
            "[benchmark]", "[pocket]", "[stage]", "[wave]", "[pools]",
            "[CENTER]", "[CONTROL-LOCK]", "[recenter]", "[shrink]",
            "[Checkpoint]", "[forced-control]", "[metabolite→parent]"
        )
        DROP_SUBSTR = (
            "phenix.pdbtools completed successfully",
            "appear unbonded and will be treated as a chain break",
            "H atom too close",
            "Open Babel hydrogenation succeeded",
            "Open Babel used to add hydrogens.",
            "[Elem histogram", "[H-Scan before]", "[H-Scan after ]",
            "Skipping invalid chain",
        )
        def __init__(self, quiet: bool, silent: bool):
            super().__init__()
            self.quiet = quiet or silent
            self.silent = silent
        def filter(self, record: logging.LogRecord) -> bool:
            msg = str(record.msg)
            lvl = record.levelno
            if lvl >= logging.ERROR:
                return True
            if any(x in msg for x in self.DROP_SUBSTR):
                return False
            if msg.startswith(self.ALLOW_PREFIX):
                return True
            if not self.quiet:
                return True
            if lvl == logging.WARNING:
                return (not self.silent)
            return False

    # Attach filters/levels to handlers: keep files verbose, console filtered
    for h in logger.handlers:
        h.addFilter(_AsciiFilter())
        if isinstance(h, logging.FileHandler):
            h.setLevel(logging.DEBUG)
        else:
            # Treat non-file as console-like
            h.setLevel(logging.INFO)
            quiet = bool(cfg.get("_QUIET_CONSOLE", True))
            silent = bool(cfg.get("_SILENT_CONSOLE", False))
            h.addFilter(_ConsoleFilter(quiet=quiet, silent=silent))

    logger.info(f"[benchmark] Starting {pdb_id}")

    paths = make_paths(cfg, base_id, pdb_file)
    logger.debug(f"[paths] pdb_path={paths.pdb_path}")
    logger.debug(f"[paths] ligand_output_dir={paths.ligand_output_dir}")
    logger.debug(f"[paths] prepped_ligands_dir={paths.prepped_ligands_dir}")

    # 1) Extract controls & build nolig
    _, control_stems = extract_ligands_to_nolig(paths, logger)
    robust_prepare_controls(paths, cfg, logger)

    # 1b) If co-crystal is a metabolite, add parent drug(s) from library
    extra_parent_ids: List[str] = []
    try:
        extra_parent_ids = ensure_parent_drugs_for_controls(
            pdb_code=pdb_id,
            ligands_raw_dir=str(paths.ligand_output_dir),
            fda_index=fda_index,
            max_additions=3,
        ) or []
        if extra_parent_ids:
            logger.info(f"[metabolite→parent] parents: {', '.join(extra_parent_ids)}")
    except Exception as e:
        logger.warning(f"[metabolite→parent] resolver failed: {e}")

    # 2) Prepare receptor
    cleaned_pdb, receptor_pdbqt = prepare_receptor(cfg, paths, logger)
    if not cleaned_pdb or not receptor_pdbqt:
        logger.warning("[benchmark] receptor prep failed; skipping protein")
        return

    # 3) Detect pocket via main logic (with fallback)
    center, detected_box, src = detect_pocket(cleaned_pdb, paths.ligand_output_dir, logger)
    logger.info(f"[pocket] source={src} center={center} box={detected_box}")
    if center:
        pockets = [("pocket1", center)]
        pocket_box_size = tuple(min(28.0, float(s)) for s in (detected_box or (24.0, 24.0, 24.0)))
    else:
        logger.error("[benchmark] no pocket could be detected; skipping protein")
        return
    if not pockets:
        c, _, _ = detect_active_site(cleaned_pdb)
        if c:
            pockets = [("pocket1", c)]
        else:
            logger.error("[benchmark] no pocket could be detected; skipping protein")
            return

    # 4) Build hints (HET/HETNAM + aliases + manual + metabolite parents)
    het_ids, het_names = parse_pdb_het_hints(Path(paths.pdb_path))
    hints: List[str] = []
    hints.extend(het_names)
    hints.extend(het_ids)
    for het in het_ids:
        hints.extend(CHEMCOMP_ALIAS.get(het.upper(), []))
    if manual_hints:
        hints.extend(manual_hints)
    hints.extend(extra_parent_ids)

    seen = set()
    final_hints: List[str] = []
    for h in hints:
        hn = _norm(h)
        if hn and hn not in seen:
            final_hints.append(h)
            seen.add(hn)
    final_hints = _merge_per_pdb_hints(pdb_id, final_hints)
    logger.debug(f"[hints] {final_hints}")

    # 5) Candidate selection (restricted to prepped_dir)
    cand_rows = select_candidates_for_protein(
        mapping=mapping,
        hints=final_hints,
        prepped_dir=prepped_dir,
        extra_parent_ids=extra_parent_ids,
        max_candidates=max_candidates,
    )
    cand_rows.sort(key=lambda t: t[1], reverse=True)
    logger.info(f"[candidates] {len(cand_rows)} selected (top score={cand_rows[0][1] if cand_rows else 'NA'})")

    # Audit CSV
    out_dir = Path(out_root) / pdb_id
    fh, cand_csv_path = _safe_open_csv_for_write(out_dir / f"benchmark_candidates_{pdb_id}.csv")
    with fh:
        w = csv.writer(fh)
        w.writerow(
            ["rank", "score", "reason", "path", "display_name", "generic_name", "brand_names", "inchikey", "corresponding_name"]
        )
        for i, (row, sc, why) in enumerate(cand_rows, start=1):
            corr = _resolve_name_for_path_or_text(row.path, fda_index) or row.generic_name or row.display_name
            w.writerow([i, sc, why, row.path, row.display_name, row.generic_name, row.brand_names, row.inchikey, corr])
    if str(cand_csv_path) != str(out_dir / f"benchmark_candidates_{pdb_id}.csv"):
        logger.info(f"[benchmark] candidates CSV was locked; wrote to {cand_csv_path.name}")

    whitelist = [row.path for (row, _, _) in cand_rows]
    if not whitelist:
        logger.warning("[benchmark] no candidate ligands matched mapping; skipping protein")
        return

    # Heavy-atom counts (needed by CenterSelector / LE)
    heavy_atom_counts: Dict[str, int] = {}
    for p in whitelist:
        try:
            heavy_atom_counts[p] = int(_count_heavy_atoms_from_pdbqt(Path(p)))
        except Exception:
            heavy_atom_counts[p] = 0

    # Control lookup (RMSD gates)
    control_lookup = build_control_lookup(paths)

    # Strictly collect THIS protein’s prepped controls and lowercase set
    prepped_control_pdbqts, control_stems_lower = collect_prepped_controls_for_protein(
        paths=paths,
        control_stems=control_stems,
        logger=logger,
    )

    # Promote hard-coded FDA controls (if present) to CONTROL set
    u = pdb_id.upper()
    if u in HARD_FDA_CONTROL_BY_PDB:
        forced_names = HARD_FDA_CONTROL_BY_PDB[u]
        forced_paths, forced_stems_lower = _promote_forced_controls_to_ctrls(
            pdb_id=pdb_id,
            mapping=mapping,
            prepped_dir=prepped_dir,
            names=forced_names,
            logger=logger,
        )
        if forced_paths:
            prepped_control_pdbqts = _dedupe_str(prepped_control_pdbqts + forced_paths)
            control_stems_lower = set(control_stems_lower) | set(forced_stems_lower)

    # Ensure HA counts include controls too
    for p in prepped_control_pdbqts:
        try:
            heavy_atom_counts[p] = int(_count_heavy_atoms_from_pdbqt(Path(p)))
        except Exception:
            heavy_atom_counts[p] = heavy_atom_counts.get(p, 0)

    # -------- single-stage per pocket --------
    pocket_strength_rows: List[List[str]] = []

    for pocket_name, center in pockets[:2]:
        stage = {
            "name": f"bench_{pocket_name}_single",
            "exhaustiveness": int(exhaustiveness),
            "num_modes": int(num_modes),
            "energy_range": int(cfg.get("ENERGY_RANGE", 9)),
        }

        params = get_recenter_params(cfg)
        guard = GlobalCenterGuard(max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2)))
        selector = CenterSelector(
            cfg,
            logger,
            control_stems=set(control_stems),
            heavy_atom_counts=heavy_atom_counts,
            initial_center=center,
        )

        try:
            cx, cy, cz = float(center[0]), float(center[1]), float(center[2])
        except Exception:
            cx, cy, cz = 0.0, 0.0, 0.0  # fallback

        box_size = tuple(min(28.0, float(s)) for s in (cfg.get("BOX_SIZE_OVERRIDE", None) or (24.0, 24.0, 24.0)))
        logger.info(f"[stage] {pdb_id}/{pocket_name} center={center} box={box_size}")

        # Ligand pool (controls first, then whitelist)
        ctrls, non_ctrls = split_controls_and_whitelist(
            whitelist_paths=whitelist,
            prepped_control_pdbqts=prepped_control_pdbqts,
            control_stems_lower=control_stems_lower,
            heavy_atom_counts=heavy_atom_counts,
            cfg=cfg,
            logger=logger,
        )
        ligands = _dedupe_str(ctrls + non_ctrls)
        stage1_original = ligands[:]
        logger.info(f"[pools] ctrls={len(ctrls)} whitelist={len(non_ctrls)} total={len(ligands)}")

        score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
        validated_ligands_last: List[str] = []
        recenter_attempts = 0

        # Optional checkpoint skip
        if bool(cfg.get("CHECKPOINT_ENABLE", False)):
            fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
            if checkpoint_should_skip(cfg, pdb_id, stage["name"], fp):
                logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                continue

        # Allow early-recenter restarts within this single-stage loop
        restarting = True
        while restarting:
            restarting = False
            guard.reset_stage()

            if not ligands:
                logger.warning(f"No ligands to dock at {stage['name']}; skipping pocket.")
                break

            logger.info(f"[benchmark] {pdb_id} | {pocket_name} | starting {stage['name']} with {len(ligands)} ligands")

            # Two-wave (controls then whitelist) if both present
            if ctrls and non_ctrls:
                logger.info(f"[wave] controls={len(ctrls)} then whitelist={len(non_ctrls)}")

                s1, v1, d1, rd1, inv1 = run_one_stage(
                    cfg=cfg,
                    pdb_id=pdb_id,
                    receptor_pdbqt=receptor_pdbqt,
                    center=center,
                    box_size=box_size,
                    stage=stage,
                    ligands=ctrls,
                    logger=logger,
                    retry_mgr=RetryManager(),
                    control_lookup=control_lookup,
                )

                # CenterSelector promotion based on controls
                try:
                    dec = selector.consider_switch(stage["name"], s1, v1, rd1, receptor_pdbqt, center, guard)
                    if dec.promoted and dec.new_center is not None:
                        old = center
                        center = dec.new_center
                        guard.mark_switch()
                        if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                            checkpoint_invalidate_from(cfg, pdb_id, [stage], start_index=0)
                        logger.info(f"[CENTER] switched (controls wave): {old} -> {center} ({dec.reason}) [global switch]")
                except Exception as e:
                    logger.warning(f"CenterSelector (controls) failed: {e}")

                # Early control lock
                lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", -6.0))
                lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
                lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))
                qualified_controls = []
                for lig in v1:
                    stem = Path(lig).stem.split("_stage")[0].lower()
                    if stem not in control_stems_lower:
                        continue
                    sc = s1.get(lig)
                    if sc is None or not isinstance(sc, (int, float)):
                        continue
                    if sc > lock_score_max:
                        continue
                    pose_path = rd1.get(lig)
                    if not pose_path:
                        continue
                    c = CenterSelector._pdbqt_centroid(pose_path)
                    if c is None:
                        continue
                    if np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                        qualified_controls.append(lig)
                if len(qualified_controls) >= lock_min_hits and not guard.locked:
                    guard.lock()
                    logger.info(f"[CONTROL-LOCK] n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Å")

                # Wave B — whitelist
                s2, v2, d2, rd2, inv2 = run_one_stage(
                    cfg=cfg,
                    pdb_id=pdb_id,
                    receptor_pdbqt=receptor_pdbqt,
                    center=center,
                    box_size=box_size,
                    stage=stage,
                    ligands=non_ctrls,
                    logger=logger,
                    retry_mgr=RetryManager(),
                    control_lookup=control_lookup,
                )

                scores, validated, distances = ({**s1, **s2}, v1 + v2, d1 + d2)
                raw_docked = {**rd1, **rd2}
                invalids = {**inv1, **inv2}
            else:
                # Single wave
                scores, validated, distances, raw_docked, invalids = run_one_stage(
                    cfg=cfg,
                    pdb_id=pdb_id,
                    receptor_pdbqt=receptor_pdbqt,
                    center=center,
                    box_size=box_size,
                    stage=stage,
                    ligands=ligands,
                    logger=logger,
                    retry_mgr=RetryManager(),
                    control_lookup=control_lookup,
                )

            validated_ligands_last = validated

            # Record (for LE calc, etc.)
            for lig, sc in scores.items():
                record_score(score_history, stage["name"], lig, sc, True)
                record_le(score_history, stage["name"], lig, sc, heavy_atom_counts)
            for lig, (sc, reason) in invalids.items():
                record_score(score_history, stage["name"], lig, sc, False, reason=str(reason))
                record_le(score_history, stage["name"], lig, sc, heavy_atom_counts)

            # CenterSelector consideration after full run
            try:
                decision = selector.consider_switch(
                    stage["name"], scores, validated, raw_docked, receptor_pdbqt, center, guard
                )
                if decision.promoted and decision.new_center is not None:
                    old = center
                    center = decision.new_center
                    guard.mark_switch()
                    if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                        checkpoint_invalidate_from(cfg, pdb_id, [stage], start_index=0)
                    logger.info(f"[CENTER] {pdb_id} {pocket_name}: {old} -> {center} ({decision.reason}) [global switch]")
            except Exception as e:
                logger.warning(f"CenterSelector failed gracefully: {e}")

            # Early recenter/expand (single-stage flavor)
            restart, center, box_size, redo_ligands, recenter_attempts = early_recenter_decision(
                0,
                scores,
                distances,
                box_size,
                center,
                stage1_original,
                recenter_attempts,
                params,
                cfg,
                pdb_id,
                receptor_pdbqt,
                logger,
                raw_docked,
                guard,
                control_anchor_hit=False,
            )
            if restart:
                ligands = redo_ligands
                if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                    checkpoint_invalidate_from(cfg, pdb_id, [stage], start_index=0)
                logger.info(f"[recenter] restarting stage with {len(ligands)} ligands; new center={center}, box={box_size}")
                restarting = True
                continue

            # Adaptive shrink (optional; harmless in single-stage)
            try:
                if cfg.get("ADAPTIVE_SHRINK_ENABLE", True) and validated:
                    med = float(np.median([d for d in distances if isinstance(d, (int, float))])) if distances else None
                    if (med is not None) and (med < float(cfg.get("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0))):
                        dec = float(cfg.get("ADAPTIVE_SHRINK_DEC", 4.0))
                        min_box = float(cfg.get("ADAPTIVE_SHRINK_MIN_BOX", 14.0))
                        new_box = tuple(max(min_box, s - dec) for s in box_size)
                        if new_box != box_size:
                            logger.info(f"[shrink] median {med:.2f} Å -> box {box_size} -> {new_box}")
                            box_size = new_box
            except Exception as _e:
                logger.warning(f"Adaptive shrink skipped: {_e}")

            # Mark checkpoint done for this single stage
            if bool(cfg.get("CHECKPOINT_ENABLE", False)):
                try:
                    fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                    checkpoint_mark_done(cfg, pdb_id, stage["name"], fp)
                except Exception:
                    pass

            # ---- outputs for this pocket ----
            out_dir.mkdir(parents=True, exist_ok=True)
            long_csv = out_dir / f"docking_score_long_{pocket_name}.csv"
            results = score_history.get(stage["name"], {})
            ordered = sorted(
                results.items(),
                key=lambda kv: (
                    kv[1].get("score", 1e9) if isinstance(kv[1].get("score"), (int, float)) else 1e9,
                    os.path.basename(kv[0]),
                ),
            )
            with open(long_csv, "w", newline="", encoding="utf-8") as f:
                f.write(f"# pocket={pocket_name}, center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f})\n")
                w = csv.writer(f)
                w.writerow(["stage", "ligand", "library_id", "corresponding_name", "score", "valid", "reason"])
                for lig, rec in ordered:
                    sc = rec.get("score")
                    sc_str = f"{sc:.2f}" if isinstance(sc, (int, float)) else ""
                    lib_id = _extract_rdk_id(lig) or ""
                    corr = _resolve_name_for_path_or_text(lig, fda_index)
                    w.writerow([stage["name"], os.path.basename(lig), lib_id, corr, sc_str, int(bool(rec.get("valid", False))), str(rec.get("reason", ""))])

            # Final validation/screenshots
            if validated_ligands_last:
                final_pose_validation_and_screenshots(
                    cfg,
                    pdb_id,
                    [stage],
                    receptor_pdbqt,
                    center,
                    validated_ligands_last,
                    score_history,
                    cleaned_pdb,
                    cfg.get("DOCKING_MODE", "discovery"),
                    logger,
                )

            # ---- Benchmark-only figure capture to DOCKED/<PDB>/<stage>/ ----
            if str(cfg.get("DOCKING_MODE", "")).lower() == "benchmark":
                with _RENDER_LOCK:
                    stage_name = stage["name"]
                    root_project = Path(cfg.get("OVERALL_DIR", str(Path(cfg["OUTPUT_DIR"]).parent)))
                    stage_dir_target = root_project / "docked" / pdb_id / stage_name
                    stage_dir_target.mkdir(parents=True, exist_ok=True)

                    results_for_stage = score_history.get(stage_name, {})
                    best_ctrl_lig, nearest_rdk_lig = pick_control_and_nearest_rdk(
                        results_for_stage, raw_docked, control_stems_lower
                    )

                    cleaned_pdb_path = str(cleaned_pdb)
                    original_pdb_path = str(Path(paths.pdb_path))
                    ctrl_pose_path = raw_docked.get(best_ctrl_lig) if best_ctrl_lig else None
                    rdk_pose_path = raw_docked.get(nearest_rdk_lig) if nearest_rdk_lig else None

                    _render_native_on_original_pdb(
                        original_pdb=original_pdb_path,
                        outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-native__NATIVE"),
                        exclude_resns=sorted(list(EXCLUDE_HET_IDS)),
                    )

                    if ctrl_pose_path and Path(ctrl_pose_path).is_file():
                        _render_three_views_with_pymol(
                            receptor_path=cleaned_pdb_path,
                            ligand_paths_and_colors=[(ctrl_pose_path, "control", "green")],
                            outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL"),
                        )

                    if rdk_pose_path and Path(rdk_pose_path).is_file():
                        _render_three_views_with_pymol(
                            receptor_path=cleaned_pdb_path,
                            ligand_paths_and_colors=[(rdk_pose_path, "rdk_closest", "magenta")],
                            outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__RDKclosest"),
                        )

                        _render_three_views_with_pymol(
                            receptor_path=original_pdb_path,
                            ligand_paths_and_colors=[(rdk_pose_path, "rdk_closest", "magenta")],
                            outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-rdk__RDKclosest"),
                        )

                        if ctrl_pose_path and Path(ctrl_pose_path).is_file():
                            _render_three_views_with_pymol(
                                receptor_path=cleaned_pdb_path,
                                ligand_paths_and_colors=[
                                    (ctrl_pose_path, "control", "green"),
                                    (rdk_pose_path, "rdk_closest", "magenta"),
                                ],
                                outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL+RDKclosest"),
                                label_top_n_res=5,
                                label_cutoff=5.0,
                            )

            # Pocket strength row
            scores_only = [rec.get("score") for rec in results.values() if isinstance(rec.get("score"), (int, float))]
            best = min(scores_only) if scores_only else None
            n_valid = sum(1 for rec in results.values() if rec.get("valid"))
            n_tested = len(results)
            pocket_strength_rows.append([
                pocket_name,
                f"{center[0]:.3f}",
                f"{center[1]:.3f}",
                f"{center[2]:.3f}",
                (f"{best:.2f}" if isinstance(best, (int, float)) else ""),
                n_valid,
                n_tested,
            ])

    if pocket_strength_rows:
        summary_csv = out_dir / "benchmark_pocket_strength.csv"
        with open(summary_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["pocket", "center_x", "center_y", "center_z", "best_score", "n_valid", "n_tested"])
            for row in pocket_strength_rows:
                w.writerow(row)

    logger.info(f"[benchmark] Finished {pdb_id} in {time.time()-t0:.1f}s")


# =============================================================================
# Auto-analysis integration
# =============================================================================

def _run_auto_analysis(
    docked_root: Path,
    mapping_csv: Path,
    only_pdb: Optional[str],
    score_tol: float,
    center_tol: float,
    rmsd_tol: float,
    include_identity: bool,
    out_dir: Optional[Path] = None,
) -> Tuple[Optional[Path], Optional[Path]]:
    try:
        import benchmark_auto_analysis as ana
    except Exception as e:
        print(f"[analysis] ERROR: could not import benchmark_auto_analysis.py: {e}")
        return None, None

    if not hasattr(ana, "run_analysis"):
        print("[analysis] ERROR: benchmark_auto_analysis.py is missing run_analysis().")
        return None, None

    try:
        return ana.run_analysis(
            docked_root=docked_root,
            mapping_csv=mapping_csv,
            only_pdb=only_pdb,
            score_tol=score_tol,
            center_tol=center_tol,
            rmsd_tol=rmsd_tol,
            include_identity=include_identity,
            out_dir=out_dir, 
        )
    except Exception as e:
        print(f"[analysis] ERROR running run_analysis(): {e}")
        return None, None


# =============================================================================
# CLI
# =============================================================================

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Benchmark mode: single ultra-stage docking of likely co-crystal FDA ligands, per pocket."
    )
    p.add_argument("--input-dir", default=DEFAULT_INPUT_DIR, help="Directory of input .pdb files.")
    p.add_argument("--out-root", default=DEFAULT_OUT_ROOT, help="Output root for benchmark results.")
    p.add_argument("--prepped", default=DEFAULT_PREPPED_DIR, help="Directory of prepped ligand .pdbqt files.")
    p.add_argument("--mapping", default=DEFAULT_MAPPING_CSV, help="Mapping CSV (fda_mapping_from_pdbqt.csv).")
    p.add_argument("--max-candidates", type=int, default=3, help="Top-N candidate ligands per protein.")
    p.add_argument("--exhaustiveness", type=int, default=24, help="Docking exhaustiveness.")
    p.add_argument("--num-modes", type=int, default=20, help="Docking num_modes.")
    p.add_argument("--hints", help="Optional manual comma-separated hints (e.g., 'imatinib,STI571').")
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() - 2 or 4)), help="Number of proteins to process in parallel.")
    p.add_argument("--verbose", action="store_true", help="Verbose logging.")
    p.add_argument("--dry-run", action="store_true", help="List resolved paths and proteins then exit.")
    #  Harvest aliases from benchmark_analysis_details.csv
    p.add_argument(
        "--alias-from-details",
        nargs="*",
        default=[],
        help=("Path(s) or glob(s) to benchmark_analysis_details.csv to harvest aliases "
              "(e.g., E:\\PythonProject\\protein_automation\\benchmarks\\*\\_analysis\\benchmark_analysis_details.csv)"),
    )
    # --- Auto-Analysis options ---
    p.add_argument("--skip-analysis", dest="run_analysis", action="store_false", help="Skip the post-run Auto-Analysis.")
    p.set_defaults(run_analysis=True)
    p.add_argument("--analysis-only-pdb", default=None, help="Restrict analysis to one PDB (e.g., 5MO4).")
    p.add_argument("--analysis-score-tol", type=float, default=1.0, help="Score tolerance (kcal/mol).")
    p.add_argument("--analysis-center-tol", type=float, default=1.0, help="Centroid distance tol (Å).")
    p.add_argument("--analysis-rmsd-tol", type=float, default=3.0, help="RMSD tolerance (Å).")
    p.add_argument("--analysis-out", default=None, help="Optional output dir for analysis CSVs. Default: <DOCKED>\\_analysis")
    p.add_argument("--quiet", action="store_true",
               help="Quieter console: keep progress/docking notes; hide prep spam. Full logs still in files.")
    p.add_argument("--silent", action="store_true",
               help="Like --quiet but also hides non-allowlisted WARNING lines (ERRORs still show).")
    p.add_argument("--analysis-no-identity",
                   action="store_true",
                   help="Skip identity analysis (sets args.analysis_no_identity=True)")
    return p

def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_argparser().parse_args(argv)

    # Console logging setup here (protein loggers are created later)
    logging.basicConfig(
        level=(logging.DEBUG if args.verbose else logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    # ---- global console level clamp (INFO default, WARNING for --quiet, ERROR for --silent)
    root = logging.getLogger()
    if args.silent:
        root.setLevel(logging.ERROR)
        for h in root.handlers:
            h.setLevel(logging.ERROR)
    elif args.quiet:
        root.setLevel(logging.WARNING)
        for h in root.handlers:
            h.setLevel(logging.WARNING)
    else:
        root.setLevel(logging.INFO)
        for h in root.handlers:
            h.setLevel(logging.INFO)

    # ---- altLoc noise filter (hide per-atom chatter; keep summaries)
    class _AltlocNoiseFilter(logging.Filter):
        """Hide per-atom altLoc chatter while keeping summary lines."""
        def __init__(self, suppress=True):
            super().__init__()
            self.suppress = suppress
            self._drop_prefixes = ("AltLocs found for", "Keeping altLoc ")

        def filter(self, record: logging.LogRecord) -> bool:
            if not self.suppress:
                return True
            try:
                msg = record.getMessage()
            except Exception:
                msg = str(record.msg)
            return not any(msg.startswith(p) for p in self._drop_prefixes)

    suppress_altloc = (
        bool(getattr(args, "silent", False)) or
        bool(getattr(args, "quiet", False)) or
        os.environ.get("ALTLOC_QUIET_DETAILS", "1").lower() not in ("0", "false")
    )
    for _h in root.handlers:
        _h.addFilter(_AltlocNoiseFilter(suppress=suppress_altloc))

    # ---- third-party noise gates (RDKit deprecations, urllib3, etc.)
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")           # kills "DEPRECATION WARNING: please use MorganGenerator"
        os.environ.setdefault("RDKIT_LOG_LEVEL", "ERROR")
    except Exception:
        pass
    for name, level in [
        ("rdkit", logging.ERROR),
        ("urllib3", logging.ERROR),
        ("PIL", logging.ERROR),
        ("matplotlib", logging.ERROR),
        ("asyncio", logging.ERROR),
    ]:
        logging.getLogger(name).setLevel(level)

    # Resolve/normalize paths (support Windows paths on Linux)
    input_dir = _resolve_existing_path(args.input_dir)
    out_root = _resolve_existing_path(args.out_root)
    prepped = _resolve_existing_path(args.prepped)
    mapping_path = _resolve_existing_path(args.mapping)

    if not args.silent:
        print("\n[paths] Resolved CLI paths:")
        _echo_path("input_dir", input_dir)
        _echo_path("out_root", out_root)
        _echo_path("prepped", prepped)
        _echo_path("mapping", mapping_path)
        print("")

    try:
        mapping_index = MappingIndex(mapping_path)
    except Exception as e:
        logging.error(f"[fatal] Could not initialize MappingIndex: {e}")
        return

    out_root.mkdir(parents=True, exist_ok=True)
    if not input_dir.exists():
        logging.error(f"[benchmark] INPUT_DIR does not exist: {input_dir}")
        return
    if not prepped.exists():
        logging.error(f"[benchmark] PREPPED library does not exist: {prepped}")
        return

    cfg = load_inputs()
    validate_config(cfg)
    cfg = dict(cfg)  # shallow copy
    cfg["INPUT_DIR"] = str(input_dir)
    cfg["DOCKED_DIR"] = str(out_root)
    cfg["OUTPUT_DIR"] = str(out_root)
    cfg["OUTPUT_LIGANDS_DIR"] = str(prepped)
    cfg["DOCKING_MODE"] = "benchmark"
    cfg.setdefault("OVERALL_DIR", str(Path(out_root).parent))

    # feed quiet/silent CLI into cfg so the per-protein logger can read it
    cfg["_QUIET_CONSOLE"] = bool(cfg.get("_QUIET_CONSOLE", False)
                             or cfg.get("quiet_console", False)
                             or (args.quiet or args.silent))
    cfg["_SILENT_CONSOLE"] = bool(cfg.get("_SILENT_CONSOLE", False)
                              or cfg.get("silent_console", False)
                              or args.silent)

    pdb_files = [f for f in os.listdir(cfg["INPUT_DIR"]) if f.lower().endswith(".pdb") and "_nolig" not in f.lower()]
    logging.info(f"[benchmark] proteins queued: {len(pdb_files)} from {cfg['INPUT_DIR']}")
    if args.dry_run:
        if not args.silent:
            for f in pdb_files[:50]:
                print("  -", f)
            if len(pdb_files) > 50:
                print(f"  ... ({len(pdb_files)-50} more)")
            print("[dry-run] Exiting without docking.")
        return

    # Build resolver index once
    fda_index = load_library_index(mapping_path)
    manual_hints = [h.strip() for h in (args.hints or "").split(",") if h.strip()] or None

    jobs = max(1, int(args.jobs))
    base_inner = int(cfg.get("MAX_PARALLEL_JOBS", 4))
    inner_for_each = max(1, base_inner // jobs)

    logging.info(f"[benchmark] jobs={jobs} inner_workers={inner_for_each}")

    # Worker
    def _work(pdb_file: str):
        cfg_local = dict(cfg)
        cfg_local["MAX_PARALLEL_JOBS"] = inner_for_each
        run_benchmark_for_protein(
            cfg=cfg_local,
            mapping=mapping_index,
            pdb_file=pdb_file,
            prepped_dir=prepped,
            out_root=out_root,
            exhaustiveness=int(args.exhaustiveness),
            num_modes=int(args.num_modes),
            max_candidates=int(args.max_candidates),
            manual_hints=manual_hints,
            fda_index=fda_index,
        )
        return pdb_file

    # Thread pool (inner docking uses subprocesses/IO; threads avoid heavy pickling)
    t0 = time.time()
    done = failed = 0
    try:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(_work, pdb): pdb for pdb in pdb_files}
            for fut in as_completed(futures):
                pdb = futures[fut]
                try:
                    fut.result()
                    done += 1
                    logging.info(f"[benchmark]  finished {pdb} ({done}/{len(pdb_files)})")
                except Exception as e:
                    failed += 1
                    logging.error(f"[benchmark]  error on {pdb}: {e} ({done+failed}/{len(pdb_files)})")
    except KeyboardInterrupt:
        logging.error("[benchmark] Interrupted by user.")

    elapsed = time.time() - t0
    logging.info(f"[benchmark] All done: finished={done}, failed={failed}, elapsed={elapsed:.1f}s")

    # Auto-Analysis pass (optional)
    if args.run_analysis:
        docked_root = Path(cfg.get("OVERALL_DIR", str(Path(out_root).parent))) / "docked"
        if not docked_root.is_dir():
            logging.warning(f"[analysis] Docked root not found at {docked_root} — skipping.")
        else:
            details, summary = _run_auto_analysis(
                docked_root=docked_root,
                mapping_csv=mapping_path,
                only_pdb=(args.analysis_only_pdb or None),
                score_tol=float(args.analysis_score_tol),
                center_tol=float(args.analysis_center_tol),
                rmsd_tol=float(args.analysis_rmsd_tol),
                include_identity=not bool(args.analysis_no_identity),
                out_dir=(Path(args.analysis_out) if args.analysis_out else None),
            )
            if details and summary:
                logging.info(f"[analysis] Details: {details}")
                logging.info(f"[analysis] Summary: {summary}")
            else:
                logging.warning("[analysis]  Analysis did not produce outputs.")


if __name__ == "__main__":
    # Example:
    # python "benchmark mode.py" --mapping "/mnt/e/PythonProject/protein_automation/fda_mapping_from_pdbqt.csv" \
    #   --input-dir "/mnt/e/PythonProject/protein_automation/input_pdbs" \
    #   --prepped "/mnt/e/PythonProject/protein_automation/prepped_ligands" \
    #   --out-root "/mnt/e/PythonProject/protein_automation/benchmarks" \
    #   --jobs 4 --verbose
    main()
