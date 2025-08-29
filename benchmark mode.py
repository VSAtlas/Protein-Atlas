# --- snip header / your docstring stays the same ---

from __future__ import annotations

import argparse
import csv
import os
import re
import time
import glob
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from capture_pose import (
    pick_control_and_nearest_rdk,
    _render_native_on_original_pdb,
    _render_three_views_with_pymol,
    _safe_open_csv_for_write,
)
from metabolite_resolver import (
    load_library_index,
    ensure_parent_drugs_for_controls,
    resolve_corresponding_name_for_rdk,
    resolve_corresponding_name_from_text,
)
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Limit simultaneous PyMOL renders (1 by default; override with env PYMOL_PARALLEL)
_RENDER_LOCK = threading.Semaphore(int(os.environ.get("PYMOL_PARALLEL", "1")))

# ---- import pipeline pieces from your main so we reuse logic verbatim ----
from main import (
    run_one_stage, make_protein_logger, make_paths, extract_ligands_to_nolig,
    robust_prepare_controls, prepare_receptor, _count_heavy_atoms_from_pdbqt,
    CenterSelector, GlobalCenterGuard, get_recenter_params, early_recenter_decision,
    fallback_recentering_if_empty, record_score, record_le,
    final_pose_validation_and_screenshots, RetryManager, build_control_lookup,
    _fingerprint_stage, checkpoint_should_skip, checkpoint_mark_done,
    checkpoint_invalidate_from,
    detect_pocket,
)

from input_and_export_functions import load_inputs, validate_config
from protein_functions import detect_active_site

# ------------------------
# Defaults (your Windows paths)
# ------------------------
DEFAULT_INPUT_DIR   = r"E:\PythonProject\protein_automation\input_pdbs"
DEFAULT_OUT_ROOT    = r"E:\PythonProject\protein_automation\benchmarks"
DEFAULT_PREPPED_DIR = r"E:\PythonProject\protein_automation\prepped_ligands"
DEFAULT_MAPPING_CSV = r"E:\PythonProject\protein_automation\fda_mapping_from_pdbqt.csv"

# ------------------------
# Util: normalization + PDB hint parsing
# ------------------------

def _norm(s: Optional[str]) -> str:
    """
    Normalize a free-text string for fuzzy matching.
    Why: reduces punctuation/whitespace/case noise so hint matching is stable.
    """
    if s is None:
        return ""
    s = str(s).strip().lower()
    s = s.replace("\u00A0", " ")
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"[-_/\\,;:\|\[\]\(\)\{\}\.\+\*'\"]+", " ", s)
    return s

def _tokenize(s: str) -> List[str]:
    """
    Tokenize into alnum chunks (plus '+').
    Why: enables token-overlap scoring between hints and mapping name fields.
    """
    if not s:
        return []
    toks = re.split(r"[^a-z0-9\+]+", s.lower())
    return [t for t in toks if t]

def _extract_rdk_id(text: str) -> Optional[str]:
    """
    Extract an 'rdk_######' id from a filename/path if present.
    Why: direct ID matches are the strongest link to mapping rows.
    """
    if not text:
        return None
    m = re.search(r"(rdk_\d{6,8})", Path(text).stem.lower())
    return m.group(1) if m else None
def parse_pdb_het_hints(pdb_path: Path) -> Tuple[List[str], List[str]]:
    """
    Collect ligand hint signals from a PDB file.

    Returns:
      het_ids   : sorted list of unique HET/HETATM residue codes (3–5 letters),
                  filtered against EXCLUDE_HET_IDS and standard residues
      het_names : deduped human-readable names from HETNAM, lightly cleaned and
                  filtered by EXCLUDE_HET_NAME_KEYWORDS
    """
    # tolerate str input
    pdb_path = Path(pdb_path)

    het_ids: set = set()
    # het -> list of (continuation_index, text)
    hetnam_map: Dict[str, List[Tuple[int, str]]] = defaultdict(list)

    if not pdb_path.is_file():
        return [], []

    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for raw in fh:
                rec = raw[:6].strip().upper()

                if rec == "HET":
                    # HET <id> ...
                    het = (raw[7:10].strip() if len(raw) >= 10 else "").upper()
                    if not het:
                        toks = raw.split()
                        if len(toks) >= 2:
                            het = toks[1].upper()
                    if het:
                        het_ids.add(het)

                elif rec == "HETNAM":
                    # Formats:
                    #   HETNAM     XXX  Name...
                    #   HETNAM   2 XXX  Continuation...
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
                    # Residue name columns 18–20 in PDB (0-based 17:20)
                    resn = raw[17:20].strip().upper()
                    if resn:
                        het_ids.add(resn)
    except Exception:
        return [], []

    # Filter: standard residues/junk that we never want as “controls”
    STANDARD_RES = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET",
        "PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL"
    }
    het_ids = {
        h for h in het_ids
        if h not in EXCLUDE_HET_IDS and h not in STANDARD_RES and 2 <= len(h) <= 5
    }

    # Build cleaned names for the remaining HETs
    het_names: List[str] = []
    for het in sorted(het_ids):
        if het in hetnam_map:
            # stitch continuation lines in order
            parts = [t for _, t in sorted(hetnam_map[het], key=lambda x: x[0])]
            nm = " ".join(parts).strip()
            if not nm:
                continue
            nm_up = nm.upper()
            if any(kw in nm_up for kw in EXCLUDE_HET_NAME_KEYWORDS):
                continue
            cleaned = _clean_alias(nm)
            het_names.append(cleaned if cleaned else nm)

    # Deduplicate names by normalized form
    seen = set()
    dedup_names: List[str] = []
    for n in het_names:
        nn = _norm(n)
        if nn and nn not in seen:
            seen.add(nn)
            dedup_names.append(n)

    return sorted(het_ids), dedup_names

# ------------------------
# Alias augmentation from prior runs (details CSV)
# ------------------------

def _norm_text(s: Optional[str]) -> str:
    """
    Aggressive normalization (alphanum only).
    Why: used in simple heuristics that don't want punctuation differences.
    """
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())

def _looks_like_brand_or_generic(name: str) -> bool:
    """
    Heuristic: keep human-style drug names; drop long IUPAC-like strings.
    Why: we only want aliases that help fuzzy matching to real drug names.
    """
    n = (name or "").strip()
    if len(n) < 3 or len(n) > 64:
        return False
    if not re.search(r"[a-zA-Z]", n):
        return False
    noisy = len(re.findall(r"[0-9\[\]\(\)\/\.\-\,;:]", n))
    return (noisy / max(1, len(n))) < 0.40

def _clean_alias(name: str) -> Optional[str]:
    """
    Clean a proposed alias and return a lowercase name if it looks usable.
    Why: normalizes punctuation/trademarks and guards against junk entries.
    """
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

def _expand_globs(paths: Sequence[str]) -> List[Path]:
    """
    Expand any glob patterns into Paths.
    Why: allows '...\\*\\_analysis\\benchmark_analysis_details.csv' inputs.
    """
    out: List[Path] = []
    for p in paths:
        if any(ch in p for ch in "*?[]"):
            out.extend(Path(x) for x in glob.glob(p))
        else:
            out.append(Path(p))
    return out

def load_aliases_from_details_csv(paths: Sequence[str]) -> Dict[str, List[str]]:
    """
    Harvest aliases from one or more benchmark_analysis_details.csv files.
    Returns: { HET_code: [aliases...] }
    Why: learn new synonyms from previous runs to improve future matching.
    """
    exp: Dict[str, set] = {}
    for csv_path in _expand_globs(paths):
        if not csv_path.is_file():
            continue
        try:
            with open(csv_path, "r", encoding="utf-8", errors="ignore", newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    het = (row.get("control_het") or "").strip().upper()
                    nm  = _clean_alias(row.get("rdk_name") or "")
                    if not het or not nm:
                        continue
                    exp.setdefault(het, set()).add(nm)
        except Exception as e:
            print(f"[alias] WARN: failed to parse {csv_path}: {e}")
    return {k: sorted(v) for k, v in exp.items()}

def merge_aliases_into_chemcomp(base: Dict[str, List[str]], extra: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """
    Merge harvested aliases into the existing CHEMCOMP_ALIAS map (lowercased).
    Why: extend our per-HET synonym coverage without losing existing entries.
    """
    out: Dict[str, List[str]] = {k: sorted(set(vv.lower() for vv in vs)) for k, vs in base.items()}
    for het, names in (extra or {}).items():
        merged = set(out.get(het, []))
        for nm in names:
            nm_norm = (nm or "").strip().lower()
            if nm_norm:
                merged.add(nm_norm)
        out[het] = sorted(merged)
    return out

# ------------------------
# Exclusions (expanded)
# ------------------------

EXCLUDE_HET_IDS = {
    # waters
    "HOH","WAT","DOD",

    # monoatomic / simple ions (halides, alkali, alkaline earth, transition, lanthanides, etc.)
    "F","CL","BR","I",
    "LI","NA","K","RB","CS","FR",
    "MG","CA","SR","BA","RA",
    "SC","TI","V","CR","MN","FE","CO","NI","CU","ZN",
    "Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","HGA","HG","PT","AU",
    "AL","GA","IN","TL","PB","BI","SN","SB","AS","SE",
    "LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU",
    # common oxyanions / salts / crystallization additives
    "SO4","SUL","PO4","HPO","DPO","CO3","NO3","SCN","CYN","BCT","CAC","IOD","BR","NCO","CLO",
    # buffers / pH agents / precipitants / scavengers / reducing agents
    "TRS","TES","BES","HEP","HEZ","MES","MOP","PIP","ADA","CHES","CAPS","BTP","POP","MOPS","PIPES",
    "DTT","BME","TCEP","IMD","IPR","IPA","ACN","ACT","ACE","FMT","CIT","TAR","TLA","GLY","HIS","ARG","LYS",
    # cryo / solvents / polyols
    "GOL","EDO","PGO","MPD","DMS","DMF","EOH","ETOH","MEO","PEO","BU3","TBA","2MO","1PE","PE8","PEU","PG4","PG5","PG6","PGE","P33",
    # PEG fragments & glymes (very common false "controls")
    "PEG","1PG","2PG","3PG","TEG","P4G","P6G","P8G","P10","P20",
    # sugars & glycans (frequent cryo/osmolytes; usually not small-molecule drugs)
    "NAG","NDG","BMA","MAN","GLC","GAL","FUC","SIA","SLB","BGC","BOG","TRE","SUC","MAL","LMT","XYP","ARA","RIB","FRU",
    # lipids / fatty acids / sterols (often membrane additives)
    "CLR","CHL","OLA","OLE","STE","PLM","MYR","PAL","LDA","LDAO","DOD","C8E","C10E","C12","C14",
    # nucleotides / energy carriers / CoA / common cofactors
    "ATP","ADP","AMP","GTP","GDP","GMP","CTP","CDP","CMP","UTP","UDP","UMP","TTP","TDP","TMP",
    "NAD","NAP","NADH","NADP","FAD","FMN","PLP","TPP","COA","ACP","SAM","SAH","THF","FOL",
    # heme and variants / porphyrins / quinones / retinal etc.
    "HEM","HEC","HEB","HEA","UQ1","UQ2","UQ5","MEN","RET","BCR",
    # phospholipid headgroups / detergents
    "CHO","CHX","CHD","TRI","CHAP","CHAPSO","DDM","DM","UDM","MNG","LMG","NG",
}

# Optional: simple name-based screen to catch variants without fixed 3-letter IDs
EXCLUDE_HET_NAME_KEYWORDS = {
    # salts / buffers / small inorganic
    "SULFATE","PHOSPHATE","CHLORIDE","BROMIDE","IODIDE","NITRATE","CARBONATE","BICARBONATE","THIOCYANATE","CACODYLATE",
    "ACETATE","FORMATE","CITRATE","TARTRATE","IMIDAZOLE","AMMONIUM",
    # solvents / cryo
    "ETHYLENE GLYCOL","DIETHYLENE GLYCOL","TRIETHYLENE GLYCOL","GLYCEROL","MPD","ISOPROPANOL","ETHANOL","METHANOL","DMSO","DMF",
    "PEG","POLYETHYLENE GLYCOL","GLYME",
    # sugars & glycans
    "GLUCOSE","GALACTOSE","MANNOSE","FRUCTOSE","MALTOSE","TREHALOSE","SUCROSE","N-ACETYLGLUCOSAMINE",
    # detergents / lipids
    "DODECYL","MALTOSIDE","NEOPENTYL GLYCOL","DIGITONIN","TWEEN","TRITON","OCTYLGLUCOSIDE","LAURYLDIMETHYLAMINE-OXIDE","LDAO",
    "CHOLESTEROL","OLEATE","PALMITATE","STEARATE",
    # nucleotides / cofactors
    "ATP","ADP","AMP","GTP","GDP","GMP","NAD","NADP","FAD","FMN","COENZYME A","S-ADENOSYLMETHIONINE","PYRIDOXAL PHOSPHATE","THIAMINE PYROPHOSPHATE",
    "HEME","UBIQUINONE","RETINAL",
}

# --- Per-PDB manual hints and hard-coded FDA controls -----------------
# These get injected into the search hints (PER_PDB_HINTS) and, if found in your
# prepped library, are PROMOTED to the control pool (HARD_FDA_CONTROL_BY_PDB).

PER_PDB_HINTS: Dict[str, List[str]] = {
    # Helps candidate selection even if native HETs are unhelpful
    "1M17": ["erlotinib", "gefitinib", "afatinib", "osimertinib"],                          # EGFR
    "2E2B": ["zidovudine", "azt", "acyclovir", "ganciclovir"],                               # viral TK/nucleoside analogs
    "2GQG": ["sunitinib", "sorafenib", "pazopanib"],                                         # broad kinase fallbacks
    "2ITN": ["sunitinib", "pazopanib", "regorafenib"],                                       # ATP-site kinase fallback
    "2RGC": ["sotorasib", "adagrasib"],                                                      # KRAS G12C (FDA)
    "3OG7": ["dasatinib", "ponatinib", "sorafenib"],                                         # broad TKIs
    "3QX3": ["etoposide", "doxorubicin", "topotecan", "irinotecan"],                         # TOP targets
    "4RT7": ["quizartinib", "gilteritinib", "midostaurin"],                                  # FLT3
    "4XUF": ["quizartinib", "gilteritinib", "midostaurin"],                                  # FLT3
    "5C7X": ["dinoprostone", "misoprostol"],                                                 # prostaglandin pathway (best-effort)
    "6ADQ": ["atovaquone"],                                                                  # bc1 Qo-site (FDA)
    "6GQO": ["axitinib", "sunitinib", "pazopanib", "cabozantinib", "regorafenib"],           # RTK set
    "6O0L": ["venetoclax", "navitoclax"],                                                    # BCL-2 family
}

# Controls here are *promoted* to the control pool (run first, anchor pocket)
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
    "5C7X": ["dinoprostone"],   # if present in library; otherwise will log a warning and skip
    "6ADQ": ["atovaquone"],
    "6GQO": ["axitinib", "sunitinib"],
    "6O0L": ["venetoclax"],
}

# ------------------------
# ChemComp alias map (expanded)
# Keys are common PDB ligand IDs; values are synonyms (lowercased ok)
# ------------------------
CHEMCOMP_ALIAS: Dict[str, List[str]] = {
    # ABL / BCR-ABL / KIT / VEGFR TKIs
    "STI": ["imatinib","gleevec","sti571"],
    "NIL": ["nilotinib","tasigna"],
    "ABL": ["asciminib","abl001","scemblix"],
    "DAS": ["dasatinib","sprycel"],
    "BOS": ["bosutinib","bosulif"],
    "PON": ["ponatinib","iclusig"],
    "AXI": ["axitinib","inlyta"],
    "SFB": ["sorafenib","nexavar","bay 43-9006","bay439006"],
    "SU1": ["sunitinib","sutent","su11248"],
    "PAZ": ["pazopanib","votrient"],
    "CAB": ["cabozantinib","cabometyx","cometriq"],
    "REG": ["regorafenib","stivarga"],
    "VAN": ["vandetanib","zd6474","caprelsa"],

    # EGFR / ERBB2
    "ERL": ["erlotinib","tarceva"],
    "GEF": ["gefitinib","iressa"],
    "AFN": ["afatinib","gilotrif"],
    "OSM": ["osimertinib","tagrisso","azd9291"],
    "LAP": ["lapatinib","tykerb"],

    # ALK/ROS1/MET/RET
    "CRZ": ["crizotinib","xalkori"],
    "CER": ["ceritinib","zykadia"],
    "ALE": ["alectinib","alecenza"],
    "LOR": ["lorlatinib","lorbrena","lorviqua"],
    "BRG": ["brigatinib","ap26113"],
    "ENT": ["entrectinib","rxdx-101"],
    # NOTE: 'CAP' merged to include both capmatinib (oncology) and capecitabine (antimetabolite).
    # This is intentionally broad so HET=CAP in PDBs still yields helpful hints.
    "CAP": ["capmatinib","tabrecta","capecitabine","xeloda"],
    "SELr": ["selpercatinib","rxdx-105","ret inhibitor"],  # avoid clash with 'SEL' (selumetinib)
    "PRT": ["pralsetinib","blud-667","gavripranib","gprc"],

    # RAS/RAF/MEK
    "VEM": ["vemurafenib","zelboraf"],
    "DAB": ["dabrafenib","tafinlar"],
    "ENC": ["encorafenib","braftovi"],
    "COB": ["cobimetinib","cotellic"],
    "BIN": ["binimetinib","mektovi"],
    "TRM": ["trametinib","mekinist"],
    "SEL": ["selumetinib","koselugo"],

    # JAK
    "RUX": ["ruxolitinib","jakafi"],
    "TOF": ["tofacitinib","xeljanz"],
    "BAR": ["baricitinib","olumiant"],
    "UPA": ["upadacitinib","rinvoq"],
    "FED": ["fedratinib","inoma","indra","indra-280"],

    # PI3K/mTOR
    "IDA": ["idelalisib","zydelig"],
    "DUV": ["duvelisib","copiktra"],
    "COP": ["copanlisib","aliqopa"],
    "API": ["alpelisib","piqray","byl719"],
    "EVR": ["everolimus","afinitor"],
    "TMS": ["temsirolimus","torisel"],
    "RAP": ["rapamycin","sirolimus"],

    # CDK4/6
    "P31": ["palbociclib","pd-0332991","ibrance"],
    "RIB": ["ribociclib","lee011","kiskali"],
    "ABE": ["abemaciclib","ly2835219","verzenio"],

    # BCL2 / apoptosis
    "ABT": ["venetoclax","abt-199","venclexta"],
    "NAV": ["navitoclax","abt-263"],

    # BTK
    "IBR": ["ibrutinib","imbruvica"],
    "ACB": ["acalabrutinib","calquence"],
    "ZAN": ["zanubrutinib","brukinsa"],

    # FLT3 (AML)
    "QUI": ["quizartinib","ac220","vantictumab"],
    "GIL": ["gilteritinib","asp2215","xospata"],
    "CRE": ["crenolanib","cp-868596"],
    "MID": ["midostaurin","pkc412","rydapt"],
    "LST": ["lestaurtinib","cep-701"],

    # IDH (AML)
    "ENA": ["enasidenib","ag-221","idhifa"],
    "IVO": ["ivosidenib","ag-120","tibsovo"],

    # Hedgehog (AML)
    "GLB": ["glasdegib","pf-04449913","daurismo"],
    "VIS": ["vismodegib","erivedge","gdc-0449"],
    "SON": ["sonidegib","odenzo","lde225"],

    # HMAs / cytotoxics used in AML
    "AZA": ["azacitidine","vidaza"],
    "DAC": ["decitabine","dacogen"],
    "ATO": ["arsenic trioxide","trisenox"],
    "DNR": ["daunorubicin"],
    "IDR": ["idarubicin"],
    "DXR": ["doxorubicin","adriamycin"],
    "ETO": ["etoposide","vp-16"],
    "TPT": ["topotecan"],
    "IRI": ["irinotecan","cpt-11"],

    # PARP inhibitors
    "OLP": ["olaparib","lynparza"],
    "NIR": ["niraparib","zejula"],
    "RUC": ["rucaparib","rubraca"],
    "TLZ": ["talazoparib","talzenna"],
    "VLP": ["veliparib","abt-888"],

    # HDAC inhibitors
    "48D": ["vorinostat","saha","zolinza"],
    "PNB": ["panobinostat","farydak"],
    "BEL": ["belinostat","beleodaq"],
    "ROM": ["romidepsin","istodax"],

    # HRT / ER / AR axis
    "TAM": ["tamoxifen"],
    "OHT": ["4-hydroxytamoxifen","hydroxytamoxifen","endoxifen","tamoxifen"],
    "BAX": ["bazedoxifene","conbriza","duavive"],
    "FUL": ["fulvestrant","faslodex"],
    "E2" : ["estradiol","17beta-estradiol","estrogen"],
    "EST": ["estradiol","estrogen"],
    "E1" : ["estrone"],
    "DHT": ["dihydrotestosterone","androstanolone"],
    "TES": ["testosterone"],
    "PRG": ["progesterone"],
    "LET": ["letrozole","femara"],
    "ANA": ["anastrozole","arimidex"],
    "EXE": ["exemestane","aromasin"],
    "BIC": ["bicalutamide","casodex"],
    "ENZ": ["enzalutamide","xtandi"],
    "APA": ["apalutamide","erleada"],
    "DAR": ["darolutamide","nubeqa"],

    # antimetabolites
    "MTX": ["methotrexate"],
    "5FU": ["5-fluorouracil","fluorouracil"],
    # "CAP" merged above
    "GEM": ["gemcitabine","gemzar"],
    "FLUa": ["fludarabine","f-ara-a"],
    "CLD": ["cladribine","2-cda"],

    # proteasome
    "BOR": ["bortezomib","velcade"],
    "CFZ": ["carfilzomib","kyprolis"],
    "IXA": ["ixazomib","ninlaro"],

    # HIV antivirals (kept because they often show up as co-crystals and in FDA libraries)
    "RTV": ["ritonavir"],
    "LPV": ["lopinavir"],
    "ATV": ["atazanavir"],
    "EFV": ["efavirenz"],

    # others seen frequently
    "MET": ["metformin"],
    "DXN": ["dexamethasone"],
    "CPT": ["camptothecin","topotecan","irinotecan"],

    # endocrine/other
    "EVE": ["everolimus","afinitor"],  # alt key to EVR
}

# --- Benchmark-derived expansions (from benchmark_analysis_details.csv) ---
CHEMCOMP_ALIAS.update({
    "0LI": ["ponatinib", "iclusig"],                # seen in 3ZOS
    "1BQ": ["ensartinib"],                          # 4I4E
    "69Q": ["enasidenib"],                          # 5I96
    "B49": ["xenazine"],                            # 3G0E
    "BRL": ["orbenin"],                             # 3DZY
    "C6F": ["lampren"],                             # 6JQR
    "CXS": ["lampren"],                             # 6JQR
    "FLC": ["pempidine"],                           # 5TQH
    "LQQ": ["palbociclib isethiolate", "palbociclib"],  # 5L2I
    "MI1": ["xeljanz", "tofacitinib"],              # 3LXK
    "P06": ["valtrex"],                             # 4XV2
    "RXT": ["jakafi", "ruxolitinib"],               # 4U5J, 6WTN
    "VGH": ["crizotinib"],                          # 2WGJ, 2XP2, 3ZBF
})
CHEMCOMP_ALIAS.update({
    # VEGFR2
    "LEV": ["lenvatinib", "lenvima"],

    # bc1 complex (respiratory chain) controls
    "ANT": ["antimycin a"],         # cytochrome bc1 Qi-site poison
    "STG": ["stigmatellin"],        # bc1 Qo-site binder
    "ATO": ["atovaquone"],          # bc1 Qo-site antimalarial

    # RAS pocket binder (for site sanity checks vs HRAS/KRAS Switch I/II pocket)
    "F0K": ["bi-2852"],             # ligand code in 6GJ8
})
# Add a few missing alias keys (must be after CHEMCOMP_ALIAS is defined)
CHEMCOMP_ALIAS.update({
    "SOT": ["sotorasib", "lumakras"],
    "ADA": ["adagrasib", "krazati"],
    "ACV": ["acyclovir", "zovirax"],
    "ZDV": ["zidovudine", "azt"],
})

# A couple of extra “junk” HET IDs observed in the CSV that are safe to skip.
EXCLUDE_HET_IDS |= {
    "2PE",   # often phenylethanol/PEG fragment in co-crystals
}

# Name-patterns that showed up as long IUPAC-like descriptors in the CSV.
EXCLUDE_HET_NAME_KEYWORDS |= {
    "SCHEMBL",
    "CARBOXYLATE",
    "CARBOXAMIDE",
    "PROPANOATE",
    "BENZAMIDE",
    "GUANIDIN",      # catches GUANIDINE/GUANIDINIUM
    "IMINIUM",
    "PIPERIDINIUM",
    "PYRIDINIUM",
    "AZONIA",
    "OXAN",          # e.g., oxan-2-yl … (sugar-like ring fragments)
}

# Re-normalize alias values to lowercase & dedupe — single pass after all updates.
for _k, _vals in list(CHEMCOMP_ALIAS.items()):
    CHEMCOMP_ALIAS[_k] = sorted(set(v.lower() for v in _vals))

# ------------------------
# Mapping index (reads your fda_mapping_from_pdbqt.csv)
# ------------------------
@dataclass
class MappingRow:
    """
    A row from the mapping CSV describing one prepped ligand file and its names.
    Why: we score and select candidate ligands by fuzzy matching these fields.
    """
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
        """Return [(field_name, value)] for all searchable name fields."""
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
    """
    In-memory index over the mapping CSV.
    Why: efficient scoring/lookups for candidate ligand selection.
    """
    def __init__(self, csv_path: Path):
        import pandas as pd
        self.csv_path = Path(csv_path)
        if not self.csv_path.exists():
            raise FileNotFoundError(f"Mapping CSV not found: {self.csv_path}")
        self.df = pd.read_csv(self.csv_path)
        if "path" not in self.df.columns:
            raise ValueError("Mapping CSV must include a 'path' column to .pdbqt files")
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
        """
        Return all rows whose filename stem contains the given rdk_id.
        Why: strong signal when we know the library ID from file names.
        """
        rid = (rdk_id or "").strip().lower()
        if not rid:
            return []
        out = []
        for row in self.rows:
            stem = Path(row.path).stem.lower()
            if rid in stem:
                out.append(row)
        return out

    def search(self, hints: Sequence[str], inchikey: Optional[str] = None, max_results: int = 6) -> List[
        Tuple["MappingRow", int, str]]:
        """
        Score all rows against text hints (and optional InChIKey).
        Why: select the best candidates per protein by path/name overlap.
        """
        hints_norm = [_norm(h) for h in hints if h]
        hint_tokens = set(t for h in hints_norm for t in _tokenize(h))

        # Also track raw uppercase HET codes (3-5 chars) for direct filename boosts
        raw_het_codes = {h.strip().upper() for h in hints if h and 2 <= len(h.strip()) <= 5}

        out: List[Tuple[MappingRow, int, str]] = []
        for row in self.rows:
            if not row.path or not Path(row.path).exists():
                continue

            best = 0
            why = ""
            p = Path(row.path)
            stem_upper = p.stem.upper()
            base_upper = p.name.upper()
            parents_upper = " ".join([pp.name.upper() for pp in p.parents])

            # 0) hard match on InChIKey if given
            if inchikey and row.inchikey and row.inchikey.strip().upper() == inchikey.strip().upper():
                best, why = 100, "inchikey_exact"

            # 1) Boost if filename/parents clearly encode a 3–5 letter HET code from hints
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
                out.append((row, best, why))

        out.sort(key=lambda t: t[1], reverse=True)
        return out[:max_results]

# ------------------------
# Pocket detection helpers
# ------------------------
# (unchanged)

# ------------------------
# Safe CSV writer (handles Excel-lock on Windows)
# ------------------------
# (unchanged)

# ------------------------
# Helper: resolve best human name for an RDK path or a control-like name
# ------------------------
def _resolve_name_for_path_or_text(p_or_text: str, fda_index) -> str:
    """
    Resolve a displayable human name for a ligand path or free-text tag.
    Why: reports & CSVs look better with real drug names than raw ids.
    """
    rdk = _extract_rdk_id(p_or_text or "")
    if rdk:
        name = resolve_corresponding_name_for_rdk(rdk, fda_index) or ""
        if name:
            return name
    # fallback: try raw text (e.g., "OHT_A600", "NIL_A601", HETNAM text)
    stem = Path(p_or_text).stem if os.path.exists(p_or_text) else str(p_or_text)
    name = resolve_corresponding_name_from_text(stem, fda_index) or ""
    return name

def select_candidates_for_protein(
    mapping: "MappingIndex",
    hints: Sequence[str],
    prepped_dir: Path,
    extra_parent_ids: Sequence[str],
    max_candidates: int,
) -> List[Tuple["MappingRow", int, str]]:
    """
    From mapping, select top-N candidate ligands limited to this protein's prepped dir.
    Why: keeps search fast & relevant to the current target's folder.
    """
    # A) text-based search
    candidates = mapping.search(hints, inchikey=None, max_results=max(6, max_candidates * 4))

    # B) parent additions (robust to id or free-text)
    parent_rows: List[Tuple[MappingRow, int, str]] = []
    for pid in extra_parent_ids or []:
        rid = _extract_rdk_id(pid or "")
        if rid:
            for r in mapping.rows_by_rdk_id(rid):
                parent_rows.append((r, 99, "parent_from_metabolite"))
        else:
            for (r, sc, why) in mapping.search([pid], inchikey=None, max_results=2):
                parent_rows.append((r, max(sc, 92), f"{why}|parent_from_metabolite"))

    # C) keep best per absolute path and filter to this protein's prepped dir
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
    """
    Split the ligand pool into controls vs. whitelist with gating/dedupe/cap.
    Why: controls run first to anchor pocket selection and early locking.
    """
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

    MAX_CONTROLS = int(cfg.get("BENCH_MAX_CONTROLS", 8))
    if len(ctrls) > MAX_CONTROLS:
        if logger:
            logger.warning(f"[SANITY] controls={len(ctrls)} looks high; capping to {MAX_CONTROLS}. Offenders will be skipped.")
            for p in ctrls[MAX_CONTROLS:]:
                logger.warning(f"[SANITY-OFFENDER] {p}")
        ctrls = ctrls[:MAX_CONTROLS]

    return ctrls, non_ctrls

def _merge_per_pdb_hints(pdb_id: str, hints: List[str]) -> List[str]:
    """
    Blend global hints with per-PDB overrides and hard-coded controls.
    Why: per-target nudges improve candidate selection when native HETs are weak.
    """
    out = list(hints or [])
    add = []
    if pdb_id.upper() in PER_PDB_HINTS:
        add.extend(PER_PDB_HINTS[pdb_id.upper()])
    if pdb_id.upper() in HARD_FDA_CONTROL_BY_PDB:
        add.extend(HARD_FDA_CONTROL_BY_PDB[pdb_id.upper()])
    # re-dedupe by normalized tokenization
    seen = set(); merged: List[str] = []
    for h in (out + add):
        hn = _norm(h)
        if hn and hn not in seen:
            seen.add(hn); merged.append(h)
    return merged

def _promote_forced_controls_to_ctrls(
    pdb_id: str,
    mapping: "MappingIndex",
    prepped_dir: Path,
    names: Sequence[str],
    logger=None,
) -> Tuple[List[str], List[str]]:
    """
    Find prepped ligands matching 'names' and promote them to *controls*.
    Why: ensures gold-standard controls run first and can lock pocket selection.
    """
    if not names:
        return [], []
    rows = mapping.search(names, inchikey=None, max_results=max(8, len(names)*3))
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
        logger.warning(f"[forced-control] {pdb_id}: requested {list(names)} but none were found in prepped library '{prepped_dir}'.")
    return promoted_paths, promoted_stems_lower

def collect_prepped_controls_for_protein(
    paths,
    control_stems: Sequence[str],
    logger=None,
) -> Tuple[List[str], set]:
    """
    Locate this protein's control ligand .pdbqt files inside its prepped dir only.
    Why: avoids cross-target leakage while still honoring co-crystal ligands.
    """
    control_stems_lower = {s.lower() for s in control_stems if s}
    prepped_control_pdbqts: List[str] = []
    prepped_dir_resolved = paths.prepped_ligands_dir.resolve()

    if paths.prepped_ligands_dir.exists():
        for p in paths.prepped_ligands_dir.glob("*.pdbqt"):
            stem0 = p.stem.split("_stage")[0].lower()
            if stem0 in control_stems_lower and p.parent.resolve() == prepped_dir_resolved:
                prepped_control_pdbqts.append(str(p.resolve()))

    prepped_control_pdbqts = sorted(set(prepped_control_pdbqts))

    # Optional: warn on anything weird (shouldn't happen with strict parent check)
    if logger:
        leaks = [pp for pp in prepped_control_pdbqts
                 if Path(pp).resolve().parent != prepped_dir_resolved]
        if leaks:
            logger.warning(f"[CONTROL-LEAK] Found {len(leaks)} controls outside {prepped_dir_resolved}")
            for pp in leaks[:10]:
                logger.warning(f"  leak -> {pp}")

    return prepped_control_pdbqts, control_stems_lower

def _dedupe_str(seq: Sequence[str]) -> List[str]:
    """
    Stable dedupe for lists of strings (keeps first occurrence).
    Why: prevents redundant work and noisy logs while preserving order.
    """
    seen = set(); out: List[str] = []
    for s in seq:
        if s not in seen:
            seen.add(s); out.append(s)
    return out

def _path_is_within(child: Path, parent: Path) -> bool:
    """
    True if 'child' path is inside 'parent' path (Windows-safe).
    Why: avoids accidental cross-folder selection.
    """
    try:
        child = child.resolve(strict=False)
        parent = parent.resolve(strict=False)
        child.relative_to(parent)
        return True
    except Exception:
        return False

# ------------------------
# Benchmark driver — single ultra-stage per pocket
# ------------------------
def run_benchmark_for_protein(
    cfg: Dict,
    mapping: "MappingIndex",
    pdb_file: str,
    prepped_dir: Path,
    out_root: Path,
    exhaustiveness: int,
    num_modes: int,
    max_candidates: int,
    manual_hints: Optional[List[str]] = None,
    fda_index=None,  # <—  pass the metabolite/library index
) -> None:
    """
    Orchestrate a single ultra-stage docking for likely FDA ligands per pocket.
    Why: quick, targeted benchmark that anchors on controls and best candidates.
    """
    base_id = os.path.splitext(pdb_file)[0]
    pdb_id = base_id.replace("_cleaned", "").upper()

    # Logger + ASCII filter for Windows consoles
    logger = make_protein_logger(cfg["DOCKED_DIR"], pdb_id, cfg)
    import logging, sys
    def _sanitize_msg(s: str) -> str:
        return (s.replace("≤", "<=").replace("≥", ">=").replace("Å", " Angstrom")
                 .replace("µ", "u").replace("°", " deg"))
    class _AsciiFilter(logging.Filter):
        def filter(self, record):
            if isinstance(record.msg, str):
                record.msg = _sanitize_msg(record.msg)
            return True
    for h in logger.handlers:
        h.addFilter(_AsciiFilter())
    for stream_name in ("stdout", "stderr"):
        try:
            getattr(sys, stream_name).reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    paths = make_paths(cfg, base_id, pdb_file)

    # 1) Extract controls & build nolig
    _, control_stems = extract_ligands_to_nolig(paths, logger)
    robust_prepare_controls(paths, cfg, logger)

    # 1b) If the co-crystal is a **metabolite**, add the **parent drug** from your library
    extra_parent_ids: List[str] = []
    try:
        extra_parent_ids = ensure_parent_drugs_for_controls(
            pdb_code=pdb_id,
            ligands_raw_dir=str(paths.ligand_output_dir),
            fda_index=fda_index,
            max_additions=3,
        ) or []
        if extra_parent_ids:
            logger.info(f"[metabolite→parent] Will also consider parents: {', '.join(extra_parent_ids)}")
    except Exception as e:
        logger.warning(f"[metabolite→parent] resolver failed: {e}")

    # 2) Prepare receptor
    cleaned_pdb, receptor_pdbqt = prepare_receptor(cfg, paths, logger)
    if not cleaned_pdb or not receptor_pdbqt:
        logger.warning("[benchmark] receptor prep failed; skipping protein")
        return

    # 3) Detect pocket using the exact logic from main
    center, detected_box, src = detect_pocket(cleaned_pdb, paths.ligand_output_dir, logger)
    if center:
        pockets = [("pocket1", center)]
        # use main’s box (clamped like main), with a safe fallback
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

    # 4) Build matching hints (HET/HETNAM + aliases + manual + metabolite parents)
    het_ids, het_names = parse_pdb_het_hints(Path(paths.pdb_path))
    hints: List[str] = []
    hints.extend(het_names)
    hints.extend(het_ids)
    for het in het_ids:
        hints.extend(CHEMCOMP_ALIAS.get(het.upper(), []))
    if manual_hints:
        hints.extend(manual_hints)
    # add names/ids returned by metabolite->parent resolver too
    hints.extend(extra_parent_ids)

    # normalize & de-dupe by meaning
    seen = set()
    final_hints: List[str] = []
    for h in hints:
        hn = _norm(h)
        if hn and hn not in seen:
            final_hints.append(h)
            seen.add(hn)
    # Merge per-PDB overrides (manual hints + hard-coded FDA controls as search terms)
    final_hints = _merge_per_pdb_hints(pdb_id, final_hints)

    # 5) Candidate selection from mapping file (must be under prepped_dir)
    cand_rows = select_candidates_for_protein(
        mapping=mapping,
        hints=final_hints,
        prepped_dir=prepped_dir,
        extra_parent_ids=extra_parent_ids,
        max_candidates=max_candidates,
    )
    cand_rows.sort(key=lambda t: t[1], reverse=True)

    # Audit (robust to file being open)
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

    # Pre-compute heavy-atom counts (needed by CenterSelector / LE)
    heavy_atom_counts: Dict[str, int] = {}
    for p in whitelist:
        try:
            heavy_atom_counts[p] = int(_count_heavy_atoms_from_pdbqt(Path(p)))
        except Exception:
            heavy_atom_counts[p] = 0

    # Prepare control lookup (RMSD gates)
    control_lookup = build_control_lookup(paths)

    # Strictly collect THIS protein’s prepped controls and lowercase set
    prepped_control_pdbqts, control_stems_lower = collect_prepped_controls_for_protein(
        paths=paths,
        control_stems=control_stems,
        logger=logger,
    )
    # Promote hard-coded FDA controls (if present in prepped library) to CONTROL set
    if pdb_id.upper() in HARD_FDA_CONTROL_BY_PDB:
        forced_names = HARD_FDA_CONTROL_BY_PDB[pdb_id.upper()]
        forced_paths, forced_stems_lower = _promote_forced_controls_to_ctrls(
            pdb_id=pdb_id,
            mapping=mapping,
            prepped_dir=prepped_dir,
            names=forced_names,
            logger=logger,
        )
        if forced_paths:
            # Inject into the control pool
            prepped_control_pdbqts = _dedupe_str(prepped_control_pdbqts + forced_paths)
            control_stems_lower = set(control_stems_lower) | set(forced_stems_lower)

    # Ensure heavy-atom counts include controls too
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

        # Dynamics state
        params = get_recenter_params(cfg)
        guard = GlobalCenterGuard(max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2)))
        selector = CenterSelector(
            cfg,
            logger,
            control_stems=set(control_stems),
            heavy_atom_counts=heavy_atom_counts,
            initial_center=center,
        )
        box_size: Tuple[float, float, float] = pocket_box_size

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
        logger.info(f"[DEBUG] pools: ctrls={len(ctrls)}, whitelist={len(non_ctrls)}, total={len(ligands)}")

        score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
        validated_ligands_last: List[str] = []
        recenter_attempts = 0

        # Optional checkpoint skip
        if bool(cfg.get("CHECKPOINT_ENABLE", False)):
            fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
            if checkpoint_should_skip(cfg, pdb_id, stage["name"], fp):
                logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                continue

        # We allow early-recenter restarts within this single-stage loop
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
                logger.info(f"Single-stage two-wave: {len(ctrls)} controls first, then {len(non_ctrls)} whitelist.")

                # Wave A — controls
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
                        logger.info(
                            f"[CENTER] switched (controls wave): {old} -> {center} ({dec.reason}) [global switch]"
                        )
                except Exception as e:
                    logger.warning(f"CenterSelector (controls) failed: {e}")

                # Early control lock (score + proximity gates)
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
                    import numpy as np
                    if np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                        qualified_controls.append(lig)
                if len(qualified_controls) >= lock_min_hits and not guard.locked:
                    guard.lock()
                    logger.info(
                        f"[CONTROL-LOCK] n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Angstrom"
                    )

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

                # Merge
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
                    logger.info(
                        f"[CENTER] {pdb_id} {pocket_name}: {old} -> {center} ({decision.reason}) [global switch]"
                    )
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
                restarting = True
                continue

            # Adaptive shrink (optional, harmless in single-stage)
            try:
                if cfg.get("ADAPTIVE_SHRINK_ENABLE", True) and validated:
                    import numpy as np
                    med = float(np.median([d for d in distances if isinstance(d, (int, float))])) if distances else None
                    if (med is not None) and (med < float(cfg.get("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0))):
                        dec = float(cfg.get("ADAPTIVE_SHRINK_DEC", 4.0))
                        min_box = float(cfg.get("ADAPTIVE_SHRINK_MIN_BOX", 14.0))
                        new_box = tuple(max(min_box, s - dec) for s in box_size)
                        if new_box != box_size:
                            logger.info(f"Adaptive shrink: median {med:.2f} Angstrom -> box {box_size} -> {new_box}")
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
                f.write(
                    f"# pocket={pocket_name}, center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f})\n"
                )
                w = csv.writer(f)
                w.writerow(["stage", "ligand", "library_id", "corresponding_name", "score", "valid", "reason"])
                for lig, rec in ordered:
                    sc = rec.get("score")
                    sc_str = f"{sc:.2f}" if isinstance(sc, (int, float)) else ""
                    lib_id = _extract_rdk_id(lig) or ""
                    corr = _resolve_name_for_path_or_text(lig, fda_index)
                    w.writerow(
                        [
                            stage["name"],
                            os.path.basename(lig),
                            lib_id,
                            corr,
                            sc_str,
                            int(bool(rec.get("valid", False))),
                            str(rec.get("reason", "")),
                        ]
                    )

            # Final validation/screenshots (unchanged call)
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

                    # (A) Original PDB with native ligand
                    _render_native_on_original_pdb(
                        original_pdb=original_pdb_path,
                        outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-native__NATIVE"),
                        exclude_resns=sorted(list(EXCLUDE_HET_IDS)),
                    )

                    # (B) Cleaned receptor + CONTROL
                    if ctrl_pose_path and Path(ctrl_pose_path).is_file():
                        _render_three_views_with_pymol(
                            receptor_path=cleaned_pdb_path,
                            ligand_paths_and_colors=[(ctrl_pose_path, "control", "green")],
                            outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL"),
                        )

                    # (C) Cleaned receptor + RDK-closest
                    if rdk_pose_path and Path(rdk_pose_path).is_file():
                        _render_three_views_with_pymol(
                            receptor_path=cleaned_pdb_path,
                            ligand_paths_and_colors=[(rdk_pose_path, "rdk_closest", "magenta")],
                            outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__RDKclosest"),
                        )

                        # (D) Original PDB + RDK-closest
                        _render_three_views_with_pymol(
                            receptor_path=original_pdb_path,
                            ligand_paths_and_colors=[(rdk_pose_path, "rdk_closest", "magenta")],
                            outprefix=str(stage_dir_target / f"{pdb_id}_orig-with-rdk__RDKclosest"),
                        )

                        # (E) Cleaned receptor + CONTROL + RDK-closest (overlay)
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
                        # (F) Interactive PyMOL session (.pml scenes) + optional launch
                        try:
                            from capture_pose import write_multiview_pml, launch_pymol_with_pml
                            pml_path = stage_dir_target / f"{pdb_id}_{stage_name}__CONTROL+RDKclosest.pml"
                            write_multiview_pml(
                                receptor_path=cleaned_pdb_path,
                                control_path=(ctrl_pose_path or ""),
                                rdk_path=(rdk_pose_path or ""),
                                out_pml=pml_path,
                                label_top_n_res=int(cfg.get("LABEL_TOP_N_RES", 5)),
                                label_cutoff=float(cfg.get("LABEL_CUTOFF_ANG", 5.0)),
                            )
                            if bool(cfg.get("OPEN_PYMOL_INTERACTIVE", False)):
                                launch_pymol_with_pml(pml_path, cfg.get("PYMOL_EXE"))
                            else:
                                logger.info(
                                    f"[PyMOL] Interactive .pml written: {pml_path.name} (double-click to open in PyMOL)")
                        except Exception as e:
                            logger.warning(f"[PyMOL] Failed to create/open interactive .pml: {e}")
                        # (G) CONTROL vs RDK-only (no receptor shown)
                        try:
                            if ctrl_pose_path and rdk_pose_path and Path(ctrl_pose_path).is_file() and Path(
                                    rdk_pose_path).is_file():
                                rdk_id_short = _extract_rdk_id(Path(rdk_pose_path).stem) or "RDKclosest"
                                _render_three_views_with_pymol(
                                    receptor_path=cleaned_pdb_path,  # load but hide it (see patch B)
                                    ligand_paths_and_colors=[
                                        (ctrl_pose_path, "control", "green"),
                                        (rdk_pose_path, "rdk_closest", "magenta"),
                                    ],
                                    outprefix=str(stage_dir_target / f"{pdb_id}_cleaned__CONTROL+{rdk_id_short}_PAIR"),
                                    hide_receptor=True,  # <— new kwarg added in patch B
                                    label_top_n_res=0,  # off for ligand-only shots
                                )
                        except Exception as e:
                            logger.warning(f"[PyMOL] pair-only render failed: {e}")

            # Pocket strength row
            scores_only = [rec.get("score") for rec in results.values() if isinstance(rec.get("score"), (int, float))]
            best = min(scores_only) if scores_only else None
            n_valid = sum(1 for rec in results.values() if rec.get("valid"))
            n_tested = len(results)
            pocket_strength_rows.append(
                [
                    pocket_name,
                    f"{center[0]:.3f}",
                    f"{center[1]:.3f}",
                    f"{center[2]:.3f}",
                    (f"{best:.2f}" if isinstance(best, (int, float)) else ""),
                    n_valid,
                    n_tested,
                ]
            )

    # Write per-protein pocket summary
    if pocket_strength_rows:
        summary_csv = out_dir / "benchmark_pocket_strength.csv"
        with open(summary_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["pocket", "center_x", "center_y", "center_z", "best_score", "n_valid", "n_tested"])
            for row in pocket_strength_rows:
                w.writerow(row)

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
    """
    Import benchmark_auto_analysis.py and run its summary scorer.
    Why: produce cross-target CSVs of match quality after docking finishes.
    """
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

# ------------------------
# CLI
# ------------------------
def build_argparser() -> argparse.ArgumentParser:
    """
    Define CLI for benchmark driver and optional auto-analysis.
    Why: make the workflow reproducible and batch-friendly.
    """
    p = argparse.ArgumentParser(description="Benchmark mode: single ultra-stage docking of likely co-crystal FDA ligands, per pocket.")
    p.add_argument("--input-dir", default=DEFAULT_INPUT_DIR)
    p.add_argument("--out-root",  default=DEFAULT_OUT_ROOT)
    p.add_argument("--prepped",   default=DEFAULT_PREPPED_DIR)
    p.add_argument("--mapping",   default=DEFAULT_MAPPING_CSV)
    p.add_argument("--max-candidates", type=int, default=3)
    p.add_argument("--exhaustiveness", type=int, default=24)
    p.add_argument("--num-modes", type=int, default=20)
    p.add_argument("--hints", help="Optional manual comma-separated hints (e.g., 'imatinib,STI571')")
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4)//2),
                   help="Number of proteins to process in parallel")
    # NEW: harvest aliases from benchmark_analysis_details.csv
    p.add_argument(
        "--alias-from-details",
        nargs="*",
        default=[],
        help="Path(s) or glob(s) to benchmark_analysis_details.csv to harvest aliases (e.g., E:\\PythonProject\\protein_automation\\benchmarks\\*\\_analysis\\benchmark_analysis_details.csv)"
    )
    # --- Auto-Analysis options ---
    p.add_argument("--skip-analysis", dest="run_analysis", action="store_false",
                   help="Skip the post-run Auto-Analysis.")
    p.set_defaults(run_analysis=True)
    p.add_argument("--analysis-only-pdb", default=None,
                   help="Optional: restrict analysis to one PDB (e.g., 5MO4).")
    p.add_argument("--analysis-score-tol", type=float, default=1.0,
                   help="Score tolerance (kcal/mol) for a point (default: 1.0).")
    p.add_argument("--analysis-center-tol", type=float, default=1.0,
                   help="Centroid distance tolerance (Å) for a point (default: 1.0).")
    p.add_argument("--analysis-rmsd-tol", type=float, default=3.0,
                   help="RMSD tolerance (Å) for a point (default: 3.0).")
    p.add_argument("--analysis-no-identity", action="store_true",
                   help="Exclude identity-match from scoring (max 3 points instead of 4).")
    p.add_argument("--analysis-out", default=None,
                   help="Optional output dir for analysis CSVs. Default: <DOCKED>\\_analysis")

    return p

def main(argv: Optional[Sequence[str]] = None) -> None:
    """
    Entry point: parse args, build indices/config, run benchmark(s), run analysis.
    Why: glue everything together for a single command execution.
    """
    args = build_argparser().parse_args(argv)

    input_dir = Path(args.input_dir)
    out_root  = Path(args.out_root)
    prepped   = Path(args.prepped)
    mapping   = MappingIndex(Path(args.mapping))

    # Harvest and merge aliases from prior analysis CSV(s) — expands CHEMCOMP_ALIAS
    alias_sources = args.alias_from_details or []
    if alias_sources:
        try:
            extra_aliases = load_aliases_from_details_csv(alias_sources)
            if extra_aliases:
                global CHEMCOMP_ALIAS
                CHEMCOMP_ALIAS = merge_aliases_into_chemcomp(CHEMCOMP_ALIAS, extra_aliases)
                # normalize again to be safe
                for _k, _vals in list(CHEMCOMP_ALIAS.items()):
                    CHEMCOMP_ALIAS[_k] = sorted(set(v.lower() for v in _vals))
                total_added = sum(len(v) for v in extra_aliases.values())
                print(f"[alias] Expanded CHEMCOMP_ALIAS with {total_added} aliases across {len(extra_aliases)} HET codes.")
        except Exception as e:
            print(f"[alias] WARN: could not augment aliases: {e}")

    out_root.mkdir(parents=True, exist_ok=True)
    if not input_dir.exists():
        print(f"[benchmark] INPUT_DIR does not exist: {input_dir}")
        return
    if not prepped.exists():
        print(f"[benchmark] PREPPED library does not exist: {prepped}")
        return

    cfg = load_inputs(); validate_config(cfg)
    cfg = dict(cfg)
    cfg["INPUT_DIR"] = str(input_dir)
    cfg["DOCKED_DIR"] = str(out_root)
    cfg["OUTPUT_DIR"] = str(out_root)
    cfg["OUTPUT_LIGANDS_DIR"] = str(prepped)
    cfg["DOCKING_MODE"] = "benchmark"
    cfg.setdefault("OVERALL_DIR", str(Path(out_root).parent))

    pdb_files = [f for f in os.listdir(cfg["INPUT_DIR"]) if f.lower().endswith(".pdb") and "_nolig" not in f.lower()]
    print(f"[benchmark] proteins queued: {len(pdb_files)} from {cfg['INPUT_DIR']}")

    # Build resolver index once
    fda_index = load_library_index(Path(args.mapping))
    manual_hints = [h.strip() for h in (args.hints or "").split(",") if h.strip()] or None

    jobs = max(1, int(args.jobs))
    # Avoid oversubscription: shrink per-protein ligand workers
    base_inner = int(cfg.get("MAX_PARALLEL_JOBS", 4))
    inner_for_each = max(1, base_inner // jobs)

    # submit in a pool
    def _work(pdb_file: str):
        """
        Worker wrapper to scale inner parallelism and run a single protein.
        Why: keeps overall CPU utilization balanced across targets.
        """
        cfg_local = dict(cfg)
        cfg_local["MAX_PARALLEL_JOBS"] = inner_for_each
        run_benchmark_for_protein(
            cfg=cfg_local,
            mapping=mapping,
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

    # Use threads: inner docking already uses subprocesses/IO; avoids heavy pickling of pandas mapping
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(_work, pdb): pdb for pdb in pdb_files}
        done, failed = 0, 0
        for fut in as_completed(futures):
            pdb = futures[fut]
            try:
                fut.result()
                done += 1
                print(f"[benchmark] ✔ finished {pdb} ({done}/{len(pdb_files)})")
            except Exception as e:
                failed += 1
                print(f"[benchmark] ✖ error on {pdb}: {e} ({done+failed}/{len(pdb_files)})")

    # --- Auto-Analysis pass (optional) ---
    if args.run_analysis:
        # Your project writes images/poses under <OVERALL_DIR>/docked
        docked_root = Path(cfg.get("OVERALL_DIR", str(Path(out_root).parent))) / "docked"
        if not docked_root.is_dir():
            print(f"[analysis] Docked root not found at {docked_root} — skipping.")
        else:
            details, summary = _run_auto_analysis(
                docked_root=docked_root,
                mapping_csv=Path(args.mapping),
                only_pdb=(args.analysis_only_pdb or None),
                score_tol=float(args.analysis_score_tol),
                center_tol=float(args.analysis_center_tol),
                rmsd_tol=float(args.analysis_rmsd_tol),
                include_identity=not bool(args.analysis_no_identity),
                out_dir=(Path(args.analysis_out) if args.analysis_out else None),
            )
            if details and summary:
                print(f"[analysis] ✅ Details: {details}")
                print(f"[analysis] ✅ Summary: {summary}")
            else:
                print("[analysis] ❌ Analysis did not produce outputs.")

if __name__ == "__main__":
    # DEV: uncomment while iterating in IDE
    # main([
    #   "--input-dir", r"E:\PythonProject\protein_automation\input_pdbs",
    #   "--out-root",  r"E:\PythonProject\protein_automation\benchmarks",
    #   "--prepped",   r"E:\PythonProject\protein_automation\prepped_ligands",
    #   "--mapping",   r"E:\PythonProject\protein_automation\fda_mapping_from_pdbqt.csv",
    #   "--alias-from-details", r"E:\PythonProject\protein_automation\benchmarks\docked\_analysis\benchmark_analysis_details.csv",
    #   "--jobs", "4",
    # ])
    main()
