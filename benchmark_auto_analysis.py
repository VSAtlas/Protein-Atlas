# --- snip header / your docstring stays the same ---
"""
Benchmark Auto-Analysis — IDE-friendly, self-contained script (v0.3)

What this does
- Crawls: DOCKED/<PDB>/bench_pocketX_single/ (e.g., bench_pocket1_single, bench_pocket2_single)
- Finds controls like: NIL_A601_bench_pocket2_single.best.pdb
- Finds RDK library poses like: rdk_0003792_bench_pocket2_single.pdbqt
- Parses Vina scores from REMARK lines, computes pose centroids, and an approximate RMSD.
- Pairs each CONTROL to the best RDK match (prefer identity when mapping name is available, else centroid-nearest).
- Scores each pair up to **4 points**:
    1) |Δ score| ≤ 1.0 kcal/mol
    2) Centroid distance ≤ 1.0 Å
    3) Pose–pose RMSD ≤ 3.0 Å
    4) Identity match (RDKit ligand resolves to the *same drug* as the control via mapping or alias)
- Writes two CSVs under DOCKED/_analysis/ by default:
    • benchmark_analysis_details.csv  (one row per control pairing; now includes pass_category)
    • benchmark_analysis_summary.csv  (one row per PDB; now includes passed/near/fail flags)

Notes
- This file is **standalone**: no project-internal imports are required. Only dependency is NumPy.
- Identity match (#4) improves when you provide your FDA mapping CSV. If omitted, #4 usually = 0 (fallback to nearest).
- RMSD uses a robust, order-agnostic approximation (nearest-neighbor pairing + Kabsch alignment) for heavy atoms only.
- “Pass / Near / Fail” in summary is computed from **best_total_points** per PDB:
    • Pass  = best_total_points == max_points (4 with identity, 3 without)
    • Near  = 1 < best_total_points < max_points
    • Fail  = best_total_points == 0
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import base64
import html as _html
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# -----------------------------
# Defaults (edit for IDE usage)
# -----------------------------
DEFAULT_DOCKED_ROOT = r"E:\PythonProject\protein_automation\docked"
DEFAULT_MAPPING_CSV = r"E:\PythonProject\protein_automation\fda_mapping_from_pdbqt.csv"  # optional
DEFAULT_ONLY_PDB = None  # e.g., "5MO4"
DEFAULT_SCORE_TOL = 1.0   # kcal/mol
DEFAULT_CENTER_TOL = 1.0  # Å
DEFAULT_RMSD_TOL = 3.0    # Å
DEFAULT_OUTDIR = None     # None -> <DOCKED>\_analysis
DEFAULT_INCLUDE_IDENTITY = True  # toggle 4th criterion
# -----------------------------
# Filename patterns
# -----------------------------
SCORE_PAT = re.compile(r"REMARK\s+VINA\s+RESULT[:\s]+(-?\d+\.\d+)")
CONTROL_PAT = re.compile(r"^(?P<het>[A-Za-z0-9]{3})_\w\d+_bench_pocket\d+_single\.best\.(?:pdb|pdbqt)$", re.I)
RDK_PAT = re.compile(r"^(rdk_\d{6,8})_bench_pocket\d+_single\.(?:pdbqt|pdb)$", re.I)

# 3-letter HET aliases → lowercase synonyms (extend as needed)
# Notes:
# - Keys are common 3–5 letter PDB ligand IDs you’re likely to see as co-crystal “controls”.
# - Values are lowercase synonyms (generic name, brand names, dev codes) to help identity matching.
# - This isn’t exhaustive, but it’s intentionally wide to improve hit rates across oncology kinases + frequent controls.
CHEMCOMP_ALIAS: Dict[str, List[str]] = {
    # ---------------------
    # ABL / BCR-ABL / KIT / VEGFR TKIs
    # ---------------------
    "STI": ["imatinib", "gleevec", "glivec", "sti571"],
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
    "AVA": ["avapritinib", "bluestone", "ayvakit", "bluestone-285", "bnd-285", "bnd285"],
    "RIP": ["ripretinib", "qinlock", "dcc-2618", "dcc2618"],
    "TIV": ["tivozanib", "fotivda"],
    "LENv": ["lenvatinib", "lenvima", "e7080"],  # use 'LENv' to avoid clash with lenalidomide

    # ---------------------
    # EGFR / ERBB2
    # ---------------------
    "ERL": ["erlotinib", "tarceva"],
    "GEF": ["gefitinib", "iressa"],
    "AFN": ["afatinib", "gilotrif"],
    "OSM": ["osimertinib", "tagrisso", "azd9291"],
    "LAP": ["lapatinib", "tykerb", "tyverb"],
    "NER": ["neratinib", "nerlynx", "hkI-272", "hki272"],
    "PYR": ["pyrotinib", "iqi", "ih-901"],

    # ---------------------
    # ALK / ROS1 / MET / RET / NTRK
    # ---------------------
    "CRZ": ["crizotinib", "xalkori"],
    "CER": ["ceritinib", "zykadia"],
    "ALE": ["alectinib", "alecenza"],
    "LOR": ["lorlatinib", "lorbrena", "lorviqua"],
    "BRG": ["brigatinib", "ap26113"],
    "ENT": ["entrectinib", "rxdx-101", "rxdx101"],
    "CAP": ["capmatinib", "tabrecta"],
    "TEP": ["tepotinib", "tepmeko", "msc2156119j"],  # 'tepmeko' (JP), add dev code
    "PRT": ["pralsetinib", "blud-667", "gavripranib", "gprc", "blud667"],
    "SELr": ["selpercatinib", "rxdx-105", "loxo-292", "loxO292", "ret inhibitor", "rxdx105", "loxo292"],  # avoid clash w/ SEL (selumetinib)
    "LAR": ["larotrectinib", "vitrakvi", "loxo-101", "loxo101"],
    "SLT": ["selitrectinib", "loxo-195", "loxo195"],

    # ---------------------
    # RAS / RAF / MEK / ERK
    # ---------------------
    "VEM": ["vemurafenib", "zelboraf", "plx4032", "ro5185426"],
    "DAB": ["dabrafenib", "tafinlar", "gsk2118436"],
    "ENC": ["encorafenib", "braftovi", "lgx818"],
    "COB": ["cobimetinib", "cotellic", "gdc-0973", "gdc0973"],
    "BIN": ["binimetinib", "mektovi", "meK162", "mek162"],
    "TRM": ["trametinib", "mekinist", "gsk1120212"],
    "SEL": ["selumetinib", "koselugo", "azd6244"],
    "ULX": ["ulixertinib", "bvd-523", "bvd523"],
    "LY4": ["ly3214996", "ly-3214996"],
    "SOT": ["sotorasib", "amg510", "lumakras", "lumykras", "amg-510"],
    "ADA": ["adagrasib", "mrtx849", "krazati", "mrtx-849"],

    # ---------------------
    # JAK
    # ---------------------
    "RUX": ["ruxolitinib", "jakafi", "jakavi"],
    "TOF": ["tofacitinib", "xeljanz"],
    "BAR": ["baricitinib", "olumiant"],
    "UPA": ["upadacitinib", "rinvoq"],
    "FED": ["fedratinib", "indra", "indra-280", "indra280", "inoma"],

    # ---------------------
    # PI3K / AKT / mTOR
    # ---------------------
    "IDA": ["idelalisib", "zydelig"],
    "DUV": ["duvelisib", "copiktra"],
    "COP": ["copanlisib", "aliqopa"],
    "API": ["alpelisib", "piqray", "byl719"],
    "IPI": ["ipatasertib", "gdc-0068", "gdc0068"],
    "CAPV": ["capivasertib", "truqap", "azu-010", "azu010", "azd5363", "azd-5363"],
    "EVR": ["everolimus", "afinitor", "certican"],
    "TMS": ["temsirolimus", "torisel"],
    "RAP": ["rapamycin", "sirolimus"],

    # ---------------------
    # CDK
    # ---------------------
    "P31": ["palbociclib", "pd-0332991", "ibrance", "pd0332991"],
    "RIB": ["ribociclib", "lee011", "kiskali", "lee-011"],
    "ABE": ["abemaciclib", "ly2835219", "verzenio", "ly-2835219"],
    "ALV": ["alvocidib", "flavopiridol", "hmds-698", "hmds698"],

    # ---------------------
    # BCL2 / apoptosis
    # ---------------------
    "ABT": ["venetoclax", "abt-199", "venclexta", "abt199"],
    "NAV": ["navitoclax", "abt-263", "abt263"],
    "NUT": ["nutlin-3", "rg7112", "mdm2 inhibitor", "nutlin3"],

    # ---------------------
    # BTK
    # ---------------------
    "IBR": ["ibrutinib", "imbruvica", "pci-32765", "pci32765"],
    "ACB": ["acalabrutinib", "calquence", "acp-196", "acp196"],
    "ZAN": ["zanubrutinib", "brukinsa", "bgb-3111", "bgb3111"],
    "PYR3": ["pirtobrutinib", "loxO-305", "loxo305"],

    # ---------------------
    # FLT3 (AML)
    # ---------------------
    "QUI": ["quizartinib", "ac220", "vantictumab"],  # keep quizartinib primary
    "GIL": ["gilteritinib", "asp2215", "xospata"],
    "CRE": ["crenolanib", "cp-868596", "cp868596"],
    "MID": ["midostaurin", "pkc412", "rydapt"],
    "LST": ["lestaurtinib", "cep-701", "cep701"],

    # ---------------------
    # IDH (AML)
    # ---------------------
    "ENA": ["enasidenib", "ag-221", "idhifa", "ag221"],
    "IVO": ["ivosidenib", "ag-120", "tibsovo", "ag120"],

    # ---------------------
    # Hedgehog (AML)
    # ---------------------
    "GLB": ["glasdegib", "pf-04449913", "daurismo", "pf04449913"],
    "VIS": ["vismodegib", "erivedge", "gdc-0449", "gdc0449"],
    "SON": ["sonidegib", "odenzo", "lde225"],

    # ---------------------
    # FGFR
    # ---------------------
    "ERD": ["erdafitinib", "balversa", "jnj-42756493", "jnj42756493"],
    "PMG": ["pemigatinib", "pemazyre", "in-109", "in109", "incb-054828", "incb054828"],
    "INF": ["infigratinib", "truseltiq", "bgj398"],
    "FUT": ["futibatinib", "lytgobi", "tas-120", "tas120"],

    # ---------------------
    # PARP inhibitors
    # ---------------------
    "OLP": ["olaparib", "lynparza"],
    "NIR": ["niraparib", "zejula"],
    "RUC": ["rucaparib", "rubraca"],
    "TLZ": ["talazoparib", "talzenna"],
    "VLP": ["veliparib", "abt-888", "abt888"],

    # ---------------------
    # HDAC inhibitors
    # ---------------------
    "48D": ["vorinostat", "saha", "zolinza"],
    "PNB": ["panobinostat", "farydak"],
    "BEL": ["belinostat", "beleodaq"],
    "ROM": ["romidepsin", "istodax"],
    "TSA": ["trichostatin a", "ts a", "trichostatin"],

    # ---------------------
    # HRT / ER / AR axis (SERMs, SERDs, AIs, AR antagonists)
    # ---------------------
    "TAM": ["tamoxifen"],
    "OHT": ["4-hydroxytamoxifen", "hydroxytamoxifen", "endoxifen", "tamoxifen"],
    "RAX": ["raloxifene", "evista"],
    "BAX": ["bazedoxifene", "conbriza", "duavive"],
    "FUL": ["fulvestrant", "faslodex"],
    "E2":  ["estradiol", "17beta-estradiol", "estrogen"],
    "EST": ["estradiol", "estrogen"],
    "E1":  ["estrone"],
    "DHT": ["dihydrotestosterone", "androstanolone"],
    "TES": ["testosterone"],
    "PRG": ["progesterone"],
    "LET": ["letrozole", "femara"],
    "ANA": ["anastrozole", "arimidex"],
    "EXE": ["exemestane", "aromasin"],
    "BIC": ["bicalutamide", "casodex"],
    "ENZ": ["enzalutamide", "xtandi", "mdv3100", "mdv-3100"],
    "APA": ["apalutamide", "erleada", "arn-509", "arn509"],
    "DAR": ["darolutamide", "nubeqa", "baye-2235", "baye2235"],

    # ---------------------
    # Antimetabolites / cytotoxics
    # ---------------------
    "MTX": ["methotrexate"],
    "5FU": ["5-fluorouracil", "fluorouracil"],
    "CAPe": ["capecitabine", "xeloda"],  # 'CAP' used by capmatinib; use 'CAPe' for capecitabine
    "GEM": ["gemcitabine", "gemzar"],
    "FLUa": ["fludarabine", "f-ara-a", "farada"],
    "CLD": ["cladribine", "2-cda"],
    "CPT": ["camptothecin"],
    "IRI": ["irinotecan", "cpt-11", "cpt11", "camptosar"],
    "TPT": ["topotecan", "hycamtin"],
    "ETO": ["etoposide", "vp-16", "vp16"],
    "DXR": ["doxorubicin", "adriamycin"],
    "DNR": ["daunorubicin"],
    "IDR": ["idarubicin"],

    # Microtubules (taxanes & vincas)
    "PTX": ["paclitaxel", "taxol"],
    "TAX": ["paclitaxel", "taxol"],
    "DTX": ["docetaxel", "taxotere"],
    "DOC": ["docetaxel", "taxotere"],
    "VCR": ["vincristine", "oncovin"],
    "VBL": ["vinblastine", "velban"],
    "VRB": ["vinorelbine", "navelbine"],

    # IMiDs
    "LEN": ["lenalidomide", "revlimid"],
    "POM": ["pomalidomide", "pomalyst"],
    "THD": ["thalidomide", "thalomid"],

    # ---------------------
    # Proteasome
    # ---------------------
    "BOR": ["bortezomib", "velcade", "ps-341", "ps341"],
    "CFZ": ["carfilzomib", "kyprolis", "pr-171", "pr171"],
    "IXA": ["ixazomib", "ninlaro", "mln9708"],

    # ---------------------
    # DNA repair / other targeted
    # ---------------------
    "ATRi": ["berzosertib", "m6620", "vx-970", "vx970", "berzosertib"],
    "ATM": ["azd0156", "azd-0156"],
    "CHK": ["prexasertib", "ly2606368", "ly-2606368"],

    # ---------------------
    # HIV antivirals (seen as crystallization controls)
    # ---------------------
    "RTV": ["ritonavir"],
    "LPV": ["lopinavir"],
    "ATV": ["atazanavir"],
    "EFV": ["efavirenz"],

    # ---------------------
    # HSP90 & bromodomains
    # ---------------------
    "GAN": ["ganetespib", "sta-9090", "sta9090"],
    "LUM": ["luminespib", "nvp-auy922", "auy922", "nvpauy922"],
    "JQ1": ["jq1", "bromodomain inhibitor", "brd4 inhibitor"],
    "AZ5": ["azd5153", "azd-5153"],

    # ---------------------
    # Common endocrine/other small molecules seen often
    # ---------------------
    "MET": ["metformin"],
    "DXN": ["dexamethasone"],
    "EVE": ["everolimus", "afinitor"],  # alt key to EVR
}
# normalize alias values to lowercase & unique
for k in list(CHEMCOMP_ALIAS.keys()):
    CHEMCOMP_ALIAS[k] = sorted(set(v.lower() for v in CHEMCOMP_ALIAS[k]))

# name normalization for robust identity checks
from typing import Optional as _Optional

def _norm_text(s: _Optional[str]) -> str:
    return re.sub(r'[^a-z0-9]+', '', (s or '').lower())

# Debug flag (env or CLI --debug)
DEBUG = bool(int(os.environ.get("BENCH_DEBUG", "0")))


def alignment_metrics(ctrl: np.ndarray, rdk: np.ndarray, nn_cap: int = 128):
    """
    Returns (rmsd, raw_centroid, aligned_centroid, n_ctrl, n_rdk).
    Aligns RDK -> CTRL using nearest-neighbor pairing + Kabsch.
    """
    n_ctrl = int(ctrl.shape[0])
    n_rdk  = int(rdk.shape[0])
    raw_centroid = None
    aligned_centroid = None
    rmsd = None

    if n_ctrl >= 1 and n_rdk >= 1:
        c_ctr = ctrl.mean(axis=0)
        r_ctr = rdk.mean(axis=0)
        raw_centroid = float(np.linalg.norm(c_ctr - r_ctr))

    if n_ctrl < 3 or n_rdk < 3:
        return rmsd, raw_centroid, aligned_centroid, n_ctrl, n_rdk

    # build NN pairs (A=ctrl reference, B=rdk moved to A)
    A = ctrl
    B = rdk
    if A.shape[0] > nn_cap:
        idx = np.linspace(0, A.shape[0] - 1, nn_cap, dtype=int)
        A = A[idx]

    used = set()
    pairs_a, pairs_b = [], []
    for a in A:
        best_j, best_d2 = None, 1e99
        for j, b in enumerate(B):
            if j in used:
                continue
            d2 = ((a - b) ** 2).sum()
            if d2 < best_d2:
                best_d2, best_j = d2, j
        if best_j is not None:
            used.add(best_j)
            pairs_a.append(a)
            pairs_b.append(B[best_j])

    if len(pairs_a) < 3:
        return rmsd, raw_centroid, aligned_centroid, n_ctrl, n_rdk

    P = np.vstack(pairs_a)  # ctrl
    Q = np.vstack(pairs_b)  # rdk

    # R that maps Q -> P (note the reversed order):
    R, _ = _kabsch(Q, P)
    Qc = Q - Q.mean(axis=0)
    Pc = P - P.mean(axis=0)
    # RMSD after alignment
    diff2 = ((Qc @ R) - Pc) ** 2
    rmsd = float(np.sqrt(diff2.sum(axis=1).mean()))

    # Apply transform to ALL rdk atoms for centroid-after-alignment
    rdk_aligned = (rdk - Q.mean(axis=0)) @ R + P.mean(axis=0)
    aligned_centroid = float(np.linalg.norm(ctrl.mean(axis=0) - rdk_aligned.mean(axis=0)))
    return rmsd, raw_centroid, aligned_centroid, n_ctrl, n_rdk


# -----------------------------
# Mapping CSV (optional, no pandas needed)
# -----------------------------
@dataclass
class MapRow:
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

    def best_name(self) -> str:
        for f in (
            self.display_name,
            self.generic_name,
            self.remark_name,
            self.sdf_title,
            self.brand_names,
            self.pubchem_name,
            self.pubchem_record_title,
        ):
            if f and str(f).strip():
                return str(f).strip()
        return ""


class MappingIndex:
    def __init__(self, csv_path: Optional[Path]):
        self.rows: List[MapRow] = []
        if not csv_path or not Path(csv_path).is_file():
            return
        try:
            with open(csv_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
                reader = csv.DictReader(f)
                for r in reader:
                    self.rows.append(
                        MapRow(
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
                        )
                    )
        except Exception as e:
            print(f"[analysis] WARNING: failed to parse mapping CSV: {e}")

    def resolve_name_for_rdk(self, rdk_id: str) -> Optional[str]:
        rid = (rdk_id or "").strip().lower()
        if not rid:
            return None
        best: Optional[str] = None
        for row in self.rows:
            p = str(row.path or "").lower()
            if rid in p:
                n = row.best_name()
                if n:
                    best = n
                    break
        return best


# -----------------------------
# PDB/PDBQT parsing helpers
# -----------------------------

def parse_vina_score(path: Path) -> Optional[float]:
    try:
        with open(path, "r", errors="ignore") as f:
            for line in f:
                m = SCORE_PAT.search(line)
                if m:
                    return float(m.group(1))
    except Exception:
        pass
    return None


def parse_coords(path: Path) -> np.ndarray:
    """Return Nx3 coords for heavy atoms from PDB/PDBQT (order-agnostic)."""
    coords: List[Tuple[float, float, float]] = []
    try:
        with open(path, "r", errors="ignore") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                atom_name = line[12:16].strip() if len(line) >= 16 else ""
                if atom_name.upper().startswith("H"):
                    continue
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except Exception:
                    parts = line.split()
                    flts = []
                    for p in parts:
                        try:
                            flts.append(float(p))
                        except Exception:
                            pass
                    if len(flts) >= 3:
                        x, y, z = flts[-3:]
                    else:
                        continue
                coords.append((x, y, z))
    except Exception:
        pass
    if not coords:
        return np.zeros((0, 3), dtype=float)
    return np.asarray(coords, dtype=float)


def centroid(pts: np.ndarray) -> Optional[np.ndarray]:
    if pts is None or pts.size == 0:
        return None
    return pts.mean(axis=0)


# -----------------------------
# Approximate RMSD (robust, order-agnostic)
# -----------------------------

def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, float]:
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    P_rot = Pc @ R
    diff2 = ((P_rot - Qc) ** 2).sum(axis=1)
    rmsd = float(np.sqrt(diff2.mean())) if len(diff2) else float("inf")
    return R, rmsd


def approximate_rmsd(coords_a: np.ndarray, coords_b: np.ndarray, nn_cap: int = 128) -> Optional[float]:
    if coords_a.shape[0] < 3 or coords_b.shape[0] < 3:
        return None
    # Use smaller set as reference; cap to nn_cap for speed
    A, B = (coords_a, coords_b) if coords_a.shape[0] <= coords_b.shape[0] else (coords_b, coords_a)
    if A.shape[0] > nn_cap:
        idx = np.linspace(0, A.shape[0] - 1, nn_cap, dtype=int)
        A = A[idx]
    used = set()
    pairs_a: List[np.ndarray] = []
    pairs_b: List[np.ndarray] = []
    for a in A:
        best_j = None
        best_d2 = 1e99
        for j, b in enumerate(B):
            if j in used:
                continue
            d2 = ((a - b) ** 2).sum()
            if d2 < best_d2:
                best_d2 = d2
                best_j = j
        if best_j is None:
            continue
        used.add(best_j)
        pairs_a.append(a)
        pairs_b.append(B[best_j])
    if len(pairs_a) < 3:
        return None
    P = np.vstack(pairs_a)
    Q = np.vstack(pairs_b)
    _, rmsd = _kabsch(P, Q)
    return rmsd


# -----------------------------
# Identity helpers
# -----------------------------

def _rdk_id_from_stem(stem: str) -> Optional[str]:
    m = re.search(r"(rdk_\d{6,8})", stem.lower())
    return m.group(1) if m else None


def resolve_rdk_name_from_mapping(rdk_id: str, mapping: MappingIndex) -> Optional[str]:
    if not mapping:
        return None
    try:
        return mapping.resolve_name_for_rdk(rdk_id)
    except Exception:
        return None


def expected_names_for_het(het: str) -> List[str]:
    return CHEMCOMP_ALIAS.get((het or "").upper().strip(), [])


# -----------------------------
# Data classes
# -----------------------------
@dataclass
class Pose:
    path: Path
    score: Optional[float]
    coords: np.ndarray
    centroid: Optional[np.ndarray]
    rdk_name: Optional[str] = None  # resolved drug name for RDKs


# -----------------------------
# Screenshot helpers (PyMOL renders)
# -----------------------------

def find_pymol_screenshots(pdb_id: str, pocket_dir: Path) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """
    Looks for three PNGs in the given pocket directory:
      <PDB>_cleaned__CONTROL+RDKclosest_side.png
      <PDB>_cleaned__CONTROL+RDKclosest_front.png
      <PDB>_cleaned__CONTROL+RDKclosest_top.png
    Returns (side, front, top) Paths if found, else None entries.
    """
    base = f"{pdb_id}_cleaned__CONTROL+RDKclosest"
    views = ["side", "front", "top"]
    out: List[Optional[Path]] = []
    for v in views:
        p = pocket_dir / f"{base}_{v}.png"
        if p.is_file():
            out.append(p)
            continue
        # Fallback: looser glob in case of minor naming differences
        cands = sorted(pocket_dir.glob(f"*{pdb_id}*CONTROL+RDKclosest*{v}*.png"))
        out.append(cands[0] if cands else None)
    return out[0], out[1], out[2]
def find_pair_only_screenshots(pdb_id: str, pocket_dir: Path, rdk_id: Optional[str]) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """
    Looks for three PNGs named like:
      <PDB>_cleaned__CONTROL+<RDKID>_PAIR_side.png
      <PDB>_cleaned__CONTROL+<RDKID>_PAIR_front.png
      <PDB>_cleaned__CONTROL+<RDKID>_PAIR_top.png
    Falls back to a loose glob if names vary.
    """
    rid = (rdk_id or "RDKclosest")
    base = f"{pdb_id}_cleaned__CONTROL+{rid}_PAIR"
    views = ["side", "front", "top"]
    out: List[Optional[Path]] = []
    for v in views:
        p = pocket_dir / f"{base}_{v}.png"
        if p.is_file():
            out.append(p)
            continue
        # Fallback: allow slight naming differences, but keep rid in the filename if possible
        cands = [q for q in pocket_dir.glob(f"*{pdb_id}*CONTROL+*PAIR*{v}*.png") if rid.lower() in q.name.lower()]
        if not cands:
            cands = list(pocket_dir.glob(f"*{pdb_id}*CONTROL+*PAIR*{v}*.png"))
        out.append(sorted(cands)[0] if cands else None)
    return out[0], out[1], out[2]


@dataclass
class PairEval:
    pdb_id: str
    pocket: str
    control_id: str          # e.g., NIL_A601
    control_het: str         # e.g., NIL
    control_file: Path
    rdk_file: Path
    rdk_id: str
    rdk_name: Optional[str]
    control_score: Optional[float]
    rdk_score: Optional[float]
    delta_kcal: Optional[float]
    centroid_dist: Optional[float]   # value used for scoring (aligned if possible)
    rmsd: Optional[float]
    flag_score: int
    flag_center: int
    flag_rmsd: int
    flag_identity: int
    total_points: int
    n_atoms_ctrl: int
    n_atoms_rdk: int
    centroid_raw: Optional[float]
    centroid_aligned: Optional[float]
    pick_reason: str  # "identity" or "centroid"
    png_side: Optional[Path] = None
    png_front: Optional[Path] = None
    png_top: Optional[Path] = None
    # NEW: interpretability/readability fields
    control_display_name: Optional[str] = None  # guessed control drug name
    flags_sum: int = 0                          # how many flags passed (incl. identity when enabled)
    confidence_label: str = ""                 # confident / plausible / weak
    confident_match: int = 0                    # 1 if confident
    good_control: int = 0                       # 1 if control looks sane (recognized + has score + size)
    png_pair_side: Optional[Path] = None
    png_pair_front: Optional[Path] = None
    png_pair_top: Optional[Path] = None



# -----------------------------
# Crawling & pairing
# -----------------------------

def find_pocket_dirs(pdb_dir: Path) -> List[Path]:
    out: List[Path] = []
    if not pdb_dir.is_dir():
        return out
    for child in pdb_dir.iterdir():
        if child.is_dir() and child.name.startswith("bench_pocket") and child.name.endswith("_single"):
            out.append(child)
    return sorted(out)


def load_controls_and_rdks(pocket_dir: Path, mapping: MappingIndex) -> Tuple[List[Tuple[str, Pose]], List[Pose]]:
    controls: List[Tuple[str, Pose]] = []
    rdks: List[Pose] = []
    for p in pocket_dir.iterdir():
        if not p.is_file():
            continue
        name = p.name
        if CONTROL_PAT.match(name):
            control_id = name.split("_bench_")[0]  # NIL_A601
            coords = parse_coords(p)
            ctr = centroid(coords)
            sc = parse_vina_score(p)
            controls.append((control_id, Pose(path=p, score=sc, coords=coords, centroid=ctr)))
        elif RDK_PAT.match(name):
            coords = parse_coords(p)
            ctr = centroid(coords)
            sc = parse_vina_score(p)
            rdk_id = _rdk_id_from_stem(Path(name).stem)

            rname = resolve_rdk_name_from_mapping(rdk_id or "", mapping)
            rdks.append(Pose(path=p, score=sc, coords=coords, centroid=ctr, rdk_name=rname))
    return controls, rdks


def _closest_by_centroid(ctrl: Pose, rdks: Sequence[Pose]) -> Optional[Pose]:
    if ctrl.centroid is None:
        return None
    best = None
    best_d = 1e99
    for r in rdks:
        if r.centroid is None:
            continue
        d = float(np.linalg.norm(ctrl.centroid - r.centroid))
        if d < best_d:
            best_d = d
            best = r
    return best


def _match_by_identity(ctrl_het: str, rdks: Sequence[Pose]) -> List[Pose]:
    expected_norm = [_norm_text(s) for s in expected_names_for_het(ctrl_het)]
    if not expected_norm:
        return []
    hits: List[Pose] = []
    for r in rdks:
        nm_norm = _norm_text(r.rdk_name)
        if nm_norm and any(en and (en in nm_norm or nm_norm in en) for en in expected_norm):
            hits.append(r)
    return hits


def evaluate_pairs(
    pdb_id: str,
    pocket_name: str,
    pocket_dir: Path,
    controls: Sequence[Tuple[str, Pose]],
    rdks: Sequence[Pose],
    score_tol: float,
    center_tol: float,
    rmsd_tol: float,
    include_identity: bool,
) -> List[PairEval]:
    # Locate the PyMOL screenshots once per pocket
    png_side, png_front, png_top = find_pymol_screenshots(pdb_id, pocket_dir)

    rows: List[PairEval] = []
    for control_id, ctrl in controls:
        m = re.match(r"^(?P<het>[A-Za-z0-9]{3})_(?P<chain>\w)(?P<res>\d+)$", control_id)
        het = (m.group("het") if m else control_id.split("_")[0]).upper()

        # Prefer identity match; fallback to centroid-nearest
        identity_flag = 0
        candidates = _match_by_identity(het, rdks)
        if candidates:
            identity_flag = 1
            picked = (min(candidates, key=lambda r: abs((r.score or 9e9) - (ctrl.score or 0.0)))
                      if ctrl.score is not None else
                      min(candidates, key=lambda r: (r.score or 9e9)))
        else:
            picked = _closest_by_centroid(ctrl, rdks)
        if picked is None:
            continue

        pick_reason = "identity" if identity_flag else "centroid"

        # alias for consistency (fixes NameError in confidence calc)
        flag_identity = identity_flag

        # Score delta
        delta_kcal: Optional[float] = None
        flag_score = 0
        if ctrl.score is not None and picked.score is not None:
            delta_kcal = abs(picked.score - ctrl.score)
            flag_score = int(delta_kcal <= score_tol)

        # Geometry metrics (use ALIGNED centroid for flag_center)
        rmsd_val, raw_ctr, aligned_ctr, n_ctrl, n_rdk = alignment_metrics(ctrl.coords, picked.coords)
        centroid_to_score = aligned_ctr if aligned_ctr is not None else raw_ctr
        flag_center = int(centroid_to_score is not None and centroid_to_score <= center_tol)
        flag_rmsd = int((rmsd_val is not None) and (rmsd_val <= rmsd_tol))

        if DEBUG and rmsd_val is not None and raw_ctr is not None and raw_ctr > (center_tol + 0.25) and flag_rmsd:
            print(f"[debug] {pdb_id}/{pocket_name}/{control_id}: RMSD ok ({rmsd_val:.2f} Å) "
                  f"but RAW centroid {raw_ctr:.2f} Å; ALIGNED centroid {aligned_ctr if aligned_ctr is None else f'{aligned_ctr:.2f}'} Å. "
                  f"Using aligned for scoring. pick={picked.path.name} via {pick_reason}")

        rdk_id = _rdk_id_from_stem(picked.path.stem) or ""
        pair_side, pair_front, pair_top = find_pair_only_screenshots(pdb_id, pocket_dir, rdk_id or "RDKclosest")

        # interpretability / labels
        ctrl_guess_list = expected_names_for_het(het)
        control_display_name = (ctrl_guess_list[0] if ctrl_guess_list else "")
        tri = flag_score + flag_center + flag_rmsd
        confident = int((flag_identity and flag_rmsd and (flag_score or flag_center)) or (tri == 3))
        confidence_label = "confident" if confident else ("plausible" if (tri >= 2 or (flag_identity and (flag_score or flag_center))) else "weak")
        good_control = int((ctrl.score is not None) and (int(ctrl.coords.shape[0]) >= 10) and bool(control_display_name))
        flags_sum = tri + (identity_flag if include_identity else 0)

        total_points = flag_score + flag_center + flag_rmsd + (identity_flag if include_identity else 0)

        rows.append(
            PairEval(
                pdb_id=pdb_id,
                pocket=pocket_name,
                control_id=control_id,
                control_het=het,
                control_file=ctrl.path,
                rdk_file=picked.path,
                rdk_id=rdk_id,
                rdk_name=picked.rdk_name,
                control_score=ctrl.score,
                rdk_score=picked.score,
                delta_kcal=delta_kcal,
                centroid_dist=centroid_to_score,   # value used for scoring
                rmsd=rmsd_val,
                flag_score=flag_score,
                flag_center=flag_center,
                flag_rmsd=flag_rmsd,
                flag_identity=identity_flag,
                total_points=total_points,
                n_atoms_ctrl=int(ctrl.coords.shape[0]),
                n_atoms_rdk=int(picked.coords.shape[0]),
                centroid_raw=raw_ctr,
                centroid_aligned=aligned_ctr,
                pick_reason=pick_reason,
                png_side=png_side,
                png_front=png_front,
                png_top=png_top,
                control_display_name=control_display_name,
                flags_sum=flags_sum,
                confidence_label=confidence_label,
                confident_match=confident,
                good_control=good_control,
                png_pair_side=pair_side,
                png_pair_front=pair_front,
                png_pair_top=pair_top,

            )
        )
    return rows


# -----------------------------
# CSV writing
# -----------------------------

def _fmt(x: Optional[float]) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return ""
    try:
        return f"{float(x):.3f}"
    except Exception:
        return ""


def _pass_category(total_points: int, include_identity: bool) -> str:
    max_points = 4 if include_identity else 3
    if total_points == max_points:
        return "pass"
    if 1 < total_points < max_points:
        return "near"
    if total_points == 0:
        return "fail"
    return ""  # exactly 1 -> no category requested


def write_details_csv(out_dir: Path, rows: Sequence[PairEval], include_identity: bool) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "benchmark_analysis_details.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "pdb_id", "pocket", "control_id", "control_het", "control_file", "rdk_file", "rdk_id",
            "control_drug_guess", "rdk_suspected_fda_name",
            "control_score_kcal", "rdk_score_kcal", "delta_kcal",
            "centroid_raw_A", "centroid_aligned_A", "centroid_dist_A",
            "rmsd_ctrl_to_rdk_A", "rmsd_rdk_to_ctrl_A",
            "score_within_tol", "centroid_within_tol", "rmsd_within_tol", "identity_match",
            "n_atoms_ctrl", "n_atoms_rdk", "pick_reason",
            "png_side", "png_front", "png_top",
            "pair_png_side", "pair_png_front", "pair_png_top",
            "flags_sum", "confident_match", "confidence_label", "good_control",
            "total_points", "pass_category",
        ])
        for r in rows:
            w.writerow([
                r.pdb_id, r.pocket, r.control_id, r.control_het,
                str(r.control_file), str(r.rdk_file), r.rdk_id,
                (r.control_display_name or ""), (r.rdk_name or ""),
                _fmt(r.control_score), _fmt(r.rdk_score), _fmt(r.delta_kcal),
                _fmt(r.centroid_raw), _fmt(r.centroid_aligned), _fmt(r.centroid_dist),
                _fmt(r.rmsd), _fmt(r.rmsd),
                r.flag_score, r.flag_center, r.flag_rmsd, r.flag_identity,
                r.n_atoms_ctrl, r.n_atoms_rdk, r.pick_reason,
                (str(r.png_side) if r.png_side else ""),
                (str(r.png_front) if r.png_front else ""),
                (str(r.png_top) if r.png_top else ""),
                (str(r.png_pair_side) if r.png_pair_side else ""),
                (str(r.png_pair_front) if r.png_pair_front else ""),
                (str(r.png_pair_top) if r.png_pair_top else ""),
                r.flags_sum, r.confident_match, r.confidence_label, r.good_control,
                r.total_points, _pass_category(r.total_points, include_identity),
            ])
    return path


# =========================
# Summary writers (drop-in)
# =========================

def write_summary_csv(out_dir: Path, rows: Sequence[PairEval], include_identity: bool) -> Path:
    """
    Per-PDB rollup with flags and counts, including 'good_control' and 'confident_match'.
    Writes: <out_dir>/benchmark_analysis_summary.csv
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    max_points = 4 if include_identity else 3
    path = out_dir / "benchmark_analysis_summary.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "pdb_id",
            "best_total_points",
            "num_controls",
            "num_pockets",
            "any_score_within_tol",
            "any_centroid_within_tol",
            "any_rmsd_within_tol",
            "any_identity_match",
            "num_good_controls",
            "num_confident_matches",
            "passed_4_or_max",
            "near_pass_2to(max-1)",
            "fail_0",
        ])
        by_pdb: Dict[str, List[PairEval]] = {}
        for r in rows:
            by_pdb.setdefault(r.pdb_id, []).append(r)

        for pdb_id in sorted(by_pdb.keys()):
            items = by_pdb[pdb_id]
            best_total = max((it.total_points for it in items), default=0)
            num_controls = len({it.control_id for it in items})
            num_pockets = len({it.pocket for it in items})
            any_score = int(any(it.flag_score for it in items))
            any_center = int(any(it.flag_center for it in items))
            any_rmsd = int(any(it.flag_rmsd for it in items))
            any_ident = int(any(it.flag_identity for it in items))
            num_good_ctrl = sum(int(it.good_control) for it in items)
            num_confident = sum(int(it.confident_match) for it in items)
            passed = int(best_total == max_points)
            near = int(1 < best_total < max_points)
            fail0 = int(best_total == 0)

            w.writerow([
                pdb_id, best_total, num_controls, num_pockets,
                any_score, any_center, any_rmsd, any_ident,
                num_good_ctrl, num_confident,
                passed, near, fail0,
            ])
    return path


# -----------------------------
# Visual reports (images embedded)
# -----------------------------

def _html_escape(s: Optional[str]) -> str:
    return _html.escape("" if s is None else str(s))


def _img_to_data_uri_or_link(p: Optional[Path], max_width_px: int = 280) -> str:
    if not p or not Path(p).is_file():
        return ""
    try:
        data = Path(p).read_bytes()
        b64 = base64.b64encode(data).decode("ascii")
        return f'<img src="data:image/png;base64,{b64}" style="max-width:{max_width_px}px; height:auto; border:1px solid #ddd;" />'
    except Exception:
        try:
            return f'<img src="{Path(p).as_uri()}" style="max-width:{max_width_px}px; height:auto; border:1px solid #ddd;" />'
        except Exception:
            return _html_escape(str(p))


def write_details_html(out_dir: Path, rows: Sequence[PairEval], include_identity: bool) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "benchmark_analysis_details_with_images.html"
    def _badge(label: str) -> str:
        cls = "bad"
        if label == "confident":
            cls = "good"
        elif label == "plausible":
            cls = "ok"
        return f"<span class='badge {cls}'>{_html_escape(label)}</span>"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("""
<!doctype html>
<html>
<head>
<meta charset="utf-8" />
<title>Benchmark Analysis — Details (with images)</title>
<style>
  body{font-family:Segoe UI,Arial,sans-serif;margin:24px}
  table{border-collapse:collapse;width:100%}
  th,td{border:1px solid #e5e5e5;padding:6px 8px;vertical-align:top;font-size:13px}
  th{background:#fafafa;position:sticky;top:0;z-index:1}
  .num{text-align:right;white-space:nowrap}
  .imgcell{white-space:nowrap}
  .muted{color:#777}
  .badge{display:inline-block;padding:2px 8px;border-radius:12px;border:1px solid #ddd;font-size:12px}
  .good{background:#e7f6ec;border-color:#c9e8d0;color:#137333}
  .ok{background:#fff8e1;border-color:#f7e0a3;color:#8a6d3b}
  .bad{background:#fdecea;border-color:#f5c6cb;color:#a61b1b}
</style>
</head>
<body>
<h2>Benchmark Analysis — Details (with images)</h2>
<p class="muted">Columns mirror <code>benchmark_analysis_details.csv</code>, with extra readability fields (scores for control & RDK, suspected FDA name, confidence, good-control) and <b>six image columns</b>. Generated by the script.</p>
<table>
<thead><tr>
  <th>pdb_id</th><th>pocket</th><th>control_id</th><th>rdk_id</th>
  <th>control_drug</th><th>rdk_suspected_fda</th>
  <th class="num">ctrl_score</th><th class="num">rdk_score</th><th class="num">Δkcal</th>
  <th class="num">centroid_A</th><th class="num">RMSD_A</th>
  <th>pick</th><th>confidence</th><th>good control</th>
  <th class="num">score✔</th><th class="num">center✔</th><th class="num">rmsd✔</th><th class="num">ident✔</th>
  <th class="imgcell">side</th><th class="imgcell">front</th><th class="imgcell">top</th>
  <th class="imgcell">pair_side</th><th class="imgcell">pair_front</th><th class="imgcell">pair_top</th>
</tr></thead>
<tbody>
""")
        for r in rows:
            fh.write("<tr>")
            fh.write(f"<td>{_html_escape(r.pdb_id)}</td>")
            fh.write(f"<td>{_html_escape(r.pocket)}</td>")
            fh.write(f"<td>{_html_escape(r.control_id)}</td>")
            fh.write(f"<td>{_html_escape(r.rdk_id)}</td>")
            fh.write(f"<td>{_html_escape(r.control_display_name or '')}</td>")
            fh.write(f"<td>{_html_escape(r.rdk_name or '')}</td>")
            fh.write(f"<td class='num'>{_html_escape(_fmt(r.control_score))}</td>")
            fh.write(f"<td class='num'>{_html_escape(_fmt(r.rdk_score))}</td>")
            fh.write(f"<td class='num'>{_html_escape(_fmt(r.delta_kcal))}</td>")
            fh.write(f"<td class='num'>{_html_escape(_fmt(r.centroid_dist))}</td>")
            fh.write(f"<td class='num'>{_html_escape(_fmt(r.rmsd))}</td>")
            fh.write(f"<td>{_html_escape(r.pick_reason)}</td>")
            fh.write(f"<td>{_badge(r.confidence_label)}</td>")
            fh.write(f"<td class='num'>{'✔' if r.good_control else '—'}</td>")
            fh.write(f"<td class='num'>{r.flag_score}</td>")
            fh.write(f"<td class='num'>{r.flag_center}</td>")
            fh.write(f"<td class='num'>{r.flag_rmsd}</td>")
            fh.write(f"<td class='num'>{r.flag_identity}</td>")
            # overlay images
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(r.png_side)}</td>")
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(r.png_front)}</td>")
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(r.png_top)}</td>")
            # pair-only images
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(r.png_pair_side)}</td>")
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(r.png_pair_front)}</td>")
            fh.write(f"<td class='imgcell'>{_img_to_data_uri_or_link(r.png_pair_top)}</td>")
            fh.write("</tr>")
        fh.write("""
</tbody>
</table>
</body>
</html>
""")
    return path

def write_details_xlsx(out_dir: Path, rows: Sequence[PairEval], include_identity: bool) -> Optional[Path]:
    """Try to write an .xlsx with embedded images. Uses xlsxwriter if available, else openpyxl. Returns path or None."""
    # First try xlsxwriter (no Pillow dependency required)
    try:
        import xlsxwriter  # type: ignore
        xlsx_path = out_dir / "benchmark_analysis_details.xlsx"
        wb = xlsxwriter.Workbook(str(xlsx_path))
        ws = wb.add_worksheet("details")
        # Column headers
        headers = [
            "pdb_id","pocket","control_id","control_het","rdk_id","rdk_suspected_fda","control_drug_guess",
            "ctrl_score","rdk_score","delta_kcal","centroid_A","rmsd_A","pick_reason",
            "confidence","good_control","score✔","center✔","rmsd✔","ident✔",
            "image_side","image_front","image_top","pair_side","pair_front","pair_top",
        ]
        for c,h in enumerate(headers):
            ws.write(0, c, h)
        # Set widths
        ws.set_column(0, 1, 10)
        ws.set_column(2, 3, 16)
        ws.set_column(4, 6, 20)
        ws.set_column(7, 12, 12)
        ws.set_column(13, 18, 11)
        ws.set_column(19, 24, 32)  # cover all 6 image columns
        # Data rows
        row = 1
        for r in rows:
            ws.write_row(row, 0, [
                r.pdb_id, r.pocket, r.control_id, r.control_het, r.rdk_id, r.rdk_name or "", r.control_display_name or "",
                _fmt(r.control_score), _fmt(r.rdk_score), _fmt(r.delta_kcal), _fmt(r.centroid_dist), _fmt(r.rmsd), r.pick_reason,
                r.confidence_label, r.good_control, r.flag_score, r.flag_center, r.flag_rmsd, r.flag_identity
            ])
            # Insert overlay images (scaled down)
            for offset, p in enumerate([r.png_side, r.png_front, r.png_top]):
                if p and Path(p).is_file():
                    try:
                        ws.insert_image(row, 19 + offset, str(p), {"x_scale":0.35, "y_scale":0.35})
                    except Exception:
                        pass
            # Insert pair-only images (scaled down)
            for offset, p in enumerate([r.png_pair_side, r.png_pair_front, r.png_pair_top]):
                if p and Path(p).is_file():
                    try:
                        ws.insert_image(row, 22 + offset, str(p), {"x_scale":0.35, "y_scale":0.35})
                    except Exception:
                        pass
            row += 1
        wb.close()
        return xlsx_path
    except Exception:
        pass

    # Fallback: openpyxl (requires Pillow for PNG sizing)
    try:
        from openpyxl import Workbook  # type: ignore
        from openpyxl.drawing.image import Image as XLImage  # type: ignore
        from openpyxl.utils import get_column_letter  # type: ignore
        xlsx_path = out_dir / "benchmark_analysis_details.xlsx"
        wb = Workbook()
        ws = wb.active
        ws.title = "details"
        headers = [
            "pdb_id","pocket","control_id","control_het","rdk_id","rdk_suspected_fda","control_drug_guess",
            "ctrl_score","rdk_score","delta_kcal","centroid_A","rmsd_A","pick_reason",
            "confidence","good_control","score✔","center✔","rmsd✔","ident✔",
            "image_side","image_front","image_top","pair_side","pair_front","pair_top",
        ]
        ws.append(headers)
        for r in rows:
            ws.append([
                r.pdb_id, r.pocket, r.control_id, r.control_het, r.rdk_id, r.rdk_name or "", r.control_display_name or "",
                _fmt(r.control_score), _fmt(r.rdk_score), _fmt(r.delta_kcal), _fmt(r.centroid_dist), _fmt(r.rmsd), r.pick_reason,
                r.confidence_label, r.good_control, r.flag_score, r.flag_center, r.flag_rmsd, r.flag_identity,
                "", "", "", "", "", ""  # placeholders for 6 images
            ])
            row_idx = ws.max_row
            # overlay
            for i, p in enumerate([r.png_side, r.png_front, r.png_top], start=20):
                try:
                    if p and Path(p).is_file():
                        img = XLImage(str(p))
                        ws.add_image(img, f"{get_column_letter(i)}{row_idx}")
                except Exception:
                    continue
            # pair-only
            for i, p in enumerate([r.png_pair_side, r.png_pair_front, r.png_pair_top], start=23):
                try:
                    if p and Path(p).is_file():
                        img = XLImage(str(p))
                        ws.add_image(img, f"{get_column_letter(i)}{row_idx}")
                except Exception:
                    continue
        wb.save(str(xlsx_path))
        return xlsx_path
    except Exception:
        return None



def write_details_visual_report(out_dir: Path, rows: Sequence[PairEval], include_identity: bool) -> Path:
    """Create a visual report attempting XLSX (with images) first, else HTML with inline images."""
    xlsx = write_details_xlsx(out_dir, rows, include_identity)
    if xlsx is not None:
        return xlsx
    return write_details_html(out_dir, rows, include_identity)


# -----------------------------
# Main
# -----------------------------

def find_pdb_dirs(docked_root: Path, only_pdb: Optional[str]) -> List[Path]:
    pdb_dirs: List[Path] = []
    for child in docked_root.iterdir():
        if not child.is_dir():
            continue
        if only_pdb and child.name.upper() != only_pdb.upper():
            continue
        pdb_dirs.append(child)
    pdb_dirs.sort()
    return pdb_dirs


def run_analysis(
    docked_root: Path,
    mapping_csv: Optional[Path],
    only_pdb: Optional[str],
    score_tol: float,
    center_tol: float,
    rmsd_tol: float,
    include_identity: bool,
    out_dir: Optional[Path] = None,
) -> Tuple[Optional[Path], Optional[Path]]:
    if not docked_root.is_dir():
        print(f"[analysis] DOCKED root not found: {docked_root}")
        return None, None

    mapping = MappingIndex(mapping_csv if mapping_csv and Path(mapping_csv).is_file() else None)

    all_rows: List[PairEval] = []
    for pdb_dir in find_pdb_dirs(docked_root, only_pdb):
        pdb_id = pdb_dir.name.upper()
        for pocket_dir in find_pocket_dirs(pdb_dir):
            controls, rdks = load_controls_and_rdks(pocket_dir, mapping)
            if not controls or not rdks:
                if DEBUG:
                    print(f"[debug] Skipping {pdb_id}/{pocket_dir.name} (controls={len(controls)} rdks={len(rdks)})")
                continue
            rows = evaluate_pairs(
                pdb_id,
                pocket_dir.name,
                pocket_dir,  # pass pocket path for PNG discovery
                controls,
                rdks,
                score_tol=score_tol,
                center_tol=center_tol,
                rmsd_tol=rmsd_tol,
                include_identity=include_identity,
            )
            all_rows.extend(rows)

    if not all_rows:
        print("[analysis] No matches found. Are the bench_pocketX_single folders populated?")
        return None, None

    out_root = out_dir if out_dir else (docked_root / "_analysis")
    details_path = write_details_csv(out_root, all_rows, include_identity)
    summary_path = write_summary_csv(out_root, all_rows, include_identity)
    # Also write a *visual* details report that actually shows images (XLSX if possible, else HTML)
    visual_path = write_details_visual_report(out_root, all_rows, include_identity)

    # ---- Nice console rollup: list PDBs by category + global best ----
    max_points = 4 if include_identity else 3
    by_pdb_best: Dict[str, int] = {}
    for r in all_rows:
        by_pdb_best[r.pdb_id] = max(by_pdb_best.get(r.pdb_id, 0), r.total_points)

    passed_pdbs = sorted([p for p, v in by_pdb_best.items() if v == max_points])
    near_pdbs   = sorted([p for p, v in by_pdb_best.items() if 1 < v < max_points])
    fail0_pdbs  = sorted([p for p, v in by_pdb_best.items() if v == 0])

    global_best = max(by_pdb_best.values()) if by_pdb_best else 0
    best_pdbs   = sorted([p for p, v in by_pdb_best.items() if v == global_best])

    print(f"[analysis] Wrote details CSV: {details_path}")
    print(f"[analysis] Wrote summary CSV: {summary_path}")
    print(f"[analysis] Wrote visual details: {visual_path}")
    print(f"[analysis] Passed (=={max_points}): {len(passed_pdbs)} -> {', '.join(passed_pdbs) if passed_pdbs else '-'}")
    print(f"[analysis] Near (2..{max_points-1}): {len(near_pdbs)} -> {', '.join(near_pdbs) if near_pdbs else '-'}")
    print(f"[analysis] Fail (==0): {len(fail0_pdbs)} -> {', '.join(fail0_pdbs) if fail0_pdbs else '-'}")
    print(f"[analysis] PDBs with BEST score = {global_best}: {', '.join(best_pdbs) if best_pdbs else '-'}")

    return details_path, summary_path


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Auto-interpret benchmark results across proteins/pockets.")
    p.add_argument("--docked-root", default=DEFAULT_DOCKED_ROOT)
    p.add_argument("--mapping", default=DEFAULT_MAPPING_CSV)
    p.add_argument("--only-pdb", default=DEFAULT_ONLY_PDB)
    p.add_argument("--score-tol", type=float, default=DEFAULT_SCORE_TOL)
    p.add_argument("--center-tol", type=float, default=DEFAULT_CENTER_TOL)
    p.add_argument("--rmsd-tol", type=float, default=DEFAULT_RMSD_TOL)
    p.add_argument("--no-identity", action="store_true", help="Exclude identity-match from scoring (max 3 points)")
    p.add_argument("--out", default=DEFAULT_OUTDIR)
    p.add_argument("--debug", action="store_true", help="Verbose per-pair diagnostics (or set BENCH_DEBUG=1)")
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_argparser().parse_args(argv)
    docked_root = Path(args.docked_root)
    mapping_csv = Path(args.mapping) if args.mapping else None
    out = Path(args.out) if args.out else None

    # honor CLI --debug
    global DEBUG
    DEBUG = DEBUG or bool(args.debug)

    run_analysis(
        docked_root=docked_root,
        mapping_csv=mapping_csv,
        only_pdb=(args.only_pdb or None),
        score_tol=float(args.score_tol),
        center_tol=float(args.center_tol),
        rmsd_tol=float(args.rmsd_tol),
        include_identity=not bool(args.no_identity),
        out_dir=out,
    )


if __name__ == "__main__":
    # Allow running directly in an IDE without supplying args
    try:
        main(None)
    except SystemExit:
        run_analysis(
            docked_root=Path(DEFAULT_DOCKED_ROOT),
            mapping_csv=Path(DEFAULT_MAPPING_CSV) if DEFAULT_MAPPING_CSV else None,
            only_pdb=DEFAULT_ONLY_PDB,
            score_tol=DEFAULT_SCORE_TOL,
            center_tol=DEFAULT_CENTER_TOL,
            rmsd_tol=DEFAULT_RMSD_TOL,
            include_identity=DEFAULT_INCLUDE_IDENTITY,
            out_dir=Path(DEFAULT_OUTDIR) if DEFAULT_OUTDIR else None,
        )
