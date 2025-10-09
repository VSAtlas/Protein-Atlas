import sys
import ctypes
import os
from ctypes import wintypes, create_unicode_buffer
from pathlib import Path
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Tuple, Dict
import re
from datetime import datetime
import shutil
from collections import defaultdict

from activesite import (
    fix_pdb_elements,
    fix_element_columns_in_file,  # unify element rewrite for columns 77–78
    _fix_ligand_element_columns_in_memory,
    derive_element,
    scan_helium_counts,
    rules_version,
    assert_no_helium_in_hydrogen_names,
    assert_no_helium_in_pdbqt,
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

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "SEC", "PYL", "MSE"
}

# Exclude common crystallization additives/buffers/metals (+ porphyrins / modified residues)
EXCLUDE_CRYSTAL_ADDITIVES = {
    "HOH", "CIT", "TAR", "SO4", "PO4", "CA", "NA", "K", "MG", "MN", "ZN",
    "GOL", "EDO", "PEG", "MPD", "TRS", "MES", "HEPES", "ACET", "ACT", "FMT",
    "MAL", "DMS", "IPA", "CLU", "NAG", "BOG",
    # common counter-ion/salt codes
    "TOS", "BES", "PTS", "OTF", "TRF", "TFA", "BF4", "PF6", "CL", "BR", "I",
    # porphyrins / heme family often not intended as small-mol ligands
    "HEM", "HEC", "HEA", "HEB", "HEO", "HEG", "HEF", "HEH",
    # phosphorylated residues frequently side-chain mods, not ligands to dock
    "PTR", "TPO", "SEP"
}

# --- Salvage / logging config ---
RUN_TAG = datetime.now().strftime("%Y%m%d_%H%M%S")
MALFORMED_LOG = Path(f"malformed_ligands_{RUN_TAG}.txt")
QUARANTINE_DIRNAME = "quarantine"

# --- Resume mode toggle ---
RESUME_SKIP = False

ALLOWED_ELEMENTS = {
    "H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "B", "Si", "Se", "Zn", "Mg", "Ca", "Mn", "Fe", "K", "Na"
}
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
# Utility / logging helpers
# =========================
MONOATOMIC_IONS = {"Cl", "Br", "I", "F", "Na", "K", "Ca", "Mg", "Zn", "Mn", "Fe", "Cu", "Co", "Ni", "Al", "Ag", "Au",
                   "Pt", "Li", "Ba", "Sr", "Cs", "Rb"}


def _looks_like_monoatomic_ion_pdbqt(lines: List[str]) -> bool:
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    if len(atom_lines) != 1:
        return False
    parts = atom_lines[0].split()
    elem = _element_from_adt(parts[-1]) if parts else ""
    return elem in MONOATOMIC_IONS


def _log_malformed(path: Path, reason: str):
    try:
        MALFORMED_LOG.parent.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    with open(MALFORMED_LOG, "a", encoding="utf-8") as fh:
        fh.write(f"{path}\t{reason}\n")


def _quick_filters(mol: Chem.Mol) -> Tuple[bool, str]:
    heavy = mol.GetNumHeavyAtoms()
    if heavy == 0:
        return False, "no_heavy_atoms"
    if heavy > MAX_HEAVY_ATOMS:
        return False, f"too_large({heavy})"
    if mol.GetNumAtoms() < MIN_ATOMS_FOR_DOCKING:
        return False, f"too_few_atoms({mol.GetNumAtoms()})"
    for a in mol.GetAtoms():
        if a.GetSymbol() not in ALLOWED_ELEMENTS:
            return False, f"disallowed_element:{a.GetSymbol()}"
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
        Chem.SanitizeMol(m)
        # Force RDKit model and keep aromatic flags on output
        Chem.SetAromaticity(m, Chem.AromaticityModel.AROMATICITY_RDKIT)
        w = Chem.SDWriter(str(out_path))
        try:
            w.SetKekulize(False)  # keep aromatic bonds marked
        except Exception:
            pass
        w.write(m);
        w.close()
        return True
    except Exception:
        return False


def _re_aromatize_mol2_in_place(mol2_path: Path, obabel_exe_short: str) -> bool:
    """Read MOL2 → impose RDKit aromaticity → write SDF (aromatic) → back to MOL2 (overwrite)."""
    try:
        m = Chem.MolFromMol2File(str(mol2_path), sanitize=True, removeHs=False)
    except Exception:
        m = None
    if m is None:
        return False

    tmp_sdf = mol2_path.with_suffix(".arom.sdf")
    if not _write_aromatic_sdf(m, tmp_sdf):
        return False

    tmp_mol2 = mol2_path.with_suffix(".arom.mol2")
    cmd = [obabel_exe_short, "-isdf", str(tmp_sdf), "-omol2", "-O", str(tmp_mol2)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
    except subprocess.CalledProcessError:
        return False
    finally:
        try:
            tmp_sdf.unlink(missing_ok=True)
        except Exception:
            pass

    if tmp_mol2.exists() and tmp_mol2.stat().st_size > 100:
        try:
            mol2_path.unlink(missing_ok=True)
            tmp_mol2.replace(mol2_path)
            return True
        except Exception:
            return False
    return False


# --- NEW: filename and H helpers ---
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
        rules_ver: str = ""
):
    """
    Backward-compatible TSV: old 4 columns + appended provenance columns.
    """
    header = (
        "ligand\tstatus\treason\tpdbqt_rel\t"
        "stage\tfailure_code\tfailure_detail\tfixes_count\trules_version\n"
    )
    if not status_log_path.exists():
        status_log_path.write_text(header, encoding="utf-8")
    with status_log_path.open("a", encoding="utf-8") as fh:
        row = [
            ligand_name, status, reason, relpath,
            stage, failure_code, failure_detail, str(fixes_count), rules_ver
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
    """
    Ensure no ADT He remains in ligand PDBQT.
    Returns (ok, fixes_count, quarantine_reason).
    """
    try:
        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
        new_lines, fixes, q_reason = assert_no_helium_in_pdbqt(lines, ligand_name)
        logging.info("[helium] ligand=%s fixes=%d quarantine=%s rules=%s",
                     mol2_file.name, fixes, str(bool(q_reason)), rules_version())

        if fixes > 0:
            with open(pdbqt_path, "w", encoding="utf-8") as out:
                out.writelines(new_lines)

        ok = (q_reason == "")
        print(f"[helium] ligand={ligand_name} fixes={fixes} quarantine={not ok} rules={_rules_version()}")

        if not ok and "inconsistent" in q_reason:
            _write_status(ligand_name, failure_code="adt_helium_inconsistent", reason=q_reason)

        return ok, fixes, q_reason

    except Exception as e:
        print(
            f"[helium] ligand={ligand_name} fixes=0 quarantine=True rules={_rules_version()} note=postwrite_exception:{e}")
        return False, 0, f"postwrite_exception:{e}"


# =========================
# OBabel-friendly SDF writer
# =========================

def _write_obabel_friendly_sdf(mol: Chem.Mol, out_path: Path) -> bool:
    try:
        m = Chem.Mol(mol)
        Chem.SanitizeMol(m)
        Chem.Kekulize(m, clearAromaticFlags=True)
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
        _log_malformed(target_mol2, "rdkit_sdf_write_fail:kekulize_or_aromaticity")
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
        ok = _attempt_obabel_series(
            base_cmd, timeout_sec=timeout_sec, threads_list=threads_list, use_fast_first=True
        )
        if not ok:
            print(f"Chunk {ci} failed entirely; moving on.")
            continue

        out_files = sorted(mol2_output_dir.glob(f"mol2_chunk{ci:04d}_*.mol2"))
        print(f"Chunk {ci}: wrote {len(out_files)} mol2 files")
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
        try:
            Chem.SanitizeMol(m)
        except Exception:
            return None

        orig_heavy = m.GetNumHeavyAtoms()
        p = _to_parent_mol(m)
        if p is None:
            _log_malformed(Path(f"rdk_{i:07d}"), "no_parent_or_too_small_after_desalting")
            return None
        if p is not m and p.GetNumHeavyAtoms() != orig_heavy:
            logging.info(f"[parent-pick] rdk_{i:07d}: {orig_heavy}?{p.GetNumHeavyAtoms()} heavy atoms")

        m = p

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
        # BUGFIX: check the actual output file (not 'pdbqt_path')
        if out.exists() and out.stat().st_size > 100:
            out_files.append(out)
            continue
        cmd = [obabel_exe, "-isdf", str(pth), "-omol2", "-O", str(out)]
        if _run_obabel(cmd, timeout_sec=120):
            out_files.append(out)
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


# === AROMATICITY AUDIT & RESCUE ===

def _count_aromatic_atoms_in_mol2(mol2_path: Path) -> int:
    m = Chem.MolFromMol2File(str(mol2_path), sanitize=True, removeHs=False)
    if m is None:
        return -1
    return sum(int(a.GetIsAromatic()) for a in m.GetAtoms())


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
    tmp = out_pdbqt.with_suffix(".obabel_tmp.pdbqt")
    cmd = [obabel_exe_short, "-imol2", str(mol2_file), "-opdbqt", "-O", str(tmp)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
    except subprocess.CalledProcessError:
        return False
    if tmp.exists() and tmp.stat().st_size > 100:
        try:
            tmp.replace(out_pdbqt)
            return True
        except Exception:
            return False
    return False


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
    lig_id = mol2_file.stem  # stable ID for logs/TSV
    try:
        if pdbqt_path.exists():
            pdbqt_path.unlink()
    except Exception:
        pass
    # Ensure explicit hydrogens: mol2_in -> mol2_h
    mol2_in = Path(mol2_file)
    mol2_h = mol2_in.with_name(mol2_in.stem + ".withH.mol2")
    okH, hstderr = add_hydrogens_mol2(mol2_in, mol2_h, obabel_exe_short)
    logging.info("[ligprep] addHs stage=primary in=%s ok=%s", mol2_in.name, okH)
    if (hstderr or "").strip():
        logging.warning("[ligprep] addHs stderr (primary) %s", hstderr.splitlines()[-1][:200])

    # Use mol2_h as input to MGLTools if AddHs succeeded
    mol2_for_mgl = mol2_h if okH else mol2_in

    # aromatic baseline from source MOL2 (keep original for comparison)
    src_arom = _count_aromatic_atoms_in_mol2(mol2_file)

    # use the H-enriched file for MGLTools (cwd is mol2_file.parent)
    mol2_short = mol2_for_mgl.name
    pdbqt_short = get_short_path_name(str(pdbqt_path.resolve()))
    cmd = [
        mgltools_python_short,
        prepare_script_short,
        "-l", mol2_short,
        "-o", pdbqt_short,
        "-U", "nphs_lps",
        "-A", "hydrogens",
    ]

    try:
        _ = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            cwd=str(mol2_file.parent),
            timeout=600
        )

        # Compare aromatic counts (MOL2 vs PDBQT)
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
        try:
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
                lines = fh.readlines()
            new_lines, fixes, q_reason = assert_no_helium_in_pdbqt(lines, mol2_file.name)
            if q_reason:
                # Quarantine and TSV provenance
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
                    reason="adt_helium_inconsistent",
                    relpath=str(qpath.name),
                    stage="validate_pdbqt",
                    failure_code="adt_helium_inconsistent",
                    failure_detail=q_reason,
                    fixes_count=0,
                    rules_ver=rules_version(),
                )
                logging.warning(
                    f"[elem-validate] file={pdbqt_path.name} stage=postwrite quarantine ligand={mol2_file.name} code=adt_helium_inconsistent")
                return (mol2_file.name, "postcheck_fail")
            if fixes > 0:
                # Overwrite file with corrections
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
                    f"[elem-validate] file={pdbqt_path.name} stage=postwrite ok ligand={mol2_file.name} fixes={fixes}")
            else:
                # still append a line so provenance is present on success
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
            logging.warning(f"[elem-validate] failed to post-check PDBQT for {mol2_file.name}: {e}")

        # --- lightweight post-write PDBQT validity gate (primary) ---
        okV, n_atoms, whyV = quick_pdbqt_validate(pdbqt_path)
        logging.info("[ligprep] validate_pdbqt result=%s atoms=%d ligand=%s",
                     "ok" if okV else "fail", n_atoms, mol2_file.stem)
        if not okV:
            # quarantine + TSV with failure_code
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
        return (mol2_file.name, "timeout")
    except subprocess.CalledProcessError as e:
        return (mol2_file.name, f"prepare_fail:{e.returncode}")


def load_mol2_lenient(path, logger):
    mol = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
    if mol is None:
        _log_malformed(Path(path), "rdkit_read_fail")
        return None
    try:
        Chem.SanitizeMol(mol)
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
    if win_probe.exists():
        return win_probe
    return lin_probe


def prep_ligands_from_pdb(ligand_output_dir: Path, ligands_mol2_dir: Path, prepped_ligands_dir: Path):
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

    pdb_files = list(ligand_output_dir.glob("*.pdb"))
    logging.info(f"Found {len(pdb_files)} PDB ligand file(s)")

    prepped_ligands_dir.mkdir(parents=True, exist_ok=True)
    status_log = prepped_ligands_dir / "ligand_prep_status.tsv"

    # make sanitized/MOL2 visible next to outputs for easier debugging
    # Create intermediates dir and a symlink to the source ligands (no per-ligand vars here)
    try:
        inter_dir = prepped_ligands_dir / "intermediates"
        inter_dir.mkdir(parents=True, exist_ok=True)
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
        logging.warning("Could not create intermediates symlink: %s", e)

    def _sanitize_pdb(in_pdb: Path, out_pdb: Path) -> bool:
        """
        Keep ONLY the target ligand residue (by resName/chain/resSeq) and drop all protein ATOM lines.
        Also retains headers/remarks/END lines.
        """
        # discover the primary ligand residue triplet from the first HETATM
        target = None  # (resn, chain, resi)
        wrote_any = False

        def _triplet(ln: str):
            # PDB fixed columns: resName 17-20, chainID 21, resSeq 22-26
            resn = (ln[17:20] if len(ln) >= 20 else "").strip()
            chain = (ln[21] if len(ln) >= 22 else " ").strip()
            resi = (ln[22:26] if len(ln) >= 26 else "").strip()
            return resn, chain, resi

        with open(in_pdb, "r", encoding="utf-8", errors="ignore") as fin:
            lines = fin.readlines()

        # first pass: find the first HETATM triplet to define the ligand
        for ln in lines:
            if ln.startswith("HETATM"):
                target = _triplet(ln)
                break

        if target:
            logging.info(f"[sanitize] {in_pdb.name}: keeping ligand {target[0]} chain {target[1]} resi {target[2]}")
        else:
            logging.warning(f"[sanitize] {in_pdb.name}: no HETATM found; nothing to keep")

        with open(out_pdb, "w", encoding="utf-8") as fout:
            for ln in lines:
                if not ln.startswith(("ATOM", "HETATM")):
                    fout.write(ln)
                    continue

                # drop all protein ATOM lines unconditionally
                if ln.startswith("ATOM"):
                    continue

                # HETATM: keep only if it matches the target ligand residue
                if target is None:
                    # if we didn't find a HETATM earlier, be conservative: drop
                    continue

                resn, chain, resi = _triplet(ln)
                if (resn, chain, resi) != target:
                    continue

                # basic altLoc/coordinate sanity (your original checks)
                altloc = ln[16] if len(ln) > 16 else " "
                keep_alt = {"", " ", "A"} | set("0123456789")
                a = altloc.strip()
                if a and (a not in keep_alt and a.upper() not in keep_alt):
                    continue
                try:
                    x = float(ln[30:38])
                    y = float(ln[38:46])
                    z = float(ln[46:54])
                    if any([(x != x), (y != y), (z != z)]) or max(abs(x), abs(y), abs(z)) > 1e6:
                        continue
                except Exception:
                    continue

                fout.write(ln)
                wrote_any = True

        if not wrote_any:
            logging.warning(f"[control-prep] {in_pdb.name}: no HETATM kept after residue filtering.")
        # require that we actually kept >=1 atom (not just headers/remarks)
        return wrote_any and out_pdb.exists() and out_pdb.stat().st_size > 0

    for pdb_file in pdb_files:
        logging.info(f"Processing: {pdb_file.name}")

        residue_name = pdb_file.stem.split("_")[0].upper()
        if residue_name in STANDARD_AMINO_ACIDS:
            logging.info(f"Skipping standard amino acid residue: {pdb_file.name}")
            continue
        if residue_name in EXCLUDE_CRYSTAL_ADDITIVES:
            _log_malformed(pdb_file, f"excluded_crystal_additive:{residue_name}")
            logging.info(f"Skipping crystallization additive: {pdb_file.name}")
            continue

        with open(pdb_file, "r", encoding="utf-8", errors="ignore") as f:
            atom_lines = [line for line in f if line.startswith(("HETATM", "ATOM"))]
        if len(atom_lines) < MIN_ATOMS_FOR_DOCKING:
            _log_malformed(pdb_file, f"tiny_ligand_fewer_than_{MIN_ATOMS_FOR_DOCKING}_atoms")
            logging.info(f"Skipping tiny ligand: {pdb_file.name}")
            continue

        sanitized = pdb_file.with_suffix(".sanitized.pdb")
        ok = _sanitize_pdb(pdb_file, sanitized)
        if not ok:
            _log_malformed(pdb_file, "sanitize_kept_no_atoms")
            logging.warning(f"Sanitization produced no atoms: {pdb_file.name}")
            continue

        # run element/element-symbol fixer on the exact file we will convert
        try:
            pre_txt = sanitized.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            pre_txt = None
        try:
            # Column-accurate rewrite for HETATM/optional ATOM → element cols (77–78)
            fix_element_columns_in_file(str(sanitized), rewrite_atoms=False)
            # Secondary invariant: if any hydrogen-named atom still carries 'He', fix & log
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

        # quick buffer/counter-ion heuristics (unchanged)
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
                    residue_name in {"UNL", "LIG"} and _polyacidic_by_counts_from_pdbfile(sanitized)):
                flag_buffer, reason = True, "counts_only_polyacidic"

        if flag_buffer:
            _log_malformed(pdb_file, f"counterion_or_buffer:{reason or 'unknown'}")
            logging.info(f"Skipping likely counter-ion/buffer ({reason or 'unknown'}): {pdb_file.name}")
            continue

        try:
            if 'm_chk' in locals() and m_chk is not None and residue_name in {"UNL", "LIG"}:
                from rdkit.Chem import rdMolDescriptors as rdmd
                rings = rdmd.CalcNumRings(m_chk)
                arom = rdmd.CalcNumAromaticRings(m_chk)
                o = sum(1 for a in m_chk.GetAtoms() if a.GetSymbol() == "O")
                hac = m_chk.GetNumHeavyAtoms()
                if rings == 0 and arom == 0 and hac >= 12 and (o / float(hac)) >= 0.40:
                    _log_malformed(pdb_file, "UNL_O_rich_ringless")
                    logging.info(f"Skipping UNL O-rich ringless fragment: {pdb_file.name}")
                    continue
        except Exception:
            pass

        # write pristine SDF reference (unchanged behavior)
        ref_dir = sanitized.parent.parent / "reference"
        ref_dir.mkdir(parents=True, exist_ok=True)
        ref_sdf = ref_dir / (pdb_file.stem + ".sdf")
        try:
            m = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
            if m is not None:
                w = Chem.SDWriter(str(ref_sdf))
                try:
                    w.SetKekulize(False)
                except Exception:
                    pass
                w.write(m);
                w.close()
        except Exception:
            with open(ref_sdf, "w") as out:
                out.write(pdb_file.stem + "\n$$$$\n")

        pdbqt_path = prepped_ligands_dir / f"{pdb_file.stem}.pdbqt"

        if pdbqt_path.exists() and pdbqt_path.stat().st_size > 100 and is_valid_ligand(pdbqt_path,
                                                                                       log_dir=prepped_ligands_dir):
            logging.info(f"[resume] Valid PDBQT already exists, skipping: {pdbqt_path.name}")
            continue

        tmp_mol2 = sanitized.with_suffix(".mol2")
        mgl_ok = False
        ob_warn_no_H = False

        ob_cmd = [
            obabel_exe_short,
            "-ipdb", get_short_path_name(str(sanitized.resolve())),
            "-omol2", "-O", get_short_path_name(str(tmp_mol2.resolve()))
        ]
        logging.info("Converting sanitized PDB -> MOL2 (no gen3d): " + " ".join(map(str, ob_cmd)))
        try:
            res_ob = subprocess.run(ob_cmd, check=True, capture_output=True, text=True, timeout=300)
            # log stderr too (helps diagnose missing hydrogens, valence issues)
            if res_ob.stderr:
                logging.warning("[obabel primary stderr] %s", res_ob.stderr.strip())

            # hint for fallback if OBabel mentions hydrogen issue
            msg = ((res_ob.stdout or "") + (res_ob.stderr or "")).lower()
            if "no explicit hydrogens" in msg or "hydrogen" in msg and "warning" in msg:
                logging.info("[hint] OBabel reported hydrogen issue -> prefer fallback route for crystal control.")
                ob_warn_no_H = True

        except subprocess.CalledProcessError as e:
            logging.warning(f"OBabel PDB->MOL2 failed for {sanitized.name}:\n{e.stderr}")

        if tmp_mol2.exists() and tmp_mol2.stat().st_size > 100:
            # enforce RDKit aromaticity onto MOL2 so ADT sees stable A/NA types
            _re_aromatize_mol2_in_place(tmp_mol2, obabel_exe_short)

            prepare_cmd = [
                mgltools_python_short, str(prepare_script),
                "-l", get_short_path_name(str(tmp_mol2.resolve())),
                "-o", get_short_path_name(str(pdbqt_path.resolve())),
                "-U", "nphs_lps",
                "-A", "hydrogens",  # <-- fixed
            ]
            logging.info("Preparing ligand (MOL2): " + " ".join(map(str, prepare_cmd)))
            try:
                result = subprocess.run(
                    [str(x) for x in prepare_cmd],
                    check=True,
                    capture_output=True,
                    text=True,
                    cwd=str(sanitized.parent),
                    timeout=600
                )
                if result.stderr:
                    logging.warning("[mgltools primary stderr] %s", result.stderr.strip())
                logging.info(result.stdout)
                mgl_ok = True
                # --- helium post-write summary (primary writer) ---
                try:
                    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
                        lines = fh.readlines()
                    new_lines, fixes, q_reason = assert_no_helium_in_pdbqt(lines, pdb_file.stem)
                    # overwrite if we auto-fixed H-named atoms labeled as He
                    if fixes > 0:
                        with open(pdbqt_path, "w", encoding="utf-8") as out:
                            out.writelines(new_lines)
                    logging.info("[helium] ligand=%s fixes=%d quarantine=%s rules=%s",
                                 pdb_file.stem, fixes, str(bool(q_reason)), rules_version())
                    if q_reason:
                        # quarantine and record provenance in TSV
                        quarantine = prepped_ligands_dir / "quarantine"
                        quarantine.mkdir(exist_ok=True)
                        qpath = quarantine / pdbqt_path.name
                        try:
                            if pdbqt_path.exists():
                                pdbqt_path.replace(qpath)
                        except Exception:
                            pass
                        _append_prep_status(
                            status_log_path=prepped_ligands_dir / "ligand_prep_status.tsv",
                            ligand_name=pdb_file.stem,
                            status="quarantined",
                            reason="adt_helium_inconsistent",
                            relpath=str(qpath.name),
                            stage="validate_pdbqt",
                            failure_code="adt_helium_inconsistent",
                            failure_detail=q_reason,
                            fixes_count=0,
                            rules_ver=rules_version(),
                        )
                        # continue to next ligand
                        continue
                    else:
                        # normal success path still gets a TSV row with fixes_count (0 or >0)
                        _append_prep_status(
                            status_log_path=prepped_ligands_dir / "ligand_prep_status.tsv",
                            ligand_name=pdb_file.stem,
                            status="ok",
                            reason="",
                            relpath=str(pdbqt_path.name),
                            stage="write_pdbqt",
                            failure_code="",
                            failure_detail="",
                            fixes_count=fixes,
                            rules_ver=rules_version(),
                        )
                except Exception as e:
                    logging.warning("[helium] postwrite summary failed for %s: %s", pdb_file.stem, e)



            except subprocess.TimeoutExpired:
                logging.error(f"Timeout preparing {tmp_mol2.name}")
            except subprocess.CalledProcessError as e:
                logging.warning(f"MGLTools failed for {tmp_mol2.name} (will try fallback):\n{e.stderr}")

        needs_fallback = (
                ob_warn_no_H
                or not mgl_ok
                or not pdbqt_path.exists()
                or pdbqt_path.stat().st_size < 100
                or not is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir)
        )
        logging.info(
            "[prep decision] mgl_ok=%s exists=%s size=%sB valid=%s ob_warn_no_H=%s -> needs_fallback=%s | %s",
            mgl_ok,
            pdbqt_path.exists(),
            (pdbqt_path.stat().st_size if pdbqt_path.exists() else 0),
            is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir) if pdbqt_path.exists() else False,
            ob_warn_no_H,
            needs_fallback,
            sanitized.name,
        )

        if needs_fallback:
            logging.info(
                f"Primary prep failed/invalid for {sanitized.name}; trying RDKit→SDF (kekulize off)→OBabel→MOL2→MGLTools.")
            tmp_sdf = sanitized.with_suffix(".tmp.sdf")
            try:
                mol = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
                if mol is None:
                    raise RuntimeError("RDKit failed to read sanitized PDB")

                try:
                    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
                except Exception:
                    pass
                try:
                    mol = Chem.AddHs(mol, addCoords=True)
                except Exception:
                    pass

                if not _write_obabel_friendly_sdf(mol, tmp_sdf):
                    raise RuntimeError("Failed to write OBabel-friendly SDF in fallback")

                ob_cmd2 = [
                    obabel_exe_short,
                    "-isdf", get_short_path_name(str(tmp_sdf.resolve())),
                    "-omol2", "-O", get_short_path_name(str(tmp_mol2.resolve()))
                ]
                logging.info("Converting SDF -> MOL2 (no gen3d): " + " ".join(map(str, ob_cmd2)))
                try:
                    res_ob2 = subprocess.run(ob_cmd2, check=True, capture_output=True, text=True, timeout=300)
                    if res_ob2.stderr:
                        logging.warning("[obabel fallback stderr] %s", res_ob2.stderr.strip())
                    logging.info(res_ob2.stdout)
                    if "no explicit hydrogens" in ((res_ob2.stdout or "") + (res_ob2.stderr or "")).lower():
                        logging.info("[hint] Fallback OBabel also reports 'no explicit hydrogens'.")
                except subprocess.CalledProcessError as e:
                    logging.error(f"OBabel SDF->MOL2 fallback failed for {sanitized.name}:\n{e.stderr}")

                if tmp_mol2.exists() and tmp_mol2.stat().st_size > 100:
                    # --- ensure explicit hydrogens before fallback MGLTools ---
                    withH = tmp_mol2.with_name(tmp_mol2.stem + ".withH.mol2")
                    okH2, hstderr2 = add_hydrogens_mol2(tmp_mol2, withH, obabel_exe_short)
                    mol2_for_mgl2 = withH if okH2 else tmp_mol2
                    logging.info("[ligprep] addHs stage=fallback in=%s ok=%s", tmp_mol2.name, okH2)
                    if (hstderr2 or "").strip():
                        logging.warning("[ligprep] addHs stderr (fallback) %s", hstderr2.splitlines()[-1][:200])

                    # build the fallback prepare command using mol2_for_mgl2
                    prepare_cmd2 = [
                        mgltools_python_short, str(prepare_script),
                        "-l", get_short_path_name(str(mol2_for_mgl2.resolve())),
                        "-o", get_short_path_name(str(pdbqt_path.resolve())),
                        "-U", "nphs_lps",
                        "-A", "hydrogens",
                    ]

                    logging.info("Preparing ligand (fallback MOL2): " + " ".join(map(str, prepare_cmd2)))
                    try:
                        res2 = subprocess.run(
                            [str(x) for x in prepare_cmd2],
                            check=True,
                            capture_output=True,
                            text=True,
                            timeout=600
                        )
                        # --- helium post-write summary (fallback writer) ---
                        try:
                            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
                                lines = fh.readlines()
                            new_lines, fixes, q_reason = assert_no_helium_in_pdbqt(lines, pdb_file.stem)
                            if fixes > 0:
                                with open(pdbqt_path, "w", encoding="utf-8") as out:
                                    out.writelines(new_lines)
                            logging.info("[helium] ligand=%s fixes=%d quarantine=%s rules=%s",
                                         pdb_file.stem, fixes, str(bool(q_reason)), rules_version())
                            if q_reason:
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
                                    ligand_name=pdb_file.stem,
                                    status="quarantined",
                                    reason="adt_helium_inconsistent",
                                    relpath=str(qpath.name),
                                    stage="validate_pdbqt",
                                    failure_code="adt_helium_inconsistent",
                                    failure_detail=q_reason,
                                    fixes_count=0,
                                    rules_ver=rules_version(),
                                )
                                return (mol2_file.name, "postcheck_fail")
                            else:
                                _append_prep_status(
                                    status_log_path=status_log_dir / "ligand_prep_status.tsv",
                                    ligand_name=mol2_file.name,
                                    status="ok",
                                    reason="",
                                    relpath=str(pdbqt_path.name),
                                    stage="write_pdbqt",
                                    failure_code="",
                                    failure_detail="",
                                    fixes_count=fixes,
                                    rules_ver=rules_version(),
                                )
                        except Exception as e:
                            logging.warning("[helium] fallback postwrite summary failed for %s: %s", mol2_file.name, e)

                        if res2.stderr:
                            logging.warning("[mgltools fallback stderr] %s", res2.stderr.strip())

                        logging.info(res2.stdout)
                    except subprocess.CalledProcessError as e:
                        logging.error(f"Fallback MGLTools failed for {sanitized.name}:\n{e.stderr}")
                    except subprocess.TimeoutExpired:
                        logging.error(f"Fallback timeout for {sanitized.name}")

            finally:
                try:
                    tmp_sdf.unlink(missing_ok=True)
                except Exception:
                    pass
        # --- extra postcheck: ensure no *protein residues* slipped into ligand PDBQT
        try:
            has_protein_res = False
            with open(pdbqt_path, "r", errors="ignore") as fh:
                for ln in fh:
                    if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                        continue
                    # PDB/PDBQT fixed columns: residue name is cols 18-20 (0-based 17:20)
                    resn = (ln[17:20] if len(ln) >= 20 else "").strip().upper()
                    if resn in STANDARD_AMINO_ACIDS:
                        has_protein_res = True
                        break
            if has_protein_res:
                _log_malformed(pdb_file, "protein_residue_in_pdbqt")
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
            pass
            continue
        except Exception:
            # don't crash prep if the check itself fails; log and proceed to other guards
            logging.warning(f"[postcheck] unable to scan for protein ATOM in {pdbqt_path.name}")

        try:
            if pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
                tmp_mol2.unlink(missing_ok=True)
        except Exception:
            pass

        # final validation & quarantine with correct log_dir
        if not pdbqt_path.exists() or pdbqt_path.stat().st_size < 100 or not is_valid_ligand(pdbqt_path,
                                                                                             log_dir=prepped_ligands_dir):
            _log_malformed(pdb_file, "pdbqt_postcheck_fail_or_small")
            quarantine = prepped_ligands_dir / QUARANTINE_DIRNAME
            quarantine.mkdir(exist_ok=True)
            try:
                if pdbqt_path.exists():
                    pdbqt_path.replace(quarantine / pdbqt_path.name)
            except Exception:
                pass
            _append_prep_status(status_log, sanitized.name, "FAIL", "postcheck_fail_or_small",
                                str((quarantine / pdbqt_path.name).relative_to(prepped_ligands_dir)) if (
                                            quarantine / pdbqt_path.name).exists() else "")
            continue

        logging.info(f"Created PDBQT: {pdbqt_path.name}")
        _append_prep_status(status_log, sanitized.name, "OK", "",
                            str(pdbqt_path.relative_to(prepped_ligands_dir)))


# =========================
# Main SDF ? MOL2 ? PDBQT pipeline
# =========================

def _valid_pdbqt(path: Path, log_dir: Path) -> bool:
    return path.exists() and path.stat().st_size > 100 and is_valid_ligand(path, log_dir=log_dir)


def prep_ligands_with_mgltools():
    print("Starting ligand preparation")

    cfg = read_config()

    ligand_extracted_dir = Path(cfg["LIGAND_EXTRACTED_DIR"]).resolve()
    ligands_mol2_dir = Path(cfg["LIGANDS_MOL2_DIR"]).resolve()
    output_ligands_dir = Path(cfg["OUTPUT_LIGANDS_DIR"]).resolve()
    prepped_ligands_dir = output_ligands_dir
    output_ligands_dir.mkdir(parents=True, exist_ok=True)

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
    status_log = output_ligands_dir / "ligand_prep_status.tsv"

    sdf_files = list(ligand_extracted_dir.glob("*.sdf"))
    print(f"Found {len(sdf_files)} SDF file(s)")
    if not sdf_files:
        return

    for sdf_file in sdf_files:
        print(f"\n=== Processing SDF: {sdf_file.name} ===")
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
            print("No MOL2 files produced; skipping this SDF.")
            continue

        print(f"Preparing {len(mol2_files)} MOL2 files with MGLTools (parallel)")
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
                pdbqt_path = output_ligands_dir / f"{mol2_file.stem}.pdbqt"
                # cache check: skip if a valid PDBQT already exists
                if pdbqt_path.exists() and pdbqt_path.stat().st_size > 100 and is_valid_ligand(pdbqt_path,
                                                                                               log_dir=output_ligands_dir):
                    logging.info(f"[resume] Valid PDBQT already exists, skipping: {pdbqt_path.name}")
                    continue

                futures.append(ex.submit(
                    _prepare_one,
                    mgltools_python_short, prepare_script_short, mol2_file, pdbqt_path, obabel_exe_short,
                    status_log_dir=prepped_ligands_dir
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

        # CLI helpers so you can repair a whole tree in one command:
        if __name__ == "__main__":
            import argparse
            ap = argparse.ArgumentParser(description="Atlas path/layout doctor")
            ap.add_argument("--migrate-legacy", metavar="PROCESSED_ROOT",
                            help="Scan processed_pdbs and migrate *_NOLIG/*_CLEANED_LIGANDS into canonical layout.")
            ap.add_argument("--expose-intermediates", nargs=2, metavar=("PROCESSED_ROOT", "PREPPED_LIGANDS_ROOT"),
                            help="Create/refresh prepped_ligands/<PDB>/intermediates symlinks for all PDBs.")
            args, _ = ap.parse_known_args()

            if args.migrate_legacy:
                root = _as_path(args.migrate_legacy).resolve()
                for d in sorted(root.iterdir()):
                    if not d.is_dir():
                        continue
                    m = re.match(r"^([A-Za-z0-9]{4})(?:_.+)?$", d.name)
                    if not m:
                        continue
                    pdb_id = m.group(1).upper()
                    try:
                        fold_legacy_layout(pdb_id, root)
                    except Exception as e:
                        logging.warning("migrate-legacy skip %s: %s", pdb_id, e)

            if args.expose_intermediates:
                proc_root = _as_path(args.expose_intermediates[0]).resolve()
                prepped_root = _as_path(args.expose_intermediates[1]).resolve()
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

        print(f"Done: {ok_count} ok, {fail_count} failed/quarantined for {sdf_file.name}")


if __name__ == "__main__":
    prep_ligands_with_mgltools()
