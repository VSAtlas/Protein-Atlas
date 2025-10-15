import os
import csv
import subprocess
import shutil
import logging
from installation import load_config
from collections import defaultdict
from logger_setup import setup_logger
from pathlib import Path
import hashlib
from typing import Iterable, Tuple, List


# Load config
config = load_config()
P2RANK_DIR = config.get("P2RANK_PATH")
import yaml
from types import SimpleNamespace

ALIASES_PATH = (
    config.get("ALIASES_PATH")
    or os.environ.get("ALIASES_YAML")
    or os.path.join(os.path.dirname(__file__), "aliases.yaml")
)

# if not found, try ./chemdb/aliases.yaml automatically
if not os.path.exists(ALIASES_PATH):
    probe = os.path.join(os.path.dirname(__file__), "chemdb", "aliases.yaml")
    if os.path.exists(probe):
        ALIASES_PATH = probe

_aliases_cache = None
_rules_cache = None

# ---  YAML loader with encoding fallbacks & punctuation cleanup ---
def _load_aliases_yaml():
    import yaml, unicodedata
    global _aliases_cache
    if _aliases_cache is not None:
        return _aliases_cache

    path = ALIASES_PATH
    try:
        # try UTF-8 first (fast path)
        with open(path, "r", encoding="utf-8") as fh:
            _aliases_cache = yaml.safe_load(fh) or {}
            return _aliases_cache
    except Exception as e_utf8:
        logging.warning("[aliases] UTF-8 load failed for %s: %s", path, e_utf8)

    try:
        # fallback: cp1252 decode + punctuation normalization → UTF-8
        raw = Path(path).read_bytes()
        text = raw.decode("cp1252", errors="strict")

        # normalize common Windows punctuation to ASCII
        repl = {
            "\u2018": "'", "\u2019": "'", "\u201c": '"', "\u201d": '"',
            "\u2013": "-", "\u2014": "-", "\u2026": "...",
            "\u00A0": " ", "\u200B": "", "\uFEFF": "",
        }
        for k, v in repl.items():
            text = text.replace(k, v)

        # fix common mojibake 'â†’' if present
        text = text.replace("â†’", "->")

        text = text.replace("\r\n", "\n").replace("\r", "\n")
        _aliases_cache = yaml.safe_load(text) or {}
        logging.info("[aliases] loaded with cp1252 fallback and sanitized punctuation")
        return _aliases_cache
    except Exception as e_cp:
        logging.warning("[aliases] cp1252 fallback failed for %s: %s", path, e_cp)
        _aliases_cache = {}
        return _aliases_cache



def get_atom_rules():
    """
    Return a SimpleNamespace of normalized sets/maps used by element fixing,
    plus compatibility views so legacy code can still do RULES["element_sets"].
    """
    global _rules_cache
    if _rules_cache is not None:
        return _rules_cache

    a  = _load_aliases_yaml() or {}
    es = a.get("element_sets", {}) or {}
    ligand_sets = a.get("ligand_sets", {}) or {}
    meeko_cfg   = a.get("meeko", {}) or {}

    # helpers
    def _flatten(items):
        out = []
        for x in (items or []):
            if isinstance(x, (list, tuple, set)):
                out.extend(_flatten(x))
            else:
                parts = [p.strip() for p in str(x).split(";")]
                out.extend([p for p in parts if p])
        return out

    def _as_set(items, up=True):
        return { (v.upper() if up else v) for v in _flatten(items) }

    # element/name logic
    peptide_like      = _as_set(es.get("peptide_like_names"))
    one_letter        = _as_set(es.get("one_letter_elements"))
    two_letter        = _as_set(es.get("two_letter_elements"))
    halide_resnames   = _as_set(es.get("halide_resnames"))
    default_element   = (es.get("default_element") or "C").upper()
    treat_backbone_ca = bool(es.get("treat_backbone_CA_as_C", True))

    # AD4 types (support both top-level and nested)
    ad4_types_yaml = a.get("ad4_types") or (a.get("pdbqt_types") or {}).get("ad4_types")
    ad4_types = _as_set(ad4_types_yaml) if ad4_types_yaml else set()

    # NEW: ligand & Meeko lists from YAML
    nucleotide_like_resnames = _as_set(ligand_sets.get("nucleotide_like_resnames"))
    meeko_drop_free_ions     = _as_set(meeko_cfg.get("drop_free_ions"))

    # retain list (and a compat copy)
    retain_raw = a.get("retain_in_receptor_resnames", []) or []
    retain_res = _as_set(retain_raw)
    compat_element_sets = dict(es)      # keep YAML shape
    compat_retain_list  = list(retain_raw)

    from types import SimpleNamespace
    _rules_cache = SimpleNamespace(
        # normalized sets / maps
        peptide_like=peptide_like,
        one_letter=one_letter,
        two_letter=two_letter,
        halide_resnames=halide_resnames,
        default_element=default_element,
        treat_backbone_ca=treat_backbone_ca,
        retain_resnames=retain_res,

        # NEW exports used elsewhere
        nucleotide_like_resnames=sorted(nucleotide_like_resnames),
        meeko_drop_free_ions=sorted(meeko_drop_free_ions),

        # name/alias maps (uppercased keys/values)
        prefix_map={ (k or "").upper(): (v or "").upper()
                     for k, v in (es.get("derive_prefix_map") or {}).items() },
        special_names={ (k or "").upper(): (v or "").upper()
                        for k, v in (es.get("special_atom_names") or {}).items() },
        halide_aliases={ (k or "").upper(): (v or "").upper()
                         for k, v in (es.get("halide_resname_aliases") or {}).items() },
        cation_aliases={ (k or "").upper(): (v or "").upper()
                         for k, v in (es.get("cation_resname_aliases") or {}).items() },

        # compatibility views for older call sites
        element_sets=compat_element_sets,
        retain_in_receptor_resnames=compat_retain_list,
        ad4_types=ad4_types,
    )
    return _rules_cache





def scan_helium_counts(text_or_lines) -> int:
    """
    Count occurrences that look like element/ADT 'He' in PDB or PDBQT.
    - PDB: fixed-width columns 77–78
    - PDBQT: trailing ADT token
    """
    n = 0
    if isinstance(text_or_lines, str):
        lines = text_or_lines.splitlines()
    else:
        lines = list(text_or_lines)
    for ln in lines:
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        # PDB fixed-width first
        if len(ln) >= 78 and ln[76:78].strip().upper() == "HE":
            n += 1
            continue
        # PDBQT trailing token
        parts = ln.split()
        if parts and parts[-1].strip().upper() == "HE":
            n += 1
    return n

def rules_version() -> str:
    """
    Return a short identifier for the current YAML/alias ruleset if available,
    else a static placeholder. Keep this stable for log grepping.
    """
    try:
        r = get_atom_rules()
        for k in ("version", "rules_version", "hash", "source_hash"):
            v = getattr(r, k, None) if hasattr(r, k) else r.get(k) if isinstance(r, dict) else None
            if v:
                return str(v)
    except Exception:
        pass
    return "rules:v1"


def format_pdb_atom_debug(line: str, line_no: int | None = None) -> str:
    """
    Render a fixed-width debug string exposing key PDB fields + element cols 77–78.
    """
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return line.strip()
    resn = line[17:20].strip()
    chain = line[21].strip() or "_"
    resi = line[22:26].strip()
    icode = line[26].strip() or " "
    name = line[12:16].strip()
    alt = line[16] if len(line) > 16 else " "
    elem = (line[76:78] if len(line) >= 78 else "  ").strip() or "?"
    ln = f"line{line_no}" if line_no is not None else "line?"
    return f"{ln} {chain}:{resi}{icode}:{resn} {name} alt={alt} elem={elem}  cols77-78='{(line[76:78] if len(line)>=78 else '  ')}'"


def scan_helium_counts_with_hits(text: str, max_hits: int = 5) -> tuple[int, list[str]]:
    """
    Wrapper around scan_helium_counts that also returns up to `max_hits` formatted offenders.
    """
    n = scan_helium_counts(text)
    if n == 0:
        return 0, []
    hits = []
    for i, ln in enumerate(text.splitlines(), start=1):
        if ln.startswith(("ATOM", "HETATM")) and (ln[76:78].strip().upper() == "HE"):
            hits.append(format_pdb_atom_debug(ln, i))
            if len(hits) >= max_hits:
                break
    return n, hits


def assert_no_helium_in_hydrogen_names(pdb_text: str) -> Tuple[str, int]:
    """
    If an ATOM/HETATM line has an atom *name* beginning with 'H' but element/ADT is 'He',
    rewrite the element to 'H'. Returns (new_text, n_fixes).
    """
    out: List[str] = []
    fixes = 0
    for ln in pdb_text.splitlines(True):
        if ln.startswith(("ATOM  ","HETATM")) and len(ln) >= 78:
            aname = ln[12:16].strip().upper()
            # ADT usually in the last whitespace token for your PDBQT; for PDB use columns 77-78
            adt_or_elem = ln.split()[-1].strip().upper()
            if aname.startswith("H") and adt_or_elem == "HE":
                # fix element columns if fixed-width PDB, else re-map the trailing ADT token
                # Prefer fixed-width column 77-78 if present:
                if ln[76:78].strip().upper() in {"HE", "H", ""}:
                    ln = ln[:76] + f"{'H':>2}" + ln[78:]
                else:
                    # trailing token path: reassemble line by replacing the last token
                    parts = ln.rstrip("\n").split()
                    parts[-1] = "H"
                    ln = " ".join(parts) + ("\n" if ln.endswith("\n") else "")
                fixes += 1
        out.append(ln)
    return "".join(out), fixes

def assert_no_helium_in_pdbqt(lines: Iterable[str], ligand_name: str) -> Tuple[List[str], int, str]:
    """
    Post-write guardrail for PDBQT:
      - If ADT='He' appears where atom NAME begins with 'H' → auto-correct to 'H', log fixes.
      - If ADT='He' appears on a non-H-named atom → quarantine is required; return reason.
    Returns (new_lines, fixes_count, quarantine_reason_or_empty).
    """
    new_lines: List[str] = []
    fixes = 0
    quarantine_reason = ""
    for ln in lines:
        if ln.startswith(("ATOM", "HETATM")):
            parts = ln.split()
            if parts:
                adt = parts[-1].strip().upper()
                if adt == "HE":
                    aname = ln[12:16].strip().upper() if len(ln) >= 16 else ""
                    if aname.startswith("H"):
                        parts[-1] = "H"
                        ln = " ".join(parts) + ("\n" if not ln.endswith("\n") else "")
                        fixes += 1
                    else:
                        quarantine_reason = f"adt_helium_inconsistent: non-H name '{aname}' has ADT=He (ligand={ligand_name})"
        new_lines.append(ln)
    return new_lines, fixes, quarantine_reason



def derive_element(aname: str, resname: str, is_het: bool, rules=None) -> str:
    """Infer element symbol from atom/residue context using YAML-driven rules."""
    if rules is None:
        rules = get_atom_rules()

    an = (aname or "").strip().upper()
    rn = (resname or "").strip().upper()

    # special atom names (e.g., OXT)
    if an in rules.special_names:
        logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                      an, rn, is_het, "special_names", rules.special_names.get(an, "?"))
        return rules.special_names[an]

    # normalize residue-name aliases (e.g., IOD->I, CL- -> CL)
    if rn in rules.halide_aliases:
        logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                      an, rn, is_het, "halide_alias", rules.halide_aliases[rn])
        rn = rules.halide_aliases[rn]
    if rn in rules.cation_aliases:
        logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                      an, rn, is_het, "cation_alias", rules.cation_aliases[rn])
        rn = rules.cation_aliases[rn]

    # derive by leading functional prefix (OE1, NE2, OD1, ND2, SD, ...)
    pref = an[:2]
    if pref in rules.prefix_map:
        logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                      an, rn, is_het, "prefix_map", rules.prefix_map.get(pref, "?"))
        return rules.prefix_map[pref]
    # halide ions by residue name for single-atom HETATMs
    if is_het and rn in rules.halide_resnames:
        logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                      an, rn, is_het, "halide_resname", rn)
        if rn in {"CL","BR"}:
            return rn[0] + rn[1].lower()
        return rn  # I, F
    
    
    
    # CA special-case (avoid backbone CA -> Calcium)
    if len(an) >= 2 and an[:2] == "CA":
        if not is_het and rules.treat_backbone_ca:
            logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                          an, rn, is_het, "backbone_CA_guard", "C")
            return "C"
        if is_het and rn in {"CA","CAL"}:
            logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                          an, rn, is_het, "Ca_ion", "Ca")
            return "Ca"

    # two-letter elements at name start (Cl, Br, Na, Mg, ...) — HETs only, and only when the
    # atom name itself looks like a stand-alone element token (length==2 or 3rd char not alpha).
    if is_het:
        two = an[:2].upper()
        looks_like_standalone = (len(an) == 2) or (len(an) >= 3 and not an[2].isalpha())
        if two in rules.two_letter and looks_like_standalone:
            t = two
            logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                          an, rn, is_het, "two_letter", t[0] + t[1].lower())
            return t[0] + t[1].lower()

        # --- Guard: avoid mislabeling organic HET atom names like "CAK","NAA" as Ca/Na ---
        # If this is a HET but NOT a known simple ion residue, and the atom name
        # continues with another alpha character (e.g., "CAK", "NAA"), prefer a
        # one-letter element guess (C/N/...) rather than a two-letter metal.
    if is_het:
        _ion_res = {
            "LI", "NA", "K", "RB", "CS",
            "MG", "CA", "SR", "BA",
            "ZN", "CU", "NI", "CO", "FE", "MN", "CD", "AL", "HG", "AG", "PB",
            "PT", "PD", "AU", "RU", "IR", "OS"
        }
        if rn not in _ion_res:
            if len(an) >= 3 and an[0].isalpha() and an[1].isalpha() and an[2].isalpha():
                c = an[0].upper()
                if c in rules.one_letter:
                    logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                                  an, rn, is_het, "het_three_letters_guard", c)
                    return c
    # hydrogens (H, 1H, 2H...)
    if an.startswith("H") or (an[:1].isdigit() and len(an) >= 2 and an[1] == "H"):
        logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                      an, rn, is_het, "hydrogen_name", "H")
        return "H"

    # one-letter defaults by first alpha
    if an and an[0].isalpha():
        c = an[0].upper()
        if c in rules.one_letter:
            logging.debug("[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                          an, rn, is_het, "one_letter", c)
            return c

    for ch in an:
        if ch.isalpha():
            c = ch.upper()
            return c if c in rules.one_letter else rules.default_element

    return rules.default_element




def ensure_model_records(pdb_input_path: str, pdb_output_path: str):
    with open(pdb_input_path, 'r') as f:
        lines = f.readlines()

    has_model = any(line.startswith("MODEL") for line in lines)

    if has_model:
        with open(pdb_output_path, 'w') as f:
            f.writelines(lines)
        logging.info(f"MODEL record found in {pdb_input_path}. File copied unchanged.")
        return

    atom_start_idx = None
    atom_end_idx = None

    for i, line in enumerate(lines):
        if line.startswith(("ATOM  ", "HETATM")):
            atom_start_idx = i
            break

    for i in range(len(lines) - 1, -1, -1):
        if lines[i].startswith(("ATOM  ", "HETATM")):
            atom_end_idx = i
            break

    if atom_start_idx is None or atom_end_idx is None:
        raise ValueError(f"No ATOM or HETATM lines found in {pdb_input_path}.")

    lines.insert(atom_start_idx, "MODEL        1\n")
    lines.insert(atom_end_idx + 2, "ENDMDL\n")

    with open(pdb_output_path, 'w') as f:
        f.writelines(lines)

    logging.info(f"MODEL/ENDMDL added in {pdb_output_path} between lines {atom_start_idx+1} and {atom_end_idx+3}.")


def remove_unparsable_hetatms(pdb_path):
    cleaned_lines = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                if 'UNK' in atom_name or 'UNX' in line:
                    continue
            cleaned_lines.append(line)
    with open(pdb_path, 'w') as f:
        f.writelines(cleaned_lines)
    logging.info(f"Unparsable HETATM entries removed from {pdb_path}.")


from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
import os

class ElementFixer(Select):
    """Rewrite element symbols safely, with peptide-like HET awareness (YAML-driven)."""
    def __init__(self):
        super().__init__()
        self.rules = get_atom_rules()

    def _is_peptidic_like(self, residue) -> bool:
        try:
            names = [a.get_name().strip().upper() for a in residue.get_atoms()]
        except Exception:
            return False
        if not names:
            return False
        hits = sum((n in self.rules.peptide_like) or (n[:2] in {"OE","NE","OD","ND","SD"}) for n in names)
        return hits >= max(4, int(0.6 * len(names)))

    def get_atom_element(self, atom):
        aname = atom.get_name()
        res   = atom.get_parent()
        rname = getattr(res, "get_resname", lambda: "")()

        try:
            hetflag = res.get_id()[0] if res is not None else " "
        except Exception:
            hetflag = " "
        is_het = (hetflag != " ")

        # if residue looks peptide-like, still derive by name (prevents ion mislabels)
        if self._is_peptidic_like(res):
            return derive_element(aname, rname, is_het, self.rules)

        return derive_element(aname, rname, is_het, self.rules)

    def accept_atom(self, atom):
        atom.element = self.get_atom_element(atom)
        return True



def fix_pdb_elements(input_path, output_path=None):
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    structure = parser.get_structure("structure", input_path)
    io.set_structure(structure)

    fixer = ElementFixer()

    if output_path is None:
        output_path = input_path

    with open(output_path, "w") as fh:
        io.save(fh, fixer)
    print(f"Saved fixed PDB to {output_path}")

from Bio.PDB import PDBParser, NeighborSearch, Selection

def calculate_ligand_protein_contacts(protein_pdb_path, ligand_lines, distance_cutoff=4.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', protein_pdb_path)
    model = structure[0]

    # Collect all protein atoms
    protein_atoms = [atom for atom in model.get_atoms() if atom.get_parent().get_id()[0] == ' ']

    # Parse ligand atoms from lines
    ligand_atoms = []
    for line in ligand_lines:
        if line.startswith(('HETATM', 'ATOM')):
            # Simple PDB atom line parsing for coords
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            # Create a dummy atom-like object
            class DummyAtom:
                def __init__(self, coord): self.coord = coord
                def get_coord(self): return self.coord
            ligand_atoms.append(DummyAtom((x,y,z)))

    # Use NeighborSearch to find protein atoms near ligand atoms
    ns = NeighborSearch(protein_atoms)
    contact_count = 0
    for latom in ligand_atoms:
        close_atoms = ns.search(latom.get_coord(), distance_cutoff)
        contact_count += len(close_atoms)

    return contact_count



def rank_ligands_by_atom_count(ligands_dict):
    return sorted(ligands_dict.items(), key=lambda item: len(item[1]), reverse=True)

def _fix_ligand_element_columns_in_memory(lines):
    """Rewrite element cols (77–78) for HET ligands using YAML rules."""
    rules = get_atom_rules()
    out = []
    for line in lines:
        if line.startswith(("ATOM  ","HETATM")) and len(line) >= 78:
            aname = line[12:16]
            resn  = line[17:20]
            is_het = line.startswith("HETATM")
            el = derive_element(aname, resn, is_het, rules)
            line = line[:76] + f"{el:>2}" + line[78:]
        out.append(line)
    return out






def fix_element_columns_in_file(src_path, dst_path=None, rewrite_atoms=False):
    rules = get_atom_rules()
    src_path = str(src_path)
    dst_path = src_path if dst_path is None else str(dst_path)
    logging.info("[elemfix] in=%s rewrite_atoms=%s", src_path, rewrite_atoms)
    out_lines = []
    with open(src_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            is_atom = line.startswith("ATOM  ")
            is_het  = line.startswith("HETATM")
            if len(line) >= 78 and (is_het or (rewrite_atoms and is_atom)):
                aname = line[12:16]
                resn  = line[17:20]
                el = derive_element(aname, resn, is_het, rules)
                line = line[:76] + f"{el:>2}" + line[78:]
            out_lines.append(line)
    with open(dst_path, "w", encoding="utf-8") as out:
        out.writelines(out_lines)
    logging.info("[elemfix] done=%s", dst_path or src_path)
    return dst_path




def extract_and_remove_ligands(pdb_path, output_cleaned_pdb, ligands_dir):
    os.makedirs(ligands_dir, exist_ok=True)
    ligands = defaultdict(list)
    ligand_coords = []  # collect all ligand atom coords for box calculation
    rules = get_atom_rules()
    retained_resnames = rules.retain_resnames  # already uppercased

    # Do NOT retain common cryos/buffers: GOL/EDO/PG4/MPD/ACT/TRS/PO4/PEG → they’ll be extracted

    with open(pdb_path, 'r') as infile, open(output_cleaned_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith('HETATM'):
                resname_raw = line[17:20]
                resname = resname_raw.strip().upper()  # normalize for set membership & filenames
                chain = line[21]
                resnum = line[22:26].strip()
                if resname not in retained_resnames:
                    ligands[(chain, resname, resnum)].append(line)
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_coords.append((x, y, z))
                    except ValueError:
                        logging.warning(f"Invalid ligand coordinates in line: {line.strip()}")
                    continue
            outfile.write(line)

    # Write separate ligand files (with element repair)
    for (chain, resname, resnum), lines in ligands.items():
        ligand_fname = os.path.join(ligands_dir, f"{resname}_{chain}{resnum}.pdb")
        fixed = _fix_ligand_element_columns_in_memory(lines)  # <<< ADD
        with open(ligand_fname, 'w') as lf:
            lf.writelines(fixed)
        logging.info(f"Saved ligand {resname} {chain}{resnum} to {ligand_fname} (elements repaired)")


    logging.info(f"Ligands extracted and removed from {pdb_path}.")
    return ligands, ligand_coords

def compute_box_from_ligand_coords(coords):
    if not coords:
        return None, None

    x_vals, y_vals, z_vals = zip(*coords)
    center = (sum(x_vals)/len(x_vals), sum(y_vals)/len(y_vals), sum(z_vals)/len(z_vals))

    buffer = 5.0
    x_range = max(x_vals) - min(x_vals) + buffer
    y_range = max(y_vals) - min(y_vals) + buffer
    z_range = max(z_vals) - min(z_vals) + buffer

    box_size = (x_range, y_range, z_range)
    return center, box_size


def get_box_from_p2rank_csv(pdb_file):
    pdb_file = os.path.abspath(pdb_file)
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    pred_dir = os.path.join(P2RANK_DIR, "test_output", f"predict_{pdb_name}")
    pred_file = os.path.join(pred_dir, f"{pdb_name}.pdb_predictions.csv")
    
    logging.info(f"Running P2Rank for: {pdb_file}")

    # Resolve P2Rank execution: accept either a dir, a jar path, or rely on PATH ("prank")
    pr_cfg = P2RANK_DIR or ""
    pr_base = None
    pr_exe = None
    pr_jar = None

    # If config points to a file and ends with .jar, treat it as the jar
    if pr_cfg and os.path.isfile(pr_cfg) and pr_cfg.lower().endswith(".jar"):
        pr_jar = pr_cfg
        pr_base = os.path.dirname(os.path.dirname(pr_cfg))  # .../bin/p2rank.jar -> parent of bin
    elif pr_cfg and os.path.isdir(pr_cfg):
        # config is an install root (directory)
        jar_candidate = os.path.join(pr_cfg, "bin", "p2rank.jar")
        if os.path.isfile(jar_candidate):
            pr_jar = jar_candidate
        pr_base = pr_cfg

    # If prank is on PATH, prefer it (simplest + portable)
    from shutil import which
    pr_on_path = which("prank")

    # Pick command: (1) prank on PATH, else (2) java -jar p2rank.jar, else error
    if pr_on_path:
        run_cmd = ["prank", "predict"]
    elif pr_jar and os.path.isfile(pr_jar):
        run_cmd = ["java", "-Xmx4G", "-jar", pr_jar, "predict"]
    else:
        logging.error("P2Rank not found. Set P2RANK_PATH to either the install directory or the p2rank.jar.")
        return None, None

    # Always write to a known output folder next to the input PDB
    out_dir = os.path.join(os.path.dirname(os.path.abspath(pdb_file)), "_p2rank")
    os.makedirs(out_dir, exist_ok=True)

    try:
        # Explicit -o ensures we know exactly where predictions land
        cmd = run_cmd + ["-f", pdb_file, "-o", out_dir]
        logging.info("P2Rank cmd: %s", " ".join(cmd))
        subprocess.run(cmd, check=True, shell=False)
        logging.info(f"P2Rank ran successfully for {pdb_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"P2Rank failed: {e}")
        return None, None

    # Prefer the explicit output location; fall back to legacy locations if needed
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    pred_file = os.path.join(out_dir, f"{pdb_name}.pdb_predictions.csv")
    if not os.path.isfile(pred_file):
        # legacy fallback: current working dir default
        legacy = os.path.join("test_output", f"predict_{pdb_name}", f"{pdb_name}.pdb_predictions.csv")
        alt = os.path.join(pr_base or "", "test_output", f"predict_{pdb_name}", f"{pdb_name}.pdb_predictions.csv")
        for probe in (legacy, alt):
            if probe and os.path.isfile(probe):
                pred_file = probe
                break

    if not os.path.isfile(pred_file):
        logging.warning(f"Prediction file not created: {pred_file}")
        return None, None
    with open(pred_file, "r", newline='') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [field.strip() for field in reader.fieldnames]

        for i, row in enumerate(reader):
            print(f"Row {i}: {row}")
            break  # for debugging

        # Rewind the file and parse again to get the top pocket
        f.seek(0)
        reader = csv.DictReader(f)
        reader.fieldnames = [field.strip() for field in reader.fieldnames]

        top_pocket = None
        for row in reader:
            try:
                rank_val = int(row.get("rank", "").strip())
                if rank_val == 1:
                    top_pocket = row
                    break
            except ValueError:
                logging.warning(f"Could not convert rank to int: {row.get('rank')}")

    if not top_pocket:
        logging.warning("No pocket found with rank 1")
        return None, None

    try:
        center = (
            float(top_pocket.get("center_x", "").strip()),
            float(top_pocket.get("center_y", "").strip()),
            float(top_pocket.get("center_z", "").strip()),
        )
        box_dim = float(top_pocket.get("surf_atoms", "20.0").strip())
        MAX_BOX_SIZE = 40.0
        box_size = (box_dim, box_dim, box_dim)
        clamped_box_size = tuple(min(dim, MAX_BOX_SIZE) for dim in box_size)
        print(f"[DEBUG] Returning center={center}, box_size={clamped_box_size}")
        return center, clamped_box_size
    except (KeyError, ValueError) as e:
        logging.error(f"Error parsing P2Rank pocket fields: {e}")
        return None, None


def detect_pocket(cleaned_pdb, logger):
    """
    Detect pocket center and box size from cleaned PDB.
    Returns (center, box_size) or (None, None) if failed.
    """
    center, box_size = get_box_from_p2rank_csv(cleaned_pdb)
    if center is None:
        logger.warning("Active-site detection failed.")
    return center, box_size

def prepare_receptor(cfg, paths, logger):
    """
    Run or reuse protein preparation to produce:
      - cleaned PDB without ligands
      - receptor PDBQT
    Returns (cleaned_pdb_path_str, receptor_pdbqt_path_str) or (None, None) on failure.
    """
    from distutils.util import strtobool
    force_reprocess = bool(strtobool(str(cfg.get("FORCE_REPROCESS", False))))
    logger.info(f"FORCE_REPROCESS={force_reprocess} | "
                f"exists(cleaned)={paths['cleaned_pdb_path'].exists()} "
                f"exists(receptor)={paths['receptor_pdbqt_path'].exists()}")

    if paths["cleaned_pdb_path"].exists() and paths["receptor_pdbqt_path"].exists() and not force_reprocess:
        logger.info("Reusing existing cleaned PDB and receptor PDBQT.")
        return norm(paths["cleaned_pdb_path"]), norm(paths["receptor_pdbqt_path"])
    import automate_protein_prep
    result = automate_protein_prep.main(str(paths["nolig_pdb_path"]))
    if not result or not isinstance(result, tuple) or len(result) != 2:
        logger.warning("Protein prep failed.")
        return None, None

    cleaned_pdb, receptor_pdbqt = result

    # Ensure receptor lives in canonical PDBQT_DIR
    try:
        if Path(receptor_pdbqt).resolve() != paths["receptor_pdbqt_path"].resolve():
            from shutil import copy2
            paths["receptor_pdbqt_path"].parent.mkdir(parents=True, exist_ok=True)
            copy2(receptor_pdbqt, paths["receptor_pdbqt_path"])
            receptor_pdbqt = str(paths["receptor_pdbqt_path"])
    except Exception as e:
        logger.warning(f"Could not relocate receptor PDBQT: {e}")

    return norm(cleaned_pdb), norm(receptor_pdbqt)


# ---------- high-level pipeline steps ----------

def extract_ligands(cfg, paths, logger):
    """
    Extract and strip ligands from input PDB into clean PDB without ligands.
    Returns the ligand count.
    """
    malformed_log = paths["ligands_mol2_dir"] / "malformed_ligands.txt"
    if malformed_log.exists():
        malformed_log.unlink()

    ligands_dict, _ = extract_and_remove_ligands(
        paths["pdb_path"], paths["nolig_pdb_path"], str(paths["ligand_output_dir"])
    )
    logger.info(f"Extracted {len(ligands_dict)} ligands → {paths['ligand_output_dir']}")
    return len(ligands_dict)


# ---------- small utils ----------

def norm(p):
    """Normalize a path to forward slashes for stable logging/keys."""
    return os.path.abspath(str(p)).replace("\\", "/")

def get_recenter_params(cfg):
    """Read early/fallback recentering knobs from config with safe defaults."""
    return {
        "EARLY_RECENTER_RATIO": float(cfg.get("EARLY_RECENTER_RATIO", 0.70)),
        "EARLY_RECENTER_MIN_EVAL": int(cfg.get("EARLY_RECENTER_MIN_EVAL", 10)),
        "EARLY_RECENTER_FAR_A": float(cfg.get("EARLY_RECENTER_FAR_A", 15.0)),
        "EARLY_RECENTER_MEDIAN_A": float(cfg.get("EARLY_RECENTER_MEDIAN_A", 10.0)),
        "ALLOW_BOX_EXPAND": bool(cfg.get("ALLOW_BOX_EXPAND", True)),
        "MAX_RECENTER_ATTEMPTS": int(cfg.get("MAX_RECENTER_ATTEMPTS", 3)),
    }
def load_aliases():
    # back-compat shim
    return get_atom_rules()

def main(pdb_file):
    base = os.path.splitext(pdb_file)[0]
    if base.endswith("_cleaned"):
        pdb_cleaned = base + ".pdb"  # already cleaned
    else:
        pdb_cleaned = base + "_cleaned.pdb"
    temp_fixed_pdb = base + "_fixed.pdb"
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    ligands_dir = os.path.join("processed_pdbs", f"{pdb_id}_ligands")
    try:
        shutil.copyfile(pdb_file, temp_fixed_pdb)
        logging.info(f"Copied PDB for fixing: {temp_fixed_pdb}")
        fix_pdb_elements(temp_fixed_pdb)

        ligands, ligand_coords = extract_and_remove_ligands(temp_fixed_pdb, pdb_cleaned, ligands_dir)

        if ligands:
            logging.info(f"Ligands removed for {pdb_file}. Ranking ligands by contacts...")
            for key, lines in ligands.items():
                logging.info(f"Ligand {key} has {len(lines)} atoms.")
            ranked_ligands = rank_ligands_by_atom_count(ligands)
            if ranked_ligands:
                # Select top ligand by contact count
                top_ligand_key, top_lines = ranked_ligands[0]
                logging.info(f"Top ligand by size: {top_ligand_key} with {len(top_lines)} atoms.")
                top_ligand_lines = top_lines

                ligand_output_dir = os.path.join(os.path.dirname(pdb_cleaned), f"{pdb_id}_ligands")
                os.makedirs(ligand_output_dir, exist_ok=True)
                ligand_path = os.path.join(ligand_output_dir, f"{pdb_id}_top_ligand.pdb")
                with open(ligand_path, 'w') as f:
                    f.writelines(ligands[top_ligand_key])
                logging.info(f"Top ligand saved to {ligand_path}")

                # Use coordinates of the top ligand only to compute box
                top_ligand_lines = ligands[top_ligand_key]
                top_ligand_coords = []
                for line in top_ligand_lines:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        top_ligand_coords.append((x, y, z))
                    except ValueError:
                        logging.warning(f"Invalid coordinates in ligand line: {line.strip()}")

                if not top_ligand_coords:
                    logging.warning("Top ligand has no valid coordinates. Falling back to P2Rank.")
                    center, box_size = get_box_from_p2rank_csv(pdb_cleaned)
                else:
                    center, box_size = compute_box_from_ligand_coords(top_ligand_coords)

            else:
                logging.warning("No ligands ranked, fallback to P2Rank.")
                center, box_size = get_box_from_p2rank_csv(pdb_cleaned)
        else:
            logging.info(f"No ligands in {pdb_cleaned}. Using P2Rank instead.")
            center, box_size = get_box_from_p2rank_csv(pdb_cleaned)

        if center and box_size:
            logging.info(f"{pdb_cleaned}: center={center}, box_size={box_size}")
            return center, box_size
        else:
            logging.warning(f"Box not determined for {pdb_cleaned}")
            return None, None

    finally:
        try:
            os.remove(temp_fixed_pdb)
            logging.info(f"Temporary file removed: {temp_fixed_pdb}")
        except OSError:
            logging.warning(f"Could not delete temp file: {temp_fixed_pdb}")
