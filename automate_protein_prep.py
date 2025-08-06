import sys
import os
import subprocess
from installation import load_config
from activesite import fix_pdb_elements
from activesite import ensure_model_records
import logging
import shutil
from logger_setup import setup_logger

# Load config once
config = load_config()

# === Configurable paths from config.txt ===
PHENIX_DIR = config["PHENIX_DIR"]
PHENIX_LIB_PATH = config["phenix_lib_path"]
PHENIX_CLEAN_SCRIPT = config["PHENIX_CLEAN_SCRIPT"]
INPUT_DIR = config["INPUT_DIR"]

# Add Phenix Python paths
sys.path.insert(0, PHENIX_LIB_PATH)
sys.path.insert(0, PHENIX_DIR)
from pathlib import Path

# PHENIX_DIR is absolute path from config, keep as is
PHENIX_DIR_PATH = Path(PHENIX_DIR)

# Use Path to join but do not normalize across drives — Path handles this correctly
REDUCE_EXE = Path(config.get("REDUCE_EXE")) if "REDUCE_EXE" in config else PHENIX_DIR_PATH / "reduce.exe"
print("[DEBUG] REDUCE_EXE =", REDUCE_EXE)

PDBTOOLS_BAT = PHENIX_DIR_PATH / "phenix.pdbtools.bat"
MOLPROBITY_BAT = PHENIX_DIR_PATH / "phenix.molprobity.bat"
PHENIX_PYTHON_BAT = PHENIX_DIR_PATH / "phenix.python.bat"
print("[DEBUG] PHENIX_DIR =", config["PHENIX_DIR"])
print("[DEBUG] PDBTOOLS_BAT =", PDBTOOLS_BAT)
print("[DEBUG] File exists:", os.path.exists(PDBTOOLS_BAT))

# Convert to strings with forward slashes (safe for subprocess on Windows)
REDUCE_EXE = str(REDUCE_EXE).replace("\\", "/")
PDBTOOLS_BAT = str(PDBTOOLS_BAT).replace("\\", "/")
MOLPROBITY_BAT = str(MOLPROBITY_BAT).replace("\\", "/")
PHENIX_PYTHON_BAT = str(PHENIX_PYTHON_BAT).replace("\\", "/")
print("[DEBUG] PHENIX_DIR =", config["PHENIX_DIR"])
print("[DEBUG] PDBTOOLS_BAT =", PDBTOOLS_BAT)
print("[DEBUG] File exists:", os.path.exists(PDBTOOLS_BAT))

def file_contains_hydrogens(pdb_path):
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")) and line[76:78].strip() == "H":
                    return True
    except Exception as e:
        logging.warning(f"Could not read {pdb_path} to check for H atoms: {e}")
    return False
def file_was_reduced(pdb_path):
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if "Reduce" in line and "protonation" in line.lower():
                    return True
    except Exception as e:
        logging.warning(f"Could not read {pdb_path} to check for Reduce header: {e}")
    return False

def filter_altlocs(pdb_input_path: str, pdb_output_path: str):
    """
    Filters alternate conformations in a PDB file.
    Keeps only one altLoc per atom: prefers ' ' (blank), then 'A', then first encountered.
    Logs all decisions.
    """
    lines = []
    atoms = {}  # key: (chain, resSeq, iCode, atom_name), value: dict altLoc->line
    removed_count = 0

    with open(pdb_input_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                atom_name = line[12:16]
                altLoc = line[16]
                chainID = line[21]
                resSeq = line[22:26].strip()
                iCode = line[26]

                key = (chainID, resSeq, iCode, atom_name.strip())
                if key not in atoms:
                    atoms[key] = {}
                atoms[key][altLoc] = line
            else:
                lines.append(line)

    # Now select one altLoc per atom and log decisions
    filtered_atoms = []
    for key, altloc_dict in atoms.items():
        altLocs = list(altloc_dict.keys())
        if len(altLocs) > 1:
            logging.info(f"AltLocs found for {key}: {altLocs}")

        if ' ' in altloc_dict:
            selected = ' '
        elif 'A' in altloc_dict:
            selected = 'A'
        else:
            selected = sorted(altloc_dict.keys())[0]

        if len(altloc_dict) > 1:
            removed = [alt for alt in altLocs if alt != selected]
            removed_count += len(removed)
            logging.info(f"Keeping altLoc '{selected}' for atom {key}, removed {removed}")

        filtered_atoms.append(altloc_dict[selected])

    filtered_atoms.sort(key=lambda l: (
        l[21],  # chain
        int(l[22:26]),  # residue number
        l[12:16].strip()  # atom name
    ))

    with open(pdb_output_path, 'w') as f:
        for line in lines:
            f.write(line)
        for atom_line in filtered_atoms:
            f.write(atom_line)

    logging.info(f"Filtered altLocs in {pdb_input_path}, output to {pdb_output_path}")
    logging.info(f"Total alternate conformers removed: {removed_count}")



def build_missing_loops(input_pdb, output_dir):
    """
    Uses MODELLER to fill in missing loops and residues in a PDB file.
    Returns the path to the filled structure or the original if filling fails.
    """
    from modeller import environ, log
    from modeller.scripts import complete_pdb
    import os

    output_pdb = os.path.join(output_dir, "modeller_filled.pdb")
    log.none()  # Suppress MODELLER log output
    logging.info(f"Running MODELLER to complete missing parts of {input_pdb}")

    try:
        env = environ()
        env.io.hetatm = True
        env.io.water = True
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')

        mdl = complete_pdb(env, input_pdb)  # this already fills missing loops and residues

        mdl.write(file=output_pdb)

        if os.path.exists(output_pdb):
            logging.info(f"MODELLER filled PDB saved to {output_pdb}")
            return output_pdb
        else:
            logging.warning(f"MODELLER did not produce expected output: {output_pdb}")
            return input_pdb

    except Exception as e:
        logging.warning(f"MODELLER failed on {input_pdb}: {e}")
        return input_pdb

def find_invalid_atoms(pdb_path):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_path)

    flagged = []
    for model in structure:
          for chain in model:
              for residue in chain:
                atoms = list(residue.get_atoms())
                atom_names = [a.get_name() for a in atoms]
                if 'CA' in atom_names and len(atoms) <= 2:
                    flagged.append((chain.id, residue.get_resname(), residue.id[1]))
    return flagged


def filter_invalid_chains(pdb_path, output_path):
    """
    Removes entire chains that do not contain at least one CA atom.
    """
    from collections import defaultdict

    chains = defaultdict(list)
    valid_chains = set()

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                chain_id = line[21]
                atom_name = line[12:16].strip()
                chains[chain_id].append(line)
                if atom_name in {"CA", "N", "C", "O"}:
                    valid_chains.add(chain_id)
            else:
                chains["HEADER"].append(line)  # Save non-atom lines separately

    with open(output_path, 'w') as f:
        for chain_id in chains:
            if chain_id == "HEADER" or chain_id in valid_chains:
                f.writelines(chains[chain_id])
            else:
                logging.warning(f"Skipping invalid chain '{chain_id}' (no CA atoms)")

def remove_unbonded_atoms(pdb_path):
    bonded_atoms = set()
    all_atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("CONECT"):
                parts = line.split()
                for atom_serial in parts[1:]:
                    bonded_atoms.add(atom_serial)
            elif line.startswith(("ATOM", "HETATM")):
                all_atoms.append(line)

    filtered = []
    for line in all_atoms:
        atom_serial = line[6:11].strip()
        if atom_serial in bonded_atoms or line[76:78].strip() != 'H':
            filtered.append(line)
        else:
            logging.info(f"Removed unbonded hydrogen: {line.strip()}")

    with open(pdb_path, 'w') as f:
        f.writelines(filtered)

def has_valid_chain(pdb_path):
    valid = False
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                valid = True
                break
    return valid


def run_prepare_receptor(input_pdb, output_pdbqt, config):
    mgltools_python = config["MGLTOOLS_PYTHON"]
    prepare_script = config["PREPARE_RECEPTOR_SCRIPT"]

    cmd = [
        mgltools_python,
        prepare_script,
        "-r", input_pdb,
        "-o", output_pdbqt,
        "-A", "none",  # disables hydrogens
        "-U", "nphs_lps_nonstdres"
    ]


    logging.info(f"Running prepare_receptor4.py: {' '.join(cmd)}")

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logging.error(f"prepare_receptor4 failed:\n{result.stderr}")
        return False

    logging.info(f"Receptor prepared: {output_pdbqt}")
    return True


def run_phenix_pdbtools(input_pdb, output_pdb, remove_waters=True):
    input_pdb = os.path.abspath(input_pdb).replace("\\", "/")
    output_pdb = os.path.abspath(output_pdb).replace("\\", "/")

    cmd = f'"{PDBTOOLS_BAT}" "{input_pdb}" output.file_name="{output_pdb}" remove="resname HOH"'
    print(f"Running: {cmd}")
    logging.info(f"Running: {cmd}")

    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print("ERROR: phenix.pdbtools failed")
        print(result.stderr)
        logging.warning(f"{result.stderr}\nERROR: phenix.pdbtools failed")
        return False

    print("phenix.pdbtools completed successfully")
    logging.info("phenix.pdbtools completed successfully")
    return True

def run_openbabel_add_h(input_pdb, output_pdb):
    cmd = ["obabel", input_pdb, "-O", output_pdb, "-h"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Open Babel failed:\n{result.stderr}")
    else:
        print("open babel completed successfully")


def assign_protonation_states(input_pdb, output_pdb, reduce_exe=None):
    """
    Uses Reduce to assign protonation states and add hydrogens.
    Falls back to Open Babel if Reduce fails.
    """
    import shutil

    input_pdb = os.path.abspath(input_pdb).replace("\\", "/")
    output_pdb = os.path.abspath(output_pdb).replace("\\", "/")

    already_reduced = file_was_reduced(input_pdb)
    has_hydrogens = file_contains_hydrogens(input_pdb)

    reduce_flags = ["-BUILD", "-quiet"]
    if has_hydrogens and not already_reduced:
        logging.info(f"Hydrogens found but Reduce not detected. Running Reduce with -nohyd.")
        reduce_flags += ["-nohyd"]

    # Attempt 1: Direct Reduce
    try:
        with open(output_pdb, "w") as out:
            result = subprocess.run(
                [REDUCE_EXE] + reduce_flags + [input_pdb],
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                cwd=os.path.dirname(REDUCE_EXE)
            )
        if result.stderr:
            logging.warning(f"Reduce stderr: {result.stderr}")
        logging.info("Reduce completed successfully on first attempt.")
        return output_pdb

    except subprocess.CalledProcessError as e:
        logging.warning(f"Reduce failed: {e.stderr}")
        fix_pdb_elements(input_pdb)

    # Attempt 2: phenix.python.bat -m phenix.reduce
    try:
        with open(output_pdb, "w") as out:
            result = subprocess.run(
                [PHENIX_PYTHON_BAT, "-m", "phenix.reduce"] + reduce_flags + [input_pdb],
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        if result.stderr:
            logging.warning(f"Retry Reduce stderr: {result.stderr}")
        logging.info("Reduce completed successfully after retry.")
        return output_pdb

    except subprocess.CalledProcessError as e:
        logging.warning("Reduce retry also failed. Trying a third Reduce attempt...")

    # Attempt 3: Reduce again with renamed input
    try:
        temp_input = input_pdb.replace(".pdb", "_retry3.pdb")
        shutil.copy(input_pdb, temp_input)

        with open(output_pdb, "w") as out:
            result = subprocess.run(
                [REDUCE_EXE] + reduce_flags + [temp_input],
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                cwd=os.path.dirname(REDUCE_EXE)
            )

        if result.stderr:
            logging.warning(f"Third Reduce attempt stderr: {result.stderr}")
        logging.info("Reduce completed successfully on third attempt.")
        return output_pdb

    except subprocess.CalledProcessError as e:
        logging.error("All Reduce attempts failed. Trying Open Babel hydrogenation fallback...")

    # Final fallback: Open Babel
    try:
        run_openbabel_add_h(input_pdb, output_pdb)
        if os.path.exists(output_pdb):
            logging.info("Hydrogens added using Open Babel fallback.")
        else:
            raise FileNotFoundError("Open Babel failed to write file.")
    except Exception as babel_error:
        logging.error(f"Open Babel hydrogenation failed: {babel_error}")
        shutil.copy(input_pdb, output_pdb)
        logging.warning("Protonation skipped entirely.")

    return output_pdb



def guess_element_from_atom_name(name: str) -> str:
    name = name.strip()
    if len(name) == 4:
        # If atom name is right-justified (columns 13-16), e.g. " CA "
        name = name[1:] if name[0].isdigit() or name[0] == ' ' else name
    name = name.upper()

    # Two-letter elements first
    two_letter = {"ZN", "FE", "MG", "MN", "CL", "NA", "CA", "CU", "CO"}
    if name[:2] in two_letter:
        return name[:2]
    elif name[0] in {"C", "H", "O", "N", "S", "P"}:
        return name[0]
    else:
        return "C"  # fallback

def strip_incomplete_residues(pdb_path, output_path, min_atoms=3):
    from Bio.PDB import PDBParser, PDBIO, Select
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_path)

    class GoodResidueSelect(Select):
        def accept_residue(self, residue):
            return len(list(residue.get_atoms())) >= min_atoms

    # === 🔧 Fix elements before writing
    for atom in structure.get_atoms():
        raw_name = atom.get_name()
        atom.element = guess_element_from_atom_name(raw_name)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, GoodResidueSelect())


def strip_nonstandard_residues(input_pdb: str, output_pdb: str):
    standard_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "HID", "HIE", "HIP", "SEC", "PYL", "MSE"
    }

    important_hetatm_residues = {
        "ZN", "MG", "NA", "CA", "MN", "CO", "FE", "CU", "HEM", "FAD", "FMN", "COA", "SAM"
    }
    known_peptides = {"CIR", "NAG", "MAN"}  # or dynamically extract from ligand list

    removed = set()
    kept_lines = []

    with open(input_pdb, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                if resname in standard_residues:
                    kept_lines.append(line)
                else:
                    removed.add(resname)
            elif line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname in important_hetatm_residues:
                    kept_lines.append(line)
                else:
                    removed.add(resname)
            else:
                kept_lines.append(line)

    with open(output_pdb, 'w') as f:
        f.writelines(kept_lines)

    logging.info(f"Removed nonstandard residues: {sorted(removed)}")
    return len(removed), output_pdb
def clean_pdb(pdb_file, output_root):
    if not os.access(output_root, os.W_OK):
        raise PermissionError(f"Cannot write to output directory: {output_root}")
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    if pdb_id.endswith("_cleaned"):
        pdb_id = pdb_id[:-8]  # Remove trailing '_cleaned'
    output_dir = os.path.join(output_root, pdb_id)
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: copy original pdb to a working file
    working_pdb = os.path.join(output_dir, f"{pdb_id}_working.pdb")
    shutil.copyfile(pdb_file, working_pdb)

    # Step 2: filter alternate locations
    filtered_pdb = os.path.join(output_dir, f"{pdb_id}_filtered.pdb")
    filter_altlocs(working_pdb, filtered_pdb)

    # Step 2.5: remove nonstandard residues
    stripped_pdb = os.path.join(output_dir, f"{pdb_id}_stripped.pdb")
    removed_count, _ = strip_nonstandard_residues(filtered_pdb, stripped_pdb)
    logging.info(f"Removed {removed_count} nonstandard residue lines.")

    # Step 3: fix element columns
    # Step 3: fix incomplete residues AND element columns
    strip_incomplete_residues(stripped_pdb, stripped_pdb)

    # Step 3.5: fill missing loops
    loop_fixed_pdb = build_missing_loops(stripped_pdb, output_dir)
    # Step 3.6: detect incomplete residues with very few atoms
    invalid_residues = find_invalid_atoms(loop_fixed_pdb)
    if invalid_residues:
        for chain_id, resname, resnum in invalid_residues:
            logging.warning(f"Incomplete residue: {resname} {chain_id}{resnum} has <=2 atoms and includes CA")
    else:
        logging.info("No incomplete CA residues found.")

    subprocess.run([
        "cmd.exe", "/c",
        PHENIX_PYTHON_BAT,
        PHENIX_CLEAN_SCRIPT,
        loop_fixed_pdb,
        output_root
    ])

    output_dir = os.path.abspath(output_dir)
    # Step 4: remove unbonded atoms + filter broken chains
    debulked_pdb = os.path.join(output_dir, f"{pdb_id}_debulked.pdb")
    shutil.copyfile(loop_fixed_pdb, debulked_pdb)
    remove_unbonded_atoms(debulked_pdb)

    chain_validated_pdb = os.path.join(output_dir, f"{pdb_id}_validated.pdb")
    filter_invalid_chains(debulked_pdb, chain_validated_pdb)

    # Step 5: protonation + final cleanup
    reduced_pdb = os.path.join(output_dir, f"{pdb_id}_reduced.pdb")
    assign_protonation_states(chain_validated_pdb, reduced_pdb)
    if not file_contains_hydrogens(reduced_pdb):
        logging.error(f"[FATAL] Reduced file missing hydrogens: {reduced_pdb}")
        return None

    fix_pdb_elements(reduced_pdb)  # Final element fix

    cleaned_pdb = os.path.abspath(os.path.join(output_dir, f"{pdb_id}_cleaned.pdb")).replace("\\", "/")

    if not run_phenix_pdbtools(reduced_pdb, cleaned_pdb):
        return None

    # Fix any broken element columns *after* Phenix touches the file
    fix_pdb_elements(cleaned_pdb)
    assert file_contains_hydrogens(cleaned_pdb), \
        f"[FATAL] Cleaned file lost hydrogens: {cleaned_pdb}"

    molprobity_log = os.path.join(output_dir, f"{pdb_id}_molprobity.log")
    with open(molprobity_log, "w") as out:
        subprocess.run([MOLPROBITY_BAT, cleaned_pdb], stdout=out)

    print(f"Cleaned: {cleaned_pdb}")
    print(f"Input PDB: {pdb_file}")
    print(f"Cleaned PDB will be written to: {cleaned_pdb}")
    logging.info(f"Cleaned: {cleaned_pdb}")
    logging.info(f"Input PDB: {pdb_file}")
    logging.info(f"Cleaned PDB will be written to: {cleaned_pdb}")

    #final check that pdb is cleaned
    from Bio.PDB.PDBParser import PDBParser
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("validate", cleaned_pdb)
        atoms = list(structure.get_atoms())
        assert len(atoms) > 0, f"[FATAL] Final cleaned PDB has no atoms: {cleaned_pdb}"
    except Exception as e:
        logging.error(f"[FATAL] Could not parse cleaned PDB: {e}")

    return cleaned_pdb


def main(pdb_filename, output_dir=r"./processed_pdbs"):
    try:
        if os.path.isabs(pdb_filename) and os.path.isfile(pdb_filename):
            pdb_path = pdb_filename
        else:
            pdb_path = os.path.join(INPUT_DIR, pdb_filename)
        pdb_id = os.path.splitext(pdb_filename)[0].upper()

        if not os.path.isfile(pdb_path):
            logging.error(f"ERROR: File does not exist: {pdb_path}")
            return None

        os.makedirs(output_dir, exist_ok=True)
        cleaned_pdb = clean_pdb(pdb_path, output_dir)

        if not cleaned_pdb:
            logging.error(f"ERROR: Cleaning failed for {pdb_filename}")
            return None

        # ===prepare receptor PDBQT ===
        output_pdbqt = os.path.join(output_dir, pdb_id, f"{pdb_id}.pdbqt")
        if not run_prepare_receptor(cleaned_pdb, output_pdbqt, config):
            logging.error(f"ERROR: Failed to prepare receptor PDBQT for {pdb_id}")
            return None

        logging.info(f"Prepared receptor PDBQT: {output_pdbqt}")
        return cleaned_pdb, output_pdbqt

    except Exception as e:
        logging.exception(f"[FATAL] automate_protein_prep.main() failed: {e}")
        return None

# Example usage
# if __name__ == "__main__":
#     test_pdb = "1a3n.pdb"
#     result = main(test_pdb)
#     if result:
#         print("Final cleaned PDB:", result)
