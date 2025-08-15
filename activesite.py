import os
import csv
import subprocess
import shutil
import logging
from installation import load_config
from collections import defaultdict
from logger_setup import setup_logger



# Load config
config = load_config()
P2RANK_DIR = config.get("P2RANK_PATH")


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
    def __init__(self):
        super().__init__()

    def get_atom_element(self, atom):
        name = atom.get_name().strip().upper()

        # Common 2-letter elements first
        if name.startswith(("ZN", "FE", "CA", "MG", "MN", "CL", "CU", "NI", "CO")):
            return name[:2]

        # Broad hydrogen matching: names starting with H or digit+H
        if name.startswith("H") or name[0] in "123456789" and name[1] == "H":
            return "H"

        # Common heavy atoms
        if name.startswith("C"):
            return "C"
        if name.startswith("O"):
            return "O"
        if name.startswith("N"):
            return "N"
        if name.startswith("S"):
            return "S"
        if name.startswith("P"):
            return "P"

        # Fallback — assume carbon to prevent crash but warn
        logging.warning(f"Unknown atom name '{name}', defaulting to 'C'")
        return "C"

    def accept_atom(self, atom):
        # Fix element symbol in the atom object itself
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

    io.save(output_path, fixer)
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


def extract_and_remove_ligands(pdb_path, output_cleaned_pdb, ligands_dir):
    os.makedirs(ligands_dir, exist_ok=True)
    ligands = defaultdict(list)
    ligand_coords = []  # collect all ligand atom coords for box calculation
    retained_resnames = {'HOH', 'SO4', 'MG', 'ZN', 'HEM', 'FAD', 'NAD', 'FE', 'CU', 'MN', 'CO', 'CA', 'CL'}

    with open(pdb_path, 'r') as infile, open(output_cleaned_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith('HETATM'):
                resname = line[17:20].strip()
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

    # Write separate ligand files
    for (chain, resname, resnum), lines in ligands.items():
        ligand_fname = os.path.join(ligands_dir, f"{resname}_{chain}{resnum}.pdb")
        with open(ligand_fname, 'w') as lf:
            lf.writelines(lines)
        logging.info(f"Saved ligand {resname} {chain}{resnum} to {ligand_fname}")

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

    jar_path = os.path.join(P2RANK_DIR, "bin", "p2rank.jar")
    if not os.path.isfile(jar_path):
        logging.error(f"Missing p2rank.jar at: {jar_path}")
        return None, None

    p2rank_executable = os.path.join(P2RANK_DIR, "prank.bat")
    cmd = [p2rank_executable, "predict", "-f", pdb_file]

    try:
        subprocess.run(cmd, check=True, shell=True)
        logging.info(f"P2Rank ran successfully for {pdb_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"P2Rank failed: {e}")
        return None, None

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
    center, box_size = detect_active_site(cleaned_pdb)
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
                top_ligand_key, top_score = ranked_ligands[0]
                logging.info(f"Top ligand by size: {top_ligand_key} with {top_score} atoms.")

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
