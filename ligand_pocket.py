import os
import logging
from collections import defaultdict
from typing import List

from Bio.PDB import NeighborSearch, PDBParser

import pdb_fixer

logger = logging.getLogger(__name__)


def extract_and_remove_ligands(pdb_path, output_cleaned_pdb, ligands_dir):
    os.makedirs(ligands_dir, exist_ok=True)
    logging.info(
        "[extract] in=%s out=%s ldir=%s", pdb_path, output_cleaned_pdb, ligands_dir
    )

    ligands = defaultdict(list)
    ligand_coords = []  # collect all ligand atom coords for box calculation
    rules = pdb_fixer.get_atom_rules()
    retained_resnames = rules.retain_resnames  # already uppercased
    # [ions] retention counters
    metal_tokens = {"ZN", "MG", "MN", "FE", "CA", "CU", "CO", "NI", "HG"}
    variant_label = (os.environ.get("APO_HOLO_MODE") or "").strip().upper() or "legacy"
    kept_metals = 0
    stripped_metals = 0
    stripped_detail: List[str] = []

    # Do NOT retain common cryos/buffers: GOL/EDO/PG4/MPD/ACT/TRS/PO4/PEG → they’ll be extracted

    with open(pdb_path, "r") as infile, open(output_cleaned_pdb, "w") as outfile:
        for line in infile:
            if line.startswith("HETATM"):
                resname_raw = line[17:20]
                resname = (
                    resname_raw.strip().upper()
                )  # normalize for set membership & filenames
                chain = line[21]
                resnum = line[22:26].strip()
                elem = (line[76:78].strip() or resname).upper()
                if resname not in retained_resnames:
                    ligands[(chain, resname, resnum)].append(line)
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_coords.append((x, y, z))
                    except ValueError:
                        logging.warning(
                            f"Invalid ligand coordinates in line: {line.strip()}"
                        )
                    if elem in metal_tokens or resname in metal_tokens:
                        stripped_metals += 1
                        stripped_detail.append(
                            f"{resname}:{chain or '-'}:{resnum or '?'}"
                        )
                        logging.info(
                            "[ions.drop] reason=not_in_retain resname=%s chain=%s resSeq=%s element=%s variant=%s",
                            resname,
                            chain.strip() or "-",
                            resnum or "?",
                            elem or "?",
                            variant_label,
                        )
                    continue
                if elem in metal_tokens or resname in metal_tokens:
                    kept_metals += 1
            outfile.write(line)

    # Write separate ligand files (with element repair)
    for (chain, resname, resnum), lines in ligands.items():
        ligand_fname = os.path.join(ligands_dir, f"{resname}_{chain}{resnum}.pdb")
        fixed = pdb_fixer.fix_ligand_element_columns(lines)
        with open(ligand_fname, "w") as lf:
            lf.writelines(fixed)
        logging.info(
            f"Saved ligand {resname} {chain}{resnum} to {ligand_fname} (elements repaired)"
        )

    logging.info(f"Ligands extracted and removed from {pdb_path}.")
    logging.info(
        "[ions.summary] stage=extract_and_remove_ligands variant=%s kept=%d stripped=%d",
        variant_label,
        kept_metals,
        stripped_metals,
    )
    logging.info(
        "[ions.diff.input→cleaned] variant=%s lost=%s",
        variant_label,
        ",".join(stripped_detail) if stripped_detail else "none",
    )
    return ligands, ligand_coords


def rank_ligands_by_atom_count(ligands_dict):
    return sorted(ligands_dict.items(), key=lambda item: len(item[1]), reverse=True)


def calculate_ligand_protein_contacts(
    protein_pdb_path,
    ligand_lines,
    distance_cutoff=4.0,
):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_pdb_path)
    model = structure[0]

    # Collect all protein atoms
    protein_atoms = [
        atom for atom in model.get_atoms() if atom.get_parent().get_id()[0] == " "
    ]

    # Parse ligand atoms from lines
    ligand_atoms = []
    for line in ligand_lines:
        if line.startswith(("HETATM", "ATOM")):
            # Simple PDB atom line parsing for coords
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            # Create a dummy atom-like object
            class DummyAtom:
                def __init__(self, coord):
                    self.coord = coord

                def get_coord(self):
                    return self.coord

            ligand_atoms.append(DummyAtom((x, y, z)))

    # Use NeighborSearch to find protein atoms near ligand atoms
    ns = NeighborSearch(protein_atoms)
    contact_count = 0
    for latom in ligand_atoms:
        close_atoms = ns.search(latom.get_coord(), distance_cutoff)
        contact_count += len(close_atoms)

    return contact_count


def compute_box_from_ligand_coords(coords):
    if not coords:
        return None, None

    x_vals, y_vals, z_vals = zip(*coords)
    center = (
        sum(x_vals) / len(x_vals),
        sum(y_vals) / len(y_vals),
        sum(z_vals) / len(z_vals),
    )

    buffer = 5.0
    x_range = max(x_vals) - min(x_vals) + buffer
    y_range = max(y_vals) - min(y_vals) + buffer
    z_range = max(z_vals) - min(z_vals) + buffer

    box_size = (x_range, y_range, z_range)
    return center, box_size
