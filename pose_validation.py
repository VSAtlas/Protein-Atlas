import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from pathlib import Path
import os
def parse_pdbqt_coordinates(pdbqt_file):
    coords = []
    elements = []
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                try:
                    x = float(parts[5])
                    y = float(parts[6])
                    z = float(parts[7])
                except (IndexError, ValueError):
                    continue
                element = parts[-1]
                coords.append([x, y, z])
                elements.append(element)
    return np.array(coords), elements



def detect_hydrogen_bonds(protein_coords, protein_elements, ligand_coords, ligand_elements, max_dist=3.5):
    donors = {'N', 'O'}
    acceptors = {'O', 'N'}

    hbonds = 0
    for lig_atom, lig_elem in zip(ligand_coords, ligand_elements):
        if lig_elem not in donors and lig_elem not in acceptors:
            continue
        for prot_atom, prot_elem in zip(protein_coords, protein_elements):
            if prot_elem not in donors and prot_elem not in acceptors:
                continue
            dist = np.linalg.norm(lig_atom - prot_atom)
            if dist < max_dist:
                hbonds += 1
                break
    return hbonds

def evaluate_lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)

    return {
        "molecular_weight": mw,
        "logP": logp,
        "num_h_donors": h_donors,
        "num_h_acceptors": h_acceptors,
        "flags": {
            "mw": mw > 500,
            "logP": logp > 5,
            "h_donors": h_donors > 5,
            "h_acceptors": h_acceptors > 10
        }
    }
from Bio.PDB import PDBParser
import numpy as np

def extract_surface_atoms(pdbqt_path, max_dist_from_center=12.0, center=None):
    coords, _ = parse_pdbqt_coordinates(pdbqt_path)
    if center is not None:
        center = np.array(center)
        coords = np.array([c for c in coords if np.linalg.norm(c - center) <= max_dist_from_center])
    return coords

def attempt_fallback_recenter(fallback_ligands, receptor_pdbqt, docking_dir, stage_name, pocket_center, logger):
    if not fallback_ligands:
        logger.warning("No fallback ligands available.")
        return None, pocket_center, None

    sorted_ligands = sorted(fallback_ligands.items(), key=lambda x: x[1])  # sort by score
    best_lig_path = None
    best_center = pocket_center

    for ligand_path, score in sorted_ligands:
        fallback_pose_path = Path(docking_dir) / stage_name / f"{Path(ligand_path).stem}_{stage_name}.pdbqt"
        if not fallback_pose_path.exists():
            logger.warning(f"Fallback pose path missing: {fallback_pose_path}")
            continue

        coords, elements = parse_pdbqt_coordinates(fallback_pose_path)
        if len(coords) == 0:
            logger.warning(f"No coordinates found in fallback pose: {fallback_pose_path}")
            continue

        new_center = np.mean(coords, axis=0).tolist()
        surface_coords = extract_surface_atoms(receptor_pdbqt, center=new_center)

        result = validate_pose_pdbqt(
            protein_pdbqt=receptor_pdbqt,
            ligand_pdbqt=fallback_pose_path,
            pocket_center=new_center,
            clash_threshold=2.0,
            DIST_THRESHOLD_SURFACE=6.0,
            DIST_THRESHOLD_CENTROID=4.5,
            surface_atom_coords=surface_coords,
        )

        if result.get("valid", False):
            logger.info(f"Recentered using ligand {Path(ligand_path).stem} — validation passed.")
            best_score = score
            return str(fallback_pose_path), new_center, best_score

        # If not valid, keep this as best option to recenter anyway
        if best_lig_path is None:
            best_lig_path = fallback_pose_path
            best_center = new_center

    logger.warning("No fallback ligand passed validation — recentering using best-scoring ligand anyway.")
    return str(best_lig_path), best_center, sorted_ligands[0][1]


def validate_pose_pdbqt(
    protein_pdbqt,
    ligand_pdbqt,
    pocket_center,
    clash_threshold=2.0,
    CLASH_TOLERANCE=3,
    DIST_THRESHOLD_SURFACE=8.0,
    DIST_THRESHOLD_CENTROID=6.0,
    ligand_smiles=None,
    surface_atom_coords=None
):
    protein_coords, protein_elements = parse_pdbqt_coordinates(protein_pdbqt)
    ligand_coords, ligand_elements = parse_pdbqt_coordinates(ligand_pdbqt)

    if ligand_coords.size == 0:
        return {
            "valid": False,
            "reason": "Ligand has no heavy atoms"
        }

    # Step 1: Surface-based or centroid-based pocket distance
    used_fallback = False
    distance_to_surface = None
    distance_to_pocket = None

    if surface_atom_coords is not None and len(surface_atom_coords) > 0:
        dists = np.linalg.norm(
            ligand_coords[:, np.newaxis, :] - surface_atom_coords[np.newaxis, :, :], axis=2
        )
        min_dist = np.min(dists)
        distance_to_surface = float(min_dist)

        if min_dist > DIST_THRESHOLD_SURFACE:
            return {
                "valid": False,
                "distance_to_surface": distance_to_surface,
                "distance_to_pocket": None,
                "clash_count": None,
                "reason": f"Too far from pocket surface (> {DIST_THRESHOLD_SURFACE} Å)"
            }

        distance_to_pocket = distance_to_surface

    else:
        # Fallback: centroid-based distance
        used_fallback = True
        ligand_centroid = np.mean(ligand_coords, axis=0)
        distance_to_pocket = np.linalg.norm(ligand_centroid - np.array(pocket_center))

        if distance_to_pocket > DIST_THRESHOLD_CENTROID:
            return {
                "valid": False,
                "distance_to_surface": None,
                "distance_to_pocket": float(distance_to_pocket),
                "clash_count": None,
                "reason": f"Too far from pocket centroid (> {DIST_THRESHOLD_CENTROID} Å)"
            }

    # Step 2: Clash detection
    clash_count = 0
    for lig_atom in ligand_coords:
        dists = np.linalg.norm(protein_coords - lig_atom, axis=1)
        clash_count += np.sum(dists < clash_threshold)

    if clash_count > CLASH_TOLERANCE:
        return {
            "valid": False,
            "distance_to_surface": distance_to_surface,
            "distance_to_pocket": float(distance_to_pocket),
            "clash_count": int(clash_count),
            "reason": f"{clash_count} clashes detected (>{CLASH_TOLERANCE})"
        }

    # Step 3: Hydrogen bond detection
    hbond_count = detect_hydrogen_bonds(
        protein_coords, protein_elements,
        ligand_coords, ligand_elements
    )

    # Step 4: Lipinski rule evaluation
    lipinski_data = evaluate_lipinski(ligand_smiles) if ligand_smiles else None

    return {
        "valid": True,
        "distance_to_surface": distance_to_surface,
        "distance_to_pocket": float(distance_to_pocket),
        "clash_count": int(clash_count),
        "hydrogen_bonds": hbond_count,
        "lipinski": lipinski_data,
        "reason": "Valid pose (fallback used)" if used_fallback else "Valid pose"
    }

