from rdkit import Chem
from rdkit.Chem import Descriptors
from input_and_export_functions import extract_best_score
import numpy as np

def parse_pdbqt_coordinates(pdbqt_file):
    coords, elements = [], []
    with open(pdbqt_file, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    continue
                elem = (line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[0])
                coords.append([x, y, z]); elements.append(elem)
    return np.array(coords), elements


# --- Minimal readers and self-docking RMSD ---

def parse_pdb_coordinates(pdb_file):
    """
    Heavy-atom PDB reader using fixed columns (like parse_pdbqt_coordinates but for PDB).
    Returns (N,3) float array and element list (H filtered out).
    """
    coords, elements = [], []
    with open(pdb_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    continue
                elem = (line[76:78].strip().upper() if len(line) >= 78 else line[12:16].strip()[0].upper())
                if elem == "H":
                    continue
                coords.append([x, y, z]); elements.append(elem)
    return np.array(coords, dtype=float), elements

def compute_redock_rmsd(crystal_lig_path: str, docked_pdbqt_path: str):
    """
    Compute heavy-atom Kabsch RMSD between the extracted crystal ligand (PDB/MOL2/SDF-as-PDB formatted)
    and the best-scoring docked pose (PDBQT). Assumes atom ordering is consistent; truncates to min length.
    Returns float RMSD in Å, or None on failure.
    """
    try:
        # Prefer PDB; if the extracted file is MOL2/SDF but written with PDB-like columns, parse still works.
        if crystal_lig_path.lower().endswith(".pdb"):
            ref_xyz, _ = parse_pdb_coordinates(crystal_lig_path)
        else:
            # best-effort parse using PDB columns (your extractor often writes PDB too)
            ref_xyz, _ = parse_pdb_coordinates(crystal_lig_path)

        prd_xyz, _ = parse_pdbqt_coordinates(docked_pdbqt_path)
        if ref_xyz.size == 0 or prd_xyz.size == 0:
            return None

        n = min(ref_xyz.shape[0], prd_xyz.shape[0])
        if n < 5:
            return None

        return _kabsch(ref_xyz[:n], prd_xyz[:n])
    except Exception:
        return None
def compute_self_rmsd(pdbqt_path: str):
    blocks = _split_pdbqt_models(pdbqt_path)
    if len(blocks) < 2:
        return None
    A = _coords_from_block(blocks[0])
    B = _coords_from_block(blocks[1])
    if A.shape[0] == 0 or B.shape[0] == 0:
        return None
    n = min(A.shape[0], B.shape[0])
    if n < 5:
        return None
    return _kabsch(A[:n], B[:n])

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


def attempt_fallback_recenter(
    fallback_ligands: dict,             # {orig_ligand_path_norm -> pose_out_path_norm}
    receptor_pdbqt: str,
    docking_dir: str,
    stage_name: str,
    pocket_center: tuple,
    logger,
    exclude_basenames: set = None,      # NEW
    max_candidates: int = 5,            # try a handful
):
    """
    Returns (fb_pose_path, new_center, best_score, chosen_ligand_path) where chosen_ligand_path
    is the original ligand we used to propose the new center. If no center found, returns (None, None, None, None).
    """
    import os
    exclude_basenames = exclude_basenames or set()

    candidates = []
    # Build a ranked list of candidates by (is_valid_first, score_asc, distance_asc)
    for lig_path, pose_path in fallback_ligands.items():
        try:
            # Your existing helper; if not available here, import from input_and_export_functions
            score = extract_best_score(pose_path) if os.path.exists(pose_path) else None
        except Exception:
            score = None

        # Validate pose relative to current (possibly wrong) center to get a distance
        res = validate_pose_pdbqt(
            protein_pdbqt=receptor_pdbqt,
            ligand_pdbqt=pose_path,
            pocket_center=pocket_center,
            clash_threshold=2.0,
            CLASH_TOLERANCE=3,
            DIST_THRESHOLD_SURFACE=6.0,
            DIST_THRESHOLD_CENTROID=4.5,
            surface_atom_coords=None
        )
        dist = res.get("distance_to_pocket", float("inf"))
        valid = bool(res.get("valid", False))
        base = os.path.basename(lig_path)

        if base in exclude_basenames:
            continue

        # sort key: valid first, then best score (lower), then closer distance
        score_key = (0 if valid else 1, (score if isinstance(score, (int, float)) else float("inf")), dist)
        candidates.append((score_key, lig_path, pose_path, score, dist))

    if not candidates:
        logger.warning("Fallback recenter: no eligible candidates after exclusions.")
        return None, None, None, None

    candidates.sort(key=lambda x: x[0])
    tried = 0
    for _, lig_path, pose_path, score, dist in candidates:
        tried += 1
        # Propose new center from this pose (ligand centroid)
        try:
            new_center = _pose_centroid_from_pdbqt(pose_path)  # helper below
        except Exception as e:
            logger.warning(f"Fallback recenter: failed to compute centroid for {os.path.basename(pose_path)}: {e}")
            new_center = None

        logger.info(f"Recenter candidate {tried}: {os.path.basename(lig_path)} | score={score} | dist={dist:.2f} | center={new_center}")
        if new_center is not None:
            logger.info(f"Recentered using ligand {os.path.basename(lig_path)} — proposed new center {new_center}.")
            return pose_path, tuple(map(float, new_center)), score, lig_path

        if tried >= max_candidates:
            break

    logger.warning("Fallback recenter: no workable candidate produced a new center.")
    return None, None, None, None


def _pose_centroid_from_pdbqt(pdbqt_path: str):
    """Compute centroid of HETATM/ATOM coordinates in a PDBQT pose."""
    xs, ys, zs = [], [], []
    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                # columns 31-38, 39-46, 47-54 (1-based PDB), but PDBQT is similar
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    xs.append(x); ys.append(y); zs.append(z)
                except ValueError:
                    continue
    if not xs:
        return None
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))


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

# =========================
# RMSD pose filtering utils
# =========================
import re
from typing import List, Tuple, Optional

def _parse_vina_score_from_lines(lines: List[str]) -> Optional[float]:
    """
    Parse the 'REMARK VINA RESULT:   -8.7  0.0  0.0' score in a block.
    Returns float or None if not found.
    """
    pat = re.compile(r"REMARK\s+VINA\s+RESULT:\s*([-\d\.]+)")
    for ln in lines:
        m = pat.search(ln)
        if m:
            try:
                return float(m.group(1))
            except Exception:
                return None
    return None

def _split_pdbqt_models(pdbqt_path: str) -> List[List[str]]:
    """
    Split a multi-model PDBQT into a list of model blocks (list-of-lines per model).
    If no MODEL/ENDMDL tags exist, treat entire file as a single-model block.
    """
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        all_lines = f.readlines()

    blocks: List[List[str]] = []
    cur: List[str] = []
    in_model = False
    for ln in all_lines:
        if ln.startswith("MODEL"):
            if cur:
                blocks.append(cur)
                cur = []
            in_model = True
            cur.append(ln)
        elif ln.startswith("ENDMDL"):
            cur.append(ln)
            blocks.append(cur)
            cur = []
            in_model = False
        else:
            cur.append(ln)

    if in_model:
        # dangling model (no ENDMDL): finalize
        blocks.append(cur)
    elif not blocks:
        # no model tags: whole file is one "block"
        blocks = [all_lines]

    return blocks

def _coords_from_block(block: List[str]) -> np.ndarray:
    xs, ys, zs = [], [], []
    for ln in block:
        if (ln.startswith("ATOM") or ln.startswith("HETATM")) and len(ln) >= 54:
            try:
                x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            except Exception:
                parts = ln.split()
                if len(parts) >= 8:
                    try:
                        x = float(parts[5]); y = float(parts[6]); z = float(parts[7])
                    except Exception:
                        continue
                else:
                    continue
            xs.append(x); ys.append(y); zs.append(z)
    if not xs:
        return np.zeros((0,3), dtype=float)
    return np.vstack([xs, ys, zs]).T.astype(float)

def _kabsch(P: np.ndarray, Q: np.ndarray) -> float:
    """
    Return RMSD between two Nx3 point sets after optimal superposition.
    Assumes P and Q have same shape and atom order.
    """
    if P.shape != Q.shape or P.size == 0:
        return float("inf")

    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)

    # covariance and SVD
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    # proper rotation (determinant check)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    U = V @ D @ Wt

    P_rot = Pc @ U
    diff = P_rot - Qc
    return float(np.sqrt((diff * diff).sum() / P.shape[0]))

def cluster_models_by_rmsd(blocks: List[List[str]], rmsd_tol: float = 2.0) -> List[int]:
    """
    Greedy clustering: keep the first model, then keep any model whose RMSD
    to all kept models is >= rmsd_tol. Returns indices of kept models.
    """
    if not blocks:
        return []

    kept = []
    kept_coords = []

    for i, b in enumerate(blocks):
        C = _coords_from_block(b)
        # skip empty models
        if C.shape[0] == 0:
            continue
        if not kept:
            kept.append(i)
            kept_coords.append(C)
            continue
        # compare to all previously kept coords
        min_r = float("inf")
        for KC in kept_coords:
            if KC.shape[0] != C.shape[0]:
                # different atom counts: fallback to centroid distance
                r = float(np.linalg.norm(C.mean(axis=0) - KC.mean(axis=0)))
            else:
                r = _kabsch(C, KC)
            if r < min_r:
                min_r = r
        if min_r >= rmsd_tol:
            kept.append(i)
            kept_coords.append(C)

    return kept

def _sort_model_indices_by_score(blocks: List[List[str]]) -> List[int]:
    """
    Sort model indices ascending by score (more negative first).
    Falls back to original order if no scores found.
    """
    scored = []
    for i, b in enumerate(blocks):
        sc = _parse_vina_score_from_lines(b)
        scored.append((i, sc))
    if all(sc is None for _, sc in scored):
        return list(range(len(blocks)))
    # sort by score (lower/more negative first), None at end
    return [i for (i, sc) in sorted(scored, key=lambda t: (t[1] is None, t[1]))]

def filter_and_rewrite_poses_by_rmsd(pdbqt_path: str,
                                     rmsd_tol: float = 2.0,
                                     max_models: Optional[int] = None) -> Tuple[int, int]:
    """
    In-place filter: keep only unique poses >= rmsd_tol apart (Kabsch RMSD).
    Optionally cap to top-N by score after dedup.
    Returns (kept_count, removed_count).
    No-ops if the file has <2 models.
    """
    blocks = _split_pdbqt_models(pdbqt_path)
    if len(blocks) <= 1:
        return (len(blocks), 0)

    # first, sort by score so greedy takes best first
    order = _sort_model_indices_by_score(blocks)
    blocks_sorted = [blocks[i] for i in order]

    kept_sorted_idx = cluster_models_by_rmsd(blocks_sorted, rmsd_tol=rmsd_tol)

    # map back to original indices
    kept_idx = [order[i] for i in kept_sorted_idx]
    # optional cap
    if max_models is not None and len(kept_idx) > max_models:
        # re-sort kept by score (ascending) and take top-N
        kept_idx = [i for i in _sort_model_indices_by_score([blocks[j] for j in kept_idx])][:max_models]

    # write back (preserve MODEL/ENDMDL if present)
    kept_blocks = [blocks[i] for i in kept_idx]
    with open(pdbqt_path, "w", encoding="utf-8", errors="ignore") as f:
        for bi, blk in enumerate(kept_blocks):
            # Ensure MODEL/ENDMDL wrappers
            has_model = any(ln.startswith("MODEL") for ln in blk)
            has_end   = any(ln.startswith("ENDMDL") for ln in blk)
            if not has_model:
                f.write(f"MODEL        {bi+1}\n")
            f.writelines(blk)
            if not has_end:
                f.write("ENDMDL\n")

    removed = max(0, len(blocks) - len(kept_blocks))
    return (len(kept_blocks), removed)
