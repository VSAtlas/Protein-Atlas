import subprocess
from input_and_export_functions import extract_best_score
from pathlib import Path
import os



def select_for_next_stage(docking_mode, i, stages, scores, logger):
    """
    Select a fraction of ligands to carry forward to the next stage based on mode and rank.
    Returns the list of next-stage ligands.
    """
    if not scores:
        return []

    percentages = {
        "discovery": [1.0, 0.1, 0.01, 0.001, 0.001],
        "polypharmacology": [1.0, 0.05, 0.005]
    }.get(docking_mode, [1.0] * len(stages))

    pct = percentages[i + 1] if i + 1 < len(percentages) else 0.01
    num_to_select = int(len(scores) * pct)
    if num_to_select < 1:
        logger.warning(f"Percentage {pct*100:.5f}% yielded <1 ligand. Using best-scoring ligand.")
        num_to_select = 1

    top_ligands = sorted(((l, s) for l, s in scores.items() if isinstance(s, (int, float))), key=lambda x: x[1])[:num_to_select]
    next_list = [l for l, _ in top_ligands if l]
    logger.info(f"Selected top {len(next_list)} ligands ({pct * 100:.5f}%) for next stage.")
    return next_list

def run_docking_task(vina_exe, config_path, ligand_name, out_path):
    try:
        if not Path(vina_exe).exists():
            raise FileNotFoundError(f"[!] AutoDock Vina not found at: {vina_exe}")
        subprocess.run([vina_exe, "--config", config_path], check=True)
        score = extract_best_score(out_path)
        return ligand_name, score
    except FileNotFoundError as e:
        print(f"❌ File error for {ligand_name}: {e}")
        return ligand_name, None
    except subprocess.CalledProcessError as e:
        print(f"❌ Docking failed for {ligand_name}: {e}")
        return ligand_name, None

def extract_models(pdbqt_path):
    with open(pdbqt_path) as f:
        lines = f.readlines()

    models = []
    current = []
    for line in lines:
        if line.startswith("MODEL"):
            current = [line]
        elif line.startswith("ENDMDL"):
            current.append(line)
            models.append(current[:])
        else:
            current.append(line)
    return models

def validate_all_poses(pdbqt_path, receptor_pdbqt, center, surface_coords, validate_fn):
    pdbqt_path = Path(pdbqt_path)
    receptor_pdbqt = Path(receptor_pdbqt)
    models = extract_models(pdbqt_path)
    best_valid_score = None
    best_valid_model = None
    temp_paths = []

    for i, model_lines in enumerate(models):
        temp_path = Path(str(pdbqt_path).replace(".pdbqt", f"_model{i + 1}.pdbqt"))
        with open(temp_path, "w") as f:
            f.writelines(model_lines)
        temp_paths.append(temp_path)

        result = validate_fn(
            protein_pdbqt=receptor_pdbqt,
            ligand_pdbqt=temp_path,
            pocket_center=center,
            clash_threshold=2.0,
            CLASH_TOLERANCE=3,
            DIST_THRESHOLD_SURFACE=6.0,
            DIST_THRESHOLD_CENTROID=4.5,
            surface_atom_coords=surface_coords,
        )

        if result.get("valid", False):
            score = extract_best_score(temp_path)
            if best_valid_score is None or score < best_valid_score:
                best_valid_score = score
                best_valid_model = temp_path

    # Now clean up all intermediate model files
    for temp_path in temp_paths:
        if os.path.exists(temp_path):
            os.remove(temp_path)

    return best_valid_model, best_valid_score

