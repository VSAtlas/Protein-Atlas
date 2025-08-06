import os
import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

from protein_functions import detect_active_site, convert_to_pdbqt
from input_and_export_functions import validate_config
import automate_protein_prep
from prep_ligands import prep_ligands_from_pdb, is_valid_ligand
from activesite import extract_and_remove_ligands
import logging
from pose_validation import validate_pose_pdbqt, extract_surface_atoms, attempt_fallback_recenter
from pathlib import Path

from input_and_export_functions import (
    load_inputs, define_docking_stages, write_score_summary_to_csv,
    extract_best_score, generate_config
)
from run_vina import run_docking_task
def main():
    print("MODELLER is working with license.")
    cfg = load_inputs()
    validate_config(cfg)

    print("Loaded config keys:", list(cfg.keys()))
    docking_mode = cfg.get("DOCKING_MODE", "discovery").lower()
    stages = define_docking_stages(docking_mode)

    pdb_files = [
        f for f in os.listdir(cfg["INPUT_DIR"])
        if f.endswith(".pdb") and "_nolig" not in f.lower()
    ]
    print("Working directory:", os.getcwd())

    overall_start = time.time()
    with tqdm(total=len(pdb_files), desc="Processing Proteins", unit="protein") as protein_pbar:
        for pdb_file in pdb_files:
            # Normalize PDB ID and prepare logging first
            base_id = os.path.splitext(pdb_file)[0]  # e.g., "1a3n_cleaned"
            pdb_id_raw = base_id.replace("_cleaned", "")  # normalized ID for consistent folder/log naming

            # Set up logger before doing anything else
            logger = logging.getLogger(pdb_id_raw)
            logger.setLevel(logging.DEBUG)
            log_file_path = os.path.join(cfg["DOCKED_DIR"], pdb_id_raw, "protein.log")
            os.makedirs(os.path.dirname(log_file_path), exist_ok=True)
            if logger.hasHandlers():
                logger.handlers.clear()
            file_handler = logging.FileHandler(log_file_path)
            file_handler.setLevel(logging.DEBUG)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

            # Log file and PDB being processed
            logger.info(f"Processing protein: {pdb_file}")
            logger.info(f"Normalized PDB ID: {pdb_id_raw}")

            # === Extract ligands first ===
            # === Define ligand_output_dir before using it ===
            pdb_id = base_id
            ligand_output_dir = Path(cfg["OUTPUT_DIR"]) / pdb_id / f"{pdb_id}_cleaned_ligands"
            ligand_output_dir.mkdir(parents=True, exist_ok=True)
            ligands_mol2_dir = Path(cfg["LIGANDS_MOL2_DIR"]) / pdb_id
            prepped_ligands_dir = Path(cfg["OUTPUT_LIGANDS_DIR"]) / pdb_id
            pdb_path = os.path.join(cfg["INPUT_DIR"], pdb_file)
            nolig_pdb_path = os.path.join(cfg["OUTPUT_DIR"], f"{base_id}_nolig.pdb")
            # Clear old malformed ligand log for this protein
            malformed_log = ligands_mol2_dir / "malformed_ligands.txt"
            if malformed_log.exists():
                malformed_log.unlink()

            ligands_dict, _ = extract_and_remove_ligands(
                pdb_path,
                nolig_pdb_path,
                str(ligand_output_dir)
            )
            logger.info(f"Extracted {len(ligands_dict)} ligands to {ligand_output_dir}")

            # === Now clean the ligand-free PDB ===
            cleaned_pdb_path = Path(cfg["OUTPUT_DIR"]) / f"{base_id}_nolig_cleaned.pdb"
            receptor_pdbqt_path = Path(cfg["PDBQT_DIR"]) / f"{base_id}_receptor.pdbqt"
            if cleaned_pdb_path.exists() and receptor_pdbqt_path.exists() and not cfg.get("FORCE_REPROCESS", False):
                logger.info("Skipping reprocessing: using existing cleaned PDB and receptor PDBQT")
                cleaned_pdb = str(cleaned_pdb_path).replace("\\", "/")
                receptor_pdbqt = str(receptor_pdbqt_path).replace("\\", "/")
            else:
                result = automate_protein_prep.main(nolig_pdb_path)
                if result is None or not isinstance(result, tuple) or len(result) != 2:
                    logger.warning(f"Protein prep failed for {pdb_file}")
                    protein_pbar.update(1)
                    continue
                cleaned_pdb, receptor_pdbqt = result
                if cleaned_pdb is None or receptor_pdbqt is None:
                    logger.warning(f"Protein prep failed for {pdb_file}")
                    protein_pbar.update(1)
                    continue

                cleaned_pdb = os.path.abspath(cleaned_pdb).replace("\\", "/")


            os.makedirs(os.path.dirname(receptor_pdbqt), exist_ok=True)

            logger.info(f"Preparing receptor PDBQT with charges: {receptor_pdbqt}")


            pdb_id = base_id



            # Detect active site
            center, box_size = detect_active_site(cleaned_pdb)
            if center is None:
                protein_pbar.update(1)
                continue


            prep_ligands_from_pdb(
                ligand_output_dir=ligand_output_dir,
                ligands_mol2_dir=ligands_mol2_dir,
                prepped_ligands_dir=prepped_ligands_dir
            )

            # Get all prepared ligand pdbqt files for docking

            # Use original MOL2s to apply RDKit filtering
            # Collect ligands from both per-protein and global folders
            # Find ligands from both the per-protein and global folders, including the root of /prepped_ligands
            ligands_pdbqt_paths = list(Path(cfg["OUTPUT_LIGANDS_DIR"]).rglob("*.pdbqt"))

            # Filter and deduplicate ligands, logging invalid ones
            seen_stems = set()
            ligands_to_dock = []
            for p in ligands_pdbqt_paths:
                if p.stem in seen_stems:
                    continue
                if not is_valid_ligand(p, cfg["OUTPUT_LIGANDS_DIR"]):  # Still uses RDKit if possible
                    logger.warning(f"Skipping malformed ligand: {p}")
                    continue
                ligands_to_dock.append(str(p))
                seen_stems.add(p.stem)

            print(f"Total valid .pdbqt ligands to screen: {len(ligands_to_dock)}")

            # 1. Gather unique ligands from mol2 for filtering
            # Step 1: Collect all valid .pdbqt ligands
            seen_stems = set()
            valid_pdbqt_ligands = {}

            for p in ligands_pdbqt_paths:
                if p.stem in seen_stems:
                    continue
                if not is_valid_ligand(p, cfg["OUTPUT_LIGANDS_DIR"]):
                    logger.warning(f"Skipping malformed ligand: {p}")
                    continue
                valid_pdbqt_ligands[p.stem] = p  # Save for lookup
                seen_stems.add(p.stem)

            print(f"Total valid .pdbqt ligands: {len(valid_pdbqt_ligands)}")

            # Step 2: Filter using .mol2 with RDKit
            ligands_mol2_paths = list(Path(cfg["LIGANDS_MOL2_DIR"]).rglob("*.mol2"))
            ligands_to_dock = []

            for mol2_path in ligands_mol2_paths:
                stem = mol2_path.stem
                if stem not in valid_pdbqt_ligands:
                    continue  # Skip ligands that failed earlier check

                mol = Chem.MolFromMol2File(str(mol2_path), sanitize=True)
                if mol is None:
                    logger.warning(f"RDKit failed to load: {mol2_path}")
                    continue

                mw = Descriptors.MolWt(mol)
                rot_bonds = Descriptors.NumRotatableBonds(mol)

                if 50 <= mw <= 2000:
                    ligands_to_dock.append(str(valid_pdbqt_ligands[stem]))
                else:
                    logger.warning(f"Ligand {stem} outside MW limits: {mw:.2f}")

            score_history = defaultdict(dict)

            # Docking stages
            for i, stage in enumerate(stages):
                logger.info(f"\nStarting {stage['name']} with {len(ligands_to_dock)} ligands...")
                scores = {}
                validated_ligands = []
                raw_docked_ligands = {}  # key = ligand path, value = docked PDBQT path

                with ThreadPoolExecutor(max_workers=cfg["MAX_PARALLEL_JOBS"]) as executor:
                    futures = {}
                    for ligand in ligands_to_dock:
                        conf_path, out_path = generate_config(
                            cfg["OVERALL_DIR"], pdb_id, receptor_pdbqt, center, box_size,
                            ligand, stage["name"], stage,
                            max(1, cfg["MAX_PARALLEL_JOBS"] // 4)
                        )
                        logger.info(f"Docked ligand saved to: {out_path}")
                        raw_docked_ligands[ligand] = out_path
                        futures[executor.submit(run_docking_task, cfg["VINA_EXE"], conf_path, ligand, out_path)] = (ligand, out_path)

                    with tqdm(total=len(futures), desc=f"Docking Ligands ({stage['name']})", unit="ligand") as ligand_pbar:
                        for future in as_completed(futures):
                            ligand, out_path = futures[future]
                            ligand_name = os.path.basename(ligand)
                            _, score = future.result()
                            if score is not None:
                                surface_coords = extract_surface_atoms(receptor_pdbqt, center=center)

                                result = validate_pose_pdbqt(
                                    protein_pdbqt=receptor_pdbqt,
                                    ligand_pdbqt=out_path,
                                    pocket_center=center,
                                    clash_threshold=2.0,
                                    CLASH_TOLERANCE=3,
                                    DIST_THRESHOLD_SURFACE=6.0,
                                    DIST_THRESHOLD_CENTROID=4.5,
                                    surface_atom_coords=surface_coords
                                )

                                logger.info(f"{ligand_name} validation: {result}")
                                if result.get("valid", False):
                                    scores[ligand] = score
                                    score_history[stage["name"]][ligand] = score
                                    validated_ligands.append(ligand)
                                    print(f"{ligand_name} | {stage['name']} score: {score:.2f} kcal/mol (valid)")
                                else:
                                    print(f"{ligand_name} pose invalid, discarded")
                            else:
                                print(f"No score found for {ligand_name}")
                                logger.warning(f"No score found for {ligand_name}")
                            ligand_pbar.update(1)


                #active site retry
                if i < len(stages) - 1:
                    if not scores:
                        logger.warning(f"No valid ligands in {stage['name']}. Attempting fallback recentering...")

                        fallback_ligand_path, new_center, best_score = attempt_fallback_recenter(
                            fallback_ligands=raw_docked_ligands,
                            receptor_pdbqt=receptor_pdbqt,
                            docking_dir=os.path.join(cfg["DOCKED_DIR"], pdb_id),
                            stage_name=stage["name"],
                            pocket_center=center,
                            logger=logger
                        )

                        # Try to find the fallback score from score_history or infer it
                        if fallback_ligand_path is None:
                            logger.warning("Fallback recovery failed. Ending docking for this protein.")
                            break

                        # Try to find the fallback score from previous stages
                        score = None
                        for s in reversed(score_history.values()):
                            if fallback_ligand_path in s:
                                score = s[fallback_ligand_path]
                                break

                        if score is None:
                            score = best_score

                        score_history[stage["name"]][fallback_ligand_path] = score
                        validated_ligands = [fallback_ligand_path]
                        center = new_center

                        # Preserve fallback ligand as only docking candidate and skip percentage selection
                        continue  # Skip to next docking stage

                    percentages = {
                        "discovery": [1.0, 0.1, 0.01, 0.001, 0.001],
                        "polypharmacology": [1.0, 0.05, 0.005]
                    }.get(docking_mode, [1.0] * len(stages))
                    pct = percentages[i + 1] if i + 1 < len(percentages) else 0.01
                    num_to_select = int(len(scores) * pct)

                    if num_to_select < 1:
                        logger.warning(f"Percentage {pct*100:.5f}% yielded <1 ligand. Using best-scoring ligand.")
                        num_to_select = 1

                    top_ligands = sorted(scores.items(), key=lambda x: x[1])[:num_to_select]
                    ligands_to_dock = [lig_path for lig_path, _ in top_ligands]
                    print(f"Selected top {len(ligands_to_dock)} ligands ({pct * 100:.5f}%) for {stages[i + 1]['name']}")
                    logger.info(f"Selected top {len(ligands_to_dock)} ligands ({pct * 100:.5f}%) for {stages[i + 1]['name']}")

            # Final validation on last stage ligands
            logger.info("\nFinal pose validation pass on last stage ligands")
            for ligand in validated_ligands:
                last_stage = stages[-1]["name"]
                out_path = Path(
                    cfg["OVERALL_DIR"]) / "docked" / pdb_id / last_stage / f"{Path(ligand).stem}_{last_stage}.pdbqt"

                surface_coords = extract_surface_atoms(pdbqt_path=receptor_pdbqt, center=center)
                if os.path.exists(out_path):
                    result = validate_pose_pdbqt(
                        protein_pdbqt=receptor_pdbqt,
                        ligand_pdbqt=out_path,
                        pocket_center=center,
                        clash_threshold=2.0,
                        DIST_THRESHOLD_SURFACE=6.0,
                        DIST_THRESHOLD_CENTROID=4.5,
                        surface_atom_coords=surface_coords,  # ← this was missing
                    )

                    logger.info(f"Final validation for {os.path.basename(ligand)}: {result}")
                else:
                    logger.warning(f"Pose file not found for {os.path.basename(ligand)} — likely discarded earlier.")
            # === Optional Screenshot Generation ===
            import subprocess

            if validated_ligands:
                top_ligand = validated_ligands[0]  # Top pose from final stage
                last_stage = stages[-1]["name"]
                # Final filtering for polypharmacology mode
                if docking_mode == "polypharmacology":
                    final_stage = stages[-1]["name"]
                    final_scores = score_history.get(final_stage, {})
                    top_ligands = sorted(final_scores.items(), key=lambda x: x[1])[:20]
                    validated_ligands = [lig for lig, _ in top_ligands]
                    logger.info(
                        f"[Polypharmacology Mode] Selected top {len(validated_ligands)} ligands for final output.")

                docked_pose_path = Path(
                    cfg["DOCKED_DIR"]) / pdb_id / last_stage / f"{Path(top_ligand).stem}_{last_stage}.pdbqt"

                if docked_pose_path.exists():
                    image_out_prefix = Path(cfg["DOCKED_DIR"]) / pdb_id / "top_pose"
                    command = [
                        str(Path(cfg["PYMOL_PATH"])), "capture_pose.py",
                        str(cleaned_pdb),
                        str(docked_pose_path),
                        str(image_out_prefix)
                    ]
                    print("Running PyMOL visualization command:", command)
                    result = subprocess.run(command, capture_output=True, text=True)
                    print("PyMOL stdout:", result.stdout)
                    print("PyMOL stderr:", result.stderr)
                    image_out_prefix.parent.mkdir(parents=True, exist_ok=True)
                    print("Expected images:")
                    print(f"{image_out_prefix}_front.png")
                    print(f"{image_out_prefix}_side.png")
                    print(f"{image_out_prefix}_top.png")
                else:
                    logger.warning(f"No top pose found to visualize for {pdb_id}")
            # Save docking score summary csv
            protein_dock_dir = os.path.join(cfg["DOCKED_DIR"], pdb_id)
            os.makedirs(protein_dock_dir, exist_ok=True)
            csv_output_path = os.path.join(protein_dock_dir, "docking_score_summary.csv")
            write_score_summary_to_csv(score_history, output_path=csv_output_path)

            protein_pbar.update(1)

    elapsed = time.time() - overall_start
    print(f"\nAll proteins processed in {elapsed / 60:.2f} minutes.")

if __name__ == "__main__":
    main()
