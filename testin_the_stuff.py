import os
import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import csv
from pathlib import Path

from input_and_export_functions import load_inputs
from pose_validation import validate_pose_pdbqt
import automate_protein_prep
from protein_functions import convert_to_pdbqt
from run_vina import run_docking_task  # Assuming this exists and calls the actual docking process
from protein_functions import detect_active_site
from input_and_export_functions import generate_config
import logging
from input_and_export_functions import (
    load_inputs, define_docking_stages, write_score_summary_to_csv,
    extract_best_score, generate_config
)
from run_vina import run_docking_task

# Test if the cleaned PDB exists and is correctly formatted
def test_cleaned_pdb_exists(pdb_id):
    path = Path(f"processed_pdbs/{pdb_id}/{pdb_id}_cleaned.pdb")
    assert path.exists(), f"Missing cleaned PDB at {path}"
    with open(path) as f:
        assert not any("HOH" in line for line in f), "Found water residues in cleaned PDB"
    print(f"[✓] Valid cleaned PDB: {path}")


# Test if ligands were extracted properly
def test_ligand_extraction(pdb_id):
    ligand_dir = Path(f"processed_pdbs/{pdb_id}/{pdb_id}_cleaned_ligands")
    ligands = list(ligand_dir.glob("*.pdb"))
    assert len(ligands) > 0, f"No ligands extracted in {ligand_dir}"
    print(f"[✓] Ligands extracted: {len(ligands)} ligands")


# Test ligand preparation (mol2 and pdbqt files)
def test_ligand_preparation(pdb_id):
    mol2_dir = Path(f"processed_pdbs/{pdb_id}/{pdb_id}_mol2_ligands")
    pdbqt_dir = Path(f"processed_pdbs/{pdb_id}/{pdb_id}_prepped_ligands")
    mol2 = list(mol2_dir.glob("*.mol2"))
    pdbqt = list(pdbqt_dir.glob("*.pdbqt"))
    assert len(mol2) > 0, f"No MOL2 ligands in {mol2_dir}"
    assert len(pdbqt) > 0, f"No prepped ligands in {pdbqt_dir}"
    print(f"[✓] Ligands prepared: {len(pdbqt)} ligands")


# Test if docking results are generated for each protein
def test_docking_outputs(pdb_id):
    stage1_dir = Path(f"docked/{pdb_id}/stage1")
    results = list(stage1_dir.glob("*.pdbqt"))
    assert len(results) > 0, f"No docking results in {stage1_dir}"
    print(f"[✓] Docking results found: {len(results)}")


# Test if docking score summary CSV is created
def test_score_csv(pdb_id):
    summary_path = Path(f"docked/{pdb_id}/docking_score_summary.csv")
    assert summary_path.exists(), f"Missing score summary at {summary_path}"
    with open(summary_path, mode="r") as file:
        reader = csv.reader(file)
        rows = list(reader)
        assert len(rows) > 1, f"{summary_path} is empty or has no valid data."
        print(f"[✓] Docking score CSV valid with {len(rows) - 1} entries")


# Check if the configuration file is valid
def test_config(cfg):
    assert isinstance(cfg, dict), "Config should be a dictionary"
    required_keys = ["input_dir", "output_dir", "pdbqt_dir", "vina_exe"]
    for key in required_keys:
        assert key in cfg, f"Missing {key} in config"
    print(f"[✓] Configuration valid: {list(cfg.keys())}")


def test_pdbqt_file_validity(pdb_id, stage="stage1"):
    stage_dir = Path(f"docked/{pdb_id}/{stage}")
    pdbqt_files = list(stage_dir.glob("*.pdbqt"))
    assert len(pdbqt_files) > 0, f"No PDBQT files in {stage_dir}"

    for pdbqt in pdbqt_files:
        with open(pdbqt, "r") as file:
            lines = file.readlines()
            assert len(lines) > 0, f"{pdbqt} is empty."
            assert any("ATOM" in line for line in lines), f"{pdbqt} does not contain any atom data."
    print(f"[✓] Valid PDBQT files in {stage_dir}")


def test_ligand_validity(pdb_id):
    ligand_dir = Path(f"processed_pdbs/{pdb_id}/{pdb_id}_prepped_ligands")
    ligands = list(ligand_dir.glob("*.pdbqt"))
    for ligand in ligands:
        with open(ligand, "r") as file:
            lines = file.readlines()
            assert len(lines) > 0, f"{ligand} is empty."
            assert any("ATOM" in line for line in lines), f"{ligand} does not contain atom data."
    print(f"[✓] Valid ligands: {len(ligands)} ligands")


def test_docking_results_validity(pdb_id):
    stage1_dir = Path(f"docked/{pdb_id}/stage1")
    docking_results = list(stage1_dir.glob("*.pdbqt"))
    assert len(docking_results) > 0, f"No docking results found in {stage1_dir}"

    for result in docking_results:
        with open(result, "r") as file:
            lines = file.readlines()
            assert len(lines) > 0, f"{result} is empty."
            assert any("REMARK" in line for line in lines), f"{result} does not contain REMARK."
            assert any("VINA" in line for line in lines), f"{result} does not contain VINA."
    print(f"[✓] Valid docking results in {stage1_dir}")


def test_docking_score_csv(pdb_id):
    score_csv_path = Path(f"docked/{pdb_id}/docking_score_summary.csv")
    assert score_csv_path.exists(), f"Missing score summary at {score_csv_path}"
    with open(score_csv_path, mode="r") as file:
        reader = csv.reader(file)
        rows = list(reader)
        assert len(rows) > 1, f"{score_csv_path} is empty or has no valid data."
        for row in rows[1:]:  # Skip the header
            score = row[1]
            try:
                assert float(score) != 0, f"Invalid docking score {score} in {score_csv_path}"
            except ValueError:
                raise AssertionError(f"Non-numeric score {score} in {score_csv_path}")
    print(f"[✓] Valid score entries in {score_csv_path}")
def test_missing_loop_completion(pdb_id):
    filtered_path = Path(f"processed_pdbs/{pdb_id}/{pdb_id}_filtered.pdb")
    cleaned_path = Path(f"processed_pdbs/{pdb_id}/{pdb_id}_cleaned.pdb")

    assert filtered_path.exists(), f"Filtered PDB not found: {filtered_path}"
    assert cleaned_path.exists(), f"Cleaned PDB not found: {cleaned_path}"

    def count_atoms(path):
        with open(path, "r") as f:
            return sum(1 for line in f if line.startswith("ATOM"))

    filtered_atoms = count_atoms(filtered_path)
    cleaned_atoms = count_atoms(cleaned_path)

    assert cleaned_atoms >= filtered_atoms, (
        f"Cleaned PDB has fewer atoms ({cleaned_atoms}) than filtered ({filtered_atoms})"
    )

    print(f"[✓] Loop/side-chain completion: {cleaned_atoms} atoms (was {filtered_atoms})")
def test_prepare_receptor_output(pdb_id):
    pdbqt_path = Path(f"{pdb_id}_receptor.pdbqt")
    assert pdbqt_path.exists(), f"Receptor PDBQT not found at {pdbqt_path}"

    with open(pdbqt_path, "r") as file:
        contents = file.read()
        assert "ROOT" in contents and "ENDROOT" in contents, "Missing ROOT/ENDROOT in receptor file"
        assert "ATOM" in contents, "No ATOM records in receptor PDBQT"

    print(f"[✓] Valid receptor PDBQT: {pdbqt_path}")


def test_assign_protonation_states(input_pdb: str, output_pdb: str):
    """
    Tests whether assign_protonation_states successfully adds hydrogens
    and produces a non-empty valid output file.
    """
    output_path = Path(output_pdb)

    # Check that output file was created
    assert output_path.exists(), f"Protonated file not created: {output_pdb}"

    with open(output_path, "r") as f:
        lines = f.readlines()

    assert len(lines) > 0, f"Output file is empty: {output_pdb}"
    assert any(" H" in line[12:16] for line in lines if line.startswith("ATOM")), \
        f"No hydrogen atoms detected in {output_pdb}"

    print(f"[✓] Protonation successful: {output_pdb} contains {sum(' H' in line[12:16] for line in lines if line.startswith('ATOM'))} hydrogen atoms")

# Ensure cleaned_pdb is correctly assigned
def run_all_tests(pdb_id, cfg):
    logger = logging.getLogger(pdb_id)
    logger.setLevel(logging.DEBUG)
    print(f"\n--- Testing pipeline for PDB: {pdb_id} ---")

    # Retrieve mgltools_python and prepare_script from cfg
    mgltools_python = cfg.get("MGLTOOLS_PYTHON")
    prepare_script = cfg.get("prepare_receptor_script")

    # Debug: Print the paths to ensure they are correct
    print(f"mgltools_python: {mgltools_python}")
    print(f"prepare_script: {prepare_script}")

    # Check if both paths are valid
    if not mgltools_python or not os.path.exists(mgltools_python):
        raise ValueError("Invalid mgltools_python path.")
    if not os.path.exists(prepare_script):
        raise ValueError(f"prepare_receptor_script not found at {prepare_script}")

    # Ensure the cleaned pdb file exists
    cleaned_pdb_path = os.path.join(cfg["input_dir"], f"{pdb_id}.pdb").replace("\\", "/")

    if not os.path.exists(cleaned_pdb_path):
        # If cleaned pdb is not found, generate it
        print(f"Cleaned PDB not found. Preparing {pdb_id}...")
        cleaned_pdb = automate_protein_prep.main(pdb_id)  # Call the function that prepares the cleaned PDB
        if cleaned_pdb is None:
            raise FileNotFoundError(f"Protein preparation failed for {pdb_id}")
        print(f"Cleaned PDB generated: {cleaned_pdb}")
        cleaned_pdb_path = cleaned_pdb  # Now set cleaned_pdb_path to the generated cleaned pdb path
    else:
        print(f"Found cleaned PDB: {cleaned_pdb_path}")

    # Step 1: Convert PDB to PDBQT (for receptor)
    receptor_pdbqt = f"{pdb_id}_receptor.pdbqt"  # Output file path for PDBQT
    success = convert_to_pdbqt(cleaned_pdb_path, receptor_pdbqt, mgltools_python, prepare_script)

    if not success:
        raise AssertionError("Failed to convert PDB to PDBQT.")

    # Step 2: Run the docking (adjust based on your pipeline)
    stages = define_docking_stages()  # Define docking stages
    # Get prepared ligands for docking from both directories
    prepped_ligands_dir = Path(cfg["output_ligands_dir"])
    control_ligands_dir = Path(f"processed_pdbs/{pdb_id}/{pdb_id}_prepped_ligands")

    # List all PDBQT ligands from both directories
    pdbqt_ligands = list(prepped_ligands_dir.glob("*.pdbqt"))
    control_ligands = list(control_ligands_dir.glob("*.pdbqt"))

    ligands_to_dock = [str(p) for p in pdbqt_ligands]  # Ensure ligands are correctly loaded
    ligands_to_dock.extend([str(p) for p in control_ligands])  # Add control ligands to the docking list

    # Detect active site (center and box size)
    center, box_size = detect_active_site(cleaned_pdb_path)  # Ensure this is done before using center and box_size

    # Loop through the ligands to dock and set config_path and out_path for each
    score_history = defaultdict(dict)  # Initialize score history for each stage
    for stage in stages:
        stage_name = stage["name"]  # Use stage["name"] for accessing the stage's name
        print(f"Starting docking for stage: {stage_name} with {len(ligands_to_dock)} ligands...")

        scores = {}
        validated_ligands = []

        with ThreadPoolExecutor(max_workers=cfg["MAX_PARALLEL_JOBS"]) as executor:
            futures = {}
            for ligand in ligands_to_dock:
                conf_path, out_path = generate_config(
                    cfg["overall_dir"], pdb_id, receptor_pdbqt, center, box_size,
                    ligand, stage_name, stage,  # Corrected here
                    max(1, cfg["MAX_PARALLEL_JOBS"] // 4)
                )
                print(f"Docked ligand saved to: {out_path}")
                futures[executor.submit(run_docking_task, cfg["vina_exe"], conf_path, ligand, out_path)] = (ligand, out_path)

            for future in as_completed(futures):
                ligand, out_path = futures[future]
                ligand_name = os.path.basename(ligand)
                _, score = future.result()
                if score is not None:
                    scores[ligand] = score
                    score_history[stage_name][ligand] = score
                    validated_ligands.append(ligand)
                    print(f"{ligand_name} | {stage_name} score: {score:.2f} kcal/mol")
                else:
                    print(f"No score found for {ligand_name}")
                    validated_ligands.append(ligand)

    # Save docking score summary CSV
    protein_dock_dir = os.path.join(cfg["docked_dir"], pdb_id)
    os.makedirs(protein_dock_dir, exist_ok=True)
    csv_output_path = os.path.join(protein_dock_dir, "docking_score_summary.csv")
    write_score_summary_to_csv(score_history, output_path=csv_output_path)

    # Run all tests
    test_cleaned_pdb_exists(pdb_id)
    test_ligand_extraction(pdb_id)
    test_ligand_preparation(pdb_id)
    test_docking_outputs(pdb_id)
    test_score_csv(pdb_id)
    test_config(cfg)

    test_pdbqt_file_validity(pdb_id)
    test_ligand_validity(pdb_id)
    test_docking_results_validity(pdb_id)
    test_docking_score_csv(pdb_id)
    test_prepare_receptor_output(pdb_id)
    test_missing_loop_completion(pdb_id)
    print(f"--- All tests passed for {pdb_id} ---\n")

def run_all_tests_no_docking(pdb_id, cfg):
    logger = logging.getLogger(pdb_id)
    logger.setLevel(logging.DEBUG)
    print(f"\n--- Testing pipeline for PDB: {pdb_id} (no docking) ---")

    # Check MGLTools paths
    mgltools_python = cfg.get("MGLTOOLS_PYTHON")
    prepare_script = cfg.get("prepare_receptor_script")

    if not mgltools_python or not os.path.exists(mgltools_python):
        raise ValueError("Invalid mgltools_python path.")
    if not os.path.exists(prepare_script):
        raise ValueError(f"prepare_receptor_script not found at {prepare_script}")

    # Ensure cleaned pdb exists or prepare it
    cleaned_pdb_path = os.path.join(cfg["input_dir"], f"{pdb_id}.pdb").replace("\\", "/")
    if not os.path.exists(cleaned_pdb_path):
        print(f"Cleaned PDB not found. Preparing {pdb_id}...")
        cleaned_pdb = automate_protein_prep.main(pdb_id)
        if cleaned_pdb is None:
            raise FileNotFoundError(f"Protein preparation failed for {pdb_id}")
        cleaned_pdb_path = cleaned_pdb
    else:
        print(f"Found cleaned PDB: {cleaned_pdb_path}")

    # Convert PDB to PDBQT (receptor only)
    receptor_pdbqt = f"{pdb_id}_receptor.pdbqt"
    success = convert_to_pdbqt(cleaned_pdb_path, receptor_pdbqt, mgltools_python, prepare_script)
    if not success:
        raise AssertionError("Failed to convert PDB to PDBQT.")

    # Run tests (skip all docking-related ones)
    test_cleaned_pdb_exists(pdb_id)
    test_ligand_extraction(pdb_id)
    test_ligand_preparation(pdb_id)
    test_config(cfg)
    test_pdbqt_file_validity(pdb_id)
    test_ligand_validity(pdb_id)
    test_prepare_receptor_output(pdb_id)
    test_missing_loop_completion(pdb_id)

    print(f"--- Non-docking tests passed for {pdb_id} ---\n")

# Example of running the tests for a given protein (pdb_id)
if __name__ == "__main__":
    cfg = load_inputs()  # Load your config
    pdb_ids = ["6T5U"]  # Add your PDB IDs here for testing
    for pdb_id in pdb_ids:
        run_all_tests(pdb_id, cfg)
