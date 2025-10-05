import os
import subprocess

def run_gnina(receptor_path, ligands_folder, center, size, output_folder, exhaustiveness=8, num_modes=3):
    """
    Runs GNINA docking on a set of ligands against a receptor.

    Args:
        receptor_path (str): Path to the receptor .pdbqt file.
        ligands_folder (str): F older containing ligand .pdbqt files.
        center (tuple): (x, y, z) center of the docking box.
        size (tuple): (x, y, z) size of the box in Å.
        output_folder (str): Where to save the docking results.
    """
    os.makedirs(output_folder, exist_ok=True)

    for ligand_file in os.listdir(ligands_folder):
        if not ligand_file.endswith(".pdbqt"):
            continue

        ligand_path = os.path.join(ligands_folder, ligand_file)
        output_path = os.path.join(output_folder, f"docked_{ligand_file}")

        cmd = [
            "gnina",
            "--receptor", receptor_path,
            "--ligand", ligand_path,
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--out", output_path,
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--seed", "0",  # Optional for reproducibility
        ]

        print(f"Running GNINA for {ligand_file}...")
        subprocess.run(cmd, check=True)

    print("Docking completed.")

