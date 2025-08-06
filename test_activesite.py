import os
import logging
from activesite import main  # Replace with actual import if modularized
import requests

logging.basicConfig(level=logging.INFO)


def download_pdb(pdb_id, output_dir="input_pdbs"):
    os.makedirs(output_dir, exist_ok=True)
    pdb_file = os.path.join(output_dir, f"{pdb_id}.pdb")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(pdb_file, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {pdb_id} to {pdb_file}")
    except Exception as e:
        print(f"Failed to download {pdb_id}: {e}")

if __name__ == "__main__":
    for pdb_id in ["1A3N", "2HYY"]:
        download_pdb(pdb_id)

def test_ligand_extraction():
    # Replace with your actual PDB test file path
    test_pdb_path = os.path.abspath("input_pdbs/1a3n.pdb")

    # Directory for saving cleaned PDB and extracted ligands
    output_dir = "processed_pdbs_test"
    ligands_dir = os.path.join(output_dir, "ligands_removed")
    os.makedirs(ligands_dir, exist_ok=True)

    # Run main function - adjust if needed
    center, box_size = main(test_pdb_path)

    print(f"Predicted binding box center: {center}")
    print(f"Predicted binding box size: {box_size}")

    # Check ligand files created
    if os.path.exists(ligands_dir):
        ligands = os.listdir(ligands_dir)
        if ligands:
            print(f"Ligands extracted and saved ({len(ligands)} files):")
            for ligand_file in ligands:
                print(" -", ligand_file)
        else:
            print("No ligands were extracted.")
    else:
        print(f"Ligand output directory {ligands_dir} does not exist.")

if __name__ == "__main__":
    test_ligand_extraction()
