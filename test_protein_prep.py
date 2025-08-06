import os
import sys
import logging
import urllib.request
import re
from automate_protein_prep import main  # <-- Replace protein_prep with your actual module filename

logging.basicConfig(level=logging.INFO)

PDB_ID = "1a3n"  # Change to any PDB you want to test
INPUT_DIR = "./input_pdbs"
OUTPUT_DIR = "./processed_pdbs_test"

def download_pdb(pdb_id, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    pdb_path = os.path.join(out_dir, f"{pdb_id}.pdb")
    if os.path.isfile(pdb_path):
        logging.info(f"PDB file already exists: {pdb_path}")
        return pdb_path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    logging.info(f"Downloading PDB {pdb_id} from {url} ...")
    try:
        urllib.request.urlretrieve(url, pdb_path)
        logging.info(f"Downloaded to {pdb_path}")
        return pdb_path
    except Exception as e:
        logging.error(f"Failed to download PDB: {e}")
        return None

def count_hydrogens(pdb_path):
    count = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip().startswith("H"):
                count += 1
    return count

def parse_molprobity_log(log_path):
    results = {}
    if not os.path.isfile(log_path):
        logging.warning(f"MolProbity log not found: {log_path}")
        return results

    with open(log_path) as f:
        text = f.read()

    # Example regex for clashscore: "Clashscore, all atoms: 3.45"
    clashscore_match = re.search(r"Clashscore, all atoms:\s*([\d\.]+)", text)
    if clashscore_match:
        results['clashscore'] = float(clashscore_match.group(1))

    # Example regex for rotamer outliers: "Rotamer outliers: 1.23%"
    rotamer_match = re.search(r"Rotamer outliers:\s*([\d\.]+)%", text)
    if rotamer_match:
        results['rotamer_outliers_%'] = float(rotamer_match.group(1))

    return results

def test_protein_prep():
    pdb_path = download_pdb(PDB_ID, INPUT_DIR)
    if not pdb_path:
        print("Failed to get PDB file. Exiting test.")
        return

    cleaned_pdb = main(f"{PDB_ID}.pdb", output_dir=OUTPUT_DIR)
    if not cleaned_pdb or not os.path.isfile(cleaned_pdb):
        print("Protein prep failed or cleaned PDB not found.")
        return

    print(f"Protein prep successful. Cleaned PDB: {cleaned_pdb}")

    num_h = count_hydrogens(cleaned_pdb)
    print(f"Number of hydrogens in cleaned PDB: {num_h}")

    molprobity_log = "C:\\Users\\micha\\PycharmProjects\\PythonProject\\protein_automation\\processed_pdbs_test\\1a3n\\1a3n_molprobity.log"
    mp_results = parse_molprobity_log(molprobity_log)
    if mp_results:
        print(f"MolProbity results:\n  Clashscore: {mp_results.get('clashscore', 'N/A')}\n"
              f"  Rotamer outliers (%): {mp_results.get('rotamer_outliers_%', 'N/A')}")
    else:
        print("No MolProbity results found or failed to parse.")

if __name__ == "__main__":
    test_protein_prep()
