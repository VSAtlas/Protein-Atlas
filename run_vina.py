import subprocess
from input_and_export_functions import extract_best_score
from pathlib import Path

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
