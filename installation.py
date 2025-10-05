import os
from pathlib import Path

def get_default_config():
    user_home = Path.home()
    base_dir = Path(__file__).resolve().parent

    return {
        "overall_dir": str(base_dir),
        "input_dir": str(base_dir / "input_pdbs"),
        "protein_dir": str(base_dir / "pdbqts"),
        "ligand_dir": str(base_dir / "input_ligands"),
        "ligand_extracted_dir": str(base_dir / "extracted_ligands"),
        "ligands_mol2_dir": str(base_dir / "ligands_mol2"),
        "output_ligands_dir": str(base_dir / "prepped_ligands"),
        "output_dir": str(base_dir / "processed_pdbs"),
        "pdbqt_dir": str(base_dir / "pdbqts"),
        "docked_dir": str(base_dir / "docked"),
        "phenix_clean_script": str(base_dir / "phenix_clean.py"),
        "vina_path": str(user_home / "OneDrive/Desktop/AutoDock-Vina-1.2.7/vina_1.2.7_win.exe"),
        "phenix_dir": str(user_home / "phenix-1.21.2-5419/Library/bin"),
        "phenix_lib_path": str(user_home / "phenix-1.21.2-5419/Lib/site-packages"),
        "p2rank_path": str(user_home / "OneDrive/Desktop/p2rank_2.5"),
        "mgltools_path": "C:/Program Files (x86)/MGLTools-1.5.7",
        "mgltools_python": "C:/Program Files (x86)/MGLTools-1.5.7/python.exe",
        "openbabel_path": "C:/Program Files (x86)/OpenBabel-3.1.1/obabel.exe",
        "cpu_only": True,
        "cpu":28
    }

def load_config(config_path="config.txt"):
    config = {}
    base_dir = Path(__file__).resolve().parent

    if Path(config_path).exists():
        with open(config_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                key, value = line.split("=", 1)
                key, value = key.strip(), value.strip()

                if value.lower() == "true":
                    value = True
                elif value.lower() == "false":
                    value = False
                elif value.startswith("C:/Users/") or value.startswith("C:\\Users\\"):
                    value = str(base_dir / Path(value).name)

                config[key] = value

    # Fill in missing keys with default values
    defaults = get_default_config()
    for key, default_value in defaults.items():
        if key not in config:
            config[key] = default_value

    return config
