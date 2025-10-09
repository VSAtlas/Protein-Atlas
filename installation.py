import os
from pathlib import Path

def get_default_config():
    user_home = Path.home()
    base_dir = Path(__file__).resolve().parent

    return {
        "OVERALL_DIR": str(base_dir),
        "INPUT_DIR": str(base_dir / "input_pdbs"),
        "PROTEIN_DIR": str(base_dir / "pdbqts"),
        "LIGAND_DIR": str(base_dir / "input_ligands"),
        "LIGAND_EXTRACTED_DIR": str(base_dir / "extracted_ligands"),
        "LIGANDS_MOL2_DIR": str(base_dir / "ligands_mol2"),
        "OUTPUT_LIGANDS_DIR": str(base_dir / "prepped_ligands"),
        # --- important separation ---
        "OUTPUT_DIR": str(base_dir / "processed_pdbs"),   # processed structures
        "DOCKED_DIR": str(base_dir / "docked"),           # docking results
        # ------------------------------------------------
        "PDBQT_DIR": str(base_dir / "pdbqts"),
        "PHENIX_CLEAN_SCRIPT": str(base_dir / "phenix_clean.py"),
        "VINA_PATH": str(user_home / "OneDrive/Desktop/AutoDock-Vina-1.2.7/vina_1.2.7_win.exe"),
        "PHENIX_DIR": str(user_home / "phenix-1.21.2-5419/Library/bin"),
        "PHENIX_LIB_PATH": str(user_home / "phenix-1.21.2-5419/Lib/site-packages"),
        "P2RANK_PATH": str(user_home / "OneDrive/Desktop/p2rank_2.5"),
        "MGLTOOLS_PATH": "C:/Program Files (x86)/MGLTools-1.5.7",
        "MGLTOOLS_PYTHON": "C:/Program Files (x86)/MGLTools-1.5.7/python.exe",
        "OPENBABEL_PATH": "C:/Program Files (x86)/OpenBabel-3.1.1/obabel.exe",
        "CPU_ONLY": True,
        "CPU": 28,
    }

def load_config(config_path="config.txt"):
    """
    Load key=value pairs from a config file, normalize keys to uppercase,
    and fill in defaults. Preserves OUTPUT_DIR ? processed_pdbs,
    DOCKED_DIR ? docked, etc.
    """
    config = {}
    base_dir = Path(__file__).resolve().parent

    if Path(config_path).exists():
        with open(config_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue

                key, value = line.split("=", 1)
                key = key.strip().upper()
                value = value.strip()

                if value.lower() == "true":
                    value = True
                elif value.lower() == "false":
                    value = False
                elif value.startswith("C:/Users/") or value.startswith("C:\\Users\\"):
                    # Map Windows absolute paths to relative equivalents under this repo
                    value = str(base_dir / Path(value).name)

                config[key] = value

    defaults = get_default_config()
    for key, default_value in defaults.items():
        config.setdefault(key, default_value)

    return config
