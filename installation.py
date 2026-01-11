from pathlib import Path
from input_and_export_functions import load_config as _load_config


def get_default_config():
    base_dir = Path(__file__).resolve().parent
    user_home = Path.home()
    return {
        "OVERALL_DIR": str(base_dir),
        "INPUT_DIR": str(base_dir / "input_pdbs"),
        "PROTEIN_DIR": str(base_dir / "pdbqts"),
        "LIGAND_DIR": str(base_dir / "input_ligands"),
        "LIGAND_EXTRACTED_DIR": str(base_dir / "extracted_ligands"),
        "LIGANDS_MOL2_DIR": str(base_dir / "ligands_mol2"),
        "OUTPUT_LIGANDS_DIR": str(base_dir / "prepped_ligands"),
        "OUTPUT_DIR": str(base_dir / "processed_pdbs"),
        "DOCKED_DIR": str(base_dir / "docked"),
        "PDBQT_DIR": str(base_dir / "pdbqts"),
    }


def load_config(config_path: str = "config.txt"):
    """
    Thin wrapper to the shared loader. Keeps existing call sites working.
    """
    return _load_config(config_path)
