# dry_run_paths.py
from pathlib import Path
from input_and_export_functions import load_config
from path_router import Paths
import os

cfg = load_config()
def make(cfg, pdb_id="1T46", pdb_file="1T46.pdb"):
    return Paths(
        pdb_id=pdb_id,
        pdb_file=pdb_file,
        over_root=Path(cfg["OVERALL_DIR"]),
        input_root=Path(cfg["INPUT_DIR"]),
        processed_root=Path(cfg["OUTPUT_DIR"]),
        docked_root=Path(cfg["DOCKED_DIR"]),
        prepped_root=Path(cfg.get("PREPPED_LIGANDS_DIR", cfg.get("OUTPUT_LIGANDS_DIR", cfg.get("PREPPED_LIGANDS_ROOT", "")))),
        ligands_mol2_root=Path(cfg["LIGANDS_MOL2_DIR"]),
        configs_root=Path(cfg.get("CONFIGS_DIR", Path(cfg["OVERALL_DIR"]) / "configs")),
        root_pdb_dir=Path(cfg["OUTPUT_DIR"]) / pdb_id,
        raw_dir=Path(cfg["OUTPUT_DIR"]) / pdb_id / "raw",
        ligands_raw_dir=Path(cfg["OUTPUT_DIR"]) / pdb_id / "ligands_raw",
        nolig_dir=Path(cfg["OUTPUT_DIR"]) / pdb_id / "nolig",
        work_dir=Path(cfg["OUTPUT_DIR"]) / pdb_id / "work",
        input_pdb_path=Path(cfg["INPUT_DIR"]) / pdb_file,
        nolig_pdb_path=Path(cfg["OUTPUT_DIR"]) / pdb_id / "nolig" / f"{pdb_id}_nolig.pdb",
    )

paths = make(cfg)

def show(tag, p: Path): print(f"{tag}: {p}")

print("\n=== Legacy (no variant, no pH) ===")
show("receptor_dir", paths.receptor_dir(None))                             # processed_pdbs/1T46/receptor/
show("cleaned_pdb", paths.receptor_cleaned_pdb(None))                      # .../1T46_cleaned.pdb
show("receptor_pdbqt", paths.receptor_pdbqt(None))                         # .../1T46.pdbqt
show("ligands_raw_dir", paths.ligand_output_dir)                           # processed_pdbs/1T46/ligands_raw/
show("prepped_lig_dir", paths.prepped_ligands_dir)                         # prepped_ligands/1T46/
show("ligands_mol2_dir", paths.ligands_mol2_dir)                           # ligands_mol2/1T46/
print("(Docked CSVs) docking_score_summary.csv expected in docked/1T46/")

print("\n=== HOLO, no pH ===")
show("receptor_dir", paths.receptor_dir("HOLO"))                            # processed_pdbs/1T46/HOLO/receptor/
show("cleaned_pdb", paths.receptor_cleaned_pdb("HOLO"))                     # .../1T46_cleaned.pdb
show("receptor_pdbqt", paths.receptor_pdbqt("HOLO"))                        # .../1T46.pdbqt

print("\n=== APO, with pH token 'pH8_0+8_1+8_2' ===")
show("receptor_dir", paths.receptor_dir("APO"))                             # processed_pdbs/1T46/APO/receptor/
show("receptor_pdbqt", paths.receptor_pdbqt("APO","pH8_0+8_1+8_2"))         # .../1T46_pH8_0+8_1+8_2.pdbqt
