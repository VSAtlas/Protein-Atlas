from __future__ import annotations

from .path_router import make_paths

def build_paths_for_protein(cfg, base_id, pdb_file, variant=None, ph_token=None):
    p = make_paths(cfg, base_id=base_id, pdb_file=pdb_file)
    return {
        "pdb_id": base_id.upper(),
        "pdb_path": str(p.input_pdb_path),
        "nolig_pdb_path": str(p.nolig_pdb_path),
        "ligand_output_dir": str(p.ligand_output_dir),
        "ligands_mol2_dir": str(p.ligands_mol2_dir),
        "prepped_ligands_dir": str(p.prepped_ligands_dir),
        "cleaned_pdb_path": str(p.receptor_cleaned_pdb(variant)),
        "receptor_pdbqt_path": str(p.receptor_pdbqt(variant, ph_token=ph_token)),
    }


# >>> BUILD_PATHS SHIM END
# -------------------------
# Backward-compat shims
# -------------------------

