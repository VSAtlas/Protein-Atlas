import os
import gzip
from pathlib import Path
from input_and_export_functions import load_config, validate_config
ligand_dir = Path(cfg.get("LIGAND_DIR", cfg.get("ligand_dir", ""))).resolve()
extracted_dir = Path(
    cfg.get(
        "EXTRACTED_LIGANDS_DIR",
        cfg.get("LIGAND_EXTRACTED_DIR", ligand_dir.parent / "extracted_ligands")
    )
).resolve()

def decompress_sdf_gz_files():
    cfg = load_config("config.txt")
    validate_config(cfg)

    # Prefer uppercase; fall back to legacy lowercase if present
    ligand_dir_str = cfg.get("LIGAND_DIR", cfg.get("ligand_dir"))
    if not ligand_dir_str:
        raise KeyError("Missing LIGAND_DIR in config.txt (or legacy ligand_dir).")
    ligand_dir = Path(ligand_dir_str).expanduser().resolve()

    # Destination for extracted SDFs (configurable)
    output_dir = Path(
        cfg.get("EXTRACTED_LIGANDS_DIR", str(ligand_dir.parent / "extracted_ligands"))
    ).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for file in ligand_dir.glob("*.sdf.gz"):
        output_file = output_dir / file.with_suffix('').name
        with gzip.open(file, "rb") as f_in, open(output_file, "wb") as f_out:
            f_out.write(f_in.read())
        print(f"Decompressed: {file.name} -> {output_file.name}")

if __name__ == "__main__":
    decompress_sdf_gz_files()
