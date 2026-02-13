import gzip
from pathlib import Path

from input_and_export_functions import load_config, validate_config


def decompress_sdf_gz_files():
    cfg = load_config()
    validate_config(cfg)

    ligand_dir_str = cfg.get("LIGAND_DIR") or cfg.get("ligand_dir")
    if not ligand_dir_str:
        raise KeyError("Missing LIGAND_DIR in config.txt (or legacy ligand_dir).")
    ligand_dir = Path(ligand_dir_str).expanduser().resolve()

    output_dir_str = (
        cfg.get("EXTRACTED_LIGANDS_DIR") or str(ligand_dir.parent / "extracted_ligands")
    )
    output_dir = Path(output_dir_str).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for file in ligand_dir.glob("*.sdf.gz"):
        output_file = output_dir / file.with_suffix("").name
        with gzip.open(file, "rb") as f_in, open(output_file, "wb") as f_out:
            f_out.write(f_in.read())
        print(f"Decompressed: {file.name} -> {output_file.name}")


if __name__ == "__main__":
    decompress_sdf_gz_files()
