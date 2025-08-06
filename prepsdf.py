import os
import gzip
from pathlib import Path

def read_config(path="config.txt"):
    config = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):  # skip empty lines and comments
                key, value = line.split("=", 1)
                config[key.strip()] = value.strip()
    return config

def decompress_sdf_gz_files():
    config = read_config()

    ligand_dir = Path(config["ligand_dir"]).expanduser().resolve()

    # Define output folder for extracted sdf files:
    # Using the same parent directory as ligand_dir, but a different subfolder
    output_dir = ligand_dir.parent / "extracted_ligands"
    output_dir.mkdir(parents=True, exist_ok=True)

    for file in ligand_dir.glob("*.sdf.gz"):
        output_file = output_dir / file.with_suffix('').name  # Remove .gz only, keep .sdf name
        with gzip.open(file, "rb") as f_in, open(output_file, "wb") as f_out:
            f_out.write(f_in.read())
        print(f"Decompressed: {file.name} -> {output_file.name}")

if __name__ == "__main__":
    decompress_sdf_gz_files()
