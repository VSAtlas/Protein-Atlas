import sys
import os
from pathlib import Path
from input_and_export_functions import load_config, validate_config
_cfg = load_config("config.txt")
validate_config(_cfg)


def clean_pdb(pdb_file, output_dir):
    pdb_id = os.path.basename(pdb_file).split(".")[0]

    # allow config to set a default outdir if caller passes "" or None
    if not output_dir:
        output_dir = _cfg.get("PHENIX_CLEAN_OUTDIR", os.path.join(_cfg["OUTPUT_DIR"], "phenix_clean"))

    out_path = os.path.join(output_dir, pdb_id)
    os.makedirs(out_path, exist_ok=True)

    suffix = _cfg.get("PHENIX_CLEAN_SUFFIX", "_cleaned.pdb")  
    cleaned_pdb_path = os.path.join(out_path, f"{pdb_id}{suffix}")

    # ===  actual cleaning logic would go here, currently placeholder ===
    # For now, just touch the file so something is written
    with open(cleaned_pdb_path, "w") as f:
        f.write("")  # Placeholder content

    cleaned_pdb_path = os.path.abspath(cleaned_pdb_path)
    print(f"Phenix cleaned: {cleaned_pdb_path}")
    return cleaned_pdb_path


if __name__ == "__main__":
    pdb_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else _cfg.get("PHENIX_CLEAN_OUTDIR", os.path.join(_cfg["OUTPUT_DIR"], "phenix_clean"))
    clean_pdb(pdb_file, output_dir)
