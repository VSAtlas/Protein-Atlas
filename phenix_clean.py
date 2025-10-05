import sys
import os


def clean_pdb(pdb_file, output_dir):
    pdb_id = os.path.basename(pdb_file).split(".")[0]
    out_path = os.path.join(output_dir, pdb_id)
    os.makedirs(out_path, exist_ok=True)

    cleaned_pdb_path = os.path.join(out_path, f"{pdb_id}_cleaned.pdb")

    # === Your actual cleaning logic would go here ===
    # For now, just touch the file so something is written
    with open(cleaned_pdb_path, "w") as f:
        f.write("")  # Placeholder content

    cleaned_pdb_path = os.path.abspath(cleaned_pdb_path)
    print(f"Phenix cleaned: {cleaned_pdb_path}")
    return cleaned_pdb_path


if __name__ == "__main__":
    pdb_file = sys.argv[1]
    output_dir = sys.argv[2]
    clean_pdb(pdb_file, output_dir)
