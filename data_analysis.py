import os
import csv
import re
from pathlib import Path
import pandas as pd

def extract_fda_scores(docked_dir):
    FDA_LIGAND_PATTERN = re.compile(r"(fda_\d+)(?:_stage\d+)?\.pdbqt", re.IGNORECASE)
    ligand_protein_scores = []

    for protein_dir in Path(docked_dir).iterdir():
        if not protein_dir.is_dir():
            continue

        pdb_id = protein_dir.name.upper()
        summary_csv = protein_dir / "docking_score_summary.csv"
        if not summary_csv.exists():
            continue

        with open(summary_csv, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                ligand_file = row.get("Ligand", "").strip()
                print(f"Checking ligand: {ligand_file}")  # <--- DEBUG

                match = FDA_LIGAND_PATTERN.search(ligand_file)
                if not match:
                    print("  ⛔ No match for FDA ligand pattern")  # <--- DEBUG
                    continue

                ligand_id = match.group(1).lower()
                print(f"  ✅ Matched FDA ligand: {ligand_id}")  # <--- DEBUG

                scores = []
                for key, val in row.items():
                    if key.lower().startswith("stage") and val and not val.endswith(".pdbqt"):
                        print(f"    → Score from {key}: {val}")  # <--- DEBUG
                        try:
                            score = float(val)
                            scores.append(score)
                        except ValueError:
                            print(f"      ⚠️ Could not convert score: {val}")  # <--- DEBUG

                if scores:
                    best_score = min(scores)
                    print(f"  🏆 Best score: {best_score}")  # <--- DEBUG
                    ligand_protein_scores.append((ligand_id, pdb_id, best_score))

    df = pd.DataFrame(ligand_protein_scores, columns=["Ligand ID", "Protein PDB ID", "Best Score"])
    df.sort_values(["Ligand ID", "Best Score"], ascending=[True, True], inplace=True)
    return df

if __name__ == "__main__":
    DOCKED_DIR = "E:/PythonProject/protein_automation/docked"
    df = extract_fda_scores(DOCKED_DIR)
    output_path = Path(DOCKED_DIR) / "fda_ligand_protein_scores.csv"
    df.to_csv(output_path, index=False)
    print(f"[✓] Scores saved to {output_path}")
