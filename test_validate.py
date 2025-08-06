from pose_validation import validate_pose_pdbqt

# Use minimal dummy files you've created
protein_file = "test_protein.pdbqt"
ligand_file = "test_ligand.pdbqt"
pocket_center = [0.0, 0.0, 0.0]  # Adjust if needed

result = validate_pose_pdbqt(
    protein_pdbqt=protein_file,
    ligand_pdbqt=ligand_file,
    pocket_center=pocket_center,
    clash_threshold=2.0,
    dist_threshold=8.0,
    ligand_smiles="CCO"  # Ethanol
)

print("Validation result:")
for k, v in result.items():
    print(f"{k}: {v}")
