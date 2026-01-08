
import os
import csv
import sys
import yaml
import shutil
import pytest
from pathlib import Path
from analysis import master_schema_export

@pytest.fixture
def mock_run_env(tmp_path):
    # Setup directories
    repo_root = tmp_path / "repo"
    repo_root.mkdir()
    
    run_id = "RUN_TEST"
    pdb = "1PDB"
    variant = "APO"
    ph = "7.4"
    
    # 1. Post Docked
    post_docked = repo_root / "post_docked" / run_id / pdb / variant / ph
    post_docked.mkdir(parents=True)
    
    # Consensus CSV
    consensus_csv = post_docked / "consensus_reranked_scorch.csv"
    with consensus_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_id", "variant", "ph_label", "ligand", "ligand_file", "library", "t_vs_decoys_consensus", "t_vs_decoys_blend", "t_vs_decoys_consensus_pre", "run_mode"])
        # Control, Valid
        writer.writerow([pdb, variant, ph, "LIG_CTRL.sanitized.pdbqt", "LIG_CTRL.sanitized.pdbqt", "LIB1", "1.0", "1.5", "0.8", "fda"])
        # FDA, Invalid
        writer.writerow([pdb, variant, ph, "LIG_FDA.pdbqt", "LIG_FDA.pdbqt", "LIB1", "2.0", "2.5", "1.8", "fda"])
        # Decoy
        writer.writerow([pdb, variant, ph, "decoys_LIG_FDA.pdbqt", "decoys_LIG_FDA.pdbqt", "LIB1", "-1.0", "-1.5", "-0.8", "dud"])

    # Posebusters
    pb_csv = post_docked / "posebusters_all_stages.csv"
    with pb_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ligand_file", "posebusters_pass", "posebusters_reason"])
        writer.writerow(["LIG_CTRL.sanitized.pdbqt", "1", "None"])
        writer.writerow(["LIG_FDA.pdbqt", "0", "Clash"])
        
    # 2. Processed PDBs (Controls)
    processed = repo_root / "processed_pdbs" / pdb / "ligands_raw"
    processed.mkdir(parents=True)
    (processed / "LIG_CTRL.sanitized.pdbqt").touch()
    
    # 3. Manifest
    manifest_path = repo_root / "post_docked" / run_id / "run_manifest.yaml"
    manifest_data = {
        "proteins": {
            pdb: {
                "stages": {
                    "pocket_detection": {
                        "details": {
                            "pocket_method": "fpocket",
                            "center_x": 10.0,
                            "center_y": 20.0,
                            "center_z": 30.0
                        }
                    }
                }
            }
        }
    }
    with manifest_path.open("w") as f:
        yaml.dump(manifest_data, f)
        
    # 4. DUD Eval Pretty
    dud_dir = repo_root / "analysis" / "dud_eval" / run_id / "post_docked"
    dud_dir.mkdir(parents=True)
    dud_file = dud_dir / f"consensus_reranked_scorch_summary_{run_id}_pretty.txt"
    with dud_file.open("w") as f:
        f.write("Header line\n")
        # variant ph runid pdb ... bedroc ef1 ef2 ef5 ef10 reason
        # We need to match the parse logic: 
        # Left: variant, ph, runid. 
        # Right: reason, ef10, ef5, ef2, ef1, bedroc?, roc_auc? 
        # The parser expects: ... BEDROC_alpha_20 EF@1% EF@2% EF@5% EF@10% Reason
        # And looks for numeric values before BEDROC for ROC.
        # Let's try to construct a line that passes:
        # APO 7.4 RUN_TEST 1PDB 0.9 0.8 20.0 5.0 4.0 3.0 2.0 OK
        # 0.9=ROC_AUC_adj, 0.8=ROC_AUC, 20.0=BEDROC, 5.0=EF1, 4.0=EF2, 3.0=EF5, 2.0=EF10
        f.write(f"{variant} {ph} {run_id} {pdb} 0.9 0.8 20.0 5.0 4.0 3.0 2.0 OK\n")

    return repo_root, run_id

def test_master_export(mock_run_env):
    repo_root, run_id = mock_run_env
    
    # Run the exporter
    # We call main directly but need to mock sys.argv
    sys.argv = ["master_schema_export.py", "--run-id", run_id, "--repo-root", str(repo_root)]
    
    assert master_schema_export.main() == 0
    
    # Verify output
    out_csv = repo_root / "data" / run_id / "master_rows.csv"
    assert out_csv.exists()
    
    with out_csv.open("r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        
    assert len(rows) == 3
    
    # Check LIG_CTRL
    ctrl = next(r for r in rows if r["ligand_base"] == "LIG_CTRL")
    assert ctrl["is_control"] == "1"
    assert ctrl["is_decoy"] == "0"
    assert ctrl["pose_valid_any"] == "1"
    assert ctrl["t_selected"] == "1.5" # Stage 2 because valid
    assert ctrl["t_selected_source"] == "stage2"
    assert ctrl["pocket_method"] == "fpocket"
    assert ctrl["center_x"] == "10.0"
    assert ctrl["ef1"] == "5.0"
    
    # Check LIG_FDA
    fda = next(r for r in rows if r["ligand_base"] == "LIG_FDA")
    assert fda["is_control"] == "0"
    assert fda["is_decoy"] == "0"
    assert fda["pose_valid_any"] == "0"
    assert fda["pose_invalid_reason_top"] == "Clash"
    assert fda["t_selected"] == "1.8" # Stage 1 because invalid (t_pre=1.8)
    assert fda["t_selected_source"] == "stage1"
    
    # Check Decoy
    decoy = next(r for r in rows if r["ligand_base"] == "decoys_LIG_FDA" and r["is_decoy"] == "1")
    assert decoy["is_decoy"] == "1"
    assert decoy["t_stage2"] == "-1.5"

    # Verify Optional Artifacts
    top_targets_csv = repo_root / "data" / run_id / "ligand_top_targets.csv"
    assert top_targets_csv.exists()
    with top_targets_csv.open("r", newline="") as f:
        top_rows = list(csv.DictReader(f))
        assert len(top_rows) == 3
        # LIG_FDA should have rank 1
        fda_top = next(r for r in top_rows if r["ligand_base"] == "LIG_FDA")
        assert fda_top["rank"] == "1"

    qc_csv = repo_root / "data" / run_id / "target_qc_summary.csv"
    assert qc_csv.exists()
    with qc_csv.open("r", newline="") as f:
        qc_rows = list(csv.DictReader(f))
        assert len(qc_rows) == 1
        assert qc_rows[0]["pdb_id"] == "1PDB"
        assert qc_rows[0]["n_rows"] == "3"
        assert qc_rows[0]["n_pose_valid"] == "1"
    
