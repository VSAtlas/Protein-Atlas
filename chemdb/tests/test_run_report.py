import csv
import sys
import yaml
import pytest
from pathlib import Path
from analysis import run_report

@pytest.fixture
def mock_report_env(tmp_path):
    repo_root = tmp_path / "repo"
    run_id = "REPORT_TEST"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)
    
    # 1. Master CSV
    # Columns: pdb_id, variant, ph_label, ligand_base, is_decoy, is_control, t_selected, t_selected_source, ef1, roc_auc
    # Added: consensus_mu_decoy, consensus_sigma_decoy, consensus_n_decoys, blend_mu_decoy, blend_sigma_decoy, blend_n_decoys
    rows = [
        # Target 1: T1|V1|P1
        {
            "pdb_id": "T1", "variant": "V1", "ph_label": "P1", "ligand_base": "L1", "t_selected": "2.0", "t_selected_source": "stage2", "is_decoy": "0", "is_control": "0", "ef1": "10.0",
            "consensus_mu_decoy": "0.5", "consensus_sigma_decoy": "0.1", "consensus_n_decoys": "100",
            "blend_mu_decoy": "0.6", "blend_sigma_decoy": "0.2", "blend_n_decoys": "100"
        },
        {
            "pdb_id": "T1", "variant": "V1", "ph_label": "P1", "ligand_base": "L2", "t_selected": "1.5", "t_selected_source": "stage2", "is_decoy": "0", "is_control": "0", "ef1": "10.0",
            "consensus_mu_decoy": "0.5", "consensus_sigma_decoy": "0.1", "consensus_n_decoys": "100",
            "blend_mu_decoy": "0.6", "blend_sigma_decoy": "0.2", "blend_n_decoys": "100"
        },
        {
            "pdb_id": "T1", "variant": "V1", "ph_label": "P1", "ligand_base": "L3", "t_selected": "1.0", "t_selected_source": "stage1", "is_decoy": "0", "is_control": "0", "ef1": "10.0",
            "consensus_mu_decoy": "0.5", "consensus_sigma_decoy": "0.1", "consensus_n_decoys": "100",
            "blend_mu_decoy": "0.6", "blend_sigma_decoy": "0.2", "blend_n_decoys": "100"
        },
        {
            "pdb_id": "T1", "variant": "V1", "ph_label": "P1", "ligand_base": "L4", "t_selected": "0.5", "t_selected_source": "stage1", "is_decoy": "0", "is_control": "0", "ef1": "10.0",
            "consensus_mu_decoy": "0.5", "consensus_sigma_decoy": "0.1", "consensus_n_decoys": "100",
            "blend_mu_decoy": "0.6", "blend_sigma_decoy": "0.2", "blend_n_decoys": "100"
        },
        {
            "pdb_id": "T1", "variant": "V1", "ph_label": "P1", "ligand_base": "L5", "t_selected": "0.0", "t_selected_source": "stage1", "is_decoy": "0", "is_control": "0", "ef1": "10.0",
            "consensus_mu_decoy": "0.5", "consensus_sigma_decoy": "0.1", "consensus_n_decoys": "100",
            "blend_mu_decoy": "0.6", "blend_sigma_decoy": "0.2", "blend_n_decoys": "100"
        },
        {
            "pdb_id": "T1", "variant": "V1", "ph_label": "P1", "ligand_base": "L6", "t_selected": "-1.0", "t_selected_source": "stage1", "is_decoy": "0", "is_control": "0", "ef1": "10.0",
            "consensus_mu_decoy": "0.5", "consensus_sigma_decoy": "0.1", "consensus_n_decoys": "100",
            "blend_mu_decoy": "0.6", "blend_sigma_decoy": "0.2", "blend_n_decoys": "100"
        },
        {
            "pdb_id": "T1", "variant": "V1", "ph_label": "P1", "ligand_base": "D1", "t_selected": "-2.0", "t_selected_source": "stage1", "is_decoy": "1", "is_control": "0", "ef1": "10.0",
            "consensus_mu_decoy": "0.5", "consensus_sigma_decoy": "0.1", "consensus_n_decoys": "100",
            "blend_mu_decoy": "0.6", "blend_sigma_decoy": "0.2", "blend_n_decoys": "100"
        },
        
        # Target 2: T2|V1|P1
        {
            "pdb_id": "T2", "variant": "V1", "ph_label": "P1", "ligand_base": "L1", "t_selected": "3.0", "t_selected_source": "stage2", "is_decoy": "0", "is_control": "0", "ef1": "5.0",
            "consensus_mu_decoy": "0.4", "consensus_sigma_decoy": "0.05", "consensus_n_decoys": "50",
            "blend_mu_decoy": "nan", "blend_sigma_decoy": "", "blend_n_decoys": "0"
        },
        {
            "pdb_id": "T2", "variant": "V1", "ph_label": "P1", "ligand_base": "D2", "t_selected": "0.0", "t_selected_source": "stage1", "is_decoy": "1", "is_control": "0", "ef1": "5.0",
            "consensus_mu_decoy": "0.4", "consensus_sigma_decoy": "0.05", "consensus_n_decoys": "50",
            "blend_mu_decoy": "nan", "blend_sigma_decoy": "", "blend_n_decoys": "0"
        },
        {
            "pdb_id": "T2", "variant": "V1", "ph_label": "P1", "ligand_base": "D3", "t_selected": "0.0", "t_selected_source": "stage1", "is_decoy": "1", "is_control": "0", "ef1": "5.0",
            "consensus_mu_decoy": "0.4", "consensus_sigma_decoy": "0.05", "consensus_n_decoys": "50",
            "blend_mu_decoy": "nan", "blend_sigma_decoy": "", "blend_n_decoys": "0"
        },
    ]
    
    with (data_dir / "master_rows.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # 2. Mock SCORCH CSV for Target 1
    scorch_dir = repo_root / "post_docked" / run_id / "T1" / "V1" / "P1"
    scorch_dir.mkdir(parents=True)
    scorch_rows = [
        {"scorch_mu_decoy": "0.7", "scorch_sigma_decoy": "0.3", "scorch_n_decoys": "120"}
    ]
    with (scorch_dir / "scorch_scores_all.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=scorch_rows[0].keys())
        writer.writeheader()
        writer.writerows(scorch_rows)
        
    return repo_root, run_id

def test_report_generation(mock_report_env):
    repo_root, run_id = mock_report_env
    
    # Run
    report = run_report.build_report(run_id, repo_root, top_n=5)
    
    # Validate structure
    assert report["run_id"] == run_id
    assert "summary" in report
    assert "targets" in report
    assert "ligands" in report
    
    # Validate Targets
    t1_key = "T1|V1|P1"
    assert t1_key in report["targets"]
    t1 = report["targets"][t1_key]
    
    # QC
    assert t1["qc"]["ef1"] == 10.0
    assert t1["qc"]["n_decoys"] == 1
    assert t1["qc"]["n_rows"] == 7
    
    # Decoy Stats T1
    ds1 = t1["qc"]["decoy_stats"]
    assert ds1["consensus"]["mu"] == 0.5
    assert ds1["consensus"]["sigma"] == 0.1
    assert ds1["consensus"]["n"] == 100
    assert ds1["blend"]["mu"] == 0.6
    assert ds1["blend"]["sigma"] == 0.2
    assert ds1["blend"]["n"] == 100
    assert ds1["scorch"]["mu"] == 0.7
    assert ds1["scorch"]["sigma"] == 0.3
    assert ds1["scorch"]["n"] == 120
    
    # Top 5
    top5 = t1["top5_ligands"]
    assert len(top5) == 5
    assert top5[0]["ligand_base"] == "L1"
    assert top5[0]["t_selected"] == 2.0
    assert top5[4]["ligand_base"] == "L5"
    
    # Check exclusion of L6 (rank 6)
    bases = [x["ligand_base"] for x in top5]
    assert "L6" not in bases
    
    # Target 2
    t2_key = "T2|V1|P1"
    t2 = report["targets"][t2_key]
    assert t2["qc"]["n_decoys"] == 2
    assert len(t2["top5_ligands"]) == 3 # Only 3 rows
    
    # Decoy Stats T2
    ds2 = t2["qc"]["decoy_stats"]
    assert ds2["consensus"]["mu"] == 0.4
    assert ds2["consensus"]["n"] == 50
    # Blend was empty/nan
    assert ds2["blend"]["mu"] is None
    assert ds2["blend"]["sigma"] is None
    # Scorch file doesn't exist for T2
    assert ds2["scorch"]["mu"] is None
    
    # Validate Ligands Section
    # L1 should be present (in top 5 of T1 and T2)
    assert "L1" in report["ligands"]
    l1_hits = report["ligands"]["L1"]["targets_in_top5"]
    assert len(l1_hits) == 2
    assert l1_hits[0]["target"] == "T1|V1|P1"
    assert l1_hits[0]["t_selected"] == 2.0
    assert l1_hits[1]["target"] == "T2|V1|P1" 
    assert l1_hits[1]["t_selected"] == 3.0
    
    # L6 should NOT be present (not in any top 5)
    assert "L6" not in report["ligands"]
    
    # Validate Exclusions
    # "alerts" should not be at top level
    assert "alerts" not in report
    
    # Ligand summaries shouldn't have forbidden keys
    # Actually ligand summary structure is simplified to just targets_in_top5 list
    assert "n_variants_seen" not in report["ligands"]["L1"]
    
    # Top hit entries exclusions
    hit = top5[0]
    assert "best_engine" not in hit
    assert "final_rank" not in hit
    
    # Global n_controls exclusion
    assert "n_controls" not in report["summary"]

