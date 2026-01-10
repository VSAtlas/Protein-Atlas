import csv
import sys
from pathlib import Path

import pytest

from analysis import master_schema_export


def _write_rows(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_target_fdr(monkeypatch, tmp_path):
    repo_root = tmp_path / "repo"
    repo_root.mkdir()
    run_id = "RUN_FDR"
    pdb_id = "2PDB"
    variant = "APO"
    ph7 = "pH7_0"
    ph6 = "pH6_0"

    # Lower threshold for faster test
    monkeypatch.setattr(master_schema_export, "MIN_DECOYS_FOR_FDR", 20)

    consensus_fields = [
        "pdb_id",
        "variant",
        "ph_label",
        "ligand",
        "ligand_file",
        "library",
        "consensus_score",
        "consensus_score_pre",
        "run_mode",
    ]

    # Target 1: uses consensus_score
    consensus_rows_ph7 = [
        {
            "pdb_id": pdb_id,
            "variant": variant,
            "ph_label": ph7,
            "ligand": "LIG1.pdbqt",
            "ligand_file": "LIG1.pdbqt",
            "library": "FDA",
            "consensus_score": "21.0",
            "consensus_score_pre": "20.0",
            "run_mode": "fda",
        },
        {
            "pdb_id": pdb_id,
            "variant": variant,
            "ph_label": ph7,
            "ligand": "LIG2.pdbqt",
            "ligand_file": "LIG2.pdbqt",
            "library": "FDA",
            "consensus_score": "19.0",
            "consensus_score_pre": "18.0",
            "run_mode": "fda",
        },
        {
            "pdb_id": pdb_id,
            "variant": variant,
            "ph_label": ph7,
            "ligand": "decoys_mock1.pdbqt",
            "ligand_file": "decoys_mock1.pdbqt",
            "library": "DECOY",
            "consensus_score": "5.0",
            "consensus_score_pre": "4.0",
            "run_mode": "dud",
        },
    ]
    ph7_dir = repo_root / "post_docked" / run_id / pdb_id / variant / ph7
    _write_rows(ph7_dir / "consensus_reranked_scorch.csv", consensus_fields, consensus_rows_ph7)

    dud_fields = [
        "pdb_id",
        "variant",
        "ph_label",
        "ligand_file",
        "library",
        "consensus_score",
        "consensus_score_pre",
    ]
    dud_rows_ph7 = [
        {
            "pdb_id": pdb_id,
            "variant": variant,
            "ph_label": ph7,
            "ligand_file": f"decoy_{i}.pdbqt",
            "library": "DECOY",
            "consensus_score": str(i + 1),
            "consensus_score_pre": str(i),
        }
        for i in range(20)
    ]
    _write_rows(ph7_dir / "dud_consensus_reranked_scorch.csv", dud_fields, dud_rows_ph7)

    # Target 2: forces fallback to consensus_score_pre
    consensus_rows_ph6 = [
        {
            "pdb_id": pdb_id,
            "variant": variant,
            "ph_label": ph6,
            "ligand": "LIG_PRE.pdbqt",
            "ligand_file": "LIG_PRE.pdbqt",
            "library": "FDA",
            "consensus_score": "",
            "consensus_score_pre": "9.0",
            "run_mode": "fda",
        },
    ]
    ph6_dir = repo_root / "post_docked" / run_id / pdb_id / variant / ph6
    _write_rows(ph6_dir / "consensus_reranked_scorch.csv", consensus_fields, consensus_rows_ph6)

    dud_rows_ph6 = [
        {
            "pdb_id": pdb_id,
            "variant": variant,
            "ph_label": ph6,
            "ligand_file": f"pre_decoy_{i}.pdbqt",
            "library": "DECOY",
            "consensus_score": "",
            "consensus_score_pre": str(i),
        }
        for i in range(20)
    ]
    _write_rows(ph6_dir / "dud_consensus_reranked_scorch.csv", dud_fields, dud_rows_ph6)

    sys.argv = [
        "master_schema_export.py",
        "--run-id",
        run_id,
        "--repo-root",
        str(repo_root),
        "--overwrite",
    ]
    assert master_schema_export.main() == 0

    out_csv = repo_root / "data" / run_id / "master_rows.csv"
    with out_csv.open("r", newline="") as f:
        rows = {r["ligand_base"]: r for r in csv.DictReader(f)}

    # Combo using consensus_score
    lig1 = rows["LIG1"]
    lig2 = rows["LIG2"]
    decoy = rows["decoys_mock1"]

    assert lig1["fdr_score_field"] == "consensus_score"
    assert lig1["fdr_n_decoys"] == "20"
    assert float(lig1["fdr_p_empirical"]) == pytest.approx(1 / 21)
    assert float(lig1["fdr_q_target"]) == pytest.approx(0.095238, rel=1e-5)
    assert lig1["fdr_hit_q10"] == "1"
    assert lig1["fdr_hit_q05"] == "0"

    assert float(lig2["fdr_p_empirical"]) == pytest.approx(3 / 21)
    assert float(lig2["fdr_q_target"]) == pytest.approx(0.142857, rel=1e-5)
    assert lig2["fdr_hit_q10"] == "0"
    assert lig2["fdr_hit_q05"] == "0"

    # Decoy rows should stay blank
    assert decoy["fdr_score_field"] == ""
    assert decoy["fdr_p_empirical"] == ""
    assert decoy["fdr_q_target"] == ""
    assert decoy["fdr_hit_q05"] == ""
    assert decoy["fdr_hit_q10"] == ""

    # Combo using consensus_score_pre
    lig_pre = rows["LIG_PRE"]
    assert lig_pre["fdr_score_field"] == "consensus_score_pre"
    assert lig_pre["fdr_n_decoys"] == "20"
    assert float(lig_pre["fdr_p_empirical"]) == pytest.approx(12 / 21)
    assert float(lig_pre["fdr_q_target"]) == pytest.approx(12 / 21)

    summary_csv = repo_root / "data" / run_id / "target_fdr_summary.csv"
    assert summary_csv.exists()
    with summary_csv.open("r", newline="") as f:
        summary_rows = list(csv.DictReader(f))
    assert len(summary_rows) == 2
    summary_map = {(r["pdb_id"], r["variant"], r["ph_label"]): r for r in summary_rows}
    assert summary_map[(pdb_id, variant, ph7)]["fdr_score_field"] == "consensus_score"
    assert summary_map[(pdb_id, variant, ph6)]["fdr_score_field"] == "consensus_score_pre"
