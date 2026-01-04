from pathlib import Path
import sys
import csv
import logging
import json

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from record_data import write_scores_csv
from rescoring_scorch import _load_consensus_top_bases, discover_stage3_roots
from rescore_reranker import rerank_consensus_with_scorch

RUN_ID = "test_run"


def _base_cfg(root: Path) -> dict:
    root.mkdir(parents=True, exist_ok=True)
    cfg = {
        "OVERALL_DIR": str(root),
        "INPUT_DIR": str(root / "input_pdbs"),
        "OUTPUT_DIR": str(root / "processed_pdbs"),
        "DOCKED_DIR": str(root / "docked"),
        "OUTPUT_LIGANDS_DIR": str(root / "prepped_ligands"),
        "USE_GNINA": "false",
        "RUN_ID": RUN_ID,
    }
    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR", "OUTPUT_LIGANDS_DIR"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)
    return cfg


def _score_history() -> dict:
    return {
        "stage1": {
            "ligand1.pdbqt": {
                "score": -7.5,
                "valid": True,
                "reason": "",
                "heavy_atoms": 12,
                "le": None,
                "pains_flag": False,
            }
        }
    }


def test_write_scores_csv_with_prefix(tmp_path):
    cfg = _base_cfg(tmp_path / "prefixed")
    csv_path = write_scores_csv(cfg, "TEST1", _score_history(), csv_prefix="dud_")
    csv_file = Path(csv_path)

    assert csv_file.name == "dud_docking_score_summary.csv"
    assert csv_file.exists()

    long_path = csv_file.with_name("dud_docking_score_long.csv")
    assert long_path.exists()
    expected_parent = Path(cfg["DOCKED_DIR"]) / RUN_ID / "TEST1"
    assert csv_file.parent.resolve() == expected_parent.resolve()


def test_write_scores_csv_default_prefix(tmp_path):
    cfg = _base_cfg(tmp_path / "default")
    csv_path = write_scores_csv(cfg, "TEST2", _score_history())
    csv_file = Path(csv_path)

    assert csv_file.name == "docking_score_summary.csv"
    assert csv_file.exists()

    long_path = csv_file.with_name("docking_score_long.csv")
    assert long_path.exists()
    expected_parent = Path(cfg["DOCKED_DIR"]) / RUN_ID / "TEST2"
    assert csv_file.parent.resolve() == expected_parent.resolve()


def test_discover_stage3_roots(tmp_path):
    variant_root = tmp_path / "variant"
    fda_root = variant_root / "stage3"
    dud_root = variant_root / "dud_stage3"
    fda_root.mkdir(parents=True)
    dud_root.mkdir(parents=True)

    roots = discover_stage3_roots(variant_root)
    assert roots["fda"] == fda_root
    assert roots["dud"] == dud_root


def test_stratified_top_fraction_includes_decoys(tmp_path):
    consensus_path = tmp_path / "consensus_docking_scores.csv"
    rows = [
        {"run_id": "r1", "pdb_id": "P1", "variant": "APO", "ph_label": "", "ligand": "fda_a_stage1.pdbqt", "consensus_score": "0.9", "library": "FDA"},
        {"run_id": "r1", "pdb_id": "P1", "variant": "APO", "ph_label": "", "ligand": "fda_b_stage1.pdbqt", "consensus_score": "0.8", "library": "FDA"},
        {"run_id": "r1", "pdb_id": "P1", "variant": "APO", "ph_label": "", "ligand": "fda_c_stage1.pdbqt", "consensus_score": "0.7", "library": "FDA"},
        {"run_id": "r1", "pdb_id": "P1", "variant": "APO", "ph_label": "", "ligand": "dud_a_stage1.pdbqt", "consensus_score": "0.6", "library": "DECOY"},
        {"run_id": "r1", "pdb_id": "P1", "variant": "APO", "ph_label": "", "ligand": "dud_b_stage1.pdbqt", "consensus_score": "0.5", "library": "DECOY"},
        {"run_id": "r1", "pdb_id": "P1", "variant": "APO", "ph_label": "", "ligand": "dud_c_stage1.pdbqt", "consensus_score": "0.4", "library": "DECOY"},
    ]
    with consensus_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    allowed_fda, allowed_dud, total_rows, k_selected, non_controls, _, _ = _load_consensus_top_bases(
        consensus_path, 0.5, logging.getLogger("test_top_frac")
    )
    assert total_rows == len(rows)
    assert non_controls == len(rows)
    assert k_selected >= 2  # stratified selection yields at least one per bucket
    assert any(base.startswith("fda_") for base in allowed_fda)
    assert any(base.startswith("dud_") for base in allowed_dud)


def test_rerank_outputs_t_scores(tmp_path):
    logger = logging.getLogger("test_rerank")
    consensus_path = tmp_path / "consensus_docking_scores.csv"
    scorch_path = tmp_path / "scorch_scores_all.csv"
    cnn_path = scorch_path.parent / "cnn_rescoring.csv"

    cons_rows = [
        {
            "run_id": "r2",
            "pdb_id": "P2",
            "variant": "HOLO",
            "ph_label": "phys",
            "ligand": "fda_a.pdbqt",
            "consensus_score": "0.9",
            "library": "FDA",
            "t_vs_decoys_consensus": "0.1",
        },
        {
            "run_id": "r2",
            "pdb_id": "P2",
            "variant": "HOLO",
            "ph_label": "phys",
            "ligand": "fda_b.pdbqt",
            "consensus_score": "0.7",
            "library": "FDA",
            "t_vs_decoys_consensus": "-0.2",
        },
        {
            "run_id": "r2",
            "pdb_id": "P2",
            "variant": "HOLO",
            "ph_label": "phys",
            "ligand": "dud_a.pdbqt",
            "consensus_score": "0.5",
            "library": "DECOY",
            "t_vs_decoys_consensus": "0.0",
        },
        {
            "run_id": "r2",
            "pdb_id": "P2",
            "variant": "HOLO",
            "ph_label": "phys",
            "ligand": "dud_b.pdbqt",
            "consensus_score": "0.3",
            "library": "DECOY",
            "t_vs_decoys_consensus": "-0.1",
        },
    ]
    with consensus_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cons_rows[0].keys()))
        writer.writeheader()
        writer.writerows(cons_rows)

    scorch_rows = [
        {"pdb_id": "P2", "variant": "HOLO", "ph": "phys", "source": "vina", "Ligand_ID": "fda_a_stage1", "SCORCH_score": "2.0", "SCORCH_certainty": "1.0", "run_mode": "fda"},
        {"pdb_id": "P2", "variant": "HOLO", "ph": "phys", "source": "vina", "Ligand_ID": "fda_b_stage1", "SCORCH_score": "1.5", "SCORCH_certainty": "1.0", "run_mode": "fda"},
    ]
    with scorch_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scorch_rows[0].keys()))
        writer.writeheader()
        writer.writerows(scorch_rows)

    dud_scorch_rows = [
        {"pdb_id": "P2", "variant": "HOLO", "ph": "phys", "source": "vina", "Ligand_ID": "dud_a_stage1", "SCORCH_score": "1.0", "SCORCH_certainty": "1.0", "run_mode": "dud"},
        {"pdb_id": "P2", "variant": "HOLO", "ph": "phys", "source": "vina", "Ligand_ID": "dud_b_stage1", "SCORCH_score": "0.5", "SCORCH_certainty": "1.0", "run_mode": "dud"},
    ]
    dud_scorch_path = scorch_path.parent / "dud_scorch_scores_all.csv"
    with dud_scorch_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dud_scorch_rows[0].keys()))
        writer.writeheader()
        writer.writerows(dud_scorch_rows)

    cnn_rows = [
        {"ligand": "fda_a_stage1", "cnn_score": "1.0", "cnn_affinity": "1.0", "stage_priority": "stage1"},
        {"ligand": "dud_a_stage1", "cnn_score": "0.8", "cnn_affinity": "1.0", "stage_priority": "stage1"},
    ]
    cnn_path.parent.mkdir(parents=True, exist_ok=True)
    with cnn_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cnn_rows[0].keys()))
        writer.writeheader()
        writer.writerows(cnn_rows)

    out_csv = scorch_path.parent / "consensus_reranked_scorch.csv"
    ok = rerank_consensus_with_scorch(
        consensus_path,
        scorch_path,
        out_csv,
        logger,
        overwrite=True,
    )
    assert ok
    with out_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        assert "t_vs_decoys_consensus" in fieldnames
        assert "t_vs_decoys_blend" in fieldnames
        rows = list(reader)
    assert any(r.get("t_vs_decoys_blend") not in ("", None) and r.get("library") == "FDA" for r in rows)
    stats_path = out_csv.parent / "decoy_stats_blend.json"
    assert stats_path.exists()
    stats = json.loads(stats_path.read_text())
    assert any(entry.get("n_decoys", 0) > 0 for entry in stats.values())
