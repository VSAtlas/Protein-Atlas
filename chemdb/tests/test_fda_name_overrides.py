import csv
import subprocess
import sys
from pathlib import Path

from analysis.reporting import fda_name_map
from analysis.reporting import run_report_core as run_report

REPO_ROOT = Path(__file__).resolve().parents[2]


def _write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def _setup_repo(tmp_path, run_id, master_rows, mapping_rows):
    repo_root = tmp_path / "repo"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)

    mapping_csv = data_dir / "mapping.csv"
    _write_csv(mapping_csv, mapping_rows)

    config_path = data_dir / "config.txt"
    config_path.write_text(f"FDA_MAPPING_CSV={mapping_csv}\n", encoding="utf-8")

    master_csv = data_dir / "master_rows.csv"
    _write_csv(master_csv, master_rows)

    return repo_root


def test_fda_mapping_loads_even_with_ligand_display_present(tmp_path):
    run_id = "FDA_DISPLAY_PRESENT"
    master_rows = [
        {
            "pdb_id": "P1",
            "variant": "V1",
            "ph_label": "P1",
            "ligand_base": "rdk_1",
            "ligand_file": "rdk_1.sdf",
            "ligand_display": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "library": "fda",
            "z_selected": "1.0",
            "is_decoy": "0",
            "is_control": "0",
            "pose_valid_any": "1",
        }
    ]
    mapping_rows = [{"rdk_id": "rdk_0000001", "generic_name": "FriendlyOne"}]
    repo_root = _setup_repo(tmp_path, run_id, master_rows, mapping_rows)

    report = run_report.build_report(run_id, repo_root, top_n=5)
    top = report["targets"]["P1|V1|P1"]["top5_ligands"]
    assert top[0]["ligand_display"] == "FriendlyOne"


def test_unfriendly_display_overridden_but_controls_preserved(tmp_path):
    run_id = "FDA_CONTROL_PRESERVE"
    master_rows = [
        {
            "pdb_id": "P1",
            "variant": "V1",
            "ph_label": "P1",
            "ligand_base": "rdk_1",
            "ligand_file": "rdk_1.sdf",
            "ligand_display": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
            "library": "fda",
            "z_selected": "1.0",
            "is_decoy": "0",
            "is_control": "0",
            "pose_valid_any": "1",
        },
        {
            "pdb_id": "P1",
            "variant": "V1",
            "ph_label": "P1",
            "ligand_base": "rdk_2",
            "ligand_file": "rdk_2.sdf",
            "ligand_display": "rdk_2",
            "library": "fda",
            "z_selected": "0.5",
            "is_decoy": "0",
            "is_control": "1",
            "pose_valid_any": "1",
        },
    ]
    mapping_rows = [
        {"rdk_id": "rdk_0000001", "generic_name": "FriendlyOne"},
        {"rdk_id": "rdk_0000002", "generic_name": "FriendlyControl"},
    ]
    repo_root = _setup_repo(tmp_path, run_id, master_rows, mapping_rows)

    report = run_report.build_report(run_id, repo_root, top_n=5)
    top = report["targets"]["P1|V1|P1"]["top5_ligands"]
    by_base = {entry["ligand_base"]: entry["ligand_display"] for entry in top}
    assert by_base["rdk_1"] == "FriendlyOne"
    assert by_base["rdk_2"] == "rdk_2"


def test_mapping_fallback_columns_generic_or_brand(tmp_path):
    mapping_csv = tmp_path / "mapping.csv"
    mapping_rows = [{"rdk_id": "rdk_0000003", "brand_names": "BrandOnly"}]
    _write_csv(mapping_csv, mapping_rows)

    fda_index = fda_name_map.try_load_fda_index(mapping_csv)
    assert fda_index is not None
    resolved = fda_name_map.resolve_ligand_display_name("rdk_3", "", fda_index)
    assert resolved == "BrandOnly"


def test_run_report_help_has_no_rdkit_numpy_import_noise():
    proc = subprocess.run(
        [sys.executable, "-m", "analysis.cli.run_report", "-h"],
        cwd=str(REPO_ROOT),
        text=True,
        capture_output=True,
        check=False,
    )
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    assert proc.returncode == 0
    assert "_ARRAY_API not found" not in combined
    assert "No module named 'prep_ligands'" not in combined
    assert "A module that was compiled using NumPy 1.x" not in combined
    assert "Generate run report YAML" in combined
