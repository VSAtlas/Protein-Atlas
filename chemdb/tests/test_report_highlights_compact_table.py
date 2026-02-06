import csv
from pathlib import Path

from analysis import run_report


def _write_pdb(path: Path, molecule_name: str) -> None:
    lines = [
        "HEADER    TEST TARGET\n",
        "COMPND    MOL_ID: 1;\n",
        f"COMPND   2 MOLECULE: {molecule_name};\n",
        "ATOM      1  N   GLY A   1      11.104  13.207   8.668  1.00 20.00           N\n",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="utf-8")


def _write_master_rows(path: Path) -> None:
    ligands = ["LIGA", "LIGB", "LIGC", "LIGD", "LIGE", "LIGF", "LIGG", "LIGH", "LIGI", "LIGJ"]
    targets = ["PDB1", "PDB2", "PDB3"]
    rows = []
    for target in targets:
        for idx, ligand in enumerate(ligands):
            t_selected: str = str(80 - idx)
            if ligand == "LIGA":
                t_selected = (
                    "100" if target == "PDB1" else "95" if target == "PDB2" else "93"
                )
            elif ligand == "LIGB":
                if target == "PDB1":
                    t_selected = "90"
                else:
                    t_selected = ""
            rows.append(
                {
                    "run_id": "testrun",
                    "pdb_id": target,
                    "variant": "V1",
                    "ph_label": "7.4",
                    "ligand_base": ligand,
                    "ligand": f"{ligand}.pdbqt",
                    "ligand_file": f"{ligand}.pdbqt",
                    "ligand_display": ligand,
                    "library": "FDA",
                    "is_decoy": "false",
                    "is_control": "false",
                    "pose_valid_any": "true",
                    "t_selected": t_selected,
                    "t_selected_source": "consensus",
                }
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def test_report_highlights_compact_table(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    run_id = "testrun"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True, exist_ok=True)

    _write_pdb(repo_root / "input_pdbs" / "PDB1.pdb", "Alpha kinase")
    _write_pdb(repo_root / "input_pdbs" / "PDB2.pdb", "Beta kinase")
    _write_pdb(repo_root / "input_pdbs" / "PDB3.pdb", "Gamma kinase")
    _write_master_rows(data_dir / "master_rows.csv")

    report = run_report.build_report(
        run_id=run_id,
        repo_root=repo_root,
        highlight_ligands=None,
        highlight_max_per_target=0,
        top_targets_n=2,
        breadth_pct_threshold=0.15,
        pct_display_decimals=1,
        tail_median_decimals=0,
        highlights_top_pct=1.0,
        highlights_max_ligands=100,
    )

    highlights = report.get("ligand_highlights") or []
    assert highlights, "Expected compact ligand highlights in report output"

    liga = next((row for row in highlights if row.get("ligand_base") == "LIGA"), None)
    assert liga is not None
    assert {"ligand", "best_target", "best_score", "best_percentile", "top_targets", "coverage", "breadth", "tail_summary"} <= set(liga.keys())
    assert "Alpha kinase (PDB1)" in str(liga.get("best_target"))
    assert str(liga.get("top_targets")).startswith("1)")
    assert "2)" in str(liga.get("top_targets"))
    assert "3)" not in str(liga.get("top_targets"))
    assert liga.get("coverage") == "3 / 3 targets"
    assert liga.get("breadth") == "3"
    assert liga.get("tail_summary") == "+1 more (median p=10%)"

    ligb = next((row for row in highlights if row.get("ligand_base") == "LIGB"), None)
    assert ligb is not None
    assert ligb.get("coverage") == "1 / 3 targets"
    assert ligb.get("breadth") == "0"
