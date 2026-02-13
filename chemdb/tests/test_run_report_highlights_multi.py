import csv

from analysis import run_report


def _write_master_rows(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def test_report_highlights_imatinib(tmp_path):
    repo_root = tmp_path / "repo"
    run_id = "HIGHLIGHT_TEST"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)

    rows = []
    for i in range(15):
        ligand_base = f"LIG_{i + 1}"
        ligand_display = f"Ligand {i + 1}"
        if i == 14:
            ligand_base = "LIG_IMATINIB"
            ligand_display = "Imatinib"
        rows.append(
            {
                "pdb_id": "T1",
                "variant": "V1",
                "ph_label": "P1",
                "ligand_base": ligand_base,
                "ligand_display": ligand_display,
                "t_selected": str(100 - i),
                "t_selected_source": "stage1",
                "is_decoy": "0",
                "is_control": "0",
            }
        )

    _write_master_rows(data_dir / "master_rows.csv", rows)

    report = run_report.build_report(run_id, repo_root)
    target = report["targets"]["T1|V1|P1"]
    highlights = target.get("highlights", [])
    hit = next(h for h in highlights if h["query"] == "imatinib")

    assert hit["found"] is True
    assert hit["ligand_display"] == "Imatinib"
    assert hit["rank"] == 15
    assert abs(hit["pct_rank"] - 1.0) < 1e-6


def test_report_multi_target_hits_capped_and_sorted(tmp_path):
    repo_root = tmp_path / "repo"
    run_id = "MULTI_TARGET_TEST"
    data_dir = repo_root / "data" / run_id
    data_dir.mkdir(parents=True)

    rows = []
    for target_id in ("T1", "T2"):
        rows.append(
            {
                "pdb_id": target_id,
                "variant": "V1",
                "ph_label": "P1",
                "ligand_base": "LIGX",
                "ligand_display": "Ligand X",
                "t_selected": "100.0",
                "t_selected_source": "stage2",
                "is_decoy": "0",
                "is_control": "0",
            }
        )
        rows.append(
            {
                "pdb_id": target_id,
                "variant": "V1",
                "ph_label": "P1",
                "ligand_base": "LIGY",
                "ligand_display": "Ligand Y",
                "t_selected": "99.0",
                "t_selected_source": "stage2",
                "is_decoy": "0",
                "is_control": "0",
            }
        )
        for i in range(38):
            rows.append(
                {
                    "pdb_id": target_id,
                    "variant": "V1",
                    "ph_label": "P1",
                    "ligand_base": f"{target_id}_LIG_{i}",
                    "ligand_display": f"{target_id} Ligand {i}",
                    "t_selected": str(98 - i),
                    "t_selected_source": "stage1",
                    "is_decoy": "0",
                    "is_control": "0",
                }
            )

    _write_master_rows(data_dir / "master_rows.csv", rows)

    report = run_report.build_report(
        run_id,
        repo_root,
        multi_target_max_pct=0.05,
        multi_target_max_hits=1,
    )
    hits = report.get("multi_target_hits", [])

    assert len(hits) == 1
    hit = hits[0]
    assert hit["ligand_base"] == "LIGX"
    assert hit["targets_qualified"] == 2
    assert len(hit["targets"]) == 2
