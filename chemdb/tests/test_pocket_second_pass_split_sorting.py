from pathlib import Path
import csv
import json

from analysis.pocket_second_pass import run_second_pass


def _write_scores_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ligand_id", "score", "label"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_performance_json(path: Path, pockets):
    payload = {
        "pdb_id": "6LU7",
        "mode": "pocket_eval",
        "pockets": pockets,
        "selected_pocket_id": "P1",
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def test_pocket_second_pass_split_sorting(tmp_path):
    eval_dir = tmp_path / "docked" / "RUNTEST" / "6LU7" / "pocket_eval"

    pockets = []
    for pocket_id in ("P1", "P2"):
        csv_path = eval_dir / pocket_id / "calibration_scores.csv"
        rows = []
        for idx in range(1, 11):
            if idx <= 5:
                score = -9.0 if pocket_id == "P1" else -4.0
            else:
                score = -4.0 if pocket_id == "P1" else -9.0
            rows.append({"ligand_id": f"S{idx}", "score": score, "label": "strong"})
        for idx in range(1, 11):
            rows.append({"ligand_id": f"N{idx}", "score": -6.0, "label": "non"})
        _write_scores_csv(csv_path, rows)
        pockets.append(
            {
                "pocket_id": pocket_id,
                "metrics": {"auc_mean": 0.50, "ef1_mean": 0.50, "auc_std": 0.0},
                "artifacts": {"scores_path": str(csv_path)},
            }
        )

    perf_path = eval_dir / "pocket_performance.json"
    _write_performance_json(perf_path, pockets)

    result = run_second_pass(perf_path)

    assert result["multi_pocket_detected"] is True
    assert result["selected_pockets"] == ["P1", "P2"]

    assignments_path = Path(result["paths"]["assignments"])
    performance_path = Path(result["paths"]["performance_sorted"])
    assert assignments_path.exists()
    assert performance_path.exists()

    assignments_payload = json.loads(assignments_path.read_text(encoding="utf-8"))
    assert assignments_payload["multi_pocket_detected"] is True
    assert assignments_payload["selected_pockets"] == ["P1", "P2"]

    assignments = assignments_payload.get("assignments", [])
    ambiguous = assignments_payload.get("ambiguous", [])
    assert all(entry.get("assigned_pocket_id") is None for entry in ambiguous)

    total_strong = 10
    counts = {"P1": 0, "P2": 0}
    for entry in assignments:
        assigned = entry["assigned_pocket_id"]
        assert assigned in counts
        counts[assigned] += 1

    assert counts["P1"] >= 0.30 * total_strong
    assert counts["P2"] >= 0.30 * total_strong

    performance_payload = json.loads(performance_path.read_text(encoding="utf-8"))
    assert performance_payload["multi_pocket_detected"] is True
    assert performance_payload["selected_pockets"] == ["P1", "P2"]

    by_id = {entry["pocket_id"]: entry for entry in performance_payload["pockets"]}
    for pocket_id in ("P1", "P2"):
        unsorted_auc = by_id[pocket_id]["unsorted_metrics"]["auc"]
        sorted_auc = by_id[pocket_id]["sorted_metrics"]["auc"]
        assert sorted_auc >= unsorted_auc + 0.10
