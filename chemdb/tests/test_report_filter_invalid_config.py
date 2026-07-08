import csv
import subprocess
import sys
from pathlib import Path

import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parents[2]


def _write_csv(path: Path, fieldnames, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _setup_repo(tmp_path: Path, run_id: str, filter_invalid: bool) -> Path:
    repo_root = tmp_path / "repo"
    repo_root.mkdir(parents=True, exist_ok=True)

    config_lines = [
        f"FILTER_INVALID={'true' if filter_invalid else 'false'}",
        "TEST_MODE_ENABLE=fda",
    ]
    _write_text(repo_root / "config.txt", "\n".join(config_lines) + "\n")

    _write_text(
        repo_root / "manifests" / run_id / "run_manifest.yaml",
        "{}\n",
    )

    combo_dir = (
        repo_root / "post_docked" / run_id / "1T46" / "HOLO" / "pH8_5"
    )

    _write_csv(
        combo_dir / "consensus_reranked_scorch.csv",
        [
            "ligand",
            "z_vs_decoys_blend",
            "z_vs_decoys_consensus",
            "z_vs_decoys_consensus_pre",
            "library",
            "run_mode",
            "final_rank",
        ],
        [
            {
                "ligand": "rdk_0002911.pdbqt",
                "z_vs_decoys_blend": "1.96319",
                "z_vs_decoys_consensus": "",
                "z_vs_decoys_consensus_pre": "",
                "library": "UNKNOWN",
                "run_mode": "fda",
                "final_rank": "17",
            },
            {
                "ligand": "other_ligand.pdbqt",
                "z_vs_decoys_blend": "1.5",
                "z_vs_decoys_consensus": "1.5",
                "z_vs_decoys_consensus_pre": "1.5",
                "library": "UNKNOWN",
                "run_mode": "fda",
                "final_rank": "18",
            },
        ],
    )

    _write_csv(
        combo_dir / "posebusters_all_stages.csv",
        ["ligand_file", "posebusters_pass", "posebusters_reason"],
        [
            {
                "ligand_file": "rdk_0002911.pdbqt",
                "posebusters_pass": "FALSE",
                "posebusters_reason": "too_far_from_pocket",
            },
            {
                "ligand_file": "other_ligand.pdbqt",
                "posebusters_pass": "TRUE",
                "posebusters_reason": "",
            },
        ],
    )

    return repo_root


def _run_master_export(repo_root: Path, run_id: str) -> None:
    script = REPO_ROOT / "analysis" / "master_schema_export.py"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "--run-id",
            run_id,
            "--repo-root",
            str(repo_root),
            "--overwrite",
        ],
        check=True,
    )


def _run_report(repo_root: Path, run_id: str) -> None:
    script = REPO_ROOT / "analysis" / "run_report.py"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "--run-id",
            run_id,
            "--repo-root",
            str(repo_root),
            "--overwrite",
            "--emit-heatmap-csv",
            "--highlight-ligands",
            "rdk_0002911",
            "--highlight-match",
            "any_contains",
        ],
        check=True,
    )


def _load_report(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _load_heatmap_rows(path: Path):
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _find_highlight(report, target_id, query):
    highlights = report["targets"][target_id]["highlights"]
    return next(h for h in highlights if h.get("query") == query)


def test_pose_invalid_kept_when_filter_false(tmp_path):
    run_id = "POSE_INVALID_KEEP"
    repo_root = _setup_repo(tmp_path, run_id, filter_invalid=False)

    _run_master_export(repo_root, run_id)
    _run_report(repo_root, run_id)

    report_path = repo_root / "data" / run_id / "report.yaml"
    heatmap_path = repo_root / "data" / run_id / "heatmap_input.csv"

    assert report_path.exists()
    report = _load_report(report_path)
    target_id = "1T46|HOLO|pH8_5"
    assert target_id in report["targets"]

    highlight = _find_highlight(report, target_id, "rdk_0002911")
    assert highlight["found"] is True
    assert highlight["pose_valid_any"] is False
    assert highlight["pose_invalid_reason_top"] == "too_far_from_pocket"

    heatmap_rows = _load_heatmap_rows(heatmap_path)
    assert any(
        r.get("target_id") == target_id and r.get("ligand_base") == "rdk_0002911"
        for r in heatmap_rows
    )


def test_pose_invalid_excluded_when_filter_true(tmp_path):
    run_id = "POSE_INVALID_DROP"
    repo_root = _setup_repo(tmp_path, run_id, filter_invalid=True)

    _run_master_export(repo_root, run_id)
    _run_report(repo_root, run_id)

    report_path = repo_root / "data" / run_id / "report.yaml"
    heatmap_path = repo_root / "data" / run_id / "heatmap_input.csv"

    report = _load_report(report_path)
    target_id = "1T46|HOLO|pH8_5"
    highlight = _find_highlight(report, target_id, "rdk_0002911")
    assert highlight["found"] is False

    heatmap_rows = _load_heatmap_rows(heatmap_path)
    assert not any(
        r.get("target_id") == target_id and r.get("ligand_base") == "rdk_0002911"
        for r in heatmap_rows
    )
