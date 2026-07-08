from __future__ import annotations

import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.reporting import run_report_core as run_report  # noqa: E402


def test_load_scorch_stats_uses_decoy_prefix(tmp_path: Path) -> None:
    decoy_prefix = "fda_dud"
    run_id = "unitrun"
    pdb_id = "1ABC"
    variant = "HOLO"
    ph = "pH7_0"
    combo_dir = tmp_path / "post_docked" / run_id / pdb_id / variant / ph
    combo_dir.mkdir(parents=True, exist_ok=True)

    path = combo_dir / f"{decoy_prefix}_scorch_scores_all.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["scorch_mu_decoy", "scorch_sigma_decoy", "scorch_n_decoys"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "scorch_mu_decoy": "1.2",
                "scorch_sigma_decoy": "0.3",
                "scorch_n_decoys": "10",
            }
        )

    stats = run_report._load_scorch_stats(
        tmp_path, run_id, pdb_id, variant, ph, decoy_prefix=decoy_prefix
    )

    assert stats["mu"] == 1.2
    assert stats["sigma"] == 0.3
    assert stats["n"] == 10
