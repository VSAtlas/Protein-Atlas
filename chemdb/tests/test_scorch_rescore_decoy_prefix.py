from __future__ import annotations

import sys
from pathlib import Path
import csv
import logging

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import post_docking.rescoring.rescoring_scorch as rescoring_scorch
import post_docking.rescoring.rescore_reranker as rescore_reranker


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def test_decoy_prefix_discovery(tmp_path: Path) -> None:
    previous = rescoring_scorch.DECOY_PREFIX_VALUE
    try:
        rescoring_scorch._set_decoy_prefix("fda_dud")
        run_id = "unitrun"
        pdb_id = "1ABC"
        variant = "HOLO"
        ph = "pH7_0"
        run_combo_root = tmp_path / "docked" / run_id / pdb_id / variant / ph
        _ensure_dir(run_combo_root / "fda_dud_stage3")
        _ensure_dir(run_combo_root / "gnina_fda_dud_stage3")

        post_combo_root = tmp_path / "post_docked" / run_id / pdb_id / variant / ph
        _ensure_dir(post_combo_root / "fda_dud_ledock_pdbqt")

        vina_candidates = rescoring_scorch.stage_dir_candidates("vina", "dud")
        assert vina_candidates[0] == "fda_dud_stage3"

        gnina_candidates = rescoring_scorch.stage_dir_candidates(
            "gnina", "dud", run_combo_root
        )
        assert "gnina_fda_dud_stage3" in gnina_candidates

        combos = rescoring_scorch.discover_combos(
            tmp_path / "docked" / run_id, tmp_path / "post_docked" / run_id
        )
        assert (pdb_id, variant, ph) in combos
    finally:
        rescoring_scorch._set_decoy_prefix(previous)


def test_multi_prefix_from_test_mode_and_reranker(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setenv("TEST_MODE_ENABLE", "dud_fda+fda+inhouse")
    cfg = {"TEST_MODE_ENABLE": "off", "DUD_PREFIX": "dud"}
    prefixes = rescoring_scorch._decoy_prefixes_from_test_mode(cfg, override=None)
    assert prefixes == ["dud_fda", "inhouse"]

    consensus_rows = [
        {
            "run_id": "r1",
            "pdb_id": "P1",
            "variant": "APO",
            "ph_label": "",
            "ligand": "fda_a_stage1.pdbqt",
            "consensus_score": "1.0",
            "library": "FDA",
        },
        {
            "run_id": "r1",
            "pdb_id": "P1",
            "variant": "APO",
            "ph_label": "",
            "ligand": "dud_fda_a_stage1.pdbqt",
            "consensus_score": "0.6",
            "library": "DECOY",
        },
        {
            "run_id": "r1",
            "pdb_id": "P1",
            "variant": "APO",
            "ph_label": "",
            "ligand": "inhouse_a_stage1.pdbqt",
            "consensus_score": "0.5",
            "library": "DECOY",
        },
    ]
    consensus_path = tmp_path / "consensus_docking_scores.csv"
    _write_csv(consensus_path, consensus_rows)

    fda_rows = [
        {
            "source": "vina",
            "Ligand_ID": "fda_a_stage1",
            "SCORCH_score": "1.0",
            "SCORCH_certainty": "1.0",
            "run_mode": "fda",
        }
    ]
    scorch_path = tmp_path / "scorch_scores_all.csv"
    _write_csv(scorch_path, fda_rows)

    previous = rescore_reranker.DECOY_PREFIX_VALUE
    try:
        for prefix in prefixes:
            decoy_rows = [
                {
                    "source": "vina",
                    "Ligand_ID": f"{prefix}_a_stage1",
                    "SCORCH_score": "0.5",
                    "SCORCH_certainty": "1.0",
                    "run_mode": "dud",
                }
            ]
            decoy_path = tmp_path / f"{prefix}_scorch_scores_all.csv"
            _write_csv(decoy_path, decoy_rows)
            out_csv = tmp_path / "consensus_reranked_scorch.csv"
            ok = rescore_reranker.rerank_consensus_with_scorch(
                consensus_path,
                [scorch_path, decoy_path],
                out_csv,
                logging.getLogger("test_multi_prefix"),
                overwrite=True,
                decoy_prefix=prefix,
            )
            assert ok
            assert (tmp_path / f"{prefix}_consensus_reranked_scorch.csv").exists()
    finally:
        rescore_reranker._set_decoy_prefix(previous)
