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
from post_docking.rescoring import scorch_orchestration
from post_docking.rescoring.scorch_types import StageSpec


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


def test_decoy_prefix_alias_discovers_legacy_dud_stage_dirs(tmp_path: Path) -> None:
    previous = rescoring_scorch.DECOY_PREFIX_VALUE
    try:
        rescoring_scorch._set_decoy_prefix("decoys")
        run_id = "unitrun"
        pdb_id = "1ABC"
        run_combo_root = tmp_path / "docked" / run_id / pdb_id
        _ensure_dir(run_combo_root / "dud_stage3")

        candidates = rescoring_scorch.stage_dir_candidates(
            "vina",
            "dud",
            run_combo_root,
        )
        assert "dud_stage3" in candidates

        mode_dirs = rescoring_scorch._discover_mode_dirs(
            (pdb_id, "", ""),
            tmp_path / "docked" / run_id,
            tmp_path / "post_docked" / run_id,
            [StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")],
            logging.getLogger("test"),
        )
        assert run_combo_root / "dud_stage3" in mode_dirs["dud"]
    finally:
        rescoring_scorch._set_decoy_prefix(previous)


def test_generic_dud_consensus_uses_generic_stage_dirs_as_dud(tmp_path: Path) -> None:
    previous = rescoring_scorch.DECOY_PREFIX_VALUE
    try:
        rescoring_scorch._set_decoy_prefix("decoys")
        run_id = "unitrun"
        pdb_id = "1ABC"
        combo = (pdb_id, "", "")
        run_root = tmp_path / "docked" / run_id
        post_root = tmp_path / "post_docked" / run_id
        combo_root = run_root / pdb_id
        _ensure_dir(combo_root / "stage1")
        _write_csv(
            combo_root / "consensus_docking_scores.csv",
            [
                {
                    "ligand": "actives_final_00001.pdbqt",
                    "run_mode": "dud",
                    "library": "dud",
                    "is_decoy": "0",
                    "consensus_score": "1.0",
                }
            ],
        )
        spec = StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv")

        mode_dirs = rescoring_scorch._discover_mode_dirs(
            combo,
            run_root,
            post_root,
            [spec],
            logging.getLogger("test"),
        )
        score_csv, score_cols, higher_is_better = rescoring_scorch._score_csv_for_spec(
            run_root,
            combo,
            spec,
            "dud",
        )

        assert combo_root / "stage1" in mode_dirs["dud"]
        assert combo_root / "stage1" not in mode_dirs["fda"]
        assert score_csv == combo_root / "consensus_docking_scores.csv"
        assert score_cols == ["consensus_score"]
        assert higher_is_better is True
    finally:
        rescoring_scorch._set_decoy_prefix(previous)


def test_flat_combo_discovery_does_not_treat_stage_dirs_as_ph(tmp_path: Path) -> None:
    run_id = "unitrun"
    pdb_id = "TEST"
    variant = "HOLO"
    variant_root = tmp_path / "docked" / run_id / pdb_id / variant
    _ensure_dir(variant_root / "stage1")
    _ensure_dir(variant_root / "stage2")
    _ensure_dir(variant_root / "stage3")

    combos = rescoring_scorch.discover_combos(
        tmp_path / "docked" / run_id,
        tmp_path / "post_docked" / run_id,
    )

    assert (pdb_id, variant, "") in combos
    assert (pdb_id, variant, "stage1") not in combos
    assert (pdb_id, variant, "stage2") not in combos
    assert (pdb_id, variant, "stage3") not in combos


def test_scoped_combo_discovery_only_walks_requested_pdb(
    tmp_path: Path, monkeypatch
) -> None:
    run_id = "unitrun"
    wanted = tmp_path / "docked" / run_id / "1ABC"
    other = tmp_path / "docked" / run_id / "9ZZZ"
    _ensure_dir(wanted / "stage1")
    _ensure_dir(other / "stage1")
    _ensure_dir(tmp_path / "post_docked" / run_id / "9ZZZ" / "scorch" / "done")
    walked: list[Path] = []
    original_rglob = Path.rglob

    def _spy_rglob(self: Path, pattern: str):
        walked.append(self)
        return original_rglob(self, pattern)

    monkeypatch.setattr(Path, "rglob", _spy_rglob)

    combos = rescoring_scorch.discover_combos(
        tmp_path / "docked" / run_id,
        tmp_path / "post_docked" / run_id,
        pdb_id="1ABC",
    )

    assert combos == {("1ABC", "", "")}
    assert wanted in walked
    assert other not in walked


def test_legacy_flat_combo_discovery_without_variant(tmp_path: Path) -> None:
    run_id = "unitrun"
    pdb_id = "5WIU"
    pdb_root = tmp_path / "docked" / run_id / pdb_id
    _ensure_dir(pdb_root / "stage1")
    _ensure_dir(pdb_root / "stage2")
    _ensure_dir(pdb_root / "stage3")

    combos = rescoring_scorch.discover_combos(
        tmp_path / "docked" / run_id,
        tmp_path / "post_docked" / run_id,
    )

    assert (pdb_id, "", "") in combos
    assert (pdb_id, "stage1", "") not in combos
    assert (pdb_id, "stage2", "") not in combos
    assert (pdb_id, "stage3", "") not in combos


def test_legacy_ph_combo_discovery_without_variant(tmp_path: Path) -> None:
    run_id = "unitrun"
    pdb_id = "5WIU"
    pdb_root = tmp_path / "docked" / run_id / pdb_id / "pH7_0"
    _ensure_dir(pdb_root / "stage1")

    combos = rescoring_scorch.discover_combos(
        tmp_path / "docked" / run_id,
        tmp_path / "post_docked" / run_id,
    )

    assert (pdb_id, "", "pH7_0") in combos
    assert (pdb_id, "pH7_0", "") not in combos


def test_scorch_receptor_path_uses_flat_receptor_without_ph() -> None:
    processed_root = Path("/tmp/processed")

    assert scorch_orchestration.receptor_path_for_combo(
        processed_root,
        ("TEST", "HOLO", ""),
    ) == processed_root / "TEST" / "HOLO" / "receptor" / "TEST.pdbqt"

    assert scorch_orchestration.receptor_path_for_combo(
        processed_root,
        ("TEST", "HOLO", "pH7_0"),
    ) == (
        processed_root
        / "TEST"
        / "HOLO"
        / "receptor"
        / "ph_ensemble"
        / "TEST_pH7_0.pdbqt"
    )


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
