from __future__ import annotations

import json
from pathlib import Path

import pytest

from cli import distributed_chunk_planner as chunk_planner
from cli import planner_combo_plan
from cli import planner_validation
from docking.completion_markers import completion_marker_path


def _write_pdbqt(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Keep file comfortably above the planner size threshold.
    path.write_text("MODEL\n" + ("ATOM\n" * 64) + "ENDMDL\n", encoding="utf-8")


def _base_cfg(prepped_root: Path, *, ph_ligand_mode: str = "off") -> dict:
    return {
        "PREPPED_LIGANDS_DIR": str(prepped_root),
        "LIBRARY_SUBDIR_DEFAULT": "fda_library",
        "TEST_LIBRARY_MAP": {"BNJS": "testlib"},
        "PH_LIGAND_MODE": ph_ligand_mode,
    }


class _FakeDistContext:
    task_id = 0
    task_count = 10

    def expected_task_ids(self) -> list[int]:
        return list(range(int(self.task_count)))


class _ThreeTaskDistContext:
    task_id = 0
    task_count = 3

    def expected_task_ids(self) -> list[int]:
        return list(range(int(self.task_count)))


def test_distributed_combo_work_items_string_false_is_non_ph_mode(
    tmp_path: Path,
) -> None:
    cfg = {"PH_ENSEMBLE": "false", "INPUT_DIR": str(tmp_path)}
    items = chunk_planner._build_combo_work_items(
        cfg, ["QRST.pdb"], variant_token="HOLO"
    )
    assert items == [("QRST.pdb", None)]


def test_resolve_combo_ligands_prefers_files_and_clamps_target(tmp_path: Path) -> None:
    prepped = tmp_path / "prepped_ligands"
    root = prepped / "testlib"
    for idx in range(1, 7):
        _write_pdbqt(root / f"decoys_final_{idx:05d}.pdbqt")
    _write_pdbqt(root / "microstates" / "decoys_final_00001__ms_deadbeef.pdbqt")
    _write_pdbqt(root / "7.0" / "decoysph7_0_final_00099.pdbqt")

    # Mismatched manifest: includes stale + filtered entries and omits valid files.
    (root / "_manifest.json").write_text(
        json.dumps(
            {
                "entries": {
                    "stale_only": "decoys_final_99999.pdbqt",
                    "microstate_only": "microstates/decoys_final_00001__ms_deadbeef.pdbqt",
                    "ph_only": "7.0/decoysph7_0_final_00099.pdbqt",
                }
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    cfg = _base_cfg(prepped)
    library, ligands = chunk_planner._resolve_combo_ligand_bases(
        cfg,
        pdb_id="BNJS",
        run_tokens=["dud"],
        target_ligands=96,
        sample_seed=1337,
    )
    assert library == "testlib"
    assert len(ligands) == 6
    assert all("__ms_" not in x for x in ligands)
    assert all("decoysph" not in x for x in ligands)


def test_resolve_combo_ligands_applies_deterministic_cap(tmp_path: Path) -> None:
    prepped = tmp_path / "prepped_ligands"
    root = prepped / "testlib"
    for idx in range(1, 13):
        _write_pdbqt(root / f"decoys_final_{idx:05d}.pdbqt")

    cfg = _base_cfg(prepped)
    _, first = chunk_planner._resolve_combo_ligand_bases(
        cfg,
        pdb_id="BNJS",
        run_tokens=["dud"],
        target_ligands=4,
        sample_seed=2026,
    )
    _, second = chunk_planner._resolve_combo_ligand_bases(
        cfg,
        pdb_id="BNJS",
        run_tokens=["dud"],
        target_ligands=4,
        sample_seed=2026,
    )
    assert len(first) == 4
    assert first == second


def test_resolve_combo_ligands_keeps_microstates_when_ph_mode_on(tmp_path: Path) -> None:
    prepped = tmp_path / "prepped_ligands"
    root = prepped / "testlib"
    _write_pdbqt(root / "decoys_final_00001.pdbqt")
    _write_pdbqt(root / "microstates" / "decoys_final_00002__ms_abc123.pdbqt")
    _write_pdbqt(root / "7.0" / "decoysph7_0_final_00003.pdbqt")

    cfg = _base_cfg(prepped, ph_ligand_mode="on")
    _, ligands = chunk_planner._resolve_combo_ligand_bases(
        cfg,
        pdb_id="BNJS",
        run_tokens=["dud"],
        target_ligands=0,
        sample_seed=7,
    )
    assert len(ligands) == 3
    assert any("__ms_" in x for x in ligands)
    assert any("decoysph7_0" in x for x in ligands)


def test_resolve_combo_ligand_streams_separates_dud_and_fda_roots(
    tmp_path: Path,
) -> None:
    prepped = tmp_path / "prepped_ligands"
    dud_root = prepped / "testlib"
    fda_root = prepped / "fda_library"
    for idx in range(1, 5):
        _write_pdbqt(dud_root / f"decoys_final_{idx:05d}.pdbqt")
        _write_pdbqt(fda_root / f"rdk_{idx:07d}.pdbqt")

    cfg = _base_cfg(prepped)
    streams = chunk_planner._resolve_combo_ligand_streams(
        cfg,
        pdb_id="BNJS",
        run_tokens=["dud", "fda"],
        target_ligands=0,
        sample_seed=1337,
    )

    by_mode = {stream["run_mode"]: stream for stream in streams}
    assert set(by_mode) == {"dud", "fda"}
    assert by_mode["dud"]["library_name"] == "testlib"
    assert by_mode["fda"]["library_name"] == "fda_library"
    assert all(str(base).startswith("decoys") for base in by_mode["dud"]["ligand_bases"])
    assert all(str(base).startswith("rdk") for base in by_mode["fda"]["ligand_bases"])


def test_build_lpt_combo_chunks_uses_stream_identity_and_production_chunk_floor(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("ATLAS_CHUNK_MIN_SIZE", raising=False)
    monkeypatch.setenv("ATLAS_CHUNK_PREFLIGHT_VALIDATE", "1")
    prepped = tmp_path / "prepped_ligands"
    dud_root = prepped / "testlib"
    fda_root = prepped / "fda_library"
    for idx in range(1, 129):
        _write_pdbqt(dud_root / f"decoys_final_{idx:05d}.pdbqt")
        _write_pdbqt(fda_root / f"rdk_{idx:07d}.pdbqt")

    chunks = planner_combo_plan._build_lpt_combo_chunks(
        _base_cfg(prepped) | {"CPU": 112},
        run_id="rid",
        dist_ctx=_FakeDistContext(),
        variant_label="LEGACY",
        combo_items=[("BNJS.pdb", None)],
        run_tokens=["dud", "fda"],
        target_ligands=0,
        sample_seed=1337,
    )

    assert {chunk["run_mode"] for chunk in chunks} == {"dud", "fda"}
    assert {chunk["library_name"] for chunk in chunks} == {"testlib", "fda_library"}
    assert all(chunk["library_root"] for chunk in chunks)
    assert min(int(chunk["ligand_count"]) for chunk in chunks) > 4
    for chunk in chunks:
        bases = [str(base) for base in chunk["ligand_bases"]]
        if chunk["run_mode"] == "dud":
            assert all(base.startswith("decoys") for base in bases)
        if chunk["run_mode"] == "fda":
            assert all(base.startswith("rdk") for base in bases)


def test_build_lpt_combo_chunks_canary_worker_scaling_splits_tiny_streams(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("ATLAS_CHUNK_MIN_SIZE", "8")
    monkeypatch.setenv("ATLAS_CHUNK_SCALE_WITH_WORKERS", "1")
    prepped = tmp_path / "prepped_ligands"
    dud_root = prepped / "testlib"
    fda_root = prepped / "fda_library"
    for idx in range(1, 56):
        _write_pdbqt(dud_root / f"decoys_final_{idx:05d}.pdbqt")
    for idx in range(1, 11):
        _write_pdbqt(fda_root / f"rdk_{idx:07d}.pdbqt")

    chunks = planner_combo_plan._build_lpt_combo_chunks(
        _base_cfg(prepped) | {"CPU": 8},
        run_id="rid",
        dist_ctx=_FakeDistContext(),
        variant_label="LEGACY",
        combo_items=[("BNJS.pdb", None)],
        run_tokens=["dud", "fda"],
        target_ligands=0,
        sample_seed=1337,
    )

    dud_chunks = [chunk for chunk in chunks if chunk["run_mode"] == "dud"]
    fda_chunks = [chunk for chunk in chunks if chunk["run_mode"] == "fda"]
    assert len(dud_chunks) > 1
    assert sum(int(chunk["ligand_count"]) for chunk in dud_chunks) == 55
    assert sum(int(chunk["ligand_count"]) for chunk in fda_chunks) == 10
    assert {chunk["library_name"] for chunk in chunks} == {"testlib", "fda_library"}


def test_build_lpt_combo_chunks_scales_production_by_workers_by_default(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("ATLAS_CHUNK_MIN_SIZE", raising=False)
    monkeypatch.delenv("ATLAS_CHUNK_SCALE_WITH_WORKERS", raising=False)
    prepped = tmp_path / "prepped_ligands"
    root = prepped / "testlib"
    for idx in range(1, 385):
        _write_pdbqt(root / f"decoys_final_{idx:05d}.pdbqt")

    chunks = planner_combo_plan._build_lpt_combo_chunks(
        _base_cfg(prepped) | {"CPU": 48},
        run_id="rid",
        dist_ctx=_ThreeTaskDistContext(),
        variant_label="HOLO",
        combo_items=[("BNJS.pdb", None)],
        run_tokens=["dud"],
        target_ligands=0,
        sample_seed=1337,
    )

    assert len(chunks) >= 10
    assert max(int(chunk["ligand_count"]) for chunk in chunks) <= 40
    assert sum(int(chunk["ligand_count"]) for chunk in chunks) == 384


def test_chunk_stream_resolution_preflight_rejects_mixed_root(
    tmp_path: Path,
) -> None:
    prepped = tmp_path / "prepped_ligands"
    dud_root = prepped / "testlib"
    fda_root = prepped / "fda_library"
    _write_pdbqt(dud_root / "decoys_final_00001.pdbqt")
    _write_pdbqt(fda_root / "rdk_0000001.pdbqt")

    cfg = _base_cfg(prepped)
    with pytest.raises(RuntimeError, match="chunk_resolution_preflight_failed"):
        chunk_planner._validate_chunk_stream_resolution(
            cfg,
            chunks=[
                {
                    "chunk_id": "bad",
                    "run_mode": "dud",
                    "library_name": "testlib",
                    "library_root": str(dud_root.resolve()),
                    "ligand_bases": ["rdk_0000001"],
                }
            ],
            sample_chunks_per_stream=1,
            max_loss_fraction=0.05,
        )


def test_resolve_combo_post_csv_prefers_canonical_over_prefixed(tmp_path: Path) -> None:
    combo_post_dir = tmp_path / "post_docked" / "runx" / "BNJS" / "HOLO" / "pH7_0"
    combo_post_dir.mkdir(parents=True, exist_ok=True)
    canonical = combo_post_dir / "consensus_reranked_scorch.csv"
    prefixed = combo_post_dir / "test_library_20_consensus_reranked_scorch.csv"
    canonical.write_text("ligand,final_score\nactives_final_00001.pdbqt,1.0\n", encoding="utf-8")
    prefixed.write_text("ligand,final_score\ndecoys_final_00001.pdbqt,1.0\n", encoding="utf-8")

    resolved = planner_validation._resolve_combo_post_consensus_csv(
        combo_post_dir,
        "test_library_20",
    )
    assert resolved == canonical


def test_verify_chunk_combo_outputs_uses_run_mode_prefixed_long_csv(
    tmp_path: Path,
) -> None:
    combo_dir = tmp_path / "docked" / "rid" / "BNJS"
    combo_dir.mkdir(parents=True, exist_ok=True)
    (combo_dir / "docking_score_summary.csv").write_text(
        "run_id,Ligand,stage1\nrid,rdk_0000001.pdbqt,-1.0\n",
        encoding="utf-8",
    )
    (combo_dir / "dud_docking_score_long.csv").write_text(
        "run_id,stage,ligand,score,valid\n"
        "rid,dud_stage1,decoys_test_library_10_2.pdbqt,-3.2,1\n",
        encoding="utf-8",
    )
    marker = completion_marker_path(combo_dir / "dud_stage1", engine="vina", chunk_id="c1")
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text("{}", encoding="utf-8")

    ok, reason = planner_validation._verify_chunk_combo_outputs(
        {"DOCKED_DIR": str(tmp_path / "docked")},
        run_id="rid",
        chunk_id="c1",
        pdb_id="BNJS",
        variant_label="LEGACY",
        ph_tag=None,
        library_name="test_library_10",
        chunk_ligand_bases=["decoys_test_library_10_2"],
        run_mode="dud",
    )

    assert ok
    assert reason == "ok"
