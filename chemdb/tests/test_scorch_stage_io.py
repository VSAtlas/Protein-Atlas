from __future__ import annotations

import logging
from pathlib import Path

from post_docking.rescoring.scorch_stage_io import materialize_inputs


def _write_ligand(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("REMARK test ligand\n", encoding="utf-8")


def test_materialize_inputs_rebuilds_stale_chunk_directory(tmp_path: Path) -> None:
    ph_root = tmp_path / "docked" / "PDB1"
    post_root = tmp_path / "post" / "PDB1"
    old_ligand = ph_root / "stage1" / "old_ligand.pdbqt"
    new_a = ph_root / "stage1" / "new_a.pdbqt"
    new_b = ph_root / "stage1" / "new_b.pdbqt"
    for ligand in (old_ligand, new_a, new_b):
        _write_ligand(ligand)

    first = materialize_inputs(
        ph_root,
        post_root,
        "vina_best",
        [old_ligand],
        overwrite=False,
        logger=logging.getLogger("test.scorch.materialize"),
        run_mode="fda",
        chunk_tag="part000",
    )
    assert first.materialized_count == 1
    assert sorted(path.name for path in first.input_dir.glob("*.pdbqt")) == [
        "old_ligand.pdbqt"
    ]

    second = materialize_inputs(
        ph_root,
        post_root,
        "vina_best",
        [new_a, new_b],
        overwrite=False,
        logger=logging.getLogger("test.scorch.materialize"),
        run_mode="fda",
        chunk_tag="part000",
    )

    assert second.input_dir == first.input_dir
    assert second.materialized_count == 2
    assert sorted(path.name for path in second.input_dir.glob("*.pdbqt")) == [
        "new_a.pdbqt",
        "new_b.pdbqt",
    ]


def test_materialize_inputs_reuses_matching_chunk_directory(tmp_path: Path) -> None:
    ph_root = tmp_path / "docked" / "PDB1"
    post_root = tmp_path / "post" / "PDB1"
    ligand = ph_root / "stage1" / "ligand.pdbqt"
    _write_ligand(ligand)

    first = materialize_inputs(
        ph_root,
        post_root,
        "vina_best",
        [ligand],
        overwrite=False,
        logger=logging.getLogger("test.scorch.materialize"),
        run_mode="fda",
        chunk_tag="part000",
    )
    marker = first.input_dir / "ligand.pdbqt"
    before = marker.stat().st_mtime_ns

    second = materialize_inputs(
        ph_root,
        post_root,
        "vina_best",
        [ligand],
        overwrite=False,
        logger=logging.getLogger("test.scorch.materialize"),
        run_mode="fda",
        chunk_tag="part000",
    )

    assert second.materialized_count == 1
    assert marker.stat().st_mtime_ns == before
