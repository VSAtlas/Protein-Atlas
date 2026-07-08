from __future__ import annotations

import logging
from pathlib import Path

from post_docking.rescoring.scorch_stage_io import materialize_inputs as _materialize_inputs


def test_materialize_inputs_skips_missing_source(tmp_path: Path) -> None:
    ph_root = tmp_path / "docked"
    combo_post_root = tmp_path / "post_docked"
    ph_root.mkdir(parents=True, exist_ok=True)
    combo_post_root.mkdir(parents=True, exist_ok=True)

    existing = ph_root / "decoys_final_00001_stage3.pdbqt"
    existing.write_text("REMARK pose\n", encoding="utf-8")
    missing = ph_root / "decoys_final_00002_stage3.pdbqt"

    result = _materialize_inputs(
        ph_root,
        combo_post_root,
        "vina_best",
        [existing, missing],
        overwrite=True,
        logger=logging.getLogger("test.scorch.materialize"),
        run_mode="fda",
    )

    assert result.input_dir is not None
    assert result.materialized_count == 1
    assert result.failed_count == 1
    assert str(missing) in result.failed_examples
    assert (result.input_dir / existing.name).exists()
    assert not (result.input_dir / missing.name).exists()


def test_materialize_inputs_hardlink_survives_source_unlink(tmp_path: Path) -> None:
    ph_root = tmp_path / "docked"
    combo_post_root = tmp_path / "post_docked"
    ph_root.mkdir(parents=True, exist_ok=True)
    combo_post_root.mkdir(parents=True, exist_ok=True)

    source = ph_root / "decoys_final_00003_stage3.pdbqt"
    source.write_text("ATOM\n", encoding="utf-8")

    result = _materialize_inputs(
        ph_root,
        combo_post_root,
        "vina_best",
        [source],
        overwrite=True,
        logger=logging.getLogger("test.scorch.materialize"),
        run_mode="fda",
    )

    assert result.input_dir is not None
    staged = result.input_dir / source.name
    assert staged.exists()

    source.unlink()
    assert staged.exists()
    assert staged.read_text(encoding="utf-8") == "ATOM\n"


def test_materialize_inputs_uses_chunk_specific_dirs(tmp_path: Path) -> None:
    ph_root = tmp_path / "docked"
    combo_post_root = tmp_path / "post_docked"
    ph_root.mkdir(parents=True, exist_ok=True)
    combo_post_root.mkdir(parents=True, exist_ok=True)

    lig_a = ph_root / "decoys_final_00011_stage3.pdbqt"
    lig_b = ph_root / "decoys_final_00022_stage3.pdbqt"
    lig_a.write_text("POSE_A\n", encoding="utf-8")
    lig_b.write_text("POSE_B\n", encoding="utf-8")

    first = _materialize_inputs(
        ph_root,
        combo_post_root,
        "vina_best",
        [lig_a],
        overwrite=False,
        logger=logging.getLogger("test.scorch.materialize.chunk_a"),
        run_mode="fda",
        chunk_tag="part000",
    )
    second = _materialize_inputs(
        ph_root,
        combo_post_root,
        "vina_best",
        [lig_b],
        overwrite=False,
        logger=logging.getLogger("test.scorch.materialize.chunk_b"),
        run_mode="fda",
        chunk_tag="part001",
    )

    assert first.input_dir is not None
    assert second.input_dir is not None
    assert first.input_dir != second.input_dir
    assert first.materialized_count == 1
    assert second.materialized_count == 1
    assert (first.input_dir / lig_a.name).exists()
    assert not (first.input_dir / lig_b.name).exists()
    assert (second.input_dir / lig_b.name).exists()
    assert not (second.input_dir / lig_a.name).exists()
