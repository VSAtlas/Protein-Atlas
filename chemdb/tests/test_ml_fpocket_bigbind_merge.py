from __future__ import annotations

import pandas as pd
import pytest

from ml.fpocket_bigbind import merge_fpocket_metrics_on_pocket, precompute_fpocket_for_bigbind_df


def _atom_line(
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:<4}{resname:>3} {chain:1}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2}\n"
    )


def test_ml_fpocket_bigbind_merge(tmp_path):
    bigbind_root = tmp_path / "BigBindV1.5"
    target_dir = bigbind_root / "SomeTarget"
    target_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdb = target_dir / "1abc_A_rec.pdb"
    receptor_pdb.write_text(
        "".join(
            [
                _atom_line(1, "N", "ALA", "A", 10, 0.0, 0.0, 0.0, "N"),
                _atom_line(2, "CA", "ALA", "A", 10, 1.0, 0.0, 0.0, "C"),
                _atom_line(3, "C", "ALA", "A", 11, 1.0, 1.0, 0.0, "C"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    pocket_pdb = target_dir / "1abc_A_rec_pocket.pdb"
    pocket_pdb.write_text(
        "".join(
            [
                _atom_line(1, "N", "ALA", "A", 10, 0.0, 0.0, 0.0, "N"),
                _atom_line(2, "CA", "ALA", "A", 10, 2.0, 0.0, 0.0, "C"),
                _atom_line(3, "C", "ALA", "A", 11, 2.0, 2.0, 0.0, "C"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    fpocket_out_root = tmp_path / "fpocket_out"
    info_dir = fpocket_out_root / "1abc_A_rec_out"
    info_dir.mkdir(parents=True, exist_ok=True)
    info_file = info_dir / "1abc_A_rec_info.txt"
    info_file.write_text(
        "\n".join(
            [
                "Pocket 1",
                "Druggability Score : 0.55",
                "Volume : 450.0",
                "Total SASA : 1000.0",
                "Polar SASA : 300.0",
                "Mean alp. sph. solvent access : 0.50",
                "",
            ]
        ),
        encoding="utf-8",
    )

    df = pd.DataFrame(
        [
            {
                "pocket": "SomeTarget/1abc_A_rec_pocket.pdb",
                "lig_smiles": "CCO",
                "active": 1,
                "ex_rec_pdb": "1ABC",
            }
        ]
    )
    run_dir = tmp_path / "run"
    metrics_df = precompute_fpocket_for_bigbind_df(
        df=df,
        bigbind_root=bigbind_root,
        atlas_cfg={
            "FPOCKET_OUTPUT_ROOT": str(fpocket_out_root),
            "FPOCKET_EXE": "/bin/false",
        },
        run_dir=run_dir,
    )

    assert not metrics_df.empty
    row = metrics_df.iloc[0]
    assert float(row["fpocket_druggability"]) == pytest.approx(0.55, abs=1e-6)
    assert float(row["fpocket_volume"]) == pytest.approx(450.0, abs=1e-6)
    assert float(row["fpocket_openness"]) == pytest.approx(0.50, abs=1e-6)
    assert float(row["fpocket_polar_fraction"]) == pytest.approx(0.30, abs=1e-6)
    assert pd.notna(row["pocket_center_x"])
    assert pd.notna(row["pocket_center_y"])
    assert pd.notna(row["pocket_center_z"])

    merged = merge_fpocket_metrics_on_pocket(df, metrics_df)
    for col in (
        "fpocket_druggability",
        "fpocket_volume",
        "fpocket_openness",
        "fpocket_polar_fraction",
        "fpocket_has_metal",
        "fpocket_tier",
        "pocket_center_x",
        "pocket_center_y",
        "pocket_center_z",
    ):
        assert pd.notna(merged.loc[0, col])
