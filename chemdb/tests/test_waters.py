from pathlib import Path

from protein_prep import waters


def _pdb_line(
    record: str,
    serial: int,
    atom: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"{record:<6}{serial:>5} {atom:<4} {resname:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{20.00:>6.2f}          {element:>2}\n"
    )


def test_compute_control_centroids_ignores_hydrogens(tmp_path: Path) -> None:
    lig_dir = tmp_path / "ligands"
    lig_dir.mkdir(parents=True, exist_ok=True)
    (lig_dir / "lig1.pdb").write_text(
        "".join(
            [
                _pdb_line("HETATM", 1, "C1", "LIG", "A", 1, 0.0, 0.0, 0.0, "C"),
                _pdb_line("HETATM", 2, "O1", "LIG", "A", 1, 2.0, 0.0, 0.0, "O"),
                _pdb_line("HETATM", 3, "H1", "LIG", "A", 1, 10.0, 0.0, 0.0, "H"),
            ]
        ),
        encoding="utf-8",
    )

    pts = waters.compute_control_centroids(lig_dir)
    assert len(pts) == 1
    x, y, z = pts[0]
    assert round(x, 3) == 1.0
    assert round(y, 3) == 0.0
    assert round(z, 3) == 0.0


def test_filter_waters_near_points(tmp_path: Path) -> None:
    src = tmp_path / "src.pdb"
    dst = tmp_path / "dst.pdb"

    src.write_text(
        "".join(
            [
                _pdb_line("ATOM", 1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
                _pdb_line("HETATM", 2, "O", "HOH", "A", 10, 1.0, 0.0, 0.0, "O"),
                _pdb_line("HETATM", 3, "O", "HOH", "A", 11, 10.0, 0.0, 0.0, "O"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    kept = waters.filter_waters_near_points(src, dst, [(0.0, 0.0, 0.0)], 2.0)
    out = dst.read_text(encoding="utf-8")

    assert kept == 1
    assert " A  10" in out
    assert " A  11" not in out
