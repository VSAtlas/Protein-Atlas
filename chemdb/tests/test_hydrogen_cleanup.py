from pathlib import Path

from protein_prep import hydrogen_cleanup


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


def test_remove_unbonded_atoms_removes_only_unbonded_hydrogen(tmp_path: Path) -> None:
    pdb = tmp_path / "input.pdb"
    pdb.write_text(
        "".join(
            [
                _pdb_line("ATOM", 1, "C1", "LIG", "A", 1, 0.0, 0.0, 0.0, "C"),
                _pdb_line("ATOM", 2, "H1", "LIG", "A", 1, 0.9, 0.0, 0.0, "H"),
                _pdb_line("ATOM", 3, "H2", "LIG", "A", 1, 8.0, 0.0, 0.0, "H"),
                "CONECT    1    2\n",
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    hydrogen_cleanup.remove_unbonded_atoms(pdb)

    text = pdb.read_text(encoding="utf-8")
    assert " H1 " in text
    assert " H2 " not in text
    assert " C1 " in text


def test_clean_hydrogens_removes_implausible_hydrogen(monkeypatch, tmp_path: Path) -> None:
    pdb = tmp_path / "geometry.pdb"
    pdb.write_text(
        "".join(
            [
                _pdb_line("ATOM", 1, "C1", "LIG", "A", 1, 0.0, 0.0, 0.0, "C"),
                _pdb_line("ATOM", 2, "H1", "LIG", "A", 1, 0.9, 0.0, 0.0, "H"),
                _pdb_line("ATOM", 3, "H2", "LIG", "A", 1, 10.0, 0.0, 0.0, "H"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(hydrogen_cleanup, "conect_coverage", lambda _path: 0.0)
    monkeypatch.setattr(
        hydrogen_cleanup,
        "_post_write_element_guard",
        lambda *_args, **_kwargs: None,
    )

    hydrogen_cleanup.clean_hydrogens(
        pdb,
        use_conect_if_reliable=True,
        conect_min_cov=0.6,
    )

    text = pdb.read_text(encoding="utf-8")
    assert " H1 " in text
    assert " H2 " not in text
