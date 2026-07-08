from pathlib import Path

from protein_prep.altloc_filter import filter_altlocs


def _atom_line(
    serial: int,
    atom_name: str,
    altloc: str,
    *,
    resname: str = "GLY",
    chain: str = "A",
    resseq: int = 1,
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0,
    occ: float = 1.0,
    element: str = "C",
) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:<4}{altloc:1}{resname:>3} {chain:1}"
        f"{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{20.00:6.2f}"
        f"          {element:>2}\n"
    )


def _atom_lines(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln.startswith(("ATOM  ", "HETATM"))]


def _pick_atom(lines: list[str], atom_name: str, chain: str, resseq: str) -> str:
    matches = [
        ln
        for ln in lines
        if ln[12:16].strip() == atom_name
        and ln[21] == chain
        and ln[22:26].strip() == resseq
    ]
    assert len(matches) == 1
    return matches[0]


def test_filter_altlocs_prefers_blank_then_a(tmp_path: Path):
    in_pdb = tmp_path / "input_altloc.pdb"
    out_pdb = tmp_path / "output_altloc.pdb"

    in_pdb.write_text(
        "".join(
            [
                "HEADER    ALTLOC TEST\n",
                _atom_line(1, "CA", "B", occ=0.40, x=1.0),
                _atom_line(2, "CA", "A", occ=0.50, x=2.0),
                _atom_line(3, "CA", " ", occ=0.60, x=3.0),
                _atom_line(4, "CB", "B", occ=0.40, x=4.0),
                _atom_line(5, "CB", "A", occ=0.50, x=5.0),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    filter_altlocs(in_pdb, out_pdb)

    assert out_pdb.exists()
    output_text = out_pdb.read_text(encoding="utf-8")
    atoms = _atom_lines(output_text)

    ca = _pick_atom(atoms, atom_name="CA", chain="A", resseq="1")
    cb = _pick_atom(atoms, atom_name="CB", chain="A", resseq="1")
    assert ca[16] == " "
    assert cb[16] == "A"
    assert " CA BGLY A   1" not in output_text
    assert " CB BGLY A   1" not in output_text


def test_filter_altlocs_no_alternates_keeps_content(tmp_path: Path):
    in_pdb = tmp_path / "input_no_altloc.pdb"
    out_pdb = tmp_path / "output_no_altloc.pdb"

    input_text = "".join(
        [
            _atom_line(1, "CA", " ", x=1.0, element="C"),
            _atom_line(2, "N", " ", x=2.0, element="N"),
        ]
    )
    in_pdb.write_text(input_text, encoding="utf-8")

    filter_altlocs(in_pdb, out_pdb)

    assert out_pdb.exists()
    assert out_pdb.read_text(encoding="utf-8") == input_text
