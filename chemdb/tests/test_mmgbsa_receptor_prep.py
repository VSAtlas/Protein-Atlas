import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _pdb_line(record, serial, name, resname, chain, resseq, x, y, z, element):
    return (
        f"{record:<6}{serial:>5} {name:<4} {resname:>3} {chain:1}"
        f"{resseq:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 20.00          {element:>2}\n"
    )


def test_mmgbsa_receptor_prep_alias_waters_metals(tmp_path, monkeypatch):
    aliases_path = tmp_path / "aliases.yaml"
    aliases_path.write_text(
        "canonical_waters: [HOH, WAT]\n"
        "canonical_metals: [ZN]\n"
        "canonical_cofactors: []\n",
        encoding="utf-8",
    )

    monkeypatch.setenv("ALIASES_YAML", str(aliases_path))

    sys.modules.pop("activesite", None)
    sys.modules.pop("protein_prep_mmgbsa", None)

    protein_prep_mmgbsa = importlib.import_module("post_docking.mmgbsa.protein_prep_mmgbsa")

    base_dir = tmp_path / "post_docked" / "RUNTEST" / "TESTP" / "HOLO" / "pH7_0" / "stage1"
    base_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = base_dir / "receptor.pdb"

    lines = [
        _pdb_line("ATOM", 1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _pdb_line("ATOM", 2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, "C"),
        _pdb_line("HETATM", 3, "ZN", "ZN", "A", 200, 1.0, 1.0, 1.0, "ZN"),
        _pdb_line("HETATM", 4, "O", "HOH", "A", 201, 1.0, 1.0, 1.0, "O"),
        _pdb_line("HETATM", 5, "O", "HOH", "A", 202, 10.0, 10.0, 10.0, "O"),
        "END\n",
    ]
    pdb_path.write_text("".join(lines), encoding="utf-8")

    result = protein_prep_mmgbsa.prep_mmgbsa_receptor(
        pdb_path=str(pdb_path),
        runid="RUNTEST",
        center=(0.0, 0.0, 0.0),
        radius=5.0,
        force=True,
    )

    output_path = Path(result["output_path"])
    assert output_path.exists()

    out_lines = output_path.read_text(encoding="utf-8").splitlines()

    assert any(
        line.startswith(("ATOM", "HETATM"))
        and line[17:20].strip() == "HOH"
        and line[22:26].strip() == "201"
        for line in out_lines
    )
    assert not any(
        line.startswith(("ATOM", "HETATM"))
        and line[17:20].strip() == "HOH"
        and line[22:26].strip() == "202"
        for line in out_lines
    )
    assert not any(
        line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "ZN" for line in out_lines
    )
