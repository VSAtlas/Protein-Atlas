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


def test_mmgbsa_receptor_sanitize_nterm_and_chainbreak():
    sys.modules.pop("protein_prep_mmgbsa", None)
    protein_prep_mmgbsa = importlib.import_module("protein_prep_mmgbsa")

    lines = [
        _pdb_line("ATOM", 1, "N", "PHE", "A", 17, 0.0, 0.0, 0.0, "N"),
        _pdb_line("ATOM", 2, "H", "PHE", "A", 17, 0.1, 0.0, 0.0, "H"),
        _pdb_line("ATOM", 3, "CA", "PHE", "A", 17, 1.2, 0.0, 0.0, "C"),
        _pdb_line("ATOM", 4, "HA", "PHE", "A", 17, 1.3, 0.0, 0.0, "H"),
        _pdb_line("ATOM", 5, "C", "PHE", "A", 17, 2.0, 0.0, 0.0, "C"),
        _pdb_line("ATOM", 6, "N", "GLY", "A", 18, 2.5, 0.0, 0.0, "N"),
        _pdb_line("ATOM", 7, "C", "PHE", "A", 329, 0.0, 0.0, 0.0, "C"),
        _pdb_line("ATOM", 8, "N", "GLY", "A", 330, 7.0, 0.0, 0.0, "N"),
        "END\n",
    ]

    sanitized, info = protein_prep_mmgbsa._sanitize_receptor_lines(
        kept_lines=lines,
        strip_nterm_h=True,
        insert_ter_on_chainbreak=True,
        chainbreak_cn_max_a=2.2,
    )

    assert lines[2] in sanitized
    assert not any(
        line.startswith("ATOM")
        and line[22:26].strip() == "17"
        and line[12:16].strip() == "H"
        for line in sanitized
    )
    assert any(
        line.startswith("ATOM")
        and line[22:26].strip() == "17"
        and line[12:16].strip() == "HA"
        for line in sanitized
    )

    idx_first_330 = next(
        idx for idx, line in enumerate(sanitized) if line.startswith("ATOM") and line[22:26].strip() == "330"
    )
    idx_last_329 = max(
        idx for idx, line in enumerate(sanitized) if line.startswith("ATOM") and line[22:26].strip() == "329"
    )
    assert any(line.startswith("TER") for line in sanitized[idx_last_329 + 1 : idx_first_330])

    assert info["inserted_TER_count"] == 1
    assert info["nterm_h_stripped_count"] == 1
