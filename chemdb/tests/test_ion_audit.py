import logging
from pathlib import Path

from protein_prep import ion_audit


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


def test_audit_ions_and_diff(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setenv("ION_AUDIT", "1")

    pdb_before = tmp_path / "before.pdb"
    pdb_after = tmp_path / "after.pdb"
    pdb_before.write_text(
        "".join(
            [
                _pdb_line("HETATM", 1, "ZN", "ZN", "A", 1, 0.0, 0.0, 0.0, "ZN"),
                _pdb_line("HETATM", 2, "NA", "NA", "A", 2, 1.0, 0.0, 0.0, "NA"),
                _pdb_line("HETATM", 3, "O", "HOH", "A", 3, 2.0, 0.0, 0.0, "O"),
                _pdb_line("HETATM", 4, "C1", "LIG", "A", 4, 3.0, 0.0, 0.0, "C"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )
    pdb_after.write_text(
        "".join(
            [
                _pdb_line("HETATM", 1, "ZN", "ZN", "A", 1, 0.0, 0.0, 0.0, "ZN"),
                _pdb_line("HETATM", 3, "O", "HOH", "A", 3, 2.0, 0.0, 0.0, "O"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    logger = logging.getLogger("test.ion_audit")
    before = ion_audit.audit_ions(pdb_before, "TEST", "input", "HOLO", logger=logger)
    after = ion_audit.audit_ions(
        pdb_after,
        "TEST",
        "final_cleaned",
        "HOLO",
        logger=logger,
    )

    assert before is not None
    assert after is not None
    assert before["counts"]["metals"]["ZN"] == 1
    assert before["counts"]["simple_ions"]["NA"] == 1
    assert before["counts"]["waters"] == 1
    assert before["counts"]["other_het"] == 1

    delta = ion_audit.diff_ions(before, after, logger=logger)
    assert delta["simple_ions"]["NA"] == -1
    assert delta["other_het"] == -1


def test_ion_audit_manager_writes_audit_json(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setenv("ION_AUDIT", "1")

    pdb_file = tmp_path / "input.pdb"
    pdb_file.write_text(
        "".join(
            [
                _pdb_line("HETATM", 1, "ZN", "ZN", "A", 1, 0.0, 0.0, 0.0, "ZN"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    mgr = ion_audit._IonAuditManager("TEST", tmp_path, "HOLO")
    assert mgr.enabled is True

    mgr.probe("input", pdb_file)
    out_json = tmp_path / "audits" / "input_ions.json"
    assert out_json.exists()
