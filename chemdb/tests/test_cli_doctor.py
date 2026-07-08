from __future__ import annotations

import io
import importlib
import json
from pathlib import Path

from cli import doctor


def _write_valid_pdbqt(path: Path, *, atom_count: int = 5) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"REMARK fake {path.name}"]
    for idx in range(1, atom_count + 1):
        lines.append(
            f"ATOM  {idx:5d}  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _stub_doctor_baseline(monkeypatch) -> None:
    monkeypatch.setattr(
        importlib,
        "import_module",
        lambda _name: object(),
    )
    monkeypatch.setattr(
        doctor,
        "resolve_tool",
        lambda _cfg, _key, _fallback: {"resolved_path": "/bin/echo", "source": "path"},
    )


def _prep_combo_layout(
    tmp_path: Path,
    *,
    pdb_id: str = "1ABC",
    receptor_atoms: int = 5,
    ligand_count: int = 2,
    include_receptor: bool = True,
    include_ligands: bool = True,
    include_input_pdb: bool = True,
) -> dict[str, str]:
    processed_root = tmp_path / "processed_pdbs"
    prepped_root = tmp_path / "prepped_ligands"
    input_dir = tmp_path / "input_pdbs"
    if include_input_pdb:
        input_dir.mkdir(parents=True, exist_ok=True)
        (input_dir / f"{pdb_id}.pdb").write_text("ATOM      1  CA  ALA A   1\n", encoding="utf-8")
    if include_receptor:
        _write_valid_pdbqt(
            processed_root / pdb_id / "receptor" / f"{pdb_id}.pdbqt",
            atom_count=receptor_atoms,
        )
    if include_ligands:
        ligand_root = prepped_root / "fda_library"
        for idx in range(ligand_count):
            _write_valid_pdbqt(ligand_root / f"lig_{idx:02d}.pdbqt")
    return {
        "INPUT_DIR": str(input_dir),
        "OUTPUT_DIR": str(processed_root),
        "PREPPED_LIGANDS_DIR": str(prepped_root),
    }


def test_run_doctor_optional_import_missing_does_not_fail(monkeypatch) -> None:
    def _fake_import(name: str):
        if name == "openbabel":
            raise ModuleNotFoundError("no module named openbabel")
        return object()

    monkeypatch.setattr(importlib, "import_module", _fake_import)
    monkeypatch.setattr(
        doctor,
        "resolve_tool",
        lambda _cfg, _key, _fallback: {"resolved_path": "/bin/echo", "source": "path"},
    )

    buf = io.StringIO()
    rc = doctor.run_doctor({}, ["atlas", "--doctor"], stream=buf)
    output = buf.getvalue()

    assert rc == 0
    assert "import:openbabel" in output
    assert "WARN" in output
    assert "doctor result: PASS" in output


def test_run_doctor_requires_scorch_when_enabled(monkeypatch) -> None:
    def _fake_resolve(_cfg, key: str, _fallback: str):
        if key == "SCORCH":
            return {"resolved_path": "", "source": "missing"}
        return {"resolved_path": "/bin/echo", "source": "path"}

    monkeypatch.setattr(
        importlib,
        "import_module",
        lambda _name: object(),
    )
    monkeypatch.setattr(doctor, "resolve_tool", _fake_resolve)

    buf = io.StringIO()
    rc = doctor.run_doctor({"USE_SCORCH": True}, ["atlas", "--doctor"], stream=buf)
    output = buf.getvalue()

    assert rc == 2
    assert "tool:SCORCH" in output
    assert "FAIL" in output
    assert "doctor result: FAIL" in output


def test_run_doctor_pdb_and_ligands_happy_path(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    cfg = _prep_combo_layout(tmp_path)

    buf = io.StringIO()
    rc = doctor.run_doctor(
        cfg,
        ["atlas", "doctor", "--pdb", "1abc", "--ligands", "fda"],
        stream=buf,
    )
    output = buf.getvalue()

    assert rc == 0
    assert "target/ligand preflight" in output
    assert "pdb_id: 1ABC" in output
    assert "receptor_pdbqt:" in output
    assert "ligand_root:" in output
    assert "ligands_checked: 2/2" in output
    assert "PASS receptor:" in output
    assert "PASS ligand:" in output
    assert "doctor result: PASS" in output


def test_run_doctor_pdb_shorthand_target_arg(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    cfg = _prep_combo_layout(tmp_path, pdb_id="P04637")

    buf = io.StringIO()
    rc = doctor.run_doctor(
        cfg,
        ["atlas", "doctor", "P04637", "--ligands", "fda"],
        stream=buf,
    )
    output = buf.getvalue()

    assert rc == 0
    assert "pdb_id: P04637" in output
    assert "doctor result: PASS" in output


def test_run_doctor_missing_receptor_pdbqt_reports_fix(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    cfg = _prep_combo_layout(tmp_path, include_receptor=False, include_ligands=False)

    buf = io.StringIO()
    rc = doctor.run_doctor(cfg, ["atlas", "doctor", "--pdb", "1ABC"], stream=buf)
    output = buf.getvalue()

    assert rc == 2
    assert "ERROR receptor:" in output
    assert "receptor PDBQT does not exist" in output
    assert "prepare the target first" in output
    assert "doctor result: FAIL (target/ligand preflight)" in output


def test_run_doctor_missing_ligand_library_dir(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    cfg = _prep_combo_layout(tmp_path, include_ligands=False)

    buf = io.StringIO()
    rc = doctor.run_doctor(cfg, ["atlas", "doctor", "--pdb", "1ABC", "--ligands", "fda"], stream=buf)
    output = buf.getvalue()

    assert rc == 2
    assert "ERROR missing: ligand library directory not found:" in output
    assert "fda_library" in output
    assert "doctor result: FAIL (target/ligand preflight)" in output


def test_run_doctor_empty_ligand_library(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    cfg = _prep_combo_layout(tmp_path, include_ligands=False)
    ligand_root = Path(cfg["PREPPED_LIGANDS_DIR"]) / "fda_library"
    ligand_root.mkdir(parents=True, exist_ok=True)

    buf = io.StringIO()
    rc = doctor.run_doctor(cfg, ["atlas", "doctor", "--pdb", "1ABC", "--ligands", "fda"], stream=buf)
    output = buf.getvalue()

    assert rc == 2
    assert "ERROR missing: no ligand PDBQT files found under:" in output
    assert "fda_library" in output


def test_run_doctor_one_atom_receptor_fails(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    cfg = _prep_combo_layout(tmp_path, receptor_atoms=1, include_ligands=False)

    buf = io.StringIO()
    rc = doctor.run_doctor(cfg, ["atlas", "doctor", "--pdb", "1ABC"], stream=buf)
    output = buf.getvalue()

    assert rc == 2
    assert "ERROR receptor:" in output
    assert "one-atom ligand/receptor is likely a failed conversion" in output


def test_run_doctor_explicit_receptor_and_ligand_file(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    receptor = tmp_path / "custom" / "receptor.pdbqt"
    ligand = tmp_path / "custom" / "ligand.pdbqt"
    _write_valid_pdbqt(receptor)
    _write_valid_pdbqt(ligand)

    buf = io.StringIO()
    rc = doctor.run_doctor(
        {},
        [
            "atlas",
            "doctor",
            "--receptor-pdbqt",
            str(receptor),
            "--ligand-file",
            str(ligand),
        ],
        stream=buf,
    )
    output = buf.getvalue()

    assert rc == 0
    assert "PASS receptor:" in output
    assert "PASS ligand:" in output
    assert "custom/receptor.pdbqt" in output or "receptor.pdbqt" in output


def test_run_doctor_ligand_dir_path(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    ligand_root = tmp_path / "my_ligands"
    _write_valid_pdbqt(ligand_root / "nested" / "one.pdbqt")

    buf = io.StringIO()
    rc = doctor.run_doctor(
        {},
        ["atlas", "doctor", "--ligand-dir", str(ligand_root)],
        stream=buf,
    )
    output = buf.getvalue()

    assert rc == 0
    assert "ligand_root:" in output
    assert "ligands_checked: 1/1" in output
    assert "PASS ligand:" in output


def test_run_doctor_missing_input_pdb(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    cfg = _prep_combo_layout(tmp_path, include_input_pdb=False, include_ligands=False)

    buf = io.StringIO()
    rc = doctor.run_doctor(cfg, ["atlas", "doctor", "--pdb", "1ABC"], stream=buf)
    output = buf.getvalue()

    assert rc == 2
    assert "ERROR missing: input PDB not found:" in output
    assert "1ABC.pdb" in output
    assert "PASS receptor:" in output


def test_run_doctor_json_includes_combo_report(tmp_path: Path, monkeypatch) -> None:
    _stub_doctor_baseline(monkeypatch)
    cfg = _prep_combo_layout(tmp_path)

    buf = io.StringIO()
    rc = doctor.run_doctor(
        cfg,
        ["atlas", "doctor", "--pdb", "1ABC", "--ligands", "fda", "--json"],
        stream=buf,
    )
    payload = json.loads(buf.getvalue())

    assert rc == 0
    assert payload["ok"] is True
    assert payload["combo"]["pdb_id"] == "1ABC"
    assert payload["combo"]["ligand_files_checked"] == 2
    assert payload["combo"]["ligand_files_total"] == 2
    assert payload["combo"]["missing"] == []
    assert any(check["kind"] == "receptor" and check["ok"] for check in payload["combo"]["checks"])
    assert any(check["kind"] == "ligand" and check["ok"] for check in payload["combo"]["checks"])
