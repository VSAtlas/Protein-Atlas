from pathlib import Path
import subprocess

import protein_prep.receptor_prep as receptor_prep


def test_clean_receptor_pdbqt_removes_user_mod_and_suffixes(tmp_path):
    pdbqt = tmp_path / "receptor.pdbqt"
    pdbqt.write_text(
        "USER  MOD reduce marker\n"
        "ATOM      1  C   LIG A   1      10.0  11.0  12.0  0.00  0.00    C  std\n"
        "HETATM    2  ZN  ZN  A   2      20.0  21.0  22.0  0.00  0.00   ZN new\n"
        "REMARK keep this\n",
        encoding="utf-8",
    )

    receptor_prep._clean_receptor_pdbqt(pdbqt)

    cleaned = pdbqt.read_text(encoding="utf-8")
    assert "USER  MOD" not in cleaned
    assert " std" not in cleaned
    assert " new" not in cleaned
    assert "REMARK keep this" in cleaned


def test_run_prepare_receptor_uses_adt_fallback_after_meeko_failure(tmp_path, monkeypatch):
    input_pdb = tmp_path / "input.pdb"
    output_pdbqt = tmp_path / "out" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    input_pdb.write_text(
        "ATOM      1  N   HIS A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  HIS A   1      11.000  10.000  10.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    fake_python = tmp_path / "pythonsh"
    fake_script = tmp_path / "prepare_receptor4.py"
    fake_python.write_text("#!/bin/sh\n", encoding="utf-8")
    fake_script.write_text("# placeholder\n", encoding="utf-8")
    fake_python.chmod(0o755)

    monkeypatch.setattr(receptor_prep, "MGLTOOLS_PYTHON", str(fake_python))
    monkeypatch.setattr(receptor_prep, "PREPARE_RECEPTOR_SCRIPT", str(fake_script))
    monkeypatch.setattr(
        receptor_prep,
        "_meeko_preflight_or_fail",
        lambda pdb, _work: Path(pdb),
    )

    build_calls = {"modern": 0, "legacy": 0, "adt": 0}

    monkeypatch.setattr(receptor_prep, "build_meeko_base_cmd", lambda: ["MEEKO_BASE"])

    def fake_build_modern(base, pdb_path, out_path):
        build_calls["modern"] += 1
        return ["MEEKO_MODERN"]

    def fake_build_legacy(base, pdb_path, out_path):
        build_calls["legacy"] += 1
        return ["MEEKO_LEGACY"]

    def fake_build_adt(py, script, pdb_path, out_path):
        build_calls["adt"] += 1
        return ["ADT"]

    monkeypatch.setattr(receptor_prep, "build_meeko_modern_cmd", fake_build_modern)
    monkeypatch.setattr(receptor_prep, "build_meeko_legacy_cmd", fake_build_legacy)
    monkeypatch.setattr(receptor_prep, "build_adt_prepare_receptor_cmd", fake_build_adt)

    calls = []
    adt_calls = {"count": 0}

    def fake_run(cmd):
        calls.append(tuple(cmd))
        tag = cmd[0]
        if tag == "ADT":
            adt_calls["count"] += 1
            if adt_calls["count"] == 2:
                output_pdbqt.write_text(
                    "ATOM      1  C   LIG A   1      10.0  11.0  12.0  0.00  0.00    C\n",
                    encoding="utf-8",
                )
                return subprocess.CompletedProcess(cmd, 0, "", "")
            return subprocess.CompletedProcess(cmd, 1, "", "adt fail")
        return subprocess.CompletedProcess(cmd, 1, "", "meeko fail")

    monkeypatch.setattr(receptor_prep, "run_subprocess_capture", fake_run)

    cfg = {
        "HIS_DEFAULT": "HIE",
        "MEEKO_ALLOW_BAD_RES": "true",
        "MGLTOOLS_PYTHON": str(fake_python),
        "PREPARE_RECEPTOR_SCRIPT": str(fake_script),
    }
    ok = receptor_prep.run_prepare_receptor(str(input_pdb), str(output_pdbqt), cfg)

    assert ok is True
    assert build_calls["modern"] == 1
    assert build_calls["legacy"] == 1
    assert build_calls["adt"] == 1
    assert calls.count(("ADT",)) == 2
    assert ("MEEKO_MODERN",) in calls
    assert output_pdbqt.exists()
