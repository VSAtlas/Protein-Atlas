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


def test_run_prepare_receptor_fails_cleanly_after_meeko_failure(tmp_path, monkeypatch):
    input_pdb = tmp_path / "input.pdb"
    output_pdbqt = tmp_path / "out" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    input_pdb.write_text(
        "ATOM      1  N   HIS A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  HIS A   1      11.000  10.000  10.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(
        receptor_prep,
        "_meeko_preflight_or_fail",
        lambda pdb, _work: Path(pdb),
    )

    build_calls = {"modern": 0, "legacy": 0}

    monkeypatch.setattr(receptor_prep, "build_meeko_base_cmd", lambda: ["MEEKO_BASE"])

    def fake_build_modern(base, pdb_path, out_path):
        build_calls["modern"] += 1
        return ["MEEKO_MODERN"]

    def fake_build_legacy(base, pdb_path, out_path):
        build_calls["legacy"] += 1
        return ["MEEKO_LEGACY"]

    monkeypatch.setattr(receptor_prep, "build_meeko_modern_cmd", fake_build_modern)
    monkeypatch.setattr(receptor_prep, "build_meeko_legacy_cmd", fake_build_legacy)

    calls = []

    def fake_run(cmd):
        calls.append(tuple(cmd))
        return subprocess.CompletedProcess(cmd, 1, "", "meeko fail")

    monkeypatch.setattr(receptor_prep, "run_subprocess_capture", fake_run)

    cfg = {
        "HIS_DEFAULT": "HIE",
        "MEEKO_ALLOW_BAD_RES": "true",
    }
    ok = receptor_prep.run_prepare_receptor(str(input_pdb), str(output_pdbqt), cfg)

    assert ok is False
    assert build_calls["modern"] == 1
    assert build_calls["legacy"] == 1
    assert ("MEEKO_MODERN",) in calls
    assert ("MEEKO_LEGACY",) in calls
