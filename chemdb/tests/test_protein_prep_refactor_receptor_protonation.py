from pathlib import Path
import subprocess

import protein_prep.protonation as protonation
import protein_prep.receptor_prep as receptor_prep


def _atom_line(element: str = "C") -> str:
    return (
        "ATOM      1  CA  GLY A   1      11.000  12.000  13.000"
        f"  1.00 20.00          {element:>2}\n"
    )


def test_assign_protonation_reduce_success(tmp_path, monkeypatch):
    in_pdb = tmp_path / "in.pdb"
    out_pdb = tmp_path / "out.pdb"
    in_pdb.write_text(_atom_line("C"), encoding="utf-8")

    monkeypatch.setattr(protonation, "fix_pdb_elements", lambda *_a, **_k: None)

    def fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(cmd, 0, _atom_line("H"), "")

    monkeypatch.setattr(protonation, "run_subprocess_capture", fake_run)

    result = protonation.assign_protonation_states(
        in_pdb, out_pdb, reduce_exe="/usr/bin/reduce"
    )

    assert result == str(out_pdb)
    assert out_pdb.exists()


def test_assign_protonation_reduce_fail_fallback_openbabel(tmp_path, monkeypatch):
    in_pdb = tmp_path / "in_fail.pdb"
    out_pdb = tmp_path / "out_fail.pdb"
    in_pdb.write_text(_atom_line("C"), encoding="utf-8")

    monkeypatch.setattr(protonation, "fix_pdb_elements", lambda *_a, **_k: None)

    def fail_run(cmd, **kwargs):
        raise RuntimeError("reduce failed")

    fallback_calls = []

    def fake_openbabel(input_pdb, output_pdb, openbabel_path=""):
        fallback_calls.append((str(input_pdb), str(output_pdb), openbabel_path))
        Path(output_pdb).write_text(_atom_line("H"), encoding="utf-8")

    monkeypatch.setattr(protonation, "run_subprocess_capture", fail_run)
    monkeypatch.setattr(protonation, "run_openbabel_add_h", fake_openbabel)

    result = protonation.assign_protonation_states(
        in_pdb, out_pdb, reduce_exe="/usr/bin/reduce"
    )

    assert result == str(out_pdb)
    assert out_pdb.exists()
    assert fallback_calls


def test_assign_protonation_retry_toggle_respects_env(tmp_path, monkeypatch):
    in_pdb = tmp_path / "in_retry.pdb"
    out_pdb = tmp_path / "out_retry.pdb"
    in_pdb.write_text(_atom_line("C"), encoding="utf-8")

    monkeypatch.setenv("REDUCE_RETRY", "0")
    monkeypatch.setenv("FALLBACK_ADDH", "1")
    monkeypatch.setattr(protonation, "fix_pdb_elements", lambda *_a, **_k: None)

    reduce_calls = []

    def fail_run(cmd, **kwargs):
        reduce_calls.append(list(cmd))
        raise RuntimeError("reduce failed")

    fallback_calls = []

    def fake_openbabel(input_pdb, output_pdb, openbabel_path=""):
        fallback_calls.append((str(input_pdb), str(output_pdb), openbabel_path))
        Path(output_pdb).write_text(_atom_line("H"), encoding="utf-8")

    monkeypatch.setattr(protonation, "run_subprocess_capture", fail_run)
    monkeypatch.setattr(protonation, "run_openbabel_add_h", fake_openbabel)

    result = protonation.assign_protonation_states(
        in_pdb, out_pdb, reduce_exe="/usr/bin/reduce"
    )

    assert result == str(out_pdb)
    assert len(reduce_calls) == 1
    assert fallback_calls


def test_run_prepare_receptor_writes_pdbqt_on_success(tmp_path, monkeypatch):
    input_pdb = tmp_path / "receptor.pdb"
    output_pdbqt = tmp_path / "out" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    input_pdb.write_text(_atom_line("C") + "END\n", encoding="utf-8")

    monkeypatch.setattr(receptor_prep, "MGLTOOLS_PYTHON", "")
    monkeypatch.setattr(receptor_prep, "PREPARE_RECEPTOR_SCRIPT", "")
    monkeypatch.setattr(receptor_prep, "_meeko_preflight_or_fail", lambda p, _w: Path(p))
    monkeypatch.setattr(receptor_prep, "build_meeko_base_cmd", lambda: ["MEEKO_BASE"])
    monkeypatch.setattr(
        receptor_prep,
        "build_meeko_modern_cmd",
        lambda base, inp, out: ["MEEKO_MODERN", "--read_pdb", str(inp), "-p", str(out)],
    )
    monkeypatch.setattr(
        receptor_prep,
        "build_meeko_legacy_cmd",
        lambda base, inp, out: ["MEEKO_LEGACY", "-r", str(inp), "-o", str(out)],
    )

    def fake_run(cmd, **kwargs):
        if cmd and cmd[0] == "MEEKO_MODERN":
            output_pdbqt.write_text(_atom_line("C"), encoding="utf-8")
            return subprocess.CompletedProcess(cmd, 0, "", "")
        return subprocess.CompletedProcess(cmd, 1, "", "failed")

    monkeypatch.setattr(receptor_prep, "run_subprocess_capture", fake_run)

    ok = receptor_prep.run_prepare_receptor(str(input_pdb), str(output_pdbqt), {})

    assert ok is True
    assert output_pdbqt.exists()


def test_run_prepare_receptor_passes_default_altloc_when_cfg_set(
    tmp_path, monkeypatch
):
    input_pdb = tmp_path / "receptor_altloc.pdb"
    output_pdbqt = tmp_path / "out2" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    input_pdb.write_text(_atom_line("C") + "END\n", encoding="utf-8")

    monkeypatch.setattr(receptor_prep, "MGLTOOLS_PYTHON", "")
    monkeypatch.setattr(receptor_prep, "PREPARE_RECEPTOR_SCRIPT", "")
    monkeypatch.setattr(receptor_prep, "_meeko_preflight_or_fail", lambda p, _w: Path(p))
    monkeypatch.setattr(receptor_prep, "build_meeko_base_cmd", lambda: ["MEEKO_BASE"])
    monkeypatch.setattr(
        receptor_prep,
        "build_meeko_modern_cmd",
        lambda base, inp, out: ["MEEKO_MODERN", "--read_pdb", str(inp), "-p", str(out)],
    )
    monkeypatch.setattr(
        receptor_prep,
        "build_meeko_legacy_cmd",
        lambda base, inp, out: ["MEEKO_LEGACY", "-r", str(inp), "-o", str(out)],
    )

    commands = []
    call_index = {"n": 0}

    def fake_run(cmd, **kwargs):
        commands.append(list(cmd))
        call_index["n"] += 1
        if call_index["n"] == 1:
            return subprocess.CompletedProcess(cmd, 1, "", "template matching failed")
        output_pdbqt.write_text(_atom_line("C"), encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(receptor_prep, "run_subprocess_capture", fake_run)

    cfg = {"MEEKO_ALLOW_BAD_RES": "true", "MEEKO_DEFAULT_ALTLOC": "A"}
    ok = receptor_prep.run_prepare_receptor(str(input_pdb), str(output_pdbqt), cfg)

    assert ok is True
    assert any(
        "--default_altloc" in cmd and "A" in cmd
        for cmd in commands
    )
