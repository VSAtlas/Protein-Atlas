from pathlib import Path
from types import SimpleNamespace

from protein_prep import phenix_tools


def test_run_phenix_pdbtools_builds_expected_command(monkeypatch, tmp_path: Path) -> None:
    inp = tmp_path / "in.pdb"
    out = tmp_path / "out.pdb"
    inp.write_text("END\n", encoding="utf-8")

    calls = []

    def _fake_run(cmd, env=None, check=False):
        calls.append((cmd, env, check))
        for token in cmd:
            if token.startswith("output.file_name="):
                out_path = Path(token.split("=", 1)[1])
                out_path.write_text("ATOM\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(
        phenix_tools, "_phenix_detect", lambda: ("linux-path", ["phenix.pdbtools"], {})
    )
    monkeypatch.setattr(phenix_tools, "_run_and_log", _fake_run)

    ok = phenix_tools.run_phenix_pdbtools(inp, out, remove_waters=True)
    assert ok is True
    assert out.exists()
    cmd = calls[0][0]
    assert cmd[0] == "phenix.pdbtools"
    assert str(inp.resolve()).replace("\\", "/") in cmd
    assert f"output.file_name={str(out.resolve()).replace(chr(92), '/')}" in cmd
    assert 'remove="resname HOH"' in cmd


def test_run_phenix_pdbtools_failure_returns_false(monkeypatch, tmp_path: Path) -> None:
    inp = tmp_path / "in.pdb"
    out = tmp_path / "out.pdb"
    inp.write_text("END\n", encoding="utf-8")

    monkeypatch.setattr(
        phenix_tools, "_phenix_detect", lambda: ("linux-path", ["phenix.pdbtools"], {})
    )
    monkeypatch.setattr(
        phenix_tools,
        "_run_and_log",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
    )
    monkeypatch.setattr(phenix_tools, "PHENIX_DIR", "", raising=False)
    monkeypatch.setattr(phenix_tools, "PHENIX_PYTHON_BAT", "", raising=False)
    monkeypatch.setattr(phenix_tools.shutil, "which", lambda *_a, **_k: None)

    assert phenix_tools.run_phenix_pdbtools(inp, out) is False


def test_run_windows_phenix_clean_script_invokes_powershell(
    monkeypatch,
    tmp_path: Path,
) -> None:
    in_pdb = tmp_path / "loop_fixed.pdb"
    out_dir = tmp_path / "nolig"
    py_bat = tmp_path / "phenix.python.bat"
    clean_script = tmp_path / "clean.py"
    in_pdb.write_text("END\n", encoding="utf-8")
    py_bat.write_text("@echo off\n", encoding="utf-8")
    clean_script.write_text("print('ok')\n", encoding="utf-8")

    payload = {}

    monkeypatch.setattr(phenix_tools, "PHENIX_PYTHON_BAT", str(py_bat), raising=False)
    monkeypatch.setattr(
        phenix_tools,
        "PHENIX_CLEAN_SCRIPT",
        str(clean_script),
        raising=False,
    )
    monkeypatch.setattr(phenix_tools, "_win_path", lambda p: str(p))

    def _fake_powershell(script: str):
        payload["script"] = script
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(phenix_tools, "_powershell", _fake_powershell)

    rc = phenix_tools.run_windows_phenix_clean_script(in_pdb, out_dir)
    assert rc == 0
    script = payload["script"]
    assert str(py_bat) in script
    assert str(clean_script) in script
    assert str(out_dir) in script


def test_run_molprobity_validate_returns_subprocess_code(
    monkeypatch,
    tmp_path: Path,
) -> None:
    pdb = tmp_path / "receptor.pdb"
    pdb.write_text("END\n", encoding="utf-8")

    seen = {}

    def _fake_run(cmd, cwd=None, stdout=None, stderr=None, text=None):
        seen["cmd"] = cmd
        seen["cwd"] = cwd
        return SimpleNamespace(returncode=2, stdout="", stderr="failed")

    monkeypatch.setattr(
        phenix_tools.shutil,
        "which",
        lambda name: "/usr/bin/phenix.molprobity" if name == "phenix.molprobity" else None,
    )
    monkeypatch.setattr(phenix_tools.subprocess, "run", _fake_run)

    rc = phenix_tools.run_molprobity_validate(pdb)
    assert rc == 2
    assert seen["cmd"][0] == "/usr/bin/phenix.molprobity"
    assert seen["cmd"][1] == str(pdb.resolve())
