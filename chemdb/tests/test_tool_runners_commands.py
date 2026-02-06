import subprocess

import pytest

import protein_prep.tool_runners as tool_runners


def test_build_openbabel_command_matches_existing_shape():
    cmd = tool_runners.build_openbabel_add_h_cmd(
        "/usr/bin/obabel", "/tmp/in.pdb", "/tmp/out.pdb"
    )
    assert cmd == ["/usr/bin/obabel", "/tmp/in.pdb", "-O", "/tmp/out.pdb", "-h"]


def test_build_meeko_and_adt_commands_match_existing_shape():
    base = ["/usr/bin/python", "-m", "meeko.cli.mk_prepare_receptor"]
    modern = tool_runners.build_meeko_modern_cmd(base, "/tmp/receptor.pdb", "/tmp/r.pdbqt")
    legacy = tool_runners.build_meeko_legacy_cmd(base, "/tmp/receptor.pdb", "/tmp/r.pdbqt")
    adt = tool_runners.build_adt_prepare_receptor_cmd(
        "/opt/mgltools/bin/pythonsh",
        "/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py",
        "/tmp/receptor.pdb",
        "/tmp/r.pdbqt",
    )

    assert modern == base + ["--read_pdb", "/tmp/receptor.pdb", "-p", "/tmp/r.pdbqt"]
    assert legacy == base + ["-r", "/tmp/receptor.pdb", "-o", "/tmp/r.pdbqt"]
    assert adt == [
        "/opt/mgltools/bin/pythonsh",
        "/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py",
        "-r",
        "/tmp/receptor.pdb",
        "-o",
        "/tmp/r.pdbqt",
        "-A",
        "none",
        "-U",
        "nphs_lps",
    ]


def test_run_openbabel_add_h_invokes_subprocess_with_expected_args(monkeypatch):
    captured = {}

    monkeypatch.setattr(tool_runners.shutil, "which", lambda name: "/usr/bin/obabel")

    def fake_run(cmd):
        captured["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(tool_runners, "run_subprocess_capture", fake_run)

    tool_runners.run_openbabel_add_h(
        "/tmp/in.pdb", "/tmp/out.pdb", openbabel_path="/usr/bin/obabel"
    )

    assert captured["cmd"] == [
        "/usr/bin/obabel",
        "/tmp/in.pdb",
        "-O",
        "/tmp/out.pdb",
        "-h",
    ]


def test_run_with_fallback_invokes_fallback_on_primary_failure():
    calls = []

    def fake_run(cmd):
        calls.append(list(cmd))
        rc = 1 if len(calls) == 1 else 0
        return subprocess.CompletedProcess(cmd, rc, "", "err" if rc else "")

    result, used_fallback = tool_runners.run_with_fallback(
        ["primary", "--flag"], ["fallback", "--flag"], run_cmd=fake_run
    )

    assert used_fallback is True
    assert result.returncode == 0
    assert calls == [["primary", "--flag"], ["fallback", "--flag"]]


def test_run_openbabel_add_h_raises_when_tool_fails(monkeypatch):
    monkeypatch.setattr(tool_runners.shutil, "which", lambda name: "/usr/bin/obabel")

    def fake_run(cmd):
        return subprocess.CompletedProcess(cmd, 2, "", "simulated failure")

    monkeypatch.setattr(tool_runners, "run_subprocess_capture", fake_run)

    with pytest.raises(RuntimeError, match="Open Babel failed"):
        tool_runners.run_openbabel_add_h(
            "/tmp/in.pdb", "/tmp/out.pdb", openbabel_path="/usr/bin/obabel"
        )
