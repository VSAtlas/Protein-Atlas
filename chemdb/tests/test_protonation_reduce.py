from pathlib import Path

import protein_prep.protonation_reduce as protonation_reduce


def _atom_line(element: str = "C") -> str:
    return (
        "ATOM      1  CA  GLY A   1      11.000  12.000  13.000"
        f"  1.00 20.00          {element:>2}\n"
    )


def test_pick_reduce_exe_prefers_env_path(tmp_path, monkeypatch):
    fake_reduce = tmp_path / "reduce"
    fake_reduce.write_text("#!/bin/sh\n", encoding="utf-8")
    monkeypatch.setenv("REDUCE_EXE", str(fake_reduce))
    monkeypatch.delenv("REDUCE_BIN", raising=False)

    picked = protonation_reduce._pick_reduce_exe()

    assert picked == str(fake_reduce)


def test_pick_reduce_exe_falls_back_to_path_lookup(tmp_path, monkeypatch):
    monkeypatch.delenv("REDUCE_EXE", raising=False)
    monkeypatch.delenv("REDUCE_BIN", raising=False)
    monkeypatch.setattr(protonation_reduce, "_REPO_ROOT", tmp_path)
    monkeypatch.setattr(protonation_reduce, "PHENIX_DIR", "")
    monkeypatch.setattr(protonation_reduce.shutil, "which", lambda name: "/usr/bin/reduce")

    picked = protonation_reduce._pick_reduce_exe()

    assert picked == "/usr/bin/reduce"


def test_assign_protonation_states_calls_openbabel_on_reduce_failure(
    tmp_path, monkeypatch
):
    in_pdb = tmp_path / "in.pdb"
    out_pdb = tmp_path / "out.pdb"
    in_pdb.write_text(_atom_line("C"), encoding="utf-8")

    monkeypatch.delenv("REDUCE_RETRY", raising=False)
    monkeypatch.delenv("FALLBACK_ADDH", raising=False)

    def fail_reduce(*args, **kwargs):
        raise RuntimeError("reduce failed in test")

    calls = []

    def fake_openbabel(input_pdb, output_pdb, openbabel_path=""):
        calls.append((str(input_pdb), str(output_pdb), openbabel_path))
        Path(output_pdb).write_text(_atom_line("H"), encoding="utf-8")

    monkeypatch.setattr(protonation_reduce.subprocess, "run", fail_reduce)
    monkeypatch.setattr(protonation_reduce, "run_openbabel_add_h", fake_openbabel)
    monkeypatch.setattr(protonation_reduce, "fix_pdb_elements", lambda *_a, **_k: None)

    result = protonation_reduce.assign_protonation_states(
        in_pdb, out_pdb, reduce_exe="/usr/bin/reduce"
    )

    assert result == str(out_pdb.resolve())
    assert calls
    assert out_pdb.exists()


def test_assign_protonation_states_skip_reduce_uses_openbabel(
    tmp_path, monkeypatch
):
    in_pdb = tmp_path / "in_skip.pdb"
    out_pdb = tmp_path / "out_skip.pdb"
    in_pdb.write_text(_atom_line("C"), encoding="utf-8")

    monkeypatch.setenv("SKIP_REDUCE", "1")

    calls = []

    def fake_openbabel(input_pdb, output_pdb, openbabel_path=""):
        calls.append((str(input_pdb), str(output_pdb), openbabel_path))
        Path(output_pdb).write_text(_atom_line("H"), encoding="utf-8")

    monkeypatch.setattr(protonation_reduce, "run_openbabel_add_h", fake_openbabel)

    result = protonation_reduce.assign_protonation_states(in_pdb, out_pdb, reduce_exe=None)

    assert result == str(out_pdb.resolve())
    assert calls
    assert out_pdb.exists()
