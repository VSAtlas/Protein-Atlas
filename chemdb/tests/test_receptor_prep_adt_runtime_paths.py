from pathlib import Path
import subprocess

import protein_prep.receptor_prep as rr


def _write_histidine_pdb(path: Path) -> None:
    path.write_text(
        "ATOM      1  N   HIS A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  HIS A   1      11.000  10.000  10.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )


def _patch_common(monkeypatch, input_pdb: Path, output_pdbqt: Path) -> None:
    monkeypatch.setattr(
        rr,
        "_cfg",
        lambda key, default="", legacy_key=None: default,
    )
    monkeypatch.setattr(rr, "_persist_subproc", lambda *args, **kwargs: None)
    monkeypatch.setattr(rr, "_clean_receptor_pdbqt", lambda *args, **kwargs: None)
    monkeypatch.setattr(rr, "_drop_free_ions_for_meeko", lambda *args, **kwargs: None)
    monkeypatch.setattr(rr, "_load_retain_allowlist", lambda cfg: (set(), "test"))
    monkeypatch.setattr(rr, "_ion_pairs_from_pdb", lambda *_args, **_kwargs: set())
    monkeypatch.setattr(rr, "_ion_pairs_from_pdbqt", lambda *_args, **_kwargs: set())
    monkeypatch.setattr(rr, "_log_ion_diff", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        rr,
        "_log_pdb_pdbqt_counts_diff",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(rr, "_ok_receptor_file", lambda _path: True)
    monkeypatch.setattr(
        rr,
        "_meeko_preflight_or_fail",
        lambda pdb_path, _work_dir: Path(pdb_path),
    )
    monkeypatch.setattr(
        rr,
        "_classify_and_rename_histidines",
        lambda _src, dst: Path(dst).write_text(
            input_pdb.read_text(encoding="utf-8"),
            encoding="utf-8",
        ),
    )
    monkeypatch.setattr(rr, "build_meeko_base_cmd", lambda: ["MEEKO_BASE"])
    monkeypatch.setattr(
        rr,
        "build_meeko_modern_cmd",
        lambda _base, _pdb, _out: ["MEEKO_MODERN"],
    )
    monkeypatch.setattr(
        rr,
        "build_meeko_legacy_cmd",
        lambda _base, _pdb, _out: ["MEEKO_LEGACY"],
    )
    monkeypatch.setattr(
        rr,
        "build_adt_prepare_receptor_cmd",
        lambda py, script, pdb, out: [py, script, "-r", pdb, "-o", out],
    )


def test_runtime_cfg_paths_enable_adt_even_if_module_globals_stale(
    tmp_path, monkeypatch
):
    input_pdb = tmp_path / "input.pdb"
    output_pdbqt = tmp_path / "out" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    _write_histidine_pdb(input_pdb)

    _patch_common(monkeypatch, input_pdb, output_pdbqt)
    monkeypatch.setattr(rr, "MGLTOOLS_PYTHON", "")
    monkeypatch.setattr(rr, "PREPARE_RECEPTOR_SCRIPT", "")
    monkeypatch.delenv("MGLTOOLS_PYTHON", raising=False)
    monkeypatch.delenv("PREPARE_RECEPTOR_SCRIPT", raising=False)

    fake_python = "/fake/mgl/pythonsh"
    fake_script = "/fake/mgl/prepare_receptor4.py"
    cfg = {
        "MGLTOOLS_PYTHON": fake_python,
        "PREPARE_RECEPTOR_SCRIPT": fake_script,
        "MEEKO_ALLOW_BAD_RES": "true",
    }

    original_is_file = Path.is_file

    def fake_is_file(path_obj: Path) -> bool:
        if str(path_obj) in {fake_python, fake_script, str(output_pdbqt)}:
            return True
        return original_is_file(path_obj)

    monkeypatch.setattr(Path, "is_file", fake_is_file)

    original_access = rr.os.access

    def fake_access(path, mode):
        if str(path) == fake_python:
            return True
        return original_access(path, mode)

    monkeypatch.setattr(rr.os, "access", fake_access)

    commands: list[list[str]] = []
    adt_calls = {"count": 0}

    def fake_run(cmd):
        commands.append(list(cmd))
        if cmd and cmd[0] == fake_python:
            adt_calls["count"] += 1
            if adt_calls["count"] == 2:
                return subprocess.CompletedProcess(cmd, 0, "", "")
            return subprocess.CompletedProcess(cmd, 1, "", "adt fail")
        return subprocess.CompletedProcess(cmd, 1, "", "meeko fail")

    monkeypatch.setattr(rr, "run_subprocess_capture", fake_run)

    ok = rr.run_prepare_receptor(str(input_pdb), str(output_pdbqt), cfg)

    assert ok is True
    adt_cmds = [cmd for cmd in commands if cmd and cmd[0] == fake_python]
    assert adt_cmds
    assert any(fake_script in cmd for cmd in adt_cmds)


def test_invalid_runtime_paths_disable_adt(tmp_path, monkeypatch):
    input_pdb = tmp_path / "input.pdb"
    output_pdbqt = tmp_path / "out" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    _write_histidine_pdb(input_pdb)

    _patch_common(monkeypatch, input_pdb, output_pdbqt)
    monkeypatch.setattr(rr, "MGLTOOLS_PYTHON", "")
    monkeypatch.setattr(rr, "PREPARE_RECEPTOR_SCRIPT", "")
    monkeypatch.delenv("MGLTOOLS_PYTHON", raising=False)
    monkeypatch.delenv("PREPARE_RECEPTOR_SCRIPT", raising=False)
    monkeypatch.setattr(rr, "_ok_receptor_file", lambda _path: False)

    invalid_python = "/fake/invalid/pythonsh"
    invalid_script = "/fake/invalid/prepare_receptor4.py"
    cfg = {
        "MGLTOOLS_PYTHON": invalid_python,
        "PREPARE_RECEPTOR_SCRIPT": invalid_script,
        "MEEKO_ALLOW_BAD_RES": "true",
    }

    original_is_file = Path.is_file

    def fake_is_file(path_obj: Path) -> bool:
        if str(path_obj) in {invalid_python, invalid_script}:
            return False
        return original_is_file(path_obj)

    monkeypatch.setattr(Path, "is_file", fake_is_file)
    monkeypatch.setattr(rr.os, "access", lambda _path, _mode: False)

    commands: list[list[str]] = []

    def fake_run(cmd):
        commands.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 1, "", "meeko fail")

    monkeypatch.setattr(rr, "run_subprocess_capture", fake_run)

    ok = rr.run_prepare_receptor(str(input_pdb), str(output_pdbqt), cfg)

    assert ok is False
    assert not any(cmd and cmd[0] == invalid_python for cmd in commands)


def test_env_overrides_cfg_for_adt_paths(tmp_path, monkeypatch):
    input_pdb = tmp_path / "input.pdb"
    output_pdbqt = tmp_path / "out" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    _write_histidine_pdb(input_pdb)

    _patch_common(monkeypatch, input_pdb, output_pdbqt)
    monkeypatch.setattr(rr, "MGLTOOLS_PYTHON", "")
    monkeypatch.setattr(rr, "PREPARE_RECEPTOR_SCRIPT", "")

    cfg_python = "/fake/cfg/pythonsh"
    cfg_script = "/fake/cfg/prepare_receptor4.py"
    env_python = "/fake/env/pythonsh"
    env_script = "/fake/env/prepare_receptor4.py"
    cfg = {
        "MGLTOOLS_PYTHON": cfg_python,
        "PREPARE_RECEPTOR_SCRIPT": cfg_script,
        "MEEKO_ALLOW_BAD_RES": "true",
    }
    monkeypatch.setenv("MGLTOOLS_PYTHON", env_python)
    monkeypatch.setenv("PREPARE_RECEPTOR_SCRIPT", env_script)

    original_is_file = Path.is_file

    def fake_is_file(path_obj: Path) -> bool:
        if str(path_obj) in {cfg_python, cfg_script, env_python, env_script}:
            return True
        return original_is_file(path_obj)

    monkeypatch.setattr(Path, "is_file", fake_is_file)

    original_access = rr.os.access

    def fake_access(path, mode):
        if str(path) in {env_python, cfg_python}:
            return True
        return original_access(path, mode)

    monkeypatch.setattr(rr.os, "access", fake_access)

    commands: list[list[str]] = []
    adt_calls = {"count": 0}

    def fake_run(cmd):
        commands.append(list(cmd))
        if cmd and cmd[0] == env_python:
            adt_calls["count"] += 1
            if adt_calls["count"] == 2:
                return subprocess.CompletedProcess(cmd, 0, "", "")
            return subprocess.CompletedProcess(cmd, 1, "", "adt fail")
        return subprocess.CompletedProcess(cmd, 1, "", "meeko fail")

    monkeypatch.setattr(rr, "run_subprocess_capture", fake_run)

    ok = rr.run_prepare_receptor(str(input_pdb), str(output_pdbqt), cfg)

    assert ok is True
    adt_cmds = [cmd for cmd in commands if cmd and cmd[0] == env_python]
    assert adt_cmds
    assert any(env_script in cmd for cmd in adt_cmds)
    assert not any(cmd and cmd[0] == cfg_python for cmd in commands)
