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
    monkeypatch.setattr(rr, "_ok_receptor_file", lambda path: Path(path).is_file())
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


def test_runtime_cfg_adt_paths_do_not_enable_legacy_fallback(
    tmp_path, monkeypatch
):
    input_pdb = tmp_path / "input.pdb"
    output_pdbqt = tmp_path / "out" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    _write_histidine_pdb(input_pdb)

    _patch_common(monkeypatch, input_pdb, output_pdbqt)
    monkeypatch.delenv("MGLTOOLS_PYTHON", raising=False)
    monkeypatch.delenv("PREPARE_RECEPTOR_SCRIPT", raising=False)

    fake_python = "/fake/mgl/pythonsh"
    fake_script = "/fake/mgl/prepare_receptor4.py"
    cfg = {
        "MGLTOOLS_PYTHON": fake_python,
        "PREPARE_RECEPTOR_SCRIPT": fake_script,
        "MEEKO_ALLOW_BAD_RES": "true",
    }

    commands: list[list[str]] = []

    def fake_run(cmd):
        commands.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 1, "", "meeko fail")

    monkeypatch.setattr(rr, "run_subprocess_capture", fake_run)

    ok = rr.run_prepare_receptor(str(input_pdb), str(output_pdbqt), cfg)

    assert ok is False
    assert commands == [["MEEKO_MODERN"], ["MEEKO_LEGACY"]]
    assert not any(cmd and cmd[0] == fake_python for cmd in commands)
    assert not hasattr(rr, "MGLTOOLS_PYTHON")
    assert not hasattr(rr, "PREPARE_RECEPTOR_SCRIPT")


def test_invalid_runtime_paths_disable_adt(tmp_path, monkeypatch):
    input_pdb = tmp_path / "input.pdb"
    output_pdbqt = tmp_path / "out" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    _write_histidine_pdb(input_pdb)

    _patch_common(monkeypatch, input_pdb, output_pdbqt)
    monkeypatch.delenv("MGLTOOLS_PYTHON", raising=False)
    monkeypatch.delenv("PREPARE_RECEPTOR_SCRIPT", raising=False)

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


def test_env_adt_paths_are_ignored_by_meeko_receptor_prep(tmp_path, monkeypatch):
    input_pdb = tmp_path / "input.pdb"
    output_pdbqt = tmp_path / "out" / "receptor.pdbqt"
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    _write_histidine_pdb(input_pdb)

    _patch_common(monkeypatch, input_pdb, output_pdbqt)

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

    commands: list[list[str]] = []

    def fake_run(cmd):
        commands.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 1, "", "meeko fail")

    monkeypatch.setattr(rr, "run_subprocess_capture", fake_run)

    ok = rr.run_prepare_receptor(str(input_pdb), str(output_pdbqt), cfg)

    assert ok is False
    assert commands == [["MEEKO_MODERN"], ["MEEKO_LEGACY"]]
    assert not any(cmd and cmd[0] == env_python for cmd in commands)
    assert not any(cmd and cmd[0] == cfg_python for cmd in commands)
