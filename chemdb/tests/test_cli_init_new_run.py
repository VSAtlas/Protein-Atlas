from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from cli import qol_cli


def _run_atlas_cli(repo_root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "cli.atlas_main_cli", *args],
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
    )


@pytest.fixture
def fake_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    (tmp_path / "config.example.txt").write_text(
        "OVERALL_DIR=.\nINPUT_DIR=input_pdbs\n",
        encoding="utf-8",
    )
    input_dir = tmp_path / "input_pdbs"

    def _repo_root() -> Path:
        return tmp_path

    def _load_effective_config() -> dict[str, str]:
        return {"INPUT_DIR": str(input_dir)}

    monkeypatch.setattr(qol_cli, "_repo_root", _repo_root)
    monkeypatch.setattr(qol_cli, "load_effective_config", _load_effective_config)
    return tmp_path


def test_atlas_init_dry_run_hermetic(
    fake_repo: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    forwarded: list[list[str]] = []

    def _fake_configure_local_tools(argv: list[str]) -> int:
        forwarded.append(list(argv))
        return 0

    monkeypatch.setattr(
        "tools.installers.configure_local_tools.main",
        _fake_configure_local_tools,
    )

    rc = qol_cli.dispatch(["init", "--dry-run", "--no-detect"])
    output = capsys.readouterr().out

    assert rc == 0
    assert "Would create" in output
    assert str(fake_repo / "config.txt") in output
    assert "atlas doctor" in output
    assert forwarded == [
        [
            "--config",
            str((fake_repo / "config.txt").resolve()),
            "--example",
            str((fake_repo / "config.example.txt").resolve()),
            "--dry-run",
            "--no-detect",
        ]
    ]
    assert not (fake_repo / "config.txt").exists()


def test_atlas_init_creates_config_from_example(
    fake_repo: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(
        "tools.installers.configure_local_tools.main",
        lambda _argv: 0,
    )

    rc = qol_cli.dispatch(["init", "--no-detect"])
    output = capsys.readouterr().out

    assert rc == 0
    assert "Created" in output
    assert (fake_repo / "config.txt").exists()
    assert (fake_repo / "config.txt").read_text(encoding="utf-8") == (
        fake_repo / "config.example.txt"
    ).read_text(encoding="utf-8")


def test_atlas_init_write_example_inputs_dry_run(
    fake_repo: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    (fake_repo / "config.txt").write_text("OVERALL_DIR=.\n", encoding="utf-8")
    monkeypatch.setattr(
        "tools.installers.configure_local_tools.main",
        lambda _argv: 0,
    )

    rc = qol_cli.dispatch(
        ["init", "--dry-run", "--no-detect", "--write-example-inputs"]
    )
    output = capsys.readouterr().out

    assert rc == 0
    assert "Would write example guide:" in output
    assert not (fake_repo / "docs" / "examples" / "new_user_inputs.md").exists()


def test_atlas_init_write_example_inputs_creates_guide(
    fake_repo: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    (fake_repo / "config.txt").write_text("OVERALL_DIR=.\n", encoding="utf-8")
    monkeypatch.setattr(
        "tools.installers.configure_local_tools.main",
        lambda _argv: 0,
    )

    rc = qol_cli.dispatch(["init", "--no-detect", "--write-example-inputs"])
    output = capsys.readouterr().out
    guide = fake_repo / "docs" / "examples" / "new_user_inputs.md"

    assert rc == 0
    assert "Wrote example guide:" in output
    assert guide.exists()
    assert "atlas new-run --pdb" in guide.read_text(encoding="utf-8")


def test_atlas_new_run_dry_run_input_missing(
    fake_repo: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    (fake_repo / "config.txt").write_text("INPUT_DIR=input_pdbs\n", encoding="utf-8")
    (fake_repo / "input_pdbs").mkdir()

    rc = qol_cli.dispatch(["new-run", "--pdb", "1abc", "--dry-run"])
    output = capsys.readouterr().out

    assert rc == 1
    assert "target: 1ABC" in output
    assert "status: input_missing" in output
    assert str(fake_repo / "input_pdbs" / "1ABC.pdb") in output


def test_atlas_new_run_dry_run_ready(
    fake_repo: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    (fake_repo / "config.txt").write_text("INPUT_DIR=input_pdbs\n", encoding="utf-8")
    input_dir = fake_repo / "input_pdbs"
    input_dir.mkdir()
    (input_dir / "1ABC.pdb").write_text("END\n", encoding="utf-8")

    monkeypatch.setattr(
        qol_cli,
        "_build_setup_report",
        lambda: {
            "ready": {
                "vina_docking": True,
                "full_receptor_prep": True,
            }
        },
    )

    rc = qol_cli.dispatch(
        ["new-run", "--pdb", "1abc", "--small-library", "--fast", "--dry-run"]
    )
    output = capsys.readouterr().out

    assert rc == 0
    assert "target: 1ABC" in output
    assert "status: ready" in output
    assert "command: atlas --pdb 1ABC --fast -test-fda" in output
    assert "--yes --no-dry-run" in output


def test_atlas_init_and_new_run_subprocess_help() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    for args in (("init", "--help"), ("new-run", "--help")):
        proc = _run_atlas_cli(repo_root, *args)
        assert proc.returncode == 0, proc.stderr
        assert "--dry-run" in proc.stdout


def test_run_panel_keeps_ligand_library_separate_from_target_ligand_filter(
    fake_repo: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    captured: list[dict[str, object]] = []

    def _fake_install_targets(_cfg: dict[str, str], **kwargs: object) -> SimpleNamespace:
        captured.append(dict(kwargs))
        return SimpleNamespace(
            panel=kwargs.get("panel_name"),
            dry_run=kwargs.get("dry_run"),
            selected_ids=["1ABC"],
            installed_ids=[],
            failures=[],
            genes_file=None,
            selected_out=fake_repo / "selected.csv",
            output_dir=fake_repo / "input_pdbs",
            manifest_path=fake_repo / "manifest.yaml",
        )

    monkeypatch.setattr(
        "cli.qol.run_panel._bindings.load_effective_config",
        lambda: {"INPUT_DIR": str(fake_repo / "input_pdbs")},
    )
    monkeypatch.setattr("cli.qol.targets._bindings.repo_root", lambda: fake_repo)
    monkeypatch.setattr("cli.target_install.install_targets", _fake_install_targets)
    monkeypatch.setattr(
        "prep_ligands.ligand_library_manager.resolve_source",
        lambda name: SimpleNamespace(name=name),
    )

    rc = qol_cli.dispatch(
        ["run-panel", "kinases", "--ligands", "chembl", "--dry-run"]
    )
    output = capsys.readouterr().out

    assert rc == 0
    assert captured[-1]["ligand_comp_ids"] == ()
    assert "running_with_library: chembl" in output
    assert "--pdb 1ABC" in output

    rc = qol_cli.dispatch(
        [
            "run-panel",
            "kinases",
            "--ligands",
            "chembl",
            "--ligand",
            "ATP",
            "--dry-run",
        ]
    )
    capsys.readouterr()

    assert rc == 0
    assert captured[-1]["ligand_comp_ids"] == ("ATP",)
