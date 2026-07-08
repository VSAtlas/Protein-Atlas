from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

# Ensure we can import project modules when this test file lives in chemdb/tests
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import prep_docking.prep_dock6 as prep_dock6
import prep_docking.prep_for_ledock as prep_for_ledock


MINIMAL_WITHH_PDB = """\
ATOM      1  N   ALA A   1      11.104  13.207   9.597  1.00 20.00           N
ATOM      2  H   ALA A   1      11.104  13.207   9.597  1.00 20.00           H
TER
END
"""


def _write_withh(ensemble_dir: Path, names: list[str]) -> None:
    for name in names:
        (ensemble_dir / name).write_text(MINIMAL_WITHH_PDB)


@pytest.fixture()
def ensemble_dirs(tmp_path: Path) -> dict[str, Path]:
    members = [
        {
            "tag": "2OJG_pH7_2",
            "withH": "2OJG_pH7_2.withH.pdb",
            "pdbqt": "2OJG_pH7_2+7_7-dup5.pdbqt",
        },
        {
            "pdbqt": "2OJG_pH6_5+6_0-dup1.pdbqt",
            "withH": "2OJG_pH6_5.withH.pdb",
        },
        {
            "label": "pH8_0+8_5-dup3",
            "withH": "2OJG_pH8_0.withH.pdb",
        },
    ]
    withh_names = [
        "2OJG_pH7_2.withH.pdb",
        "2OJG_pH6_5.withH.pdb",
        "2OJG_pH8_0.withH.pdb",
    ]
    out: dict[str, Path] = {}
    for variant in ("HOLO", "APO"):
        ensemble_dir = tmp_path / variant / "receptor" / "ph_ensemble"
        ensemble_dir.mkdir(parents=True, exist_ok=True)
        (ensemble_dir / "ensemble.json").write_text(json.dumps({"members": members}))
        _write_withh(ensemble_dir, withh_names)
        out[variant] = ensemble_dir
    return out


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("pH7_2+7_7-dup5", "pH7_2"),
        ("2OJG_pH7_2+7_7-dup5", "pH7_2"),
        ("2OJG_pH7_2.withH.pdb", "pH7_2"),
        ("2OJG_pH6_5+6_0-dup1.pdbqt", "pH6_5"),
        ("pH8_0+8_5-dup3", "pH8_0"),
        ("pH7_2", "pH7_2"),
    ],
)
def test_canonicalize_ph_key_variants(raw: str, expected: str) -> None:
    assert prep_for_ledock._canonicalize_ph_key(raw, "2OJG") == expected


@pytest.mark.parametrize("variant", ["HOLO", "APO"])
@pytest.mark.parametrize(
    "query, expected_name",
    [
        ("pH7_2", "2OJG_pH7_2.withH.pdb"),
        ("pH7_2+7_7", "2OJG_pH7_2.withH.pdb"),
        ("pH7_2+7_7-dup5", "2OJG_pH7_2.withH.pdb"),
        ("2OJG_pH7_2+7_7-dup5", "2OJG_pH7_2.withH.pdb"),
        ("2OJG_pH7_2.withH.pdb", "2OJG_pH7_2.withH.pdb"),
        ("2OJG_pH7_2+7_7-dup5.pdbqt", "2OJG_pH7_2.withH.pdb"),
        ("pH6_5+6_0-dup1", "2OJG_pH6_5.withH.pdb"),
        ("2OJG_pH6_5+6_0-dup1.pdbqt", "2OJG_pH6_5.withH.pdb"),
        ("pH8_0+8_5-dup3", "2OJG_pH8_0.withH.pdb"),
    ],
)
def test_resolve_withh_from_manifest_accepts_degenerate_labels(
    ensemble_dirs: dict[str, Path], variant: str, query: str, expected_name: str
) -> None:
    ensemble_dir = ensemble_dirs[variant]
    logger = logging.getLogger("test.ensemble")
    resolved = prep_for_ledock._resolve_withh_from_manifest(
        ensemble_dir, "2OJG", query, logger
    )
    assert resolved == ensemble_dir / expected_name
    assert resolved.is_absolute()


@pytest.mark.parametrize("variant", ["HOLO", "APO"])
def test_ensure_ledock_receptor_uses_canonical_withh(
    ensemble_dirs: dict[str, Path], monkeypatch: pytest.MonkeyPatch, variant: str
) -> None:
    calls: list[dict[str, object]] = []

    def fake_ph_ensemble_dir(
        pdb_id: str, variant: str | None = None, legacy: bool = False
    ) -> Path:
        return ensemble_dirs[variant or "HOLO"]

    def fake_run(cmd, check=False, cwd=None, **_kwargs):
        calls.append({"cmd": cmd, "cwd": cwd})
        if cwd:
            Path(cwd, "pro.pdb").write_text(MINIMAL_WITHH_PDB)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(prep_for_ledock, "ph_ensemble_dir", fake_ph_ensemble_dir)
    monkeypatch.setattr(prep_for_ledock.subprocess, "run", fake_run)

    logger = logging.getLogger("test.ledock")
    out = prep_for_ledock.ensure_ledock_receptor(
        cfg={},
        pdb_id="2OJG",
        variant=variant,
        ph_label="2OJG_pH7_2+7_7-dup5",
        logger=logger,
    )

    assert out is not None
    assert out.name == "pro.pdb"
    assert out.exists()
    assert calls
    assert calls[0]["cmd"][1] == "2OJG_pH7_2.withH.pdb"


def test_ensure_dock6_surface_selects_base_withh(
    ensemble_dirs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[dict[str, object]] = []

    def fake_ph_ensemble_dir(
        pdb_id: str, variant: str | None = None, legacy: bool = False
    ) -> Path:
        return ensemble_dirs[variant or "HOLO"]

    def fake_run(cmd, check=False, cwd=None, **_kwargs):
        calls.append({"cmd": cmd, "cwd": cwd})
        if "-o" in cmd:
            out_path = Path(cmd[cmd.index("-o") + 1])
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text("MS")
        return SimpleNamespace(returncode=0)

    def fake_ensure_dock6_grids(
        cfg,
        dock6_root: Path,
        receptor_pdb: Path,
        logger,
        grid_prefix: str,
        grid_spacing: float,
    ) -> Path:
        grid_path = dock6_root / f"{grid_prefix}.nrg"
        grid_path.write_text("GRID")
        return grid_path

    monkeypatch.setattr(prep_dock6, "ph_ensemble_dir", fake_ph_ensemble_dir)
    monkeypatch.setattr(prep_dock6.subprocess, "run", fake_run)
    monkeypatch.setattr(prep_dock6, "_ensure_sphgen_spheres", lambda **_kwargs: None)
    monkeypatch.setattr(prep_dock6, "_ensure_selected_spheres", lambda **_kwargs: None)
    monkeypatch.setattr(prep_dock6, "_ensure_site_box", lambda **_kwargs: None)
    monkeypatch.setattr(prep_dock6, "_log_site_alignment", lambda **_kwargs: None)
    monkeypatch.setattr(prep_dock6, "_ensure_dock6_grids", fake_ensure_dock6_grids)

    logger = logging.getLogger("test.dock6")
    out = prep_dock6.ensure_dock6_surface(
        cfg={"FAST_MODE": True},
        pdb_id="2OJG",
        variant="HOLO",
        ph_label="pH7_2+7_7-dup5",
        logger=logger,
    )

    assert out is not None
    assert out.name == "rec.ms"
    assert out.exists()

    dock6_root = ensemble_dirs["HOLO"] / "dock6"
    expected_noh = dock6_root / "2OJG_pH7_2_noH.pdb"
    assert expected_noh.exists()

    assert calls
    assert calls[0]["cmd"][1] == str(expected_noh)
