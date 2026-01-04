from __future__ import annotations

import logging
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import docking
import path_router
from docking_subruns import ProteinDockingContext, resolve_center_box_for_ph
from fallback_recenter import RecenterParams


class _StubPaths:
    def __init__(self, pdb_id: str) -> None:
        self.pdb_id = pdb_id


def test_control_centers_by_ph_calls_per_tag(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ph_tags = ["pH3_0", "pH6_0", "pH9_0"]
    monkeypatch.setattr(path_router, "load_ph_tags", lambda _pdb_id, variant=None: ph_tags)

    rec_paths: dict[str, Path] = {}

    def fake_receptor_file(pdb_id: str, variant=None, ph_tag=None, legacy: bool = False) -> Path:
        path = tmp_path / f"{pdb_id}_{ph_tag}.pdbqt"
        path.write_text("RECEPTOR\n", encoding="utf-8")
        rec_paths[str(ph_tag)] = path
        return path

    monkeypatch.setattr(path_router, "receptor_file", fake_receptor_file)

    calls: list[tuple[str, str | None]] = []
    centers = {
        "pH3_0": (1.0, 2.0, 3.0),
        "pH6_0": (4.0, 5.0, 6.0),
        "pH9_0": (7.0, 8.0, 9.0),
    }

    def fake_select_center(cfg, paths, receptor_pdbqt, logger, *, variant=None, ph_token=None, legacy=False):
        calls.append((receptor_pdbqt, ph_token))
        return centers[str(ph_token)], (24.0, 24.0, 24.0)

    monkeypatch.setattr(docking, "select_center_via_control_redock", fake_select_center)

    cfg = {
        "PH_ENSEMBLE": True,
        "CONTROL_BOX_A": 24.0,
        "BOX_SIZE_MAX_A": 28.0,
        "RUN_ID": "",
    }
    paths = _StubPaths("TEST")
    logger = logging.getLogger("test.control_redock_per_ph")

    center_by_ph, box_by_ph, source_by_ph = docking._control_centers_by_ph(
        cfg,
        paths,
        logger,
        variant_token="HOLO",
        legacy_mode=False,
        cleaned_pdb="dummy",
        fallback_center=(0.0, 0.0, 0.0),
        fallback_box=(10.0, 10.0, 10.0),
    )

    assert len(calls) == len(ph_tags)
    for tag in ph_tags:
        assert (str(rec_paths[tag]), tag) in calls
        assert center_by_ph[tag] == centers[tag]
        assert box_by_ph[tag] == (24.0, 24.0, 24.0)
        assert source_by_ph[tag] == "control"


def test_resolve_center_box_for_ph_uses_per_ph_map() -> None:
    ctx = ProteinDockingContext(
        cfg={},
        paths=None,  # type: ignore[arg-type]
        logger=logging.getLogger("test.control_resolve"),
        pdb_id="TEST",
        variant_env="",
        variant_token=None,
        variant_label="legacy",
        legacy_mode=False,
        cleaned_pdb=None,
        receptor_pdbqt=None,
        center=(0.0, 0.0, 0.0),
        box_size=(10.0, 10.0, 10.0),
        center_by_ph={"pH6_0": (1.0, 1.0, 1.0)},
        box_by_ph={"pH6_0": (20.0, 20.0, 20.0)},
        center_source_by_ph={"pH6_0": "control"},
        stages=[],
        recenter_params=RecenterParams(),
        control_stems=[],
        control_lookup={},
    )

    center, box, src = resolve_center_box_for_ph(ctx, "pH6_0")
    assert center == (1.0, 1.0, 1.0)
    assert box == (20.0, 20.0, 20.0)
    assert src == "control_ph"

    center_fallback, box_fallback, src_fallback = resolve_center_box_for_ph(ctx, "pH8_0")
    assert center_fallback == (0.0, 0.0, 0.0)
    assert box_fallback == (10.0, 10.0, 10.0)
    assert src_fallback == "control_missing"
