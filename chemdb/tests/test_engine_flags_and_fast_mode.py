import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from docking.docking_gnina import _resolve_gnina_stage_params, should_run_gnina_for_target
from docking.docking_ledock import resolve_ledock_stage_params, should_run_ledock_for_target


def _hard_td():
    return types.SimpleNamespace(difficulty="hard")


def _degenerate_td():
    return types.SimpleNamespace(difficulty="degenerate")


def _easy_td():
    return types.SimpleNamespace(difficulty="easy")


def _base_ledock_cfg(**overrides):
    cfg = {
        "USE_LEDOCK": "true",
        "PH_ENSEMBLE": True,
        "_PH_ENSEMBLE_CANONICAL": {"TEST": {"pH7_0": {}}},
    }
    cfg.update(overrides)
    return cfg


def test_gnina_stage_params_respect_fast_mode_off():
    base_stage = {"name": "stage1", "exhaustiveness": 8, "num_modes": 9}
    cfg = {"FAST_MODE": False, "VINA_VERBOSITY": 0}
    resolved = _resolve_gnina_stage_params(cfg, base_stage)
    assert resolved["exhaustiveness"] == 8
    assert resolved["num_modes"] == 9
    assert resolved["verbosity"] == 0


def test_gnina_stage_params_respect_fast_mode_on():
    base_stage = {"name": "stage1", "exhaustiveness": 8, "num_modes": 9}
    cfg = {"FAST_MODE": True, "VINA_VERBOSITY": 2}
    resolved = _resolve_gnina_stage_params(cfg, base_stage)
    assert resolved["exhaustiveness"] == 1
    assert resolved["num_modes"] == 1
    assert resolved["verbosity"] == 2


def test_should_run_gnina_disabled_by_use_gnina_false():
    cfg = {
        "USE_GNINA": "false",
        "ENABLE_GNINA": True,
        "GNINA_EXE": "/bin/true",
        "_FILE_CFG": {},
    }
    assert not should_run_gnina_for_target(_hard_td(), cfg)


def test_should_run_gnina_enabled_for_hard_or_degenerate():
    cfg = {
        "USE_GNINA": "true",
        "ENABLE_GNINA": True,
        "GNINA_EXE": "/bin/true",
        "_FILE_CFG": {},
    }
    assert should_run_gnina_for_target(_hard_td(), cfg)
    assert should_run_gnina_for_target(_degenerate_td(), cfg)


def test_should_run_gnina_not_run_for_easy_even_if_enabled():
    cfg = {
        "USE_GNINA": "true",
        "ENABLE_GNINA": True,
        "GNINA_EXE": "/bin/true",
        "_FILE_CFG": {},
    }
    assert not should_run_gnina_for_target(_easy_td(), cfg)


def test_should_run_gnina_requires_gnina_exe_and_enable_flag():
    td = _hard_td()
    cfg_no_exe = {
        "USE_GNINA": "true",
        "ENABLE_GNINA": True,
        "GNINA_EXE": "",
        "_FILE_CFG": {},
    }
    assert not should_run_gnina_for_target(td, cfg_no_exe)

    cfg_disabled = {
        "USE_GNINA": "true",
        "ENABLE_GNINA": False,
        "GNINA_EXE": "/bin/true",
        "_FILE_CFG": {},
    }
    assert not should_run_gnina_for_target(td, cfg_disabled)


def test_ledock_stage_params_default_vs_fast_mode():
    cfg_slow = {"FAST_MODE": False}
    key_slow, params_slow = resolve_ledock_stage_params(cfg_slow, "stage1")
    assert key_slow == "stage1"
    assert params_slow["n_poses"] > 1

    cfg_fast = {"FAST_MODE": True}
    key_fast, params_fast = resolve_ledock_stage_params(cfg_fast, "stage1")
    assert key_fast == "stage1"
    assert params_fast["n_poses"] == 1
    assert pytest.approx(params_fast["rmsd"], rel=0, abs=1e-9) == 1.5


def test_should_run_ledock_disabled_by_use_ledock_false():
    cfg = _base_ledock_cfg(USE_LEDOCK="false")
    assert not should_run_ledock_for_target(cfg)


def test_should_run_ledock_requires_ph_ensemble_flag():
    cfg = _base_ledock_cfg(PH_ENSEMBLE=False)
    assert not should_run_ledock_for_target(cfg)


def test_should_run_ledock_requires_non_empty_manifest():
    cfg = {
        "USE_LEDOCK": "true",
        "PH_ENSEMBLE": True,
        "_PH_ENSEMBLE_CANONICAL": {},
    }
    assert not should_run_ledock_for_target(cfg)


def test_should_run_ledock_true_when_all_flags_set():
    cfg = _base_ledock_cfg()
    assert should_run_ledock_for_target(cfg)
