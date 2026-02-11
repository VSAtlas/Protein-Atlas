from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def test_select_control_engines_defaults_to_vina_only():
    import docking.docking_control_redock as docking

    cfg = {"CONTROL_CONSENSUS": False}
    engines = docking._select_control_engines(
        cfg, use_gnina=True, use_ledock=True, use_dock6=True
    )
    assert engines == ["vina"]


def test_select_control_engines_adds_engines_when_enabled():
    import docking.docking_control_redock as docking

    cfg = {"CONTROL_CONSENSUS": True}
    engines = docking._select_control_engines(
        cfg, use_gnina=True, use_ledock=True, use_dock6=True
    )
    assert engines == ["vina", "gnina", "ledock", "dock6"]


def test_select_control_engines_respects_disabled_capabilities():
    import docking.docking_control_redock as docking

    cfg = {"CONTROL_CONSENSUS": True}
    engines = docking._select_control_engines(
        cfg, use_gnina=False, use_ledock=True, use_dock6=False
    )
    assert engines == ["vina", "ledock"]
