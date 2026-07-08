#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ruff: noqa: F401,E402,F403

"""Compatibility wrapper for DUD evaluation.

This keeps `analysis/dud_eval.py` executable for existing workflows while
routing imports and behavior to the canonical package layout.
"""

from pathlib import Path
import sys
import warnings

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
if str(_REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "src"))

warnings.warn(
    "analysis.dud_eval is deprecated; use analysis.cli.dud_eval or analysis.dud_eval_core.*. Removal after 2 release cycles.",
    DeprecationWarning,
    stacklevel=2,
)

from analysis.dud_eval_core.api import *
from analysis.dud_eval_core.orchestrate import main as _main


if __name__ == "__main__":
    raise SystemExit(_main())
