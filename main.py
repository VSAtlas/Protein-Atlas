#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compatibility bridge for the Atlas docking pipeline.

No new implementation logic belongs in this root module. The maintained
pipeline implementation lives in :mod:`cli.main_pipeline`; this file exists so
``python main.py`` and historical ``import main`` call sites keep working.
"""

from __future__ import annotations

import sys
from pathlib import Path


_REPO_ROOT = Path(__file__).resolve().parent
_SRC_ROOT = _REPO_ROOT / "src"
if str(_SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(_SRC_ROOT))
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from cli import main_pipeline as _pipeline  # noqa: E402

# Make `import main` behave like the historical implementation module, including
# private helper access used by compatibility tests and developer scripts.
sys.modules[__name__] = _pipeline


if __name__ == "__main__":
    _pipeline._run_with_email_notification()
