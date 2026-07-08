"""Compatibility wrapper for Python's top-level ``sitecustomize`` hook."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"

src_s = str(SRC_ROOT)
if src_s not in sys.path:
    sys.path.insert(0, src_s)

from config.sitecustomize_support import ensure_sitecustomize  # noqa: E402

ensure_sitecustomize(REPO_ROOT)
