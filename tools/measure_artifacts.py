#!/usr/bin/env python3
"""Compatibility shim for root-level measure_artifacts.py."""

from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    script = Path(__file__).resolve().parents[1] / "measure_artifacts.py"
    runpy.run_path(str(script), run_name="__main__")
