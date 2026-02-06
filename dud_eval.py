"""
Shim module to expose analysis/dud_eval.py as a top-level import.
"""

from analysis.dud_eval import *  # type: ignore  # noqa: F401,F403
from analysis.dud_eval import main as _main  # type: ignore[attr-defined]


if __name__ == "__main__":
    _main()
