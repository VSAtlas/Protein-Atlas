from __future__ import annotations

import sys
from typing import Sequence


_PAIR_FLAGS = {
    "--pairs",
    "--selected-pairs",
    "--pair-manifest",
    "--library-manifest",
}
_COMPARE_FLAGS = {"--compare-run", "--compare-runid"}


def _has_any_flag(argv: Sequence[str], flags: set[str]) -> bool:
    return any(argument.split("=", 1)[0] in flags for argument in argv)


def _is_compare_addon_run(argv: Sequence[str]) -> bool:
    return _has_any_flag(argv, _PAIR_FLAGS) and _has_any_flag(argv, _COMPARE_FLAGS)


def _cmd_run(argv: Sequence[str]) -> int:
    if _is_compare_addon_run(argv):
        from cli.qol.ml import _cmd_ml

        return _cmd_ml(["score-addons", *argv])

    from main import main as root_main

    old_argv = sys.argv[:]
    try:
        sys.argv = [old_argv[0], *argv]
        root_main()
    finally:
        sys.argv = old_argv
    return 0
