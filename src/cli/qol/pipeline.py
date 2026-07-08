from __future__ import annotations

import sys
from typing import Sequence



def _cmd_run(argv: Sequence[str]) -> int:
    from main import main as root_main

    old_argv = sys.argv[:]
    try:
        sys.argv = [old_argv[0], *argv]
        root_main()
    finally:
        sys.argv = old_argv
    return 0
