from __future__ import annotations

import sys
from typing import Sequence


from cli.qol.onboarding import (
    _cmd_demo,
    _cmd_help,
    _cmd_init,
    _cmd_setup_report,
    _cmd_smoke,
)
from cli.qol.new_run import _cmd_first_run, _cmd_new_run
from cli.qol.status import _cmd_runs, _cmd_status
from cli.qol.reporting import _cmd_analysis, _cmd_report
from cli.qol.maintenance import (
    _cmd_artifacts,
    _cmd_debug,
    _cmd_dev,
    _cmd_reproduce,
    _cmd_screenshot,
    _cmd_throughput,
)
from cli.qol.slurm import _cmd_slurm
from cli.qol.targets import _cmd_targets
from cli.qol.run_panel import _cmd_run_panel
from cli.qol.ligands import _cmd_ligands, _cmd_library_alias
from cli.qol.ml import _cmd_ml
from cli.qol.pipeline import _cmd_run
from cli.qol.publish import _cmd_publish


def dispatch(argv: Sequence[str] | None = None) -> int | None:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        return None
    command = args[0].strip().lower()
    handlers = {
        "help": _cmd_help,
        "init": _cmd_init,
        "demo": _cmd_demo,
        "smoke": _cmd_smoke,
        "new-run": _cmd_new_run,
        "first-run": _cmd_first_run,
        "setup-report": _cmd_setup_report,
        "status": _cmd_status,
        "runs": _cmd_runs,
        "report": _cmd_report,
        "analysis": _cmd_analysis,
        "debug": _cmd_debug,
        "artifacts": _cmd_artifacts,
        "throughput": _cmd_throughput,
        "screenshot": _cmd_screenshot,
        "dev": _cmd_dev,
        "reproduce": _cmd_reproduce,
        "slurm": _cmd_slurm,
        "targets": _cmd_targets,
        "run-panel": _cmd_run_panel,
        "ligands": _cmd_ligands,
        "ml": _cmd_ml,
        "hmdb": lambda rest: _cmd_library_alias("hmdb", rest),
        "fda": lambda rest: _cmd_library_alias("fda", rest),
        "chembl": lambda rest: _cmd_library_alias("chembl", rest),
        "chebi": lambda rest: _cmd_library_alias("chebi", rest),
        "coconut": lambda rest: _cmd_library_alias("coconut", rest),
        "run": _cmd_run,
        "publish": _cmd_publish,
    }
    handler = handlers.get(command)
    if handler is None:
        return None
    return handler(args[1:])
