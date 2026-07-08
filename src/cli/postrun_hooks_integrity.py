# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
import os
import shlex
import subprocess
import sys
from typing import Any, Mapping

from cli.postrun_hooks_common import _is_no_library_docking
from cli.postrun_hooks_support import _resolve_hook_roots, _with_repo_src_on_pythonpath


def _maybe_run_throughput_integrity(
    cfg: Mapping[str, Any], run_id: str, *, strict: bool = True
) -> bool:
    """
    Run throughput integrity accounting and return True on pass.
    """
    logger = logging.getLogger("throughput-integrity-hook")
    if not run_id:
        logger.info("[throughput-integrity.skip] reason=missing_run_id")
        return True
    if _is_no_library_docking(cfg):
        logger.info("[throughput-integrity.skip] reason=no_library_docking run_id=%s", run_id)
        return True

    roots = _resolve_hook_roots(cfg)
    cmd = [
        sys.executable,
        "-m",
        "analysis.cli.throughput_integrity",
        "--run-id",
        str(run_id),
        "--repo-root",
        str(roots.analysis_root),
    ]
    if strict:
        cmd.append("--strict")
    else:
        cmd.append("--no-strict")

    logger.info(
        "[throughput-integrity.run] run_id=%s strict=%s cmd=%s",
        run_id,
        str(strict).lower(),
        shlex.join(cmd),
    )
    try:
        env = _with_repo_src_on_pythonpath(dict(os.environ), roots.code_root)
        result = subprocess.run(cmd, cwd=str(roots.code_root), check=False, env=env)
    except Exception:
        logger.warning(
            "[throughput-integrity.fail] run_id=%s reason=unexpected_exception",
            run_id,
            exc_info=True,
        )
        return False

    if result.returncode != 0:
        logger.warning(
            "[throughput-integrity.fail] run_id=%s returncode=%s strict=%s",
            run_id,
            result.returncode,
            str(strict).lower(),
        )
        return False

    logger.info(
        "[throughput-integrity.done] run_id=%s returncode=%s strict=%s",
        run_id,
        result.returncode,
        str(strict).lower(),
    )
    return True
