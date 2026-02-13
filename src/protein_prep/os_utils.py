from __future__ import annotations

import logging
import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import Union


def _is_wsl() -> bool:
    try:
        return "microsoft" in platform.release().lower() or "WSL_INTEROP" in os.environ
    except Exception:
        return False


def _bin_on_path(name: str) -> bool:
    return shutil.which(name) is not None


def _win_path(p: Union[str, Path]) -> str:
    s = str(p)
    if not _is_wsl():
        return s
    try:
        out = subprocess.check_output(["wslpath", "-w", s], text=True).strip()
        return out
    except Exception:
        return s


def _powershell(cmd: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def _as_path(p) -> Path:
    return p if isinstance(p, Path) else Path(p)


def _run_and_log(
    cmd: list[str], *, env: dict, check: bool = False
) -> subprocess.CompletedProcess:
    logging.info("[phenix] exec: %s", " ".join(cmd))
    cp = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env
    )
    if cp.stdout:
        logging.debug("[phenix][stdout]\n%s", cp.stdout)
    if cp.stderr:
        logging.debug("[phenix][stderr]\n%s", cp.stderr)
    if check and cp.returncode != 0:
        raise subprocess.CalledProcessError(cp.returncode, cmd, cp.stdout, cp.stderr)
    return cp
