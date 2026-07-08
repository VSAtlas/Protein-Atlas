from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any

from config.runtime_config import load_config
from cli.qol import _bindings


def _display_path(value: Any, *, root: Path | None = None) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    path = Path(text).expanduser()
    if not path.is_absolute():
        return text
    root = root or _bindings.repo_root()
    for base, label in (
        (root, "."),
        (Path(sys.prefix), "$PYTHON_PREFIX"),
        (Path(os.environ.get("CONDA_PREFIX", "") or ""), "$CONDA_PREFIX"),
        (Path(os.environ.get("VIRTUAL_ENV", "") or ""), "$VIRTUAL_ENV"),
        (Path.home(), "~"),
    ):
        if not str(base):
            continue
        try:
            rel = path.relative_to(base.expanduser())
        except ValueError:
            continue
        rel_text = rel.as_posix()
        if label == ".":
            return rel_text or "."
        return f"{label}/{rel_text}" if rel_text else label
    if path.parts[:2] in {("/", "bin"), ("/", "usr")}:
        return path.as_posix()
    tail = "/".join(path.parts[-3:])
    return f".../{tail}" if tail else path.name


def load_effective_config() -> dict[str, Any]:
    return load_config("config.txt", base_dir=_bindings.repo_root())
