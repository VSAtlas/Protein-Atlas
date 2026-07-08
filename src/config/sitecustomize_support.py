"""Support for the repository-level :mod:`sitecustomize` compatibility hook."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _is_writable(path: Path) -> bool:
    try:
        path.mkdir(parents=True, exist_ok=True)
        probe = path / ".sitecustomize_probe"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()
        return True
    except Exception:
        return False


def ensure_sitecustomize(repo_root: Path) -> None:
    """Keep local test/script imports and micromamba cache writes predictable."""
    home = Path(os.environ.get("HOME", "") or "").expanduser()
    cache_root = home / ".cache" / "mamba"

    if not home or not _is_writable(cache_root):
        fallback = repo_root / ".home"
        fallback_cache = fallback / ".cache" / "mamba"
        if _is_writable(fallback_cache):
            os.environ["HOME"] = str(fallback)

    src_root = repo_root / "src"
    for path in (src_root, repo_root):
        path_s = str(path)
        if path_s not in sys.path:
            sys.path.insert(0, path_s)
