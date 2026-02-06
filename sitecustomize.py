"""
Project-level sitecustomize to keep sandboxed micromamba runs working during tests.
Really this is just so that codex agent can run tests

If the default HOME is not writable (common in restricted sandboxes), micromamba
fails to acquire ~/.cache/mamba/proc/proc.lock. We redirect HOME to a repo-local
fallback that we can write to, but only when the existing HOME looks unusable.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _is_writable(path: Path) -> bool:
    try:
        path.mkdir(parents=True, exist_ok=True)
        probe = path / ".sitecustomize_probe"
        probe.write_text("ok")
        probe.unlink()
        return True
    except Exception:
        return False


def _ensure_writable_home() -> None:
    home = Path(os.environ.get("HOME", "") or "").expanduser()
    cache_root = home / ".cache" / "mamba"

    if home and _is_writable(cache_root):
        return

    fallback = Path(__file__).resolve().parent / ".home"
    fallback_cache = fallback / ".cache" / "mamba"
    if _is_writable(fallback_cache):
        os.environ["HOME"] = str(fallback)


_ensure_writable_home()

# Ensure the repository root and src directory are always importable for tests and scripts.
REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))
