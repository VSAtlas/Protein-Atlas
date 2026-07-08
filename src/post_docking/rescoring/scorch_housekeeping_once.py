from __future__ import annotations

import json
import logging
import os
import time
from pathlib import Path
from typing import Callable

from config.output_paths import output_root


def _marker_dir(repo_root: Path, run_id: str) -> Path:
    legacy_root = Path(repo_root).expanduser() / "manifests"
    preferred_root = output_root(repo_root, "manifests")
    if legacy_root.exists() and not preferred_root.exists():
        root = legacy_root
    elif not (Path(repo_root).expanduser() / "outputs").exists():
        root = legacy_root
    else:
        root = preferred_root
    return root / str(run_id) / "scorch_housekeeping"


def _write_done(done_path: Path) -> None:
    done_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = done_path.with_suffix(done_path.suffix + ".tmp")
    tmp.write_text(json.dumps({"done_at": time.time()}, sort_keys=True), encoding="utf-8")
    tmp.replace(done_path)


def run_housekeeping_once(
    *,
    repo_root: Path,
    run_id: str,
    step: str,
    action: Callable[[], None],
    logger: logging.Logger,
    wait_timeout_sec: float = 900.0,
    poll_sec: float = 2.0,
) -> bool:
    marker_dir = _marker_dir(repo_root, run_id)
    marker_dir.mkdir(parents=True, exist_ok=True)
    token = str(step).strip().lower().replace(" ", "_")
    done_path = marker_dir / f"{token}.done.json"
    claim_path = marker_dir / f"{token}.claim"
    if done_path.exists():
        return True

    try:
        fd = os.open(str(claim_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        deadline = time.time() + max(5.0, float(wait_timeout_sec))
        while time.time() < deadline:
            if done_path.exists():
                return True
            time.sleep(max(0.2, float(poll_sec)))
        logger.warning(
            "[scorch.housekeeping.once] action=wait_timeout step=%s run_id=%s wait_timeout_sec=%.1f",
            token,
            run_id,
            float(wait_timeout_sec),
        )
        return False

    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(f"{os.getpid()} {time.time():.3f}\n")
            handle.flush()
            os.fsync(handle.fileno())
        action()
        _write_done(done_path)
        return True
    finally:
        try:
            claim_path.unlink()
        except FileNotFoundError:
            pass
