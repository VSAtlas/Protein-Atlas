from __future__ import annotations

import os
import signal
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any, Callable, Mapping, Optional, Sequence

from config.normalize import _to_bool

_LOCK_DISABLE_ENV = "ATLAS_SCORCH_DISABLE_LOCK_FALLBACK"
_RUNTIME_DISABLE_ENV = "ATLAS_SCORCH_DISABLE_RUNTIME_RETRY"
_TRUE_VALUES = {"1", "true", "yes", "on", "y"}


def lock_fallback_enabled() -> bool:
    return not _to_bool(os.environ.get(_LOCK_DISABLE_ENV), False)


def runtime_retry_enabled() -> bool:
    return not _to_bool(os.environ.get(_RUNTIME_DISABLE_ENV), False)


def is_mamba_lock_error(stderr: str, stdout: str) -> bool:
    text = f"{stderr}\n{stdout}".lower()
    has_mamba = "libmamba" in text or "micromamba" in text or "mamba" in text
    if not has_mamba or "lock" not in text:
        return False
    signatures = (
        "could not set lock",
        "lockfile acquisition failed",
        "failed to lock",
        "resource temporarily unavailable",
    )
    return any(sig in text for sig in signatures)


def is_runtime_pool_error(stderr: str, stdout: str) -> bool:
    text = f"{stderr}\n{stdout}".lower()
    signatures = (
        "brokenprocesspool",
        "a task has failed to un-serialize",
        "memoryerror",
        "cannot recover from memoryerrors",
    )
    return any(sig in text for sig in signatures)


def build_env_python_retry_cmd(
    cmd: Sequence[str],
    *,
    scorch_python: Optional[Path],
    scorch_script: Optional[Path],
) -> Optional[list[str]]:
    parts = [str(x) for x in cmd]
    if len(parts) < 3 or parts[0] != "micromamba" or parts[1] != "run":
        return None
    script = str(scorch_script) if scorch_script is not None else ""
    script_idx = -1
    if script:
        try:
            script_idx = parts.index(script)
        except ValueError:
            script_idx = -1
    if script_idx < 0:
        for i, token in enumerate(parts):
            if token.endswith("scorch.py"):
                script_idx = i
                script = token
                break
    if script_idx < 0 or not script:
        return None
    python_path = scorch_python
    if python_path is None:
        for flag in ("-p", "--prefix"):
            if flag in parts:
                idx = parts.index(flag)
                if idx + 1 < len(parts):
                    python_path = Path(parts[idx + 1]) / "bin" / "python"
                    break
    python_bin = str(
        python_path if python_path is not None else Path(sys.executable).resolve()
    )
    return [python_bin, script, *parts[script_idx + 1 :]]


def _rewrite_threads(cmd: Sequence[str], threads: int) -> list[str]:
    parts = [str(x) for x in cmd]
    try:
        idx = parts.index("--threads")
    except ValueError:
        return parts
    if idx + 1 >= len(parts):
        return parts
    out = list(parts)
    out[idx + 1] = str(max(1, int(threads)))
    return out


def _extract_threads(cmd: Sequence[str]) -> Optional[int]:
    parts = [str(x) for x in cmd]
    try:
        idx = parts.index("--threads")
    except ValueError:
        return None
    if idx + 1 >= len(parts):
        return None
    try:
        return int(parts[idx + 1])
    except Exception:
        return None


def run_with_lock_fallback(
    *,
    cmd: Sequence[str],
    run_kwargs: Mapping[str, Any],
    logger: Any,
    component: str,
    source: str,
    stage: str,
    scorch_python: Optional[Path],
    scorch_script: Optional[Path],
) -> tuple[subprocess.CompletedProcess[str], bool]:
    kwargs = dict(run_kwargs)
    timeout_sec = kwargs.pop("timeout", None)
    heartbeat_fn = kwargs.pop("heartbeat_fn", None)
    heartbeat_sec = kwargs.pop("heartbeat_sec", None)
    used_retry = False
    current_cmd = [str(x) for x in cmd]
    current = run_timed_subprocess(
        current_cmd,
        run_kwargs=kwargs,
        timeout_sec=timeout_sec,
        heartbeat_fn=heartbeat_fn,
        heartbeat_sec=heartbeat_sec,
    )
    if current.returncode == 0:
        return current, False
    if not lock_fallback_enabled():
        pass
    elif is_mamba_lock_error(current.stderr or "", current.stdout or ""):
        retry_cmd = build_env_python_retry_cmd(
            current_cmd,
            scorch_python=scorch_python,
            scorch_script=scorch_script,
        )
        if retry_cmd is not None:
            logger.warning(
                "%s action=score lock_retry=true source=%s stage=%s reason=mamba_lock_detected returncode=%s",
                component,
                source,
                stage,
                current.returncode,
            )
            current = run_timed_subprocess(
                retry_cmd,
                run_kwargs=kwargs,
                timeout_sec=timeout_sec,
                heartbeat_fn=heartbeat_fn,
                heartbeat_sec=heartbeat_sec,
            )
            current_cmd = retry_cmd
            used_retry = True
            logger.info(
                "%s action=score lock_retry=true source=%s stage=%s mode=env_python retry_returncode=%s",
                component,
                source,
                stage,
                current.returncode,
            )
        else:
            logger.warning(
                "%s action=score lock_retry=false source=%s stage=%s reason=fallback_unavailable",
                component,
                source,
                stage,
            )
    if (
        current.returncode != 0
        and runtime_retry_enabled()
        and is_runtime_pool_error(current.stderr or "", current.stdout or "")
    ):
        prev_threads = _extract_threads(current_cmd)
        if prev_threads is None or int(prev_threads) > 1:
            runtime_cmd = _rewrite_threads(current_cmd, 1)
            logger.warning(
                "%s action=score runtime_retry=true source=%s stage=%s reason=runtime_pool_error threads_prev=%s threads_retry=1 returncode=%s",
                component,
                source,
                stage,
                str(prev_threads),
                current.returncode,
            )
            current = run_timed_subprocess(
                runtime_cmd,
                run_kwargs=kwargs,
                timeout_sec=timeout_sec,
                heartbeat_fn=heartbeat_fn,
                heartbeat_sec=heartbeat_sec,
            )
            used_retry = True
            logger.info(
                "%s action=score runtime_retry=true source=%s stage=%s retry_returncode=%s",
                component,
                source,
                stage,
                current.returncode,
            )
    return current, used_retry


def run_timed_subprocess(
    cmd: Sequence[str],
    *,
    run_kwargs: Mapping[str, Any],
    timeout_sec: Optional[float] = None,
    heartbeat_fn: Optional[Callable[[], None]] = None,
    heartbeat_sec: Optional[float] = None,
) -> subprocess.CompletedProcess[str]:
    """
    Run a subprocess with optional process-group timeout and heartbeat support.

    With no timeout/heartbeat this intentionally falls back to subprocess.run so
    existing tests and monkeypatches keep their old behavior.
    """
    kwargs = dict(run_kwargs)
    if timeout_sec is None and heartbeat_fn is None:
        return subprocess.run([str(x) for x in cmd], **kwargs)
    timeout_value: Optional[float]
    try:
        timeout_value = float(timeout_sec) if timeout_sec is not None else None
    except Exception:
        timeout_value = None
    if timeout_value is not None and timeout_value <= 0:
        timeout_value = None
    try:
        heartbeat_interval = float(heartbeat_sec) if heartbeat_sec is not None else 60.0
    except Exception:
        heartbeat_interval = 60.0
    heartbeat_interval = max(1.0, float(heartbeat_interval))

    text_mode = bool(kwargs.get("text") or kwargs.get("universal_newlines"))
    kwargs.pop("check", None)
    capture_output = bool(kwargs.pop("capture_output", False))
    if capture_output:
        kwargs.setdefault("stdout", subprocess.PIPE)
        kwargs.setdefault("stderr", subprocess.PIPE)
    kwargs.setdefault("start_new_session", True)

    stop_event = threading.Event()

    def _heartbeat_loop() -> None:
        if heartbeat_fn is None:
            return
        while not stop_event.wait(heartbeat_interval):
            try:
                heartbeat_fn()
            except Exception:
                pass

    hb_thread: Optional[threading.Thread] = None
    if heartbeat_fn is not None:
        hb_thread = threading.Thread(
            target=_heartbeat_loop,
            name="scorch-subprocess-heartbeat",
            daemon=True,
        )
        hb_thread.start()

    proc = subprocess.Popen([str(x) for x in cmd], **kwargs)
    timed_out = False
    try:
        stdout, stderr = proc.communicate(timeout=timeout_value)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except Exception:
            try:
                proc.terminate()
            except Exception:
                pass
        try:
            stdout, stderr = proc.communicate(timeout=10.0)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except Exception:
                try:
                    proc.kill()
                except Exception:
                    pass
            stdout, stderr = proc.communicate()
    finally:
        stop_event.set()
        if hb_thread is not None:
            hb_thread.join(timeout=1.0)

    rc = int(proc.returncode if proc.returncode is not None else 124)
    if timed_out:
        rc = 124
        timeout_msg = (
            f"SCORCH subprocess timeout after {float(timeout_value or 0.0):.1f}s"
        )
        if text_mode:
            stderr = (stderr or "") + ("\n" if stderr else "") + timeout_msg
        elif isinstance(stderr, (bytes, bytearray)):
            stderr = bytes(stderr) + (b"\n" if stderr else b"") + timeout_msg.encode()
    return subprocess.CompletedProcess([str(x) for x in cmd], rc, stdout, stderr)
