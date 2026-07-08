from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Callable


def _env_truthy(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def send_run_email(status: int, start_time: str, end_time: str) -> None:
    """Best-effort notification for direct `python main.py` runs."""
    try:
        try:
            host = os.uname().nodename
        except AttributeError:
            host = "unknown-host"

        subject = f"Atlas run exited with status {status} on {host}"
        invoked = Path(sys.argv[0]).name or "atlas"
        cmd_line = " ".join([invoked, *sys.argv[1:]])
        body = "\n".join(
            [
                f"Atlas run finished on host: {host}",
                f"Start time: {start_time}",
                f"End time:   {end_time}",
                f"Exit status: {status}",
                "",
                "Command:",
                cmd_line,
                "",
            ]
        )

        mail_bin = shutil.which("mail")
        if not mail_bin:
            logging.info("[notify] mail command unavailable; skipping email")
            return

        try:
            proc = subprocess.Popen(
                [mail_bin, "-s", subject, "mpg2352@utexas.edu"],
                stdin=subprocess.PIPE,
                text=True,
            )
        except Exception as exc:
            logging.warning("[notify] failed to spawn mail command err=%s", exc)
            return

        try:
            proc.communicate(body, timeout=30)
        except Exception as exc:
            logging.warning("[notify] mail command failed err=%s", exc)
            try:
                proc.kill()
            except Exception:
                pass
    except Exception:
        try:
            logging.exception("[notify] unexpected error while building notification email")
        except Exception:
            pass


def run_with_email_notification(
    *,
    bootstrap_runtime: Callable[[], None],
    smoke_emit_config_demo: Callable[[], None],
    main_func: Callable[[], None],
) -> None:
    """
    Compatibility wrapper used only by direct `python main.py` execution.

    Keep production CLI routing in `src/cli/atlas_main_cli.py` and pipeline
    behavior in owner modules; this wrapper only preserves script-exit behavior.
    """
    start_time = time.strftime("%Y-%m-%d %H:%M:%S")
    status = 0

    try:
        bootstrap_runtime()
        if _env_truthy(os.environ.get("ATLAS_SMOKE_EMIT_CONFIG_DEMO")):
            smoke_emit_config_demo()
        main_func()
    except SystemExit as exc:
        code = exc.code
        status = code if isinstance(code, int) else 1
        raise
    except BaseException:
        status = 1
        raise
    finally:
        end_time = time.strftime("%Y-%m-%d %H:%M:%S")
        try:
            send_run_email(status=status, start_time=start_time, end_time=end_time)
        except Exception:
            try:
                logging.exception("[notify] email wrapper raised unexpectedly")
            except Exception:
                pass
