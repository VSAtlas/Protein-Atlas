from __future__ import annotations

import logging
import os
from collections.abc import Iterator
from contextlib import contextmanager


def _system_exit_code(code: object) -> int:
    if code is None:
        return 0
    if isinstance(code, int):
        return code
    return 1


@contextmanager
def process_lifecycle(*, run_id: str, is_resume: bool) -> Iterator[None]:
    """Emit one start record and one truthful terminal record for a pipeline run."""
    logging.info(
        "[run.lifecycle] action=process_start run_id=%s pid=%d ppid=%d resume=%s",
        run_id,
        os.getpid(),
        os.getppid(),
        is_resume,
    )
    try:
        yield
    except SystemExit as exc:
        exit_code = _system_exit_code(exc.code)
        status = "clean" if exit_code == 0 else "failed"
        log = logging.info if exit_code == 0 else logging.error
        log(
            "[run.lifecycle] action=process_exit run_id=%s pid=%d "
            "status=%s reason=system_exit exit_code=%d",
            run_id,
            os.getpid(),
            status,
            exit_code,
        )
        raise
    except KeyboardInterrupt:
        logging.warning(
            "[run.lifecycle] action=process_exit run_id=%s pid=%d "
            "status=interrupted reason=keyboard_interrupt exit_code=130",
            run_id,
            os.getpid(),
        )
        raise
    except Exception as exc:
        logging.exception(
            "[run.lifecycle] action=unhandled_exception run_id=%s pid=%d "
            "resume=%s exception_type=%s",
            run_id,
            os.getpid(),
            is_resume,
            type(exc).__name__,
        )
        raise
    else:
        logging.info(
            "[run.lifecycle] action=process_exit run_id=%s pid=%d status=clean",
            run_id,
            os.getpid(),
        )
