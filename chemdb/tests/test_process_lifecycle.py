from __future__ import annotations

import logging

import pytest

from cli.process_lifecycle import process_lifecycle


def _messages(caplog: pytest.LogCaptureFixture) -> list[str]:
    return [record.getMessage() for record in caplog.records]


def test_process_lifecycle_logs_clean_return(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.INFO):
        with process_lifecycle(run_id="clean", is_resume=False):
            pass

    messages = _messages(caplog)
    assert any("action=process_start run_id=clean" in message for message in messages)
    assert any("action=process_exit run_id=clean" in message for message in messages)
    assert any("status=clean" in message for message in messages)
    assert not any("action=unhandled_exception" in message for message in messages)


def test_process_lifecycle_treats_system_exit_zero_as_clean(
    caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level(logging.INFO), pytest.raises(SystemExit) as exc_info:
        with process_lifecycle(run_id="plan", is_resume=False):
            raise SystemExit(0)

    assert exc_info.value.code == 0
    messages = _messages(caplog)
    assert any("status=clean reason=system_exit exit_code=0" in m for m in messages)
    assert not any("action=unhandled_exception" in message for message in messages)


def test_process_lifecycle_treats_nonzero_system_exit_as_failed(
    caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level(logging.INFO), pytest.raises(SystemExit) as exc_info:
        with process_lifecycle(run_id="failed", is_resume=True):
            raise SystemExit(1)

    assert exc_info.value.code == 1
    messages = _messages(caplog)
    assert any("status=failed reason=system_exit exit_code=1" in m for m in messages)
    assert not any("action=unhandled_exception" in message for message in messages)


def test_process_lifecycle_logs_keyboard_interrupt_separately(
    caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level(logging.INFO), pytest.raises(KeyboardInterrupt):
        with process_lifecycle(run_id="interrupted", is_resume=True):
            raise KeyboardInterrupt

    messages = _messages(caplog)
    assert any("status=interrupted reason=keyboard_interrupt" in m for m in messages)
    assert not any("action=unhandled_exception" in message for message in messages)


def test_process_lifecycle_logs_true_unhandled_exception(
    caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level(logging.INFO), pytest.raises(RuntimeError, match="boom"):
        with process_lifecycle(run_id="crash", is_resume=False):
            raise RuntimeError("boom")

    messages = _messages(caplog)
    assert any(
        "action=unhandled_exception run_id=crash" in message
        and "exception_type=RuntimeError" in message
        for message in messages
    )
    assert not any("action=process_exit run_id=crash" in message for message in messages)
