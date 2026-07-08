"""Runtime indirection for cli.qol_cli monkeypatch compatibility."""

from __future__ import annotations

from typing import Any


def repo_root():
    from cli import qol_cli

    return qol_cli._repo_root()


def load_effective_config() -> dict[str, Any]:
    from cli import qol_cli

    return qol_cli.load_effective_config()


def build_setup_report() -> dict[str, Any]:
    from cli import qol_cli

    return qol_cli._build_setup_report()


def subprocess():
    from cli import qol_cli

    return qol_cli.subprocess
