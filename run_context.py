# -*- coding: utf-8 -*-
from __future__ import annotations

import datetime
import json
import logging
import os
from typing import Any, Dict, Mapping, Optional

from run_manifest import load_run_manifest
from cli_utils import _cli_val


def _resolve_run_id(argv: list[str]) -> str:
    cli_run_id = _cli_val(argv, "--run-id") or _cli_val(argv, "-run-id")
    env_run_id = (os.environ.get("ATLAS_RUN_ID") or "").strip()
    if cli_run_id:
        return cli_run_id
    if env_run_id and env_run_id.lower() not in {"smoke_demo", "main_smoke_demo"}:
        return env_run_id
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def _prepare_run_logfile(run_id: str) -> str:
    logs_dir = os.path.join(os.getcwd(), "logs")
    os.makedirs(logs_dir, exist_ok=True)
    return os.path.join(logs_dir, f"main_{run_id}.log")


class ConfigDict(dict):
    __slots__ = ()
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as exc:
            raise AttributeError(key) from exc
    def __setattr__(self, key, value):
        self[key] = value
    def copy(self):
        return ConfigDict(super().copy())


def _apply_resume_config_from_snapshot(
    cfg: Mapping[str, Any],
    run_id: str,
) -> dict:
    """
    For resume mode:

    - Load run_manifest.yaml for run_id.
    - Find the per-run config snapshot (paths.run_dir + paths.config_file).
    - Load it as a dict.
    - Return that dict so the caller can use it as the canonical cfg.

    If anything goes wrong, this logs a WARNING and returns the original cfg as a plain dict.
    """
    logger = logging.getLogger("resume.config")
    try:
        manifest = load_run_manifest(cfg, run_id)
        if not manifest:
            logger.warning(
                "[resume.config] manifest missing or empty for run_id=%s; "
                "falling back to current cfg.",
                run_id,
            )
            return dict(cfg)

        paths = manifest.get("paths") or {}
        run_dir = paths.get("run_dir")
        cfg_file_name = paths.get("config_file", "run_config.yaml")

        if not run_dir:
            logger.warning(
                "[resume.config] manifest missing run_dir for run_id=%s; "
                "falling back to current cfg.",
                run_id,
            )
            return dict(cfg)

        cfg_path = os.path.join(run_dir, cfg_file_name)
        if not os.path.exists(cfg_path):
            logger.warning(
                "[resume.config] config snapshot missing for run_id=%s path=%s; "
                "falling back to current cfg.",
                run_id,
                cfg_path,
            )
            return dict(cfg)

        try:
            with open(cfg_path, "r", encoding="utf-8") as handle:
                snapshot = json.load(handle)
            if not isinstance(snapshot, dict):
                raise ValueError("snapshot is not a dict")
            logger.info(
                "[resume.config] using snapshot for run_id=%s path=%s", run_id, cfg_path
            )
            return snapshot
        except Exception as exc:
            logger.warning(
                "[resume.config] failed to load snapshot for run_id=%s path=%s err=%s; "
                "falling back to current cfg.",
                run_id,
                cfg_path,
                exc,
            )
            return dict(cfg)
    except Exception:
        logger.warning(
            "[resume.config] unexpected error for run_id=%s; falling back to current cfg.",
            run_id,
            exc_info=True,
        )
        return dict(cfg)
