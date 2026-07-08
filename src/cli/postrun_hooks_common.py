# -*- coding: utf-8 -*-
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from config.value_access import to_bool


def _is_no_library_docking(cfg: Mapping[str, Any]) -> bool:
    return bool(to_bool(cfg.get("NO_LIBRARY_DOCKING"), default=False))


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)

