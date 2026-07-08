"""Shared config value access and scalar coercion helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, cast


def cfg_get(cfg: object | None, key: str, default: object = None) -> object:
    if cfg is None:
        return default
    if isinstance(cfg, Mapping):
        return cfg.get(key, default)
    getter = getattr(cfg, "get", None)
    if callable(getter):
        try:
            return getter(key, default)
        except (TypeError, AttributeError):
            pass
    return getattr(cfg, key, default)


def cfg_first(cfg: object | None, keys: Sequence[str], default: object = None) -> object:
    for key in keys:
        value = cfg_get(cfg, key, None)
        if value is None:
            continue
        if str(value).strip() == "":
            continue
        return value
    return default


def cfg_has_value(cfg: object | None, key: str) -> bool:
    value = cfg_get(cfg, key, None)
    return value is not None and str(value).strip() != ""


def to_bool(value: object, default: bool | None = False) -> bool | None:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    token = str(value).strip().lower()
    if token in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if token in {"0", "false", "f", "no", "n", "off"}:
        return False
    return default


def to_int(value: object, default: int | None = None) -> int | None:
    try:
        return int(cast(Any, value))
    except (TypeError, ValueError):
        return default


def to_float(value: object, default: float | None = None) -> float | None:
    try:
        return float(cast(Any, value))
    except (TypeError, ValueError):
        return default
