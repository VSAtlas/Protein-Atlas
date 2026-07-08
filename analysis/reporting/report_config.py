"""Shared report config file readers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

from analysis.reporting.value_utils import parse_bool, strip_quotes
from config.output_paths import run_output_dir

_INLINE_COMMENT_RE = re.compile(r"\s+#.*$")


def report_config_candidates(repo_root: Path, run_id: str) -> list[Path]:
    return [
        run_output_dir(repo_root, "data", run_id) / "config.txt",
        run_output_dir(repo_root, "data", run_id) / "config_snapshot.txt",
        run_output_dir(repo_root, "manifests", run_id) / "config.txt",
        Path(repo_root) / "config.txt",
    ]


def read_config_key(
    path: Path,
    key_names: str | Iterable[str],
    *,
    strip_inline_comments: bool = False,
    keys_case_insensitive: bool = False,
) -> str | None:
    keys = {key_names} if isinstance(key_names, str) else set(key_names)
    keys_upper = {k.upper() for k in keys} if keys_case_insensitive else None
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if strip_inline_comments:
                    stripped = _INLINE_COMMENT_RE.sub("", stripped).strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                key, raw_value = stripped.split("=", 1)
                file_key = key.strip()
                if keys_case_insensitive:
                    keys_upper_set = keys_upper or {k.upper() for k in keys}
                    if file_key.upper() not in keys_upper_set:
                        continue
                elif file_key not in keys:
                    continue
                return strip_quotes(raw_value)
    except OSError:
        return None
    return None


def read_config_keys(
    path: Path,
    key_names: Iterable[str],
    *,
    strip_inline_comments: bool = False,
    keys_case_insensitive: bool = False,
) -> dict[str, str | None]:
    """Read many keys from one pass over a config file (same semantics as :func:`read_config_key`)."""
    keys_list = list(dict.fromkeys(key_names))
    result: dict[str, str | None] = {k: None for k in keys_list}
    if not path.exists():
        return result
    if keys_case_insensitive:
        upper_to_canon: dict[str, str] = {}
        for k in keys_list:
            upper_to_canon.setdefault(k.upper(), k)
    else:
        upper_to_canon = {}
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if strip_inline_comments:
                    stripped = _INLINE_COMMENT_RE.sub("", stripped).strip()
                if not stripped or stripped.startswith("#") or "=" not in stripped:
                    continue
                key, raw_value = stripped.split("=", 1)
                file_key = key.strip()
                if keys_case_insensitive:
                    canon = upper_to_canon.get(file_key.upper())
                    if canon is None:
                        continue
                elif file_key not in result:
                    continue
                else:
                    canon = file_key
                if result[canon] is None:
                    result[canon] = strip_quotes(raw_value)
    except OSError:
        return result
    return result


def resolve_config_key(
    repo_root: Path,
    run_id: str,
    key_names: str | Iterable[str],
) -> str | None:
    for path in report_config_candidates(repo_root, run_id):
        value = read_config_key(path, key_names)
        if value is not None:
            return value
    return None


def resolve_config_int(repo_root: Path, run_id: str, key: str) -> int | None:
    value = resolve_config_key(repo_root, run_id, key)
    if value is None:
        return None
    try:
        parsed = int(value)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def resolve_config_float(
    repo_root: Path,
    run_id: str,
    key: str,
    *,
    allow_percent: bool = False,
) -> float | None:
    value = resolve_config_key(repo_root, run_id, key)
    if value is None:
        return None
    text = value.strip()
    scale = 1.0
    if allow_percent and text.endswith("%"):
        text = text[:-1].strip()
        scale = 0.01
    try:
        return float(text) * scale
    except ValueError:
        return None


def resolve_config_bool(repo_root: Path, run_id: str, key: str) -> bool | None:
    return parse_bool(resolve_config_key(repo_root, run_id, key), default=None)
