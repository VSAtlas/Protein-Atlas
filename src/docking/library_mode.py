from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Dict, Optional


_RESERVED_TOKENS = ("dud", "hmdb", "fda")
_OFF_VALUES = {"", "0", "false", "no", "off", "none", "null"}
_ON_VALUES = {"true", "yes", "on", "1"}


def is_truthy(cfg: dict, key: str, default: bool = True) -> bool:
    if key not in cfg:
        return default
    value = cfg.get(key)
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    lowered = str(value).strip().lower()
    if lowered in _ON_VALUES:
        return True
    if lowered in _OFF_VALUES:
        return False
    return bool(value)


def _parse_test_libraries_value(raw: object) -> list[str]:
    if isinstance(raw, bool):
        return ["dud", "fda"] if raw else ["fda"]

    s = str(raw).strip()
    if not s:
        return ["fda"]
    lowered = s.lower()

    if lowered in _OFF_VALUES:
        return ["fda"]
    if lowered in _ON_VALUES:
        return ["dud", "fda"]
    if lowered == "default":
        return ["fda"]
    if lowered in {"both", "fda_dud", "dud_fda", "fda+dud", "dud+fda"}:
        return ["dud", "fda"]
    if lowered in {"fda-dud", "dud-fda"}:
        return ["dud", "fda"]

    tokens = [tok for tok in re.split(r"[+,\s]+", lowered) if tok]
    normalized: list[str] = []
    for tok in tokens:
        if tok in {"and", "off", "none", "null"}:
            continue
        if tok == "default":
            tok = "fda"
        normalized.append(tok)

    if not normalized:
        return ["fda"]

    deduped: list[str] = []
    seen: set[str] = set()
    for tok in normalized:
        if tok not in seen:
            seen.add(tok)
            deduped.append(tok)

    return deduped


def parse_test_libraries(cfg: dict) -> list[str]:
    raw = (
        os.environ.get("TEST_MODE_ENABLE")
        if "TEST_MODE_ENABLE" in os.environ
        else cfg.get("TEST_MODE_ENABLE", "off")
    )
    return _parse_test_libraries_value(raw)


def _coerce_test_map(m) -> Dict[str, str]:
    import ast as _ast
    import json as _json

    if isinstance(m, dict):
        return {str(k).upper(): str(v) for k, v in m.items()}
    s = str(m).strip()
    if not s:
        return {}
    parsed = None
    try:
        parsed = _json.loads(s)
    except Exception:
        try:
            parsed = _ast.literal_eval(s)
        except Exception:
            parsed = {}
    return {
        str(k).upper(): str(v)
        for k, v in (parsed if isinstance(parsed, dict) else {}).items()
    }


def _is_pytest_context(cfg: dict) -> bool:
    try:
        selection_mode = str(cfg.get("PDB_SELECTION_MODE", "")).strip().lower()
    except Exception:
        selection_mode = ""
    env_selection = str(os.environ.get("PDB_SELECTION_MODE", "")).strip().lower()
    return (
        "pytest" in os.environ.get("PYTEST_CURRENT_TEST", "")
        or selection_mode == "test_library_map"
        or env_selection == "test_library_map"
    )


def _dedup_paths(seq: list[Path]) -> list[Path]:
    deduped: list[Path] = []
    seen: set[str] = set()
    for candidate in seq:
        if not candidate:
            continue
        path_obj = Path(candidate)
        try:
            key = str(path_obj.resolve())
        except Exception:
            key = str(path_obj)
        if key not in seen:
            seen.add(key)
            deduped.append(path_obj)
    return deduped


def compute_allowed_library_roots(
    cfg: dict,
    pdb_id: str,
    logger: logging.Logger,
    *,
    tokens_override: Optional[list[str]] = None,
) -> list[Path]:
    tokens = tokens_override or parse_test_libraries(cfg)

    base_root = Path(
        cfg.get("OUTPUT_LIGANDS_DIR")
        or cfg.get("PREPPED_LIGANDS_ROOT")
        or "prepped_ligands"
    )
    subdir_default = str(cfg.get("LIBRARY_SUBDIR_DEFAULT", "fda_library"))
    hmdb_subdir_raw = str(cfg.get("HMDB_LIBRARY_SUBDIR", "hmdb"))
    hmdb_test_subdir = str(cfg.get("HMDB_TEST_LIBRARY_SUBDIR", "hmdb_test_library_10"))

    pytest_mode = _is_pytest_context(cfg)
    hmdb_use_test = pytest_mode and ("hmdb" in tokens)
    hmdb_subdir = hmdb_test_subdir if hmdb_use_test else hmdb_subdir_raw
    if hmdb_use_test and hmdb_subdir != hmdb_subdir_raw:
        logger.info(
            "[ligands.hmdb-test] pytest=%s subdir=%s raw=%s",
            pytest_mode,
            hmdb_subdir,
            hmdb_subdir_raw,
        )

    maybe_map = cfg.get("TEST_LIBRARY_MAP", {})
    test_map = _coerce_test_map(maybe_map)
    logger.info(
        "[lib-roots.map] raw_type=%s keys=%d",
        type(maybe_map).__name__,
        len(test_map),
    )
    mapped_value = test_map.get(str(pdb_id).upper())

    if "dud" in tokens and not mapped_value:
        logger.warning(
            "[test-mode] PDB %s missing from TEST_LIBRARY_MAP; using default library=%s",
            pdb_id,
            subdir_default,
        )

    roots: list[Path] = []
    reserved_only = all(token in _RESERVED_TOKENS for token in tokens)
    if reserved_only:
        has_dud = "dud" in tokens
        has_hmdb = "hmdb" in tokens
        has_fda = "fda" in tokens
        if has_dud and has_hmdb and has_fda:
            tokens = ["dud", "hmdb", "fda"]
        elif has_hmdb and has_dud and not has_fda:
            tokens = ["hmdb", "dud"]
        elif has_hmdb and has_fda and not has_dud:
            tokens = ["hmdb", "fda"]
        elif has_dud and has_fda and not has_hmdb:
            tokens = ["dud", "fda"]
        elif has_dud:
            tokens = ["dud"]
        elif has_hmdb:
            tokens = ["hmdb"]
        elif has_fda:
            tokens = ["fda"]

    cfg["_TEST_MODE_EFFECTIVE_TOKENS"] = list(tokens)

    for token in tokens:
        if token == "fda":
            roots.append(base_root / subdir_default)
        elif token == "hmdb":
            roots.append(base_root / hmdb_subdir)
        elif token == "dud":
            if mapped_value:
                roots.append(base_root / mapped_value)
            else:
                roots.append(base_root / subdir_default)
        else:
            roots.append(base_root / token)

    roots = _dedup_paths(roots)

    allowed_noncontrol_roots: list[Path] = []
    for root in roots:
        if root.exists():
            allowed_noncontrol_roots.append(root)
        else:
            logger.warning("[ligands.test-roots] missing=%s", root)

    extra_paths_cfg: list[str] = cfg.get("EXTRA_LIGAND_ROOTS", []) or []
    for d in extra_paths_cfg:
        p = Path(d)
        if p.exists():
            allowed_noncontrol_roots.append(p)
        else:
            logger.warning("[ligands.extra-roots] missing=%s", p)

    extra_dirs = str(cfg.get("LIBRARY_EXTRA_DIRS", "")).strip()
    if extra_dirs:
        for d in extra_dirs.split(";"):
            d = d.strip()
            if not d:
                continue
            p = Path(d)
            if p.exists():
                allowed_noncontrol_roots.append(p)

    allowed_noncontrol_roots = _dedup_paths(allowed_noncontrol_roots)

    logger.info(
        "[ligands.allowed-roots] tokens=%s pdb=%s roots=%d",
        "+".join(tokens),
        pdb_id,
        len(allowed_noncontrol_roots),
    )
    cfg["_ALLOWED_NONCONTROL_ROOTS"] = [str(p) for p in allowed_noncontrol_roots]
    cfg["_TEST_MODE_EFFECTIVE"] = "+".join(tokens)

    primary_root_str = (
        str(allowed_noncontrol_roots[0]) if allowed_noncontrol_roots else None
    )
    logger.info(
        "[ph_ligand.roots] tokens=%s pdb=%s ph_root=%s noncontrol_roots=%s",
        "+".join(tokens),
        pdb_id,
        primary_root_str,
        cfg.get("_ALLOWED_NONCONTROL_ROOTS"),
    )
    return allowed_noncontrol_roots


__all__ = [
    "compute_allowed_library_roots",
    "is_truthy",
    "parse_test_libraries",
]
