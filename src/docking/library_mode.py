from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Any, Dict, Mapping, Optional


_RESERVED_TOKENS = ("dud", "hmdb", "fda")
_OFF_VALUES = {"", "0", "false", "no", "off", "none", "null"}
_ON_VALUES = {"true", "yes", "on", "1"}
_REPO_ROOT = Path(__file__).resolve().parents[2]


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


def test_libraries_value_with_token(raw: object, required_token: str) -> str:
    tokens = _parse_test_libraries_value(raw)
    token = str(required_token).strip().lower()
    if token and token not in tokens:
        tokens.insert(0, token)
    return "+".join(tokens)


def parse_test_libraries(cfg: Mapping[str, Any]) -> list[str]:
    raw = cfg.get("_TEST_MODE_ENABLE_OVERRIDE")
    if raw is None:
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


def _coerce_test_map_cached(cfg: dict, m: object) -> Dict[str, str]:
    cache_key: tuple[str, str]
    if isinstance(m, dict):
        # Dict parse cost is low; avoid stale cache risk from in-place dict mutation.
        return _coerce_test_map(m)
    cache_key = (type(m).__name__, str(m))
    cached_key = cfg.get("_TEST_LIBRARY_MAP_CACHE_KEY")
    cached_val = cfg.get("_TEST_LIBRARY_MAP_CACHE_VAL")
    if cached_key == cache_key and isinstance(cached_val, dict):
        return {str(k): str(v) for k, v in cached_val.items()}
    parsed = _coerce_test_map(m)
    cfg["_TEST_LIBRARY_MAP_CACHE_KEY"] = cache_key
    cfg["_TEST_LIBRARY_MAP_CACHE_VAL"] = dict(parsed)
    return parsed


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


def _fixture_prepped_root(cfg: dict) -> Path:
    raw = cfg.get("TEST_FIXTURE_PREPPED_LIGANDS_DIR")
    if raw:
        return Path(raw)
    return _REPO_ROOT / "chemdb" / "tests" / "fixtures" / "prepped_ligands"


def _allow_fixture_fallback(cfg: dict, pdb_id: str, pytest_mode: bool) -> bool:
    if pytest_mode:
        return True
    if str(pdb_id).upper() == "TEST":
        return True
    return is_truthy(cfg, "USE_TEST_FIXTURES", default=False)


def _log_allowed_roots(
    cfg: dict,
    logger: logging.Logger,
    *,
    tokens: list[str],
    pdb_id: str,
    roots: list[Path],
) -> None:
    log_key = "|".join(
        [
            "+".join(tokens),
            str(pdb_id).upper(),
            ";".join(str(p) for p in roots),
        ]
    )
    seen_raw = cfg.get("_ALLOWED_ROOTS_LOG_SEEN")
    if isinstance(seen_raw, set):
        seen = seen_raw
    else:
        seen = set()
        cfg["_ALLOWED_ROOTS_LOG_SEEN"] = seen
    first_emit = log_key not in seen
    if first_emit:
        seen.add(log_key)
    log_fn = logger.info if first_emit else logger.debug
    log_fn(
        "[ligands.allowed-roots] tokens=%s pdb=%s roots=%d",
        "+".join(tokens),
        pdb_id,
        len(roots),
    )
    primary_root_str = str(roots[0]) if roots else None
    log_fn(
        "[ph_ligand.roots] tokens=%s pdb=%s ph_root=%s noncontrol_roots=%s",
        "+".join(tokens),
        pdb_id,
        primary_root_str,
        cfg.get("_ALLOWED_NONCONTROL_ROOTS"),
    )


def compute_allowed_library_roots(
    cfg: dict,
    pdb_id: str,
    logger: logging.Logger,
    *,
    tokens_override: Optional[list[str]] = None,
) -> list[Path]:
    tokens = tokens_override or parse_test_libraries(cfg)

    base_root = Path(
        cfg.get("PREPPED_LIGANDS_DIR")
        or cfg.get("OUTPUT_LIGANDS_DIR")
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
    test_map = _coerce_test_map_cached(cfg, maybe_map)
    map_fp = f"{type(maybe_map).__name__}:{len(test_map)}"
    if str(cfg.get("_LIB_ROOTS_MAP_LOG_FINGERPRINT", "")) != map_fp:
        logger.info(
            "[lib-roots.map] raw_type=%s keys=%d",
            type(maybe_map).__name__,
            len(test_map),
        )
        cfg["_LIB_ROOTS_MAP_LOG_FINGERPRINT"] = map_fp
    mapped_value = test_map.get(str(pdb_id).upper())

    if "dud" in tokens and not mapped_value:
        logger.warning(
            "[test-mode] PDB %s missing from TEST_LIBRARY_MAP; using default library=%s",
            pdb_id,
            subdir_default,
        )

    roots_cache_raw = cfg.get("_ALLOWED_LIBRARY_ROOTS_CACHE")
    roots_cache: Dict[str, list[str]]
    if isinstance(roots_cache_raw, dict):
        roots_cache = roots_cache_raw
    else:
        roots_cache = {}
        cfg["_ALLOWED_LIBRARY_ROOTS_CACHE"] = roots_cache
    cache_key = "|".join(
        [
            str(pdb_id).upper(),
            ",".join(tokens),
            str(base_root),
            str(_fixture_prepped_root(cfg)),
            str(subdir_default),
            str(hmdb_subdir),
            str(mapped_value or ""),
            str(cfg.get("LIBRARY_EXTRA_DIRS", "") or ""),
            ";".join(str(x) for x in (cfg.get("EXTRA_LIGAND_ROOTS", []) or [])),
        ]
    )
    cached_paths = roots_cache.get(cache_key)
    if isinstance(cached_paths, list):
        cached_roots = [Path(p) for p in cached_paths]
        cfg["_ALLOWED_NONCONTROL_ROOTS"] = [str(p) for p in cached_roots]
        cfg["_TEST_MODE_EFFECTIVE"] = "+".join(tokens)
        _log_allowed_roots(
            cfg,
            logger,
            tokens=tokens,
            pdb_id=pdb_id,
            roots=cached_roots,
        )
        return cached_roots

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
    fixture_root = _fixture_prepped_root(cfg)
    allow_fixture_fallback = _allow_fixture_fallback(cfg, pdb_id, pytest_mode)
    for root in roots:
        if root.exists():
            allowed_noncontrol_roots.append(root)
        elif allow_fixture_fallback and (fixture_root / root.name).exists():
            fixture_path = fixture_root / root.name
            allowed_noncontrol_roots.append(fixture_path)
            logger.info(
                "[ligands.test-fixture-root] missing=%s using=%s",
                root,
                fixture_path,
            )
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

    cfg["_ALLOWED_NONCONTROL_ROOTS"] = [str(p) for p in allowed_noncontrol_roots]
    cfg["_TEST_MODE_EFFECTIVE"] = "+".join(tokens)
    _log_allowed_roots(
        cfg,
        logger,
        tokens=tokens,
        pdb_id=pdb_id,
        roots=allowed_noncontrol_roots,
    )
    roots_cache[cache_key] = [str(p) for p in allowed_noncontrol_roots]
    return allowed_noncontrol_roots


__all__ = [
    "compute_allowed_library_roots",
    "is_truthy",
    "parse_test_libraries",
]
