"""Canonical residue loader helpers for active-site processing."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
import re
from typing import Any, Set

from . import ALIASES_PATH, yaml as _yaml

_DEFAULT_CFG_CACHE: Mapping[str, Any] | None = None

_TOKEN_SPLIT = re.compile(r"[;\s,]+")


def _iter_values(payload: Any) -> Iterable[Any]:
    if payload is None:
        return []
    if isinstance(payload, Mapping):
        items: list[Any] = []
        for key, value in payload.items():
            items.append(key)
            items.extend(_iter_values(value))
        return items
    if isinstance(payload, (str, bytes)):
        return [payload]
    if isinstance(payload, Iterable):
        items: list[Any] = []
        for entry in payload:
            items.extend(_iter_values(entry))
        return items
    return [payload]


def _load_default_cfg() -> Mapping[str, Any]:
    global _DEFAULT_CFG_CACHE
    if _DEFAULT_CFG_CACHE is not None:
        return _DEFAULT_CFG_CACHE
    try:
        with open(ALIASES_PATH, "r", encoding="utf-8") as handle:
            data = _yaml.safe_load(handle) or {}
    except Exception:
        data = {}
    if not isinstance(data, Mapping):
        data = {}
    _DEFAULT_CFG_CACHE = data
    return data


def _normalize_tokens(raw_tokens: Iterable[Any]) -> Set[str]:
    tokens: Set[str] = set()
    for raw in raw_tokens:
        if raw is None:
            continue
        text = str(raw).strip()
        if not text:
            continue
        for piece in _TOKEN_SPLIT.split(text):
            token = piece.strip()
            if token:
                tokens.add(token.upper())
    return tokens


def _collect_alias_synonyms(
    alias_map: Any,
    canonical: Set[str],
) -> Set[str]:
    if not isinstance(alias_map, Mapping):
        return set()
    synonyms: Set[str] = set()
    for alias, target in alias_map.items():
        alias_tokens = _normalize_tokens([alias])
        target_tokens = _normalize_tokens([target])
        if canonical & target_tokens:
            synonyms |= alias_tokens
    return synonyms


def _resolve_cfg(cfg: Mapping[str, Any] | None) -> Mapping[str, Any]:
    if cfg:
        return cfg
    return _load_default_cfg()


def load_canonical_metals(cfg: Mapping[str, Any] | None) -> Set[str]:
    """Load canonical metal residue names and their aliases from configuration."""

    cfg = _resolve_cfg(cfg)
    canonical = _normalize_tokens(_iter_values(cfg.get("canonical_metals")))

    element_sets = cfg.get("element_sets")
    if isinstance(element_sets, Mapping):
        canonical |= _collect_alias_synonyms(
            element_sets.get("cation_resname_aliases"), canonical
        )
        canonical |= _collect_alias_synonyms(
            element_sets.get("halide_resname_aliases"), canonical
        )

    canonical |= _collect_alias_synonyms(cfg.get("retain_element_alias_map"), canonical)
    return {token.strip().upper() for token in canonical if token.strip()}


def load_canonical_cofactors(cfg: Mapping[str, Any] | None) -> Set[str]:
    """Load canonical cofactor residue names from configuration."""

    cfg = _resolve_cfg(cfg)
    return _normalize_tokens(_iter_values(cfg.get("canonical_cofactors")))


def load_canonical_waters(cfg: Mapping[str, Any] | None) -> Set[str]:
    """Load canonical water residue names from configuration."""

    cfg = _resolve_cfg(cfg)
    return _normalize_tokens(_iter_values(cfg.get("canonical_waters")))


__all__ = [
    "load_canonical_metals",
    "load_canonical_cofactors",
    "load_canonical_waters",
]
