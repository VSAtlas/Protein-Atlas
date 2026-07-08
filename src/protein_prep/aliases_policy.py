"""Shared alias/policy constants used across protein-prep modules."""

from __future__ import annotations

from typing import Iterable

from protein_prep.pdb_fixer_runtime import get_atom_rules


def load_aliases():
    return get_atom_rules()


def _flatten_semicolons(items):
    out = []
    for item in items or []:
        parts = [p.strip().upper() for p in str(item).split(";") if p.strip()]
        out.extend(parts)
    return out


def _canonize_two_letter(xs) -> set[str]:
    out = set()
    for line in xs or []:
        for tok in str(line).split(";"):
            t = tok.strip()
            if t:
                out.add(t.upper())
    return out


def _to_upper_set(values: Iterable[str] | None) -> set[str]:
    out: set[str] = set()
    for val in values or []:
        if val is None:
            continue
        text = str(val).strip()
        if text:
            out.add(text.upper())
    return out


ALIASES = get_atom_rules()
RULES = ALIASES.__dict__ if hasattr(ALIASES, "__dict__") else dict(ALIASES)
_ES = RULES.get("element_sets", {}) or {}
_ONE = {str(s).upper() for s in (_ES.get("one_letter_elements") or [])}
_TWORAW = _ES.get("two_letter_elements", []) or []
_TWO = _canonize_two_letter(_TWORAW)
_ALIAS_SETS = getattr(ALIASES, "alias_sets", None)
_ELEMENT_ALIAS_MAP = {
    str(k).strip().upper(): str(v).strip().upper()
    for k, v in (getattr(_ALIAS_SETS, "element_alias", {}) or {}).items()
    if str(k).strip() and str(v).strip()
}
_NORMALIZE_RESNAME_FN = getattr(ALIASES, "normalize_resname", None)


def _normalize_resname(token: str) -> str:
    base = (token or "").strip().upper()
    if not base:
        return ""
    if callable(_NORMALIZE_RESNAME_FN):
        try:
            normalized = _NORMALIZE_RESNAME_FN(token)
            if normalized:
                text = str(normalized).strip().upper()
                if text:
                    return text
        except Exception:
            pass
    return _ELEMENT_ALIAS_MAP.get(base, base)


_POLICY_MODE = getattr(ALIASES, "policy_mode", "LEGACY")
_RETAIN_VARIANT = _to_upper_set(getattr(ALIASES, "retain_resnames", []))
_RETAIN_VARIANT_CANONICAL = {
    _normalize_resname(tok)
    for tok in getattr(ALIASES, "retain_resnames", [])
    if _normalize_resname(tok)
}

_WATER_NAMES = _to_upper_set(
    getattr(ALIASES, "waters", getattr(_ALIAS_SETS, "waters", set()))
)

_COFACTOR_NAMES = _to_upper_set(getattr(ALIASES, "cofactors", set()))
_COFACTOR_CANONICAL = {
    _normalize_resname(tok)
    for tok in getattr(ALIASES, "cofactors", set())
    if _normalize_resname(tok)
}
_COFACTOR_RAW_ALL = _to_upper_set(
    getattr(ALIASES, "cofactors_all", getattr(_ALIAS_SETS, "cofactors", set()))
)

_ELEMENT_TOKENS_RAW = _to_upper_set(
    getattr(ALIASES, "element_tokens", getattr(_ALIAS_SETS, "element_tokens", set()))
)
_ELEM_CANON = _to_upper_set(
    getattr(
        ALIASES,
        "elem_tokens_canonical",
        getattr(_ALIAS_SETS, "elem_tokens_canonical", set()),
    )
)

_SALT_RESNAMES = {"NA", "K", "CL", "BR", "I"}
_METAL_RESNAMES = {"MG", "MN", "FE", "ZN", "CU", "CO", "NI", "CA"}

_ION_AUDIT_METALS = {
    "ZN",
    "MG",
    "MN",
    "FE",
    "CO",
    "NI",
    "CU",
    "CD",
    "HG",
    "CA",
}
_ION_AUDIT_SIMPLE_IONS = {"NA", "K", "CL", "BR", "I"}
_ION_AUDIT_WATERS = {"HOH", "WAT"}
_ION_AUDIT_ALIAS_MAP = {
    "ZN1": "ZN",
    "ZN2": "ZN",
    "ZN3": "ZN",
    "ZN+": "ZN",
    "ZN+2": "ZN",
    "MG1": "MG",
    "MG2": "MG",
    "MG+": "MG",
    "MN2": "MN",
    "MN3": "MN",
    "FE2": "FE",
    "FE3": "FE",
    "CO2": "CO",
    "NI2": "NI",
    "CU1": "CU",
    "CU2": "CU",
    "CD2": "CD",
    "HG2": "HG",
    "CA1": "CA",
    "CA2": "CA",
    "NA1": "NA",
    "K1": "K",
    "CL-": "CL",
    "BR-": "BR",
    "I-": "I",
}
_ION_AUDIT_ENABLED_VALUES = {"1", "true", "yes"}
_ION_AUDIT_DISABLED_VALUES = {"0", "false", "no"}

_ION_BREADCRUMB_METAL_ORDER = (
    "ZN",
    "HG",
    "MG",
    "FE",
    "MN",
    "CO",
    "NI",
    "CU",
    "CD",
    "CA",
)
_ION_BREADCRUMB_SIMPLE_ORDER = ("NA", "K", "CL", "BR", "I")

# Logging-state globals intentionally kept in one place for legacy hydration users.
_ALIASES_BIND_LOGGED = False
_ION_KEEP_LOGGED: set[str] = set()
_COFACTOR_DROP_LOGGED: set[str] = set()
_IONS_CFG_CACHE: tuple[set[str], str] | None = None
_IONS_CFG_LOGGED = False
