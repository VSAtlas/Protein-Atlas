"""Shared trust predicates for versioned FDA mapping identity rows."""

from __future__ import annotations

from typing import Mapping


def clean_mapping_value(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "null", "<na>"} else text


def is_resolved_identity_status(value: object) -> bool:
    status = clean_mapping_value(value).casefold()
    return status.startswith("manifest_resolved_") or status.startswith("terminal_")


def resolved_preferred_name(row: Mapping[str, object]) -> str:
    if not is_resolved_identity_status(row.get("identity_resolution_status")):
        return ""
    return clean_mapping_value(row.get("resolved_preferred_name"))


def authoritative_drugcentral_id(row: Mapping[str, object]) -> str:
    """Return the identity-status-authorized DrugCentral ID, if one exists."""

    if is_resolved_identity_status(row.get("identity_resolution_status")):
        return clean_mapping_value(row.get("resolved_drugcentral_id"))
    return clean_mapping_value(row.get("drugcentral_id"))


def has_resolved_identity(row: Mapping[str, object]) -> bool:
    return bool(resolved_preferred_name(row))


def suppresses_legacy_aliases(row: Mapping[str, object]) -> bool:
    """Resolved and terminal identities admit only their preferred name and exact ID."""

    return has_resolved_identity(row)


__all__ = [
    "authoritative_drugcentral_id",
    "clean_mapping_value",
    "has_resolved_identity",
    "is_resolved_identity_status",
    "resolved_preferred_name",
    "suppresses_legacy_aliases",
]
