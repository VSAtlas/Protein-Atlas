from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

TRUTHY_TEXT = {"1", "true", "yes", "y", "on"}
DECOY_ROLE_TOKENS = {"dud", "decoy", "decoys"}


def truthy_text(value: Any) -> bool:
    return str(value or "").strip().lower() in TRUTHY_TEXT


def decoy_role_token(value: Any) -> bool:
    return str(value or "").strip().lower() in DECOY_ROLE_TOKENS


def row_is_control(row: Mapping[str, Any], key: str = "is_control") -> bool:
    return truthy_text(row.get(key))


def row_has_explicit_decoy(row: Mapping[str, Any], key: str = "is_decoy") -> bool:
    return truthy_text(row.get(key))


def row_is_decoy_by_role(
    row: Mapping[str, Any], *, run_mode_key: str = "run_mode", library_key: str = "library"
) -> bool:
    return decoy_role_token(row.get(run_mode_key)) or decoy_role_token(
        row.get(library_key)
    )


def ligand_filename_from_row(row: Mapping[str, Any]) -> str:
    ligand_file = row.get("ligand_file") or row.get("ligand") or row.get("Ligand_ID") or ""
    try:
        return Path(str(ligand_file)).name
    except Exception:
        return str(ligand_file)
