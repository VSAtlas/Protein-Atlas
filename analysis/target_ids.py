"""Helpers for report target identifiers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


def _clean(value: Any) -> str:
    return str(value or "").strip()


@dataclass(frozen=True)
class TargetId:
    pdb_id: str
    variant: str
    ph_label: str

    def format(self) -> str:
        return build_target_id(self.pdb_id, self.variant, self.ph_label)


def build_target_id(pdb_id: Any, variant: Any, ph_label: Any) -> str:
    return "|".join([_clean(pdb_id), _clean(variant), _clean(ph_label)])


def parse_target_id(target_id: Any) -> TargetId:
    parts = _clean(target_id).split("|")
    return TargetId(
        pdb_id=(parts[0].strip().upper() if len(parts) > 0 else ""),
        variant=(parts[1].strip() if len(parts) > 1 else ""),
        ph_label=(parts[2].strip() if len(parts) > 2 else ""),
    )


def pdb_id_from_target_id(target_id: Any) -> str:
    return parse_target_id(target_id).pdb_id
