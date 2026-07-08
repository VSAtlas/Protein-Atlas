from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Mapping, Union

ReviewRecord = dict[str, Any]
ValidationResult = dict[str, list[str] | bool]

REQUIRED_REVIEW_FIELDS = (
    "ligand_protonation",
    "tautomer",
    "stereochemistry",
    "net_charge",
    "parameter_review",
    "input_chemistry_authority",
    "water_policy",
    "metal_policy",
    "apo_holo",
    "ph",
)

APPROVED_STATUSES = frozenset({"reviewed", "approved", "ok"})
APO_HOLO_VALUES = frozenset({"APO", "HOLO"})


def review_record_path(root: Union[str, Path], target_id: str, ligand_id: str) -> Path:
    """Return the canonical per-target/per-ligand MMGBSA review JSON path."""
    target = _safe_path_token(target_id, "target_id")
    ligand = _safe_path_token(ligand_id, "ligand_id")
    return Path(root) / target / f"{ligand}.mmgbsa_review.json"


def default_review_record(target_id: str, ligand_id: str) -> ReviewRecord:
    """Create a minimal JSON-serializable MMGBSA chemistry review record."""
    return {
        "schema": "mmgbsa_chemistry_review_v1",
        "target_id": str(target_id),
        "ligand_id": str(ligand_id),
        "ligand_protonation": _blank_review_section(),
        "tautomer": _blank_review_section(),
        "stereochemistry": _blank_review_section(),
        "net_charge": {
            "reviewed": False,
            "status": None,
            "value": None,
            "source": None,
            "notes": "",
        },
        "parameter_review": _blank_review_section(),
        "input_chemistry_authority": _blank_review_section(),
        "water_policy": {
            "reviewed": False,
            "status": None,
            "policy": None,
            "retained_waters": [],
            "source": None,
            "notes": "",
        },
        "metal_policy": {
            "reviewed": False,
            "status": None,
            "policy": None,
            "retained_metals": [],
            "source": None,
            "notes": "",
        },
        "apo_holo": {
            "reviewed": False,
            "status": None,
            "value": None,
            "source": None,
            "notes": "",
        },
        "ph": {
            "reviewed": False,
            "status": None,
            "value": None,
            "source": None,
            "notes": "",
        },
        "notes": "",
    }


def write_review_record(path: Union[str, Path], record: Mapping[str, Any]) -> Path:
    """Write a review record JSON atomically and return the written path."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    _assert_json_object(record)
    tmp_path = out_path.with_name(f".{out_path.name}.{os.getpid()}.tmp")
    try:
        with tmp_path.open("w", encoding="utf-8") as handle:
            json.dump(dict(record), handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp_path, out_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    return out_path


def validate_review_record(
    record: Mapping[str, Any],
    *,
    production: bool = False,
) -> ValidationResult:
    """Return validation problems/warnings for an MMGBSA chemistry review record."""
    problems: list[str] = []
    warnings: list[str] = []

    for key in ("target_id", "ligand_id"):
        if not _nonempty_text(record.get(key)):
            problems.append(f"{key} is required")

    for field in REQUIRED_REVIEW_FIELDS:
        section = record.get(field)
        if not isinstance(section, Mapping):
            problems.append(f"{field} must be an object")
            continue
        _validate_review_section(field, section, problems, warnings, production)

    _validate_net_charge(record.get("net_charge"), problems)
    _validate_apo_holo(record.get("apo_holo"), problems, warnings)
    _validate_ph(record.get("ph"), problems, warnings)
    _validate_policy_list(
        record.get("water_policy"),
        "water_policy",
        "retained_waters",
        problems,
    )
    _validate_policy_list(
        record.get("metal_policy"),
        "metal_policy",
        "retained_metals",
        problems,
    )
    _validate_authority(record.get("input_chemistry_authority"), problems, production)
    _validate_suspicious_chemistry(
        record.get("suspicious_chemistry"), problems, warnings, production
    )

    notes = record.get("notes")
    if notes is not None and not isinstance(notes, str):
        warnings.append("notes should be a string")

    return {
        "valid": not problems,
        "production_ready": production and not problems,
        "problems": problems,
        "warnings": warnings,
    }


def _blank_review_section() -> ReviewRecord:
    return {
        "reviewed": False,
        "status": None,
        "source": None,
        "notes": "",
    }


def _validate_review_section(
    field: str,
    section: Mapping[str, Any],
    problems: list[str],
    warnings: list[str],
    production: bool,
) -> None:
    reviewed = section.get("reviewed")
    status = section.get("status")
    if not isinstance(reviewed, bool):
        problems.append(f"{field}.reviewed must be true or false")
    if status is not None and not isinstance(status, str):
        problems.append(f"{field}.status must be a string or null")
    if production:
        _validate_production_section(field, reviewed, status, problems)
        return
    if reviewed is True and str(status or "").lower() not in APPROVED_STATUSES:
        warnings.append(f"{field} is reviewed but not reviewed/approved/ok")


def _validate_production_section(
    field: str,
    reviewed: object,
    status: object,
    problems: list[str],
) -> None:
    if reviewed is not True:
        problems.append(f"{field} must be reviewed for production")
    if str(status or "").lower() not in APPROVED_STATUSES:
        problems.append(f"{field}.status must be reviewed, approved, or ok for production")


def _validate_net_charge(section: object, problems: list[str]) -> None:
    if not isinstance(section, Mapping):
        return
    value = section.get("value")
    if value is not None and not isinstance(value, int):
        problems.append("net_charge.value must be an integer or null")


def _validate_apo_holo(
    section: object,
    problems: list[str],
    warnings: list[str],
) -> None:
    if not isinstance(section, Mapping):
        return
    value = section.get("value")
    if value is None:
        warnings.append("apo_holo.value is not set")
        return
    if not isinstance(value, str) or value.upper() not in APO_HOLO_VALUES:
        problems.append("apo_holo.value must be APO or HOLO")


def _validate_ph(
    section: object,
    problems: list[str],
    warnings: list[str],
) -> None:
    if not isinstance(section, Mapping):
        return
    value = section.get("value")
    if value is None:
        warnings.append("ph.value is not set")
        return
    if not isinstance(value, (int, float)):
        problems.append("ph.value must be numeric or null")
        return
    if not 0.0 <= float(value) <= 14.0:
        problems.append("ph.value must be between 0 and 14")


def _validate_policy_list(
    section: object,
    section_name: str,
    list_name: str,
    problems: list[str],
) -> None:
    if not isinstance(section, Mapping):
        return
    values = section.get(list_name)
    if values is None:
        return
    if not isinstance(values, list):
        problems.append(f"{section_name}.{list_name} must be a list")
        return
    for idx, value in enumerate(values):
        if not isinstance(value, str):
            problems.append(f"{section_name}.{list_name}[{idx}] must be a string")


def _validate_authority(
    section: object, problems: list[str], production: bool
) -> None:
    if not production or not isinstance(section, Mapping):
        return
    if section.get("authoritative") is not True:
        problems.append("input_chemistry_authority.authoritative must be true for production")


def _validate_suspicious_chemistry(
    section: object,
    problems: list[str],
    warnings: list[str],
    production: bool,
) -> None:
    if not isinstance(section, Mapping):
        return
    severity = str(section.get("severity", "none") or "none").lower()
    reasons = section.get("reasons", [])
    if severity in {"warning", "blocker", "critical"}:
        warnings.append(f"suspicious_chemistry severity={severity}")
    if production and severity in {"blocker", "critical"}:
        if isinstance(reasons, list):
            reason_text = ",".join(str(reason) for reason in reasons[:8])
        else:
            reason_text = str(reasons)
        problems.append(f"suspicious_chemistry blocks production: {reason_text}")


def _assert_json_object(record: Mapping[str, Any]) -> None:
    try:
        json.dumps(dict(record), sort_keys=True)
    except TypeError as exc:
        raise TypeError("review record must be JSON-serializable") from exc


def _safe_path_token(value: str, field: str) -> str:
    token = str(value).strip()
    if not token:
        raise ValueError(f"{field} is required")
    if token in {".", ".."} or "/" in token or "\\" in token:
        raise ValueError(f"{field} must be a single path component")
    return token


def _nonempty_text(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())
