"""Fail-closed typed artifact contract for Atlas release archive indexes."""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Mapping


ARTIFACT_ROLE_SCOPES: dict[str, str] = {
    "source_receptor": "receptor",
    "prepared_receptor": "receptor",
    "source_ligand": "ligand",
    "prepared_ligand": "ligand",
    "docking_pose": "pair",
    "all_docking_poses": "pair",
    "all_pose_scores": "pair",
    "native_ligand": "receptor",
    "native_redock_pose": "receptor",
    "pose_image": "pair",
    "receptor_validation_report": "receptor",
    "native_redock_validation_report": "receptor",
    "pose_validation_report": "pair",
    "search_box_definition": "receptor",
    "run_configuration": "run",
    "release_manifest": "run",
    "raw_stdout": "pair",
    "raw_stderr": "pair",
    "software_environment": "run",
}
ARTIFACT_ROLES = tuple(sorted(ARTIFACT_ROLE_SCOPES))
ARTIFACT_PUBLICATION_POLICIES: dict[str, str] = {
    "source_receptor": "private",
    "prepared_receptor": "public_approved",
    "source_ligand": "private",
    "prepared_ligand": "public_approved",
    "docking_pose": "public_if_selected_for_release",
    "all_docking_poses": "private",
    "all_pose_scores": "private",
    "native_ligand": "public_approved",
    "native_redock_pose": "public_approved",
    "pose_image": "public_approved",
    "receptor_validation_report": "public_after_sanitization",
    "native_redock_validation_report": "public_after_sanitization",
    "pose_validation_report": "public_after_sanitization",
    "search_box_definition": "public_approved",
    "run_configuration": "public_after_sanitization",
    "release_manifest": "public_after_sanitization",
    "raw_stdout": "private",
    "raw_stderr": "private",
    "software_environment": "public_after_sanitization",
}
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")


class ArtifactContractError(ValueError):
    """Raised when an artifact archive index is unsafe or ambiguous."""


def artifact_scope(role: str) -> str:
    """Return the immutable association scope for one controlled role."""
    try:
        return ARTIFACT_ROLE_SCOPES[role]
    except KeyError as exc:
        raise ArtifactContractError(f"unknown artifact_role: {role!r}") from exc


def artifact_publication_policy(role: str) -> str:
    """Return the approved public-release policy for one controlled role."""
    try:
        return ARTIFACT_PUBLICATION_POLICIES[role]
    except KeyError as exc:
        raise ArtifactContractError(f"unknown artifact_role: {role!r}") from exc


def artifact_public_metadata_allowed(role: str, entry: Mapping[str, Any]) -> bool:
    """Return whether one artifact may appear in a sanitized public projection."""
    policy = artifact_publication_policy(role)
    if policy == "private":
        return False
    if policy == "public_if_selected_for_release":
        return entry.get("selected_for_release") is True
    return True


def artifact_public_record_allowed(
    role: str,
    entry: Mapping[str, Any],
    *,
    verified: Any,
    sha256: Any,
) -> bool:
    """Require public policy plus verified, content-addressed artifact identity."""
    if verified != 1 or not _SHA256.fullmatch(str(sha256 or "").strip()):
        return False
    return artifact_public_metadata_allowed(role, entry)


def validate_artifact_index(index: Any) -> list[str]:
    """Return structural errors for a typed release archive index."""
    if not isinstance(index, Mapping):
        return ["artifact index top-level value must be an object"]
    groups = index.get("groups")
    if not isinstance(groups, list):
        return ["artifact index groups must be a list"]

    errors: list[str] = []
    identities: set[tuple[str, str]] = set()
    for group_index, group in enumerate(groups):
        prefix = f"groups[{group_index}]"
        if not isinstance(group, Mapping):
            errors.append(f"{prefix} must be an object")
            continue
        group_key = group.get("group")
        if not isinstance(group_key, Mapping):
            errors.append(f"{prefix}.group must be an object")
        else:
            if not str(group_key.get("pdb_id") or "").strip():
                errors.append(f"{prefix}.group.pdb_id must be non-empty")
            for field in ("variant", "ph"):
                if field not in group_key:
                    errors.append(f"{prefix}.group.{field} must be present")
        if not isinstance(group.get("verified"), bool):
            errors.append(f"{prefix}.verified must be a boolean")
        entries = group.get("entries")
        if not isinstance(entries, list):
            errors.append(f"{prefix}.entries must be a list")
            continue
        for entry_index, entry in enumerate(entries):
            entry_prefix = f"{prefix}.entries[{entry_index}]"
            if not isinstance(entry, Mapping):
                errors.append(f"{entry_prefix} must be an object")
                continue
            role = str(entry.get("artifact_role") or "").strip()
            if role not in ARTIFACT_ROLE_SCOPES:
                allowed = ", ".join(ARTIFACT_ROLES)
                errors.append(f"{entry_prefix}.artifact_role must be one of: {allowed}")
                scope = ""
            else:
                scope = ARTIFACT_ROLE_SCOPES[role]
            if "selected_for_release" in entry and not isinstance(
                entry.get("selected_for_release"), bool
            ):
                errors.append(f"{entry_prefix}.selected_for_release must be a boolean")
            ligand = str(entry.get("ligand_canonical_id") or "").strip()
            if scope in {"ligand", "pair"} and not ligand:
                errors.append(
                    f"{entry_prefix}.ligand_canonical_id is required for "
                    f"artifact_role {role!r}"
                )
            if scope in {"run", "receptor"} and ligand:
                errors.append(
                    f"{entry_prefix}.ligand_canonical_id is not allowed for "
                    f"artifact_role {role!r}"
                )

            archive_path = str(entry.get("archive_path") or "").strip()
            member_name = str(entry.get("member_name") or "").strip()
            if not archive_path:
                errors.append(f"{entry_prefix}.archive_path must be non-empty")
            if not member_name:
                errors.append(f"{entry_prefix}.member_name must be non-empty")
            identity = (archive_path, member_name)
            if archive_path and member_name:
                if identity in identities:
                    errors.append(
                        f"{entry_prefix} duplicates archive_path/member_name: "
                        f"{archive_path!r}, {member_name!r}"
                    )
                identities.add(identity)

            digest = str(entry.get("sha256") or "").strip()
            if not _SHA256.fullmatch(digest):
                errors.append(f"{entry_prefix}.sha256 must be 64 hexadecimal digits")
            size = entry.get("size_bytes")
            if isinstance(size, bool) or not isinstance(size, int) or size < 0:
                errors.append(
                    f"{entry_prefix}.size_bytes must be a non-negative integer"
                )
            if not str(entry.get("file_type") or "").strip():
                errors.append(f"{entry_prefix}.file_type must be non-empty")
    return errors


def load_artifact_index(path: Path) -> dict[str, Any]:
    """Load and validate one release artifact index without silent fallback."""
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ArtifactContractError(
            f"cannot parse artifact index JSON ({exc})"
        ) from exc
    errors = validate_artifact_index(value)
    if errors:
        raise ArtifactContractError("; ".join(errors))
    return dict(value)


def artifact_role_counts(index: Mapping[str, Any]) -> dict[str, int]:
    """Summarize controlled roles in an already validated index."""
    counts: Counter[str] = Counter()
    for group in index.get("groups", []):
        for entry in group.get("entries", []):
            counts[str(entry["artifact_role"])] += 1
    return dict(sorted(counts.items()))
