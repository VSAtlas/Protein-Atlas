"""Release-manifest loading and structural validation."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Mapping

import yaml  # type: ignore[import-untyped]

REQUIRED_POLICY_NAMES = (
    "receptor_selection",
    "native_redocking",
    "failure_handling",
    "normalization",
)
_SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


class ReleaseManifestError(ValueError):
    """Raised when an Atlas release manifest is invalid."""


def load_release_manifest(path: Path) -> dict[str, Any]:
    """Load a JSON or YAML release manifest and validate its structure."""
    try:
        with path.open("r", encoding="utf-8") as handle:
            if path.suffix.lower() == ".json":
                loaded = json.load(handle)
            else:
                loaded = yaml.safe_load(handle)
    except (OSError, json.JSONDecodeError, yaml.YAMLError) as exc:
        raise ReleaseManifestError(f"could not read release manifest {path}: {exc}") from exc
    if not isinstance(loaded, dict):
        raise ReleaseManifestError("release manifest must be a mapping")
    errors = validate_release_manifest(loaded)
    if errors:
        raise ReleaseManifestError("; ".join(errors))
    return loaded


def _run_id(entry: Any) -> str:
    if isinstance(entry, str):
        return entry.strip()
    if isinstance(entry, Mapping):
        return str(entry.get("run_id") or "").strip()
    return ""


def validate_release_manifest(manifest: Mapping[str, Any]) -> list[str]:
    """Return all structural errors without touching runtime artifacts."""
    errors: list[str] = []
    release_id = str(manifest.get("release_id") or "").strip()
    if not release_id:
        errors.append("release_id is required")
    elif not _SAFE_ID.fullmatch(release_id):
        errors.append("release_id must use only letters, numbers, '.', '_', or '-'")

    schema_version = manifest.get("schema_version")
    if schema_version not in (1, "1", "1.0"):
        errors.append("schema_version must be 1")

    runs = manifest.get("runs")
    if not isinstance(runs, list) or not runs:
        errors.append("runs must be a non-empty list")
    else:
        seen: set[str] = set()
        for index, entry in enumerate(runs):
            run_id = _run_id(entry)
            if not run_id:
                errors.append(f"runs[{index}] must be a run ID or mapping with run_id")
            elif run_id in seen:
                errors.append(f"duplicate run_id: {run_id}")
            else:
                seen.add(run_id)

    policies = manifest.get("scientific_policies")
    if not isinstance(policies, Mapping):
        errors.append("scientific_policies must be a mapping")
    else:
        for name in REQUIRED_POLICY_NAMES:
            value = policies.get(name)
            if not isinstance(value, Mapping) or not value:
                errors.append(f"scientific_policies.{name} must be a non-empty mapping")

    return errors


def iter_run_entries(manifest: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Normalize string and mapping run declarations."""
    result: list[dict[str, Any]] = []
    for value in manifest.get("runs", []):
        if isinstance(value, str):
            result.append({"run_id": value})
        else:
            result.append(dict(value))
    return result
