"""Release-manifest loading and structural validation."""

from __future__ import annotations

import json
import re
from pathlib import PurePosixPath
from pathlib import Path
from typing import Any, Mapping

import yaml  # type: ignore[import-untyped]

from analysis.atlas_database.annotations import ANNOTATION_SOURCE_NAMES

REQUIRED_POLICY_NAMES = (
    "receptor_selection",
    "native_redocking",
    "failure_handling",
    "normalization",
)
_SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_ATTEMPT_SELECTOR_FIELDS = {
    "pdb_id",
    "variant",
    "ph_label",
    "ligand_canonical_id",
    "completion_relpath",
    "completion_sha256",
    "result_row_number",
    "result_sha256",
}


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
        raise ReleaseManifestError(
            f"could not read release manifest {path}: {exc}"
        ) from exc
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


def _selector_context_token(value: Any) -> str:
    token = str(value or "").strip()
    if token.lower() in {"base", "none", "null"} or token.lower().endswith(".csv"):
        return ""
    return token


def _validate_attempt_selections(
    entry: Mapping[str, Any], run_index: int, errors: list[str]
) -> None:
    selectors = entry.get("attempt_selections")
    if selectors is None:
        return
    prefix = f"runs[{run_index}].attempt_selections"
    if not isinstance(selectors, list):
        errors.append(f"{prefix} must be a list")
        return
    seen: set[tuple[str, str, str, str]] = set()
    for selector_index, selector in enumerate(selectors):
        item_prefix = f"{prefix}[{selector_index}]"
        if not isinstance(selector, Mapping):
            errors.append(f"{item_prefix} must be a mapping")
            continue
        unknown = sorted(set(selector) - _ATTEMPT_SELECTOR_FIELDS)
        for field in unknown:
            errors.append(f"{item_prefix} has unknown field: {field}")
        for field in ("pdb_id", "variant", "ph_label", "ligand_canonical_id"):
            if field not in selector:
                errors.append(f"{item_prefix}.{field} must be present")
        pdb_id = str(selector.get("pdb_id") or "").strip().upper()
        ligand = str(selector.get("ligand_canonical_id") or "").strip()
        if not pdb_id:
            errors.append(f"{item_prefix}.pdb_id must be non-empty")
        if not ligand:
            errors.append(f"{item_prefix}.ligand_canonical_id must be non-empty")
        key = (
            pdb_id,
            str(selector.get("variant") or "").strip().upper(),
            _selector_context_token(selector.get("ph_label")),
            ligand,
        )
        if key in seen:
            errors.append(f"{item_prefix} duplicates selector key {key!r}")
        seen.add(key)

        completion_relpath = str(selector.get("completion_relpath") or "").strip()
        completion_sha256 = str(selector.get("completion_sha256") or "").strip()
        has_completion = bool(completion_relpath or completion_sha256)
        if has_completion and not completion_relpath:
            errors.append(f"{item_prefix}.completion_relpath is required")
        if has_completion and not _SHA256.fullmatch(completion_sha256):
            errors.append(
                f"{item_prefix}.completion_sha256 must be 64 hexadecimal digits"
            )
        if completion_relpath:
            path = PurePosixPath(completion_relpath)
            if (
                path.is_absolute()
                or ".." in path.parts
                or "\\" in completion_relpath
                or completion_relpath.endswith("/")
            ):
                errors.append(
                    f"{item_prefix}.completion_relpath must be a safe POSIX path "
                    "relative to paths.docked"
                )

        result_row = selector.get("result_row_number")
        result_sha256 = str(selector.get("result_sha256") or "").strip()
        has_result = result_row is not None or bool(result_sha256)
        if has_result and (
            isinstance(result_row, bool)
            or not isinstance(result_row, int)
            or result_row < 2
        ):
            errors.append(f"{item_prefix}.result_row_number must be an integer >= 2")
        if has_result and not _SHA256.fullmatch(result_sha256):
            errors.append(f"{item_prefix}.result_sha256 must be 64 hexadecimal digits")
        if not has_completion and not has_result:
            errors.append(
                f"{item_prefix} must select a completion attempt, a result attempt, or both"
            )


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
            if isinstance(entry, Mapping):
                _validate_attempt_selections(entry, index, errors)

    policies = manifest.get("scientific_policies")
    if not isinstance(policies, Mapping):
        errors.append("scientific_policies must be a mapping")
    else:
        for name in REQUIRED_POLICY_NAMES:
            value = policies.get(name)
            if not isinstance(value, Mapping) or not value:
                errors.append(f"scientific_policies.{name} must be a non-empty mapping")

    annotations = manifest.get("annotations")
    if annotations is not None:
        if not isinstance(annotations, Mapping):
            errors.append("annotations must be a mapping")
        else:
            unknown = sorted(set(annotations) - set(ANNOTATION_SOURCE_NAMES))
            for name in unknown:
                errors.append(f"unknown annotations source: {name}")
            for name in ANNOTATION_SOURCE_NAMES:
                if name not in annotations:
                    continue
                spec = annotations[name]
                if isinstance(spec, str):
                    valid = bool(spec.strip())
                elif isinstance(spec, Mapping):
                    valid = bool(str(spec.get("path") or "").strip())
                else:
                    valid = False
                if not valid:
                    errors.append(
                        f"annotations.{name} must be a path string or mapping with path"
                    )

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
