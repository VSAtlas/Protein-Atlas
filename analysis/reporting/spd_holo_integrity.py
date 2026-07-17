"""Post-run, structured-evidence integrity audit for the frozen SPD protein panel."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from analysis.reporting.manifest_utils import load_run_manifest
from config.output_paths import run_dir_candidates, run_output_dir


SCHEMA_VERSION = "atlas.spd_holo_integrity.v1"
EXPECTED_COHORT_COUNT = 93
EXPECTED_STRICT_HOLO_COUNT = 76
EXPECTED_APO_FALLBACK_COUNT = 17
DEFAULT_SELECTION = Path("analysis/gene_list/spd_targets_final_93_v2_selected.csv")
DEFAULT_STRICT_MANIFEST = Path("analysis/gene_list/spd_targets_strict_manifest_final.csv")
_RUN_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
_PREP_COMPLETE = {"complete", "completed", "success", "succeeded"}
_DIRECT_ARTIFACT_DIRS = {"ligands_raw", "receptor", "work"}

CSV_FIELDS = (
    "schema_version",
    "run_id",
    "gene",
    "pdb_id",
    "expected_input_class",
    "expected_nontrivial_comp_ids",
    "manifest_entry_keys",
    "manifest_entry_count",
    "manifest_prep_statuses",
    "prep_complete",
    "processed_target_dir",
    "artifact_layouts",
    "ligand_component_paths",
    "extracted_component_ids",
    "expected_component_matches",
    "receptor_pdbqt_paths",
    "receptor_pdbqt_count",
    "meeko_sidecar_paths",
    "meeko_status",
    "retention_sidecar_path",
    "retention_status",
    "retained_after_cofactors",
    "retained_after_metals",
    "retention_missing_after_count",
    "retention_candidate_ids",
    "prep_integrity_certified",
    "experimentally_holo_certified",
    "cohort_requirement_certified",
    "certification_status",
    "reason_codes",
    "evidence_sha256s",
)


class HoloIntegrityInputError(ValueError):
    """Raised before outputs are written when the audit input is not authoritative."""


@dataclass(frozen=True)
class HoloIntegrityResult:
    csv_path: Path
    json_path: Path
    counts: Mapping[str, int]
    exit_code: int


@dataclass(frozen=True)
class _CohortRow:
    gene: str
    pdb_id: str
    expected_input_class: str
    expected_components: tuple[str, ...]


@dataclass(frozen=True)
class _ArtifactLineage:
    selected_layout: Path | None
    discovered_layouts: tuple[Path, ...]
    manifest_receptors: tuple[Path, ...]
    issues: tuple[str, ...]


def audit_spd_holo_integrity(
    *,
    repo_root: Path,
    run_id: str,
    selection_csv: Path,
    strict_manifest_csv: Path,
    output_dir: Path | None = None,
    strict: bool = True,
) -> HoloIntegrityResult:
    """Audit prepared SPD targets using only run-scoped structured evidence."""

    root = repo_root.expanduser().resolve()
    token = str(run_id).strip()
    if not _RUN_ID_RE.fullmatch(token):
        raise HoloIntegrityInputError(f"invalid run ID: {run_id!r}")
    selection = _resolve_input(root, selection_csv)
    strict_manifest = _resolve_input(root, strict_manifest_csv)
    cohort = _load_frozen_cohort(selection, strict_manifest)

    manifest, manifest_path = load_run_manifest(root, token)
    if not isinstance(manifest, dict) or manifest_path is None:
        raise HoloIntegrityInputError(f"run manifest is missing or unreadable: {token}")
    manifest_run_id = str(manifest.get("run_id") or "").strip()
    if manifest_run_id != token:
        raise HoloIntegrityInputError(
            f"run manifest ID mismatch: expected={token} actual={manifest_run_id or '<empty>'}"
        )
    proteins = manifest.get("proteins")
    if not isinstance(proteins, dict):
        raise HoloIntegrityInputError("run manifest proteins section is not a mapping")
    processed_root = _resolve_processed_root(root, token, manifest)
    if processed_root is None:
        raise HoloIntegrityInputError(
            f"no run-scoped processed-protein root exists for run {token}"
        )

    entries_by_pdb = _group_manifest_entries(proteins)
    rows = [
        _audit_target(
            cohort_row,
            repo_root=root,
            run_id=token,
            processed_root=processed_root,
            manifest_entries=entries_by_pdb.get(cohort_row.pdb_id, ()),
        )
        for cohort_row in cohort
    ]
    evidence_count = sum(int(row["structured_evidence_count"]) for row in rows)
    if evidence_count == 0:
        raise HoloIntegrityInputError(
            "no cohort-matching structured manifest or prepared artifacts were found"
        )

    counts = _counts(rows)
    out_dir = (
        output_dir.expanduser().resolve()
        if output_dir is not None
        else run_output_dir(root, "data", token) / "holo_integrity"
    )
    if out_dir.exists() and not out_dir.is_dir():
        raise HoloIntegrityInputError(f"holo-integrity output is not a directory: {out_dir}")
    csv_bytes = _csv_bytes(rows)
    csv_path = out_dir / "holo_integrity.csv"
    json_path = out_dir / "holo_integrity.json"
    payload = {
        "schema_version": SCHEMA_VERSION,
        "run_id": token,
        "inputs": {
            "selection_csv": str(selection),
            "selection_sha256": _sha256(selection),
            "strict_manifest_csv": str(strict_manifest),
            "strict_manifest_sha256": _sha256(strict_manifest),
            "run_manifest": str(manifest_path.resolve()),
            "run_manifest_sha256": _sha256(manifest_path),
            "processed_pdb_dir": str(processed_root),
        },
        "outputs": {
            "csv": str(csv_path),
            "csv_sha256": hashlib.sha256(csv_bytes).hexdigest(),
            "json": str(json_path),
        },
        "counts": counts,
        "policy": {
            "cohort_is_frozen": True,
            "expected_cohort_count": EXPECTED_COHORT_COUNT,
            "expected_strict_holo_count": EXPECTED_STRICT_HOLO_COUNT,
            "expected_documented_apo_fallback_count": EXPECTED_APO_FALLBACK_COUNT,
            "strict_holo_requires_expected_component_extraction": True,
            "prep_complete_and_nonempty_receptor_required": True,
            "retention_gap_inconsistent_or_unreadable_blocks": True,
            "absent_retention_sidecar_is_nonblocking_not_reported": True,
            "documented_apo_fallback_is_never_experimentally_holo_certified": True,
            "raw_logs_read": False,
        },
        "rows": rows,
    }
    json_bytes = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    _write_atomic_pair(out_dir, csv_path, csv_bytes, json_path, json_bytes)
    exit_code = 1 if strict and counts["cohort_requirement_uncertified"] else 0
    return HoloIntegrityResult(csv_path, json_path, counts, exit_code)


def _resolve_input(repo_root: Path, path: Path) -> Path:
    candidate = path.expanduser()
    if not candidate.is_absolute():
        candidate = repo_root / candidate
    candidate = candidate.resolve()
    if not candidate.is_file() or candidate.stat().st_size <= 0:
        raise HoloIntegrityInputError(f"required SPD cohort input is missing: {candidate}")
    return candidate


def _load_frozen_cohort(selection: Path, strict_manifest: Path) -> list[_CohortRow]:
    selected = _read_csv(selection)
    strict_rows = _read_csv(strict_manifest)
    selected_by_pair: dict[tuple[str, str], Mapping[str, str]] = {}
    for row in selected:
        pair = (_clean(row.get("gene")).upper(), _clean(row.get("pdb_id")).upper())
        if not all(pair) or pair in selected_by_pair:
            raise HoloIntegrityInputError("selection must contain unique nonempty gene/PDB pairs")
        selected_by_pair[pair] = row
    strict_by_pair: dict[tuple[str, str], Mapping[str, str]] = {}
    for row in strict_rows:
        if _clean(row.get("status")) != "strict_holo_downloaded":
            continue
        if _clean(row.get("downloaded")).casefold() != "yes":
            continue
        pair = (_clean(row.get("gene")).upper(), _clean(row.get("pdb_id")).upper())
        if not all(pair) or pair in strict_by_pair:
            raise HoloIntegrityInputError(
                "strict manifest must contain unique nonempty downloaded gene/PDB pairs"
            )
        strict_by_pair[pair] = row
    selected_pairs = set(selected_by_pair)
    strict_pairs = set(strict_by_pair)
    fallback_pairs = selected_pairs - strict_pairs
    if (
        len(selected_pairs) != EXPECTED_COHORT_COUNT
        or len(strict_pairs) != EXPECTED_STRICT_HOLO_COUNT
        or len(fallback_pairs) != EXPECTED_APO_FALLBACK_COUNT
        or not strict_pairs.issubset(selected_pairs)
    ):
        raise HoloIntegrityInputError(
            "frozen SPD cohort drifted: "
            f"selection={len(selected_pairs)} strict={len(strict_pairs)} "
            f"fallback={len(fallback_pairs)}"
        )
    cohort: list[_CohortRow] = []
    for pair in sorted(selected_pairs):
        selected_row = selected_by_pair[pair]
        if pair in strict_by_pair:
            source = strict_by_pair[pair]
            input_class = "strict_holo"
        else:
            source = selected_row
            input_class = "documented_apo_fallback"
        components = _tokens(source.get("nontrivial_comp_ids"))
        if input_class == "strict_holo" and not components:
            raise HoloIntegrityInputError(
                f"strict holo pair lacks expected component IDs: {pair[0]}/{pair[1]}"
            )
        cohort.append(_CohortRow(pair[0], pair[1], input_class, components))
    return cohort


def _resolve_processed_root(
    repo_root: Path, run_id: str, manifest: Mapping[str, Any]
) -> Path | None:
    raw_paths = manifest.get("paths")
    candidates: list[Path] = []
    if isinstance(raw_paths, dict) and raw_paths.get("processed_pdb_dir"):
        raw = Path(str(raw_paths["processed_pdb_dir"])).expanduser()
        candidate = raw if raw.is_absolute() else repo_root / raw
        if candidate.name == run_id:
            candidates.append(candidate)
    candidates.extend(run_dir_candidates(repo_root, "processed_pdbs", run_id))
    seen: set[str] = set()
    for raw in candidates:
        candidate = raw.expanduser().resolve()
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.name == run_id and candidate.is_dir():
            return candidate
    return None


def _group_manifest_entries(
    proteins: Mapping[str, Any],
) -> dict[str, tuple[tuple[str, Mapping[str, Any]], ...]]:
    grouped: defaultdict[str, list[tuple[str, Mapping[str, Any]]]] = defaultdict(list)
    for key, raw_entry in proteins.items():
        if not isinstance(raw_entry, dict):
            continue
        key_pdb = str(key).split("|", 1)[0].strip().upper()
        field_pdb = _clean(raw_entry.get("pdb_id")).upper()
        pdb_id = field_pdb or key_pdb
        if pdb_id:
            grouped[pdb_id].append((str(key), raw_entry))
    return {
        pdb_id: tuple(sorted(entries, key=lambda item: item[0]))
        for pdb_id, entries in grouped.items()
    }


def _audit_target(
    cohort: _CohortRow,
    *,
    repo_root: Path,
    run_id: str,
    processed_root: Path,
    manifest_entries: Sequence[tuple[str, Mapping[str, Any]]],
) -> dict[str, Any]:
    reasons: list[str] = []
    entry_keys = [key for key, _entry in manifest_entries]
    prep_statuses = [_prep_status(entry) for _key, entry in manifest_entries]
    prep_complete = bool(prep_statuses) and all(
        status in _PREP_COMPLETE for status in prep_statuses
    )
    if not manifest_entries:
        reasons.append("manifest_entry_missing")
    elif not prep_complete:
        reasons.append("prep_not_complete")

    target_root = processed_root / cohort.pdb_id
    lineage = _resolve_artifact_lineage(
        repo_root, processed_root, target_root, manifest_entries
    )
    reasons.extend(lineage.issues)
    selected_layouts = (
        [lineage.selected_layout] if lineage.selected_layout is not None else []
    )
    reported_layouts = selected_layouts or list(lineage.discovered_layouts)
    ligand_paths, ligand_issues = _artifact_files(
        selected_layouts, "ligands_raw", "*.pdb"
    )
    ligand_paths = [
        path for path in ligand_paths if "sanitized" not in path.stem.casefold()
    ]
    receptor_paths, receptor_issues = _artifact_files(
        selected_layouts, "receptor", "*.pdbqt", recursive=True
    )
    meeko_paths, meeko_path_issues = _artifact_files(
        selected_layouts, "receptor", "*.meeko_input_stage.json", recursive=True
    )
    retention_paths, retention_path_issues = _artifact_files(
        selected_layouts, "work", "*.retained_hets_audit.json"
    )
    reasons.extend(ligand_issues)
    reasons.extend(receptor_issues)
    reasons.extend(meeko_path_issues)
    reasons.extend(retention_path_issues)

    component_ids: list[str] = []
    for path in ligand_paths:
        component = _first_component_id(path)
        if component is None:
            reasons.append("ligand_component_unreadable")
        else:
            component_ids.append(component)
    component_ids = sorted(set(component_ids))
    expected_matches = sorted(set(cohort.expected_components) & set(component_ids))
    if cohort.expected_input_class == "strict_holo" and not expected_matches:
        reasons.append("strict_expected_component_not_extracted")

    receptor_paths = [path for path in receptor_paths if path.stat().st_size > 0]
    if not receptor_paths:
        reasons.append("nonempty_receptor_missing")
    if lineage.manifest_receptors and not set(lineage.manifest_receptors).issubset(
        receptor_paths
    ):
        reasons.append("manifest_receptor_artifact_missing")
    lineage_root = lineage.selected_layout or target_root
    meeko_status, meeko_issues = _meeko_status(
        meeko_paths, receptor_paths, lineage_root
    )
    reasons.extend(meeko_issues)
    retention = _retention_status(
        retention_paths, layout_root=lineage.selected_layout
    )
    reasons.extend(retention["issues"])

    blocking = sorted(set(reasons))
    prep_integrity = prep_complete and bool(receptor_paths) and not blocking
    experimental_holo = bool(
        cohort.expected_input_class == "strict_holo"
        and prep_integrity
        and expected_matches
    )
    requirement_certified = (
        experimental_holo
        if cohort.expected_input_class == "strict_holo"
        else prep_integrity
    )
    if requirement_certified and cohort.expected_input_class == "strict_holo":
        certification = "certified_strict_holo"
    elif requirement_certified:
        certification = "certified_documented_apo_fallback_prep_only"
    else:
        certification = "uncertified"

    evidence_paths = sorted(
        {*ligand_paths, *receptor_paths, *meeko_paths, *retention["evidence_paths"]},
        key=lambda path: str(path),
    )
    evidence_hashes = {str(path): _sha256(path) for path in evidence_paths}
    return {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "gene": cohort.gene,
        "pdb_id": cohort.pdb_id,
        "expected_input_class": cohort.expected_input_class,
        "expected_nontrivial_comp_ids": list(cohort.expected_components),
        "manifest_entry_keys": entry_keys,
        "manifest_entry_count": len(entry_keys),
        "manifest_prep_statuses": sorted(set(prep_statuses)),
        "prep_complete": prep_complete,
        "processed_target_dir": str(target_root),
        "artifact_layouts": [str(path) for path in reported_layouts],
        "ligand_component_paths": [str(path) for path in ligand_paths],
        "extracted_component_ids": component_ids,
        "expected_component_matches": expected_matches,
        "receptor_pdbqt_paths": [str(path) for path in receptor_paths],
        "receptor_pdbqt_count": len(receptor_paths),
        "meeko_sidecar_paths": [str(path) for path in meeko_paths],
        "meeko_status": meeko_status,
        "retention_sidecar_path": retention["path"],
        "retention_status": retention["status"],
        "retained_after_cofactors": retention["after_cofactors"],
        "retained_after_metals": retention["after_metals"],
        "retention_missing_after_count": retention["missing_after_count"],
        "retention_candidate_ids": retention["candidate_ids"],
        "prep_integrity_certified": prep_integrity,
        "experimentally_holo_certified": experimental_holo,
        "cohort_requirement_certified": requirement_certified,
        "certification_status": certification,
        "reason_codes": blocking,
        "evidence_sha256s": evidence_hashes,
        "structured_evidence_count": len(entry_keys) + len(evidence_paths),
    }


def _prep_status(entry: Mapping[str, Any]) -> str:
    stages = entry.get("stages")
    prep = stages.get("prep") if isinstance(stages, dict) else None
    return _clean(prep.get("status") if isinstance(prep, dict) else None).casefold()


def _resolve_artifact_lineage(
    repo_root: Path,
    processed_root: Path,
    target_root: Path,
    manifest_entries: Sequence[tuple[str, Mapping[str, Any]]],
) -> _ArtifactLineage:
    layouts, layout_issues = _artifact_layouts(target_root)
    issues = list(layout_issues)
    if not manifest_entries:
        return _ArtifactLineage(None, tuple(layouts), (), tuple(issues))

    bindings, binding_issues = _manifest_layout_bindings(
        repo_root,
        processed_root,
        target_root,
        layouts,
        manifest_entries,
    )
    issues.extend(binding_issues)
    manifest_receptors = tuple(
        sorted({receptor for receptor, _layout in bindings}, key=str)
    )
    if binding_issues:
        return _ArtifactLineage(
            None, tuple(layouts), manifest_receptors, tuple(issues)
        )

    bound_layouts = {layout for _receptor, layout in bindings}
    if len(bound_layouts) != 1:
        issue = (
            "artifact_layout_ambiguous"
            if len(bound_layouts) > 1
            else "manifest_receptor_layout_unresolved"
        )
        issues.append(issue)
        return _ArtifactLineage(
            None, tuple(layouts), manifest_receptors, tuple(issues)
        )
    selected = next(iter(bound_layouts))
    return _ArtifactLineage(
        selected, tuple(layouts), manifest_receptors, tuple(issues)
    )


def _manifest_layout_bindings(
    repo_root: Path,
    processed_root: Path,
    target_root: Path,
    layouts: Sequence[Path],
    manifest_entries: Sequence[tuple[str, Mapping[str, Any]]],
) -> tuple[list[tuple[Path, Path]], list[str]]:
    bindings: list[tuple[Path, Path]] = []
    issues: list[str] = []
    for _key, entry in manifest_entries:
        raw_path = _manifest_receptor_path(entry)
        if not raw_path:
            issues.append("manifest_receptor_path_missing")
            continue
        matches: dict[tuple[str, str], tuple[Path, Path]] = {}
        for candidate in _manifest_path_candidates(
            raw_path, repo_root, processed_root, target_root
        ):
            for layout in layouts:
                if _is_within(candidate, layout / "receptor"):
                    matches[(str(candidate), str(layout))] = (candidate, layout)
        if not matches:
            issues.append("manifest_receptor_layout_unresolved")
            continue
        if len(matches) > 1:
            issues.append("manifest_receptor_layout_ambiguous")
            continue
        bindings.append(next(iter(matches.values())))
    return bindings, issues


def _manifest_receptor_path(entry: Mapping[str, Any]) -> str:
    stages = entry.get("stages")
    prep = stages.get("prep") if isinstance(stages, Mapping) else None
    details = prep.get("details") if isinstance(prep, Mapping) else None
    if not isinstance(details, Mapping):
        return ""
    return _clean(details.get("receptor_pdbqt"))


def _manifest_path_candidates(
    raw_path: str,
    repo_root: Path,
    processed_root: Path,
    target_root: Path,
) -> tuple[Path, ...]:
    raw = Path(raw_path).expanduser()
    candidates = (
        (raw,)
        if raw.is_absolute()
        else (
            target_root / raw,
            processed_root / raw,
            repo_root / raw,
        )
    )
    resolved: dict[str, Path] = {}
    for candidate in candidates:
        try:
            path = candidate.resolve()
        except (OSError, RuntimeError):
            continue
        resolved[str(path)] = path
    return tuple(resolved.values())


def _artifact_layouts(target_root: Path) -> tuple[list[Path], list[str]]:
    if not target_root.is_dir():
        return [], []
    layouts: list[Path] = [target_root]
    issues: list[str] = []
    try:
        children = sorted(target_root.iterdir(), key=lambda path: path.name.casefold())
    except OSError:
        return [], ["processed_target_unreadable"]
    for child in children:
        if not child.is_dir() or child.name in _DIRECT_ARTIFACT_DIRS:
            continue
        if child.is_symlink() or not _is_within(child, target_root):
            if any((child / name).exists() for name in _DIRECT_ARTIFACT_DIRS):
                issues.append("artifact_layout_symlink_or_escape")
            continue
        if any((child / name).is_dir() for name in _DIRECT_ARTIFACT_DIRS):
            layouts.append(child)
    return layouts, issues


def _artifact_files(
    layouts: Sequence[Path],
    subdir: str,
    pattern: str,
    *,
    recursive: bool = False,
) -> tuple[list[Path], list[str]]:
    paths: dict[str, Path] = {}
    issues: list[str] = []
    for layout in layouts:
        directory = layout / subdir
        if not directory.exists():
            continue
        if directory.is_symlink() or not directory.is_dir() or not _is_within(directory, layout):
            issues.append("artifact_directory_inconsistent")
            continue
        try:
            candidates = directory.rglob(pattern) if recursive else directory.glob(pattern)
            for path in candidates:
                if path.is_symlink() or not path.is_file() or not _is_within(path, directory):
                    issues.append("artifact_file_symlink_or_escape")
                    continue
                paths[str(path.resolve())] = path.resolve()
        except OSError:
            issues.append("artifact_directory_unreadable")
    return sorted(paths.values(), key=lambda path: str(path)), issues


def _first_component_id(path: Path) -> str | None:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(("ATOM  ", "HETATM")):
                    component = line[17:20].strip().upper()
                    return component or None
    except OSError:
        return None
    return None


def _meeko_status(
    sidecars: Sequence[Path], receptors: Sequence[Path], layout_root: Path
) -> tuple[str, list[str]]:
    if not sidecars:
        return "not_reported", []
    receptor_set = {path.resolve() for path in receptors}
    issues: list[str] = []
    unreadable = False
    for path in sidecars:
        try:
            payload = _read_json(path)
        except ValueError:
            unreadable = True
            continue
        if payload.get("ok") is not True:
            issues.append("meeko_sidecar_not_ok")
        output_raw = _clean(payload.get("output_pdbqt"))
        if not output_raw:
            issues.append("meeko_sidecar_output_missing")
            continue
        output = Path(output_raw).expanduser()
        if not output.is_absolute():
            output = layout_root / output
        output = output.resolve()
        if output not in receptor_set or not _is_within(output, layout_root):
            issues.append("meeko_sidecar_output_inconsistent")
    if unreadable:
        return "unreadable", [*issues, "meeko_sidecar_unreadable"]
    if issues:
        return "inconsistent", issues
    return "ok", []


def _retention_status(
    paths: Sequence[Path], *, layout_root: Path | None = None
) -> dict[str, Any]:
    empty: dict[str, Any] = {
        "path": "",
        "status": "not_reported",
        "after_cofactors": None,
        "after_metals": None,
        "missing_after_count": None,
        "candidate_ids": [],
        "issues": [],
        "evidence_paths": [],
    }
    if not paths:
        return empty
    latest = max(paths, key=lambda path: (path.stat().st_mtime_ns, str(path)))
    try:
        payload = _read_json(latest)
    except ValueError:
        return {
            **empty,
            "path": str(latest),
            "status": "unreadable",
            "issues": ["retention_sidecar_unreadable"],
            "evidence_paths": [latest],
        }
    if layout_root is not None and not _retention_payload_matches_layout(
        payload, layout_root
    ):
        return {
            **empty,
            "path": str(latest),
            "status": "inconsistent",
            "issues": ["retention_sidecar_layout_inconsistent"],
            "evidence_paths": [latest],
        }
    after = payload.get("after")
    candidates = payload.get("candidates")
    missing = payload.get("missing_after_count")
    if (
        not isinstance(after, dict)
        or not _nonnegative_int(after.get("cofactors"))
        or not _nonnegative_int(after.get("metals"))
        or not _nonnegative_int(missing)
        or not isinstance(candidates, list)
    ):
        return {
            **empty,
            "path": str(latest),
            "status": "inconsistent",
            "issues": ["retention_sidecar_inconsistent"],
            "evidence_paths": [latest],
        }
    candidate_ids = sorted(
        {
            _clean(item.get("token") or item.get("resname")).upper()
            for item in candidates
            if isinstance(item, dict) and _clean(item.get("token") or item.get("resname"))
        }
    )
    missing_count = int(str(missing))
    status = "retention_gap" if missing_count else "ok"
    issues = ["retention_gap"] if missing_count else []
    return {
        "path": str(latest),
        "status": status,
        "after_cofactors": int(str(after["cofactors"])),
        "after_metals": int(str(after["metals"])),
        "missing_after_count": missing_count,
        "candidate_ids": candidate_ids,
        "issues": issues,
        "evidence_paths": [latest],
    }


def _retention_payload_matches_layout(
    payload: Mapping[str, Any], layout_root: Path
) -> bool:
    for field in ("source_pdb", "target_pdb"):
        raw = _clean(payload.get(field))
        if not raw:
            return False
        path = Path(raw).expanduser()
        if not path.is_absolute():
            path = layout_root / path
        if not _is_within(path, layout_root):
            return False
    return True


def _counts(rows: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    classes = Counter(str(row["expected_input_class"]) for row in rows)
    prep_certified = sum(bool(row["prep_integrity_certified"]) for row in rows)
    experimental = sum(bool(row["experimentally_holo_certified"]) for row in rows)
    requirement = sum(bool(row["cohort_requirement_certified"]) for row in rows)
    return {
        "cohort_rows": len(rows),
        "strict_holo_expected": classes["strict_holo"],
        "documented_apo_fallback_expected": classes["documented_apo_fallback"],
        "prep_integrity_certified": prep_certified,
        "prep_integrity_uncertified": len(rows) - prep_certified,
        "experimentally_holo_certified": experimental,
        "documented_apo_fallback_experimentally_holo_certified": sum(
            bool(row["experimentally_holo_certified"])
            for row in rows
            if row["expected_input_class"] == "documented_apo_fallback"
        ),
        "cohort_requirement_certified": requirement,
        "cohort_requirement_uncertified": len(rows) - requirement,
    }


def _csv_bytes(rows: Sequence[Mapping[str, Any]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(CSV_FIELDS), extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        flat = dict(row)
        for field in (
            "expected_nontrivial_comp_ids",
            "manifest_entry_keys",
            "manifest_prep_statuses",
            "artifact_layouts",
            "ligand_component_paths",
            "extracted_component_ids",
            "expected_component_matches",
            "receptor_pdbqt_paths",
            "meeko_sidecar_paths",
            "retention_candidate_ids",
            "reason_codes",
        ):
            flat[field] = ";".join(str(value) for value in row.get(field, []))
        hashes = row.get("evidence_sha256s") or {}
        flat["evidence_sha256s"] = ";".join(
            f"{path}={digest}" for path, digest in sorted(hashes.items())
        )
        for field in (
            "prep_complete",
            "prep_integrity_certified",
            "experimentally_holo_certified",
            "cohort_requirement_certified",
        ):
            flat[field] = "true" if row.get(field) else "false"
        writer.writerow(flat)
    return buffer.getvalue().encode("utf-8")


def _write_atomic_pair(
    output_dir: Path,
    csv_path: Path,
    csv_bytes: bytes,
    json_path: Path,
    json_bytes: bytes,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_tmp = csv_path.with_suffix(csv_path.suffix + ".part")
    json_tmp = json_path.with_suffix(json_path.suffix + ".part")
    if csv_tmp.exists() or json_tmp.exists():
        raise HoloIntegrityInputError(
            f"stale partial holo-integrity output exists under {output_dir}"
        )
    try:
        _write_fsynced(csv_tmp, csv_bytes)
        _write_fsynced(json_tmp, json_bytes)
        os.replace(csv_tmp, csv_path)
        os.replace(json_tmp, json_path)
    except Exception:
        for path in (csv_tmp, json_tmp):
            if path.exists():
                path.unlink()
        raise


def _write_fsynced(path: Path, content: bytes) -> None:
    with path.open("xb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise HoloIntegrityInputError(f"CSV has no header: {path}")
        return [dict(row) for row in reader]


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"invalid JSON: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON root is not an object: {path}")
    return value


def _tokens(value: object) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                token.strip().upper()
                for token in re.split(r"[;,|\s]+", _clean(value))
                if token.strip()
            }
        )
    )


def _nonnegative_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _is_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
    except (OSError, ValueError):
        return False
    return True


def _clean(value: object) -> str:
    text = str(value or "").strip()
    return "" if text.casefold() == "nan" else text


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = [
    "DEFAULT_SELECTION",
    "DEFAULT_STRICT_MANIFEST",
    "HoloIntegrityInputError",
    "HoloIntegrityResult",
    "audit_spd_holo_integrity",
]
