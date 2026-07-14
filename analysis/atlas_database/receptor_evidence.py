"""Evidence-only APO/HOLO and metal-audit extraction for Atlas releases."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from analysis.atlas_database.importer import (
    _context_values,
    _load_mapping,
    _run_paths,
)
from analysis.atlas_database.manifest import iter_run_entries, load_release_manifest
from config.output_paths import run_output_dir
from protein_prep.pdb_fixer_runtime import (
    get_atom_rules,
    load_canonical_cofactors,
    load_canonical_metals,
)
from protein_prep.pdb_records import (
    atom_name,
    chain_id,
    insertion_code,
    is_atom_record,
    line_element,
    line_xyz,
    pdbqt_line_element,
    pdbqt_line_xyz,
    record_type,
    resseq,
    residue_name,
)

CHEMISTRY_EVIDENCE_SCHEMA_VERSION = 2
CHEMISTRY_CLASSIFICATION_METHOD = "atlas_prepared_receptor_retained_chemistry_v2"
CONTROLLED_RECEPTOR_CLASSIFICATIONS = frozenset({"APO", "HOLO"})
_AUDIT_FILENAMES = (
    "metal_site_audit.json",
    "autodock4zn_plan.json",
    "metal_parameterization_plan.json",
    "metal_bound_het_audit.json",
    "metal_chemistry_audit.json",
)
_AUDIT_PATTERNS = (
    "**/metal_site_audit.json",
    "**/*retained_hets_audit.json",
    "**/audits/*_ions.json",
    "**/metal_bound_het_audit.json",
    "**/metal_chemistry_audit.json",
    "**/autodock4zn_plan.json",
    "**/metal_parameterization_plan.json",
)


class ReceptorEvidenceError(ValueError):
    """Raised when receptor evidence cannot be generated unambiguously."""


def _text(value: Any) -> str:
    return str(value).strip() if value is not None else ""


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _token_set(
    values: Iterable[str],
    normalize: Callable[[str], str] | None = None,
) -> set[str]:
    tokens: set[str] = set()
    for value in values:
        token = _text(value).upper()
        if not token:
            continue
        normalized = normalize(token) if normalize is not None else token
        if normalized:
            tokens.add(_text(normalized).upper())
    return tokens


def _policy_sha256(
    *,
    metals: Iterable[str],
    cofactors: Iterable[str],
    halides: Iterable[str],
    ambiguous_cations: Iterable[str],
    nucleotide_like: Iterable[str],
) -> str:
    payload = json.dumps(
        {
            "canonical_metals": sorted(set(metals)),
            "canonical_cofactors": sorted(set(cofactors)),
            "canonical_halides": sorted(set(halides)),
            "ambiguous_cations": sorted(set(ambiguous_cations)),
            "nucleotide_like": sorted(set(nucleotide_like)),
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _chemistry_sets(
    *,
    canonical_metals: Iterable[str] | None,
    canonical_cofactors: Iterable[str] | None,
) -> tuple[set[str], set[str], set[str], set[str], set[str], Callable[[str], str]]:
    rules = get_atom_rules()
    normalizer = getattr(rules, "normalize_resname", None)

    def normalize(value: str) -> str:
        token = _text(value).upper()
        if not token:
            return ""
        if callable(normalizer):
            try:
                normalized = _text(normalizer(token)).upper()
                if normalized:
                    return normalized
            except Exception:
                pass
        return token

    metals = _token_set(
        canonical_metals
        if canonical_metals is not None
        else load_canonical_metals(None),
        normalize,
    )
    cofactors = _token_set(
        canonical_cofactors
        if canonical_cofactors is not None
        else load_canonical_cofactors(None),
        normalize,
    )
    halide_aliases = _mapping(getattr(rules, "halide_aliases", {}))
    halides = _token_set(
        [
            *list(getattr(rules, "halide_resnames", []) or []),
            *list(halide_aliases.keys()),
            *list(halide_aliases.values()),
        ],
        normalize,
    )
    free_cations = _token_set(
        getattr(rules, "meeko_drop_free_ions", []) or [],
        normalize,
    )
    cation_aliases = _mapping(getattr(rules, "cation_aliases", {}))
    free_cations |= {
        normalize(alias)
        for alias, target in cation_aliases.items()
        if normalize(_text(target)) in free_cations
    }
    nucleotide_like = _token_set(
        getattr(rules, "nucleotide_like_resnames", []) or [],
        normalize,
    )
    metals -= halides
    return metals, cofactors, halides, free_cations, nucleotide_like, normalize


def _category(
    identities: set[tuple[str, ...]],
    tokens: set[str],
) -> dict[str, Any]:
    return {"count": len(identities), "tokens": sorted(tokens)}


def audit_prepared_receptor_chemistry(
    path: Path,
    *,
    canonical_metals: Iterable[str] | None = None,
    canonical_cofactors: Iterable[str] | None = None,
    file_role: str = "prepared_receptor",
) -> dict[str, Any]:
    """Extract retained-chemistry evidence without qualifying a receptor.

    Existing but empty, malformed, unsupported, or protein-free files remain
    unresolved. Ambiguous free-ion, halide, and nucleotide-like evidence is
    reported separately and is not converted into an APO/HOLO classification.
    """

    receptor = Path(path).expanduser().resolve()
    if not receptor.is_file():
        raise ReceptorEvidenceError(f"prepared receptor does not exist: {receptor}")
    (
        metals,
        cofactors,
        halides,
        ambiguous_cations,
        nucleotide_like,
        normalize,
    ) = _chemistry_sets(
        canonical_metals=canonical_metals,
        canonical_cofactors=canonical_cofactors,
    )
    lower_name = receptor.name.lower()
    if lower_name.endswith(".pdbqt"):
        file_format = "pdbqt"
        coordinate_parser = pdbqt_line_xyz
        element_parser = pdbqt_line_element
    elif lower_name.endswith((".pdb", ".ent")):
        file_format = "pdb"
        coordinate_parser = line_xyz
        element_parser = line_element
    else:
        file_format = "unsupported"
        coordinate_parser = line_xyz
        element_parser = line_element

    category_identities: dict[str, set[tuple[str, ...]]] = {
        "metal": set(),
        "ambiguous_cation": set(),
        "halide": set(),
        "cofactor": set(),
        "nucleotide_like": set(),
        "other_heterogen": set(),
    }
    category_tokens: dict[str, set[str]] = {
        name: set() for name in category_identities
    }
    atom_record_count = 0
    valid_coordinate_count = 0
    invalid_coordinate_count = 0
    protein_coordinate_count = 0
    nonzero_coordinate_count = 0
    with receptor.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not is_atom_record(line):
                continue
            atom_record_count += 1
            xyz = coordinate_parser(line)
            if xyz is None:
                invalid_coordinate_count += 1
                continue
            valid_coordinate_count += 1
            if any(abs(value) > 0.0 for value in xyz):
                nonzero_coordinate_count += 1
            raw_residue = residue_name(line)
            residue = normalize(raw_residue)
            raw_element = _text(element_parser(line)).upper()
            element = normalize(raw_element)
            record = record_type(line)
            identity = (
                record,
                chain_id(line),
                resseq(line),
                insertion_code(line),
                residue,
                atom_name(line),
                element,
            )
            residue_identity = identity[:5]
            token = residue or element
            category_name = ""
            if residue in halides or (
                record == "HETATM" and element in halides
            ):
                category_name = "halide"
            elif residue in ambiguous_cations or (
                record == "HETATM" and element in ambiguous_cations
            ):
                category_name = "ambiguous_cation"
            elif residue in metals or (record == "HETATM" and element in metals):
                category_name = "metal"
            elif residue in cofactors:
                category_name = "cofactor"
            elif residue in nucleotide_like:
                category_name = "nucleotide_like"
            elif record == "HETATM":
                category_name = "other_heterogen"
            else:
                protein_coordinate_count += 1
            if residue in nucleotide_like:
                category_tokens["nucleotide_like"].add(residue)
                category_identities["nucleotide_like"].add(residue_identity)
            if category_name:
                category_tokens[category_name].add(token)
                if category_name in {"cofactor", "nucleotide_like", "other_heterogen"}:
                    category_identities[category_name].add(residue_identity)
                else:
                    category_identities[category_name].add(identity)

    validation_errors: list[str] = []
    size_bytes = receptor.stat().st_size
    if size_bytes == 0:
        validation_errors.append("empty_file")
    if file_format == "unsupported":
        validation_errors.append("unsupported_receptor_format")
    if atom_record_count == 0:
        validation_errors.append("no_atom_records")
    if invalid_coordinate_count:
        validation_errors.append("malformed_atom_coordinates")
    if valid_coordinate_count == 0:
        validation_errors.append("no_valid_coordinates")
    if protein_coordinate_count == 0:
        validation_errors.append("no_protein_coordinates")
    if valid_coordinate_count and nonzero_coordinate_count == 0:
        validation_errors.append("all_coordinates_at_origin")

    unambiguous_retained = bool(
        category_identities["metal"] or category_identities["cofactor"]
    )
    ambiguous_retained = bool(
        category_identities["ambiguous_cation"]
        or category_identities["halide"]
        or category_identities["nucleotide_like"]
    )
    if validation_errors:
        evidence_status = "unresolved"
        classification = ""
        classification_candidate = None
    elif unambiguous_retained:
        evidence_status = "observed"
        classification = "HOLO"
        classification_candidate = "HOLO"
    elif ambiguous_retained:
        evidence_status = "ambiguous"
        classification = ""
        classification_candidate = "HOLO"
    else:
        evidence_status = "observed"
        classification = "APO"
        classification_candidate = "APO"

    categories = {
        name: _category(category_identities[name], category_tokens[name])
        for name in category_identities
    }
    return {
        "evidence_schema_version": CHEMISTRY_EVIDENCE_SCHEMA_VERSION,
        "classification_method": CHEMISTRY_CLASSIFICATION_METHOD,
        "chemistry_evidence_status": evidence_status,
        "receptor_classification": classification,
        "classification_candidate": classification_candidate,
        "file_role": _text(file_role) or "prepared_receptor",
        "file_format": file_format,
        "prepared_receptor_path": str(receptor),
        "prepared_receptor_sha256": _sha256(receptor),
        "prepared_receptor_size_bytes": size_bytes,
        "atom_record_count": atom_record_count,
        "valid_coordinate_count": valid_coordinate_count,
        "invalid_coordinate_count": invalid_coordinate_count,
        "protein_coordinate_count": protein_coordinate_count,
        "validation_status": "invalid" if validation_errors else "valid",
        "validation_errors": validation_errors,
        "ambiguous_retained_chemistry": ambiguous_retained,
        "retained_categories": categories,
        "retained_metal_atom_count": categories["metal"]["count"],
        "retained_metal_tokens": categories["metal"]["tokens"],
        "retained_cofactor_residue_count": categories["cofactor"]["count"],
        "retained_cofactor_tokens": categories["cofactor"]["tokens"],
        "retained_ambiguous_cation_count": categories["ambiguous_cation"]["count"],
        "retained_ambiguous_cation_tokens": categories["ambiguous_cation"]["tokens"],
        "retained_halide_count": categories["halide"]["count"],
        "retained_halide_tokens": categories["halide"]["tokens"],
        "retained_nucleotide_like_residue_count": categories["nucleotide_like"]["count"],
        "retained_nucleotide_like_tokens": categories["nucleotide_like"]["tokens"],
        "canonical_chemistry_policy_sha256": _policy_sha256(
            metals=metals,
            cofactors=cofactors,
            halides=halides,
            ambiguous_cations=ambiguous_cations,
            nucleotide_like=nucleotide_like,
        ),
        "canonical_metals": sorted(metals),
        "canonical_cofactors": sorted(cofactors),
        "canonical_halides": sorted(halides),
        "ambiguous_cations": sorted(ambiguous_cations),
        "nucleotide_like_resnames": sorted(nucleotide_like),
    }


def _current_run_root(run_manifest_path: Path) -> Path | None:
    for parent in run_manifest_path.resolve().parents:
        if parent.name == "outputs":
            return parent.parent
    return None


def _resolve_recorded_path(
    value: Any,
    *,
    run_manifest_path: Path,
    run_manifest: Mapping[str, Any],
) -> Path | None:
    token = _text(value)
    if not token:
        return None
    recorded = Path(token).expanduser()
    candidates = [
        recorded,
        run_manifest_path.parent / recorded,
    ]
    recorded_root_text = _text(
        _mapping(run_manifest.get("paths")).get("all_dirs")
        or _mapping(run_manifest.get("command")).get("ATLAS_ALL_DIRS")
    )
    recorded_root = Path(recorded_root_text).expanduser()
    current_root = _current_run_root(run_manifest_path)
    if recorded.is_absolute() and recorded_root_text and current_root is not None:
        try:
            candidates.append(current_root / recorded.relative_to(recorded_root))
        except ValueError:
            pass
    for candidate in candidates:
        try:
            resolved = candidate.resolve()
        except OSError:
            continue
        if resolved.is_file():
            return resolved
    return None


def _audit_kind(path: Path) -> str:
    name = path.name.lower()
    if name == "metal_site_audit.json":
        return "metal_site_interaction"
    if name.endswith("retained_hets_audit.json"):
        return "retained_chemistry"
    if name.endswith("_ions.json"):
        return "ion_stage"
    if name == "metal_bound_het_audit.json":
        return "metal_bound_het_interaction"
    if name == "metal_chemistry_audit.json":
        return "metal_chemistry"
    if name == "autodock4zn_plan.json":
        return "autodock4zn_plan"
    if name == "metal_parameterization_plan.json":
        return "metal_parameterization_plan"
    return "other_structured_receptor_audit"


def _safe_path_text(path: Path) -> str:
    try:
        return str(path.resolve())
    except OSError:
        return str(path)


def _safe_sha256(path: Path) -> tuple[str | None, str | None]:
    try:
        return _sha256(path), None
    except OSError as exc:
        return None, str(exc)


def _load_structured_audit(
    path: Path,
) -> tuple[dict[str, Any], Mapping[str, Any] | None]:
    source_sha256, hash_error = _safe_sha256(path)
    record: dict[str, Any] = {
        "audit_kind": _audit_kind(path),
        "source_path": _safe_path_text(path),
        "source_sha256": source_sha256,
    }
    if hash_error:
        record["source_error"] = hash_error
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        record.update(
            {
                "parse_status": "invalid_json",
                "parse_error": str(exc),
            }
        )
        return record, None
    if not isinstance(payload, Mapping):
        record.update(
            {
                "parse_status": "unexpected_shape",
                "payload": payload,
            }
        )
        return record, None
    record.update({"parse_status": "parsed", "payload": payload})
    return record, payload


def _normalized_variant(value: Any) -> str:
    token = _text(value).upper()
    return "" if token in {"", "NONE", "NULL", "LEGACY"} else token


def _normalized_ph(value: Any) -> str:
    token = _text(value)
    return "" if token.lower() in {"", "base", "none", "null"} else token


def _context_docked_root(
    docked_run_root: Path,
    pdb_id: str,
    variant: str,
    ph_label: str,
) -> Path:
    """Mirror the canonical path-router docked context beneath a resolved run root."""

    root = Path(docked_run_root) / pdb_id
    normalized_variant = _normalized_variant(variant)
    if normalized_variant:
        root = root / normalized_variant
    normalized_ph = _normalized_ph(ph_label)
    if normalized_ph:
        root = root / normalized_ph
    return root


def _canonical_receptor_candidate(
    processed_run_root: Path,
    pdb_id: str,
    variant: str,
    ph_label: str,
    *,
    cleaned: bool,
) -> Path:
    """Mirror the path-router receptor layout beneath a resolved run root."""

    root = Path(processed_run_root) / pdb_id
    normalized_variant = _normalized_variant(variant)
    if normalized_variant:
        root = root / normalized_variant
    receptor_root = root / "receptor"
    if cleaned:
        return receptor_root / f"{pdb_id}_cleaned.pdb"
    normalized_ph = _normalized_ph(ph_label)
    if normalized_ph:
        return receptor_root / "ph_ensemble" / f"{pdb_id}_{normalized_ph}.pdbqt"
    return receptor_root / f"{pdb_id}.pdbqt"


def _context_details(context: Mapping[str, Any]) -> Mapping[str, Any]:
    stages = _mapping(context.get("stages"))
    prep = _mapping(stages.get("prep"))
    return _mapping(prep.get("details"))


def _resolve_context_receptors(
    *,
    raw_context: Mapping[str, Any],
    sibling_contexts: Sequence[tuple[str, Mapping[str, Any], str, str, str]],
    pdb_id: str,
    variant: str,
    ph_label: str,
    processed_run_root: Path,
    run_manifest_path: Path,
    run_manifest: Mapping[str, Any],
) -> tuple[Path | None, Path | None, dict[str, Any]]:
    details = _context_details(raw_context)
    cleaned = _resolve_recorded_path(
        details.get("cleaned_pdb"),
        run_manifest_path=run_manifest_path,
        run_manifest=run_manifest,
    )
    docking = _resolve_recorded_path(
        details.get("receptor_pdbqt"),
        run_manifest_path=run_manifest_path,
        run_manifest=run_manifest,
    )
    cleaned_source = "context_manifest" if cleaned is not None else ""
    docking_source = "context_manifest" if docking is not None else ""

    for _, sibling, sibling_pdb, sibling_variant, _ in sibling_contexts:
        if sibling_pdb != pdb_id or sibling_variant != variant:
            continue
        sibling_details = _context_details(sibling)
        if cleaned is None:
            cleaned = _resolve_recorded_path(
                sibling_details.get("cleaned_pdb"),
                run_manifest_path=run_manifest_path,
                run_manifest=run_manifest,
            )
            if cleaned is not None:
                cleaned_source = "sibling_context_manifest"
        if cleaned is not None:
            break

    if cleaned is None:
        candidate = _canonical_receptor_candidate(
            processed_run_root,
            pdb_id,
            variant,
            ph_label,
            cleaned=True,
        )
        if candidate.is_file():
            cleaned = candidate.resolve()
            cleaned_source = "canonical_output_path"
    if docking is None:
        candidate = _canonical_receptor_candidate(
            processed_run_root,
            pdb_id,
            variant,
            ph_label,
            cleaned=False,
        )
        if candidate.is_file():
            docking = candidate.resolve()
            docking_source = "canonical_output_path"

    return cleaned, docking, {
        "cleaned_receptor_resolution": cleaned_source or "unavailable",
        "exact_docking_receptor_resolution": docking_source or "unavailable",
    }


def _audit_context_matches(
    payload: Mapping[str, Any],
    *,
    pdb_id: str,
    variant: str,
    ph_label: str,
) -> tuple[bool, str | None]:
    declared_pdb = _text(payload.get("pdb_id") or payload.get("pdb")).upper()
    if declared_pdb and declared_pdb != pdb_id:
        return False, "declared_pdb_mismatch"
    if "variant" in payload:
        declared_variant = _normalized_variant(payload.get("variant"))
        if declared_variant != _normalized_variant(variant):
            return False, "declared_variant_mismatch"
    if "ph_label" in payload or "ph" in payload:
        declared_ph = _normalized_ph(payload.get("ph_label") or payload.get("ph"))
        if declared_ph != _normalized_ph(ph_label):
            return False, "declared_ph_mismatch"
    return True, None


_RECEPTOR_REFERENCE_KEYS = frozenset(
    {
        "file",
        "target_pdb",
        "receptor_pdb",
        "receptor_pdbqt",
        "receptor_path",
        "prepared_receptor_path",
    }
)


def _audit_receptor_references(payload: Any) -> list[Any]:
    references: list[Any] = []
    if isinstance(payload, Mapping):
        for key, value in payload.items():
            if _text(key).lower() in _RECEPTOR_REFERENCE_KEYS:
                references.append(value)
            references.extend(_audit_receptor_references(value))
    elif isinstance(payload, list):
        for value in payload:
            references.extend(_audit_receptor_references(value))
    return references


def _resolved_reference_matches(
    values: Iterable[Any],
    expected_paths: set[Path],
    *,
    run_manifest_path: Path,
    run_manifest: Mapping[str, Any],
) -> bool:
    for value in values:
        resolved = _resolve_recorded_path(
            value,
            run_manifest_path=run_manifest_path,
            run_manifest=run_manifest,
        )
        if resolved is not None and resolved in expected_paths:
            return True
    return False


def _audit_candidates(
    *,
    processed_run_root: Path,
    docked_context_root: Path,
    pdb_id: str,
) -> list[Path]:
    candidates: set[Path] = set()
    if docked_context_root.is_dir():
        for filename in _AUDIT_FILENAMES:
            path = docked_context_root / filename
            if path.is_file():
                candidates.add(path)
        candidates.update(
            path
            for path in docked_context_root.glob("*retained_hets_audit.json")
            if path.is_file()
        )
        candidates.update(
            path
            for path in (docked_context_root / "audits").glob("*_ions.json")
            if path.is_file()
        )
    pdb_root = processed_run_root / pdb_id
    if pdb_root.is_dir():
        for pattern in _AUDIT_PATTERNS:
            candidates.update(path for path in pdb_root.glob(pattern) if path.is_file())
    return sorted(candidates)


def _bound_structured_audits(
    *,
    processed_run_root: Path,
    docked_context_root: Path,
    pdb_id: str,
    variant: str,
    ph_label: str,
    cleaned_path: Path | None,
    docking_path: Path | None,
    run_manifest_path: Path,
    run_manifest: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    expected_paths = {
        path.resolve()
        for path in (cleaned_path, docking_path)
        if path is not None and path.is_file()
    }
    exact_root = docked_context_root.resolve()
    attached: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = []
    for path in _audit_candidates(
        processed_run_root=processed_run_root,
        docked_context_root=docked_context_root,
        pdb_id=pdb_id,
    ):
        record, payload = _load_structured_audit(path)
        exact_context_path = path.parent.resolve() == exact_root
        if payload is not None:
            context_matches, mismatch_reason = _audit_context_matches(
                payload,
                pdb_id=pdb_id,
                variant=variant,
                ph_label=ph_label,
            )
            if not context_matches:
                rejected.append(
                    {
                        "audit_kind": record["audit_kind"],
                        "source_path": record["source_path"],
                        "source_sha256": record.get("source_sha256"),
                        "rejection_reason": mismatch_reason,
                    }
                )
                continue
            references = _audit_receptor_references(payload)
            receptor_reference_matches = _resolved_reference_matches(
                references,
                expected_paths,
                run_manifest_path=run_manifest_path,
                run_manifest=run_manifest,
            )
            if references and expected_paths and not receptor_reference_matches:
                rejected.append(
                    {
                        "audit_kind": record["audit_kind"],
                        "source_path": record["source_path"],
                        "source_sha256": record.get("source_sha256"),
                        "rejection_reason": "receptor_path_mismatch",
                    }
                )
                continue
        else:
            receptor_reference_matches = False

        if exact_context_path:
            binding_status = (
                "exact_context_path_and_receptor"
                if receptor_reference_matches
                else "exact_context_path"
            )
        elif receptor_reference_matches:
            binding_status = "exact_receptor_reference"
        else:
            rejected.append(
                {
                    "audit_kind": record["audit_kind"],
                    "source_path": record["source_path"],
                    "source_sha256": record.get("source_sha256"),
                    "rejection_reason": "not_bound_to_exact_context",
                }
            )
            continue
        record["binding_status"] = binding_status
        record["bound_context"] = {
            "pdb_id": pdb_id,
            "variant": variant,
            "ph_label": ph_label,
            "cleaned_receptor_sha256": (
                _sha256(cleaned_path) if cleaned_path is not None else None
            ),
            "exact_docking_receptor_sha256": (
                _sha256(docking_path) if docking_path is not None else None
            ),
        }
        attached.append(record)
    return attached, rejected


def _chemistry_signature(evidence: Mapping[str, Any]) -> str:
    categories = _mapping(evidence.get("retained_categories"))
    retained = {
        name: categories.get(name)
        for name in (
            "metal",
            "ambiguous_cation",
            "halide",
            "cofactor",
            "nucleotide_like",
        )
    }
    return json.dumps(retained, sort_keys=True, separators=(",", ":"))


def _combined_chemistry(
    cleaned: Mapping[str, Any] | None,
    docking: Mapping[str, Any] | None,
) -> dict[str, Any]:
    if docking is None:
        return {
            "chemistry_evidence_status": "unresolved",
            "receptor_classification": "",
            "preparation_comparison_status": "exact_docking_receptor_unavailable",
            "chemistry_evidence_error": "exact docking receptor file is unavailable",
        }
    docking_status = _text(docking.get("chemistry_evidence_status"))
    docking_classification = _text(docking.get("receptor_classification"))
    if docking_status != "observed":
        return {
            "chemistry_evidence_status": docking_status or "unresolved",
            "receptor_classification": "",
            "preparation_comparison_status": (
                "exact_docking_receptor_unresolved"
                if docking_status != "ambiguous"
                else "exact_docking_receptor_ambiguous"
            ),
            "chemistry_evidence_error": "; ".join(
                _text(value) for value in docking.get("validation_errors", []) if _text(value)
            )
            or None,
        }
    if cleaned is None:
        return {
            "chemistry_evidence_status": "observed",
            "receptor_classification": docking_classification,
            "preparation_comparison_status": "cleaned_receptor_unavailable",
            "chemistry_evidence_error": None,
        }
    cleaned_status = _text(cleaned.get("chemistry_evidence_status"))
    if cleaned_status != "observed":
        return {
            "chemistry_evidence_status": "unresolved",
            "receptor_classification": "",
            "preparation_comparison_status": (
                "cleaned_receptor_ambiguous"
                if cleaned_status == "ambiguous"
                else "cleaned_receptor_unresolved"
            ),
            "chemistry_evidence_error": "cleaned receptor chemistry could not be compared",
        }
    cleaned_classification = _text(cleaned.get("receptor_classification"))
    if cleaned_classification != docking_classification:
        return {
            "chemistry_evidence_status": "unresolved",
            "receptor_classification": "",
            "preparation_comparison_status": "classification_disagreement",
            "chemistry_evidence_error": (
                f"cleaned receptor is {cleaned_classification}; "
                f"exact docking receptor is {docking_classification}"
            ),
        }
    comparison_status = (
        "matched"
        if _chemistry_signature(cleaned) == _chemistry_signature(docking)
        else "retained_chemistry_disagreement"
    )
    return {
        "chemistry_evidence_status": "observed",
        "receptor_classification": docking_classification,
        "preparation_comparison_status": comparison_status,
        "chemistry_evidence_error": None,
    }


def _context_key_text(key: tuple[str, str, str, str]) -> str:
    run_id, pdb_id, variant, ph_label = key
    return "|".join((run_id, pdb_id, variant, ph_label or "base"))


def _parse_context_key(value: str) -> tuple[str, str, str, str]:
    parts = _text(value).split("|")
    if len(parts) != 4 or any(not part.strip() for part in parts[:3]):
        raise ReceptorEvidenceError(
            "expected context keys must use run_id|PDB|VARIANT|pH_or_base"
        )
    return (
        parts[0].strip(),
        parts[1].strip().upper(),
        parts[2].strip().upper(),
        _normalized_ph(parts[3]),
    )


def _inventory_sha256(keys: Iterable[tuple[str, str, str, str]]) -> str:
    payload = json.dumps(
        [list(key) for key in sorted(keys)],
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_inventory(
    actual_keys: Sequence[tuple[str, str, str, str]],
    *,
    expected_context_count: int | None,
    expected_context_keys: Iterable[str] | None,
    expected_inventory_sha256: str | None,
) -> tuple[list[str] | None, str]:
    if expected_context_count is not None and len(actual_keys) != expected_context_count:
        raise ReceptorEvidenceError(
            f"frozen receptor inventory requires {expected_context_count} contexts; "
            f"observed {len(actual_keys)}"
        )
    expected_key_texts: list[str] | None = None
    if expected_context_keys is not None:
        parsed = [_parse_context_key(value) for value in expected_context_keys]
        if len(set(parsed)) != len(parsed):
            raise ReceptorEvidenceError("frozen receptor inventory contains duplicate keys")
        if set(parsed) != set(actual_keys):
            missing = sorted(set(parsed) - set(actual_keys))
            unexpected = sorted(set(actual_keys) - set(parsed))
            raise ReceptorEvidenceError(
                "frozen receptor inventory key mismatch: "
                f"missing={[ _context_key_text(key) for key in missing ]}; "
                f"unexpected={[ _context_key_text(key) for key in unexpected ]}"
            )
        expected_key_texts = [_context_key_text(key) for key in sorted(parsed)]
    inventory_sha256 = _inventory_sha256(actual_keys)
    if expected_inventory_sha256 is not None:
        expected_hash = _text(expected_inventory_sha256).lower()
        if expected_hash.startswith("sha256:"):
            expected_hash = expected_hash.removeprefix("sha256:")
        if len(expected_hash) != 64 or any(
            char not in "0123456789abcdef" for char in expected_hash
        ):
            raise ReceptorEvidenceError(
                "expected inventory SHA-256 must be 64 hexadecimal characters"
            )
        if inventory_sha256 != expected_hash:
            raise ReceptorEvidenceError(
                "frozen receptor inventory SHA-256 mismatch: "
                f"expected {expected_hash}; observed {inventory_sha256}"
            )
    return expected_key_texts, inventory_sha256


def build_receptor_evidence_records(
    manifest_path: Path,
    repo_root: Path,
    *,
    expected_context_count: int | None = None,
    expected_context_keys: Iterable[str] | None = None,
    expected_inventory_sha256: str | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build exact-context receptor annotations without qualifying receptors."""

    release_manifest_path = Path(manifest_path).expanduser().resolve()
    manifest = load_release_manifest(release_manifest_path)
    records: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str, str]] = set()
    canonical_metals = load_canonical_metals(None)
    canonical_cofactors = load_canonical_cofactors(None)
    for run_entry in iter_run_entries(manifest):
        run_id = _text(run_entry.get("run_id"))
        if not run_id:
            raise ReceptorEvidenceError("release run entry requires a non-empty run_id")
        run_paths = _run_paths(Path(repo_root), run_entry)
        run_manifest_path = run_paths["manifest"].expanduser().resolve()
        if not run_manifest_path.is_file():
            raise ReceptorEvidenceError(
                f"run manifest does not exist for {run_id}: {run_manifest_path}"
            )
        run_manifest = _load_mapping(run_manifest_path)
        proteins = _mapping(run_manifest.get("proteins"))
        contexts: list[tuple[str, Mapping[str, Any], str, str, str]] = []
        for raw_key, raw_context in proteins.items():
            if not isinstance(raw_context, Mapping):
                continue
            pdb_id, variant, ph_label = _context_values(str(raw_key), raw_context)
            if not pdb_id:
                raise ReceptorEvidenceError(
                    f"run {run_id} contains a receptor context without a PDB ID"
                )
            contexts.append((str(raw_key), raw_context, pdb_id, variant, ph_label))
        processed_run_root = run_output_dir(
            Path(repo_root), "processed_pdbs", run_id
        ).resolve()
        docked_run_root = run_paths["docked"].expanduser().resolve()
        run_manifest_sha256 = _sha256(run_manifest_path)
        for raw_key, raw_context, pdb_id, variant, ph_label in contexts:
            key = (run_id, pdb_id, variant, ph_label)
            if key in seen:
                raise ReceptorEvidenceError(f"duplicate receptor context: {key!r}")
            seen.add(key)
            cleaned_path, docking_path, path_resolution = _resolve_context_receptors(
                raw_context=raw_context,
                sibling_contexts=contexts,
                pdb_id=pdb_id,
                variant=variant,
                ph_label=ph_label,
                processed_run_root=processed_run_root,
                run_manifest_path=run_manifest_path,
                run_manifest=run_manifest,
            )
            cleaned_chemistry: dict[str, Any] | None = None
            docking_chemistry: dict[str, Any] | None = None
            if cleaned_path is not None:
                cleaned_chemistry = audit_prepared_receptor_chemistry(
                    cleaned_path,
                    canonical_metals=canonical_metals,
                    canonical_cofactors=canonical_cofactors,
                    file_role="cleaned_prepared_receptor",
                )
            if docking_path is not None:
                docking_chemistry = audit_prepared_receptor_chemistry(
                    docking_path,
                    canonical_metals=canonical_metals,
                    canonical_cofactors=canonical_cofactors,
                    file_role="exact_docking_receptor",
                )
            combined = _combined_chemistry(cleaned_chemistry, docking_chemistry)
            observed = _text(combined["receptor_classification"])
            chemistry_status = _text(combined["chemistry_evidence_status"])
            conflict = bool(
                observed
                and variant in CONTROLLED_RECEPTOR_CLASSIFICATIONS
                and observed != variant
            )
            docked_context_root = _context_docked_root(
                docked_run_root,
                pdb_id,
                variant,
                ph_label,
            )
            audits, audit_rejections = _bound_structured_audits(
                processed_run_root=processed_run_root,
                docked_context_root=docked_context_root,
                pdb_id=pdb_id,
                variant=variant,
                ph_label=ph_label,
                cleaned_path=cleaned_path,
                docking_path=docking_path,
                run_manifest_path=run_manifest_path,
                run_manifest=run_manifest,
            )
            primary_chemistry = docking_chemistry or cleaned_chemistry or {}
            record = {
                "run_id": run_id,
                "pdb_id": pdb_id,
                "variant": variant,
                "ph_label": ph_label,
                "receptor_classification": observed,
                "classification_method": CHEMISTRY_CLASSIFICATION_METHOD,
                "chemistry_evidence_status": chemistry_status,
                "prepared_receptor_sha256": (
                    docking_chemistry.get("prepared_receptor_sha256")
                    if docking_chemistry
                    else None
                ),
                "canonical_chemistry_policy_sha256": primary_chemistry.get(
                    "canonical_chemistry_policy_sha256"
                ),
                "retained_metal_atom_count": (
                    docking_chemistry.get("retained_metal_atom_count")
                    if docking_chemistry
                    else None
                ),
                "retained_cofactor_residue_count": (
                    docking_chemistry.get("retained_cofactor_residue_count")
                    if docking_chemistry
                    else None
                ),
                "requested_observed_conflict": conflict,
                "qualification_status": "pending_receptor_quality_review",
                "qualification_reason": (
                    "evidence only; receptor quality and native redocking "
                    "qualification have not been applied"
                ),
                "provenance": {
                    "evidence_schema_version": CHEMISTRY_EVIDENCE_SCHEMA_VERSION,
                    "run_manifest_path": str(run_manifest_path),
                    "run_manifest_sha256": run_manifest_sha256,
                    "run_manifest_context_key": raw_key,
                    "requested_variant": variant,
                    "observed_receptor_classification": observed or None,
                    "requested_observed_conflict": conflict,
                    "classification_method": CHEMISTRY_CLASSIFICATION_METHOD,
                    "chemistry_evidence_status": chemistry_status,
                    "chemistry_evidence_error": combined.get(
                        "chemistry_evidence_error"
                    ),
                    "preparation_comparison_status": combined.get(
                        "preparation_comparison_status"
                    ),
                    "prepared_receptor_sha256": (
                        docking_chemistry.get("prepared_receptor_sha256")
                        if docking_chemistry
                        else None
                    ),
                    "canonical_chemistry_policy_sha256": primary_chemistry.get(
                        "canonical_chemistry_policy_sha256"
                    ),
                    "retained_metal_atom_count": (
                        docking_chemistry.get("retained_metal_atom_count")
                        if docking_chemistry
                        else None
                    ),
                    "retained_cofactor_residue_count": (
                        docking_chemistry.get("retained_cofactor_residue_count")
                        if docking_chemistry
                        else None
                    ),
                    **path_resolution,
                    "cleaned_receptor_chemistry": cleaned_chemistry,
                    "exact_docking_receptor_chemistry": docking_chemistry,
                    "prepared_receptor_chemistry": docking_chemistry,
                    "structured_receptor_audits": audits,
                    "structured_receptor_audit_rejections": audit_rejections,
                },
            }
            records.append(record)
    records.sort(
        key=lambda row: (
            row["run_id"],
            row["pdb_id"],
            row["variant"],
            row["ph_label"],
        )
    )
    actual_keys = [
        (
            _text(row["run_id"]),
            _text(row["pdb_id"]).upper(),
            _text(row["variant"]).upper(),
            _normalized_ph(row["ph_label"]),
        )
        for row in records
    ]
    expected_key_texts, inventory_sha256 = _validate_inventory(
        actual_keys,
        expected_context_count=expected_context_count,
        expected_context_keys=expected_context_keys,
        expected_inventory_sha256=expected_inventory_sha256,
    )
    observed_count = sum(bool(row["receptor_classification"]) for row in records)
    conflict_count = sum(bool(row["requested_observed_conflict"]) for row in records)
    ambiguous_count = sum(
        row["chemistry_evidence_status"] == "ambiguous" for row in records
    )
    audit_rows = [
        audit
        for row in records
        for audit in _mapping(row["provenance"]).get(
            "structured_receptor_audits", []
        )
        if isinstance(audit, Mapping)
    ]
    audit_kinds = [_text(audit.get("audit_kind")) for audit in audit_rows]
    audit_rejection_count = sum(
        len(
            _mapping(row["provenance"]).get(
                "structured_receptor_audit_rejections", []
            )
        )
        for row in records
    )
    summary = {
        "evidence_schema_version": CHEMISTRY_EVIDENCE_SCHEMA_VERSION,
        "release_id": manifest["release_id"],
        "release_manifest_path": str(release_manifest_path),
        "release_manifest_sha256": _sha256(release_manifest_path),
        "expected_context_count": expected_context_count,
        "expected_context_keys": expected_key_texts,
        "expected_inventory_sha256": expected_inventory_sha256,
        "context_count": len(records),
        "context_keys": [_context_key_text(key) for key in actual_keys],
        "context_inventory_sha256": inventory_sha256,
        "chemistry_observed_count": observed_count,
        "chemistry_ambiguous_count": ambiguous_count,
        "chemistry_unresolved_count": len(records) - observed_count,
        "requested_observed_conflict_count": conflict_count,
        "structured_receptor_audit_count": len(audit_kinds),
        "structured_receptor_audit_parsed_count": sum(
            audit.get("parse_status") == "parsed" for audit in audit_rows
        ),
        "structured_receptor_audit_rejection_count": audit_rejection_count,
        "structured_receptor_audit_kind_counts": {
            kind: audit_kinds.count(kind) for kind in sorted(set(audit_kinds))
        },
        "classification_counts": {
            token: sum(row["receptor_classification"] == token for row in records)
            for token in sorted(CONTROLLED_RECEPTOR_CLASSIFICATIONS)
        },
    }
    return records, summary


def write_receptor_evidence(
    manifest_path: Path,
    repo_root: Path,
    output_path: Path,
    *,
    expected_context_count: int | None = None,
    expected_context_keys: Iterable[str] | None = None,
    expected_inventory_sha256: str | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Write hashed receptor annotation evidence and its compact summary."""

    output = Path(output_path).expanduser().resolve()
    summary_path = output.with_suffix(".summary.json")
    if not overwrite:
        for path in (output, summary_path):
            if path.exists():
                raise ReceptorEvidenceError(f"output already exists: {path}")
    records, summary = build_receptor_evidence_records(
        manifest_path,
        repo_root,
        expected_context_count=expected_context_count,
        expected_context_keys=expected_context_keys,
        expected_inventory_sha256=expected_inventory_sha256,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        yaml.safe_dump({"records": records}, sort_keys=False), encoding="utf-8"
    )
    summary = {
        **summary,
        "annotation_path": str(output),
        "annotation_sha256": _sha256(output),
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


__all__ = [
    "CHEMISTRY_CLASSIFICATION_METHOD",
    "CHEMISTRY_EVIDENCE_SCHEMA_VERSION",
    "CONTROLLED_RECEPTOR_CLASSIFICATIONS",
    "ReceptorEvidenceError",
    "audit_prepared_receptor_chemistry",
    "build_receptor_evidence_records",
    "write_receptor_evidence",
]
