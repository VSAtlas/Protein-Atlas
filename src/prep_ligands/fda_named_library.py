"""Build a manifest-backed, human-readable FDA canonical-parent library.

The builder never mutates legacy PDBQTs. Compatible inputs are copied byte for
byte, missing canonical parents are prepared with the shared Meeko backend, and
chemistry unsuitable for conventional small-molecule docking is reported rather
than silently coerced into a PDBQT.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from rdkit import Chem, rdBase

from prep_ligands.fda_identity_audit import (
    ACTIVE,
    ADDITIVE,
    AMBIGUOUS_MIXTURE,
    COMBINATION,
    METAL_COMPLEX,
    SALT_OR_SOLVATE,
    STANDARDIZATION_VERSION,
    _classify_substance,
    _parent_mol,
    _read_source_sdf,
    _safe_mol,
    _structure_from_mol,
)
from prep_ligands.library_index import LibraryIndex
from prep_ligands.prep_ligands_meeko import prepare_mol_to_pdbqt_with_meeko


SCHEMA_VERSION = 1
SLUG_VERSION = "nfkd-ascii-lower-underscore-v1"
PREP_POLICY = "audit_standardized_parent_charge_stereo_meeko_etkdg_seed_0xA7A52"
_UNSAFE_CLASSES = {ADDITIVE, COMBINATION, METAL_COMPLEX, AMBIGUOUS_MIXTURE}
_REUSABLE_RELATIONS = {"approved_full_form_docked", "approved_parent_docked"}
_SUPPORTED_ATOMIC_NUMBERS = {
    5,  # B
    6,  # C
    7,  # N
    8,  # O
    9,  # F
    15,  # P
    16,  # S
    17,  # Cl
    35,  # Br
    53,  # I
}

_MANIFEST_FIELDS = [
    "canonical_parent_id",
    "preferred_name",
    "filename",
    "slug_version",
    "parent_smiles",
    "parent_inchikey",
    "drugcentral_ids",
    "source_record_indices",
    "approved_full_form_names",
    "approved_full_form_exact_inchikeys",
    "approved_source_classes",
    "approved_form_relation",
    "materialization_action",
    "materialization_status",
    "status_reason_codes",
    "output_pdbqt_path",
    "output_pdbqt_sha256",
    "output_pdbqt_bytes",
    "sidecar_path",
    "sidecar_sha256",
    "legacy_selected_rdk_id",
    "legacy_selected_mapping_row",
    "legacy_selected_path",
    "legacy_selected_sha256",
    "legacy_candidate_rdk_ids",
    "legacy_candidate_paths",
    "legacy_candidate_sha256s",
    "legacy_selection_policy",
    "trigger_rdk_ids",
    "prep_policy",
    "standardization_version",
]


def build_named_fda_library(
    *,
    repaired_mapping_csv: Path,
    redock_delta_csv: Path,
    source_sdf: Path,
    approval_manifest: Path,
    output_dir: Path,
    workers: int = 8,
    plan_only: bool = False,
) -> dict[str, Any]:
    """Materialize the repaired FDA parent library and deterministic manifests."""

    repaired_mapping_csv = repaired_mapping_csv.expanduser().resolve()
    redock_delta_csv = redock_delta_csv.expanduser().resolve()
    source_sdf = source_sdf.expanduser().resolve()
    approval_manifest = approval_manifest.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    _require_input_files(
        repaired_mapping_csv, redock_delta_csv, source_sdf, approval_manifest
    )
    source_index = _read_source_sdf(source_sdf, approval_manifest)
    if not source_index.approval_provenance_verified:
        raise ValueError(
            "FDA source SDF is not backed by a valid exact-filter approval manifest: "
            f"{source_index.approval_provenance_reason}"
        )

    source_groups, unsafe_records = _source_parent_groups(source_index)
    parent_mols = _source_parent_mols(source_sdf, source_groups)
    repaired = _read_csv(repaired_mapping_csv)
    delta = _read_csv(redock_delta_csv)
    legacy_candidates = _legacy_candidates(repaired, set(source_groups))
    covered = set(legacy_candidates)
    delta_keys = _validate_delta(delta, source_groups)
    _validate_partition(set(source_groups), covered, delta_keys)

    input_hashes = {
        "repaired_mapping_csv": _sha256(repaired_mapping_csv),
        "redock_delta_csv": _sha256(redock_delta_csv),
        "source_sdf": _sha256(source_sdf),
        "approval_manifest": _sha256(approval_manifest),
    }
    input_paths = {
        "repaired_mapping_csv": str(repaired_mapping_csv),
        "redock_delta_csv": str(redock_delta_csv),
        "source_sdf": str(source_sdf),
        "approval_manifest": str(approval_manifest),
    }
    _validate_output_scope(output_dir, input_paths.values(), repaired)
    names = _assign_names(source_groups)
    triggers = _trigger_rdk_ids(repaired)
    rows = _build_plan_rows(
        source_groups,
        names,
        legacy_candidates,
        delta_keys,
        triggers,
        output_dir,
    )
    _annotate_delta_suitability(rows, parent_mols, source_groups)
    if plan_only:
        return _plan_summary(
            rows, unsafe_records, input_paths, input_hashes, output_dir
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    _validate_resume_inputs(output_dir / "fda_named_library_summary.json", input_hashes)
    parent_sdf = output_dir / "fda_canonical_parent_source.sdf"
    _write_parent_sdf(parent_sdf, rows, parent_mols)
    prep_context = {
        "input_paths": input_paths,
        "input_sha256": input_hashes,
        "source_parent_sdf": str(parent_sdf),
        "source_parent_sdf_sha256": _sha256(parent_sdf),
        "rdkit_version": rdBase.rdkitVersion,
        "meeko_version": _package_version("meeko"),
    }
    _materialize_reused(rows, output_dir, prep_context)
    _materialize_delta(
        rows,
        output_dir,
        prep_context,
        parent_mols,
        workers=max(1, min(32, int(workers))),
    )

    rows.sort(key=lambda row: str(row["parent_inchikey"]))
    manifest_path = output_dir / "fda_named_library_manifest.csv"
    quarantine_path = output_dir / "fda_named_library_quarantine.csv"
    _write_csv(manifest_path, rows, _MANIFEST_FIELDS)
    _write_csv(
        quarantine_path,
        [row for row in rows if row["materialization_status"] != "ready"],
        _MANIFEST_FIELDS,
    )
    ready_paths = [
        Path(str(row["output_pdbqt_path"]))
        for row in rows
        if row["materialization_status"] == "ready"
    ]
    _validate_ready_outputs(rows, ready_paths, output_dir)
    library_index = LibraryIndex().write_manifest_for_root(
        output_dir,
        [path.relative_to(output_dir) for path in ready_paths],
    )
    summary = _summary(
        rows=rows,
        unsafe_records=unsafe_records,
        input_paths=input_paths,
        input_hashes=input_hashes,
        prep_context=prep_context,
        manifest_path=manifest_path,
        quarantine_path=quarantine_path,
        library_index=library_index,
        output_dir=output_dir,
    )
    _write_json(output_dir / "fda_named_library_summary.json", summary)
    return summary


def _require_input_files(*paths: Path) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise ValueError("missing FDA named-library input(s): " + ", ".join(missing))


def _validate_output_scope(
    output_dir: Path,
    input_paths: Iterable[str],
    repaired: Sequence[Mapping[str, str]],
) -> None:
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"FDA named-library output is not a directory: {output_dir}")
    _validate_input_output_scope(output_dir, input_paths)
    _validate_legacy_output_scope(output_dir, repaired)


def _validate_input_output_scope(
    output_dir: Path, input_paths: Iterable[str]
) -> None:
    resolved_inputs = {Path(path).resolve() for path in input_paths}
    if output_dir in resolved_inputs or any(
        output_dir in path.parents for path in resolved_inputs
    ):
        raise ValueError("FDA named-library output contains or replaces an input")


def _validate_legacy_output_scope(
    output_dir: Path, repaired: Sequence[Mapping[str, str]]
) -> None:
    legacy_roots = {
        Path(row["selected_pdbqt_path"]).expanduser().resolve().parent
        for row in repaired
        if row.get("selected_pdbqt_path")
    }
    for root in legacy_roots:
        if (
            output_dir == root
            or output_dir in root.parents
            or root in output_dir.parents
        ):
            raise ValueError(
                "FDA named-library output must be separate from the legacy library: "
                f"{root}"
            )


def _source_parent_groups(source_index: Any) -> tuple[dict[str, list[Any]], list[Any]]:
    groups: dict[str, list[Any]] = defaultdict(list)
    unsafe: list[Any] = []
    for record in source_index.by_index.values():
        source_class = _classify_substance(record.structure, record.name, False)
        key = record.structure.parent_inchikey
        if source_class in _UNSAFE_CLASSES or not key:
            unsafe.append(record)
            continue
        groups[key].append(record)
    return dict(sorted(groups.items())), unsafe


def _source_parent_mols(
    source_sdf: Path, source_groups: Mapping[str, Sequence[Any]]
) -> dict[str, Chem.Mol]:
    representative_indices = {
        int(_representative_record(records).index): key
        for key, records in source_groups.items()
    }
    parents: dict[str, Chem.Mol] = {}
    supplier = Chem.SDMolSupplier(str(source_sdf), removeHs=False, sanitize=False)
    for record_index, raw_mol in enumerate(supplier, start=1):
        key = representative_indices.get(record_index)
        if key is None:
            continue
        safe_mol = _safe_mol(raw_mol)
        if safe_mol is None:
            raise ValueError(f"FDA source SDF record is unreadable: {record_index}")
        parent = _parent_mol(safe_mol)
        evidence = _structure_from_mol(parent, f"fda_named_parent:{record_index}")
        if evidence.parent_inchikey != key:
            raise ValueError(
                "FDA source parent does not reproduce audited parent identity: "
                f"record={record_index} expected={key} actual={evidence.parent_inchikey}"
            )
        parents[key] = parent
    missing = set(source_groups) - set(parents)
    if missing:
        raise ValueError(
            f"FDA source parent molecules are missing: count={len(missing)}"
        )
    return parents


def _legacy_candidates(
    repaired: Sequence[Mapping[str, str]], safe_parent_keys: set[str]
) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in repaired:
        key = row.get("approved_parent_inchikey", "")
        if key not in safe_parent_keys:
            continue
        if row.get("approved_docked_relation") not in _REUSABLE_RELATIONS:
            continue
        path = Path(row.get("selected_pdbqt_path", "")).expanduser()
        stated_hash = row.get("selected_pdbqt_sha256", "")
        if not path.is_file() or path.stat().st_size <= 0 or not stated_hash:
            continue
        actual_hash = _sha256(path)
        if actual_hash != stated_hash:
            raise ValueError(f"legacy PDBQT checksum changed: {path}")
        grouped[key].append(dict(row))
    return dict(grouped)


def _validate_delta(
    delta: Sequence[Mapping[str, str]], source_groups: Mapping[str, Sequence[Any]]
) -> set[str]:
    keys = [row.get("desired_parent_inchikey", "") for row in delta]
    if not all(keys) or len(keys) != len(set(keys)):
        raise ValueError("FDA redock delta parent keys are empty or duplicated")
    extra = set(keys) - set(source_groups)
    if extra:
        raise ValueError(f"FDA redock delta contains non-approved parent keys: {sorted(extra)}")
    for row in delta:
        key = row["desired_parent_inchikey"]
        expected = source_groups[key][0].structure.parent_smiles
        if row.get("desired_parent_smiles") != expected:
            raise ValueError(f"FDA redock delta parent SMILES mismatch: {key}")
    return set(keys)


def _validate_partition(safe: set[str], covered: set[str], delta: set[str]) -> None:
    if covered & delta:
        raise ValueError("covered and redock FDA parent sets overlap")
    if covered | delta != safe:
        missing = safe - covered - delta
        extra = covered | delta - safe
        raise ValueError(
            f"FDA parent partition is incomplete: missing={len(missing)} extra={len(extra)}"
        )


def _assign_names(source_groups: Mapping[str, Sequence[Any]]) -> dict[str, str]:
    bases: dict[str, str] = {}
    for key, records in source_groups.items():
        representative = _representative_record(records)
        base = _slug(representative.name)
        if not base:
            base = f"drugcentral_{_record_id(representative) or representative.index}"
        bases[key] = base
    collisions = Counter(value.casefold() for value in bases.values())
    assigned: dict[str, str] = {}
    used: set[str] = set()
    for key in sorted(bases):
        records = source_groups[key]
        name = bases[key]
        if collisions[name.casefold()] > 1:
            record_id = _record_id(_representative_record(records))
            name = f"{name}__dc{record_id or _representative_record(records).index}"
        if name.casefold() in used:
            name = f"{name}__{hashlib.sha256(key.encode()).hexdigest()[:10]}"
        if name.casefold() in used:
            raise ValueError(f"FDA output filename collision: {name}")
        used.add(name.casefold())
        assigned[key] = name
    return assigned


def _representative_record(records: Sequence[Any]) -> Any:
    return sorted(
        records,
        key=lambda record: (
            _classify_substance(record.structure, record.name, False) != ACTIVE,
            record.structure.fragment_count != 1,
            len(str(record.name)),
            int(record.index),
        ),
    )[0]


def _build_plan_rows(
    source_groups: Mapping[str, Sequence[Any]],
    names: Mapping[str, str],
    legacy_candidates: Mapping[str, Sequence[Mapping[str, str]]],
    delta_keys: set[str],
    triggers: Mapping[str, Sequence[str]],
    output_dir: Path,
) -> list[dict[str, Any]]:
    return [
        _build_plan_row(
            key=key,
            records=records,
            name=names[key],
            candidates=list(legacy_candidates.get(key, ())),
            in_delta=key in delta_keys,
            trigger_rdk_ids=triggers.get(key, ()),
            output_dir=output_dir,
        )
        for key, records in source_groups.items()
    ]


def _build_plan_row(
    *,
    key: str,
    records: Sequence[Any],
    name: str,
    candidates: Sequence[Mapping[str, str]],
    in_delta: bool,
    trigger_rdk_ids: Sequence[str],
    output_dir: Path,
) -> dict[str, Any]:
    representative = _representative_record(records)
    selected, action = _plan_action(key, candidates, in_delta)
    filename = f"{name}.pdbqt"
    row = {
        "canonical_parent_id": f"fda_parent_{key}",
        "preferred_name": representative.name,
        "filename": filename,
        "slug_version": SLUG_VERSION,
        "parent_smiles": representative.structure.parent_smiles,
        "parent_inchikey": key,
        "approved_form_relation": _form_relation(records),
        "materialization_action": action,
        "materialization_status": "planned",
        "status_reason_codes": "",
        "output_pdbqt_path": str(output_dir / filename),
        "output_pdbqt_sha256": "",
        "output_pdbqt_bytes": "",
        "sidecar_path": str(
            (output_dir / filename).with_suffix(".ligprep_source.json")
        ),
        "sidecar_sha256": "",
        "trigger_rdk_ids": _joined(trigger_rdk_ids),
        "prep_policy": PREP_POLICY,
        "standardization_version": STANDARDIZATION_VERSION,
    }
    row.update(_record_lineage_fields(records))
    row.update(_candidate_lineage_fields(candidates))
    row.update(_selected_legacy_fields(selected))
    return row


def _plan_action(
    key: str,
    candidates: Sequence[Mapping[str, str]],
    in_delta: bool,
) -> tuple[Mapping[str, str] | None, str]:
    selected = _select_legacy_candidate(candidates) if candidates else None
    action = "legacy_byte_copy" if selected else "prepare_parent"
    if in_delta == bool(selected):
        raise ValueError(f"FDA parent action disagrees with repair partition: {key}")
    return selected, action


def _record_lineage_fields(records: Sequence[Any]) -> dict[str, str]:
    return {
        "drugcentral_ids": _joined(_record_id(record) for record in records),
        "source_record_indices": _joined(str(record.index) for record in records),
        "approved_full_form_names": _joined(record.name for record in records),
        "approved_full_form_exact_inchikeys": _joined(
            record.structure.exact_inchikey for record in records
        ),
        "approved_source_classes": _joined(
            _classify_substance(record.structure, record.name, False)
            for record in records
        ),
    }


def _candidate_lineage_fields(
    candidates: Sequence[Mapping[str, str]],
) -> dict[str, str]:
    return {
        "legacy_candidate_rdk_ids": _joined(
            candidate.get("rdk_id", "") for candidate in candidates
        ),
        "legacy_candidate_paths": _joined(
            candidate.get("selected_pdbqt_path", "") for candidate in candidates
        ),
        "legacy_candidate_sha256s": _joined(
            candidate.get("selected_pdbqt_sha256", "") for candidate in candidates
        ),
    }


def _selected_legacy_fields(
    selected: Mapping[str, str] | None,
) -> dict[str, str]:
    if selected is None:
        return {
            "legacy_selected_rdk_id": "",
            "legacy_selected_mapping_row": "",
            "legacy_selected_path": "",
            "legacy_selected_sha256": "",
            "legacy_selection_policy": "",
        }
    return {
        "legacy_selected_rdk_id": selected.get("rdk_id", ""),
        "legacy_selected_mapping_row": selected.get("mapping_row_number", ""),
        "legacy_selected_path": selected.get("selected_pdbqt_path", ""),
        "legacy_selected_sha256": selected.get("selected_pdbqt_sha256", ""),
        "legacy_selection_policy": "prepared_state_exact_then_mapping_row_rdk_path_v1",
    }


def _select_legacy_candidate(
    candidates: Sequence[Mapping[str, str]],
) -> Mapping[str, str]:
    state_rank = {
        "high_fidelity_exact": 0,
        "coordinate_connectivity_exact": 1,
        "high_fidelity_parent": 2,
        "coordinate_connectivity_parent": 3,
    }
    return sorted(
        candidates,
        key=lambda row: (
            state_rank.get(row.get("prepared_connectivity_state", ""), 9),
            int(row.get("mapping_row_number", "0") or 0),
            row.get("rdk_id", ""),
            row.get("selected_pdbqt_path", ""),
        ),
    )[0]


def _materialize_reused(
    rows: Sequence[dict[str, Any]],
    output_dir: Path,
    prep_context: Mapping[str, Any],
) -> None:
    for row in rows:
        if row["materialization_action"] != "legacy_byte_copy":
            continue
        src = Path(str(row["legacy_selected_path"]))
        dst = Path(str(row["output_pdbqt_path"]))
        expected = str(row["legacy_selected_sha256"])
        if _ready_existing(
            dst,
            row["parent_inchikey"],
            expected,
            prep_context,
        ):
            _finish_ready_row(row, dst)
            continue
        if dst.exists():
            raise ValueError(f"refusing to overwrite changed FDA named PDBQT: {dst}")
        shutil.copy2(src, dst)
        if _sha256(dst) != expected:
            raise RuntimeError(f"byte-preserving FDA PDBQT copy failed: {src}")
        _write_sidecar(row, dst, prep_context, writer="legacy_byte_copy")
        _finish_ready_row(row, dst)


def _materialize_delta(
    rows: Sequence[dict[str, Any]],
    output_dir: Path,
    prep_context: Mapping[str, Any],
    parent_mols: Mapping[str, Chem.Mol],
    *,
    workers: int,
) -> None:
    tasks = _pending_delta_rows(rows, prep_context)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                _prepare_delta_row, row, parent_mols, output_dir
            )
            for row in tasks
        ]
        for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            row, result = future.result()
            _record_delta_result(row, result, output_dir, prep_context)
            if index % 100 == 0:
                print(f"[fda-named-library] prepared {index}/{len(tasks)}")


def _pending_delta_rows(
    rows: Sequence[dict[str, Any]],
    prep_context: Mapping[str, Any],
) -> list[dict[str, Any]]:
    tasks = []
    for row in rows:
        if row["materialization_action"] != "prepare_parent":
            continue
        if row["materialization_status"] == "quarantined_non_dockable":
            continue
        dst = Path(str(row["output_pdbqt_path"]))
        if _ready_existing(
            dst,
            row["parent_inchikey"],
            "",
            prep_context,
        ):
            _finish_ready_row(row, dst)
            continue
        if dst.exists():
            raise ValueError(f"refusing to overwrite changed FDA named PDBQT: {dst}")
        tasks.append(row)
    return tasks


def _prepare_delta_row(
    row: dict[str, Any],
    parent_mols: Mapping[str, Chem.Mol],
    output_dir: Path,
) -> tuple[dict[str, Any], Any]:
    mol = Chem.Mol(parent_mols[str(row["parent_inchikey"])])
    result = prepare_mol_to_pdbqt_with_meeko(
        mol,
        Path(str(row["output_pdbqt_path"])),
        ligand_name=str(row["preferred_name"]),
        log_dir=output_dir,
    )
    return row, result


def _record_delta_result(
    row: dict[str, Any],
    result: Any,
    output_dir: Path,
    prep_context: Mapping[str, Any],
) -> None:
    dst = Path(str(row["output_pdbqt_path"]))
    if result is None or not result.ok or not dst.is_file():
        row["materialization_status"] = "prep_failed"
        row["status_reason_codes"] = (
            result.status if result is not None else "parent_source_mol_unavailable"
        )
        return
    _write_sidecar(row, dst, prep_context, writer="meeko")
    _finish_ready_row(row, dst)


def _annotate_delta_suitability(
    rows: Sequence[dict[str, Any]],
    parent_mols: Mapping[str, Chem.Mol],
    source_groups: Mapping[str, Sequence[Any]],
) -> None:
    planned = 0
    quarantined = 0
    for row in rows:
        if row["materialization_action"] != "prepare_parent":
            continue
        planned += 1
        key = str(row["parent_inchikey"])
        reasons = _dockability_reasons(parent_mols[key], key, source_groups[key])
        if not reasons:
            continue
        quarantined += 1
        row["materialization_status"] = "quarantined_non_dockable"
        row["status_reason_codes"] = ";".join(reasons)
    if planned == 599 and quarantined != 13:
        raise RuntimeError(
            "FDA delta suitability policy drifted: "
            f"planned={planned} expected_quarantined=13 actual={quarantined}"
        )


def _dockability_reasons(
    mol: Chem.Mol, parent_inchikey: str, source_records: Sequence[Any]
) -> list[str]:
    if mol.GetNumAtoms() == 0:
        return ["degenerate_parent_graph"]
    evidence = _structure_from_mol(mol, "named_fda_parent")
    if evidence.parent_inchikey != parent_inchikey:
        return ["parent_inchikey_recanonicalization_mismatch"]
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if _heavy_graph_is_degenerate(mol, len(heavy_atoms)):
        return _degenerate_parent_reasons(parent_inchikey, source_records)
    if _is_tiny_carbon_free_or_radical(heavy_atoms):
        return ["tiny_carbon_free_inorganic_or_radical_parent"]
    unsupported = _unsupported_parent_elements(heavy_atoms)
    if unsupported:
        return ["unsupported_element_capability:" + ",".join(unsupported)]
    return []


def _heavy_graph_is_degenerate(mol: Chem.Mol, heavy_atom_count: int) -> bool:
    if heavy_atom_count <= 1:
        return True
    return not any(
        bond.GetBeginAtom().GetAtomicNum() > 1
        and bond.GetEndAtom().GetAtomicNum() > 1
        for bond in mol.GetBonds()
    )


def _is_tiny_carbon_free_or_radical(heavy_atoms: Sequence[Any]) -> bool:
    if len(heavy_atoms) > 2:
        return False
    carbon_free = all(atom.GetAtomicNum() != 6 for atom in heavy_atoms)
    has_radical = any(atom.GetNumRadicalElectrons() > 0 for atom in heavy_atoms)
    return carbon_free or has_radical


def _degenerate_parent_reasons(
    parent_inchikey: str, source_records: Sequence[Any]
) -> list[str]:
    reasons = ["degenerate_parent_graph"]
    if any(
        record.structure.exact_inchikey != parent_inchikey
        for record in source_records
    ):
        reasons.append("unsafe_parent_collapse")
    return reasons


def _unsupported_parent_elements(heavy_atoms: Sequence[Any]) -> list[str]:
    return sorted(
        {
            atom.GetSymbol()
            for atom in heavy_atoms
            if atom.GetAtomicNum() not in _SUPPORTED_ATOMIC_NUMBERS
        }
    )


def _write_sidecar(
    row: Mapping[str, Any],
    pdbqt_path: Path,
    prep_context: Mapping[str, Any],
    *,
    writer: str,
) -> None:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "writer": writer,
        "canonical_parent_id": row["canonical_parent_id"],
        "preferred_name": row["preferred_name"],
        "parent_smiles": row["parent_smiles"],
        "parent_inchikey": row["parent_inchikey"],
        "drugcentral_ids": str(row["drugcentral_ids"]).split(";") if row["drugcentral_ids"] else [],
        "source_record_indices": str(row["source_record_indices"]).split(";"),
        "approved_full_form_names": str(row["approved_full_form_names"]).split(";"),
        "approved_full_form_exact_inchikeys": str(
            row["approved_full_form_exact_inchikeys"]
        ).split(";"),
        "approved_form_relation": row["approved_form_relation"],
        "pdbqt_path": str(pdbqt_path),
        "pdbqt_sha256": _sha256(pdbqt_path),
        "pdbqt_bytes": pdbqt_path.stat().st_size,
        "prep_policy": PREP_POLICY,
        "ph_policy": "preserve_standardized_parent_charge; no synthetic pH microstate",
        "rdkit_version": prep_context["rdkit_version"],
        "meeko_version": prep_context["meeko_version"],
        "standardization_version": STANDARDIZATION_VERSION,
        "input_paths": dict(prep_context["input_paths"]),
        "input_sha256": dict(prep_context["input_sha256"]),
        "source_parent_sdf": prep_context["source_parent_sdf"],
        "source_parent_sdf_sha256": prep_context["source_parent_sdf_sha256"],
        "legacy_selected_rdk_id": row["legacy_selected_rdk_id"],
        "legacy_selected_path": row["legacy_selected_path"],
        "legacy_selected_sha256": row["legacy_selected_sha256"],
        "trigger_rdk_ids": str(row["trigger_rdk_ids"]).split(";")
        if row["trigger_rdk_ids"]
        else [],
    }
    _write_json(Path(str(row["sidecar_path"])), payload)


def _finish_ready_row(row: dict[str, Any], path: Path) -> None:
    row["materialization_status"] = "ready"
    row["status_reason_codes"] = ""
    row["output_pdbqt_sha256"] = _sha256(path)
    row["output_pdbqt_bytes"] = path.stat().st_size
    sidecar = Path(str(row["sidecar_path"]))
    row["sidecar_sha256"] = _sha256(sidecar) if sidecar.is_file() else ""


def _ready_existing(
    path: Path,
    parent_inchikey: object,
    expected_hash: str,
    prep_context: Mapping[str, Any],
) -> bool:
    if not _existing_pdbqt_matches(path, expected_hash):
        return False
    sidecar = path.with_suffix(".ligprep_source.json")
    payload = _read_sidecar_payload(sidecar)
    return bool(
        payload
        and _sidecar_context_matches(payload, parent_inchikey, prep_context)
    )


def _existing_pdbqt_matches(path: Path, expected_hash: str) -> bool:
    if not path.is_file() or path.stat().st_size <= 0:
        return False
    return not expected_hash or _sha256(path) == expected_hash


def _read_sidecar_payload(path: Path) -> Mapping[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def _sidecar_context_matches(
    payload: Mapping[str, Any],
    parent_inchikey: object,
    prep_context: Mapping[str, Any],
) -> bool:
    checks = (
        str(payload.get("parent_inchikey") or "") == str(parent_inchikey),
        payload.get("input_sha256") == dict(prep_context["input_sha256"]),
        payload.get("source_parent_sdf_sha256")
        == prep_context["source_parent_sdf_sha256"],
        payload.get("rdkit_version") == prep_context["rdkit_version"],
        payload.get("meeko_version") == prep_context["meeko_version"],
    )
    return all(checks)


def _write_parent_sdf(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    parent_mols: Mapping[str, Chem.Mol],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    writer = Chem.SDWriter(str(tmp))
    try:
        for row in sorted(rows, key=lambda item: str(item["parent_inchikey"])):
            mol = Chem.Mol(parent_mols[str(row["parent_inchikey"])])
            mol.SetProp("_Name", str(row["preferred_name"]))
            mol.SetProp("PARENT_INCHIKEY", str(row["parent_inchikey"]))
            mol.SetProp("DRUGCENTRAL_IDS", str(row["drugcentral_ids"]))
            mol.SetProp("SOURCE_RECORD_INDICES", str(row["source_record_indices"]))
            mol.SetProp("APPROVED_FULL_FORM_NAMES", str(row["approved_full_form_names"]))
            writer.write(mol)
    finally:
        writer.close()
    os.replace(tmp, path)


def _validate_ready_outputs(
    rows: Sequence[Mapping[str, Any]], ready_paths: Sequence[Path], output_dir: Path
) -> None:
    expected = _expected_ready_paths(rows)
    actual = {
        path.resolve()
        for path in output_dir.glob("*.pdbqt")
        if path.is_file() and not path.is_symlink()
    }
    if expected != actual or expected != {path.resolve() for path in ready_paths}:
        raise RuntimeError("FDA named library contains unmanifested or missing PDBQTs")
    _validate_unique_ready_hashes(rows)


def _expected_ready_paths(rows: Sequence[Mapping[str, Any]]) -> set[Path]:
    return {
        Path(str(row["output_pdbqt_path"])).resolve()
        for row in rows
        if row["materialization_status"] == "ready"
    }


def _validate_unique_ready_hashes(rows: Sequence[Mapping[str, Any]]) -> None:
    hashes: dict[str, str] = {}
    for row in rows:
        if row["materialization_status"] != "ready":
            continue
        digest = str(row["output_pdbqt_sha256"])
        key = str(row["parent_inchikey"])
        other = hashes.get(digest)
        if other and other != key:
            raise RuntimeError(f"PDBQT bytes shared by distinct FDA parents: {other}, {key}")
        hashes[digest] = key


def _plan_summary(
    rows: Sequence[Mapping[str, Any]],
    unsafe_records: Sequence[Any],
    input_paths: Mapping[str, str],
    input_hashes: Mapping[str, str],
    output_dir: Path,
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "plan_only": True,
        "slug_version": SLUG_VERSION,
        "prep_policy": PREP_POLICY,
        "standardization_version": STANDARDIZATION_VERSION,
        "rdkit_version": rdBase.rdkitVersion,
        "meeko_version": _package_version("meeko"),
        "input_paths": dict(input_paths),
        "input_sha256": dict(input_hashes),
        "output_dir": str(output_dir),
        "canonical_parent_rows": len(rows),
        "chemistry_trigger_rdk_ids": _chemistry_trigger_ids(rows),
        "unsafe_approved_source_records_excluded": len(unsafe_records),
        "action_counts": dict(
            sorted(Counter(str(row["materialization_action"]) for row in rows).items())
        ),
        "status_counts": dict(
            sorted(Counter(str(row["materialization_status"]) for row in rows).items())
        ),
        "ready_pdbqt_count": 0,
        "policy": {
            "legacy_inputs_mutated": False,
            "docking_launched": False,
            "plan_wrote_outputs": False,
            "filename_is_identity_authority": False,
        },
    }


def _summary(
    *,
    rows: Sequence[Mapping[str, Any]],
    unsafe_records: Sequence[Any],
    input_paths: Mapping[str, str],
    input_hashes: Mapping[str, str],
    prep_context: Mapping[str, Any],
    manifest_path: Path,
    quarantine_path: Path,
    library_index: Path,
    output_dir: Path,
) -> dict[str, Any]:
    status_counts = Counter(str(row["materialization_status"]) for row in rows)
    action_counts = Counter(str(row["materialization_action"]) for row in rows)
    ready = [row for row in rows if row["materialization_status"] == "ready"]
    return {
        "schema_version": SCHEMA_VERSION,
        "slug_version": SLUG_VERSION,
        "prep_policy": PREP_POLICY,
        "standardization_version": STANDARDIZATION_VERSION,
        "rdkit_version": rdBase.rdkitVersion,
        "meeko_version": prep_context["meeko_version"],
        "input_paths": dict(input_paths),
        "input_sha256": dict(input_hashes),
        "source_parent_sdf": prep_context["source_parent_sdf"],
        "source_parent_sdf_sha256": prep_context["source_parent_sdf_sha256"],
        "canonical_parent_rows": len(rows),
        "chemistry_trigger_rdk_ids": _chemistry_trigger_ids(rows),
        "unsafe_approved_source_records_excluded": len(unsafe_records),
        "action_counts": dict(sorted(action_counts.items())),
        "status_counts": dict(sorted(status_counts.items())),
        "ready_pdbqt_count": len(ready),
        "manifest_csv": str(manifest_path),
        "manifest_sha256": _sha256(manifest_path),
        "quarantine_csv": str(quarantine_path),
        "quarantine_sha256": _sha256(quarantine_path),
        "library_index": str(library_index),
        "library_index_sha256": _sha256(library_index),
        "library_inventory": _inventory(output_dir),
        "policy": {
            "legacy_inputs_mutated": False,
            "docking_launched": False,
            "filename_is_identity_authority": False,
            "approved_full_form_and_docked_parent_are_separate": True,
            "non_dockable_approved_substances_are_quarantined": True,
        },
    }


def _inventory(output_dir: Path) -> dict[str, Any]:
    paths = sorted(path for path in output_dir.glob("*.pdbqt") if path.is_file())
    digest = hashlib.sha256()
    for path in paths:
        relative = path.relative_to(output_dir).as_posix()
        digest.update(relative.encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return {"count": len(paths), "sha256": digest.hexdigest()}


def _chemistry_trigger_ids(rows: Sequence[Mapping[str, Any]]) -> list[str]:
    values: set[str] = set()
    for row in rows:
        values.update(
            token
            for token in str(row.get("trigger_rdk_ids") or "").split(";")
            if token
        )
    return sorted(values)


def _trigger_rdk_ids(
    repaired: Sequence[Mapping[str, str]],
) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = defaultdict(list)
    for row in repaired:
        if row.get("chemistry_redock_required") != "true":
            continue
        key = row.get("approved_parent_inchikey", "")
        if key:
            grouped[key].append(row.get("rdk_id", ""))
    return dict(grouped)


def _form_relation(records: Sequence[Any]) -> str:
    if all(record.structure.exact_inchikey == record.structure.parent_inchikey for record in records):
        return "full_equals_parent"
    classes = {
        _classify_substance(record.structure, record.name, False) for record in records
    }
    if SALT_OR_SOLVATE in classes:
        return "counterion_or_solvate_removed"
    return "charge_or_metal_normalized"


def _record_id(record: Any) -> str:
    return str(record.identifiers.get("id") or "").strip()


def _slug(value: object) -> str:
    text = unicodedata.normalize("NFKD", str(value or ""))
    ascii_text = text.encode("ascii", "ignore").decode("ascii").casefold()
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", ascii_text)).strip("_")


def _joined(values: Iterable[object]) -> str:
    return ";".join(sorted({str(value).strip() for value in values if str(value).strip()}))


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _write_csv(
    path: Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str]
) -> None:
    tmp = path.with_suffix(path.suffix + ".part")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _validate_resume_inputs(summary_path: Path, input_hashes: Mapping[str, str]) -> None:
    if not summary_path.is_file():
        return
    try:
        old = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"invalid existing FDA named-library summary: {summary_path}") from exc
    if old.get("input_sha256") != dict(input_hashes):
        raise ValueError("FDA named-library inputs changed; choose a new output directory")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a named, manifest-backed FDA canonical-parent PDBQT library."
    )
    parser.add_argument("--repaired-mapping-csv", type=Path, required=True)
    parser.add_argument("--redock-delta-csv", type=Path, required=True)
    parser.add_argument("--source-sdf", type=Path, required=True)
    parser.add_argument("--approval-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--plan-only", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    summary = build_named_fda_library(
        repaired_mapping_csv=args.repaired_mapping_csv,
        redock_delta_csv=args.redock_delta_csv,
        source_sdf=args.source_sdf,
        approval_manifest=args.approval_manifest,
        output_dir=args.output_dir,
        workers=args.workers,
        plan_only=args.plan_only,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = ["build_named_fda_library", "main"]
