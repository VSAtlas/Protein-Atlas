"""Identity-safe staging of measured BindingDB FDA inactive pairs.

This module deliberately permits only exact, full InChIKey joins to the active
FDA mapping. Absence of an assay, parent-InChIKey similarity, and drug-name
similarity never create a negative label or an identity match.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd


DIRECT_ENDPOINTS = frozenset({"KI", "KD", "IC50"})
SAFE_NEGATIVE_RELATIONS = frozenset({"=", ">", ">="})
SAFE_ACTIVE_RELATIONS = frozenset({"=", "<", "<="})
ACCEPTED_FDA_REGULATORY_STATUSES = frozenset(
    {"drugcentral_fda_approved", "fda_approved_current_or_historical"}
)
CONFIRMED_TERMINAL_DISPOSITIONS = frozenset(
    {
        "confirmed_fda_active_ingredient",
        "confirmed_fda_salt_parent_docked",
        "incorrect_name_structure_mapping_resolved",
    }
)
PROJECTED_LABEL_STATUS = "strict_measured_matched_inactive"
IDENTITY_MATCH_METHOD = "exact_full_inchikey"

_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_SAFE_TOKEN = re.compile(r"^[A-Za-z0-9_.-]+$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")

BINDINGDB_REQUIRED_COLUMNS = frozenset(
    {
        "drug_id",
        "drug_name",
        "target_id",
        "uniprot",
        "assay_id",
        "activity_nM",
        "activity_relation",
        "activity_type",
        "activity_units",
        "source",
        "source_specific_activity_label",
        "bindingdb_reactant_set_id",
        "bindingdb_monomerid",
        "bindingdb_ligand_name",
        "target_name",
        "target_organism",
        "document_ids",
        "pubchem_aid",
        "publication_year",
        "inchikey",
        "smiles",
    }
)
MAPPING_REQUIRED_COLUMNS = frozenset(
    {
        "rdk_id",
        "identity_structure_validated",
        "identity_exact_inchikey",
        "regulatory_status",
        "terminal_disposition",
        "canonical_parent_id",
        "canonical_pdbqt_path",
        "canonical_pdbqt_sha256",
    }
)

SELECTED_COLUMNS = (
    "pdb_id",
    "ligand_base",
    "ligand_file_stem",
    "drug_id",
    "generic_name",
    "display_name",
    "target_gene",
    "target_uniprot",
    "target_id",
    "projected_label",
    "projected_label_status",
    "relation_domain",
    "canonical_pair_key",
    "source_name",
    "source_family",
    "representative_source",
    "upstream_source",
    "evidence_source_names",
    "evidence_sources",
    "source_label_policy",
    "evidence_type",
    "negative_evidence_type",
    "negative_source",
    "negative_confidence",
    "production_truth_allowed",
    "training_allowed",
    "benchmark_only",
    "probe_sensitivity_only",
    "activity_nM",
    "activity_relation",
    "activity_type",
    "activity_units",
    "assay_id",
    "assay_count",
    "source_row_count",
    "duplicate_assay_rows_collapsed",
    "bindingdb_assay_ids",
    "bindingdb_reactant_set_ids",
    "bindingdb_monomerids",
    "bindingdb_document_ids",
    "bindingdb_publication_years",
    "bindingdb_activity_types",
    "bindingdb_activity_relations",
    "bindingdb_activity_values_nm",
    "bindingdb_source_row_numbers",
    "bindingdb_assay_provenance_json",
    "inchikey",
    "identity_exact_inchikey",
    "identity_match_method",
    "identity_structure_validated",
    "mapping_match_count",
    "mapping_eligible_count",
    "mapping_selection_reason",
    "regulatory_status",
    "terminal_disposition",
    "canonical_parent_id",
    "canonical_pdbqt_path",
    "canonical_pdbqt_sha256",
    "pdbqt_path",
    "pdbqt_available",
)


@dataclass(frozen=True)
class TargetSpec:
    """A strict BindingDB target-to-Atlas-PDB mapping supplied by the caller."""

    uniprot: str
    gene: str
    pdb_id: str

    def normalized(self) -> "TargetSpec":
        uniprot = _clean(self.uniprot).upper()
        gene = _clean(self.gene).upper()
        pdb_id = _clean(self.pdb_id).upper()
        if not uniprot or not gene or not pdb_id:
            raise ValueError("target UniProt, gene, and PDB ID must all be nonempty")
        if not _SAFE_TOKEN.fullmatch(uniprot) or not _SAFE_TOKEN.fullmatch(gene):
            raise ValueError(f"unsafe target identifier: {self!r}")
        if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
            raise ValueError(f"PDB ID must contain exactly four letters/digits: {pdb_id}")
        return TargetSpec(uniprot=uniprot, gene=gene, pdb_id=pdb_id)


def build_target_specs(
    target_uniprots: Sequence[str],
    target_genes: Sequence[str],
    pdb_ids: Sequence[str],
) -> list[TargetSpec]:
    """Build aligned target specifications from repeatable CLI-style values."""

    lengths = {len(target_uniprots), len(target_genes), len(pdb_ids)}
    if lengths != {len(target_uniprots)} or not target_uniprots:
        raise ValueError(
            "target_uniprot, target_gene, and pdb_id must be nonempty aligned lists"
        )
    specs = [
        TargetSpec(uniprot, gene, pdb_id).normalized()
        for uniprot, gene, pdb_id in zip(
            target_uniprots, target_genes, pdb_ids, strict=True
        )
    ]
    if len(set(specs)) != len(specs):
        raise ValueError("duplicate target UniProt/gene/PDB specifications are not allowed")
    return specs


def stage_bindingdb_fda_inactives(
    *,
    bindingdb_path: str | Path,
    mapping_path: str | Path,
    targets: Sequence[TargetSpec],
    out_dir: str | Path,
    library_name: str,
    inactive_nm: float = 10_000.0,
    active_nm: float = 1_000.0,
    config_path: str | Path = "config.txt",
    prepped_root: str | Path | None = None,
    chunksize: int = 100_000,
) -> dict[str, Any]:
    """Stage strict BindingDB measured inactives backed by validated FDA PDBQTs.

    Only target-scoped measured rows can become labels. Untested pairs remain
    unknown, and identity resolution never falls back to parent InChIKey or name.
    """

    source_path = Path(bindingdb_path)
    active_mapping_path = Path(mapping_path)
    output = Path(out_dir)
    normalized_targets = [target.normalized() for target in targets]
    _validate_inputs(
        source_path=source_path,
        mapping_path=active_mapping_path,
        targets=normalized_targets,
        inactive_nm=inactive_nm,
        active_nm=active_nm,
        chunksize=chunksize,
        library_name=library_name,
    )
    output.mkdir(parents=True, exist_ok=True)

    source_rows, source_total = _read_target_rows(
        source_path, normalized_targets, chunksize=chunksize
    )
    expanded = _expand_target_rows(source_rows, normalized_targets)
    classified, exclusions = _classify_assay_rows(
        expanded, inactive_nm=inactive_nm, active_nm=active_nm
    )

    mapping = pd.read_csv(
        active_mapping_path, dtype=str, keep_default_na=False, low_memory=False
    )
    _require_columns(mapping, MAPPING_REQUIRED_COLUMNS, "FDA mapping")
    mapping = mapping.copy()
    mapping["_identity_exact_inchikey"] = mapping[
        "identity_exact_inchikey"
    ].map(_normalize_inchikey)
    mapping_groups = {
        key: group.copy()
        for key, group in mapping.loc[
            mapping["_identity_exact_inchikey"].map(_valid_inchikey)
        ].groupby("_identity_exact_inchikey", sort=False)
    }

    accepted, identity_exclusions = _resolve_candidate_identities(
        classified,
        mapping_groups=mapping_groups,
        mapping_path=active_mapping_path,
    )
    exclusions.extend(identity_exclusions)
    selected, duplicate_exclusions = _collapse_assay_rows(accepted)
    exclusions.extend(duplicate_exclusions)

    selected, copy_manifest, library_map = _copy_target_libraries(
        selected,
        library_name=library_name,
        config_path=Path(config_path),
        prepped_root=Path(prepped_root) if prepped_root is not None else None,
    )

    selected_path = output / "selected_matched_inactive_pairs.csv"
    exclusions_path = output / "bindingdb_fda_inactive_exclusions.csv"
    copy_manifest_path = output / "canonical_pdbqt_copy_manifest.csv"
    library_map_path = output / "per_pdb_library_map.json"
    manifest_path = output / "bindingdb_fda_inactive_manifest.json"

    selected.reindex(columns=SELECTED_COLUMNS).to_csv(selected_path, index=False)
    _exclusions_frame(exclusions).to_csv(exclusions_path, index=False)
    copy_manifest.to_csv(copy_manifest_path, index=False)
    library_map_path.write_text(
        json.dumps(library_map, indent=2, sort_keys=True), encoding="utf-8"
    )

    manifest: dict[str, Any] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "bindingdb": str(source_path),
        "bindingdb_file_size": source_path.stat().st_size,
        "bindingdb_file_mtime_ns": source_path.stat().st_mtime_ns,
        "mapping": str(active_mapping_path),
        "mapping_sha256": _sha256(active_mapping_path),
        "targets": [target.__dict__ for target in normalized_targets],
        "thresholds_nm": {
            "active_conflict_at_or_below": float(active_nm),
            "inactive_at_or_above": float(inactive_nm),
        },
        "counts": {
            "bindingdb_rows_total": int(source_total),
            "target_scoped_source_rows": int(len(source_rows)),
            "target_pdb_expanded_rows": int(len(expanded)),
            "safe_negative_source_rows": int(len(classified)),
            "identity_accepted_source_rows": int(len(accepted)),
            "selected_pairs": int(len(selected)),
            "selected_ligands": int(selected["ligand_base"].nunique())
            if not selected.empty
            else 0,
            "excluded_source_rows": int(len(exclusions)),
            "copied_pdbqts": int(copy_manifest["copy_verified"].sum())
            if not copy_manifest.empty
            else 0,
        },
        "policy": {
            "label": (
                "human, direct Ki/Kd/IC50 with a safe exact/lower-bound value "
                f">= {float(inactive_nm):g} nM"
            ),
            "active_conflict": (
                "same target UniProt and exact full InChIKey with a safe exact/"
                f"upper-bound value <= {float(active_nm):g} nM excludes the pair"
            ),
            "identity_join": IDENTITY_MATCH_METHOD,
            "parent_inchikey_matching": "prohibited",
            "name_only_matching": "prohibited",
            "absence_as_negative": "prohibited; untested/absent pairs remain unknown",
            "accepted_regulatory_statuses": sorted(
                ACCEPTED_FDA_REGULATORY_STATUSES
            ),
            "confirmed_terminal_dispositions": sorted(
                CONFIRMED_TERMINAL_DISPOSITIONS
            ),
            "canonical_pdbqt": (
                "source file must exist and both pre-copy and post-copy SHA256 "
                "must equal canonical_pdbqt_sha256"
            ),
            "duplicate_assays": (
                "collapsed by assay/source/endpoint/value identity while preserving "
                "all source-row provenance in JSON"
            ),
        },
        "outputs": {
            "selected_pairs": str(selected_path),
            "exclusions": str(exclusions_path),
            "copy_manifest": str(copy_manifest_path),
            "per_pdb_library_map": str(library_map_path),
            "manifest": str(manifest_path),
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return manifest


def _validate_inputs(
    *,
    source_path: Path,
    mapping_path: Path,
    targets: Sequence[TargetSpec],
    inactive_nm: float,
    active_nm: float,
    chunksize: int,
    library_name: str,
) -> None:
    if not source_path.is_file():
        raise FileNotFoundError(f"BindingDB table does not exist: {source_path}")
    if not mapping_path.is_file():
        raise FileNotFoundError(f"FDA mapping does not exist: {mapping_path}")
    if not targets:
        raise ValueError("at least one target specification is required")
    if not math.isfinite(active_nm) or active_nm <= 0:
        raise ValueError("active_nm must be finite and positive")
    if not math.isfinite(inactive_nm) or inactive_nm <= active_nm:
        raise ValueError("inactive_nm must be finite and greater than active_nm")
    if chunksize < 1:
        raise ValueError("chunksize must be positive")
    if not _SAFE_TOKEN.fullmatch(library_name):
        raise ValueError("library_name may contain only letters, digits, ., _, and -")


def _read_target_rows(
    path: Path, targets: Sequence[TargetSpec], *, chunksize: int
) -> tuple[pd.DataFrame, int]:
    wanted = {target.uniprot for target in targets}
    frames: list[pd.DataFrame] = []
    total = 0
    for chunk in pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
        chunksize=chunksize,
    ):
        _require_columns(chunk, BINDINGDB_REQUIRED_COLUMNS, "BindingDB table")
        chunk = chunk.copy()
        chunk["_source_row_number"] = range(total + 2, total + len(chunk) + 2)
        total += len(chunk)
        mask = chunk["uniprot"].map(_clean).str.upper().isin(wanted)
        if mask.any():
            frames.append(chunk.loc[mask].copy())
    if frames:
        return pd.concat(frames, ignore_index=True), total
    columns = list(BINDINGDB_REQUIRED_COLUMNS) + ["_source_row_number"]
    return pd.DataFrame(columns=columns), total


def _expand_target_rows(
    rows: pd.DataFrame, targets: Sequence[TargetSpec]
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for target in targets:
        subset = rows.loc[
            rows["uniprot"].map(_clean).str.upper().eq(target.uniprot)
        ].copy()
        subset["target_uniprot"] = target.uniprot
        subset["target_gene"] = target.gene
        subset["pdb_id"] = target.pdb_id
        frames.append(subset)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return rows.assign(target_uniprot="", target_gene="", pdb_id="").iloc[0:0]


def _classify_assay_rows(
    rows: pd.DataFrame, *, inactive_nm: float, active_nm: float
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if rows.empty:
        return rows.copy(), []
    classified = rows.copy()
    classified["_endpoint"] = classified["activity_type"].map(_normalize_endpoint)
    classified["_relation"] = classified["activity_relation"].map(
        _normalize_relation
    )
    classified["_activity_nm"] = pd.to_numeric(
        classified["activity_nM"], errors="coerce"
    )
    classified["_inchikey"] = classified["inchikey"].map(_normalize_inchikey)
    classified["_human"] = classified["target_organism"].map(_is_human)

    active_mask = (
        classified["_human"]
        & classified["_endpoint"].isin(DIRECT_ENDPOINTS)
        & classified["_relation"].isin(SAFE_ACTIVE_RELATIONS)
        & classified["_activity_nm"].notna()
        & classified["_activity_nm"].le(active_nm)
        & classified["_inchikey"].map(_valid_inchikey)
    )
    active_conflicts = set(
        classified.loc[active_mask, ["target_uniprot", "_inchikey"]].itertuples(
            index=False, name=None
        )
    )

    accepted_positions: list[int] = []
    exclusions: list[dict[str, Any]] = []
    for position, row in classified.iterrows():
        reasons: list[str] = []
        activity = row["_activity_nm"]
        if not bool(row["_human"]):
            reasons.append("non_human_target")
        if row["_endpoint"] not in DIRECT_ENDPOINTS:
            reasons.append("not_direct_ki_kd_ic50")
        if pd.isna(activity) or not math.isfinite(float(activity)) or float(activity) <= 0:
            reasons.append("invalid_activity_nm")
        elif float(activity) < inactive_nm:
            reasons.append("below_inactive_threshold")
        if row["_relation"] not in SAFE_NEGATIVE_RELATIONS:
            reasons.append("relation_does_not_prove_inactivity")
        if not _valid_inchikey(row["_inchikey"]):
            reasons.append("missing_or_invalid_full_inchikey")
        conflict_key = (row["target_uniprot"], row["_inchikey"])
        if _valid_inchikey(row["_inchikey"]) and conflict_key in active_conflicts:
            reasons.append("same_target_measured_active_conflict")
        if reasons:
            exclusions.append(_exclusion_record(row, "assay_filter", reasons))
        else:
            accepted_positions.append(position)
    return classified.loc[accepted_positions].copy(), exclusions


def _resolve_candidate_identities(
    candidates: pd.DataFrame,
    *,
    mapping_groups: Mapping[str, pd.DataFrame],
    mapping_path: Path,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if candidates.empty:
        return candidates.copy(), []
    accepted: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    file_hash_cache: dict[Path, str] = {}
    selection_cache: dict[tuple[str, str], tuple[dict[str, Any] | None, str]] = {}
    for _, row in candidates.iterrows():
        inchikey = _clean(row["_inchikey"])
        source_drug_id = _clean(row.get("drug_id")).lower()
        cache_key = (inchikey, source_drug_id)
        if cache_key not in selection_cache:
            selection_cache[cache_key] = _select_mapping_row(
                mapping_groups.get(inchikey),
                source_drug_id=source_drug_id,
                mapping_path=mapping_path,
                file_hash_cache=file_hash_cache,
            )
        selected_mapping, selection_reason = selection_cache[cache_key]
        if selected_mapping is None:
            exclusions.append(
                _exclusion_record(row, "identity_filter", [selection_reason])
            )
            continue
        payload = row.to_dict()
        payload.update(selected_mapping)
        payload["_mapping_selection_reason"] = selection_reason
        accepted.append(payload)
    return pd.DataFrame(accepted), exclusions


def _select_mapping_row(
    rows: pd.DataFrame | None,
    *,
    source_drug_id: str,
    mapping_path: Path,
    file_hash_cache: dict[Path, str],
) -> tuple[dict[str, Any] | None, str]:
    if rows is None or rows.empty:
        return None, "no_exact_full_inchikey_mapping"
    eligible: list[dict[str, Any]] = []
    failure_reasons: set[str] = set()
    for _, row in rows.iterrows():
        reasons: list[str] = []
        if not _truthy(row.get("identity_structure_validated")):
            reasons.append("identity_structure_not_validated")
        if _clean(row.get("regulatory_status")) not in ACCEPTED_FDA_REGULATORY_STATUSES:
            reasons.append("regulatory_status_not_accepted")
        if _clean(row.get("terminal_disposition")) not in CONFIRMED_TERMINAL_DISPOSITIONS:
            reasons.append("terminal_disposition_not_confirmed")
        if not _clean(row.get("rdk_id")):
            reasons.append("missing_rdk_id")
        expected = _clean(row.get("canonical_pdbqt_sha256")).lower()
        if not _SHA256.fullmatch(expected):
            reasons.append("missing_or_invalid_canonical_pdbqt_sha256")
        canonical = _resolve_canonical_path(
            _clean(row.get("canonical_pdbqt_path")), mapping_path
        )
        if canonical is None or not canonical.is_file():
            reasons.append("canonical_pdbqt_missing")
        actual = ""
        if canonical is not None and canonical.is_file() and _SHA256.fullmatch(expected):
            actual = file_hash_cache.setdefault(canonical, _sha256(canonical))
            if actual.lower() != expected:
                reasons.append("canonical_pdbqt_sha256_mismatch")
        if reasons:
            failure_reasons.update(reasons)
            continue
        payload = row.to_dict()
        payload["_canonical_pdbqt_resolved"] = str(canonical)
        payload["_canonical_pdbqt_actual_sha256"] = actual
        eligible.append(payload)
    if not eligible:
        suffix = ",".join(sorted(failure_reasons)) or "no_eligible_mapping_row"
        return None, f"exact_mapping_failed:{suffix}"

    frame = pd.DataFrame(eligible)
    if frame["_canonical_pdbqt_actual_sha256"].nunique() != 1:
        return None, "ambiguous_exact_identity_multiple_canonical_structures"
    frame["_source_id_match"] = frame["rdk_id"].map(_clean).str.lower().eq(
        source_drug_id
    )
    terminal_order = {
        "confirmed_fda_active_ingredient": 0,
        "confirmed_fda_salt_parent_docked": 1,
        "incorrect_name_structure_mapping_resolved": 2,
    }
    regulatory_order = {
        "drugcentral_fda_approved": 0,
        "fda_approved_current_or_historical": 1,
    }
    frame["_terminal_order"] = frame["terminal_disposition"].map(
        terminal_order
    ).fillna(99)
    frame["_regulatory_order"] = frame["regulatory_status"].map(
        regulatory_order
    ).fillna(99)
    frame = frame.sort_values(
        ["_source_id_match", "_terminal_order", "_regulatory_order", "rdk_id"],
        ascending=[False, True, True, True],
        kind="stable",
    )
    selected = frame.iloc[0].to_dict()
    selected["_mapping_match_count"] = int(len(rows))
    selected["_mapping_eligible_count"] = int(len(frame))
    reason = (
        "source_drug_id_and_exact_inchikey"
        if bool(selected["_source_id_match"])
        else "deterministic_exact_inchikey_representative"
    )
    return selected, reason


def _collapse_assay_rows(
    accepted: pd.DataFrame,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if accepted.empty:
        return pd.DataFrame(columns=SELECTED_COLUMNS), []
    selected_rows: list[dict[str, Any]] = []
    duplicate_exclusions: list[dict[str, Any]] = []
    for (pdb_id, ligand_base), group in accepted.groupby(
        ["pdb_id", "rdk_id"], sort=True, dropna=False
    ):
        if group["_inchikey"].nunique() != 1:
            for _, row in group.iterrows():
                duplicate_exclusions.append(
                    _exclusion_record(
                        row, "pair_collapse", ["pair_maps_multiple_exact_inchikeys"]
                    )
                )
            continue
        evidence_key_columns = [
            "assay_id",
            "bindingdb_reactant_set_id",
            "bindingdb_monomerid",
            "_endpoint",
            "_relation",
            "_activity_nm",
            "source",
            "document_ids",
        ]
        evidence: list[dict[str, Any]] = []
        duplicate_count = 0
        for _, evidence_group in group.groupby(
            evidence_key_columns, sort=False, dropna=False
        ):
            first = evidence_group.iloc[0]
            source_rows = sorted(
                int(value) for value in evidence_group["_source_row_number"].tolist()
            )
            duplicate_count += max(0, len(evidence_group) - 1)
            evidence.append(
                {
                    "source_row_numbers": source_rows,
                    "assay_id": _clean(first.get("assay_id")),
                    "reactant_set_id": _clean(
                        first.get("bindingdb_reactant_set_id")
                    ),
                    "monomer_id": _clean(first.get("bindingdb_monomerid")),
                    "endpoint": _clean(first.get("_endpoint")),
                    "relation": _clean(first.get("_relation")),
                    "activity_nm": float(first["_activity_nm"]),
                    "source": _clean(first.get("source")),
                    "source_call": _clean(
                        first.get("source_specific_activity_label")
                    ),
                    "document_ids": _clean(first.get("document_ids")),
                    "publication_year": _clean(first.get("publication_year")),
                    "target_name": _clean(first.get("target_name")),
                }
            )
            for _, duplicate in evidence_group.iloc[1:].iterrows():
                duplicate_exclusions.append(
                    _exclusion_record(
                        duplicate,
                        "pair_collapse",
                        ["duplicate_assay_row_collapsed_with_provenance"],
                    )
                )
        representative = group.assign(
            _endpoint_order=group["_endpoint"].map({"KD": 0, "KI": 1, "IC50": 2})
        ).sort_values(
            ["_endpoint_order", "_activity_nm", "_source_row_number"],
            kind="stable",
        ).iloc[0]
        source_names = _join_unique(group["source"])
        preferred_name = _first_nonempty(
            representative,
            "preferred_identity",
            "resolved_preferred_name",
            "display_name",
            "generic_name",
            "bindingdb_ligand_name",
            "drug_name",
        )
        exact_inchikey = _clean(representative["_inchikey"])
        source_row_numbers = sorted(
            int(value) for value in group["_source_row_number"].tolist()
        )
        selected_rows.append(
            {
                "pdb_id": _clean(pdb_id).upper(),
                "ligand_base": _clean(ligand_base).lower(),
                "ligand_file_stem": _clean(ligand_base).lower(),
                "drug_id": preferred_name,
                "generic_name": preferred_name,
                "display_name": preferred_name,
                "target_gene": _clean(representative["target_gene"]).upper(),
                "target_uniprot": _clean(
                    representative["target_uniprot"]
                ).upper(),
                "target_id": _clean(representative["target_gene"]).upper(),
                "projected_label": 0,
                "projected_label_status": PROJECTED_LABEL_STATUS,
                "relation_domain": "drug_target_activity",
                "canonical_pair_key": (
                    f"drug_target_activity::{_clean(ligand_base).lower()}::"
                    f"{_clean(representative['target_uniprot']).upper()}"
                ),
                "source_name": "BindingDB",
                "source_family": "BindingDB",
                "representative_source": "BindingDB",
                "upstream_source": "BindingDB",
                "evidence_source_names": source_names or "BindingDB",
                "evidence_sources": source_names or "BindingDB",
                "source_label_policy": PROJECTED_LABEL_STATUS,
                "evidence_type": "measured_inactive_direct_binding",
                "negative_evidence_type": "measured_inactive_direct_binding",
                "negative_source": "BindingDB",
                "negative_confidence": "strict",
                "production_truth_allowed": True,
                "training_allowed": True,
                "benchmark_only": False,
                "probe_sensitivity_only": False,
                "activity_nM": float(representative["_activity_nm"]),
                "activity_relation": _clean(representative["_relation"]),
                "activity_type": _clean(representative["_endpoint"]),
                "activity_units": "nM",
                "assay_id": _clean(representative.get("assay_id")),
                "assay_count": len(evidence),
                "source_row_count": len(group),
                "duplicate_assay_rows_collapsed": duplicate_count,
                "bindingdb_assay_ids": _join_unique(group["assay_id"]),
                "bindingdb_reactant_set_ids": _join_unique(
                    group["bindingdb_reactant_set_id"]
                ),
                "bindingdb_monomerids": _join_unique(
                    group["bindingdb_monomerid"]
                ),
                "bindingdb_document_ids": _join_unique(group["document_ids"]),
                "bindingdb_publication_years": _join_unique(
                    group["publication_year"]
                ),
                "bindingdb_activity_types": _join_unique(group["_endpoint"]),
                "bindingdb_activity_relations": _join_unique(group["_relation"]),
                "bindingdb_activity_values_nm": _join_unique(
                    group["_activity_nm"].map(lambda value: f"{float(value):g}")
                ),
                "bindingdb_source_row_numbers": "|".join(
                    str(value) for value in source_row_numbers
                ),
                "bindingdb_assay_provenance_json": json.dumps(
                    evidence, sort_keys=True, separators=(",", ":")
                ),
                "inchikey": exact_inchikey,
                "identity_exact_inchikey": exact_inchikey,
                "identity_match_method": IDENTITY_MATCH_METHOD,
                "identity_structure_validated": True,
                "mapping_match_count": int(
                    representative.get("_mapping_match_count", 0)
                ),
                "mapping_eligible_count": int(
                    representative.get("_mapping_eligible_count", 0)
                ),
                "mapping_selection_reason": _clean(
                    representative.get("_mapping_selection_reason")
                ),
                "regulatory_status": _clean(
                    representative.get("regulatory_status")
                ),
                "terminal_disposition": _clean(
                    representative.get("terminal_disposition")
                ),
                "canonical_parent_id": _clean(
                    representative.get("canonical_parent_id")
                ),
                "canonical_pdbqt_path": _clean(
                    representative.get("_canonical_pdbqt_resolved")
                ),
                "canonical_pdbqt_sha256": _clean(
                    representative.get("_canonical_pdbqt_actual_sha256")
                ),
                "pdbqt_path": "",
                "pdbqt_available": False,
            }
        )
    selected = pd.DataFrame(selected_rows, columns=SELECTED_COLUMNS)
    if not selected.empty and selected.duplicated(["pdb_id", "ligand_base"]).any():
        raise RuntimeError("internal error: selected pairs are not unique")
    return selected, duplicate_exclusions


def _copy_target_libraries(
    selected: pd.DataFrame,
    *,
    library_name: str,
    config_path: Path,
    prepped_root: Path | None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    copy_columns = [
        "pdb_id",
        "ligand_base",
        "library",
        "canonical_source_path",
        "destination_path",
        "expected_sha256",
        "source_sha256",
        "destination_sha256",
        "copy_status",
        "copy_verified",
    ]
    if selected.empty:
        return selected.copy(), pd.DataFrame(columns=copy_columns), {}
    root = prepped_root or _prepped_root(config_path, _clean(selected.iloc[0]["pdb_id"]))
    root.mkdir(parents=True, exist_ok=True)
    root_resolved = root.resolve()
    output = selected.copy()
    copy_rows: list[dict[str, Any]] = []
    library_map: dict[str, str] = {}
    target_dirs: set[Path] = set()
    for index, row in output.iterrows():
        pdb_id = _clean(row["pdb_id"]).upper()
        ligand_base = _clean(row["ligand_base"]).lower()
        if not _SAFE_TOKEN.fullmatch(ligand_base):
            raise RuntimeError(f"unsafe ligand_base from FDA mapping: {ligand_base!r}")
        library = f"{library_name}_{pdb_id}"
        target_dir = root / library
        target_dir.resolve().relative_to(root_resolved)
        target_dir.mkdir(parents=True, exist_ok=True)
        target_dirs.add(target_dir)
        library_map[pdb_id] = library
        source = Path(_clean(row["canonical_pdbqt_path"]))
        expected = _clean(row["canonical_pdbqt_sha256"]).lower()
        source_hash = _sha256(source)
        if source_hash != expected:
            raise RuntimeError(
                f"canonical PDBQT changed after validation: {source} "
                f"expected={expected} actual={source_hash}"
            )
        destination = target_dir / f"{ligand_base}.pdbqt"
        if destination.exists():
            existing_hash = _sha256(destination)
            if existing_hash != expected:
                raise FileExistsError(
                    f"refusing to overwrite mismatched staged PDBQT: {destination}"
                )
            copy_status = "reused_identical"
        else:
            partial = destination.with_suffix(".pdbqt.partial")
            if partial.exists():
                partial.unlink()
            shutil.copy2(source, partial)
            partial_hash = _sha256(partial)
            if partial_hash != expected:
                partial.unlink(missing_ok=True)
                raise RuntimeError(
                    f"copied PDBQT checksum mismatch: {source} -> {destination}"
                )
            os.replace(partial, destination)
            copy_status = "copied_verified"
        destination_hash = _sha256(destination)
        verified = destination_hash == expected
        if not verified:
            raise RuntimeError(f"post-copy checksum mismatch: {destination}")
        output.at[index, "pdbqt_path"] = str(destination)
        output.at[index, "pdbqt_available"] = True
        copy_rows.append(
            {
                "pdb_id": pdb_id,
                "ligand_base": ligand_base,
                "library": library,
                "canonical_source_path": str(source),
                "destination_path": str(destination),
                "expected_sha256": expected,
                "source_sha256": source_hash,
                "destination_sha256": destination_hash,
                "copy_status": copy_status,
                "copy_verified": verified,
            }
        )
    from prep_ligands.library_index import LibraryIndex

    for target_dir in sorted(target_dirs):
        LibraryIndex().load([target_dir])
    return output, pd.DataFrame(copy_rows, columns=copy_columns), library_map


def _prepped_root(config_path: Path, pdb_id: str) -> Path:
    from config.runtime_config import load_config
    from path_router import make_paths

    repo_root = Path(__file__).resolve().parents[2]
    config = config_path if config_path.is_absolute() else repo_root / config_path
    cfg = load_config(str(config), base_dir=repo_root)
    return make_paths(cfg, pdb_id, f"{pdb_id}.pdb").prepped_root


def _resolve_canonical_path(value: str, mapping_path: Path) -> Path | None:
    if not value:
        return None
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    mapping_resolved = mapping_path.resolve()
    repo_root = (
        mapping_resolved.parent.parent.parent
        if mapping_resolved.parent.name == "data"
        and mapping_resolved.parent.parent.name == "chemdb"
        else Path.cwd().resolve()
    )
    candidates = [repo_root / path, mapping_resolved.parent / path, Path.cwd() / path]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return candidates[0].resolve()


def _exclusion_record(
    row: Mapping[str, Any] | pd.Series, stage: str, reasons: Iterable[str]
) -> dict[str, Any]:
    payload = {
        column: row.get(column, "")
        for column in sorted(BINDINGDB_REQUIRED_COLUMNS)
    }
    payload.update(
        {
            "source_row_number": row.get("_source_row_number", ""),
            "pdb_id": row.get("pdb_id", ""),
            "target_gene": row.get("target_gene", ""),
            "target_uniprot": row.get("target_uniprot", row.get("uniprot", "")),
            "normalized_endpoint": row.get("_endpoint", ""),
            "normalized_relation": row.get("_relation", ""),
            "normalized_activity_nm": row.get("_activity_nm", ""),
            "normalized_exact_inchikey": row.get("_inchikey", ""),
            "exclusion_stage": stage,
            "exclude_reason": "|".join(sorted(set(reasons))),
            "identity_policy": IDENTITY_MATCH_METHOD,
            "absence_as_negative": False,
        }
    )
    return payload


def _exclusions_frame(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    if rows:
        return pd.DataFrame(rows)
    columns = sorted(BINDINGDB_REQUIRED_COLUMNS) + [
        "source_row_number",
        "pdb_id",
        "target_gene",
        "target_uniprot",
        "normalized_endpoint",
        "normalized_relation",
        "normalized_activity_nm",
        "normalized_exact_inchikey",
        "exclusion_stage",
        "exclude_reason",
        "identity_policy",
        "absence_as_negative",
    ]
    return pd.DataFrame(columns=columns)


def _require_columns(frame: pd.DataFrame, required: Iterable[str], name: str) -> None:
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{name} is missing required columns: {', '.join(missing)}")


def _normalize_endpoint(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]", "", _clean(value).upper())


def _normalize_relation(value: Any) -> str:
    relation = _clean(value).replace("≤", "<=").replace("≥", ">=").upper()
    aliases = {"==": "=", "EQ": "=", "GT": ">", "GE": ">=", "LT": "<", "LE": "<="}
    return aliases.get(relation, relation)


def _normalize_inchikey(value: Any) -> str:
    return _clean(value).upper()


def _valid_inchikey(value: Any) -> bool:
    return bool(_INCHIKEY.fullmatch(_clean(value).upper()))


def _is_human(value: Any) -> bool:
    organism = re.sub(r"\s+", " ", _clean(value).casefold())
    return organism in {"human", "homo sapiens", "9606"}


def _truthy(value: Any) -> bool:
    return _clean(value).casefold() in {"1", "true", "t", "yes", "y"}


def _clean(value: Any) -> str:
    if value is None:
        return ""
    try:
        if bool(pd.isna(value)):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _first_nonempty(row: Mapping[str, Any] | pd.Series, *columns: str) -> str:
    for column in columns:
        value = _clean(row.get(column, ""))
        if value:
            return value
    return ""


def _join_unique(values: Iterable[Any]) -> str:
    unique = sorted({_clean(value) for value in values if _clean(value)})
    return "|".join(unique)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = ["TargetSpec", "build_target_specs", "stage_bindingdb_fda_inactives"]
