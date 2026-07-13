from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.source_tables import read_source_table


DEFAULT_LABEL_CANDIDATES = (
    "spd_binding_label",
    "combined_activity_ml_label",
    "ml_binary_label",
)
DIRECT_ENDPOINTS = {"KI", "KD", "IC50", "EC50", "AC50"}
ACTIVE_TEXT = {
    "1",
    "1.0",
    "true",
    "yes",
    "active",
    "hit",
    "positive",
    "strict_positive",
}
INACTIVE_TEXT = {
    "0",
    "0.0",
    "false",
    "no",
    "inactive",
    "non-hit",
    "nonhit",
    "negative",
    "measured_or_reliable_negative",
}
INVALID_TEXT = (
    "outside typical range",
    "potential transcription error",
    "non standard unit",
)
DESCRIPTOR_COLUMNS = ("rdkit_mol_wt", "rdkit_mol_logp", "rdkit_tpsa")
SMILES_COLUMNS = ("canonical_smiles", "smiles", "ligand_smiles", "smiles_neutral")
INCHIKEY_COLUMNS = (
    "inchikey",
    "standard_inchi_key",
    "standard_inchikey",
    "ligand_inchikey",
)
NAME_COLUMNS = ("generic_name", "display_name", "drug_name", "mapped_drug_name")


@dataclass(frozen=True)
class SourceSpec:
    name: str
    path: Path


def parse_source_spec(value: str | Path) -> SourceSpec:
    """Parse a repeatable ``[SOURCE=]PATH`` source argument."""

    text = str(value).strip()
    if not text:
        raise ValueError("source path must not be empty")
    if "=" in text:
        name, raw_path = text.split("=", 1)
        if name.strip() and raw_path.strip():
            return SourceSpec(name=_canonical_source_name(name), path=Path(raw_path))
    path = Path(text)
    return SourceSpec(name=_canonical_source_name(path.stem), path=path)


def stage_matched_inactive_additions(
    *,
    spd_dataset_path: str | Path,
    source_paths: Sequence[str | Path | SourceSpec] | Mapping[str, str | Path],
    out_dir: str | Path,
    max_per_target: int = 10,
    label_col: str | None = None,
    active_nm: float = 1_000.0,
    inactive_nm: float = 10_000.0,
    min_tanimoto: float = 0.0,
    max_descriptor_distance: float | None = None,
    descriptor_weight: float = 0.25,
    include_existing: bool = False,
    library_name: str = "matched_inactive_additions",
    prepare: bool = False,
    force: bool = False,
    config_path: str | Path = "config.txt",
) -> dict[str, Any]:
    """Stage local, measured matched-inactive additions for an SPD table.

    Candidate compounds must resolve to an existing SPD ``ligand_base``. Source
    evidence is never downloaded, the input SPD table is never mutated, and
    this workflow does not launch docking.
    """

    _validate_thresholds(
        max_per_target=max_per_target,
        active_nm=active_nm,
        inactive_nm=inactive_nm,
        min_tanimoto=min_tanimoto,
        max_descriptor_distance=max_descriptor_distance,
        descriptor_weight=descriptor_weight,
    )
    dataset_path = Path(spd_dataset_path)
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    specs = _source_specs(source_paths)
    if not specs:
        raise ValueError("at least one local external source path is required")
    missing = [str(spec.path) for spec in specs if not spec.path.is_file()]
    if missing:
        raise FileNotFoundError(f"external source paths do not exist: {missing}")

    spd = read_source_table(dataset_path)
    selected_label = _select_label_column(spd, label_col)
    target_meta = _target_metadata(spd)
    compound_meta, compound_lookups = _compound_metadata(spd)
    descriptor_scales = _descriptor_scales(compound_meta)

    source_rows = _load_source_rows(
        specs,
        active_nm=active_nm,
        inactive_nm=inactive_nm,
    )
    source_rows = _resolve_targets(source_rows, target_meta)
    source_rows = _resolve_ligands(source_rows, compound_meta, compound_lookups)

    existing_pairs = _pair_keys(spd)
    labels = pd.to_numeric(spd[selected_label], errors="coerce")
    positive_pairs = _pair_keys(spd.loc[labels.eq(1)])
    positives = _positive_reference_rows(spd, labels, compound_meta)
    candidates = _collapse_candidate_pairs(
        source_rows,
        target_meta=target_meta,
        existing_pairs=existing_pairs,
        positive_pairs=positive_pairs,
        include_existing=include_existing,
    )
    candidates = _match_candidates(
        candidates,
        positives=positives,
        descriptor_scales=descriptor_scales,
        min_tanimoto=min_tanimoto,
        max_descriptor_distance=max_descriptor_distance,
        descriptor_weight=descriptor_weight,
    )
    candidates, selected_pairs = _select_candidates(
        candidates,
        target_meta=target_meta,
        max_per_target=max_per_target,
    )

    source_rows_path = output / "matched_inactive_source_row_classification.csv"
    all_candidates_path = output / "all_matched_inactive_candidates.csv"
    selected_path = output / "selected_matched_inactive_pairs.csv"
    unique_ligands_path = output / "selected_unique_ligands.csv"
    summary_path = output / "matched_inactive_per_target_source_summary.csv"
    sdf_path = output / f"{library_name}.sdf"
    sdf_failures_path = output / "sdf_generation_failures.csv"
    availability_path = output / "matched_inactive_preparation_availability.csv"
    forced_pdb_manifest_path = output / "forced_pdb_manifest.csv"

    prepared_selected, unique_ligands, preparation, sdf_failures, availability = (
        _preparation_outputs(
            selected_pairs,
            output=output,
            sdf_path=sdf_path,
            library_name=library_name,
            prepare=prepare,
            force=force,
            config_path=Path(config_path),
        )
    )
    summary = _per_target_source_summary(candidates, prepared_selected, source_rows)
    forced_pdb_manifest = _forced_pdb_manifest(prepared_selected)

    source_rows.to_csv(source_rows_path, index=False)
    candidates.to_csv(all_candidates_path, index=False)
    prepared_selected.to_csv(selected_path, index=False)
    unique_ligands.to_csv(unique_ligands_path, index=False)
    summary.to_csv(summary_path, index=False)
    sdf_failures.to_csv(sdf_failures_path, index=False)
    availability.to_csv(availability_path, index=False)
    forced_pdb_manifest.to_csv(forced_pdb_manifest_path, index=False)

    manifest: dict[str, Any] = {
        "spd_dataset": str(dataset_path),
        "out_dir": str(output),
        "label_col": selected_label,
        "source_paths": [{"name": spec.name, "path": str(spec.path)} for spec in specs],
        "max_per_target": int(max_per_target),
        "thresholds": {
            "active_nm": float(active_nm),
            "inactive_nm": float(inactive_nm),
            "min_tanimoto": float(min_tanimoto),
            "max_descriptor_distance": (
                float(max_descriptor_distance)
                if max_descriptor_distance is not None
                else None
            ),
            "descriptor_weight": float(descriptor_weight),
        },
        "include_existing": bool(include_existing),
        "library_name": library_name,
        "counts": {
            "spd_rows": int(len(spd)),
            "spd_targets": int(target_meta["target_uniprot"].nunique()),
            "spd_ligands": int(compound_meta["ligand_base"].nunique()),
            "spd_positive_pairs": int(len(positive_pairs)),
            "source_rows": int(len(source_rows)),
            "source_measured_inactive_rows": int(
                source_rows["evidence_inactive"].sum()
            ),
            "source_active_rows": int(source_rows["evidence_active"].sum()),
            "source_target_mapped_rows": int(
                source_rows["target_mapping_status"].eq("matched_spd_uniprot").sum()
            ),
            "source_ligand_resolved_rows": int(
                source_rows["ligand_resolution_status"].eq("resolved_spd_ligand").sum()
            ),
            "candidate_pairs": int(len(candidates)),
            "eligible_pairs": int(candidates["eligible"].sum())
            if not candidates.empty
            else 0,
            "selected_unique_pairs": int(
                prepared_selected["candidate_pair_key"].nunique()
            )
            if not prepared_selected.empty
            else 0,
            "selected_pair_pdb_rows": int(len(prepared_selected)),
            "selected_unique_ligands": int(len(unique_ligands)),
            "forced_pdb_manifest_rows": int(len(forced_pdb_manifest)),
        },
        "descriptor_standardization": descriptor_scales,
        "preparation": preparation,
        "outputs": {
            "source_row_classification": str(source_rows_path),
            "all_candidates": str(all_candidates_path),
            "selected_pairs": str(selected_path),
            "selected_unique_ligands": str(unique_ligands_path),
            "per_target_source_summary": str(summary_path),
            "sdf": str(sdf_path),
            "sdf_failures": str(sdf_failures_path),
            "preparation_availability": str(availability_path),
            "forced_pdb_manifest": str(forced_pdb_manifest_path),
            "forced_pdb_manifest_rows": int(len(forced_pdb_manifest)),
        },
        "policy": {
            "generic_measured_inactive": (
                "exact or lower-bound direct activity >= inactive_nm; upper bounds cannot establish inactivity"
            ),
            "generic_measured_active": (
                "exact or upper-bound direct activity <= active_nm, or usable explicit active call"
            ),
            "toxcast_inactive": (
                "explicit non-hit with available assay/chemical/overall usability checks passing"
            ),
            "active_conflict": "any usable active evidence excludes the resolved drug-target pair",
            "compound_domain": "existing SPD/FDA-RDK ligand_base identities only",
            "target_domain": "exact SPD target_uniprot mapping and SPD PDB expansion only",
            "existing_pairs_excluded_by_default": True,
            "nearest_positive": (
                "same-target SPD positive maximizing Morgan Tanimoto minus descriptor_weight times "
                "standardized MW/logP/TPSA distance"
            ),
            "input_spd_mutated": False,
            "network_downloads": False,
            "docking_launched": False,
        },
    }
    manifest_path = output / "matched_inactive_addition_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    return manifest


def _validate_thresholds(
    *,
    max_per_target: int,
    active_nm: float,
    inactive_nm: float,
    min_tanimoto: float,
    max_descriptor_distance: float | None,
    descriptor_weight: float,
) -> None:
    if max_per_target < 1:
        raise ValueError("max_per_target must be at least 1")
    if active_nm <= 0 or inactive_nm <= active_nm:
        raise ValueError("inactive_nm must be greater than a positive active_nm")
    if not 0.0 <= min_tanimoto <= 1.0:
        raise ValueError("min_tanimoto must be between 0 and 1")
    if max_descriptor_distance is not None and max_descriptor_distance < 0:
        raise ValueError("max_descriptor_distance must be nonnegative")
    if descriptor_weight < 0:
        raise ValueError("descriptor_weight must be nonnegative")


def _source_specs(
    source_paths: Sequence[str | Path | SourceSpec] | Mapping[str, str | Path],
) -> list[SourceSpec]:
    if isinstance(source_paths, Mapping):
        return [
            SourceSpec(_canonical_source_name(name), Path(path))
            for name, path in source_paths.items()
        ]
    specs: list[SourceSpec] = []
    for value in source_paths:
        specs.append(
            value if isinstance(value, SourceSpec) else parse_source_spec(value)
        )
    return specs


def _select_label_column(frame: pd.DataFrame, requested: str | None) -> str:
    if requested:
        if requested not in frame.columns:
            raise ValueError(f"SPD dataset lacks requested label column {requested!r}")
        return requested
    for column in DEFAULT_LABEL_CANDIDATES:
        if column in frame.columns:
            return column
    raise ValueError(
        "SPD dataset lacks a supported positive label column; pass label_col explicitly"
    )


def _target_metadata(frame: pd.DataFrame) -> pd.DataFrame:
    if "target_uniprot" not in frame.columns:
        raise ValueError(
            "SPD dataset must contain target_uniprot for target/PDB restriction"
        )
    columns = [
        column
        for column in (
            "target_uniprot",
            "target_gene",
            "gene_symbol",
            "target_id",
            "pdb_id",
            "target_family",
            "protein_class",
        )
        if column in frame.columns
    ]
    meta = frame[columns].copy()
    meta["target_uniprot"] = meta["target_uniprot"].map(_upper)
    if "target_gene" not in meta.columns:
        meta["target_gene"] = _first_nonempty_series(meta, ("gene_symbol", "target_id"))
    for column in ("target_gene", "pdb_id", "target_family", "protein_class"):
        if column not in meta.columns:
            meta[column] = ""
        meta[column] = meta[column].map(_clean)
    meta = meta[meta["target_uniprot"].ne("")].copy()
    if meta.empty:
        raise ValueError("SPD dataset has no nonempty target_uniprot values")
    meta = meta[
        ["target_uniprot", "target_gene", "pdb_id", "target_family", "protein_class"]
    ].drop_duplicates()
    return (
        meta.groupby(["target_uniprot", "pdb_id"], dropna=False, sort=True)
        .agg(
            target_gene=("target_gene", _first_nonempty_value),
            target_family=("target_family", _first_nonempty_value),
            protein_class=("protein_class", _first_nonempty_value),
        )
        .reset_index()
    )


def _compound_metadata(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, dict[str, set[str]]]]:
    if "ligand_base" not in frame.columns:
        raise ValueError(
            "SPD dataset must contain ligand_base for existing FDA/RDK resolution"
        )
    work = frame.copy()
    work["ligand_base"] = work["ligand_base"].map(_clean)
    work = work[work["ligand_base"].ne("")].copy()
    if work.empty:
        raise ValueError("SPD dataset has no nonempty ligand_base values")

    rows: list[dict[str, Any]] = []
    for ligand_base, group in work.groupby("ligand_base", sort=True):
        smiles, smiles_column = _first_group_value(group, SMILES_COLUMNS)
        inchikey, inchikey_column = _first_group_value(group, INCHIKEY_COLUMNS)
        descriptor = {
            column: _first_numeric(group, column) for column in DESCRIPTOR_COLUMNS
        }
        chemistry = _chemistry(smiles) if smiles else None
        for column in DESCRIPTOR_COLUMNS:
            if pd.isna(descriptor[column]) and chemistry is not None:
                descriptor[column] = chemistry[column]
        scaffold = _first_group_value(group, ("scaffold_key", "murcko_scaffold"))[0]
        if not scaffold and chemistry is not None:
            scaffold = str(chemistry["scaffold_key"])
        pdbqt_path, _ = _first_group_value(
            group,
            ("pdbqt_path", "ligand_pdbqt_path", "ligand_path", "path"),
        )
        rows.append(
            {
                "ligand_base": ligand_base,
                "display_name": _first_group_value(group, NAME_COLUMNS)[0]
                or ligand_base,
                "generic_name": _first_group_value(
                    group, ("generic_name", "display_name")
                )[0]
                or ligand_base,
                "smiles": smiles,
                "canonical_smiles": str(chemistry["canonical_smiles"])
                if chemistry is not None
                else smiles,
                "inchikey": inchikey,
                "scaffold_key": scaffold,
                "spd_structure_source_column": smiles_column,
                "spd_inchikey_source_column": inchikey_column,
                "pdbqt_path": pdbqt_path,
                **descriptor,
            }
        )
    compounds = pd.DataFrame(rows)
    lookups: dict[str, dict[str, set[str]]] = {
        "base": {},
        "inchikey": {},
        "smiles": {},
        "name": {},
    }
    for _, row in compounds.iterrows():
        base = _clean(row["ligand_base"])
        _add_lookup(lookups["base"], _norm(base), base)
        _add_lookup(
            lookups["inchikey"], _inchikey_connectivity(row.get("inchikey")), base
        )
        _add_lookup(
            lookups["smiles"], _canonical_smiles(row.get("canonical_smiles")), base
        )
        for field in ("display_name", "generic_name"):
            _add_lookup(lookups["name"], _norm(row.get(field)), base)
    return compounds, lookups


def _load_source_rows(
    specs: Sequence[SourceSpec],
    *,
    active_nm: float,
    inactive_nm: float,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for source_index, spec in enumerate(specs):
        raw = read_source_table(spec.path)
        normalized = _normalize_source(raw, spec=spec, source_index=source_index)
        if normalized.empty:
            classifications = pd.DataFrame(
                index=normalized.index,
                columns=[
                    "evidence_state",
                    "evidence_active",
                    "evidence_inactive",
                    "evidence_row_conflict",
                    "evidence_classification_basis",
                    "toxcast_usability_basis",
                    "active_nm_threshold",
                    "inactive_nm_threshold",
                ],
            )
        else:
            classifications = normalized.apply(
                lambda row: _classify_evidence(
                    row, active_nm=active_nm, inactive_nm=inactive_nm
                ),
                axis=1,
                result_type="expand",
            )
        frames.append(pd.concat([normalized, classifications], axis=1))
    return pd.concat(frames, ignore_index=True, sort=False)


def _normalize_source(
    raw: pd.DataFrame, *, spec: SourceSpec, source_index: int
) -> pd.DataFrame:
    out = pd.DataFrame(index=raw.index)
    stated_source = _first_nonempty_series(
        raw, ("source", "source_name", "dataset", "label_source")
    )
    out["source_family"] = stated_source.map(
        lambda value: _canonical_source_name(value or spec.name)
    )
    out["source_name"] = stated_source.where(stated_source.ne(""), spec.name)
    out["source_path"] = str(spec.path)
    out["source_row_id"] = [f"{source_index}:{index}" for index in raw.index]
    out["source_target_uniprot"] = _first_nonempty_series(
        raw,
        ("target_uniprot", "uniprot", "uniprot_id", "uniprot_accession", "target_id"),
    ).map(_upper)
    out["source_target_gene"] = _first_nonempty_series(
        raw,
        ("target_gene", "gene_symbol", "gene", "cardiac_gene_symbol"),
    )
    out["source_ligand_base"] = _first_nonempty_series(
        raw,
        ("ligand_base", "drug_id", "rdk_id", "compound_id"),
    )
    out["source_drug_name"] = _first_nonempty_series(
        raw, NAME_COLUMNS + ("ligand_name",)
    )
    out["source_inchikey"] = _first_nonempty_series(raw, INCHIKEY_COLUMNS)
    out["source_smiles"] = _first_nonempty_series(raw, SMILES_COLUMNS)
    out["activity_nM"] = _activity_nm(raw)
    out["activity_relation"] = _first_nonempty_series(
        raw,
        ("activity_relation", "standard_relation", "relation", "operator"),
    ).map(_normalize_relation)
    out["activity_type"] = _first_nonempty_series(
        raw,
        ("activity_type", "standard_type", "endpoint_type", "endpoint", "type"),
    )
    out["activity_outcome"] = _first_nonempty_series(
        raw,
        (
            "source_specific_activity_label",
            "activity_outcome",
            "activity_class",
            "toxcast_hit_call",
            "hit_call",
            "hitc",
            "label",
            "toxcast_hts_label",
        ),
    )
    out["toxcast_hit_call"] = _first_nonempty_series(
        raw,
        ("toxcast_hit_call", "hit_call", "hitc"),
    )
    out["assay_id"] = _first_nonempty_series(
        raw,
        ("assay_id", "assay_chembl_id", "bindingdb_reactant_set_id", "aeid"),
    )
    out["publication_year"] = _first_nonempty_series(
        raw,
        (
            "publication_year",
            "activity_publication_year",
            "document_year",
            "database_release_year",
        ),
    )
    out["document_ids"] = _first_nonempty_series(
        raw,
        (
            "document_ids",
            "activity_document_ids",
            "document_chembl_id",
            "pmid",
            "pubmed_id",
        ),
    )
    out["data_validity_comment"] = _first_nonempty_series(
        raw,
        (
            "data_validity_comment",
            "chembl_data_validity_comment",
            "standard_text_value",
        ),
    )
    out["potential_duplicate"] = _first_nonempty_series(
        raw,
        ("potential_duplicate", "chembl_potential_duplicate"),
    )
    for column, aliases in {
        "toxcast_label_usable": ("toxcast_label_usable", "label_usable"),
        "toxcast_assay_quality_pass": (
            "toxcast_assay_quality_pass",
            "assay_quality_pass",
        ),
        "toxcast_chemical_qc_pass": ("toxcast_chemical_qc_pass", "chemical_qc_pass"),
        "toxcast_cytotoxicity_confounded": (
            "toxcast_cytotoxicity_confounded",
            "cytotoxicity_confounded",
        ),
    }.items():
        out[column] = _first_nonempty_series(raw, aliases)
        out[f"{column}_available"] = any(alias in raw.columns for alias in aliases)
    out["source_provenance_json"] = [
        json.dumps(
            {str(key): value for key, value in row.items() if pd.notna(value)},
            sort_keys=True,
            default=str,
        )
        for row in raw.to_dict("records")
    ]
    return out


def _activity_nm(raw: pd.DataFrame) -> pd.Series:
    values = pd.to_numeric(
        _first_nonempty_series(
            raw,
            (
                "activity_nM",
                "activity_nm",
                "activity_value_nM",
                "standard_value_nm",
                "standard_value",
                "value_nM",
                "value",
                "ac50_nM",
                "ac50_nm",
            ),
        ),
        errors="coerce",
    )
    units = _first_nonempty_series(
        raw, ("activity_units", "standard_units", "units", "unit")
    ).str.lower()
    values = values.where(~units.isin({"um", "µm", "micromolar"}), values * 1_000.0)
    values = values.where(~units.isin({"mm", "millimolar"}), values * 1_000_000.0)
    values = values.where(~units.isin({"pm", "picomolar"}), values / 1_000.0)
    pchembl = pd.to_numeric(
        _first_nonempty_series(raw, ("pchembl_value", "pchembl", "pxc50")),
        errors="coerce",
    )
    return values.fillna((10.0 ** (-pchembl)) * 1_000_000_000.0)


def _classify_evidence(
    row: pd.Series,
    *,
    active_nm: float,
    inactive_nm: float,
) -> dict[str, Any]:
    relation = _normalize_relation(row.get("activity_relation"))
    endpoint = _clean(row.get("activity_type")).upper()
    endpoint_root = re.sub(r"[^A-Z0-9].*$", "", endpoint)
    outcome = _norm(row.get("activity_outcome"))
    toxcast_outcome = _norm(row.get("toxcast_hit_call"))
    value = pd.to_numeric(pd.Series([row.get("activity_nM")]), errors="coerce").iloc[0]
    duplicate = _truth_value(row.get("potential_duplicate")) is True
    validity = _norm(row.get("data_validity_comment"))
    invalid = duplicate or any(token in validity for token in INVALID_TEXT)
    exact = relation in {"", "="}
    upper_bound = relation in {"<", "<="}
    lower_bound = relation in {">", ">="}
    direct_endpoint = endpoint_root in DIRECT_ENDPOINTS or endpoint_root == ""
    explicit_active = outcome in ACTIVE_TEXT or toxcast_outcome in ACTIVE_TEXT
    explicit_inactive = outcome in INACTIVE_TEXT or toxcast_outcome in INACTIVE_TEXT
    is_toxcast = _canonical_source_name(row.get("source_family")) == "ToxCast"

    toxcast_usable = True
    toxcast_checks: list[str] = []
    if is_toxcast:
        for column, expected in (
            ("toxcast_label_usable", True),
            ("toxcast_assay_quality_pass", True),
            ("toxcast_chemical_qc_pass", True),
            ("toxcast_cytotoxicity_confounded", False),
        ):
            if bool(row.get(f"{column}_available", False)):
                observed = _truth_value(row.get(column))
                passed = observed is expected
                toxcast_usable &= passed
                toxcast_checks.append(f"{column}={'pass' if passed else 'fail'}")
        if not toxcast_checks:
            toxcast_checks.append("explicit_hit_call_no_qc_columns")

    numeric_active = bool(
        direct_endpoint
        and pd.notna(value)
        and (exact or upper_bound)
        and float(value) <= active_nm
    )
    numeric_inactive = bool(
        direct_endpoint
        and pd.notna(value)
        and (exact or lower_bound)
        and float(value) >= inactive_nm
    )
    if is_toxcast:
        active = explicit_active and toxcast_usable and not invalid
        inactive = explicit_inactive and toxcast_usable and not invalid
        basis = (
            "conflicting_toxcast_hit_calls"
            if active and inactive
            else "toxcast_usable_hit"
            if active
            else "toxcast_usable_non_hit"
            if inactive
            else "toxcast_unusable_or_unknown"
        )
    else:
        active = (
            numeric_active or (explicit_active and direct_endpoint)
        ) and not invalid
        inactive = numeric_inactive and not invalid
        if active and inactive:
            basis = "conflicting_active_inactive_row"
        elif inactive:
            basis = "censor_aware_measured_inactive"
        elif active:
            basis = "censor_aware_measured_or_explicit_active"
        elif explicit_inactive and not numeric_inactive:
            basis = "explicit_inactive_without_threshold_support"
        elif (
            pd.notna(value) and float(value) > active_nm and float(value) < inactive_nm
        ):
            basis = "gray_zone"
        elif upper_bound and pd.notna(value) and float(value) >= inactive_nm:
            basis = "upper_bound_cannot_establish_inactive"
        else:
            basis = "unknown_or_unusable"
    conflict = bool(active and inactive)
    state = (
        "conflicting_row"
        if conflict
        else "active"
        if active
        else "measured_inactive"
        if inactive
        else "unknown"
    )
    return {
        "evidence_state": state,
        "evidence_active": bool(active),
        "evidence_inactive": bool(inactive),
        "evidence_row_conflict": conflict,
        "evidence_classification_basis": basis,
        "toxcast_usability_basis": ";".join(toxcast_checks),
        "active_nm_threshold": float(active_nm),
        "inactive_nm_threshold": float(inactive_nm),
    }


def _resolve_targets(
    source_rows: pd.DataFrame, target_meta: pd.DataFrame
) -> pd.DataFrame:
    work = source_rows.copy()
    canonical = {
        _upper(value): _upper(value)
        for value in target_meta["target_uniprot"].dropna().unique()
        if _upper(value)
    }
    work["target_uniprot"] = work["source_target_uniprot"].map(
        lambda value: canonical.get(_upper(value), "")
    )
    work["target_mapping_status"] = work["target_uniprot"].map(
        lambda value: "matched_spd_uniprot" if value else "no_spd_uniprot_match"
    )
    target_summary = (
        target_meta.sort_values(["target_uniprot", "pdb_id"])
        .groupby("target_uniprot", sort=False)
        .agg(
            target_gene=("target_gene", _first_nonempty_value),
            target_family=("target_family", _first_nonempty_value),
            protein_class=("protein_class", _first_nonempty_value),
            pdb_ids=(
                "pdb_id",
                lambda values: ";".join(
                    sorted({_clean(value).upper() for value in values if _clean(value)})
                ),
            ),
        )
        .reset_index()
    )
    return work.merge(target_summary, on="target_uniprot", how="left")


def _resolve_ligands(
    source_rows: pd.DataFrame,
    compounds: pd.DataFrame,
    lookups: dict[str, dict[str, set[str]]],
) -> pd.DataFrame:
    if source_rows.empty:
        work = source_rows.copy()
        work["ligand_base"] = pd.Series(dtype="object")
        work["ligand_resolution_status"] = pd.Series(dtype="object")
        work["ligand_match_type"] = pd.Series(dtype="object")
        work["source_spd_inchikey_conflict"] = pd.Series(dtype=bool)
        return work
    resolutions = source_rows.apply(
        lambda row: _resolve_ligand_row(row, lookups),
        axis=1,
        result_type="expand",
    )
    work = pd.concat([source_rows, resolutions], axis=1)
    compound_columns = {
        column: f"spd_{column}" if column != "ligand_base" else column
        for column in compounds.columns
    }
    work = work.merge(
        compounds.rename(columns=compound_columns),
        on="ligand_base",
        how="left",
        validate="many_to_one",
    )
    source_key = work["source_inchikey"].map(_inchikey_connectivity)
    spd_key = work.get("spd_inchikey", pd.Series("", index=work.index)).map(
        _inchikey_connectivity
    )
    work["source_spd_inchikey_conflict"] = (
        source_key.ne("") & spd_key.ne("") & source_key.ne(spd_key)
    )
    work.loc[work["source_spd_inchikey_conflict"], "ligand_resolution_status"] = (
        "source_spd_identity_conflict"
    )
    return work


def _resolve_ligand_row(
    row: pd.Series,
    lookups: dict[str, dict[str, set[str]]],
) -> dict[str, str]:
    queries = (
        ("ligand_base", "base", _norm(row.get("source_ligand_base"))),
        ("inchikey", "inchikey", _inchikey_connectivity(row.get("source_inchikey"))),
        ("smiles", "smiles", _canonical_smiles(row.get("source_smiles"))),
        ("drug_name", "name", _norm(row.get("source_drug_name"))),
        ("drug_id_as_name", "name", _norm(row.get("source_ligand_base"))),
    )
    ambiguous = False
    for match_type, lookup_name, query in queries:
        if not query:
            continue
        matches = lookups[lookup_name].get(query, set())
        if len(matches) == 1:
            return {
                "ligand_base": next(iter(matches)),
                "ligand_resolution_status": "resolved_spd_ligand",
                "ligand_match_type": match_type,
            }
        ambiguous |= len(matches) > 1
    return {
        "ligand_base": "",
        "ligand_resolution_status": "ambiguous_spd_ligand"
        if ambiguous
        else "no_spd_ligand_match",
        "ligand_match_type": "",
    }


def _collapse_candidate_pairs(
    source_rows: pd.DataFrame,
    *,
    target_meta: pd.DataFrame,
    existing_pairs: set[str],
    positive_pairs: set[str],
    include_existing: bool,
) -> pd.DataFrame:
    mapped = source_rows[
        source_rows["target_mapping_status"].eq("matched_spd_uniprot")
        & source_rows["ligand_resolution_status"].eq("resolved_spd_ligand")
    ].copy()
    if mapped.empty:
        return _empty_candidates()
    mapped["candidate_pair_key"] = mapped.apply(
        lambda row: _pair_key(row.get("target_uniprot"), row.get("ligand_base")), axis=1
    )
    rows: list[dict[str, Any]] = []
    for pair_key, group in mapped.groupby("candidate_pair_key", sort=True):
        inactive = group[group["evidence_inactive"]].copy()
        if inactive.empty:
            continue
        inactive = inactive.sort_values(
            ["activity_nM", "source_family", "source_row_id"],
            ascending=[False, True, True],
            na_position="last",
        )
        representative = inactive.iloc[0]
        active_conflict = bool(
            group["evidence_active"].any() or group["evidence_row_conflict"].any()
        )
        already_present = pair_key in existing_pairs
        spd_positive_conflict = pair_key in positive_pairs
        target_uniprot = _clean(representative["target_uniprot"])
        pdb_ids = ";".join(
            sorted(
                {
                    _clean(value).upper()
                    for value in target_meta.loc[
                        target_meta["target_uniprot"].eq(target_uniprot), "pdb_id"
                    ]
                    if _clean(value)
                }
            )
        )
        reasons: list[str] = []
        if active_conflict:
            reasons.append("source_conflicted_active_pair")
        if spd_positive_conflict:
            reasons.append("spd_positive_pair_conflict")
        if already_present and not include_existing:
            reasons.append("pair_already_present")
        if not pdb_ids:
            reasons.append("no_spd_pdb")
        if not _clean(representative.get("spd_canonical_smiles")):
            reasons.append("missing_spd_structure")
        provenance = [json.loads(value) for value in group["source_provenance_json"]]
        rows.append(
            {
                "candidate_pair_key": pair_key,
                "target_uniprot": target_uniprot,
                "target_gene": _clean(representative.get("target_gene")),
                "target_family": _clean(representative.get("target_family")),
                "protein_class": _clean(representative.get("protein_class")),
                "pdb_ids": pdb_ids,
                "ligand_base": _clean(representative.get("ligand_base")),
                "display_name": _clean(representative.get("spd_display_name")),
                "generic_name": _clean(representative.get("spd_generic_name")),
                "smiles": _clean(representative.get("spd_canonical_smiles")),
                "inchikey": _clean(representative.get("spd_inchikey")),
                "scaffold_key": _clean(representative.get("spd_scaffold_key")),
                "rdkit_mol_wt": representative.get("spd_rdkit_mol_wt"),
                "rdkit_mol_logp": representative.get("spd_rdkit_mol_logp"),
                "rdkit_tpsa": representative.get("spd_rdkit_tpsa"),
                "spd_pdbqt_path": _clean(representative.get("spd_pdbqt_path")),
                "representative_source": _clean(representative.get("source_family")),
                "representative_source_name": _clean(representative.get("source_name")),
                "representative_source_path": _clean(representative.get("source_path")),
                "representative_source_row_id": _clean(
                    representative.get("source_row_id")
                ),
                "activity_nM": representative.get("activity_nM"),
                "activity_relation": _clean(representative.get("activity_relation")),
                "activity_type": _clean(representative.get("activity_type")),
                "assay_id": _clean(representative.get("assay_id")),
                "publication_year": _clean(representative.get("publication_year")),
                "document_ids": _clean(representative.get("document_ids")),
                "evidence_sources": ";".join(
                    sorted(set(group["source_family"].map(_clean)) - {""})
                ),
                "evidence_source_names": ";".join(
                    sorted(set(group["source_name"].map(_clean)) - {""})
                ),
                "evidence_source_paths": ";".join(
                    sorted(set(group["source_path"].map(_clean)) - {""})
                ),
                "evidence_source_row_ids": ";".join(
                    sorted(set(group["source_row_id"].map(_clean)) - {""})
                ),
                "n_source_evidence_rows": int(len(group)),
                "n_measured_inactive_rows": int(group["evidence_inactive"].sum()),
                "n_active_rows": int(group["evidence_active"].sum()),
                "n_conflicting_rows": int(group["evidence_row_conflict"].sum()),
                "source_conflicted_active_pair": active_conflict,
                "pair_already_present": already_present,
                "spd_positive_pair_conflict": spd_positive_conflict,
                "include_existing_policy": bool(include_existing),
                "evidence_provenance_json": json.dumps(
                    provenance, sort_keys=True, default=str
                ),
                "projected_label": 0,
                "projected_label_status": "strict_measured_matched_inactive",
                "negative_evidence_type": "censor_aware_measured_inactive",
                "eligible": not reasons,
                "exclusion_reasons": ";".join(reasons),
            }
        )
    return pd.DataFrame(rows) if rows else _empty_candidates()


def _positive_reference_rows(
    spd: pd.DataFrame,
    labels: pd.Series,
    compounds: pd.DataFrame,
) -> pd.DataFrame:
    positive = spd.loc[labels.eq(1)].copy()
    if positive.empty:
        return pd.DataFrame()
    positive["target_uniprot"] = positive["target_uniprot"].map(_upper)
    positive["ligand_base"] = positive["ligand_base"].map(_clean)
    positive = positive[
        positive["target_uniprot"].ne("") & positive["ligand_base"].ne("")
    ].copy()
    keep = [
        column
        for column in ("target_uniprot", "target_gene", "pdb_id", "ligand_base")
        if column in positive
    ]
    positive = positive[keep].drop_duplicates(["target_uniprot", "ligand_base"])
    return positive.merge(
        compounds, on="ligand_base", how="left", validate="many_to_one"
    )


def _match_candidates(
    candidates: pd.DataFrame,
    *,
    positives: pd.DataFrame,
    descriptor_scales: dict[str, float],
    min_tanimoto: float,
    max_descriptor_distance: float | None,
    descriptor_weight: float,
) -> pd.DataFrame:
    if candidates.empty:
        return candidates
    rows: list[dict[str, Any]] = []
    for _, candidate in candidates.iterrows():
        row = candidate.to_dict()
        target_positives = (
            positives[
                positives["target_uniprot"].eq(candidate["target_uniprot"])
                & positives["ligand_base"].ne(candidate["ligand_base"])
            ]
            if not positives.empty
            else pd.DataFrame()
        )
        matches: list[dict[str, Any]] = []
        candidate_chemistry = _chemistry(candidate.get("smiles"))
        if candidate_chemistry is not None:
            for _, positive in target_positives.iterrows():
                positive_chemistry = _chemistry(positive.get("canonical_smiles"))
                if positive_chemistry is None:
                    continue
                tanimoto = _tanimoto(
                    candidate_chemistry["fingerprint"],
                    positive_chemistry["fingerprint"],
                )
                distance = _standardized_distance(
                    candidate, positive, descriptor_scales
                )
                if math.isnan(distance):
                    continue
                score = float(tanimoto) - descriptor_weight * distance
                matches.append(
                    {
                        "nearest_positive_ligand_base": _clean(
                            positive.get("ligand_base")
                        ),
                        "nearest_positive_display_name": _clean(
                            positive.get("display_name")
                        ),
                        "nearest_positive_pdb_id": _clean(
                            positive.get("pdb_id")
                        ).upper(),
                        "nearest_positive_smiles": _clean(
                            positive.get("canonical_smiles")
                        ),
                        "nearest_positive_scaffold_key": _clean(
                            positive.get("scaffold_key")
                        ),
                        "nearest_positive_tanimoto": float(tanimoto),
                        "nearest_positive_descriptor_distance": float(distance),
                        "nearest_positive_match_score": float(score),
                    }
                )
        if matches:
            best = sorted(
                matches,
                key=lambda item: (
                    -item["nearest_positive_match_score"],
                    -item["nearest_positive_tanimoto"],
                    item["nearest_positive_descriptor_distance"],
                    item["nearest_positive_ligand_base"],
                ),
            )[0]
            row.update(best)
            candidate_scaffold = _clean(candidate.get("scaffold_key"))
            positive_scaffold = _clean(best["nearest_positive_scaffold_key"])
            row["scaffold_match"] = bool(
                candidate_scaffold
                and positive_scaffold
                and candidate_scaffold == positive_scaffold
            )
            reasons = _reason_list(row.get("exclusion_reasons"))
            if best["nearest_positive_tanimoto"] < min_tanimoto:
                reasons.append("below_min_tanimoto")
            if (
                max_descriptor_distance is not None
                and best["nearest_positive_descriptor_distance"]
                > max_descriptor_distance
            ):
                reasons.append("above_max_descriptor_distance")
        else:
            row.update(
                {
                    "nearest_positive_ligand_base": "",
                    "nearest_positive_display_name": "",
                    "nearest_positive_pdb_id": "",
                    "nearest_positive_smiles": "",
                    "nearest_positive_scaffold_key": "",
                    "nearest_positive_tanimoto": math.nan,
                    "nearest_positive_descriptor_distance": math.nan,
                    "nearest_positive_match_score": math.nan,
                    "scaffold_match": False,
                }
            )
            reasons = _reason_list(row.get("exclusion_reasons"))
            reasons.append("no_same_target_spd_positive_match")
        row["exclusion_reasons"] = ";".join(dict.fromkeys(reasons))
        row["eligible"] = not reasons
        rows.append(row)
    return pd.DataFrame(rows)


def _select_candidates(
    candidates: pd.DataFrame,
    *,
    target_meta: pd.DataFrame,
    max_per_target: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if candidates.empty:
        return candidates, pd.DataFrame()
    ranked = candidates.copy()
    ranked["selected_pair"] = False
    ranked["selection_rank_within_target"] = pd.Series(
        pd.NA, index=ranked.index, dtype="Int64"
    )
    eligible = ranked[ranked["eligible"]].copy()
    eligible = eligible.sort_values(
        [
            "target_uniprot",
            "nearest_positive_match_score",
            "nearest_positive_tanimoto",
            "nearest_positive_descriptor_distance",
            "n_measured_inactive_rows",
            "ligand_base",
        ],
        ascending=[True, False, False, True, False, True],
    )
    selected_indices: list[int] = []
    for _, group in eligible.groupby("target_uniprot", sort=True):
        chosen = group.head(max_per_target)
        selected_indices.extend(chosen.index.tolist())
        for rank, index in enumerate(chosen.index, start=1):
            ranked.loc[index, "selection_rank_within_target"] = rank
    ranked.loc[selected_indices, "selected_pair"] = True
    selected_unique = ranked.loc[selected_indices].copy()
    if selected_unique.empty:
        return ranked, pd.DataFrame()
    expanded = selected_unique.drop(
        columns=["pdb_ids", "target_gene", "target_family", "protein_class"]
    ).merge(
        target_meta,
        on="target_uniprot",
        how="inner",
        validate="many_to_many",
    )
    expanded["pdb_id"] = expanded["pdb_id"].map(lambda value: _clean(value).upper())
    expanded = expanded[expanded["pdb_id"].ne("")].copy()
    expanded["selection_reason"] = (
        "censor_aware_measured_inactive;nearest_same_target_spd_positive;"
        "morgan_tanimoto_and_standardized_physchem_match"
    )
    expanded["ligand_domain"] = "spd_existing_fda_rdk"
    expanded["training_domain"] = "fda_atlas_current_ligand"
    expanded["production_truth_allowed"] = True
    expanded["benchmark_only"] = False
    return ranked, expanded.sort_values(
        ["target_uniprot", "selection_rank_within_target", "pdb_id", "ligand_base"]
    ).reset_index(drop=True)


def _preparation_outputs(
    selected: pd.DataFrame,
    *,
    output: Path,
    sdf_path: Path,
    library_name: str,
    prepare: bool,
    force: bool,
    config_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any], pd.DataFrame, pd.DataFrame]:
    failures = pd.DataFrame(columns=["ligand_file_stem", "reason"])
    if selected.empty:
        selected = selected.reindex(
            columns=[
                "candidate_pair_key",
                "target_uniprot",
                "target_gene",
                "pdb_id",
                "ligand_base",
                "display_name",
                "smiles",
                "evidence_sources",
                "projected_label",
                "projected_label_status",
                "nearest_positive_ligand_base",
                "nearest_positive_tanimoto",
                "nearest_positive_descriptor_distance",
                "scaffold_match",
                "selection_rank_within_target",
            ]
        )
        unique = pd.DataFrame(
            columns=[
                "ligand_base",
                "ligand_file_stem",
                "display_name",
                "generic_name",
                "smiles",
                "inchikey",
                "scaffold_key",
                "target_uniprots",
                "pdb_ids",
                "evidence_sources",
                "n_selected_target_pairs",
                "sdf_available",
                "pdbqt_available",
                "evidence_provenance_json",
            ]
        )
        availability = pd.DataFrame(
            columns=[
                "pdb_id",
                "target_uniprot",
                "ligand_base",
                "ligand_file_stem",
                "sdf_available",
                "pdbqt_path",
                "pdbqt_available",
                "preparation_status",
            ]
        )
        return (
            selected,
            unique,
            {
                "requested": bool(prepare),
                "prepared": False,
                "status": "no_selected_pairs",
                "sdf_available": False,
                "pdbqt_available_rows": int(availability["pdbqt_available"].sum()),
            },
            failures,
            availability,
        )

    from analysis.ml.target_positive_addition_prep import (
        _finalize_selected,
        _prepare_libraries,
        _write_sdf,
    )

    finalized = _finalize_selected(selected, library_name=library_name)
    failures = _write_sdf(finalized, sdf_path)
    failed_stems = set(
        failures.get("ligand_file_stem", pd.Series(dtype=str)).astype(str)
    )
    availability = finalized[
        ["pdb_id", "target_uniprot", "target_gene", "ligand_base", "ligand_file_stem"]
    ].copy()
    availability["sdf_path"] = str(sdf_path)
    availability["sdf_available"] = sdf_path.is_file() & ~availability[
        "ligand_file_stem"
    ].isin(failed_stems)
    availability["pdbqt_path"] = finalized.get(
        "spd_pdbqt_path",
        pd.Series("", index=finalized.index),
    ).map(_clean)
    availability["pdbqt_available"] = availability["pdbqt_path"].map(_path_is_file)
    availability["preparation_status"] = availability["pdbqt_available"].map(
        {True: "existing_spd_pdbqt", False: "not_requested"}
    )
    prep_manifest: dict[str, Any] = {
        "requested": bool(prepare),
        "prepared": False,
        "status": "not_requested",
        "sdf_available": bool(sdf_path.is_file()),
        "sdf_failures": int(len(failures)),
        "pdbqt_available_rows": int(availability["pdbqt_available"].sum()),
    }
    if prepare and sdf_path.is_file():
        try:
            prep_result = _prepare_libraries(
                finalized,
                sdf_path=sdf_path,
                output=output,
                library_name=library_name,
                force=force,
                config_path=config_path,
            )
            prepared_availability = pd.read_csv(
                prep_result["availability_manifest"],
                low_memory=False,
            )
            availability = availability.drop(
                columns=["pdbqt_path", "pdbqt_available"]
            ).merge(
                prepared_availability[
                    [
                        "pdb_id",
                        "target_gene",
                        "ligand_base",
                        "ligand_file_stem",
                        "pdbqt_path",
                        "pdbqt_available",
                    ]
                ],
                on=["pdb_id", "target_gene", "ligand_base", "ligand_file_stem"],
                how="left",
                validate="one_to_one",
            )
            availability["pdbqt_available"] = (
                availability["pdbqt_available"].fillna(False).astype(bool)
            )
            availability["preparation_status"] = availability["pdbqt_available"].map(
                {True: "prepared", False: "prep_missing"}
            )
            prep_manifest = {"requested": True, "status": "completed", **prep_result}
        except (
            Exception
        ) as exc:  # pragma: no cover - depends on optional local Meeko installation.
            availability["preparation_status"] = availability["pdbqt_available"].map(
                {True: "prep_failed_existing_spd_pdbqt", False: "prep_failed"}
            )
            availability["preparation_error"] = str(exc)
            prep_manifest = {
                "requested": True,
                "prepared": False,
                "status": "failed",
                "error": str(exc),
                "sdf_available": bool(sdf_path.is_file()),
                "sdf_failures": int(len(failures)),
                "pdbqt_available_rows": int(availability["pdbqt_available"].sum()),
            }
    finalized = finalized.merge(
        availability[
            [
                "pdb_id",
                "target_gene",
                "ligand_base",
                "ligand_file_stem",
                "sdf_available",
                "pdbqt_path",
                "pdbqt_available",
                "preparation_status",
            ]
        ],
        on=["pdb_id", "target_gene", "ligand_base", "ligand_file_stem"],
        how="left",
        validate="one_to_one",
    )
    unique = _unique_selected_ligands(finalized)
    return finalized, unique, prep_manifest, failures, availability


def _unique_selected_ligands(selected: pd.DataFrame) -> pd.DataFrame:
    if selected.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for ligand_base, group in selected.groupby("ligand_base", sort=True):
        first = group.iloc[0]
        rows.append(
            {
                "ligand_base": ligand_base,
                "ligand_file_stem": _clean(first.get("ligand_file_stem")),
                "display_name": _clean(first.get("display_name")),
                "generic_name": _clean(first.get("generic_name")),
                "smiles": _clean(first.get("smiles")),
                "inchikey": _clean(first.get("inchikey")),
                "scaffold_key": _clean(first.get("scaffold_key")),
                "rdkit_mol_wt": first.get("rdkit_mol_wt"),
                "rdkit_mol_logp": first.get("rdkit_mol_logp"),
                "rdkit_tpsa": first.get("rdkit_tpsa"),
                "target_uniprots": ";".join(
                    sorted(set(group["target_uniprot"].map(_clean)) - {""})
                ),
                "pdb_ids": ";".join(sorted(set(group["pdb_id"].map(_clean)) - {""})),
                "evidence_sources": ";".join(
                    sorted(
                        {
                            source
                            for value in group["evidence_sources"].map(_clean)
                            for source in value.split(";")
                            if source
                        }
                    )
                ),
                "n_selected_target_pairs": int(group["candidate_pair_key"].nunique()),
                "sdf_available": bool(group["sdf_available"].fillna(False).any()),
                "pdbqt_available": bool(group["pdbqt_available"].fillna(False).any()),
                "evidence_provenance_json": json.dumps(
                    [
                        json.loads(value)
                        for value in group["evidence_provenance_json"].drop_duplicates()
                    ],
                    sort_keys=True,
                    default=str,
                ),
            }
        )
    return pd.DataFrame(rows)


def _forced_pdb_manifest(selected: pd.DataFrame) -> pd.DataFrame:
    columns = ["gene", "pdb_id"]
    if selected.empty:
        return pd.DataFrame(columns=columns)
    required = {"target_gene", "pdb_id"}
    if not required.issubset(selected.columns):
        raise ValueError("selected matched-inactive pairs lack target_gene or pdb_id")
    manifest = selected[["target_gene", "pdb_id"]].rename(
        columns={"target_gene": "gene"}
    )
    manifest["gene"] = manifest["gene"].map(_clean)
    manifest["pdb_id"] = manifest["pdb_id"].map(_upper)
    manifest = manifest[manifest["gene"].ne("") & manifest["pdb_id"].ne("")]
    return manifest.drop_duplicates(columns).sort_values(columns).reset_index(drop=True)


def _per_target_source_summary(
    candidates: pd.DataFrame,
    selected: pd.DataFrame,
    source_rows: pd.DataFrame,
) -> pd.DataFrame:
    if candidates.empty:
        return pd.DataFrame(
            columns=[
                "target_uniprot",
                "target_gene",
                "source_family",
                "n_source_rows",
                "n_candidate_pairs",
                "n_eligible_pairs",
                "n_selected_pairs",
            ]
        )
    selected_keys = set(
        selected.get("candidate_pair_key", pd.Series(dtype=str)).astype(str)
    )
    exploded = candidates.copy()
    exploded["source_family"] = exploded["evidence_sources"].str.split(";")
    exploded = exploded.explode("source_family")
    exploded["selected_for_source"] = exploded["candidate_pair_key"].isin(selected_keys)
    exploded["selected_tanimoto"] = pd.to_numeric(
        exploded["nearest_positive_tanimoto"],
        errors="coerce",
    ).where(exploded["selected_for_source"])
    source_counts = (
        source_rows[
            source_rows["target_mapping_status"].eq("matched_spd_uniprot")
            & source_rows["evidence_inactive"]
        ]
        .groupby(["target_uniprot", "source_family"], dropna=False)
        .size()
        .rename("n_source_rows")
    )
    summary = (
        exploded.groupby(
            ["target_uniprot", "target_gene", "source_family"],
            dropna=False,
        )
        .agg(
            n_candidate_pairs=("candidate_pair_key", "nunique"),
            n_eligible_pairs=("eligible", "sum"),
            n_selected_pairs=("selected_for_source", "sum"),
            n_active_conflict_pairs=("source_conflicted_active_pair", "sum"),
            n_existing_pair_exclusions=("pair_already_present", "sum"),
            mean_selected_tanimoto=("selected_tanimoto", "mean"),
        )
        .reset_index()
    )
    summary = summary.merge(
        source_counts.reset_index(),
        on=["target_uniprot", "source_family"],
        how="left",
    )
    summary["n_source_rows"] = summary["n_source_rows"].fillna(0).astype(int)
    return (
        summary[
            [
                "target_uniprot",
                "target_gene",
                "source_family",
                "n_source_rows",
                "n_candidate_pairs",
                "n_eligible_pairs",
                "n_selected_pairs",
                "n_active_conflict_pairs",
                "n_existing_pair_exclusions",
                "mean_selected_tanimoto",
            ]
        ]
        .sort_values(["target_uniprot", "source_family"])
        .reset_index(drop=True)
    )


def _descriptor_scales(compounds: pd.DataFrame) -> dict[str, float]:
    scales: dict[str, float] = {}
    for column in DESCRIPTOR_COLUMNS:
        values = pd.to_numeric(compounds[column], errors="coerce")
        scale = float(values.std(ddof=0)) if values.notna().any() else math.nan
        scales[column] = scale if math.isfinite(scale) and scale > 0 else 1.0
    return scales


def _standardized_distance(
    candidate: Mapping[str, Any] | pd.Series,
    positive: Mapping[str, Any] | pd.Series,
    scales: Mapping[str, float],
) -> float:
    squared: list[float] = []
    for column in DESCRIPTOR_COLUMNS:
        left = pd.to_numeric(pd.Series([candidate.get(column)]), errors="coerce").iloc[
            0
        ]
        right = pd.to_numeric(pd.Series([positive.get(column)]), errors="coerce").iloc[
            0
        ]
        if pd.isna(left) or pd.isna(right):
            return math.nan
        squared.append(((float(left) - float(right)) / float(scales[column])) ** 2)
    return math.sqrt(sum(squared) / len(squared))


@lru_cache(maxsize=100_000)
def _chemistry(smiles: Any) -> dict[str, Any] | None:
    text = _clean(smiles)
    if not text:
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import (
            Crippen,
            Descriptors,
            rdFingerprintGenerator,
            rdMolDescriptors,
        )
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError as exc:  # pragma: no cover - environment contract.
        raise RuntimeError(
            "RDKit is required for matched-inactive chemical matching"
        ) from exc
    molecule = Chem.MolFromSmiles(text)
    if molecule is None:
        return None
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mol_wt = getattr(Descriptors, "MolWt")
    mol_logp = getattr(Crippen, "MolLogP")
    return {
        "canonical_smiles": Chem.MolToSmiles(molecule, isomericSmiles=True),
        "fingerprint": generator.GetFingerprint(molecule),
        "scaffold_key": MurckoScaffold.MurckoScaffoldSmiles(
            mol=molecule,
            includeChirality=False,
        ),
        "rdkit_mol_wt": float(mol_wt(molecule)),
        "rdkit_mol_logp": float(mol_logp(molecule)),
        "rdkit_tpsa": float(rdMolDescriptors.CalcTPSA(molecule)),
    }


def _tanimoto(left: Any, right: Any) -> float:
    from rdkit import DataStructs

    return float(DataStructs.TanimotoSimilarity(left, right))


@lru_cache(maxsize=100_000)
def _canonical_smiles(value: Any) -> str:
    chemistry = _chemistry(value)
    return str(chemistry["canonical_smiles"]) if chemistry is not None else ""


def _pair_keys(frame: pd.DataFrame) -> set[str]:
    if not {"target_uniprot", "ligand_base"}.issubset(frame.columns):
        return set()
    return {
        _pair_key(target, ligand)
        for target, ligand in zip(
            frame["target_uniprot"], frame["ligand_base"], strict=False
        )
        if _upper(target) and _clean(ligand)
    }


def _pair_key(target_uniprot: Any, ligand_base: Any) -> str:
    return f"{_upper(target_uniprot)}||{_norm(ligand_base)}"


def _empty_candidates() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "candidate_pair_key",
            "target_uniprot",
            "target_gene",
            "pdb_ids",
            "ligand_base",
            "evidence_sources",
            "eligible",
            "exclusion_reasons",
        ]
    )


def _first_nonempty_series(frame: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    result = pd.Series("", index=frame.index, dtype="object")
    for column in columns:
        if column not in frame.columns:
            continue
        values = frame[column].map(_clean)
        result = result.where(result.ne(""), values)
    return result


def _first_group_value(frame: pd.DataFrame, columns: Sequence[str]) -> tuple[str, str]:
    for column in columns:
        if column not in frame.columns:
            continue
        values = frame[column].map(_clean)
        values = values[values.ne("")]
        if not values.empty:
            return values.iloc[0], column
    return "", ""


def _first_numeric(frame: pd.DataFrame, column: str) -> float:
    if column not in frame.columns:
        return math.nan
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    return float(values.iloc[0]) if not values.empty else math.nan


def _first_nonempty_value(values: pd.Series) -> str:
    cleaned = values.map(_clean)
    cleaned = cleaned[cleaned.ne("")]
    return cleaned.iloc[0] if not cleaned.empty else ""


def _add_lookup(lookup: dict[str, set[str]], key: str, ligand_base: str) -> None:
    if key:
        lookup.setdefault(key, set()).add(ligand_base)


def _reason_list(value: Any) -> list[str]:
    return [item for item in _clean(value).split(";") if item]


def _normalize_relation(value: Any) -> str:
    text = _clean(value).replace("≤", "<=").replace("≥", ">=")
    return "" if text.lower() in {"nan", "none", "exact"} else text


def _canonical_source_name(value: Any) -> str:
    text = _norm(value)
    if "toxcast" in text or "tox21" in text or "invitrodb" in text:
        return "ToxCast"
    if "bindingdb" in text:
        return "BindingDB"
    if "papyrus" in text:
        return "Papyrus"
    if "chembl" in text:
        return "ChEMBL"
    return _clean(value) or "external"


def _inchikey_connectivity(value: Any) -> str:
    text = _clean(value).upper()
    return text.split("-", 1)[0] if text else ""


def _truth_value(value: Any) -> bool | None:
    text = _norm(value)
    if text in {"1", "1.0", "true", "yes", "y", "pass", "passed"}:
        return True
    if text in {"0", "0.0", "false", "no", "n", "fail", "failed"}:
        return False
    return None


def _path_is_file(value: Any) -> bool:
    text = _clean(value)
    return bool(text and Path(text).is_file())


def _clean(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).strip().split())


def _norm(value: Any) -> str:
    return _clean(value).lower()


def _upper(value: Any) -> str:
    return _clean(value).upper()
