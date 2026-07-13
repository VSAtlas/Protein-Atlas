from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from analysis.external.source_tables import read_source_table
from analysis.ml.chemical_clusters import add_butina_chemical_clusters
from analysis.ml.target_positive_addition_prep import (
    _finalize_selected,
    _prepare_libraries,
    _resolve_structures,
    _write_sdf,
)

DIRECT_ENDPOINTS = {"IC50", "KI", "KD"}
SAFE_INACTIVE_RELATIONS = {"", "=", ">", ">="}
SAFE_ACTIVE_RELATIONS = {"", "=", "<", "<="}


def stage_ttd_strict_nonbinders(
    *,
    ttd_path: str | Path,
    positive_pairs_path: str | Path,
    out_dir: str | Path,
    library_name: str = "ttd_strict_nonbinders",
    nonbinder_nm: float = 200_000.0,
    active_nm: float = 1_000.0,
    negatives_per_positive: int = 2,
    max_per_target: int = 10,
    resolve_pubchem: bool = True,
    prepare: bool = False,
    force: bool = False,
    config_path: str | Path = "config.txt",
) -> dict[str, Any]:
    """Stage strict TTD nonbinders matched to locally selected TTD positives.

    TTD defines nonbinders above 200 uM. Rows are retained only when a direct
    Ki/Kd/IC50 value safely proves that bound and no same-target TTD active
    measurement conflicts. These controls remain external sensitivity rows.
    """

    if nonbinder_nm <= active_nm:
        raise ValueError("nonbinder_nm must be greater than active_nm")
    if negatives_per_positive < 1 or max_per_target < 1:
        raise ValueError("selection counts must be positive")

    source_path = Path(ttd_path)
    positive_path = Path(positive_pairs_path)
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)

    source = read_source_table(source_path)
    positives = pd.read_csv(positive_path, low_memory=False)
    _require_columns(
        source,
        {
            "drug_id",
            "gene_symbol",
            "activity_nM",
            "activity_relation",
            "activity_type",
            "pubchem_cid",
        },
        source_path,
    )
    _require_columns(
        positives,
        {"source", "target_gene", "target_uniprot", "pdb_id", "smiles"},
        positive_path,
    )

    classified = _classify_ttd_rows(
        source,
        active_nm=active_nm,
        nonbinder_nm=nonbinder_nm,
    )
    ttd_positives = positives[
        positives["source"].fillna("").astype(str).str.casefold().eq("ttd")
    ].copy()
    target_meta = (
        ttd_positives[
            [
                "target_gene",
                "target_uniprot",
                "pdb_id",
                "target_family",
                "protein_class",
            ]
        ]
        .drop_duplicates("target_gene")
        .copy()
    )
    positive_counts = ttd_positives.groupby("target_gene").size().to_dict()
    positive_ids = set(
        zip(
            ttd_positives["target_gene"].map(_clean),
            ttd_positives.get("drug_id", pd.Series("", index=ttd_positives.index)).map(_clean),
            strict=False,
        )
    )

    candidates = _collapse_candidates(
        classified,
        target_meta=target_meta,
        positive_ids=positive_ids,
        source_path=source_path,
    )
    candidates, resolution = _resolve_structures(
        candidates,
        cache_path=output / "pubchem_cid_properties.csv",
        fetch=resolve_pubchem,
    )
    candidates = candidates[
        candidates["structure_resolution_status"].isin(
            {"source_structure", "pubchem_resolved"}
        )
    ].copy()
    candidates = _standardize_parent_fragments(candidates)
    candidates = _add_nearest_positive_similarity(candidates, ttd_positives)
    candidates, cluster_summary = add_butina_chemical_clusters(
        candidates,
        distance_threshold=0.35,
    )
    selected = _select_diverse_controls(
        candidates,
        positive_counts=positive_counts,
        negatives_per_positive=negatives_per_positive,
        max_per_target=max_per_target,
    )
    selected = _finalize_selected(selected, library_name=library_name)

    classified_path = output / "ttd_row_classification.csv"
    candidates_path = output / "ttd_strict_nonbinder_candidates.csv"
    selected_path = output / "selected_ttd_strict_nonbinders.csv"
    summary_path = output / "ttd_strict_nonbinder_target_summary.csv"
    sdf_path = output / f"{library_name}.sdf"
    failures_path = output / "sdf_generation_failures.csv"

    failures = _write_sdf(selected, sdf_path)
    preparation: dict[str, Any] = {
        "requested": bool(prepare),
        "prepared": False,
        "status": "not_requested",
    }
    if prepare and not selected.empty and sdf_path.is_file():
        preparation = {
            "requested": True,
            "status": "completed",
            **_prepare_libraries(
                selected,
                sdf_path=sdf_path,
                output=output,
                library_name=library_name,
                force=force,
                config_path=Path(config_path),
            ),
        }
        availability = pd.read_csv(
            preparation["availability_manifest"],
            low_memory=False,
        )
        selected = selected.merge(
            availability,
            on=["pdb_id", "target_gene", "ligand_base", "ligand_file_stem"],
            how="left",
            validate="one_to_one",
        )

    summary = _target_summary(
        classified,
        candidates,
        selected,
        positive_counts=positive_counts,
        negatives_per_positive=negatives_per_positive,
        max_per_target=max_per_target,
    )
    classified.to_csv(classified_path, index=False)
    candidates.to_csv(candidates_path, index=False)
    selected.to_csv(selected_path, index=False)
    summary.to_csv(summary_path, index=False)
    failures.to_csv(failures_path, index=False)

    manifest = {
        "ttd_path": str(source_path),
        "positive_pairs_path": str(positive_path),
        "out_dir": str(output),
        "library_name": library_name,
        "thresholds": {
            "active_nm": float(active_nm),
            "strict_nonbinder_nm": float(nonbinder_nm),
        },
        "selection": {
            "negatives_per_positive": int(negatives_per_positive),
            "max_per_target": int(max_per_target),
        },
        "counts": {
            "ttd_rows": int(len(classified)),
            "strict_nonbinder_rows": int(classified["strict_nonbinder"].sum()),
            "active_rows": int(classified["strict_active"].sum()),
            "conflicting_pair_rows": int(classified["pair_has_active_conflict"].sum()),
            "resolved_candidates": int(len(candidates)),
            "multifragment_candidates_parent_standardized": int(
                candidates["structure_standardization_status"]
                .eq("largest_fragment_parent")
                .sum()
            ),
            "selected_pairs": int(len(selected)),
            "selected_unique_ligands": int(
                selected.get("ligand_file_stem", pd.Series(dtype=str)).nunique()
            ),
        },
        "pubchem_resolution": resolution,
        "chemical_clustering": cluster_summary,
        "preparation": preparation,
        "outputs": {
            "row_classification": str(classified_path),
            "candidates": str(candidates_path),
            "selected_pairs": str(selected_path),
            "target_summary": str(summary_path),
            "sdf": str(sdf_path),
            "sdf_failures": str(failures_path),
        },
        "policy": {
            "source_definition": "TTD strict nonbinder >200 uM",
            "comparable_active_conflicts_excluded": True,
            "absence_is_negative": False,
            "spd_binding_label_modified": False,
            "spd_exposure_label_modified": False,
            "production_fda_truth": False,
            "training_domain": "external_probe_sensitivity",
            "multifragment_policy": (
                "retain raw_smiles and use deterministic largest-heavy-atom "
                "parent fragment for docking"
            ),
        },
    }
    manifest_path = output / "ttd_strict_nonbinder_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    return manifest


def _classify_ttd_rows(
    frame: pd.DataFrame,
    *,
    active_nm: float,
    nonbinder_nm: float,
) -> pd.DataFrame:
    work = frame.copy()
    work["target_gene"] = work["gene_symbol"].map(_clean).str.upper()
    work["drug_id"] = work["drug_id"].map(_clean)
    work["activity_nM"] = pd.to_numeric(work["activity_nM"], errors="coerce")
    work["activity_relation"] = work["activity_relation"].map(_clean)
    work["activity_type"] = work["activity_type"].map(_clean).str.upper()
    direct = work["activity_type"].isin(DIRECT_ENDPOINTS)
    work["strict_nonbinder"] = (
        direct
        & work["activity_relation"].isin(SAFE_INACTIVE_RELATIONS)
        & work["activity_nM"].gt(nonbinder_nm)
    )
    work["strict_active"] = (
        direct
        & work["activity_relation"].isin(SAFE_ACTIVE_RELATIONS)
        & work["activity_nM"].le(active_nm)
    )
    work["source_pair_key"] = work["target_gene"] + "|" + work["drug_id"]
    conflict_keys = set(work.loc[work["strict_active"], "source_pair_key"])
    work["pair_has_active_conflict"] = work["source_pair_key"].isin(conflict_keys)
    return work


def _collapse_candidates(
    classified: pd.DataFrame,
    *,
    target_meta: pd.DataFrame,
    positive_ids: set[tuple[str, str]],
    source_path: Path,
) -> pd.DataFrame:
    strict = classified[
        classified["strict_nonbinder"]
        & ~classified["pair_has_active_conflict"]
    ].copy()
    strict = strict[
        strict["target_gene"].isin(set(target_meta["target_gene"].map(_clean)))
    ]
    strict = strict[
        ~strict.apply(
            lambda row: (_clean(row["target_gene"]), _clean(row["drug_id"]))
            in positive_ids,
            axis=1,
        )
    ]
    strict = strict.sort_values(
        ["target_gene", "drug_id", "activity_nM"],
        ascending=[True, True, False],
    ).drop_duplicates(["target_gene", "drug_id"], keep="first")
    strict = strict.merge(
        target_meta,
        on="target_gene",
        how="inner",
        validate="many_to_one",
    )
    strict["source"] = "TTD"
    strict["source_path"] = str(source_path)
    strict["source_row_id"] = strict.get("assay_id", strict["source_pair_key"])
    strict["activity_uM"] = strict["activity_nM"] / 1000.0
    strict["drug_name"] = strict.get("drug_name", pd.Series("", index=strict.index))
    strict["display_name"] = strict["drug_name"]
    strict["generic_name"] = strict["drug_name"]
    strict["smiles"] = strict.get("smiles", pd.Series("", index=strict.index))
    strict["inchikey"] = strict.get("inchikey", pd.Series("", index=strict.index))
    strict["ligand_base"] = strict["drug_id"]
    strict["projected_label"] = 0
    strict["projected_label_status"] = "ttd_strict_nonbinder_gt200um"
    strict["negative_evidence_type"] = "ttd_strict_nonbinder_gt200um"
    strict["training_allowed"] = True
    strict["production_truth_allowed"] = False
    strict["benchmark_only"] = False
    strict["ligand_domain"] = "external_probe_sensitivity"
    strict["training_domain"] = "external_probe_sensitivity"
    strict["selection_reason"] = (
        "same_source_ttd_strict_nonbinder;active_conflict_excluded"
    )
    return strict.reset_index(drop=True)


def _add_nearest_positive_similarity(
    candidates: pd.DataFrame,
    positives: pd.DataFrame,
) -> pd.DataFrame:
    work = candidates.copy()
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    positive_fingerprints: dict[str, list[tuple[str, Any]]] = {}
    for _, row in positives.iterrows():
        mol = Chem.MolFromSmiles(_clean(row.get("smiles")))
        if mol is None:
            continue
        positive_fingerprints.setdefault(_clean(row.get("target_gene")), []).append(
            (_clean(row.get("ligand_base")), generator.GetFingerprint(mol))
        )

    nearest_ids: list[str] = []
    similarities: list[float] = []
    for _, row in work.iterrows():
        mol = Chem.MolFromSmiles(_clean(row.get("smiles")))
        references = positive_fingerprints.get(_clean(row.get("target_gene")), [])
        if mol is None or not references:
            nearest_ids.append("")
            similarities.append(float("nan"))
            continue
        fingerprint = generator.GetFingerprint(mol)
        scored = [
            (float(DataStructs.TanimotoSimilarity(fingerprint, reference)), ligand)
            for ligand, reference in references
        ]
        similarity, ligand = max(scored)
        nearest_ids.append(ligand)
        similarities.append(similarity)
    work["nearest_ttd_positive_ligand_base"] = nearest_ids
    work["nearest_ttd_positive_tanimoto"] = similarities
    return work


def _select_diverse_controls(
    candidates: pd.DataFrame,
    *,
    positive_counts: dict[str, int],
    negatives_per_positive: int,
    max_per_target: int,
) -> pd.DataFrame:
    selected_indices: list[int] = []
    eligible = candidates[
        candidates.get(
            "structure_parent_valid",
            pd.Series(True, index=candidates.index),
        )
    ].copy()
    for target, group in eligible.groupby("target_gene", sort=True):
        requested = min(
            max_per_target,
            max(1, int(positive_counts.get(target, 0)) * negatives_per_positive),
        )
        ranked = group.sort_values(
            ["nearest_ttd_positive_tanimoto", "activity_nM", "drug_id"],
            ascending=[False, False, True],
            na_position="last",
        )
        diverse = ranked.drop_duplicates("chemical_cluster", keep="first")
        chosen = list(diverse.head(requested).index)
        if len(chosen) < requested:
            chosen.extend(
                index
                for index in ranked.index
                if index not in chosen
            )
        selected_indices.extend(chosen[:requested])
    selected = candidates.loc[selected_indices].copy()
    selected["selection_rank_within_target"] = (
        selected.groupby("target_gene").cumcount() + 1
    )
    selected["requested_controls_for_target"] = selected["target_gene"].map(
        lambda target: min(
            max_per_target,
            max(1, int(positive_counts.get(target, 0)) * negatives_per_positive),
        )
    )
    return selected.reset_index(drop=True)


def _target_summary(
    classified: pd.DataFrame,
    candidates: pd.DataFrame,
    selected: pd.DataFrame,
    *,
    positive_counts: dict[str, int],
    negatives_per_positive: int,
    max_per_target: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    targets = sorted(positive_counts)
    for target in targets:
        requested = min(
            max_per_target,
            max(1, int(positive_counts[target]) * negatives_per_positive),
        )
        rows.append(
            {
                "target_gene": target,
                "ttd_positive_rows": int(positive_counts[target]),
                "strict_nonbinder_source_rows": int(
                    (
                        classified["target_gene"].eq(target)
                        & classified["strict_nonbinder"]
                    ).sum()
                ),
                "active_conflict_source_rows": int(
                    (
                        classified["target_gene"].eq(target)
                        & classified["pair_has_active_conflict"]
                    ).sum()
                ),
                "resolved_candidate_pairs": int(
                    candidates["target_gene"].eq(target).sum()
                ),
                "requested_controls": requested,
                "selected_controls": int(
                    selected.get(
                        "target_gene",
                        pd.Series("", index=selected.index),
                    ).eq(target).sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def _standardize_parent_fragments(frame: pd.DataFrame) -> pd.DataFrame:
    work = frame.copy()
    raw_smiles: list[str] = []
    parent_smiles: list[str] = []
    fragment_counts: list[int] = []
    parent_fractions: list[float] = []
    statuses: list[str] = []
    valid: list[bool] = []
    for value in work["smiles"]:
        raw = _clean(value)
        molecule = Chem.MolFromSmiles(raw) if raw else None
        raw_smiles.append(raw)
        if molecule is None:
            parent_smiles.append("")
            fragment_counts.append(0)
            parent_fractions.append(float("nan"))
            statuses.append("invalid_structure")
            valid.append(False)
            continue
        fragments = Chem.GetMolFrags(
            molecule,
            asMols=True,
            sanitizeFrags=True,
        )
        parent = max(
            fragments,
            key=lambda fragment: (
                fragment.GetNumHeavyAtoms(),
                fragment.GetNumAtoms(),
                Chem.MolToSmiles(fragment, isomericSmiles=True),
            ),
        )
        total_heavy = max(1, molecule.GetNumHeavyAtoms())
        parent_smiles.append(Chem.MolToSmiles(parent, isomericSmiles=True))
        fragment_counts.append(len(fragments))
        parent_fractions.append(parent.GetNumHeavyAtoms() / total_heavy)
        statuses.append(
            "single_fragment" if len(fragments) == 1 else "largest_fragment_parent"
        )
        valid.append(parent.GetNumHeavyAtoms() > 0)
    work["raw_smiles"] = raw_smiles
    work["raw_structure_fragment_count"] = fragment_counts
    work["parent_heavy_atom_fraction"] = parent_fractions
    work["structure_standardization_status"] = statuses
    work["structure_parent_valid"] = valid
    work["smiles"] = parent_smiles
    return work


def _require_columns(
    frame: pd.DataFrame,
    required: set[str],
    path: Path,
) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")


def _clean(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()
