from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from analysis.ml.feature_sets import PHYSICHEM_DESCRIPTOR_FEATURES
from analysis.ml.labels import binary_label_series


# These resolutions are deliberately coarser than RDKit's floating-point output.
# A signature that remains unique here is usable as a practical identity proxy.
DESCRIPTOR_RESOLUTIONS = {
    "rdkit_mol_wt": 0.1,
    "rdkit_mol_logp": 0.01,
    "rdkit_tpsa": 0.1,
    "rdkit_hbd": 1.0,
    "rdkit_hba": 1.0,
    "rdkit_rotatable_bonds": 1.0,
    "rdkit_formal_charge": 1.0,
    "rdkit_aromatic_rings": 1.0,
    "rdkit_fraction_csp3": 0.01,
    "rdkit_qed": 0.01,
}


def _stable_signature(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    text = frame[columns].astype("string").fillna("<NA>").agg("|".join, axis=1)
    return text.map(lambda value: hashlib.sha256(value.encode("utf-8")).hexdigest()[:16])


def _rounded_descriptors(frame: pd.DataFrame, descriptors: list[str]) -> pd.DataFrame:
    rounded = pd.DataFrame(index=frame.index)
    for descriptor in descriptors:
        values = pd.to_numeric(frame[descriptor], errors="coerce")
        resolution = DESCRIPTOR_RESOLUTIONS[descriptor]
        rounded[descriptor] = (values / resolution).round().astype("Int64")
    return rounded


def _signature_capacity(
    frame: pd.DataFrame,
    *,
    signature: pd.Series,
    signature_type: str,
    descriptors: str,
) -> dict[str, Any]:
    counts = signature.value_counts(dropna=False)
    n_drugs = int(len(frame))
    singleton_drugs = int(counts.eq(1).sum())
    return {
        "signature_type": signature_type,
        "descriptors": descriptors,
        "n_descriptors": len(descriptors.split("|")),
        "n_drugs": n_drugs,
        "n_unique_signatures": int(len(counts)),
        "unique_signature_fraction": float(len(counts) / n_drugs),
        "singleton_drug_fraction": float(singleton_drugs / n_drugs),
        "largest_collision_group": int(counts.max()),
        "median_drugs_per_signature": float(counts.median()),
    }


def _drug_level_table(
    frame: pd.DataFrame,
    *,
    label_col: str,
    descriptors: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    work = frame.copy()
    work[label_col] = binary_label_series(work[label_col])
    work = work.loc[work[label_col].notna()].copy()
    for descriptor in descriptors:
        work[descriptor] = pd.to_numeric(work[descriptor], errors="coerce")
    rounded = _rounded_descriptors(work, descriptors)
    work["_practical_signature"] = _stable_signature(rounded, descriptors)

    identity_cols = [
        col
        for col in (
            "display_name",
            "generic_name",
            "canonical_smiles",
            "smiles",
            "inchikey",
            "scaffold_key",
            "chemical_cluster",
        )
        if col in work.columns
    ]
    invariance_rows: list[dict[str, Any]] = []
    for descriptor in descriptors:
        counts = work.groupby("drug_id", dropna=False)[descriptor].nunique(dropna=False)
        invariance_rows.append(
            {
                "descriptor": descriptor,
                "n_drugs": int(len(counts)),
                "n_drugs_with_multiple_values": int(counts.gt(1).sum()),
                "max_values_within_drug": int(counts.max()),
            }
        )

    signature_counts = (
        work.groupby(["drug_id", "_practical_signature"], dropna=False)
        .size()
        .rename("n_rows_for_signature")
        .reset_index()
    )
    modal = (
        signature_counts.sort_values(
            ["drug_id", "n_rows_for_signature", "_practical_signature"],
            ascending=[True, False, True],
            na_position="last",
        )
        .drop_duplicates("drug_id")
        .rename(columns={"_practical_signature": "modal_practical_signature"})
    )
    modal_rows = work.merge(
        modal[["drug_id", "modal_practical_signature"]],
        on="drug_id",
        how="inner",
        validate="many_to_one",
    )
    modal_rows = modal_rows.loc[
        modal_rows["_practical_signature"].eq(
            modal_rows["modal_practical_signature"]
        )
    ]
    descriptor_profile = modal_rows.groupby("drug_id", dropna=False).agg(
        {
            **{descriptor: "median" for descriptor in descriptors},
            **{column: "first" for column in identity_cols},
            "modal_practical_signature": "first",
        }
    )
    label_profile = work.groupby("drug_id", dropna=False)[label_col].agg(
        ["sum", "count", "mean"]
    )
    label_profile.columns = ["n_positive", "n_labeled_rows", "positive_rate"]
    drug = label_profile.join(descriptor_profile, how="left")

    conflict_drugs = signature_counts.groupby("drug_id", dropna=False).size()
    conflict_ids = conflict_drugs.loc[conflict_drugs.gt(1)].index
    conflicts = work.loc[work["drug_id"].isin(conflict_ids)].copy()
    if conflicts.empty:
        conflict_detail = pd.DataFrame()
    else:
        conflict_aggregations: dict[str, Any] = {
            label_col: ["sum", "count", "mean"],
            **{descriptor: "median" for descriptor in descriptors},
            **{
                column: lambda values: "|".join(
                    sorted(set(values.dropna().astype(str)))[:10]
                )
                for column in identity_cols
            },
        }
        conflict_detail = conflicts.groupby(
            ["drug_id", "_practical_signature"], dropna=False
        ).agg(conflict_aggregations)
        conflict_detail.columns = [
            "n_positive",
            "n_rows",
            "positive_rate",
            *descriptors,
            *identity_cols,
        ]
        conflict_detail = conflict_detail.reset_index().rename(
            columns={"_practical_signature": "practical_signature"}
        )
    return drug.reset_index(), pd.DataFrame(invariance_rows), conflict_detail


def _same_nonmissing(left: Any, right: Any) -> bool:
    return bool(pd.notna(left) and pd.notna(right) and left == right)


def _nearest_neighbors(
    drug: pd.DataFrame,
    *,
    descriptors: list[str],
) -> pd.DataFrame:
    numeric = drug[descriptors].apply(pd.to_numeric, errors="coerce")
    complete = numeric.notna().all(axis=1)
    eligible = drug.loc[complete].reset_index(drop=True)
    values = numeric.loc[complete].to_numpy(dtype=float)
    scale = values.std(axis=0, ddof=0)
    scale[~np.isfinite(scale) | (scale == 0.0)] = 1.0
    standardized = (values - values.mean(axis=0)) / scale
    distances = np.sqrt(
        np.square(standardized[:, None, :] - standardized[None, :, :]).sum(axis=2)
    )
    np.fill_diagonal(distances, np.inf)

    rows: list[dict[str, Any]] = []
    cluster_values = (
        eligible["chemical_cluster"].astype("string")
        if "chemical_cluster" in eligible.columns
        else pd.Series(pd.NA, index=eligible.index, dtype="string")
    )
    scaffold_values = (
        eligible["scaffold_key"].astype("string")
        if "scaffold_key" in eligible.columns
        else pd.Series(pd.NA, index=eligible.index, dtype="string")
    )
    for index, row in eligible.iterrows():
        candidates = [
            ("nearest_any", np.ones(len(eligible), dtype=bool)),
            (
                "nearest_outside_cluster",
                cluster_values.ne(cluster_values.iloc[index]).fillna(True).to_numpy(),
            ),
            (
                "nearest_outside_scaffold",
                scaffold_values.ne(scaffold_values.iloc[index]).fillna(True).to_numpy(),
            ),
        ]
        for neighbor_type, mask in candidates:
            mask[index] = False
            candidate_indices = np.flatnonzero(mask)
            if not len(candidate_indices):
                continue
            neighbor_index = int(candidate_indices[np.argmin(distances[index, candidate_indices])])
            neighbor = eligible.iloc[neighbor_index]
            rows.append(
                {
                    "drug_id": row["drug_id"],
                    "drug_name": row.get("generic_name", row.get("display_name")),
                    "neighbor_type": neighbor_type,
                    "neighbor_drug_id": neighbor["drug_id"],
                    "neighbor_drug_name": neighbor.get(
                        "generic_name", neighbor.get("display_name")
                    ),
                    "standardized_euclidean_distance": float(
                        distances[index, neighbor_index]
                    ),
                    "positive_rate": float(row["positive_rate"]),
                    "neighbor_positive_rate": float(neighbor["positive_rate"]),
                    "absolute_positive_rate_difference": float(
                        abs(row["positive_rate"] - neighbor["positive_rate"])
                    ),
                    "same_chemical_cluster": _same_nonmissing(
                        cluster_values.iloc[index],
                        cluster_values.iloc[neighbor_index],
                    ),
                    "same_scaffold": _same_nonmissing(
                        scaffold_values.iloc[index],
                        scaffold_values.iloc[neighbor_index],
                    ),
                }
            )
    return pd.DataFrame(rows)


def run_descriptor_constellation_audit(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    label_col: str = "spd_binding_label",
) -> dict[str, Any]:
    dataset = Path(dataset_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    header = pd.read_csv(dataset, nrows=0)
    descriptors = list(PHYSICHEM_DESCRIPTOR_FEATURES)
    required = {"drug_id", label_col, *descriptors}
    missing = sorted(required - set(header.columns))
    if missing:
        raise ValueError(f"missing required columns: {', '.join(missing)}")
    optional = [
        col
        for col in (
            "display_name",
            "generic_name",
            "canonical_smiles",
            "smiles",
            "inchikey",
            "scaffold_key",
            "chemical_cluster",
        )
        if col in header.columns
    ]
    frame = pd.read_csv(
        dataset,
        usecols=["drug_id", label_col, *descriptors, *optional],
        low_memory=False,
    )
    drug, invariance, within_drug_conflicts = _drug_level_table(
        frame,
        label_col=label_col,
        descriptors=descriptors,
    )
    rounded = _rounded_descriptors(drug, descriptors)
    exact = drug[descriptors].apply(pd.to_numeric, errors="coerce").round(12)

    capacities: list[dict[str, Any]] = []
    for signature_type, values in (("exact_12dp", exact), ("practical_resolution", rounded)):
        full_signature = _stable_signature(values, descriptors)
        capacities.append(
            _signature_capacity(
                drug,
                signature=full_signature,
                signature_type=signature_type,
                descriptors="|".join(descriptors),
            )
        )
        for first, second in combinations(descriptors, 2):
            pair = [first, second]
            capacities.append(
                _signature_capacity(
                    drug,
                    signature=_stable_signature(values, pair),
                    signature_type=signature_type,
                    descriptors="|".join(pair),
                )
            )

    practical_signature = _stable_signature(rounded, descriptors)
    profile = drug.copy()
    profile["practical_constellation_signature"] = practical_signature
    profile["practical_signature_drug_count"] = practical_signature.map(
        practical_signature.value_counts()
    )
    collisions = profile.loc[profile["practical_signature_drug_count"].gt(1)].copy()
    neighbors = _nearest_neighbors(drug, descriptors=descriptors)

    capacity = pd.DataFrame(capacities).sort_values(
        ["n_descriptors", "singleton_drug_fraction", "unique_signature_fraction"],
        ascending=[False, False, False],
    )
    profile.to_csv(out / "descriptor_drug_profiles.csv", index=False)
    invariance.to_csv(out / "descriptor_within_drug_invariance.csv", index=False)
    within_drug_conflicts.to_csv(
        out / "descriptor_within_drug_conflicts.csv", index=False
    )
    capacity.to_csv(out / "descriptor_identity_capacity.csv", index=False)
    collisions.to_csv(out / "descriptor_constellation_collisions.csv", index=False)
    neighbors.to_csv(out / "descriptor_nearest_neighbors.csv", index=False)

    full_practical = capacity.loc[
        capacity["signature_type"].eq("practical_resolution")
        & capacity["n_descriptors"].eq(len(descriptors))
    ].iloc[0]
    nearest_any = neighbors.loc[neighbors["neighbor_type"].eq("nearest_any")]
    conflict_counts = (
        within_drug_conflicts.groupby("drug_id", dropna=False).size()
        if not within_drug_conflicts.empty
        else pd.Series(dtype=int)
    )
    manifest = {
        "dataset": str(dataset),
        "label_col": label_col,
        "n_labeled_rows": int(drug["n_labeled_rows"].sum()),
        "n_drugs": int(len(drug)),
        "n_drugs_with_multiple_practical_signatures": int(len(conflict_counts)),
        "max_practical_signatures_per_drug": (
            int(conflict_counts.max()) if len(conflict_counts) else 1
        ),
        "status": "failed" if len(conflict_counts) else "passed",
        "descriptors": descriptors,
        "resolutions": DESCRIPTOR_RESOLUTIONS,
        "practical_full_signature_unique_fraction": float(
            full_practical["unique_signature_fraction"]
        ),
        "practical_full_signature_singleton_fraction": float(
            full_practical["singleton_drug_fraction"]
        ),
        "nearest_neighbor_same_cluster_fraction": float(
            nearest_any["same_chemical_cluster"].mean()
        ),
        "nearest_neighbor_positive_rate_spearman": float(
            nearest_any["positive_rate"].corr(
                nearest_any["neighbor_positive_rate"], method="spearman"
            )
        ),
        "interpretation": (
            "Signature uniqueness measures identity capacity, not causal feature use. "
            "Use it with grouped OOF transfer gaps, SHAP interactions, and nearest-neighbor "
            "failure review before attributing errors to descriptor memorization."
        ),
        "outputs": {
            "drug_profiles": str(out / "descriptor_drug_profiles.csv"),
            "invariance": str(out / "descriptor_within_drug_invariance.csv"),
            "within_drug_conflicts": str(
                out / "descriptor_within_drug_conflicts.csv"
            ),
            "identity_capacity": str(out / "descriptor_identity_capacity.csv"),
            "collisions": str(out / "descriptor_constellation_collisions.csv"),
            "nearest_neighbors": str(out / "descriptor_nearest_neighbors.csv"),
        },
    }
    (out / "descriptor_constellation_audit_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
