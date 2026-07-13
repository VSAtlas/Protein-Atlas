from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class ChemicalClusterSummary:
    method: str
    input_rows: int
    unique_smiles: int
    valid_smiles: int
    invalid_smiles: int
    n_clusters: int
    distance_threshold: float

    def as_dict(self) -> dict[str, Any]:
        return {
            "chemical_cluster_method": self.method,
            "chemical_cluster_input_rows": self.input_rows,
            "chemical_cluster_unique_smiles": self.unique_smiles,
            "chemical_cluster_valid_smiles": self.valid_smiles,
            "chemical_cluster_invalid_smiles": self.invalid_smiles,
            "chemical_cluster_n_clusters": self.n_clusters,
            "chemical_cluster_distance_threshold": self.distance_threshold,
        }


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _first_smiles_column(df: pd.DataFrame) -> str | None:
    for field in ("smiles", "canonical_smiles", "ligand_smiles"):
        if field in df.columns and df[field].notna().any():
            return field
    return None


def add_butina_chemical_clusters(
    df: pd.DataFrame,
    *,
    distance_threshold: float = 0.35,
    radius: int = 2,
    n_bits: int = 2048,
    max_unique_smiles: int = 20000,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add ECFP/Butina chemical clusters for harder chemical holdouts.

    The split code can already consume a `chemical_cluster` column. This helper
    creates one from molecular fingerprints so publication stress tests are not
    limited to Murcko scaffold labels. The threshold is Tanimoto distance, so
    0.35 roughly groups molecules with similarity >= 0.65.
    """

    out = df.copy()
    smiles_col = _first_smiles_column(out)
    if smiles_col is None:
        out["butina_cluster"] = "butina_unassigned"
        out["chemical_cluster"] = "butina_unassigned"
        return out, ChemicalClusterSummary(
            method="butina",
            input_rows=len(out),
            unique_smiles=0,
            valid_smiles=0,
            invalid_smiles=0,
            n_clusters=0,
            distance_threshold=distance_threshold,
        ).as_dict()

    unique_smiles = sorted({_clean(value) for value in out[smiles_col].dropna() if _clean(value)})
    if len(unique_smiles) > max_unique_smiles:
        raise ValueError(
            "Butina clustering is O(n^2) in the number of unique SMILES; "
            f"got {len(unique_smiles)} unique molecules, limit is {max_unique_smiles}. "
            "Use scaffold/auto clustering or precompute clusters for this table."
        )

    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import rdFingerprintGenerator
        from rdkit.ML.Cluster import Butina
    except Exception as exc:  # pragma: no cover - depends on optional RDKit env
        raise RuntimeError("RDKit is required for --chemical-cluster butina") from exc

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    valid_smiles: list[str] = []
    fps = []
    for smiles in unique_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        valid_smiles.append(smiles)
        fps.append(generator.GetFingerprint(mol))

    distances: list[float] = []
    for i in range(1, len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        distances.extend([1.0 - sim for sim in sims])
    clusters = Butina.ClusterData(distances, len(fps), distance_threshold, isDistData=True)

    lookup: dict[str, str] = {}
    for cluster_index, member_indices in enumerate(clusters):
        cluster_id = f"butina:{cluster_index:05d}"
        for member_index in member_indices:
            lookup[valid_smiles[member_index]] = cluster_id

    cluster_values = out[smiles_col].map(lambda value: lookup.get(_clean(value), "butina_unassigned"))
    out["butina_cluster"] = cluster_values
    out["chemical_cluster"] = cluster_values
    out["chemical_cluster_source"] = "butina_ecfp"
    summary = ChemicalClusterSummary(
        method="butina",
        input_rows=len(out),
        unique_smiles=len(unique_smiles),
        valid_smiles=len(valid_smiles),
        invalid_smiles=len(unique_smiles) - len(valid_smiles),
        n_clusters=len(clusters),
        distance_threshold=distance_threshold,
    )
    return out, summary.as_dict()
