from __future__ import annotations

from typing import Tuple

import pandas as pd

from analysis.ml.audit_utils import split_source_tokens


CHEMICAL_CLUSTER_FIELDS = [
    "chemical_cluster",
    "chemical_cluster_id",
    "ecfp_cluster",
    "butina_cluster",
    "umap_cluster",
    "scaffold_cluster",
    "scaffold_key",
    "ligand_chemotype",
]

TARGET_FAMILY_FIELDS = [
    "target_family",
    "protein_family",
    "target_class",
    "protein_class",
]

SOURCE_HOLDOUT_FIELDS = [
    "label_source",
    "source_family",
    "upstream_source",
]


def _entity_holdout(df: pd.DataFrame, field: str, seed: int, test_fraction: float) -> tuple[pd.Index, pd.Index]:
    entities = pd.Series(df[field].dropna().unique()).sample(frac=1.0, random_state=seed).tolist()
    n_test = max(1, int(round(len(entities) * test_fraction))) if entities else 0
    held_out = set(entities[:n_test])
    test_mask = df[field].isin(held_out)
    return df.index[~test_mask], df.index[test_mask]


def _first_available_field(df: pd.DataFrame, fields: list[str]) -> str | None:
    for field in fields:
        if field in df.columns and df[field].notna().any():
            return field
    return None


def _murcko_scaffold(smiles: str) -> str:
    if not str(smiles or "").strip():
        return ""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        mol = Chem.MolFromSmiles(str(smiles).strip())
        if mol is None:
            return ""
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        return scaffold or Chem.MolToSmiles(mol, isomericSmiles=False)
    except Exception:
        return ""


def _derived_chemical_cluster(df: pd.DataFrame) -> pd.Series:
    smiles_field = _first_available_field(df, ["smiles", "canonical_smiles", "ligand_smiles"])
    if smiles_field is None:
        raise ValueError(
            "chemical_cluster_holdout requires one of "
            f"{CHEMICAL_CLUSTER_FIELDS} or a SMILES column"
        )
    scaffolds = df[smiles_field].map(_murcko_scaffold).fillna("")
    if scaffolds.str.len().gt(0).any():
        return "murcko|" + scaffolds
    raise ValueError(
        "chemical_cluster_holdout could not derive Murcko scaffolds from SMILES; "
        "provide chemical_cluster, butina_cluster, scaffold_key, or ligand_chemotype"
    )


def make_split(df: pd.DataFrame, split_mode: str = "random", seed: int = 42, test_fraction: float = 0.2) -> Tuple[pd.Index, pd.Index]:
    if split_mode == "random":
        test = df.sample(frac=test_fraction, random_state=seed).index
        return df.index.difference(test), test
    if split_mode == "drug_holdout":
        return _entity_holdout(df, "drug_id", seed, test_fraction)
    if split_mode == "target_holdout":
        return _entity_holdout(df, "target_id", seed, test_fraction)
    if split_mode == "chemical_cluster_holdout":
        field = _first_available_field(df, CHEMICAL_CLUSTER_FIELDS)
        work = df
        if field is None:
            work = df.copy()
            field = "_derived_chemical_cluster"
            work[field] = _derived_chemical_cluster(df)
        return _entity_holdout(work, field, seed, test_fraction)
    if split_mode == "target_family_holdout":
        field = _first_available_field(df, TARGET_FAMILY_FIELDS)
        if field is None:
            raise ValueError(f"target_family_holdout requires one of {TARGET_FAMILY_FIELDS}")
        return _entity_holdout(df, field, seed, test_fraction)
    if split_mode == "source_holdout":
        field = _first_available_field(df, SOURCE_HOLDOUT_FIELDS)
        if field is None:
            raise ValueError(f"source_holdout requires one of {SOURCE_HOLDOUT_FIELDS}")
        exploded = df[field].map(split_source_tokens)
        sources = sorted({source for tokens in exploded for source in tokens})
        if len(sources) < 2:
            raise ValueError(f"source_holdout requires at least two distinct {field} values")
        held_out = set(
            pd.Series(sources)
            .sample(frac=1.0, random_state=seed)
            .tolist()[: max(1, int(round(len(sources) * test_fraction)))]
        )
        test_mask = exploded.map(lambda tokens: bool(set(tokens) & held_out))
        train_mask = ~test_mask
        if not train_mask.any() or not test_mask.any():
            raise ValueError(
                f"source_holdout produced empty train/test for {field}: "
                f"train={int(train_mask.sum())} test={int(test_mask.sum())}"
            )
        return df.index[train_mask], df.index[test_mask]
    if split_mode in {"scaffold_holdout", "protein_class_holdout"}:
        field = (
            "scaffold_key"
            if split_mode == "scaffold_holdout" and "scaffold_key" in df.columns
            else "ligand_chemotype"
            if split_mode == "scaffold_holdout"
            else "protein_class"
        )
        if field not in df.columns:
            raise ValueError(f"{split_mode} requires {field}")
        return _entity_holdout(df, field, seed, test_fraction)
    if split_mode == "temporal_holdout":
        for field in ("label_publication_year", "evidence_publication_year", "activity_publication_year", "database_release_year"):
            if field in df.columns:
                years = pd.to_numeric(df[field], errors="coerce")
                observed = years.notna()
                missing_fraction = float((~observed).mean()) if len(years) else 1.0
                if missing_fraction > 0.50:
                    raise ValueError(
                        f"temporal_holdout year coverage too low for {field}: "
                        f"missing_fraction={missing_fraction:.3f} > 0.500"
                    )
                cutoff = years.quantile(0.8)
                return df.index[years <= cutoff], df.index[years > cutoff]
        raise ValueError("temporal_holdout requires an evidence/date year column")
    raise ValueError(f"unsupported split_mode: {split_mode}")


def split_overlap_summary(
    train: pd.DataFrame,
    test: pd.DataFrame,
    split_mode: str,
) -> dict[str, object]:
    fields_by_mode = {
        "drug_holdout": ["drug_id"],
        "target_holdout": ["target_id"],
        "scaffold_holdout": ["scaffold_key" if "scaffold_key" in train.columns and "scaffold_key" in test.columns else "ligand_chemotype"],
        "protein_class_holdout": ["protein_class"],
        "chemical_cluster_holdout": [
            _first_available_field(train, CHEMICAL_CLUSTER_FIELDS)
            or _first_available_field(test, CHEMICAL_CLUSTER_FIELDS)
            or "scaffold_key"
        ],
        "target_family_holdout": [
            _first_available_field(train, TARGET_FAMILY_FIELDS)
            or _first_available_field(test, TARGET_FAMILY_FIELDS)
            or "protein_class"
        ],
        "source_holdout": [
            _first_available_field(train, SOURCE_HOLDOUT_FIELDS)
            or _first_available_field(test, SOURCE_HOLDOUT_FIELDS)
            or "label_source"
        ],
    }
    fields = fields_by_mode.get(split_mode, [])
    overlaps: dict[str, int] = {}
    for field in fields:
        if field not in train.columns or field not in test.columns:
            overlaps[field] = -1
            continue
        overlaps[field] = len(set(train[field].dropna()) & set(test[field].dropna()))
    return {
        "split_mode": split_mode,
        "n_train": len(train),
        "n_test": len(test),
        "overlaps": overlaps,
        "passes_holdout": all(value == 0 for value in overlaps.values()),
    }
