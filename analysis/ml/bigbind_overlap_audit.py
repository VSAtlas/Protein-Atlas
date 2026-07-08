from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

ATLAS_SMILES_COLS = ["canonical_smiles", "smiles", "ligand_smiles", "banana_smiles"]
ATLAS_TARGET_COLS = ["target_uniprot", "uniprot", "target_id"]
BIGBIND_SPLIT_FILES = {
    "train": "activities_train.csv",
    "test": "activities_test.csv",
    "val": "activities_val.csv",
    "sna_1_train": "activities_sna_1_train.csv",
    "sna_1_test": "activities_sna_1_test.csv",
}


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _upper(value: Any) -> str:
    return _clean(value).upper()


def _lower(value: Any) -> str:
    return _clean(value).lower()


def _first_nonempty_frame(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in cols:
        if col not in df.columns:
            continue
        values = df[col].map(_clean)
        out = out.where(out.astype(str).str.len().gt(0), values)
    return out


def _canonical_smiles_and_inchikey(smiles: Any) -> tuple[str, str, str]:
    text = _clean(smiles)
    if not text:
        return "", "", "missing_smiles"
    try:
        from rdkit import Chem
    except Exception:
        return text, "", "rdkit_unavailable"
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return text, "", "rdkit_parse_failed"
    canonical = Chem.MolToSmiles(mol, isomericSmiles=False)
    inchikey = Chem.MolToInchiKey(mol)
    return canonical, inchikey, "ok"


def _add_ligand_keys(df: pd.DataFrame, smiles_col: str, *, canonicalize: bool = True) -> pd.DataFrame:
    out = df.copy()
    cleaned = out[smiles_col].map(_clean)
    if not canonicalize:
        out["_canonical_smiles_key"] = cleaned
        out["_inchikey_key"] = cleaned
        out["_ligand_key_status"] = "raw_smiles_fast_mode"
        return out
    unique = pd.Series(cleaned.dropna().unique(), dtype="object")
    cache = {value: _canonical_smiles_and_inchikey(value) for value in unique.tolist()}
    out["_canonical_smiles_key"] = cleaned.map(lambda value: cache.get(value, (value, "", "missing_smiles"))[0])
    out["_inchikey_key"] = cleaned.map(lambda value: cache.get(value, (value, "", "missing_smiles"))[1])
    out["_ligand_key_status"] = cleaned.map(lambda value: cache.get(value, (value, "", "missing_smiles"))[2])
    return out


def _split_files(bigbind_dir: Path) -> list[tuple[str, Path]]:
    root = bigbind_dir
    if not any((root / name).exists() for name in BIGBIND_SPLIT_FILES.values()):
        for child_name in ["BigBindV1", "BigBindV1.5"]:
            child = bigbind_dir / child_name
            if any((child / name).exists() for name in BIGBIND_SPLIT_FILES.values()):
                root = child
                break
    found: list[tuple[str, Path]] = []
    for split, filename in BIGBIND_SPLIT_FILES.items():
        path = root / filename
        if path.exists():
            found.append((split, path))
    return found


def _load_bigbind_splits(bigbind_dir: Path, *, canonicalize_ligands: bool = True) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for split, path in _split_files(bigbind_dir):
        header = pd.read_csv(path, nrows=0).columns.tolist()
        usecols = [col for col in ["lig_smiles", "uniprot", "pocket", "ex_rec_pdb", "active"] if col in header]
        if "lig_smiles" not in usecols:
            continue
        frame = pd.read_csv(path, usecols=usecols, low_memory=False)
        frame["bigbind_split"] = split
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    combined = _add_ligand_keys(
        combined.rename(columns={"lig_smiles": "_bigbind_smiles"}),
        "_bigbind_smiles",
        canonicalize=canonicalize_ligands,
    )
    combined["_target_key"] = combined.get("uniprot", pd.Series("", index=combined.index)).map(_upper)
    combined["_pdb_key"] = combined.get("ex_rec_pdb", pd.Series("", index=combined.index)).map(_lower)
    combined["_ligand_target_key"] = combined["_inchikey_key"] + "|" + combined["_target_key"]
    combined["_ligand_pdb_key"] = combined["_inchikey_key"] + "|" + combined["_pdb_key"]
    return combined


def _load_atlas(atlas_table: Path) -> pd.DataFrame:
    atlas = pd.read_csv(atlas_table, low_memory=False)
    atlas = atlas.copy()
    atlas["_atlas_smiles"] = _first_nonempty_frame(atlas, ATLAS_SMILES_COLS)
    atlas = _add_ligand_keys(atlas, "_atlas_smiles")
    atlas["_target_key"] = _first_nonempty_frame(atlas, ATLAS_TARGET_COLS).map(_upper)
    atlas["_pdb_key"] = atlas.get("pdb_id", pd.Series("", index=atlas.index)).map(_lower)
    atlas["_ligand_target_key"] = atlas["_inchikey_key"] + "|" + atlas["_target_key"]
    atlas["_ligand_pdb_key"] = atlas["_inchikey_key"] + "|" + atlas["_pdb_key"]
    return atlas


def _join_splits(values: set[str]) -> str:
    return ";".join(sorted(v for v in values if v))


def _split_lookup(bigbind: pd.DataFrame, key_col: str) -> dict[str, set[str]]:
    work = bigbind[[key_col, "bigbind_split"]].dropna().copy()
    work = work[work[key_col].astype(str).str.len().gt(0)]
    lookup: dict[str, set[str]] = {}
    for key, split in work.itertuples(index=False):
        lookup.setdefault(str(key), set()).add(str(split))
    return lookup


def _mark_membership(atlas: pd.DataFrame, bigbind: pd.DataFrame) -> pd.DataFrame:
    out = atlas.copy()
    for name, key_col in {
        "ligand": "_inchikey_key",
        "target": "_target_key",
        "pdb": "_pdb_key",
        "ligand_target": "_ligand_target_key",
        "ligand_pdb": "_ligand_pdb_key",
    }.items():
        lookup = _split_lookup(bigbind, key_col)
        out[f"bigbind_{name}_overlap"] = out[key_col].map(lambda key: bool(key and str(key) in lookup))
        out[f"bigbind_{name}_splits"] = out[key_col].map(lambda key: _join_splits(lookup.get(str(key), set())) if key else "")
    return out


def _summary(atlas: pd.DataFrame, bigbind: pd.DataFrame, membership: pd.DataFrame) -> dict[str, Any]:
    split_rows = []
    for split, group in bigbind.groupby("bigbind_split", dropna=False):
        split_rows.append(
            {
                "split": str(split),
                "rows": int(len(group)),
                "unique_ligands_inchikey": int(group["_inchikey_key"].replace("", pd.NA).nunique(dropna=True)),
                "unique_targets": int(group["_target_key"].replace("", pd.NA).nunique(dropna=True)),
                "unique_pdbs": int(group["_pdb_key"].replace("", pd.NA).nunique(dropna=True)),
            }
        )
    overlap_counts = {
        name: int(membership[f"bigbind_{name}_overlap"].sum())
        for name in ["ligand", "target", "pdb", "ligand_target", "ligand_pdb"]
        if f"bigbind_{name}_overlap" in membership.columns
    }
    return {
        "atlas_rows": int(len(atlas)),
        "atlas_unique_ligands_inchikey": int(atlas["_inchikey_key"].replace("", pd.NA).nunique(dropna=True)),
        "atlas_unique_targets": int(atlas["_target_key"].replace("", pd.NA).nunique(dropna=True)),
        "atlas_unique_pdbs": int(atlas["_pdb_key"].replace("", pd.NA).nunique(dropna=True)),
        "bigbind_rows": int(len(bigbind)),
        "bigbind_ligand_key_mode": "canonical_inchikey" if bool(bigbind.get("_ligand_key_status", pd.Series()).ne("raw_smiles_fast_mode").any()) else "raw_smiles_fast_mode",
        "bigbind_splits": split_rows,
        "atlas_overlap_row_counts": overlap_counts,
        "policy": "PBAS/BANANA independence checks should flag exact ligand, target, PDB, ligand-target, and ligand-PDB overlaps before claiming BigBind-independent evaluation.",
    }


def audit_bigbind_overlap(
    atlas_table: str | Path,
    bigbind_dir: str | Path,
    out_dir: str | Path,
    *,
    canonicalize_bigbind_ligands: bool = True,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    atlas = _load_atlas(Path(atlas_table))
    bigbind = _load_bigbind_splits(Path(bigbind_dir), canonicalize_ligands=canonicalize_bigbind_ligands)
    if bigbind.empty:
        raise ValueError(f"no BigBind activity split files found under {bigbind_dir}")
    membership = _mark_membership(atlas, bigbind)
    summary = _summary(atlas, bigbind, membership)
    public_cols = [
        col
        for col in [
            "drug_id",
            "target_id",
            "target_uniprot",
            "pdb_id",
            "_inchikey_key",
            "_target_key",
            "_pdb_key",
            "bigbind_ligand_overlap",
            "bigbind_ligand_splits",
            "bigbind_target_overlap",
            "bigbind_target_splits",
            "bigbind_pdb_overlap",
            "bigbind_pdb_splits",
            "bigbind_ligand_target_overlap",
            "bigbind_ligand_target_splits",
            "bigbind_ligand_pdb_overlap",
            "bigbind_ligand_pdb_splits",
        ]
        if col in membership.columns
    ]
    membership[public_cols].to_csv(out / "bigbind_atlas_overlap_rows.csv", index=False)
    pd.DataFrame(summary["bigbind_splits"]).to_csv(out / "bigbind_split_summary.csv", index=False)
    (out / "bigbind_overlap_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary
