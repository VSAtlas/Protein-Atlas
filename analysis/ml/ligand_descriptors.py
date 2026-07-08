from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


DESCRIPTOR_COLUMNS = [
    "rdkit_mol_wt",
    "rdkit_mol_logp",
    "rdkit_tpsa",
    "rdkit_hbd",
    "rdkit_hba",
    "rdkit_rotatable_bonds",
    "rdkit_formal_charge",
    "rdkit_aromatic_rings",
    "rdkit_fraction_csp3",
    "rdkit_qed",
]
DESCRIPTOR_SOURCE_COL = "rdkit_descriptor_source"
DESCRIPTOR_STATUS_COL = "rdkit_descriptor_status"
SMILES_FIELDS = ["smiles", "canonical_smiles", "ligand_smiles", "smiles_neutral", "isomeric_smiles"]
FDA_MAPPING_PATH = Path("chemdb/data/fda_mapping_from_pdbqt.csv")
NAME_FIELDS = [
    "drug_id",
    "ligand_base",
    "display_name",
    "generic_name",
    "pubchem_name",
    "drugcentral_generic_name",
]


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _norm(value: Any) -> str:
    return _clean(value).lower()


def _row_first(row: pd.Series, fields: list[str]) -> tuple[str, str]:
    for field in fields:
        if field in row.index:
            value = _clean(row.get(field))
            if value:
                return value, field
    return "", ""


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_mapping_smiles(root: Path | None = None) -> dict[str, tuple[str, str]]:
    base = root or _repo_root()
    path = base / FDA_MAPPING_PATH
    if not path.exists():
        return {}
    wanted = {
        "path",
        "generic_name",
        "display_name",
        "pubchem_name",
        "drugcentral_generic_name",
        "smiles_neutral",
        "smiles",
        "canonical_smiles",
    }
    try:
        mapping = pd.read_csv(path, low_memory=False, usecols=lambda col: col in wanted)
    except Exception:
        return {}
    lookup: dict[str, tuple[str, str]] = {}
    for row in mapping.to_dict("records"):
        smiles = _clean(row.get("smiles_neutral") or row.get("canonical_smiles") or row.get("smiles"))
        if not smiles:
            continue
        source = str(path)
        pdbqt_path = _clean(row.get("path"))
        if pdbqt_path:
            lookup[f"base:{Path(pdbqt_path).stem.lower()}"] = (smiles, source)
        for field in ("generic_name", "display_name", "pubchem_name", "drugcentral_generic_name"):
            key = _norm(row.get(field))
            if key:
                lookup[f"name:{key}"] = (smiles, source)
    return lookup


def _mapping_smiles(row: pd.Series, lookup: dict[str, tuple[str, str]]) -> tuple[str, str]:
    for field in NAME_FIELDS:
        if field not in row.index:
            continue
        value = _clean(row.get(field))
        if not value:
            continue
        key_type = "base" if field == "ligand_base" else "name"
        found = lookup.get(f"{key_type}:{value.lower()}")
        if found:
            return found
    return "", ""


def _descriptor_record(smiles: str) -> dict[str, float] | None:
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    except Exception:
        return None
    try:
        mol = Chem.MolFromSmiles(str(smiles).strip())
    except Exception:
        mol = None
    if mol is None:
        return None
    mol_wt = getattr(Descriptors, "MolWt")
    mol_logp = getattr(Crippen, "MolLogP")
    hbd = getattr(Lipinski, "NumHDonors")
    hba = getattr(Lipinski, "NumHAcceptors")
    rotatable = getattr(Lipinski, "NumRotatableBonds")
    aromatic = getattr(Lipinski, "NumAromaticRings")
    return {
        "rdkit_mol_wt": float(mol_wt(mol)),
        "rdkit_mol_logp": float(mol_logp(mol)),
        "rdkit_tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "rdkit_hbd": float(hbd(mol)),
        "rdkit_hba": float(hba(mol)),
        "rdkit_rotatable_bonds": float(rotatable(mol)),
        "rdkit_formal_charge": float(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
        "rdkit_aromatic_rings": float(aromatic(mol)),
        "rdkit_fraction_csp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "rdkit_qed": float(QED.qed(mol)),
    }


def add_ligand_physchem_descriptors(
    df: pd.DataFrame,
    *,
    repo_root: Path | None = None,
    overwrite: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add RDKit ligand descriptors from row SMILES or the FDA PDBQT mapping.

    Missing descriptors remain NaN. The function never treats an unmapped ligand as
    chemically negative evidence; it records a status/source instead.
    """
    out = df.copy()
    for col in DESCRIPTOR_COLUMNS:
        if col not in out.columns or overwrite:
            out[col] = pd.NA
    if DESCRIPTOR_SOURCE_COL not in out.columns or overwrite:
        out[DESCRIPTOR_SOURCE_COL] = ""
    if DESCRIPTOR_STATUS_COL not in out.columns or overwrite:
        out[DESCRIPTOR_STATUS_COL] = ""

    lookup = _load_mapping_smiles(repo_root)
    cache: dict[str, dict[str, float] | None] = {}
    filled = 0
    mapped = 0
    row_smiles = 0
    parse_failed = 0
    missing = 0
    skipped_existing = 0

    for idx, row in out.iterrows():
        if not overwrite and all(pd.notna(out.at[idx, col]) for col in DESCRIPTOR_COLUMNS):
            skipped_existing += 1
            continue
        smiles, source = _row_first(row, SMILES_FIELDS)
        if smiles:
            row_smiles += 1
            source = source or "row_smiles"
        else:
            smiles, source = _mapping_smiles(row, lookup)
            if smiles:
                mapped += 1
        if not smiles:
            out.at[idx, DESCRIPTOR_STATUS_COL] = "missing_smiles"
            missing += 1
            continue
        if smiles not in cache:
            cache[smiles] = _descriptor_record(smiles)
        desc = cache[smiles]
        if desc is None:
            out.at[idx, DESCRIPTOR_STATUS_COL] = "rdkit_parse_failed"
            out.at[idx, DESCRIPTOR_SOURCE_COL] = source
            parse_failed += 1
            continue
        for col, value in desc.items():
            out.at[idx, col] = value
        out.at[idx, DESCRIPTOR_SOURCE_COL] = source
        out.at[idx, DESCRIPTOR_STATUS_COL] = "ok"
        filled += 1

    summary = {
        "ligand_descriptor_rows": int(len(out)),
        "ligand_descriptor_filled_rows": int(filled),
        "ligand_descriptor_row_smiles_rows": int(row_smiles),
        "ligand_descriptor_mapping_rows": int(mapped),
        "ligand_descriptor_missing_smiles_rows": int(missing),
        "ligand_descriptor_parse_failed_rows": int(parse_failed),
        "ligand_descriptor_skipped_existing_rows": int(skipped_existing),
        "ligand_descriptor_columns": DESCRIPTOR_COLUMNS,
    }
    return out, summary
