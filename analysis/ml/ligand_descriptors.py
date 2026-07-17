from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from analysis.ml.ligand_functional_topology import (
    LIGAND_STRUCTURE_DESCRIPTOR_GROUPS,
)


PHYSICHEM_DESCRIPTOR_COLUMNS = [
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
SMILES_FIELDS = [
    "smiles",
    "canonical_smiles",
    "ligand_smiles",
    "smiles_neutral",
    "isomeric_smiles",
]
FDA_MAPPING_PATH = Path("chemdb/data/fda_mapping_from_pdbqt.csv")
IDENTITY_FIELDS = ["rdk_id", "_join_rdk", "ligand_rdk_id", "rdk", "ligand_base"]
NAME_FIELDS = [
    "drug_id",
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
        "rdk_id",
        "generic_name",
        "display_name",
        "pubchem_name",
        "drugcentral_generic_name",
        "smiles_neutral",
        "smiles",
        "canonical_smiles",
    }
    try:
        mapping = pd.read_csv(
            path,
            low_memory=False,
            usecols=lambda col: col in wanted,
        )
    except Exception:
        return {}
    lookup: dict[str, tuple[str, str]] = {}
    for row in mapping.to_dict("records"):
        smiles = _clean(
            row.get("smiles_neutral")
            or row.get("canonical_smiles")
            or row.get("smiles")
        )
        if not smiles:
            continue
        source = str(path)
        rdk_id = _norm(row.get("rdk_id"))
        if rdk_id:
            lookup[f"base:{rdk_id}"] = (smiles, source)
        pdbqt_path = _clean(row.get("path"))
        if pdbqt_path:
            lookup[f"base:{Path(pdbqt_path).stem.lower()}"] = (smiles, source)
        for field in (
            "generic_name",
            "display_name",
            "pubchem_name",
            "drugcentral_generic_name",
        ):
            key = _norm(row.get(field))
            if key:
                lookup[f"name:{key}"] = (smiles, source)
    return lookup


def _mapping_smiles(
    row: pd.Series,
    lookup: dict[str, tuple[str, str]],
) -> tuple[str, str]:
    for field in IDENTITY_FIELDS:
        if field not in row.index:
            continue
        value = _norm(row.get(field))
        if value and (found := lookup.get(f"base:{value}")):
            return found
    for field in NAME_FIELDS:
        if field not in row.index:
            continue
        value = _norm(row.get(field))
        if not value:
            continue
        found = lookup.get(f"name:{value}")
        if found:
            return found
    return "", ""


DescriptorRecord = dict[str, float | int]
DescriptorCalculator = Callable[[Any], DescriptorRecord]


@dataclass(frozen=True)
class DescriptorBlock:
    name: str
    columns: tuple[str, ...]
    calculate: DescriptorCalculator


def _physchem_descriptor_record(mol: Any) -> DescriptorRecord:
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

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
        "rdkit_formal_charge": float(
            sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        ),
        "rdkit_aromatic_rings": float(aromatic(mol)),
        "rdkit_fraction_csp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "rdkit_qed": float(QED.qed(mol)),
    }


DESCRIPTOR_BLOCKS = (
    DescriptorBlock(
        name="physchem",
        columns=tuple(PHYSICHEM_DESCRIPTOR_COLUMNS),
        calculate=_physchem_descriptor_record,
    ),
) + tuple(
    DescriptorBlock(
        name=group.name,
        columns=group.columns,
        calculate=group.calculate,
    )
    for group in LIGAND_STRUCTURE_DESCRIPTOR_GROUPS
)
DESCRIPTOR_GROUP_COLUMNS = {block.name: block.columns for block in DESCRIPTOR_BLOCKS}
DESCRIPTOR_COLUMNS = [column for block in DESCRIPTOR_BLOCKS for column in block.columns]
if len(DESCRIPTOR_COLUMNS) != len(set(DESCRIPTOR_COLUMNS)):
    raise RuntimeError("registered ligand descriptor columns must be unique")


def _calculate_registered_descriptors(mol: Any) -> DescriptorRecord:
    record: DescriptorRecord = {}
    for block in DESCRIPTOR_BLOCKS:
        values = block.calculate(mol)
        if set(values) != set(block.columns):
            raise RuntimeError(
                f"descriptor block {block.name!r} violated its column contract"
            )
        record.update(values)
    return record


def _descriptor_record(smiles: str) -> DescriptorRecord | None:
    try:
        from rdkit import Chem
    except Exception:
        return None
    try:
        mol = Chem.MolFromSmiles(str(smiles).strip())
    except Exception:
        mol = None
    if mol is None:
        return None
    return _calculate_registered_descriptors(mol)


def add_ligand_physchem_descriptors(
    df: pd.DataFrame,
    *,
    repo_root: Path | None = None,
    overwrite: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add registered RDKit descriptors from row SMILES or the FDA mapping.

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
    cache: dict[str, DescriptorRecord | None] = {}
    filled = 0
    mapped = 0
    row_smiles = 0
    parse_failed = 0
    missing = 0
    skipped_existing = 0

    for idx, row in out.iterrows():
        if not overwrite and all(
            pd.notna(out.at[idx, col]) for col in DESCRIPTOR_COLUMNS
        ):
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
            if overwrite or pd.isna(out.at[idx, col]):
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
        "ligand_descriptor_groups": {
            name: list(columns) for name, columns in DESCRIPTOR_GROUP_COLUMNS.items()
        },
    }
    return out, summary
