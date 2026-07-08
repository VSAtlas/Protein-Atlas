from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.ligand_descriptors import SMILES_FIELDS, _load_mapping_smiles, _mapping_smiles, _row_first

GENERIC_CHEMOTYPES = {
    "",
    "nan",
    "none",
    "null",
    "other",
    "unknown",
    "unassigned",
    "other / unassigned",
    "other / structural scaffold assigned",
    "structural scaffold assigned",
}

CHEMOTYPE_PATTERNS: tuple[tuple[str, str], ...] = (
    ("beta_lactam", "C1C(=O)N1"),
    ("sulfonamide", "S(=O)(=O)N"),
    ("carboxylic_acid", "C(=O)[O;H1,-1]"),
    ("phosphate_or_phosphonate", "P(=O)(O)O"),
    ("quinolone_like", "c1ccc2c(c1)ncc(=O)o2"),
    ("benzodiazepine_like", "c1ccc2c(c1)NCN=C2"),
    ("phenothiazine_like", "S1c2ccccc2Nc2ccccc12"),
    ("xanthine_like", "O=c1[nH]cnc2[nH]cnc12"),
    ("azole", "[nH0,nH1]1cccc1"),
    ("nucleoside_like", "n1cnc2c(N)ncnc12"),
    ("lactone", "C(=O)O[C;R]"),
    ("urea_or_carbamate", "N-C(=O)-N"),
)


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _norm(value: Any) -> str:
    return _clean(value).lower()


def _is_generic(value: Any) -> bool:
    return _norm(value) in GENERIC_CHEMOTYPES


def _safe_token(value: str) -> str:
    cleaned = _norm(value).replace("-", "_").replace("/", "_")
    cleaned = re.sub(r"[^a-z0-9_]+", "_", cleaned).strip("_")
    return cleaned or "structural_other"


def _safe_chemotype_token(value: str) -> str:
    if value.startswith("murcko:"):
        digest = hashlib.sha1(value.encode("utf-8")).hexdigest()[:12]
        return f"murcko_{digest}"
    return _safe_token(value)


def _smiles_for_row(row: pd.Series, lookup: dict[str, tuple[str, str]]) -> tuple[str, str]:
    smiles, source = _row_first(row, SMILES_FIELDS)
    if smiles:
        return smiles, source or "row_smiles"
    return _mapping_smiles(row, lookup)


def _pattern_cache() -> list[tuple[str, Any]]:
    from rdkit import Chem

    cache: list[tuple[str, Any]] = []
    for name, smarts in CHEMOTYPE_PATTERNS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            cache.append((name, patt))
    return cache


def _chemotype_from_smiles(smiles: str, patterns: list[tuple[str, Any]]) -> tuple[str, str]:
    from rdkit import Chem
    from rdkit.Chem import Lipinski, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return "", "rdkit_parse_failed"
    for name, patt in patterns:
        if mol.HasSubstructMatch(patt):
            return name, "rdkit_smarts"

    atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    heavy = max(1, mol.GetNumHeavyAtoms())
    oxygen_fraction = atom_symbols.count("O") / heavy
    nitrogen_fraction = atom_symbols.count("N") / heavy
    aromatic_rings = Lipinski.NumAromaticRings(mol)
    aliphatic_rings = Lipinski.NumAliphaticRings(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    fraction_csp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot = Lipinski.NumRotatableBonds(mol)

    if ring_count >= 4 and aliphatic_rings >= 3 and fraction_csp3 >= 0.45:
        return "steroid_like", "rdkit_shape_rule"
    if oxygen_fraction >= 0.30 and hbd >= 3 and hba >= 4:
        return "glycoside_or_polyhydroxy", "rdkit_composition_rule"
    if aromatic_rings >= 2 and nitrogen_fraction >= 0.08:
        return "aromatic_n_heterocycle", "rdkit_composition_rule"
    if aromatic_rings >= 2:
        return "polyaromatic", "rdkit_composition_rule"
    if aromatic_rings == 1:
        return "monoaromatic", "rdkit_composition_rule"
    if rot >= 8 and heavy >= 30:
        return "flexible_large_lipophile", "rdkit_composition_rule"
    if fraction_csp3 >= 0.65:
        return "saturated_aliphatic", "rdkit_composition_rule"

    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    if scaffold:
        return f"murcko:{scaffold}", "rdkit_murcko_fallback"
    return "structural_other", "rdkit_fallback"


def add_structural_ligand_chemotypes(
    df: pd.DataFrame,
    *,
    repo_root: Path | None = None,
    overwrite_generic_only: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Backfill ligand chemotype from SMILES/RDKit without dropping original columns."""

    out = df.copy()
    original = out["ligand_chemotype"] if "ligand_chemotype" in out.columns else pd.Series("", index=out.index)
    if "ligand_chemotype_original" not in out.columns:
        out["ligand_chemotype_original"] = original
    for col in [
        "ligand_chemotype_structural",
        "ligand_chemotype_structural_raw",
        "ligand_chemotype_broad",
        "ligand_chemotype_source",
        "ligand_chemotype_status",
    ]:
        if col not in out.columns:
            out[col] = ""

    lookup = _load_mapping_smiles(repo_root)
    patterns: list[tuple[str, Any]] = []
    try:
        patterns = _pattern_cache()
    except Exception:
        pass

    cache: dict[str, tuple[str, str]] = {}
    filled = 0
    preserved_existing = 0
    missing_smiles = 0
    parse_failed = 0
    generic_before = int(original.map(_is_generic).sum())

    for idx, row in out.iterrows():
        current = original.loc[idx] if idx in original.index else ""
        should_fill = not overwrite_generic_only or _is_generic(current)
        if not should_fill:
            preserved_existing += 1
            safe_current = _safe_chemotype_token(_clean(current))
            out.at[idx, "ligand_chemotype_structural_raw"] = _clean(current)
            out.at[idx, "ligand_chemotype_structural"] = safe_current
            out.at[idx, "ligand_chemotype_broad"] = _safe_token(_clean(current))
            out.at[idx, "ligand_chemotype_source"] = "existing_ligand_chemotype"
            out.at[idx, "ligand_chemotype_status"] = "preserved_existing"
            continue
        smiles, smiles_source = _smiles_for_row(row, lookup)
        if not smiles:
            missing_smiles += 1
            out.at[idx, "ligand_chemotype_status"] = "missing_smiles"
            continue
        if smiles not in cache:
            try:
                cache[smiles] = _chemotype_from_smiles(smiles, patterns)
            except Exception:
                cache[smiles] = ("", "rdkit_error")
        chemotype, source = cache[smiles]
        if not chemotype:
            parse_failed += 1
            out.at[idx, "ligand_chemotype_status"] = source or "rdkit_parse_failed"
            continue
        safe_chemotype = _safe_chemotype_token(chemotype)
        broad = _safe_token(chemotype.split(":", 1)[0])
        out.at[idx, "ligand_chemotype_structural_raw"] = chemotype
        out.at[idx, "ligand_chemotype_structural"] = safe_chemotype
        out.at[idx, "ligand_chemotype_broad"] = broad
        out.at[idx, "ligand_chemotype_source"] = f"{source}:{smiles_source}"
        out.at[idx, "ligand_chemotype_status"] = "filled_structural"
        out.at[idx, "ligand_chemotype"] = safe_chemotype
        filled += 1

    summary = {
        "ligand_chemotype_rows": int(len(out)),
        "ligand_chemotype_generic_before": int(generic_before),
        "ligand_chemotype_filled_structural_rows": int(filled),
        "ligand_chemotype_preserved_existing_rows": int(preserved_existing),
        "ligand_chemotype_missing_smiles_rows": int(missing_smiles),
        "ligand_chemotype_parse_failed_rows": int(parse_failed),
        "ligand_chemotype_policy": "Existing non-generic chemotypes are preserved; generic/missing values are backfilled from SMILES/RDKit with provenance.",
    }
    return out, summary
