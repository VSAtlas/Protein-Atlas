"""Extracted-ligand MOL2 candidate generation for crystal prep."""

from dataclasses import dataclass
import logging
import os
from pathlib import Path
from typing import List, Optional

from rdkit import Chem

from prep_ligands.prep_ligands_common_conversion import (
    _log_std_diff,
    _pick_best_existing_mol2,
    _re_aromatize_mol2_in_place,
    _run_obabel,
    standardize_mol_with_activesite,
)
from prep_ligands.prep_ligands_reporting import _audit_protonation_metrics


@dataclass
class ExtractedMol2Selection:
    selected_mol2: Path
    strategy_labels: List[str]
    fallbacks_used: List[str]


def prepare_extracted_mol2_selection(
    *,
    tmp_pdb: Path,
    sanitized: Path,
    tmp_mol2: Path,
    lig_id: str,
    prepped_ligands_dir: Path,
    obabel_exe_short: Optional[str],
) -> ExtractedMol2Selection:
    _audit_protonation_metrics(tag=lig_id, mol_or_path=tmp_pdb, context="pre")
    strategy_labels: List[str] = []
    fallbacks_used: List[str] = []
    proto_a_path = sanitized.with_suffix(".protoA.mol2")
    proto_b_path = sanitized.with_suffix(".protoB.mol2")

    try:
        mol_h = Chem.MolFromPDBFile(str(tmp_pdb), sanitize=False, removeHs=False)
        if mol_h is None:
            raise ValueError("rdkit_from_pdb_failed")
        mol_h = Chem.AddHs(mol_h, addCoords=True)
        try:
            Chem.Kekulize(mol_h, clearAromaticFlags=True)
        except Exception:
            pass
        writer = getattr(Chem, "MolToMol2File", None)
        if not callable(writer):
            raise AttributeError("RDKit MolToMol2File is unavailable in this build")
        writer(mol_h, str(proto_a_path))
        if proto_a_path.exists() and proto_a_path.stat().st_size > 100:
            strategy_labels.append("RDKit-AddHs")
        else:
            fallbacks_used.append("RDKit-AddHs-fail")
    except Exception:
        fallbacks_used.append("RDKit-AddHs-fail")
    if obabel_exe_short:
        _run_obabel(
            [
                obabel_exe_short,
                "-ipdb",
                str(sanitized),
                "-omol2",
                "-O",
                str(proto_b_path),
                "-p",
                os.environ.get("LIGPREP_PH", "7.4"),
                "--partialcharge",
                "gasteiger",
            ],
            timeout_sec=600,
        )
    selected_mol2 = _pick_best_existing_mol2(
        (proto_a_path, proto_b_path),
        tmp_mol2,
    )
    if selected_mol2 == tmp_mol2:
        strategy_labels.append("no-extra-H (ADT adds)")
        fallbacks_used.append("ADT/Meeko-rescue-possible")

    if obabel_exe_short:
        try:
            _re_aromatize_mol2_in_place(selected_mol2, obabel_exe_short)
        except Exception as exc:
            logging.error("[ligprep] re_arom crash for %s: %s", selected_mol2.name, exc)

    try:
        mol_raw = Chem.MolFromMol2File(str(selected_mol2), sanitize=False, removeHs=False)
        old_smiles = ""
        try:
            mol_old = Chem.MolFromMol2File(str(selected_mol2), sanitize=True, removeHs=False)
            old_smiles = (
                Chem.MolToSmiles(mol_old, isomericSmiles=True) if mol_old is not None else ""
            )
        except Exception:
            pass
        mol_std = standardize_mol_with_activesite(mol_raw)
        new_smiles = Chem.MolToSmiles(mol_std, isomericSmiles=True) if mol_std is not None else ""
        if old_smiles and new_smiles and old_smiles != new_smiles:
            logging.warning(
                "[std:audit][extracted] %s: SMILES changed %s -> %s",
                selected_mol2.name,
                old_smiles,
                new_smiles,
            )
            _log_std_diff(prepped_ligands_dir, lig_id, "extracted", old_smiles, new_smiles)
    except Exception as exc:
        logging.warning("[std:extracted] standardization skipped due to error: %s", exc)
    return ExtractedMol2Selection(
        selected_mol2=selected_mol2,
        strategy_labels=strategy_labels,
        fallbacks_used=fallbacks_used,
    )


__all__ = ["ExtractedMol2Selection", "prepare_extracted_mol2_selection"]
