"""Meeko-first ligand preparation helpers.

This module keeps the publication-facing ligand-prep path centered on SDF
chemistry and Meeko PDBQT writing. Extracted PDB ligands are treated as
coordinate carriers: when possible, bond orders are restored from CCD/template
chemistry before the same SDF->PDBQT backend is used.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
import os
from pathlib import Path
import urllib.error
import urllib.request

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

from prep_ligands.prep_ligands_common_validation import (
    QUARANTINE_DIRNAME,
    assert_no_helium_in_pdbqt,
    is_valid_ligand,
    quick_pdbqt_validate,
    validate_pdbqt_invariants,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class MeekoPrepResult:
    ok: bool
    status: str
    writer: str
    detail: str = ""


@dataclass(frozen=True)
class RestoredLigandSdf:
    sdf_path: Path
    template_source: str
    restored_bond_orders: bool
    detail: str = ""


def _load_meeko_api():
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except Exception as exc:  # pragma: no cover - environment failure path
        raise RuntimeError(
            "Meeko ligand prep is required for the default ligand pipeline. "
            "Install `meeko` in the docking environment."
        ) from exc
    return MoleculePreparation, PDBQTWriterLegacy


def _mol_has_3d(mol: Chem.Mol) -> bool:
    return mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()


def _prepare_mol_for_meeko(mol: Chem.Mol, ligand_name: str) -> Chem.Mol:
    prepared = Chem.Mol(mol)
    if prepared.GetNumAtoms() == 0:
        raise ValueError("empty_molecule")
    try:
        Chem.SanitizeMol(prepared)
    except Exception:
        prepared.UpdatePropertyCache(strict=False)

    # Add only missing implicit hydrogens even when explicit isotope H/D atoms
    # are already present. RDKit preserves those existing atoms, isotopes, and
    # conformer coordinates; skipping AddHs entirely would leave the remaining
    # heavy atoms under-hydrogenated (for example, deutetrabenazine).
    prepared = Chem.AddHs(prepared, addCoords=_mol_has_3d(prepared))

    if not _mol_has_3d(prepared):
        params = AllChem.ETKDGv3()  # type: ignore[attr-defined]
        params.randomSeed = 0xA7A52
        rc = AllChem.EmbedMolecule(prepared, params)  # type: ignore[attr-defined]
        if rc != 0:
            rc = AllChem.EmbedMolecule(  # type: ignore[attr-defined]
                prepared, randomSeed=0xA7A52, useRandomCoords=True
            )
        if rc != 0:
            raise ValueError("rdkit_embed_failed")
        try:
            AllChem.UFFOptimizeMolecule(prepared, maxIters=200)  # type: ignore[attr-defined]
        except Exception:
            logger.debug("[meeko] UFF optimize skipped ligand=%s", ligand_name)

    if not prepared.HasProp("_Name"):
        prepared.SetProp("_Name", ligand_name)
    return prepared


def prepare_mol_to_pdbqt_with_meeko(
    mol: Chem.Mol,
    pdbqt_path: Path,
    *,
    ligand_name: str,
    log_dir: Path,
) -> MeekoPrepResult:
    """Write one RDKit molecule to PDBQT using Meeko and validate the result."""
    try:
        MoleculePreparation, PDBQTWriterLegacy = _load_meeko_api()
        prepped = _prepare_mol_for_meeko(mol, ligand_name)
        last_result = MeekoPrepResult(False, "prepare_fail:not_attempted", "meeko")
        for rigid_macrocycles in (False, True):
            preparator = MoleculePreparation(rigid_macrocycles=rigid_macrocycles)
            setups = (
                preparator.prepare(prepped)
                if hasattr(preparator, "prepare")
                else preparator(prepped)
            )
            if not setups:
                last_result = MeekoPrepResult(
                    False, "prepare_fail:no_setups", "meeko"
                )
                continue

            written = PDBQTWriterLegacy.write_string(setups[0])
            if isinstance(written, tuple):
                pdbqt_text, ok, err = written
                if not ok:
                    last_result = MeekoPrepResult(
                        False, "prepare_fail:writer", "meeko", str(err)
                    )
                    continue
            else:
                pdbqt_text = str(written)

            pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
            pdbqt_path.write_text(pdbqt_text, encoding="utf-8")
            lines = pdbqt_text.splitlines(keepends=True)
            new_lines, fixes, q_reason = assert_no_helium_in_pdbqt(lines, ligand_name)
            if q_reason:
                last_result = MeekoPrepResult(
                    False, "postcheck_fail:helium", "meeko", q_reason
                )
                if not rigid_macrocycles:
                    continue
                _quarantine(pdbqt_path)
                return last_result
            if fixes:
                pdbqt_path.write_text("".join(new_lines), encoding="utf-8")

            ok_quick, _n_atoms, reason = quick_pdbqt_validate(pdbqt_path)
            if not ok_quick:
                last_result = MeekoPrepResult(
                    False, f"postcheck_fail:{reason}", "meeko"
                )
                if not rigid_macrocycles:
                    continue
                _quarantine(pdbqt_path)
                return last_result

            inv_ok, _metrics, inv_fail = validate_pdbqt_invariants(pdbqt_path)
            if not inv_ok:
                last_result = MeekoPrepResult(
                    False,
                    "postcheck_fail:" + ",".join(inv_fail),
                    "meeko",
                )
                if not rigid_macrocycles:
                    continue
                _quarantine(pdbqt_path)
                return last_result

            if not is_valid_ligand(pdbqt_path, log_dir=log_dir):
                last_result = MeekoPrepResult(
                    False, "postcheck_fail:is_valid_ligand", "meeko"
                )
                if not rigid_macrocycles:
                    continue
                _quarantine(pdbqt_path)
                return last_result
            status = "ok_rigid_macrocycle" if rigid_macrocycles else "ok"
            return MeekoPrepResult(True, status, "meeko")
        if pdbqt_path.exists():
            _quarantine(pdbqt_path)
        return last_result
    except Exception as exc:
        logger.warning("[meeko] ligand=%s failed: %s", ligand_name, exc)
        return MeekoPrepResult(False, "prepare_fail:meeko_exception", "meeko", str(exc))


def prepare_sdf_to_pdbqt_with_meeko(
    sdf_path: Path,
    pdbqt_path: Path,
    *,
    ligand_name: str | None = None,
    log_dir: Path,
) -> MeekoPrepResult:
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = next((candidate for candidate in supplier if candidate is not None), None)
    if mol is None:
        return MeekoPrepResult(False, "prepare_fail:sdf_read", "meeko")
    name = ligand_name or (mol.GetProp("_Name").strip() if mol.HasProp("_Name") else "")
    return prepare_mol_to_pdbqt_with_meeko(
        mol,
        pdbqt_path,
        ligand_name=name or sdf_path.stem,
        log_dir=log_dir,
    )


def restore_extracted_pdb_to_sdf(
    pdb_path: Path,
    sdf_path: Path,
    *,
    residue_name: str,
    template_cache_dir: Path,
) -> RestoredLigandSdf:
    """Restore extracted PDB ligand chemistry and write coordinate-preserving SDF."""
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    template = _load_template_mol(residue_name, template_cache_dir)
    pdb_mol = Chem.MolFromPDBFile(str(pdb_path), sanitize=False, removeHs=True)
    if pdb_mol is None or pdb_mol.GetNumAtoms() == 0:
        raise ValueError(f"could_not_read_extracted_pdb:{pdb_path}")

    restored = None
    template_source = "none"
    detail = ""
    if template is not None:
        try:
            template_source = template.GetProp("_Name") if template.HasProp("_Name") else residue_name
            template_no_h = Chem.RemoveHs(template, sanitize=False)
            restored = AllChem.AssignBondOrdersFromTemplate(template_no_h, pdb_mol)
            Chem.SanitizeMol(restored)
        except Exception as exc:
            detail = f"template_assign_failed:{exc}"
            logger.warning(
                "[ligprep-extracted] CCD/template bond-order restore failed residue=%s pdb=%s err=%s",
                residue_name,
                pdb_path.name,
                exc,
            )
            restored = _restore_template_coordinates_by_mcs(template, pdb_mol, pdb_path)
            if restored is not None:
                detail += ";template_mcs_coordinate_transfer"

    if restored is None:
        restored = Chem.Mol(pdb_mol)
        try:
            Chem.SanitizeMol(restored)
        except Exception:
            restored.UpdatePropertyCache(strict=False)
        template_source = "pdb_proximity_fallback"

    restored.SetProp("_Name", pdb_path.stem)
    restored_h = Chem.AddHs(restored, addCoords=True)
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(restored_h)
    writer.close()
    return RestoredLigandSdf(
        sdf_path=sdf_path,
        template_source=template_source,
        restored_bond_orders=template is not None and restored is not None and "pdb_proximity_fallback" not in template_source,
        detail=detail,
    )


def _remove_all_hydrogens(mol: Chem.Mol) -> Chem.Mol:
    editable = Chem.RWMol(mol)
    for atom_idx in sorted(
        [atom.GetIdx() for atom in editable.GetAtoms() if atom.GetAtomicNum() == 1],
        reverse=True,
    ):
        editable.RemoveAtom(int(atom_idx))
    return editable.GetMol()


def _restore_template_coordinates_by_mcs(
    template: Chem.Mol,
    pdb_mol: Chem.Mol,
    pdb_path: Path,
) -> Chem.Mol | None:
    """Use CCD/template chemistry with extracted PDB coordinates when MCS is complete."""
    try:
        template_heavy = _remove_all_hydrogens(template)
        pdb_heavy = _remove_all_hydrogens(pdb_mol)
        if (
            template_heavy.GetNumAtoms() == 0
            or pdb_heavy.GetNumAtoms() == 0
            or pdb_heavy.GetNumConformers() == 0
        ):
            return None
        mcs = rdFMCS.FindMCS(
            [template_heavy, pdb_heavy],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
            timeout=10,
        )
        if mcs.canceled or mcs.numAtoms != pdb_heavy.GetNumAtoms():
            return None
        pattern = Chem.MolFromSmarts(mcs.smartsString)
        if pattern is None:
            return None
        template_match = template_heavy.GetSubstructMatch(pattern)
        pdb_match = pdb_heavy.GetSubstructMatch(pattern)
        if len(template_match) != pdb_heavy.GetNumAtoms() or len(pdb_match) != pdb_heavy.GetNumAtoms():
            return None

        source_conf = pdb_heavy.GetConformer()
        conf = Chem.Conformer(template_heavy.GetNumAtoms())
        mapped_template = set()
        for template_idx, pdb_idx in zip(template_match, pdb_match):
            pos = source_conf.GetAtomPosition(int(pdb_idx))
            conf.SetAtomPosition(int(template_idx), pos)
            mapped_template.add(int(template_idx))
        if len(mapped_template) != template_heavy.GetNumAtoms():
            return None

        restored = Chem.Mol(template_heavy)
        restored.RemoveAllConformers()
        conf.Set3D(True)
        restored.AddConformer(conf)
        Chem.SanitizeMol(restored)
        logger.info(
            "[ligprep-extracted] CCD/template MCS coordinate transfer succeeded pdb=%s atoms=%d",
            pdb_path.name,
            restored.GetNumAtoms(),
        )
        return restored
    except Exception as exc:
        logger.warning(
            "[ligprep-extracted] CCD/template MCS coordinate transfer failed pdb=%s err=%s",
            pdb_path.name,
            exc,
        )
        return None


def _load_template_mol(residue_name: str, template_cache_dir: Path) -> Chem.Mol | None:
    ccd_id = residue_name.strip().upper()
    if not ccd_id or ccd_id in {"UNL", "LIG"}:
        return None

    for candidate in _template_candidates(ccd_id, template_cache_dir):
        mol = _read_template(candidate)
        if mol is not None:
            mol.SetProp("_Name", str(candidate))
            return mol

    downloaded = _download_ccd_ideal_sdf(ccd_id, template_cache_dir)
    if downloaded is not None:
        mol = _read_template(downloaded)
        if mol is not None:
            mol.SetProp("_Name", str(downloaded))
            return mol
    return None


def _template_candidates(ccd_id: str, template_cache_dir: Path) -> list[Path]:
    roots = [template_cache_dir]
    extra = os.environ.get("LIGPREP_TEMPLATE_SDF_DIR") or os.environ.get("CCD_TEMPLATE_DIR")
    if extra:
        roots.insert(0, Path(extra).expanduser())
    names = [
        f"{ccd_id}_ideal.sdf",
        f"{ccd_id}.sdf",
        f"{ccd_id}_model.sdf",
        f"{ccd_id}.mol",
    ]
    return [root / name for root in roots for name in names]


def _read_template(path: Path) -> Chem.Mol | None:
    if not path.is_file():
        return None
    if path.suffix.lower() == ".sdf":
        supplier = Chem.SDMolSupplier(str(path), removeHs=True)
        return next((mol for mol in supplier if mol is not None), None)
    if path.suffix.lower() == ".mol":
        return Chem.MolFromMolFile(str(path), removeHs=True)
    return None


def _download_ccd_ideal_sdf(ccd_id: str, template_cache_dir: Path) -> Path | None:
    if str(os.environ.get("LIGPREP_ALLOW_CCD_DOWNLOAD", "1")).lower() in {
        "0",
        "false",
        "no",
    }:
        return None
    template_cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = template_cache_dir / f"{ccd_id}_ideal.sdf"
    url = f"https://files.rcsb.org/ligands/download/{ccd_id}_ideal.sdf"
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            payload = response.read()
        if payload:
            out_path.write_bytes(payload)
            return out_path
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        logger.info("[ligprep-extracted] CCD download unavailable ccd=%s err=%s", ccd_id, exc)
    return None


def _quarantine(pdbqt_path: Path) -> None:
    try:
        if not pdbqt_path.exists():
            return
        quarantine = pdbqt_path.parent / QUARANTINE_DIRNAME
        quarantine.mkdir(exist_ok=True)
        pdbqt_path.replace(quarantine / pdbqt_path.name)
    except Exception:
        pass


__all__ = [
    "MeekoPrepResult",
    "RestoredLigandSdf",
    "prepare_mol_to_pdbqt_with_meeko",
    "prepare_sdf_to_pdbqt_with_meeko",
    "restore_extracted_pdb_to_sdf",
]
