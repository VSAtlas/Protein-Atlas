from __future__ import annotations

from dataclasses import dataclass
import importlib
import json
import logging
import os
from pathlib import Path
from typing import Any, Mapping

from protein_prep.pdb_records import pdbqt_line_element, pdbqt_line_xyz, record_type
from post_docking.mmgbsa.rdkit_sdf_io import write_rdkit_mol_to_sdf
from post_docking.mmgbsa.operation_exceptions import (
    INT_COERCE_ERRORS,
    MMGBSA_POSE_EXPORT_ERRORS,
    OPTIONAL_IMPORT_ERRORS,
    PATH_RESOLVE_ERRORS,
    RDKIT_PROPERTY_ERRORS,
    TEXT_FILE_READ_ERRORS,
)


@dataclass(frozen=True)
class PoseSdfResult:
    ok: bool
    status: str
    sdf_path: Path | None = None
    metadata_path: Path | None = None


def _source_sidecar_path(ligand_pdbqt: Path) -> Path:
    return ligand_pdbqt.with_suffix(".ligprep_source.json")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        text = path.read_text(encoding="utf-8")
    except TEXT_FILE_READ_ERRORS:
        return {}
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        return {}
    return payload if isinstance(payload, dict) else {}


def _ligand_base(path: Path) -> str:
    stem = path.stem
    for suffix in (".sanitized", "_ctrl_redock", "_ctrl_redock_vina"):
        stem = stem.replace(suffix, "")
    return stem


def _fallback_crystal_source_meta(ligand_path: Path) -> tuple[dict[str, Any], Path | None]:
    """Find extracted crystal-control SDF provenance for legacy PDBQTs without sidecars."""
    try:
        parts = ligand_path.resolve(strict=False).parts
    except PATH_RESOLVE_ERRORS:
        parts = ligand_path.parts
    try:
        idx = parts.index("prepped_ligands")
    except ValueError:
        return {}, None
    if len(parts) <= idx + 1:
        return {}, None

    pdb_id = parts[idx + 1]
    base = _ligand_base(ligand_path)
    repo_root = Path(*parts[:idx]) if idx > 0 else Path.cwd()
    run_id = (os.environ.get("ATLAS_RUN_ID") or "").strip()
    if not run_id:
        return {}, None
    ref_dir = repo_root / "outputs" / "processed_pdbs" / run_id / pdb_id / "reference"
    candidates = [
        ref_dir / f"{base}.restored.sdf",
        ref_dir / f"{base}.sanitized.sdf",
        ref_dir / f"{base}.sdf",
    ]
    for source_sdf in candidates:
        if not source_sdf.is_file() or source_sdf.stat().st_size <= 0:
            continue
        return (
            {
                "schema_version": 1,
                "writer": "atlas_mmgbsa_pose_artifacts",
                "source_format": "sdf",
                "source_sdf": str(source_sdf),
                "source_record_index": 1,
                "source_record_name": base,
                "source_kind": "legacy_extracted_reference_fallback",
                "chemistry_authoritative": False,
                "notes": "Recovered source SDF provenance from extracted crystal-control reference artifacts.",
            },
            None,
        )
    return {}, None


_pdbqt_atom_element = pdbqt_line_element


def _first_model_pdbqt_coords(pdbqt_path: Path) -> list[tuple[str, float, float, float]]:
    coords: list[tuple[str, float, float, float]] = []
    in_model = False
    saw_model = False
    with pdbqt_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            rec = record_type(line)
            if rec == "MODEL":
                if saw_model:
                    break
                saw_model = True
                in_model = True
                continue
            if rec == "ENDMDL" and in_model:
                break
            if saw_model and not in_model:
                continue
            if rec not in {"ATOM", "HETATM"}:
                continue
            xyz = pdbqt_line_xyz(line)
            if xyz is None:
                continue
            x, y, z = xyz
            coords.append((_pdbqt_atom_element(line), x, y, z))
    return coords


def _source_mol_from_sdf(source_sdf: Path, record_index: int):
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(source_sdf), removeHs=False)
    wanted = max(1, int(record_index))
    for idx, mol in enumerate(supplier, start=1):
        if idx == wanted:
            return mol
    return None


def _infer_hydrogen_position(mol, atom_idx: int, conf, new_positions: dict[int, tuple[float, float, float]]):
    atom = mol.GetAtomWithIdx(atom_idx)
    for neighbor in atom.GetNeighbors():
        nbr_idx = int(neighbor.GetIdx())
        if nbr_idx not in new_positions:
            continue
        base = new_positions[nbr_idx]
        if conf is None:
            return base
        try:
            old_h = conf.GetAtomPosition(atom_idx)
            old_n = conf.GetAtomPosition(nbr_idx)
            return (
                base[0] + float(old_h.x - old_n.x),
                base[1] + float(old_h.y - old_n.y),
                base[2] + float(old_h.z - old_n.z),
            )
        except RDKIT_PROPERTY_ERRORS:
            return base
    return (0.0, 0.0, 0.0)


def _pose_mol_with_docked_coords(source_mol, pdbqt_coords: list[tuple[str, float, float, float]]):
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    mol = Chem.Mol(source_mol)
    if mol.GetNumConformers() == 0:
        mol = Chem.AddHs(mol, addCoords=False)
    old_conf = mol.GetConformer() if mol.GetNumConformers() else None

    all_xyz = [(x, y, z) for _elem, x, y, z in pdbqt_coords]
    heavy_xyz = [(x, y, z) for elem, x, y, z in pdbqt_coords if elem.upper() != "H"]
    heavy_atom_indexes = [
        int(atom.GetIdx()) for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1
    ]

    inferred_h = False
    new_positions: dict[int, tuple[float, float, float]] = {}
    if len(all_xyz) == mol.GetNumAtoms():
        new_positions = {idx: xyz for idx, xyz in enumerate(all_xyz)}
    elif len(heavy_xyz) == len(heavy_atom_indexes):
        new_positions = {
            atom_idx: xyz for atom_idx, xyz in zip(heavy_atom_indexes, heavy_xyz)
        }
        inferred_h = True
        for atom in mol.GetAtoms():
            atom_idx = int(atom.GetIdx())
            if atom_idx in new_positions:
                continue
            new_positions[atom_idx] = _infer_hydrogen_position(
                mol, atom_idx, old_conf, new_positions
            )
    else:
        raise ValueError(
            f"atom_count_mismatch:pdbqt_atoms={len(all_xyz)} "
            f"pdbqt_heavy={len(heavy_xyz)} source_atoms={mol.GetNumAtoms()} "
            f"source_heavy={len(heavy_atom_indexes)}"
        )

    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom_idx in range(mol.GetNumAtoms()):
        x, y, z = new_positions[atom_idx]
        conf.SetAtomPosition(atom_idx, Point3D(float(x), float(y), float(z)))
    conf.Set3D(True)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol, inferred_h


def preserve_authoritative_pose_sdf(
    *,
    ligand_pdbqt: str | Path,
    docked_pdbqt: str | Path,
    cfg: Mapping[str, Any],
    logger: logging.Logger | None = None,
) -> PoseSdfResult:
    log = logger or logging.getLogger("mmgbsa.pose_artifacts")
    ligand_path = Path(ligand_pdbqt)
    docked_path = Path(docked_pdbqt)
    sidecar = _source_sidecar_path(ligand_path)
    if not sidecar.is_file():
        source_meta, sidecar_for_meta = _fallback_crystal_source_meta(ligand_path)
        if not source_meta:
            return PoseSdfResult(False, "missing_ligprep_source_sidecar")
    else:
        source_meta = _read_json(sidecar)
        sidecar_for_meta = sidecar
    source_sdf = Path(str(source_meta.get("source_sdf", "") or "")).expanduser()
    if not source_sdf.is_file():
        return PoseSdfResult(False, "missing_source_sdf")
    if not docked_path.is_file():
        return PoseSdfResult(False, "missing_docked_pdbqt")

    out_sdf = docked_path.with_suffix(".sdf")
    out_meta = docked_path.with_suffix(".mmgbsa_pose.json")
    if (
        out_sdf.is_file()
        and out_sdf.stat().st_size > 0
        and not bool(cfg.get("_MMGBSA_POSE_SDF_FORCE", False))
    ):
        return PoseSdfResult(True, "exists", out_sdf, out_meta if out_meta.exists() else None)

    try:
        record_index = int(source_meta.get("source_record_index", 1) or 1)
    except INT_COERCE_ERRORS:
        record_index = 1
    try:
        importlib.import_module("rdkit")
    except OPTIONAL_IMPORT_ERRORS as exc:
        log.warning(
            "[mmgbsa.pose_sdf] action=rdkit_import_failed err=%s",
            exc,
            exc_info=True,
        )
        return PoseSdfResult(False, "rdkit_unavailable")
    try:
        source_mol = _source_mol_from_sdf(source_sdf, record_index)
        if source_mol is None:
            return PoseSdfResult(False, "source_sdf_record_unreadable")
        coords = _first_model_pdbqt_coords(docked_path)
        if not coords:
            return PoseSdfResult(False, "docked_pdbqt_has_no_coordinates")
        pose_mol, inferred_h = _pose_mol_with_docked_coords(source_mol, coords)
        authoritative = source_meta.get("chemistry_authoritative") is not False
        pose_mol.SetProp("_Name", docked_path.stem)
        pose_mol.SetProp("MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE", str(authoritative).lower())
        pose_mol.SetProp(
            "MMGBSA_INPUT_CHEMISTRY_SOURCE",
            str(sidecar_for_meta) if sidecar_for_meta else str(source_sdf),
        )
        pose_mol.SetProp("MMGBSA_DOCKED_PDBQT", str(docked_path))
        pose_mol.SetProp("MMGBSA_SOURCE_SDF", str(source_sdf))

        out_sdf.parent.mkdir(parents=True, exist_ok=True)
        if not write_rdkit_mol_to_sdf(pose_mol, out_sdf, set_kekulize_false=False):
            return PoseSdfResult(False, "sdf_write_failed")

        payload = {
            "schema_version": 1,
            "created_by": "atlas_mmgbsa_pose_artifacts",
            "docked_pdbqt": str(docked_path),
            "ligand_pdbqt": str(ligand_path),
            "source_sidecar": str(sidecar_for_meta) if sidecar_for_meta else "",
            "source_sdf": str(source_sdf),
            "source_record_index": record_index,
            "source_record_name": source_meta.get("source_record_name", ""),
            "output_sdf": str(out_sdf),
            "chemistry_authoritative": authoritative,
            "source_kind": source_meta.get("source_kind", ""),
            "bond_orders_restored": source_meta.get("bond_orders_restored"),
            "coordinates_from": "docked_pdbqt_first_model",
            "bond_orders_charges_from": "source_sdf",
            "hydrogen_coordinates_inferred": bool(inferred_h),
            "pdbqt_to_sdf_reconstruction": False,
        }
        tmp_meta = out_meta.with_suffix(out_meta.suffix + ".part")
        tmp_meta.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(tmp_meta, out_meta)
        log.info(
            "[mmgbsa.pose_sdf] status=ok docked=%s sdf=%s source=%s inferred_h=%s",
            docked_path,
            out_sdf,
            source_sdf,
            inferred_h,
        )
        return PoseSdfResult(True, "ok", out_sdf, out_meta)
    except MMGBSA_POSE_EXPORT_ERRORS as exc:
        log.warning(
            "[mmgbsa.pose_sdf] action=pose_export_failed status=failed docked=%s ligand=%s reason=%s",
            docked_path,
            ligand_path,
            exc,
            exc_info=True,
        )
        return PoseSdfResult(False, f"failed:{exc}")


__all__ = ["PoseSdfResult", "preserve_authoritative_pose_sdf"]
