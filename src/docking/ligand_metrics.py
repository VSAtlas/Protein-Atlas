from __future__ import annotations

import contextlib
import logging
import math
import os
import tempfile
from pathlib import Path
from typing import Dict, Optional

from rdkit import Chem, rdBase
from rdkit.Chem import rdFMCS, rdMolAlign

from docking.pose_validation import compute_redock_rmsd
from path_router.path_router import make_paths


def _rdkit_quiet_logs():
    try:
        return rdBase.BlockLogs()
    except Exception:
        return contextlib.nullcontext()


def compute_ligand_efficiency(
    score: Optional[float], heavy_atoms: Optional[int]
) -> Optional[float]:
    try:
        if score is None or heavy_atoms is None or int(heavy_atoms) <= 0:
            return None
        return float(-float(score) / int(heavy_atoms))
    except Exception:
        return None


def record_le(
    score_history: Dict[str, Dict[str, Dict]],
    stage_name: str,
    lig_path: str,
    score: Optional[float],
    heavy_atom_counts: Dict[str, int],
) -> Optional[float]:
    ha = heavy_atom_counts.get(lig_path)
    le = compute_ligand_efficiency(score, ha)
    stage_map = score_history.setdefault(stage_name, {})
    rec = stage_map.setdefault(lig_path, {})
    rec["heavy_atoms"] = int(ha) if isinstance(ha, (int, float)) else None
    rec["le"] = le
    return le


def write_gnina_scores_csv(*args, **kwargs):
    from docking.docking_gnina_support import write_gnina_scores_csv as _impl

    return _impl(*args, **kwargs)


def _pose_path_for(
    csv_cfg: Dict,
    pdb_id: str,
    stage_name: str,
    lig_path: str,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
) -> str:
    """Build the expected pose path for a ligand at a given stage."""
    from pathlib import Path

    paths = make_paths(csv_cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_token = (
        variant or os.environ.get("APO_HOLO_VARIANT", "") or ""
    ).strip().upper() or None
    stage_dir = paths.docked_stage_dir(variant_token, stage_name, ph_token)
    return str(stage_dir / f"{Path(lig_path).stem}_{stage_name}.pdbqt")


def write_scores_csv(
    cfg: Dict,
    pdb_id: str,
    score_history: Dict[str, Dict[str, Dict]],
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
    *,
    csv_prefix: str = "",
    summary_basename: str = "docking_score_summary.csv",
    long_basename: str = "docking_score_long.csv",
) -> str:
    from docking.docking_vina import write_scores_csv as _impl

    return _impl(
        cfg,
        pdb_id,
        score_history,
        ph_label=ph_label,
        variant=variant,
        csv_prefix=csv_prefix,
        summary_basename=summary_basename,
        long_basename=long_basename,
    )


def _read_any_lig(path: str):
    mol = None
    loader = "unknown"
    sanitize = True
    sanitize_failed = False
    load_err: object | None = None
    ext = os.path.splitext(path)[1].lower()
    rmsd_logger = logging.getLogger("rmsd")

    try:
        with _rdkit_quiet_logs():
            if ext in (".sdf", ".sd"):
                loader = "SDMolSupplier"
                suppl = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
                mol = next((m for m in suppl if m is not None), None)
            elif ext == ".mol2":
                loader = "MolFromMol2File"
                mol = Chem.MolFromMol2File(path, sanitize=True, removeHs=False)
            elif ext == ".pdbqt":
                loader = "MolFromPDBFile(pdbqt)"
                sanitize = False
                mol = Chem.MolFromPDBFile(
                    path, sanitize=False, removeHs=False, proximityBonding=True
                )
                if mol is not None:
                    try:
                        Chem.SanitizeMol(mol)
                    except Exception as exc:
                        sanitize_failed = True
                        rmsd_logger.warning(
                            "[read_any] sanitize failed for path='%s' err=%r",
                            path,
                            exc,
                        )
            elif ext == ".pdb":
                loader = "MolFromPDBFile"
                try:
                    mol = Chem.MolFromPDBFile(
                        path, sanitize=True, removeHs=False, proximityBonding=True
                    )
                except Exception as exc:
                    load_err = exc
                    mol = None
                if mol is None:
                    try:
                        with open(path, "rt", errors="ignore") as fh:
                            head = fh.read(1024)
                    except Exception:
                        head = ""
                    looks_like_pdbqt = any(
                        token in head for token in ("REMARK VINA", "TORSDOF", "ROOT")
                    )
                    if looks_like_pdbqt:
                        loader = "MolFromPDBFile(pdbqt-fallback)"
                        sanitize = False
                        try:
                            mol = Chem.MolFromPDBFile(
                                path,
                                sanitize=False,
                                removeHs=False,
                                proximityBonding=True,
                            )
                        except Exception as exc:
                            load_err = exc
                            mol = None
                        if mol is not None:
                            try:
                                Chem.SanitizeMol(mol)
                            except Exception as exc:
                                sanitize_failed = True
                                rmsd_logger.warning(
                                    "[read_any] sanitize failed for path='%s' err=%r",
                                    path,
                                    exc,
                                )
            else:
                loader = "auto"
                mol = Chem.MolFromMolFile(path, sanitize=True, removeHs=False)
    except Exception as exc:
        load_err = exc if load_err is None else load_err
        mol = None
    if (mol is None) and (load_err is None):
        load_err = "load_returned_None"

    try:
        if mol is not None:
            from rdkit.Chem import rdMolDescriptors

            try:
                from rdkit.Chem import inchi

                inchikey = inchi.MolToInchiKey(mol)
            except Exception:
                inchikey = "NA"
            extra = " sanitize_failed=True" if sanitize_failed else ""
            rmsd_logger.info(
                "[read_any] loader=%s sanitize=%s path='%s'%s atoms=%s heavy=%s formula=%s inchikey=%s",
                loader,
                sanitize,
                path,
                extra,
                mol.GetNumAtoms(),
                mol.GetNumHeavyAtoms(),
                rdMolDescriptors.CalcMolFormula(mol),
                inchikey,
            )
        else:
            rmsd_logger.info(
                "[read_any] loader=%s sanitize=%s path='%s' mol=None err=%r",
                loader,
                sanitize,
                path,
                load_err,
            )
    except Exception:
        pass
    return mol


def _is_readable_ref(pth: Path) -> bool:
    try:
        mol = _read_any_lig(str(pth))
        return (mol is not None) and (mol.GetNumHeavyAtoms() > 0)
    except Exception:
        return False


def compute_rmsd(ref_path: str, docked_path: str) -> float:
    rmsd_logger = logging.getLogger("rmsd")

    def _clamp(val: float) -> float:
        try:
            if 0 < val < 0.01:
                return 0.01
        except Exception:
            pass
        return val

    coord_rmsd = None
    tmp_to_cleanup: list[Path] = []
    try:
        coord_ref = ref_path
        if not str(ref_path).lower().endswith(".pdb"):
            mol = _read_any_lig(ref_path)
            if mol is not None:
                tmp_ref = Path(tempfile.mkstemp(suffix=".pdb")[1])
                Chem.MolToPDBFile(mol, str(tmp_ref))
                coord_ref = str(tmp_ref)
                tmp_to_cleanup.append(tmp_ref)

        coord_rmsd = compute_redock_rmsd(coord_ref, docked_path)
        if coord_rmsd is not None and math.isfinite(coord_rmsd):
            rmsd_logger.info(
                "[rmsd.coord] ref='%s' dock='%s' rmsd=%.3fA (kabsch)",
                coord_ref,
                docked_path,
                coord_rmsd,
            )
            return _clamp(coord_rmsd)
    except Exception as exc:
        rmsd_logger.warning(
            "[rmsd.coord] failed ref='%s' dock='%s' err=%r; falling back to RDKit/MCS",
            ref_path,
            docked_path,
            exc,
        )
    finally:
        for tmp in tmp_to_cleanup:
            try:
                tmp.unlink()
            except Exception:
                pass

    ref = _read_any_lig(ref_path)
    dock = _read_any_lig(docked_path)
    if not ref or not dock:
        return float("inf")

    try:
        rmsd = float(rdMolAlign.GetBestRMS(ref, dock))
        rmsd_logger.info(
            "[rmsd.best] ref='%s' dock='%s' rmsd=%.3fA",
            ref_path,
            docked_path,
            rmsd,
        )
        return _clamp(rmsd)
    except Exception as exc:
        rmsd_logger.warning(
            "[rmsd.best] failed ref='%s' dock='%s' err=%r; falling back to MCS",
            ref_path,
            docked_path,
            exc,
        )

    try:
        Chem.FastFindRings(ref)
        Chem.FastFindRings(dock)
        mcs = rdFMCS.FindMCS(
            [ref, dock],
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            matchValences=True,
        )
        patt = Chem.MolFromSmarts(mcs.smartsString)
        if patt is None:
            return float("inf")
        ref_match = ref.GetSubstructMatch(patt)
        dock_match = dock.GetSubstructMatch(patt)
        if not ref_match or not dock_match or (len(ref_match) != len(dock_match)):
            return float("inf")
        atom_map = list(zip(dock_match, ref_match))
        rmsd = float(rdMolAlign.AlignMol(dock, ref, atomMap=atom_map))
        rmsd_logger.info(
            "[rmsd.mcs] ref='%s' dock='%s' rmsd=%.3fA atoms=%d",
            ref_path,
            docked_path,
            rmsd,
            len(atom_map),
        )
        return _clamp(rmsd)
    except Exception as exc:
        rmsd_logger.warning(
            "[rmsd.mcs] failed ref='%s' dock='%s' err=%r",
            ref_path,
            docked_path,
            exc,
        )
        return float("inf")
