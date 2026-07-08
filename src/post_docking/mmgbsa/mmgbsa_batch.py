from __future__ import annotations

import csv
import json
import logging
import os
import shutil
import stat
import subprocess
import tarfile
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

from post_docking.mmgbsa.mmgbsa_pipeline import _maybe_run_mmgbsa_for_pdb
from post_docking.mmgbsa.rdkit_sdf_io import write_rdkit_mol_to_sdf as _write_rdkit_mol_to_sdf

LOGGER = logging.getLogger("mmgbsa.batch")


@dataclass(frozen=True)
class MMGBSAHit:
    run_id: str
    pdb_id: str
    pdb_file: str
    variant: Optional[str]
    ph_label: Optional[str]
    stage_dir: str
    ligand_id: str
    source_sdf_path: str = ""
    source_mol2_path: str = ""


@dataclass(frozen=True)
class MMGBSAJob:
    run_id: str
    pdb_id: str
    pdb_file: str
    variant: Optional[str]
    ph_label: Optional[str]
    stage_dir: str
    ligand_ids: tuple[str, ...]
    ligand_sources: tuple[tuple[str, str, str], ...] = ()


def _first_value(row: Mapping[str, Any], names: Sequence[str]) -> str:
    lowered = {str(k).strip().lower(): v for k, v in row.items()}
    for name in names:
        value = lowered.get(name.lower())
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return text
    return ""


def _as_bool(value: Any) -> Optional[bool]:
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "significant", "hit"}:
        return True
    if text in {"0", "false", "no", "n", "not_significant"}:
        return False
    return None


def _as_float(value: Any) -> Optional[float]:
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def _normalize_variant(raw: str) -> Optional[str]:
    token = str(raw or "").strip()
    if not token or token.lower() in {"base", "legacy", "none", "null", "-"}:
        return None
    return token.upper()


def _normalize_ph(raw: str) -> Optional[str]:
    token = str(raw or "").strip()
    if not token or token.lower() in {"base", "none", "null", "-"}:
        return None
    return token


def _row_is_significant(
    row: Mapping[str, Any], require_significant: bool, max_q_value: Optional[float]
) -> bool:
    if not require_significant and max_q_value is None:
        return True

    sig_raw = _first_value(
        row,
        (
            "significant",
            "is_significant",
            "heatmap_significant",
            "pair_significant",
            "selected",
        ),
    )
    sig_val = _as_bool(sig_raw) if sig_raw else None
    if require_significant and sig_val is False:
        return False
    if require_significant and sig_val is True and max_q_value is None:
        return True

    q_raw = _first_value(row, ("q_value", "q", "fdr", "adj_p_value", "p_value", "p"))
    q_val = _as_float(q_raw)
    if max_q_value is not None:
        return q_val is not None and q_val <= max_q_value
    return sig_val is True


def load_report_hits(
    report_csv: Path,
    *,
    run_id: str,
    default_stage_dir: str = "stage1",
    default_variant: Optional[str] = None,
    default_ph_label: Optional[str] = None,
    require_significant: bool = False,
    max_q_value: Optional[float] = None,
    score_column: str = "",
    score_direction: str = "lower",
    top_n: int = 0,
) -> list[MMGBSAHit]:
    rows: list[dict[str, str]] = []
    with report_csv.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if _row_is_significant(row, require_significant, max_q_value):
                rows.append(dict(row))

    if score_column:
        reverse = score_direction.strip().lower() in {"higher", "desc", "descending"}
        missing_sentinel = float("-inf") if reverse else float("inf")

        def _score_key(row: Mapping[str, Any]) -> float:
            value = _as_float(row.get(score_column))
            return value if value is not None else missing_sentinel

        rows.sort(key=_score_key, reverse=reverse)
    if top_n > 0:
        rows = rows[:top_n]

    hits: list[MMGBSAHit] = []
    seen: set[tuple[str, str, Optional[str], Optional[str], str, str]] = set()
    for row in rows:
        pdb_id = _first_value(row, ("pdb_id", "pdb", "target_pdb", "structure_id"))
        ligand_id = _first_value(
            row,
            (
                "ligand_base",
                "ligand_id",
                "ligand",
                "drug_id",
                "compound_id",
                "molecule_id",
            ),
        )
        if not pdb_id or not ligand_id:
            LOGGER.warning("[mmgbsa.batch] row=skip reason=missing_pdb_or_ligand")
            continue

        row_run_id = _first_value(row, ("run_id", "run", "atlas_run_id")) or run_id
        variant = _normalize_variant(
            _first_value(row, ("variant", "variant_label", "apo_holo", "structure_variant"))
            or (default_variant or "")
        )
        ph_label = _normalize_ph(
            _first_value(row, ("ph_label", "ph_tag", "ph", "ph_dir"))
            or (default_ph_label or "")
        )
        stage_dir = (
            _first_value(row, ("stage_dir", "stage", "input_stage_dir"))
            or default_stage_dir
        )
        pdb_file = _first_value(row, ("pdb_file", "pdb_path")) or f"{pdb_id}.pdb"
        source_sdf = _first_value(
            row, ("sdf_path", "ligand_sdf", "source_sdf", "pdbbind_sdf")
        )
        source_mol2 = _first_value(
            row, ("mol2_path", "ligand_mol2", "source_mol2", "pdbbind_mol2")
        )
        key = (
            row_run_id,
            pdb_id.upper(),
            variant,
            ph_label,
            stage_dir,
            ligand_id,
        )
        if key in seen:
            continue
        seen.add(key)
        hits.append(
            MMGBSAHit(
                run_id=row_run_id,
                pdb_id=pdb_id.upper(),
                pdb_file=pdb_file,
                variant=variant,
                ph_label=ph_label,
                stage_dir=stage_dir,
                ligand_id=ligand_id,
                source_sdf_path=source_sdf,
                source_mol2_path=source_mol2,
            )
        )
    return hits


def group_hits(hits: Iterable[MMGBSAHit]) -> list[MMGBSAJob]:
    grouped: dict[tuple[str, str, str, Optional[str], Optional[str], str], set[str]] = {}
    sources: dict[
        tuple[str, str, str, Optional[str], Optional[str], str], dict[str, tuple[str, str]]
    ] = {}
    for hit in hits:
        key = (
            hit.run_id,
            hit.pdb_id,
            hit.pdb_file,
            hit.variant,
            hit.ph_label,
            hit.stage_dir,
        )
        grouped.setdefault(key, set()).add(hit.ligand_id)
        if hit.source_sdf_path or hit.source_mol2_path:
            sources.setdefault(key, {})[hit.ligand_id] = (
                hit.source_sdf_path,
                hit.source_mol2_path,
            )

    jobs: list[MMGBSAJob] = []
    def _group_sort_key(
        item: tuple[tuple[str, str, str, Optional[str], Optional[str], str], set[str]]
    ) -> tuple[str, str, str, str, str, str]:
        run_id, pdb_id, pdb_file, variant, ph_label, stage_dir = item[0]
        return (
            run_id,
            pdb_id,
            pdb_file,
            variant or "",
            ph_label or "",
            stage_dir,
        )

    for (run_id, pdb_id, pdb_file, variant, ph_label, stage_dir), ligand_ids in sorted(
        grouped.items(), key=_group_sort_key
    ):
        group_key = (run_id, pdb_id, pdb_file, variant, ph_label, stage_dir)
        source_map = sources.get(group_key, {})
        jobs.append(
            MMGBSAJob(
                run_id=run_id,
                pdb_id=pdb_id,
                pdb_file=pdb_file,
                variant=variant,
                ph_label=ph_label,
                stage_dir=stage_dir,
                ligand_ids=tuple(sorted(ligand_ids)),
                ligand_sources=tuple(
                    (ligand_id, source_sdf, source_mol2)
                    for ligand_id, (source_sdf, source_mol2) in sorted(
                        source_map.items()
                    )
                ),
            )
        )
    return jobs


def _safe_member_path(dest_dir: Path, member_name: str) -> Path:
    target = (dest_dir / member_name).resolve()
    dest_resolved = dest_dir.resolve()
    if target != dest_resolved and dest_resolved not in target.parents:
        raise ValueError(f"archive member escapes destination: {member_name}")
    return target


def _normalize_artifact_member_name(member_name: str, dest_dir: Path) -> str:
    normalized = str(member_name or "").replace("\\", "/")
    if (
        normalized.startswith("/")
        or normalized == ".."
        or normalized.startswith("../")
        or "/../" in normalized
    ):
        return normalized
    parts = tuple(part for part in normalized.split("/") if part not in {"", "."})
    if len(parts) >= 3 and parts[0] in {"post_docked", "docked"}:
        run_id = dest_dir.name
        if parts[1] == run_id:
            return str(Path(*parts[2:]))
    return str(Path(*parts)) if parts else ""


def _extract_tar_members(archive: tarfile.TarFile, dest_dir: Path) -> None:
    for member in archive:
        member_name = _normalize_artifact_member_name(member.name, dest_dir)
        if not member_name:
            continue
        target = _safe_member_path(dest_dir, member_name)
        if member.isdir():
            target.mkdir(parents=True, exist_ok=True)
            continue
        if not member.isfile():
            raise ValueError(f"unsupported archive member type: {member.name}")
        source = archive.extractfile(member)
        if source is None:
            raise ValueError(f"unable to read archive member: {member.name}")
        target.parent.mkdir(parents=True, exist_ok=True)
        tmp_target = target.with_suffix(target.suffix + ".part")
        with source, tmp_target.open("wb") as handle:
            shutil.copyfileobj(source, handle)
        os.replace(tmp_target, target)


def extract_artifact(archive_path: Path, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    suffixes = "".join(archive_path.suffixes).lower()
    if suffixes.endswith(".tar.zst"):
        proc = subprocess.Popen(
            ["zstd", "-dc", str(archive_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert proc.stdout is not None
        try:
            with tarfile.open(fileobj=proc.stdout, mode="r|") as archive:
                _extract_tar_members(archive, dest_dir)
        finally:
            proc.stdout.close()
        _, stderr = proc.communicate()
        if proc.returncode != 0:
            detail = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(f"zstd extraction failed for {archive_path}: {detail}")
    elif suffixes.endswith((".tar.gz", ".tgz", ".tar.bz2", ".tar.xz", ".tar")):
        with tarfile.open(archive_path, "r:*") as archive:
            _extract_tar_members(archive, dest_dir)
    elif archive_path.suffix.lower() == ".zip":
        with zipfile.ZipFile(archive_path) as archive:
            for info in archive.infolist():
                member_name = _normalize_artifact_member_name(info.filename, dest_dir)
                if not member_name:
                    continue
                target = _safe_member_path(dest_dir, member_name)
                if info.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                file_type = (info.external_attr >> 16) & 0o170000
                if file_type and file_type != stat.S_IFREG:
                    raise ValueError(
                        f"unsupported archive member type: {info.filename}"
                    )
                target.parent.mkdir(parents=True, exist_ok=True)
                tmp_target = target.with_suffix(target.suffix + ".part")
                with archive.open(info) as source, tmp_target.open("wb") as handle:
                    shutil.copyfileobj(source, handle)
                os.replace(tmp_target, target)
    else:
        raise ValueError(f"unsupported artifact archive: {archive_path}")
    return dest_dir


def _safe_ligand_stem(ligand_id: str) -> str:
    text = Path(str(ligand_id or "ligand")).stem
    safe = "".join(ch if ch.isalnum() or ch in {"_", "-", "."} else "_" for ch in text)
    return safe.strip("._") or "ligand"


def _post_docked_base(cfg: Mapping[str, Any], job: MMGBSAJob) -> Path:
    root = Path(str(cfg.get("POST_DOCKED_DIR", "post_docked") or "post_docked")).expanduser()
    if not root.is_absolute():
        root = Path(str(cfg.get("OVERALL_DIR", "."))) / root
    base = root / job.run_id / job.pdb_id
    if job.variant:
        base = base / job.variant.upper()
    return base


def _sdf_sanitizes(path: Path) -> bool:
    try:
        from rdkit import Chem

        supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
        return bool(supplier and len(supplier) and supplier[0] is not None)
    except (ImportError, OSError, RuntimeError, ValueError):
        return False


def _write_sdf_from_mol2(mol2_path: Path, sdf_path: Path) -> bool:
    try:
        from rdkit import Chem

        mol = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
        if mol is None:
            return False
        return _write_rdkit_mol_to_sdf(mol, sdf_path)
    except (ImportError, OSError, RuntimeError, ValueError):
        return False


def _write_sdf_from_mol2_with_sdf_coords(
    mol2_path: Path,
    coord_sdf_path: Path,
    sdf_path: Path,
) -> bool:
    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D

        template = Chem.MolFromMol2File(str(mol2_path), sanitize=True, removeHs=False)
        if template is None:
            template = Chem.MolFromMol2File(
                str(mol2_path), sanitize=False, removeHs=False
            )
        if template is None:
            return False
        supplier = Chem.SDMolSupplier(
            str(coord_sdf_path), sanitize=False, removeHs=False
        )
        coord_mol = supplier[0] if supplier and len(supplier) else None
        if coord_mol is None or coord_mol.GetNumConformers() == 0:
            return False

        coord_conf = coord_mol.GetConformer()
        coord_heavy = [
            coord_conf.GetAtomPosition(atom.GetIdx())
            for atom in coord_mol.GetAtoms()
            if atom.GetAtomicNum() > 1
        ]
        template_heavy = [
            int(atom.GetIdx()) for atom in template.GetAtoms() if atom.GetAtomicNum() > 1
        ]
        if len(coord_heavy) != len(template_heavy):
            return False

        old_conf = template.GetConformer() if template.GetNumConformers() else None
        new_positions: dict[int, tuple[float, float, float]] = {
            atom_idx: (float(pos.x), float(pos.y), float(pos.z))
            for atom_idx, pos in zip(template_heavy, coord_heavy)
        }
        for atom in template.GetAtoms():
            atom_idx = int(atom.GetIdx())
            if atom_idx in new_positions:
                continue
            inferred = None
            for neighbor in atom.GetNeighbors():
                nbr_idx = int(neighbor.GetIdx())
                if nbr_idx not in new_positions:
                    continue
                base = new_positions[nbr_idx]
                if old_conf is None:
                    inferred = base
                    break
                old_atom = old_conf.GetAtomPosition(atom_idx)
                old_nbr = old_conf.GetAtomPosition(nbr_idx)
                inferred = (
                    base[0] + float(old_atom.x - old_nbr.x),
                    base[1] + float(old_atom.y - old_nbr.y),
                    base[2] + float(old_atom.z - old_nbr.z),
                )
                break
            new_positions[atom_idx] = inferred or (0.0, 0.0, 0.0)

        conf = Chem.Conformer(template.GetNumAtoms())
        for atom_idx in range(template.GetNumAtoms()):
            x, y, z = new_positions[atom_idx]
            conf.SetAtomPosition(atom_idx, Point3D(x, y, z))
        conf.Set3D(True)
        template.RemoveAllConformers()
        template.AddConformer(conf, assignId=True)
        return _write_rdkit_mol_to_sdf(template, sdf_path)
    except (ImportError, OSError, RuntimeError, ValueError, IndexError):
        return False


def _write_source_metadata(
    sdf_path: Path,
    *,
    ligand_id: str,
    source_sdf: Path | None,
    source_mol2: Path | None,
    source_used: str,
) -> None:
    payload = {
        "schema": "mmgbsa_report_source_pose_v1",
        "ligand_id": ligand_id,
        "input_chemistry_authoritative": True,
        "chemistry_authoritative": True,
        "source_kind": "report_supplied_pdbbind_or_curated_ligand",
        "source_sdf": str(source_sdf) if source_sdf else "",
        "source_mol2": str(source_mol2) if source_mol2 else "",
        "source_used": source_used,
        "pdbqt_to_sdf_reconstruction": False,
        "notes": "Report-selected ligand source was materialized for MMGBSA.",
    }
    meta_path = sdf_path.with_suffix(".mmgbsa_pose.json")
    tmp_path = meta_path.with_suffix(meta_path.suffix + ".part")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp_path, meta_path)


def _existing_source_paths(
    source_sdf_raw: str, source_mol2_raw: str
) -> tuple[Path | None, Path | None]:
    source_sdf = Path(source_sdf_raw).expanduser() if source_sdf_raw else None
    source_mol2 = Path(source_mol2_raw).expanduser() if source_mol2_raw else None
    if source_sdf is not None and not source_sdf.is_file():
        source_sdf = None
    if source_mol2 is not None and not source_mol2.is_file():
        source_mol2 = None
    return source_sdf, source_mol2


def _copy_report_ligand_sources(
    stage_dir: Path, stem: str, source_sdf: Path | None, source_mol2: Path | None
) -> tuple[Path, str]:
    dest_sdf = stage_dir / f"{stem}.sdf"
    source_used = "sdf"
    if source_sdf is not None and source_mol2 is not None:
        if _write_sdf_from_mol2_with_sdf_coords(source_mol2, source_sdf, dest_sdf):
            source_used = "mol2_template_sdf_coordinates"
        elif source_sdf.resolve() != dest_sdf.resolve():
            shutil.copy2(source_sdf, dest_sdf)
    elif source_sdf is not None:
        if source_sdf.resolve() != dest_sdf.resolve():
            shutil.copy2(source_sdf, dest_sdf)
    if source_mol2 is not None:
        dest_mol2 = stage_dir / f"{stem}.source.mol2"
        if source_mol2.resolve() != dest_mol2.resolve():
            shutil.copy2(source_mol2, dest_mol2)
    if (source_sdf is None or not _sdf_sanitizes(dest_sdf)) and source_mol2 is not None:
        if _write_sdf_from_mol2(source_mol2, dest_sdf):
            source_used = "mol2_to_sdf"
    return dest_sdf, source_used


def _materialize_one_report_ligand(
    stage_dir: Path,
    ligand_id: str,
    source_sdf: Path | None,
    source_mol2: Path | None,
) -> bool:
    if source_sdf is None and source_mol2 is None:
        return False
    dest_sdf, source_used = _copy_report_ligand_sources(
        stage_dir, _safe_ligand_stem(ligand_id), source_sdf, source_mol2
    )
    if not _materialized_sdf_ready(dest_sdf):
        return False
    _write_source_metadata(
        dest_sdf,
        ligand_id=ligand_id,
        source_sdf=source_sdf,
        source_mol2=source_mol2,
        source_used=source_used,
    )
    return True


def _materialized_sdf_ready(dest_sdf: Path) -> bool:
    return dest_sdf.exists() and dest_sdf.stat().st_size > 0


def _materialize_report_ligand_sources(job: MMGBSAJob, cfg: Mapping[str, Any]) -> bool:
    if not job.ligand_sources:
        return False
    ph_label = job.ph_label or "pH7_0"
    stage_dir = _post_docked_base(cfg, job) / ph_label / job.stage_dir
    stage_dir.mkdir(parents=True, exist_ok=True)
    materialized = [
        _materialize_one_report_ligand(
            stage_dir,
            ligand_id,
            *_existing_source_paths(source_sdf_raw, source_mol2_raw),
        )
        for ligand_id, source_sdf_raw, source_mol2_raw in job.ligand_sources
    ]
    return any(materialized)


def run_jobs(
    jobs: Sequence[MMGBSAJob],
    *,
    cfg: Mapping[str, Any],
    test_mode: str = "off",
    dry_run: bool = False,
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for job in jobs:
        job_cfg = dict(cfg)
        job_cfg.update(
            {
                "MMGBSA_ENABLED": True,
                "MMGBSA_INPUT_STAGE_DIR": job.stage_dir,
                "MMGBSA_SELECTED_LIGANDS": ",".join(job.ligand_ids),
                "MMGBSA_MAX_LIGANDS": 0,
            }
        )
        if job.ph_label:
            job_cfg["MMGBSA_INPUT_PH_LABELS"] = job.ph_label
        materialized_sources = False
        if not dry_run:
            materialized_sources = _materialize_report_ligand_sources(job, job_cfg)
        if materialized_sources:
            job_cfg["MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE"] = True
            job_cfg["MMGBSA_INPUT_CHEMISTRY_SOURCE"] = "report_supplied_pdbbind_or_curated_ligand"
        payload = {
            "job": asdict(job),
            "dry_run": dry_run,
            "planned_source_materialization": bool(job.ligand_sources),
            "materialized_sources": materialized_sources,
            "status": "planned" if dry_run else "started",
        }
        results.append(payload)
        if dry_run:
            continue
        _maybe_run_mmgbsa_for_pdb(
            job_cfg,
            job.pdb_file,
            job.pdb_id,
            job.variant,
            job.run_id,
            test_mode,
            job.variant is None,
        )
        payload["status"] = "submitted"
    return results


def write_plan(path: Path, results: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".part")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(list(results), handle, indent=2, sort_keys=True)
    os.replace(tmp_path, path)
