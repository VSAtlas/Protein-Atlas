"""Bulk SDF ligand preparation through the shared Meeko backend."""

from __future__ import annotations

import concurrent.futures
from dataclasses import dataclass
from datetime import datetime
import json
import logging
import os
from pathlib import Path
import shutil
from typing import Dict, Iterable, Optional, Set

from rdkit import Chem

from prep_ligands.prep_ligands_meeko import prepare_mol_to_pdbqt_with_meeko
from prep_ligands.prep_ligands_microstates import compute_microstate_id
from prep_ligands.prep_ligands_reporting import _append_prep_status
from prep_ligands.prep_ligands_runtime import (
    BulkContext,
    BulkInitArgs,
    init_bulk_context,
    ph_label,
    relative_to_output,
)
from prep_ligands.prep_ligands_common import is_valid_ligand


logger = logging.getLogger(__name__)


@dataclass(slots=True)
class _MeekoPrepareTask:
    lig_stem: str
    ph_label: str
    ph_value: float
    pdbqt_path: Path
    pdbqt_rel: str
    source_sdf: Path
    source_record_index: int
    source_record_name: str


def _record_name(mol: Chem.Mol, sdf_path: Path, record_index: int, total_hint: int) -> str:
    raw_name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
    if raw_name:
        return raw_name
    if total_hint <= 1:
        return sdf_path.stem
    return f"{sdf_path.stem}_{record_index:05d}"


def _safe_stem(ctx: BulkContext, raw_name: str, seq: int) -> str:
    if ctx.rename_active:
        return f"{ctx.rename_prefix}{seq:0{ctx.rename_pad}d}"
    from prep_ligands.prep_ligands_runtime import sanitize_ligand_name_for_filename

    return sanitize_ligand_name_for_filename(raw_name)


def _matches_only(
    only_set: Optional[Set[str]],
    *,
    raw_name: str,
    safe_stem: str,
    sdf_path: Path,
    record_index: int,
) -> bool:
    if not only_set:
        return True
    candidates = {
        raw_name,
        safe_stem,
        sdf_path.stem,
        f"{sdf_path.stem}_{record_index:05d}",
    }
    return bool(candidates & only_set)


def _discover_sdf_files(ctx: BulkContext) -> list[Path]:
    sdf_files = (
        [Path(ctx.in_sdf_env).resolve()]
        if ctx.in_sdf_env
        else list(ctx.ligand_extracted_dir.glob("*.sdf"))
    )
    initial_count = len(sdf_files)
    sdf_files = [path for path in sdf_files if "docked" not in path.name.lower()]
    skipped = initial_count - len(sdf_files)
    if skipped > 0:
        logger.info(
            "[ligprep.skip] skipped %d SDF(s) with 'docked' in name under %s",
            skipped,
            ctx.ligand_extracted_dir,
        )
    print(f"Found {len(sdf_files)} SDF file(s)")
    return sdf_files


def _count_sdf_records(sdf_path: Path) -> int:
    try:
        count = 0
        with open(sdf_path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("$$$$"):
                    count += 1
        return max(1, count)
    except Exception:
        return 1


def _iter_sdf_molecules(
    sdf_files: Iterable[Path],
    ctx: BulkContext,
    only_set: Optional[Set[str]],
) -> list[tuple[Chem.Mol, Path, str, str, int]]:
    records: list[tuple[Chem.Mol, Path, str, str, int]] = []
    seq = ctx.rename_start
    for sdf_path in sdf_files:
        total_hint = _count_sdf_records(sdf_path)
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        for record_index, mol in enumerate(supplier, start=1):
            if mol is None:
                logger.warning(
                    "[ligprep.meeko] skipped unreadable SDF record file=%s index=%d",
                    sdf_path,
                    record_index,
                )
                continue
            raw_name = _record_name(mol, sdf_path, record_index, total_hint)
            safe_stem = _safe_stem(ctx, raw_name, seq)
            if not _matches_only(
                only_set,
                raw_name=raw_name,
                safe_stem=safe_stem,
                sdf_path=sdf_path,
                record_index=record_index,
            ):
                continue
            seq += 1
            mol.SetProp("_Name", raw_name)
            records.append((mol, sdf_path, raw_name, safe_stem, record_index))
    return records


def _pdbqt_name_for_task(ctx: BulkContext, lig_stem: str, current_ph: float) -> tuple[Path, str, str]:
    current_ph_label = ph_label(current_ph) if ctx.use_ph_subdirs else ""
    out_dir = ctx.output_ligands_dir / current_ph_label if current_ph_label else ctx.output_ligands_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if ctx.use_ph_subdirs and current_ph_label:
        base, _, rest = lig_stem.partition("_")
        lig_stem_ph = f"{base}{current_ph_label}_{rest}" if rest else f"{lig_stem}{current_ph_label}"
        pdbqt_name = f"{lig_stem_ph}.pdbqt"
    else:
        pdbqt_name = f"{lig_stem}.pdbqt"
    pdbqt_path = out_dir / pdbqt_name
    return pdbqt_path, relative_to_output(pdbqt_path, ctx.output_ligands_dir), current_ph_label


def _update_microstate_registry(ctx: BulkContext, task: _MeekoPrepareTask) -> None:
    if (
        not ctx.microstate_dedup
        or ctx.ph_values is None
        or ctx.microstates.directory is None
        or ctx.microstates.registry is None
        or ctx.microstates.index is None
        or not task.pdbqt_path.is_file()
    ):
        return

    microstate_id = compute_microstate_id(task.pdbqt_path)
    if not microstate_id:
        logger.debug(
            "[microstate] skip_dedup_empty_id ligand=%s ph_label=%s ph_value=%.2f path=%s",
            task.lig_stem,
            task.ph_label,
            task.ph_value,
            task.pdbqt_path,
        )
        return

    entry = ctx.microstates.index.get(microstate_id)
    if entry is None:
        canonical_name = f"{task.lig_stem}__ms_{microstate_id}.pdbqt"
        canonical_rel_path = f"microstates/{canonical_name}"
        canonical_path = ctx.microstates.directory / canonical_name
        try:
            canonical_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(task.pdbqt_path, canonical_path)
            _copy_source_sidecar(task.pdbqt_path, canonical_path)
        except Exception as exc:
            logger.warning(
                "[microstate] copy_to_canonical_failed ligand=%s src=%s dst=%s err=%s",
                task.lig_stem,
                task.pdbqt_path,
                canonical_path,
                exc,
            )
        entry = {
            "microstate_id": microstate_id,
            "pdbqt_path": canonical_rel_path,
            "aliases": [],
        }
        ctx.microstates.registry["microstates"].append(entry)
        ctx.microstates.index[microstate_id] = entry
        ctx.microstates.dirty = True
        logger.debug(
            "[microstate] new_microstate library=%s ligand=%s microstate_id=%s canonical=%s",
            ctx.library_name,
            task.lig_stem,
            microstate_id,
            canonical_name,
        )
    else:
        canonical_rel_path = entry.get("pdbqt_path")
        if canonical_rel_path:
            canonical_path = ctx.prepped_ligands_dir / canonical_rel_path
            if not canonical_path.exists():
                try:
                    canonical_path.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(task.pdbqt_path, canonical_path)
                    _copy_source_sidecar(task.pdbqt_path, canonical_path)
                except Exception as exc:
                    logger.warning(
                        "[microstate] canonical_missing_copy_failed ligand=%s src=%s dst=%s err=%s",
                        task.lig_stem,
                        task.pdbqt_path,
                        canonical_path,
                        exc,
                    )
        logger.debug(
            "[microstate] reuse_microstate library=%s ligand=%s microstate_id=%s canonical=%s",
            ctx.library_name,
            task.lig_stem,
            microstate_id,
            canonical_rel_path,
        )

    alias = {
        "ligand_stem": task.lig_stem,
        "ph_label": task.ph_label,
        "ph_value": task.ph_value,
    }
    if alias not in entry.get("aliases", []):
        entry.setdefault("aliases", []).append(alias)
        ctx.microstates.dirty = True


def _ligprep_source_sidecar(pdbqt_path: Path) -> Path:
    return pdbqt_path.with_suffix(".ligprep_source.json")


def _copy_source_sidecar(src_pdbqt: Path, dst_pdbqt: Path) -> None:
    src = _ligprep_source_sidecar(src_pdbqt)
    if src.is_file():
        shutil.copy2(src, _ligprep_source_sidecar(dst_pdbqt))


def _write_source_sidecar(ctx: BulkContext, task: _MeekoPrepareTask) -> None:
    sidecar = _ligprep_source_sidecar(task.pdbqt_path)
    payload = {
        "schema_version": 1,
        "created_at": datetime.utcnow().isoformat() + "Z",
        "writer": "meeko",
        "pdbqt_path": str(task.pdbqt_path),
        "pdbqt_rel": task.pdbqt_rel,
        "ligand_stem": task.lig_stem,
        "ph_label": task.ph_label,
        "ph_value": task.ph_value,
        "source_format": "sdf",
        "source_sdf": str(task.source_sdf),
        "source_record_index": task.source_record_index,
        "source_record_name": task.source_record_name,
        "library_name": getattr(ctx, "library_name", ""),
        "chemistry_authoritative": True,
        "notes": "PDBQT was generated from this SDF record; downstream docked pose SDFs should reuse this chemistry and only replace coordinates.",
    }
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    tmp = sidecar.with_suffix(sidecar.suffix + ".part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, sidecar)


def bulk_prepare_sdfs_with_meeko(
    ctx: BulkContext,
    only_set: Optional[Set[str]],
) -> list[Path]:
    sdf_files = _discover_sdf_files(ctx)
    if not sdf_files:
        unit_sdf_dir = ctx.ligands_mol2_dir / "_rdkit_embedded_sdf"
        if unit_sdf_dir.is_dir():
            sdf_files = sorted(unit_sdf_dir.glob("*.sdf"))
            if sdf_files:
                print(f"[ligprep.meeko] using cached per-ligand SDFs from {unit_sdf_dir}")
    if not sdf_files:
        ctx.save_microstates(force=True)
        return []

    if ctx.ph_values is not None and len(ctx.eff_ph_values) > 1:
        logger.warning(
            "[ligprep.meeko] multiple pH labels requested, but Meeko backend does not fake pH microstates. "
            "Provide real protonation/tautomer SDF records from Molscrub/Scrubber for publication pH ensembles."
        )

    records = _iter_sdf_molecules(sdf_files, ctx, only_set)
    if only_set and not records:
        print("[test-mode] No requested ligands were found. Nothing to do; exiting cleanly.")
        ctx.save_microstates(force=True)
        return []

    tasks: list[tuple[_MeekoPrepareTask, Chem.Mol]] = []
    resume_skips = 0
    resume_examples: list[str] = []
    for mol, sdf_path, raw_name, lig_stem, record_index in records:
        for current_ph in ctx.eff_ph_values:
            pdbqt_path, pdbqt_rel, current_ph_label = _pdbqt_name_for_task(
                ctx, lig_stem, float(current_ph)
            )
            resume_skip = (
                not ctx.force
                and pdbqt_path.exists()
                and pdbqt_path.stat().st_size > 100
                and is_valid_ligand(pdbqt_path, log_dir=ctx.output_ligands_dir)
            )
            if resume_skip:
                resume_skips += 1
                if len(resume_examples) < 10:
                    resume_examples.append(pdbqt_rel)
                continue
            task = _MeekoPrepareTask(
                lig_stem=lig_stem,
                ph_label=current_ph_label,
                ph_value=float(current_ph),
                pdbqt_path=pdbqt_path,
                pdbqt_rel=pdbqt_rel,
                source_sdf=sdf_path,
                source_record_index=record_index,
                source_record_name=raw_name,
            )
            tasks.append((task, mol))

    workers = max(1, min(32, int(os.environ.get("CPU", "8")), max(1, len(tasks))))
    futures: Dict[concurrent.futures.Future, _MeekoPrepareTask] = {}
    result_paths: list[Path] = []
    print(
        f"Preparing {len(tasks)} SDF ligand record(s) with Meeko into {ctx.output_ligands_dir}"
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        for task, mol in tasks:
            futures[
                executor.submit(
                    prepare_mol_to_pdbqt_with_meeko,
                    Chem.Mol(mol),
                    task.pdbqt_path,
                    ligand_name=task.lig_stem,
                    log_dir=ctx.output_ligands_dir,
                )
            ] = task

        total = len(futures)
        for i, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            task = futures[future]
            result = future.result()
            result_paths.append(ctx.output_ligands_dir / task.pdbqt_rel)
            try:
                _append_prep_status(
                    ctx.status_log,
                    task.lig_stem,
                    "OK" if result.ok else "FAIL",
                    "" if result.ok else result.status,
                    task.pdbqt_rel,
                )
            except Exception:
                pass
            if result.ok:
                _write_source_sidecar(ctx, task)
                _update_microstate_registry(ctx, task)
            if i % 100 == 0 or not result.ok:
                print(f"[{i}/{total}] {task.lig_stem}: {result.status}")

    print(
        f"[ligprep.meeko] scheduled={len(tasks)} resume_skips={resume_skips} "
        f"force={ctx.force} examples_skipped={resume_examples[:5]}"
    )
    ctx.save_microstates(force=True)
    return result_paths


def prep_ligands_with_meeko(
    *,
    force: bool = False,
    only: Optional[Set[str]] = None,
    ph_values: Optional[list[float]] = None,
    microstate_dedup: bool = False,
    in_sdf: Optional[Path] = None,
    in_sdf_dir: Optional[Path] = None,
    in_pdb_dir: Optional[Path] = None,
    mol2_dir: Optional[Path] = None,
    out_pdbqt_dir: Optional[Path] = None,
    status_log: Optional[Path] = None,
    root_dir: Optional[Path] = None,
) -> object:
    ctx = init_bulk_context(
        BulkInitArgs(
            force=force,
            microstate_dedup=microstate_dedup,
            ph_values=ph_values,
            in_sdf=in_sdf,
            in_sdf_dir=in_sdf_dir,
            in_pdb_dir=in_pdb_dir,
            mol2_dir=mol2_dir,
            out_pdbqt_dir=out_pdbqt_dir,
            status_log=status_log,
            root_dir=root_dir,
        )
    )

    if ctx["in_pdb_dir_env"] and not ctx["in_sdf_env"]:
        from prep_ligands.prep_ligands_crystal import prep_ligands_from_pdb

        result = prep_ligands_from_pdb(
            Path(ctx["in_pdb_dir_env"]).resolve(),
            ctx["ligands_mol2_dir"],
            ctx["output_ligands_dir"],
        )
        ctx["maybe_save_microstate_registry"](force=True)
        return result
    return bulk_prepare_sdfs_with_meeko(ctx, only)


__all__ = ["bulk_prepare_sdfs_with_meeko", "prep_ligands_with_meeko"]
