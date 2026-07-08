"""Phase-4 bulk ligand-prep scheduling and result handling."""

from __future__ import annotations

import concurrent.futures
import logging
import os
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

from prep_ligands.prep_ligands_bulk_sdf import rename_mol2_with_prefix
from prep_ligands.prep_ligands_common import _prepare_one, collapse_sanitized_path, is_valid_ligand
from prep_ligands.prep_ligands_microstates import compute_microstate_id
from prep_ligands.prep_ligands_reporting import _append_prep_status
from prep_ligands.prep_ligands_runtime import (
    BulkContext,
    ph_label,
    relative_to_output,
    sanitize_ligand_name_for_filename,
)


logger = logging.getLogger(__name__)


@dataclass(slots=True)
class _BulkPrepareTask:
    lig_stem: str
    ph_label: str
    ph_value: float
    pdbqt_path: Path
    pdbqt_rel: str
    mol2_input: Path


def _normalize_bulk_mol2(mol2_file: Path) -> Path:
    normalized = collapse_sanitized_path(mol2_file)
    if normalized == mol2_file:
        return mol2_file
    try:
        normalized.write_bytes(mol2_file.read_bytes())
        logging.info("[tidy] normalized double-sanitized -> %s", normalized.name)
        return normalized
    except Exception as exc:
        logging.warning("[tidy] unable to normalize %s: %s", mol2_file, exc)
        return mol2_file


def _resolve_ligand_stem(ctx: BulkContext, mol2_file: Path, seq: int) -> tuple[Path, str]:
    mol2_input = _normalize_bulk_mol2(mol2_file)
    lig_stem = mol2_input.stem
    name_path = mol2_input.with_suffix(".name")
    parent_name_raw = ""
    try:
        if name_path.exists():
            parent_name_raw = name_path.read_text(encoding="utf-8", errors="ignore").strip()
    except Exception as exc:
        logging.warning("[ligprep] unable to read .name for %s: %s", mol2_input.name, exc)

    if ctx.rename_active:
        lig_stem = f"{ctx.rename_prefix}{seq:0{ctx.rename_pad}d}"
        if ctx.rename_force or mol2_input.stem != lig_stem:
            try:
                mol2_input = rename_mol2_with_prefix(mol2_input, lig_stem)
            except Exception as exc:
                logging.warning("[ligprep] rename failed for %s: %s", mol2_input, exc)
    elif not ctx.is_fda_library and parent_name_raw:
        lig_stem = sanitize_ligand_name_for_filename(parent_name_raw)
    return mol2_input, lig_stem


def _iter_prepare_tasks(ctx: BulkContext, mol2_files: Iterable[Path]) -> tuple[list[_BulkPrepareTask], int, list[str]]:
    resume_skips = 0
    resume_examples: list[str] = []
    ligprep_debug_seen = 0
    ordered = sorted(mol2_files, key=lambda path: path.name) if ctx.rename_active else list(mol2_files)
    tasks: list[_BulkPrepareTask] = []

    for seq, mol2_file in enumerate(ordered, start=ctx.rename_start):
        mol2_input, lig_stem = _resolve_ligand_stem(ctx, Path(mol2_file), seq)
        name_path = mol2_input.with_suffix(".name")
        parent_name_raw = ""
        try:
            if name_path.exists():
                parent_name_raw = name_path.read_text(encoding="utf-8", errors="ignore").strip()
        except Exception:
            pass
        if ligprep_debug_seen < 5:
            print(
                f"[ligprep-debug] mol2={mol2_input.name} rename_active={ctx.rename_active} "
                f"name_path_exists={name_path.exists()} parent_name_raw={parent_name_raw!r} "
                f"lig_stem={lig_stem}"
            )
            ligprep_debug_seen += 1

        for current_ph in ctx.eff_ph_values:
            current_ph_label = ph_label(current_ph) if ctx.use_ph_subdirs else None
            out_dir = (
                ctx.output_ligands_dir / current_ph_label
                if current_ph_label
                else ctx.output_ligands_dir
            )
            out_dir.mkdir(parents=True, exist_ok=True)

            if ctx.use_ph_subdirs and current_ph_label:
                base, _, rest = lig_stem.partition("_")
                lig_stem_ph = f"{base}{current_ph_label}_{rest}" if rest else f"{lig_stem}{current_ph_label}"
                pdbqt_name = f"{lig_stem_ph}.pdbqt"
            else:
                pdbqt_name = f"{lig_stem}.pdbqt"

            pdbqt_path = out_dir / pdbqt_name
            pdbqt_rel = relative_to_output(pdbqt_path, ctx.output_ligands_dir)

            resume_skip = (
                not ctx.force
                and pdbqt_path.exists()
                and pdbqt_path.stat().st_size > 100
                and is_valid_ligand(pdbqt_path, log_dir=ctx.output_ligands_dir)
            )
            if resume_skip:
                logging.info(
                    "[resume] Valid PDBQT exists, skipping: %s ph=%.2f",
                    pdbqt_path,
                    current_ph,
                )
                resume_skips += 1
                if len(resume_examples) < 10:
                    resume_examples.append(pdbqt_rel)
                continue

            if ctx.force and pdbqt_path.exists():
                logging.info("[force] Overwriting existing PDBQT: %s", pdbqt_rel)

            tasks.append(
                _BulkPrepareTask(
                    lig_stem=lig_stem,
                    ph_label=current_ph_label or "",
                    ph_value=float(current_ph),
                    pdbqt_path=pdbqt_path,
                    pdbqt_rel=pdbqt_rel,
                    mol2_input=mol2_input,
                )
            )
    return tasks, resume_skips, resume_examples


def _update_microstate_registry(ctx: BulkContext, task: _BulkPrepareTask) -> None:
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
        logger.debug(
            "microstate_dedup: add_alias stem=%s ph_label=%s ph_value=%.2f microstate_id=%s",
            task.lig_stem,
            task.ph_label,
            task.ph_value,
            entry["microstate_id"],
        )
    if ctx.microstates.alias_index is not None:
        ctx.microstates.alias_index[(task.lig_stem, task.ph_label, task.ph_value)] = entry


def bulk_prepare_pdbqts(ctx: BulkContext, mol2_files: List[Path]) -> List[Path]:
    if not mol2_files:
        return []

    print(
        f"Preparing {len(mol2_files)} MOL2 files with MGLTools (parallel) into {ctx.output_ligands_dir}"
    )
    tasks, resume_skips, resume_examples = _iter_prepare_tasks(ctx, mol2_files)
    futures: Dict[concurrent.futures.Future[tuple[str, str]], _BulkPrepareTask] = {}

    with ThreadPoolExecutor(max_workers=max(1, int(os.environ.get("CPU", "8")))) as executor:
        for task in tasks:
            logger.debug(
                "prep_ligands: scheduling ligand prep stem=%s ph=%.2f",
                task.lig_stem,
                task.ph_value,
            )
            futures[
                executor.submit(
                    _prepare_one,
                    ctx.mgltools_python_short,
                    ctx.prepare_script_short,
                    task.mol2_input,
                    task.pdbqt_path,
                    obabel_exe_short=ctx.obabel_exe_short,
                    status_log_dir=ctx.output_ligands_dir,
                    ph=task.ph_value,
                )
            ] = task

        total = len(futures)
        result_paths: List[Path] = []
        for i, future in enumerate(as_completed(futures), start=1):
            task = futures[future]
            name, status = future.result()
            lig_name = Path(name).stem
            result_paths.append(ctx.output_ligands_dir / task.pdbqt_rel)
            try:
                _append_prep_status(
                    ctx.status_log,
                    lig_name,
                    "OK" if status == "ok" else "FAIL",
                    "" if status == "ok" else status,
                    task.pdbqt_rel,
                )
            except Exception:
                pass

            if status == "ok":
                _update_microstate_registry(ctx, task)
            if i % 100 == 0 or status != "ok":
                print(f"[{i}/{total}] {name}: {status}")

    print(
        f"[ligprep] scheduled={len(tasks)} resume_skips={resume_skips} force={ctx.force} "
        f"examples_skipped={resume_examples[:5]}"
    )
    return result_paths


__all__ = ["bulk_prepare_pdbqts"]
