"""Phase-3 bulk SDF to MOL2 orchestration helpers."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Set

from prep_ligands.prep_ligands_common import collapse_sanitized_path
from prep_ligands.prep_ligands_runtime import BulkContext

logger = logging.getLogger(__name__)


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
        logging.info(
            "[ligprep.skip] skipped %d SDF(s) with 'docked' in name under %s",
            skipped,
            ctx.ligand_extracted_dir,
        )
    print(f"Found {len(sdf_files)} SDF file(s)")
    return sdf_files


def _normalize_sdf_path(sdf_file: Path) -> Path:
    normalized = collapse_sanitized_path(sdf_file)
    if normalized == sdf_file:
        return sdf_file
    try:
        normalized.write_bytes(sdf_file.read_bytes())
        logging.info("[tidy] normalized double-sanitized SDF -> %s", normalized.name)
        return normalized
    except Exception as exc:
        logging.warning("[tidy] unable to normalize SDF %s: %s", sdf_file, exc)
        return sdf_file


def _recover_preexisting_mol2s(ctx: BulkContext) -> list[Path]:
    for root in (ctx.ligands_mol2_dir, ctx.output_ligands_dir):
        probe = sorted(root.glob("*.mol2"))
        if probe:
            logging.warning(
                "[compat] Found pre-existing MOL2s under %s (DEPRECATED layout); continuing.",
                root,
            )
            return probe
    return []


def _filter_requested_mol2s(
    mol2_files: Sequence[Path],
    only_set: Optional[Set[str]],
    *,
    maybe_save: Callable[..., None],
) -> list[Path]:
    if not only_set:
        return list(mol2_files)

    stems = {path.stem for path in mol2_files}
    requested = set(only_set)
    filtered = [path for path in mol2_files if path.stem in requested]
    sample = sorted(requested)[:10]
    logging.info(
        "[test-mode] enabled count=%d remain=%d sample=%s",
        len(requested),
        len(filtered),
        ",".join(sample),
    )
    print(
        f"[test-mode] enabled count={len(requested)} remain={len(filtered)} sample={','.join(sample)}"
    )
    missing = sorted(requested - stems)
    if missing:
        logging.warning(
            "[test-mode] requested ligands not found among MOL2s: %s",
            ",".join(missing[:20]) + ("..." if len(missing) > 20 else ""),
        )
    if not filtered:
        print(
            "[test-mode] No requested ligands were found. Nothing to do; exiting cleanly."
        )
        maybe_save(force=True)
    return filtered


def _convert_single_sdf(
    ctx: BulkContext,
    sdf_file: Path,
    *,
    only_set: Optional[Set[str]],
    rdkit_embed_fn: Callable[..., List[Path]],
    write_parent_sdf_fn: Callable[[Path, Path], int],
    convert_parallel_fn: Callable[..., List[Path]],
) -> list[Path]:
    print(f"=== Processing SDF: {sdf_file.name} ===")
    sdf_abs = _normalize_sdf_path(sdf_file).resolve()
    max_workers = max(1, int(os.environ.get("CPU", "8")))

    if ctx.use_rdkit_for_3d:
        print(
            "Using RDKit ETKDG for 3D with parent-picking; OBabel only for format conversion "
        )
        return rdkit_embed_fn(
            sdf_abs,
            ctx.ligands_mol2_dir,
            obabel_exe=ctx.obabel_exe_short,
            max_workers=max_workers,
            only_set=only_set,
        )

    print("Using OBabel --gen3d; pre-cleaning SDF to parent-only ")
    cleaned_sdf = ctx.ligands_mol2_dir / f"{sdf_abs.stem}_parents.sdf"
    n_kept = write_parent_sdf_fn(sdf_abs, cleaned_sdf)
    print(f"Parent-only SDF kept {n_kept} records")
    if n_kept == 0:
        print("No parent molecules survived desalting; skipping.")
        return []
    return convert_parallel_fn(
        cleaned_sdf,
        ctx.ligands_mol2_dir,
        ctx.obabel_exe_short,
        threads=ctx.obabel_threads,
        timeout_sec=ctx.obabel_timeout_s,
        chunk_size=ctx.chunk_size,
    )


def bulk_generate_mol2_files(
    ctx: BulkContext,
    only_set: Optional[Set[str]],
    *,
    rdkit_embed_fn: Callable[..., List[Path]],
    write_parent_sdf_fn: Callable[[Path, Path], int],
    convert_parallel_fn: Callable[..., List[Path]],
) -> List[Path]:
    sdf_files = _discover_sdf_files(ctx)
    if not sdf_files:
        ctx.save_microstates(force=True)
        return []

    all_mol2_files: list[Path] = []
    for sdf_file in sdf_files:
        mol2_files = _convert_single_sdf(
            ctx,
            sdf_file,
            only_set=only_set,
            rdkit_embed_fn=rdkit_embed_fn,
            write_parent_sdf_fn=write_parent_sdf_fn,
            convert_parallel_fn=convert_parallel_fn,
        )
        if not mol2_files:
            mol2_files = _recover_preexisting_mol2s(ctx)
        if not mol2_files:
            print("No MOL2 files produced; skipping this SDF.")
            continue

        filtered = _filter_requested_mol2s(
            mol2_files,
            only_set,
            maybe_save=ctx.save_microstates,
        )
        if only_set and not filtered:
            return []
        all_mol2_files.extend(filtered)
    return all_mol2_files


__all__ = ["bulk_generate_mol2_files"]
