from __future__ import annotations

import logging
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple


def _convert_one_pdbqt_to_mol2(src: Path, dst: Path, logger: logging.Logger) -> bool:
    """
    Convert a single PDBQT file to MOL2 using Open Babel.
    Returns True if conversion was performed, False if skipped (already exists).
    """
    if dst.exists() and dst.stat().st_size > 0:
        return False

    dst.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "obabel",
        "-ipdbqt",
        str(src),
        "-omol2",
        "-O",
        str(dst),
    ]
    logger.debug("[ledock.mol2.convert] cmd=%s", " ".join(cmd))
    subprocess.run(cmd, check=True)

    return True


def _convert_one_pdbqt_to_mol2_worker(arg: Tuple[Path, Path]) -> Tuple[Path, Path, bool, str | None]:
    src, dst = arg
    worker_logger = logging.getLogger("ledock.mol2")
    try:
        converted = _convert_one_pdbqt_to_mol2(src, dst, worker_logger)
        return src, dst, converted, None
    except Exception as e:
        return src, dst, False, str(e)


def mirror_library_to_mol2(
    library_root: Path,
    logger: logging.Logger,
    max_workers: int | None = None,
) -> None:
    """
    Mirror all .pdbqt ligands under `library_root` into a
    `library_root / f\"{library_root.name}_mol2\"` subtree, preserving
    relative directory structure, converting with Open Babel.
    """
    if not library_root.exists():
        logger.warning("[ledock.mol2] library_root_missing=%s", library_root)
        return

    library_name = library_root.name
    mol2_root = library_root / f"{library_name}_mol2"

    tasks: List[Tuple[Path, Path]] = []

    for src in library_root.rglob("*.pdbqt"):
        if mol2_root in src.parents:
            continue

        rel = src.relative_to(library_root)
        dst = mol2_root / rel.with_suffix(".mol2")
        tasks.append((src, dst))

    if not tasks:
        logger.info("[ledock.mol2] library=%s no_pdbqt_files_found", library_root)
        return

    logger.info(
        "[ledock.mol2] library=%s mol2_root=%s n_pdbqt=%d",
        library_root,
        mol2_root,
        len(tasks),
    )

    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        for src, dst, _converted, err in ex.map(_convert_one_pdbqt_to_mol2_worker, tasks):
            if err:
                logger.warning(
                    "[ledock.mol2.error] src=%s dst=%s reason=%s",
                    src,
                    dst,
                    err,
                )


def _group_library_roots_from_ligands(ligands: Iterable[Path | str]) -> List[Path]:
    """
    Given a set of ligand PDBQT paths, infer distinct library roots of the form
    .../prepped_ligands/<library>, based on the first path component after 'prepped_ligands'.
    """
    roots: Dict[Path, bool] = {}
    for p in ligands:
        path = Path(p)
        try:
            parts = path.resolve(strict=False).parts
        except Exception:
            parts = path.parts
        try:
            idx = parts.index("prepped_ligands")
        except ValueError:
            continue
        if idx + 1 >= len(parts):
            continue
        library_root = Path(*parts[: idx + 2])
        roots[library_root] = True
    return list(roots.keys())


def ensure_mol2_for_ledock(
    cfg: Dict[str, Any],
    ligands: Iterable[Path | str],
    logger: logging.Logger,
    max_workers: int | None = None,
) -> None:
    """
    For the given set of ligand PDBQT paths, ensure their libraries under
    `prepped_ligands/<library>` have mirrored MOL2 trees under
    `prepped_ligands/<library>/<library>_mol2`.
    """
    ligands = list(ligands)
    if not ligands:
        logger.info("[ledock.mol2] no_ligands_for_run; nothing to do")
        return

    library_roots = _group_library_roots_from_ligands(ligands)
    if not library_roots:
        logger.info("[ledock.mol2] no_library_roots_detected")
        return

    logger.info(
        "[ledock.mol2] preparing_mol2_for_libraries=%s",
        ", ".join(str(r) for r in library_roots),
    )

    for root in library_roots:
        mirror_library_to_mol2(root, logger=logger, max_workers=max_workers)
