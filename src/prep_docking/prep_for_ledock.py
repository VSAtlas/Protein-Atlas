from __future__ import annotations

import json
import logging
import os
import re
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from path_router import ph_ensemble_dir


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


def _convert_one_pdbqt_to_mol2_worker(
    arg: Tuple[Path, Path],
) -> Tuple[Path, Path, bool, str | None]:
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
        for src, dst, _converted, err in ex.map(
            _convert_one_pdbqt_to_mol2_worker, tasks
        ):
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


def map_pdbqt_to_mol2_path(
    pdbqt_path: Path | str,
    logger: Optional[logging.Logger] = None,
) -> Optional[Path]:
    """
    Map a PDBQT ligand path under prepped_ligands/<library>/... to the
    mirrored MOL2 path under prepped_ligands/<library>/<library>_mol2/... .
    """
    path = Path(pdbqt_path)
    if path.suffix.lower() == ".mol2":
        return path

    try:
        resolved = path.resolve(strict=False)
    except Exception:
        resolved = path

    parts = resolved.parts
    try:
        idx = parts.index("prepped_ligands")
    except ValueError:
        if logger:
            logger.debug("[ledock.mol2.map] outside_prepped_ligands=%s", path)
        return None

    if idx + 1 >= len(parts):
        return None

    library_root = Path(*parts[: idx + 2])
    library_name = library_root.name
    mol2_root = library_root / f"{library_name}_mol2"

    try:
        rel = resolved.relative_to(mol2_root)
        return (mol2_root / rel).with_suffix(".mol2")
    except Exception:
        pass

    try:
        rel = resolved.relative_to(library_root)
    except Exception:
        return None

    return (mol2_root / rel).with_suffix(".mol2")


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


def _normalize_ph_label(ph_label: Optional[str]) -> Optional[str]:
    token = str(ph_label).strip() if ph_label is not None else ""
    return token or None


def _canonicalize_ph_key(raw: str | None, pdb_id: str | None = None) -> Optional[str]:
    token = str(raw).strip() if raw is not None else ""
    if not token:
        return None

    # Normalize degenerate tokens to a stable base key for withH selection.
    looks_like_path = (
        any(sep in token for sep in ("/", "\\")) or "." in Path(token).name
    )
    if looks_like_path:
        token = Path(token).name
        token = Path(token).stem
        if token.endswith(".withH"):
            token = token[: -len(".withH")]

    if pdb_id:
        prefix = f"{str(pdb_id).upper()}_"
        if token.upper().startswith(prefix):
            token = token[len(prefix) :]

    if "+" in token:
        token = token.split("+", 1)[0]

    token = re.sub(r"-dup\d*$", "", token, flags=re.IGNORECASE)
    token = token.strip()
    return token or None


def _variant_for_ph(variant: Optional[str], legacy_mode: bool) -> Optional[str]:
    v = (str(variant).strip().upper() or None) if variant is not None else None
    if v:
        return v
    if legacy_mode:
        return None
    return "HOLO"


def _label_from_manifest_entry(entry: Dict[str, Any], pdb_id: str) -> Optional[str]:
    if not isinstance(entry, dict):
        return None
    for key in ("label", "ph_label"):
        if entry.get(key):
            return str(entry[key])
    tag = entry.get("tag")
    if tag:
        tag_str = str(tag)
        prefix = f"{pdb_id}_"
        return tag_str[len(prefix) :] if tag_str.startswith(prefix) else tag_str
    for key in ("withH", "pdbqt", "receptor_pdbqt", "output_pdbqt", "path", "receptor"):
        raw = entry.get(key)
        if not raw:
            continue
        stem = Path(str(raw)).stem
        if stem.endswith(".withH"):
            stem = stem[: -len(".withH")]
        prefix = f"{pdb_id}_"
        return stem[len(prefix) :] if stem.startswith(prefix) else stem
    return None


def get_lepro_exe(cfg: Dict[str, Any]) -> str:
    """
    Resolve the LePro executable path, honoring cfg and file_cfg overrides.
    """
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("LEPRO_EXE") or None
        if not candidate:
            file_cfg = cfg.get("_FILE_CFG")
            if isinstance(file_cfg, dict):
                candidate = file_cfg.get("LEPRO_EXE") or None
    if not candidate:
        candidate = shutil.which("lepro") or "lepro"
    logger = logging.getLogger("ledock.lepro")
    logger.info("[ledock.lepro.exe] path=%s", candidate)
    return str(candidate)


def _resolve_withh_from_manifest(
    ensemble_dir: Path,
    pdb_id: str,
    ph_label: str,
    logger: logging.Logger,
) -> Optional[Path]:
    manifest_path = ensemble_dir / "ensemble.json"
    if not manifest_path.exists():
        return None
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8", errors="ignore"))
    except Exception as exc:
        logger.warning("[ledock.manifest.error] path=%s reason=%s", manifest_path, exc)
        return None
    members = payload.get("members") if isinstance(payload, dict) else None
    if not isinstance(members, list):
        return None
    want = _canonicalize_ph_key(ph_label, pdb_id)
    if not want:
        logger.info(
            "[ledock.manifest.miss] pdb=%s raw=%s want=%s path=%s",
            pdb_id,
            ph_label,
            want,
            manifest_path,
        )
        return None
    for entry in members:
        entry_label_raw = _label_from_manifest_entry(entry, str(pdb_id).upper())
        if not entry_label_raw:
            continue
        entry_key = _canonicalize_ph_key(entry_label_raw, pdb_id)
        if not entry_key:
            continue
        if entry_key.strip().lower() != want.lower():
            continue
        withh = entry.get("withH") or entry.get("withH_pdb") or entry.get("withH_path")
        if not withh:
            continue
        path = Path(str(withh))
        if not path.is_absolute():
            path = manifest_path.parent / path
        return path
    logger.info(
        "[ledock.manifest.miss] pdb=%s raw=%s want=%s path=%s",
        pdb_id,
        ph_label,
        want,
        manifest_path,
    )
    return None


def ensure_ledock_receptor(
    cfg: Dict[str, Any],
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    logger: logging.Logger,
) -> Optional[Path]:
    """
    Ensure a LePro-processed receptor (pro.pdb) exists for the given target.
    """
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )
    variant_for_ph = _variant_for_ph(variant_token, legacy_mode)
    ph_token_raw = _normalize_ph_label(ph_label)
    if not ph_token_raw:
        logger.warning("[ledock.receptor.skip] reason=missing_ph_label pdb=%s", pdb_id)
        return None
    ph_key = _canonicalize_ph_key(ph_token_raw, pdb_id)
    if not ph_key:
        logger.warning(
            "[ledock.receptor.skip] reason=invalid_ph_label pdb=%s raw=%s",
            pdb_id,
            ph_token_raw,
        )
        return None

    ensemble_dir = ph_ensemble_dir(pdb_id, variant=variant_for_ph, legacy=legacy_mode)
    pro_path = ensemble_dir / "pro.pdb"
    if pro_path.exists() and pro_path.stat().st_size > 0:
        logger.info(
            "[ledock.receptor] pdb=%s variant=%s ph=%s receptor_pdb=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token_raw,
            pro_path,
        )
        return pro_path

    withh_path = _resolve_withh_from_manifest(ensemble_dir, pdb_id, ph_key, logger)
    if withh_path is None:
        prefix = f"{str(pdb_id).upper()}_"
        fallback_label = ph_key
        if not fallback_label.startswith(prefix):
            fallback_label = f"{prefix}{fallback_label}"
        candidate = ensemble_dir / f"{fallback_label}.withH.pdb"
        if not candidate.exists():
            logger.warning(
                "[ledock.receptor.missing] pdb=%s variant=%s ph=%s path=%s",
                pdb_id,
                variant_for_ph or "HOLO",
                ph_token_raw,
                candidate,
            )
            return None
        withh_path = candidate

    lepro_exe = get_lepro_exe(cfg)
    try:
        logger.info(
            "[ledock.lepro] exe=%s input=%s cwd=%s",
            lepro_exe,
            withh_path.name,
            ensemble_dir,
        )
        subprocess.run([lepro_exe, withh_path.name], check=True, cwd=ensemble_dir)
    except Exception as exc:
        logger.warning(
            "[ledock.lepro.error] pdb=%s variant=%s ph=%s reason=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token_raw,
            exc,
        )
        return None

    if pro_path.exists() and pro_path.stat().st_size > 0:
        logger.info(
            "[ledock.receptor] pdb=%s variant=%s ph=%s receptor_pdb=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token_raw,
            pro_path,
        )
        return pro_path

    logger.warning(
        "[ledock.receptor.missing] pdb=%s variant=%s ph=%s path=%s",
        pdb_id,
        variant_for_ph or "HOLO",
        ph_token_raw,
        pro_path,
    )
    return None


def symlink_mol2_for_stage(
    ligand_pairs: Iterable[Tuple[Path, Path]],
    stage_root: Path,
    logger: logging.Logger,
) -> Dict[Path, Path]:
    """
    Create symlinks for MOL2 ligands under stage_root using the original ligand
    stem to keep downstream .dok names aligned, and return a mapping from
    original PDBQT path to alias MOL2 path.
    """
    alias_map: Dict[Path, Path] = {}
    for idx, (lig_path, mol2_path) in enumerate(ligand_pairs, start=1):
        alias_stem = Path(lig_path).stem
        alias_ext = Path(mol2_path).suffix
        alias_name = f"{alias_stem}{alias_ext}"
        alias_path = stage_root / alias_name
        try:
            if alias_path.exists() or alias_path.is_symlink():
                alias_path.unlink()
        except FileNotFoundError:
            pass
        os.symlink(mol2_path, alias_path)
        logger.info("[ledock.symlink] src=%s alias=%s", mol2_path, alias_path)
        alias_map[Path(lig_path)] = alias_path

    return alias_map
