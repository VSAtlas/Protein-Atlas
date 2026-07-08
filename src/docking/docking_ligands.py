from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.MolStandardize import rdMolStandardize

from docking.docking_ligand_selection import (
    _is_under_ph_subdir,
    _iter_pdbqt_dirfirst,
    select_ligands_for_next,
)
from docking.active_preserving_cap import is_dud_active_key, stable_active_preserving_limit
from docking.library_mode import (
    _coerce_test_map,
    _parse_test_libraries_value,
    compute_allowed_library_roots,
    parse_test_libraries,
)
from docking.pose_validation import compute_redock_rmsd
from docking.ligand_metrics import _is_readable_ref, _read_any_lig, compute_rmsd
from prep_ligands.prep_ligands_crystal import prep_ligands_from_pdb
from path_router.path_router import Paths
from prep_ligands.library_index import LibraryIndex

def norm(p: str | Path) -> str:
    """Normalize path to a clean, forward-slash string for logs & keys."""
    return os.path.abspath(str(p)).replace("\\", "/")


def _chunk_ligand_key(name: str | Path) -> str:
    base = Path(str(name or "")).name.strip().lower()
    for suffix in (".pdbqt", ".mol2", ".sdf", ".pdb"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    if base.endswith(".sanitized"):
        base = base[: -len(".sanitized")]
    return base.strip()


def _count_heavy_atoms_from_pdbqt(pdbqt_path: Path) -> int:
    heavy = 0
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            element = line[76:78].strip() if len(line) >= 78 else ""
            if element:
                if element.upper() != "H":
                    heavy += 1
            else:
                atom_name = line[12:16].strip()
                if not atom_name.upper().startswith("H"):
                    heavy += 1
    return heavy


def _load_mol_any(pdbqt_path: Path, obabel_exe: str | None) -> Chem.Mol | None:
    base = pdbqt_path.with_suffix("")
    # prefer SDF, then MOL2, then PDB
    sdf = base.with_suffix(".sdf")
    mol2 = base.with_suffix(".mol2")
    pdb = base.with_suffix(".pdb")
    if sdf.exists():
        supp = Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True)
        for m in supp:
            if m:
                return m
    if mol2.exists():
        m = Chem.MolFromMol2File(str(mol2), sanitize=True, removeHs=False)
        if m:
            return m
    if pdb.exists():
        m = Chem.MolFromPDBFile(str(pdb), sanitize=True, removeHs=False)
        if m:
            return m
    # fallback: PDBQT -> SDF via obabel (Linux-friendly)
    if obabel_exe:
        with tempfile.TemporaryDirectory() as td:
            out = Path(td, "tmp.sdf")
            try:
                subprocess.check_call(
                    [
                        obabel_exe,
                        "-ipdbqt",
                        str(pdbqt_path),
                        "-osdf",
                        "-O",
                        str(out),
                        "--retype",
                        "--addh",
                    ]
                )
                supp = Chem.SDMolSupplier(str(out), removeHs=False, sanitize=True)
                for m in supp:
                    if m:
                        return m
            except Exception:
                return None
    return None


def _standardize(m: Chem.Mol) -> Chem.Mol:
    parent = rdMolStandardize.ChargeParent(m)  # neutralize/parent
    rdMolStandardize.Normalize(parent)  # FG normalization
    Chem.SanitizeMol(parent)
    return parent


# Build catalog with PAINS A/B/C
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
pains_catalog = FilterCatalog.FilterCatalog(params)


def _resolve_test_mode(cfg) -> str:
    """
    Compatibility helper: return legacy mode strings for reserved-only configs,
    or a "+"-joined token string when custom library tokens are used.
    """
    raw = (
        os.environ.get("TEST_MODE_ENABLE")
        if "TEST_MODE_ENABLE" in os.environ
        else cfg.get("TEST_MODE_ENABLE", "off")
    )
    if isinstance(raw, bool):
        if not raw:
            return "off"
    s = str(raw).strip()
    if not s or s.lower() in {"", "0", "false", "no", "off", "none", "null"}:
        return "off"

    tokens = parse_test_libraries(cfg)
    reserved = {"dud", "fda", "hmdb"}
    if tokens and all(tok in reserved for tok in tokens):
        token_set = set(tokens)
        if token_set == {"dud"}:
            return "dud"
        if token_set == {"fda"}:
            return "fda"
        if token_set == {"hmdb"}:
            return "hmdb"
        if token_set == {"dud", "fda"}:
            return "fda+dud"
        if token_set == {"hmdb", "dud"}:
            return "hmdb+dud"
        if token_set == {"hmdb", "fda"}:
            return "hmdb+fda"
        if token_set == {"dud", "fda", "hmdb"}:
            return "fda+dud+hmdb"
    return "+".join(tokens) if tokens else "off"


def _dedup_index_roots(seq: list[Path]) -> list[Path]:
    deduped: list[Path] = []
    seen: set[str] = set()
    for candidate in seq:
        if not candidate:
            continue
        path_obj = Path(candidate)
        try:
            key = str(path_obj.resolve())
        except Exception:
            key = str(path_obj)
        if key not in seen:
            seen.add(key)
            deduped.append(path_obj)
    return deduped


def _lib_roots_for_pdb(
    cfg: Dict,
    pdb_id: str,
    paths: Paths,
    logger: logging.Logger,
    *,
    test_mode_override: Optional[str] = None,
) -> tuple[list[Path], list[Path]]:
    del paths
    tokens_override = (
        _parse_test_libraries_value(test_mode_override)
        if test_mode_override is not None
        else None
    )
    allowed_noncontrol_roots = compute_allowed_library_roots(
        cfg,
        pdb_id,
        logger,
        tokens_override=tokens_override,
    )
    return [], allowed_noncontrol_roots


def prepare_and_filter_ligands(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    *,
    run_mode: Optional[str] = None,
) -> Tuple[List[str], Dict[str, int], Dict[str, bool]]:
    """
    Gathers candidate ligands, keeps existing validation/PAINS logic, and
    filters the *non-control* pool to allowed library roots derived from
    TEST_MODE_ENABLE tokens.

    run_mode:
      - None: use parsed TEST_MODE_ENABLE tokens
      - "<token>": treat run_mode as a single token subrun (reserved or custom)

    Controls are *never* filtered out here.
    LIBRARY_EXTRA_DIRS remain included (unchanged).
    """
    overall_tokens = parse_test_libraries(cfg)
    effective_tokens = overall_tokens
    if run_mode is not None:
        token = str(run_mode).strip()
        if token:
            lowered = token.lower()
            if lowered in {"default", "off", "none", "null"}:
                lowered = "fda"
            effective_tokens = [lowered]

    chunk_keys_raw = cfg.get("_CHUNK_LIGAND_KEYS")
    if isinstance(chunk_keys_raw, (list, tuple, set)):
        chunk_keys = {_chunk_ligand_key(str(x)) for x in chunk_keys_raw if str(x).strip()}
    else:
        chunk_keys = set()

    def _record_soft_failure(reason: str) -> None:
        token = str(reason or "unknown_soft_failure").strip() or "unknown_soft_failure"
        try:
            current = int(cfg.get("_SOFTFAIL_LIGAND_COUNT", 0) or 0)
        except Exception:
            current = 0
        cfg["_SOFTFAIL_LIGAND_COUNT"] = int(current + 1)
        reasons_raw = cfg.get("_SOFTFAIL_LIGAND_REASONS")
        if isinstance(reasons_raw, list):
            reasons = [str(x) for x in reasons_raw if str(x).strip()]
        else:
            reasons = []
        if token not in reasons:
            reasons.append(token)
        cfg["_SOFTFAIL_LIGAND_REASONS"] = reasons

    # Chunk workers dock only manifest-selected library ligands. Extracted-control
    # prep/redock is owned by the combo prep task, so repeating it per chunk turns
    # startup into a filesystem race.
    if chunk_keys:
        logger.info(
            "[ligprep.extracted] skip reason=distributed_chunk chunk_id=%s selected=%d",
            str(cfg.get("_CHUNK_ID", "")),
            len(chunk_keys),
        )
    else:
        try:
            prep_ligands_from_pdb(
                ligand_output_dir=paths.ligand_output_dir,
                ligands_mol2_dir=paths.ligands_mol2_dir,
                prepped_ligands_dir=paths.prepped_ligands_dir,
            )
        except Exception as exc:
            _record_soft_failure("extracted_ligprep_exception")
            logger.warning(
                "[ligprep.extracted.softfail] pdb=%s run_mode=%s error=%s",
                paths.pdb_id,
                str(run_mode or "default"),
                exc,
                exc_info=True,
            )

    chunk_library_root = str(cfg.get("_CHUNK_LIBRARY_ROOT") or "").strip()
    if chunk_library_root:
        root_path = Path(chunk_library_root)
        allowed_noncontrol_roots = [root_path] if root_path.exists() else []
        logger.info(
            "[ligands.chunk-root] chunk_id=%s run_mode=%s root=%s exists=%s",
            str(cfg.get("_CHUNK_ID", "")),
            str(run_mode or ""),
            chunk_library_root,
            int(root_path.exists()),
        )
    else:
        allowed_noncontrol_roots = compute_allowed_library_roots(
            cfg,
            paths.pdb_id.upper(),
            logger,
            tokens_override=effective_tokens,
        )
    per_index_roots: list[Path] = []
    if paths.prepped_ligands_dir:
        per_index_roots.append(paths.prepped_ligands_dir)

    index_library_roots: list[Path] = _dedup_index_roots(allowed_noncontrol_roots)
    per_index_roots = _dedup_index_roots(per_index_roots)
    index_roots = _dedup_index_roots(per_index_roots + index_library_roots)

    manifest_filename = str(cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json"))
    lib_index = cfg.get("_LIB_INDEX")
    if not isinstance(lib_index, LibraryIndex):
        lib_index = LibraryIndex(manifest_filename=manifest_filename, logger=logger)
        cfg["_LIB_INDEX"] = lib_index
    if index_roots:
        lib_index.load(index_roots)
    cfg["_LIB_INDEX_PER_ROOTS"] = [str(p) for p in per_index_roots]
    cfg["_LIB_INDEX_LIBRARY_ROOTS"] = [str(p) for p in index_library_roots]

    def _finalize_selected_paths(final_paths: list[Path], *, controls_count: int) -> Tuple[List[str], Dict[str, int], Dict[str, bool]]:
        cfg["ALL_LIGAND_PATHS"] = [str(p) for p in final_paths]
        cfg["ALL_LIGAND_PATHS_VALID"] = [str(p) for p in final_paths]
        if chunk_keys:
            cfg["_CHUNK_EFFECTIVE_LIGAND_KEYS"] = sorted(
                {
                    key
                    for key in (_chunk_ligand_key(p.name) for p in final_paths)
                    if key
                }
            )
        else:
            cfg.pop("_CHUNK_EFFECTIVE_LIGAND_KEYS", None)

        pains_flags: Dict[str, bool] = {}
        heavy_atom_counts: Dict[str, int] = {}
        for p in final_paths:
            pains_flags[str(p)] = False
            pains_flags[p.stem] = False
            try:
                ha = _count_heavy_atoms_from_pdbqt(p)
            except Exception:
                ha = 0
            heavy_atom_counts[str(p)] = ha
            heavy_atom_counts[p.stem] = ha

        ligands = [str(p) for p in final_paths]
        non_control_count = max(0, len(final_paths) - int(controls_count))
        logger.info(
            f"Selected ligands -> controls={int(controls_count)} + non-controls={non_control_count} = total={len(ligands)}"
        )
        bad = [p for p in ligands if not str(p).lower().endswith(".pdbqt")]
        if bad:
            raise ValueError(f"Ligand is not a .pdbqt file: {bad[0]}")
        return ligands, heavy_atom_counts, pains_flags

    if chunk_keys:
        resolved_chunk_paths: list[Path] = []
        missing_chunk_keys: list[str] = []
        lookup_roots = list(index_library_roots)
        for key in sorted(chunk_keys):
            found = lib_index.lookup(key, lookup_roots, allow_prefix=False)
            if found is None:
                filename = key if key.lower().endswith(".pdbqt") else f"{key}.pdbqt"
                found = lib_index.lookup_filename(filename, lookup_roots)
            if found is None:
                missing_chunk_keys.append(key)
                continue
            try:
                if not found.exists() or found.stat().st_size <= 100:
                    missing_chunk_keys.append(key)
                    continue
            except Exception:
                missing_chunk_keys.append(key)
                continue
            resolved_chunk_paths.append(Path(found))

        if resolved_chunk_paths and not missing_chunk_keys:
            logger.info(
                "[ligands.chunk-direct] chunk_id=%s selected=%d roots=%d",
                str(cfg.get("_CHUNK_ID", "")),
                len(resolved_chunk_paths),
                len(lookup_roots),
            )
            return _finalize_selected_paths(resolved_chunk_paths, controls_count=0)
        logger.warning(
            "[ligands.chunk-direct] action=fallback_scan chunk_id=%s resolved=%d missing=%d sample_missing=%s",
            str(cfg.get("_CHUNK_ID", "")),
            len(resolved_chunk_paths),
            len(missing_chunk_keys),
            ",".join(missing_chunk_keys[:5]) or "none",
        )

    def _under(p: Path, root: Path) -> bool:
        try:
            p.resolve().relative_to(root.resolve())
            return True
        except Exception:
            pass
        try:
            p.absolute().relative_to(root.absolute())
            return True
        except Exception:
            return False

    def _enumerate_noncontrol_candidates_via_index(
        cfg: Dict,
        allowed_roots: list[Path],
        logger: logging.Logger,
    ) -> list[Path]:
        from prep_ligands.prep_ligands_microstates import enumerate_ligands_for_docking

        resolved_roots = _dedup_index_roots([Path(r) for r in allowed_roots if r])
        resolved_existing = [r for r in resolved_roots if r.exists()]
        logger.info(
            "[lib-roots] non-control roots = %s",
            [str(p) for p in resolved_existing],
        )

        candidates: list[Path] = []
        microstate_roots = [
            r for r in resolved_existing if (r / "microstates.json").exists()
        ]
        if microstate_roots:
            try:
                all_ms: list[Path] = []
                for root in microstate_roots:
                    ms_paths = enumerate_ligands_for_docking(
                        requested_ph_values=None,
                        cfg=cfg,
                        pdb_id=paths.pdb_id.upper(),
                        root_dir=root,
                        microstate_dedup=True,
                        force=False,
                    )
                    all_ms.extend(ms_paths)
                ms_filtered = [
                    p
                    for p in all_ms
                    if any(_under(p, root) for root in resolved_existing)
                ]
                candidates.extend(ms_filtered)
                logger.info(
                    "[lib-index.microstate] roots=%d ligands=%d",
                    len(microstate_roots),
                    len(ms_filtered),
                )
            except Exception as exc:
                logger.warning("[lib-index.microstate] error=%s", exc)

        manifest_candidates: list[Path] = []
        manifest_count_by_root: dict[Path, int] = {}
        if isinstance(cfg.get("_LIB_INDEX"), LibraryIndex):
            index_obj: LibraryIndex = cfg["_LIB_INDEX"]
            for root in resolved_existing:
                manifest = getattr(index_obj, "_cache", {}).get(Path(root))
                if not manifest:
                    continue
                rel_values: list[str] = []
                # filenames is the canonical full inventory; entries is lookup-oriented.
                if getattr(manifest, "filenames", None):
                    rel_values.extend(str(v) for v in manifest.filenames.values())
                elif getattr(manifest, "entries", None):
                    rel_values.extend(str(v) for v in manifest.entries.values())
                if rel_values:
                    manifest_count_by_root[Path(root)] = len(rel_values)
                    for rel in rel_values:
                        manifest_candidates.append(Path(root) / rel)
        if manifest_candidates:
            candidates.extend(manifest_candidates)
            logger.info(
                "[lib-index.manifest] roots=%d ligands=%d",
                len(resolved_existing),
                len(manifest_candidates),
            )

        # Guard against stale/partial manifests: compare against a lightweight
        # filesystem enumeration and augment when the manifest undercounts.
        fs_fallback_added = 0
        for root in resolved_existing:
            fs_paths = [Path(p) for p in _iter_pdbqt_dirfirst(root)]
            if not fs_paths:
                continue
            fs_count = len(fs_paths)
            manifest_count = int(manifest_count_by_root.get(Path(root), 0))
            if manifest_count >= fs_count:
                continue
            logger.warning(
                "[lib-index.manifest.undercount] root=%s manifest=%d filesystem=%d action=augment_with_filesystem",
                str(root),
                manifest_count,
                fs_count,
            )
            candidates.extend(fs_paths)
            fs_fallback_added += fs_count
        if fs_fallback_added:
            logger.info(
                "[lib-index.filesystem] roots=%d added=%d",
                len(resolved_existing),
                fs_fallback_added,
            )

        if not candidates:
            logger.warning(
                "[lib-index] no usable index detected; falling back to filesystem scan under %d roots (this is slow)",
                len(resolved_existing),
            )
            seen: set[str] = set()
            for root in resolved_existing:
                for p in _iter_pdbqt_dirfirst(root):
                    pn = norm(p)
                    if pn not in seen:
                        seen.add(pn)
                        candidates.append(p)

        deduped: list[Path] = []
        seen_keys: set[str] = set()
        for p in candidates:
            key = norm(p)
            if key not in seen_keys:
                seen_keys.add(key)
                deduped.append(Path(p))
        deduped.sort(key=lambda p: norm(p))
        return deduped

    scan_roots: list[Path] = []
    prepped_lig_root = paths.prepped_ligands_dir
    if prepped_lig_root and prepped_lig_root.exists():
        scan_roots.append(prepped_lig_root)

    scan_roots.extend(allowed_noncontrol_roots)
    scan_roots = _dedup_index_roots(scan_roots)

    single_selector = cfg.get("_EFFECTIVE_SINGLE_LIGAND")
    crawl_allowed = not bool(single_selector)
    logger.info(
        "[ligands.scan.guard] roots=%d crawl_allowed=%s selector=%s",
        len(scan_roots),
        crawl_allowed,
        single_selector,
    )

    seen: set[str] = set()
    controls: list[Path] = []
    per_protein_allow = {"controls", "reference"}
    if prepped_lig_root and prepped_lig_root.exists():
        logger.info(
            "[ligands.scan.root] root=%s allowed_subdirs=%s crawl=%s",
            prepped_lig_root,
            per_protein_allow,
            crawl_allowed,
        )
        for p in _iter_pdbqt_dirfirst(
            prepped_lig_root, allowed_subdirs=per_protein_allow
        ):
            pn = norm(p)
            if pn not in seen:
                seen.add(pn)
                controls.append(p)

    noncontrol_candidates = _enumerate_noncontrol_candidates_via_index(
        cfg, allowed_noncontrol_roots, logger
    )

    all_pdbqt_paths: list[Path] = []
    all_pdbqt_paths.extend(controls)
    all_pdbqt_paths.extend(noncontrol_candidates)

    logger.info(
        "[ligands.scan] roots=%d found=%d", len(scan_roots), len(all_pdbqt_paths)
    )
    cfg["ALL_LIGAND_PATHS"] = [str(p) for p in all_pdbqt_paths]

    if not all_pdbqt_paths:
        logger.warning("No .pdbqt ligands were found under the configured roots.")
        return [], {}, {}

    valid_pdbqt: Dict[str, Path] = {}
    for p in all_pdbqt_paths:
        try:
            if p.exists() and p.stat().st_size > 100:
                valid_pdbqt[norm(p)] = p
            else:
                logger.debug(f"Excluded malformed ligand (missing or too small): {p}")
        except Exception:
            logger.debug(f"Excluded malformed ligand (exception): {p}")

    cfg["ALL_LIGAND_PATHS_VALID"] = [str(p) for p in valid_pdbqt.values()]
    logger.info("[ligands.valid] count=%d", len(valid_pdbqt))

    controls_valid: list[Path] = []
    noncontrols: list[Path] = []
    for p in valid_pdbqt.values():
        if prepped_lig_root and _under(p, prepped_lig_root):
            controls_valid.append(p)
        else:
            noncontrols.append(p)

    filtered_noncontrols: list[Path] = []
    for p in noncontrols:
        keep = False
        for root in allowed_noncontrol_roots:
            if root.exists() and _under(p, root):
                keep = True
                break
        if keep:
            filtered_noncontrols.append(p)

    logger.info(
        "[ligands.coverage] roots=%d noncontrol_candidates=%d filtered_noncontrols=%d controls=%d",
        len(allowed_noncontrol_roots),
        len(noncontrols),
        len(filtered_noncontrols),
        len(controls_valid),
    )

    # When pH ligand mode is disabled, drop pH-annotated microstate ligands
    # from the non-control pool. We detect these by 'pH' in the filename stem.
    ph_mode_raw = str(cfg.get("PH_LIGAND_MODE", "off")).strip().lower()
    ph_ligand_mode_on = ph_mode_raw not in ("", "off", "none", "false", "0")
    if not ph_ligand_mode_on:
        micro_removed: list[Path] = []
        ph_dir_removed: list[Path] = []
        stem_removed: list[Path] = []
        filtered: list[Path] = []
        for p in filtered_noncontrols:
            is_microstates = any(part.lower() == "microstates" for part in p.parts)
            is_ph_subdir = _is_under_ph_subdir(p, allowed_noncontrol_roots)
            has_ph_stem = "pH" in p.stem

            if is_microstates:
                micro_removed.append(p)
                continue
            if is_ph_subdir:
                ph_dir_removed.append(p)
                continue
            if has_ph_stem:
                stem_removed.append(p)
                continue
            filtered.append(p)

        filtered_noncontrols = filtered
        logger.info(
            "[ligands.ph-filter] mode=off microstates_removed=%d ph_subdir_removed=%d stem_removed=%d",
            len(micro_removed),
            len(ph_dir_removed),
            len(stem_removed),
        )
        if micro_removed:
            logger.info(
                "[ligands.ph-filter.examples] category=microstates sample=%s",
                [str(p) for p in micro_removed[:3]],
            )
        if ph_dir_removed:
            logger.info(
                "[ligands.ph-filter.examples] category=ph_subdir sample=%s",
                [str(p) for p in ph_dir_removed[:3]],
            )

    if chunk_keys:
        before = len(filtered_noncontrols)
        filtered_noncontrols = [
            p for p in filtered_noncontrols if _chunk_ligand_key(p.name) in chunk_keys
        ]
        logger.info(
            "[ligands.chunk-filter] chunk_id=%s selected=%d before=%d after=%d",
            str(cfg.get("_CHUNK_ID", "")),
            len(chunk_keys),
            before,
            len(filtered_noncontrols),
        )
    elif cfg.get("ENABLE_DISTRIBUTED_COMBO_CHUNKS") and isinstance(
        cfg.get("_BENCH_SMALL_TARGET_LIGANDS"), int
    ):
        target_n = int(cfg.get("_BENCH_SMALL_TARGET_LIGANDS", 0) or 0)
        if target_n > 0 and len(filtered_noncontrols) > target_n:
            seed = int(cfg.get("_BENCH_SMALL_SEED", 1337) or 1337)
            active_before = sum(
                1 for path in filtered_noncontrols if is_dud_active_key(path.name)
            )
            before = len(filtered_noncontrols)
            filtered_noncontrols = stable_active_preserving_limit(
                filtered_noncontrols,
                limit=target_n,
                seed=seed,
                scope=str(paths.pdb_id),
                key=lambda p: _chunk_ligand_key(p.name),
            )
            active_after = sum(
                1 for path in filtered_noncontrols if is_dud_active_key(path.name)
            )
            logger.info(
                "[ligands.small-profile-cap] pdb=%s before=%d after=%d active_before=%d active_after=%d seed=%d",
                paths.pdb_id,
                before,
                len(filtered_noncontrols),
                active_before,
                active_after,
                seed,
            )

    # --- Optional: build library manifests from scan results ---
    if cfg.get("LIBRARY_MANIFEST_BUILD_ON_SCAN"):
        try:
            manifest_filename = str(
                cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json")
            )
            by_root: Dict[Path, list[Path]] = {}

            for root in allowed_noncontrol_roots:
                root = Path(root)
                if not root.exists():
                    continue
                for lig in filtered_noncontrols:
                    lig_path = Path(lig)
                    if not lig_path.exists():
                        continue
                    if not _under(lig_path, root):
                        continue
                    by_root.setdefault(root, []).append(lig_path)

            for root, ligs in by_root.items():
                manifest_path = root / manifest_filename
                if manifest_path.exists():
                    logger.info(
                        "[lib-manifest.scan.skip] root=%s reason=exists path=%s",
                        str(root),
                        str(manifest_path),
                    )
                    continue

                if not ligs:
                    continue

                logger.info(
                    "[lib-manifest.scan.build] root=%s ligands=%d manifest=%s",
                    str(root),
                    len(ligs),
                    str(manifest_path),
                )

                entries: list[str] = []
                for lig in ligs:
                    try:
                        rel = Path(lig).resolve().relative_to(root.resolve())
                        entries.append(rel.as_posix())
                    except Exception:
                        entries.append(os.path.relpath(str(lig), str(root)))

                tmp_path = manifest_path.with_suffix(manifest_path.suffix + ".tmp")

                try:
                    lib_index = LibraryIndex(
                        manifest_filename=manifest_filename,
                        logger=logger,
                    )
                    if hasattr(lib_index, "write_manifest_for_root"):
                        lib_index.write_manifest_for_root(root, entries, tmp_path)
                    else:
                        import json

                        tmp_path.parent.mkdir(parents=True, exist_ok=True)
                        data = {"root": str(root), "entries": entries}
                        tmp_path.write_text(json.dumps(data, indent=2, sort_keys=True))

                    tmp_path.replace(manifest_path)
                except Exception:
                    logger.exception(
                        "[lib-manifest.scan.error] root=%s manifest=%s",
                        str(root),
                        str(manifest_path),
                    )
        except Exception:
            logger.exception(
                "[lib-manifest.scan.error] unexpected failure during build_on_scan"
            )
    # ------------------------------------------------------------

    # Merge back: controls (unaltered) + filtered non-controls
    if cfg.get("_EFFECTIVE_SINGLE_LIGAND") and cfg.get("_SINGLE_RESOLVED_PATH"):
        resolved_path = Path(cfg["_SINGLE_RESOLVED_PATH"])
        final_paths = controls_valid + [resolved_path]
        filtered_noncontrols = [resolved_path]
        logger.info(
            "[single.fuel] resolved=%s controls=%d (blocking non-control pool)",
            cfg["_SINGLE_RESOLVED_PATH"],
            len(controls_valid),
        )
    else:
        final_paths = controls_valid + filtered_noncontrols

    if chunk_keys:
        cfg["_CHUNK_EFFECTIVE_LIGAND_KEYS"] = sorted(
            {
                key
                for key in (_chunk_ligand_key(p.name) for p in filtered_noncontrols)
                if key
            }
        )
    else:
        cfg.pop("_CHUNK_EFFECTIVE_LIGAND_KEYS", None)

    return _finalize_selected_paths(final_paths, controls_count=len(controls))

def validate_ligand(
    ligand_name: str,
    docked_path: str,
    crystal_path: Optional[str] = None,
    rmsd_thresh: float = 2.0,
    self_rmsd: Optional[float] = None,
    logger=None,
) -> bool:
    """
    Validate ligand docking.
      * If crystal structure available ? use redocking RMSD (hard gate).
      * Otherwise (non-controls) ? self-RMSD is *log-only* (never reject).
    """
    if crystal_path and Path(crystal_path).exists():
        # Prefer coordinate-based Kabsch RMSD for control redock when a crystal ligand is available.
        # This uses pose_validation.compute_redock_rmsd, which operates directly on coordinates,
        # and falls back to the RDKit/MCS-based compute_rmsd if needed.
        redock_rmsd = None
        try:
            redock_rmsd = compute_redock_rmsd(crystal_path, docked_path)
        except Exception as e:
            if logger:
                logger.warning(
                    f"[validate] {ligand_name}: compute_redock_rmsd failed for "
                    f"crystal='{crystal_path}' docked='{docked_path}'; "
                    f"falling back to RDKit/MCS RMSD; err={e!r}"
                )

        # If the coordinate-based RMSD could not be computed (None), fall back to the
        # original RDKit/MCS RMSD implementation to preserve behavior.
        if redock_rmsd is None:
            redock_rmsd = compute_rmsd(crystal_path, docked_path)
        if logger:
            sr = f"{self_rmsd:.2f}" if isinstance(self_rmsd, (int, float)) else "n/a"
            logger.info(
                f"[validate] {ligand_name}: redock_RMSD={redock_rmsd:.2f} A, self_RMSD={sr}"
            )
        if redock_rmsd <= rmsd_thresh:
            return True
        else:
            if logger:
                logger.warning(
                    f"[validate] {ligand_name}: redocking failed (RMSD {redock_rmsd:.2f} A > {rmsd_thresh:.2f})"
                )
            return False

    # Non-controls: log self-RMSD but do not gate on it
    try:
        sr_val = float(self_rmsd) if self_rmsd is not None else None
    except Exception:
        sr_val = None
    if logger:
        sr_txt = f"{sr_val:.2f}" if isinstance(sr_val, (int, float)) else "n/a"
        logger.info(f"[validate] {ligand_name}: self_RMSD={sr_txt} A (LOG-ONLY)")
    return True


__all__ = [
    "select_ligands_for_next",
    "prepare_and_filter_ligands",
    "compute_rmsd",
    "validate_ligand",
    "_coerce_test_map",
    "_resolve_test_mode",
    "_lib_roots_for_pdb",
    "_dedup_index_roots",
    "_count_heavy_atoms_from_pdbqt",
    "_read_any_lig",
    "_is_readable_ref",
]
