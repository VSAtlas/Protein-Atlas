# single_ligand_index.py
# Helpers for single-ligand selection and manifest-backed library indexing.

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Dict, Optional

from docking.pockets.active_site_runtime import norm
from prep_ligands.library_index import LibraryIndex
from prep_ligands.fda_mapping_identity import (
    resolved_preferred_name,
    suppresses_legacy_aliases,
)
from path_router import Paths

_SINGLE_ALLOW_PREFIX = False
_SINGLE_SKIP_GLOBAL = True
_SINGLE_MANIFEST_ONLY = True
_SINGLE_ORDER = ["fda_library", "per_protein"]


# [single-index] Deduplicate manifest roots and retain stable ordering.
def _dedupe_manifest_roots(seq) -> list[Path]:
    deduped: list[Path] = []
    seen: set[str] = set()
    for entry in seq:
        if not entry:
            continue
        path_obj = Path(entry)
        try:
            key = str(path_obj.resolve())
        except Exception:
            key = str(path_obj)
        if key not in seen:
            seen.add(key)
            deduped.append(path_obj)
    return deduped


# [single-index] Prime the manifest-backed index before single-ligand lookups.
def _ensure_single_ligand_index(
    cfg: Dict, paths: Paths, logger: logging.Logger
) -> None:
    cfg.setdefault("PREPPED_LIGANDS_DIR", str(paths.prepped_ligands_dir.parent))

    per_roots = _dedupe_manifest_roots(
        [
            getattr(paths, "prepped_ligands_dir", None),
        ]
    )
    library_roots = _dedupe_manifest_roots(
        [
            cfg.get("PREPPED_LIGANDS_DIR"),
        ]
    )

    cfg["_LIB_INDEX_PER_ROOTS"] = [str(p) for p in per_roots]
    cfg["_LIB_INDEX_LIBRARY_ROOTS"] = [str(p) for p in library_roots]

    if not per_roots and not library_roots:
        cfg["_LIB_INDEX"] = None
        return

    manifest_filename = str(cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json"))
    strict_mtime = cfg.get("LIBRARY_MANIFEST_STRICT_MTIME_CHECK", True)
    index = cfg.get("_LIB_INDEX")
    if not isinstance(index, LibraryIndex):
        index = LibraryIndex(
            manifest_filename=manifest_filename,
            logger=logger,
            strict_mtime_check=bool(strict_mtime),
        )
        cfg["_LIB_INDEX"] = index

    index_roots = _dedupe_manifest_roots([*per_roots, *library_roots])
    if index_roots:
        status_parts = []
        for root in index_roots:
            manifest_path = Path(root) / manifest_filename
            manifest_state = "present" if manifest_path.exists() else "missing"
            status_parts.append(f"{manifest_path}={manifest_state}")
        index.load(index_roots)
        logger.info(
            "[single.index] per_roots=%s lib_roots=%s manifests=%s",
            [str(p) for p in per_roots],
            [str(p) for p in library_roots],
            ";".join(status_parts) or "none",
        )
        logger.info(
            "[single.index] Using %d library roots for single-ligand selection",
            len(index_roots),
        )
    else:
        logger.info(
            "[single.index] per_roots=%s lib_roots=%s manifests=none",
            [str(p) for p in per_roots],
            [str(p) for p in library_roots],
        )


def _resolve_single_ligand(
    selector: str, pdb_id: str, cfg: Dict, logger: logging.Logger
) -> Optional[Path]:
    """Resolve SINGLE_LIGAND selector using manifest-backed lookups (if enabled)."""
    if not selector:
        return None

    allow_prefix = _SINGLE_ALLOW_PREFIX
    skip_global = _SINGLE_SKIP_GLOBAL
    manifest_only = _SINGLE_MANIFEST_ONLY
    suggestions_cap = int(cfg.get("SINGLE_LIGAND_SUGGESTIONS", 5) or 5)

    per_roots = _dedupe_manifest_roots(cfg.get("_LIB_INDEX_PER_ROOTS", []))
    library_roots = _dedupe_manifest_roots(cfg.get("_LIB_INDEX_LIBRARY_ROOTS", []))

    lib_index = cfg.get("_LIB_INDEX")
    has_index = isinstance(lib_index, LibraryIndex)

    prepped_ligands_dir = cfg.get("PREPPED_LIGANDS_DIR")
    fda_root = (
        Path(prepped_ligands_dir).joinpath("fda_library")
        if prepped_ligands_dir
        else None
    )
    name_map: Optional[dict[str, set[str]]] = None
    fda_logged = False
    key = _norm_name_key(selector)

    effective_order = list(_SINGLE_ORDER)
    use_manifest = has_index and (
        bool(per_roots) or (not skip_global and bool(library_roots))
    )

    logger.info(
        "[single.debug] selector=%s raw_order=%s order=%s skip_global=%s manifest_only=%s has_index=%s use_manifest=%s per_roots=%s lib_roots=%s allow_prefix=%s",
        selector,
        "fixed",
        effective_order,
        skip_global,
        manifest_only,
        str(has_index).lower(),
        str(use_manifest).lower(),
        [str(p) for p in per_roots],
        [str(p) for p in library_roots],
        str(allow_prefix).lower(),
    )

    def _ensure_name_map() -> dict[str, set[str]]:
        nonlocal name_map
        if name_map is None:
            name_map = _load_fda_name_map(cfg, logger)
        return name_map

    def _resolve_fda_scope() -> Optional[Path]:
        nonlocal fda_logged
        if fda_root is None or not fda_root.exists():
            return None
        mapping = _ensure_name_map()
        basenames: list[str] = []
        if key:
            if key in mapping:
                basenames.extend(sorted(mapping.get(key, set())))
            if not basenames and allow_prefix:
                for map_key, values in mapping.items():
                    if map_key.startswith(key):
                        basenames.extend(sorted(values))
        basenames = list(dict.fromkeys(basenames))
        if len(basenames) > 1:
            logger.error(
                "[single.name.ambiguous] key=%s candidates=%s action=fail_closed",
                key,
                basenames,
            )
            fda_logged = True
            return None
        probe_target: Path = fda_root if not basenames else fda_root / basenames[0]
        if not fda_logged:
            logger.info(
                "[single.name] key=%s basenames=%s probe=%s",
                key,
                basenames,
                norm(probe_target),
            )
            fda_logged = True
        for basename in basenames:
            candidate = fda_root / basename
            if candidate.exists():
                logger.info(
                    "[single.fda.hit] key=%s basename=%s path=%s",
                    key,
                    basename,
                    norm(candidate),
                )
                logger.info(
                    "[single.lookup.hit] source=fda selector=%s path=%s",
                    selector,
                    norm(candidate),
                )
                return candidate
        return None

    def _suggest_from_index() -> list[str]:
        if isinstance(lib_index, LibraryIndex):
            ordered: list[Path] = []
            seen: set[str] = set()
            for candidate in [*per_roots, *library_roots]:
                try:
                    key_str = str(candidate.resolve())
                except Exception:
                    key_str = str(candidate)
                if key_str not in seen:
                    seen.add(key_str)
                    ordered.append(candidate)
            return lib_index.suggest(selector, ordered, suggestions_cap)
        return []

    def _log_miss(suggestions: list[str]) -> None:
        logger.error(
            "[single.lookup.miss] selector=%s suggestions=[%s]",
            selector,
            ",".join(suggestions),
        )

    hit = _resolve_fda_scope()
    if hit:
        return hit

    if manifest_only and not use_manifest:
        _log_miss(_suggest_from_index())
        return None

    if not has_index:
        _log_miss(_suggest_from_index())
        return None
    assert isinstance(lib_index, LibraryIndex)

    if not use_manifest:
        _log_miss(_suggest_from_index())
        return None

    if per_roots:
        hit = lib_index.lookup(selector, per_roots, allow_prefix=allow_prefix)
        if hit:
            logger.info(
                "[single.lookup.hit] source=manifest scope=per_protein selector=%s path=%s",
                selector,
                norm(hit),
            )
            return hit

    if not skip_global and library_roots:
        hit = lib_index.lookup(selector, library_roots, allow_prefix=allow_prefix)
        if hit:
            logger.info(
                "[single.lookup.hit] source=manifest scope=global selector=%s path=%s",
                selector,
                norm(hit),
            )
            return hit

    _log_miss(_suggest_from_index())
    return None


# --- FDA name mapping (CSV) ---------------------------------------------------
# Lets SINGLE_LIGAND resolve by generic/brand/synonym (e.g., "imatinib", "Gleevec").
_FDA_NAME_MAP_CACHE: dict[str, set[str]] | None = None


def _norm_name_key(s: str) -> str:
    import re

    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def _split_multi_names(v: str) -> list[str]:
    import re

    parts = re.split(r"[|;,/]", v or "")
    return [p.strip() for p in parts if p and p.strip()]


def _load_fda_name_map(cfg: Dict, logger: logging.Logger) -> dict[str, set[str]]:
    """
    Build dict: normalized_name -> {pdbqt_basename, ...}
    CSV must have at least: column 'path' pointing to a *.pdbqt, plus name columns.
    """
    global _FDA_NAME_MAP_CACHE
    if isinstance(_FDA_NAME_MAP_CACHE, dict):
        return _FDA_NAME_MAP_CACHE

    import csv
    from pathlib import Path

    csv_path = (
        os.environ.get("FDA_MAPPING_CSV", "").strip()
        or str(cfg.get("FDA_MAPPING_CSV", "")).strip()
    )
    mapping: dict[str, set[str]] = {}
    if not csv_path:
        _FDA_NAME_MAP_CACHE = {}
        return _FDA_NAME_MAP_CACHE

    p = Path(csv_path)
    if not p.exists():
        logger.info(
            f"[single:name] FDA_MAPPING_CSV not found at {csv_path} (name lookup disabled)."
        )
        _FDA_NAME_MAP_CACHE = {}
        return _FDA_NAME_MAP_CACHE

    try:
        with open(p, newline="", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f)
            for row in reader:
                path = (row.get("path") or "").strip()
                if not path.endswith(".pdbqt"):
                    continue
                base = os.path.basename(path)
                resolved_name = resolved_preferred_name(row)

                # "single" name fields
                singles = [resolved_name, Path(base).stem]
                if not suppresses_legacy_aliases(row):
                    singles.extend(
                        [
                            row.get("display_name", ""),
                            row.get("generic_name", ""),
                            row.get("rxnorm_generic_name", ""),
                            row.get("drugcentral_generic_name", ""),
                            row.get("pubchem_name", ""),
                            row.get("pubchem_record_title", ""),
                        ]
                    )
                # multi-value fields (split)
                multis = []
                if not suppresses_legacy_aliases(row):
                    for col in (
                        "brand_names",
                        "rxnorm_brand_names",
                        "drugcentral_brand_names",
                        "pubchem_synonyms",
                    ):
                        v = row.get(col, "")
                        if v:
                            multis.extend(_split_multi_names(v))

                for nm in [*singles, *multis]:
                    key = _norm_name_key(nm)
                    if key:
                        mapping.setdefault(key, set()).add(base)

        logger.info(f"[single:name] Loaded FDA name map ({len(mapping)} keys) from {p}")
    except Exception as e:
        logger.warning(f"[single:name] Failed to load name map: {e}")
        mapping = {}

    _FDA_NAME_MAP_CACHE = mapping
    return mapping
