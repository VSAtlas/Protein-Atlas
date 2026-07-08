from __future__ import annotations

import json
import logging
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Mapping, Sequence

from cli.cli_utils import _cli_has, _cli_val
from config.value_access import to_bool
from prep_ligands.library_index import LibraryIndex

_MANIFEST_PREFLIGHT_CHECK_WORKERS = 8
_MANIFEST_PREFLIGHT_PROGRESS_EVERY = 5


def _manifest_preflight_root_status(
    root: Path,
    *,
    manifest_filename: str,
) -> tuple[Path, str]:
    path = Path(root)
    if not path.exists():
        return path, "missing_root"
    manifest_path = path / manifest_filename
    if not manifest_path.exists():
        return path, "missing_manifest"
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8")) or {}
    except Exception:
        return path, "manifest_unreadable"
    if not isinstance(payload, dict):
        return path, "manifest_invalid"
    entries = payload.get("entries")
    filenames = payload.get("filenames")
    if not isinstance(entries, dict) or not isinstance(filenames, dict):
        return path, "manifest_invalid"
    return path, "healthy"


def _evaluate_manifest_preflight_roots(
    cfg: Mapping[str, Any],
    *,
    roots: Sequence[Path],
    logger: logging.Logger,
) -> tuple[list[Path], dict[str, int]]:
    manifest_filename = str(cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json"))
    deduped: list[Path] = []
    seen: set[str] = set()
    for root in roots:
        path = Path(root)
        try:
            key = str(path.resolve())
        except Exception:
            key = str(path)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(path)
    summary = {
        "roots_total": len(deduped),
        "healthy": 0,
        "missing_root": 0,
        "missing_manifest": 0,
        "manifest_unreadable": 0,
        "manifest_invalid": 0,
        "broken_total": 0,
    }
    if not deduped:
        return [], summary

    started = time.perf_counter()
    workers = max(
        1,
        min(
            int(_MANIFEST_PREFLIGHT_CHECK_WORKERS),
            len(deduped),
            int(os.cpu_count() or 1),
        ),
    )
    logger.info(
        "[distributed.chunk.manifest-preflight.check.start] roots=%d workers=%d",
        len(deduped),
        workers,
    )
    statuses: list[tuple[Path, str]] = []
    if workers <= 1:
        for idx, root in enumerate(deduped, start=1):
            statuses.append(
                _manifest_preflight_root_status(
                    root,
                    manifest_filename=manifest_filename,
                )
            )
            if (
                idx % int(_MANIFEST_PREFLIGHT_PROGRESS_EVERY) == 0
                or idx == len(deduped)
            ):
                elapsed = time.perf_counter() - started
                logger.info(
                    "[distributed.chunk.manifest-preflight.check.progress] checked=%d total=%d elapsed_s=%.2f",
                    idx,
                    len(deduped),
                    elapsed,
                )
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(
                    _manifest_preflight_root_status,
                    root,
                    manifest_filename=manifest_filename,
                )
                for root in deduped
            ]
            for idx, fut in enumerate(futures, start=1):
                try:
                    statuses.append(fut.result())
                except Exception:
                    statuses.append((deduped[idx - 1], "manifest_unreadable"))
                if (
                    idx % int(_MANIFEST_PREFLIGHT_PROGRESS_EVERY) == 0
                    or idx == len(deduped)
                ):
                    elapsed = time.perf_counter() - started
                    logger.info(
                        "[distributed.chunk.manifest-preflight.check.progress] checked=%d total=%d elapsed_s=%.2f",
                        idx,
                        len(deduped),
                        elapsed,
                    )

    broken_roots: list[Path] = []
    for root, status in statuses:
        if status == "healthy":
            summary["healthy"] += 1
            continue
        summary[status] = int(summary.get(status, 0)) + 1
        if status != "missing_root":
            broken_roots.append(root)
    summary["broken_total"] = len(broken_roots)
    elapsed = time.perf_counter() - started
    logger.info(
        "[distributed.chunk.manifest-preflight.check.done] roots=%d healthy=%d broken=%d missing_root=%d elapsed_s=%.2f",
        int(summary.get("roots_total", 0)),
        int(summary.get("healthy", 0)),
        int(summary.get("broken_total", 0)),
        int(summary.get("missing_root", 0)),
        elapsed,
    )
    return broken_roots, summary


def _rebuild_library_manifests(
    cfg: Mapping[str, Any],
    *,
    roots: Sequence[Path],
    logger: logging.Logger,
) -> dict[str, int]:
    manifest_filename = str(cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json"))
    summary = {
        "roots_total": 0,
        "built": 0,
        "rebuilt": 0,
        "unchanged": 0,
        "missing_root": 0,
        "failed": 0,
    }
    if not roots:
        return summary
    deduped: list[Path] = []
    seen: set[str] = set()
    for root in roots:
        path = Path(root)
        try:
            key = str(path.resolve())
        except Exception:
            key = str(path)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(path)
    pre_state: dict[Path, tuple[bool, float]] = {}
    usable_roots: list[Path] = []
    for root in deduped:
        if not root.exists():
            summary["missing_root"] += 1
            continue
        usable_roots.append(root)
        manifest_path = root / manifest_filename
        exists = manifest_path.exists()
        mtime = 0.0
        if exists:
            try:
                mtime = float(manifest_path.stat().st_mtime)
            except Exception:
                mtime = 0.0
        pre_state[root] = (exists, mtime)
    summary["roots_total"] = len(usable_roots)
    if not usable_roots:
        return summary
    try:
        idx = LibraryIndex(manifest_filename=manifest_filename, logger=logger)
        idx.load(usable_roots)
    except Exception:
        logger.exception("[lib-manifest.rebuild] unexpected load failure")
        summary["failed"] = len(usable_roots)
        return summary

    for root in usable_roots:
        manifest_path = root / manifest_filename
        pre_exists, pre_mtime = pre_state.get(root, (False, 0.0))
        if not manifest_path.exists():
            summary["failed"] += 1
            continue
        try:
            post_mtime = float(manifest_path.stat().st_mtime)
        except Exception:
            post_mtime = pre_mtime
        if not pre_exists:
            summary["built"] += 1
        elif abs(post_mtime - pre_mtime) > 1e-6:
            summary["rebuilt"] += 1
        else:
            summary["unchanged"] += 1
    return summary


def _resolve_manifest_rebuild_roots(
    cfg: Mapping[str, Any], argv: Sequence[str], logger: logging.Logger
) -> list[Path]:
    rebuild_root_raw = _cli_val(list(argv), "--rebuild-root")
    if rebuild_root_raw:
        rebuild_root = Path(str(rebuild_root_raw)).expanduser().resolve()
        if not rebuild_root.exists() or not rebuild_root.is_dir():
            print(
                f"ERROR: --rebuild-root path does not exist or is not a directory: {rebuild_root}",
                file=sys.stderr,
            )
            sys.exit(2)
        base_prepped = rebuild_root / "prepped_ligands"
    else:
        base_prepped = Path(
            str(
                cfg.get("PREPPED_LIGANDS_DIR")
                or cfg.get("OUTPUT_LIGANDS_DIR")
                or "prepped_ligands"
            )
        )
    base_prepped = base_prepped.expanduser().resolve()
    if not base_prepped.exists() or not base_prepped.is_dir():
        logger.warning(
            "[lib-manifest.rebuild] prepped_root_missing=%s action=none",
            str(base_prepped),
        )
        return []

    roots: list[Path] = []
    seen: set[str] = set()
    for child in sorted(base_prepped.iterdir()):
        if not child.is_dir():
            continue
        name = child.name.strip()
        if not name or name.startswith(".") or name == "__pycache__":
            continue
        key = str(child.resolve())
        if key in seen:
            continue
        seen.add(key)
        roots.append(child)
    logger.info(
        "[lib-manifest.rebuild] prepped_root=%s discovered_roots=%d",
        str(base_prepped),
        len(roots),
    )
    return roots


def maybe_run_manifest_rebuild_only(
    cfg: Mapping[str, Any], argv: Sequence[str]
) -> bool:
    if not _cli_has(list(argv), "--rebuild"):
        return False
    logger = logging.getLogger("lib-manifest")
    started = time.perf_counter()
    roots = _resolve_manifest_rebuild_roots(cfg, argv, logger)
    summary = _rebuild_library_manifests(
        cfg,
        roots=roots,
        logger=logger,
    )
    elapsed = time.perf_counter() - started
    print(
        "[lib-manifest.rebuild] "
        f"roots={summary.get('roots_total', 0)} "
        f"built={summary.get('built', 0)} "
        f"rebuilt={summary.get('rebuilt', 0)} "
        f"unchanged={summary.get('unchanged', 0)} "
        f"missing_root={summary.get('missing_root', 0)} "
        f"failed={summary.get('failed', 0)} "
        f"elapsed_s={elapsed:.2f}"
    )
    return True


def _manifest_prebuild_enabled(cfg: Mapping[str, Any], *, distributed_enabled: bool) -> bool:
    default = bool(distributed_enabled)
    raw = cfg.get("CHUNK_PLANNER_PREBUILD_MANIFESTS", default)
    env = os.environ.get("ATLAS_CHUNK_PLANNER_PREBUILD_MANIFESTS")
    if env is not None:
        raw = env
    return bool(to_bool(raw, default=None))
