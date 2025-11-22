"""Helper for manifest-backed ligand library lookups."""
from __future__ import annotations

import json
import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


@dataclass
class _Manifest:
    entries: dict[str, str]  # normalized key -> relative path
    filenames: dict[str, str]  # lowercase filename -> relative path


class LibraryIndex:
    """Manifest-backed resolver for ligand libraries."""

    def __init__(
        self,
        *,
        manifest_filename: str = "_manifest.json",
        logger: logging.Logger | None = None,
    ) -> None:
        self._manifest_filename = manifest_filename
        self._logger = logger or logging.getLogger(__name__)
        self._cache: dict[Path, _Manifest] = {}
        self._roots: list[Path] = []

    @staticmethod
    def _normalize_exact(selector: str) -> str:
        stem = Path(selector).stem.lower()
        return stem

    @staticmethod
    def _normalize_base(selector: str) -> str:
        stem = LibraryIndex._normalize_exact(selector)
        base = stem.split("_stage")[0]
        return base

    @staticmethod
    def _manifest_path(root: Path, filename: str) -> Path:
        return root / filename

    def load(self, roots: Sequence[Path]) -> None:
        seen: set[str] = set()
        ordered: list[Path] = []
        for root in roots:
            if not root:
                continue
            key = str(Path(root).resolve())
            if key in seen:
                continue
            seen.add(key)
            ordered.append(Path(root))
        self._roots = ordered
        status_entries: list[str] = []
        for root in ordered:
            manifest_path = self._manifest_path(Path(root), self._manifest_filename)
            manifest_state = "present" if manifest_path.exists() else "missing"
            status_entries.append(f"{Path(root)}:{manifest_state}")
        roots_summary = ",".join(status_entries) if status_entries else "none"
        self._logger.info("[lib-index.load] roots=%s", roots_summary)
        self._cache = {}
        for root in ordered:
            self._cache[root] = self._load_one(root)

    def _load_one(self, root: Path) -> _Manifest:
        root = Path(root)
        manifest_path = self._manifest_path(root, self._manifest_filename)
        entries: dict[str, str] = {}
        filenames: dict[str, str] = {}
        try:
            root_stat = root.stat()
        except FileNotFoundError:
            root_stat = None

        current_mtime = getattr(root_stat, "st_mtime", None)
        data: dict | None = None
        if manifest_path.exists():
            try:
                with open(manifest_path, "r", encoding="utf-8") as fh:
                    data = json.load(fh)
            except Exception as exc:
                self._logger.warning(
                    "[lib-index.error] root=%s path=%s reason=%s",
                    str(root),
                    str(manifest_path),
                    exc,
                )
                data = None

        needs_build = True
        if data and isinstance(data, dict):
            recorded = data.get("mtime")
            if recorded is not None and recorded == current_mtime:
                raw_entries = data.get("entries")
                raw_files = data.get("filenames")
                if isinstance(raw_entries, dict) and isinstance(raw_files, dict):
                    entries = {str(k): str(v) for k, v in raw_entries.items()}
                    filenames = {str(k): str(v) for k, v in raw_files.items()}
                    needs_build = False
        if needs_build and current_mtime is not None:
            entries, filenames = self._build_manifest(root, manifest_path, current_mtime)
        return _Manifest(entries=entries, filenames=filenames)

    def _build_manifest(self, root: Path, manifest_path: Path, mtime: float) -> tuple[dict[str, str], dict[str, str]]:
        entries: dict[str, str] = {}
        filenames: dict[str, str] = {}
        count = 0
        start = time.perf_counter()
        for dirpath, _dirnames, filenames_list in os.walk(root):
            for name in filenames_list:
                if not name.lower().endswith(".pdbqt"):
                    continue
                abs_path = Path(dirpath) / name
                rel_path = os.path.relpath(abs_path, root)
                rel_key = rel_path.replace("\\", "/")
                filename_key = name.lower()
                if filename_key not in filenames:
                    filenames[filename_key] = rel_key
                exact_key = self._normalize_exact(name)
                if exact_key not in entries:
                    entries[exact_key] = rel_key
                base_key = self._normalize_base(name)
                if base_key and base_key not in entries:
                    entries[base_key] = rel_key
                count += 1
        try:
            manifest_path.parent.mkdir(parents=True, exist_ok=True)
            with open(manifest_path, "w", encoding="utf-8") as fh:
                json.dump({"mtime": mtime, "entries": entries, "filenames": filenames}, fh, indent=2, sort_keys=True)
        except Exception as exc:
            self._logger.warning(
                "[lib-index.error] root=%s path=%s reason=%s",
                str(root),
                str(manifest_path),
                exc,
            )
        elapsed = time.perf_counter() - start
        self._logger.info("[lib-index.build] root=%s count=%d time_sec=%.3f", str(root), count, elapsed)
        return entries, filenames

    def write_manifest_for_root(
        self,
        root: Path,
        relative_entries: Iterable[str | Path],
        tmp_path: Path | None = None,
    ) -> Path:
        """Write a manifest for ``root`` using already-resolved relative entries.

        Intended for callers that have just scanned a library tree and want to
        persist a manifest without re-crawling the filesystem. ``relative_entries``
        should already be relative to ``root``.
        """

        root = Path(root)
        manifest_path = tmp_path or self._manifest_path(root, self._manifest_filename)
        try:
            root_stat = root.stat()
            mtime = getattr(root_stat, "st_mtime", time.time())
        except Exception:
            mtime = time.time()

        entries: dict[str, str] = {}
        filenames: dict[str, str] = {}
        count = 0
        start = time.perf_counter()
        for rel in relative_entries:
            rel_path = Path(rel)
            rel_key = rel_path.as_posix()
            name = rel_path.name
            filename_key = name.lower()
            if filename_key not in filenames:
                filenames[filename_key] = rel_key
            exact_key = self._normalize_exact(name)
            if exact_key not in entries:
                entries[exact_key] = rel_key
            base_key = self._normalize_base(name)
            if base_key and base_key not in entries:
                entries[base_key] = rel_key
            count += 1

        try:
            manifest_path.parent.mkdir(parents=True, exist_ok=True)
            with open(manifest_path, "w", encoding="utf-8") as fh:
                json.dump({"mtime": mtime, "entries": entries, "filenames": filenames}, fh, indent=2, sort_keys=True)
        except Exception as exc:
            self._logger.warning(
                "[lib-index.error] root=%s path=%s reason=%s",
                str(root),
                str(manifest_path),
                exc,
            )

        elapsed = time.perf_counter() - start
        self._logger.info("[lib-index.build.scan] root=%s count=%d time_sec=%.3f", str(root), count, elapsed)
        return manifest_path

    def lookup(
        self,
        selector: str,
        roots: Sequence[Path] | None = None,
        *,
        allow_prefix: bool = False,
    ) -> Path | None:
        search_roots = list(roots) if roots is not None else self._roots
        exact_key = self._normalize_exact(selector)
        base_key = self._normalize_base(selector)
        keys = [k for k in {exact_key, base_key} if k]
        for root in search_roots:
            manifest = self._cache.get(Path(root))
            if not manifest:
                continue
            for key in keys:
                rel = manifest.entries.get(key)
                if rel:
                    return Path(root) / rel
            if allow_prefix and base_key:
                for key, rel in manifest.entries.items():
                    if key.startswith(base_key):
                        return Path(root) / rel
        return None

    def lookup_filename(self, filename: str, roots: Sequence[Path]) -> Path | None:
        if not filename:
            return None
        key = filename.lower()
        for root in roots:
            manifest = self._cache.get(Path(root))
            if not manifest:
                continue
            rel = manifest.filenames.get(key)
            if rel:
                return Path(root) / rel
        return None

    def suggest(self, selector: str, roots: Iterable[Path], limit: int) -> list[str]:
        import difflib

        selector_key = self._normalize_base(selector)
        candidates: list[str] = []
        seen: set[str] = set()
        for root in roots:
            manifest = self._cache.get(Path(root))
            if not manifest:
                continue
            for key in manifest.entries.keys():
                if key not in seen:
                    seen.add(key)
                    candidates.append(key)
        if not candidates:
            return []
        limit = max(1, int(limit))
        matches = difflib.get_close_matches(selector_key, candidates, n=limit)
        return matches

    @property
    def roots(self) -> Sequence[Path]:
        return list(self._roots)
