#!/usr/bin/env python3
"""Measure file artifacts under one or more roots (stdlib only).

Modes:
- Tree mode (default): measure current filesystem tree(s)
- Compare mode: compare two roots (--compare A B)
- Run snapshot mode: --run-id + --capture baseline|final
- Run diff mode: --run-id + --diff

Run mode stores all state and reports under size_artifacts/<run_id>/ unless --out is provided
for report files.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import heapq
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple, cast


def _norm_rel_path(path: str) -> str:
    return path.replace("\\", "/").lstrip("/")


def _norm_pattern(pattern: str) -> str:
    return pattern.replace("\\", "/")


def _match_excluded(rel_path: str, patterns: Sequence[str], is_dir: bool = False) -> bool:
    if not patterns:
        return False
    rel_norm = _norm_rel_path(rel_path)
    candidates = [rel_norm, f"/{rel_norm}"]
    if is_dir:
        candidates.append(f"{rel_norm}/")
        candidates.append(f"/{rel_norm}/")
    for raw_pattern in patterns:
        pattern = _norm_pattern(raw_pattern)
        for candidate in candidates:
            if fnmatch.fnmatch(candidate, pattern):
                return True
    return False


def _match_excluded_with_root(root_label: str, rel_path: str, patterns: Sequence[str], is_dir: bool = False) -> bool:
    if _match_excluded(rel_path, patterns, is_dir=is_dir):
        return True
    root_norm = _norm_rel_path(root_label)
    if not root_norm:
        return False
    combined = f"{root_norm}/{_norm_rel_path(rel_path)}" if rel_path else root_norm
    return _match_excluded(combined, patterns, is_dir=is_dir)


def walk_tree(
    root: str | os.PathLike[str],
    follow_symlinks: bool = False,
    exclude_patterns: Sequence[str] | None = None,
    quiet: bool = False,
) -> Iterator[Tuple[str, str, int, bool]]:
    """Iteratively walk root using os.scandir.

    Yields tuples: (rel_path, abs_path, size_bytes, is_file).
    """

    root_abs = os.path.abspath(os.fspath(root))
    patterns = list(exclude_patterns or [])

    def warn(message: str) -> None:
        if not quiet:
            print(f"[measure_artifacts] warning: {message}", file=sys.stderr)

    stack: List[Tuple[str, str]] = [("", root_abs)]
    seen_dirs: set[Tuple[int, int]] = set()

    if follow_symlinks:
        try:
            st_root = os.stat(root_abs, follow_symlinks=True)
            seen_dirs.add((st_root.st_dev, st_root.st_ino))
        except OSError as exc:
            warn(f"cannot stat root '{root_abs}': {exc}")
            return

    while stack:
        rel_dir, abs_dir = stack.pop()
        try:
            with os.scandir(abs_dir) as entries:
                for entry in entries:
                    rel_path = entry.name if not rel_dir else f"{rel_dir}/{entry.name}"

                    try:
                        is_dir = entry.is_dir(follow_symlinks=follow_symlinks)
                    except OSError as exc:
                        warn(f"cannot inspect '{entry.path}': {exc}")
                        continue

                    if _match_excluded(rel_path, patterns, is_dir=is_dir):
                        continue

                    if is_dir:
                        if follow_symlinks:
                            try:
                                st_dir = entry.stat(follow_symlinks=True)
                            except OSError as exc:
                                warn(f"cannot stat directory '{entry.path}': {exc}")
                                continue
                            key = (st_dir.st_dev, st_dir.st_ino)
                            if key in seen_dirs:
                                continue
                            seen_dirs.add(key)
                        stack.append((rel_path, entry.path))
                        yield (rel_path, entry.path, 0, False)
                        continue

                    try:
                        is_file = entry.is_file(follow_symlinks=follow_symlinks)
                    except OSError as exc:
                        warn(f"cannot inspect file type '{entry.path}': {exc}")
                        continue

                    if not is_file:
                        continue

                    try:
                        st_file = entry.stat(follow_symlinks=follow_symlinks)
                    except OSError as exc:
                        warn(f"cannot stat file '{entry.path}': {exc}")
                        continue

                    yield (rel_path, entry.path, int(st_file.st_size), True)
        except OSError as exc:
            warn(f"cannot read directory '{abs_dir}': {exc}")


def _ext_for_path(path_value: str) -> str:
    _, ext = os.path.splitext(path_value)
    return ext.lower() if ext else "(none)"


def _depth_bucket(rel_path: str, max_depth: int) -> str:
    parts = _norm_rel_path(rel_path).split("/")
    dirs = parts[:-1]
    if max_depth <= 0 or not dirs:
        return "."
    depth = min(len(dirs), max_depth)
    if depth <= 0:
        return "."
    return "/".join(dirs[:depth])


def _file_parent_dir(rel_path: str) -> str:
    parts = _norm_rel_path(rel_path).split("/")
    dirs = parts[:-1]
    if not dirs:
        return "."
    return "/".join(dirs)


def _iter_ancestors(rel_dir: str) -> Iterable[str]:
    if rel_dir == ".":
        return []
    parts = rel_dir.split("/")
    return ["/".join(parts[: idx + 1]) for idx in range(len(parts))]


def _write_csv(path: Path, headers: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(list(headers))
        writer.writerows(rows)


def _as_int(value: object) -> int:
    return cast(int, value)


def _sorted_stat_items(stats: Dict[str, Dict[str, int]]) -> List[Tuple[str, Dict[str, int]]]:
    return sorted(stats.items(), key=lambda kv: (-kv[1]["bytes"], -kv[1]["count"], kv[0]))


def _join_root_rel(root_label: str, rel_path: str) -> str:
    root_norm = _norm_rel_path(root_label)
    rel_norm = _norm_rel_path(rel_path)
    if not root_norm:
        return rel_norm
    if not rel_norm:
        return root_norm
    return f"{root_norm}/{rel_norm}"


def _build_analysis_from_records(
    records: Iterable[Tuple[str, str, int]],
    *,
    dir_depth: int,
    top_n: int,
    dir_depth_full: bool,
    verbose: bool,
) -> Dict[str, object]:
    total_files = 0
    total_bytes = 0
    ext_stats: Dict[str, Dict[str, int]] = {}
    dir_stats: Dict[str, Dict[str, int]] = {".": {"count": 0, "bytes": 0}}
    full_dir_stats: Dict[str, Dict[str, int]] = {}
    file_heap: List[Tuple[int, str, str]] = []

    for root_label, rel_path, size in records:
        total_files += 1
        total_bytes += size

        display_path = _join_root_rel(root_label, rel_path)
        ext = _ext_for_path(display_path)
        e = ext_stats.setdefault(ext, {"count": 0, "bytes": 0})
        e["count"] += 1
        e["bytes"] += size

        dir_stats["."]["count"] += 1
        dir_stats["."]["bytes"] += size
        bucket = _depth_bucket(display_path, dir_depth)
        if bucket != ".":
            d = dir_stats.setdefault(bucket, {"count": 0, "bytes": 0})
            d["count"] += 1
            d["bytes"] += size

        parent_dir = _file_parent_dir(display_path)
        if parent_dir != ".":
            if dir_depth_full:
                ancestors = _iter_ancestors(parent_dir)
            else:
                ancestors = [bucket] if bucket != "." else []
            for anc in ancestors:
                fd = full_dir_stats.setdefault(anc, {"count": 0, "bytes": 0})
                fd["count"] += 1
                fd["bytes"] += size

        file_item = (size, display_path, ext)
        if len(file_heap) < top_n:
            heapq.heappush(file_heap, file_item)
        elif size > file_heap[0][0]:
            heapq.heapreplace(file_heap, file_item)

        if verbose and total_files % 100000 == 0:
            print(
                f"[measure_artifacts] processed {total_files:,} files ({total_bytes:,} bytes)",
                file=sys.stderr,
            )

    top_files = [
        {"path": path, "size": size, "extension": ext}
        for size, path, ext in sorted(file_heap, key=lambda item: (-item[0], item[1]))
    ]

    if dir_depth_full:
        top_dir_source = full_dir_stats
        dir_scope = "full_ancestors"
    else:
        top_dir_source = {k: v for k, v in dir_stats.items() if k != "."}
        dir_scope = f"depth_limited_{dir_depth}"

    top_dirs = [
        {"directory": d, "bytes": vals["bytes"], "count": vals["count"]}
        for d, vals in sorted(
            top_dir_source.items(),
            key=lambda kv: (-kv[1]["bytes"], -kv[1]["count"], kv[0]),
        )[:top_n]
    ]

    return {
        "root": "multi-root",
        "roots": [],
        "total_files": total_files,
        "total_bytes": total_bytes,
        "extension_stats": ext_stats,
        "directory_stats": dir_stats,
        "top_files": top_files,
        "top_directories": top_dirs,
        "top_directories_scope": dir_scope,
    }


def _collect_tree_records(
    roots: Sequence[str],
    *,
    follow_symlinks: bool,
    exclude_patterns: Sequence[str],
    quiet: bool,
) -> List[Tuple[str, str, int]]:
    records: List[Tuple[str, str, int]] = []
    for root in roots:
        root_abs = os.path.abspath(root)
        for rel_path, abs_path, size, is_file in walk_tree(
            root_abs,
            follow_symlinks=follow_symlinks,
            exclude_patterns=exclude_patterns,
            quiet=quiet,
        ):
            if not is_file:
                continue
            if _match_excluded_with_root(root, rel_path, exclude_patterns, is_dir=False):
                continue
            records.append((root, rel_path, size))
    return records


def analyze_roots(
    roots: Sequence[str],
    *,
    dir_depth: int,
    top_n: int,
    follow_symlinks: bool,
    exclude_patterns: Sequence[str],
    quiet: bool,
    verbose: bool,
    dir_depth_full: bool,
) -> Dict[str, object]:
    records = _collect_tree_records(
        roots,
        follow_symlinks=follow_symlinks,
        exclude_patterns=exclude_patterns,
        quiet=quiet,
    )
    result = _build_analysis_from_records(
        records,
        dir_depth=dir_depth,
        top_n=top_n,
        dir_depth_full=dir_depth_full,
        verbose=verbose,
    )
    result["roots"] = [os.path.abspath(r) for r in roots]
    if len(roots) == 1:
        result["root"] = os.path.abspath(roots[0])
    return result


def analyze_root(
    root: str | os.PathLike[str],
    *,
    dir_depth: int,
    top_n: int,
    follow_symlinks: bool,
    exclude_patterns: Sequence[str],
    quiet: bool,
    verbose: bool,
    dir_depth_full: bool,
) -> Dict[str, object]:
    root_str = os.fspath(root)
    records: List[Tuple[str, str, int]] = []
    for rel_path, _abs_path, size, is_file in walk_tree(
        root_str,
        follow_symlinks=follow_symlinks,
        exclude_patterns=exclude_patterns,
        quiet=quiet,
    ):
        if not is_file:
            continue
        records.append(("", rel_path, size))
    result = _build_analysis_from_records(
        records,
        dir_depth=dir_depth,
        top_n=top_n,
        dir_depth_full=dir_depth_full,
        verbose=verbose,
    )
    result["roots"] = [os.path.abspath(root_str)]
    result["root"] = os.path.abspath(root_str)
    return result


def _print_summary(result: Dict[str, object], top_print: int = 10) -> None:
    roots = cast(List[str], result.get("roots", []))
    if roots:
        if len(roots) == 1:
            print(f"Root: {roots[0]}")
        else:
            print("Roots:")
            for root in roots:
                print(f"  - {root}")
    else:
        print(f"Root: {result['root']}")

    print(f"Total files: {_as_int(result['total_files']):,}")
    print(f"Total bytes: {_as_int(result['total_bytes']):,}")

    ext_stats = result["extension_stats"]
    assert isinstance(ext_stats, dict)
    print("\nTop extensions by bytes:")
    for ext, vals in _sorted_stat_items(ext_stats)[:top_print]:
        print(f"  {ext:>10}  files={vals['count']:,}  bytes={vals['bytes']:,}")

    dir_stats = result["directory_stats"]
    assert isinstance(dir_stats, dict)
    print("\nTop directories by bytes:")
    for d, vals in _sorted_stat_items({k: v for k, v in dir_stats.items() if k != "."})[:top_print]:
        print(f"  {d:<40} files={vals['count']:,}  bytes={vals['bytes']:,}")

    print("\nTop files by size:")
    top_files = result["top_files"]
    assert isinstance(top_files, list)
    for item in top_files[:top_print]:
        print(f"  {item['path']:<60} size={int(item['size']):,}  ext={item['extension']}")


def _write_single_reports(out_dir: Path, result: Dict[str, object], metadata: Dict[str, object]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_payload = {
        "metadata": metadata,
        "totals": {
            "files": result["total_files"],
            "bytes": result["total_bytes"],
        },
        "roots": result.get("roots", []),
        "top_directories_scope": result["top_directories_scope"],
        "top_files": result["top_files"],
        "top_directories": result["top_directories"],
    }
    (out_dir / "summary.json").write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")

    ext_stats = result["extension_stats"]
    assert isinstance(ext_stats, dict)
    ext_rows = [[ext, vals["count"], vals["bytes"]] for ext, vals in _sorted_stat_items(ext_stats)]
    _write_csv(out_dir / "by_extension.csv", ["extension", "count", "bytes"], ext_rows)

    dir_stats = result["directory_stats"]
    assert isinstance(dir_stats, dict)
    dir_rows = [[d, vals["count"], vals["bytes"]] for d, vals in _sorted_stat_items(dir_stats)]
    _write_csv(out_dir / "by_dir.csv", ["directory", "count", "bytes"], dir_rows)

    top_files = result["top_files"]
    assert isinstance(top_files, list)
    top_file_rows = [[item["path"], item["size"], item["extension"]] for item in top_files]
    _write_csv(out_dir / "top_files.csv", ["path", "size", "extension"], top_file_rows)


def _compute_delta_rows(
    stats_a: Dict[str, Dict[str, int]], stats_b: Dict[str, Dict[str, int]]
) -> List[Tuple[str, int, int, int, int, int, int]]:
    keys = sorted(set(stats_a) | set(stats_b))
    rows: List[Tuple[str, int, int, int, int, int, int]] = []
    for key in keys:
        a = stats_a.get(key, {"count": 0, "bytes": 0})
        b = stats_b.get(key, {"count": 0, "bytes": 0})
        rows.append((
            key,
            a["count"],
            a["bytes"],
            b["count"],
            b["bytes"],
            b["count"] - a["count"],
            b["bytes"] - a["bytes"],
        ))
    rows.sort(key=lambda row: (-abs(row[6]), -abs(row[5]), row[0]))
    return rows


def _write_compare_reports(
    out_dir: Path,
    result_a: Dict[str, object],
    result_b: Dict[str, object],
    metadata: Dict[str, object],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    ext_a = cast(Dict[str, Dict[str, int]], result_a["extension_stats"])
    ext_b = cast(Dict[str, Dict[str, int]], result_b["extension_stats"])
    dir_a = cast(Dict[str, Dict[str, int]], result_a["directory_stats"])
    dir_b = cast(Dict[str, Dict[str, int]], result_b["directory_stats"])

    rows_ext = _compute_delta_rows(ext_a, ext_b)
    rows_dir = _compute_delta_rows(dir_a, dir_b)

    compare_payload = {
        "metadata": metadata,
        "delta": {
            "files": _as_int(result_b["total_files"]) - _as_int(result_a["total_files"]),
            "bytes": _as_int(result_b["total_bytes"]) - _as_int(result_a["total_bytes"]),
        },
        "roots": {
            "a": result_a["root"],
            "b": result_b["root"],
        },
        "top_extension_deltas": rows_ext[:10],
        "top_directory_deltas": rows_dir[:10],
    }
    (out_dir / "compare.json").write_text(json.dumps(compare_payload, indent=2), encoding="utf-8")

    _write_csv(
        out_dir / "compare_by_extension.csv",
        [
            "extension",
            "count_a",
            "bytes_a",
            "count_b",
            "bytes_b",
            "delta_count",
            "delta_bytes",
        ],
        rows_ext,
    )
    _write_csv(
        out_dir / "compare_by_dir.csv",
        [
            "directory",
            "count_a",
            "bytes_a",
            "count_b",
            "bytes_b",
            "delta_count",
            "delta_bytes",
        ],
        rows_dir,
    )


def _print_compare_summary(result_a: Dict[str, object], result_b: Dict[str, object], top_print: int = 10) -> None:
    delta_files = _as_int(result_b["total_files"]) - _as_int(result_a["total_files"])
    delta_bytes = _as_int(result_b["total_bytes"]) - _as_int(result_a["total_bytes"])

    print(f"Compare (B - A):\n  A={result_a['root']}\n  B={result_b['root']}")
    print(f"File count delta: {delta_files:+,}")
    print(f"Bytes delta:      {delta_bytes:+,}")

    ext_a = cast(Dict[str, Dict[str, int]], result_a["extension_stats"])
    ext_b = cast(Dict[str, Dict[str, int]], result_b["extension_stats"])
    dir_a = cast(Dict[str, Dict[str, int]], result_a["directory_stats"])
    dir_b = cast(Dict[str, Dict[str, int]], result_b["directory_stats"])

    ext_rows = _compute_delta_rows(ext_a, ext_b)
    dir_rows = _compute_delta_rows(dir_a, dir_b)

    print("\nTop extension deltas by absolute bytes:")
    for row in ext_rows[:top_print]:
        print(f"  {row[0]:>10}  d_count={row[5]:+}  d_bytes={row[6]:+}")

    print("\nTop directory deltas by absolute bytes:")
    for row in dir_rows[:top_print]:
        print(f"  {row[0]:<40} d_count={row[5]:+}  d_bytes={row[6]:+}")


def _run_state_dir(run_id: str) -> Path:
    return Path("size_artifacts") / run_id


def _snapshot_path(run_id: str, kind: str) -> Path:
    return _run_state_dir(run_id) / f"{kind}.jsonl"


def _capture_snapshot(
    run_id: str,
    kind: str,
    roots: Sequence[str],
    *,
    follow_symlinks: bool,
    exclude_patterns: Sequence[str],
    quiet: bool,
    verbose: bool,
    max_records: int | None,
) -> Path:
    out_path = _snapshot_path(run_id, kind)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    warned_cap = False
    with out_path.open("w", encoding="utf-8") as handle:
        for root in roots:
            root_abs = os.path.abspath(root)
            for rel_path, abs_path, size, is_file in walk_tree(
                root_abs,
                follow_symlinks=follow_symlinks,
                exclude_patterns=exclude_patterns,
                quiet=quiet,
            ):
                if not is_file:
                    continue
                if _match_excluded_with_root(root, rel_path, exclude_patterns, is_dir=False):
                    continue

                if max_records is not None and count >= max_records:
                    if not warned_cap and not quiet:
                        print(
                            f"[measure_artifacts] warning: max records reached ({max_records}); truncating snapshot",
                            file=sys.stderr,
                        )
                        warned_cap = True
                    continue

                try:
                    st = os.stat(abs_path, follow_symlinks=follow_symlinks)
                except OSError as exc:
                    if not quiet:
                        print(f"[measure_artifacts] warning: cannot stat file '{abs_path}': {exc}", file=sys.stderr)
                    continue

                rec = {
                    "root": root,
                    "rel": _norm_rel_path(rel_path),
                    "size": int(size),
                    "mtime_ns": int(st.st_mtime_ns),
                }
                handle.write(json.dumps(rec) + "\n")
                count += 1
                if verbose and count % 100000 == 0:
                    print(f"[measure_artifacts] snapshot {kind}: {count:,} files", file=sys.stderr)

    if verbose:
        print(f"[measure_artifacts] snapshot {kind} wrote {count:,} file records", file=sys.stderr)
    return out_path


def _load_snapshot(path: Path, roots_filter: set[str] | None = None) -> Dict[Tuple[str, str], Tuple[int, int]]:
    result: Dict[Tuple[str, str], Tuple[int, int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            root = str(rec["root"])
            if roots_filter is not None and root not in roots_filter:
                continue
            rel = _norm_rel_path(str(rec["rel"]))
            size = int(rec["size"])
            mtime_ns = int(rec["mtime_ns"])
            result[(root, rel)] = (size, mtime_ns)
    return result


def _diff_snapshots(
    baseline: Dict[Tuple[str, str], Tuple[int, int]],
    final: Dict[Tuple[str, str], Tuple[int, int]],
) -> Tuple[
    List[Tuple[str, str, int]],
    List[Tuple[str, str, int, int]],
    List[Tuple[str, str, int]],
]:
    created: List[Tuple[str, str, int]] = []
    modified: List[Tuple[str, str, int, int]] = []
    deleted: List[Tuple[str, str, int]] = []

    for key, (new_size, new_mtime) in final.items():
        old = baseline.get(key)
        if old is None:
            created.append((key[0], key[1], new_size))
            continue
        old_size, old_mtime = old
        if old_size != new_size or old_mtime != new_mtime:
            modified.append((key[0], key[1], old_size, new_size))

    for key, (old_size, _) in baseline.items():
        if key not in final:
            deleted.append((key[0], key[1], old_size))

    created.sort(key=lambda x: (-x[2], x[0], x[1]))
    modified.sort(key=lambda x: (-x[3], x[0], x[1]))
    deleted.sort(key=lambda x: (-x[2], x[0], x[1]))
    return created, modified, deleted


def _write_run_diff_reports(
    out_dir: Path,
    *,
    result: Dict[str, object],
    metadata: Dict[str, object],
    created: List[Tuple[str, str, int]],
    modified: List[Tuple[str, str, int, int]],
    deleted: List[Tuple[str, str, int]],
    top_n: int,
) -> None:
    _write_single_reports(out_dir, result, metadata)

    created_rows = [[_norm_rel_path(rel), size, root] for root, rel, size in created]
    modified_rows = [[_norm_rel_path(rel), old_size, new_size, root] for root, rel, old_size, new_size in modified]
    _write_csv(out_dir / "created_files.csv", ["path", "size", "root"], created_rows)
    _write_csv(
        out_dir / "modified_files.csv",
        ["path", "old_size", "new_size", "root"],
        modified_rows,
    )

    deleted_rows_all = [[_norm_rel_path(rel), old_size, root] for root, rel, old_size in deleted]
    deleted_rows = deleted_rows_all[:top_n]
    _write_csv(out_dir / "deleted_files.csv", ["path", "old_size", "root"], deleted_rows)
    if len(deleted_rows_all) > len(deleted_rows):
        deleted_meta = {
            "deleted_total": len(deleted_rows_all),
            "deleted_rows_written": len(deleted_rows),
            "truncated": True,
        }
        (out_dir / "deleted_count.json").write_text(json.dumps(deleted_meta, indent=2), encoding="utf-8")


def _print_run_diff_summary(
    run_id: str,
    result: Dict[str, object],
    created: List[Tuple[str, str, int]],
    modified: List[Tuple[str, str, int, int]],
    deleted: List[Tuple[str, str, int]],
) -> None:
    print(f"Run ID: {run_id}")
    print(f"Created files: {len(created):,}")
    print(f"Modified files: {len(modified):,}")
    print(f"Deleted files: {len(deleted):,}")
    print(f"Run footprint files (created+modified): {_as_int(result['total_files']):,}")
    print(f"Run footprint bytes (new size): {_as_int(result['total_bytes']):,}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Measure artifact files under a root directory")

    parser.add_argument("--root", action="append", default=[], help="Root directory to analyze; repeatable")
    parser.add_argument("--roots-from", help="File with newline-separated roots to analyze")
    parser.add_argument("--compare", nargs=2, metavar=("ROOT_A", "ROOT_B"), help="Compare ROOT_B against ROOT_A")

    parser.add_argument("--run-id", help="Run identifier for run-scoped snapshots/diff")
    parser.add_argument("--capture", choices=["baseline", "final"], help="Capture run snapshot")
    parser.add_argument("--diff", action="store_true", help="Diff baseline/final snapshots for run-id")
    parser.add_argument("--max-records", type=int, help="Optional cap for snapshot records")

    parser.add_argument("--out", help="Output report directory")
    parser.add_argument("--dir-depth", type=int, default=3, help="Directory aggregation depth (default: 3)")
    parser.add_argument("--dir-depth-full", action="store_true", help="Use full ancestor accumulation for largest directories")
    parser.add_argument("--top", type=int, default=50, help="Top-N files/directories (default: 50)")
    parser.add_argument("--follow-symlinks", action="store_true", help="Follow symlinks during traversal")
    parser.add_argument("--exclude", action="append", default=[], help="Glob-like exclude pattern; repeatable")
    parser.add_argument("--quiet", action="store_true", help="Suppress warnings")
    parser.add_argument("--verbose", action="store_true", help="Verbose progress output")
    return parser


def _resolve_roots(root_args: Sequence[str], roots_from: str | None) -> List[str]:
    roots: List[str] = list(root_args)
    if roots_from:
        path = Path(roots_from)
        for raw in path.read_text(encoding="utf-8").splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            roots.append(line)

    if not roots:
        return ["docked"]

    seen: set[str] = set()
    ordered: List[str] = []
    for root in roots:
        key = root.strip()
        if not key or key in seen:
            continue
        seen.add(key)
        ordered.append(key)
    return ordered


def _metadata(args: argparse.Namespace, roots: Sequence[str]) -> Dict[str, object]:
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(),
        "python": sys.version,
        "roots": [os.path.abspath(r) for r in roots],
        "args": vars(args),
        "exclude_matching": "patterns are matched against both rel and root/rel normalized paths",
    }


def _validate_roots_exist(parser: argparse.ArgumentParser, roots: Sequence[str]) -> None:
    for root in roots:
        if not os.path.isdir(root):
            parser.error(f"root is not a directory: {root}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.dir_depth < 0:
        parser.error("--dir-depth must be >= 0")
    if args.top <= 0:
        parser.error("--top must be > 0")
    if args.max_records is not None and args.max_records <= 0:
        parser.error("--max-records must be > 0")

    roots = _resolve_roots(args.root, args.roots_from)

    if args.compare:
        if args.run_id or args.capture or args.diff:
            parser.error("--compare cannot be combined with --run-id/--capture/--diff")
        root_a, root_b = args.compare
        if not os.path.isdir(root_a):
            parser.error(f"ROOT_A is not a directory: {root_a}")
        if not os.path.isdir(root_b):
            parser.error(f"ROOT_B is not a directory: {root_b}")

        out_dir = Path(args.out) if args.out else Path(root_a) / "_artifact_reports"

        result_a = analyze_roots(
            [root_a],
            dir_depth=args.dir_depth,
            top_n=args.top,
            follow_symlinks=args.follow_symlinks,
            exclude_patterns=args.exclude,
            quiet=args.quiet,
            verbose=args.verbose,
            dir_depth_full=args.dir_depth_full,
        )
        result_b = analyze_roots(
            [root_b],
            dir_depth=args.dir_depth,
            top_n=args.top,
            follow_symlinks=args.follow_symlinks,
            exclude_patterns=args.exclude,
            quiet=args.quiet,
            verbose=args.verbose,
            dir_depth_full=args.dir_depth_full,
        )

        _print_compare_summary(result_a, result_b)
        metadata = _metadata(args, [root_a, root_b])
        _write_compare_reports(out_dir, result_a, result_b, metadata)
        print(f"\nCompare reports written to: {out_dir}")
        return 0

    run_mode = bool(args.run_id)
    if args.capture or args.diff:
        if not args.run_id:
            parser.error("--capture/--diff requires --run-id")
        run_mode = True

    if run_mode and args.run_id is None:
        parser.error("run mode requires --run-id")

    _validate_roots_exist(parser, roots)

    if args.capture:
        snap_path = _capture_snapshot(
            args.run_id,
            args.capture,
            roots,
            follow_symlinks=args.follow_symlinks,
            exclude_patterns=args.exclude,
            quiet=args.quiet,
            verbose=args.verbose,
            max_records=args.max_records,
        )
        print(f"Snapshot written: {snap_path}")
        return 0

    if args.diff:
        assert args.run_id is not None
        baseline_path = _snapshot_path(args.run_id, "baseline")
        final_path = _snapshot_path(args.run_id, "final")
        if not baseline_path.exists():
            parser.error(f"missing baseline snapshot: {baseline_path}")
        if not final_path.exists():
            parser.error(f"missing final snapshot: {final_path}")

        roots_filter = set(roots)
        baseline = _load_snapshot(baseline_path, roots_filter=roots_filter)
        final = _load_snapshot(final_path, roots_filter=roots_filter)
        created, modified, deleted = _diff_snapshots(baseline, final)

        footprint_records: List[Tuple[str, str, int]] = []
        for root, rel, size in created:
            footprint_records.append((root, rel, size))
        for root, rel, _, new_size in modified:
            footprint_records.append((root, rel, new_size))

        result = _build_analysis_from_records(
            footprint_records,
            dir_depth=args.dir_depth,
            top_n=args.top,
            dir_depth_full=args.dir_depth_full,
            verbose=args.verbose,
        )
        result["roots"] = [os.path.abspath(r) for r in roots]

        out_dir = Path(args.out) if args.out else _run_state_dir(args.run_id) / "reports"
        _print_run_diff_summary(args.run_id, result, created, modified, deleted)
        _print_summary(result)
        metadata = _metadata(args, roots)
        metadata["run_id"] = args.run_id
        metadata["diff_counts"] = {
            "created": len(created),
            "modified": len(modified),
            "deleted": len(deleted),
            "footprint": len(footprint_records),
        }
        _write_run_diff_reports(
            out_dir,
            result=result,
            metadata=metadata,
            created=created,
            modified=modified,
            deleted=deleted,
            top_n=args.top,
        )
        print(f"\nRun diff reports written to: {out_dir}")
        return 0

    if args.out:
        out_dir = Path(args.out)
    elif run_mode and args.run_id:
        out_dir = _run_state_dir(args.run_id) / "reports"
    elif len(roots) == 1:
        out_dir = Path(roots[0]) / "_artifact_reports"
    else:
        out_dir = Path("_artifact_reports")

    result = analyze_roots(
        roots,
        dir_depth=args.dir_depth,
        top_n=args.top,
        follow_symlinks=args.follow_symlinks,
        exclude_patterns=args.exclude,
        quiet=args.quiet,
        verbose=args.verbose,
        dir_depth_full=args.dir_depth_full,
    )
    _print_summary(result)
    metadata = _metadata(args, roots)
    _write_single_reports(out_dir, result, metadata)
    print(f"\nReports written to: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
