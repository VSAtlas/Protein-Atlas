#!/usr/bin/env python3
"""
Summarize run-scoped errors for a specific run_id.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

PDB_RE = re.compile(r"(?<![A-Za-z0-9])[0-9][A-Za-z0-9]{3}(?![A-Za-z0-9])")
TS_RE = re.compile(r"^\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}")
TRACEBACK_START = "Traceback (most recent call last):"
LOG_SUMMARY_DIRNAME = "logs_summary"

LOG_EXTS = {".log", ".err", ".out", ".txt"}
DEEP_LOG_PREFIXES = ("stderr", "stdout")
SKIP_EXTS = {
    ".pdb",
    ".pdbqt",
    ".sdf",
    ".mol2",
    ".gz",
    ".zip",
    ".tar",
    ".tgz",
    ".bz2",
}

FAILURE_MARKERS = [
    "CalledProcessError",
    "non-zero exit status",
    "return code",
    "segmentation fault",
    "killed",
    "aborted",
]


@dataclass(frozen=True)
class FileMeta:
    path: Path
    source: str
    pdb_id: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize run-scoped error logs for a given run_id."
    )
    parser.add_argument("run_id", help="Run ID to summarize.")
    parser.add_argument("--root", default=".", help="Repo root (default: .)")
    parser.add_argument("--deep", action="store_true", help="Enable deep scans.")
    parser.add_argument(
        "--recent-hours",
        type=float,
        default=168,
        help="Only include logs modified within this window (0 disables).",
    )
    parser.add_argument(
        "--tail-kb",
        type=int,
        default=512,
        help="Read only the tail of each file, in KB.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Show top N signatures by frequency.",
    )
    parser.add_argument(
        "--max-examples",
        type=int,
        default=3,
        help="Max examples per signature.",
    )
    parser.add_argument("--json-out", help="Optional JSON output path.")
    return parser.parse_args()


def is_recent(path: Path, recent_hours: float, now: float) -> bool:
    if recent_hours <= 0:
        return True
    try:
        mtime = path.stat().st_mtime
    except FileNotFoundError:
        return False
    cutoff = now - recent_hours * 3600
    return mtime >= cutoff


def tail_text(path: Path, tail_kb: int) -> str:
    max_bytes = max(tail_kb, 1) * 1024
    with path.open("rb") as fh:
        try:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            offset = max(size - max_bytes, 0)
            fh.seek(offset)
        except OSError:
            fh.seek(0)
        data = fh.read()
    return data.decode("utf-8", errors="replace")


def normalize_signature(text: str) -> str:
    text = re.sub(r"\b\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}[.,]?\d*\b", "<TS>", text)
    text = re.sub(r"\b\d{2}:\d{2}:\d{2}[.,]?\d*\b", "<TS>", text)
    text = re.sub(r"0x[0-9a-fA-F]+", "<HEX>", text)
    text = re.sub(r"/[^ \n\t]+", "<PATH>", text)
    text = re.sub(r"\b\d+\b", "<INT>", text)
    return text


def safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
    return cleaned or "run"


def extract_tracebacks(lines: list[str]) -> list[str]:
    snippets = []
    i = 0
    while i < len(lines):
        if lines[i].startswith(TRACEBACK_START):
            block = [lines[i]]
            i += 1
            while i < len(lines):
                block.append(lines[i])
                if re.match(r"^[A-Za-z_][\w\.]*:", lines[i]):
                    i += 1
                    break
                i += 1
            snippets.append("\n".join(block))
            continue
        i += 1
    return snippets


def extract_log_errors(lines: list[str]) -> list[str]:
    snippets = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if TS_RE.match(line) and (" ERROR " in line or " CRITICAL " in line):
            block = [line]
            i += 1
            while i < len(lines) and not TS_RE.match(lines[i]):
                block.append(lines[i])
                i += 1
            snippets.append("\n".join(block))
            continue
        i += 1
    return snippets


def extract_failure_markers(lines: list[str]) -> list[str]:
    snippets = []
    for line in lines:
        for marker in FAILURE_MARKERS:
            if marker.lower() in line.lower():
                snippets.append(line)
                break
    return snippets


def extract_snippets(text: str) -> list[str]:
    lines = text.splitlines()
    snippets: list[str] = []
    snippets.extend(extract_tracebacks(lines))
    snippets.extend(extract_log_errors(lines))
    snippets.extend(extract_failure_markers(lines))
    return snippets


def categorize_snippet(snippet: str) -> str:
    if snippet.startswith(TRACEBACK_START):
        return "traceback"
    lines = snippet.splitlines()
    if lines:
        first = lines[0]
        if TS_RE.match(first) and (" ERROR " in first or " CRITICAL " in first):
            return "log_error"
    for marker in FAILURE_MARKERS:
        if marker.lower() in snippet.lower():
            return "failure_marker"
    return "other"


def infer_pdb_id(text: str) -> str | None:
    match = PDB_RE.search(text)
    if match:
        return match.group(0).upper()
    return None


def path_has_skip_ext(path: Path) -> bool:
    suffixes = [s.lower() for s in path.suffixes]
    if not suffixes:
        return False
    if suffixes[-1] in SKIP_EXTS:
        return True
    if len(suffixes) >= 2 and "".join(suffixes[-2:]) in SKIP_EXTS:
        return True
    return False


def should_include_failed(
    path: Path, run_id: str, recent_hours: float, now: float, tail_kb: int
) -> tuple[bool, str]:
    if path.suffix.lower() not in LOG_EXTS:
        return False, ""
    if not is_recent(path, recent_hours, now):
        return False, ""
    name_match = run_id in path.name
    try:
        tail = tail_text(path, tail_kb)
    except OSError:
        return False, ""
    if name_match:
        return True, tail
    if f"run_id={run_id}" in tail or f"docked/{run_id}" in tail or run_id in tail:
        return True, tail
    return False, ""


def list_run_protein_logs(run_dir: Path) -> Iterable[Path]:
    if not run_dir.is_dir():
        return []
    candidates = []
    for entry in run_dir.iterdir():
        if entry.is_dir():
            protein_log = entry / "protein.log"
            if protein_log.is_file():
                candidates.append(protein_log)
    return candidates


def list_top_level_logs(run_dir: Path) -> Iterable[Path]:
    if not run_dir.is_dir():
        return []
    common = {
        "main.log",
        "run.log",
        "pipeline.log",
        "workflow.log",
        "error.log",
    }
    logs = []
    for entry in run_dir.iterdir():
        if not entry.is_file():
            continue
        if entry.name in common or entry.suffix.lower() == ".log":
            logs.append(entry)
    return logs


def list_deep_logs(root: Path, run_id: str) -> list[FileMeta]:
    deep_logs: list[FileMeta] = []
    docked_root = root / "docked" / run_id
    candidate_pdbs = {p.upper() for p in PDB_RE.findall(run_id)}
    if docked_root.is_dir():
        for entry in docked_root.iterdir():
            if entry.is_dir() and len(entry.name) == 4 and entry.name.isalnum():
                candidate_pdbs.add(entry.name.upper())
    deepcoy_root = root / "DeepCoy_duds" / "deepcoy_work"
    for pdb_id in sorted(candidate_pdbs):
        pdb_root = deepcoy_root / pdb_id
        if not pdb_root.is_dir():
            continue
        for dirpath, _, filenames in os.walk(pdb_root):
            for filename in filenames:
                lower = filename.lower()
                if not (
                    lower.startswith(DEEP_LOG_PREFIXES)
                    or Path(filename).suffix.lower() in LOG_EXTS
                ):
                    continue
                path = Path(dirpath) / filename
                if path_has_skip_ext(path):
                    continue
                deep_logs.append(FileMeta(path=path, source="deepcoy", pdb_id=pdb_id))
    return deep_logs


def list_deep_run_logs(run_dir: Path, source: str) -> list[FileMeta]:
    entries: list[FileMeta] = []
    if not run_dir.is_dir():
        return entries
    for dirpath, _, filenames in os.walk(run_dir):
        for filename in filenames:
            path = Path(dirpath) / filename
            if path_has_skip_ext(path):
                continue
            if path.suffix.lower() not in LOG_EXTS:
                continue
            pdb_id = infer_pdb_id(path.as_posix())
            entries.append(FileMeta(path=path, source=source, pdb_id=pdb_id))
    return entries


def collect_files(args: argparse.Namespace) -> list[FileMeta]:
    root = Path(args.root).resolve()
    run_id = args.run_id
    now = time.time()

    files: dict[Path, FileMeta] = {}

    docked_root = root / "docked" / run_id
    post_docked_root = root / "post_docked" / run_id
    failed_root = root / "failed"

    for path in list_run_protein_logs(docked_root):
        if is_recent(path, args.recent_hours, now):
            pdb_id = infer_pdb_id(path.parent.name)
            files[path] = FileMeta(path=path, source="docked", pdb_id=pdb_id)

    for path in list_top_level_logs(docked_root):
        if is_recent(path, args.recent_hours, now):
            files[path] = FileMeta(path=path, source="docked")

    for path in list_top_level_logs(post_docked_root):
        if is_recent(path, args.recent_hours, now):
            files[path] = FileMeta(path=path, source="post_docked")

    if failed_root.is_dir():
        for entry in failed_root.iterdir():
            if not entry.is_file():
                continue
            include, tail = should_include_failed(
                entry, run_id, args.recent_hours, now, args.tail_kb
            )
            if not include:
                continue
            pdb_id = infer_pdb_id(entry.name) or infer_pdb_id(tail)
            files[entry] = FileMeta(path=entry, source="failed", pdb_id=pdb_id)

    if args.deep:
        for meta in list_deep_run_logs(docked_root, "docked"):
            if is_recent(meta.path, args.recent_hours, now):
                files.setdefault(meta.path, meta)
        for meta in list_deep_run_logs(post_docked_root, "post_docked"):
            if is_recent(meta.path, args.recent_hours, now):
                files.setdefault(meta.path, meta)
        for meta in list_deep_logs(root, run_id):
            if is_recent(meta.path, args.recent_hours, now):
                files.setdefault(meta.path, meta)

    return list(files.values())


def summarize_files(
    files: list[FileMeta], tail_kb: int, max_examples: int
) -> tuple[dict[str, dict[str, object]], int, int]:
    signatures: dict[str, dict[str, object]] = {}
    raw_hits = 0
    scanned_files = 0

    for meta in files:
        try:
            text = tail_text(meta.path, tail_kb)
        except OSError:
            continue
        scanned_files += 1
        snippets = extract_snippets(text)
        if not snippets:
            continue
        for snippet in snippets:
            raw_hits += 1
            category = categorize_snippet(snippet)
            signature = normalize_signature(snippet)
            entry = signatures.setdefault(signature, {"count": 0, "examples": []})
            entry["count"] = int(entry["count"]) + 1
            if len(entry["examples"]) < max_examples:
                entry["examples"].append(
                    {
                        "path": str(meta.path.resolve()),
                        "source": meta.source,
                        "pdb_id": meta.pdb_id,
                        "snippet": snippet,
                        "category": category,
                    }
                )
    return signatures, scanned_files, raw_hits


def format_example(example: dict[str, object]) -> str:
    parts = [
        f"path={example['path']}",
        f"source={example['source']}",
    ]
    if example.get("pdb_id"):
        parts.append(f"pdb={example['pdb_id']}")
    header = " ".join(parts)
    snippet = example["snippet"]
    return f"{header}\n{snippet}"


def format_report(
    run_id: str,
    signatures: dict[str, dict[str, object]],
    scanned_files: int,
    raw_hits: int,
    top: int,
) -> str:
    sorted_sigs = sorted(
        signatures.items(),
        key=lambda item: int(item[1]["count"]),
        reverse=True,
    )
    lines = [
        f"run_id: {run_id}",
        f"scanned_files: {scanned_files}",
        f"raw_hits: {raw_hits}",
        f"unique_signatures: {len(signatures)}",
    ]
    if not sorted_sigs:
        return "\n".join(lines) + "\n"
    lines.append("")
    lines.append("Signatures:")
    for idx, (signature, payload) in enumerate(sorted_sigs[:top], start=1):
        lines.append("")
        lines.append(f"[{idx}] count={payload['count']}")
        lines.append(signature)
        for example in payload["examples"]:
            lines.append(format_example(example))
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    files = collect_files(args)
    signatures, scanned_files, raw_hits = summarize_files(
        files, args.tail_kb, args.max_examples
    )

    if args.json_out:
        json_payload = {
            "run_id": args.run_id,
            "scanned_files": scanned_files,
            "raw_hits": raw_hits,
            "unique_signatures": len(signatures),
            "signatures": [
                {
                    "signature": signature,
                    "count": payload["count"],
                    "examples": payload["examples"],
                }
                for signature, payload in sorted(
                    signatures.items(),
                    key=lambda item: int(item[1]["count"]),
                    reverse=True,
                )
            ],
        }
        Path(args.json_out).write_text(
            json.dumps(json_payload, indent=2), encoding="utf-8"
        )

    report_text = format_report(
        args.run_id, signatures, scanned_files, raw_hits, args.top
    )
    summary_root = Path(args.root).resolve() / LOG_SUMMARY_DIRNAME
    summary_root.mkdir(parents=True, exist_ok=True)
    summary_path = summary_root / f"{safe_filename(args.run_id)}_summary.log"
    summary_path.write_text(report_text, encoding="utf-8")
    print(report_text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
