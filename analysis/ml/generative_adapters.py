from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
from collections.abc import Iterator, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SMILES_COLUMNS = (
    "smiles",
    "SMILES",
    "generated_smiles",
    "canonical_smiles",
    "mol",
    "molecule",
)
SCORE_COLUMNS = (
    "score",
    "total_score",
    "reward",
    "likelihood",
    "prior_likelihood",
    "agent_likelihood",
)
ID_COLUMNS = ("molecule_id", "compound_id", "sample_id", "id", "name", "identifier")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def atlas_tools_dir() -> Path:
    default_tools = repo_root().parent.parent / "tools"
    return Path(os.environ.get("ATLAS_TOOLS_DIR") or default_tools).expanduser()


def default_generative_root() -> Path:
    return Path(os.environ.get("ATLAS_GENERATIVE_ROOT") or (atlas_tools_dir() / "ml" / "generative")).expanduser()


def default_reinvent_root() -> Path:
    return Path(os.environ.get("ATLAS_REINVENT4_ROOT") or (default_generative_root() / "REINVENT4")).expanduser()


def default_guacamol_root() -> Path:
    return Path(os.environ.get("ATLAS_GUACAMOL_ROOT") or (default_generative_root() / "guacamol")).expanduser()


def default_moses_root() -> Path:
    return Path(os.environ.get("ATLAS_MOSES_ROOT") or (default_generative_root() / "moses")).expanduser()


def default_generative_python() -> Path | None:
    explicit = os.environ.get("ATLAS_GENERATIVE_PYTHON")
    if explicit:
        return Path(explicit).expanduser()
    candidate = atlas_tools_dir() / "envs" / "atlas-ml-generative" / "bin" / "python"
    return candidate if candidate.exists() else None


def get_git_sha(root: Path | None = None) -> str:
    try:
        proc = subprocess.run(
            ["git", "-C", str(root or repo_root()), "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
        )
    except Exception:
        return "unknown"
    value = proc.stdout.strip()
    return value if proc.returncode == 0 and value else "unknown"


def environment_snapshot() -> dict[str, Any]:
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: str | Path | None) -> dict[str, Any]:
    if path is None:
        return {"path": None, "exists": False}
    resolved = Path(path).expanduser()
    if not resolved.exists():
        return {"path": str(resolved), "exists": False}
    if resolved.is_dir():
        return {
            "path": str(resolved),
            "exists": True,
            "type": "directory",
        }
    return {
        "path": str(resolved),
        "exists": True,
        "type": "file",
        "size_bytes": int(resolved.stat().st_size),
        "sha256": sha256_file(resolved),
    }


def write_json(path: str | Path, payload: dict[str, Any]) -> Path:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return out_path


def normalize_command_remainder(command: Sequence[str] | None) -> list[str]:
    if not command:
        return []
    values = list(command)
    if values and values[0] == "--":
        values = values[1:]
    return values


def run_external_command(
    command: Sequence[str],
    *,
    out_dir: str | Path,
    stem: str,
    cwd: str | Path | None = None,
    timeout_sec: int | None = None,
) -> dict[str, Any]:
    cmd = list(command)
    if not cmd:
        raise ValueError("external command is empty")

    run_dir = Path(out_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = run_dir / f"{stem}.stdout.log"
    stderr_path = run_dir / f"{stem}.stderr.log"
    started = utc_now()
    started_wall = datetime.now(timezone.utc)
    status = "completed"
    returncode = 1
    error = ""

    try:
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
            proc = subprocess.run(
                cmd,
                cwd=str(cwd) if cwd else None,
                stdout=stdout,
                stderr=stderr,
                text=True,
                timeout=timeout_sec if timeout_sec and timeout_sec > 0 else None,
                check=False,
            )
        returncode = int(proc.returncode)
        if returncode != 0:
            status = "failed"
    except FileNotFoundError as exc:
        status = "missing_executable"
        returncode = 127
        error = str(exc)
        stderr_path.write_text(error + "\n", encoding="utf-8")
        stdout_path.touch()
    except subprocess.TimeoutExpired as exc:
        status = "timeout"
        returncode = 124
        error = str(exc)
    ended = utc_now()
    duration = (datetime.now(timezone.utc) - started_wall).total_seconds()
    return {
        "status": status,
        "command": cmd,
        "cwd": str(Path(cwd).expanduser()) if cwd else None,
        "returncode": returncode,
        "started_at": started,
        "ended_at": ended,
        "duration_seconds": round(duration, 3),
        "stdout": file_record(stdout_path),
        "stderr": file_record(stderr_path),
        "error": error,
    }


def iter_smiles_records(
    path: str | Path,
    *,
    smiles_column: str | None = None,
    score_column: str | None = None,
) -> Iterator[dict[str, Any]]:
    source_path = Path(path)
    suffix = source_path.suffix.lower()
    if suffix == ".jsonl":
        yield from _iter_jsonl(source_path, smiles_column=smiles_column, score_column=score_column)
        return
    if suffix == ".json":
        yield from _iter_json(source_path, smiles_column=smiles_column, score_column=score_column)
        return
    if suffix in {".csv", ".tsv"} or _looks_delimited(source_path):
        yield from _iter_delimited(source_path, smiles_column=smiles_column, score_column=score_column)
        return
    yield from _iter_text_smiles(source_path)


def standardize_smiles_file(
    input_path: str | Path,
    out_path: str | Path,
    *,
    source_name: str,
    smiles_column: str | None = None,
    score_column: str | None = None,
) -> dict[str, Any]:
    source = Path(input_path)
    if not source.exists():
        raise FileNotFoundError(f"SMILES source does not exist: {source}")

    output = Path(out_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["molecule_id", "smiles", "source", "rank", "score", "raw_record_json"]
    rows = 0
    seen: set[str] = set()
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for rank, record in enumerate(
            iter_smiles_records(source, smiles_column=smiles_column, score_column=score_column),
            start=1,
        ):
            smiles = _clean(record.get("smiles"))
            if not smiles:
                continue
            molecule_id = _clean(record.get("molecule_id")) or f"generated_{rank:06d}"
            seen.add(smiles)
            rows += 1
            writer.writerow(
                {
                    "molecule_id": molecule_id,
                    "smiles": smiles,
                    "source": source_name,
                    "rank": rank,
                    "score": _clean(record.get("score")),
                    "raw_record_json": json.dumps(record.get("raw", record), sort_keys=True, default=str),
                }
            )

    return {
        "status": "captured",
        "source": file_record(source),
        "standardized_smiles": file_record(output),
        "rows": rows,
        "unique_smiles": len(seen),
        "source_name": source_name,
        "smiles_column": smiles_column,
        "score_column": score_column,
    }


def write_benchmark_smiles(
    input_path: str | Path,
    out_path: str | Path,
    *,
    smiles_column: str | None = None,
) -> dict[str, Any]:
    source = Path(input_path)
    if not source.exists():
        raise FileNotFoundError(f"generated SMILES file does not exist: {source}")
    output = Path(out_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    seen: set[str] = set()
    with output.open("w", encoding="utf-8") as handle:
        for record in iter_smiles_records(source, smiles_column=smiles_column):
            smiles = _clean(record.get("smiles"))
            if not smiles:
                continue
            rows += 1
            seen.add(smiles)
            handle.write(smiles + "\n")
    return {
        "status": "prepared",
        "source": file_record(source),
        "benchmark_smiles": file_record(output),
        "rows": rows,
        "unique_smiles": len(seen),
        "smiles_column": smiles_column,
    }


def collect_result_files(
    paths: Sequence[str | Path],
    *,
    copy_to: str | Path | None = None,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    copy_root = Path(copy_to) if copy_to else None
    if copy_root:
        copy_root.mkdir(parents=True, exist_ok=True)
    for raw_path in paths:
        source = Path(raw_path).expanduser()
        record = file_record(source)
        if copy_root and source.exists() and source.is_file():
            target = _unique_copy_path(copy_root / source.name)
            shutil.copy2(source, target)
            record["copied_to"] = file_record(target)
        records.append(record)
    return records


def extract_metrics(paths: Sequence[str | Path], *, suite: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for raw_path in paths:
        path = Path(raw_path)
        if not path.exists() or not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix == ".json":
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                continue
            rows.extend(_metrics_from_json(payload, suite=suite, source_file=path, prefix=""))
        elif suffix in {".csv", ".tsv"}:
            rows.extend(_metrics_from_delimited(path, suite=suite))
    return rows


def write_metrics_csv(path: str | Path, metrics: Sequence[dict[str, Any]]) -> Path:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["suite", "source_file", "metric", "value"]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in metrics:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    return out_path


def _clean(value: Any) -> str:
    if value is None:
        return ""
    return " ".join(str(value).strip().split())


def _looks_delimited(path: Path) -> bool:
    if not path.exists() or not path.is_file():
        return False
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            first = handle.readline()
    except OSError:
        return False
    lowered = first.lower()
    return "smiles" in lowered and ("," in first or "\t" in first)


def _iter_text_smiles(path: Path) -> Iterator[dict[str, Any]]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            yield {
                "molecule_id": parts[1] if len(parts) > 1 else "",
                "smiles": parts[0],
                "score": "",
                "raw": {"line": line_no, "text": stripped},
            }


def _iter_delimited(
    path: Path,
    *,
    smiles_column: str | None = None,
    score_column: str | None = None,
) -> Iterator[dict[str, Any]]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else _sniff_delimiter(path)
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            yield _record_from_mapping(row, smiles_column=smiles_column, score_column=score_column)


def _iter_jsonl(
    path: Path,
    *,
    smiles_column: str | None = None,
    score_column: str | None = None,
) -> Iterator[dict[str, Any]]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            try:
                payload = json.loads(text)
            except json.JSONDecodeError:
                payload = {"smiles": text}
            if isinstance(payload, dict):
                record = _record_from_mapping(payload, smiles_column=smiles_column, score_column=score_column)
                record["raw"] = {"line": line_no, "payload": payload}
                yield record
            elif isinstance(payload, str):
                yield {"molecule_id": "", "smiles": payload, "score": "", "raw": {"line": line_no, "payload": payload}}


def _iter_json(
    path: Path,
    *,
    smiles_column: str | None = None,
    score_column: str | None = None,
) -> Iterator[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    for item in _json_items(payload):
        if isinstance(item, dict):
            yield _record_from_mapping(item, smiles_column=smiles_column, score_column=score_column)
        elif isinstance(item, str):
            yield {"molecule_id": "", "smiles": item, "score": "", "raw": item}


def _json_items(payload: Any) -> Iterator[Any]:
    if isinstance(payload, list):
        yield from payload
        return
    if isinstance(payload, dict):
        for key in ("molecules", "smiles", "generated", "results", "samples", "data"):
            value = payload.get(key)
            if isinstance(value, list):
                yield from value
                return
        yield payload


def _record_from_mapping(
    row: dict[str, Any],
    *,
    smiles_column: str | None = None,
    score_column: str | None = None,
) -> dict[str, Any]:
    smiles = _first_value(row, [smiles_column] if smiles_column else SMILES_COLUMNS)
    score = _first_value(row, [score_column] if score_column else SCORE_COLUMNS)
    molecule_id = _first_value(row, ID_COLUMNS)
    return {
        "molecule_id": molecule_id,
        "smiles": smiles,
        "score": score,
        "raw": row,
    }


def _first_value(row: dict[str, Any], columns: Sequence[str | None]) -> str:
    by_lower = {str(key).lower(): key for key in row.keys()}
    for column in columns:
        if not column:
            continue
        key = column if column in row else by_lower.get(str(column).lower())
        if key is None:
            continue
        value = _clean(row.get(key))
        if value:
            return value
    return ""


def _sniff_delimiter(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            sample = handle.read(4096)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return str(dialect.delimiter)
    except Exception:
        return ","


def _unique_copy_path(path: Path) -> Path:
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    for idx in range(1, 10_000):
        candidate = path.with_name(f"{stem}_{idx}{suffix}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"could not allocate unique copy path under {path.parent}")


def _metrics_from_json(payload: Any, *, suite: str, source_file: Path, prefix: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if isinstance(payload, dict):
        label = _clean(payload.get("benchmark_name") or payload.get("name") or payload.get("metric"))
        if label and _is_number(payload.get("score")):
            rows.append(_metric_row(suite, source_file, label, payload["score"]))
        if label and _is_number(payload.get("value")):
            rows.append(_metric_row(suite, source_file, label, payload["value"]))
        for key, value in payload.items():
            next_prefix = f"{prefix}.{key}" if prefix else str(key)
            if _is_number(value):
                rows.append(_metric_row(suite, source_file, next_prefix, value))
            elif isinstance(value, (dict, list)):
                rows.extend(_metrics_from_json(value, suite=suite, source_file=source_file, prefix=next_prefix))
    elif isinstance(payload, list):
        for idx, item in enumerate(payload):
            rows.extend(_metrics_from_json(item, suite=suite, source_file=source_file, prefix=f"{prefix}[{idx}]"))
    return rows


def _metrics_from_delimited(path: Path, *, suite: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    delimiter = "\t" if path.suffix.lower() == ".tsv" else _sniff_delimiter(path)
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for idx, row in enumerate(reader):
            if idx >= 2000:
                break
            metric = _first_value(row, ("metric", "benchmark", "benchmark_name", "name"))
            value = _first_value(row, ("value", "score", "mean"))
            if metric and _is_number(value):
                rows.append(_metric_row(suite, path, metric, value))
                continue
            for key, raw_value in row.items():
                if _is_number(raw_value):
                    label = f"{metric}.{key}" if metric else str(key)
                    rows.append(_metric_row(suite, path, label, raw_value))
    return rows


def _is_number(value: Any) -> bool:
    if isinstance(value, bool) or value is None:
        return False
    try:
        float(value)
    except (TypeError, ValueError):
        return False
    return True


def _metric_row(suite: str, source_file: Path, metric: str, value: Any) -> dict[str, Any]:
    return {
        "suite": suite,
        "source_file": str(source_file),
        "metric": metric,
        "value": float(value),
    }
