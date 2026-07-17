"""Create the ready FDA parent delta without rerunning ligand preparation."""

from __future__ import annotations

import argparse
import csv
import errno
import hashlib
import json
import os
import shutil
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence

from prep_ligands.library_index import LibraryIndex


SCHEMA_VERSION = 1
SOURCE_MANIFEST = "fda_named_library_manifest.csv"
SOURCE_QUARANTINE = "fda_named_library_quarantine.csv"
SOURCE_SUMMARY = "fda_named_library_summary.json"
OUTPUT_MANIFEST = "fda_delta_manifest.csv"
OUTPUT_QUARANTINE = "fda_delta_quarantine.csv"
OUTPUT_SUMMARY = "fda_delta_summary.json"
EXTRA_FIELDS = (
    "delta_subset_schema_version",
    "source_named_library_dir",
    "source_named_manifest_sha256",
    "source_output_pdbqt_path",
    "source_output_pdbqt_sha256",
    "source_sidecar_path",
    "source_sidecar_sha256",
    "pdbqt_materialization_method",
    "sidecar_materialization_method",
)


def materialize_fda_delta_library(
    *, named_library_dir: Path, output_dir: Path, prefer_hardlinks: bool = True
) -> dict[str, Any]:
    """Materialize exactly the 586 ``prepare_parent`` + ``ready`` entries."""

    source = named_library_dir.expanduser().resolve()
    output = output_dir.expanduser().resolve()
    _validate_roots(source, output)
    source_hashes, delta, fields = _load_delta_source(source)

    output.mkdir(parents=True, exist_ok=True)
    _validate_resume(output, source, source_hashes)
    materialized = _materialize_delta_rows(
        delta,
        source=source,
        output=output,
        source_manifest_hash=source_hashes["manifest"],
        prefer_hardlinks=prefer_hardlinks,
    )
    ready, quarantined = _partition_materialized(materialized)
    _validate_inventory(output, ready)

    output_manifest, output_quarantine = _write_delta_manifests(
        output,
        materialized=materialized,
        quarantined=quarantined,
        source_delta=delta,
        source_fields=fields,
    )
    index_path, relative_ready = _ensure_delta_index(output, ready)
    result = _summary(
        source,
        output,
        source_hashes,
        materialized,
        output_manifest,
        output_quarantine,
        index_path,
    )
    _write_json_once(output / OUTPUT_SUMMARY, result)
    return result


def _load_delta_source(
    source: Path,
) -> tuple[dict[str, str], list[dict[str, str]], list[str]]:
    manifest = source / SOURCE_MANIFEST
    quarantine = source / SOURCE_QUARANTINE
    summary_path = source / SOURCE_SUMMARY
    _require_files(manifest, quarantine, summary_path)
    source_hashes = {
        "manifest": _sha256(manifest),
        "quarantine": _sha256(quarantine),
        "summary": _sha256(summary_path),
    }
    rows, fields = _read_csv(manifest)
    quarantine_rows, _ = _read_csv(quarantine)
    source_summary = _read_json(summary_path)
    delta = _validate_source(
        source, rows, quarantine_rows, source_summary, source_hashes
    )
    return source_hashes, delta, fields


def _materialize_delta_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    source: Path,
    output: Path,
    source_manifest_hash: str,
    prefer_hardlinks: bool,
) -> list[dict[str, str]]:
    return [
        _materialize_row(
            row,
            source=source,
            output=output,
            source_manifest_hash=source_manifest_hash,
            prefer_hardlinks=prefer_hardlinks,
        )
        for row in rows
    ]


def _partition_materialized(
    materialized: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    materialized.sort(key=lambda row: row["parent_inchikey"])
    ready = [row for row in materialized if row["materialization_status"] == "ready"]
    quarantined = [
        row
        for row in materialized
        if row["materialization_status"] == "quarantined_non_dockable"
    ]
    if (len(materialized), len(ready), len(quarantined)) != (599, 586, 13):
        raise RuntimeError("FDA delta output partition drifted")
    return ready, quarantined


def _write_delta_manifests(
    output: Path,
    *,
    materialized: Sequence[Mapping[str, str]],
    quarantined: Sequence[Mapping[str, str]],
    source_delta: Sequence[Mapping[str, str]],
    source_fields: Sequence[str],
) -> tuple[Path, Path]:
    output_fields = [
        *source_fields,
        *[field for field in EXTRA_FIELDS if field not in source_fields],
    ]
    output_manifest = output / OUTPUT_MANIFEST
    output_quarantine = output / OUTPUT_QUARANTINE
    _write_csv_once(
        output_manifest,
        materialized,
        output_fields,
        adopt_rows=source_delta,
        adopt_fields=source_fields,
    )
    source_delta_quarantine = [
        row
        for row in source_delta
        if row["materialization_status"] == "quarantined_non_dockable"
    ]
    _write_csv_once(
        output_quarantine,
        quarantined,
        output_fields,
        adopt_rows=source_delta_quarantine,
        adopt_fields=source_fields,
    )
    return output_manifest, output_quarantine


def _ensure_delta_index(
    output: Path, ready: Sequence[Mapping[str, str]]
) -> tuple[Path, list[Path]]:
    index_path = output / "_manifest.json"
    relative_ready = [Path(row["filename"]) for row in ready]
    if not index_path.exists():
        LibraryIndex().write_manifest_for_root(output, relative_ready)
    _validate_index(index_path, relative_ready, output)
    return index_path, relative_ready


def _validate_source(
    source: Path,
    rows: Sequence[Mapping[str, str]],
    quarantine_rows: Sequence[Mapping[str, str]],
    summary: Mapping[str, Any],
    source_hashes: Mapping[str, str],
) -> list[dict[str, str]]:
    _validate_source_partition(rows)
    _validate_source_summary(summary, source_hashes)
    delta = [
        dict(row) for row in rows if row["materialization_action"] == "prepare_parent"
    ]
    _validate_delta_partition(delta)
    _validate_delta_quarantine(delta, quarantine_rows)
    _validate_delta_filenames(delta)
    for row in delta:
        _validate_source_row(source, row)
    return delta


def _validate_source_partition(rows: Sequence[Mapping[str, str]]) -> None:
    actions = dict(Counter(row.get("materialization_action", "") for row in rows))
    statuses = dict(Counter(row.get("materialization_status", "") for row in rows))
    expected_actions = {"legacy_byte_copy": 1264, "prepare_parent": 599}
    expected_statuses = {"ready": 1850, "quarantined_non_dockable": 13}
    if (
        len(rows) != 1863
        or actions != expected_actions
        or statuses != expected_statuses
    ):
        raise ValueError(
            f"finalized FDA manifest partition drifted: rows={len(rows)} "
            f"actions={actions} statuses={statuses}"
        )


def _validate_source_summary(
    summary: Mapping[str, Any], source_hashes: Mapping[str, str]
) -> None:
    expected_values = {
        "manifest_sha256": source_hashes["manifest"],
        "quarantine_sha256": source_hashes["quarantine"],
        "action_counts": {"legacy_byte_copy": 1264, "prepare_parent": 599},
        "status_counts": {"ready": 1850, "quarantined_non_dockable": 13},
    }
    for key, expected_value in expected_values.items():
        if summary.get(key) != expected_value:
            raise ValueError("finalized FDA summary disagrees with source artifacts")
    expected_counts = {"canonical_parent_rows": 1863, "ready_pdbqt_count": 1850}
    for key, expected_count in expected_counts.items():
        if int(summary.get(key, -1)) != expected_count:
            raise ValueError("finalized FDA summary disagrees with source artifacts")


def _validate_delta_partition(delta: Sequence[Mapping[str, str]]) -> None:
    delta_statuses = dict(Counter(row["materialization_status"] for row in delta))
    if delta_statuses != {"ready": 586, "quarantined_non_dockable": 13}:
        raise ValueError(f"FDA delta status partition drifted: {delta_statuses}")


def _validate_delta_quarantine(
    delta: Sequence[Mapping[str, str]],
    quarantine_rows: Sequence[Mapping[str, str]],
) -> None:
    quarantined_keys = {
        row["parent_inchikey"]
        for row in delta
        if row["materialization_status"] == "quarantined_non_dockable"
    }
    source_quarantine_keys = {row.get("parent_inchikey", "") for row in quarantine_rows}
    if len(quarantine_rows) != 13 or source_quarantine_keys != quarantined_keys:
        raise ValueError("finalized FDA quarantine does not match delta quarantine")


def _validate_delta_filenames(delta: Sequence[Mapping[str, str]]) -> None:
    names = [row["filename"] for row in delta]
    if len({name.casefold() for name in names}) != 599:
        raise ValueError("FDA delta filenames are not case-insensitively unique")


def _validate_source_row(source: Path, row: Mapping[str, str]) -> None:
    filename = row["filename"]
    _validate_delta_filename(filename)
    if row["materialization_status"] != "ready":
        _validate_quarantined_source_row(row, filename)
        return
    pdbqt, sidecar = _source_artifact_paths(row)
    _validate_source_artifact_scope(source, filename, pdbqt, sidecar)
    _validate_file(pdbqt, row["output_pdbqt_sha256"], row["output_pdbqt_bytes"])
    _validate_file(sidecar, row["sidecar_sha256"], None)


def _validate_delta_filename(filename: str) -> None:
    if Path(filename).name != filename or not filename.casefold().endswith(".pdbqt"):
        raise ValueError(f"unsafe FDA delta filename: {filename!r}")


def _validate_quarantined_source_row(row: Mapping[str, str], filename: str) -> None:
    if row.get("output_pdbqt_sha256") or row.get("sidecar_sha256"):
        raise ValueError(f"quarantined FDA row has output bytes: {filename}")


def _source_artifact_paths(row: Mapping[str, str]) -> tuple[Path, Path]:
    return (
        Path(row["output_pdbqt_path"]).expanduser().resolve(),
        Path(row["sidecar_path"]).expanduser().resolve(),
    )


def _validate_source_artifact_scope(
    source: Path, filename: str, pdbqt: Path, sidecar: Path
) -> None:
    if (
        pdbqt.parent != source
        or sidecar.parent != source
        or pdbqt.name != filename
        or sidecar.name != Path(filename).with_suffix(".ligprep_source.json").name
        or pdbqt.is_symlink()
        or sidecar.is_symlink()
    ):
        raise ValueError(f"FDA delta source artifact scope mismatch: {filename}")


def _materialize_row(
    row: Mapping[str, str],
    *,
    source: Path,
    output: Path,
    source_manifest_hash: str,
    prefer_hardlinks: bool,
) -> dict[str, str]:
    result = dict(row)
    filename = row["filename"]
    source_pdbqt = Path(row["output_pdbqt_path"]).resolve()
    source_sidecar = Path(row["sidecar_path"]).resolve()
    output_pdbqt = output / filename
    output_sidecar = output_pdbqt.with_suffix(".ligprep_source.json")
    ready = row["materialization_status"] == "ready"
    result.update(
        {
            "delta_subset_schema_version": str(SCHEMA_VERSION),
            "source_named_library_dir": str(source),
            "source_named_manifest_sha256": source_manifest_hash,
            "source_output_pdbqt_path": str(source_pdbqt) if ready else "",
            "source_output_pdbqt_sha256": row.get("output_pdbqt_sha256", ""),
            "source_sidecar_path": str(source_sidecar) if ready else "",
            "source_sidecar_sha256": row.get("sidecar_sha256", ""),
            "output_pdbqt_path": str(output_pdbqt),
            "sidecar_path": str(output_sidecar),
            "pdbqt_materialization_method": "",
            "sidecar_materialization_method": "",
        }
    )
    if ready:
        result["pdbqt_materialization_method"] = _link_or_copy(
            source_pdbqt, output_pdbqt, row["output_pdbqt_sha256"], prefer_hardlinks
        )
        result["sidecar_materialization_method"] = _link_or_copy(
            source_sidecar, output_sidecar, row["sidecar_sha256"], prefer_hardlinks
        )
    return result


def _link_or_copy(source: Path, target: Path, digest: str, hardlink: bool) -> str:
    existing_method = _existing_materialization_method(source, target, digest)
    if existing_method is not None:
        return existing_method
    if hardlink:
        hardlink_method = _try_hardlink(source, target, digest)
        if hardlink_method is not None:
            return hardlink_method
    return _copy_artifact(source, target, digest)


def _path_exists(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def _existing_materialization_method(
    source: Path, target: Path, digest: str
) -> str | None:
    if not _path_exists(target):
        return None
    if target.is_symlink() or not target.is_file() or _sha256(target) != digest:
        raise ValueError(f"refusing to overwrite changed FDA delta artifact: {target}")
    return _storage_method(source, target)


def _try_hardlink(source: Path, target: Path, digest: str) -> str | None:
    try:
        os.link(source, target)
    except OSError as exc:
        concurrent_method = _concurrent_materialization_method(
            source, target, digest, exc
        )
        if concurrent_method is not None:
            return concurrent_method
        if exc.errno not in {
            errno.EXDEV,
            errno.EPERM,
            errno.EACCES,
            errno.EOPNOTSUPP,
            errno.ENOTSUP,
        }:
            raise
        return None
    if _sha256(target) != digest:
        raise RuntimeError(f"FDA delta hardlink checksum mismatch: {target}")
    return "hardlink"


def _concurrent_materialization_method(
    source: Path, target: Path, digest: str, cause: OSError
) -> str | None:
    if not _path_exists(target):
        return None
    if target.is_symlink() or not target.is_file() or _sha256(target) != digest:
        raise ValueError(
            f"concurrent FDA delta artifact disagrees with source: {target}"
        ) from cause
    return _storage_method(source, target)


def _copy_artifact(source: Path, target: Path, digest: str) -> str:
    part = target.with_name(target.name + ".part")
    if part.exists():
        raise ValueError(f"stale FDA delta partial blocks resume: {part}")
    try:
        with source.open("rb") as src, part.open("xb") as dst:
            shutil.copyfileobj(src, dst, 1024 * 1024)
            dst.flush()
            os.fsync(dst.fileno())
        if _sha256(part) != digest:
            raise RuntimeError(f"FDA delta copy checksum mismatch: {part}")
        os.replace(part, target)
    except Exception:
        _remove_partial(part)
        raise
    return "copy"


def _remove_partial(path: Path) -> None:
    if path.exists():
        path.unlink()


def _storage_method(source: Path, target: Path) -> str:
    left, right = source.stat(), target.stat()
    return (
        "hardlink"
        if (left.st_dev, left.st_ino) == (right.st_dev, right.st_ino)
        else "copy"
    )


def _validate_inventory(output: Path, rows: Sequence[Mapping[str, str]]) -> None:
    expected_pdbqt, expected_sidecar = _expected_inventory(rows)
    _validate_inventory_names(output, expected_pdbqt, expected_sidecar)
    for row in rows:
        _validate_inventory_row(output, row)


def _expected_inventory(
    rows: Sequence[Mapping[str, str]],
) -> tuple[set[str], set[str]]:
    expected_pdbqt = {row["filename"] for row in rows}
    expected_sidecar = {
        Path(name).with_suffix(".ligprep_source.json").name for name in expected_pdbqt
    }
    return expected_pdbqt, expected_sidecar


def _validate_inventory_names(
    output: Path, expected_pdbqt: set[str], expected_sidecar: set[str]
) -> None:
    actual_pdbqt = _regular_inventory_names(output, "*.pdbqt")
    actual_sidecar = _regular_inventory_names(output, "*.ligprep_source.json")
    if actual_pdbqt != expected_pdbqt or actual_sidecar != expected_sidecar:
        raise RuntimeError(
            "FDA delta output has missing, extra, or symlinked artifacts"
        )


def _regular_inventory_names(output: Path, pattern: str) -> set[str]:
    return {
        path.name
        for path in output.glob(pattern)
        if path.is_file() and not path.is_symlink()
    }


def _validate_inventory_row(output: Path, row: Mapping[str, str]) -> None:
    output_pdbqt = output / row["filename"]
    _validate_file(
        output_pdbqt,
        row["output_pdbqt_sha256"],
        row["output_pdbqt_bytes"],
    )
    _validate_file(
        output_pdbqt.with_suffix(".ligprep_source.json"),
        row["sidecar_sha256"],
        None,
    )


def _validate_index(path: Path, entries: Sequence[Path], root: Path) -> None:
    payload = _read_json(path)
    expected = {entry.name.casefold(): entry.as_posix() for entry in entries}
    filenames = payload.get("filenames")
    raw_entries = payload.get("entries")
    if (
        not isinstance(filenames, dict)
        or filenames != expected
        or not isinstance(raw_entries, dict)
    ):
        raise ValueError("FDA delta LibraryIndex does not exactly match ready files")
    if not LibraryIndex._manifest_references_valid(
        root, entries=raw_entries, filenames=filenames
    ):
        raise ValueError("FDA delta LibraryIndex contains an invalid reference")


def _validate_resume(output: Path, source: Path, hashes: Mapping[str, str]) -> None:
    path = output / OUTPUT_SUMMARY
    if path.exists():
        old = _read_json(path)
        if old.get("source_named_library_dir") != str(source) or old.get(
            "source_sha256"
        ) != dict(hashes):
            raise ValueError(
                "FDA delta resume source changed; choose a new output directory"
            )


def _summary(
    source: Path,
    output: Path,
    hashes: Mapping[str, str],
    rows: Sequence[Mapping[str, str]],
    manifest: Path,
    quarantine: Path,
    index: Path,
) -> dict[str, Any]:
    ready = [row for row in rows if row["materialization_status"] == "ready"]
    inventory = hashlib.sha256()
    for path in sorted(output.glob("*.pdbqt"), key=lambda item: item.name.casefold()):
        inventory.update(path.name.encode())
        inventory.update(b"\0")
        inventory.update(path.read_bytes())
        inventory.update(b"\0")
    return {
        "schema_version": SCHEMA_VERSION,
        "source_named_library_dir": str(source),
        "source_sha256": dict(hashes),
        "output_dir": str(output),
        "selection": {
            "materialization_action": "prepare_parent",
            "materialization_status": "ready",
        },
        "delta_manifest_rows": len(rows),
        "ready_pdbqt_count": len(ready),
        "quarantine_count": len(rows) - len(ready),
        "status_counts": dict(
            sorted(Counter(row["materialization_status"] for row in rows).items())
        ),
        "pdbqt_materialization_methods": dict(
            sorted(
                Counter(row["pdbqt_materialization_method"] for row in ready).items()
            )
        ),
        "manifest_csv": str(manifest),
        "manifest_sha256": _sha256(manifest),
        "quarantine_csv": str(quarantine),
        "quarantine_sha256": _sha256(quarantine),
        "library_index": str(index),
        "library_index_sha256": _sha256(index),
        "library_inventory": {"count": len(ready), "sha256": inventory.hexdigest()},
        "policy": {
            "meeko_rerun": False,
            "source_bytes_mutated": False,
            "quarantine_rows_emit_no_pdbqt": True,
            "checksums_verified_before_and_after_materialization": True,
        },
    }


def _validate_roots(source: Path, output: Path) -> None:
    if not source.is_dir():
        raise ValueError(f"finalized FDA named-library directory not found: {source}")
    if output.exists() and not output.is_dir():
        raise ValueError(f"FDA delta output is not a directory: {output}")
    if source == output or source in output.parents or output in source.parents:
        raise ValueError("FDA delta output must be separate from the named library")


def _require_files(*paths: Path) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise ValueError(
            "missing finalized FDA named-library inputs: " + ", ".join(missing)
        )


def _validate_file(path: Path, digest: str, size: object) -> None:
    if len(digest) != 64 or not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"missing or invalid FDA delta artifact: {path}")
    if size not in (None, "") and path.stat().st_size != int(str(size)):
        raise ValueError(f"FDA delta artifact size changed: {path}")
    if _sha256(path) != digest.casefold():
        raise ValueError(f"FDA delta artifact checksum changed: {path}")


def _read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        return [dict(row) for row in reader], list(reader.fieldnames)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"invalid JSON: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON root is not an object: {path}")
    return value


def _write_csv_once(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
    *,
    adopt_rows: Sequence[Mapping[str, Any]] = (),
    adopt_fields: Sequence[str] = (),
) -> None:
    import io

    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(fields), extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    content = buffer.getvalue().encode()
    if path.exists() and path.read_bytes() != content:
        existing_rows, existing_fields = _read_csv(path)
        normalized_adopt = [
            {field: str(row.get(field, "")) for field in adopt_fields}
            for row in adopt_rows
        ]
        if existing_fields != list(adopt_fields) or existing_rows != normalized_adopt:
            raise ValueError(
                f"refusing to overwrite changed FDA delta metadata: {path}"
            )
        _replace_atomic(path, content)
        return
    _write_once(path, content)


def _write_json_once(path: Path, value: Mapping[str, Any]) -> None:
    _write_once(path, (json.dumps(value, indent=2, sort_keys=True) + "\n").encode())


def _write_once(path: Path, content: bytes) -> None:
    if path.exists():
        if not path.is_file() or path.read_bytes() != content:
            raise ValueError(
                f"refusing to overwrite changed FDA delta metadata: {path}"
            )
        return
    part = path.with_name(path.name + ".part")
    if part.exists():
        raise ValueError(f"stale FDA delta metadata partial blocks resume: {part}")
    with part.open("xb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(part, path)


def _replace_atomic(path: Path, content: bytes) -> None:
    part = path.with_name(path.name + ".part")
    if part.exists():
        raise ValueError(f"stale FDA delta metadata partial blocks adoption: {part}")
    with part.open("xb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(part, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--named-library-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--copy", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    result = materialize_fda_delta_library(
        named_library_dir=args.named_library_dir,
        output_dir=args.output_dir,
        prefer_hardlinks=not args.copy,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = ["main", "materialize_fda_delta_library"]
