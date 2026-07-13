"""Deterministic deployment assets and offline integrity checks for Atlas sites."""

from __future__ import annotations

import hashlib
import json
import mimetypes
import sqlite3
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Iterable, Mapping
from urllib.parse import unquote, urlsplit

BUNDLE_SCHEMA_VERSION = 1
INVENTORY_NAME = "release_inventory.json"
CHECKSUM_NAME = "release_checksums.sha256"
REPORT_NAME = "release_verification.json"
_EXCLUDED = {INVENTORY_NAME, CHECKSUM_NAME, REPORT_NAME}
_FORBIDDEN_PATHS = ("/stor/", "/home/", "/tmp/")
_TEXT_SUFFIXES = {".css", ".csv", ".html", ".js", ".json", ".sha256", ".txt"}
_PRIVATE_NAMES = (".private.", "build_summary.private")
_OUTSIDE_BUNDLE = Path("__outside_bundle__")

PAGES_HEADERS = """/*
  Cache-Control: public, max-age=0, must-revalidate
  Content-Security-Policy: default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'none'
  Permissions-Policy: camera=(), geolocation=(), microphone=(), payment=(), usb=()
  Referrer-Policy: no-referrer
  X-Content-Type-Options: nosniff
  X-Frame-Options: DENY

/assets/*
  Cache-Control: public, max-age=3600, must-revalidate

/downloads/*
  Cache-Control: public, max-age=31536000, immutable
"""


class _LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        wanted = (
            "href"
            if tag in {"a", "link"}
            else "src"
            if tag in {"img", "script"}
            else ""
        )
        if not wanted:
            return
        for name, value in attrs:
            if name == wanted and value:
                self.links.append(value)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _bundle_files(site_dir: Path) -> list[Path]:
    return sorted(
        (
            path
            for path in site_dir.rglob("*")
            if path.is_file()
            and not path.is_symlink()
            and path.relative_to(site_dir).as_posix() not in _EXCLUDED
        ),
        key=lambda path: path.relative_to(site_dir).as_posix(),
    )


def _entry(site_dir: Path, path: Path) -> dict[str, Any]:
    relative = path.relative_to(site_dir).as_posix()
    media_type = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
    return {
        "path": relative,
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
        "media_type": media_type,
    }


def write_release_inventory(site_dir: Path) -> dict[str, Any]:
    """Write deterministic inventory, checksums, and Cloudflare Pages headers."""
    site_dir = Path(site_dir).resolve()
    if not site_dir.is_dir():
        raise FileNotFoundError(f"site directory does not exist: {site_dir}")
    (site_dir / "_headers").write_text(PAGES_HEADERS, encoding="utf-8")
    manifest_path = site_dir / "site_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    entries = [_entry(site_dir, path) for path in _bundle_files(site_dir)]
    inventory = {
        "bundle_schema_version": BUNDLE_SCHEMA_VERSION,
        "release_id": str(manifest.get("release_id") or ""),
        "file_count": len(entries),
        "total_size_bytes": sum(int(entry["size_bytes"]) for entry in entries),
        "files": entries,
    }
    (site_dir / INVENTORY_NAME).write_text(
        json.dumps(inventory, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    checksum_lines = [f"{entry['sha256']}  {entry['path']}" for entry in entries]
    (site_dir / CHECKSUM_NAME).write_text(
        "\n".join(checksum_lines) + "\n", encoding="utf-8"
    )
    return inventory


def _json_links(value: Any) -> Iterable[str]:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if str(key).lower() in {"url", "href", "src", "entrypoint"} and isinstance(
                item, str
            ):
                yield item
            else:
                yield from _json_links(item)
    elif isinstance(value, list):
        for item in value:
            yield from _json_links(item)


def _local_target(site_dir: Path, source: Path, value: str) -> Path | None:
    parsed = urlsplit(value)
    if (
        parsed.scheme
        or parsed.netloc
        or not parsed.path
        or parsed.path.startswith("data:")
    ):
        return None
    decoded = unquote(parsed.path)
    candidate = (
        site_dir / decoded.lstrip("/")
        if decoded.startswith("/")
        else source.parent / decoded
    )
    try:
        resolved = candidate.resolve()
        resolved.relative_to(site_dir)
    except (OSError, ValueError):
        return _OUTSIDE_BUNDLE
    return resolved


def _link_errors(site_dir: Path) -> list[dict[str, str]]:
    errors: list[dict[str, str]] = []
    for source in _bundle_files(site_dir):
        links: Iterable[str] = ()
        try:
            if source.suffix.lower() == ".html":
                parser = _LinkParser()
                parser.feed(source.read_text(encoding="utf-8"))
                links = parser.links
            elif source.suffix.lower() == ".json":
                links = _json_links(json.loads(source.read_text(encoding="utf-8")))
        except (OSError, UnicodeDecodeError, ValueError):
            continue
        link_source = (
            site_dir / "index.html" if source.suffix.lower() == ".json" else source
        )
        for link in links:
            target = _local_target(site_dir, link_source, link)
            if target == _OUTSIDE_BUNDLE or (
                target is not None and not target.exists()
            ):
                errors.append(
                    {
                        "code": "missing_local_link",
                        "source": source.relative_to(site_dir).as_posix(),
                        "target": link,
                    }
                )
    return errors


def _sqlite_leaks(path: Path) -> list[str]:
    leaks: list[str] = []
    try:
        with sqlite3.connect(path) as connection:
            tables = [
                row[0]
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
                )
            ]
            for table in tables:
                columns = [
                    row[1]
                    for row in connection.execute(f"PRAGMA table_info({table})")
                    if "TEXT" in str(row[2]).upper()
                ]
                for column in columns:
                    for (value,) in connection.execute(
                        f"SELECT {column} FROM {table} WHERE {column} IS NOT NULL"
                    ):
                        if any(prefix in str(value) for prefix in _FORBIDDEN_PATHS):
                            leaks.append(f"{table}.{column}")
                            break
    except sqlite3.DatabaseError:
        leaks.append("invalid_sqlite_database")
    return sorted(set(leaks))


def _leak_errors(site_dir: Path) -> list[dict[str, str]]:
    errors: list[dict[str, str]] = []
    for path in site_dir.rglob("*"):
        relative = path.relative_to(site_dir).as_posix()
        if path.is_symlink():
            errors.append({"code": "symlink_not_allowed", "path": relative})
            continue
        if not path.is_file():
            continue
        if any(token in path.name for token in _PRIVATE_NAMES):
            errors.append({"code": "private_file_in_bundle", "path": relative})
        if path.suffix.lower() in {".sqlite", ".db"}:
            for location in _sqlite_leaks(path):
                errors.append(
                    {
                        "code": "private_path_in_sqlite",
                        "path": relative,
                        "location": location,
                    }
                )
        elif path.suffix.lower() in _TEXT_SUFFIXES or path.name == "_headers":
            try:
                text = path.read_text(encoding="utf-8")
            except (OSError, UnicodeDecodeError):
                continue
            if any(prefix in text for prefix in _FORBIDDEN_PATHS):
                errors.append({"code": "private_path_in_text", "path": relative})
    return errors


def _checksum_errors(
    site_dir: Path, actual: Mapping[str, Mapping[str, Any]]
) -> list[dict[str, str]]:
    path = site_dir / CHECKSUM_NAME
    if not path.is_file():
        return []
    declared: dict[str, str] = {}
    try:
        for line in path.read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            if len(digest) != 64 or name in declared:
                raise ValueError
            declared[name] = digest
    except (OSError, ValueError):
        return [{"code": "checksum_file_invalid", "path": CHECKSUM_NAME}]
    errors: list[dict[str, str]] = []
    for name in sorted(set(declared) | set(actual)):
        if name not in declared:
            errors.append({"code": "file_missing_from_checksums", "path": name})
        elif name not in actual:
            errors.append({"code": "checksum_file_missing", "path": name})
        elif declared[name] != actual[name]["sha256"]:
            errors.append({"code": "declared_checksum_mismatch", "path": name})
    return errors


def verify_release_bundle(
    site_dir: Path, *, report_path: Path | None = None
) -> dict[str, Any]:
    """Verify inventory, links, and public-safety invariants without network access."""
    site_dir = Path(site_dir).resolve()
    canonical_report = site_dir / REPORT_NAME
    if canonical_report.is_symlink():
        raise ValueError(
            f"refusing symlinked canonical verification report: {canonical_report}"
        )
    destination = (
        Path(report_path).resolve() if report_path is not None else canonical_report
    )
    if destination != canonical_report and destination.is_relative_to(site_dir):
        raise ValueError(
            "custom verification reports must be written outside the deployable "
            f"site directory: {destination}"
        )
    errors: list[dict[str, Any]] = []
    required = {
        "index.html",
        "site_manifest.json",
        "_headers",
        INVENTORY_NAME,
        CHECKSUM_NAME,
    }
    for name in sorted(required):
        if not (site_dir / name).is_file():
            errors.append({"code": "required_file_missing", "path": name})
    inventory_path = site_dir / INVENTORY_NAME
    expected: dict[str, dict[str, Any]] = {}
    if inventory_path.is_file():
        try:
            inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
            rows = inventory["files"]
            if not isinstance(rows, list):
                raise TypeError
            for row in rows:
                name = str(row["path"])
                if name in expected:
                    raise ValueError
                expected[name] = dict(row)
            manifest = json.loads(
                (site_dir / "site_manifest.json").read_text(encoding="utf-8")
            )
            if (
                inventory.get("bundle_schema_version") != BUNDLE_SCHEMA_VERSION
                or inventory.get("release_id") != str(manifest.get("release_id") or "")
                or inventory.get("file_count") != len(rows)
                or inventory.get("total_size_bytes")
                != sum(int(row["size_bytes"]) for row in rows)
            ):
                raise ValueError
        except (OSError, ValueError, KeyError, TypeError):
            errors.append({"code": "inventory_invalid", "path": INVENTORY_NAME})
    actual = {
        path.relative_to(site_dir).as_posix(): _entry(site_dir, path)
        for path in _bundle_files(site_dir)
    }
    for name in sorted(set(expected) | set(actual)):
        if name not in expected:
            errors.append({"code": "file_missing_from_inventory", "path": name})
        elif name not in actual:
            errors.append({"code": "inventory_file_missing", "path": name})
        elif expected[name] != actual[name]:
            errors.append({"code": "checksum_or_size_mismatch", "path": name})
    errors.extend(_checksum_errors(site_dir, actual))
    errors.extend(_link_errors(site_dir))
    errors.extend(_leak_errors(site_dir))
    report = {
        "verification_schema_version": BUNDLE_SCHEMA_VERSION,
        "site_directory": site_dir.name,
        "status": "passed" if not errors else "failed",
        "file_count": len(actual),
        "error_count": len(errors),
        "errors": errors,
    }
    destination.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def prepare_release_bundle(site_dir: Path) -> dict[str, Any]:
    """Write deployment assets and return the canonical verification report."""
    write_release_inventory(site_dir)
    return verify_release_bundle(site_dir)
