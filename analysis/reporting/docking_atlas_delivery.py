"""Provider-neutral download projection and offline deployment preflight."""

from __future__ import annotations

import hashlib
import json
import mimetypes
import re
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.parse import quote, urlsplit, urlunsplit

from analysis.reporting.docking_atlas_edge_assets import (
    APP_CSS,
    APP_JS,
    INDEX_HTML,
    PUBLIC_HEADERS,
    WORKER_SOURCE,
    WRANGLER_TEMPLATE,
)

DELIVERY_SCHEMA_VERSION = 1
DOWNLOAD_MANIFEST_NAME = "download_projection.json"
PREFLIGHT_REPORT_NAME = "deployment_preflight.json"
MIB = 1024 * 1024
_CLOUDFLARE_ASSET_BYTES = 25 * MIB
_CLOUDFLARE_ASSET_COUNTS = {"free": 20_000, "paid": 100_000}
_R2_BUCKET = re.compile(r"^[a-z0-9](?:[a-z0-9-]{1,61}[a-z0-9])$")
_WORKER = re.compile(r"^[A-Za-z0-9](?:[A-Za-z0-9-]{0,61}[A-Za-z0-9])?$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_mapping(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{label} is missing or symlinked: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{label} is not valid JSON: {path}") from exc
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must contain a JSON object: {path}")
    return dict(value)


def _public_base_url(value: str) -> str:
    parsed = urlsplit(value.strip())
    if (
        parsed.scheme != "https"
        or not parsed.netloc
        or parsed.username is not None
        or parsed.password is not None
        or parsed.query
        or parsed.fragment
    ):
        raise ValueError(
            "download base URL must be an HTTPS origin/path without credentials, "
            "query parameters, or a fragment"
        )
    segments = [segment for segment in parsed.path.split("/") if segment]
    if any(segment in {".", ".."} for segment in segments):
        raise ValueError("download base URL path must not contain dot segments")
    normalized_path = "/" + "/".join(
        quote(segment, safe="-._~") for segment in segments
    )
    normalized_path = normalized_path.rstrip("/")
    return urlunsplit(("https", parsed.netloc.lower(), normalized_path, "", ""))


def _declared_downloads(site_dir: Path) -> list[dict[str, Any]]:
    downloads_root = site_dir / "downloads"
    if downloads_root.is_symlink() or not downloads_root.is_dir():
        raise ValueError(
            f"downloads directory is missing or symlinked: {downloads_root}"
        )
    manifest = _read_mapping(site_dir / "site_manifest.json", "site manifest")
    value = manifest.get("downloads", [])
    if not isinstance(value, list):
        raise ValueError("site manifest downloads must be a JSON array")
    rows: list[dict[str, Any]] = []
    for index, item in enumerate(value):
        if not isinstance(item, Mapping):
            raise ValueError(f"site manifest downloads[{index}] must be an object")
        rows.append(dict(item))
    return rows


def build_download_projection(
    site_dir: Path,
    *,
    release_id: str,
    release_token: str,
    download_base_url: str,
) -> dict[str, Any]:
    """Map verified local downloads to immutable, provider-neutral public URLs."""
    site_dir = Path(site_dir).resolve()
    base_url = _public_base_url(download_base_url)
    seen_paths: set[str] = set()
    seen_urls: set[str] = set()
    uploads: list[dict[str, Any]] = []
    public_entries: list[dict[str, Any]] = []
    for index, row in enumerate(_declared_downloads(site_dir)):
        local_url = str(row.get("url") or row.get("href") or "").strip()
        parsed = urlsplit(local_url)
        if (
            parsed.scheme
            or parsed.netloc
            or parsed.query
            or parsed.fragment
            or not parsed.path.startswith("downloads/")
        ):
            raise ValueError(
                f"declared download {index} must be a local downloads/ path"
            )
        local_path = parsed.path
        candidate = site_dir / local_path
        try:
            resolved = candidate.resolve(strict=True)
            relative = resolved.relative_to((site_dir / "downloads").resolve())
        except (OSError, ValueError) as exc:
            raise ValueError(
                f"declared download escapes or is missing from downloads/: {local_path}"
            ) from exc
        if candidate.is_symlink() or not resolved.is_file():
            raise ValueError(f"declared download is missing or symlinked: {local_path}")
        relative_name = relative.as_posix()
        if relative_name in seen_paths:
            raise ValueError(f"duplicate declared download path: {local_path}")
        seen_paths.add(relative_name)
        digest = _sha256(resolved)
        declared_hash = str(row.get("content_hash") or "")
        if declared_hash != f"sha256:{digest}":
            raise ValueError(
                f"declared download hash does not match file bytes: {local_path}"
            )
        object_key = f"releases/{release_token}/downloads/{relative_name}"
        encoded_key = "/".join(
            quote(part, safe="-._~") for part in object_key.split("/")
        )
        public_url = f"{base_url}/{encoded_key}"
        if public_url in seen_urls:
            raise ValueError(f"duplicate projected public URL: {public_url}")
        seen_urls.add(public_url)
        media_type = (
            mimetypes.guess_type(resolved.name)[0] or "application/octet-stream"
        )
        upload = {
            "source_path": local_path,
            "object_key": object_key,
            "public_url": public_url,
            "sha256": digest,
            "size_bytes": resolved.stat().st_size,
            "media_type": media_type,
            "cache_control": "public, max-age=31536000, immutable",
        }
        uploads.append(upload)
        public_entries.append(
            {
                "label": str(row.get("label") or resolved.name),
                "url": public_url,
                "content_hash": f"sha256:{digest}",
                "size_bytes": upload["size_bytes"],
                "media_type": media_type,
            }
        )
    return {
        "schema_version": DELIVERY_SCHEMA_VERSION,
        "release_id": release_id,
        "release_token": release_token,
        "provider_contract": "https_object_origin",
        "download_base_url": base_url,
        "immutability": "release-qualified keys; never replace bytes at an existing key",
        "upload_count": len(uploads),
        "upload_bytes": sum(int(row["size_bytes"]) for row in uploads),
        "uploads": uploads,
        "public_entries": public_entries,
    }


def _edge_objects(edge_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    manifest = _read_mapping(edge_dir / "object_manifest.json", "edge object manifest")
    value = manifest.get("objects")
    if not isinstance(value, list):
        raise ValueError("edge object manifest objects must be a JSON array")
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, item in enumerate(value):
        if not isinstance(item, Mapping):
            raise ValueError(f"edge object manifest objects[{index}] must be an object")
        row = dict(item)
        key = str(row.get("key") or "")
        path = edge_dir / "objects" / key
        try:
            resolved = path.resolve(strict=True)
            resolved.relative_to((edge_dir / "objects").resolve())
        except (OSError, ValueError) as exc:
            raise ValueError(
                f"edge object is missing or escapes object root: {key}"
            ) from exc
        if not key or key in seen or path.is_symlink() or not resolved.is_file():
            raise ValueError(
                f"edge object key is duplicate, missing, or symlinked: {key}"
            )
        seen.add(key)
        if row.get("size_bytes") != resolved.stat().st_size or row.get(
            "sha256"
        ) != _sha256(resolved):
            raise ValueError(f"edge object integrity mismatch: {key}")
        rows.append(row)
    tree_paths = list((edge_dir / "objects").rglob("*"))
    symlinks = [path for path in tree_paths if path.is_symlink()]
    if symlinks:
        raise ValueError(f"symlinked edge object is forbidden: {symlinks[0]}")
    actual = {
        path.relative_to(edge_dir / "objects").as_posix()
        for path in tree_paths
        if path.is_file()
    }
    if actual != seen:
        raise ValueError("edge object tree does not exactly match object_manifest.json")
    return manifest, rows


def _static_assets(edge_dir: Path) -> list[dict[str, Any]]:
    public = edge_dir / "public"
    if public.is_symlink() or not public.is_dir():
        raise ValueError(f"edge public directory is missing or symlinked: {public}")
    rows: list[dict[str, Any]] = []
    for path in sorted(public.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"symlinked edge static asset is forbidden: {path}")
        if path.is_file():
            rows.append(
                {
                    "path": path.relative_to(public).as_posix(),
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
            )
    return rows


def _validate_generated_shell(
    edge_dir: Path,
    *,
    marker: Mapping[str, Any],
    object_manifest: Mapping[str, Any],
) -> None:
    release_token = str(object_manifest.get("release_token") or "")
    if (
        marker.get("schema_version") != DELIVERY_SCHEMA_VERSION
        or marker.get("release_token") != release_token
    ):
        raise ValueError("edge bundle marker does not match the object manifest")
    expected = {
        "public/index.html": INDEX_HTML.encode("utf-8"),
        "public/_headers": PUBLIC_HEADERS.encode("utf-8"),
        "public/assets/app.css": APP_CSS.encode("utf-8"),
        "public/assets/app.js": APP_JS.encode("utf-8"),
        "src/index.mjs": WORKER_SOURCE.encode("utf-8"),
        "wrangler.toml": WRANGLER_TEMPLATE.encode("utf-8"),
    }
    for relative, content in expected.items():
        path = edge_dir / relative
        if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
            raise ValueError(
                f"generated edge support file integrity mismatch: {relative}"
            )
    config = _read_mapping(
        edge_dir / "public" / "release-config.json", "edge release configuration"
    )
    if config.get("release_token") != release_token or config.get(
        "release_id"
    ) != object_manifest.get("release_id"):
        raise ValueError("edge release configuration does not match object manifest")


def _safe_report_path(edge_dir: Path, report_path: Path | None) -> Path:
    destination = (
        Path(report_path).resolve() if report_path else edge_dir / PREFLIGHT_REPORT_NAME
    )
    if destination.is_symlink():
        raise ValueError(
            f"refusing symlinked deployment preflight report: {destination}"
        )
    return destination


def _cloudflare_name_errors(
    *,
    worker_name: str | None,
    production_bucket: str | None,
    preview_bucket: str | None,
) -> list[dict[str, str]]:
    blockers: list[dict[str, str]] = []
    if not worker_name:
        blockers.append(
            {"code": "cloudflare_worker_name_required", "field": "worker_name"}
        )
    elif not _WORKER.fullmatch(worker_name):
        blockers.append(
            {"code": "cloudflare_worker_name_invalid", "field": "worker_name"}
        )
    for field, value in (
        ("production_bucket", production_bucket),
        ("preview_bucket", preview_bucket),
    ):
        if not value:
            blockers.append({"code": f"cloudflare_{field}_required", "field": field})
        elif not _R2_BUCKET.fullmatch(value):
            blockers.append({"code": f"cloudflare_{field}_invalid", "field": field})
    if production_bucket and production_bucket == preview_bucket:
        blockers.append(
            {"code": "cloudflare_buckets_must_be_distinct", "field": "preview_bucket"}
        )
    return blockers


def _wrangler_config(
    worker_name: str, production_bucket: str, preview_bucket: str
) -> str:
    return f'''name = "{worker_name}"
main = "src/index.mjs"
compatibility_date = "2026-07-13"

[assets]
directory = "./public"
binding = "ASSETS"
not_found_handling = "single-page-application"
run_worker_first = ["/api/*"]

[[r2_buckets]]
binding = "ATLAS_RELEASES"
bucket_name = "{production_bucket}"
preview_bucket_name = "{preview_bucket}"
'''


def build_deployment_preflight(
    edge_dir: Path,
    *,
    site_dir: Path | None = None,
    provider: str = "generic",
    cloudflare_plan: str = "free",
    worker_name: str | None = None,
    production_bucket: str | None = None,
    preview_bucket: str | None = None,
    report_path: Path | None = None,
) -> dict[str, Any]:
    """Verify a delivery bundle offline and list explicit remaining resources."""
    edge_dir = Path(edge_dir).resolve()
    marker = _read_mapping(edge_dir / ".atlas-edge-bundle.json", "edge bundle marker")
    if marker.get("generator") != "atlas publish edge-bundle":
        raise ValueError("edge bundle marker has an unexpected generator")
    object_manifest, objects = _edge_objects(edge_dir)
    _validate_generated_shell(edge_dir, marker=marker, object_manifest=object_manifest)
    static_assets = _static_assets(edge_dir)
    projection_path = edge_dir / DOWNLOAD_MANIFEST_NAME
    projection = (
        _read_mapping(projection_path, "download projection")
        if projection_path.is_file() and not projection_path.is_symlink()
        else None
    )
    blockers: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    if projection is None:
        blockers.append(
            {
                "code": "download_projection_required",
                "action": "rebuild edge bundle with --download-base-url https://<download-host>",
            }
        )
    elif site_dir is None:
        blockers.append(
            {
                "code": "download_source_site_required",
                "action": "pass --site-dir for offline download-byte verification",
            }
        )
    else:
        resolved_site = Path(site_dir).resolve()
        downloads_root = resolved_site / "downloads"
        if downloads_root.is_symlink() or not downloads_root.is_dir():
            raise ValueError(
                f"download source directory is missing or symlinked: {downloads_root}"
            )
        uploads = projection.get("uploads")
        if not isinstance(uploads, list):
            raise ValueError("download projection uploads must be a JSON array")
        if projection.get("release_token") != object_manifest.get("release_token"):
            raise ValueError(
                "download projection release token does not match edge objects"
            )
        if projection.get("upload_count") != len(uploads):
            raise ValueError("download projection upload count is inconsistent")
        public_entries = projection.get("public_entries")
        if not isinstance(public_entries, list) or len(public_entries) != len(uploads):
            raise ValueError("download projection public entries are inconsistent")
        base_url = _public_base_url(str(projection.get("download_base_url") or ""))
        expected_bytes = 0
        seen_keys: set[str] = set()
        for index, item in enumerate(uploads):
            if not isinstance(item, Mapping):
                raise ValueError(
                    f"download projection uploads[{index}] must be an object"
                )
            source_path = str(item.get("source_path") or "")
            object_key = str(item.get("object_key") or "")
            expected_prefix = (
                f"releases/{object_manifest.get('release_token')}/downloads/"
            )
            if not object_key.startswith(expected_prefix) or object_key in seen_keys:
                raise ValueError(
                    f"invalid or duplicate projected object key: {object_key}"
                )
            seen_keys.add(object_key)
            encoded_key = "/".join(
                quote(part, safe="-._~") for part in object_key.split("/")
            )
            expected_public_url = f"{base_url}/{encoded_key}"
            public_entry = public_entries[index]
            if not isinstance(public_entry, Mapping):
                raise ValueError(
                    f"download projection public_entries[{index}] must be an object"
                )
            if (
                item.get("public_url") != expected_public_url
                or public_entry.get("url") != expected_public_url
                or public_entry.get("content_hash") != f"sha256:{item.get('sha256')}"
                or public_entry.get("size_bytes") != item.get("size_bytes")
            ):
                raise ValueError(
                    f"projected download URL/metadata is inconsistent: {object_key}"
                )
            candidate = resolved_site / source_path
            try:
                resolved = candidate.resolve(strict=True)
                resolved.relative_to(downloads_root.resolve())
            except (OSError, ValueError) as exc:
                raise ValueError(
                    "projected download source is missing or escapes downloads/: "
                    f"{source_path}"
                ) from exc
            if candidate.is_symlink() or not resolved.is_file():
                raise ValueError(
                    f"projected download source is missing or symlinked: {source_path}"
                )
            if (
                item.get("sha256") != _sha256(resolved)
                or item.get("size_bytes") != resolved.stat().st_size
            ):
                raise ValueError(
                    f"projected download integrity mismatch: {source_path}"
                )
            expected_bytes += resolved.stat().st_size
        if projection.get("upload_bytes") != expected_bytes:
            raise ValueError("download projection byte count is inconsistent")
    if provider not in {"generic", "cloudflare-workers-r2"}:
        raise ValueError(f"unsupported deployment provider: {provider}")
    provider_details: dict[str, Any] = {"name": provider}
    required_resources: list[dict[str, Any]] = []
    credential_environment_names: list[str] = []
    generated_config: str | None = None
    if provider == "cloudflare-workers-r2":
        if cloudflare_plan not in _CLOUDFLARE_ASSET_COUNTS:
            raise ValueError("Cloudflare plan must be 'free' or 'paid'")
        blockers.extend(
            _cloudflare_name_errors(
                worker_name=worker_name,
                production_bucket=production_bucket,
                preview_bucket=preview_bucket,
            )
        )
        asset_limit = _CLOUDFLARE_ASSET_COUNTS[cloudflare_plan]
        too_large = [
            row
            for row in static_assets
            if int(row["size_bytes"]) > _CLOUDFLARE_ASSET_BYTES
        ]
        if len(static_assets) > asset_limit:
            blockers.append(
                {
                    "code": "cloudflare_static_asset_count_exceeded",
                    "actual": len(static_assets),
                    "limit": asset_limit,
                }
            )
        if too_large:
            blockers.append(
                {
                    "code": "cloudflare_static_asset_size_exceeded",
                    "limit_bytes": _CLOUDFLARE_ASSET_BYTES,
                    "paths": [str(row["path"]) for row in too_large],
                }
            )
        if projection is not None and urlsplit(
            str(projection.get("download_base_url") or "")
        ).path not in {"", "/"}:
            blockers.append(
                {
                    "code": "cloudflare_r2_download_origin_must_not_have_path",
                    "field": "download_base_url",
                    "action": (
                        "use an origin-only custom domain URL for direct R2 public "
                        "bucket delivery"
                    ),
                }
            )
        warnings.append(
            {
                "code": "cloudflare_remote_state_not_verified",
                "action": (
                    "after explicit approval, verify bucket/domain ownership and "
                    "post-upload hashes remotely"
                ),
            }
        )
        required_resources = [
            {
                "kind": "worker_service",
                "name": worker_name,
                "creation": "created only by a future explicitly approved Wrangler deploy",
            },
            {
                "kind": "r2_bucket",
                "environment": "production",
                "name": production_bucket,
                "creation": "create in Cloudflare dashboard before upload",
            },
            {
                "kind": "r2_bucket",
                "environment": "preview",
                "name": preview_bucket,
                "creation": "create in Cloudflare dashboard before remote preview",
            },
            {
                "kind": "https_download_origin",
                "name": projection.get("download_base_url") if projection else None,
                "creation": "attach a production custom domain to the download bucket",
            },
        ]
        credential_environment_names = [
            "CLOUDFLARE_ACCOUNT_ID",
            "CLOUDFLARE_API_TOKEN",
        ]
        if worker_name and production_bucket and preview_bucket and not blockers:
            generated_config = _wrangler_config(
                worker_name, production_bucket, preview_bucket
            )
        provider_details.update(
            {
                "plan": cloudflare_plan,
                "static_asset_count_limit": asset_limit,
                "static_asset_size_limit_bytes": _CLOUDFLARE_ASSET_BYTES,
                "configuration": {
                    "worker_name": worker_name,
                    "production_bucket": production_bucket,
                    "preview_bucket": preview_bucket,
                    "binding": "ATLAS_RELEASES",
                },
                "official_limit_reference": "https://developers.cloudflare.com/workers/platform/limits/#static-assets",
                "official_auth_reference": "https://developers.cloudflare.com/workers/wrangler/system-environment-variables/",
            }
        )
    elif projection is not None:
        warnings.append(
            {
                "code": "generic_provider_limits_not_enforced",
                "action": "compare the recorded asset/object counts with the selected host",
            }
        )

    report = {
        "schema_version": DELIVERY_SCHEMA_VERSION,
        "status": "ready" if not blockers else "blocked",
        "provider": provider_details,
        "release_id": object_manifest.get("release_id"),
        "release_token": object_manifest.get("release_token"),
        "static_assets": {
            "count": len(static_assets),
            "bytes": sum(int(row["size_bytes"]) for row in static_assets),
            "largest_bytes": max(
                (int(row["size_bytes"]) for row in static_assets), default=0
            ),
        },
        "edge_objects": {
            "count": len(objects),
            "bytes": sum(int(row["size_bytes"]) for row in objects),
        },
        "downloads": {
            "projected": projection is not None,
            "count": int(projection.get("upload_count", 0)) if projection else 0,
            "bytes": int(projection.get("upload_bytes", 0)) if projection else 0,
            "manifest": DOWNLOAD_MANIFEST_NAME if projection else None,
        },
        "required_resources": required_resources,
        "credential_environment_names": credential_environment_names,
        "secret_values_inspected": False,
        "network_checks_performed": False,
        "upload_performed": False,
        "deployment_performed": False,
        "blockers": blockers,
        "warnings": warnings,
    }
    destination = _safe_report_path(edge_dir, report_path)
    config_path = edge_dir / "wrangler.preflight.toml"
    if generated_config is not None and config_path.is_symlink():
        raise ValueError(f"refusing symlinked Wrangler preflight config: {config_path}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if generated_config is not None:
        config_path.write_text(generated_config, encoding="utf-8")
    return report


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="atlas publish preflight",
        description=(
            "Validate an Atlas delivery bundle and report required resources; "
            "never upload or deploy."
        ),
    )
    parser.add_argument("--edge-dir", required=True, type=Path)
    parser.add_argument("--site-dir", type=Path)
    parser.add_argument(
        "--provider",
        choices=("generic", "cloudflare-workers-r2"),
        default="generic",
    )
    parser.add_argument("--cloudflare-plan", choices=("free", "paid"), default="free")
    parser.add_argument("--worker-name")
    parser.add_argument("--production-bucket")
    parser.add_argument("--preview-bucket")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    try:
        report = build_deployment_preflight(
            args.edge_dir,
            site_dir=args.site_dir,
            provider=args.provider,
            cloudflare_plan=args.cloudflare_plan,
            worker_name=args.worker_name,
            production_bucket=args.production_bucket,
            preview_bucket=args.preview_bucket,
            report_path=args.report,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "ready" else 2


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "DOWNLOAD_MANIFEST_NAME",
    "PREFLIGHT_REPORT_NAME",
    "build_deployment_preflight",
    "build_download_projection",
    "main",
]
