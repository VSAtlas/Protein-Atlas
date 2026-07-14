"""Build an optional R2-backed SPA bundle from a public Atlas release site."""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from analysis.reporting.docking_atlas_edge_assets import (
    APP_CSS,
    APP_JS,
    INDEX_HTML,
    PUBLIC_HEADERS,
    WORKER_SOURCE,
    WRANGLER_TEMPLATE,
)
from analysis.reporting.docking_atlas_delivery import (
    DOWNLOAD_MANIFEST_NAME,
    build_download_projection,
)

_PRIVATE_PATH_MARKERS = (b"/stor/", b"/home/", b"/tmp/")
_BUNDLE_MARKER = ".atlas-edge-bundle.json"
_TOKEN_LIMIT = 512
MAX_BROWSER_PAYLOAD_BYTES = 128 * 1024 * 1024
MAX_BROWSER_PAIR_COUNT = 50_000


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("utf-8")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def _write_json(path: Path, value: Any) -> None:
    _write_bytes(path, _json_bytes(value))


def _required_mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a JSON object")
    return dict(value)


def _mapping_list(value: Any, label: str) -> list[dict[str, Any]]:
    if not isinstance(value, list):
        raise ValueError(f"{label} must be a JSON array")
    return [
        _required_mapping(item, f"{label}[{index}]") for index, item in enumerate(value)
    ]


def _identifier(row: Mapping[str, Any], names: Sequence[str], label: str) -> str:
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip():
            return str(value).strip()
    raise ValueError(f"{label} lacks a stable identifier ({', '.join(names)})")


def _route_token(value: str, label: str) -> str:
    token = base64.urlsafe_b64encode(value.encode("utf-8")).decode("ascii").rstrip("=")
    if not token or len(token) > _TOKEN_LIMIT:
        raise ValueError(
            f"{label} produces a route identifier outside 1-{_TOKEN_LIMIT} characters"
        )
    return token


def _release_token(release_id: str, source_sha256: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", release_id.casefold()).strip("-")
    slug = (slug or "release")[:96].rstrip("-")
    return f"{slug}-{source_sha256[:16]}"


def _display_label(row: Mapping[str, Any], identifier: str) -> str:
    for name in ("display_name", "ligand_display_name", "receptor_label", "name"):
        value = row.get(name)
        if value is not None and str(value).strip():
            return str(value).strip()
    return identifier


def _object_entry(objects_root: Path, path: Path, data: bytes) -> dict[str, Any]:
    return {
        "key": path.relative_to(objects_root).as_posix(),
        "size_bytes": len(data),
        "sha256": _sha256(data),
        "content_type": "application/json; charset=utf-8",
        "cache_control": "public, max-age=31536000, immutable",
    }


def _write_object(
    objects_root: Path,
    relative_key: str,
    value: Any,
    entries: list[dict[str, Any]],
) -> dict[str, Any]:
    path = objects_root / relative_key
    data = _json_bytes(value)
    _write_bytes(path, data)
    entry = _object_entry(objects_root, path, data)
    entries.append(entry)
    return entry


def _rank_value(pair: Mapping[str, Any], entity_kind: str) -> Any:
    if entity_kind == "targets":
        return pair.get("rank_within_receptor", pair.get("protein_rank"))
    return pair.get("rank_across_receptors", pair.get("drug_rank"))


def _pair_summary(
    pair: Mapping[str, Any],
    *,
    pair_id: str,
    pair_route_id: str,
    target_id: str,
    target_route_id: str,
    target_label: str,
    drug_id: str,
    drug_route_id: str,
    drug_label: str,
) -> dict[str, Any]:
    return {
        "id": pair_id,
        "pair_route_id": pair_route_id,
        "target_id": target_id,
        "target_route_id": target_route_id,
        "target_label": target_label,
        "drug_id": drug_id,
        "drug_route_id": drug_route_id,
        "drug_label": drug_label,
        "final_status": pair.get("final_status"),
        "failure_code": pair.get("failure_code"),
        "failure_reason": pair.get("failure_reason"),
        "final_score": pair.get("final_score"),
        "final_score_source": pair.get("final_score_source"),
        "rank_eligible": pair.get("rank_eligible"),
        "ranking_eligibility_reason": pair.get("ranking_eligibility_reason"),
        "rank_within_receptor": pair.get(
            "rank_within_receptor", pair.get("protein_rank")
        ),
        "rank_across_receptors": pair.get(
            "rank_across_receptors", pair.get("drug_rank")
        ),
        "pose_valid": pair.get("pose_valid"),
    }


def _counterpart_summary(
    summary: Mapping[str, Any], entity_kind: str
) -> dict[str, Any]:
    result = dict(summary)
    if entity_kind == "targets":
        result.update(
            {
                "counterpart_kind": "drugs",
                "counterpart_route_id": summary["drug_route_id"],
                "counterpart_label": summary["drug_label"],
            }
        )
    else:
        result.update(
            {
                "counterpart_kind": "targets",
                "counterpart_route_id": summary["target_route_id"],
                "counterpart_label": summary["target_label"],
            }
        )
    result["rank"] = _rank_value(summary, entity_kind)
    return result


def _pair_sort_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    rank = row.get("rank")
    return (
        rank is None,
        float(rank) if isinstance(rank, (int, float)) else 0.0,
        str(row.get("counterpart_label") or "").casefold(),
        str(row.get("pair_route_id") or ""),
    )


def _prepare_output(output_dir: Path, overwrite: bool) -> None:
    if output_dir.is_symlink():
        raise ValueError(f"refusing symlinked edge output directory: {output_dir}")
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        return
    marker = output_dir / _BUNDLE_MARKER
    if not overwrite:
        raise ValueError(f"edge bundle already exists: {output_dir}; use --overwrite")
    if marker.is_symlink() or not marker.is_file():
        raise ValueError(
            f"refusing to replace unrecognized directory without {_BUNDLE_MARKER}: "
            f"{output_dir}"
        )
    try:
        marker_payload = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(
            f"refusing to replace directory with invalid {_BUNDLE_MARKER}: {output_dir}"
        ) from exc
    if not isinstance(marker_payload, Mapping) or (
        marker_payload.get("schema_version") != 1
        or marker_payload.get("generator") != "atlas publish edge-bundle"
    ):
        raise ValueError(
            "refusing to replace directory without a valid Atlas edge bundle "
            f"marker: {output_dir}"
        )
    shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)


def _validate_locations(site_dir: Path, output_dir: Path) -> None:
    resolved_site = site_dir.resolve()
    resolved_output = output_dir.resolve()
    if (
        resolved_output == resolved_site
        or resolved_output.is_relative_to(resolved_site)
        or resolved_site.is_relative_to(resolved_output)
    ):
        raise ValueError(
            "edge output and verified static site must be disjoint directories"
        )


def _source_payload(site_dir: Path) -> tuple[dict[str, Any], bytes]:
    browser_path = site_dir / "downloads" / "release_browser.json"
    if browser_path.is_symlink() or not browser_path.is_file():
        raise ValueError(
            f"public browser payload is missing or symlinked: {browser_path}"
        )
    payload_size = browser_path.stat().st_size
    if payload_size > MAX_BROWSER_PAYLOAD_BYTES:
        raise ValueError(
            "bounded edge projection refuses release_browser.json larger than "
            f"{MAX_BROWSER_PAYLOAD_BYTES} bytes (found {payload_size}); this command "
            "loads browser JSON in memory and is not the publication-scale path. "
            "Use the --database streaming SQLite projection rather than generating "
            "per-pair static HTML."
        )
    raw = browser_path.read_bytes()
    for marker in _PRIVATE_PATH_MARKERS:
        if marker in raw:
            raise ValueError(
                f"public browser payload contains machine-local path marker: "
                f"{marker.decode('ascii')}"
            )
    try:
        payload = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid public browser payload: {browser_path}") from exc
    return _required_mapping(payload, "release_browser.json"), raw


def _entity_rows(
    payload: Mapping[str, Any],
    *,
    primary_name: str,
    fallback_name: str,
    identifier_names: Sequence[str],
) -> tuple[dict[str, dict[str, Any]], dict[str, str], dict[str, str]]:
    raw = payload.get(primary_name)
    if raw is None:
        raw = payload.get(fallback_name, [])
    rows = _mapping_list(raw, primary_name)
    by_id: dict[str, dict[str, Any]] = {}
    route_ids: dict[str, str] = {}
    labels: dict[str, str] = {}
    for index, row in enumerate(rows):
        identifier = _identifier(row, identifier_names, f"{primary_name}[{index}]")
        if identifier in by_id:
            raise ValueError(f"duplicate {primary_name} identifier: {identifier}")
        route_id = _route_token(identifier, f"{primary_name} identifier")
        entity = dict(row)
        entity["id"] = identifier
        entity["route_id"] = route_id
        by_id[identifier] = entity
        route_ids[identifier] = route_id
        labels[identifier] = _display_label(row, identifier)
    return by_id, route_ids, labels


def build_edge_bundle(
    site_dir: Path,
    output_dir: Path,
    *,
    overwrite: bool = False,
    download_base_url: str | None = None,
) -> dict[str, Any]:
    """Decompose a public browser payload into immutable, R2-ready objects."""
    _validate_locations(site_dir, output_dir)
    payload, source_bytes = _source_payload(site_dir)
    source_sha256 = _sha256(source_bytes)
    release = _required_mapping(payload.get("release"), "release")
    release_id = _identifier(release, ("id", "release_id", "version"), "release")
    release_token = _release_token(release_id, source_sha256)
    download_projection = (
        build_download_projection(
            site_dir,
            release_id=release_id,
            release_token=release_token,
            download_base_url=download_base_url,
        )
        if download_base_url
        else None
    )

    targets, target_routes, target_labels = _entity_rows(
        payload,
        primary_name="targets",
        fallback_name="proteins",
        identifier_names=("id", "target_id", "target_key", "protein_id"),
    )
    drugs, drug_routes, drug_labels = _entity_rows(
        payload,
        primary_name="ligands",
        fallback_name="drugs",
        identifier_names=("id", "drug_id", "canonical_id"),
    )
    pairs = _mapping_list(payload.get("pairs", []), "pairs")
    if len(pairs) > MAX_BROWSER_PAIR_COUNT:
        raise ValueError(
            "bounded edge projection refuses more than "
            f"{MAX_BROWSER_PAIR_COUNT} pair cells (found {len(pairs)}); this command "
            "is for conference/intermediate snapshots, not the full matrix. A "
            "streaming SQLite edge projection is required for publication scale."
        )

    artifacts = _mapping_list(payload.get("artifacts", []), "artifacts")
    artifacts_by_pair: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for artifact in artifacts:
        pair_reference = artifact.get("pair_cell_id")
        if pair_reference is not None:
            artifacts_by_pair[str(pair_reference)].append(artifact)
    validated_pair_ids: set[str] = set()
    for index, pair in enumerate(pairs):
        pair_id = _identifier(pair, ("pair_cell_id", "id"), f"pairs[{index}]")
        if pair_id in validated_pair_ids:
            raise ValueError(f"duplicate pair cell identifier: {pair_id}")
        validated_pair_ids.add(pair_id)
        _route_token(pair_id, "pair cell identifier")
        target_id = _identifier(
            pair, ("target_id", "target_key", "protein_id"), f"pairs[{index}] target"
        )
        drug_id = _identifier(
            pair, ("drug_id", "ligand_canonical_id"), f"pairs[{index}] drug"
        )
        if target_id not in targets:
            raise ValueError(f"pair {pair_id} references unknown target: {target_id}")
        if drug_id not in drugs:
            raise ValueError(f"pair {pair_id} references unknown drug: {drug_id}")

    unknown_artifact_pairs = sorted(set(artifacts_by_pair) - validated_pair_ids)
    if unknown_artifact_pairs:
        preview = ", ".join(repr(value) for value in unknown_artifact_pairs[:5])
        suffix = "" if len(unknown_artifact_pairs) <= 5 else ", ..."
        raise ValueError(
            f"artifact references unknown pair cell identifier(s): {preview}{suffix}"
        )

    _prepare_output(output_dir, overwrite)
    _write_json(
        output_dir / _BUNDLE_MARKER,
        {
            "schema_version": 1,
            "generator": "atlas publish edge-bundle",
            "state": "building",
        },
    )

    objects_root = output_dir / "objects"
    public_root = output_dir / "public"
    prefix = f"releases/{release_token}"
    object_entries: list[dict[str, Any]] = []
    target_pairs: dict[str, list[dict[str, Any]]] = defaultdict(list)
    drug_pairs: dict[str, list[dict[str, Any]]] = defaultdict(list)
    pair_shards: dict[str, list[dict[str, Any]]] = defaultdict(list)
    pair_ids_seen: set[str] = set()

    for index, pair in enumerate(pairs):
        pair_id = _identifier(pair, ("pair_cell_id", "id"), f"pairs[{index}]")
        if pair_id in pair_ids_seen:
            raise ValueError(f"duplicate pair cell identifier: {pair_id}")
        pair_ids_seen.add(pair_id)
        target_id = _identifier(
            pair,
            ("target_id", "target_key", "protein_id"),
            f"pairs[{index}] target",
        )
        drug_id = _identifier(
            pair,
            ("drug_id", "ligand_canonical_id"),
            f"pairs[{index}] drug",
        )
        if target_id not in targets:
            raise ValueError(f"pair {pair_id} references unknown target: {target_id}")
        if drug_id not in drugs:
            raise ValueError(f"pair {pair_id} references unknown drug: {drug_id}")
        pair_route_id = _route_token(pair_id, "pair cell identifier")
        summary = _pair_summary(
            pair,
            pair_id=pair_id,
            pair_route_id=pair_route_id,
            target_id=target_id,
            target_route_id=target_routes[target_id],
            target_label=target_labels[target_id],
            drug_id=drug_id,
            drug_route_id=drug_routes[drug_id],
            drug_label=drug_labels[drug_id],
        )
        target_pairs[target_id].append(_counterpart_summary(summary, "targets"))
        drug_pairs[drug_id].append(_counterpart_summary(summary, "drugs"))
        shard = hashlib.sha256(pair_route_id.encode("ascii")).hexdigest()[:2]
        pair_shards[shard].append(summary)
        pair_record = {
            "schema_version": 1,
            "release_id": release_id,
            "release_token": release_token,
            "id": pair_id,
            "route_id": pair_route_id,
            "target": {
                "id": target_id,
                "route_id": target_routes[target_id],
                "label": target_labels[target_id],
            },
            "drug": {
                "id": drug_id,
                "route_id": drug_routes[drug_id],
                "label": drug_labels[drug_id],
            },
            "pair": pair,
            "artifacts": artifacts_by_pair.get(pair_id, []),
        }
        _write_object(
            objects_root,
            f"{prefix}/records/pairs/{pair_route_id}.json",
            pair_record,
            object_entries,
        )

    target_index: list[dict[str, Any]] = []
    for target_id, entity in targets.items():
        entity_pairs = sorted(target_pairs.get(target_id, []), key=_pair_sort_key)
        record = {
            "schema_version": 1,
            "release_id": release_id,
            "release_token": release_token,
            "id": target_id,
            "route_id": target_routes[target_id],
            "entity": entity,
            "pairs": entity_pairs,
        }
        _write_object(
            objects_root,
            f"{prefix}/records/targets/{target_routes[target_id]}.json",
            record,
            object_entries,
        )
        target_index.append(
            {
                "id": target_id,
                "route_id": target_routes[target_id],
                "label": target_labels[target_id],
                "pair_count": len(entity_pairs),
                "primary_score_count": sum(
                    int(row.get("final_score") is not None) for row in entity_pairs
                ),
            }
        )

    drug_index: list[dict[str, Any]] = []
    for drug_id, entity in drugs.items():
        entity_pairs = sorted(drug_pairs.get(drug_id, []), key=_pair_sort_key)
        record = {
            "schema_version": 1,
            "release_id": release_id,
            "release_token": release_token,
            "id": drug_id,
            "route_id": drug_routes[drug_id],
            "entity": entity,
            "pairs": entity_pairs,
        }
        _write_object(
            objects_root,
            f"{prefix}/records/drugs/{drug_routes[drug_id]}.json",
            record,
            object_entries,
        )
        drug_index.append(
            {
                "id": drug_id,
                "route_id": drug_routes[drug_id],
                "label": drug_labels[drug_id],
                "pair_count": len(entity_pairs),
                "primary_score_count": sum(
                    int(row.get("final_score") is not None) for row in entity_pairs
                ),
            }
        )

    target_index.sort(key=lambda row: (str(row["label"]).casefold(), str(row["id"])))
    drug_index.sort(key=lambda row: (str(row["label"]).casefold(), str(row["id"])))
    _write_object(
        objects_root,
        f"{prefix}/indexes/targets.json",
        {
            "schema_version": 1,
            "release_token": release_token,
            "kind": "targets",
            "records": target_index,
        },
        object_entries,
    )
    _write_object(
        objects_root,
        f"{prefix}/indexes/drugs.json",
        {
            "schema_version": 1,
            "release_token": release_token,
            "kind": "drugs",
            "records": drug_index,
        },
        object_entries,
    )

    shard_index: list[dict[str, Any]] = []
    for shard, shard_rows in sorted(pair_shards.items()):
        shard_rows.sort(
            key=lambda row: (
                str(row["target_label"]).casefold(),
                str(row["drug_label"]).casefold(),
                str(row["id"]),
            )
        )
        entry = _write_object(
            objects_root,
            f"{prefix}/indexes/pairs/{shard}.json",
            {
                "schema_version": 1,
                "release_token": release_token,
                "kind": "pairs",
                "shard": shard,
                "records": shard_rows,
            },
            object_entries,
        )
        shard_index.append(
            {
                "id": shard,
                "count": len(shard_rows),
                "sha256": entry["sha256"],
                "size_bytes": entry["size_bytes"],
                "api_path": f"/api/releases/{release_token}/indexes/pairs/{shard}",
            }
        )
    _write_object(
        objects_root,
        f"{prefix}/indexes/pairs.json",
        {
            "schema_version": 1,
            "release_token": release_token,
            "kind": "pair_shards",
            "count": len(pairs),
            "shards": shard_index,
        },
        object_entries,
    )

    coverage = _required_mapping(payload.get("coverage", {}), "coverage")
    release_manifest = {
        "schema_version": 1,
        "release": release,
        "release_token": release_token,
        "source_payload_sha256": source_sha256,
        "counts": {
            "targets": len(targets),
            "drugs": len(drugs),
            "pairs": len(pairs),
            "pair_index_shards": len(shard_index),
        },
        "coverage": coverage,
        "score_contract": payload.get("score_contract", {}),
        "scientific_policies": payload.get("scientific_policies", {}),
        "downloads": (
            download_projection["public_entries"] if download_projection else []
        ),
        "route_templates": {
            "release": f"/releases/{release_token}",
            "target": f"/releases/{release_token}/targets/{{target_route_id}}",
            "drug": f"/releases/{release_token}/drugs/{{drug_route_id}}",
            "pair": f"/releases/{release_token}/pairs/{{pair_route_id}}",
        },
        "api_templates": {
            "target": f"/api/releases/{release_token}/targets/{{target_route_id}}",
            "drug": f"/api/releases/{release_token}/drugs/{{drug_route_id}}",
            "pair": f"/api/releases/{release_token}/pairs/{{pair_route_id}}",
        },
    }
    _write_object(
        objects_root,
        f"{prefix}/manifest.json",
        release_manifest,
        object_entries,
    )

    _write_bytes(public_root / "index.html", INDEX_HTML.encode("utf-8"))
    _write_bytes(public_root / "_headers", PUBLIC_HEADERS.encode("utf-8"))
    _write_bytes(public_root / "assets" / "app.css", APP_CSS.encode("utf-8"))
    _write_bytes(public_root / "assets" / "app.js", APP_JS.encode("utf-8"))
    _write_json(
        public_root / "release-config.json",
        {"release_id": release_id, "release_token": release_token},
    )
    _write_bytes(output_dir / "src" / "index.mjs", WORKER_SOURCE.encode("utf-8"))
    _write_bytes(output_dir / "wrangler.toml", WRANGLER_TEMPLATE.encode("utf-8"))

    object_entries.sort(key=lambda row: str(row["key"]))
    object_manifest = {
        "schema_version": 1,
        "release_id": release_id,
        "release_token": release_token,
        "source": {
            "file": "release_browser.json",
            "sha256": source_sha256,
        },
        "object_prefix": prefix,
        "immutability": "content-qualified release prefix; never overwrite object keys",
        "objects": object_entries,
    }
    _write_json(output_dir / "object_manifest.json", object_manifest)
    if download_projection is not None:
        _write_json(output_dir / DOWNLOAD_MANIFEST_NAME, download_projection)
    summary = {
        "schema_version": 1,
        "release_id": release_id,
        "release_token": release_token,
        "source_payload_sha256": source_sha256,
        "input_mode": "bounded_browser_json",
        "publication_scale_supported": False,
        "bounded_input": {
            "max_browser_payload_bytes": MAX_BROWSER_PAYLOAD_BYTES,
            "max_pair_count": MAX_BROWSER_PAIR_COUNT,
            "memory_shape": (
                "browser payload plus pair/entity summary projections are materialized"
            ),
        },
        "static_fallback_modified": False,
        "deploy_performed": False,
        "counts": release_manifest["counts"],
        "r2_object_count": len(object_entries),
        "r2_bytes": sum(int(row["size_bytes"]) for row in object_entries),
        "entrypoint": f"/releases/{release_token}",
        "object_manifest": "object_manifest.json",
        "download_projection": (
            DOWNLOAD_MANIFEST_NAME if download_projection is not None else None
        ),
        "download_count": (
            int(download_projection["upload_count"]) if download_projection else 0
        ),
        "download_bytes": (
            int(download_projection["upload_bytes"]) if download_projection else 0
        ),
        "worker": "src/index.mjs",
        "wrangler_template": "wrangler.toml",
    }
    _write_json(output_dir / "edge_bundle_summary.json", summary)
    _write_json(
        output_dir / _BUNDLE_MARKER,
        {
            "schema_version": 1,
            "generator": "atlas publish edge-bundle",
            "release_token": release_token,
        },
    )
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas publish edge-bundle",
        description=(
            "Build an optional dependency-free SPA, read-only Worker, and immutable "
            "R2 object tree from an existing public Atlas site."
        ),
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--site-dir", type=Path)
    source.add_argument(
        "--database",
        type=Path,
        help=(
            "Stream a public Atlas SQLite snapshot into coarse pair shards "
            "without materializing the pair matrix"
        ),
    )
    parser.add_argument("--source-site-dir", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--batch-rows", type=int, default=1000)
    parser.add_argument("--coarse-shard-rows", type=int, default=2000)
    parser.add_argument("--max-pairs", type=int, default=2_000_000)
    parser.add_argument(
        "--download-base-url",
        help=(
            "HTTPS object-origin base URL used to project verified release "
            "downloads; no files are uploaded"
        ),
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    try:
        if args.database is not None:
            from analysis.reporting.docking_atlas_edge_sqlite import (
                build_streaming_edge_bundle,
            )

            database_path = Path(args.database)
            default_root = (
                database_path.parent.parent.parent
                if database_path.parent.name == "downloads"
                else database_path.parent
            )
            output_dir = args.out_dir or default_root / "edge"
            summary = build_streaming_edge_bundle(
                database_path,
                output_dir,
                source_site_dir=args.source_site_dir,
                overwrite=args.overwrite,
                download_base_url=args.download_base_url,
                batch_rows=args.batch_rows,
                coarse_shard_rows=args.coarse_shard_rows,
                max_pairs=args.max_pairs,
            )
        else:
            site_dir = Path(args.site_dir)
            output_dir = args.out_dir or site_dir.parent / "edge"
            summary = build_edge_bundle(
                site_dir,
                output_dir,
                overwrite=args.overwrite,
                download_base_url=args.download_base_url,
            )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["build_edge_bundle", "main"]
