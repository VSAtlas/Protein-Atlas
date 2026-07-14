from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from analysis.reporting.docking_atlas_delivery import build_deployment_preflight
from analysis.reporting.docking_atlas_edge import build_edge_bundle
from cli.qol.publish import _cmd_publish


def _write_public_site(site_dir: Path) -> bytes:
    payload = {
        "release": {
            "id": "Atlas v0.1 demo",
            "title": "Docking Atlas v0.1",
            "summary": "Auditable test release",
        },
        "coverage": {"pair_count": 2, "rank_eligible_count": 1},
        "score_contract": {
            "primary_field": "final_score",
            "direction": "higher_is_better",
            "no_fallback": True,
        },
        "scientific_policies": {"native_redocking": {"threshold": "frozen"}},
        "targets": [
            {
                "id": "run-a|1ABC|HOLO|pH7.4",
                "display_name": "Target alpha",
                "pair_count": 2,
            }
        ],
        "ligands": [
            {"id": "Drug A", "display_name": "Drug alpha"},
            {"id": "Drug/B", "display_name": "Drug beta"},
        ],
        "pairs": [
            {
                "pair_cell_id": 11,
                "target_id": "run-a|1ABC|HOLO|pH7.4",
                "drug_id": "Drug A",
                "final_status": "valid",
                "final_score": 4.5,
                "rank_eligible": 1,
                "rank_within_receptor": 1,
                "rank_across_receptors": 1,
            },
            {
                "pair_cell_id": 12,
                "target_id": "run-a|1ABC|HOLO|pH7.4",
                "drug_id": "Drug/B",
                "final_status": "failed",
                "failure_code": "engine_error",
                "failure_reason": "structured failure retained",
                "rank_eligible": 0,
                "ranking_eligibility_reason": "missing_final_score",
            },
        ],
        "artifacts": [
            {
                "artifact_id": 7,
                "pair_cell_id": 11,
                "member_name": "pose.sdf",
                "sha256": "a" * 64,
            }
        ],
    }
    data = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    browser = site_dir / "downloads" / "release_browser.json"
    browser.parent.mkdir(parents=True, exist_ok=True)
    browser.write_bytes(data)
    (site_dir / "index.html").write_text("static fallback", encoding="utf-8")
    return data


def _declare_download(site_dir: Path) -> Path:
    path = site_dir / "downloads" / "pairs.csv"
    path.write_text("drug,score\nA,4.5\n", encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    (site_dir / "site_manifest.json").write_text(
        json.dumps(
            {
                "release_id": "Atlas v0.1 demo",
                "downloads": [
                    {
                        "label": "Complete pair table (CSV)",
                        "url": "downloads/pairs.csv",
                        "content_hash": f"sha256:{digest}",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return path


def test_edge_bundle_decomposes_release_into_stable_r2_records(tmp_path: Path) -> None:
    site_dir = tmp_path / "site"
    source = _write_public_site(site_dir)
    static_hash = hashlib.sha256((site_dir / "index.html").read_bytes()).hexdigest()

    output_dir = tmp_path / "edge"
    summary = build_edge_bundle(site_dir, output_dir)

    expected_suffix = hashlib.sha256(source).hexdigest()[:16]
    release_token = summary["release_token"]
    assert release_token == f"atlas-v0-1-demo-{expected_suffix}"
    assert summary["static_fallback_modified"] is False
    assert summary["deploy_performed"] is False
    assert summary["input_mode"] == "bounded_browser_json"
    assert summary["publication_scale_supported"] is False
    assert summary["bounded_input"]["max_pair_count"] == 50_000
    assert summary["counts"] == {
        "targets": 1,
        "drugs": 2,
        "pairs": 2,
        "pair_index_shards": 2,
    }
    assert (
        hashlib.sha256((site_dir / "index.html").read_bytes()).hexdigest()
        == static_hash
    )

    manifest = json.loads(
        (
            output_dir / "objects" / "releases" / release_token / "manifest.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["route_templates"] == {
        "release": f"/releases/{release_token}",
        "target": f"/releases/{release_token}/targets/{{target_route_id}}",
        "drug": f"/releases/{release_token}/drugs/{{drug_route_id}}",
        "pair": f"/releases/{release_token}/pairs/{{pair_route_id}}",
    }

    target_index = json.loads(
        (
            output_dir
            / "objects"
            / "releases"
            / release_token
            / "indexes"
            / "targets.json"
        ).read_text(encoding="utf-8")
    )
    target_route = target_index["records"][0]["route_id"]
    target_record = json.loads(
        (
            output_dir
            / "objects"
            / "releases"
            / release_token
            / "records"
            / "targets"
            / f"{target_route}.json"
        ).read_text(encoding="utf-8")
    )
    assert [row["final_status"] for row in target_record["pairs"]] == [
        "valid",
        "failed",
    ]
    failed_route = target_record["pairs"][1]["pair_route_id"]
    failed_record = json.loads(
        (
            output_dir
            / "objects"
            / "releases"
            / release_token
            / "records"
            / "pairs"
            / f"{failed_route}.json"
        ).read_text(encoding="utf-8")
    )
    assert failed_record["pair"]["failure_reason"] == "structured failure retained"

    object_manifest = json.loads(
        (output_dir / "object_manifest.json").read_text(encoding="utf-8")
    )
    assert object_manifest["object_prefix"] == f"releases/{release_token}"
    assert summary["r2_object_count"] == len(object_manifest["objects"])
    assert all(
        item["key"].startswith(f"releases/{release_token}/")
        for item in object_manifest["objects"]
    )
    assert (output_dir / "public" / "index.html").is_file()
    assert (output_dir / "public" / "_headers").is_file()
    assert (output_dir / "public" / "assets" / "app.js").is_file()
    assert (output_dir / "src" / "index.mjs").is_file()
    assert (output_dir / "wrangler.toml").is_file()


def test_edge_bundle_is_deterministic_and_only_overwrites_its_own_output(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    first = tmp_path / "edge-one"
    second = tmp_path / "edge-two"
    build_edge_bundle(site_dir, first)
    build_edge_bundle(site_dir, second)

    assert (first / "object_manifest.json").read_bytes() == (
        second / "object_manifest.json"
    ).read_bytes()
    with pytest.raises(ValueError, match="use --overwrite"):
        build_edge_bundle(site_dir, first)

    unrelated = tmp_path / "unrelated"
    unrelated.mkdir()
    (unrelated / "keep.txt").write_text("keep", encoding="utf-8")
    with pytest.raises(ValueError, match="unrecognized directory"):
        build_edge_bundle(site_dir, unrelated, overwrite=True)
    assert (unrelated / "keep.txt").read_text(encoding="utf-8") == "keep"

    rebuilt = build_edge_bundle(site_dir, first, overwrite=True)
    assert rebuilt["release_token"]


def test_edge_bundle_rejects_private_paths_and_pair_reference_drift(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    browser = site_dir / "downloads" / "release_browser.json"
    payload = json.loads(browser.read_text(encoding="utf-8"))
    payload["release"]["summary"] = "/stor/private/result"
    browser.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="machine-local path"):
        build_edge_bundle(site_dir, tmp_path / "edge-private")

    _write_public_site(site_dir)
    payload = json.loads(browser.read_text(encoding="utf-8"))
    payload["pairs"][0]["drug_id"] = "unknown-drug"
    browser.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown drug"):
        build_edge_bundle(site_dir, tmp_path / "edge-drift")


def test_edge_bundle_rejects_unknown_non_null_artifact_pair_reference(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    browser = site_dir / "downloads" / "release_browser.json"
    payload = json.loads(browser.read_text(encoding="utf-8"))
    payload["artifacts"][0]["pair_cell_id"] = 999
    browser.write_text(json.dumps(payload), encoding="utf-8")
    output_dir = tmp_path / "edge-artifact-drift"

    with pytest.raises(ValueError, match="artifact references unknown pair"):
        build_edge_bundle(site_dir, output_dir)

    assert not output_dir.exists()


def test_edge_bundle_rejects_inputs_above_explicit_memory_bounds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from analysis.reporting import docking_atlas_edge

    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    too_large = tmp_path / "edge-too-large"
    monkeypatch.setattr(docking_atlas_edge, "MAX_BROWSER_PAYLOAD_BYTES", 32)
    with pytest.raises(ValueError, match="loads browser JSON in memory"):
        docking_atlas_edge.build_edge_bundle(site_dir, too_large)
    assert not too_large.exists()

    monkeypatch.setattr(docking_atlas_edge, "MAX_BROWSER_PAYLOAD_BYTES", 10**9)
    monkeypatch.setattr(docking_atlas_edge, "MAX_BROWSER_PAIR_COUNT", 1)
    too_many = tmp_path / "edge-too-many"
    with pytest.raises(ValueError, match="conference/intermediate snapshots"):
        docking_atlas_edge.build_edge_bundle(site_dir, too_many)
    assert not too_many.exists()


def test_generated_worker_is_read_only_and_path_constrained(tmp_path: Path) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    output_dir = tmp_path / "edge"
    build_edge_bundle(site_dir, output_dir)

    worker = (output_dir / "src" / "index.mjs").read_text(encoding="utf-8")
    assert 'request.method !== "GET" && request.method !== "HEAD"' in worker
    assert "ATLAS_RELEASES.get(key)" in worker
    assert "ATLAS_RELEASES.head(key)" in worker
    assert ".put(" not in worker
    assert ".list(" not in worker
    assert ".delete(" not in worker
    assert 'parts[0] !== "api"' in worker
    assert "TOKEN.test" in worker
    assert "secureStaticResponse(await env.ASSETS.fetch(request))" in worker
    assert '"Content-Security-Policy"' in worker
    assert '"X-Frame-Options": "DENY"' in worker
    assert '"X-Content-Type-Options": "nosniff"' in worker
    wrangler = (output_dir / "wrangler.toml").read_text(encoding="utf-8")
    assert 'run_worker_first = ["/api/*"]' in wrangler
    assert 'not_found_handling = "single-page-application"' in wrangler


def test_generated_static_assets_apply_security_headers_without_worker_routing(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    output_dir = tmp_path / "edge"
    build_edge_bundle(site_dir, output_dir)

    headers = (output_dir / "public" / "_headers").read_text(encoding="utf-8")
    assert headers.startswith("/*\n")
    assert "Content-Security-Policy: default-src 'self'" in headers
    assert "Permissions-Policy:" in headers
    assert "Referrer-Policy: no-referrer" in headers
    assert "X-Content-Type-Options: nosniff" in headers
    assert "X-Frame-Options: DENY" in headers
    assert "/assets/*\n  Cache-Control:" in headers


def test_publish_cli_forwards_edge_bundle_without_touching_build(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from analysis.reporting import docking_atlas_edge

    captured: list[str] = []

    def fake_main(forwarded: list[str]) -> int:
        captured.extend(forwarded)
        return 9

    monkeypatch.setattr(docking_atlas_edge, "main", fake_main)
    assert (
        _cmd_publish(
            [
                "edge-bundle",
                "--site-dir",
                "release/site",
                "--out-dir",
                "release/edge",
                "--overwrite",
                "--download-base-url",
                "https://downloads.example.org",
            ]
        )
        == 9
    )
    assert captured == [
        "--site-dir",
        "release/site",
        "--out-dir",
        "release/edge",
        "--overwrite",
        "--download-base-url",
        "https://downloads.example.org",
    ]


def test_publish_cli_forwards_streaming_sqlite_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from analysis.reporting import docking_atlas_edge

    captured: list[str] = []

    def fake_main(forwarded: list[str]) -> int:
        captured.extend(forwarded)
        return 11

    monkeypatch.setattr(docking_atlas_edge, "main", fake_main)
    assert (
        _cmd_publish(
            [
                "edge-bundle",
                "--database",
                "release/site/downloads/docking_atlas.sqlite",
                "--source-site-dir",
                "release/site",
                "--out-dir",
                "release/edge",
                "--batch-rows",
                "500",
                "--coarse-shard-rows",
                "4000",
                "--max-pairs",
                "800000",
                "--download-base-url",
                "https://downloads.example.org",
            ]
        )
        == 11
    )
    assert captured == [
        "--database",
        "release/site/downloads/docking_atlas.sqlite",
        "--source-site-dir",
        "release/site",
        "--out-dir",
        "release/edge",
        "--download-base-url",
        "https://downloads.example.org",
        "--batch-rows",
        "500",
        "--coarse-shard-rows",
        "4000",
        "--max-pairs",
        "800000",
    ]


def test_edge_bundle_rejects_source_ancestor_output_and_forged_marker(
    tmp_path: Path,
) -> None:
    release_root = tmp_path / "release-root"
    site_dir = release_root / "site"
    _write_public_site(site_dir)
    (release_root / ".atlas-edge-bundle.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "generator": "atlas publish edge-bundle",
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="disjoint directories"):
        build_edge_bundle(site_dir, release_root, overwrite=True)
    assert (site_dir / "index.html").read_text(encoding="utf-8") == "static fallback"

    forged = tmp_path / "forged-edge"
    forged.mkdir()
    (forged / ".atlas-edge-bundle.json").write_text("{}", encoding="utf-8")
    keep = forged / "keep.txt"
    keep.write_text("keep", encoding="utf-8")
    with pytest.raises(ValueError, match="valid Atlas edge bundle marker"):
        build_edge_bundle(site_dir, forged, overwrite=True)
    assert keep.read_text(encoding="utf-8") == "keep"


def test_edge_bundle_projects_verified_downloads_to_provider_neutral_urls(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    download = _declare_download(site_dir)
    output_dir = tmp_path / "edge"

    summary = build_edge_bundle(
        site_dir,
        output_dir,
        download_base_url="https://downloads.example.org/atlas",
    )

    projection = json.loads(
        (output_dir / "download_projection.json").read_text(encoding="utf-8")
    )
    expected_url = (
        "https://downloads.example.org/atlas/releases/"
        f"{summary['release_token']}/downloads/pairs.csv"
    )
    assert projection["provider_contract"] == "https_object_origin"
    assert projection["upload_count"] == 1
    assert projection["upload_bytes"] == download.stat().st_size
    assert projection["uploads"][0]["public_url"] == expected_url
    assert projection["uploads"][0]["source_path"] == "downloads/pairs.csv"
    assert not Path(projection["uploads"][0]["source_path"]).is_absolute()
    release_manifest = json.loads(
        (
            output_dir
            / "objects"
            / "releases"
            / summary["release_token"]
            / "manifest.json"
        ).read_text(encoding="utf-8")
    )
    assert release_manifest["downloads"] == projection["public_entries"]
    app = (output_dir / "public" / "assets" / "app.js").read_text(encoding="utf-8")
    assert "Download frozen release" in app
    assert 'target.protocol !== "https:"' in app


def test_download_projection_rejects_insecure_url_and_hash_drift(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    download = _declare_download(site_dir)

    with pytest.raises(ValueError, match="HTTPS origin"):
        build_edge_bundle(
            site_dir,
            tmp_path / "edge-http",
            download_base_url="http://downloads.example.org",
        )
    assert not (tmp_path / "edge-http").exists()

    download.write_text("mutated\n", encoding="utf-8")
    with pytest.raises(ValueError, match="hash does not match"):
        build_edge_bundle(
            site_dir,
            tmp_path / "edge-hash",
            download_base_url="https://downloads.example.org",
        )
    assert not (tmp_path / "edge-hash").exists()


def test_cloudflare_preflight_is_credential_free_and_reports_exact_resources(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    _declare_download(site_dir)
    edge_dir = tmp_path / "edge"
    build_edge_bundle(
        site_dir,
        edge_dir,
        download_base_url="https://downloads.example.org",
    )

    report = build_deployment_preflight(
        edge_dir,
        site_dir=site_dir,
        provider="cloudflare-workers-r2",
        cloudflare_plan="free",
        worker_name="docking-atlas-edge",
        production_bucket="atlas-docking-releases",
        preview_bucket="atlas-docking-releases-preview",
    )

    assert report["status"] == "ready"
    assert report["provider"]["static_asset_count_limit"] == 20_000
    assert report["provider"]["static_asset_size_limit_bytes"] == 25 * 1024 * 1024
    assert report["credential_environment_names"] == [
        "CLOUDFLARE_ACCOUNT_ID",
        "CLOUDFLARE_API_TOKEN",
    ]
    assert report["secret_values_inspected"] is False
    assert report["network_checks_performed"] is False
    assert report["upload_performed"] is False
    assert report["deployment_performed"] is False
    assert [row["kind"] for row in report["required_resources"]] == [
        "worker_service",
        "r2_bucket",
        "r2_bucket",
        "https_download_origin",
    ]
    config = (edge_dir / "wrangler.preflight.toml").read_text(encoding="utf-8")
    assert 'name = "docking-atlas-edge"' in config
    assert 'bucket_name = "atlas-docking-releases"' in config
    assert 'preview_bucket_name = "atlas-docking-releases-preview"' in config
    assert "CLOUDFLARE_API_TOKEN" not in config


def test_preflight_blocks_unselected_resources_and_detects_download_drift(
    tmp_path: Path,
) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    download = _declare_download(site_dir)
    edge_dir = tmp_path / "edge"
    build_edge_bundle(
        site_dir,
        edge_dir,
        download_base_url="https://downloads.example.org/atlas",
    )

    blocked = build_deployment_preflight(
        edge_dir, site_dir=site_dir, provider="cloudflare-workers-r2"
    )
    assert blocked["status"] == "blocked"
    assert {row["code"] for row in blocked["blockers"]} >= {
        "cloudflare_worker_name_required",
        "cloudflare_production_bucket_required",
        "cloudflare_preview_bucket_required",
        "cloudflare_r2_download_origin_must_not_have_path",
    }

    path_blocked = build_deployment_preflight(
        edge_dir,
        site_dir=site_dir,
        provider="cloudflare-workers-r2",
        worker_name="docking-atlas-edge",
        production_bucket="atlas-docking-releases",
        preview_bucket="atlas-docking-releases-preview",
    )
    assert {row["code"] for row in path_blocked["blockers"]} == {
        "cloudflare_r2_download_origin_must_not_have_path"
    }
    assert not (edge_dir / "wrangler.preflight.toml").exists()

    download.write_text("changed after projection\n", encoding="utf-8")
    with pytest.raises(ValueError, match="integrity mismatch"):
        build_deployment_preflight(
            edge_dir,
            site_dir=site_dir,
            provider="generic",
        )


def test_publish_cli_forwards_preflight_without_deploying(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from analysis.reporting import docking_atlas_delivery

    captured: list[str] = []

    def fake_main(forwarded: list[str]) -> int:
        captured.extend(forwarded)
        return 7

    monkeypatch.setattr(docking_atlas_delivery, "main", fake_main)
    assert (
        _cmd_publish(
            [
                "preflight",
                "--edge-dir",
                "release/edge",
                "--site-dir",
                "release/site",
                "--provider",
                "cloudflare-workers-r2",
                "--worker-name",
                "docking-atlas-edge",
                "--production-bucket",
                "atlas-docking-releases",
                "--preview-bucket",
                "atlas-docking-releases-preview",
            ]
        )
        == 7
    )
    assert captured == [
        "--edge-dir",
        "release/edge",
        "--provider",
        "cloudflare-workers-r2",
        "--cloudflare-plan",
        "free",
        "--site-dir",
        "release/site",
        "--worker-name",
        "docking-atlas-edge",
        "--production-bucket",
        "atlas-docking-releases",
        "--preview-bucket",
        "atlas-docking-releases-preview",
    ]


def test_preflight_rejects_tampered_worker_source(tmp_path: Path) -> None:
    site_dir = tmp_path / "site"
    _write_public_site(site_dir)
    _declare_download(site_dir)
    edge_dir = tmp_path / "edge"
    build_edge_bundle(
        site_dir,
        edge_dir,
        download_base_url="https://downloads.example.org",
    )
    (edge_dir / "src" / "index.mjs").write_text("tampered\n", encoding="utf-8")

    with pytest.raises(ValueError, match="support file integrity mismatch"):
        build_deployment_preflight(
            edge_dir,
            site_dir=site_dir,
            provider="generic",
        )
