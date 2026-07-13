from __future__ import annotations

import json
from pathlib import Path

import pytest

from analysis.atlas_database.release_bundle import (
    CHECKSUM_NAME,
    INVENTORY_NAME,
    prepare_release_bundle,
    verify_release_bundle,
)


def _site(tmp_path: Path) -> Path:
    site = tmp_path / "site"
    (site / "assets").mkdir(parents=True)
    (site / "downloads").mkdir()
    (site / "index.html").write_text(
        '<script src="assets/app.js"></script><a href="downloads/data.csv">data</a>',
        encoding="utf-8",
    )
    (site / "assets" / "app.js").write_text('console.log("atlas");\n', encoding="utf-8")
    (site / "downloads" / "data.csv").write_text(
        "drug,score\nA,-9.0\n", encoding="utf-8"
    )
    (site / "site_manifest.json").write_text(
        json.dumps(
            {
                "release_id": "atlas-v0.1",
                "entrypoint": "index.html",
                "downloads": [{"href": "downloads/data.csv"}],
            }
        ),
        encoding="utf-8",
    )
    return site


def _error_codes(report: dict[str, object]) -> set[str]:
    errors = report["errors"]
    assert isinstance(errors, list)
    return {
        str(error["code"])
        for error in errors
        if isinstance(error, dict) and "code" in error
    }


def test_release_inventory_and_checksums_are_deterministic(tmp_path: Path) -> None:
    site = _site(tmp_path)

    first = prepare_release_bundle(site)
    inventory = (site / INVENTORY_NAME).read_bytes()
    checksums = (site / CHECKSUM_NAME).read_bytes()
    second = prepare_release_bundle(site)

    assert first["status"] == "passed"
    assert second["status"] == "passed"
    assert (site / INVENTORY_NAME).read_bytes() == inventory
    assert (site / CHECKSUM_NAME).read_bytes() == checksums
    headers = (site / "_headers").read_text(encoding="utf-8")
    assert "script-src 'self'" in headers
    assert "style-src 'self' 'unsafe-inline'" in headers


def test_verifier_detects_payload_mutation(tmp_path: Path) -> None:
    site = _site(tmp_path)
    assert prepare_release_bundle(site)["status"] == "passed"

    (site / "assets" / "app.js").write_text("mutated\n", encoding="utf-8")
    report = verify_release_bundle(site)

    assert report["status"] == "failed"
    assert {"checksum_or_size_mismatch", "declared_checksum_mismatch"} <= _error_codes(
        report
    )


def test_verifier_fails_for_missing_and_escaping_links(tmp_path: Path) -> None:
    site = _site(tmp_path)
    (tmp_path / "outside.html").write_text("outside\n", encoding="utf-8")
    (site / "index.html").write_text(
        '<a href="missing.html">missing</a><a href="../outside.html">outside</a>',
        encoding="utf-8",
    )

    report = prepare_release_bundle(site)

    assert report["status"] == "failed"
    missing = [
        error for error in report["errors"] if error["code"] == "missing_local_link"
    ]
    assert {error["target"] for error in missing} == {"missing.html", "../outside.html"}


def test_verifier_fails_for_private_path_leak(tmp_path: Path) -> None:
    site = _site(tmp_path)
    (site / "assets" / "app.js").write_text(
        'const source = "/stor/private/run/result.csv";\n',
        encoding="utf-8",
    )

    report = prepare_release_bundle(site)

    assert report["status"] == "failed"
    assert "private_path_in_text" in _error_codes(report)


def test_inventory_does_not_follow_symlinks_and_verifier_rejects_them(
    tmp_path: Path,
) -> None:
    site = _site(tmp_path)
    outside = tmp_path / "outside.txt"
    outside.write_text("external content\n", encoding="utf-8")
    link = site / "assets" / "external.txt"
    link.symlink_to(outside)

    report = prepare_release_bundle(site)
    inventory = json.loads((site / INVENTORY_NAME).read_text(encoding="utf-8"))

    assert report["status"] == "failed"
    assert "symlink_not_allowed" in _error_codes(report)
    assert "assets/external.txt" not in {row["path"] for row in inventory["files"]}


def test_verifier_rejects_duplicate_or_tampered_inventory_metadata(
    tmp_path: Path,
) -> None:
    site = _site(tmp_path)
    assert prepare_release_bundle(site)["status"] == "passed"
    inventory_path = site / INVENTORY_NAME
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    inventory["files"].append(dict(inventory["files"][0]))
    inventory["file_count"] += 1
    inventory["total_size_bytes"] += inventory["files"][0]["size_bytes"]
    inventory_path.write_text(json.dumps(inventory), encoding="utf-8")

    report = verify_release_bundle(site)

    assert report["status"] == "failed"
    assert "inventory_invalid" in _error_codes(report)


def test_verifier_rejects_untracked_custom_report_inside_site(
    tmp_path: Path,
) -> None:
    site = _site(tmp_path)
    assert prepare_release_bundle(site)["status"] == "passed"
    custom_report = site / "custom-verification.json"

    with pytest.raises(ValueError, match="outside the deployable site"):
        verify_release_bundle(site, report_path=custom_report)

    assert not custom_report.exists()
    assert verify_release_bundle(site)["status"] == "passed"


def test_verifier_rejects_symlinked_canonical_report_without_outside_write(
    tmp_path: Path,
) -> None:
    site = _site(tmp_path)
    assert prepare_release_bundle(site)["status"] == "passed"
    outside = tmp_path / "outside-report.json"
    outside.write_text("keep outside\n", encoding="utf-8")
    canonical_report = site / "release_verification.json"
    canonical_report.unlink()
    canonical_report.symlink_to(outside)

    with pytest.raises(ValueError, match="symlinked canonical verification report"):
        verify_release_bundle(site)

    assert canonical_report.is_symlink()
    assert outside.read_text(encoding="utf-8") == "keep outside\n"
