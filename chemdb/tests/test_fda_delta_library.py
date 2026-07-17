from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from prep_ligands import fda_delta_library as delta


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _valid_summary() -> tuple[dict[str, object], dict[str, str]]:
    source_hashes = {
        "manifest": "a" * 64,
        "quarantine": "b" * 64,
    }
    summary: dict[str, object] = {
        "manifest_sha256": source_hashes["manifest"],
        "quarantine_sha256": source_hashes["quarantine"],
        "action_counts": {"legacy_byte_copy": 1264, "prepare_parent": 599},
        "status_counts": {"ready": 1850, "quarantined_non_dockable": 13},
        "canonical_parent_rows": 1863,
        "ready_pdbqt_count": 1850,
    }
    return summary, source_hashes


def test_link_or_copy_copies_and_resumes_exact_bytes(tmp_path: Path) -> None:
    source = tmp_path / "source.pdbqt"
    target = tmp_path / "target.pdbqt"
    source.write_bytes(b"REMARK exact FDA parent\n")
    digest = _sha256(source)

    assert delta._link_or_copy(source, target, digest, hardlink=False) == "copy"
    assert target.read_bytes() == source.read_bytes()
    assert delta._link_or_copy(source, target, digest, hardlink=False) == "copy"


def test_link_or_copy_fails_closed_on_mismatch_and_stale_partial(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.pdbqt"
    source.write_bytes(b"source bytes\n")
    digest = _sha256(source)

    changed = tmp_path / "changed.pdbqt"
    changed.write_bytes(b"different bytes\n")
    with pytest.raises(ValueError, match="refusing to overwrite changed"):
        delta._link_or_copy(source, changed, digest, hardlink=False)

    bad_digest = tmp_path / "bad-digest.pdbqt"
    with pytest.raises(RuntimeError, match="copy checksum mismatch"):
        delta._link_or_copy(source, bad_digest, "0" * 64, hardlink=False)
    assert not bad_digest.exists()
    assert not bad_digest.with_name("bad-digest.pdbqt.part").exists()

    stale = tmp_path / "stale.pdbqt"
    stale.with_name("stale.pdbqt.part").write_bytes(b"partial\n")
    with pytest.raises(ValueError, match="stale FDA delta partial blocks resume"):
        delta._link_or_copy(source, stale, digest, hardlink=False)
    assert not stale.exists()


@pytest.mark.parametrize(
    ("key", "invalid"),
    [
        ("manifest_sha256", "0" * 64),
        ("quarantine_sha256", "0" * 64),
        ("action_counts", {"legacy_byte_copy": 1264, "prepare_parent": 598}),
        ("status_counts", {"ready": 1849, "quarantined_non_dockable": 13}),
        ("canonical_parent_rows", 1862),
        ("ready_pdbqt_count", 1849),
    ],
)
def test_validate_source_summary_requires_exact_hashes_and_counts(
    key: str, invalid: object
) -> None:
    summary, source_hashes = _valid_summary()
    delta._validate_source_summary(summary, source_hashes)

    summary[key] = invalid
    with pytest.raises(
        ValueError, match="finalized FDA summary disagrees with source artifacts"
    ):
        delta._validate_source_summary(summary, source_hashes)


def test_validate_source_row_accepts_exact_ready_and_empty_quarantine(
    tmp_path: Path,
) -> None:
    source = tmp_path / "named-library"
    source.mkdir()
    pdbqt = source / "approved-parent.pdbqt"
    sidecar = source / "approved-parent.ligprep_source.json"
    pdbqt.write_bytes(b"REMARK approved parent\n")
    sidecar.write_bytes(b"{}\n")
    ready = {
        "filename": pdbqt.name,
        "materialization_status": "ready",
        "output_pdbqt_path": str(pdbqt),
        "output_pdbqt_sha256": _sha256(pdbqt),
        "output_pdbqt_bytes": str(pdbqt.stat().st_size),
        "sidecar_path": str(sidecar),
        "sidecar_sha256": _sha256(sidecar),
    }
    delta._validate_source_row(source.resolve(), ready)

    quarantined = {
        "filename": "unsafe-parent.pdbqt",
        "materialization_status": "quarantined_non_dockable",
        "output_pdbqt_sha256": "",
        "sidecar_sha256": "",
    }
    delta._validate_source_row(source.resolve(), quarantined)

    with pytest.raises(ValueError, match="unsafe FDA delta filename"):
        delta._validate_source_row(
            source.resolve(), {**quarantined, "filename": "../escape.pdbqt"}
        )
    with pytest.raises(ValueError, match="quarantined FDA row has output bytes"):
        delta._validate_source_row(
            source.resolve(), {**quarantined, "output_pdbqt_sha256": "f" * 64}
        )
