from __future__ import annotations

import json
from pathlib import Path

import pytest

from prep_ligands import ligand_library_manager as manager


def _sdf_record(title: str, **properties: str) -> str:
    lines = [
        title,
        "  ATLAS",
        "",
        "  0  0  0  0  0  0  0  0  0  0999 V2000",
        "M  END",
    ]
    for key, value in properties.items():
        lines.extend((f">  <{key}>", value, ""))
    lines.append("$$$$")
    return "\n".join(lines) + "\n"


def _filtered_titles(path: Path) -> list[str]:
    return [record.splitlines()[0] for record in manager._iter_sdf_record_text(path)]


def test_headerless_drugcentral_metadata_keeps_first_row_and_exact_id(
    tmp_path: Path,
) -> None:
    metadata = tmp_path / "FDA_Approved.csv"
    metadata.write_text("21,Alpha Drug\n42,Beta Drug\n", encoding="utf-8")
    source = tmp_path / "structures.sdf"
    source.write_text(
        "".join(
            (
                _sdf_record("Alpha Drug", ID="21", NOTES="approved"),
                # An exact approved name must not override a conflicting DrugCentral ID.
                _sdf_record("Alpha Drug", ID="210", NOTES="contains 21 and Alpha Drug"),
                _sdf_record("Beta Drug", ID="42"),
            )
        ),
        encoding="utf-8",
    )
    output = tmp_path / "filtered.sdf"

    approved = manager._read_drugcentral_approval_keys([metadata])
    kept = manager._filter_sdf_records(source, output, approved)

    assert approved.drugcentral_ids == frozenset({"21", "42"})
    assert approved.names == frozenset({"alpha drug", "beta drug"})
    assert kept == 2
    assert _filtered_titles(output) == ["Alpha Drug", "Beta Drug"]


def test_filter_uses_only_exact_supported_sdf_fields(tmp_path: Path) -> None:
    metadata = tmp_path / "FDA_Approved.csv"
    metadata.write_text(
        "STRUCT_ID,PREFERRED_NAME,APPLICATION_NUMBER\n7,Cat,substring bait\n",
        encoding="utf-8",
    )
    source = tmp_path / "structures.sdf"
    source.write_text(
        "".join(
            (
                _sdf_record("Cat", STRUCT_ID="7"),
                _sdf_record("Catherine", STRUCT_ID="70", NOTES="Cat 7 substring bait"),
                _sdf_record("substring bait", NOTES="7 Cat"),
            )
        ),
        encoding="utf-8",
    )
    output = tmp_path / "filtered.sdf"

    kept = manager._filter_sdf_records(
        source,
        output,
        manager._read_drugcentral_approval_keys([metadata]),
    )

    assert kept == 1
    assert _filtered_titles(output) == ["Cat"]


def test_header_and_property_aliases_allow_exact_name_and_identifier_fallbacks(
    tmp_path: Path,
) -> None:
    metadata = tmp_path / "FDA_Approved.csv"
    metadata.write_text(
        "DrugCentral ID,Generic Name,CAS Registry Number,InChI Key\n"
        "101,Example Drug,123-45-6,ABCDEFGHIJKLMN-OPQRSTUVWX-Y\n",
        encoding="utf-8",
    )
    source = tmp_path / "structures.sdf"
    source.write_text(
        "".join(
            (
                _sdf_record("by id", DRUGCENTRAL_ID="101"),
                _sdf_record("by preferred name", DRUG_NAME="  EXAMPLE   DRUG "),
                _sdf_record("by cas", CAS_RN="123-45-6"),
                _sdf_record("by key", INCHI_KEY="abcdefghijklmn-opqrstuvwx-y"),
                _sdf_record("not exact", DRUG_NAME="Example Drug hydrochloride"),
            )
        ),
        encoding="utf-8",
    )
    output = tmp_path / "filtered.sdf"

    kept = manager._filter_sdf_records(
        source,
        output,
        manager._read_drugcentral_approval_keys([metadata]),
    )

    assert kept == 4
    assert _filtered_titles(output) == [
        "by id",
        "by preferred name",
        "by cas",
        "by key",
    ]


def test_fda_filter_fails_closed_when_approval_metadata_has_no_usable_keys(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    paths = manager.LigandLibraryPaths(
        raw_dir=tmp_path,
        raw_sdf=tmp_path / "fda.sdf",
        library_dir=tmp_path / "library",
        status_log=tmp_path / "status.tsv",
        source_manifest=tmp_path / "manifest.json",
    )
    paths.raw_sdf.write_text(
        _sdf_record("Not necessarily FDA", ID="7"), encoding="utf-8"
    )
    metadata = tmp_path / "FDA_Approved.csv"
    metadata.write_text("\n", encoding="utf-8")
    monkeypatch.setattr(
        manager, "_download_metadata", lambda *_args, **_kwargs: [metadata]
    )

    with pytest.raises(RuntimeError, match="refusing to leave the unfiltered"):
        manager._filter_drugcentral_fda(paths, manager.SOURCES["fda"], force=False)

    assert not paths.raw_sdf.exists()
    quarantined = tmp_path / "fda.unvalidated.sdf"
    assert _filtered_titles(quarantined) == ["Not necessarily FDA"]


def test_cached_fda_sdf_is_revalidated_before_fetch_returns(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    paths = manager.LigandLibraryPaths(
        raw_dir=tmp_path,
        raw_sdf=tmp_path / "fda.sdf",
        library_dir=tmp_path / "library",
        status_log=tmp_path / "status.tsv",
        source_manifest=tmp_path / "manifest.json",
    )
    paths.raw_sdf.write_text(_sdf_record("cached", ID="7"), encoding="utf-8")
    called: list[bool] = []
    monkeypatch.setattr(manager, "paths_for_source", lambda *_args: paths)
    monkeypatch.setattr(
        manager,
        "_filter_drugcentral_fda",
        lambda *_args, **kwargs: called.append(bool(kwargs.get("force"))),
    )

    returned = manager.fetch_source("fda", {}, force=False)

    assert returned == paths
    assert called == [False]


def test_prepare_fda_rejects_cache_without_valid_filter_manifest(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    paths = manager.LigandLibraryPaths(
        raw_dir=tmp_path,
        raw_sdf=tmp_path / "fda.sdf",
        library_dir=tmp_path / "library",
        status_log=tmp_path / "status.tsv",
        source_manifest=tmp_path / "manifest.json",
    )
    paths.raw_sdf.write_text(_sdf_record("unvalidated", ID="7"), encoding="utf-8")
    monkeypatch.setattr(manager, "paths_for_source", lambda *_args: paths)

    with pytest.raises(RuntimeError, match="valid exact-filter manifest"):
        manager.prepare_source("fda", {})

    assert not paths.library_dir.exists()


def test_install_fda_rejects_stale_prepared_library_after_source_repair(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    paths = manager.LigandLibraryPaths(
        raw_dir=tmp_path / "raw",
        raw_sdf=tmp_path / "raw" / "fda.sdf",
        library_dir=tmp_path / "library",
        status_log=tmp_path / "library" / "status.tsv",
        source_manifest=tmp_path / "library" / "atlas_ligand_source_manifest.json",
    )
    paths.raw_dir.mkdir()
    paths.library_dir.mkdir()
    paths.raw_sdf.write_text(_sdf_record("Approved Drug", ID="7"), encoding="utf-8")
    metadata = paths.raw_dir / "FDA_Approved.csv"
    metadata.write_text("7,Approved Drug\n", encoding="utf-8")
    approved = manager._read_drugcentral_approval_keys([metadata])
    manager._write_fda_filter_manifest(paths, [metadata], approved, kept=1)
    (paths.library_dir / "legacy_unfiltered_0001.pdbqt").write_text(
        "TORSDOF 0\n",
        encoding="utf-8",
    )
    paths.source_manifest.write_text(
        json.dumps(
            {
                "source": {"name": "fda"},
                "raw_sdf_sha256": manager._sha256(paths.raw_sdf),
                "prepared": True,
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(manager, "fetch_source", lambda *_args, **_kwargs: paths)
    monkeypatch.setattr(manager, "paths_for_source", lambda *_args: paths)

    with pytest.raises(RuntimeError, match="different or unverifiable source SDF"):
        manager.install_source("fda", {}, force_prep=False)

    assert (paths.library_dir / "legacy_unfiltered_0001.pdbqt").exists()
