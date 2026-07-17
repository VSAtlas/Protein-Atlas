"""Download, prepare, and document Atlas ligand libraries."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

from config.output_paths import output_root
from prep_ligands.library_index import LibraryIndex
from prep_ligands.prep_ligands_bulk_meeko import prep_ligands_with_meeko


@dataclass(frozen=True)
class LigandSource:
    name: str
    library_subdir: str
    description: str
    url: str
    license_note: str
    citation: str
    aliases: tuple[str, ...] = ()
    metadata_urls: tuple[str, ...] = ()
    download_page: str = ""
    download_regex: str = ""


@dataclass(frozen=True)
class LigandLibraryPaths:
    raw_dir: Path
    raw_sdf: Path
    library_dir: Path
    status_log: Path
    source_manifest: Path


SOURCES: dict[str, LigandSource] = {
    "hmdb": LigandSource(
        name="hmdb",
        aliases=("metabolites",),
        library_subdir="hmdb",
        description="Human Metabolome Database metabolite structures.",
        url="https://hmdb.ca/system/downloads/current/structures.zip",
        license_note=(
            "HMDB is freely available for public use; commercial redistribution "
            "requires explicit permission and significant use should cite HMDB."
        ),
        citation="HMDB current release; cite the HMDB database paper.",
    ),
    "fda": LigandSource(
        name="fda",
        aliases=("fda-approved", "drugcentral-fda"),
        library_subdir="fda_library",
        description="FDA-approved small-molecule drugs from DrugCentral structures.",
        url="https://unmtid-dbs.net/download/DrugCentral/2023/structures.molV2.sdf.gz",
        metadata_urls=("https://drugcentral.org/static/FDA_Approved.csv",),
        license_note=(
            "DrugCentral provides public downloads; confirm local citation/license "
            "requirements for publication or redistribution."
        ),
        citation="DrugCentral downloads and FDA-approved-drug list.",
    ),
    "chembl": LigandSource(
        name="chembl",
        aliases=("bioactive", "druglike"),
        library_subdir="chembl",
        description="ChEMBL bioactive drug-like molecule structures from EMBL-EBI.",
        url="https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/",
        download_page="https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/",
        download_regex=r'href="([^"]*chembl_\d+\.sdf\.gz)"',
        license_note="ChEMBL data are provided under Creative Commons Attribution-ShareAlike 3.0.",
        citation="ChEMBL Database in 2023 and the ChEMBL release used for the run.",
    ),
    "chebi": LigandSource(
        name="chebi",
        aliases=("chebi-3star", "chebi-lite"),
        library_subdir="chebi",
        description="High-confidence ChEBI 3-star chemical structures from EMBL-EBI.",
        url="https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi_lite_3_stars.sdf.gz",
        license_note="ChEBI is distributed by EMBL-EBI; cite the ChEBI database and release used.",
        citation="ChEBI SDF downloads from EMBL-EBI.",
    ),
    "coconut": LigandSource(
        name="coconut",
        aliases=("natural-products", "natural_products", "np"),
        library_subdir="coconut",
        description="COCONUT open natural products 3D SDF structures.",
        url="https://coconut.naturalproducts.net/download",
        download_page="https://coconut.naturalproducts.net/download",
        download_regex=r'href="(https://coconut\.s3\.uni-jena\.de/[^"]*coconut_sdf_3d-[^"]*\.zip)"',
        license_note="COCONUT data are released under Creative Commons CC0.",
        citation="COCONUT 2.0 and the COCONUT release used for the run.",
    ),
}


def known_sources() -> list[LigandSource]:
    return [SOURCES[key] for key in sorted(SOURCES)]


def resolve_source(name: str) -> LigandSource:
    token = str(name).strip().lower()
    for source in SOURCES.values():
        if token == source.name or token in source.aliases:
            return source
    choices = ", ".join(sorted(SOURCES))
    raise ValueError(f"unknown ligand source '{name}'. Known sources: {choices}")


def paths_for_source(cfg: dict[str, Any], source: LigandSource) -> LigandLibraryPaths:
    repo_root = Path(__file__).resolve().parents[2]
    data_root = Path(
        str(
            cfg.get("LIGAND_SOURCE_CACHE_DIR")
            or output_root(repo_root, "data") / "ligand_sources"
        )
    )
    if not data_root.is_absolute():
        data_root = repo_root / data_root
    prepped_root = Path(
        str(
            cfg.get("PREPPED_LIGANDS_DIR")
            or cfg.get("OUTPUT_LIGANDS_DIR")
            or repo_root / "prepped_ligands"
        )
    )
    if not prepped_root.is_absolute():
        prepped_root = repo_root / prepped_root
    raw_dir = data_root / source.name
    library_dir = prepped_root / source.library_subdir
    return LigandLibraryPaths(
        raw_dir=raw_dir,
        raw_sdf=raw_dir / f"{source.name}.sdf",
        library_dir=library_dir,
        status_log=library_dir / "ligand_prep_status.tsv",
        source_manifest=library_dir / "atlas_ligand_source_manifest.json",
    )


def library_has_pdbqts(library_dir: Path) -> bool:
    try:
        return any(library_dir.rglob("*.pdbqt"))
    except Exception:
        return False


def fetch_source(
    source_name: str,
    cfg: dict[str, Any],
    *,
    source_sdf: Path | None = None,
    source_url: str | None = None,
    force: bool = False,
    limit: int | None = None,
) -> LigandLibraryPaths:
    source = resolve_source(source_name)
    paths = paths_for_source(cfg, source)
    paths.raw_dir.mkdir(parents=True, exist_ok=True)

    if paths.raw_sdf.exists() and not force:
        if source.name == "fda":
            _filter_drugcentral_fda(paths, source, force=False)
        return paths

    if source_sdf is not None:
        _materialize_sdf_from_path(Path(source_sdf), paths.raw_sdf)
    else:
        resolved_url = source_url or _resolve_download_url(source)
        downloaded = _download_to_cache(resolved_url, paths.raw_dir, force=force)
        _materialize_sdf_from_path(downloaded, paths.raw_sdf)

    if source.name == "fda":
        _filter_drugcentral_fda(paths, source, force=force)

    if limit is not None and limit > 0:
        _limit_sdf_records(paths.raw_sdf, limit)
        if source.name == "fda":
            _refresh_fda_filter_manifest(paths)
    return paths


def prepare_source(
    source_name: str,
    cfg: dict[str, Any],
    *,
    force: bool = False,
    limit: int | None = None,
) -> LigandLibraryPaths:
    source = resolve_source(source_name)
    paths = paths_for_source(cfg, source)
    if not paths.raw_sdf.exists():
        raise FileNotFoundError(
            f"source SDF not found: {paths.raw_sdf}. Run `atlas ligands fetch {source.name}` first."
        )
    if source.name == "fda":
        valid, reason, _ = validate_fda_filter_manifest(paths.raw_sdf)
        if not valid:
            raise RuntimeError(
                "FDA source cache is not backed by a valid exact-filter manifest "
                f"({reason}); rerun `atlas ligands fetch fda` before preparation."
            )
        _guard_fda_prepared_library(paths, source, force=force)
    if limit is not None and limit > 0:
        limited = paths.raw_dir / f"{source.name}.limit{limit}.sdf"
        shutil.copy2(paths.raw_sdf, limited)
        _limit_sdf_records(limited, limit)
        in_sdf = limited
    else:
        in_sdf = paths.raw_sdf

    paths.library_dir.mkdir(parents=True, exist_ok=True)
    prep_ligands_with_meeko(
        force=force,
        in_sdf=in_sdf,
        mol2_dir=paths.raw_dir / "intermediates",
        out_pdbqt_dir=paths.library_dir,
        status_log=paths.status_log,
    )
    _write_library_manifest(paths.library_dir)
    _write_source_manifest(source, paths, prepared=True)
    return paths


def install_source(
    source_name: str,
    cfg: dict[str, Any],
    *,
    source_sdf: Path | None = None,
    source_url: str | None = None,
    force_fetch: bool = False,
    force_prep: bool = False,
    limit: int | None = None,
) -> LigandLibraryPaths:
    source = resolve_source(source_name)
    paths = fetch_source(
        source_name,
        cfg,
        source_sdf=source_sdf,
        source_url=source_url,
        force=force_fetch,
    )
    if force_prep or not _prepared_library_is_current(paths, source):
        paths = prepare_source(source_name, cfg, force=force_prep, limit=limit)
    return paths


def _guard_fda_prepared_library(
    paths: LigandLibraryPaths,
    source: LigandSource,
    *,
    force: bool,
) -> None:
    if not library_has_pdbqts(paths.library_dir):
        return
    if _prepared_library_is_current(paths, source):
        return
    if not force:
        raise RuntimeError(
            "Existing FDA PDBQT library was prepared from a different or "
            "unverifiable source SDF. Rerun with --force-prep to preserve it "
            "under an .unvalidated directory and rebuild from the filtered source."
        )
    _quarantine_unvalidated_library(paths.library_dir)


def _quarantine_unvalidated_library(library_dir: Path) -> Path:
    candidate = library_dir.with_name(f"{library_dir.name}.unvalidated")
    suffix = 1
    while candidate.exists():
        candidate = library_dir.with_name(f"{library_dir.name}.unvalidated.{suffix}")
        suffix += 1
    library_dir.replace(candidate)
    return candidate


def _prepared_library_is_current(
    paths: LigandLibraryPaths,
    source: LigandSource,
) -> bool:
    if not library_has_pdbqts(paths.library_dir):
        return False
    payload = _read_source_manifest(paths.source_manifest)
    if payload is None or not payload.get("prepared"):
        return False
    if _source_manifest_name(payload) != source.name:
        return False
    if _clean_manifest_text(payload.get("raw_sdf_sha256")) != _sha256(paths.raw_sdf):
        return False
    if source.name != "fda":
        return True
    return _prepared_fda_manifest_matches(paths, payload)


def _read_source_manifest(path: Path) -> dict[str, object] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _source_manifest_name(payload: dict[str, object]) -> str:
    source = payload.get("source")
    if not isinstance(source, dict):
        return ""
    return _clean_manifest_text(source.get("name"))


def _prepared_fda_manifest_matches(
    paths: LigandLibraryPaths,
    payload: dict[str, object],
) -> bool:
    valid, _, manifest = validate_fda_filter_manifest(paths.raw_sdf)
    if not valid or manifest is None:
        return False
    stored_hash = _clean_manifest_text(payload.get("fda_filter_manifest_sha256"))
    return bool(stored_hash and stored_hash == _sha256(manifest))


def _download_to_cache(url: str, raw_dir: Path, *, force: bool) -> Path:
    parsed = urllib.parse.urlparse(url)
    name = Path(parsed.path).name or "downloaded.sdf"
    out = raw_dir / name
    if out.exists() and not force:
        return out
    request = urllib.request.Request(
        url, headers={"User-Agent": "atlas-ligand-fetch/1.0"}
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            with out.open("wb") as handle:
                shutil.copyfileobj(response, handle)
    except urllib.error.HTTPError as exc:
        if exc.code == 403:
            raise RuntimeError(
                f"download forbidden by upstream source: {url}. "
                "Some ligand providers require an interactive browser/form download. "
                "Download the SDF/SDF.GZ/ZIP manually, then rerun with "
                "`atlas ligands install <source> --source-sdf <path>`."
            ) from exc
        raise
    return out


def _resolve_download_url(source: LigandSource) -> str:
    if not source.download_page or not source.download_regex:
        return source.url
    page = _read_url_text(source.download_page)
    matches = re.findall(source.download_regex, page, flags=re.IGNORECASE)
    if not matches:
        raise RuntimeError(
            f"could not discover current download URL for {source.name} from {source.download_page}. "
            "Use `atlas ligands install <source> --source-url <url>` or "
            "`atlas ligands install <source> --source-sdf <path>`."
        )
    return urllib.parse.urljoin(source.download_page, matches[0])


def _read_url_text(url: str) -> str:
    request = urllib.request.Request(
        url, headers={"User-Agent": "atlas-ligand-fetch/1.0"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8", errors="ignore")


def _materialize_sdf_from_path(source_path: Path, out_sdf: Path) -> None:
    source_path = source_path.expanduser().resolve()
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    if source_path.suffix.lower() == ".zip":
        with zipfile.ZipFile(source_path) as archive:
            sdf_names = [
                name for name in archive.namelist() if name.lower().endswith(".sdf")
            ]
            if not sdf_names:
                raise ValueError(f"no .sdf file found in archive: {source_path}")
            with archive.open(sdf_names[0]) as src, out_sdf.open("wb") as dst:
                shutil.copyfileobj(src, dst)
        return
    if source_path.suffix.lower() == ".gz":
        with gzip.open(source_path, "rb") as src, out_sdf.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        return
    if source_path.resolve() != out_sdf.resolve():
        shutil.copy2(source_path, out_sdf)


def _download_metadata(
    source: LigandSource, raw_dir: Path, *, force: bool
) -> list[Path]:
    paths: list[Path] = []
    for url in source.metadata_urls:
        paths.append(_download_to_cache(url, raw_dir, force=force))
    return paths


def _filter_drugcentral_fda(
    paths: LigandLibraryPaths, source: LigandSource, *, force: bool
) -> None:
    if not force:
        valid, _, _ = validate_fda_filter_manifest(paths.raw_sdf)
        if valid:
            return
    metadata = _download_metadata(source, paths.raw_dir, force=force)
    approved_keys = _read_drugcentral_approval_keys(metadata)
    if not approved_keys:
        _quarantine_unvalidated_fda_sdf(paths)
        raise RuntimeError(
            "DrugCentral FDA approval metadata contained no usable exact identifiers "
            "or names; refusing to leave the unfiltered all-drug SDF labeled as FDA."
        )
    filtered = paths.raw_dir / "fda.filtered.sdf"
    kept = _filter_sdf_records(paths.raw_sdf, filtered, approved_keys)
    if kept == 0:
        filtered.unlink(missing_ok=True)
        _quarantine_unvalidated_fda_sdf(paths)
        raise RuntimeError(
            "DrugCentral FDA filtering produced zero SDF records. "
            "Use `atlas ligands fetch fda --source-sdf <file>` with a curated SDF if the upstream schema changed."
        )
    filtered.replace(paths.raw_sdf)
    _write_fda_filter_manifest(paths, metadata, approved_keys, kept)


def _quarantine_unvalidated_fda_sdf(paths: LigandLibraryPaths) -> Path | None:
    if not paths.raw_sdf.exists():
        return None
    candidate = paths.raw_dir / "fda.unvalidated.sdf"
    suffix = 1
    while candidate.exists():
        candidate = paths.raw_dir / f"fda.unvalidated.{suffix}.sdf"
        suffix += 1
    paths.raw_sdf.replace(candidate)
    return candidate


@dataclass(frozen=True)
class _DrugCentralApprovalKeys:
    drugcentral_ids: frozenset[str]
    names: frozenset[str]
    cas_numbers: frozenset[str]
    inchi_keys: frozenset[str]

    def __bool__(self) -> bool:
        return any(
            (self.drugcentral_ids, self.names, self.cas_numbers, self.inchi_keys)
        )


def _write_fda_filter_manifest(
    paths: LigandLibraryPaths,
    metadata: Sequence[Path],
    approved_keys: _DrugCentralApprovalKeys,
    kept: int,
) -> Path:
    manifest = paths.raw_dir / "fda_filter_manifest.json"
    payload = {
        "schema_version": 1,
        "filter_version": "drugcentral-fda-exact-v1",
        "filtered_sdf": str(paths.raw_sdf),
        "filtered_sdf_sha256": _sha256(paths.raw_sdf),
        "kept_records": kept,
        "approval_metadata": [
            {"path": str(path), "sha256": _sha256(path)} for path in metadata
        ],
        "approval_key_counts": {
            "drugcentral_ids": len(approved_keys.drugcentral_ids),
            "names": len(approved_keys.names),
            "cas_numbers": len(approved_keys.cas_numbers),
            "inchi_keys": len(approved_keys.inchi_keys),
        },
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def _refresh_fda_filter_manifest(paths: LigandLibraryPaths) -> None:
    manifest = paths.raw_dir / "fda_filter_manifest.json"
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        raise RuntimeError(
            "FDA source was filtered but its provenance manifest could not be refreshed."
        )
    if not isinstance(payload, dict):
        raise RuntimeError("FDA filter provenance manifest is not a JSON object.")
    payload["filtered_sdf"] = str(paths.raw_sdf)
    payload["filtered_sdf_sha256"] = _sha256(paths.raw_sdf)
    payload["kept_records"] = sum(1 for _ in _iter_sdf_record_text(paths.raw_sdf))
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def validate_fda_filter_manifest(
    source_sdf: Path,
    manifest_path: Path | None = None,
) -> tuple[bool, str, Path | None]:
    source = Path(source_sdf).expanduser()
    manifest = _resolve_fda_manifest_path(source, manifest_path)
    if not source.is_file():
        return False, "source_sdf_missing", manifest
    if not manifest.is_file():
        return False, "approval_manifest_missing", manifest
    payload, load_reason = _load_fda_filter_manifest(manifest)
    if payload is None:
        return False, load_reason, manifest
    for validator in (
        _manifest_version_reason,
        lambda value: _manifest_sdf_reason(value, source),
        lambda value: _manifest_metadata_reason(value, manifest),
    ):
        reason = validator(payload)
        if reason:
            return False, reason, manifest
    return True, "", manifest.resolve()


def _resolve_fda_manifest_path(
    source: Path,
    manifest_path: Path | None,
) -> Path:
    if manifest_path is not None:
        return Path(manifest_path).expanduser()
    return source.parent / "fda_filter_manifest.json"


def _load_fda_filter_manifest(
    manifest: Path,
) -> tuple[dict[str, object] | None, str]:
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None, "approval_manifest_invalid"
    if not isinstance(payload, dict):
        return None, "approval_manifest_invalid"
    return payload, ""


def _manifest_version_reason(payload: dict[str, object]) -> str:
    if (
        payload.get("schema_version") != 1
        or payload.get("filter_version") != "drugcentral-fda-exact-v1"
    ):
        return "approval_manifest_version_unsupported"
    return ""


def _manifest_sdf_reason(payload: dict[str, object], source: Path) -> str:
    if _clean_manifest_text(payload.get("filtered_sdf_sha256")) != _sha256(source):
        return "approval_manifest_sdf_hash_mismatch"
    expected_records = _manifest_record_count(payload)
    actual_records = sum(1 for _ in _iter_sdf_record_text(source))
    if expected_records <= 0 or actual_records != expected_records:
        return "approval_manifest_record_count_mismatch"
    key_counts = payload.get("approval_key_counts")
    if not isinstance(key_counts, dict) or not any(
        isinstance(value, int) and value > 0 for value in key_counts.values()
    ):
        return "approval_manifest_has_no_approval_keys"
    return ""


def _manifest_record_count(payload: dict[str, object]) -> int:
    try:
        return int(str(payload.get("kept_records", 0)))
    except (TypeError, ValueError):
        return 0


def _manifest_metadata_reason(payload: dict[str, object], manifest: Path) -> str:
    metadata = payload.get("approval_metadata")
    if not isinstance(metadata, list) or not metadata:
        return "approval_manifest_metadata_missing"
    for item in metadata:
        if not isinstance(item, dict):
            return "approval_manifest_metadata_invalid"
        if not _manifest_metadata_item_valid(item, manifest):
            return "approval_manifest_metadata_hash_mismatch"
    return ""


def _manifest_metadata_item_valid(
    item: dict[str, object],
    manifest: Path,
) -> bool:
    metadata_path = Path(_clean_manifest_text(item.get("path"))).expanduser()
    if not metadata_path.is_absolute():
        metadata_path = manifest.parent / metadata_path
    expected_hash = _clean_manifest_text(item.get("sha256"))
    if not metadata_path.is_file() or not expected_hash:
        return False
    return _sha256(metadata_path) == expected_hash


def _clean_manifest_text(value: object) -> str:
    return str(value or "").strip()


_APPROVAL_HEADER_KIND = {
    "id": "drugcentral_ids",
    "dcid": "drugcentral_ids",
    "drugid": "drugcentral_ids",
    "structid": "drugcentral_ids",
    "structureid": "drugcentral_ids",
    "drugcentralid": "drugcentral_ids",
    "drugcentralstructid": "drugcentral_ids",
    "name": "names",
    "inn": "names",
    "drugname": "names",
    "genericname": "names",
    "preferredname": "names",
    "cas": "cas_numbers",
    "casrn": "cas_numbers",
    "casnumber": "cas_numbers",
    "casregistrynumber": "cas_numbers",
    "inchikey": "inchi_keys",
}

_SDF_PROPERTY_KIND = {
    "id": "drugcentral_ids",
    "dcid": "drugcentral_ids",
    "drugid": "drugcentral_ids",
    "structid": "drugcentral_ids",
    "structureid": "drugcentral_ids",
    "drugcentralid": "drugcentral_ids",
    "drugcentralstructid": "drugcentral_ids",
    "name": "names",
    "inn": "names",
    "drugname": "names",
    "genericname": "names",
    "preferredname": "names",
    "cas": "cas_numbers",
    "casrn": "cas_numbers",
    "casnumber": "cas_numbers",
    "casregistrynumber": "cas_numbers",
    "inchikey": "inchi_keys",
}

_SDF_PROPERTY_RE = re.compile(r"^>\s*<([^>]+)>")


def _normalized_field_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lstrip("\ufeff").casefold())


def _normalized_exact_value(value: str) -> str:
    return " ".join(value.strip().casefold().split())


def _read_drugcentral_approval_keys(paths: Sequence[Path]) -> _DrugCentralApprovalKeys:
    values: dict[str, set[str]] = {
        "drugcentral_ids": set(),
        "names": set(),
        "cas_numbers": set(),
        "inchi_keys": set(),
    }

    def add_row(row: Sequence[str], columns: dict[str, int]) -> None:
        for kind, column in columns.items():
            if column >= len(row):
                continue
            value = _normalized_exact_value(str(row[column] or ""))
            if value:
                values[kind].add(value)

    for path in paths:
        try:
            with path.open(
                "r", encoding="utf-8", errors="ignore", newline=""
            ) as handle:
                reader = csv.reader(handle)
                first_row = next(reader, None)
                if first_row is None:
                    continue
                columns: dict[str, int] = {}
                for column, header in enumerate(first_row):
                    kind = _APPROVAL_HEADER_KIND.get(
                        _normalized_field_name(str(header or ""))
                    )
                    if kind is not None and kind not in columns:
                        columns[kind] = column
                if not columns:
                    # The current DrugCentral FDA_Approved.csv has no header and
                    # consists of ``DrugCentral ID,preferred name`` rows.
                    columns = {"drugcentral_ids": 0, "names": 1}
                    add_row(first_row, columns)
                for row in reader:
                    add_row(row, columns)
        except Exception:
            continue
    return _DrugCentralApprovalKeys(
        drugcentral_ids=frozenset(values["drugcentral_ids"]),
        names=frozenset(values["names"]),
        cas_numbers=frozenset(values["cas_numbers"]),
        inchi_keys=frozenset(values["inchi_keys"]),
    )


def _sdf_record_keys(record: str) -> _DrugCentralApprovalKeys:
    values: dict[str, set[str]] = {
        "drugcentral_ids": set(),
        "names": set(),
        "cas_numbers": set(),
        "inchi_keys": set(),
    }
    lines = record.splitlines()
    if lines:
        title = _normalized_exact_value(lines[0])
        if title:
            values["names"].add(title)

    line_number = 0
    while line_number < len(lines):
        match = _SDF_PROPERTY_RE.match(lines[line_number])
        if match is None:
            line_number += 1
            continue
        kind = _SDF_PROPERTY_KIND.get(_normalized_field_name(match.group(1)))
        line_number += 1
        property_lines: list[str] = []
        while line_number < len(lines) and lines[line_number].strip():
            if lines[line_number].startswith("$$$$"):
                break
            property_lines.append(lines[line_number])
            line_number += 1
        if kind is not None:
            value = _normalized_exact_value(" ".join(property_lines))
            if value:
                values[kind].add(value)

    return _DrugCentralApprovalKeys(
        drugcentral_ids=frozenset(values["drugcentral_ids"]),
        names=frozenset(values["names"]),
        cas_numbers=frozenset(values["cas_numbers"]),
        inchi_keys=frozenset(values["inchi_keys"]),
    )


def _record_is_approved(
    record_keys: _DrugCentralApprovalKeys,
    approved_keys: _DrugCentralApprovalKeys,
) -> bool:
    # DrugCentral's own structure ID is authoritative when both inputs expose it.
    # Do not allow an unrelated record with a coincidentally matching name to
    # override a conflicting exact ID.
    if record_keys.drugcentral_ids and approved_keys.drugcentral_ids:
        return bool(record_keys.drugcentral_ids & approved_keys.drugcentral_ids)
    for record_values, approved_values in (
        (record_keys.inchi_keys, approved_keys.inchi_keys),
        (record_keys.cas_numbers, approved_keys.cas_numbers),
        (record_keys.names, approved_keys.names),
    ):
        if record_values and approved_values and record_values & approved_values:
            return True
    return False


def _filter_sdf_records(
    in_sdf: Path,
    out_sdf: Path,
    approved_keys: _DrugCentralApprovalKeys,
) -> int:
    kept = 0
    out_sdf.unlink(missing_ok=True)
    for record in _iter_sdf_record_text(in_sdf):
        if _record_is_approved(_sdf_record_keys(record), approved_keys):
            with out_sdf.open("a", encoding="utf-8") as handle:
                handle.write(record)
                if not record.endswith("\n"):
                    handle.write("\n")
            kept += 1
    return kept


def _limit_sdf_records(path: Path, limit: int) -> None:
    tmp = path.with_suffix(path.suffix + ".limit.tmp")
    kept = 0
    with tmp.open("w", encoding="utf-8") as handle:
        for record in _iter_sdf_record_text(path):
            if kept >= limit:
                break
            handle.write(record)
            if not record.endswith("\n"):
                handle.write("\n")
            kept += 1
    tmp.replace(path)


def _iter_sdf_record_text(path: Path) -> Iterable[str]:
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            chunks.append(line)
            if line.startswith("$$$$"):
                yield "".join(chunks)
                chunks = []
    if chunks:
        yield "".join(chunks)


def _write_library_manifest(library_dir: Path) -> Path:
    rel_entries = [
        path.relative_to(library_dir)
        for path in sorted(library_dir.rglob("*.pdbqt"))
        if path.is_file()
    ]
    return LibraryIndex().write_manifest_for_root(library_dir, rel_entries)


def _write_source_manifest(
    source: LigandSource,
    paths: LigandLibraryPaths,
    *,
    prepared: bool,
) -> Path:
    pdbqts = (
        sorted(paths.library_dir.rglob("*.pdbqt")) if paths.library_dir.exists() else []
    )
    payload = {
        "source": asdict(source),
        "raw_sdf": str(paths.raw_sdf),
        "raw_sdf_sha256": _sha256(paths.raw_sdf) if paths.raw_sdf.exists() else "",
        "library_dir": str(paths.library_dir),
        "prepared": prepared,
        "prepared_count": len(pdbqts),
        "meeko_writer": "meeko",
    }
    if source.name == "fda":
        filter_manifest = paths.raw_dir / "fda_filter_manifest.json"
        payload["fda_filter_manifest"] = str(filter_manifest)
        payload["fda_filter_manifest_sha256"] = _sha256(filter_manifest)
    paths.source_manifest.parent.mkdir(parents=True, exist_ok=True)
    paths.source_manifest.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    return paths.source_manifest


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = [
    "LigandLibraryPaths",
    "LigandSource",
    "fetch_source",
    "install_source",
    "known_sources",
    "library_has_pdbqts",
    "paths_for_source",
    "prepare_source",
    "resolve_source",
    "validate_fda_filter_manifest",
]
