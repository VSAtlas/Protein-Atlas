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
        str(cfg.get("LIGAND_SOURCE_CACHE_DIR") or output_root(repo_root, "data") / "ligand_sources")
    )
    if not data_root.is_absolute():
        data_root = repo_root / data_root
    prepped_root = Path(
        str(cfg.get("PREPPED_LIGANDS_DIR") or cfg.get("OUTPUT_LIGANDS_DIR") or repo_root / "prepped_ligands")
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
    paths = fetch_source(
        source_name,
        cfg,
        source_sdf=source_sdf,
        source_url=source_url,
        force=force_fetch,
    )
    if force_prep or not library_has_pdbqts(paths.library_dir):
        paths = prepare_source(source_name, cfg, force=force_prep, limit=limit)
    return paths


def _download_to_cache(url: str, raw_dir: Path, *, force: bool) -> Path:
    parsed = urllib.parse.urlparse(url)
    name = Path(parsed.path).name or "downloaded.sdf"
    out = raw_dir / name
    if out.exists() and not force:
        return out
    request = urllib.request.Request(url, headers={"User-Agent": "atlas-ligand-fetch/1.0"})
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
    request = urllib.request.Request(url, headers={"User-Agent": "atlas-ligand-fetch/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8", errors="ignore")


def _materialize_sdf_from_path(source_path: Path, out_sdf: Path) -> None:
    source_path = source_path.expanduser().resolve()
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    if source_path.suffix.lower() == ".zip":
        with zipfile.ZipFile(source_path) as archive:
            sdf_names = [name for name in archive.namelist() if name.lower().endswith(".sdf")]
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


def _download_metadata(source: LigandSource, raw_dir: Path, *, force: bool) -> list[Path]:
    paths: list[Path] = []
    for url in source.metadata_urls:
        paths.append(_download_to_cache(url, raw_dir, force=force))
    return paths


def _filter_drugcentral_fda(paths: LigandLibraryPaths, source: LigandSource, *, force: bool) -> None:
    metadata = _download_metadata(source, paths.raw_dir, force=force)
    approved_tokens = _read_approval_tokens(metadata)
    if not approved_tokens:
        return
    filtered = paths.raw_dir / "fda.filtered.sdf"
    kept = _filter_sdf_records(paths.raw_sdf, filtered, approved_tokens)
    if kept == 0:
        raise RuntimeError(
            "DrugCentral FDA filtering produced zero SDF records. "
            "Use `atlas ligands fetch fda --source-sdf <file>` with a curated SDF if the upstream schema changed."
        )
    filtered.replace(paths.raw_sdf)


def _read_approval_tokens(paths: Sequence[Path]) -> set[str]:
    tokens: set[str] = set()
    preferred = ("drugcentral", "struct", "id", "name", "inn", "cas")
    for path in paths:
        try:
            with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    for key, value in row.items():
                        key_l = str(key or "").lower()
                        if not any(part in key_l for part in preferred):
                            continue
                        value_s = str(value or "").strip()
                        if value_s:
                            tokens.add(value_s.lower())
        except Exception:
            continue
    return tokens


def _filter_sdf_records(in_sdf: Path, out_sdf: Path, tokens: set[str]) -> int:
    kept = 0
    out_sdf.unlink(missing_ok=True)
    for record in _iter_sdf_record_text(in_sdf):
        haystack = record.lower()
        if any(token in haystack for token in tokens):
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
    pdbqts = sorted(paths.library_dir.rglob("*.pdbqt")) if paths.library_dir.exists() else []
    payload = {
        "source": asdict(source),
        "raw_sdf": str(paths.raw_sdf),
        "raw_sdf_sha256": _sha256(paths.raw_sdf) if paths.raw_sdf.exists() else "",
        "library_dir": str(paths.library_dir),
        "prepared": prepared,
        "prepared_count": len(pdbqts),
        "meeko_writer": "meeko",
    }
    paths.source_manifest.parent.mkdir(parents=True, exist_ok=True)
    paths.source_manifest.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
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
]
