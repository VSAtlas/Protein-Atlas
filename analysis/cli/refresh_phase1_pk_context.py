from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import zipfile

import pandas as pd

from analysis.external.openfda_pk import fetch_openfda_pk_context
from analysis.external.openfda_pk_review import audit_openfda_pk_cache
from analysis.external.pk_context import (
    combine_pk_context,
    load_existing_phase1_pk_context,
    load_flat_pk_context,
    load_spd_pk_context,
    write_pk_context_outputs,
)
from analysis.external.pkdb_api import probe_pkdb_api
from analysis.external.source_tables import download_to_cache


VERSION = "Atlasv0.0.03"
NCATS_FRDB_URL = "https://drugs.ncats.io/downloads-public/frdb-v2024-12-30.zip"
NCATS_FRDB_SOURCES = (
    ("2024-12-30", NCATS_FRDB_URL),
    ("2023-07-05", "https://drugs.ncats.io/downloads-public/frdb-v2023-07-05.zip"),
    ("2023-02-15", "https://drugs.ncats.io/downloads-public/frdb-v2023-02-15.zip"),
    ("2021-09-09", "https://drugs.ncats.io/downloads-public/frdb-v2021-09-09.zip"),
    ("2021-05-10", "https://drugs.ncats.io/downloads-public/frdb-v2021-05-10.zip"),
)
PKDB_BULK_URL = "https://pk-db.com/api/v1/filter/?download=true&concise=false"
DEFAULT_MODEL_TABLE = Path(
    "data/AtlasSPD_phase1/combined_activity_source_matched_20260709/"
    "model_ready/spd_binding_deduplicated.csv"
)
DEFAULT_SPD = Path(
    "data/external/spd/sutherland_2023_spd_supplementary_data_1_15.xlsx"
)


def _extract_archive(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as handle:
        handle.extractall(destination)


def _add_flat_source(
    parts: list[pd.DataFrame],
    statuses: list[dict[str, object]],
    *,
    path: Path,
    source_name: str,
    version: str,
) -> None:
    if not path.exists() or path.stat().st_size <= 3:
        statuses.append(
            {
                "source": source_name,
                "status": "unavailable",
                "path": str(path),
                "rows": 0,
                "reason": "file missing or empty",
            }
        )
        return
    try:
        context = load_flat_pk_context(path, source_name=source_name, source_version=version)
    except Exception as exc:
        statuses.append(
            {
                "source": source_name,
                "status": "parse_failed",
                "path": str(path),
                "rows": 0,
                "reason": f"{type(exc).__name__}: {exc}",
            }
        )
        return
    parts.append(context)
    statuses.append(
        {
            "source": source_name,
            "status": "ingested",
            "path": str(path),
            "rows": int(len(context)),
            "reason": "",
        }
    )

def _prepare_ncats(
    *,
    external_root: Path,
    download_sources: bool,
    parts: list[pd.DataFrame],
    statuses: list[dict[str, object]],
) -> None:
    source_dir = external_root / "ncats_inxight"
    archive = source_dir / "frdb-v2024-12-30.zip"
    extracted = source_dir / "frdb-v2024-12-30"
    selected_version = "2024-12-30"
    download_errors: list[str] = []
    archive_valid = archive.exists() and zipfile.is_zipfile(archive)
    if download_sources and not archive_valid:
        try:
            download_to_cache(
                NCATS_FRDB_URL,
                archive,
                retries=3,
                sleep_sec=5.0,
                overwrite=True,
            )
            archive_valid = zipfile.is_zipfile(archive)
        except Exception as exc:
            download_errors.append(f"2024-12-30: {type(exc).__name__}: {exc}")
    if download_sources and not archive_valid:
        for version, url in NCATS_FRDB_SOURCES[1:]:
            candidate = source_dir / f"frdb-v{version}.zip"
            try:
                download_to_cache(
                    url,
                    candidate,
                    retries=2,
                    sleep_sec=5.0,
                    overwrite=True,
                )
            except Exception as exc:
                download_errors.append(f"{version}: {type(exc).__name__}: {exc}")
                continue
            if zipfile.is_zipfile(candidate):
                archive = candidate
                extracted = source_dir / f"frdb-v{version}"
                selected_version = version
                archive_valid = True
                break
    if archive_valid and not extracted.exists():
        try:
            _extract_archive(archive, extracted)
        except Exception as exc:
            statuses.append(
                {
                    "source": "NCATS_Inxight_FRDB",
                    "status": "archive_failed",
                    "path": str(archive),
                    "rows": 0,
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )
    candidates = list(extracted.glob("**/frdb-pk.tsv")) if extracted.exists() else []
    if candidates:
        _add_flat_source(
            parts,
            statuses,
            path=candidates[0],
            source_name="NCATS_Inxight_FRDB",
            version=selected_version,
        )
    elif not any(row["source"] == "NCATS_Inxight_FRDB" for row in statuses):
        statuses.append(
            {
                "source": "NCATS_Inxight_FRDB",
                "status": "unavailable",
                "path": str(archive),
                "rows": 0,
                "reason": (
                    "all official current/archive downloads were unavailable; "
                    + (
                        " | ".join(download_errors)
                        if download_errors
                        else "network retry skipped for this refresh"
                    )
                ),
            }
        )


def _prepare_pkdb(
    *,
    external_root: Path,
    download_sources: bool,
    parts: list[pd.DataFrame],
    statuses: list[dict[str, object]],
) -> None:
    source_dir = external_root / "pkdb"
    health_path = source_dir / "pkdb_api_health.json"
    if download_sources:
        health = probe_pkdb_api(out_dir=source_dir)
    elif health_path.exists():
        health = json.loads(health_path.read_text(encoding="utf-8"))
    else:
        health = {
            "status": "not_checked",
            "reason": "network checks disabled by --skip-download",
        }
    if health.get("status") != "available":
        statuses.append(
            {
                "source": "PK-DB",
                "status": str(health.get("status") or "unavailable"),
                "path": str(source_dir / "pkdb_api_health.json"),
                "rows": 0,
                "reason": str(health.get("reason") or "PK-DB API unavailable"),
                "advertised_outputs": health.get("filter_advertised_outputs", 0),
                "retrieved_outputs": health.get("outputs_endpoint_count", 0),
                "archive_outputs_csv_bytes": health.get(
                    "archive_outputs_csv_bytes", 0
                ),
            }
        )
        return
    extracted = source_dir / "pkdb_all_full"
    archive = source_dir / "pkdb_all_full.zip"
    output_path = extracted / "outputs.csv"
    usable_output = output_path.exists() and output_path.stat().st_size > 3
    if download_sources and not usable_output:
        try:
            download_to_cache(
                PKDB_BULK_URL,
                archive,
                retries=2,
                sleep_sec=3.0,
                overwrite=True,
            )
        except Exception as exc:
            statuses.append(
                {
                    "source": "PK-DB",
                    "status": "download_failed",
                    "path": str(archive),
                    "rows": 0,
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )
            return
    if archive.exists() and zipfile.is_zipfile(archive) and not usable_output:
        try:
            _extract_archive(archive, extracted)
        except Exception as exc:
            statuses.append(
                {
                    "source": "PK-DB",
                    "status": "archive_failed",
                    "path": str(archive),
                    "rows": 0,
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )
            return
    if output_path.exists() and output_path.stat().st_size > 3:
        _add_flat_source(
            parts,
            statuses,
            path=output_path,
            source_name="PK-DB",
            version="live_api_export",
        )
    else:
        statuses.append(
            {
                "source": "PK-DB",
                "status": "api_export_incomplete",
                "path": str(archive),
                "rows": 0,
                "reason": (
                    "PK-DB health probe passed but bulk outputs.csv remained empty; "
                    "do not interpret as zero PK coverage"
                ),
            }
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Refresh contextual PK sources and join them to AtlasSPD Phase 1."
    )
    parser.add_argument("--model-table", type=Path, default=DEFAULT_MODEL_TABLE)
    parser.add_argument("--spd-workbook", type=Path, default=DEFAULT_SPD)
    parser.add_argument("--external-root", type=Path, default=Path("data/external"))
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/AtlasSPD_phase1/pk_context_v0_0_03"),
    )
    parser.add_argument("--openfda", choices=("auto", "always", "never"), default="auto")
    parser.add_argument("--openfda-max-drugs", type=int, default=0)
    parser.add_argument("--openfda-sleep-sec", type=float, default=0.25)
    parser.add_argument("--skip-openfda-source-review", action="store_true")
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--drugbank-cmax", type=Path, default=None)
    parser.add_argument("--drugbank-protein-binding", type=Path, default=None)
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    )
    if not args.model_table.exists():
        raise FileNotFoundError(f"Phase 1 model table not found: {args.model_table}")
    model_table = pd.read_csv(args.model_table, low_memory=False)
    parts: list[pd.DataFrame] = []
    statuses: list[dict[str, object]] = []

    existing = load_existing_phase1_pk_context(model_table)
    parts.append(existing)
    statuses.append(
        {
            "source": "existing_phase1",
            "status": "ingested",
            "path": str(args.model_table),
            "rows": int(len(existing)),
            "reason": "preserved before external backfill",
        }
    )

    if args.spd_workbook.exists():
        spd = load_spd_pk_context(args.spd_workbook)
        parts.append(spd)
        statuses.append(
            {
                "source": "SPD",
                "status": "ingested",
                "path": str(args.spd_workbook),
                "rows": int(len(spd)),
                "reason": "",
            }
        )
    else:
        statuses.append(
            {
                "source": "SPD",
                "status": "unavailable",
                "path": str(args.spd_workbook),
                "rows": 0,
                "reason": "workbook missing",
            }
        )

    _prepare_ncats(
        external_root=args.external_root,
        download_sources=not args.skip_download,
        parts=parts,
        statuses=statuses,
    )
    _prepare_pkdb(
        external_root=args.external_root,
        download_sources=not args.skip_download,
        parts=parts,
        statuses=statuses,
    )

    for source_name, path in (
        ("DrugBank_Cmax", args.drugbank_cmax),
        ("DrugBank_protein_binding", args.drugbank_protein_binding),
    ):
        if path is not None:
            _add_flat_source(
                parts,
                statuses,
                path=path,
                source_name=source_name,
                version="BYOL",
            )
        else:
            statuses.append(
                {
                    "source": source_name,
                    "status": "byol_required",
                    "path": "",
                    "rows": 0,
                    "reason": "licensed DrugBank export not provided",
                }
            )

    openfda_cache = args.external_root / "dailymed_spl" / "phase1_openfda"
    openfda_review_manifest: dict[str, object] = {}
    if args.openfda != "never":
        openfda_context, openfda_manifest = fetch_openfda_pk_context(
            model_table,
            cache_dir=openfda_cache / "records",
            max_drugs=max(0, args.openfda_max_drugs),
            sleep_sec=max(0.0, args.openfda_sleep_sec),
            reuse_cache=args.openfda == "auto",
        )
        if not openfda_context.empty:
            parts.append(openfda_context)
        statuses.append(
            {
                "source": "DailyMed_openFDA_SPL",
                "status": "ingested" if not openfda_context.empty else "no_numeric_rows",
                "path": str(openfda_cache),
                "rows": int(len(openfda_context)),
                "reason": "" if not openfda_context.empty else "no label context rows matched",
                **openfda_manifest,
            }
        )
        if not args.skip_openfda_source_review:
            openfda_review_manifest = audit_openfda_pk_cache(
                cache_dir=openfda_cache / "records",
                model_table=model_table,
                out_dir=args.out_dir / "openfda_source_text_review",
            )
    else:
        statuses.append(
            {
                "source": "DailyMed_openFDA_SPL",
                "status": "skipped",
                "path": str(openfda_cache),
                "rows": 0,
                "reason": "disabled by --openfda never",
            }
        )

    context = combine_pk_context(parts)
    manifest = write_pk_context_outputs(
        context=context,
        model_table=model_table,
        out_dir=args.out_dir,
    )
    status_payload = {
        "version": VERSION,
        "description": "Context-preserving PK ingestion for AtlasSPD Phase 1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "model_table": str(args.model_table),
        "sources": statuses,
        "output_manifest": manifest,
        "drugbank_policy": "optional BYOL; no licensed data bundled",
        "label_policy": "spd_exposure_label is not recomputed from external PK",
        "openfda_source_text_review": openfda_review_manifest,
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "pk_source_status.json").write_text(
        json.dumps(status_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    pd.DataFrame(statuses).to_csv(args.out_dir / "pk_source_status.csv", index=False)
    print(json.dumps(status_payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
