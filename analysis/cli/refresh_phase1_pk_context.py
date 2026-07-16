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
    load_ncats_frdb_pk_context,
    load_reviewed_openfda_pk_context,
    load_spd_pk_context,
    write_pk_context_outputs,
)
from analysis.external.pkdb_api import probe_pkdb_api
from analysis.external.pkdb_recovery import (
    build_canonical_pk_identity_table,
    recover_pkdb_context,
)
from analysis.external.source_tables import download_to_cache
from analysis.external.spl_pk_adjudication import adjudicate_spl_pk_candidates
from analysis.external.spl_pk_candidate_review import review_spl_pk_candidates
from analysis.external.spl_pk_context import load_adjudicated_spl_pk_context


VERSION = "Atlasv0.0.15"
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
DEFAULT_SPD = Path("data/external/spd/sutherland_2023_spd_supplementary_data_1_15.xlsx")


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
        context = load_flat_pk_context(
            path, source_name=source_name, source_version=version
        )
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
        try:
            context = load_ncats_frdb_pk_context(
                candidates[0], source_version=selected_version
            )
        except Exception as exc:
            statuses.append(
                {
                    "source": "NCATS_Inxight_FRDB",
                    "status": "parse_failed",
                    "path": str(candidates[0]),
                    "rows": 0,
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )
        else:
            parts.append(context)
            statuses.append(
                {
                    "source": "NCATS_Inxight_FRDB",
                    "status": "ingested",
                    "path": str(candidates[0]),
                    "rows": int(len(context)),
                    "reason": "source-specific FRDB schema adapter",
                }
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
    model_table: pd.DataFrame,
    source_rights_manifest: Path | None,
    parts: list[pd.DataFrame],
    statuses: list[dict[str, object]],
) -> None:
    source_dir = external_root / "pkdb"
    if source_rights_manifest is not None and not source_rights_manifest.is_file():
        raise FileNotFoundError(
            f"PK-DB source-rights manifest does not exist: {source_rights_manifest}"
        )
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

    extracted = source_dir / "pkdb_all_full"
    archive = source_dir / "pkdb_all_full.zip"
    output_path = extracted / "outputs.csv"
    usable_output = output_path.exists() and output_path.stat().st_size > 3
    if health.get("status") == "available":
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
                        "source": "PK-DB_documented_export",
                        "status": "download_failed",
                        "path": str(archive),
                        "rows": 0,
                        "reason": f"{type(exc).__name__}: {exc}",
                    }
                )
        if archive.exists() and zipfile.is_zipfile(archive) and not usable_output:
            try:
                _extract_archive(archive, extracted)
            except Exception as exc:
                statuses.append(
                    {
                        "source": "PK-DB_documented_export",
                        "status": "archive_failed",
                        "path": str(archive),
                        "rows": 0,
                        "reason": f"{type(exc).__name__}: {exc}",
                    }
                )
        usable_output = output_path.exists() and output_path.stat().st_size > 3
        if usable_output:
            _add_flat_source(
                parts,
                statuses,
                path=output_path,
                source_name="PK-DB",
                version="live_api_export",
            )
            return

    recovery_dir = source_dir / "recovered_phase1"
    recovery_manifest_path = recovery_dir / "pkdb_recovery_manifest.json"
    recovery_context_path = recovery_dir / "pkdb_pk_context.csv"
    recovery_manifest: dict[str, object] = {}
    if download_sources:
        try:
            source_rights = (
                pd.read_csv(source_rights_manifest, dtype="object")
                if source_rights_manifest is not None
                and source_rights_manifest.is_file()
                else None
            )
            recovery_manifest = recover_pkdb_context(
                model_table=model_table,
                out_dir=recovery_dir,
                workers=4,
                source_rights=source_rights,
            )
        except Exception as exc:
            statuses.append(
                {
                    "source": "PK-DB_open_study_recovery",
                    "status": "recovery_failed",
                    "path": str(recovery_dir),
                    "rows": 0,
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )
    elif recovery_manifest_path.exists():
        recovery_manifest = json.loads(
            recovery_manifest_path.read_text(encoding="utf-8")
        )

    recovered = pd.DataFrame()
    if recovery_context_path.exists() and recovery_context_path.stat().st_size > 0:
        try:
            recovered = pd.read_csv(recovery_context_path, low_memory=False)
        except (OSError, pd.errors.ParserError, UnicodeDecodeError) as exc:
            statuses.append(
                {
                    "source": "PK-DB_open_study_recovery",
                    "status": "parse_failed",
                    "path": str(recovery_context_path),
                    "rows": 0,
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )
    if not recovered.empty:
        parts.append(recovered)
        status = (
            "ingested_rights_verified_context"
            if int(str(recovery_manifest.get("model_ready_rows", 0) or 0)) > 0
            else "ingested_context_quarantined"
        )
    else:
        status = "no_model_ready_open_license_context"
    statuses.append(
        {
            "source": "PK-DB_open_study_recovery",
            "status": status,
            "path": str(recovery_context_path),
            "rows": int(len(recovered)),
            "reason": str(
                health.get("reason")
                or "documented output export unavailable; source TSV fallback used"
            ),
            "advertised_outputs": health.get("filter_advertised_outputs", 0),
            "retrieved_outputs": health.get("outputs_endpoint_count", 0),
            "matched_studies": recovery_manifest.get("matched_studies", 0),
            "rights_excluded_studies": recovery_manifest.get(
                "rights_excluded_studies", 0
            ),
            "access_excluded_studies": recovery_manifest.get(
                "public_access_excluded_studies", 0
            ),
            "recovered_context_rows": recovery_manifest.get("recovered_rows", 0),
            "model_ready_rows": recovery_manifest.get("model_ready_rows", 0),
        }
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Refresh contextual PK sources and join them to AtlasSPD Phase 1."
    )
    parser.add_argument("--model-table", type=Path, default=DEFAULT_MODEL_TABLE)
    parser.add_argument(
        "--fda-mapping",
        type=Path,
        default=Path("chemdb/data/fda_mapping_from_pdbqt.csv"),
        help="Validated ligand-base identity mapping used for PK-DB matching.",
    )
    parser.add_argument("--spd-workbook", type=Path, default=DEFAULT_SPD)
    parser.add_argument("--external-root", type=Path, default=Path("data/external"))
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/AtlasSPD_phase1/pk_context_v0_0_15"),
    )
    parser.add_argument(
        "--openfda", choices=("auto", "always", "never"), default="auto"
    )
    parser.add_argument("--openfda-max-drugs", type=int, default=0)
    parser.add_argument("--openfda-sleep-sec", type=float, default=0.25)
    parser.add_argument("--skip-openfda-source-review", action="store_true")
    parser.add_argument(
        "--spl-review-decisions",
        type=Path,
        default=Path(
            "data/external/dailymed_spl/phase1_openfda/spl_pk_review_decisions.csv"
        ),
        help=(
            "Explicit reviewed decisions for SPL clearance, absolute "
            "bioavailability, and maximum-dose candidates. Missing files are "
            "treated as no approvals, not as an error."
        ),
    )
    parser.add_argument(
        "--reviewed-openfda-context",
        type=Path,
        default=None,
        help=(
            "Previously source-text-adjudicated DailyMed context CSV. "
            "These rows can be reused with --openfda never."
        ),
    )
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument(
        "--pkdb-source-rights-manifest",
        type=Path,
        default=None,
        help=(
            "Audited PK-DB study-rights CSV. Without an affirmative per-study "
            "training_allowed decision and rights reference, recovered rows "
            "remain contextual and are excluded from training."
        ),
    )
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
    if not args.fda_mapping.is_file():
        raise FileNotFoundError(
            f"canonical FDA mapping does not exist: {args.fda_mapping}"
        )
    pk_identity_table = build_canonical_pk_identity_table(
        model_table,
        pd.read_csv(args.fda_mapping, low_memory=False),
    )
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
        model_table=pk_identity_table,
        source_rights_manifest=args.pkdb_source_rights_manifest,
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
    spl_candidate_manifest: dict[str, object] = {}
    spl_adjudication_manifest: dict[str, object] = {}
    if args.reviewed_openfda_context is not None:
        if not args.reviewed_openfda_context.is_file():
            raise FileNotFoundError(
                "reviewed openFDA context does not exist: "
                f"{args.reviewed_openfda_context}"
            )
        reviewed_context = load_reviewed_openfda_pk_context(
            args.reviewed_openfda_context
        )
        parts.append(reviewed_context)
        statuses.append(
            {
                "source": "DailyMed_openFDA_SPL_reviewed",
                "status": "ingested",
                "path": str(args.reviewed_openfda_context),
                "rows": int(len(reviewed_context)),
                "reason": "source-text-adjudicated contextual PK only",
            }
        )
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
                "status": "ingested"
                if not openfda_context.empty
                else "no_numeric_rows",
                "path": str(openfda_cache),
                "rows": int(len(openfda_context)),
                "reason": ""
                if not openfda_context.empty
                else "no label context rows matched",
                **openfda_manifest,
            }
        )
        if not args.skip_openfda_source_review:
            cmax_review_dir = args.out_dir / "openfda_source_text_review"
            openfda_review_manifest = audit_openfda_pk_cache(
                cache_dir=openfda_cache / "records",
                model_table=model_table,
                out_dir=cmax_review_dir,
            )
            cmax_adjudicated = cmax_review_dir / "openfda_spl_cmax_adjudicated.csv"
            reviewed_cmax = load_reviewed_openfda_pk_context(cmax_adjudicated)
            if not reviewed_cmax.empty:
                parts.append(reviewed_cmax)
            statuses.append(
                {
                    "source": "DailyMed_openFDA_SPL_Cmax_reviewed",
                    "status": "ingested"
                    if not reviewed_cmax.empty
                    else "no_model_ready_rows",
                    "path": str(cmax_adjudicated),
                    "rows": int(len(reviewed_cmax)),
                    "reason": "dose bound to the adjudicated Cmax scenario",
                }
            )

            spl_review_dir = args.out_dir / "spl_pk_source_review"
            spl_candidate_manifest = review_spl_pk_candidates(
                cache_dir=openfda_cache / "records",
                out_dir=spl_review_dir / "candidates",
            )
            spl_adjudication_manifest = adjudicate_spl_pk_candidates(
                candidates=spl_review_dir
                / "candidates"
                / "spl_pk_candidate_review.csv",
                out_dir=spl_review_dir / "adjudicated",
                review_decisions=(
                    args.spl_review_decisions
                    if args.spl_review_decisions.is_file()
                    else None
                ),
            )
            adjudicated_context = load_adjudicated_spl_pk_context(
                clearance_path=(
                    spl_review_dir / "adjudicated" / "accepted_clearance_context.csv"
                ),
                bioavailability_path=(
                    spl_review_dir
                    / "adjudicated"
                    / "accepted_absolute_bioavailability_context.csv"
                ),
            )
            if not adjudicated_context.empty:
                parts.append(adjudicated_context)
            statuses.append(
                {
                    "source": "DailyMed_openFDA_SPL_semantic_review",
                    "status": (
                        "ingested"
                        if not adjudicated_context.empty
                        else "no_model_ready_rows"
                    ),
                    "path": str(spl_review_dir / "adjudicated"),
                    "rows": int(len(adjudicated_context)),
                    "reason": (
                        "machine same-clause gate or structured human-reviewed absolute "
                        "bioavailability context; maximum recommended dose remains sensitivity-only"
                    ),
                }
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
        "pkdb_identity_mapping": {
            "path": str(args.fda_mapping),
            "source_rows": int(len(model_table)),
            "canonical_rows": int(len(pk_identity_table)),
            "canonical_identity_rows": int(
                pk_identity_table.get(
                    "pk_identity_source",
                    pd.Series("", index=pk_identity_table.index),
                )
                .fillna("")
                .eq("canonical_fda_mapping_validated_parent")
                .sum()
            ),
        },
        "sources": statuses,
        "output_manifest": manifest,
        "drugbank_policy": "optional BYOL; no licensed data bundled",
        "label_policy": "spd_exposure_label is not recomputed from external PK",
        "openfda_source_text_review": openfda_review_manifest,
        "spl_pk_candidate_review": spl_candidate_manifest,
        "spl_pk_adjudication": spl_adjudication_manifest,
        "spl_review_decisions": (
            str(args.spl_review_decisions)
            if args.spl_review_decisions.is_file()
            else "not_provided"
        ),
        "pkdb_source_rights_manifest": (
            str(args.pkdb_source_rights_manifest)
            if args.pkdb_source_rights_manifest is not None
            and args.pkdb_source_rights_manifest.is_file()
            else "not_provided"
        ),
        "dose_policy": {
            "primary": "dose matched to the adjudicated Cmax study",
            "sensitivity": "maximum labeled recommended adult dose in a separate artifact",
            "excluded_from_primary": "highest studied dose and maximum tolerated dose",
        },
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
