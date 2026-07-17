from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Sequence


from config.output_paths import output_root, run_output_dir
from cli.qol.pipeline import _cmd_run
from cli.qol.utils import _display_path
from cli.qol import _bindings


def _cmd_ligands(argv: Sequence[str]) -> int:
    tokens = list(argv)
    if tokens[:2] == ["audit", "stereo"]:
        from analysis.cli.audit_ligand_stereo import main as audit_stereo

        return audit_stereo(tokens[2:])
    parser = _build_ligands_parser()
    args, pipeline_args = parser.parse_known_args(tokens)
    if args.ligand_cmd != "run" and pipeline_args:
        parser.error("unrecognized arguments: " + " ".join(pipeline_args))
    special_result = _dispatch_special_ligand_command(args)
    if special_result is not None:
        return special_result
    return _run_ligand_action(args, pipeline_args)


def _dispatch_special_ligand_command(args: argparse.Namespace) -> int | None:
    if args.ligand_cmd == "sources":
        return _print_ligand_sources()
    if args.ligand_cmd == "audit" and args.audit_source == "fda":
        return _cmd_ligands_audit_fda(args)
    if args.ligand_cmd == "repair":
        return _dispatch_fda_repair_command(args)
    if args.ligand_cmd == "promote" and args.promotion_source == "fda-v3":
        return _cmd_ligands_promote_fda_v3(args)
    return None


def _dispatch_fda_repair_command(args: argparse.Namespace) -> int | None:
    if args.repair_source == "fda":
        return _cmd_ligands_repair_fda(args)
    if args.repair_source == "fda-delta":
        return _cmd_ligands_repair_fda_delta(args)
    if args.repair_source == "fda-legacy":
        return _cmd_ligands_repair_fda_legacy(args)
    if args.repair_source == "fda-legacy-bind-run":
        return _cmd_ligands_repair_fda_legacy_bind_run(args)
    if args.repair_source == "fda-terminal-v3":
        return _cmd_ligands_repair_fda_terminal_v3(args)
    return None


def _build_ligands_parser() -> argparse.ArgumentParser:
    from prep_ligands.ligand_library_manager import known_sources

    source_choices = sorted(
        {
            token
            for source in known_sources()
            for token in (source.name, *source.aliases)
        }
    )
    parser = argparse.ArgumentParser(
        prog="atlas ligands",
        description="Download, prepare, and run Atlas ligand libraries.",
    )
    sub = parser.add_subparsers(dest="ligand_cmd", required=True)
    sub.add_parser("sources", help="List built-in ligand library sources.")
    _add_fda_audit_parser(sub)
    _add_fda_repair_parser(sub)
    _add_fda_promotion_parser(sub)
    _add_ligand_action_parsers(sub, source_choices)
    return parser


def _add_fda_audit_parser(sub: argparse._SubParsersAction) -> None:
    audit = sub.add_parser(
        "audit",
        help="Audit a local ligand library without installing or downloading data.",
    )
    audit_sub = audit.add_subparsers(dest="audit_source", required=True)
    audit_sub.add_parser(
        "stereo",
        add_help=False,
        help="Report source-defined and unresolved ligand stereochemistry.",
    )
    fda_audit = audit_sub.add_parser(
        "fda",
        help="Cross-check FDA ligand names, structures, salts, and eligibility.",
        description=(
            "Audit the local FDA mapping against available source and prepared-library "
            "evidence. This command is offline and never installs a library."
        ),
    )
    fda_audit.add_argument(
        "--mapping-csv",
        help=(
            "FDA mapping CSV (default: <repo>/chemdb/data/fda_mapping_from_pdbqt.csv)."
        ),
    )
    fda_audit.add_argument(
        "--source-sdf",
        help=(
            "FDA source SDF used for stable-ID and structure verification. FDA "
            "approval additionally requires a validated filter manifest."
        ),
    )
    fda_audit.add_argument(
        "--approval-manifest",
        help=(
            "DrugCentral exact-filter manifest (default: fda_filter_manifest.json "
            "beside --source-sdf)."
        ),
    )
    fda_audit.add_argument(
        "--library-dir",
        help="Optional prepared FDA PDBQT library directory.",
    )
    fda_audit.add_argument(
        "--out-dir",
        help="Audit output directory (default: <repo>/outputs/data/fda_identity_audit).",
    )
    fda_audit.add_argument(
        "--classified-non-medication-csv",
        help=(
            "Curated non-medication classifications (default: "
            "<repo>/docs/fda_non_medication_substances.csv)."
        ),
    )
    fda_audit.add_argument(
        "--score-csv",
        action="append",
        default=[],
        help=(
            "Explicit historical score CSV to reconcile (repeatable). No score "
            "files are discovered automatically."
        ),
    )
    fda_audit.add_argument(
        "--score-run-id",
        action="append",
        default=[],
        help=(
            "Explicit run ID whose fda_master_rows.csv (or master_rows.csv) is "
            "reconciled; repeatable."
        ),
    )
    fda_audit.add_argument(
        "--fail-on-blockers",
        action="store_true",
        help="Return exit code 1 when the audit finds blocking identity mismatches.",
    )
    fda_audit.add_argument(
        "--fail-on-blocking",
        action="store_true",
        dest="fail_on_blockers",
        help=argparse.SUPPRESS,
    )


def _add_fda_repair_parser(sub: argparse._SubParsersAction) -> None:
    repair = sub.add_parser(
        "repair",
        help="Build a separate, provenance-backed repaired ligand library.",
    )
    repair_sub = repair.add_subparsers(dest="repair_source", required=True)
    fda_repair = repair_sub.add_parser(
        "fda",
        help="Materialize the manifest-approved FDA canonical-parent library.",
        description=(
            "Copy compatible legacy FDA parent PDBQTs byte-for-byte and prepare "
            "only the approved parent delta. The original library is unchanged "
            "and this command never launches docking."
        ),
    )
    fda_repair.add_argument("--repaired-mapping-csv", required=True)
    fda_repair.add_argument("--redock-delta-csv", required=True)
    fda_repair.add_argument("--source-sdf", required=True)
    fda_repair.add_argument("--approval-manifest", required=True)
    fda_repair.add_argument("--output-dir", required=True)
    fda_repair.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Parallel Meeko preparation workers, capped at 32 (default: 8).",
    )
    fda_repair.add_argument(
        "--plan-only",
        action="store_true",
        help="Validate and report the complete cohort without writing any files.",
    )
    fda_delta = repair_sub.add_parser(
        "fda-delta",
        help="Extract the ready FDA parent delta without rerunning preparation.",
    )
    fda_delta.add_argument("--named-library-dir", required=True)
    fda_delta.add_argument("--output-dir", required=True)
    fda_delta.add_argument(
        "--copy",
        action="store_true",
        help="Copy bytes instead of preferring same-filesystem hardlinks.",
    )
    legacy = repair_sub.add_parser(
        "fda-legacy",
        help="Write the versioned, non-destructive legacy FDA reconciliation bundle.",
        description=(
            "Preserve every usable legacy PDBQT, apply manifest-backed identity "
            "corrections only in a new mapping, attach corrected identities and "
            "checksums to explicit score tables, and retain unresolved rows for "
            "a second pass. No input or score table is overwritten."
        ),
    )
    legacy.add_argument("--original-mapping-csv", required=True)
    legacy.add_argument("--repaired-mapping-csv", required=True)
    legacy.add_argument("--redock-delta-csv", required=True)
    legacy.add_argument("--named-library-manifest-csv")
    legacy.add_argument(
        "--score-csv",
        action="append",
        default=[],
        help=(
            "Explicit score CSV (repeatable). When omitted, only the three "
            "hash-pinned SPD FDA master tables are used."
        ),
    )
    legacy.add_argument("--output-dir", required=True)
    bind_run = repair_sub.add_parser(
        "fda-legacy-bind-run",
        help="Hash-bind a legacy reconciliation bundle to an active FDA delta run.",
        description=(
            "Write one post-launch binding JSON beside the legacy reconciliation "
            "bundle. The command verifies the ready-parent hashes, frozen targets, "
            "delta summary, and live run manifest without rewriting the pre-launch "
            "redock linkage or any run output."
        ),
    )
    bind_run.add_argument("--reconciliation-dir", required=True)
    bind_run.add_argument("--run-id", required=True)
    bind_run.add_argument("--target-dir", required=True)
    bind_run.add_argument("--target-manifest-csv", required=True)
    bind_run.add_argument("--delta-summary-json", required=True)
    bind_run.add_argument("--run-manifest-yaml", required=True)
    bind_run.add_argument("--cpu-count", type=int, required=True)
    bind_run.add_argument("--observed-worker-concurrency", type=int)
    bind_run.add_argument("--expected-target-count", type=int, default=93)
    bind_run.add_argument("--expected-ready-count", type=int, default=586)
    bind_run.add_argument("--expected-apo-holo-mode", default="holo")
    bind_run.add_argument("--expected-run-status", default="running")
    terminal_v3 = repair_sub.add_parser(
        "fda-terminal-v3",
        help="Resolve every legacy FDA-map row to a terminal structure-first disposition.",
        description=(
            "Cross-check the v2 reconciliation against full DrugCentral, a "
            "validated PubChem cache, Drugs@FDA, and the canonical named library. "
            "Outputs are versioned; no input or PDBQT is modified."
        ),
    )
    terminal_v3.add_argument("--mapping-v2-csv", required=True)
    terminal_v3.add_argument("--unresolved-v2-csv", required=True)
    terminal_v3.add_argument("--repaired-mapping-csv", required=True)
    terminal_v3.add_argument("--full-drugcentral-sdf", required=True)
    terminal_v3.add_argument("--approved-source-sdf", required=True)
    terminal_v3.add_argument("--fda-approved-csv", required=True)
    terminal_v3.add_argument("--pubchem-records-jsonl", required=True)
    terminal_v3.add_argument("--pubchem-manifest-json", required=True)
    terminal_v3.add_argument("--pubchem-inchikey-records-jsonl", required=True)
    terminal_v3.add_argument("--pubchem-inchikey-manifest-json", required=True)
    terminal_v3.add_argument("--drugsfda-products-txt", required=True)
    terminal_v3.add_argument("--named-library-manifest-csv", required=True)
    terminal_v3.add_argument("--repair-quarantine-csv", required=True)
    terminal_v3.add_argument("--output-dir", required=True)


def _add_fda_promotion_parser(sub: argparse._SubParsersAction) -> None:
    promote = sub.add_parser(
        "promote",
        help="Promote a verified versioned ligand mapping without changing defaults.",
    )
    promote_sub = promote.add_subparsers(dest="promotion_source", required=True)
    fda_v3 = promote_sub.add_parser(
        "fda-v3",
        help="Write a hash-bound promotion manifest for a verified FDA v3 mapping.",
        description=(
            "Require a passing independent verifier report whose mapping SHA-256 "
            "matches the supplied v3 CSV. The mapping is not copied and no config "
            "is changed unless --update-config names an existing file."
        ),
    )
    fda_v3.add_argument(
        "--mapping-v3-csv",
        required=True,
        help="Explicit fda_legacy_mapping_v3.csv path.",
    )
    fda_v3.add_argument(
        "--verifier-json",
        required=True,
        help="Explicit passing independent verifier JSON path.",
    )
    fda_v3.add_argument("--output-dir", required=True)
    fda_v3.add_argument(
        "--score-csv",
        action="append",
        default=[],
        help=(
            "Explicit score CSV to annotate (repeatable). Joins require both an "
            "exact RDK ID and matching PDBQT SHA-256."
        ),
    )
    fda_v3.add_argument(
        "--update-config",
        metavar="CONFIG_PATH",
        help=(
            "Explicitly update FDA_MAPPING_CSV in this existing config file. "
            "No config file is changed when omitted."
        ),
    )


def _add_ligand_action_parsers(
    sub: argparse._SubParsersAction, source_choices: Sequence[str]
) -> None:
    for name in ("fetch", "prep", "install", "run"):
        cmd = sub.add_parser(name, help=f"{name} a ligand library source.")
        cmd.add_argument("source", choices=source_choices)
        cmd.add_argument(
            "--source-sdf", help="Use a local SDF/SDF.GZ/ZIP instead of downloading."
        )
        cmd.add_argument("--source-url", help="Override the built-in source URL.")
        cmd.add_argument("--force-fetch", action="store_true")
        cmd.add_argument("--force-prep", action="store_true")
        cmd.add_argument(
            "--limit", type=int, default=0, help="Prepare only the first N SDF records."
        )
        cmd.add_argument(
            "--dry-run",
            action="store_true",
            help="Show planned paths and commands without fetching, prepping, or running.",
        )
        if name == "run":
            cmd.add_argument("--no-install", action="store_true")


def _run_ligand_action(args: argparse.Namespace, pipeline_args: Sequence[str]) -> int:
    from prep_ligands.ligand_library_manager import resolve_source

    cfg = _bindings.load_effective_config()
    source = resolve_source(args.source)
    limit = int(args.limit or 0) or None
    source_sdf = Path(args.source_sdf).expanduser() if args.source_sdf else None
    if args.dry_run:
        return _dry_run_ligand_action(
            args, pipeline_args, cfg, source, source_sdf, limit
        )
    try:
        return _execute_ligand_action(
            args, pipeline_args, cfg, source, source_sdf, limit
        )
    except RuntimeError as exc:
        print(f"atlas ligands: {exc}", file=sys.stderr)
        return 2


def _dry_run_ligand_action(
    args: argparse.Namespace,
    pipeline_args: Sequence[str],
    cfg: dict[str, Any],
    source: Any,
    source_sdf: Path | None,
    limit: int | None,
) -> int:
    from prep_ligands.ligand_library_manager import paths_for_source

    paths = paths_for_source(cfg, source)
    _print_ligand_paths("planned", source.name, paths)
    print(f"action: {args.ligand_cmd}")
    if source_sdf:
        print(f"source_sdf: {_display_path(source_sdf)}")
    if args.source_url:
        print(f"source_url: {args.source_url}")
    if limit is not None:
        print(f"limit: {limit}")
    if args.ligand_cmd == "run":
        forwarded = _forwarded_pipeline_args(pipeline_args)
        print(
            "planned_pipeline_args: " + (" ".join(forwarded) if forwarded else "<none>")
        )
    return 0


def _execute_ligand_action(
    args: argparse.Namespace,
    pipeline_args: Sequence[str],
    cfg: dict[str, Any],
    source: Any,
    source_sdf: Path | None,
    limit: int | None,
) -> int:
    from prep_ligands.ligand_library_manager import (
        fetch_source,
        install_source,
        paths_for_source,
        prepare_source,
    )

    if args.ligand_cmd == "fetch":
        paths = fetch_source(
            source.name,
            cfg,
            source_sdf=source_sdf,
            source_url=args.source_url,
            force=args.force_fetch,
            limit=limit,
        )
        _print_ligand_paths("fetched", source.name, paths)
        return 0
    if args.ligand_cmd == "prep":
        if source_sdf or args.source_url:
            fetch_source(
                source.name,
                cfg,
                source_sdf=source_sdf,
                source_url=args.source_url,
                force=args.force_fetch,
                limit=limit,
            )
        paths = prepare_source(source.name, cfg, force=args.force_prep, limit=limit)
        _print_ligand_paths("prepared", source.name, paths)
        return 0
    if args.ligand_cmd == "install":
        paths = install_source(
            source.name,
            cfg,
            source_sdf=source_sdf,
            source_url=args.source_url,
            force_fetch=args.force_fetch,
            force_prep=args.force_prep,
            limit=limit,
        )
        _print_ligand_paths("installed", source.name, paths)
        return 0
    if args.ligand_cmd == "run":
        paths = paths_for_source(cfg, source)
        if not args.no_install:
            paths = install_source(
                source.name,
                cfg,
                source_sdf=source_sdf,
                source_url=args.source_url,
                force_fetch=args.force_fetch,
                force_prep=args.force_prep,
                limit=limit,
            )
            _print_ligand_paths("installed", source.name, paths)
        return _run_with_library(source.name, list(pipeline_args or []))
    return 1


def _cmd_ligands_audit_fda(args: argparse.Namespace) -> int:
    from prep_ligands.fda_identity_audit import run_fda_identity_audit

    repo_root = _bindings.repo_root()

    try:
        config = _fda_audit_config(args, repo_root)
        result = run_fda_identity_audit(config)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"atlas ligands audit fda: {exc}", file=sys.stderr)
        return 2

    _print_fda_audit_result(result, config.output_dir, repo_root)
    return 1 if args.fail_on_blockers and result.has_blockers else 0


def _cmd_ligands_repair_fda(args: argparse.Namespace) -> int:
    from prep_ligands.fda_named_library import build_named_fda_library

    try:
        summary = build_named_fda_library(
            repaired_mapping_csv=Path(args.repaired_mapping_csv),
            redock_delta_csv=Path(args.redock_delta_csv),
            source_sdf=Path(args.source_sdf),
            approval_manifest=Path(args.approval_manifest),
            output_dir=Path(args.output_dir),
            workers=int(args.workers),
            plan_only=bool(args.plan_only),
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"atlas ligands repair fda: {exc}", file=sys.stderr)
        return 2
    print("ligand_repair: fda")
    print(f"status: {'planned' if args.plan_only else 'complete'}")
    print(f"canonical_parent_rows: {summary['canonical_parent_rows']}")
    print(f"action_counts: {json.dumps(summary['action_counts'], sort_keys=True)}")
    print(f"status_counts: {json.dumps(summary['status_counts'], sort_keys=True)}")
    print(f"ready_pdbqt_count: {summary['ready_pdbqt_count']}")
    print(
        "output_dir: "
        + _display_path(Path(args.output_dir).expanduser(), root=_bindings.repo_root())
    )
    return 0


def _cmd_ligands_repair_fda_delta(args: argparse.Namespace) -> int:
    from prep_ligands.fda_delta_library import materialize_fda_delta_library

    try:
        summary = materialize_fda_delta_library(
            named_library_dir=Path(args.named_library_dir),
            output_dir=Path(args.output_dir),
            prefer_hardlinks=not bool(args.copy),
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"atlas ligands repair fda-delta: {exc}", file=sys.stderr)
        return 2
    print("ligand_repair: fda-delta")
    print("status: complete")
    print(f"delta_manifest_rows: {summary['delta_manifest_rows']}")
    print(f"ready_pdbqt_count: {summary['ready_pdbqt_count']}")
    print(f"quarantine_count: {summary['quarantine_count']}")
    print(
        "materialization_methods: "
        + json.dumps(summary["pdbqt_materialization_methods"], sort_keys=True)
    )
    print(
        "output_dir: "
        + _display_path(Path(args.output_dir).expanduser(), root=_bindings.repo_root())
    )
    return 0


def _cmd_ligands_repair_fda_legacy(args: argparse.Namespace) -> int:
    from prep_ligands.fda_legacy_reconciliation import reconcile_legacy_fda

    repo_root = _bindings.repo_root()
    score_csvs = (
        [Path(value).expanduser() for value in args.score_csv]
        if args.score_csv
        else _default_fda_legacy_score_csvs(repo_root)
    )
    try:
        outputs, summary = reconcile_legacy_fda(
            original_mapping_csv=Path(args.original_mapping_csv),
            repaired_mapping_csv=Path(args.repaired_mapping_csv),
            redock_delta_csv=Path(args.redock_delta_csv),
            score_csvs=score_csvs,
            output_dir=Path(args.output_dir),
            named_library_manifest_csv=Path(args.named_library_manifest_csv)
            if args.named_library_manifest_csv
            else None,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"atlas ligands repair fda-legacy: {exc}", file=sys.stderr)
        return 2
    print("ligand_repair: fda-legacy")
    print("status: complete")
    print(f"counts: {json.dumps(summary['counts'], sort_keys=True)}")
    for label, path in (
        ("versioned_mapping_csv", outputs.versioned_mapping_csv),
        ("identity_changes_csv", outputs.identity_changes_csv),
        ("file_inventory_csv", outputs.file_inventory_csv),
        ("score_identity_joins_csv", outputs.score_identity_joins_csv),
        ("unresolved_second_pass_csv", outputs.unresolved_second_pass_csv),
        ("redock_linkage_csv", outputs.redock_linkage_csv),
        ("summary_json", outputs.summary_json),
    ):
        print(f"{label}: {_display_path(path, root=repo_root)}")
    return 0


def _cmd_ligands_repair_fda_legacy_bind_run(args: argparse.Namespace) -> int:
    from prep_ligands.fda_run_binding import bind_legacy_reconciliation_run

    repo_root = _bindings.repo_root()
    try:
        output, payload = bind_legacy_reconciliation_run(
            reconciliation_dir=Path(args.reconciliation_dir),
            run_id=str(args.run_id),
            target_dir=Path(args.target_dir),
            target_manifest_csv=Path(args.target_manifest_csv),
            delta_summary_json=Path(args.delta_summary_json),
            run_manifest_yaml=Path(args.run_manifest_yaml),
            cpu_count=int(args.cpu_count),
            observed_worker_concurrency=args.observed_worker_concurrency,
            expected_target_count=int(args.expected_target_count),
            expected_ready_count=int(args.expected_ready_count),
            expected_apo_holo_mode=str(args.expected_apo_holo_mode),
            expected_run_status=str(args.expected_run_status),
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"atlas ligands repair fda-legacy-bind-run: {exc}", file=sys.stderr)
        return 2
    binding_sha256 = hashlib.sha256(output.binding_json.read_bytes()).hexdigest()
    print("ligand_repair: fda-legacy-bind-run")
    print("status: complete")
    print(f"run_id: {payload['run']['run_id']}")
    print(f"run_status: {payload['run']['observed_status']}")
    print(
        "frozen_target_count: "
        f"{payload['frozen_target_manifest']['target_count']}"
    )
    print(
        "ready_ligand_count: "
        f"{payload['ready_ligand_binding']['ready_row_count']}"
    )
    print(f"binding_json_sha256: {binding_sha256}")
    print(f"binding_json: {_display_path(output.binding_json, root=repo_root)}")
    return 0


def _cmd_ligands_repair_fda_terminal_v3(args: argparse.Namespace) -> int:
    from prep_ligands.fda_terminal_resolver import resolve_fda_terminal_v3

    repo_root = _bindings.repo_root()
    try:
        outputs, summary = resolve_fda_terminal_v3(
            mapping_v2_csv=Path(args.mapping_v2_csv),
            unresolved_v2_csv=Path(args.unresolved_v2_csv),
            repaired_mapping_csv=Path(args.repaired_mapping_csv),
            full_drugcentral_sdf=Path(args.full_drugcentral_sdf),
            approved_source_sdf=Path(args.approved_source_sdf),
            fda_approved_csv=Path(args.fda_approved_csv),
            pubchem_records_jsonl=Path(args.pubchem_records_jsonl),
            pubchem_manifest_json=Path(args.pubchem_manifest_json),
            pubchem_inchikey_records_jsonl=Path(
                args.pubchem_inchikey_records_jsonl
            ),
            pubchem_inchikey_manifest_json=Path(
                args.pubchem_inchikey_manifest_json
            ),
            drugsfda_products_txt=Path(args.drugsfda_products_txt),
            named_library_manifest_csv=Path(args.named_library_manifest_csv),
            repair_quarantine_csv=Path(args.repair_quarantine_csv),
            output_dir=Path(args.output_dir),
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"atlas ligands repair fda-terminal-v3: {exc}", file=sys.stderr)
        return 2
    print("ligand_repair: fda-terminal-v3")
    print(f"status: {summary['completion_status']}")
    print(f"counts: {json.dumps(summary['counts'], sort_keys=True)}")
    for label, path in (
        ("mapping_v3_csv", outputs.mapping_v3_csv),
        ("terminal_dispositions_csv", outputs.terminal_dispositions_csv),
        ("fda_subset_csv", outputs.fda_subset_csv),
        ("excluded_subset_csv", outputs.excluded_subset_csv),
        ("summary_json", outputs.summary_json),
    ):
        print(f"{label}: {_display_path(path, root=repo_root)}")
    return 0 if summary["complete"] else 1


def _cmd_ligands_promote_fda_v3(args: argparse.Namespace) -> int:
    from prep_ligands.fda_mapping_promotion import (
        promote_verified_fda_mapping_v3,
    )

    repo_root = _bindings.repo_root()
    try:
        outputs, manifest = promote_verified_fda_mapping_v3(
            mapping_v3_csv=Path(args.mapping_v3_csv),
            verifier_json=Path(args.verifier_json),
            output_dir=Path(args.output_dir),
            score_csvs=tuple(Path(value) for value in args.score_csv),
            update_config=Path(args.update_config) if args.update_config else None,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"atlas ligands promote fda-v3: {exc}", file=sys.stderr)
        return 2
    manifest_sha256 = hashlib.sha256(
        outputs.promotion_manifest_json.read_bytes()
    ).hexdigest()
    print("ligand_promotion: fda-v3")
    print("status: complete")
    print(f"mapping_sha256: {manifest['mapping']['sha256']}")
    print(f"mapping_rows: {manifest['mapping']['rows']}")
    print(f"verifier_sha256: {manifest['independent_verifier']['sha256']}")
    print(f"config_changed: {str(manifest['config_update']['changed']).lower()}")
    print(f"promotion_manifest_sha256: {manifest_sha256}")
    if outputs.score_identity_joins_csv is not None:
        print(
            "score_identity_joins_csv: "
            + _display_path(outputs.score_identity_joins_csv, root=repo_root)
        )
    print(
        "promotion_manifest_json: "
        + _display_path(outputs.promotion_manifest_json, root=repo_root)
    )
    return 0


def _default_fda_legacy_score_csvs(repo_root: Path) -> list[Path]:
    return [
        output_root(repo_root, "data") / "SPD_addons" / "fda_master_rows.csv",
        output_root(repo_root, "data") / "SPDaddon2" / "fda_master_rows.csv",
        output_root(repo_root, "data") / "SPDsparsefloor" / "fda_master_rows.csv",
    ]


def _expanded_cli_path(value: str | None) -> Path | None:
    return Path(value).expanduser() if value else None


def _fda_audit_config(args: argparse.Namespace, repo_root: Path) -> Any:
    from prep_ligands.fda_identity_audit import FDAIdentityAuditConfig

    return FDAIdentityAuditConfig(
        mapping_csv=_expanded_cli_path(args.mapping_csv)
        or repo_root / "chemdb" / "data" / "fda_mapping_from_pdbqt.csv",
        source_sdf=_expanded_cli_path(args.source_sdf),
        approval_manifest=_expanded_cli_path(args.approval_manifest),
        library_dir=_expanded_cli_path(args.library_dir),
        output_dir=_expanded_cli_path(args.out_dir)
        or output_root(repo_root, "data") / "fda_identity_audit",
        classified_nonmed_csv=_expanded_cli_path(args.classified_non_medication_csv)
        or repo_root / "docs" / "fda_non_medication_substances.csv",
        score_csvs=_fda_score_csvs(args, repo_root),
    )


_RUN_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


def _fda_score_csvs(args: argparse.Namespace, repo_root: Path) -> tuple[Path, ...]:
    paths = [Path(value).expanduser() for value in args.score_csv]
    for raw_run_id in args.score_run_id:
        run_id = str(raw_run_id).strip()
        if not _RUN_ID_RE.fullmatch(run_id):
            raise ValueError(f"invalid FDA score run ID: {raw_run_id!r}")
        run_dir = run_output_dir(repo_root, "data", run_id)
        candidates = (run_dir / "fda_master_rows.csv", run_dir / "master_rows.csv")
        score_csv = next((path for path in candidates if path.is_file()), None)
        if score_csv is None:
            raise ValueError(
                f"no fda_master_rows.csv or master_rows.csv for score run {run_id!r} "
                f"under {run_dir}"
            )
        paths.append(score_csv)
    return tuple(dict.fromkeys(path.resolve() for path in paths))


def _print_fda_audit_result(result: Any, output_dir: Path, repo_root: Path) -> None:
    print("ligand_audit: fda")
    print("status: complete")
    print(f"total_rows: {result.total_rows}")
    print(f"blocker_count: {result.blocker_count}")
    print(f"review_count: {result.review_count}")
    print(f"output_dir: {_display_path(output_dir, root=repo_root)}")
    for label, path in (
        ("row_audit_csv", result.row_audit_csv),
        ("evidence_jsonl", result.evidence_jsonl),
        ("summary_json", result.summary_json),
        ("review_csv", result.review_csv),
        ("quarantine_manifest_csv", result.quarantine_manifest_csv),
        ("eligible_manifest_csv", result.eligible_manifest_csv),
        ("source_catalog_csv", result.source_catalog_csv),
        (
            "approval_cross_reference_csv",
            getattr(result, "approval_cross_reference_csv", None),
        ),
        ("repaired_mapping_csv", getattr(result, "repaired_mapping_csv", None)),
        (
            "score_reuse_manifest_csv",
            getattr(result, "score_reuse_manifest_csv", None),
        ),
        ("redock_delta_csv", getattr(result, "redock_delta_csv", None)),
        (
            "repair_quarantine_csv",
            getattr(result, "repair_quarantine_csv", None),
        ),
    ):
        if path is None:
            continue
        print(f"{label}: {_display_path(path, root=repo_root)}")


def _print_ligand_sources() -> int:
    from prep_ligands.ligand_library_manager import known_sources

    for source in known_sources():
        aliases = f" aliases={','.join(source.aliases)}" if source.aliases else ""
        print(f"{source.name}: {source.description}{aliases}")
        print(f"  library_subdir: {source.library_subdir}")
        print(f"  url: {source.url}")
        print(f"  license: {source.license_note}")
        print(f"  citation: {source.citation}")
    print(
        "Use `atlas ligands install chembl`, `chebi`, `coconut`, `fda`, or `hmdb` to download and prepare."
    )
    return 0


def _print_ligand_paths(action: str, source_name: str, paths: Any) -> None:
    root = _bindings.repo_root()
    print(f"ligand_source: {source_name}")
    print(f"status: {action}")
    print(f"raw_sdf: {_display_path(paths.raw_sdf, root=root)}")
    print(f"library_dir: {_display_path(paths.library_dir, root=root)}")
    print(f"source_manifest: {_display_path(paths.source_manifest, root=root)}")


def _cmd_library_alias(source_name: str, argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog=f"atlas {source_name}",
        description=(
            f"Install the {source_name.upper()} ligand library if needed, then run Atlas "
            f"with TEST_MODE_ENABLE={source_name}."
        ),
    )
    parser.add_argument("--no-install", action="store_true")
    parser.add_argument("--source-sdf")
    parser.add_argument("--source-url")
    parser.add_argument("--force-fetch", action="store_true")
    parser.add_argument("--force-prep", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    args, pipeline_args = parser.parse_known_args(list(argv))

    from prep_ligands.ligand_library_manager import (
        install_source,
        paths_for_source,
        resolve_source,
    )

    cfg = _bindings.load_effective_config()
    source = resolve_source(source_name)
    if args.dry_run:
        _print_ligand_paths("planned", source.name, paths_for_source(cfg, source))
        forwarded = _forwarded_pipeline_args(pipeline_args)
        print(
            "planned_pipeline_args: " + (" ".join(forwarded) if forwarded else "<none>")
        )
        return 0
    if not args.no_install:
        try:
            paths = install_source(
                source.name,
                cfg,
                source_sdf=Path(args.source_sdf).expanduser()
                if args.source_sdf
                else None,
                source_url=args.source_url,
                force_fetch=args.force_fetch,
                force_prep=args.force_prep,
                limit=int(args.limit or 0) or None,
            )
        except RuntimeError as exc:
            print(f"atlas {source.name}: {exc}", file=sys.stderr)
            return 2
        _print_ligand_paths("installed", source.name, paths)
    else:
        _print_ligand_paths(
            "using_existing", source.name, paths_for_source(cfg, source)
        )
    return _run_with_library(source.name, list(pipeline_args or []))


def _forwarded_pipeline_args(pipeline_args: Sequence[str] | None) -> list[str]:
    """Drop an optional argparse-style ``--`` separator between known and forwarded args."""
    forwarded = list(pipeline_args or [])
    if forwarded and forwarded[0] == "--":
        forwarded.pop(0)
    return forwarded


def _run_with_library(source_name: str, pipeline_args: Sequence[str]) -> int:
    env_updates = {"TEST_MODE_ENABLE": source_name}
    old_env = {key: os.environ.get(key) for key in env_updates}
    try:
        for key, value in env_updates.items():
            os.environ[key] = value
        forwarded = _forwarded_pipeline_args(pipeline_args)
        print(f"running_with_library: {source_name}")
        print("pipeline_args: " + (" ".join(forwarded) if forwarded else "<none>"))
        return _cmd_run(forwarded)
    finally:
        for key, old_value in old_env.items():
            if old_value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = old_value
