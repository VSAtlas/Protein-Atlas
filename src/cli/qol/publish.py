from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Sequence


def _forward_receptor_evidence(args: argparse.Namespace) -> int:
    forwarded = ["--manifest", str(args.manifest), "--out", str(args.out)]
    if args.repo_root is not None:
        forwarded.extend(["--repo-root", str(args.repo_root)])
    for context_key in args.expected_context_key or []:
        forwarded.extend(["--expected-context-key", str(context_key)])
    if args.expected_inventory_sha256 is not None:
        forwarded.extend(
            ["--expected-inventory-sha256", args.expected_inventory_sha256]
        )
    if args.expected_context_count is not None:
        forwarded.extend(
            ["--expected-context-count", str(args.expected_context_count)]
        )
    if args.overwrite:
        forwarded.append("--overwrite")
    from analysis.cli import build_receptor_evidence

    return int(build_receptor_evidence.main(forwarded))


def _forward_verify(args: argparse.Namespace) -> int:
    forwarded = ["--site-dir", str(args.site_dir)]
    if args.report is not None:
        forwarded.extend(["--report", str(args.report)])
    from analysis.cli import verify_atlas_release

    return int(verify_atlas_release.main(forwarded))


def _forward_edge_bundle(args: argparse.Namespace) -> int:
    forwarded = (
        ["--database", str(args.database)]
        if args.database is not None
        else ["--site-dir", str(args.site_dir)]
    )
    optional_paths = (
        ("source-site-dir", args.source_site_dir),
        ("out-dir", args.out_dir),
    )
    for flag, value in optional_paths:
        if value is not None:
            forwarded.extend([f"--{flag}", str(value)])
    if args.overwrite:
        forwarded.append("--overwrite")
    if args.download_base_url:
        forwarded.extend(["--download-base-url", args.download_base_url])
    for name in ("batch_rows", "coarse_shard_rows", "max_pairs"):
        value = getattr(args, name)
        if value is not None:
            forwarded.extend([f"--{name.replace('_', '-')}", str(value)])
    from analysis.reporting import docking_atlas_edge

    return int(docking_atlas_edge.main(forwarded))


def _forward_preflight(args: argparse.Namespace) -> int:
    forwarded = [
        "--edge-dir",
        str(args.edge_dir),
        "--provider",
        str(args.provider),
        "--cloudflare-plan",
        str(args.cloudflare_plan),
    ]
    for name in (
        "site_dir",
        "worker_name",
        "production_bucket",
        "preview_bucket",
        "report",
    ):
        value = getattr(args, name)
        if value is not None:
            forwarded.extend([f"--{name.replace('_', '-')}", str(value)])
    from analysis.reporting import docking_atlas_delivery

    return int(docking_atlas_delivery.main(forwarded))


def _forward_database_command(args: argparse.Namespace) -> int:
    forwarded = [str(args.publish_cmd), "--manifest", str(args.manifest)]
    if args.out_dir is not None:
        forwarded.extend(["--out-dir", str(args.out_dir)])
    if args.publish_cmd == "audit" and args.strict:
        forwarded.append("--strict")
    if args.publish_cmd == "build" and args.overwrite:
        forwarded.append("--overwrite")
    if args.publish_cmd == "build" and args.no_parquet:
        forwarded.append("--no-parquet")
    from analysis.cli import build_atlas_database

    return int(build_atlas_database.main(forwarded))


def _cmd_publish(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas publish",
        description="Audit, build, verify, or project an Atlas publication bundle.",
    )
    sub = parser.add_subparsers(dest="publish_cmd", required=True)

    audit = sub.add_parser(
        "audit",
        help="Audit publication inputs without building the database.",
    )
    audit.add_argument("--manifest", required=True, type=Path)
    audit.add_argument("--out-dir", type=Path)
    audit.add_argument("--strict", action="store_true")

    build = sub.add_parser(
        "build",
        help="Build a publication database from audited Atlas inputs.",
    )
    build.add_argument("--manifest", required=True, type=Path)
    build.add_argument("--out-dir", type=Path)
    build.add_argument("--overwrite", action="store_true")
    build.add_argument("--no-parquet", action="store_true")

    receptor_evidence = sub.add_parser(
        "receptor-evidence",
        help=(
            "Extract exact prepared-receptor chemistry and structured metal "
            "audits without qualifying receptors."
        ),
    )
    receptor_evidence.add_argument("--manifest", required=True, type=Path)
    receptor_evidence.add_argument("--out", required=True, type=Path)
    receptor_evidence.add_argument("--repo-root", type=Path)
    receptor_evidence.add_argument("--expected-context-count", type=int)
    receptor_evidence.add_argument("--expected-context-key", action="append")
    receptor_evidence.add_argument(
        "--expected-inventory-sha256"
    )
    receptor_evidence.add_argument("--overwrite", action="store_true")

    verify = sub.add_parser(
        "verify",
        help="Verify a generated publication site offline.",
    )
    verify.add_argument("--site-dir", required=True, type=Path)
    verify.add_argument("--report", type=Path)

    edge_bundle = sub.add_parser(
        "edge-bundle",
        help=(
            "Build an optional object-backed SPA bundle from a public site or "
            "streaming SQLite snapshot; does not deploy it."
        ),
    )
    edge_source = edge_bundle.add_mutually_exclusive_group(required=True)
    edge_source.add_argument("--site-dir", type=Path)
    edge_source.add_argument("--database", type=Path)
    edge_bundle.add_argument("--source-site-dir", type=Path)
    edge_bundle.add_argument("--out-dir", type=Path)
    edge_bundle.add_argument("--overwrite", action="store_true")
    edge_bundle.add_argument("--download-base-url")
    edge_bundle.add_argument("--batch-rows", type=int)
    edge_bundle.add_argument("--coarse-shard-rows", type=int)
    edge_bundle.add_argument("--max-pairs", type=int)

    preflight = sub.add_parser(
        "preflight",
        help=(
            "Validate an edge delivery bundle and list required resources; "
            "does not upload or deploy."
        ),
    )
    preflight.add_argument("--edge-dir", required=True, type=Path)
    preflight.add_argument("--site-dir", type=Path)
    preflight.add_argument(
        "--provider",
        choices=("generic", "cloudflare-workers-r2"),
        default="generic",
    )
    preflight.add_argument(
        "--cloudflare-plan", choices=("free", "paid"), default="free"
    )
    preflight.add_argument("--worker-name")
    preflight.add_argument("--production-bucket")
    preflight.add_argument("--preview-bucket")
    preflight.add_argument("--report", type=Path)

    args = parser.parse_args(list(argv))
    handler: Callable[[argparse.Namespace], int] | None = {
        "receptor-evidence": _forward_receptor_evidence,
        "verify": _forward_verify,
        "edge-bundle": _forward_edge_bundle,
        "preflight": _forward_preflight,
    }.get(args.publish_cmd)
    return handler(args) if handler is not None else _forward_database_command(args)


__all__ = ["_cmd_publish"]
