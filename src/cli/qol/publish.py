from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence


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

    verify = sub.add_parser(
        "verify",
        help="Verify a generated publication site offline.",
    )
    verify.add_argument("--site-dir", required=True, type=Path)
    verify.add_argument("--report", type=Path)

    edge_bundle = sub.add_parser(
        "edge-bundle",
        help=(
            "Build an optional R2-backed SPA bundle from an existing public site; "
            "does not deploy it."
        ),
    )
    edge_bundle.add_argument("--site-dir", required=True, type=Path)
    edge_bundle.add_argument("--out-dir", type=Path)
    edge_bundle.add_argument("--overwrite", action="store_true")

    args = parser.parse_args(list(argv))
    if args.publish_cmd == "verify":
        forwarded = ["--site-dir", str(args.site_dir)]
        if args.report is not None:
            forwarded.extend(["--report", str(args.report)])
        from analysis.cli import verify_atlas_release

        return int(verify_atlas_release.main(forwarded))

    if args.publish_cmd == "edge-bundle":
        forwarded = ["--site-dir", str(args.site_dir)]
        if args.out_dir is not None:
            forwarded.extend(["--out-dir", str(args.out_dir)])
        if args.overwrite:
            forwarded.append("--overwrite")
        from analysis.reporting import docking_atlas_edge

        return int(docking_atlas_edge.main(forwarded))

    forwarded = [
        str(args.publish_cmd),
        "--manifest",
        str(args.manifest),
    ]
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


__all__ = ["_cmd_publish"]
