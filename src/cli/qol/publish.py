from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence


def _cmd_publish(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas publish",
        description="Audit or build a publication database from an Atlas manifest.",
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

    args = parser.parse_args(list(argv))
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

    from analysis.cli import build_atlas_database

    return int(build_atlas_database.main(forwarded))


__all__ = ["_cmd_publish"]
