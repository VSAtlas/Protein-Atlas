from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.descriptor_constellation_audit import (
    run_descriptor_constellation_audit,
)
from analysis.ml.feature_sets import get_feature_set


def _resolution(value: str) -> tuple[str, float]:
    feature, separator, raw = value.partition("=")
    if not separator or not feature.strip():
        raise argparse.ArgumentTypeError("resolutions must use FEATURE=VALUE")
    try:
        resolution = float(raw)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("resolution must be numeric") from exc
    if resolution <= 0:
        raise argparse.ArgumentTypeError("resolution must be positive")
    return feature.strip(), resolution


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Measure whether RDKit descriptor combinations can act as drug-identity "
            "proxies and inspect nearest chemical neighbors."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--label", default="spd_binding_label")
    parser.add_argument("--feature-set")
    parser.add_argument("--feature", action="append", default=None)
    parser.add_argument(
        "--resolution",
        action="append",
        type=_resolution,
        default=[],
        metavar="FEATURE=VALUE",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.feature_set and args.feature:
        raise SystemExit("use either --feature-set or --feature, not both")
    descriptors = get_feature_set(args.feature_set) if args.feature_set else args.feature
    resolutions = dict(args.resolution)
    if len(resolutions) != len(args.resolution):
        raise SystemExit("duplicate --resolution feature")
    result = run_descriptor_constellation_audit(
        args.dataset,
        args.out_dir,
        label_col=args.label,
        descriptors=descriptors,
        resolutions=resolutions,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
