from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from analysis.ml.pk_mechanistic_features import (
    SPD_FREE_CMAX_TARGET,
    TARGET_ENDPOINTS,
    materialize_pk_mechanistic_features,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Copy recognized contextual PK variables into explicit pk_mech_* "
            "columns without conversion or imputation."
        )
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help=(
            "Coverage manifest path. Defaults to "
            "<out stem>.pk_mechanistic_manifest.json."
        ),
    )
    parser.add_argument(
        "--target-endpoint",
        choices=TARGET_ENDPOINTS,
        default=SPD_FREE_CMAX_TARGET,
        help=(
            "Exposure target used for endpoint-alignment gates. The default "
            "keeps external pk_context administration audit-only for the SPD "
            "free_cmax_um target."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.input.is_file():
        raise FileNotFoundError(f"Input CSV not found: {args.input}")
    manifest_path = args.manifest or args.out.with_suffix(
        ".pk_mechanistic_manifest.json"
    )
    if manifest_path.resolve() == args.out.resolve():
        raise ValueError("Manifest path must differ from output CSV path")

    source = pd.read_csv(args.input, low_memory=False)
    augmented, manifest = materialize_pk_mechanistic_features(
        source,
        target_endpoint=args.target_endpoint,
    )
    manifest["input"] = str(args.input)
    manifest["output"] = str(args.out)
    manifest["manifest"] = str(manifest_path)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    augmented.to_csv(args.out, index=False)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
