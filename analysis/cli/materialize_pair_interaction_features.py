from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from analysis.ml.pair_interaction_features import (
    DEFAULT_PAIR_KEY_COLUMNS,
    GEOMETRY_PROVIDER_NAME,
    PAIR_FEATURE_PROVIDER_REGISTRY,
    POSE_POLICY_EXACT,
    POSE_POLICY_STAGE1_MODEL1,
    materialize_pair_interaction_features,
)


def _json_object(value: str) -> dict[str, Any]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise argparse.ArgumentTypeError(f"invalid JSON: {exc.msg}") from exc
    if not isinstance(parsed, dict):
        raise argparse.ArgumentTypeError("provider config must be a JSON object")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Append content-bound, PDBQT-defensible pair-interaction features. "
            "Exact pose provenance is required unless the explicit Stage-1/model-1 "
            "sensitivity policy is selected."
        )
    )
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--manifest-out",
        type=Path,
        help="Default: <out>.manifest.json",
    )
    parser.add_argument(
        "--pose-policy",
        choices=(POSE_POLICY_EXACT, POSE_POLICY_STAGE1_MODEL1),
        default=POSE_POLICY_EXACT,
        help=(
            "exact-supplied validates pose_path/pose_sha256/pose_model_index and "
            "receptor_path/receptor_sha256. stage1-model-1 is sensitivity-only."
        ),
    )
    parser.add_argument(
        "--frozen-run-id",
        help="Required only for --pose-policy stage1-model-1.",
    )
    parser.add_argument(
        "--frozen-run-root",
        type=Path,
        help=(
            "Frozen Atlas repository or relocated bundle root used by the canonical "
            "run output router; required for stage1-model-1."
        ),
    )
    parser.add_argument(
        "--pair-key-column",
        action="append",
        dest="pair_key_columns",
        help=(
            "Repeat to override the default exact key: "
            + ", ".join(DEFAULT_PAIR_KEY_COLUMNS)
        ),
    )
    parser.add_argument(
        "--provider",
        action="append",
        choices=PAIR_FEATURE_PROVIDER_REGISTRY.names(),
        dest="providers",
        help=f"Repeat to append provider blocks; default: {GEOMETRY_PROVIDER_NAME}.",
    )
    parser.add_argument(
        "--provider-config-json",
        type=_json_object,
        default={},
        help=(
            "JSON object keyed by provider name. Configuration is versioned and "
            "content-hashed into every output row and the manifest."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Explicitly replace existing output and manifest files.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest_out = args.manifest_out or Path(f"{args.out}.manifest.json")
    provider_configs = args.provider_config_json
    selected_providers = args.providers or (GEOMETRY_PROVIDER_NAME,)
    unknown_configs = sorted(set(provider_configs) - set(selected_providers))
    if unknown_configs:
        raise ValueError(
            "configuration supplied for providers that were not selected: "
            f"{unknown_configs}"
        )
    manifest = materialize_pair_interaction_features(
        args.dataset,
        args.out,
        manifest_out,
        pose_policy=args.pose_policy,
        pair_key_columns=args.pair_key_columns or DEFAULT_PAIR_KEY_COLUMNS,
        provider_names=selected_providers,
        provider_configs=provider_configs,
        frozen_run_id=args.frozen_run_id,
        frozen_run_root=args.frozen_run_root,
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=True, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
