from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.spd_censor_relabel import relabel_spd_binding_censor_aware


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Recompute SPD binding and combined activity labels with explicit "
            "AC50 censor-relation handling."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--active-um", type=float, default=1.0)
    parser.add_argument("--inactive-um", type=float, default=10.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = relabel_spd_binding_censor_aware(
        args.dataset,
        args.out,
        active_um=args.active_um,
        inactive_um=args.inactive_um,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
