from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.spd_four_expert_tables import build_spd_four_expert_tables


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build SPD binding/exposure/tissue/mechanism expert tables for Atlas ML."
    )
    parser.add_argument("--spd-table", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--active-um", type=float, default=1.0)
    parser.add_argument("--inactive-um", type=float, default=10.0)
    args = parser.parse_args(argv)
    manifest = build_spd_four_expert_tables(
        args.spd_table,
        args.out_dir,
        active_um=args.active_um,
        inactive_um=args.inactive_um,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
