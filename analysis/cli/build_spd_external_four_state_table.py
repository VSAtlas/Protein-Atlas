from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.spd_external_four_state_merge import build_spd_external_four_state_merged_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Merge an SPD model-ready table with external four-state evidence while preserving separate and combined labels."
    )
    parser.add_argument("--spd-table", required=True, type=Path)
    parser.add_argument("--four-state-table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--run-dir", type=Path, default=None)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--no-tissue-labels", action="store_true")
    parser.add_argument("--no-refresh-metadata", action="store_true")
    args = parser.parse_args(argv)
    summary = build_spd_external_four_state_merged_table(
        args.spd_table,
        args.four_state_table,
        args.out,
        run_dir=args.run_dir,
        repo_root=args.repo_root,
        add_tissue_labels=not args.no_tissue_labels,
        refresh_metadata=not args.no_refresh_metadata,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
