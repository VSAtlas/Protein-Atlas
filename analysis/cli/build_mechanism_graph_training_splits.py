from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.mechanism_graph_splits import DEFAULT_SPLIT_MODES, build_mechanism_graph_training_splits


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build KGE/hypergraph-ready mechanism graph edge manifests and negative-edge policy files."
    )
    parser.add_argument("--mechanism-pu-table", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--panel", default="cardiac_qt")
    parser.add_argument("--label-col", default="mechanism_pu_label")
    parser.add_argument("--split-modes", nargs="*", default=list(DEFAULT_SPLIT_MODES))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--validation-fraction", type=float, default=0.10)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    args = parser.parse_args(argv)
    manifest = build_mechanism_graph_training_splits(
        args.mechanism_pu_table,
        args.out_dir,
        panel=args.panel,
        split_modes=args.split_modes,
        label_col=args.label_col,
        seed=args.seed,
        validation_fraction=args.validation_fraction,
        test_fraction=args.test_fraction,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
