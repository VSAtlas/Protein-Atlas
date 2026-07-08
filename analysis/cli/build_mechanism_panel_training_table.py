from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.mechanism_panel_training import build_mechanism_panel_training_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build an ADR-panel-specific mechanism training table while preserving quarantined audit columns.",
    )
    parser.add_argument("--source-table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--panel", default="cardiac_qt")
    parser.add_argument("--source-label-col", default="mechanism_pu_label")
    parser.add_argument("--label-col", default=None)
    parser.add_argument("--topk-out", default=None, type=Path)
    parser.add_argument("--feature-source-table", default=None, type=Path)
    args = parser.parse_args(argv)
    manifest = build_mechanism_panel_training_table(
        args.source_table,
        args.out,
        panel=args.panel,
        source_label_col=args.source_label_col,
        label_col=args.label_col,
        topk_out=args.topk_out,
        feature_source_table=args.feature_source_table,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
