from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.mechanism_panel_audit import write_panel_readiness_audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit ADR/site panel readiness for mechanism PU labels.")
    parser.add_argument("--mechanism-pu-table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--label-col", default="mechanism_pu_label")
    args = parser.parse_args(argv)
    write_panel_readiness_audit(
        args.mechanism_pu_table,
        args.out,
        label_col=args.label_col,
        source_balance_path=args.out.with_name(args.out.stem + "_source_balance.csv"),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
