from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from analysis.ml.target_representation import write_target_pocket_features


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build one deterministic physicochemical/geometry descriptor row "
            "per Atlas pocket PDB."
        )
    )
    parser.add_argument("--pocket-map", type=Path, required=True, nargs="+")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--model-table", type=Path, default=None)
    parser.add_argument("--joined-out", type=Path, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    missing = [path for path in args.pocket_map if not path.exists()]
    if missing:
        raise FileNotFoundError(f"pocket maps not found: {missing}")
    pocket_map = pd.concat(
        [pd.read_csv(path, low_memory=False) for path in args.pocket_map],
        ignore_index=True,
    )
    manifest = write_target_pocket_features(
        pocket_map,
        out_dir=args.out_dir,
        cache_root=args.cache_dir,
    )
    if args.model_table is not None:
        if not args.model_table.exists():
            raise FileNotFoundError(f"model table not found: {args.model_table}")
        model_table = pd.read_csv(args.model_table, low_memory=False)
        target_table = pd.read_csv(manifest["table"], low_memory=False)
        if "pdb_id" not in model_table:
            raise ValueError("model table must include pdb_id")
        model_table["pdb_id"] = model_table["pdb_id"].astype(str).str.upper()
        target_table["pdb_id"] = target_table["pdb_id"].astype(str).str.upper()
        replace = [
            column
            for column in target_table.columns
            if column != "pdb_id" and column in model_table.columns
        ]
        joined = model_table.drop(columns=replace).merge(
            target_table,
            on="pdb_id",
            how="left",
            validate="many_to_one",
        )
        joined_out = args.joined_out or args.out_dir / "target_pocket_model_ready.csv"
        joined_out.parent.mkdir(parents=True, exist_ok=True)
        joined.to_csv(joined_out, index=False)
        manifest["joined_model_table"] = str(joined_out)
        manifest["joined_rows"] = int(len(joined))
        (args.out_dir / "target_pocket_features_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
