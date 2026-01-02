#!/usr/bin/env python
import argparse
import json
import sys
from pathlib import Path

import utils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert SMILES to a DeepCoy-compatible valid JSON file.")
    parser.add_argument("--in-smi", required=True, help="Input SMILES file.")
    parser.add_argument("--out-json", required=True, help="Output JSON file.")
    parser.add_argument("--dataset", default="zinc", help="Dataset name (default: zinc).")
    parser.add_argument("--max-n", type=int, default=None, help="Maximum number of valid rows to emit.")
    return parser.parse_args()


def main():
    args = parse_args()

    in_path = Path(args.in_smi)
    out_path = Path(args.out_json)

    if not in_path.is_file():
        print(f"ERROR: Input SMILES file not found: {in_path}", file=sys.stderr)
        return 2

    try:
        dataset_meta = utils.dataset_info(args.dataset)
    except Exception:
        print(f"ERROR: Unknown dataset: {args.dataset}", file=sys.stderr)
        return 2

    max_bucket = int(max(dataset_meta["bucket_sizes"])) - 1

    rows = []
    total = 0
    skipped = 0

    with in_path.open("r") as f:
        for line in f:
            total += 1
            line = line.strip()
            if not line:
                continue

            smiles = line.split()[0]
            nodes, edges = utils.to_graph(smiles, args.dataset)
            if not nodes or not edges:
                skipped += 1
                continue

            if len(nodes) > max_bucket:
                skipped += 1
                continue

            rows.append({
                "graph_in": edges,
                "graph_out": edges,
                "node_features_in": nodes,
                "node_features_out": nodes,
                "smiles_in": smiles,
                "smiles_out": smiles,
            })

            if args.max_n is not None and len(rows) >= args.max_n:
                break

    print(f"Processed {total} lines. Valid: {len(rows)}. Skipped: {skipped}.")

    if not rows:
        print("ERROR: No valid SMILES entries found.", file=sys.stderr)
        return 4

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(rows, f, indent=2)

    return 0


if __name__ == "__main__":
    sys.exit(main())
