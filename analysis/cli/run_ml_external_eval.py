from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.external_eval import run_external_model_eval


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate a trained Atlas ML model on an external labeled table.")
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", default=None)
    parser.add_argument("--name", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    manifest = run_external_model_eval(model_dir=args.model_dir, dataset=args.dataset, label_col=args.label, name=args.name, out_dir=args.out_dir)
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
