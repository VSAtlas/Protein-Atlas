from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.split_locking import DEFAULT_LOCK_SPLITS, lock_run_splits, validate_run_splits


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Lock or validate Atlas ML split manifests for a run.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    lock = sub.add_parser("lock")
    lock.add_argument("--run-id", required=True)
    lock.add_argument("--run-dir", required=True, type=Path)
    lock.add_argument("--out-dir", required=True, type=Path)
    lock.add_argument("--experts", nargs="*", default=None)
    lock.add_argument("--splits", nargs="*", default=DEFAULT_LOCK_SPLITS)
    lock.add_argument("--seed", type=int, default=42)
    lock.add_argument("--validation-fraction", type=float, default=0.15)
    validate = sub.add_parser("validate")
    validate.add_argument("--run-dir", required=True, type=Path)
    validate.add_argument("--split-root", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.cmd == "lock":
        manifest = lock_run_splits(
            run_id=args.run_id,
            run_dir=args.run_dir,
            out_dir=args.out_dir,
            experts=args.experts,
            split_modes=args.splits,
            seed=args.seed,
            validation_fraction=args.validation_fraction,
        )
        print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
        return 0
    result = validate_run_splits(run_dir=args.run_dir, split_root=args.split_root)
    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    return 0 if result.get("status") == "valid" else 1


if __name__ == "__main__":
    raise SystemExit(main())
