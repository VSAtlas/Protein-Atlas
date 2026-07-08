from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.banana import default_banana_root, run_banana_inference


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run BANANA inference on prepared Atlas BANANA inputs.")
    parser.add_argument("--banana-inputs", required=True, type=Path)
    parser.add_argument(
        "--banana-root",
        type=Path,
        default=None,
        help=f"BANANA checkout root (default: {default_banana_root()})",
    )
    parser.add_argument("--python-executable", type=Path, default=None)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--num-workers", type=int, default=4)
    parser.add_argument("--gpu", action="store_true", help="Use BANANA's default cuda:0 path instead of CPU.")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    run_banana_inference(
        args.banana_inputs,
        args.out,
        banana_root=args.banana_root,
        python_executable=args.python_executable,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        no_gpu=not args.gpu,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
