from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.model_run_ledger import write_model_run_record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Write or append a compact Atlas ML model-run ledger row.")
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--ledger", type=Path, default=None)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--task", default=None)
    parser.add_argument("--notes", default="")
    parser.add_argument("--precision-k", type=int, default=20)
    args = parser.parse_args(argv)
    write_model_run_record(
        args.model_dir,
        repo_root=args.repo_root,
        run_id=args.run_id,
        task=args.task,
        notes=args.notes,
        precision_k=args.precision_k,
        ledger_path=args.ledger,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
