from __future__ import annotations

import argparse
from pathlib import Path

from analysis.ml.graph_model_backends import write_graph_model_backend_audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit optional graph/KGE ML backend availability.")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    write_graph_model_backend_audit(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
