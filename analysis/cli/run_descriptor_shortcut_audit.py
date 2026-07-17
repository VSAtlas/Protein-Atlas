from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.ml.descriptor_shortcut_audit import run_descriptor_shortcut_audit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run grouped leave-one-RDKit-descriptor-out models and drug-level "
            "label-association audits."
        )
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--label", default="spd_binding_label")
    parser.add_argument(
        "--feature-set",
        default="spd_binding_nonleaky_consensus_z",
    )
    parser.add_argument("--model", default="lightgbm")
    parser.add_argument(
        "--group-cols",
        nargs="+",
        default=["drug_id", "chemical_cluster", "target_id", "target_family"],
    )
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--bootstraps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--model-n-jobs", type=int, default=8)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = run_descriptor_shortcut_audit(
        args.dataset,
        args.out_dir,
        label_col=args.label,
        feature_set=args.feature_set,
        model_type=args.model,
        group_cols=args.group_cols,
        n_splits=args.n_splits,
        n_bootstraps=args.bootstraps,
        seed=args.seed,
        model_n_jobs=args.model_n_jobs,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
