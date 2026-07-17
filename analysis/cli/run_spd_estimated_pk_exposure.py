from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.io import load_config
from analysis.ml.spd_estimated_pk_exposure import run_spd_estimated_pk_exposure_model


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Train non-memorizing SPD exposure model: predict AC50 and free Cmax, then combine them."
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--label", default="spd_exposure_label")
    parser.add_argument(
        "--potency-feature-set",
        default="spd_binding_pair_final_full_no_qed",
    )
    parser.add_argument(
        "--pk-feature-set",
        default="ligand_physchem_descriptors_no_qed",
    )
    parser.add_argument("--split", default="drug_holdout")
    parser.add_argument("--model", choices=["ridge", "random_forest", "lightgbm"], default="ridge")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--censored-policy",
        choices=["exclude", "bound"],
        default="exclude",
        help="Exclude right-censored AC50 bounds from potency fitting, or treat bounds as exact for sensitivity only.",
    )
    parser.add_argument("--config", type=Path, default=None)
    args = parser.parse_args(argv)

    config = load_config(args.config)
    seed = int(args.seed if args.seed is not None else config.get("project", {}).get("random_seed", 42))
    manifest = run_spd_estimated_pk_exposure_model(
        args.dataset,
        args.out_dir,
        label_col=args.label,
        potency_feature_set=args.potency_feature_set,
        pk_feature_set=args.pk_feature_set,
        split_mode=args.split,
        model_type=args.model,
        seed=seed,
        censored_policy=args.censored_policy,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
