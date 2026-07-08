from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.banana import build_banana_input_table, run_banana_inference
from analysis.external.banana_pockets import build_banana_pocket_map
from analysis.io import load_config
from analysis.cli.run_mechanism_pu_evaluation import run_leave_one_source_out
from analysis.ml.audit_suite import run_ml_audit_suite
from analysis.ml.dataset_sampling import materialize_sampled_dataset, materialize_split_locked_sampled_dataset
from analysis.ml.feature_metadata import enrich_ml_feature_metadata, refresh_tables
from analysis.ml.feature_sets import get_feature_set
from analysis.ml.mechanism_pu_tables import build_mechanism_pu_table
from analysis.ml.mechanism_recovery import write_mechanism_topk_recovery
from analysis.ml.spd_estimated_pk_exposure import run_spd_estimated_pk_exposure_model
from analysis.ml.split_manifest import dataframe_content_hash
from analysis.ml.spd_four_expert_tables import build_spd_four_expert_tables
from analysis.ml.train_classifier import train_ml_model
from analysis.ml.tissue_pu_recovery import run_tissue_pu_recovery_report


EXPERT_DEFAULTS: dict[str, dict[str, Any]] = {
    "binding": {
        "label": "spd_binding_label",
        "feature_set": "spd_binding_nonleaky",
        "split": "drug_holdout",
        "candidates": [
            "spd_four_experts/model_ready/ml_spd_binding_model_ready.csv",
        ],
        "audit_candidates": [
            "source_benchmark_comparison_v2/tables/model_ready/ml_toxcast_hts_model_ready.csv",
            "source_benchmark_comparison_v2/tables/model_ready/ml_spd_exposure_model_ready.csv",
            "source_benchmark_comparison_v2/tables/model_ready/ml_matched_source_transfer_model_ready.csv",
        ],
    },
    "exposure": {
        "label": "spd_exposure_relevant",
        "feature_set": "spd_exposure_nonleaky",
        "split": "drug_holdout",
        "candidates": [
            "spd_four_experts/model_ready/ml_spd_exposure_model_ready.csv",
            "source_benchmark_comparison_v2/tables/ml_spd_exposure_table.csv",
            "spd_exposure_nonleaky_model_ready.csv",
            "source_benchmark_comparison_v2/tables/model_ready/ml_spd_exposure_model_ready.csv",
        ],
        "audit_candidates": [
            "source_benchmark_comparison_v2/tables/model_ready/ml_curated_bioactivity_model_ready.csv",
            "source_benchmark_comparison_v2/tables/model_ready/ml_toxcast_hts_model_ready.csv",
        ],
    },
    "tissue": {
        "label": "tissue_site_label",
        "feature_set": "spd_tissue_site_nonleaky",
        "split": "drug_holdout",
        "candidates": [
            "spd_four_experts/model_ready/ml_spd_tissue_model_ready.csv",
            "tissue_site/model_ready/ml_tissue_site_model_ready.csv",
        ],
        "audit_candidates": [
            "source_benchmark_comparison_v2/tables/model_ready/ml_spd_exposure_model_ready.csv",
        ],
    },
    "mechanism": {
        "label": "mechanism_ml_label",
        "feature_set": "spd_mechanism_nonleaky",
        "split": "drug_holdout",
        "candidates": [
            "spd_four_experts/model_ready/ml_spd_mechanism_model_ready.csv",
            "negative_evidence/mechanism_four_state_labels_deep_sources_v7_adrecs_positive_source_balanced.csv",
            "negative_evidence/mechanism_four_state_labels_deep_sources_source_balanced.csv",
            "negative_evidence/mechanism_four_state_labels.csv",
        ],
        "audit_candidates": [
            "adr_mechanism_panel/adr_mechanism_panel_table.csv",
            "negative_evidence/mechanism_four_state_labels.csv",
        ],
    },
}



BANANA_REQUEST_FEATURES = {
    "banana_score",
    "banana_binding_probability",
    "banana_score_normalized",
    "banana_atlas_blend_score",
    "binding_expert_score",
}


def _feature_set_requests_banana(feature_set: str) -> bool:
    try:
        return bool(BANANA_REQUEST_FEATURES.intersection(get_feature_set(feature_set)))
    except Exception:
        return False


def _maybe_populate_banana_scores(
    dataset: Path,
    *,
    expert: str,
    label_col: str,
    feature_set: str,
    run_dir: Path,
    out_dir: Path,
    repo_root: Path,
    mode: str,
    scope: str,
    batch_size: int,
    num_workers: int,
) -> Path:
    if mode == "never" or not _feature_set_requests_banana(feature_set):
        return dataset
    df = pd.read_csv(dataset, low_memory=False)
    labelable = pd.to_numeric(df.get(label_col, pd.Series(pd.NA, index=df.index)), errors="coerce").notna()
    score = pd.to_numeric(df.get("banana_score", pd.Series(pd.NA, index=df.index)), errors="coerce")
    target_mask = labelable if scope == "labelable" else pd.Series(True, index=df.index)
    target_rows = int(target_mask.sum())
    missing_target_scores = int((target_mask & score.isna()).sum())
    requested_banana_features = [
        feature
        for feature in get_feature_set(feature_set)
        if feature in BANANA_REQUEST_FEATURES and feature in df.columns
    ]
    cached_banana_nonmissing = {
        feature: int(pd.to_numeric(df.loc[target_mask, feature], errors="coerce").notna().sum())
        for feature in requested_banana_features
        if feature not in {"banana_score_feature_source", "banana_binding_probability_feature_source"}
    }
    banana_dir = out_dir / "_banana_features" / expert
    banana_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = banana_dir / "banana_feature_population_manifest.json"
    manifest: dict[str, object] = {
        "dataset": str(dataset),
        "expert": expert,
        "label_col": label_col,
        "feature_set": feature_set,
        "mode": mode,
        "scope": scope,
        "target_rows": target_rows,
        "missing_target_scores": missing_target_scores,
        "requested_banana_features": requested_banana_features,
        "cached_banana_nonmissing": cached_banana_nonmissing,
        "status": "pending",
    }
    if target_rows == 0:
        manifest["status"] = "skipped_no_target_rows"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
        return dataset
    if mode == "auto" and any(count > 0 for count in cached_banana_nonmissing.values()):
        manifest["status"] = "skipped_cached_banana_features_present"
        manifest["policy"] = "auto mode treats existing BANANA-derived feature columns as cache; use --banana-scoring always to fill missing raw BANANA rows"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
        return dataset
    if mode == "auto" and missing_target_scores == 0:
        manifest["status"] = "skipped_scores_already_present"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
        return dataset

    pocket_map = banana_dir / "pocket_map.csv"
    if not pocket_map.exists():
        build_banana_pocket_map(
            dataset,
            pocket_map,
            raw_pdb_dir=banana_dir / "raw_pdbs",
            pocket_dir=banana_dir / "pockets",
            receptor_search_dirs=[repo_root / "outputs/processed_pdbs", repo_root / "input_pdbs"],
            allow_download=True,
            require_complete=False,
        )
    banana_inputs = banana_dir / "banana_inputs.csv"
    if not banana_inputs.exists():
        build_banana_input_table(
            dataset,
            banana_inputs,
            fda_mapping_path=repo_root / "chemdb/data/fda_mapping_from_pdbqt.csv",
            pocket_map_path=pocket_map,
        )
    input_df = pd.read_csv(banana_inputs, low_memory=False)
    if scope == "labelable":
        infer_inputs = banana_dir / "banana_inputs.labelable.csv"
        input_df.loc[labelable].to_csv(infer_inputs, index=False)
    else:
        infer_inputs = banana_inputs
    infer_df = pd.read_csv(infer_inputs, low_memory=False)
    ready_rows = int(infer_df.get("banana_input_status", pd.Series("", index=infer_df.index)).astype(str).eq("ready").sum())
    manifest["ready_rows"] = ready_rows
    if ready_rows == 0:
        manifest["status"] = "skipped_no_ready_banana_inputs"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
        return dataset
    scores = banana_dir / ("banana_scores.labelable.csv" if scope == "labelable" else "banana_scores.csv")
    if not scores.exists():
        run_banana_inference(
            infer_inputs,
            scores,
            batch_size=batch_size,
            num_workers=num_workers,
            no_gpu=True,
        )
    enriched, summary = enrich_ml_feature_metadata(
        df,
        run_dir=run_dir,
        repo_root=repo_root,
        extra_feature_tables=[scores],
    )
    out_path = banana_dir / f"{dataset.stem}.with_banana.csv"
    enriched.to_csv(out_path, index=False)
    manifest.update(
        {
            "status": "written",
            "pocket_map": str(pocket_map),
            "banana_inputs": str(banana_inputs),
            "inference_inputs": str(infer_inputs),
            "banana_scores": str(scores),
            "output_dataset": str(out_path),
            "feature_metadata_summary": summary,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return out_path

def _first_existing(run_dir: Path, rel_paths: list[str]) -> Path | None:
    for rel in rel_paths:
        path = run_dir / rel
        if path.exists():
            return path
    return None


def _configure_mlflow_defaults(*, out_dir: Path, run_id: str, disabled: bool) -> None:
    if disabled:
        os.environ["ATLAS_DISABLE_MLFLOW"] = "1"
        return
    os.environ.setdefault("ATLAS_ENABLE_MLFLOW", "1")
    os.environ.setdefault("ATLAS_MLFLOW_TRACKING_URI", f"file:{out_dir.parent / 'mlruns'}")
    os.environ.setdefault("ATLAS_MLFLOW_EXPERIMENT", f"atlas/{run_id}/ml")


def _hash_policy(payload: dict[str, object]) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode("utf-8")).hexdigest()[:16]


def _split_manifest_overrides(values: list[str] | None) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in values or []:
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key and value:
            out[key] = Path(value)
    return out


def _resolve_run_dir(repo_root: Path, run_id: str, run_dir: Path | None) -> Path:
    if run_dir is not None:
        return run_dir.resolve()
    candidates = [
        repo_root / "data" / run_id,
        repo_root / "outputs" / "data" / run_id,
        repo_root / "outputs" / run_id / "outputs" / "data" / run_id,
        repo_root / "outputs" / run_id,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _dataset_label_counts(dataset: Path, label_col: str) -> dict[str, Any]:
    df = pd.read_csv(dataset, usecols=lambda col: col == label_col, low_memory=False)
    labels = pd.to_numeric(df[label_col], errors="coerce")
    return {
        "rows": int(len(df)),
        "labeled_rows": int(labels.notna().sum()),
        "label_counts": {
            str(key): int(value)
            for key, value in labels.value_counts(dropna=False).items()
        },
    }


def _existing_audit_candidates(run_dir: Path, rel_paths: list[str]) -> list[str]:
    out: list[str] = []
    for rel in rel_paths:
        path = run_dir / rel
        if path.exists():
            out.append(f"{Path(rel).stem}={path}")
    return out


SPD_SOURCE_CANDIDATES = [
    "source_benchmark_comparison_v2/tables/ml_spd_exposure_table.csv",
    "spd_atlas_benchmark.csv",
    "pair_feature_table.csv",
    "master_rows.csv",
]


def _prepare_spd_expert_tables(
    *,
    run_dir: Path,
    out_dir: Path,
    mode: str,
    source_table: Path | None,
    spd_panel_path: Path | None = None,
    target_map_path: Path | None = None,
    ligand_map_path: Path | None = None,
    target_metadata_path: Path | None = None,
    mechanism_label_path: Path | None = None,
    target_expression_path: Path | None = None,
) -> dict[str, Path]:
    if mode == "never":
        return {}
    source = source_table or _first_existing(run_dir, SPD_SOURCE_CANDIDATES)
    if source is None or not source.exists():
        if mode == "always":
            raise FileNotFoundError(
                "no SPD source table found; pass --spd-expert-source-table or build source_benchmark_comparison_v2 first"
            )
        return {}
    out = out_dir / "_spd_expert_tables"
    manifest = build_spd_four_expert_tables(
        source,
        out,
        spd_panel_path=spd_panel_path,
        target_map_path=target_map_path,
        ligand_map_path=ligand_map_path,
        target_metadata_path=target_metadata_path,
        mechanism_label_path=mechanism_label_path,
        target_expression_path=target_expression_path,
        run_dir=run_dir,
    )
    prepared: dict[str, Path] = {}
    for expert, info in (manifest.get("experts") or {}).items():
        model_ready = info.get("model_ready_path")
        if model_ready and Path(model_ready).exists():
            prepared[str(expert)] = Path(model_ready)
    return prepared


def _resolve_group_reweight_cols(dataset: Path, requested: list[str] | None) -> list[str] | None:
    if not requested:
        return None
    if requested == ["none"]:
        return None
    if "auto" not in requested:
        return requested
    header = pd.read_csv(dataset, nrows=5, low_memory=False)
    candidates = [
        "external_source_family",
        "source_family",
        "label_source",
        "external_upstream_source",
        "upstream_source",
        "target_family",
        "protein_class",
    ]
    resolved = [col for col in candidates if col in header.columns]
    return resolved or None


def _metadata_refreshed_dataset(
    dataset: Path,
    *,
    run_dir: Path,
    out_dir: Path,
    mode: str,
    chemical_cluster: str,
    target_family: str,
    source_lineage: str,
    drop_columns: list[str] | None,
    use_cache: bool,
    refresh_cache: bool,
) -> Path:
    if mode == "never":
        return dataset
    policy = {
        "chemical_cluster": chemical_cluster,
        "target_family": target_family,
        "source_lineage": source_lineage,
        "drop_columns": drop_columns or [],
    }
    if use_cache:
        frame = pd.read_csv(dataset, low_memory=False)
        source_hash = dataframe_content_hash(frame)
        cache_dir = out_dir / "feature_cache" / source_hash / _hash_policy(policy)
        manifest_path = cache_dir / "feature_cache_manifest.json"
        if manifest_path.exists() and not refresh_cache:
            try:
                cached = json.loads(manifest_path.read_text(encoding="utf-8"))
                cached_output = Path(str(cached.get("output_dataset", "")))
                if cached_output.exists() and cached.get("source_dataset_hash") == source_hash:
                    return cached_output
            except Exception:
                pass
        refresh_dir = cache_dir
    else:
        source_hash = None
        refresh_dir = out_dir / "_feature_metadata"
    manifest = refresh_tables(
        [dataset],
        out_dir=refresh_dir,
        in_place=False,
        run_dir=run_dir,
        chemical_cluster=chemical_cluster,
        target_family=target_family,
        source_lineage=source_lineage,
        drop_columns=drop_columns,
    )
    outputs = manifest.get("outputs") or []
    output_dataset = Path(outputs[0]) if outputs else dataset
    if use_cache:
        cache_manifest = {
            "source_dataset": str(dataset),
            "source_dataset_hash": source_hash,
            "output_dataset": str(output_dataset),
            "policy": policy,
            "refresh_manifest": manifest,
        }
        refresh_dir.mkdir(parents=True, exist_ok=True)
        (refresh_dir / "feature_cache_manifest.json").write_text(json.dumps(cache_manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return output_dataset



def _run_estimated_pk_exposure(
    *,
    dataset: Path,
    label_col: str,
    potency_feature_set: str,
    pk_feature_set: str,
    split_mode: str,
    out_dir: Path,
    seed: int,
    model_type: str,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "expert": "exposure_estimated_pk",
        "dataset": str(dataset),
        "label": label_col,
        "potency_feature_set": potency_feature_set,
        "pk_feature_set": pk_feature_set,
        "split": split_mode,
        "status": "pending",
    }
    try:
        manifest = run_spd_estimated_pk_exposure_model(
            dataset,
            out_dir / "exposure_estimated_pk",
            label_col=label_col,
            potency_feature_set=potency_feature_set,
            pk_feature_set=pk_feature_set,
            split_mode=split_mode,
            model_type=model_type,
            seed=seed,
        )
        result.update(
            {
                "status": "trained",
                "model_dir": str(out_dir / "exposure_estimated_pk"),
                "metrics": manifest.get("metrics", {}),
                "potency_features": manifest.get("potency_features", []),
                "pk_features": manifest.get("pk_features", []),
                "policy": manifest.get("policy", []),
                "warning": manifest.get("warning"),
            }
        )
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = str(exc)
    return result

def _run_expert(
    *,
    name: str,
    dataset: Path,
    label_col: str,
    feature_set: str,
    split_mode: str,
    out_dir: Path,
    seed: int,
    model_type: str,
    bootstraps: int,
    run_audit: bool,
    class_weight: str | None,
    allow_partial_rescoring_features: bool,
    allow_label_definition_features: bool,
    nested_model_selection: bool,
    pu_mode: str,
    pu_bags: int,
    pu_unlabeled_ratio: float,
    candidates: list[str],
    source_transfer_train: list[str] | None,
    source_transfer_test: list[str] | None,
    claim_mode: str,
    repo_root: Path,
    split_manifest: Path | None = None,
    dataset_provenance: dict[str, object] | None = None,
    calibration_method: str = "none",
    strict_feature_set: bool = False,
    model_params: dict[str, object] | None = None,
    group_reweight_cols: list[str] | None = None,
    group_reweight_include_label: bool = True,
    group_reweight_max_factor: float = 5.0,
) -> dict[str, Any]:
    expert_dir = out_dir / name
    model_dir = expert_dir / "model"
    result: dict[str, Any] = {
        "expert": name,
        "dataset": str(dataset),
        "label": label_col,
        "feature_set": feature_set,
        "split": split_mode,
        "status": "pending",
    }
    try:
        result["label_balance"] = _dataset_label_counts(dataset, label_col)
        trained = train_ml_model(
            dataset,
            label_col,
            feature_set,
            model_type,
            split_mode,
            model_dir,
            seed,
            n_bootstraps=bootstraps,
            class_weight=class_weight,
            allow_partial_rescoring_features=allow_partial_rescoring_features,
            allow_label_definition_features=allow_label_definition_features,
            nested_model_selection=nested_model_selection,
            pu_mode=pu_mode,
            pu_bags=pu_bags,
            pu_unlabeled_ratio=pu_unlabeled_ratio,
            claim_mode=claim_mode,
            repo_root=repo_root,
            split_manifest=split_manifest,
            dataset_provenance=dataset_provenance,
            calibration_method=calibration_method,
            strict_feature_set=strict_feature_set,
            model_params=model_params,
            group_reweight_cols=group_reweight_cols,
            group_reweight_include_label=group_reweight_include_label,
            group_reweight_max_factor=group_reweight_max_factor,
        )
        result.update(
            {
                "status": "trained",
                "model_dir": str(model_dir),
                "metrics": trained.get("metrics", {}),
                "n_train": trained.get("n_train"),
                "n_test": trained.get("n_test"),
                "features": trained.get("features", []),
                "claim_readiness": trained.get("claim_readiness", {}),
                "split_manifest_input": trained.get("split_manifest_input"),
                "sample_manifest": (dataset_provenance or {}).get("sample_manifest"),
            }
        )
        if run_audit:
            audit = run_ml_audit_suite(
                dataset,
                label_col,
                expert_dir / "audit",
                feature_set=feature_set,
                splits=[
                    "random",
                    "drug_holdout",
                    "target_holdout",
                    "scaffold_holdout",
                    "chemical_cluster_holdout",
                    "target_family_holdout",
                    "temporal_holdout",
                    "source_holdout",
                ],
                source_col="label_source",
                seed=seed,
                allow_partial_rescoring_features=allow_partial_rescoring_features,
                allow_label_definition_features=allow_label_definition_features,
                candidates=candidates,
                strict_source_overlap=True,
                source_transfer_train=source_transfer_train,
                source_transfer_test=source_transfer_test,
                predictions_path=model_dir / "model_predictions.csv",
                claim_mode=claim_mode,
            )
            result["audit"] = audit.get("outputs", {})
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = str(exc)
    return result


def _add_override(parser: argparse.ArgumentParser, expert: str) -> None:
    parser.add_argument(f"--{expert}-dataset", type=Path, default=None)
    parser.add_argument(f"--{expert}-label", default=EXPERT_DEFAULTS[expert]["label"])
    parser.add_argument(f"--{expert}-feature-set", default=EXPERT_DEFAULTS[expert]["feature_set"])
    parser.add_argument(f"--{expert}-split", default=EXPERT_DEFAULTS[expert]["split"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run a reproducible Atlas ML pass under one output directory. "
            "Only --out-dir is required; expert datasets default to files under data/<run-id>/."
        )
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", default="pilotstudy")
    parser.add_argument("--run-dir", type=Path, default=None)
    parser.add_argument(
        "--experts",
        nargs="*",
        choices=sorted(EXPERT_DEFAULTS),
        default=["binding", "exposure", "tissue", "mechanism"],
    )
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument("--model-n-jobs", type=int, default=0, help="Threads/jobs for model families that support parallelism; 0 keeps model defaults.")
    parser.add_argument("--bootstraps", type=int, default=200)
    parser.add_argument("--class-weight", choices=["balanced", "none"], default="balanced")
    parser.add_argument(
        "--pu-mode",
        choices=[
            "standard_binary",
            "positive_unlabeled_weighted",
            "case_control_matched",
            "bagging_pu",
            "stratified_bagging_pu",
            "propensity_weighted_pu",
            "elkan_noto",
            "pulsnar_style",
        ],
        default="standard_binary",
    )
    parser.add_argument("--pu-bags", type=int, default=50)
    parser.add_argument("--pu-unlabeled-ratio", type=float, default=1.0)
    parser.add_argument("--skip-audit", action="store_true")
    parser.add_argument("--run-ablations", action="store_true")
    parser.add_argument(
        "--run-estimated-pk-exposure",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Explicit opt-in: predict free Cmax from RDKit descriptors and combine with predicted AC50 for sensitivity analysis.",
    )
    parser.add_argument("--exposure-potency-feature-set", default="spd_potency_physchem")
    parser.add_argument("--exposure-pk-feature-set", default="ligand_physchem_descriptors")
    parser.add_argument("--exposure-regressor", choices=["ridge", "random_forest"], default="ridge")
    parser.add_argument(
        "--feature-metadata",
        choices=["auto", "always", "never"],
        default="auto",
        help="Refresh chemical_cluster, target_family, and source-lineage columns before training.",
    )
    parser.add_argument("--chemical-cluster", choices=["auto", "scaffold", "smiles", "chemotype", "none"], default="auto")
    parser.add_argument("--target-family", choices=["auto", "protein_class", "gene_heuristic", "none"], default="auto")
    parser.add_argument("--source-lineage", choices=["auto", "none"], default="auto")
    parser.add_argument(
        "--drop-column",
        nargs="*",
        default=None,
        help="Optional columns to remove from refreshed training tables for sensitivity passes.",
    )
    parser.add_argument(
        "--audit-candidate",
        nargs="*",
        default=[],
        help="Optional benchmark/calibration tables as NAME=PATH or PATH.",
    )
    parser.add_argument("--source-transfer-train", nargs="*", default=None)
    parser.add_argument("--source-transfer-test", nargs="*", default=None)
    parser.add_argument(
        "--prepare-spd-expert-tables",
        choices=["auto", "always", "never"],
        default="auto",
        help="Build objective-specific SPD binding/exposure/tissue/mechanism tables before selecting expert datasets.",
    )
    parser.add_argument(
        "--spd-expert-source-table",
        type=Path,
        default=None,
        help="Optional full-SPD Atlas table used to prepare four expert tables.",
    )
    parser.add_argument(
        "--spd-panel-table",
        type=Path,
        default=None,
        help="Optional SPD assay-pair table used to enrich raw Atlas run masters; defaults to local staged SPD panel files.",
    )
    parser.add_argument(
        "--spd-target-map",
        type=Path,
        default=None,
        help="Optional PDB-to-SPD-target map; defaults to local spd_target_selected.csv files.",
    )
    parser.add_argument(
        "--fda-ligand-map",
        type=Path,
        default=None,
        help="Optional ligand_base-to-FDA-drug map; defaults to chemdb/data/fda_mapping_from_pdbqt.csv.",
    )
    parser.add_argument(
        "--spd-target-metadata",
        type=Path,
        default=None,
        help="Optional SPD target metadata with protein class; defaults to analysis/gene_list/spd_assays_all.csv.",
    )
    parser.add_argument(
        "--mechanism-label-table",
        type=Path,
        default=None,
        help="Optional mechanism four-state label table to project onto this run; defaults to local staged mechanism labels.",
    )
    parser.add_argument(
        "--target-expression-table",
        type=Path,
        default=None,
        help="Optional HPA/GTEx/Bgee/OpenTargets-style target expression table for tissue scoring.",
    )
    parser.add_argument(
        "--allow-partial-rescoring-features",
        action="store_true",
        help="Allow SCORCH/final fields only for fully rescored subset datasets.",
    )
    parser.add_argument(
        "--allow-label-definition-features",
        action="store_true",
        help="Allow label-defining evidence fields for explicitly named sensitivity models only.",
    )
    parser.add_argument(
        "--nested-model-selection",
        action="store_true",
        default=True,
        help="Use an inner validation split for supported hyperparameter selection before outer-holdout scoring.",
    )
    parser.add_argument(
        "--no-nested-model-selection",
        dest="nested_model_selection",
        action="store_false",
        help="Disable the default inner validation split for a faster smoke run.",
    )
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    parser.add_argument("--sample-rows", type=int, default=0, help="Deterministically train on at most this many rows per expert for smoke iteration.")
    parser.add_argument("--sample-seed", type=int, default=42)
    parser.add_argument("--sample-strategy", choices=["stratified", "random"], default="stratified")
    parser.add_argument("--split-manifest", type=Path, default=None, help="Locked split manifest directory/CSV/JSON to reuse for every expert.")
    parser.add_argument("--expert-split-manifest", nargs="*", default=None, help="Per-expert locked splits as expert=PATH.")
    parser.add_argument("--no-mlflow", action="store_true", help="Disable default local MLflow logging for this training pass.")
    parser.add_argument("--no-feature-cache", dest="feature_cache", action="store_false", default=True)
    parser.add_argument("--refresh-feature-cache", action="store_true")
    parser.add_argument("--banana-scoring", choices=["auto", "always", "never"], default="auto", help="Populate BANANA scores for feature sets that request BANANA. Default scores labelable rows only.")
    parser.add_argument("--banana-score-scope", choices=["labelable", "all"], default="labelable")
    parser.add_argument("--banana-batch-size", type=int, default=512)
    parser.add_argument("--banana-num-workers", type=int, default=8)
    parser.add_argument(
        "--run-mechanism-pu-evaluation",
        action="store_true",
        help="Build mechanism PU labels, top-K recovery, and optional leave-one-source-out PU reports from the raw mechanism expert table.",
    )
    parser.add_argument(
        "--run-tissue-pu-recovery",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="For the tissue expert, run PU top-K recovery instead of ordinary supervised AUROC/AUPRC.",
    )
    parser.add_argument(
        "--tissue-pu-mode",
        choices=["bagging_pu", "stratified_bagging_pu"],
        default="stratified_bagging_pu",
    )
    parser.add_argument("--tissue-pu-bags", type=int, default=50)
    parser.add_argument("--tissue-pu-unlabeled-ratio", type=float, default=1.0)
    parser.add_argument("--tissue-pu-top-k", nargs="+", type=int, default=[10, 20, 50])
    parser.add_argument("--tissue-pu-max-eval-background", type=int, default=10000)
    parser.add_argument(
        "--allow-tissue-expression-features",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Allow expression/site features in tissue PU recovery when tissue labels are externally sourced.",
    )
    parser.add_argument("--skip-mechanism-source-held-pu", action="store_true")
    parser.add_argument("--mechanism-pu-max-sources", type=int, default=0)
    parser.add_argument("--calibration", choices=["none", "raw", "sigmoid", "isotonic"], default="none")
    parser.add_argument(
        "--group-reweight-cols",
        nargs="*",
        default=None,
        help="Optional capped inverse-frequency reweighting over source/family strata; use auto for source/family defaults.",
    )
    parser.add_argument(
        "--group-reweight-include-label",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include the supervised label in group-reweight strata.",
    )
    parser.add_argument("--group-reweight-max-factor", type=float, default=5.0)
    parser.add_argument("--strict-feature-set", action="store_true", help="Fail expert training if the selected feature set is not fully present after metadata refresh.")
    for expert in sorted(EXPERT_DEFAULTS):
        _add_override(parser, expert)
    args = parser.parse_args(argv)

    repo_root = args.repo_root.resolve()
    run_dir = _resolve_run_dir(repo_root, args.run_id, args.run_dir)
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    _configure_mlflow_defaults(out_dir=out, run_id=str(args.run_id), disabled=bool(args.no_mlflow))
    expert_split_manifests = _split_manifest_overrides(args.expert_split_manifest)
    config = load_config(args.config)
    seed = int(config.get("project", {}).get("random_seed", 42))
    class_weight = None if args.class_weight == "none" else args.class_weight
    model_params: dict[str, object] | None = None
    if int(args.model_n_jobs) > 0:
        model_params = {"n_jobs": int(args.model_n_jobs), "thread_count": int(args.model_n_jobs)}

    prepared_spd_tables = _prepare_spd_expert_tables(
        run_dir=run_dir,
        out_dir=out,
        mode=args.prepare_spd_expert_tables,
        source_table=args.spd_expert_source_table,
        spd_panel_path=args.spd_panel_table,
        target_map_path=args.spd_target_map,
        ligand_map_path=args.fda_ligand_map,
        target_metadata_path=args.spd_target_metadata,
        mechanism_label_path=args.mechanism_label_table,
        target_expression_path=args.target_expression_table,
    )

    results: list[dict[str, Any]] = []
    tissue_pu_outputs: dict[str, object] = {}
    for expert in args.experts:
        defaults = EXPERT_DEFAULTS[expert]
        dataset = (
            getattr(args, f"{expert}_dataset")
            or prepared_spd_tables.get(expert)
            or _first_existing(run_dir, list(defaults["candidates"]))
        )
        if dataset is None:
            results.append(
                {
                    "expert": expert,
                    "status": "skipped",
                    "reason": "default_dataset_not_found",
                    "searched": [str(run_dir / rel) for rel in defaults["candidates"]],
                }
            )
            continue
        label_col = getattr(args, f"{expert}_label")
        feature_set = getattr(args, f"{expert}_feature_set")
        split_mode = getattr(args, f"{expert}_split")
        if args.feature_metadata in {"auto", "always"}:
            dataset = _metadata_refreshed_dataset(
                dataset,
                run_dir=run_dir,
                out_dir=out,
                mode=args.feature_metadata,
                chemical_cluster=args.chemical_cluster,
                target_family=args.target_family,
                source_lineage=args.source_lineage,
                drop_columns=args.drop_column,
                use_cache=bool(args.feature_cache),
                refresh_cache=bool(args.refresh_feature_cache),
            )
        dataset = _maybe_populate_banana_scores(
            dataset,
            expert=expert,
            label_col=label_col,
            feature_set=feature_set,
            run_dir=run_dir,
            out_dir=out,
            repo_root=repo_root,
            mode=args.banana_scoring,
            scope=args.banana_score_scope,
            batch_size=int(args.banana_batch_size),
            num_workers=int(args.banana_num_workers),
        )
        split_manifest = expert_split_manifests.get(expert) or args.split_manifest
        group_reweight_cols = _resolve_group_reweight_cols(dataset, args.group_reweight_cols)
        dataset_provenance: dict[str, object] = {}
        if args.sample_rows and int(args.sample_rows) > 0:
            if split_manifest is not None:
                sample_manifest = materialize_split_locked_sampled_dataset(
                    dataset,
                    label_col=label_col,
                    split_manifest=split_manifest,
                    out_dir=out / "smoke_samples" / expert,
                    sample_rows=int(args.sample_rows),
                    seed=int(args.sample_seed),
                    strategy=str(args.sample_strategy),
                    split_mode=split_mode,
                )
                split_manifest = Path(str(sample_manifest["sampled_split_manifest"]))
            else:
                sample_manifest = materialize_sampled_dataset(
                    dataset,
                    label_col=label_col,
                    out_dir=out / "smoke_samples" / expert,
                    sample_rows=int(args.sample_rows),
                    seed=int(args.sample_seed),
                    strategy=str(args.sample_strategy),
                )
            dataset = Path(str(sample_manifest["sampled_dataset"]))
            dataset_provenance = {"sampled": True, "sample_manifest": sample_manifest}
        audit_candidates = [
            *list(args.audit_candidate),
            *_existing_audit_candidates(run_dir, list(defaults.get("audit_candidates", []))),
        ]
        if expert == "tissue" and args.run_tissue_pu_recovery:
            tissue_model_dir = out / "tissue" / "model"
            try:
                recovery = run_tissue_pu_recovery_report(
                    dataset,
                    label_col=label_col,
                    feature_set=feature_set,
                    model_type=args.model,
                    split_mode=split_mode,
                    out_dir=tissue_model_dir,
                    seed=seed,
                    pu_mode=args.tissue_pu_mode,
                    pu_bags=int(args.tissue_pu_bags),
                    pu_unlabeled_ratio=float(args.tissue_pu_unlabeled_ratio),
                    top_k=[int(k) for k in args.tissue_pu_top_k],
                    max_eval_background=int(args.tissue_pu_max_eval_background),
                    class_weight=class_weight,
                    allow_tissue_expression_features=bool(args.allow_tissue_expression_features),
                    strict_feature_set=bool(args.strict_feature_set),
                    repo_root=repo_root,
                )
                tissue_pu_outputs = {
                    "status": "written",
                    "model_dir": str(tissue_model_dir),
                    "manifest": str(tissue_model_dir / "tissue_pu_recovery_manifest.json"),
                    "predictions": str(tissue_model_dir / "tissue_pu_predictions.csv"),
                    "topk_recovery": str(tissue_model_dir / "tissue_pu_topk_recovery.csv"),
                    "case_studies": str(tissue_model_dir / "tissue_pu_case_studies.csv"),
                    "n_predictions": recovery.get("n_predictions"),
                }
                results.append(
                    {
                        "expert": expert,
                        "dataset": str(dataset),
                        "label": label_col,
                        "feature_set": feature_set,
                        "split": split_mode,
                        "status": "trained",
                        "model_dir": str(tissue_model_dir),
                        "metrics": {
                            "precision_at_K": recovery.get("record", {}).get("precision_at_K"),
                            "enrichment_at_K": recovery.get("record", {}).get("enrichment_at_K"),
                            "heldout_positive_recovery_at_K": recovery.get("record", {}).get("heldout_positive_recovery_at_K"),
                        },
                        "n_train": recovery.get("record", {}).get("n_train"),
                        "n_test": recovery.get("record", {}).get("n_test"),
                        "features": recovery.get("manifest", {}).get("features", []),
                    }
                )
            except Exception as exc:
                tissue_pu_outputs = {"status": "failed", "error": str(exc), "dataset": str(dataset)}
                results.append(
                    {
                        "expert": expert,
                        "dataset": str(dataset),
                        "label": label_col,
                        "feature_set": feature_set,
                        "split": split_mode,
                        "status": "failed",
                        "error": str(exc),
                    }
                )
            continue
        results.append(
            _run_expert(
                name=expert,
                dataset=dataset,
                label_col=label_col,
                feature_set=feature_set,
                split_mode=split_mode,
                out_dir=out,
                seed=seed,
                model_type=args.model,
                bootstraps=max(0, args.bootstraps),
                run_audit=not bool(args.skip_audit),
                class_weight=class_weight,
                allow_partial_rescoring_features=bool(args.allow_partial_rescoring_features),
                allow_label_definition_features=bool(args.allow_label_definition_features),
                nested_model_selection=bool(args.nested_model_selection),
                pu_mode=args.pu_mode,
                pu_bags=args.pu_bags,
                pu_unlabeled_ratio=args.pu_unlabeled_ratio,
                candidates=audit_candidates,
                source_transfer_train=args.source_transfer_train,
                source_transfer_test=args.source_transfer_test,
                claim_mode=args.claim_mode,
                repo_root=repo_root,
                split_manifest=split_manifest,
                dataset_provenance=dataset_provenance,
                calibration_method=args.calibration,
                strict_feature_set=bool(args.strict_feature_set),
                model_params=model_params,
                group_reweight_cols=group_reweight_cols,
                group_reweight_include_label=bool(args.group_reweight_include_label),
                group_reweight_max_factor=float(args.group_reweight_max_factor),
            )
        )
        if expert == "exposure" and args.run_estimated_pk_exposure:
            results.append(
                _run_estimated_pk_exposure(
                    dataset=dataset,
                    label_col=label_col,
                    potency_feature_set=args.exposure_potency_feature_set,
                    pk_feature_set=args.exposure_pk_feature_set,
                    split_mode=split_mode,
                    out_dir=out,
                    seed=seed,
                    model_type=args.exposure_regressor,
                )
            )
        if args.run_ablations and expert in {"binding", "exposure"}:
            ablation_feature_set = "pilot_binding_only"
            results.append(
                _run_expert(
                    name=f"{expert}_ablation_atlas_consensus_only",
                    dataset=dataset,
                    label_col=label_col,
                    feature_set=ablation_feature_set,
                    split_mode=split_mode,
                    out_dir=out,
                    seed=seed,
                    model_type=args.model,
                    bootstraps=max(0, args.bootstraps),
                    run_audit=False,
                    class_weight=class_weight,
                    allow_partial_rescoring_features=False,
                    allow_label_definition_features=False,
                    nested_model_selection=bool(args.nested_model_selection),
                    pu_mode=args.pu_mode,
                    pu_bags=args.pu_bags,
                    pu_unlabeled_ratio=args.pu_unlabeled_ratio,
                    candidates=[],
                    source_transfer_train=None,
                    source_transfer_test=None,
                    claim_mode=args.claim_mode,
                    repo_root=repo_root,
                    calibration_method=args.calibration,
                    strict_feature_set=bool(args.strict_feature_set),
                    model_params=model_params,
                )
            )

    mechanism_pu_outputs: dict[str, object] = {}
    if args.run_mechanism_pu_evaluation:
        mechanism_raw_candidates = [
            out / "_spd_expert_tables" / "ml_spd_mechanism_table.csv",
            run_dir / "spd_four_experts" / "ml_spd_mechanism_table.csv",
            run_dir / "_spd_expert_tables" / "ml_spd_mechanism_table.csv",
        ]
        mechanism_raw = next((path for path in mechanism_raw_candidates if path.exists()), None)
        if mechanism_raw is None:
            mechanism_pu_outputs = {
                "status": "skipped",
                "reason": "raw_mechanism_table_not_found",
                "searched": [str(path) for path in mechanism_raw_candidates],
            }
        else:
            mechanism_out = out / "mechanism_pu_evaluation"
            mechanism_out.mkdir(parents=True, exist_ok=True)
            pu_table = mechanism_out / "mechanism_pu_table.csv"
            build_mechanism_pu_table(mechanism_raw, pu_table)
            full_feature_pu_table: Path | None = None
            if args.feature_metadata in {"auto", "always"}:
                full_feature_pu_table = _metadata_refreshed_dataset(
                    pu_table,
                    run_dir=run_dir,
                    out_dir=mechanism_out,
                    mode=args.feature_metadata,
                    chemical_cluster=args.chemical_cluster,
                    target_family=args.target_family,
                    source_lineage=args.source_lineage,
                    drop_columns=args.drop_column,
                    use_cache=False,
                    refresh_cache=True,
                )
                explicit_full_feature = mechanism_out / "mechanism_pu_table.full_features.csv"
                if full_feature_pu_table.exists():
                    pd.read_csv(full_feature_pu_table, low_memory=False).to_csv(explicit_full_feature, index=False)
                    pu_table = explicit_full_feature
            topk = write_mechanism_topk_recovery(
                pu_table,
                mechanism_out / "mechanism_topk_recovery.csv",
                label_col="mechanism_pu_label",
            )
            loso_rows = 0
            if not args.skip_mechanism_source_held_pu:
                loso = run_leave_one_source_out(
                    pu_table,
                    mechanism_out,
                    source_col="label_source",
                    label_col="mechanism_pu_label",
                    feature_set="spd_mechanism_nonleaky",
                    model_type=args.model,
                    pu_mode=args.pu_mode,
                    seed=seed,
                    max_sources=int(args.mechanism_pu_max_sources or 0),
                    strict_feature_set=bool(args.strict_feature_set),
                )
                loso_rows = int(len(loso))
            mechanism_pu_outputs = {
                "status": "written",
                "raw_mechanism_table": str(mechanism_raw),
                "pu_table": str(pu_table),
                "full_feature_pu_table": str(full_feature_pu_table) if full_feature_pu_table else None,
                "topk_recovery": str(mechanism_out / "mechanism_topk_recovery.csv"),
                "topk_rows": int(len(topk)),
                "leave_one_source_out_rows": loso_rows,
            }

    manifest = {
        "run_id": args.run_id,
        "run_dir": str(run_dir),
        "out_dir": str(out),
        "feature_metadata": {
            "mode": args.feature_metadata,
            "chemical_cluster": args.chemical_cluster,
            "target_family": args.target_family,
            "source_lineage": args.source_lineage,
            "drop_columns": args.drop_column or [],
            "feature_cache": bool(args.feature_cache),
            "refresh_feature_cache": bool(args.refresh_feature_cache),
        },
        "banana_scoring": {
            "mode": args.banana_scoring,
            "scope": args.banana_score_scope,
            "batch_size": int(args.banana_batch_size),
            "num_workers": int(args.banana_num_workers),
        },
        "group_reweighting": {
            "group_reweight_cols": args.group_reweight_cols or [],
            "include_label": bool(args.group_reweight_include_label),
            "max_factor": float(args.group_reweight_max_factor),
        },
        "sampling": {
            "sample_rows": int(args.sample_rows or 0),
            "sample_seed": int(args.sample_seed),
            "sample_strategy": args.sample_strategy,
        },
        "tissue_pu_recovery": tissue_pu_outputs,
        "split_manifest": str(args.split_manifest) if args.split_manifest else None,
        "expert_split_manifests": {key: str(value) for key, value in expert_split_manifests.items()},
        "prepared_spd_expert_tables": {key: str(value) for key, value in prepared_spd_tables.items()},
        "mechanism_pu_evaluation": mechanism_pu_outputs,
        "experts": results,
    }
    manifest_path = out / "ml_training_pass_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    trained = [row for row in results if row.get("status") == "trained"]
    failed = [row for row in results if row.get("status") == "failed"]
    return 1 if failed or not trained else 0


if __name__ == "__main__":
    raise SystemExit(main())
