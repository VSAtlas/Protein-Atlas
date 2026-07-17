"""Dock selected ligands in a frozen reference-run Vina score space.

This module intentionally uses raw stage-1 Vina energies. Atlas consensus scores
are within-library percentiles, so percentiles from a small add-on library are
not comparable with percentiles from a full FDA/DUD run.
"""

from __future__ import annotations

import hashlib
import json
import math
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist
from typing import Any, Mapping

import numpy as np
import pandas as pd

from analysis.ml.compare_run_consensus import (
    _read_comparison_run,
    compare_consensus_to_run_decoys,
)
from config.runtime_config import load_config
from docking.run_vina import run_docking_task
from src.config.output_paths import run_output_dir

SCORE_COLUMN = "z_vs_compare_run_vina_stage1"
SCORCH_SCORE_COLUMN = "z_vs_compare_run_scorch_composite"
MIN_SCORCH_DECOY_NULL = 200
MIN_FULL_SCORCH_NULL_FRACTION = 0.95
MIN_REFERENCE_EXHAUSTIVENESS = 2
PROVENANCE_SUFFIXES = (
    "source",
    "run_id",
    "run_sha256",
    "null_sha256",
    "decoy_n",
    "decoy_mu",
    "decoy_sigma",
    "score_space_sha256",
    "receptor_sha256",
    "ligand_pdbqt_sha256",
)
_PATH_KEYS = {"receptor", "ligand", "out", "cpu"}


@dataclass(frozen=True)
class ReferenceContext:
    comparison_run_id: str
    pdb_id: str
    receptor: Path
    receptor_sha256: str
    template: Path
    parameters: dict[str, str]
    decoys: np.ndarray
    null_path: Path
    null_sha256: str
    run_sha256: str
    score_space_sha256: str
    config_count: int
    config_excluded_count: int
    config_context_sha256: str


def score_selected_against_reference_vina(
    *,
    selected_pairs_path: str | Path,
    comparison_run_id: str,
    out_dir: str | Path,
    repo_root: str | Path = ".",
    vina_exe: str = "vina",
    workers: int = 8,
    minimum_reference_exhaustiveness: int = MIN_REFERENCE_EXHAUSTIVENESS,
    run_id: str | None = None,
) -> dict[str, Any]:
    """Dock selected pairs with frozen stage-1 settings and standardize to DUD."""

    if not 1 <= int(workers) <= 32:
        raise ValueError("workers must be between 1 and 32")
    if int(minimum_reference_exhaustiveness) < 2:
        raise ValueError("minimum_reference_exhaustiveness must be at least 2")
    if run_id is not None:
        _validate_run_id(run_id)

    selected_path = Path(selected_pairs_path).resolve()
    selected_sha256 = _sha256_file(selected_path)
    output = Path(out_dir).resolve()
    root = Path(repo_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    workers = int(workers)

    selected = pd.read_csv(selected_path, low_memory=False)
    _require_columns(selected, {"pdb_id", "ligand_base", "pdbqt_path"})
    selected = _deduplicate_selected(selected)
    run_roots = {
        name: run_output_dir(root, name, comparison_run_id)
        for name in ("configs", "docked", "processed_pdbs")
    }
    missing_roots = [name for name, path in run_roots.items() if not path.is_dir()]
    if missing_roots:
        raise FileNotFoundError(
            f"comparison run {comparison_run_id!r} is missing roots: {missing_roots}"
        )

    resolved_vina = _resolve_executable(vina_exe)
    vina_version = _vina_version(resolved_vina)
    vina_sha256 = _sha256_file(resolved_vina)
    contexts = {
        pdb_id: _load_reference_context(
            pdb_id=pdb_id,
            comparison_run_id=comparison_run_id,
            run_roots=run_roots,
            vina_version=vina_version,
            vina_sha256=vina_sha256,
            minimum_exhaustiveness=int(minimum_reference_exhaustiveness),
        )
        for pdb_id in sorted(selected["pdb_id"].map(_pdb_id).unique())
    }
    reference_settings_audit = _reference_vina_settings_audit(contexts)
    guarded_inputs = {
        selected_path: selected_sha256,
        **{
            path: _sha256_file(path)
            for context in contexts.values()
            for path in (context.template, context.receptor, context.null_path)
        },
    }

    work_items: list[tuple[int, dict[str, Any], ReferenceContext]] = []
    immediate_rows: list[dict[str, Any]] = []
    for index, row in selected.reset_index(drop=True).iterrows():
        values = row.to_dict()
        pdb_id = _pdb_id(values.get("pdb_id"))
        ligand = Path(str(values.get("pdbqt_path") or "")).expanduser()
        context = contexts[pdb_id]
        if not ligand.is_file():
            immediate_rows.append(
                _failed_row(values, context, "missing_ligand_pdbqt", ligand)
            )
            continue
        work_items.append((index, values, context))

    rows = list(immediate_rows)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                _score_one,
                index=index,
                row=row,
                context=context,
                output=output,
                vina_exe=str(resolved_vina),
            ): index
            for index, row, context in work_items
        }
        for future in as_completed(futures):
            rows.append(future.result())

    scores = pd.DataFrame(rows).sort_values(["pdb_id", "ligand_base"], kind="stable")
    scores = _attach_bh_qvalues(scores)
    score_path = output / "reference_vina_scores.csv"
    scores.to_csv(score_path, index=False)

    context_rows = [
        {
            "pdb_id": context.pdb_id,
            "comparison_run_id": comparison_run_id,
            "reference_receptor": str(context.receptor),
            "reference_receptor_sha256": context.receptor_sha256,
            "reference_template": str(context.template),
            "reference_run_sha256": context.run_sha256,
            "reference_null_sha256": context.null_sha256,
            "reference_decoy_n": int(context.decoys.size),
            "reference_decoy_mu": float(context.decoys.mean()),
            "reference_decoy_sigma": float(context.decoys.std(ddof=0)),
            "reference_decoy_median": float(np.median(context.decoys)),
            "reference_decoy_mad": float(
                np.median(np.abs(context.decoys - np.median(context.decoys)))
            ),
            "reference_decoy_unique_scores": int(np.unique(context.decoys).size),
            "reference_config_count": context.config_count,
            "reference_config_excluded_count": context.config_excluded_count,
            "reference_config_context_sha256": context.config_context_sha256,
            "score_space_sha256": context.score_space_sha256,
            "parameters_json": json.dumps(context.parameters, sort_keys=True),
        }
        for context in contexts.values()
    ]
    context_path = output / "reference_vina_contexts.csv"
    pd.DataFrame(context_rows).to_csv(context_path, index=False)

    _verify_hashes_unchanged(guarded_inputs)
    status_counts = scores["reference_vina_status"].value_counts().to_dict()
    manifest: dict[str, Any] = {
        "schema": "atlas.reference-vina-compare.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_id": str(run_id or ""),
        "selected_pairs": str(selected_path),
        "selected_pairs_sha256": selected_sha256,
        "comparison_run_id": comparison_run_id,
        "repo_root": str(root),
        "run_roots": {name: str(path) for name, path in run_roots.items()},
        "vina_executable": str(resolved_vina),
        "vina_executable_sha256": vina_sha256,
        "vina_version": vina_version,
        "workers": workers,
        "minimum_reference_exhaustiveness_policy": MIN_REFERENCE_EXHAUSTIVENESS,
        "reference_vina_settings_audit": reference_settings_audit,
        "warnings": list(reference_settings_audit["warnings"]),
        "counts": {
            "selected_pairs": int(len(selected)),
            "scored_pairs": int(scores[SCORE_COLUMN].notna().sum()),
            "failed_pairs": int(scores[SCORE_COLUMN].isna().sum()),
            "targets": int(len(contexts)),
            "q_bh_global_le_0_05": int(
                pd.to_numeric(scores["reference_vina_q_bh_global"], errors="coerce")
                .le(0.05)
                .sum()
            ),
            "q_bh_target_le_0_05": int(
                pd.to_numeric(scores["reference_vina_q_bh_target"], errors="coerce")
                .le(0.05)
                .sum()
            ),
            "status": {str(key): int(value) for key, value in status_counts.items()},
        },
        "outputs": {
            "scores": str(score_path),
            "contexts": str(context_path),
        },
        "policy": {
            "score_direction": "higher modified z/tail z means stronger binding; raw Vina is lower-better",
            "z_formula": "0.67448975 * (median_DUD_stage1 - candidate_vina_stage1) / MAD_DUD_stage1",
            "empirical_p_formula": "(1 + count(DUD_score <= candidate_score)) / (n_DUD + 1)",
            "bh_q_global_family": (
                "Benjamini-Hochberg across all successfully scored selected pairs"
            ),
            "bh_q_target_family": (
                "Benjamini-Hochberg separately within each PDB target"
            ),
            "empirical_tail_z_formula": "Phi^-1(1 - clipped_empirical_p)",
            "null_membership": (
                "exact comparison run, dud_stage1, unique decoys_fda_* ligand, "
                "finite raw score; custom pose-valid flags are ignored on both sides"
            ),
            "common_score_space": (
                "candidate uses frozen reference receptor, grid, exhaustiveness, "
                "energy range, and num_modes"
            ),
            "tiny_library_percentiles_used": False,
            "scorch_interpretation": "not computed; SCORCH remains a selected-tranche conditional analysis",
            "random_seed": "reference run and candidate configs do not establish a fixed seed",
            "historical_vina_version": (
                "not recorded in the imported run; candidate executable/version "
                "are hashed and this limitation is explicit"
            ),
        },
    }
    manifest_path = output / "reference_vina_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def plan_reference_addon_scoring(
    *,
    selected_pairs_path: str | Path,
    comparison_run_id: str,
    out_dir: str | Path,
    run_id: str,
    repo_root: str | Path = ".",
    vina_exe: str = "vina",
    workers: int = 8,
    minimum_reference_exhaustiveness: int = MIN_REFERENCE_EXHAUSTIVENESS,
    scorch: bool = False,
) -> dict[str, Any]:
    """Validate immutable inputs and return a workload-free execution plan."""

    if not 1 <= int(workers) <= 32:
        raise ValueError("workers must be between 1 and 32")
    if int(minimum_reference_exhaustiveness) < 2:
        raise ValueError("minimum_reference_exhaustiveness must be at least 2")
    _validate_run_id(run_id)
    if scorch and run_id == comparison_run_id:
        raise ValueError("SCORCH add-on run_id must differ from comparison_run_id")

    selected_path = Path(selected_pairs_path).resolve()
    selected = pd.read_csv(selected_path, low_memory=False)
    _require_columns(selected, {"pdb_id", "ligand_base", "pdbqt_path"})
    selected = _deduplicate_selected(selected)
    root = Path(repo_root).resolve()
    run_roots = {
        name: run_output_dir(root, name, comparison_run_id)
        for name in ("configs", "docked", "processed_pdbs")
    }
    missing_roots = [name for name, path in run_roots.items() if not path.is_dir()]
    if missing_roots:
        raise FileNotFoundError(
            f"comparison run {comparison_run_id!r} is missing roots: {missing_roots}"
        )

    resolved_vina = _resolve_executable(vina_exe)
    vina_version = _vina_version(resolved_vina)
    vina_sha256 = _sha256_file(resolved_vina)
    targets = sorted(selected["pdb_id"].map(_pdb_id).unique())
    contexts = {
        pdb_id: _load_reference_context(
            pdb_id=pdb_id,
            comparison_run_id=comparison_run_id,
            run_roots=run_roots,
            vina_version=vina_version,
            vina_sha256=vina_sha256,
            minimum_exhaustiveness=int(minimum_reference_exhaustiveness),
        )
        for pdb_id in targets
    }
    reference_settings_audit = _reference_vina_settings_audit(contexts)
    missing_ligands = [
        str(path)
        for path in selected["pdbqt_path"].map(
            lambda value: Path(str(value)).expanduser()
        )
        if not path.is_file()
    ]
    plan: dict[str, Any] = {
        "schema": "atlas.reference-addon-plan.v1",
        "status": "planned",
        "run_id": run_id,
        "selected_pairs": str(selected_path),
        "selected_pairs_sha256": _sha256_file(selected_path),
        "comparison_run_id": comparison_run_id,
        "out_dir": str(Path(out_dir).resolve()),
        "workers": int(workers),
        "minimum_reference_exhaustiveness_policy": MIN_REFERENCE_EXHAUSTIVENESS,
        "reference_vina_settings_audit": reference_settings_audit,
        "warnings": list(reference_settings_audit["warnings"]),
        "scorch": bool(scorch),
        "candidate_scorch_coverage": 1.0 if scorch else 0.0,
        "counts": {
            "selected_pairs": int(len(selected)),
            "targets": int(len(targets)),
            "missing_ligand_pdbqts": int(len(missing_ligands)),
        },
        "targets": targets,
        "missing_ligand_pdbqts": missing_ligands,
        "vina_executable": str(resolved_vina),
        "vina_executable_sha256": vina_sha256,
        "vina_version": vina_version,
        "reference_contexts": [
            {
                "pdb_id": context.pdb_id,
                "receptor": str(context.receptor),
                "receptor_sha256": context.receptor_sha256,
                "stage1_decoy_n": int(context.decoys.size),
                "stage1_null_sha256": context.null_sha256,
                "score_space_sha256": context.score_space_sha256,
                "reference_exhaustiveness": _parameter_int(
                    context.parameters, "exhaustiveness"
                ),
                "parameters": dict(context.parameters),
            }
            for context in contexts.values()
        ],
    }
    if scorch:
        comparison_post_root, scorch_nulls = _load_scorch_nulls(
            root, comparison_run_id, targets
        )
        plan["comparison_post_docked_root"] = str(comparison_post_root)
        scorch_rows = [
            _scorch_null_manifest_row(
                null,
                stage1_decoy_n=int(contexts[pdb_id].decoys.size),
            )
            for pdb_id, null in scorch_nulls.items()
        ]
        plan["scorch_nulls"] = scorch_rows
        truncated = [
            row
            for row in scorch_rows
            if float(row["stage1_coverage_fraction"]) < MIN_FULL_SCORCH_NULL_FRACTION
        ]
        plan["scorch_reference_ready"] = True
        plan["scorch_reference_full_coverage"] = not truncated
        plan["scorch_null_selection_conditioned"] = bool(truncated)
        plan["scorch_score_claim_status"] = (
            "selection_conditioned_sensitivity"
            if truncated
            else "full_reference_null"
        )
        if truncated:
            plan["warnings"].append(
                "SCORCH DUD nulls cover only the Stage-1-selected subset. "
                "Docking and rescoring may continue, but standardized SCORCH "
                "scores are selection-conditioned sensitivity evidence and are "
                "not claim-grade until the same-PDB DUD null is fully rescored."
            )
        plan["scorch_workspace_policy"] = (
            "isolated under out_dir/_scorch_workspace; comparison run is read-only"
        )
    return plan


def score_reference_vina_candidates_with_scorch(
    *,
    stage1_manifest: Mapping[str, Any],
    selected_pairs_path: str | Path,
    comparison_run_id: str,
    out_dir: str | Path,
    run_id: str,
    repo_root: str | Path = ".",
    workers: int = 8,
    require_full_scorch_null: bool = False,
) -> dict[str, Any]:
    """SCORCH every candidate pose and standardize against frozen DUD SCORCH."""

    if not 1 <= int(workers) <= 32:
        raise ValueError("workers must be between 1 and 32")
    _validate_run_id(run_id)
    if run_id == comparison_run_id:
        raise ValueError("SCORCH add-on run_id must differ from comparison_run_id")

    root = Path(repo_root).resolve()
    output = Path(out_dir).resolve()
    selected_path = Path(selected_pairs_path).resolve()
    expected_selected_hash = str(stage1_manifest.get("selected_pairs_sha256") or "")
    if not expected_selected_hash:
        raise ValueError("stage-1 manifest lacks selected_pairs_sha256")
    if _sha256_file(selected_path) != expected_selected_hash:
        raise RuntimeError("selected-pair manifest changed after stage-1 planning")

    stage1_outputs = stage1_manifest.get("outputs")
    if not isinstance(stage1_outputs, Mapping):
        raise ValueError("stage-1 manifest lacks outputs")
    stage1_score_path = Path(str(stage1_outputs.get("scores") or "")).resolve()
    context_path = Path(str(stage1_outputs.get("contexts") or "")).resolve()
    if not stage1_score_path.is_file() or not context_path.is_file():
        raise FileNotFoundError("stage-1 score/context artifacts are missing")

    scores = pd.read_csv(stage1_score_path, low_memory=False)
    _require_columns(
        scores,
        {
            "pdb_id",
            "ligand_base",
            "reference_vina_status",
            "reference_vina_stage1_kcal_mol",
            "reference_vina_pose",
            f"{SCORE_COLUMN}_ligand_pdbqt_sha256",
        },
    )
    scores["pdb_id"] = scores["pdb_id"].map(_pdb_id)
    scores["ligand_base"] = scores["ligand_base"].map(_ligand_base)
    if scores.duplicated(["pdb_id", "ligand_base"]).any():
        raise ValueError("stage-1 scores contain duplicate PDB-ligand keys")
    raw_vina = pd.to_numeric(scores["reference_vina_stage1_kcal_mol"], errors="coerce")
    failed = ~scores["reference_vina_status"].eq(
        "scored_common_reference_context"
    ) | ~np.isfinite(raw_vina)
    if failed.any():
        raise RuntimeError(
            f"{int(failed.sum())} candidates lack valid stage-1 Vina scores"
        )

    contexts = pd.read_csv(context_path, low_memory=False)
    _require_columns(
        contexts,
        {
            "pdb_id",
            "reference_receptor",
            "reference_receptor_sha256",
            "reference_decoy_n",
            "score_space_sha256",
        },
    )
    contexts["pdb_id"] = contexts["pdb_id"].map(_pdb_id)
    if contexts["pdb_id"].duplicated().any():
        raise ValueError("stage-1 contexts contain duplicate target rows")
    targets = sorted(scores["pdb_id"].unique())
    if set(contexts["pdb_id"]) != set(targets):
        raise ValueError("stage-1 target/context coverage is mixed or incomplete")

    comparison_post_root, scorch_nulls = _load_scorch_nulls(
        root, comparison_run_id, targets
    )
    stage1_decoy_counts = {
        str(row["pdb_id"]): int(row["reference_decoy_n"])
        for row in contexts.to_dict(orient="records")
    }
    scorch_null_rows = [
        _scorch_null_manifest_row(
            null,
            stage1_decoy_n=stage1_decoy_counts[pdb_id],
        )
        for pdb_id, null in scorch_nulls.items()
    ]
    truncated_nulls = [
        row
        for row in scorch_null_rows
        if float(row["stage1_coverage_fraction"]) < MIN_FULL_SCORCH_NULL_FRACTION
    ]
    if require_full_scorch_null:
        _require_full_scorch_nulls(scorch_nulls, stage1_decoy_counts)
    selection_conditioned = bool(truncated_nulls)
    score_claim_status = (
        "selection_conditioned_sensitivity"
        if selection_conditioned
        else "full_reference_null"
    )
    null_warning = (
        "SCORCH DUD null covers only the Stage-1-selected subset; z-scores are "
        "selection-conditioned sensitivity evidence, not claim-grade."
        if selection_conditioned
        else ""
    )
    guarded_inputs: dict[Path, str] = {
        selected_path: expected_selected_hash,
        stage1_score_path: _sha256_file(stage1_score_path),
        context_path: _sha256_file(context_path),
    }
    for path in _scorch_null_source_paths(scorch_nulls):
        guarded_inputs[path] = _sha256_file(path)
    for row in contexts.to_dict(orient="records"):
        receptor = Path(str(row["reference_receptor"])).resolve()
        expected = str(row["reference_receptor_sha256"])
        if not receptor.is_file() or _sha256_file(receptor) != expected:
            raise RuntimeError(f"reference receptor hash mismatch: {receptor}")
        guarded_inputs[receptor] = expected
    for row in scores.to_dict(orient="records"):
        pose = Path(str(row["reference_vina_pose"])).resolve()
        if not pose.is_file():
            raise FileNotFoundError(f"candidate Vina pose is missing: {pose}")
        guarded_inputs[pose] = _sha256_file(pose)

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    workspace = output / "_scorch_workspace" / f"{run_id}_{timestamp}"
    docked_root = workspace / "docked"
    post_docked_root = workspace / "post_docked"
    processed_root = workspace / "processed_pdbs"
    docked_run = docked_root / run_id
    processed_run = processed_root / run_id
    context_by_pdb = contexts.set_index("pdb_id").to_dict(orient="index")

    for pdb_id, group in scores.groupby("pdb_id", sort=True):
        target_root = docked_run / pdb_id
        stage_root = target_root / "stage1"
        stage_root.mkdir(parents=True, exist_ok=False)
        receptor_source = Path(
            str(context_by_pdb[pdb_id]["reference_receptor"])
        ).resolve()
        receptor_target = processed_run / pdb_id / "receptor" / f"{pdb_id}.pdbqt"
        _copy_verified(receptor_source, receptor_target)

        rank = raw_vina.loc[group.index].rank(method="average", ascending=True)
        denominator = max(1, len(group) - 1)
        consensus = 1.0 - (rank - 1.0) / denominator
        consensus_rows: list[dict[str, Any]] = []
        for index, row in group.iterrows():
            ligand_base = _ligand_base(row["ligand_base"])
            pose_source = Path(str(row["reference_vina_pose"])).resolve()
            pose_target = stage_root / f"{ligand_base}_stage1.pdbqt"
            _copy_verified(pose_source, pose_target)
            percentile = float(consensus.loc[index])
            consensus_rows.append(
                {
                    "run_id": run_id,
                    "run_mode": "fda",
                    "pdb_id": pdb_id,
                    "variant": "",
                    "ph_label": "",
                    "ligand": f"{ligand_base}.pdbqt",
                    "library": "ml_addon",
                    "is_decoy": 0,
                    "consensus_score": percentile,
                    "p_energy": percentile,
                    "p_ledock": 0.0,
                    "p_cnn": 0.0,
                    "p_vina": percentile,
                    "p_gnina_energy": 0.0,
                    "p_dock6": 0.0,
                    "n_engines_with_data": 1,
                    "best_engine": "vina",
                    "best_signal": "vina",
                    "z_vs_decoys_consensus": "",
                    "t_vs_decoys_consensus": "",
                }
            )
        pd.DataFrame(consensus_rows).to_csv(
            target_root / "consensus_docking_scores.csv", index=False
        )

    _verify_hashes_unchanged(guarded_inputs)
    cfg = load_config(config_path=str(root / "config.txt"), base_dir=root)
    cfg.update(
        {
            "OVERALL_DIR": str(root),
            "OUTPUT_DIR": str(processed_root),
            "USE_SCORCH": True,
            "SCORCH_TOP_FRACTION": 1.0,
            "USE_GNINA": False,
            "USE_LEDOCK": False,
            "USE_DOCK6": False,
            "CPU": int(workers),
            "GLOBAL_DOCKING_SCHEDULER": "",
        }
    )
    from post_docking.rescoring.rescoring_scorch import main as scorch_main

    scorch_return_code = int(
        scorch_main(
            [
                "--run-id",
                run_id,
                "--repo-root",
                str(root),
                "--docked-root",
                str(docked_root),
                "--post-docked-root",
                str(post_docked_root),
                "--jobs",
                str(int(workers)),
                "--threads",
                "1",
                "--overwrite",
                "--skip-autofix",
                "--decoy-prefix",
                "dud",
            ],
            cfg_override=cfg,
        )
        or 0
    )
    if scorch_return_code != 0:
        raise RuntimeError(
            f"SCORCH engine failed for isolated add-on run: {scorch_return_code}"
        )

    post_run = post_docked_root / run_id
    result_frames: list[pd.DataFrame] = []
    for pdb_id in targets:
        result_path = post_run / pdb_id / "consensus_reranked_scorch.csv"
        if not result_path.is_file():
            raise FileNotFoundError(f"missing candidate SCORCH output: {result_path}")
        frame = pd.read_csv(result_path, low_memory=False)
        _require_columns(frame, {"ligand", "scorch_composite"})
        frame["pdb_id"] = pdb_id
        frame["ligand_base"] = frame["ligand"].map(_ligand_base)
        frame["scorch_output_file"] = str(result_path)
        frame["scorch_output_sha256"] = _sha256_file(result_path)
        result_frames.append(frame)
    scorch_scores = pd.concat(result_frames, ignore_index=True)
    if scorch_scores.duplicated(["pdb_id", "ligand_base"]).any():
        raise ValueError("SCORCH output contains duplicate PDB-ligand keys")

    expected_keys = set(zip(scores["pdb_id"], scores["ligand_base"], strict=False))
    observed_keys = set(
        zip(scorch_scores["pdb_id"], scorch_scores["ligand_base"], strict=False)
    )
    if expected_keys != observed_keys:
        missing = sorted(expected_keys - observed_keys)
        extra = sorted(observed_keys - expected_keys)
        raise RuntimeError(
            "SCORCH candidate coverage mismatch: "
            f"missing={missing[:8]}, extra={extra[:8]}"
        )
    composite = pd.to_numeric(scorch_scores["scorch_composite"], errors="coerce")
    if not np.isfinite(composite).all():
        raise RuntimeError(
            f"{int((~np.isfinite(composite)).sum())} candidates lack SCORCH scores"
        )

    stage1_only = [
        column
        for column in scores.columns
        if column not in scorch_scores.columns or column in {"pdb_id", "ligand_base"}
    ]
    combined = scorch_scores.merge(
        scores[stage1_only],
        on=["pdb_id", "ligand_base"],
        how="inner",
        validate="one_to_one",
    )
    coverage_by_pdb = {
        str(row["pdb_id"]): float(row["stage1_coverage_fraction"])
        for row in scorch_null_rows
    }
    combined["scorch_null_stage1_coverage_fraction"] = combined["pdb_id"].map(
        coverage_by_pdb
    )
    combined["scorch_null_selection_conditioned"] = selection_conditioned
    combined["scorch_score_claim_status"] = score_claim_status
    combined["scorch_reference_warning"] = null_warning
    raw_path = output / "reference_scorch_scores_raw.csv"
    combined.to_csv(raw_path, index=False)
    score_path = output / "reference_scorch_scores.csv"
    comparison_manifest = compare_consensus_to_run_decoys(
        raw_path,
        comparison_post_root,
        score_path,
        comparison_run_id=comparison_run_id,
        compare_stream="post_docked_scorch",
        compare_score_column="scorch_composite",
    )
    standardized = pd.read_csv(score_path, low_memory=False)
    z_scores = pd.to_numeric(standardized[SCORCH_SCORE_COLUMN], errors="coerce")
    if len(standardized) != len(scores) or not np.isfinite(z_scores).all():
        raise RuntimeError("SCORCH comparison z-score coverage is incomplete")
    _verify_hashes_unchanged(guarded_inputs)

    manifest_path = output / "reference_scorch_manifest.json"
    manifest: dict[str, Any] = {
        "schema": "atlas.reference-scorch-compare.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "run_id": run_id,
        "comparison_run_id": comparison_run_id,
        "selected_pairs": str(selected_path),
        "selected_pairs_sha256": expected_selected_hash,
        "workers": int(workers),
        "candidate_scorch_coverage": 1.0,
        "candidate_count": int(len(scores)),
        "target_count": int(len(targets)),
        "scorch_engine_return_code": scorch_return_code,
        "isolated_workspace": str(workspace),
        "comparison_run_mutated": False,
        "scorch_null_minimum_n": MIN_SCORCH_DECOY_NULL,
        "scorch_reference_full_coverage": not selection_conditioned,
        "scorch_null_selection_conditioned": selection_conditioned,
        "scorch_score_claim_status": score_claim_status,
        "warnings": [null_warning] if null_warning else [],
        "scorch_nulls": scorch_null_rows,
        "policy": {
            "candidate_selection": "all candidate reference_vina_pose rows",
            "candidate_top_fraction": 1.0,
            "scorch_source": "Vina pose only",
            "z_formula": (
                "z_vs_compare_run_scorch_composite = "
                "(candidate scorch_composite - same-PDB DUD SCORCH mean) / "
                "same-PDB DUD SCORCH population SD"
            ),
            "comparison_run_access": "read-only with before/after SHA256 guards",
            "score_interpretation": score_claim_status,
            "strict_full_null_requested": bool(require_full_scorch_null),
        },
        "hashes": {
            "stage1_scores_sha256": _sha256_file(stage1_score_path),
            "stage1_contexts_sha256": _sha256_file(context_path),
            "raw_scorch_scores_sha256": _sha256_file(raw_path),
            "scorch_scores_sha256": _sha256_file(score_path),
        },
        "outputs": {
            "scores": str(score_path),
            "raw_scores": str(raw_path),
            "comparison_manifest": str(comparison_manifest["outputs"]["manifest"]),
            "comparison_nulls": str(comparison_manifest["outputs"]["nulls"]),
            "manifest": str(manifest_path),
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def _load_scorch_nulls(
    repo_root: Path,
    comparison_run_id: str,
    targets: list[str],
) -> tuple[Path, dict[str, dict[str, Any]]]:
    comparison_root = run_output_dir(
        repo_root, "post_docked", comparison_run_id
    ).resolve()
    nulls, _ = _read_comparison_run(
        comparison_root,
        run_id=comparison_run_id,
        compare_stream="post_docked_scorch",
        compare_score_column="scorch_composite",
    )
    missing = sorted(set(targets) - set(nulls))
    if missing:
        raise ValueError(
            "comparison run lacks DUD SCORCH nulls for: " + ", ".join(missing)
        )
    selected = {pdb_id: nulls[pdb_id] for pdb_id in targets}
    for pdb_id, null in selected.items():
        n_decoys = int(null["n_decoys"])
        if n_decoys < MIN_SCORCH_DECOY_NULL:
            raise ValueError(
                f"insufficient DUD SCORCH null for {pdb_id}: "
                f"n={n_decoys}, required>={MIN_SCORCH_DECOY_NULL}"
            )
        if (
            str(null.get("variant") or "").strip()
            or str(null.get("ph_label") or "").strip()
        ):
            raise ValueError(
                f"comparison DUD SCORCH null for {pdb_id} uses a non-flat "
                "variant/pH context that cannot be mixed with this Vina reference"
            )
    return comparison_root, selected


def _require_full_scorch_nulls(
    nulls: Mapping[str, Mapping[str, Any]],
    stage1_decoy_counts: Mapping[str, int],
) -> None:
    truncated: list[str] = []
    for pdb_id, null in nulls.items():
        stage1_n = int(stage1_decoy_counts[pdb_id])
        scorch_n = int(null["n_decoys"])
        fraction = scorch_n / stage1_n if stage1_n else 0.0
        if fraction < MIN_FULL_SCORCH_NULL_FRACTION:
            truncated.append(f"{pdb_id}:{scorch_n}/{stage1_n} ({fraction:.1%})")
    if truncated:
        raise ValueError(
            "all-candidate SCORCH standardization requires a full same-PDB DUD "
            "SCORCH null; selection-truncated references found: " + ", ".join(truncated)
        )


def _reference_paths(
    null: Mapping[str, Any],
    field: str,
    *,
    required: bool,
) -> list[Path]:
    encoded = str(null.get(field) or "[]")
    values = json.loads(encoded)
    if not isinstance(values, list) or (required and not values):
        requirement = "non-empty list" if required else "list"
        raise ValueError(f"{field} must encode a {requirement}")
    paths: list[Path] = []
    for value in values:
        path = Path(str(value)).resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        paths.append(path)
    return paths


def _scorch_null_source_paths(
    nulls: Mapping[str, Mapping[str, Any]],
) -> list[Path]:
    paths: list[Path] = []
    for null in nulls.values():
        paths.extend(
            _reference_paths(null, "decoy_reference_file", required=True)
        )
        paths.extend(
            _reference_paths(null, "score_reference_file", required=False)
        )
    return sorted(set(paths))


def _scorch_null_manifest_row(
    null: Mapping[str, Any],
    *,
    stage1_decoy_n: int | None = None,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "pdb_id": str(null["pdb_id"]),
        "n_decoys": int(null["n_decoys"]),
        "n_unique": int(null["n_unique"]),
        "mu": float(null["mu"]),
        "sigma": float(null["sigma"]),
        "null_sha256": str(null["null_sha256"]),
        "decoy_reference_file": str(null["decoy_reference_file"]),
        "decoy_reference_sha256": str(null["decoy_reference_sha256"]),
    }
    for field in ("score_reference_file", "score_reference_sha256"):
        value = null.get(field)
        if value:
            row[field] = str(value)
    if stage1_decoy_n is not None:
        stage1_n = int(stage1_decoy_n)
        row["stage1_decoy_n"] = stage1_n
        row["stage1_coverage_fraction"] = (
            int(null["n_decoys"]) / stage1_n if stage1_n else 0.0
        )
    return row


def _copy_verified(source: Path, destination: Path) -> None:
    source_hash = _sha256_file(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    if _sha256_file(destination) != source_hash:
        raise RuntimeError(f"staged file hash mismatch: {destination}")


def _verify_hashes_unchanged(expected: Mapping[Path, str]) -> None:
    changed = [
        str(path)
        for path, digest in expected.items()
        if not path.is_file() or _sha256_file(path) != digest
    ]
    if changed:
        raise RuntimeError(
            "immutable scoring inputs changed during execution: "
            + ", ".join(changed[:8])
        )


def _validate_run_id(value: str) -> str:
    run_id = str(value).strip()
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    if (
        not run_id
        or run_id in {".", ".."}
        or any(character not in allowed for character in run_id)
    ):
        raise ValueError(
            "run_id must contain only letters, digits, dot, underscore, or hyphen"
        )
    return run_id


def _parameter_int(parameters: Mapping[str, Any], key: str) -> int | None:
    normalized = {
        str(name).strip().casefold(): value for name, value in parameters.items()
    }
    try:
        return int(float(normalized[key.casefold()]))
    except (KeyError, TypeError, ValueError):
        return None


def _reference_vina_settings_audit(
    contexts: Mapping[str, ReferenceContext],
) -> dict[str, Any]:
    per_target = [
        {
            "pdb_id": context.pdb_id,
            "exhaustiveness": _parameter_int(
                context.parameters, "exhaustiveness"
            ),
            "num_modes": _parameter_int(context.parameters, "num_modes"),
            "energy_range": _parameter_int(context.parameters, "energy_range"),
            "config_count": int(context.config_count),
            "config_excluded_count": int(context.config_excluded_count),
            "config_context_sha256": context.config_context_sha256,
        }
        for context in contexts.values()
    ]
    exhaustiveness = [
        int(row["exhaustiveness"])
        for row in per_target
        if row["exhaustiveness"] is not None
    ]
    warnings: list[str] = []
    if exhaustiveness and min(exhaustiveness) <= MIN_REFERENCE_EXHAUSTIVENESS:
        warnings.append(
            "The comparison run uses the minimum accepted Vina exhaustiveness "
            "for at least one target. Settings are reproduced exactly, but this "
            "does not establish docking convergence; keep the add-on analysis "
            "exploratory until an exhaustiveness/convergence sensitivity passes."
        )
    if len(set(exhaustiveness)) > 1:
        warnings.append(
            "Reference Vina exhaustiveness differs across targets; preserve "
            "target-stratified provenance and do not interpret z-scores as one "
            "homogeneous raw-energy scale."
        )
    return {
        "status": (
            "minimum_settings_no_convergence_evidence"
            if warnings
            else "settings_recorded"
        ),
        "policy_minimum_exhaustiveness": MIN_REFERENCE_EXHAUSTIVENESS,
        "observed_exhaustiveness_min": min(exhaustiveness)
        if exhaustiveness
        else None,
        "observed_exhaustiveness_max": max(exhaustiveness)
        if exhaustiveness
        else None,
        "heterogeneous_exhaustiveness": len(set(exhaustiveness)) > 1,
        "per_target": per_target,
        "warnings": warnings,
    }


def _attach_bh_qvalues(scores: pd.DataFrame) -> pd.DataFrame:
    """Attach BH q-values without converting failed or missing scores to evidence."""

    result = scores.copy()
    p_values = pd.to_numeric(result["reference_vina_empirical_p"], errors="coerce")
    result["reference_vina_q_bh_global"] = _benjamini_hochberg(p_values)
    result["reference_vina_q_bh_target"] = pd.Series(
        np.nan, index=result.index, dtype=float
    )
    for _, indices in result.groupby("pdb_id", dropna=False, sort=False).groups.items():
        index = pd.Index(indices)
        result.loc[index, "reference_vina_q_bh_target"] = _benjamini_hochberg(
            p_values.loc[index]
        )
    return result


def _benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Return monotone Benjamini-Hochberg adjusted p-values on the input index."""

    numeric = pd.to_numeric(p_values, errors="coerce")
    valid = numeric.notna() & numeric.between(0.0, 1.0, inclusive="both")
    output = pd.Series(np.nan, index=numeric.index, dtype=float)
    if not valid.any():
        return output
    ordered = numeric.loc[valid].sort_values(kind="stable")
    count = len(ordered)
    ranks = np.arange(1, count + 1, dtype=float)
    adjusted = ordered.to_numpy(dtype=float) * float(count) / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output.loc[ordered.index] = np.clip(adjusted, 0.0, 1.0)
    return output


def _score_one(
    *,
    index: int,
    row: dict[str, Any],
    context: ReferenceContext,
    output: Path,
    vina_exe: str,
) -> dict[str, Any]:
    ligand_base = _ligand_base(row.get("ligand_base"))
    ligand_path = Path(str(row.get("pdbqt_path"))).expanduser().resolve()
    ligand_sha256 = _sha256_file(ligand_path)
    config_path = output / "configs" / context.pdb_id / f"{ligand_base}_stage1.txt"
    pose_path = output / "poses" / context.pdb_id / f"{ligand_base}_stage1.pdbqt"
    _write_config(
        config_path,
        context=context,
        ligand=ligand_path,
        out=pose_path,
    )
    pose_path.unlink(missing_ok=True)
    pose_path.with_suffix(pose_path.suffix + ".failed.txt").unlink(missing_ok=True)
    _, score = run_docking_task(
        vina_exe,
        str(config_path),
        ligand_base,
        str(pose_path),
        write_failure_marker_flag=True,
        cpu_override=1,
    )
    if score is None or not math.isfinite(float(score)):
        return _failed_row(
            row,
            context,
            "vina_failed_or_missing_score",
            ligand_path,
            ligand_sha256=ligand_sha256,
            config_path=config_path,
            pose_path=pose_path,
        )

    raw_score = float(score)
    mu = float(context.decoys.mean())
    sigma = float(context.decoys.std(ddof=0))
    median = float(np.median(context.decoys))
    mad = float(np.median(np.abs(context.decoys - median)))
    if mad <= 0:
        return _failed_row(
            row,
            context,
            "reference_null_zero_mad",
            ligand_path,
            ligand_sha256=ligand_sha256,
            config_path=config_path,
            pose_path=pose_path,
        )
    z_score = 0.67448975 * (median - raw_score) / mad
    more_extreme = int(np.count_nonzero(context.decoys <= raw_score))
    empirical_p = (more_extreme + 1.0) / (context.decoys.size + 1.0)
    tail_probability = min(
        max(empirical_p, 1.0 / (context.decoys.size + 1.0)),
        context.decoys.size / (context.decoys.size + 1.0),
    )
    empirical_tail_z = NormalDist().inv_cdf(1.0 - tail_probability)
    result = dict(row)
    result.update(
        {
            "pdb_id": context.pdb_id,
            "ligand_base": ligand_base,
            "reference_vina_status": "scored_common_reference_context",
            "reference_vina_stage1_kcal_mol": raw_score,
            SCORE_COLUMN: z_score,
            "reference_vina_empirical_p": empirical_p,
            "reference_vina_empirical_tail_z": empirical_tail_z,
            "reference_vina_top10pct": bool(empirical_p <= 0.10),
            "reference_vina_config": str(config_path),
            "reference_vina_pose": str(pose_path),
            f"{SCORE_COLUMN}_source": "frozen_reference_raw_vina_stage1_dud_null",
            f"{SCORE_COLUMN}_run_id": _run_id_from_context(context),
            f"{SCORE_COLUMN}_run_sha256": context.run_sha256,
            f"{SCORE_COLUMN}_null_sha256": context.null_sha256,
            f"{SCORE_COLUMN}_decoy_n": int(context.decoys.size),
            f"{SCORE_COLUMN}_decoy_mu": mu,
            f"{SCORE_COLUMN}_decoy_sigma": sigma,
            f"{SCORE_COLUMN}_decoy_median": median,
            f"{SCORE_COLUMN}_decoy_mad": mad,
            f"{SCORE_COLUMN}_decoy_unique_scores": int(np.unique(context.decoys).size),
            f"{SCORE_COLUMN}_score_space_sha256": context.score_space_sha256,
            f"{SCORE_COLUMN}_receptor_sha256": context.receptor_sha256,
            f"{SCORE_COLUMN}_ligand_pdbqt_sha256": ligand_sha256,
        }
    )
    return result


def _load_reference_context(
    *,
    pdb_id: str,
    comparison_run_id: str,
    run_roots: dict[str, Path],
    vina_version: str,
    vina_sha256: str,
    minimum_exhaustiveness: int,
) -> ReferenceContext:
    stage_dir = run_roots["configs"] / pdb_id / "dud_stage1"
    templates = sorted(stage_dir.glob("decoys_fda_*_dud_stage1.txt"))
    if not templates:
        raise FileNotFoundError(
            f"missing reference DUD stage-1 Vina configs for {pdb_id}"
        )
    (
        template,
        parameters,
        config_context_sha256,
        allowed_decoy_ligands,
        config_excluded_count,
    ) = _validate_reference_configs(templates, pdb_id=pdb_id)
    missing = [
        key
        for key in (
            "center_x",
            "center_y",
            "center_z",
            "size_x",
            "size_y",
            "size_z",
            "exhaustiveness",
            "energy_range",
            "num_modes",
        )
        if key not in parameters
    ]
    if missing:
        raise ValueError(f"reference Vina config for {pdb_id} misses {missing}")
    exhaustiveness = int(float(parameters["exhaustiveness"]))
    if exhaustiveness < minimum_exhaustiveness:
        raise ValueError(
            f"reference Vina config for {pdb_id} uses exhaustiveness={exhaustiveness}; "
            f"ML add-on scoring requires >= {minimum_exhaustiveness}"
        )
    receptor = run_roots["processed_pdbs"] / pdb_id / "receptor" / f"{pdb_id}.pdbqt"
    if not receptor.is_file():
        holo = (
            run_roots["processed_pdbs"]
            / pdb_id
            / "HOLO"
            / "receptor"
            / f"{pdb_id}.pdbqt"
        )
        receptor = holo if holo.is_file() else receptor
    if not receptor.is_file():
        raise FileNotFoundError(f"missing frozen reference receptor for {pdb_id}")
    null_path = run_roots["docked"] / pdb_id / "dud_docking_score_long.csv"
    null_table = pd.read_csv(null_path, low_memory=False)
    _require_columns(null_table, {"run_id", "stage", "ligand", "score"})
    stage = null_table["stage"].astype(str).str.strip().str.casefold()
    run_match = (
        null_table["run_id"].fillna("").astype(str).str.strip().eq(comparison_run_id)
    )
    ligand = null_table["ligand"].fillna("").astype(str).str.strip()
    ligand_base = ligand.map(_ligand_base)
    mask = stage.eq("dud_stage1") & run_match & ligand_base.isin(allowed_decoy_ligands)
    null_members = null_table.loc[mask].copy()
    if null_members["ligand"].duplicated().any():
        duplicate_count = int(null_members["ligand"].duplicated(keep=False).sum())
        raise ValueError(
            f"reference DUD stage-1 null for {pdb_id} has "
            f"{duplicate_count} duplicate ligand rows"
        )
    decoys = pd.to_numeric(null_table.loc[mask, "score"], errors="coerce")
    decoys = decoys[np.isfinite(decoys)].to_numpy(dtype=float)
    mad = float(np.median(np.abs(decoys - np.median(decoys)))) if decoys.size else 0.0
    if decoys.size < 200 or mad <= 0:
        raise ValueError(
            f"invalid reference DUD stage-1 null for {pdb_id}: n={decoys.size}"
        )
    receptor_sha256 = _sha256_file(receptor)
    null_sha256 = _sha256_file(null_path)
    run_payload = {
        "comparison_run_id": comparison_run_id,
        "pdb_id": pdb_id,
        "reference_template_sha256": _sha256_file(template),
        "reference_config_context_sha256": config_context_sha256,
        "reference_config_count": len(allowed_decoy_ligands),
        "reference_config_excluded_count": config_excluded_count,
        "reference_null_sha256": null_sha256,
        "reference_receptor_sha256": receptor_sha256,
    }
    run_sha256 = _json_sha256(run_payload)
    score_payload = {
        "stage": "stage1",
        "engine": "vina",
        "parameters": {
            key: value for key, value in parameters.items() if key not in _PATH_KEYS
        },
        "config_context_sha256": config_context_sha256,
        "receptor_sha256": receptor_sha256,
        "vina_version": vina_version,
        "vina_executable_sha256": vina_sha256,
        "null_sha256": null_sha256,
    }
    return ReferenceContext(
        comparison_run_id=comparison_run_id,
        pdb_id=pdb_id,
        receptor=receptor,
        receptor_sha256=receptor_sha256,
        template=template,
        parameters=parameters,
        decoys=decoys,
        null_path=null_path,
        null_sha256=null_sha256,
        run_sha256=run_sha256,
        score_space_sha256=_json_sha256(score_payload),
        config_count=len(allowed_decoy_ligands),
        config_excluded_count=config_excluded_count,
        config_context_sha256=config_context_sha256,
    )


def _validate_reference_configs(
    templates: list[Path],
    *,
    pdb_id: str,
) -> tuple[Path, dict[str, str], str, set[str], int]:
    """Select the dominant score-independent DUD stage-1 config context."""

    contexts: dict[str, dict[str, Any]] = {}
    expected_receptor_suffix = f"/{pdb_id}/receptor/{pdb_id}.pdbqt".casefold()
    for template in templates:
        parameters = _parse_config(template)
        receptor = parameters.get("receptor", "").replace("\\", "/").casefold()
        if not receptor.endswith(expected_receptor_suffix):
            raise ValueError(
                f"reference config receptor does not match {pdb_id}: {template}"
            )
        context = {
            key: value
            for key, value in parameters.items()
            if key not in _PATH_KEYS and key != "verbosity"
        }
        signature = _json_sha256(context)
        group = contexts.setdefault(
            signature,
            {"parameters": parameters, "templates": []},
        )
        group["templates"].append(template)
    if not contexts:
        raise ValueError(f"no reference configs found for {pdb_id}")
    signature, modal = max(
        contexts.items(),
        key=lambda item: len(item[1]["templates"]),
    )
    modal_templates: list[Path] = modal["templates"]
    minimum_dominance = max(200, math.ceil(0.90 * len(templates)))
    if len(modal_templates) < minimum_dominance:
        counts = {key: len(value["templates"]) for key, value in contexts.items()}
        raise ValueError(
            f"reference DUD stage-1 configs for {pdb_id} lack a dominant "
            f"score context: {counts}"
        )
    allowed_ligands = {
        path.name.removesuffix("_dud_stage1.txt") for path in modal_templates
    }
    payload = {
        "pdb_id": pdb_id,
        "config_count": len(modal_templates),
        "config_excluded_count": len(templates) - len(modal_templates),
        "score_context_sha256": signature,
    }
    return (
        modal_templates[0],
        modal["parameters"],
        _json_sha256(payload),
        allowed_ligands,
        len(templates) - len(modal_templates),
    )


def _write_config(
    path: Path,
    *,
    context: ReferenceContext,
    ligand: Path,
    out: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"receptor = {context.receptor}",
        f"ligand = {ligand}",
        *(
            f"{key} = {value}"
            for key, value in context.parameters.items()
            if key not in _PATH_KEYS
        ),
        "cpu = 1",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _failed_row(
    row: dict[str, Any],
    context: ReferenceContext,
    status: str,
    ligand_path: Path,
    *,
    ligand_sha256: str = "",
    config_path: Path | None = None,
    pose_path: Path | None = None,
) -> dict[str, Any]:
    result = dict(row)
    result.update(
        {
            "pdb_id": context.pdb_id,
            "ligand_base": _ligand_base(row.get("ligand_base")),
            "reference_vina_status": status,
            "reference_vina_stage1_kcal_mol": pd.NA,
            SCORE_COLUMN: pd.NA,
            "reference_vina_empirical_p": pd.NA,
            "reference_vina_empirical_tail_z": pd.NA,
            "reference_vina_top10pct": pd.NA,
            "reference_vina_config": str(config_path or ""),
            "reference_vina_pose": str(pose_path or ""),
            f"{SCORE_COLUMN}_source": "",
            f"{SCORE_COLUMN}_run_id": _run_id_from_context(context),
            f"{SCORE_COLUMN}_run_sha256": context.run_sha256,
            f"{SCORE_COLUMN}_null_sha256": context.null_sha256,
            f"{SCORE_COLUMN}_decoy_n": int(context.decoys.size),
            f"{SCORE_COLUMN}_decoy_mu": float(context.decoys.mean()),
            f"{SCORE_COLUMN}_decoy_sigma": float(context.decoys.std(ddof=0)),
            f"{SCORE_COLUMN}_decoy_median": float(np.median(context.decoys)),
            f"{SCORE_COLUMN}_decoy_mad": float(
                np.median(np.abs(context.decoys - np.median(context.decoys)))
            ),
            f"{SCORE_COLUMN}_decoy_unique_scores": int(np.unique(context.decoys).size),
            f"{SCORE_COLUMN}_score_space_sha256": context.score_space_sha256,
            f"{SCORE_COLUMN}_receptor_sha256": context.receptor_sha256,
            f"{SCORE_COLUMN}_ligand_pdbqt_sha256": ligand_sha256,
            "pdbqt_path": str(ligand_path),
        }
    )
    return result


def _deduplicate_selected(frame: pd.DataFrame) -> pd.DataFrame:
    work = frame.copy()
    work["pdb_id"] = work["pdb_id"].map(_pdb_id)
    if "ligand_file_stem" in work:
        file_stem = work["ligand_file_stem"].fillna("").astype(str).str.strip()
        ligand_base = file_stem.where(file_stem.ne(""), work["ligand_base"])
    else:
        ligand_base = work["ligand_base"]
    work["ligand_base"] = ligand_base.map(_ligand_base)
    work["pdbqt_path"] = work["pdbqt_path"].fillna("").astype(str).str.strip()
    invalid = work["pdb_id"].eq("") | work["ligand_base"].eq("")
    if invalid.any():
        raise ValueError(
            f"selected pairs contain {int(invalid.sum())} invalid pair keys"
        )
    conflicts = work.groupby(["pdb_id", "ligand_base"])["pdbqt_path"].nunique().gt(1)
    if conflicts.any():
        raise ValueError(
            f"selected pairs contain conflicting ligand paths for {int(conflicts.sum())} keys"
        )
    return work.drop_duplicates(["pdb_id", "ligand_base"], keep="first")


def _parse_config(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip().casefold()] = value.strip()
    return values


def _resolve_executable(value: str) -> Path:
    candidate = Path(value).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    resolved = shutil.which(value)
    if not resolved:
        raise FileNotFoundError(f"Vina executable not found: {value}")
    return Path(resolved).resolve()


def _vina_version(executable: Path) -> str:
    result = subprocess.run(
        [str(executable), "--version"],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    text = (result.stdout or result.stderr or "").strip()
    return text.splitlines()[0] if text else "unknown"


def _run_id_from_context(context: ReferenceContext) -> str:
    return context.comparison_run_id


def _require_columns(frame: pd.DataFrame, required: set[str]) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"table missing required columns: {missing}")


def _pdb_id(value: Any) -> str:
    text = str(value or "").strip().upper()
    return "" if text.casefold() == "nan" else text


def _ligand_base(value: Any) -> str:
    text = Path(str(value or "").strip()).name
    if text.casefold().endswith(".pdbqt"):
        text = text[:-6]
    return "" if text.casefold() == "nan" else text


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


__all__ = [
    "SCORE_COLUMN",
    "plan_reference_addon_scoring",
    "score_reference_vina_candidates_with_scorch",
    "score_selected_against_reference_vina",
]
