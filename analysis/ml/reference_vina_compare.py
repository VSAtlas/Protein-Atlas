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
from typing import Any

import numpy as np
import pandas as pd

from docking.run_vina import run_docking_task
from src.config.output_paths import run_output_dir

SCORE_COLUMN = "z_vs_compare_run_vina_stage1"
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
    minimum_reference_exhaustiveness: int = 2,
) -> dict[str, Any]:
    """Dock selected pairs with frozen stage-1 settings and standardize to DUD."""

    if not 1 <= int(workers) <= 32:
        raise ValueError("workers must be between 1 and 32")
    if int(minimum_reference_exhaustiveness) < 2:
        raise ValueError("minimum_reference_exhaustiveness must be at least 2")

    selected_path = Path(selected_pairs_path).resolve()
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

    scores = pd.DataFrame(rows).sort_values(
        ["pdb_id", "ligand_base"], kind="stable"
    )
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

    status_counts = scores["reference_vina_status"].value_counts().to_dict()
    manifest: dict[str, Any] = {
        "schema": "atlas.reference-vina-compare.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "selected_pairs": str(selected_path),
        "selected_pairs_sha256": _sha256_file(selected_path),
        "comparison_run_id": comparison_run_id,
        "repo_root": str(root),
        "run_roots": {name: str(path) for name, path in run_roots.items()},
        "vina_executable": str(resolved_vina),
        "vina_executable_sha256": vina_sha256,
        "vina_version": vina_version,
        "workers": workers,
        "minimum_reference_exhaustiveness": int(
            minimum_reference_exhaustiveness
        ),
        "counts": {
            "selected_pairs": int(len(selected)),
            "scored_pairs": int(scores[SCORE_COLUMN].notna().sum()),
            "failed_pairs": int(scores[SCORE_COLUMN].isna().sum()),
            "targets": int(len(contexts)),
            "q_bh_global_le_0_05": int(
                pd.to_numeric(
                    scores["reference_vina_q_bh_global"], errors="coerce"
                )
                .le(0.05)
                .sum()
            ),
            "q_bh_target_le_0_05": int(
                pd.to_numeric(
                    scores["reference_vina_q_bh_target"], errors="coerce"
                )
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
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest["outputs"]["manifest"] = str(manifest_path)
    return manifest


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
            f"{SCORE_COLUMN}_decoy_unique_scores": int(
                np.unique(context.decoys).size
            ),
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
    ) = _validate_reference_configs(
        templates, pdb_id=pdb_id
    )
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
        holo = run_roots["processed_pdbs"] / pdb_id / "HOLO" / "receptor" / f"{pdb_id}.pdbqt"
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
    mask = (
        stage.eq("dud_stage1")
        & run_match
        & ligand_base.isin(allowed_decoy_ligands)
    )
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
        raise ValueError(f"invalid reference DUD stage-1 null for {pdb_id}: n={decoys.size}")
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
        "parameters": {key: value for key, value in parameters.items() if key not in _PATH_KEYS},
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
        counts = {
            key: len(value["templates"])
            for key, value in contexts.items()
        }
        raise ValueError(
            f"reference DUD stage-1 configs for {pdb_id} lack a dominant "
            f"score context: {counts}"
        )
    allowed_ligands = {
        path.name.removesuffix("_dud_stage1.txt")
        for path in modal_templates
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
            f"{SCORE_COLUMN}_decoy_unique_scores": int(
                np.unique(context.decoys).size
            ),
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
        raise ValueError(f"selected pairs contain {int(invalid.sum())} invalid pair keys")
    conflicts = (
        work.groupby(["pdb_id", "ligand_base"])["pdbqt_path"].nunique().gt(1)
    )
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


__all__ = ["SCORE_COLUMN", "score_selected_against_reference_vina"]
