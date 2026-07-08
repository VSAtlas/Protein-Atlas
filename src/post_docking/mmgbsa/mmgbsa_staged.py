from __future__ import annotations

import json
import logging
import os
import shutil
import shlex
import stat
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from post_docking.mmgbsa.mmgbsa_batch import (
    MMGBSAJob,
    _post_docked_base,
    _safe_ligand_stem,
    extract_artifact,
    group_hits,
    load_report_hits,
    run_jobs,
)
from post_docking.mmgbsa.mmgbsa_explicit import (
    build_and_run_explicit_workflow,
    run_tleap,
    write_solvated_tleap,
)
from post_docking.mmgbsa.mmgbsa_mpi_env import apply_blas_single_thread_env_defaults
from post_docking.mmgbsa.mmgbsa_pipeline import (
    _mmgbsa_append_summary,
    _mmgbsa_frame_aggregate,
    _mmgbsa_frame_qc_ok,
    _mmgbsa_md_aggregate,
    _mmgbsa_summary_frame_fields,
)
from post_docking.mmgbsa.mmgbsa_qc import detect_outliers, summarize_values
from post_docking.mmgbsa.run_mmgbsa import parse_mmpbsa_delta_total, run_mmgbsa

LOGGER = logging.getLogger("mmgbsa.staged")

VALIDATION_TIERS = (
    "custom",
    "contract",
    "topology",
    "tiny-md",
    "canary",
    "wang-lite",
    "production",
)


def _write_json(path: Path, payload: Mapping[str, Any] | Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".part")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    os.replace(tmp_path, path)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def _abs(path: str | Path, base: Path | None = None) -> str:
    value = Path(path).expanduser()
    if not value.is_absolute() and base is not None:
        value = base / value
    return str(value.resolve())


def _cfg_paths(cfg: Mapping[str, Any], job: MMGBSAJob) -> tuple[Path, Path, Path]:
    base = _post_docked_base(cfg, job) / (job.ph_label or "pH7_0")
    mmgbsa_root = base / "mmgbsa"
    work_stage = mmgbsa_root / "work" / job.stage_dir
    return base, mmgbsa_root, work_stage


def _parse_tleap_inputs(build_leap: Path) -> dict[str, str]:
    if not build_leap.exists():
        return {}
    inputs: dict[str, str] = {}
    for line in build_leap.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if "loadmol2" in stripped:
            inputs["ligand_mol2"] = stripped.split("loadmol2", 1)[1].strip()
        elif stripped.startswith("loadamberparams "):
            inputs["ligand_frcmod"] = stripped.split(None, 1)[1].strip()
        elif "loadpdb" in stripped:
            inputs["receptor_pdb"] = stripped.split("loadpdb", 1)[1].strip()
    return {key: _abs(value, build_leap.parent) for key, value in inputs.items()}


def _fallback_inputs(mmgbsa_root: Path, stage_dir: str, ligand_stem: str) -> dict[str, str]:
    return {
        "ligand_mol2": _abs(mmgbsa_root / "mol2" / stage_dir / f"{ligand_stem}.mol2"),
        "ligand_frcmod": _abs(mmgbsa_root / "frcmod" / stage_dir / f"{ligand_stem}.frcmod"),
    }


def _find_ligand_work_dir(work_stage: Path, ligand_id: str) -> Path | None:
    stem = _safe_ligand_stem(ligand_id)
    exact = work_stage / stem
    if (exact / "complex.prmtop").exists():
        return exact
    if not work_stage.exists():
        return None
    matches = [
        child
        for child in work_stage.iterdir()
        if child.is_dir()
        and child.name.startswith(stem)
        and (child / "complex.prmtop").exists()
    ]
    return sorted(matches, key=lambda path: path.name)[0] if matches else None


def _ligand_record(cfg: Mapping[str, Any], job: MMGBSAJob, ligand_id: str) -> dict[str, Any]:
    base_dir, mmgbsa_root, work_stage = _cfg_paths(cfg, job)
    work_dir = _find_ligand_work_dir(work_stage, ligand_id)
    ligand_stem = work_dir.name if work_dir is not None else _safe_ligand_stem(ligand_id)
    record: dict[str, Any] = {
        "run_id": job.run_id,
        "pdb_id": job.pdb_id,
        "pdb_file": job.pdb_file,
        "variant": job.variant or "",
        "ph_label": job.ph_label or "pH7_0",
        "stage_dir": job.stage_dir,
        "ligand_id": ligand_id,
        "ligand_stem": ligand_stem,
        "base_dir": _abs(base_dir),
        "mmgbsa_root": _abs(mmgbsa_root),
        "work_dir": _abs(work_dir) if work_dir is not None else "",
    }
    if work_dir is None:
        record["ready"] = False
        record["missing"] = ["ligand_work_dir"]
        return record

    inputs = _fallback_inputs(mmgbsa_root, job.stage_dir, ligand_stem)
    inputs.update(_parse_tleap_inputs(work_dir / "build.leap"))
    record.update(inputs)
    record.update(
        {
            "dry_complex_prmtop": _abs(work_dir / "complex.prmtop"),
            "dry_receptor_prmtop": _abs(work_dir / "receptor.prmtop"),
            "dry_ligand_prmtop": _abs(work_dir / "ligand.prmtop"),
            "ready": True,
        }
    )
    required = (
        "receptor_pdb",
        "ligand_mol2",
        "ligand_frcmod",
        "dry_complex_prmtop",
        "dry_receptor_prmtop",
        "dry_ligand_prmtop",
    )
    missing = [key for key in required if not Path(str(record.get(key, ""))).exists()]
    record["ready"] = not missing
    record["missing"] = missing
    return record


def _gpu_list(raw: str) -> list[str]:
    tokens = [part.strip() for part in str(raw or "").replace(";", ",").split(",")]
    return [token for token in tokens if token]


def _seed_for_replicate(base_seed: int, replicate: int) -> int:
    return max(1, int(base_seed) + (int(replicate) - 1) * 10007)


def _cap(value: int, *, default: int = 1, upper: int = 32) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        parsed = default
    return max(1, min(parsed, upper))


def _write_shell(path: Path, lines: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(lines) + "\n"
    tmp_path = path.with_suffix(path.suffix + ".part")
    tmp_path.write_text(text, encoding="utf-8")
    os.replace(tmp_path, path)
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP)


def _quote_paths(paths: Sequence[Path]) -> str:
    return " ".join(shlex.quote(str(path)) for path in paths)


def _write_scripts(
    *,
    out_dir: Path,
    repo_root: Path,
    python_exe: str,
    config: str,
    md_jobs: Sequence[Path],
    aggregate_jobs: Sequence[Path],
    max_parallel: int,
) -> dict[str, str]:
    config_arg = shlex.quote(config)
    py = shlex.quote(python_exe)
    repo = shlex.quote(str(repo_root))
    md_script = out_dir / "run_md_local.sh"
    agg_script = out_dir / "run_aggregate_local.sh"
    _write_shell(
        md_script,
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f"cd {repo}",
            f"MAX_PARALLEL=\"${{MMGBSA_STAGED_MAX_PARALLEL:-{max_parallel}}}\"",
            f"JOBS=({_quote_paths(md_jobs)})",
            'for job in "${JOBS[@]}"; do',
            '  while [ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL" ]; do wait -n; done',
            f"  {py} -m post_docking.mmgbsa.cli staged-run-md --job \"$job\" --config {config_arg} &",
            "done",
            "wait",
        ],
    )
    _write_shell(
        agg_script,
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f"cd {repo}",
            f"JOBS=({_quote_paths(aggregate_jobs)})",
            'for job in "${JOBS[@]}"; do',
            f"  {py} -m post_docking.mmgbsa.cli staged-aggregate --job \"$job\" --config {config_arg}",
            "done",
        ],
    )
    return {"md_local": str(md_script), "aggregate_local": str(agg_script)}


def _make_md_job(
    ligand: Mapping[str, Any],
    *,
    job_path: Path,
    replicate: int,
    seed: int,
    gpu: str,
    cpu_cap: int,
    md_ranks: int,
    md_engine: str,
    openmm_start_stage: str,
    prod_ps: float,
    frame_stride_ps: float,
    validation_tier: str,
    cfg_overrides: Mapping[str, Any],
    fixture_out_dir: Path | None,
) -> dict[str, Any]:
    work_dir = Path(str(ligand["work_dir"]))
    tier_slug = validation_tier.strip().lower().replace("-", "_")
    if tier_slug in {"custom", "production"}:
        rep_dir = work_dir / f"rep{replicate}"
    else:
        rep_dir = work_dir / f"rep{replicate}_{tier_slug}"
    payload = {
        "schema_version": 1,
        "kind": "mmgbsa_staged_md",
        "job_path": str(job_path),
        "replicate": replicate,
        "seed": seed,
        "gpu": gpu,
        "cpu_cap": cpu_cap,
        "md_mpi_ranks": md_ranks,
        "md_engine": md_engine,
        "openmm_start_stage": openmm_start_stage,
        "prod_ps": prod_ps,
        "frame_stride_ps": frame_stride_ps,
        "validation_tier": validation_tier,
        "cfg_overrides": dict(cfg_overrides),
        "fixture_out_dir": str(fixture_out_dir) if fixture_out_dir else "",
        "ligand": dict(ligand),
        "rep_dir": str(rep_dir),
        "explicit_dir": str(rep_dir / "explicit"),
        "status_json": str(rep_dir / "explicit" / "md_job_status.json"),
    }
    _write_json(job_path, payload)
    return payload


def _make_aggregate_job(
    ligand: Mapping[str, Any],
    *,
    job_path: Path,
    replicate_jobs: Sequence[Mapping[str, Any]],
    cpu_cap: int,
    mmpbsa_ranks: int,
) -> dict[str, Any]:
    work_dir = Path(str(ligand["work_dir"]))
    validation_tier = (
        str(replicate_jobs[0].get("validation_tier", "custom"))
        if replicate_jobs
        else "custom"
    )
    if validation_tier in {"contract", "topology"}:
        summary_json = job_path.with_suffix(".summary.json")
        status_json = job_path.with_suffix(".status.json")
    elif validation_tier not in {"custom", "production"}:
        tier_slug = validation_tier.replace("-", "_")
        summary_json = work_dir / f"mmgbsa_replicate_summary_{tier_slug}.json"
        status_json = work_dir / f"aggregate_status_{tier_slug}.json"
    else:
        summary_json = work_dir / "mmgbsa_replicate_summary.json"
        status_json = work_dir / "aggregate_status.json"
    payload = {
        "schema_version": 1,
        "kind": "mmgbsa_staged_aggregate",
        "job_path": str(job_path),
        "cpu_cap": cpu_cap,
        "mmpbsa_mpi_ranks": mmpbsa_ranks,
        "ligand": dict(ligand),
        "replicates": [dict(item) for item in replicate_jobs],
        "validation_tier": validation_tier,
        "summary_json": str(summary_json),
        "status_json": str(status_json),
    }
    _write_json(job_path, payload)
    return payload


def _prepare_cfg(cfg: Mapping[str, Any]) -> dict[str, Any]:
    prepared = dict(cfg)
    prepared.update(
        {
            "MMGBSA_ENABLED": True,
            "MMGBSA_MD_ENABLED": False,
            "MMGBSA_MD_SOLVENT_MODEL": "explicit",
            "MMGBSA_MMPBSA_ENABLED": False,
            "MMGBSA_MMPBSA_RUN": False,
            "MMGBSA_TLEAP_ENABLED": True,
            "MMGBSA_TLEAP_RUN": True,
        }
    )
    return prepared


def _prepare_report_jobs(
    *,
    report: Path,
    run_id: str,
    cfg: Mapping[str, Any],
    artifacts: Sequence[Path],
    artifact_dest: Path | None,
    default_stage_dir: str,
    default_variant: str,
    default_ph_label: str,
    require_significant: bool,
    max_q_value: float | None,
    score_column: str,
    score_direction: str,
    top_n: int,
    run_prepare: bool,
    test_mode: str,
    dry_run: bool,
) -> tuple[list[MMGBSAJob], list[dict[str, Any]]]:
    artifact_root = artifact_dest or Path(str(cfg.get("POST_DOCKED_DIR", "post_docked"))) / run_id
    for archive in artifacts:
        extract_artifact(archive, artifact_root)
    hits = load_report_hits(
        report,
        run_id=run_id,
        default_stage_dir=default_stage_dir,
        default_variant=default_variant,
        default_ph_label=default_ph_label,
        require_significant=require_significant,
        max_q_value=max_q_value,
        score_column=score_column,
        score_direction=score_direction,
        top_n=max(0, top_n),
    )
    jobs = group_hits(hits)
    if not run_prepare:
        return list(jobs), []
    prepare_results = run_jobs(
        jobs,
        cfg=_prepare_cfg(cfg),
        test_mode=test_mode,
        dry_run=dry_run,
    )
    return list(jobs), prepare_results


def _staged_caps(
    *,
    cfg: Mapping[str, Any],
    gpus: str,
    cpus: int,
    md_engine: str,
    openmm_start_stage: str,
    md_ranks: int,
    mmpbsa_ranks: int,
) -> dict[str, Any]:
    gpu_tokens = _gpu_list(gpus)
    max_parallel = max(1, min(len(gpu_tokens) if gpu_tokens else 1, _cap(cpus)))
    per_job_cpu = max(1, _cap(cpus) // max_parallel)
    return {
        "base_seed": int(float(cfg.get("MMGBSA_MD_BASE_SEED", 12345) or 12345)),
        "cpus_total": _cap(cpus),
        "gpus": gpu_tokens,
        "max_parallel_md_jobs": max_parallel,
        "per_md_job_cpu_cap": per_job_cpu,
        "md_mpi_ranks": _md_rank_cap(
            md_engine=md_engine,
            openmm_start_stage=openmm_start_stage,
            requested=md_ranks,
            per_job_cpu=per_job_cpu,
        ),
        "mmpbsa_mpi_ranks": _cap(mmpbsa_ranks, upper=_cap(cpus)),
    }


def _md_rank_cap(
    *,
    md_engine: str,
    openmm_start_stage: str,
    requested: int,
    per_job_cpu: int,
) -> int:
    if md_engine != "openmm":
        return _cap(requested, upper=per_job_cpu)
    start_stage = str(openmm_start_stage or "min1").strip().lower()
    if start_stage in {"min1", "minimize", "all"}:
        return 1
    hybrid_requested = per_job_cpu if int(requested or 0) <= 1 else requested
    return _cap(hybrid_requested, default=per_job_cpu, upper=per_job_cpu)


def _tier_cfg_overrides(validation_tier: str) -> dict[str, Any]:
    tier = validation_tier.strip().lower()
    if tier == "tiny-md":
        return {
            "MMGBSA_MD_DT_PS": 0.0005,
            "MMGBSA_MD_MIN_STEPS": 500,
            "MMGBSA_MD_RESCUE_STEPS": 2000,
            "MMGBSA_MD_HEAT_PS": 2.0,
            "MMGBSA_MD_EQUIL_PS": 2.0,
            "MMGBSA_MD_PROD_PS": 2.0,
            "MMGBSA_MD_FRAME_STRIDE_PS": 0.5,
            "MMGBSA_MMPBSA_USE_TRAJ_FRAMES": True,
        }
    if tier == "canary":
        return {
            "MMGBSA_MD_DT_PS": 0.0005,
            "MMGBSA_MD_MIN_STEPS": 5000,
            "MMGBSA_MD_HEAT_PS": 2.0,
            "MMGBSA_MD_EQUIL_PS": 50.0,
            "MMGBSA_MD_PROD_PS": 500.0,
            "MMGBSA_MD_FRAME_STRIDE_PS": 10.0,
            "MMGBSA_MMPBSA_USE_TRAJ_FRAMES": True,
        }
    if tier == "wang-lite":
        return {
            "MMGBSA_MD_DT_PS": 0.0005,
            "MMGBSA_MD_MIN_STEPS": 5000,
            "MMGBSA_MD_HEAT_PS": 2.0,
            "MMGBSA_MD_EQUIL_PS": 100.0,
            "MMGBSA_MD_PROD_PS": 1000.0,
            "MMGBSA_MD_FRAME_STRIDE_PS": 10.0,
            "MMGBSA_MMPBSA_USE_TRAJ_FRAMES": True,
            "MMGBSA_KEEP_WATERS": True,
            "MMGBSA_WATER_POLICY": "ACTIVE_SITE",
            "MMGBSA_WATER_KEEP_RADIUS_A": 4.0,
            "MMGBSA_WATER_KEEP_RADIUS": 4.0,
            "MMGBSA_TLEAP_WATER_MODEL": "tip3p",
            "MMGBSA_EXPLICIT_WATER_BUFFER_A": 8.0,
            "MMGBSA_QC_MIN_FRAMES": 50,
            "MMGBSA_QC_MIN_BLOCKS": 5,
            "MMGBSA_QC_MAX_SEM_KCAL": 2.0,
            "MMGBSA_QC_MAX_BLOCK_RANGE_KCAL": 12.0,
        }
    if tier == "production":
        return {
            "MMGBSA_MD_MIN_STEPS": 10000,
            "MMGBSA_MD_HEAT_PS": 100.0,
            "MMGBSA_MD_EQUIL_PS": 1000.0,
            "MMGBSA_MD_PROD_PS": 10000.0,
            "MMGBSA_MD_FRAME_STRIDE_PS": 10.0,
            "MMGBSA_MMPBSA_USE_TRAJ_FRAMES": True,
            "MMGBSA_KEEP_WATERS": True,
            "MMGBSA_WATER_POLICY": "ACTIVE_SITE",
            "MMGBSA_WATER_KEEP_RADIUS_A": 4.0,
            "MMGBSA_WATER_KEEP_RADIUS": 4.0,
            "MMGBSA_TLEAP_WATER_MODEL": "tip3p",
            "MMGBSA_EXPLICIT_WATER_BUFFER_A": 8.0,
        }
    return {}


def _tier_runtime_values(
    validation_tier: str,
    *,
    replicates: int,
    prod_ps: float,
    frame_stride_ps: float,
) -> tuple[int, float, float, dict[str, Any]]:
    tier = validation_tier.strip().lower()
    overrides = _tier_cfg_overrides(tier)
    if tier == "tiny-md":
        return 1, 2.0, 0.5, overrides
    if tier == "canary":
        return 1, 500.0, 10.0, overrides
    if tier == "wang-lite":
        return 1, 1000.0, 10.0, overrides
    if tier == "production":
        return max(3, replicates), 10000.0, 10.0, overrides
    return max(1, replicates), prod_ps, frame_stride_ps, overrides


def _write_staged_job_files(
    *,
    out_dir: Path,
    ligands: Sequence[Mapping[str, Any]],
    caps: Mapping[str, Any],
    replicates: int,
    prod_ps: float,
    frame_stride_ps: float,
    md_engine: str,
    openmm_start_stage: str,
    validation_tier: str,
    cfg_overrides: Mapping[str, Any],
    fixture_out_dir: Path | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[Path], list[Path]]:
    md_job_paths: list[Path] = []
    aggregate_job_paths: list[Path] = []
    md_payloads: list[dict[str, Any]] = []
    aggregate_payloads: list[dict[str, Any]] = []
    gpu_tokens = list(caps.get("gpus", []))
    base_seed = int(caps.get("base_seed", 12345))
    for ligand_idx, ligand in enumerate(ligands, start=1):
        replicate_payloads: list[dict[str, Any]] = []
        for rep in range(1, max(1, replicates) + 1):
            gpu = gpu_tokens[(ligand_idx + rep - 2) % len(gpu_tokens)] if gpu_tokens else ""
            job_path = out_dir / "jobs" / "md" / f"{ligand_idx:04d}_rep{rep}.json"
            payload = _make_md_job(
                ligand,
                job_path=job_path,
                replicate=rep,
                seed=_seed_for_replicate(base_seed + ligand_idx * 100_000, rep),
                gpu=gpu,
                cpu_cap=int(caps["per_md_job_cpu_cap"]),
                md_ranks=int(caps["md_mpi_ranks"]),
                md_engine=md_engine,
                openmm_start_stage=openmm_start_stage,
                prod_ps=prod_ps,
                frame_stride_ps=frame_stride_ps,
                validation_tier=validation_tier,
                cfg_overrides=cfg_overrides,
                fixture_out_dir=fixture_out_dir,
            )
            md_payloads.append(payload)
            replicate_payloads.append(payload)
            md_job_paths.append(job_path)
        aggregate_path = out_dir / "jobs" / "aggregate" / f"{ligand_idx:04d}.json"
        aggregate_payload = _make_aggregate_job(
            ligand,
            job_path=aggregate_path,
            replicate_jobs=replicate_payloads,
            cpu_cap=int(caps["cpus_total"]),
            mmpbsa_ranks=int(caps["mmpbsa_mpi_ranks"]),
        )
        aggregate_payloads.append(aggregate_payload)
        aggregate_job_paths.append(aggregate_path)
    return md_payloads, aggregate_payloads, md_job_paths, aggregate_job_paths


def stage_from_report(
    *,
    report: Path,
    run_id: str,
    out_dir: Path,
    cfg: Mapping[str, Any],
    config_path: str,
    artifacts: Sequence[Path] = (),
    artifact_dest: Path | None = None,
    default_stage_dir: str = "stage1",
    default_variant: str = "",
    default_ph_label: str = "",
    require_significant: bool = False,
    max_q_value: float | None = None,
    score_column: str = "",
    score_direction: str = "lower",
    top_n: int = 0,
    run_prepare: bool = True,
    test_mode: str = "off",
    dry_run: bool = False,
    replicates: int = 1,
    prod_ps: float = 4000.0,
    frame_stride_ps: float = 10.0,
    md_engine: str = "openmm",
    gpus: str = "0",
    cpus: int = 32,
    md_ranks: int = 1,
    mmpbsa_ranks: int = 32,
    openmm_start_stage: str = "heat",
    python_exe: str = "",
    validation_tier: str = "custom",
    fixture_out_dir: Path | None = None,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg_work = dict(cfg)
    tier = validation_tier.strip().lower()
    if tier not in VALIDATION_TIERS:
        raise ValueError(f"unsupported validation tier: {validation_tier}")
    replicates, prod_ps, frame_stride_ps, cfg_overrides = _tier_runtime_values(
        tier,
        replicates=replicates,
        prod_ps=prod_ps,
        frame_stride_ps=frame_stride_ps,
    )
    cfg_stage = dict(cfg_work)
    cfg_stage.update(cfg_overrides)
    jobs, prepare_results = _prepare_report_jobs(
        report=report,
        run_id=run_id,
        cfg=cfg_stage,
        artifacts=artifacts,
        artifact_dest=artifact_dest,
        default_stage_dir=default_stage_dir,
        default_variant=default_variant,
        default_ph_label=default_ph_label,
        require_significant=require_significant,
        max_q_value=max_q_value,
        score_column=score_column,
        score_direction=score_direction,
        top_n=top_n,
        run_prepare=run_prepare,
        test_mode=test_mode,
        dry_run=dry_run,
    )
    ligands = [
        _ligand_record(cfg_stage, job, ligand_id)
        for job in jobs
        for ligand_id in job.ligand_ids
    ]
    ready_ligands = [ligand for ligand in ligands if ligand.get("ready")]
    caps = _staged_caps(
        cfg=cfg_stage,
        gpus=gpus,
        cpus=cpus,
        md_engine=md_engine,
        openmm_start_stage=openmm_start_stage,
        md_ranks=md_ranks,
        mmpbsa_ranks=mmpbsa_ranks,
    )
    md_payloads, aggregate_payloads, md_job_paths, aggregate_job_paths = _write_staged_job_files(
        out_dir=out_dir,
        ligands=ready_ligands,
        caps=caps,
        replicates=replicates,
        prod_ps=prod_ps,
        frame_stride_ps=frame_stride_ps,
        md_engine=md_engine,
        openmm_start_stage=openmm_start_stage,
        validation_tier=tier,
        cfg_overrides=cfg_overrides,
        fixture_out_dir=fixture_out_dir,
    )

    scripts = _write_scripts(
        out_dir=out_dir,
        repo_root=Path.cwd(),
        python_exe=python_exe or sys.executable,
        config=config_path,
        md_jobs=md_job_paths,
        aggregate_jobs=aggregate_job_paths,
        max_parallel=int(caps["max_parallel_md_jobs"]),
    )
    manifest = {
        "schema_version": 1,
        "workflow": "mmgbsa_staged_report",
        "run_id": run_id,
        "report": str(report),
        "out_dir": str(out_dir),
        "prepare_ran": run_prepare,
        "prepare_results": prepare_results,
        "jobs": [asdict(job) for job in jobs],
        "ligands": ligands,
        "md_jobs": md_payloads,
        "aggregate_jobs": aggregate_payloads,
        "scripts": scripts,
        "caps": caps,
        "validation_tier": tier,
        "cfg_overrides": cfg_overrides,
        "fixture_out_dir": str(fixture_out_dir) if fixture_out_dir else "",
    }
    _write_json(out_dir / "staged_manifest.json", manifest)
    return manifest


def _job_cpu_cap(job: Mapping[str, Any]) -> int:
    return _cap(int(job.get("cpu_cap", 1)))


def _apply_thread_env(cpu_cap: int) -> None:
    env = os.environ
    env["GLOBAL_SCHEDULER_CPUS"] = str(cpu_cap)
    env["CPU"] = str(cpu_cap)
    apply_blas_single_thread_env_defaults(env)


def _apply_gpu_env(job: Mapping[str, Any], *, mmpbsa: bool) -> None:
    gpu = str(job.get("gpu", "") or "").strip()
    if not gpu or mmpbsa:
        return
    # ROCm recommends ROCR_VISIBLE_DEVICES for Linux GPU isolation. On this
    # MI100 node, setting HIP_VISIBLE_DEVICES/GPU_DEVICE_ORDINAL together with
    # ROCR_VISIBLE_DEVICES made nonzero GPU masks invisible to OpenMM/HIP.
    os.environ["ROCR_VISIBLE_DEVICES"] = gpu
    for env_name in ("HIP_VISIBLE_DEVICES", "GPU_DEVICE_ORDINAL", "CUDA_VISIBLE_DEVICES"):
        os.environ.pop(env_name, None)


def _openmm_logical_device_index_for_mask(gpu_token: str) -> str:
    token = str(gpu_token or "").strip()
    if not token or token.lower() == "all":
        return token
    return "0"


def _apply_job_caps(
    cfg: Mapping[str, Any], job: Mapping[str, Any], *, mmpbsa: bool
) -> dict[str, Any]:
    updated = dict(cfg)
    cpu_cap = _job_cpu_cap(job)
    updated["CPU"] = cpu_cap
    updated["GLOBAL_SCHEDULER_CPUS"] = cpu_cap
    if mmpbsa:
        updated["MMGBSA_MMPBSA_MPI_RANKS"] = _cap(
            int(job.get("mmpbsa_mpi_ranks", cpu_cap)), upper=cpu_cap
        )
    else:
        updated["MMGBSA_MD_MPI_RANKS"] = _cap(int(job.get("md_mpi_ranks", 1)), upper=cpu_cap)
        updated["MMGBSA_MD_ENGINE"] = str(job.get("md_engine", "openmm") or "openmm")
        updated["MMGBSA_OPENMM_START_STAGE"] = str(job.get("openmm_start_stage", "heat") or "heat")
        updated["MMGBSA_MD_PROD_PS"] = float(job.get("prod_ps", 4000.0) or 4000.0)
        updated["MMGBSA_MD_FRAME_STRIDE_PS"] = float(job.get("frame_stride_ps", 10.0) or 10.0)
        if str(job.get("validation_tier", "") or "").strip().lower() == "production":
            updated["MMGBSA_PROTOCOL"] = "production"
        if updated["MMGBSA_MD_ENGINE"] == "openmm":
            gpu_token = str(job.get("gpu", "")).strip()
            updated["MMGBSA_OPENMM_DEVICE_INDEX"] = _openmm_logical_device_index_for_mask(gpu_token)
            updated["MMGBSA_OPENMM_CPU_THREADS"] = "1"
    overrides = job.get("cfg_overrides", {})
    if isinstance(overrides, Mapping):
        updated.update(dict(overrides))
    _apply_thread_env(cpu_cap)
    _apply_gpu_env(job, mmpbsa=mmpbsa)
    return updated


def _required_ligand_inputs(ligand: Mapping[str, Any]) -> tuple[bool, list[str]]:
    required = (
        "receptor_pdb",
        "ligand_mol2",
        "ligand_frcmod",
        "dry_complex_prmtop",
        "dry_receptor_prmtop",
        "dry_ligand_prmtop",
    )
    missing = [key for key in required if not Path(str(ligand.get(key, ""))).exists()]
    return not missing, missing


def _run_contract_job(job: Mapping[str, Any]) -> dict[str, Any]:
    ok, missing = _required_ligand_inputs(job["ligand"])
    return {
        "ok": ok,
        "skipped": True,
        "validation_tier": "contract",
        "missing": missing,
        "message": "validated staged job schema and required input paths",
    }


def _run_topology_job(
    *,
    job: Mapping[str, Any],
    cfg: Mapping[str, Any],
    force: bool,
) -> dict[str, Any]:
    ligand = job["ligand"]
    plan = write_solvated_tleap(
        receptor_pdb=str(ligand["receptor_pdb"]),
        ligand_mol2=str(ligand["ligand_mol2"]),
        ligand_frcmod=str(ligand["ligand_frcmod"]),
        dry_complex_prmtop=str(ligand["dry_complex_prmtop"]),
        dry_receptor_prmtop=str(ligand["dry_receptor_prmtop"]),
        dry_ligand_prmtop=str(ligand["dry_ligand_prmtop"]),
        out_dir=str(job["explicit_dir"]),
        force=force,
    )
    tleap = run_tleap(plan, cfg, force=force)
    return {
        "ok": bool(tleap.get("ok")),
        "validation_tier": "topology",
        "plan": plan,
        "tleap": tleap,
    }


def _copy_fixture_outputs(job: Mapping[str, Any], result: Mapping[str, Any]) -> dict[str, str]:
    fixture_root = str(job.get("fixture_out_dir", "") or "").strip()
    if not fixture_root or not result.get("ok"):
        return {}
    ligand = job["ligand"]
    label = f"{ligand.get('pdb_id', 'pdb')}_{ligand.get('ligand_stem', 'ligand')}_rep{job.get('replicate', 1)}"
    dest = Path(fixture_root) / label
    explicit_dir = Path(str(job["explicit_dir"]))
    copied: dict[str, str] = {}
    for name in ("complex_solvated.prmtop", "complex_solvated.inpcrd", "prod.nc", "prod_stripped.nc"):
        source = explicit_dir / name
        if source.exists() and source.stat().st_size > 0:
            dest.mkdir(parents=True, exist_ok=True)
            target = dest / name
            shutil.copy2(source, target)
            copied[name] = str(target)
    if copied:
        _write_json(dest / "fixture_manifest.json", {"job": dict(job), "files": copied})
    return copied


def run_staged_md_job(*, job_json: Path, cfg: Mapping[str, Any], force: bool = False) -> dict[str, Any]:
    job = _read_json(job_json)
    run_cfg = _apply_job_caps(cfg, job, mmpbsa=False)
    ligand = job["ligand"]
    result: dict[str, Any]
    try:
        tier = str(job.get("validation_tier", "custom") or "custom")
        if tier == "contract":
            result = _run_contract_job(job)
        elif tier == "topology":
            result = _run_topology_job(job=job, cfg=run_cfg, force=force)
        else:
            result = build_and_run_explicit_workflow(
                receptor_pdb=str(ligand["receptor_pdb"]),
                ligand_mol2=str(ligand["ligand_mol2"]),
                ligand_frcmod=str(ligand["ligand_frcmod"]),
                dry_complex_prmtop=str(ligand["dry_complex_prmtop"]),
                dry_receptor_prmtop=str(ligand["dry_receptor_prmtop"]),
                dry_ligand_prmtop=str(ligand["dry_ligand_prmtop"]),
                out_dir=str(job["explicit_dir"]),
                cfg=run_cfg,
                seed=int(job.get("seed", 12345)),
                force=force,
                run=True,
            )
            fixture_files = _copy_fixture_outputs(job, result)
            if fixture_files:
                result["fixture_files"] = fixture_files
        status = {"ok": bool(result.get("ok")), "job": job, "result": result}
    except Exception as exc:
        LOGGER.exception("staged MD job failed: %s", job_json)
        status = {"ok": False, "job": job, "error": f"{type(exc).__name__}: {exc}"}
    _write_json(Path(str(job["status_json"])), status)
    return status


def _aggregate_one_replicate(
    *,
    replicate: Mapping[str, Any],
    ligand: Mapping[str, Any],
    cfg: Mapping[str, Any],
    force: bool,
) -> dict[str, Any]:
    rep_idx = int(replicate.get("replicate", 1))
    rep_dir = Path(str(replicate["rep_dir"]))
    explicit_dir = Path(str(replicate["explicit_dir"]))
    traj = explicit_dir / "prod_stripped.nc"
    if not traj.exists():
        traj = explicit_dir / "stripped.nc"
    if not traj.exists() or traj.stat().st_size == 0:
        return {"replicate": rep_idx, "ok": False, "notes": "missing_stripped_traj", "trajectory": str(traj)}
    try:
        mmpbsa_result = run_mmgbsa(
            complex_prmtop=str(ligand["dry_complex_prmtop"]),
            receptor_prmtop=str(ligand["dry_receptor_prmtop"]),
            ligand_prmtop=str(ligand["dry_ligand_prmtop"]),
            trajectory_path=str(traj),
            work_dir=str(rep_dir),
            cfg=cfg,
            force=force,
            run=True,
        )
        csv_path = Path(str(mmpbsa_result.get("out_csv", rep_dir / "FINAL_RESULTS_MMPBSA.csv")))
        frame_meta = _mmgbsa_frame_aggregate(csv_path, cfg, LOGGER)
        delta = frame_meta.get("delta_total")
        if delta is None:
            delta = parse_mmpbsa_delta_total(csv_path)
        ok = delta is not None and _mmgbsa_frame_qc_ok(cfg, frame_meta)
        result = {
            "replicate": rep_idx,
            "ok": ok,
            "delta_total": float(delta) if delta is not None else None,
            "results_csv": str(csv_path),
            "trajectory": str(traj),
            "notes": "" if ok else "frame_qc_failed_or_missing_delta",
        }
        result.update(_mmgbsa_summary_frame_fields(frame_meta))
        return result
    except Exception as exc:
        LOGGER.exception("staged MMPBSA replicate failed: %s", rep_dir)
        return {"replicate": rep_idx, "ok": False, "notes": f"{type(exc).__name__}: {exc}", "trajectory": str(traj)}


def _ok_replicates(rep_results: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [
        dict(rep)
        for rep in rep_results
        if rep.get("ok") and rep.get("delta_total") is not None
    ]


def _outlier_text(ok_reps: list[dict[str, Any]]) -> str:
    outlier_indexes = detect_outliers([float(rep["delta_total"]) for rep in ok_reps])
    for idx in outlier_indexes:
        if idx < len(ok_reps):
            ok_reps[idx]["outlier"] = True
    return ";".join(
        str(ok_reps[idx].get("replicate", idx + 1))
        for idx in outlier_indexes
        if idx < len(ok_reps)
    )


def _best_replicate(
    ok_reps: Sequence[Mapping[str, Any]], agg_delta: float | None
) -> dict[str, Any] | None:
    if agg_delta is not None and ok_reps:
        return dict(
            min(
                ok_reps,
                key=lambda rep: abs(float(rep["delta_total"]) - float(agg_delta)),
            )
        )
    return dict(ok_reps[0]) if ok_reps else None


def _replicate_rollup(
    rep_results: Sequence[Mapping[str, Any]], method: str
) -> tuple[list[dict[str, Any]], float | None, dict[str, Any], str, dict[str, Any] | None]:
    ok_reps = _ok_replicates(rep_results)
    ok_deltas = [float(rep["delta_total"]) for rep in ok_reps]
    agg_delta = _mmgbsa_md_aggregate(ok_deltas, method)
    stats = summarize_values(ok_deltas)
    return ok_reps, agg_delta, stats, _outlier_text(ok_reps), _best_replicate(ok_reps, agg_delta)


def _aggregate_summary_payload(
    *,
    job: Mapping[str, Any],
    rep_results: Sequence[Mapping[str, Any]],
    ok_reps: Sequence[Mapping[str, Any]],
    agg_delta: float | None,
    stats: Mapping[str, Any],
    outlier_text: str,
    best_rep: Mapping[str, Any] | None,
    method: str,
) -> dict[str, Any]:
    return {
        "ok": agg_delta is not None,
        "n_reps_total": len(rep_results),
        "n_reps_ok": len(ok_reps),
        "agg_method": method,
        "agg_delta": agg_delta,
        "replicate_stats": dict(stats),
        "outlier_replicates": outlier_text,
        "best_replicate": dict(best_rep) if best_rep else None,
        "replicates": list(rep_results),
        "job": dict(job),
    }


def _summary_csv_row(
    *,
    ligand: Mapping[str, Any],
    rep_results: Sequence[Mapping[str, Any]],
    ok_reps: Sequence[Mapping[str, Any]],
    agg_delta: float | None,
    stats: Mapping[str, Any],
    outlier_text: str,
    best_rep: Mapping[str, Any] | None,
    method: str,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "stage_dir": ligand.get("stage_dir", ""),
        "ligand_stem": ligand.get("ligand_stem", ""),
        "delta_total": f"{float(agg_delta):.6g}" if agg_delta is not None else "",
        "results_csv": str(best_rep.get("results_csv", "")) if best_rep else "",
        "ok": agg_delta is not None,
        "notes": f"staged_md_reps_ok={len(ok_reps)}/{len(rep_results)} rep_agg={method}",
        "rep_mean": stats.get("mean", "") if stats else "",
        "rep_sd": stats.get("sd", "") if stats else "",
        "rep_sem": stats.get("sem", "") if stats else "",
        "rep_ci95": stats.get("ci95", "") if stats else "",
        "rep_outliers": outlier_text,
        "rep_agg_method": method,
        "replicates_ok": len(ok_reps),
        "replicates_total": len(rep_results),
    }
    if best_rep:
        row.update(_mmgbsa_summary_frame_fields(best_rep))
    return row


def _run_validation_aggregate_job(job: Mapping[str, Any], validation_tier: str) -> dict[str, Any]:
    replicates = list(job.get("replicates", []))
    missing_job_paths = [
        str(rep.get("job_path", ""))
        for rep in replicates
        if not Path(str(rep.get("job_path", ""))).is_file()
    ]
    result = {
        "ok": not missing_job_paths,
        "skipped": True,
        "validation_tier": validation_tier,
        "message": "validated staged aggregate schema without MMPBSA",
        "replicates_total": len(replicates),
        "missing_job_paths": missing_job_paths,
        "job": dict(job),
    }
    _write_json(Path(str(job["summary_json"])), result)
    _write_json(Path(str(job["status_json"])), result)
    return result


def run_staged_aggregate_job(*, job_json: Path, cfg: Mapping[str, Any], force: bool = False) -> dict[str, Any]:
    job = _read_json(job_json)
    run_cfg = _apply_job_caps(cfg, job, mmpbsa=True)
    for replicate in job.get("replicates", []):
        overrides = replicate.get("cfg_overrides", {})
        if isinstance(overrides, Mapping):
            run_cfg.update(overrides)
    validation_tier = str(job.get("validation_tier", "custom") or "custom").strip().lower()
    if validation_tier in {"contract", "topology"}:
        return _run_validation_aggregate_job(job, validation_tier)
    ligand = job["ligand"]
    rep_results = [
        _aggregate_one_replicate(
            replicate=replicate,
            ligand=ligand,
            cfg=run_cfg,
            force=force,
        )
        for replicate in job.get("replicates", [])
    ]
    method = str(run_cfg.get("MMGBSA_MD_REPLICATE_AGG", "mean") or "mean").strip().lower()
    ok_reps, agg_delta, stats, outlier_text, best_rep = _replicate_rollup(
        rep_results, method
    )
    summary = _aggregate_summary_payload(
        job=job,
        rep_results=rep_results,
        ok_reps=ok_reps,
        agg_delta=agg_delta,
        stats=stats,
        outlier_text=outlier_text,
        best_rep=best_rep,
        method=method,
    )
    _write_json(Path(str(job["summary_json"])), summary)
    _write_json(Path(str(job["status_json"])), summary)
    row = _summary_csv_row(
        ligand=ligand,
        rep_results=rep_results,
        ok_reps=ok_reps,
        agg_delta=agg_delta,
        stats=stats,
        outlier_text=outlier_text,
        best_rep=best_rep,
        method=method,
    )
    _mmgbsa_append_summary(Path(str(ligand["mmgbsa_root"])) / "mmgbsa_results_summary.csv", row)
    return summary


__all__ = [
    "run_staged_aggregate_job",
    "run_staged_md_job",
    "stage_from_report",
]
