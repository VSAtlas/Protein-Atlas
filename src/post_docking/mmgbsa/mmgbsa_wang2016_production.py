from __future__ import annotations

import csv
import json
import shlex
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

from post_docking.mmgbsa._atomic_io import write_json_atomic, write_text_atomic
from post_docking.mmgbsa.mmgbsa_batch import MMGBSAHit, load_report_hits
from post_docking.mmgbsa.mmgbsa_review import (
    default_review_record,
    review_record_path,
    validate_review_record,
    write_review_record,
)


def write_wang2016_production_package(
    *,
    report: Path,
    run_id: str,
    out_dir: Path,
    artifacts: Sequence[Path] = (),
    artifact_dest: Path | None = None,
    overall_dir: Path | None = None,
    post_docked_dir: Path | None = None,
    docked_dir: Path | None = None,
    configs_dir: Path | None = None,
    cpus: int = 32,
    partition: str = "",
    hours: int = 72,
    python_exe: str | None = None,
    md_prod_ps: float = 10000.0,
    md_replicates: int = 3,
    frame_stride_ps: float = 10.0,
    analysis_start_ps: float = 0.0,
    analysis_end_ps: float = 0.0,
    analysis_interval: int = 1,
    production_strict: bool = True,
    md_engine: str = "sander",
    gpus: str = "",
    openmm_start_stage: str = "prod",
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    package_report = _write_pdbbind_enriched_report(report, out_dir)
    hits = load_report_hits(package_report, run_id=run_id)
    review_manifest = _write_review_manifest(out_dir / "review_records", hits)
    script_path = out_dir / "run_wang2016_mmgbsa_production.sbatch"
    local_script_path = out_dir / "run_wang2016_mmgbsa_production_local.sh"
    plan_out = out_dir / "wang2016_mmgbsa_production_plan.json"
    log_path = out_dir / "wang2016_mmgbsa_production_%j.log"
    local_log_path = out_dir / "wang2016_mmgbsa_production_local.log"
    pid_path = out_dir / "wang2016_mmgbsa_production_local.pid"
    python_path = python_exe or sys.executable
    script_text = _slurm_script(
        report=package_report,
        run_id=run_id,
        artifacts=artifacts,
        artifact_dest=artifact_dest,
        overall_dir=overall_dir,
        post_docked_dir=post_docked_dir,
        docked_dir=docked_dir,
        configs_dir=configs_dir,
        plan_out=plan_out,
        cpus=min(max(1, int(cpus)), 32),
        partition=partition,
        hours=max(1, int(hours)),
        log_path=log_path,
        python_exe=python_path,
        md_prod_ps=md_prod_ps,
        md_replicates=md_replicates,
        frame_stride_ps=frame_stride_ps,
        analysis_start_ps=analysis_start_ps,
        analysis_end_ps=analysis_end_ps,
        analysis_interval=analysis_interval,
        production_strict=production_strict,
        md_engine=md_engine,
        gpus="all",
        openmm_start_stage=openmm_start_stage,
    )
    local_script_text = _local_script(
        report=package_report,
        run_id=run_id,
        artifacts=artifacts,
        artifact_dest=artifact_dest,
        overall_dir=overall_dir,
        post_docked_dir=post_docked_dir,
        docked_dir=docked_dir,
        configs_dir=configs_dir,
        plan_out=plan_out,
        cpus=min(max(1, int(cpus)), 32),
        python_exe=python_path,
        md_prod_ps=md_prod_ps,
        md_replicates=md_replicates,
        frame_stride_ps=frame_stride_ps,
        analysis_start_ps=analysis_start_ps,
        analysis_end_ps=analysis_end_ps,
        analysis_interval=analysis_interval,
        production_strict=production_strict,
        md_engine=md_engine,
        gpus=gpus,
        openmm_start_stage=openmm_start_stage,
    )
    write_text_atomic(script_path, script_text, require_nonempty=True)
    write_text_atomic(local_script_path, local_script_text, require_nonempty=True)
    script_path.chmod(script_path.stat().st_mode | 0o100)
    local_script_path.chmod(local_script_path.stat().st_mode | 0o100)
    result = {
        "schema_version": 1,
        "report": str(report),
        "package_report": str(package_report),
        "run_id": run_id,
        "script": str(script_path),
        "submit_command": f"sbatch {shlex.quote(str(script_path))}",
        "local_script": str(local_script_path),
        "local_submit_command": (
            f"nohup {shlex.quote(str(local_script_path))} > "
            f"{shlex.quote(str(local_log_path))} 2>&1 & echo $! > "
            f"{shlex.quote(str(pid_path))}"
        ),
        "local_log_path": str(local_log_path),
        "local_pid_path": str(pid_path),
        "plan_out": str(plan_out),
        "log_path": str(log_path),
        "cpus": min(max(1, int(cpus)), 32),
        "md_prod_ps": float(md_prod_ps),
        "md_replicates": int(md_replicates),
        "frame_stride_ps": float(frame_stride_ps),
        "analysis_start_ps": float(analysis_start_ps),
        "analysis_end_ps": float(analysis_end_ps),
        "analysis_interval": int(analysis_interval),
        "production_strict": bool(production_strict),
        "md_engine": str(md_engine),
        "gpus": str(gpus),
        "openmm_start_stage": str(openmm_start_stage),
        "review_manifest": review_manifest,
        "notes": [
            "This writes a scheduler-ready production job; it does not mark ligand chemistry as reviewed.",
            "Review records must be approved or replaced with curated/precharged inputs before publication claims.",
        ],
    }
    write_json_atomic(out_dir / "wang2016_mmgbsa_production_package.json", result)
    return result


def _write_pdbbind_enriched_report(report: Path, out_dir: Path) -> Path:
    pdbbind_root = (
        report.parent.parent / "pdbbindpp2020_targets" / "pbpp-2020"
    ).expanduser()
    if not pdbbind_root.is_dir():
        return report
    enriched = out_dir / "wang2016_mmgbsa_production_enriched_report.csv"
    with report.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        for field in ("pdbbind_mol2", "pdbbind_sdf"):
            if field not in fieldnames:
                fieldnames.append(field)
        rows: list[dict[str, str]] = []
        for row in reader:
            pdb_id = str(row.get("pdb_id", "") or "").strip().lower()
            ligand_dir = pdbbind_root / pdb_id
            mol2 = ligand_dir / f"{pdb_id}_ligand.mol2"
            sdf = ligand_dir / f"{pdb_id}_ligand.sdf"
            out_row = dict(row)
            if mol2.is_file() and not str(out_row.get("pdbbind_mol2", "")).strip():
                out_row["pdbbind_mol2"] = str(mol2)
            if sdf.is_file() and not str(out_row.get("pdbbind_sdf", "")).strip():
                out_row["pdbbind_sdf"] = str(sdf)
            rows.append(out_row)
    with enriched.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return enriched


def _write_review_manifest(review_root: Path, hits: Sequence[MMGBSAHit]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for hit in hits:
        prep_metadata = _prep_metadata_for_hit(hit)
        record = _review_record_for_hit(hit, prep_metadata)
        path = review_record_path(review_root, hit.pdb_id, hit.ligand_id)
        write_review_record(path, record)
        validation = validate_review_record(record, production=True)
        rows.append(
            {
                "pdb_id": hit.pdb_id,
                "ligand_id": hit.ligand_id,
                "review_record": str(path),
                "production_ready": validation["production_ready"],
                "problems": validation["problems"],
                "warnings": validation["warnings"],
                "prep_metadata": prep_metadata.get("metadata_path", ""),
            }
        )
    manifest = {
        "schema_version": 1,
        "curation_status": "requires_review",
        "records": rows,
    }
    write_json_atomic(review_root / "wang2016_ligand_review_manifest.json", manifest)
    return manifest


def _review_record_for_hit(hit: MMGBSAHit, prep_metadata: Mapping[str, Any]) -> dict[str, Any]:
    record = default_review_record(hit.pdb_id, hit.ligand_id)
    source = hit.source_mol2_path or hit.source_sdf_path or "report"
    for key in (
        "ligand_protonation",
        "tautomer",
        "stereochemistry",
        "parameter_review",
        "input_chemistry_authority",
    ):
        section = record[key]
        section["reviewed"] = False
        section["status"] = "requires_review"
        section["source"] = source
        section["notes"] = "Staged automatically from Wang redock report; not expert-reviewed."
    record["net_charge"].update(
        {
            "reviewed": False,
            "status": "requires_review",
            "value": prep_metadata.get("net_charge_used"),
            "source": prep_metadata.get("metadata_path", source),
            "charge_method": prep_metadata.get("charge_method", ""),
            "notes": "Confirm formal charge and charge method before publication MMGBSA.",
        }
    )
    record["water_policy"].update(
        {
            "reviewed": False,
            "status": "requires_review",
            "policy": "structural_waters_unreviewed",
            "source": "Atlas receptor prep outputs",
            "notes": "Confirm retained/displaced structural waters against Wang/PDBbind protocol.",
        }
    )
    record["metal_policy"].update(
        {
            "reviewed": False,
            "status": "requires_review",
            "policy": "metals_unreviewed",
            "source": "Atlas receptor prep outputs",
            "notes": "Confirm metal/cofactor retention and force-field treatment.",
        }
    )
    record["apo_holo"].update(
        {
            "reviewed": False,
            "status": "requires_review",
            "value": "HOLO",
            "source": "Wang crystal-control benchmark",
            "notes": "Crystal-control benchmark starts from holo complexes.",
        }
    )
    record["ph"].update(
        {
            "reviewed": False,
            "status": "requires_review",
            "value": None,
            "source": "PDB header/PROPKA context",
            "notes": "Confirm assay/crystallization pH and ligand state context.",
        }
    )
    record["notes"] = (
        "Automatic Wang benchmark review scaffold. Replace with curated ligand "
        "charge/protonation/tautomer/water/metal decisions before publication use."
    )
    record["source_sdf_path"] = hit.source_sdf_path
    record["source_mol2_path"] = hit.source_mol2_path
    return record


def _prep_metadata_for_hit(hit: MMGBSAHit) -> dict[str, Any]:
    source = Path(hit.source_sdf_path or hit.source_mol2_path or "")
    if not source.name:
        return {}
    try:
        stage_dir = source.parent.name
        ph_dir = source.parent.parent
        candidate = ph_dir / "mmgbsa" / "mol2" / stage_dir / f"{source.stem}.mmgbsa_prep.json"
    except Exception:
        return {}
    if not candidate.exists():
        return {}
    try:
        data = json.loads(candidate.read_text(encoding="utf-8"))
    except Exception:
        return {}
    if isinstance(data, dict):
        data["metadata_path"] = str(candidate)
        return data
    return {}


def _wang_runtime_arg_parser(
    *,
    cpus: int,
    md_prod_ps: float,
    md_replicates: int,
    frame_stride_ps: float,
    analysis_start_ps: float,
    analysis_end_ps: float,
    analysis_interval: int,
    production_strict: bool,
    md_engine: str,
    gpus: str,
    openmm_start_stage: str,
    include_cpu_flag: bool,
) -> list[str]:
    cpu_lines = [f'CPU="${{CPU:-{cpus}}}"'] if include_cpu_flag else []
    cpu_case = (
        [
            "    --cpus)",
            "      CPU=\"$2\"",
            "      shift 2",
            "      ;;",
        ]
        if include_cpu_flag
        else []
    )
    strict_default = "true" if production_strict else "false"
    return [
        *cpu_lines,
        f'WANG_MD_PROD_PS="${{MMGBSA_MD_PROD_PS:-{float(md_prod_ps)}}}"',
        f'WANG_MD_REPLICATES="${{MMGBSA_MD_NREPLICATES:-{int(md_replicates)}}}"',
        f'WANG_FRAME_STRIDE_PS="${{MMGBSA_MD_FRAME_STRIDE_PS:-{float(frame_stride_ps)}}}"',
        f'WANG_ANALYSIS_START_PS="${{MMGBSA_ANALYSIS_START_PS:-{float(analysis_start_ps)}}}"',
        f'WANG_ANALYSIS_END_PS="${{MMGBSA_ANALYSIS_END_PS:-{float(analysis_end_ps)}}}"',
        f'WANG_ANALYSIS_INTERVAL="${{MMGBSA_ANALYSIS_INTERVAL:-{int(analysis_interval)}}}"',
        f'WANG_PRODUCTION_STRICT="${{MMGBSA_PRODUCTION_STRICT:-{strict_default}}}"',
        f'WANG_MD_ENGINE="${{MMGBSA_MD_ENGINE:-{str(md_engine)}}}"',
        f'WANG_GPUS="${{MMGBSA_OPENMM_DEVICE_INDEX:-{str(gpus)}}}"',
        f'WANG_OPENMM_START_STAGE="${{MMGBSA_OPENMM_START_STAGE:-{str(openmm_start_stage)}}}"',
        'WANG_MD_MPI_RANKS="${MMGBSA_MD_MPI_RANKS:-auto}"',
        'WANG_MMPBSA_MPI_RANKS="${MMGBSA_MMPBSA_MPI_RANKS:-auto}"',
        "usage() {",
        '  echo "Usage: $0 [--prod-ps PS] [--replicates N] [--frame-stride-ps PS] [--analysis-start-ps PS] [--analysis-end-ps PS] [--analysis-interval N] [--md-engine sander|openmm] [--gpus IDX|all] [--openmm-start-stage STAGE] [--md-ranks N|auto] [--mmpbsa-ranks N|auto] [--validation-mode|--strict-production] [--cpus N]" >&2',
        "}",
        "while [ \"$#\" -gt 0 ]; do",
        "  case \"$1\" in",
        "    --prod-ps|--md-prod-ps)",
        "      WANG_MD_PROD_PS=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --replicates|--md-replicates)",
        "      WANG_MD_REPLICATES=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --frame-stride-ps)",
        "      WANG_FRAME_STRIDE_PS=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --analysis-start-ps)",
        "      WANG_ANALYSIS_START_PS=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --analysis-end-ps)",
        "      WANG_ANALYSIS_END_PS=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --analysis-interval)",
        "      WANG_ANALYSIS_INTERVAL=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --md-ranks)",
        "      WANG_MD_MPI_RANKS=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --md-engine)",
        "      WANG_MD_ENGINE=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --gpus|--gpu-indices)",
        "      WANG_GPUS=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --openmm-start-stage)",
        "      WANG_OPENMM_START_STAGE=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --mmpbsa-ranks)",
        "      WANG_MMPBSA_MPI_RANKS=\"$2\"",
        "      shift 2",
        "      ;;",
        "    --validation-mode)",
        "      WANG_PRODUCTION_STRICT=false",
        "      shift",
        "      ;;",
        "    --strict-production)",
        "      WANG_PRODUCTION_STRICT=true",
        "      shift",
        "      ;;",
        *cpu_case,
        "    -h|--help)",
        "      usage",
        "      exit 0",
        "      ;;",
        "    *)",
        '      echo "Unknown Wang MMGBSA option: $1" >&2',
        "      usage",
        "      exit 2",
        "      ;;",
        "  esac",
        "done",
    ]


def _slurm_script(
    *,
    report: Path,
    run_id: str,
    artifacts: Sequence[Path],
    artifact_dest: Path | None,
    overall_dir: Path | None,
    post_docked_dir: Path | None,
    docked_dir: Path | None,
    configs_dir: Path | None,
    plan_out: Path,
    cpus: int,
    partition: str,
    hours: int,
    log_path: Path,
    python_exe: str,
    md_prod_ps: float,
    md_replicates: int,
    frame_stride_ps: float,
    analysis_start_ps: float,
    analysis_end_ps: float,
    analysis_interval: int,
    production_strict: bool,
    md_engine: str,
    gpus: str,
    openmm_start_stage: str,
) -> str:
    amber_mpi_prefix = Path(__file__).resolve().parents[5] / "tools" / "envs" / "ambertools-mpi"
    artifact_args = " ".join(f"--artifact {shlex.quote(str(path))}" for path in artifacts)
    optional_args = _optional_args(
        {
            "--artifact-dest": artifact_dest,
            "--overall-dir": overall_dir,
            "--post-docked-dir": post_docked_dir,
            "--docked-dir": docked_dir,
            "--configs-dir": configs_dir,
        }
    )
    partition_line = f"#SBATCH --partition={partition}\n" if partition else ""
    command = " ".join(
        part
        for part in (
            f"{shlex.quote(python_exe)} -m post_docking.mmgbsa.cli batch-from-report",
            f"--report {shlex.quote(str(report))}",
            f"--run-id {shlex.quote(run_id)}",
            artifact_args,
            optional_args,
            f"--plan-out {shlex.quote(str(plan_out))}",
            "--log-level INFO",
        )
        if part
    )
    return "\n".join(
        [
            "#!/usr/bin/env bash",
            f"#SBATCH --job-name=wang2016-mmgbsa-{run_id}",
            f"#SBATCH --cpus-per-task={cpus}",
            f"#SBATCH --time={hours}:00:00",
            f"#SBATCH --output={log_path}",
            partition_line.rstrip(),
            "set -euo pipefail",
            *_wang_runtime_arg_parser(
                cpus=cpus,
                md_prod_ps=md_prod_ps,
                md_replicates=md_replicates,
                frame_stride_ps=frame_stride_ps,
                analysis_start_ps=analysis_start_ps,
                analysis_end_ps=analysis_end_ps,
                analysis_interval=analysis_interval,
                production_strict=production_strict,
                md_engine=md_engine,
                gpus=gpus or "all",
                openmm_start_stage=openmm_start_stage,
                include_cpu_flag=False,
            ),
            "export CPU=${SLURM_CPUS_PER_TASK:-" + str(cpus) + "}",
            "export GLOBAL_SCHEDULER_CPUS=${CPU}",
            f"export MMGBSA_AMBERTOOLS_PREFIX=${{MMGBSA_AMBERTOOLS_PREFIX:-{shlex.quote(str(amber_mpi_prefix))}}}",
            "export AMBERTOOLS_PREFIX=${AMBERTOOLS_PREFIX:-$MMGBSA_AMBERTOOLS_PREFIX}",
            "export MMGBSA_MD_MPI_ENABLED=${MMGBSA_MD_MPI_ENABLED:-true}",
            "export MMGBSA_MD_ENGINE=${WANG_MD_ENGINE}",
            "export MMGBSA_OPENMM_PLATFORM=${MMGBSA_OPENMM_PLATFORM:-auto}",
            "export MMGBSA_OPENMM_DEVICE_INDEX=${MMGBSA_OPENMM_DEVICE_INDEX:-0}",
            "export MMGBSA_OPENMM_START_STAGE=${WANG_OPENMM_START_STAGE}",
            "if [ \"${WANG_MD_ENGINE}\" = \"openmm\" ] && [ \"${WANG_GPUS}\" != \"\" ] && [ \"${WANG_GPUS}\" != \"all\" ]; then",
            "  export ROCR_VISIBLE_DEVICES=${WANG_GPUS}",
            "  unset HIP_VISIBLE_DEVICES GPU_DEVICE_ORDINAL CUDA_VISIBLE_DEVICES",
            "fi",
            "export MMGBSA_MD_MPI_RANKS=${WANG_MD_MPI_RANKS}",
            "export MMGBSA_MMPBSA_MPI_ENABLED=${MMGBSA_MMPBSA_MPI_ENABLED:-true}",
            "export MMGBSA_MMPBSA_MPI_RANKS=${WANG_MMPBSA_MPI_RANKS}",
            "export MMGBSA_PROTOCOL=production",
            "export MMGBSA_PRODUCTION_STRICT=${WANG_PRODUCTION_STRICT}",
            "export MMGBSA_MD_NREPLICATES=${WANG_MD_REPLICATES}",
            "export MMGBSA_MD_PROD_PS=${WANG_MD_PROD_PS}",
            "export MMGBSA_MD_FRAME_STRIDE_PS=${WANG_FRAME_STRIDE_PS}",
            "export MMGBSA_ANALYSIS_START_PS=${WANG_ANALYSIS_START_PS}",
            "export MMGBSA_ANALYSIS_END_PS=${WANG_ANALYSIS_END_PS}",
            "export MMGBSA_ANALYSIS_INTERVAL=${WANG_ANALYSIS_INTERVAL}",
            "export MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2=${MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2:-true}",
            "export MMGBSA_LIGAND_ANTECHAMBER_TIMEOUT_SEC=${MMGBSA_LIGAND_ANTECHAMBER_TIMEOUT_SEC:-1800}",
            ": ${MMGBSA_CHEMISTRY_REVIEW_STATUS:?Set after Wang ligand/receptor review}",
            ": ${MMGBSA_PROTONATION_REVIEW_STATUS:?Set after ligand protonation review}",
            ": ${MMGBSA_TAUTOMER_REVIEW_STATUS:?Set after ligand tautomer review}",
            ": ${MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS:?Set after ligand stereochemistry review}",
            ": ${MMGBSA_NET_CHARGE_REVIEW_STATUS:?Set after ligand net-charge review}",
            ": ${MMGBSA_PARAMETER_REVIEW_STATUS:?Set after ligand parameter review}",
            ": ${MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE:?Set true only after confirming authoritative SDF/MOL2 chemistry}",
            command,
            "",
        ]
    )


def _local_script(
    *,
    report: Path,
    run_id: str,
    artifacts: Sequence[Path],
    artifact_dest: Path | None,
    overall_dir: Path | None,
    post_docked_dir: Path | None,
    docked_dir: Path | None,
    configs_dir: Path | None,
    plan_out: Path,
    cpus: int,
    python_exe: str,
    md_prod_ps: float,
    md_replicates: int,
    frame_stride_ps: float,
    analysis_start_ps: float,
    analysis_end_ps: float,
    analysis_interval: int,
    production_strict: bool,
    md_engine: str,
    gpus: str,
    openmm_start_stage: str,
) -> str:
    amber_mpi_prefix = Path(__file__).resolve().parents[5] / "tools" / "envs" / "ambertools-mpi"
    artifact_args = " ".join(f"--artifact {shlex.quote(str(path))}" for path in artifacts)
    optional_args = _optional_args(
        {
            "--artifact-dest": artifact_dest,
            "--overall-dir": overall_dir,
            "--post-docked-dir": post_docked_dir,
            "--docked-dir": docked_dir,
            "--configs-dir": configs_dir,
        }
    )
    command = " ".join(
        part
        for part in (
            f"{shlex.quote(python_exe)} -m post_docking.mmgbsa.cli batch-from-report",
            f"--report {shlex.quote(str(report))}",
            f"--run-id {shlex.quote(run_id)}",
            artifact_args,
            optional_args,
            f"--plan-out {shlex.quote(str(plan_out))}",
            "--log-level INFO",
        )
        if part
    )
    return "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            *_wang_runtime_arg_parser(
                cpus=cpus,
                md_prod_ps=md_prod_ps,
                md_replicates=md_replicates,
                frame_stride_ps=frame_stride_ps,
                analysis_start_ps=analysis_start_ps,
                analysis_end_ps=analysis_end_ps,
                analysis_interval=analysis_interval,
                production_strict=production_strict,
                md_engine=md_engine,
                gpus=gpus,
                openmm_start_stage=openmm_start_stage,
                include_cpu_flag=True,
            ),
            f"CPU=${{CPU:-{cpus}}}",
            "if [ \"$CPU\" -gt 32 ]; then",
            "  echo \"Refusing to run Wang MMGBSA with CPU=$CPU; maximum is 32.\" >&2",
            "  exit 2",
            "fi",
            "export CPU",
            "CPU_LAST=$((CPU - 1))",
            "CPUSET=0-$CPU_LAST",
            "export GLOBAL_SCHEDULER_CPUS=${GLOBAL_SCHEDULER_CPUS:-$CPU}",
            f"export MMGBSA_AMBERTOOLS_PREFIX=${{MMGBSA_AMBERTOOLS_PREFIX:-{shlex.quote(str(amber_mpi_prefix))}}}",
            "export AMBERTOOLS_PREFIX=${AMBERTOOLS_PREFIX:-$MMGBSA_AMBERTOOLS_PREFIX}",
            "export MMGBSA_MD_MPI_ENABLED=${MMGBSA_MD_MPI_ENABLED:-true}",
            "export MMGBSA_MD_ENGINE=${WANG_MD_ENGINE}",
            "export MMGBSA_OPENMM_PLATFORM=${MMGBSA_OPENMM_PLATFORM:-auto}",
            "export MMGBSA_OPENMM_DEVICE_INDEX=${MMGBSA_OPENMM_DEVICE_INDEX:-0}",
            "export MMGBSA_OPENMM_START_STAGE=${WANG_OPENMM_START_STAGE}",
            "if [ \"${WANG_MD_ENGINE}\" = \"openmm\" ] && [ \"${WANG_GPUS}\" != \"\" ] && [ \"${WANG_GPUS}\" != \"all\" ]; then",
            "  export ROCR_VISIBLE_DEVICES=${WANG_GPUS}",
            "  unset HIP_VISIBLE_DEVICES GPU_DEVICE_ORDINAL CUDA_VISIBLE_DEVICES",
            "fi",
            "export MMGBSA_MD_MPI_RANKS=${WANG_MD_MPI_RANKS}",
            "export MMGBSA_MMPBSA_MPI_ENABLED=${MMGBSA_MMPBSA_MPI_ENABLED:-true}",
            "export MMGBSA_MMPBSA_MPI_RANKS=${WANG_MMPBSA_MPI_RANKS}",
            "export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}",
            "export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}",
            "export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}",
            "export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}",
            "export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}",
            "export MMGBSA_PROTOCOL=production",
            "export MMGBSA_FORCE=${MMGBSA_FORCE:-true}",
            "export MMGBSA_PRODUCTION_STRICT=${WANG_PRODUCTION_STRICT}",
            "export MMGBSA_MD_NREPLICATES=${WANG_MD_REPLICATES}",
            "export MMGBSA_MD_PROD_PS=${WANG_MD_PROD_PS}",
            "export MMGBSA_MD_FRAME_STRIDE_PS=${WANG_FRAME_STRIDE_PS}",
            "export MMGBSA_ANALYSIS_START_PS=${WANG_ANALYSIS_START_PS}",
            "export MMGBSA_ANALYSIS_END_PS=${WANG_ANALYSIS_END_PS}",
            "export MMGBSA_ANALYSIS_INTERVAL=${WANG_ANALYSIS_INTERVAL}",
            "export MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2=${MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2:-true}",
            "export MMGBSA_LIGAND_ANTECHAMBER_TIMEOUT_SEC=${MMGBSA_LIGAND_ANTECHAMBER_TIMEOUT_SEC:-1800}",
            ": ${MMGBSA_CHEMISTRY_REVIEW_STATUS:?Set after Wang ligand/receptor review}",
            ": ${MMGBSA_PROTONATION_REVIEW_STATUS:?Set after ligand protonation review}",
            ": ${MMGBSA_TAUTOMER_REVIEW_STATUS:?Set after ligand tautomer review}",
            ": ${MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS:?Set after ligand stereochemistry review}",
            ": ${MMGBSA_NET_CHARGE_REVIEW_STATUS:?Set after ligand net-charge review}",
            ": ${MMGBSA_PARAMETER_REVIEW_STATUS:?Set after ligand parameter review}",
            ": ${MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE:?Set true only after confirming authoritative SDF/MOL2 chemistry}",
            "if command -v taskset >/dev/null 2>&1; then",
            f"  exec taskset -c \"$CPUSET\" {command}",
            "fi",
            f"exec {command}",
            "",
        ]
    )


def _optional_args(options: Mapping[str, Path | None]) -> str:
    parts: list[str] = []
    for flag, value in options.items():
        if value is not None:
            parts.append(f"{flag} {shlex.quote(str(value))}")
    return " ".join(parts)


__all__ = ["write_wang2016_production_package"]
