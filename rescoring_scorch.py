from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import sys
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from input_and_export_functions import load_config

COMPONENT = "[scorch-rescore]"
SCORCH_SCRIPT: Path | None = None
SCORCH_ENV: str = "scorch-env"
SCORCH_ROOT: Path | None = None
VINA_STAGE_DIRS = ("stage1", "stage2", "stage3")
GNINA_STAGE_DIRS = ("gnina_stage1", "gnina_stage2", "gnina_stage3")
POST_STAGE_DIRS = ("ledock_pdbqt", "dock6_pdbqt")
OUTPUT_NAMES = {
    "stage1": "scorch_scores_stage1.csv",
    "stage2": "scorch_scores_stage2.csv",
    "stage3": "scorch_scores_stage3.csv",
    "gnina_stage1": "scorch_scores_gnina_stage1.csv",
    "gnina_stage2": "scorch_scores_gnina_stage2.csv",
    "gnina_stage3": "scorch_scores_gnina_stage3.csv",
    "ledock_pdbqt": "scorch_scores_ledock.csv",
    "dock6_pdbqt": "scorch_scores_dock6.csv",
}


@dataclass(frozen=True)
class StageSpec:
    source: str  # vina, gnina, ledock, dock6
    stage_dir: str
    output_name: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run SCORCH rescoring for a run-id")
    parser.add_argument("--run-id", required=True, help="Run identifier under docked/")
    parser.add_argument(
        "--repo-root",
        default=None,
        help="Repository root (default: directory containing rescoring_scorch.py)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Threads passed to SCORCH (default: 16)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel SCORCH processes (default: 1; keep small to avoid oversubscription)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-run scoring even if outputs already exist",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging for this script",
    )
    parser.add_argument(
        "--skip-autofix",
        action="store_true",
        help="Skip automatic pose_bust / prep_for_scorch repair steps",
    )
    return parser.parse_args()


def configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("rescoring_scorch")


def _resolve_roots(args: argparse.Namespace) -> Tuple[Path, Path, Path, Path]:
    repo_root = Path(args.repo_root).resolve() if args.repo_root else Path(__file__).resolve().parent
    docked_root = repo_root / "docked"
    post_docked_root = repo_root / "post_docked"
    processed_root = repo_root / "processed_pdbs"
    return repo_root, docked_root, post_docked_root, processed_root


def _preflight(logger: logging.Logger) -> bool:
    global SCORCH_SCRIPT, SCORCH_ENV, SCORCH_ROOT

    if shutil.which("micromamba") is None:
        logger.error("%s action=preflight status=failed reason=missing_micromamba", COMPONENT)
        return False

    try:
        cfg = load_config()
    except Exception as exc:
        logger.error("%s action=preflight status=failed reason=config_load_error error=%s", COMPONENT, exc)
        return False

    script_cfg = cfg.get("SCORCH_SCRIPT")
    env_cfg = cfg.get("SCORCH_ENV")

    if script_cfg:
        SCORCH_SCRIPT = Path(script_cfg)

    if env_cfg:
        SCORCH_ENV = str(env_cfg)

    if not SCORCH_SCRIPT:
        logger.error("%s action=preflight status=failed reason=missing_scorch_script_cfg", COMPONENT)
        return False

    SCORCH_SCRIPT = SCORCH_SCRIPT.resolve()
    if not SCORCH_SCRIPT.exists():
        logger.error("%s action=preflight status=failed reason=missing_scorch path=%s", COMPONENT, SCORCH_SCRIPT)
        return False

    SCORCH_ROOT = SCORCH_SCRIPT.parent
    logger.info(
        "%s action=preflight status=ok scorch_script=%s scorch_env=%s scorch_root=%s",
        COMPONENT,
        SCORCH_SCRIPT,
        SCORCH_ENV,
        SCORCH_ROOT,
    )
    return True


def _collect_combo_from_rel(parts: Sequence[str]) -> Optional[Tuple[str, str, str]]:
    if len(parts) < 3:
        return None
    return parts[0], parts[1], parts[2]


def discover_combos(run_root: Path, post_root: Path) -> Set[Tuple[str, str, str]]:
    combos: Set[Tuple[str, str, str]] = set()
    for stage in VINA_STAGE_DIRS + GNINA_STAGE_DIRS:
        for path in run_root.rglob(stage):
            try:
                rel = path.relative_to(run_root).parts
            except ValueError:
                continue
            combo = _collect_combo_from_rel(rel)
            if combo:
                combos.add(combo)

    for stage_dir in POST_STAGE_DIRS:
        for path in post_root.rglob(stage_dir):
            try:
                rel = path.relative_to(post_root).parts
            except ValueError:
                continue
            combo = _collect_combo_from_rel(rel)
            if combo:
                combos.add(combo)
    return combos


def _posebusters_missing(post_root: Path, combos: Set[Tuple[str, str, str]]) -> bool:
    if not post_root.exists():
        return True
    if combos:
        for pdb_id, variant, ph in combos:
            csv_path = post_root / pdb_id / variant / ph / "posebusters_all_stages.csv"
            if not csv_path.exists():
                return True
        return False
    # If no combos, fall back to presence of any posebusters output
    any_pb = any(post_root.rglob("posebusters_all_stages.csv"))
    return not any_pb


def _prep_missing(post_root: Path, combos: Set[Tuple[str, str, str]]) -> bool:
    def _dir_empty(path: Path) -> bool:
        return not path.exists() or not any(path.glob("*.pdbqt"))

    if combos:
        for pdb_id, variant, ph in combos:
            if _dir_empty(post_root / pdb_id / variant / ph / "ledock_pdbqt"):
                return True
            if _dir_empty(post_root / pdb_id / variant / ph / "dock6_pdbqt"):
                return True
        return False
    # Fallback: look for any prepared ligands at all
    any_pdbqt = any(post_root.rglob("*_pdbqt/*.pdbqt"))
    return not any_pdbqt


def _run_pose_bust(run_id: str, repo_root: Path, overwrite: bool, logger: logging.Logger) -> bool:
    cmd = [sys.executable, str(repo_root / "pose_bust.py"), "--run-id", run_id]
    if overwrite:
        cmd.append("--overwrite")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.error(
            "%s action=autofix status=failed step=pose_bust returncode=%s stderr=%s stdout=%s",
            COMPONENT,
            proc.returncode,
            proc.stderr.strip(),
            proc.stdout.strip(),
        )
        return False
    logger.info("%s action=autofix status=ok step=pose_bust run_id=%s", COMPONENT, run_id)
    return True


def _run_prep_for_scorch(run_id: str, repo_root: Path, overwrite: bool, logger: logging.Logger) -> bool:
    cmd = [sys.executable, str(repo_root / "prep_for_scorch.py"), "--run-id", run_id]
    if overwrite:
        cmd.append("--overwrite")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.error(
            "%s action=autofix status=failed step=prep_for_scorch returncode=%s stderr=%s stdout=%s",
            COMPONENT,
            proc.returncode,
            proc.stderr.strip(),
            proc.stdout.strip(),
        )
        return False
    logger.info("%s action=autofix status=ok step=prep_for_scorch run_id=%s", COMPONENT, run_id)
    return True


def _scorch_command(receptor: Path, ligands: Path, threads: int) -> List[str]:
    return [
        "micromamba",
        "run",
        "-n",
        SCORCH_ENV,
        "python",
        str(SCORCH_SCRIPT),
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligands),
        "--out",
        "{out}",
        "--threads",
        str(threads),
        "--verbose",
    ]


def _annotate_csv(csv_path: Path, metadata: Dict[str, str], logger: logging.Logger) -> bool:
    try:
        with csv_path.open() as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            rows = list(reader)
    except Exception as exc:
        logger.error(
            "%s action=annotate status=failed path=%s reason=read_error error=%s",
            COMPONENT,
            csv_path,
            exc,
        )
        return False

    if not rows:
        logger.warning("%s action=annotate status=warning reason=no_rows path=%s", COMPONENT, csv_path)

    ligand_keys = ("ligand", "Ligand", "ligand_file", "Ligand_file", "file", "filename", "name", "molecule", "Molecule")
    out_rows: List[Dict[str, str]] = []
    for row in rows:
        ligand_val = ""
        for key in ligand_keys:
            if row.get(key):
                ligand_val = row[key]
                break
        ligand_val = Path(ligand_val).name if ligand_val else ""
        row.update(metadata)
        row["ligand_file"] = ligand_val or metadata.get("ligand_file", "")
        out_rows.append(row)

    meta_fields = list(metadata.keys()) + ["ligand_file"]
    out_fields = meta_fields + [fn for fn in fieldnames if fn not in meta_fields]

    try:
        with csv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=out_fields)
            writer.writeheader()
            for row in out_rows:
                writer.writerow(row)
    except Exception as exc:
        logger.error(
            "%s action=annotate status=failed path=%s reason=write_error error=%s",
            COMPONENT,
            csv_path,
            exc,
        )
        return False
    return True


def _score_stage(
    spec: StageSpec,
    combo: Tuple[str, str, str],
    run_root: Path,
    post_root: Path,
    receptor: Path,
    threads: int,
    overwrite: bool,
    logger: logging.Logger,
) -> Tuple[bool, Optional[Path]]:
    pdb_id, variant, ph = combo
    if spec.source in {"vina", "gnina"}:
        lig_path = run_root / pdb_id / variant / ph / spec.stage_dir
    else:
        lig_path = post_root / pdb_id / variant / ph / spec.stage_dir

    out_path = post_root / pdb_id / variant / ph / spec.output_name
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not overwrite:
        logger.info(
            "%s action=score status=skip source=%s stage=%s reason=exists output=%s",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            out_path,
        )
        return True, out_path

    if not lig_path.exists():
        logger.warning(
            "%s action=score status=skip source=%s stage=%s reason=missing_ligands path=%s",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            lig_path,
        )
        return True, None

    ligands_available = list(lig_path.glob("*.pdbqt")) if lig_path.is_dir() else ([lig_path] if lig_path.suffix == ".pdbqt" else [])
    if not ligands_available:
        logger.warning(
            "%s action=score status=skip source=%s stage=%s reason=no_pdbqt path=%s",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            lig_path,
        )
        return True, None

    cmd = _scorch_command(receptor, lig_path, threads)
    cmd[cmd.index("{out}")] = str(out_path)
    run_kwargs = {"capture_output": True, "text": True}
    if SCORCH_ROOT is not None:
        run_kwargs["cwd"] = str(SCORCH_ROOT)

    proc = subprocess.run(cmd, **run_kwargs)
    if proc.returncode != 0:
        logger.error(
            "%s action=score status=failed source=%s stage=%s returncode=%s stderr=%s",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            proc.returncode,
            proc.stderr.strip(),
        )
        return False, None

    if not out_path.exists() or out_path.stat().st_size == 0:
        logger.error(
            "%s action=score status=failed source=%s stage=%s reason=empty_output output=%s",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            out_path,
        )
        return False, None

    metadata = {
        "pdb_id": pdb_id,
        "variant": variant,
        "ph": ph,
        "source": spec.source,
        "stage_dir": spec.stage_dir,
    }
    if not _annotate_csv(out_path, metadata, logger):
        return False, None

    logger.info(
        "%s action=score status=ok source=%s stage=%s ligands=%d output=%s",
        COMPONENT,
        spec.source,
        spec.stage_dir,
        len(ligands_available),
        out_path,
    )
    return True, out_path


def _aggregate_all(post_root: Path, specs: List[StageSpec], combos: Iterable[Tuple[str, str, str]], logger: logging.Logger) -> None:
    rows: List[Dict[str, str]] = []
    fields: List[str] = []
    for pdb_id, variant, ph in combos:
        for spec in specs:
            csv_path = post_root / pdb_id / variant / ph / spec.output_name
            if not csv_path.exists():
                continue
            try:
                with csv_path.open() as handle:
                    reader = csv.DictReader(handle)
                    if not reader.fieldnames:
                        continue
                    for row in reader:
                        rows.append(row)
                    for fn in reader.fieldnames:
                        if fn not in fields:
                            fields.append(fn)
            except Exception as exc:
                logger.error(
                    "%s action=aggregate status=failed path=%s reason=read_error error=%s",
                    COMPONENT,
                    csv_path,
                    exc,
                )
                continue

    if not rows:
        logger.warning("%s action=aggregate status=skip reason=no_rows", COMPONENT)
        return

    all_out = post_root / "scorch_scores_all.csv"
    try:
        with all_out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
    except Exception as exc:
        logger.error(
            "%s action=aggregate status=failed path=%s reason=write_error error=%s",
            COMPONENT,
            all_out,
            exc,
        )
        return
    logger.info("%s action=aggregate status=ok rows=%d output=%s", COMPONENT, len(rows), all_out)


def main() -> int:
    args = parse_args()
    logger = configure_logging(args.verbose)
    if not _preflight(logger):
        return 1

    repo_root, docked_root, post_root, processed_root = _resolve_roots(args)
    run_root = docked_root / args.run_id
    post_run_root = post_root / args.run_id
    if not run_root.exists():
        logger.error("%s action=preflight status=failed reason=missing_run_root path=%s", COMPONENT, run_root)
        return 1

    specs: List[StageSpec] = [
        StageSpec("vina", "stage1", OUTPUT_NAMES["stage1"]),
        StageSpec("vina", "stage2", OUTPUT_NAMES["stage2"]),
        StageSpec("vina", "stage3", OUTPUT_NAMES["stage3"]),
        StageSpec("gnina", "gnina_stage1", OUTPUT_NAMES["gnina_stage1"]),
        StageSpec("gnina", "gnina_stage2", OUTPUT_NAMES["gnina_stage2"]),
        StageSpec("gnina", "gnina_stage3", OUTPUT_NAMES["gnina_stage3"]),
        StageSpec("ledock", "ledock_pdbqt", OUTPUT_NAMES["ledock_pdbqt"]),
        StageSpec("dock6", "dock6_pdbqt", OUTPUT_NAMES["dock6_pdbqt"]),
    ]

    combos = discover_combos(run_root, post_run_root)
    logger.info("%s action=discover status=ok combos=%d", COMPONENT, len(combos))

    if not args.skip_autofix:
        if _posebusters_missing(post_run_root, combos):
            _run_pose_bust(args.run_id, repo_root, args.overwrite, logger)
        if _prep_missing(post_run_root, combos):
            _run_prep_for_scorch(args.run_id, repo_root, args.overwrite, logger)

    combos = discover_combos(run_root, post_run_root)
    if not combos:
        logger.warning("%s action=discover status=skip reason=no_combos run_id=%s", COMPONENT, args.run_id)
        return 0

    tasks: List[Tuple[StageSpec, Tuple[str, str, str], Path]] = []
    skipped_missing_receptor = 0
    for combo in sorted(combos):
        pdb_id, variant, ph = combo
        receptor = processed_root / pdb_id / variant / "receptor" / "ph_ensemble" / f"{pdb_id}_{ph}.pdbqt"
        if not receptor.exists():
            logger.error(
                "%s action=score status=skip reason=missing_receptor pdb_id=%s variant=%s ph=%s receptor=%s",
                COMPONENT,
                pdb_id,
                variant,
                ph,
                receptor,
            )
            skipped_missing_receptor += 1
            continue
        for spec in specs:
            tasks.append((spec, combo, receptor))

    total_jobs = len(tasks)
    completed = 0
    failed_jobs = 0

    if args.jobs <= 1:
        for spec, combo, receptor in tasks:
            ok, _ = _score_stage(
                spec, combo, run_root, post_run_root, receptor, args.threads, args.overwrite, logger
            )
            if ok:
                completed += 1
            else:
                failed_jobs += 1
    else:
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            future_map = {
                pool.submit(
                    _score_stage, spec, combo, run_root, post_run_root, receptor, args.threads, args.overwrite, logger
                ): (spec, combo)
                for spec, combo, receptor in tasks
            }
            for future in as_completed(future_map):
                try:
                    ok, _ = future.result()
                    if ok:
                        completed += 1
                    else:
                        failed_jobs += 1
                except Exception as exc:  # defensive
                    failed_jobs += 1
                    spec, combo = future_map[future]
                    logger.error(
                        "%s action=score status=failed reason=worker_exception source=%s stage=%s combo=%s error=%s",
                        COMPONENT,
                        spec.source,
                        spec.stage_dir,
                        combo,
                        exc,
                    )

    _aggregate_all(post_run_root, specs, combos, logger)

    logger.info(
        "%s action=summary status=%s combos=%d total_jobs=%d completed=%d failed=%d missing_receptor=%d",
        COMPONENT,
        "ok" if failed_jobs == 0 else "failed",
        len(combos),
        total_jobs,
        completed,
        failed_jobs,
        skipped_missing_receptor,
    )
    return 0 if failed_jobs == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
