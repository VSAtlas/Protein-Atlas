from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


LEDOCK_STAGE_DIRS = {"ledock_stage1", "ledock_stage2", "ledock_stage3"}
LEDOCK_DUD_STAGE_DIRS = {"ledock_dud_stage1", "ledock_dud_stage2", "ledock_dud_stage3"}
LEDOCK_DUD_STAGE_DIRS_LEGACY = {
    "dud_ledock_stage1",
    "dud_ledock_stage2",
    "dud_ledock_stage3",
}
DOCK6_STAGE_DIRS = {"dock6_stage1", "dock6_stage2", "dock6_stage3"}
DOCK6_DUD_STAGE_DIRS = {"dock6_dud_stage1", "dock6_dud_stage2", "dock6_dud_stage3"}
DOCK6_DUD_STAGE_DIRS_LEGACY = {
    "dud_dock6_stage1",
    "dud_dock6_stage2",
    "dud_dock6_stage3",
}
COMPONENT = "[prep-for-scorch]"
EXAMPLE_LIMIT = 8


@dataclass(frozen=True)
class LedockTask:
    input_path: Path
    output_dir: Path
    stage_dir: str
    pdb: str
    variant: str
    ph: str


@dataclass(frozen=True)
class Dock6Task:
    input_path: Path
    output_dir: Path
    stage_dir: str


@dataclass
class Dock6Result:
    molecules: int = 0
    converted: int = 0
    skipped: int = 0
    failed: int = 0
    mappings: Optional[List[str]] = None

    def __post_init__(self) -> None:
        if self.mappings is None:
            self.mappings = []


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert LeDock DOK and DOCK6 MOL2 outputs to PDBQT under post_docked/"
    )
    parser.add_argument("--run-id", required=True, help="Run identifier under docked/")
    parser.add_argument(
        "--repo-root",
        default=None,
        help="Repository root (default: directory containing prep_for_scorch.py)",
    )
    parser.add_argument(
        "--docked-root",
        default=None,
        help="Docked root (default: <repo-root>/docked)",
    )
    parser.add_argument(
        "--post-docked-root",
        default=None,
        help="Post-docked root (default: <repo-root>/post_docked)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-run conversions even if outputs already exist",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers for independent conversions (default: 1)",
    )
    return parser.parse_args()


def configure_logging() -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    return logging.getLogger("prep_for_scorch")


def resolve_roots(args: argparse.Namespace) -> Tuple[Path, Path, Path]:
    repo_root = (
        Path(args.repo_root).resolve()
        if args.repo_root
        else Path(__file__).resolve().parents[3]
    )
    docked_root = (
        Path(args.docked_root).resolve() if args.docked_root else repo_root / "docked"
    )
    post_docked_root = (
        Path(args.post_docked_root).resolve()
        if args.post_docked_root
        else repo_root / "post_docked"
    )
    return repo_root, docked_root, post_docked_root


def _find_stage_parts(
    path: Path, run_root: Path, stage_dirs: set[str]
) -> Optional[Tuple[str, str, str, str]]:
    rel_parts = path.relative_to(run_root).parts
    stage_idx = next(
        (idx for idx, part in enumerate(rel_parts) if part in stage_dirs), None
    )
    if stage_idx is None:
        return None
    if stage_idx < 3:
        return None
    pdb_id, variant, ph = rel_parts[stage_idx - 3 : stage_idx]
    stage_dir = rel_parts[stage_idx]
    return pdb_id, variant, ph, stage_dir


def _pose_base_from_name(name: str) -> str:
    stem = name
    stem = re.sub(r"\.(pdbqt|mol2|dok)$", "", stem, flags=re.IGNORECASE)
    stem = stem.replace(".sanitized", "")
    stem = re.sub(r"(_gnina_dud_stage\d+)$", "", stem)
    stem = re.sub(r"(_dud_gnina_stage\d+)$", "", stem)
    stem = re.sub(r"(_dock6_dud_stage\d+)$", "", stem)
    stem = re.sub(r"(_dud_dock6_stage\d+)$", "", stem)
    stem = re.sub(r"(_dud_stage\d+)$", "", stem)
    stem = re.sub(r"(__ledock_stage\d+)$", "", stem)
    stem = re.sub(r"(__dock6_stage\d+)$", "", stem)
    stem = re.sub(r"(_gnina_stage\d+)$", "", stem)
    stem = re.sub(r"(_stage\d+)$", "", stem)
    stem = stem.replace("__", "_")
    stem = re.sub(r"_+$", "", stem)
    return stem


def _load_decoy_bases(
    score_csvs: List[Path], run_root: Optional[Path], logger: logging.Logger
) -> Dict[Tuple[str, str, str], set[str]]:
    decoys: Dict[Tuple[str, str, str], set[str]] = {}
    if not score_csvs:
        return decoys
    for score_csv in score_csvs:
        if not score_csv.exists():
            continue
        combo_from_path: Optional[Tuple[str, str, str]] = None
        if run_root is not None:
            try:
                rel_parts = score_csv.relative_to(run_root).parts
            except ValueError:
                rel_parts = ()
            if len(rel_parts) >= 3:
                combo_from_path = (rel_parts[0], rel_parts[1], rel_parts[2])
        try:
            with score_csv.open("r", newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    pdb_id = str(row.get("pdb_id", "")).strip()
                    variant = str(row.get("variant", "")).strip()
                    ph = str(row.get("ph_label", "") or row.get("ph", "")).strip()
                    if combo_from_path:
                        if not pdb_id:
                            pdb_id = combo_from_path[0]
                        if not variant:
                            variant = combo_from_path[1]
                        if not ph:
                            ph = combo_from_path[2]
                    ligand = str(row.get("ligand", "")).strip()
                    if not (pdb_id and variant):
                        continue
                    base = _pose_base_from_name(Path(ligand).name)
                    if not base:
                        continue
                    decoys.setdefault((pdb_id, variant, ph), set()).add(base)
        except Exception as exc:
            logger.warning(
                "%s action=decoys status=skip reason=read_error path=%s error=%s",
                COMPONENT,
                score_csv,
                exc,
            )
    return decoys


def discover_tasks(
    run_root: Path, post_run_root: Path, logger: logging.Logger
) -> Tuple[List[LedockTask], List[Dock6Task], List[LedockTask], List[Dock6Task]]:
    ledock_tasks: List[LedockTask] = []
    dock6_tasks: List[Dock6Task] = []
    ledock_dud_tasks: List[LedockTask] = []
    dock6_dud_tasks: List[Dock6Task] = []
    ledock_dud_combos: set[Tuple[str, str, str]] = set()

    for dok_path in run_root.rglob("*.dok"):
        stage_parts = _find_stage_parts(
            dok_path,
            run_root,
            LEDOCK_STAGE_DIRS | LEDOCK_DUD_STAGE_DIRS | LEDOCK_DUD_STAGE_DIRS_LEGACY,
        )
        if not stage_parts:
            continue
        pdb_id, variant, ph, stage_dir = stage_parts
        if (
            stage_dir in LEDOCK_DUD_STAGE_DIRS
            or stage_dir in LEDOCK_DUD_STAGE_DIRS_LEGACY
        ):
            output_dir = post_run_root / pdb_id / variant / ph / "dud_ledock_pdbqt"
            ledock_dud_tasks.append(
                LedockTask(
                    input_path=dok_path,
                    output_dir=output_dir,
                    stage_dir=stage_dir,
                    pdb=pdb_id,
                    variant=variant,
                    ph=ph,
                )
            )
            ledock_dud_combos.add((pdb_id, variant, ph))
        else:
            output_dir = post_run_root / pdb_id / variant / ph / "ledock_pdbqt"
            ledock_tasks.append(
                LedockTask(
                    input_path=dok_path,
                    output_dir=output_dir,
                    stage_dir=stage_dir,
                    pdb=pdb_id,
                    variant=variant,
                    ph=ph,
                )
            )

    for mol2_path in run_root.rglob("*.mol2"):
        stage_parts = _find_stage_parts(
            mol2_path,
            run_root,
            DOCK6_STAGE_DIRS | DOCK6_DUD_STAGE_DIRS | DOCK6_DUD_STAGE_DIRS_LEGACY,
        )
        if not stage_parts:
            continue
        pdb_id, variant, ph, stage_dir = stage_parts
        if (
            stage_dir in DOCK6_DUD_STAGE_DIRS
            or stage_dir in DOCK6_DUD_STAGE_DIRS_LEGACY
        ):
            output_dir = post_run_root / pdb_id / variant / ph / "dud_dock6_pdbqt"
            dock6_dud_tasks.append(
                Dock6Task(
                    input_path=mol2_path, output_dir=output_dir, stage_dir=stage_dir
                )
            )
        else:
            output_dir = post_run_root / pdb_id / variant / ph / "dock6_pdbqt"
            dock6_tasks.append(
                Dock6Task(
                    input_path=mol2_path, output_dir=output_dir, stage_dir=stage_dir
                )
            )

    dud_score_csvs = list(run_root.rglob("dud_ledock_docking_score_long.csv"))
    dud_bases_by_combo = _load_decoy_bases(dud_score_csvs, run_root, logger)
    if dud_bases_by_combo:
        for task in ledock_tasks:
            combo = (task.pdb, task.variant, task.ph)
            if combo in ledock_dud_combos:
                continue
            decoy_bases = dud_bases_by_combo.get(combo, set())
            if not decoy_bases:
                continue
            base = _pose_base_from_name(task.input_path.name)
            if base and base in decoy_bases:
                ledock_dud_tasks.append(
                    LedockTask(
                        input_path=task.input_path,
                        output_dir=post_run_root
                        / task.pdb
                        / task.variant
                        / task.ph
                        / "dud_ledock_pdbqt",
                        stage_dir=task.stage_dir,
                        pdb=task.pdb,
                        variant=task.variant,
                        ph=task.ph,
                    )
                )

    logger.info(
        "%s action=discover ledock_found=%d dock6_found=%d dud_ledock_found=%d dud_dock6_found=%d run_root=%s",
        COMPONENT,
        len(ledock_tasks),
        len(dock6_tasks),
        len(ledock_dud_tasks),
        len(dock6_dud_tasks),
        run_root,
    )
    return ledock_tasks, dock6_tasks, ledock_dud_tasks, dock6_dud_tasks


def _ensure_obabel(logger: logging.Logger) -> bool:
    if shutil.which("obabel") is None:
        logger.error(
            "%s action=preflight status=failed reason=missing_obabel", COMPONENT
        )
        return False
    return True


def _validate_output(path: Path) -> bool:
    try:
        return path.exists() and path.stat().st_size > 0
    except OSError:
        return False


def convert_ledock(
    task: LedockTask, overwrite: bool, logger: logging.Logger
) -> Tuple[str, Optional[str]]:
    output_path = task.output_dir / f"{task.input_path.stem}__{task.stage_dir}.pdbqt"
    if output_path.exists() and not overwrite:
        logger.info(
            "%s action=convert status=skip kind=ledock stage=%s input=%s output=%s",
            COMPONENT,
            task.stage_dir,
            task.input_path,
            output_path,
        )
        return "skipped", None

    task.output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "obabel",
        "-ipdb",
        str(task.input_path),
        "-opdbqt",
        "-O",
        str(output_path),
        "-xr",
        "--partialcharge",
        "gasteiger",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.error(
            "%s action=convert status=failed kind=ledock stage=%s input=%s returncode=%s stderr=%s",
            COMPONENT,
            task.stage_dir,
            task.input_path,
            proc.returncode,
            proc.stderr.strip(),
        )
        try:
            output_path.unlink(missing_ok=True)
        except Exception:
            pass
        return "failed", None

    if not _validate_output(output_path):
        logger.error(
            "%s action=convert status=failed kind=ledock stage=%s input=%s reason=empty_output output=%s",
            COMPONENT,
            task.stage_dir,
            task.input_path,
            output_path,
        )
        try:
            output_path.unlink(missing_ok=True)
        except Exception:
            pass
        return "failed", None

    logger.info(
        "%s action=convert status=ok kind=ledock stage=%s input=%s output=%s",
        COMPONENT,
        task.stage_dir,
        task.input_path,
        output_path,
    )
    return "converted", f"{task.input_path} -> {output_path}"


def _sanitize_mol_name(name: str) -> str:
    cleaned = name.strip().replace(" ", "_").replace(os.sep, "_").replace("\\", "_")
    return cleaned or "mol"


def _strip_trailing_mol2(name: str) -> str:
    # Normalize DOCK6 mol names to avoid base mismatches across stage outputs.
    cleaned = name.strip()
    return re.sub(r"\.mol2$", "", cleaned, flags=re.IGNORECASE)


def _safe_unlink(path: Path) -> None:
    try:
        path.unlink(missing_ok=True)
    except Exception:
        pass


def _split_mol2_blocks(mol2_path: Path) -> List[Tuple[str, List[str]]]:
    blocks: List[Tuple[str, List[str]]] = []
    try:
        with mol2_path.open() as handle:
            lines = handle.readlines()
    except OSError:
        return blocks

    current_name: Optional[str] = None
    current_lines: List[str] = []
    idx = 0
    while idx < len(lines):
        line = lines[idx]
        if line.strip() == "@<TRIPOS>MOLECULE":
            if current_name and current_lines:
                blocks.append((current_name, current_lines))
            current_lines = [line]
            if idx + 1 < len(lines):
                current_name = lines[idx + 1].strip() or f"mol_{len(blocks)}"
                current_lines.append(lines[idx + 1])
                idx += 2
                continue
            current_name = f"mol_{len(blocks)}"
        else:
            current_lines.append(line)
        idx += 1

    if current_name and current_lines:
        blocks.append((current_name, current_lines))
    return blocks


def convert_dock6(
    task: Dock6Task, overwrite: bool, logger: logging.Logger
) -> Dock6Result:
    result = Dock6Result()
    task.output_dir.mkdir(parents=True, exist_ok=True)
    blocks = _split_mol2_blocks(task.input_path)
    if not blocks:
        logger.error(
            "%s action=split status=failed kind=dock6 stage=%s input=%s reason=no_blocks",
            COMPONENT,
            task.stage_dir,
            task.input_path,
        )
        result.failed += 1
        return result

    logger.info(
        "%s action=split status=ok kind=dock6 stage=%s input=%s molecules=%d",
        COMPONENT,
        task.stage_dir,
        task.input_path,
        len(blocks),
    )

    name_counts: Dict[str, int] = {}
    for mol_name, lines in blocks:
        result.molecules += 1
        base_name = _sanitize_mol_name(_strip_trailing_mol2(mol_name))
        name_counts[base_name] = name_counts.get(base_name, 0) + 1
        count = name_counts[base_name]
        prefix = f"{base_name}__{task.stage_dir}"
        out_base = prefix if count == 1 else f"{prefix}__dup{count}"
        mol2_path = task.output_dir / f"{out_base}.mol2"
        output_path = task.output_dir / f"{out_base}.pdbqt"

        if output_path.exists() and not overwrite:
            logger.info(
                "%s action=convert status=skip kind=dock6 stage=%s input=%s output=%s",
                COMPONENT,
                task.stage_dir,
                task.input_path,
                output_path,
            )
            result.skipped += 1
            _safe_unlink(mol2_path)
            continue

        try:
            with mol2_path.open("w") as handle:
                handle.writelines(lines)
        except Exception as exc:
            logger.error(
                "%s action=split status=failed kind=dock6 stage=%s input=%s mol_name=%s error=%s",
                COMPONENT,
                task.stage_dir,
                task.input_path,
                out_base,
                exc,
            )
            result.failed += 1
            _safe_unlink(mol2_path)
            continue

        cmd = [
            "obabel",
            "-imol2",
            str(mol2_path),
            "-opdbqt",
            "-O",
            str(output_path),
            "-xr",
            "--partialcharge",
            "gasteiger",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            logger.error(
                "%s action=convert status=failed kind=dock6 stage=%s input=%s mol_name=%s returncode=%s stderr=%s",
                COMPONENT,
                task.stage_dir,
                task.input_path,
                out_base,
                proc.returncode,
                proc.stderr.strip(),
            )
            try:
                output_path.unlink(missing_ok=True)
            except Exception:
                pass
            result.failed += 1
            _safe_unlink(mol2_path)
            continue

        if not _validate_output(output_path):
            logger.error(
                "%s action=convert status=failed kind=dock6 stage=%s input=%s mol_name=%s reason=empty_output",
                COMPONENT,
                task.stage_dir,
                task.input_path,
                out_base,
            )
            try:
                output_path.unlink(missing_ok=True)
            except Exception:
                pass
            result.failed += 1
            _safe_unlink(mol2_path)
            continue

        result.converted += 1
        if result.mappings is not None:
            result.mappings.append(
                f"{task.input_path}::{out_base}.mol2 -> {output_path}"
            )
        _safe_unlink(mol2_path)

    return result


def convert_ledock_via_sdf(
    task: LedockTask, post_run_root: Path, overwrite: bool, logger: logging.Logger
) -> Tuple[str, Optional[str]]:
    output_path = task.output_dir / f"{task.input_path.stem}__{task.stage_dir}.pdbqt"
    if output_path.exists() and not overwrite:
        logger.info(
            "%s action=convert status=skip kind=ledock_fallback stage=%s input=%s output=%s",
            COMPONENT,
            task.stage_dir,
            task.input_path,
            output_path,
        )
        return "skipped", None

    sdf_path = (
        post_run_root
        / task.pdb
        / task.variant
        / task.ph
        / task.stage_dir
        / f"{task.input_path.stem}.sdf"
    )
    if not sdf_path.exists() or not _validate_output(sdf_path):
        logger.error(
            "%s action=convert status=failed kind=ledock_fallback stage=%s input=%s reason=missing_sdf sdf=%s",
            COMPONENT,
            task.stage_dir,
            task.input_path,
            sdf_path,
        )
        return "failed", None

    cmd = [
        "obabel",
        "-isdf",
        str(sdf_path),
        "-opdbqt",
        "-O",
        str(output_path),
        "-xr",
        "--partialcharge",
        "gasteiger",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not _validate_output(output_path):
        logger.warning(
            "%s action=convert status=retry kind=ledock_fallback stage=%s input=%s reason=%s stderr=%s",
            COMPONENT,
            task.stage_dir,
            task.input_path,
            "gasteiger_failed" if proc.returncode != 0 else "empty_output",
            proc.stderr.strip(),
        )
        _safe_unlink(output_path)
        retry_cmd = [
            "obabel",
            "-isdf",
            str(sdf_path),
            "-opdbqt",
            "-O",
            str(output_path),
            "-xr",
        ]
        retry_proc = subprocess.run(retry_cmd, capture_output=True, text=True)
        if retry_proc.returncode != 0 or not _validate_output(output_path):
            logger.error(
                "%s action=convert status=failed kind=ledock_fallback stage=%s input=%s returncode=%s stderr=%s",
                COMPONENT,
                task.stage_dir,
                task.input_path,
                retry_proc.returncode,
                retry_proc.stderr.strip(),
            )
            _safe_unlink(output_path)
            return "failed", None

    logger.info(
        "%s action=convert status=ok kind=ledock_fallback stage=%s input=%s output=%s",
        COMPONENT,
        task.stage_dir,
        task.input_path,
        output_path,
    )
    return "converted", f"{task.input_path} -> {output_path}"


def _run_ledock_tasks(
    tasks: List[LedockTask], overwrite: bool, workers: int, logger: logging.Logger
) -> Tuple[int, int, int, List[str], List[LedockTask]]:
    converted = 0
    skipped = 0
    failed = 0
    mappings: List[str] = []
    failed_tasks: List[LedockTask] = []

    def _runner(task: LedockTask) -> Tuple[str, Optional[str]]:
        try:
            return convert_ledock(task, overwrite, logger)
        except Exception as exc:  # defensive guard
            logger.error(
                "%s action=convert status=failed kind=ledock stage=%s input=%s reason=exception error=%s",
                COMPONENT,
                task.stage_dir,
                task.input_path,
                exc,
            )
            return "failed", None

    if workers <= 1:
        results = [(_runner(task), task) for task in tasks]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            future_map = {pool.submit(_runner, task): task for task in tasks}
            results = []
            for future in as_completed(future_map):
                task = future_map[future]
                try:
                    results.append((future.result(), task))
                except Exception as exc:
                    logger.error(
                        "%s action=convert status=failed kind=ledock stage=%s input=%s reason=worker_exception error=%s",
                        COMPONENT,
                        task.stage_dir,
                        task.input_path,
                        exc,
                    )
                    results.append((("failed", None), task))

    for (status, mapping), task in results:
        if status == "converted":
            converted += 1
        elif status == "skipped":
            skipped += 1
        else:
            failed += 1
            failed_tasks.append(task)
        if mapping and len(mappings) < EXAMPLE_LIMIT:
            mappings.append(mapping)

    return converted, skipped, failed, mappings, failed_tasks


def _run_dock6_tasks(
    tasks: List[Dock6Task], overwrite: bool, workers: int, logger: logging.Logger
) -> Tuple[int, int, int, int, List[str]]:
    molecules = 0
    converted = 0
    skipped = 0
    failed = 0
    mappings: List[str] = []

    def _runner(task: Dock6Task) -> Dock6Result:
        try:
            return convert_dock6(task, overwrite, logger)
        except Exception as exc:  # defensive guard
            logger.error(
                "%s action=convert status=failed kind=dock6 stage=%s input=%s reason=exception error=%s",
                COMPONENT,
                task.stage_dir,
                task.input_path,
                exc,
            )
            res = Dock6Result()
            res.failed = 1
            return res

    if workers <= 1:
        results = [_runner(task) for task in tasks]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            future_map = {pool.submit(_runner, task): task for task in tasks}
            results = []
            for future in as_completed(future_map):
                try:
                    results.append(future.result())
                except Exception as exc:
                    task = future_map[future]
                    logger.error(
                        "%s action=convert status=failed kind=dock6 stage=%s input=%s reason=worker_exception error=%s",
                        COMPONENT,
                        task.stage_dir,
                        task.input_path,
                        exc,
                    )
                    res = Dock6Result()
                    res.failed = 1
                    results.append(res)

    for res in results:
        molecules += res.molecules
        converted += res.converted
        skipped += res.skipped
        failed += res.failed
        if res.mappings:
            for mapping in res.mappings:
                if len(mappings) < EXAMPLE_LIMIT:
                    mappings.append(mapping)

    return molecules, converted, skipped, failed, mappings


def run_pose_bust_conversion(
    run_id: str, repo_root: Path, overwrite: bool, logger: logging.Logger
) -> bool:
    pose_bust_path = repo_root / "pose_bust.py"
    if not pose_bust_path.exists():
        logger.error(
            "%s action=fallback status=failed reason=missing_pose_bust path=%s",
            COMPONENT,
            pose_bust_path,
        )
        return False

    cmd = [
        sys.executable,
        str(pose_bust_path),
        "--run-id",
        run_id,
        "--skip-posebusters",
    ]
    if overwrite:
        cmd.append("--overwrite")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.error(
            "%s action=fallback status=failed reason=pose_bust_error returncode=%s stderr=%s stdout=%s",
            COMPONENT,
            proc.returncode,
            proc.stderr.strip(),
            proc.stdout.strip(),
        )
        return False

    logger.info(
        "%s action=fallback status=pose_bust_ok run_id=%s overwrite=%s",
        COMPONENT,
        run_id,
        overwrite,
    )
    return True


def _cleanup_residual_mol2(post_run_root: Path, logger: logging.Logger) -> None:
    residuals = list(post_run_root.rglob("dock6_pdbqt/*.mol2"))
    residuals += list(post_run_root.rglob("dud_dock6_pdbqt/*.mol2"))
    if not residuals:
        return
    for path in residuals:
        _safe_unlink(path)
    logger.info(
        "%s action=cleanup kind=dock6 temp_mol2_removed=%d root=%s",
        COMPONENT,
        len(residuals),
        post_run_root,
    )


def _migrate_dock6_pdbqt_names(
    post_run_root: Path, overwrite: bool, logger: logging.Logger
) -> None:
    if not overwrite or not post_run_root.exists():
        return

    renamed = 0
    overwritten = 0
    skipped = 0
    failed = 0
    for path in list(post_run_root.rglob("dock6_pdbqt/*.pdbqt")) + list(
        post_run_root.rglob("dud_dock6_pdbqt/*.pdbqt")
    ):
        if ".mol2__dock6_stage" not in path.name.lower():
            continue
        new_name = re.sub(
            r"\.mol2(?=__dock6_stage\d+)", "", path.name, flags=re.IGNORECASE
        )
        if new_name == path.name:
            continue
        new_path = path.with_name(new_name)
        if new_path.exists():
            try:
                new_path.unlink()
            except Exception as exc:
                failed += 1
                logger.warning(
                    "%s action=migrate status=failed kind=dock6 reason=unlink_failed source=%s target=%s error=%s | remove collision",
                    COMPONENT,
                    path,
                    new_path,
                    exc,
                )
                continue
            overwritten += 1
        try:
            path.replace(new_path)
            renamed += 1
        except Exception as exc:
            failed += 1
            logger.warning(
                "%s action=migrate status=failed kind=dock6 reason=rename_failed source=%s target=%s error=%s | normalize dock6 name",
                COMPONENT,
                path,
                new_path,
                exc,
            )

    if renamed or overwritten or skipped or failed:
        logger.info(
            "%s action=migrate status=ok kind=dock6 renamed=%d overwritten=%d skipped=%d failed=%d root=%s | normalize dock6 names",
            COMPONENT,
            renamed,
            overwritten,
            skipped,
            failed,
            post_run_root,
        )


def _run_ledock_with_fallback(
    tasks: List[LedockTask],
    run_id: str,
    repo_root: Path,
    post_run_root: Path,
    overwrite: bool,
    workers: int,
    logger: logging.Logger,
    label: str,
) -> Tuple[int, int, int, List[str]]:
    (
        ledock_converted,
        ledock_skipped,
        ledock_failed,
        ledock_mappings,
        ledock_failed_tasks,
    ) = _run_ledock_tasks(tasks, overwrite, max(1, workers), logger)

    fallback_targets: List[LedockTask] = ledock_failed_tasks
    logger.info(
        "%s action=fallback status=collect kind=%s failed_tasks=%d",
        COMPONENT,
        label,
        len(fallback_targets),
    )
    if not fallback_targets:
        for task in tasks:
            out_path = (
                task.output_dir / f"{task.input_path.stem}__{task.stage_dir}.pdbqt"
            )
            if not _validate_output(out_path):
                fallback_targets.append(task)

    ledock_failed = len(fallback_targets)
    if fallback_targets:
        logger.info(
            "%s action=fallback status=precheck kind=%s ledock_failed=%d candidates=%d",
            COMPONENT,
            label,
            ledock_failed,
            len(fallback_targets),
        )

        logger.warning(
            "%s action=fallback status=triggered kind=%s failed_count=%d",
            COMPONENT,
            label,
            len(fallback_targets),
        )
        logger.info(
            "%s action=fallback status=candidates kind=%s examples=%s",
            COMPONENT,
            label,
            [task.input_path.name for task in fallback_targets[:EXAMPLE_LIMIT]],
        )
        fallback_ok = run_pose_bust_conversion(run_id, repo_root, overwrite, logger)
        if not fallback_ok:
            logger.error(
                "%s action=fallback status=failed kind=%s reason=pose_bust_run failed_count=%d",
                COMPONENT,
                label,
                ledock_failed,
            )
            logger.warning(
                "%s action=fallback status=continue_after_pose_bust_failure kind=%s candidates=%d",
                COMPONENT,
                label,
                len(fallback_targets),
            )

        if fallback_targets:
            still_failed: List[LedockTask] = []
            converted_via_fallback = 0
            for task in fallback_targets:
                status, mapping = convert_ledock_via_sdf(
                    task, post_run_root, overwrite, logger
                )
                if status == "converted":
                    ledock_converted += 1
                    converted_via_fallback += 1
                elif status == "skipped":
                    ledock_skipped += 1
                else:
                    still_failed.append(task)
                if mapping and len(ledock_mappings) < EXAMPLE_LIMIT:
                    ledock_mappings.append(mapping)

            ledock_failed = len(still_failed)
            logger.info(
                "%s action=fallback status=ok kind=%s converted=%d still_failed=%d pose_bust_ok=%s",
                COMPONENT,
                label,
                converted_via_fallback,
                ledock_failed,
                fallback_ok,
            )
    else:
        ledock_failed = 0

    return ledock_converted, ledock_skipped, ledock_failed, ledock_mappings


def main() -> int:
    args = parse_args()
    logger = configure_logging()
    repo_root, docked_root, post_docked_root = resolve_roots(args)

    if not _ensure_obabel(logger):
        return 1

    run_root = docked_root / args.run_id
    post_run_root = post_docked_root / args.run_id
    if not run_root.exists():
        logger.error(
            "%s action=preflight status=failed reason=missing_run_root run_root=%s",
            COMPONENT,
            run_root,
        )
        return 1

    _migrate_dock6_pdbqt_names(post_run_root, args.overwrite, logger)

    ledock_tasks, dock6_tasks, ledock_dud_tasks, dock6_dud_tasks = discover_tasks(
        run_root, post_run_root, logger
    )
    if ledock_dud_tasks:
        logger.info(
            "[prep_for_scorch.dud] source=ledock found_in=%s n_in=%d out_dir=%s",
            ledock_dud_tasks[0].input_path.parent,
            len(ledock_dud_tasks),
            ledock_dud_tasks[0].output_dir,
        )
    else:
        logger.info(
            "[prep_for_scorch.dud] source=ledock status=skip reason=no_inputs run_root=%s",
            run_root,
        )
    if dock6_dud_tasks:
        logger.info(
            "[prep_for_scorch.dud] source=dock6 found_in=%s n_in=%d out_dir=%s",
            dock6_dud_tasks[0].input_path.parent,
            len(dock6_dud_tasks),
            dock6_dud_tasks[0].output_dir,
        )
    else:
        logger.info(
            "[prep_for_scorch.dud] source=dock6 status=skip reason=no_inputs run_root=%s",
            run_root,
        )

    (
        ledock_dud_converted,
        ledock_dud_skipped,
        ledock_dud_failed,
        ledock_dud_mappings,
    ) = _run_ledock_with_fallback(
        ledock_dud_tasks,
        args.run_id,
        repo_root,
        post_run_root,
        args.overwrite,
        max(1, args.workers),
        logger,
        "dud_ledock",
    )
    (
        ledock_converted,
        ledock_skipped,
        ledock_failed,
        ledock_mappings,
    ) = _run_ledock_with_fallback(
        ledock_tasks,
        args.run_id,
        repo_root,
        post_run_root,
        args.overwrite,
        max(1, args.workers),
        logger,
        "ledock",
    )

    (
        dock6_dud_molecules,
        dock6_dud_converted,
        dock6_dud_skipped,
        dock6_dud_failed,
        dock6_dud_mappings,
    ) = _run_dock6_tasks(dock6_dud_tasks, args.overwrite, max(1, args.workers), logger)
    (
        dock6_molecules,
        dock6_converted,
        dock6_skipped,
        dock6_failed,
        dock6_mappings,
    ) = _run_dock6_tasks(dock6_tasks, args.overwrite, max(1, args.workers), logger)

    _cleanup_residual_mol2(post_run_root, logger)

    if ledock_mappings or dock6_mappings or ledock_dud_mappings or dock6_dud_mappings:
        print("Example conversions:")
        for mapping in (
            ledock_mappings + dock6_mappings + ledock_dud_mappings + dock6_dud_mappings
        ):
            print(f"  {mapping}")

    ledock_found = len(ledock_tasks)
    dock6_found = len(dock6_tasks)
    dud_ledock_found = len(ledock_dud_tasks)
    dud_dock6_found = len(dock6_dud_tasks)
    print(
        "ledock_found={found} ledock_converted={conv} ledock_skipped={skip} ledock_failed={fail} dud_ledock_found={dud_found} dud_ledock_converted={dud_conv} dud_ledock_skipped={dud_skip} dud_ledock_failed={dud_fail}".format(
            found=ledock_found,
            conv=ledock_converted,
            skip=ledock_skipped,
            fail=ledock_failed,
            dud_found=dud_ledock_found,
            dud_conv=ledock_dud_converted,
            dud_skip=ledock_dud_skipped,
            dud_fail=ledock_dud_failed,
        )
    )
    print(
        "dock6_mol2_found={found} dock6_molecules_split={mols} dock6_converted={conv} dock6_skipped={skip} dock6_failed={fail} dud_dock6_found={dud_found} dud_dock6_molecules_split={dud_mols} dud_dock6_converted={dud_conv} dud_dock6_skipped={dud_skip} dud_dock6_failed={dud_fail}".format(
            found=dock6_found,
            mols=dock6_molecules,
            conv=dock6_converted,
            skip=dock6_skipped,
            fail=dock6_failed,
            dud_found=dud_dock6_found,
            dud_mols=dock6_dud_molecules,
            dud_conv=dock6_dud_converted,
            dud_skip=dock6_dud_skipped,
            dud_fail=dock6_dud_failed,
        )
    )

    has_failures = (
        ledock_failed > 0
        or dock6_failed > 0
        or ledock_dud_failed > 0
        or dock6_dud_failed > 0
    )
    if has_failures:
        logger.error("%s action=summary status=failed", COMPONENT)
    else:
        logger.info("%s action=summary status=ok", COMPONENT)
    return 1 if has_failures else 0


if __name__ == "__main__":
    sys.exit(main())
