from __future__ import annotations

import argparse
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd
from input_and_export_functions import load_config


STAGE_DIRS = {
    "stage1",
    "stage2",
    "stage3",
    "gnina_stage1",
    "gnina_stage2",
    "gnina_stage3",
    "ledock_stage1",
    "ledock_stage2",
    "ledock_stage3",
    "dock6_stage1",
    "dock6_stage2",
    "dock6_stage3",
}
EXCLUDE_NAMES = {"receptor.pdbqt", "protein.pdbqt", "receptor.pdb", "protein.pdb"}
ENGINE_BY_STAGE_DIR = {
    "stage1": "vina",
    "stage2": "vina",
    "stage3": "vina",
    "gnina_stage1": "gnina",
    "gnina_stage2": "gnina",
    "gnina_stage3": "gnina",
    "ledock_stage1": "ledock",
    "ledock_stage2": "ledock",
    "ledock_stage3": "ledock",
    "dock6_stage1": "dock6",
    "dock6_stage2": "dock6",
    "dock6_stage3": "dock6",
}
STAGE_NORMALIZED = {
    "stage1": "stage1",
    "stage2": "stage2",
    "stage3": "stage3",
    "gnina_stage1": "stage1",
    "gnina_stage2": "stage2",
    "gnina_stage3": "stage3",
    "ledock_stage1": "stage1",
    "ledock_stage2": "stage2",
    "ledock_stage3": "stage3",
    "dock6_stage1": "stage1",
    "dock6_stage2": "stage2",
    "dock6_stage3": "stage3",
}
VARIANT_CANON: Dict[str, str] = {}


@dataclass(frozen=True)
class LigandTask:
    input_path: Path
    output_path: Path
    stage: str
    kind: str  # "pdbqt" or "dok"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert docked ligand PDBQT and LeDock DOK files to SDF under post_docked."
    )
    parser.add_argument("--run-id", required=True, help="Run identifier under docked/")
    parser.add_argument(
        "--repo-root",
        default=None,
        help="Repository root (default: directory containing pose_bust.py)",
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
        help="Reconvert even if output already exists",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help="Number of parallel workers for obabel conversions (0 uses config-based max)",
    )
    parser.add_argument(
        "--run-posebusters",
        dest="run_posebusters",
        action="store_true",
        default=True,
        help="Run PoseBusters on converted ligands (default: enabled)",
    )
    parser.add_argument(
        "--skip-posebusters",
        dest="run_posebusters",
        action="store_false",
        help="Skip PoseBusters validation",
    )
    parser.add_argument(
        "--posebusters-relative-cutoff",
        type=float,
        default=0.92,
        help="Relative distance cutoff for PoseBusters pass criterion",
    )
    parser.add_argument(
        "--posebusters-max-workers",
        type=int,
        default=0,
        help="Max workers to pass to PoseBusters (0 or None lets PoseBusters decide)",
    )
    parser.add_argument(
        "--posebusters-chunk",
        type=int,
        default=0,
        help="Chunk size for PoseBusters inputs (0 = no chunking; send all SDFs in a single bust call)",
    )
    parser.add_argument(
        "--include-dock6",
        dest="include_dock6",
        action="store_true",
        default=True,
        help="Include DOCK6 mol2 conversions (default: enabled)",
    )
    parser.add_argument(
        "--skip-dock6",
        dest="include_dock6",
        action="store_false",
        help="Skip DOCK6 mol2 conversions",
    )
    return parser.parse_args()


def configure_logging() -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    return logging.getLogger("pose_bust")


def _find_stage(rel_parts: Tuple[str, ...]) -> str | None:
    for part in reversed(rel_parts[:-1]):
        if part in STAGE_DIRS:
            return part
    return None


def discover_inputs(run_root: Path, post_root: Path) -> List[LigandTask]:
    tasks: List[LigandTask] = []
    for path in (
        list(run_root.rglob("*.pdbqt"))
        + list(run_root.rglob("*.dok"))
        + list(run_root.rglob("*.mol2"))
    ):
        if path.name in EXCLUDE_NAMES:
            logging.info("[pose-bust] action=skip reason=excluded input=%s", path)
            continue

        rel_parts = path.relative_to(run_root).parts
        stage = _find_stage(rel_parts)
        if stage is None:
            continue
        # Only accept dock6 mol2 files inside dock6 stage dirs
        if path.suffix.lower() == ".mol2" and not stage.startswith("dock6_stage"):
            continue

        rel_path = path.relative_to(run_root)
        out_path = (post_root / rel_path).with_suffix(".sdf")
        if path.suffix.lower() == ".dok":
            kind = "dok"
        elif path.suffix.lower() == ".mol2":
            kind = "dock6"
        else:
            kind = "pdbqt"
        tasks.append(
            LigandTask(input_path=path, output_path=out_path, stage=stage, kind=kind)
        )

    return tasks


def _parse_pose_pdb(pdb_path: Path) -> Tuple[int, List[Tuple[float, float, float]]]:
    atoms: List[Tuple[float, float, float]] = []
    with pdb_path.open() as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((x, y, z))
            except ValueError:
                continue
    return len(atoms), atoms


def _validate_ledock_identity(
    pdb_files: Sequence[Path], sdf_path: Path, tol: float = 1e-2
) -> bool:
    try:
        from rdkit import Chem  # type: ignore
    except Exception as exc:
        logging.error(
            "[pose-bust] action=validate status=failed reason=missing_rdkit sdf=%s error=%s",
            sdf_path,
            exc,
        )
        return False

    suppl = Chem.SDMolSupplier(
        str(sdf_path), removeHs=False, sanitize=False, strictParsing=False
    )
    mols = []
    for m in suppl:
        if m is None:
            logging.error(
                "[pose-bust] action=validate status=failed reason=mol_load_failed sdf=%s",
                sdf_path,
            )
            return False
        mols.append(m)

    if len(mols) != len(pdb_files):
        logging.error(
            "[pose-bust] action=validate status=failed reason=count_mismatch sdf=%s poses=%d mols=%d",
            sdf_path,
            len(pdb_files),
            len(mols),
        )
        return False

    for idx, (pdb_file, mol) in enumerate(zip(pdb_files, mols)):
        pdb_atoms, pdb_coords = _parse_pose_pdb(pdb_file)
        try:
            conf = mol.GetConformer()
        except Exception as exc:
            logging.error(
                "[pose-bust] action=validate status=failed reason=missing_conformer pose_idx=%d pdb=%s error=%s",
                idx,
                pdb_file,
                exc,
            )
            return False
        sdf_coords = [
            (
                conf.GetAtomPosition(i).x,
                conf.GetAtomPosition(i).y,
                conf.GetAtomPosition(i).z,
            )
            for i in range(mol.GetNumAtoms())
        ]
        if pdb_atoms != len(sdf_coords):
            logging.error(
                "[pose-bust] action=validate status=failed reason=atom_count_mismatch pose_idx=%d pdb=%s pdb_atoms=%d sdf_atoms=%d",
                idx,
                pdb_file,
                pdb_atoms,
                len(sdf_coords),
            )
            return False
        max_diff = 0.0
        for (px, py, pz), (sx, sy, sz) in zip(pdb_coords, sdf_coords):
            max_diff = max(max_diff, abs(px - sx), abs(py - sy), abs(pz - sz))
        if max_diff > tol:
            logging.error(
                "[pose-bust] action=validate status=failed reason=coord_mismatch pose_idx=%d pdb=%s max_diff=%.4f tol=%.4f",
                idx,
                pdb_file,
                max_diff,
                tol,
            )
            return False
    return True


def _collect_pose_pdbs(out_dir: Path, stem: str) -> List[Path]:
    pose_files = sorted(out_dir.glob(f"{stem}_dock*.pdb"))
    if not pose_files:
        pose_files = sorted(
            p
            for p in out_dir.glob(f"{stem}*.pdb")
            if p.name not in EXCLUDE_NAMES and p.is_file()
        )
    return pose_files


def _convert_dok(task: LigandTask, overwrite: bool) -> str:
    output = task.output_path
    out_dir = output.parent
    stem = task.input_path.stem

    if output.exists() and not overwrite:
        logging.info(
            "[pose-bust] action=skip reason=exists stage=%s kind=dok input=%s output=%s",
            task.stage,
            task.input_path,
            output,
        )
        return "skipped"

    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob(f"{stem}_dock*.pdb"):
        try:
            stale.unlink()
        except Exception:
            pass
    for stale in out_dir.glob(f"{stem}*.pdb"):
        try:
            stale.unlink()
        except Exception:
            pass

    local_dok = out_dir / task.input_path.name

    def _cleanup_local_dok():
        try:
            local_dok.unlink(missing_ok=True)
        except Exception:
            pass

    def _cleanup_pose_files(pose_paths: Sequence[Path] | None = None):
        paths = list(pose_paths) if pose_paths else []
        if not paths:
            paths = list(out_dir.glob(f"{stem}_dock*.pdb")) + list(
                out_dir.glob(f"{stem}*.pdb")
            )
        for p in paths:
            try:
                p.unlink(missing_ok=True)
            except Exception:
                pass

    try:
        if overwrite or not local_dok.exists():
            shutil.copy2(task.input_path, local_dok)
    except Exception as exc:
        logging.error(
            "[pose-bust] action=split status=failed stage=%s input=%s reason=copy_failed error=%s",
            task.stage,
            task.input_path,
            exc,
        )
        _cleanup_local_dok()
        return "failed"

    split_cmd = ["ledock", "-spli", local_dok.name]
    split_proc = subprocess.run(split_cmd, capture_output=True, text=True, cwd=out_dir)
    if split_proc.returncode != 0:
        logging.warning(
            "[pose-bust] action=split status=nonzero_return stage=%s input=%s returncode=%s stderr=%s stdout=%s",
            task.stage,
            task.input_path,
            split_proc.returncode,
            split_proc.stderr.strip(),
            split_proc.stdout.strip(),
        )

    pose_files = _collect_pose_pdbs(out_dir, stem)
    if not pose_files:
        logging.warning(
            "[pose-bust] action=split status=skip reason=no_pose_files stage=%s input=%s out_dir=%s",
            task.stage,
            task.input_path,
            out_dir,
        )
        _cleanup_local_dok()
        _cleanup_pose_files()
        return "skipped"

    obabel_cmd = [
        "obabel",
        "-ipdb",
        *[str(p) for p in pose_files],
        "-osdf",
        "-O",
        str(output),
    ]
    convert_proc = subprocess.run(obabel_cmd, capture_output=True, text=True)
    if convert_proc.returncode != 0:
        logging.warning(
            "[pose-bust] action=convert status=skip stage=%s kind=dok input=%s returncode=%s stderr=%s",
            task.stage,
            task.input_path,
            convert_proc.returncode,
            convert_proc.stderr.strip(),
        )
        try:
            output.unlink(missing_ok=True)
        except Exception:
            pass
        _cleanup_local_dok()
        _cleanup_pose_files(pose_files)
        return "skipped"

    if not output.exists() or output.stat().st_size == 0:
        logging.warning(
            "[pose-bust] action=convert status=skip reason=empty_output stage=%s input=%s output=%s",
            task.stage,
            task.input_path,
            output,
        )
        try:
            output.unlink(missing_ok=True)
        except Exception:
            pass
        _cleanup_local_dok()
        _cleanup_pose_files(pose_files)
        return "skipped"

    if not _validate_ledock_identity(pose_files, output, tol=1e-2):
        try:
            output.unlink(missing_ok=True)
        except Exception:
            pass
        _cleanup_local_dok()
        _cleanup_pose_files(pose_files)
        return "skipped"

    logging.info(
        "[pose-bust] action=convert status=ok stage=%s kind=dok input=%s output=%s poses=%d",
        task.stage,
        task.input_path,
        output,
        len(pose_files),
    )
    _cleanup_local_dok()
    _cleanup_pose_files(pose_files)
    return "converted"


def run_conversion(task: LigandTask, overwrite: bool) -> str:
    if task.kind == "dok":
        return _convert_dok(task, overwrite)
    if task.kind == "dock6":
        return _convert_dock6(task, overwrite)

    output = task.output_path
    if output.exists() and not overwrite:
        logging.info(
            "[pose-bust] action=skip reason=exists stage=%s kind=pdbqt input=%s output=%s",
            task.stage,
            task.input_path,
            output,
        )
        return "skipped"

    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["obabel", "-ipdbqt", str(task.input_path), "-osdf", "-O", str(output)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logging.error(
            "[pose-bust] action=convert status=failed stage=%s kind=pdbqt input=%s returncode=%s stderr=%s",
            task.stage,
            task.input_path,
            proc.returncode,
            proc.stderr.strip(),
        )
        try:
            output.unlink(missing_ok=True)
        except Exception:
            pass
        return "failed"

    try:
        if not output.exists() or output.stat().st_size == 0:
            logging.error(
                "[pose-bust] action=convert status=failed stage=%s kind=pdbqt input=%s reason=empty_output output=%s",
                task.stage,
                task.input_path,
                output,
            )
            return "failed"
    except OSError as exc:
        logging.error(
            "[pose-bust] action=convert status=failed stage=%s kind=pdbqt input=%s reason=stat_error error=%s",
            task.stage,
            task.input_path,
            exc,
        )
        return "failed"

    logging.debug(
        "[pose-bust] action=convert status=ok stage=%s kind=pdbqt input=%s output=%s",
        task.stage,
        task.input_path,
        output,
    )
    return "converted"


def process_tasks(
    tasks: Iterable[LigandTask], overwrite: bool, workers: int
) -> Dict[LigandTask, str]:
    statuses: Dict[LigandTask, str] = {}
    workers = max(1, workers)

    if workers == 1:
        for task in tasks:
            statuses[task] = run_conversion(task, overwrite)
        return statuses

    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_map = {executor.submit(run_conversion, t, overwrite): t for t in tasks}
        for future in as_completed(future_map):
            task = future_map[future]
            try:
                statuses[task] = future.result()
            except Exception as exc:  # defensive guard for unexpected worker errors
                logging.error(
                    "[pose-bust] action=convert status=failed stage=%s input=%s reason=exception error=%s",
                    task.stage,
                    task.input_path,
                    exc,
                )
                statuses[task] = "failed"
    return statuses


def validate_outputs(tasks: Iterable[LigandTask]) -> List[LigandTask]:
    missing: List[LigandTask] = []
    for task in tasks:
        if task.kind == "dock6":
            stage_dir = task.output_path.parent
            sdfs = list(stage_dir.glob("*.sdf"))
            if not sdfs:
                missing.append(task)
                logging.error(
                    "[pose-bust] status=missing stage=%s input=%s output_dir=%s exists=False size_ok=False",
                    task.stage,
                    task.input_path,
                    stage_dir,
                )
            continue
        try:
            exists = task.output_path.exists()
            size_ok = task.output_path.stat().st_size > 0 if exists else False
        except OSError:
            exists = False
            size_ok = False

        if not exists or not size_ok:
            missing.append(task)
            logging.error(
                "[pose-bust] status=missing stage=%s input=%s output=%s exists=%s size_ok=%s",
                task.stage,
                task.input_path,
                task.output_path,
                exists,
                size_ok,
            )
    return missing


def _canonical_variant(variant: str) -> str:
    return VARIANT_CANON.get(variant, variant)


def resolve_receptor(pdb_id: str, variant: str, ph: str, repo_root: Path) -> Path:
    variant_canon = _canonical_variant(variant)
    return (
        repo_root
        / "processed_pdbs"
        / pdb_id
        / variant_canon
        / "receptor"
        / "ph_ensemble"
        / f"{pdb_id}_{ph}.withH.pdb"
    )


def _as_bool(val) -> bool:
    if isinstance(val, bool):
        return val
    if val is None:
        return False
    if isinstance(val, (int, float)):
        if isinstance(val, float) and math.isnan(val):
            return False
        return val != 0
    sval = str(val).strip().lower()
    return sval in {"true", "1", "yes", "y", "t"}


def _as_float(val) -> float:
    try:
        return float(val)
    except Exception:
        return float("nan")


def _compute_posebusters_pass(row: pd.Series, cutoff: float) -> bool:
    mol_loaded = _as_bool(row.get("mol_pred_loaded"))
    sanitization = _as_bool(row.get("sanitization"))
    passes_valence_checks = _as_bool(row.get("passes_valence_checks"))
    internal_steric_clash = _as_bool(row.get("internal_steric_clash"))
    bond_lengths = _as_bool(row.get("bond_lengths"))
    bond_angles = _as_bool(row.get("bond_angles"))
    most_extreme_clash_protein = _as_bool(row.get("most_extreme_clash_protein"))
    rel_dist = _as_float(row.get("most_extreme_relative_distance_protein"))
    if math.isnan(rel_dist):
        rel_dist = float("-inf")

    return (
        mol_loaded
        and sanitization
        and passes_valence_checks
        and internal_steric_clash
        and bond_lengths
        and bond_angles
        and (most_extreme_clash_protein or rel_dist >= cutoff)
    )


def _chunk_list(items: Sequence[Path], size: int) -> List[List[Path]]:
    if size <= 0:
        return [list(items)]
    return [list(items[i : i + size]) for i in range(0, len(items), size)]


def _split_dock6_mol2(input_path: Path) -> List[Tuple[str, List[str]]]:
    blocks: List[Tuple[str, List[str]]] = []
    current: List[str] = []
    current_name: str | None = None
    with input_path.open() as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip() == "@<TRIPOS>MOLECULE":
            if current and current_name:
                blocks.append((current_name, current))
            current = [line]
            if i + 1 < len(lines):
                current_name = lines[i + 1].strip() or f"mol_{len(blocks)}"
                current.append(lines[i + 1])
                i += 2
                continue
            else:
                current_name = f"mol_{len(blocks)}"
        elif current is not None:
            current.append(line)
        i += 1
    if current and current_name:
        blocks.append((current_name, current))
    return blocks


def _sanitize_mol2_name(name: str) -> str:
    candidate = name.strip()
    for ext in (".mol2", ".sdf"):
        if candidate.lower().endswith(ext):
            candidate = candidate[: -len(ext)]
    candidate = Path(candidate).stem or "mol"
    candidate = candidate.replace(os.sep, "_").replace("\\", "_")
    return candidate


def _cleanup_stage_mol2(stage_dir: Path) -> None:
    for mol2 in stage_dir.glob("*.mol2"):
        try:
            mol2.unlink()
        except Exception:
            pass


def _convert_dock6(task: LigandTask, overwrite: bool) -> str:
    # Skip if already converted unless overwrite
    stage_dir = task.output_path.parent
    if not overwrite:
        existing = list(stage_dir.glob("*.sdf"))
        if existing:
            _cleanup_stage_mol2(stage_dir)
            logging.info(
                "[pose-bust] action=skip reason=exists stage=%s kind=dock6 input=%s sdfs=%d",
                task.stage,
                task.input_path,
                len(existing),
            )
            return "skipped"

    stage_dir.mkdir(parents=True, exist_ok=True)
    blocks = _split_dock6_mol2(task.input_path)
    if not blocks:
        logging.warning(
            "[pose-bust] action=convert status=skip reason=no_blocks stage=%s input=%s",
            task.stage,
            task.input_path,
        )
        return "skipped"

    mol_counts = 0
    name_counts: Dict[str, int] = {}
    failed_blocks = 0
    for name, lines in blocks:
        base_name = _sanitize_mol2_name(name if name else f"mol_{mol_counts}")
        count = name_counts.get(base_name, 0) + 1
        name_counts[base_name] = count
        out_base = f"{base_name}" if count == 1 else f"{base_name}__dup{count}"
        mol2_path = stage_dir / f"{out_base}.mol2"
        sdf_path = stage_dir / f"{out_base}.sdf"

        try:
            with mol2_path.open("w") as f:
                f.writelines(lines)
        except Exception as exc:
            logging.error(
                "[pose-bust] action=split status=failed stage=%s input=%s mol_name=%s error=%s",
                task.stage,
                task.input_path,
                out_base,
                exc,
            )
            failed_blocks += 1
            continue

        if sdf_path.exists() and not overwrite:
            mol_counts += 1
            continue

        cmd = ["obabel", "-imol2", str(mol2_path), "-osdf", "-O", str(sdf_path)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            logging.warning(
                "[pose-bust] action=convert status=skip stage=%s kind=dock6 input=%s mol_name=%s returncode=%s stderr=%s",
                task.stage,
                task.input_path,
                out_base,
                proc.returncode,
                proc.stderr.strip(),
            )
            try:
                sdf_path.unlink(missing_ok=True)
            except Exception:
                pass
            failed_blocks += 1
            continue

        try:
            if not sdf_path.exists() or sdf_path.stat().st_size == 0:
                logging.warning(
                    "[pose-bust] action=convert status=skip reason=empty_output stage=%s kind=dock6 input=%s mol_name=%s",
                    task.stage,
                    task.input_path,
                    out_base,
                )
                failed_blocks += 1
                continue
        except OSError as exc:
            logging.warning(
                "[pose-bust] action=convert status=skip reason=stat_error stage=%s kind=dock6 input=%s mol_name=%s error=%s",
                task.stage,
                task.input_path,
                out_base,
                exc,
            )
            failed_blocks += 1
            continue

        mol_counts += 1

    logging.info(
        "[pose-bust] action=convert status=ok stage=%s kind=dock6 input=%s molecules=%d failed_blocks=%d",
        task.stage,
        task.input_path,
        mol_counts,
        failed_blocks,
    )
    if failed_blocks == 0:
        _cleanup_stage_mol2(stage_dir)
    if mol_counts == 0:
        return "skipped"
    return "converted"


def _df_missing_receptor(
    sdfs: List[Path],
    pdb_id: str,
    variant: str,
    ph: str,
    stage_dir: str,
) -> pd.DataFrame:
    engine = ENGINE_BY_STAGE_DIR.get(stage_dir, "")
    stage = STAGE_NORMALIZED.get(stage_dir, stage_dir)
    rows = []
    for sdf in sdfs:
        rows.append(
            {
                "pdb_id": pdb_id,
                "variant": variant,
                "ph": ph,
                "engine": engine,
                "stage": stage,
                "stage_dir": stage_dir,
                "ligand_file": sdf.name,
                "posebusters_pass": False,
                "posebusters_reason": "missing_receptor",
                "mol_pred_loaded": False,
                "sanitization": False,
                "passes_valence_checks": False,
                "internal_steric_clash": False,
                "bond_lengths": False,
                "bond_angles": False,
                "most_extreme_clash_protein": False,
                "most_extreme_relative_distance_protein": float("nan"),
            }
        )
    return pd.DataFrame(rows)


def _run_bust_for_stage(
    stage_dir_path: Path,
    sdfs: List[Path],
    receptor: Path,
    stage_dir: str,
    pdb_id: str,
    variant: str,
    ph: str,
    cutoff: float,
    max_workers: int,
    chunk_size: int,
) -> Tuple[pd.DataFrame, bool]:
    engine = ENGINE_BY_STAGE_DIR.get(stage_dir, "")
    stage = STAGE_NORMALIZED.get(stage_dir, stage_dir)
    chunks = _chunk_list(sdfs, chunk_size)
    logging.info(
        "[pose-bust] action=posebusters status=info stage_dir=%s sdfs=%d chunk_size=%d chunks=%d workers=%d",
        stage_dir_path,
        len(sdfs),
        chunk_size,
        len(chunks),
        max_workers if max_workers and max_workers > 0 else 0,
    )
    chunk_frames: List[pd.DataFrame] = []
    for idx, chunk in enumerate(chunks):
        tmp_fd, tmp_path_str = tempfile.mkstemp(
            suffix=f".posebusters_chunk{idx}.csv", dir=stage_dir_path
        )
        os.close(tmp_fd)
        tmp_path = Path(tmp_path_str)
        cmd = [
            "bust",
            *[str(p) for p in chunk],
            "-p",
            str(receptor),
            "--outfmt",
            "csv",
            "--full-report",
            "--output",
            str(tmp_path),
        ]
        if max_workers and max_workers > 0:
            cmd.extend(["--max-workers", str(max_workers)])
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            logging.error(
                "[pose-bust] action=posebusters status=failed stage_dir=%s receptor=%s returncode=%s stderr=%s stdout=%s",
                stage_dir_path,
                receptor,
                proc.returncode,
                proc.stderr.strip(),
                proc.stdout.strip(),
            )
            try:
                tmp_path.unlink(missing_ok=True)
            except Exception:
                pass
            return pd.DataFrame(), False

        try:
            df_chunk = pd.read_csv(tmp_path)
        except Exception as exc:
            logging.error(
                "[pose-bust] action=posebusters status=failed reason=read_error stage_dir=%s error=%s",
                stage_dir_path,
                exc,
            )
            try:
                tmp_path.unlink(missing_ok=True)
            except Exception:
                pass
            return pd.DataFrame(), False

        chunk_frames.append(df_chunk)
        try:
            tmp_path.unlink(missing_ok=True)
        except Exception:
            pass

    if not chunk_frames:
        return pd.DataFrame(), False

    df_stage = pd.concat(chunk_frames, ignore_index=True)

    # Robust ligand_file assignment to tolerate PoseBusters row/SDF mismatches.
    if "mol_pred" in df_stage.columns:
        df_stage["ligand_file"] = (
            df_stage["mol_pred"].astype(str).apply(lambda x: Path(x).name)
        )
    else:
        sdf_names = [Path(p).name for p in sdfs]
        n_rows = len(df_stage)
        n_sdfs = len(sdf_names)

        if n_rows == n_sdfs:
            df_stage["ligand_file"] = sdf_names
        elif n_sdfs == 0:
            logging.warning(
                "[pose-bust] action=posebusters status=warning reason=no_sdfs_but_rows stage_dir=%s rows=%d",
                stage_dir_path,
                n_rows,
            )
            df_stage["ligand_file"] = [f"row_{i:05d}" for i in range(n_rows)]
        else:
            repeats = math.ceil(n_rows / n_sdfs)
            repeated = (sdf_names * repeats)[:n_rows]
            logging.warning(
                "[pose-bust] action=posebusters status=warning reason=row_sdf_mismatch stage_dir=%s rows=%d sdfs=%d",
                stage_dir_path,
                n_rows,
                n_sdfs,
            )
            df_stage["ligand_file"] = repeated

    df_stage["pdb_id"] = pdb_id
    df_stage["variant"] = variant
    df_stage["ph"] = ph
    df_stage["engine"] = engine
    df_stage["stage"] = stage
    df_stage["stage_dir"] = stage_dir

    df_stage["posebusters_pass"] = df_stage.apply(
        lambda r: _compute_posebusters_pass(r, cutoff), axis=1
    )

    dup_mask = df_stage.duplicated(subset=["ligand_file"], keep="first")
    if dup_mask.any():
        logging.warning(
            "[pose-bust] action=posebusters status=warning reason=duplicate_ligand_entries stage_dir=%s duplicates=%d",
            stage_dir_path,
            int(dup_mask.sum()),
        )
        df_stage = df_stage.loc[~dup_mask].copy()

    out_csv = stage_dir_path / f"posebusters_{stage_dir}.csv"
    try:
        df_stage.to_csv(out_csv, index=False)
    except Exception as exc:
        logging.error(
            "[pose-bust] action=posebusters status=failed reason=write_error stage_dir=%s output=%s error=%s",
            stage_dir_path,
            out_csv,
            exc,
        )
        return pd.DataFrame(), False

    return df_stage, True


def run_posebusters_for_run(
    post_root: Path,
    repo_root: Path,
    run_id: str,
    cutoff: float,
    max_workers: int,
    chunk_size: int,
) -> Dict[str, object]:
    stage_root = post_root / run_id
    if not stage_root.exists():
        logging.error(
            "[pose-bust] action=posebusters status=failed reason=missing_post_root path=%s",
            stage_root,
        )
        return {"ok": False}

    if shutil.which("bust") is None:
        logging.error(
            "[pose-bust] action=posebusters status=failed reason=missing_bust_binary"
        )
        return {"ok": False}

    stage_frames: List[pd.DataFrame] = []
    total_ligands = 0
    stages_with_receptor = 0
    stages_missing_receptor = 0
    bust_failures = 0

    for stage_dir in STAGE_DIRS:
        for stage_dir_path in stage_root.rglob(stage_dir):
            if not stage_dir_path.is_dir():
                continue
            rel_parts = stage_dir_path.relative_to(stage_root).parts
            if len(rel_parts) < 4:
                logging.warning(
                    "[pose-bust] action=posebusters status=skip reason=unexpected_layout path=%s",
                    stage_dir_path,
                )
                continue
            pdb_id, variant, ph, stage_dir_name = (
                rel_parts[0],
                rel_parts[1],
                rel_parts[2],
                rel_parts[3],
            )
            sdfs = sorted(stage_dir_path.glob("*.sdf"))
            if not sdfs:
                continue

            total_ligands += len(sdfs)
            receptor = resolve_receptor(pdb_id, variant, ph, repo_root)
            if not receptor.exists():
                logging.error(
                    "[pose-bust] action=posebusters status=failed reason=missing_receptor pdb_id=%s variant=%s ph=%s receptor=%s",
                    pdb_id,
                    variant,
                    ph,
                    receptor,
                )
                df_missing = _df_missing_receptor(
                    sdfs, pdb_id, variant, ph, stage_dir_name
                )
                try:
                    df_missing.to_csv(
                        stage_dir_path / f"posebusters_{stage_dir_name}.csv",
                        index=False,
                    )
                except Exception as exc:
                    logging.error(
                        "[pose-bust] action=posebusters status=failed reason=write_error stage_dir=%s error=%s",
                        stage_dir_path,
                        exc,
                    )
                stage_frames.append(df_missing)
                stages_missing_receptor += 1
                continue

            df_stage, ok = _run_bust_for_stage(
                stage_dir_path,
                sdfs,
                receptor,
                stage_dir_name,
                pdb_id,
                variant,
                ph,
                cutoff,
                max_workers,
                chunk_size,
            )
            if not ok:
                bust_failures += 1
            else:
                stages_with_receptor += 1
            if not df_stage.empty:
                stage_frames.append(df_stage)

    if total_ligands == 0:
        logging.info(
            "[pose-bust] action=posebusters status=skip reason=no_ligands run_id=%s",
            run_id,
        )
        return {"ok": True}

    if not stage_frames:
        logging.error(
            "[pose-bust] action=posebusters status=failed reason=no_results run_id=%s",
            run_id,
        )
        return {"ok": False}

    all_df = pd.concat(stage_frames, ignore_index=True)

    dup_mask = all_df.duplicated(
        subset=["pdb_id", "variant", "ph", "stage_dir", "ligand_file"], keep="first"
    )
    if dup_mask.any():
        logging.warning(
            "[pose-bust] action=posebusters status=warning reason=duplicate_ligands overall_duplicates=%d",
            int(dup_mask.sum()),
        )

    # Per (pdb_id, variant, ph) rollups
    group_fields = ["pdb_id", "variant", "ph"]
    for (gp_pdb, gp_var, gp_ph), df_group in all_df.groupby(group_fields, dropna=False):
        combo_dir = stage_root / gp_pdb / gp_var / gp_ph
        combo_dir.mkdir(parents=True, exist_ok=True)
        combo_all_csv = combo_dir / "posebusters_all_stages.csv"
        try:
            df_group.to_csv(combo_all_csv, index=False)
        except Exception as exc:
            logging.error(
                "[pose-bust] action=posebusters status=failed reason=write_error output=%s error=%s",
                combo_all_csv,
                exc,
            )
            return {"ok": False}

        dup_mask_combo = df_group.duplicated(
            subset=["pdb_id", "variant", "ph", "stage_dir", "ligand_file"],
            keep="first",
        )
        if dup_mask_combo.any():
            logging.warning(
                "[pose-bust] action=posebusters status=warning reason=duplicate_ligands group=%s duplicates=%d",
                (gp_pdb, gp_var, gp_ph),
                int(dup_mask_combo.sum()),
            )

        passfail_combo = df_group.loc[
            ~dup_mask_combo,
            [
                "pdb_id",
                "variant",
                "ph",
                "engine",
                "stage",
                "stage_dir",
                "ligand_file",
                "posebusters_pass",
            ],
        ]
        combo_passfail_csv = combo_dir / "posebusters_passfail.csv"
        try:
            passfail_combo.to_csv(combo_passfail_csv, index=False)
        except Exception as exc:
            logging.error(
                "[pose-bust] action=posebusters status=failed reason=write_error output=%s error=%s",
                combo_passfail_csv,
                exc,
            )
            return {"ok": False}

    all_missing_receptor = (
        stages_with_receptor == 0 and stages_missing_receptor > 0 and total_ligands > 0
    )
    ok = bust_failures == 0 and not all_missing_receptor
    if all_missing_receptor:
        logging.error(
            "[pose-bust] action=posebusters status=failed reason=all_missing_receptor run_id=%s",
            run_id,
        )
    if bust_failures > 0:
        logging.error(
            "[pose-bust] action=posebusters status=failed reason=bust_failures count=%d run_id=%s",
            bust_failures,
            run_id,
        )

    return {"ok": ok}


def main() -> int:
    args = parse_args()
    logger = configure_logging()

    try:
        cfg = load_config()
    except Exception as exc:
        logger.error(
            "[pose-bust] action=init status=failed reason=config_load_error error=%s",
            exc,
        )
        return 1

    cpu_cfg = cfg.get("CPU")
    max_jobs_cfg = cfg.get("MAX_PARALLEL_JOBS")

    try:
        cpu_limit = int(cpu_cfg) if cpu_cfg is not None else (os.cpu_count() or 1)
    except Exception:
        cpu_limit = os.cpu_count() or 1

    try:
        max_jobs_limit = int(max_jobs_cfg) if max_jobs_cfg is not None else cpu_limit
    except Exception:
        max_jobs_limit = cpu_limit

    GLOBAL_MAX_WORKERS = max(1, min(cpu_limit, max_jobs_limit, os.cpu_count() or 1))

    if args.workers is None or args.workers <= 0:
        convert_workers = GLOBAL_MAX_WORKERS
    else:
        convert_workers = max(1, min(args.workers, GLOBAL_MAX_WORKERS))

    if args.posebusters_max_workers is None or args.posebusters_max_workers <= 0:
        posebusters_workers = GLOBAL_MAX_WORKERS
    else:
        posebusters_workers = max(
            1, min(args.posebusters_max_workers, GLOBAL_MAX_WORKERS)
        )

    if shutil.which("obabel") is None:
        logger.error("[pose-bust] action=init status=failed reason=missing_obabel")
        return 1

    repo_root = (
        Path(args.repo_root).resolve()
        if args.repo_root
        else Path(__file__).resolve().parent
    )
    docked_root = (
        Path(args.docked_root).resolve() if args.docked_root else repo_root / "docked"
    )
    post_docked_root = (
        Path(args.post_docked_root).resolve()
        if args.post_docked_root
        else repo_root / "post_docked"
    )

    run_root = docked_root / args.run_id
    output_root = post_docked_root / args.run_id

    if not run_root.exists():
        logger.error(
            "[pose-bust] action=init status=failed reason=missing_run_root path=%s",
            run_root,
        )
        return 1

    tasks = discover_inputs(run_root, output_root)
    if not args.include_dock6:
        tasks = [t for t in tasks if t.kind != "dock6"]
    if not tasks:
        logger.error(
            "[pose-bust] action=scan status=failed reason=no_ligands_found run_root=%s",
            run_root,
        )
        return 1

    pdbqt_tasks = [t for t in tasks if t.kind == "pdbqt"]
    dok_tasks = [t for t in tasks if t.kind == "dok"]
    dock6_tasks = [t for t in tasks if t.kind == "dock6"]
    if dok_tasks and shutil.which("ledock") is None:
        logger.error(
            "[pose-bust] action=init status=failed reason=missing_ledock dok_tasks=%d",
            len(dok_tasks),
        )
        return 1

    logger.info(
        "[pose-bust] action=scan status=ok run_id=%s found=%d pdbqt=%d dok=%d dock6=%d docked_root=%s post_docked_root=%s",
        args.run_id,
        len(tasks),
        len(pdbqt_tasks),
        len(dok_tasks),
        len(dock6_tasks),
        docked_root,
        post_docked_root,
    )

    for sample in tasks[:3]:
        logger.info(
            "[pose-bust] action=example stage=%s input=%s output=%s",
            sample.stage,
            sample.input_path,
            sample.output_path,
        )

    statuses: Dict[LigandTask, str] = {}
    if pdbqt_tasks:
        logger.info(
            "[pose-bust] action=convert status=info note=parallel_pdbqt tasks=%d workers=%d",
            len(pdbqt_tasks),
            convert_workers,
        )
        statuses.update(process_tasks(pdbqt_tasks, args.overwrite, convert_workers))
    if dok_tasks:
        logger.info(
            "[pose-bust] action=convert status=info note=parallel_ledock dok_tasks=%d workers=%d",
            len(dok_tasks),
            convert_workers,
        )
        statuses.update(
            process_tasks(dok_tasks, args.overwrite, workers=convert_workers)
        )
    if dock6_tasks:
        logger.info(
            "[pose-bust] action=convert status=info note=parallel_dock6 dock6_tasks=%d workers=%d",
            len(dock6_tasks),
            convert_workers,
        )
        statuses.update(
            process_tasks(dock6_tasks, args.overwrite, workers=convert_workers)
        )

    converted = sum(1 for s in statuses.values() if s == "converted")
    skipped = sum(1 for s in statuses.values() if s == "skipped")
    failed = sum(1 for s in statuses.values() if s == "failed")

    missing = validate_outputs(tasks)
    failed_total = failed + len(missing)

    converted_pdbqt = sum(
        1 for t, s in statuses.items() if t.kind == "pdbqt" and s == "converted"
    )
    converted_dok = sum(
        1 for t, s in statuses.items() if t.kind == "dok" and s == "converted"
    )
    converted_dock6 = sum(
        1 for t, s in statuses.items() if t.kind == "dock6" and s == "converted"
    )
    logger.info(
        "[pose-bust] action=summary run_id=%s found=%d pdbqt=%d dok=%d dock6=%d converted=%d converted_pdbqt=%d converted_dok=%d converted_dock6=%d skipped=%d failed=%d missing=%d",
        args.run_id,
        len(tasks),
        len(pdbqt_tasks),
        len(dok_tasks),
        len(dock6_tasks),
        converted,
        converted_pdbqt,
        converted_dok,
        converted_dock6,
        skipped,
        failed,
        len(missing),
    )

    if failed_total > 0:
        logger.error(
            "[pose-bust] action=complete status=failed run_id=%s failed=%d",
            args.run_id,
            failed_total,
        )
        return 1

    pb_ok = True
    if args.run_posebusters:
        logger.info(
            "[pose-bust] action=posebusters status=info note=run workers=%d chunk=%d",
            posebusters_workers,
            args.posebusters_chunk,
        )
        pb_result = run_posebusters_for_run(
            post_docked_root,
            repo_root,
            args.run_id,
            args.posebusters_relative_cutoff,
            posebusters_workers,
            args.posebusters_chunk,
        )
        pb_ok = bool(pb_result.get("ok", False))
    else:
        logger.info("[pose-bust] action=posebusters status=skip reason=disabled")

    if not pb_ok:
        logger.error(
            "[pose-bust] action=complete status=failed run_id=%s failed=posebusters",
            args.run_id,
        )
        return 1

    logger.info("[pose-bust] action=complete status=ok run_id=%s", args.run_id)
    return 0


if __name__ == "__main__":
    sys.exit(main())
