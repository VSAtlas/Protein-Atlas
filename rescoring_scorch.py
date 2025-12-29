from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import math
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from input_and_export_functions import load_config

try:
    from rescore_reranker import find_consensus_csv, rerank_consensus_with_scorch
except Exception:  # pragma: no cover - optional dependency
    find_consensus_csv = None
    rerank_consensus_with_scorch = None

COMPONENT = "[scorch-rescore]"
SCORCH_SCRIPT: Path | None = None
SCORCH_ENV: str = "scorch-env"
SCORCH_ROOT: Path | None = None
SCORCH_TOP_FRACTION_DEFAULT = 0.15
SCORCH_TOP_FRACTION_KEY = "SCORCH_TOP_FRACTION"
VINA_STAGE_DIRS = ("stage1", "stage2", "stage3")
GNINA_STAGE_DIRS = ("gnina_stage1", "gnina_stage2", "gnina_stage3")
POST_STAGE_DIRS = ("ledock_pdbqt", "dock6_pdbqt")


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


def _collect_stage_pdbqts(ph_root: Path, stage_dir: str) -> List[Path]:
    candidates = [
        ph_root / stage_dir,
        ph_root / "retries" / stage_dir,
        ph_root / stage_dir / "retries",
        ph_root / f"{stage_dir}_retries",
        ph_root / f"{stage_dir}_retry",
    ]

    for child in ph_root.iterdir() if ph_root.exists() else []:
        if child.is_dir() and child.name.startswith(stage_dir) and "retry" in child.name:
            candidates.append(child)

    collected: Set[Path] = set()
    for candidate in candidates:
        if not candidate.exists():
            continue
        if candidate.is_file() and candidate.suffix == ".pdbqt":
            collected.add(candidate.resolve())
            continue
        if candidate.is_dir():
            for pdbqt in candidate.rglob("*.pdbqt"):
                collected.add(pdbqt.resolve())

    return sorted(collected)


def _pose_base_from_path(p: Path) -> str:
    stem = p.stem
    stem = stem.replace(".sanitized", "")
    stem = re.sub(r"(__ledock_stage\d+)$", "", stem)
    stem = re.sub(r"(__dock6_stage\d+)$", "", stem)
    stem = re.sub(r"(_gnina_stage\d+)$", "", stem)
    stem = re.sub(r"(_stage\d+)$", "", stem)
    stem = stem.replace("__", "_")
    stem = re.sub(r"_+$", "", stem)
    return stem


def _control_base_from_path(p: Path) -> str:
    # Match docking_controls canonicalization for extracted control ligands.
    stem = p.stem.split("_stage")[0]
    return stem.split(".sanitized")[0]


def _load_control_bases(processed_root: Path, pdb_id: str, logger: logging.Logger) -> Set[str]:
    roots = [
        processed_root / pdb_id / "ligands_raw",
        processed_root / f"{pdb_id}_NOLIG" / "ligands_raw",
    ]
    bases: Set[str] = set()
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            base = _control_base_from_path(path)
            if base:
                bases.add(base)
    logger.debug(
        "%s action=controls status=ok pdb_id=%s bases=%d roots=%s",
        COMPONENT,
        pdb_id,
        len(bases),
        ",".join(str(r) for r in roots),
    )
    return bases


def _stage_priority(name: str) -> int:
    name_lower = name.lower()
    if name_lower.endswith("stage3") or "stage3" in name_lower:
        return 3
    if name_lower.endswith("stage2") or "stage2" in name_lower:
        return 2
    if name_lower.endswith("stage1") or "stage1" in name_lower:
        return 1
    return 0


def _collect_best_pose_per_base(
    ph_root: Path, stage_dirs: Sequence[str], allowed_bases: Optional[Set[str]], logger: logging.Logger
) -> Tuple[List[Path], Dict[int, int], int]:
    best: Dict[str, Tuple[int, Path]] = {}
    total_candidates = 0
    for stage_dir in stage_dirs:
        stage_pdbqts = _collect_stage_pdbqts(ph_root, stage_dir)
        total_candidates += len(stage_pdbqts)
        for pdbqt in stage_pdbqts:
            base = _pose_base_from_path(pdbqt)
            if allowed_bases is not None and base not in allowed_bases:
                continue
            priority = _stage_priority(stage_dir)
            current = best.get(base)
            if current is None or priority > current[0]:
                best[base] = (priority, pdbqt)

    stage_counts: Dict[int, int] = {}
    for priority, _ in best.values():
        stage_counts[priority] = stage_counts.get(priority, 0) + 1

    ligands = [entry[1] for entry in sorted(best.values(), key=lambda t: (-t[0], str(t[1])))]
    logger.debug(
        "%s action=collect_best status=ok stage_dirs=%s candidates=%d selected=%d stage3=%d stage2=%d stage1=%d stage0=%d",
        COMPONENT,
        ",".join(stage_dirs),
        total_candidates,
        len(ligands),
        stage_counts.get(3, 0),
        stage_counts.get(2, 0),
        stage_counts.get(1, 0),
        stage_counts.get(0, 0),
    )
    return ligands, stage_counts, total_candidates


def _load_consensus_top_bases(
    consensus_csv: Path, frac: float, logger: logging.Logger, control_bases: Optional[Set[str]] = None
) -> Tuple[Set[str], int, int, int, int]:
    rows: List[Dict[str, str]] = []
    try:
        with consensus_csv.open() as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            missing = {"ligand", "consensus_score"} - set(fieldnames)
            if missing:
                logger.error(
                    "%s action=select status=failed reason=missing_fields path=%s missing=%s",
                    COMPONENT,
                    consensus_csv,
                    ",".join(sorted(missing)),
                )
                return set(), 0, 0, 0, 0
            for row in reader:
                rows.append(row)
    except Exception as exc:
        logger.error(
            "%s action=select status=failed reason=read_error path=%s error=%s",
            COMPONENT,
            consensus_csv,
            exc,
        )
        return set(), 0, 0, 0, 0

    if not rows:
        return set(), 0, 0, 0, 0

    def _score(row: Dict[str, str]) -> float:
        try:
            return float(row.get("consensus_score", ""))
        except Exception:
            return float("-inf")

    control_set = control_bases or set()
    scored_rows: List[Tuple[float, str, bool]] = []
    for row in rows:
        lig = str(row.get("ligand", "")).strip()
        if not lig:
            continue
        base = _pose_base_from_path(Path(lig))
        if not base:
            continue
        scored_rows.append((_score(row), base, base in control_set))

    if not scored_rows:
        return set(), len(rows), 0, 0, 0

    controls_in_consensus = {base for _, base, is_ctrl in scored_rows if is_ctrl}
    non_controls = [(score, base) for score, base, is_ctrl in scored_rows if not is_ctrl]
    non_controls_sorted = sorted(non_controls, key=lambda item: item[0], reverse=True)
    k = max(1, math.ceil(len(non_controls_sorted) * frac)) if non_controls_sorted else 0

    allowed: Set[str] = set(controls_in_consensus)
    for _, base in non_controls_sorted[:k]:
        allowed.add(base)

    return allowed, len(rows), k, len(non_controls_sorted), len(controls_in_consensus)


def _materialize_inputs(
    ph_root: Path, combo_post_root: Path, stage_dir: str, ligands: List[Path], overwrite: bool, logger: logging.Logger
) -> Optional[Path]:
    if not ligands:
        return None

    input_dir = combo_post_root / ".scorch_inputs" / stage_dir
    if input_dir.exists() and not overwrite:
        existing = list(input_dir.glob("*.pdbqt"))
        if existing:
            return input_dir
        shutil.rmtree(input_dir)
    elif input_dir.exists():
        shutil.rmtree(input_dir)

    input_dir.mkdir(parents=True, exist_ok=True)

    existing_names: Set[str] = set()
    for src in ligands:
        name = src.name
        if name in existing_names:
            try:
                rel = src.relative_to(ph_root)
                rel_str = str(rel)
            except ValueError:
                rel_str = str(src)
            digest = hashlib.sha1(rel_str.encode("utf-8")).hexdigest()[:8]
            name = f"{digest}_{name}"
        existing_names.add(name)
        dst = input_dir / name
        try:
            dst.symlink_to(src)
        except OSError:
            try:
                shutil.copy2(src, dst)
            except Exception as exc:
                logger.error(
                    "%s action=link status=failed stage=%s source=%s reason=copy_error error=%s",
                    COMPONENT,
                    stage_dir,
                    src,
                    exc,
                )
                return None

    return input_dir


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
    allowed_bases: Optional[Set[str]] = None,
) -> Tuple[bool, Optional[Path]]:
    pdb_id, variant, ph = combo
    combo_post_root = post_root / pdb_id / variant / ph
    out_path = combo_post_root / spec.output_name
    combo_post_root.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not overwrite:
        logger.info(
            "%s action=score status=skip source=%s stage=%s reason=exists output=%s",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            out_path,
        )
        return True, out_path

    if spec.source in {"vina", "gnina"}:
        ph_root = run_root / pdb_id / variant / ph
        stage_counts: Optional[Dict[int, int]] = None
        stage_dirs: Sequence[str] = (spec.stage_dir,)
        if spec.stage_dir == "vina_best":
            stage_dirs = ("stage3", "stage2", "stage1")
            ligands_available, stage_counts, before_count = _collect_best_pose_per_base(
                ph_root, stage_dirs, allowed_bases, logger
            )
        elif spec.stage_dir == "gnina_best":
            stage_dirs = ("gnina_stage3", "gnina_stage2", "gnina_stage1")
            ligands_available, stage_counts, before_count = _collect_best_pose_per_base(
                ph_root, stage_dirs, allowed_bases, logger
            )
        else:
            ligands_available = _collect_stage_pdbqts(ph_root, spec.stage_dir)
            before_count = len(ligands_available)
            if allowed_bases is not None:
                ligands_available = [p for p in ligands_available if _pose_base_from_path(p) in allowed_bases]
        if not ligands_available:
            logger.warning(
                "%s action=score status=skip source=%s stage=%s reason=no_pdbqt stage_dirs=%s path=%s",
                COMPONENT,
                spec.source,
                spec.stage_dir,
                ",".join(stage_dirs),
                ph_root,
            )
            return True, None
        if stage_counts is not None:
            logger.debug(
                "%s action=score source=%s stage=%s candidates=%d ligands_unique=%d stage3=%d stage2=%d stage1=%d stage0=%d",
                COMPONENT,
                spec.source,
                spec.stage_dir,
                before_count,
                len(ligands_available),
                stage_counts.get(3, 0),
                stage_counts.get(2, 0),
                stage_counts.get(1, 0),
                stage_counts.get(0, 0),
            )
        else:
            logger.debug(
                "%s action=score source=%s stage=%s ligands_before=%d ligands_after=%d",
                COMPONENT,
                spec.source,
                spec.stage_dir,
                before_count,
                len(ligands_available),
            )
        lig_path = _materialize_inputs(ph_root, combo_post_root, spec.stage_dir, ligands_available, overwrite, logger)
        if lig_path is None:
            logger.error(
                "%s action=score status=failed source=%s stage=%s reason=materialize_failed",
                COMPONENT,
                spec.source,
                spec.stage_dir,
            )
            return False, None
    else:
        lig_path = combo_post_root / spec.stage_dir
        if not lig_path.exists():
            logger.warning(
                "%s action=score status=skip source=%s stage=%s reason=missing_ligands path=%s",
                COMPONENT,
                spec.source,
                spec.stage_dir,
                lig_path,
            )
            return True, None
        ligands_available = list(lig_path.rglob("*.pdbqt")) if lig_path.is_dir() else ([lig_path] if lig_path.suffix == ".pdbqt" else [])
        before_count = len(ligands_available)
        if allowed_bases is not None:
            ligands_available = [p for p in ligands_available if _pose_base_from_path(p) in allowed_bases]
        if not ligands_available:
            logger.warning(
                "%s action=score status=skip source=%s stage=%s reason=no_selected_pdbqt path=%s",
                COMPONENT,
                spec.source,
                spec.stage_dir,
                lig_path,
            )
            return True, None
        logger.debug(
            "%s action=score source=%s stage=%s ligands_before=%d ligands_after=%d",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            before_count,
            len(ligands_available),
        )

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


def _aggregate_combo(
    post_root: Path, specs: List[StageSpec], combo: Tuple[str, str, str], logger: logging.Logger
) -> Optional[Path]:
    pdb_id, variant, ph = combo
    combo_dir = post_root / pdb_id / variant / ph
    rows: List[Dict[str, str]] = []
    fields: List[str] = []

    for spec in specs:
        csv_path = combo_dir / spec.output_name
        if not csv_path.exists() or csv_path.stat().st_size == 0:
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
        logger.warning(
            "%s action=aggregate status=skip reason=no_rows pdb_id=%s variant=%s ph=%s",
            COMPONENT,
            pdb_id,
            variant,
            ph,
        )
        return None

    out_path = combo_dir / "scorch_scores_all.csv"
    try:
        with out_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
    except Exception as exc:
        logger.error(
            "%s action=aggregate status=failed path=%s reason=write_error error=%s",
            COMPONENT,
            out_path,
            exc,
        )
        return None

    logger.info(
        "%s action=aggregate status=ok rows=%d output=%s pdb_id=%s variant=%s ph=%s",
        COMPONENT,
        len(rows),
        out_path,
        pdb_id,
        variant,
        ph,
    )
    return out_path


def main() -> int:
    args = parse_args()
    logger = configure_logging(args.verbose)

    try:
        cfg = load_config()
    except Exception as exc:
        logger.error(
            "%s action=preflight status=failed reason=config_load_error error=%s",
            COMPONENT,
            exc,
        )
        return 1

    raw_use = cfg.get("USE_SCORCH", False)
    use_scorch = False
    if isinstance(raw_use, bool):
        use_scorch = raw_use
    else:
        s = str(raw_use).strip().lower()
        use_scorch = s in ("1", "true", "yes", "on")

    top_fraction = SCORCH_TOP_FRACTION_DEFAULT
    if SCORCH_TOP_FRACTION_KEY in cfg:
        try:
            candidate = float(cfg.get(SCORCH_TOP_FRACTION_KEY, SCORCH_TOP_FRACTION_DEFAULT))
            if 0 < candidate <= 1:
                top_fraction = candidate
            else:
                logger.warning(
                    "%s action=preflight status=warn reason=invalid_top_fraction value=%s using_default=%.2f",
                    COMPONENT,
                    candidate,
                    SCORCH_TOP_FRACTION_DEFAULT,
                )
        except Exception as exc:
            logger.warning(
                "%s action=preflight status=warn reason=parse_top_fraction_failed error=%s using_default=%.2f",
                COMPONENT,
                exc,
                SCORCH_TOP_FRACTION_DEFAULT,
            )

    if not use_scorch:
        logger.info(
            "%s action=skip status=ok reason=USE_SCORCH_false run_id=%s",
            COMPONENT,
            args.run_id,
        )
        return 0

    if not _preflight(logger):
        return 1

    repo_root, docked_root, post_root, processed_root = _resolve_roots(args)
    run_root = docked_root / args.run_id
    post_run_root = post_root / args.run_id
    if not run_root.exists():
        logger.error("%s action=preflight status=failed reason=missing_run_root path=%s", COMPONENT, run_root)
        return 1

    specs: List[StageSpec] = [
        StageSpec("vina", "vina_best", "scorch_scores_vina_best.csv"),
        StageSpec("gnina", "gnina_best", "scorch_scores_gnina_best.csv"),
        StageSpec("ledock", "ledock_pdbqt", "scorch_scores_ledock.csv"),
        StageSpec("dock6", "dock6_pdbqt", "scorch_scores_dock6.csv"),
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

    tasks: List[Tuple[StageSpec, Tuple[str, str, str], Path, Set[str]]] = []
    skipped_missing_receptor = 0
    skipped_missing_consensus = 0
    combos_with_tasks: Set[Tuple[str, str, str]] = set()
    control_cache: Dict[str, Set[str]] = {}
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
        consensus_csv = run_root / pdb_id / variant / ph / "consensus_docking_scores.csv"
        if pdb_id not in control_cache:
            control_cache[pdb_id] = _load_control_bases(processed_root, pdb_id, logger)
        allowed_bases, total_rows, selected_rows, noncontrol_rows, control_rows = _load_consensus_top_bases(
            consensus_csv, top_fraction, logger, control_cache[pdb_id]
        )
        if total_rows == 0:
            logger.warning(
                "%s action=select status=skip reason=missing_consensus pdb_id=%s variant=%s ph=%s path=%s",
                COMPONENT,
                pdb_id,
                variant,
                ph,
                consensus_csv,
            )
            skipped_missing_consensus += 1
            continue
        logger.info(
            "%s action=select status=ok pdb_id=%s variant=%s ph=%s consensus_rows=%d noncontrol_rows=%d selected_noncontrol=%d controls_in_consensus=%d frac=%.2f path=%s",
            COMPONENT,
            pdb_id,
            variant,
            ph,
            total_rows,
            noncontrol_rows,
            selected_rows,
            control_rows,
            top_fraction,
            consensus_csv,
        )
        for spec in specs:
            tasks.append((spec, combo, receptor, allowed_bases))
        combos_with_tasks.add(combo)

    total_jobs = len(tasks)
    completed = 0
    failed_jobs = 0
    combos_attempted = len(combos_with_tasks)

    if args.jobs <= 1:
        for spec, combo, receptor, allowed_bases in tasks:
            ok, _ = _score_stage(
                spec, combo, run_root, post_run_root, receptor, args.threads, args.overwrite, logger, allowed_bases
            )
            if ok:
                completed += 1
            else:
                failed_jobs += 1
    else:
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            future_map = {
                pool.submit(
                    _score_stage,
                    spec,
                    combo,
                    run_root,
                    post_run_root,
                    receptor,
                    args.threads,
                    args.overwrite,
                    logger,
                    allowed_bases,
                ): (spec, combo)
                for spec, combo, receptor, allowed_bases in tasks
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

    for combo in sorted(combos_with_tasks):
        all_path = _aggregate_combo(post_run_root, specs, combo, logger)
        if all_path and rerank_consensus_with_scorch and find_consensus_csv:
            try:
                pdb_id, variant, ph = combo
                dock_combo_dir = run_root / pdb_id / variant / ph
                consensus_csv = find_consensus_csv(dock_combo_dir)
                if not consensus_csv:
                    logger.warning(
                        "%s action=rerank status=skip reason=missing_consensus pdb_id=%s variant=%s ph=%s dock_dir=%s",
                        COMPONENT,
                        pdb_id,
                        variant,
                        ph,
                        dock_combo_dir,
                    )
                else:
                    out_csv = all_path.parent / "consensus_reranked_scorch.csv"
                    rerank_consensus_with_scorch(consensus_csv, all_path, out_csv, logger, overwrite=args.overwrite)
            except Exception as exc:
                logger.warning(
                    "%s action=rerank status=skip reason=exception combo=%s error=%s",
                    COMPONENT,
                    combo,
                    exc,
                )
        elif all_path and not rerank_consensus_with_scorch:
            logger.warning("%s action=rerank status=skip reason=reranker_import_failed combo=%s", COMPONENT, combo)

    logger.info(
        "%s action=summary status=%s combos=%d combos_attempted=%d total_jobs=%d completed=%d failed=%d missing_receptor=%d missing_consensus=%d",
        COMPONENT,
        "ok" if failed_jobs == 0 else "failed",
        len(combos),
        combos_attempted,
        total_jobs,
        completed,
        failed_jobs,
        skipped_missing_receptor,
        skipped_missing_consensus,
    )
    return 0 if failed_jobs == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
