from __future__ import annotations

import argparse
import csv
import logging
import math
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

# Ensure repository root is in sys.path for root-level imports
REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from input_and_export_functions import load_config

COMPONENT = "[cnn-rescore]"
CNN_TOP_FRACTION_DEFAULT = 0.10
CNN_TOP_FRACTION_KEY = "CNN_TOP_FRACTION"
CNN_MODEL_DEFAULT = "crossdock_default2018"
CNN_MODEL_KEY = "CNN_MODEL"
GNINA_EXE_KEY = "GNINA_EXE"

DISCOVERY_STAGE_DIRS = (
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
)
# Ignore gnina_stage* during autofix so non-gnina SDF gaps still trigger pose_bust.
AUTOFIX_STAGE_DIRS = (
    "stage1",
    "stage2",
    "stage3",
    "ledock_stage1",
    "ledock_stage2",
    "ledock_stage3",
    "dock6_stage1",
    "dock6_stage2",
    "dock6_stage3",
)


@dataclass(frozen=True)
class StageSpec:
    source: str
    stage_dirs: Tuple[str, ...]


@dataclass(frozen=True)
class ComboContext:
    run_id: str
    pdb_id: str
    variant: str
    ph: str
    receptor: Path
    combo_post_root: Path
    allowed_bases: Set[str]
    logger: logging.Logger


@dataclass(frozen=True)
class RescoreTask:
    context: ComboContext
    source: str
    selected_stage_dir: str
    stage_priority: int
    base: str
    input_sdf: Path
    output_sdf: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run GNINA CNN rescoring for a run-id")
    parser.add_argument("--run-id", required=True, help="Run identifier under docked/")
    parser.add_argument(
        "--repo-root",
        default=None,
        help="Repository root (default: directory containing rescoring_cnn.py)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel CNN rescoring jobs (default: 1)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs",
    )
    parser.add_argument(
        "--skip-autofix",
        action="store_true",
        help="Skip automatic pose_bust repair step",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging for this script",
    )
    return parser.parse_args()


def configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("rescoring_cnn")


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


def _collect_combo_from_rel(parts: Sequence[str]) -> Optional[Tuple[str, str, str]]:
    if len(parts) < 3:
        return None
    return parts[0], parts[1], parts[2]


def discover_combos(run_root: Path, post_root: Path) -> Set[Tuple[str, str, str]]:
    combos: Set[Tuple[str, str, str]] = set()
    for stage in DISCOVERY_STAGE_DIRS:
        for path in run_root.rglob(stage):
            try:
                rel = path.relative_to(run_root).parts
            except ValueError:
                continue
            combo = _collect_combo_from_rel(rel)
            if combo:
                combos.add(combo)

    for stage in DISCOVERY_STAGE_DIRS:
        for path in post_root.rglob(stage):
            try:
                rel = path.relative_to(post_root).parts
            except ValueError:
                continue
            combo = _collect_combo_from_rel(rel)
            if combo:
                combos.add(combo)
    return combos


def _collect_stage_sdfs(ph_root: Path, stage_dir: str) -> List[Path]:
    candidates = [
        ph_root / stage_dir,
        ph_root / "retries" / stage_dir,
        ph_root / stage_dir / "retries",
        ph_root / f"{stage_dir}_retries",
        ph_root / f"{stage_dir}_retry",
    ]

    for child in ph_root.iterdir() if ph_root.exists() else []:
        if (
            child.is_dir()
            and child.name.startswith(stage_dir)
            and "retry" in child.name
        ):
            candidates.append(child)

    collected: Set[Path] = set()
    for candidate in candidates:
        if not candidate.exists():
            continue
        if candidate.is_file() and candidate.suffix.lower() == ".sdf":
            collected.add(candidate.resolve())
            continue
        if candidate.is_dir():
            for sdf in candidate.rglob("*"):
                if sdf.is_file() and sdf.suffix.lower() == ".sdf":
                    collected.add(sdf.resolve())

    return sorted(collected, key=lambda p: str(p))


def _has_any_sdf_in_stage_dirs(ph_root: Path, stage_dirs: Sequence[str]) -> bool:
    if not ph_root.exists():
        return False
    for stage_dir in stage_dirs:
        candidates = [
            ph_root / stage_dir,
            ph_root / "retries" / stage_dir,
            ph_root / stage_dir / "retries",
            ph_root / f"{stage_dir}_retries",
            ph_root / f"{stage_dir}_retry",
        ]
        for child in ph_root.iterdir() if ph_root.exists() else []:
            if (
                child.is_dir()
                and child.name.startswith(stage_dir)
                and "retry" in child.name
            ):
                candidates.append(child)
        for candidate in candidates:
            if not candidate.exists():
                continue
            if candidate.is_file():
                if candidate.suffix.lower() == ".sdf":
                    return True
                continue
            for sdf in candidate.rglob("*.sdf"):
                if sdf.is_file():
                    return True
    return False


def _pose_base_from_path(p: Path) -> str:
    # Keep normalization identical to rescoring_scorch for consensus compatibility.
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


def _load_control_bases(
    processed_root: Path, pdb_id: str, logger: logging.Logger
) -> Set[str]:
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


def _collect_best_sdf_per_base(
    ph_root: Path,
    stage_dirs: Sequence[str],
    allowed_bases: Optional[Set[str]],
    logger: logging.Logger,
) -> Tuple[Dict[str, Tuple[int, str, Path]], Dict[int, int], int]:
    best: Dict[str, Tuple[int, str, Path]] = {}
    total_candidates = 0
    for stage_dir in stage_dirs:
        stage_sdfs = _collect_stage_sdfs(ph_root, stage_dir)
        total_candidates += len(stage_sdfs)
        for sdf in stage_sdfs:
            base = _pose_base_from_path(sdf)
            if allowed_bases is not None and base not in allowed_bases:
                continue
            priority = _stage_priority(stage_dir)
            current = best.get(base)
            if current is None:
                best[base] = (priority, stage_dir, sdf)
                continue
            if priority > current[0]:
                best[base] = (priority, stage_dir, sdf)
                continue
            if priority == current[0] and str(sdf) < str(current[2]):
                best[base] = (priority, stage_dir, sdf)

    stage_counts: Dict[int, int] = {}
    for priority, _, _ in best.values():
        stage_counts[priority] = stage_counts.get(priority, 0) + 1

    logger.debug(
        "%s action=collect_best status=ok stage_dirs=%s candidates=%d selected=%d stage3=%d stage2=%d stage1=%d stage0=%d",
        COMPONENT,
        ",".join(stage_dirs),
        total_candidates,
        len(best),
        stage_counts.get(3, 0),
        stage_counts.get(2, 0),
        stage_counts.get(1, 0),
        stage_counts.get(0, 0),
    )
    return best, stage_counts, total_candidates


def _load_consensus_top_bases(
    consensus_csv: Path,
    frac: float,
    logger: logging.Logger,
    control_bases: Optional[Set[str]] = None,
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
    non_controls = [
        (score, base) for score, base, is_ctrl in scored_rows if not is_ctrl
    ]
    non_controls_sorted = sorted(non_controls, key=lambda item: item[0], reverse=True)
    k = max(1, math.ceil(len(non_controls_sorted) * frac)) if non_controls_sorted else 0

    allowed: Set[str] = set(controls_in_consensus)
    for _, base in non_controls_sorted[:k]:
        allowed.add(base)
    return allowed, len(rows), k, len(non_controls_sorted), len(controls_in_consensus)


def _has_sdf(path: Path) -> bool:
    if not path.exists():
        return False
    if path.is_file():
        return path.suffix.lower() == ".sdf"
    for sdf in path.rglob("*"):
        if sdf.is_file() and sdf.suffix.lower() == ".sdf":
            return True
    return False


def _sdfs_missing(post_root: Path, combos: Set[Tuple[str, str, str]]) -> bool:
    if not post_root.exists():
        return True
    if combos:
        for pdb_id, variant, ph in combos:
            combo_root = post_root / pdb_id / variant / ph
            if not _has_any_sdf_in_stage_dirs(combo_root, AUTOFIX_STAGE_DIRS):
                return True
        return False
    for sdf in post_root.rglob("*.sdf"):
        if any(part.startswith("gnina_stage") for part in sdf.parts):
            continue
        return False
    return True


def _run_pose_bust(
    run_id: str, repo_root: Path, overwrite: bool, logger: logging.Logger
) -> bool:
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
    logger.info(
        "%s action=autofix status=ok step=pose_bust run_id=%s", COMPONENT, run_id
    )
    return True


def _resolve_gnina_exe(
    cfg: Dict[str, object], repo_root: Path, logger: logging.Logger
) -> Optional[Path]:
    exe_cfg = cfg.get(GNINA_EXE_KEY) or cfg.get(GNINA_EXE_KEY.lower())
    if exe_cfg:
        exe_path = Path(str(exe_cfg)).expanduser()
        if exe_path.exists():
            return exe_path
        logger.error(
            "%s action=preflight status=failed reason=missing_gnina_exe_cfg path=%s",
            COMPONENT,
            exe_path,
        )
        return None

    bundled = repo_root / "tools" / "gnina"
    if bundled.exists():
        return bundled

    which = shutil.which("gnina")
    if which:
        return Path(which)

    logger.error(
        "%s action=preflight status=failed reason=missing_gnina_exe", COMPONENT
    )
    return None


def _safe_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(str(value).strip())
    except Exception:
        return None


def _extract_sdf_tags(block: str) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    lines = block.splitlines()
    for idx, line in enumerate(lines):
        s = line.strip()
        if not s.startswith(">") or "<" not in s or ">" not in s:
            continue
        start = s.find("<") + 1
        end = s.find(">", start)
        if end <= start:
            continue
        tag = s[start:end].strip()
        value = ""
        for j in range(idx + 1, len(lines)):
            candidate = lines[j].strip()
            if not candidate:
                continue
            if candidate.startswith(">") and "<" in candidate and ">" in candidate:
                break
            if candidate == "$$$$":
                break
            value = candidate
            break
        tags[tag] = value
    return tags


def _select_best_pose(
    poses: List[Dict[str, Optional[float]]],
) -> Optional[Dict[str, Optional[float]]]:
    if not poses:
        return None
    if any(p.get("cnn_score") is not None for p in poses):
        max_score = max(p["cnn_score"] for p in poses if p.get("cnn_score") is not None)
        top = [p for p in poses if p.get("cnn_score") == max_score]
        if len(top) == 1:
            return top[0]
        with_affinity = [p for p in top if p.get("cnn_affinity") is not None]
        if with_affinity:
            min_affinity = min(
                p["cnn_affinity"]
                for p in with_affinity
                if p.get("cnn_affinity") is not None
            )
            for pose in top:
                if pose.get("cnn_affinity") == min_affinity:
                    return pose
        return top[0]
    return poses[0]


def _parse_sdf_scores(
    path: Path, logger: logging.Logger
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[str]]:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except Exception as exc:
        return None, None, None, f"read_error: {exc}"

    if not text.strip():
        return None, None, None, "empty_output"

    poses: List[Dict[str, Optional[float]]] = []
    for block in text.split("$$$$"):
        if not block.strip():
            continue
        tags = _extract_sdf_tags(block)
        tags_lower = {k.lower(): v for k, v in tags.items()}
        affinity_val = tags.get("Affinity") or tags_lower.get("affinity")
        if affinity_val is None:
            affinity_val = tags.get("minimizedAffinity") or tags_lower.get(
                "minimizedaffinity"
            )
        poses.append(
            {
                "cnn_score": _safe_float(tags_lower.get("cnnscore")),
                "cnn_affinity": _safe_float(tags_lower.get("cnnaffinity")),
                "affinity": _safe_float(affinity_val),
            }
        )

    best = _select_best_pose(poses)
    if best is None:
        logger.warning(
            "%s action=parse status=warning reason=no_poses path=%s", COMPONENT, path
        )
        return None, None, None, "no_poses"
    return best.get("cnn_score"), best.get("cnn_affinity"), best.get("affinity"), None


def _compact_error(text: str, limit: int = 400) -> str:
    if not text:
        return ""
    cleaned = " ".join(text.strip().split())
    return cleaned[:limit]


def _gnina_command(
    gnina_exe: Path, receptor: Path, input_sdf: Path, output_sdf: Path, model: str
) -> List[str]:
    return [
        str(gnina_exe),
        "-r",
        str(receptor),
        "-l",
        str(input_sdf),
        "--score_only",
        "--cnn_scoring",
        "rescore",
        "--cnn",
        model,
        "--autobox_ligand",
        str(input_sdf),
        "-o",
        str(output_sdf),
    ]


def _empty_row(context: ComboContext, base: str) -> Dict[str, str]:
    return {
        "run_id": context.run_id,
        "pdb_id": context.pdb_id,
        "variant": context.variant,
        "ph": context.ph,
        "source": "",
        "selected_stage_dir": "",
        "stage_priority": "",
        "ligand": base,
        "ligand_file": "",
        "input_sdf": "",
        "output_sdf": "",
        "cnn_score": "",
        "cnn_affinity": "",
        "affinity": "",
        "status": "missing_pose",
        "error": "no_pose_for_allowed_base",
    }


def _stage_label(source: str, stage_dir: str) -> str:
    if stage_dir.startswith(("gnina_", "ledock_", "dock6_")):
        return stage_dir
    if stage_dir.startswith("stage"):
        return f"vina_{stage_dir}"
    return stage_dir or source


def _rescore_task(
    task: RescoreTask, gnina_exe: Path, model: str, overwrite: bool
) -> Dict[str, str]:
    context = task.context
    logger = context.logger
    stage_label = _stage_label(task.source, task.selected_stage_dir)
    row = {
        "run_id": context.run_id,
        "pdb_id": context.pdb_id,
        "variant": context.variant,
        "ph": context.ph,
        "source": task.source,
        "selected_stage_dir": task.selected_stage_dir,
        "stage_priority": stage_label,
        "ligand": task.base,
        "ligand_file": task.input_sdf.name,
        "input_sdf": str(task.input_sdf),
        "output_sdf": str(task.output_sdf),
        "cnn_score": "",
        "cnn_affinity": "",
        "affinity": "",
        "status": "",
        "error": "",
    }

    if not task.input_sdf.exists():
        row["status"] = "missing_pose"
        row["error"] = "input_sdf_missing"
        logger.warning(
            "%s action=rescore status=skip reason=input_missing source=%s base=%s input=%s",
            COMPONENT,
            task.source,
            task.base,
            task.input_sdf,
        )
        return row

    task.output_sdf.parent.mkdir(parents=True, exist_ok=True)

    old_output = task.output_sdf.parent / f"{task.input_sdf.name}.cnn.sdf"
    if not task.output_sdf.exists() and old_output.exists():
        try:
            old_output.rename(task.output_sdf)
        except Exception:
            try:
                shutil.copy2(old_output, task.output_sdf)
            except Exception as exc:
                row["status"] = "failed"
                row["error"] = f"rename_output_failed: {exc}"
                logger.error(
                    "%s action=rescore status=failed reason=rename_failed source=%s base=%s output=%s",
                    COMPONENT,
                    task.source,
                    task.base,
                    task.output_sdf,
                )
                return row

    if task.output_sdf.exists() and not overwrite:
        logger.info(
            "%s action=rescore status=skip reason=output_exists source=%s base=%s output=%s",
            COMPONENT,
            task.source,
            task.base,
            task.output_sdf,
        )
    else:
        cmd = _gnina_command(
            gnina_exe, context.receptor, task.input_sdf, task.output_sdf, model
        )
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            row["status"] = "failed"
            row["error"] = _compact_error(
                f"gnina_returncode={proc.returncode} stderr={proc.stderr} stdout={proc.stdout}"
            )
            logger.error(
                "%s action=rescore status=failed reason=gnina_failed source=%s base=%s returncode=%s stderr=%s",
                COMPONENT,
                task.source,
                task.base,
                proc.returncode,
                _compact_error(proc.stderr),
            )
            return row

    if not task.output_sdf.exists() or task.output_sdf.stat().st_size == 0:
        row["status"] = "failed"
        row["error"] = "missing_output"
        logger.error(
            "%s action=rescore status=failed reason=missing_output source=%s base=%s output=%s",
            COMPONENT,
            task.source,
            task.base,
            task.output_sdf,
        )
        return row

    # Parse GNINA output without RDKit to keep dependencies minimal.
    cnn_score, cnn_affinity, affinity, parse_error = _parse_sdf_scores(
        task.output_sdf, logger
    )
    if parse_error:
        row["status"] = "failed"
        row["error"] = parse_error
        logger.error(
            "%s action=rescore status=failed reason=parse_error source=%s base=%s error=%s",
            COMPONENT,
            task.source,
            task.base,
            parse_error,
        )
        return row

    if cnn_score is not None:
        row["cnn_score"] = str(cnn_score)
    if cnn_affinity is not None:
        row["cnn_affinity"] = str(cnn_affinity)
    if affinity is not None:
        row["affinity"] = str(affinity)
    row["status"] = "ok"
    logger.info(
        "%s action=rescore status=ok source=%s base=%s stage=%s output=%s",
        COMPONENT,
        task.source,
        task.base,
        task.selected_stage_dir,
        task.output_sdf,
    )
    return row


def _ensure_combo_logger(
    base_logger: logging.Logger, combo: Tuple[str, str, str], post_run_root: Path
) -> logging.Logger:
    pdb_id, variant, ph = combo
    name = f"rescoring_cnn.{pdb_id}.{variant}.{ph}"
    logger = logging.getLogger(name)
    logger.setLevel(base_logger.level)
    logger.propagate = True
    combo_dir = post_run_root / pdb_id / variant / ph
    combo_dir.mkdir(parents=True, exist_ok=True)
    log_path = combo_dir / "cnn_rescore.log"
    if not any(
        isinstance(handler, logging.FileHandler)
        and getattr(handler, "baseFilename", "") == str(log_path)
        for handler in logger.handlers
    ):
        handler = logging.FileHandler(log_path)
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        logger.addHandler(handler)
    return logger


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

    top_fraction = CNN_TOP_FRACTION_DEFAULT
    if CNN_TOP_FRACTION_KEY in cfg:
        try:
            candidate = float(cfg.get(CNN_TOP_FRACTION_KEY, CNN_TOP_FRACTION_DEFAULT))
            if 0 < candidate <= 1:
                top_fraction = candidate
            else:
                logger.warning(
                    "%s action=preflight status=warn reason=invalid_top_fraction value=%s using_default=%.2f",
                    COMPONENT,
                    candidate,
                    CNN_TOP_FRACTION_DEFAULT,
                )
        except Exception as exc:
            logger.warning(
                "%s action=preflight status=warn reason=parse_top_fraction_failed error=%s using_default=%.2f",
                COMPONENT,
                exc,
                CNN_TOP_FRACTION_DEFAULT,
            )

    cnn_model = (
        str(cfg.get(CNN_MODEL_KEY, CNN_MODEL_DEFAULT)).strip() or CNN_MODEL_DEFAULT
    )

    repo_root, docked_root, post_root, processed_root = _resolve_roots(args)
    run_root = docked_root / args.run_id
    post_run_root = post_root / args.run_id
    if not run_root.exists():
        logger.error(
            "%s action=preflight status=failed reason=missing_run_root path=%s",
            COMPONENT,
            run_root,
        )
        return 1

    gnina_exe = _resolve_gnina_exe(cfg, repo_root, logger)
    if gnina_exe is None:
        return 1

    logger.info(
        "%s action=preflight status=ok run_id=%s gnina_exe=%s cnn_model=%s top_fraction=%.2f cnn_skip_gnina=true",
        COMPONENT,
        args.run_id,
        gnina_exe,
        cnn_model,
        top_fraction,
    )

    stage_specs = [
        StageSpec("vina_best", ("stage3", "stage2", "stage1")),
        StageSpec("ledock_best", ("ledock_stage3", "ledock_stage2", "ledock_stage1")),
        StageSpec("dock6_best", ("dock6_stage3", "dock6_stage2", "dock6_stage1")),
    ]

    combos = discover_combos(run_root, post_run_root)
    logger.info("%s action=discover status=ok combos=%d", COMPONENT, len(combos))

    if not args.skip_autofix and _sdfs_missing(post_run_root, combos):
        _run_pose_bust(args.run_id, repo_root, args.overwrite, logger)

    combos = discover_combos(run_root, post_run_root)
    if not combos:
        logger.warning(
            "%s action=discover status=skip reason=no_combos run_id=%s",
            COMPONENT,
            args.run_id,
        )
        return 0

    combo_contexts: List[ComboContext] = []
    skipped_missing_receptor = 0
    skipped_missing_consensus = 0
    skipped_existing = 0

    control_cache: Dict[str, Set[str]] = {}
    for combo in sorted(combos):
        combo_logger = _ensure_combo_logger(logger, combo, post_run_root)
        pdb_id, variant, ph = combo
        combo_post_root = post_run_root / pdb_id / variant / ph
        output_csv = combo_post_root / "cnn_rescoring.csv"
        if output_csv.exists() and not args.overwrite:
            combo_logger.info(
                "%s action=combo status=skip reason=exists pdb_id=%s variant=%s ph=%s path=%s",
                COMPONENT,
                pdb_id,
                variant,
                ph,
                output_csv,
            )
            skipped_existing += 1
            continue

        receptor = (
            processed_root
            / pdb_id
            / variant
            / "receptor"
            / "ph_ensemble"
            / f"{pdb_id}_{ph}.pdbqt"
        )
        if not receptor.exists():
            combo_logger.error(
                "%s action=combo status=skip reason=missing_receptor pdb_id=%s variant=%s ph=%s receptor=%s",
                COMPONENT,
                pdb_id,
                variant,
                ph,
                receptor,
            )
            skipped_missing_receptor += 1
            continue

        consensus_csv = (
            run_root / pdb_id / variant / ph / "consensus_docking_scores.csv"
        )
        if pdb_id not in control_cache:
            control_cache[pdb_id] = _load_control_bases(
                processed_root, pdb_id, combo_logger
            )
        (
            allowed_bases,
            total_rows,
            selected_rows,
            noncontrol_rows,
            control_rows,
        ) = _load_consensus_top_bases(
            consensus_csv, top_fraction, combo_logger, control_cache[pdb_id]
        )
        if total_rows == 0:
            combo_logger.warning(
                "%s action=select status=skip reason=missing_consensus pdb_id=%s variant=%s ph=%s path=%s",
                COMPONENT,
                pdb_id,
                variant,
                ph,
                consensus_csv,
            )
            skipped_missing_consensus += 1
            continue

        combo_logger.info(
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

        combo_contexts.append(
            ComboContext(
                run_id=args.run_id,
                pdb_id=pdb_id,
                variant=variant,
                ph=ph,
                receptor=receptor,
                combo_post_root=combo_post_root,
                allowed_bases=allowed_bases,
                logger=combo_logger,
            )
        )

    if not combo_contexts:
        logger.warning(
            "%s action=summary status=skip reason=no_valid_combos run_id=%s",
            COMPONENT,
            args.run_id,
        )
        return 0

    tasks: List[RescoreTask] = []
    results_by_combo: Dict[Tuple[str, str, str], List[Dict[str, str]]] = {
        (ctx.pdb_id, ctx.variant, ctx.ph): [] for ctx in combo_contexts
    }
    allowed_bases_by_combo: Dict[Tuple[str, str, str], Set[str]] = {
        (ctx.pdb_id, ctx.variant, ctx.ph): ctx.allowed_bases for ctx in combo_contexts
    }

    for ctx in combo_contexts:
        ph_root = post_run_root / ctx.pdb_id / ctx.variant / ctx.ph
        for spec in stage_specs:
            best_map, stage_counts, total_candidates = _collect_best_sdf_per_base(
                ph_root, spec.stage_dirs, ctx.allowed_bases, ctx.logger
            )
            if not best_map:
                ctx.logger.warning(
                    "%s action=select status=skip reason=no_sdf source=%s pdb_id=%s variant=%s ph=%s stage_dirs=%s",
                    COMPONENT,
                    spec.source,
                    ctx.pdb_id,
                    ctx.variant,
                    ctx.ph,
                    ",".join(spec.stage_dirs),
                )
                continue
            ctx.logger.debug(
                "%s action=select status=ok source=%s pdb_id=%s variant=%s ph=%s candidates=%d selected=%d stage3=%d stage2=%d stage1=%d stage0=%d",
                COMPONENT,
                spec.source,
                ctx.pdb_id,
                ctx.variant,
                ctx.ph,
                total_candidates,
                len(best_map),
                stage_counts.get(3, 0),
                stage_counts.get(2, 0),
                stage_counts.get(1, 0),
                stage_counts.get(0, 0),
            )
            for base, (priority, stage_dir, sdf_path) in best_map.items():
                output_sdf = (
                    ctx.combo_post_root
                    / "cnn_outputs"
                    / spec.source
                    / f"{sdf_path.stem}.cnn.sdf"
                )
                tasks.append(
                    RescoreTask(
                        context=ctx,
                        source=spec.source,
                        selected_stage_dir=stage_dir,
                        stage_priority=priority,
                        base=base,
                        input_sdf=sdf_path,
                        output_sdf=output_sdf,
                    )
                )

    tasks = sorted(
        tasks,
        key=lambda t: (
            t.context.pdb_id,
            t.context.variant,
            t.context.ph,
            t.source,
            t.base,
            str(t.input_sdf),
        ),
    )

    total_jobs = len(tasks)
    completed = 0
    failed_jobs = 0

    if args.jobs <= 1:
        for task in tasks:
            row = _rescore_task(task, gnina_exe, cnn_model, args.overwrite)
            results_by_combo[
                (task.context.pdb_id, task.context.variant, task.context.ph)
            ].append(row)
            if row.get("status") == "failed":
                failed_jobs += 1
            else:
                completed += 1
    else:
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            future_map = {
                pool.submit(
                    _rescore_task, task, gnina_exe, cnn_model, args.overwrite
                ): task
                for task in tasks
            }
            for future in as_completed(future_map):
                task = future_map[future]
                try:
                    row = future.result()
                except Exception as exc:
                    row = _empty_row(task.context, task.base)
                    row["status"] = "failed"
                    row["error"] = _compact_error(str(exc))
                    task.context.logger.error(
                        "%s action=rescore status=failed reason=worker_exception source=%s base=%s error=%s",
                        COMPONENT,
                        task.source,
                        task.base,
                        exc,
                    )
                results_by_combo[
                    (task.context.pdb_id, task.context.variant, task.context.ph)
                ].append(row)
                if row.get("status") == "failed":
                    failed_jobs += 1
                else:
                    completed += 1

    output_fields = [
        "run_id",
        "pdb_id",
        "variant",
        "ph",
        "stage_priority",
        "ligand",
        "cnn_score",
        "cnn_affinity",
        "affinity",
    ]

    for ctx in combo_contexts:
        combo_key = (ctx.pdb_id, ctx.variant, ctx.ph)
        rows = results_by_combo.get(combo_key, [])
        allowed_bases = allowed_bases_by_combo.get(combo_key, set())
        seen_bases = {row.get("ligand") for row in rows if row.get("ligand")}
        missing_bases = sorted(base for base in allowed_bases if base not in seen_bases)
        for base in missing_bases:
            rows.append(_empty_row(ctx, base))

        if not rows:
            ctx.logger.warning(
                "%s action=write status=skip reason=no_rows pdb_id=%s variant=%s ph=%s",
                COMPONENT,
                ctx.pdb_id,
                ctx.variant,
                ctx.ph,
            )
            continue

        rows_sorted = sorted(
            rows,
            key=lambda r: (
                r.get("ligand", ""),
                r.get("source", ""),
                r.get("selected_stage_dir", ""),
                r.get("input_sdf", ""),
            ),
        )

        out_csv = ctx.combo_post_root / "cnn_rescoring.csv"
        ctx.combo_post_root.mkdir(parents=True, exist_ok=True)
        try:
            with out_csv.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=output_fields)
                writer.writeheader()
                for row in rows_sorted:
                    writer.writerow({key: row.get(key, "") for key in output_fields})
        except Exception as exc:
            ctx.logger.error(
                "%s action=write status=failed reason=write_error pdb_id=%s variant=%s ph=%s error=%s",
                COMPONENT,
                ctx.pdb_id,
                ctx.variant,
                ctx.ph,
                exc,
            )
            failed_jobs += 1
            continue

        ctx.logger.info(
            "%s action=write status=ok rows=%d output=%s pdb_id=%s variant=%s ph=%s",
            COMPONENT,
            len(rows_sorted),
            out_csv,
            ctx.pdb_id,
            ctx.variant,
            ctx.ph,
        )

    logger.info(
        "%s action=summary status=%s combos=%d combos_attempted=%d total_jobs=%d completed=%d failed=%d missing_receptor=%d missing_consensus=%d skipped_existing=%d",
        COMPONENT,
        "ok" if failed_jobs == 0 else "failed",
        len(combos),
        len(combo_contexts),
        total_jobs,
        completed,
        failed_jobs,
        skipped_missing_receptor,
        skipped_missing_consensus,
        skipped_existing,
    )
    return 0 if failed_jobs == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
