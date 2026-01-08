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

# Ensure repository root is in sys.path for root-level imports
REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from input_and_export_functions import load_config

try:
    try:
        from .rescore_reranker import find_consensus_csv, rerank_consensus_with_scorch
    except (ImportError, ValueError, SystemError):
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
DUD_STAGE_DIRS = ("dud_stage1", "dud_stage2", "dud_stage3")
GNINA_DUD_STAGE_DIRS = ("gnina_dud_stage1", "gnina_dud_stage2", "gnina_dud_stage3")
GNINA_DUD_STAGE_DIRS_LEGACY = ("dud_gnina_stage1", "dud_gnina_stage2", "dud_gnina_stage3")
DOCK6_DUD_STAGE_DIRS = ("dock6_dud_stage1", "dock6_dud_stage2", "dock6_dud_stage3")
DOCK6_DUD_STAGE_DIRS_LEGACY = ("dud_dock6_stage1", "dud_dock6_stage2", "dud_dock6_stage3")
LEDOCK_DUD_STAGE_DIRS = ("ledock_dud_stage1", "ledock_dud_stage2", "ledock_dud_stage3")
LEDOCK_DUD_STAGE_DIRS_LEGACY = ("dud_ledock_stage1", "dud_ledock_stage2", "dud_ledock_stage3")
DUD_POST_STAGE_DIRS = ("dud_ledock_pdbqt", "dud_dock6_pdbqt")
VINA_STAGE_DIRS = ("stage1", "stage2", "stage3") + DUD_STAGE_DIRS
GNINA_STAGE_DIRS = ("gnina_stage1", "gnina_stage2", "gnina_stage3") + GNINA_DUD_STAGE_DIRS + GNINA_DUD_STAGE_DIRS_LEGACY
DOCK6_STAGE_DIRS = ("dock6_stage1", "dock6_stage2", "dock6_stage3") + DOCK6_DUD_STAGE_DIRS + DOCK6_DUD_STAGE_DIRS_LEGACY
LEDOCK_STAGE_DIRS = ("ledock_stage1", "ledock_stage2", "ledock_stage3") + LEDOCK_DUD_STAGE_DIRS + LEDOCK_DUD_STAGE_DIRS_LEGACY
POST_STAGE_DIRS = ("ledock_pdbqt", "dock6_pdbqt") + DUD_POST_STAGE_DIRS


@dataclass(frozen=True)
class StageSpec:
    source: str  # vina, gnina, ledock, dock6
    stage_dir: str
    output_name: str


def discover_stage3_roots(variant_root: Path) -> Dict[str, Path]:
    roots: Dict[str, Path] = {}
    fda = variant_root / "stage3"
    dud = variant_root / "dud_stage3"
    if fda.exists():
        roots["fda"] = fda
    if dud.exists():
        roots["dud"] = dud
    logger = logging.getLogger("rescoring_scorch")
    logger.info(
        "[scorch.discover] variant_root=%s fda=%s dud=%s",
        str(variant_root),
        str(fda) if fda.exists() else "None",
        str(dud) if dud.exists() else "None",
    )
    return roots


def _discover_mode_dirs(
    combo: Tuple[str, str, str],
    run_root: Path,
    post_root: Path,
    specs: Sequence[StageSpec],
    logger: logging.Logger,
) -> Dict[str, List[Path]]:
    pdb_id, variant, ph = combo
    mode_dirs: Dict[str, Set[Path]] = {"fda": set(), "dud": set()}
    for spec in specs:
        root = run_root / pdb_id / variant / ph if spec.source in {"vina", "gnina"} else post_root / pdb_id / variant / ph
        for mode in ("fda", "dud"):
            for stage_dir in stage_dir_candidates(spec.source, mode, root):
                candidate = root / stage_dir
                if candidate.exists():
                    mode_dirs[mode].add(candidate)
    fda_list = sorted(mode_dirs["fda"])
    dud_list = sorted(mode_dirs["dud"])
    logger.info(
        "[scorch.discover] pdb_id=%s variant=%s ph=%s fda=%s dud=%s",
        pdb_id,
        variant,
        ph,
        ",".join(str(p) for p in fda_list) if fda_list else "None",
        ",".join(str(p) for p in dud_list) if dud_list else "None",
    )
    return {"fda": fda_list, "dud": dud_list}


def _prefer_existing(
    root: Optional[Path], primary: Sequence[str], legacy: Sequence[str]
) -> Tuple[str, ...]:
    if root is None:
        return tuple(primary)
    if any((root / name).exists() for name in primary):
        return tuple(primary)
    if any((root / name).exists() for name in legacy):
        return tuple(legacy)
    return tuple(primary)


def stage_dir_candidates(source: str, mode: str, root: Optional[Path] = None) -> Tuple[str, ...]:
    src = (source or "").lower()
    mode_norm = (mode or "").lower()
    if src == "vina":
        return ("stage3", "stage2", "stage1") if mode_norm != "dud" else ("dud_stage3", "dud_stage2", "dud_stage1")
    if src == "gnina":
        if mode_norm != "dud":
            return ("gnina_stage3", "gnina_stage2", "gnina_stage1")
        return _prefer_existing(root, GNINA_DUD_STAGE_DIRS, GNINA_DUD_STAGE_DIRS_LEGACY)
    if src == "dock6":
        if mode_norm != "dud":
            return ("dock6_pdbqt",)
        return _prefer_existing(root, ("dud_dock6_pdbqt",), ("dud_dock6_pdbqt",))
    if src == "ledock":
        if mode_norm != "dud":
            return ("ledock_pdbqt",)
        return _prefer_existing(root, ("dud_ledock_pdbqt",), ("dud_ledock_pdbqt",))
    return ()


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
    repo_root = Path(args.repo_root).resolve() if args.repo_root else Path(__file__).resolve().parents[3]
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
    for stage in VINA_STAGE_DIRS + GNINA_STAGE_DIRS + DOCK6_STAGE_DIRS + LEDOCK_STAGE_DIRS:
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
    stem = re.sub(r"(_dud_gnina_stage\d+)$", "", stem)
    stem = re.sub(r"(_gnina_dud_stage\d+)$", "", stem)
    stem = re.sub(r"(_dock6_dud_stage\d+)$", "", stem)
    stem = re.sub(r"(_dud_dock6_stage\d+)$", "", stem)
    stem = re.sub(r"(_dud_stage\d+)$", "", stem)
    stem = re.sub(r"(__dud_ledock_stage\d+)$", "", stem)
    stem = re.sub(r"(__dud_dock6_stage\d+)$", "", stem)
    stem = re.sub(r"(__ledock_stage\d+)$", "", stem)
    stem = re.sub(r"(__dock6_stage\d+)$", "", stem)
    stem = re.sub(r"(__dock6_dud_stage\d+)$", "", stem)
    stem = re.sub(r"(_gnina_stage\d+)$", "", stem)
    stem = re.sub(r"(_stage\d+)$", "", stem)
    stem = re.sub(r"\.(mol2|pdbqt)$", "", stem, flags=re.IGNORECASE)
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


def _stage_priority(stage_dir: str, pdbqt: Optional[Path] = None) -> int:
    name_lower = stage_dir.lower()
    if name_lower.endswith("stage3") or "stage3" in name_lower:
        return 3
    if name_lower.endswith("stage2") or "stage2" in name_lower:
        return 2
    if name_lower.endswith("stage1") or "stage1" in name_lower:
        return 1
    if pdbqt is not None:
        stem = pdbqt.stem
        match = re.search(r"__dock6_stage(\d+)", stem, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
        match = re.search(r"__ledock_stage(\d+)", stem, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    return 0


def _collect_best_pose_per_base(
    ph_root: Path, stage_dirs: Sequence[str], allowed_bases: Optional[Set[str]], logger: logging.Logger
) -> Tuple[List[Path], Dict[int, int], int, Set[str]]:
    best: Dict[str, Tuple[int, Path]] = {}
    total_candidates = 0
    for stage_dir in stage_dirs:
        stage_pdbqts = _collect_stage_pdbqts(ph_root, stage_dir)
        total_candidates += len(stage_pdbqts)
        for pdbqt in stage_pdbqts:
            base = _pose_base_from_path(pdbqt)
            priority = _stage_priority(stage_dir, pdbqt)
            current = best.get(base)
            if current is None or priority > current[0]:
                best[base] = (priority, pdbqt)

    available_bases = set(best.keys())
    if allowed_bases is not None:
        filtered = {base: entry for base, entry in best.items() if base in allowed_bases}
    else:
        filtered = best

    stage_counts: Dict[int, int] = {}
    for priority, _ in filtered.values():
        stage_counts[priority] = stage_counts.get(priority, 0) + 1

    ligands = [entry[1] for entry in sorted(filtered.values(), key=lambda t: (-t[0], str(t[1])))]
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
    return ligands, stage_counts, total_candidates, available_bases


def _load_consensus_top_bases(
    consensus_csv: Path, frac: float, logger: logging.Logger, control_bases: Optional[Set[str]] = None
) -> Tuple[Set[str], Set[str], int, int, int, int, int]:
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
                return set(), set(), 0, 0, 0, 0, 0
            for row in reader:
                rows.append(row)
    except Exception as exc:
        logger.error(
            "%s action=select status=failed reason=read_error path=%s error=%s",
            COMPONENT,
            consensus_csv,
            exc,
        )
        return set(), set(), 0, 0, 0, 0, 0

    if not rows:
        return set(), set(), 0, 0, 0, 0, 0

    def _score(row: Dict[str, str]) -> float:
        try:
            return float(row.get("consensus_score", ""))
        except Exception:
            return float("-inf")

    control_set = control_bases or set()
    scored_rows: List[Tuple[float, str, bool, str]] = []
    controls_in_consensus: Set[str] = set()
    has_library = "library" in (rows[0].keys() if rows else [])
    unknown_count = 0
    for row in rows:
        lig = str(row.get("ligand", "")).strip()
        if not lig:
            continue
        base = _pose_base_from_path(Path(lig))
        if not base:
            continue
        is_ctrl = base in control_set
        if is_ctrl:
            controls_in_consensus.add(base)
        lib = str(row.get("library", "")).strip().upper() if has_library else ""
        if lib == "":
            lib = "UNKNOWN"
        scored_rows.append((_score(row), base, is_ctrl, lib))
        if lib == "UNKNOWN" and not is_ctrl:
            unknown_count += 1

    if not scored_rows:
        empty = set(control_set)
        return empty, empty, len(rows), 0, 0, len(control_set), len(controls_in_consensus)

    if not has_library:
        non_controls = [(score, base) for score, base, is_ctrl, _ in scored_rows if not is_ctrl]
        non_controls_sorted = sorted(non_controls, key=lambda item: item[0], reverse=True)
        k = max(1, math.ceil(len(non_controls_sorted) * frac)) if non_controls_sorted else 0

        allowed: Set[str] = set(control_set)
        for _, base in non_controls_sorted[:k]:
            allowed.add(base)

        logger.info(
            "[select.top_frac] library=absent total_rows=%d non_controls=%d k=%d controls=%d controls_in_consensus=%d frac=%.3f",
            len(rows),
            len(non_controls_sorted),
            k,
            len(control_set),
            len(controls_in_consensus),
            frac,
        )
        return allowed, allowed, len(rows), k, len(non_controls_sorted), len(control_set), len(controls_in_consensus)

    buckets = {"FDA": [], "DECOY": []}
    non_controls_total = 0
    for score, base, is_ctrl, lib in scored_rows:
        if is_ctrl:
            continue
        non_controls_total += 1
        if lib in buckets:
            buckets[lib].append((score, base))
    allowed_fda: Set[str] = set(control_set)
    allowed_decoy: Set[str] = set(control_set)
    k_fda = 0
    k_decoy = 0
    for lib, bucket in buckets.items():
        sorted_bucket = sorted(bucket, key=lambda item: item[0], reverse=True)
        k = max(1, math.ceil(len(sorted_bucket) * frac)) if sorted_bucket else 0
        if lib == "FDA":
            k_fda = k
        elif lib == "DECOY":
            k_decoy = k
        for _, base in sorted_bucket[:k]:
            if lib == "FDA":
                allowed_fda.add(base)
            elif lib == "DECOY":
                allowed_decoy.add(base)

    logger.info(
        "[select.top_frac] library=stratified total_rows=%d non_controls=%d non_controls_fda=%d k_fda=%d non_controls_decoy=%d k_decoy=%d controls=%d controls_in_consensus=%d unknown=%d frac=%.3f",
        len(rows),
        non_controls_total,
        len(buckets["FDA"]),
        k_fda,
        len(buckets["DECOY"]),
        k_decoy,
        len(control_set),
        len(controls_in_consensus),
        unknown_count,
        frac,
    )
    return (
        allowed_fda,
        allowed_decoy,
        len(rows),
        k_fda + k_decoy,
        non_controls_total,
        len(control_set),
        len(controls_in_consensus),
    )



def _select_top_bases_from_score_csv(
    score_csv: Path,
    frac: float,
    control_bases: Optional[Set[str]],
    *,
    higher_is_better: bool,
    logger: logging.Logger,
    score_cols: Optional[Sequence[str]] = None,
) -> Tuple[Set[str], int, int, int]:
    control_set = control_bases or set()
    if not score_csv.exists():
        logger.info(
            "%s action=select status=skip reason=missing_score_csv path=%s",
            COMPONENT,
            score_csv,
        )
        return set(control_set), 0, 0, len(control_set)

    ligand_keys = (
        "ligand",
        "Ligand",
        "ligand_file",
        "Ligand_file",
        "Ligand_ID",
        "ligand_id",
        "ligand_name",
        "name",
        "molecule",
        "Molecule",
    )
    best_by_base: Dict[str, float] = {}
    noncontrol_seen: Set[str] = set()
    n_rows = 0
    try:
        with score_csv.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            for row in reader:
                n_rows += 1
                lig_val = ""
                for key in ligand_keys:
                    if row.get(key):
                        lig_val = str(row.get(key, "")).strip()
                        break
                if not lig_val:
                    continue
                base = _pose_base_from_path(Path(lig_val))
                if not base:
                    continue
                if base in control_set:
                    continue
                noncontrol_seen.add(base)

                candidates: List[float] = []
                if score_cols:
                    for col in score_cols:
                        raw = row.get(col, "")
                        try:
                            val = float(str(raw).strip())
                        except Exception:
                            continue
                        if math.isfinite(val):
                            candidates.append(val)
                else:
                    for col in fields:
                        if col in ligand_keys or col in {"run_id", "variant", "stage", "valid", "reason"}:
                            continue
                        raw = row.get(col, "")
                        try:
                            val = float(str(raw).strip())
                        except Exception:
                            continue
                        if math.isfinite(val):
                            candidates.append(val)
                if not candidates:
                    continue
                score_val = max(candidates) if higher_is_better else min(candidates)
                current = best_by_base.get(base)
                if current is None:
                    best_by_base[base] = score_val
                else:
                    if higher_is_better and score_val > current:
                        best_by_base[base] = score_val
                    elif (not higher_is_better) and score_val < current:
                        best_by_base[base] = score_val
    except Exception as exc:
        logger.warning(
            "%s action=select status=skip reason=read_error path=%s error=%s",
            COMPONENT,
            score_csv,
            exc,
        )
        return set(control_set), 0, 0, len(control_set)

    pool_bases = sorted(noncontrol_seen)
    k = max(1, math.ceil(len(pool_bases) * frac)) if pool_bases else 0

    ordered_scored = sorted(
        best_by_base.items(),
        key=lambda item: item[1],
        reverse=higher_is_better,
    )
    ordered_bases = [base for base, _ in ordered_scored]
    if len(ordered_bases) < len(pool_bases):
        remaining = [base for base in pool_bases if base not in best_by_base]
        ordered_bases.extend(remaining)

    allowed: Set[str] = set(control_set)
    for base in ordered_bases[:k]:
        allowed.add(base)

    logger.info(
        "%s action=select status=ok path=%s rows=%d n_pool=%d k=%d controls=%d allowed_total=%d scored=%d",
        COMPONENT,
        score_csv,
        n_rows,
        len(pool_bases),
        k,
        len(control_set),
        len(allowed),
        len(best_by_base),
    )
    return allowed, len(pool_bases), k, len(control_set)


def _score_csv_for_spec(
    run_root: Path,
    combo: Tuple[str, str, str],
    spec: StageSpec,
    run_mode: str,
) -> Tuple[Optional[Path], List[str], bool]:
    pdb_id, variant, ph = combo
    combo_root = run_root / pdb_id / variant / ph
    dud = run_mode == "dud"
    higher_is_better = False

    if spec.source == "vina":
        summary = combo_root / ("dud_docking_score_summary.csv" if dud else "docking_score_summary.csv")
        long_csv = combo_root / ("dud_docking_score_long.csv" if dud else "docking_score_long.csv")
        if summary.exists():
            cols = [f"{'dud_' if dud else ''}stage3", f"{'dud_' if dud else ''}stage2", f"{'dud_' if dud else ''}stage1"]
            return summary, cols, higher_is_better
        return long_csv, ["score"], higher_is_better

    if spec.source == "gnina":
        summary = combo_root / ("dud_gnina_docking_score_summary.csv" if dud else "gnina_docking_score_summary.csv")
        long_csv = combo_root / ("dud_gnina_docking_score_long.csv" if dud else "gnina_docking_score_long.csv")
        if summary.exists():
            prefix = "gnina_dud_stage" if dud else "gnina_stage"
            cols = [f"{prefix}3", f"{prefix}2", f"{prefix}1"]
            return summary, cols, higher_is_better
        return long_csv, ["gnina_primary_score", "gnina_minimized_affinity_kcal"], higher_is_better

    if spec.source == "ledock":
        summary = combo_root / ("dud_ledock_docking_score_summary.csv" if dud else "ledock_docking_score_summary.csv")
        long_csv = combo_root / ("dud_ledock_docking_score_long.csv" if dud else "ledock_docking_score_long.csv")
        if summary.exists():
            cols = [f"{'dud_' if dud else ''}stage3", f"{'dud_' if dud else ''}stage2", f"{'dud_' if dud else ''}stage1"]
            return summary, cols, higher_is_better
        return long_csv, ["ledock_best_score_kcal"], higher_is_better

    if spec.source == "dock6":
        summary = combo_root / ("dud_dock6_docking_score_summary.csv" if dud else "dock6_docking_score_summary.csv")
        long_csv = combo_root / ("dud_dock6_docking_score_long.csv" if dud else "dock6_docking_score_long.csv")
        if summary.exists():
            cols = [f"{'dud_' if dud else ''}stage3", f"{'dud_' if dud else ''}stage2", f"{'dud_' if dud else ''}stage1"]
            return summary, cols, higher_is_better
        return long_csv, ["dock6_grid_score"], higher_is_better

    return None, [], higher_is_better



def _materialize_inputs(
    ph_root: Path,
    combo_post_root: Path,
    stage_dir: str,
    ligands: List[Path],
    overwrite: bool,
    logger: logging.Logger,
    *,
    run_mode: str,
) -> Optional[Path]:
    if not ligands:
        return None

    input_dir = combo_post_root / ".scorch_inputs" / run_mode / stage_dir
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


def _log_score_csv_coverage(
    source: str, combo: Tuple[str, str, str], score_csv: Path, allowed_bases: Set[str], logger: logging.Logger
) -> None:
    if not score_csv.exists():
        logger.debug(
            "%s action=select source=%s pdb_id=%s variant=%s ph=%s status=skip reason=missing_score_csv path=%s",
            COMPONENT,
            source,
            combo[0],
            combo[1],
            combo[2],
            score_csv,
        )
        return
    try:
        with score_csv.open() as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "ligand" not in reader.fieldnames:
                logger.debug(
                    "%s action=select source=%s pdb_id=%s variant=%s ph=%s status=skip reason=missing_ligand_column path=%s",
                    COMPONENT,
                    source,
                    combo[0],
                    combo[1],
                    combo[2],
                    score_csv,
                )
                return
            bases_in_csv = {
                base
                for row in reader
                for base in [_pose_base_from_path(Path(str(row.get("ligand", "")).strip()))]
                if row.get("ligand") and base
            }
    except Exception as exc:
        logger.debug(
            "%s action=select source=%s pdb_id=%s variant=%s ph=%s status=skip reason=read_error path=%s error=%s",
            COMPONENT,
            source,
            combo[0],
            combo[1],
            combo[2],
            score_csv,
            exc,
        )
        return

    covered = len(allowed_bases & bases_in_csv)
    logger.info(
        "%s action=select source=%s pdb_id=%s variant=%s ph=%s allowed=%d score_rows=%d covered_by_scores=%d path=%s | score CSV coverage only",
        COMPONENT,
        source,
        combo[0],
        combo[1],
        combo[2],
        len(allowed_bases),
        len(bases_in_csv),
        covered,
        score_csv,
    )


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


def _prep_missing(run_root: Path, post_root: Path, combos: Set[Tuple[str, str, str]]) -> bool:
    def _dir_empty(path: Path) -> bool:
        return not path.exists() or not any(path.glob("*.pdbqt"))

    if combos:
        for pdb_id, variant, ph in combos:
            combo_run = run_root / pdb_id / variant / ph
            if _dir_empty(post_root / pdb_id / variant / ph / "ledock_pdbqt"):
                return True
            if _dir_empty(post_root / pdb_id / variant / ph / "dock6_pdbqt"):
                return True
            if any((combo_run / d).exists() for d in LEDOCK_DUD_STAGE_DIRS + LEDOCK_DUD_STAGE_DIRS_LEGACY):
                if _dir_empty(post_root / pdb_id / variant / ph / "dud_ledock_pdbqt"):
                    return True
            if (combo_run / "dud_ledock_docking_score_long.csv").exists():
                if _dir_empty(post_root / pdb_id / variant / ph / "dud_ledock_pdbqt"):
                    return True
            if any((combo_run / d).exists() for d in DOCK6_DUD_STAGE_DIRS + DOCK6_DUD_STAGE_DIRS_LEGACY):
                if _dir_empty(post_root / pdb_id / variant / ph / "dud_dock6_pdbqt"):
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
    cmd = [sys.executable, str(repo_root / "src/post_docking/rescoring/prep_for_scorch.py"), "--run-id", run_id]
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
    control_bases: Optional[Set[str]] = None,
    run_mode: str = "fda",
    stage_dirs_override: Optional[Sequence[str]] = None,
    score_csv: Optional[Path] = None,
    ) -> Tuple[bool, Optional[Path]]:
    pdb_id, variant, ph = combo
    combo_post_root = post_root / pdb_id / variant / ph
    output_name = spec.output_name if run_mode == "fda" else f"dud_{spec.output_name}"
    out_path = combo_post_root / output_name
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
        stage_dirs: Sequence[str] = stage_dirs_override or (spec.stage_dir,)
        use_best = spec.stage_dir in {"vina_best", "gnina_best"}
        if use_best:
            if stage_dirs_override is None and spec.stage_dir == "vina_best":
                stage_dirs = ("stage3", "stage2", "stage1")
            elif stage_dirs_override is None and spec.stage_dir == "gnina_best":
                stage_dirs = ("gnina_stage3", "gnina_stage2", "gnina_stage1")
            ligands_available, stage_counts, total_candidates, available_bases = _collect_best_pose_per_base(
                ph_root, stage_dirs, allowed_bases, logger
            )
        else:
            ligands_available = []
            total_candidates = 0
            available_bases = set()
            for sd in stage_dirs:
                stage_pdbqts = _collect_stage_pdbqts(ph_root, sd)
                total_candidates += len(stage_pdbqts)
                ligands_available.extend(stage_pdbqts)
                available_bases.update(_pose_base_from_path(p) for p in stage_pdbqts)
            if allowed_bases is not None:
                ligands_available = [p for p in ligands_available if _pose_base_from_path(p) in allowed_bases]
            stage_counts = {}
        if allowed_bases is not None and not ligands_available and total_candidates > 0:
            logger.debug(
                "%s action=select status=debug source=%s stage=%s reason=filtered_empty candidate_bases=%s allowed_sample=%s",
                COMPONENT,
                spec.source,
                spec.stage_dir,
                sorted(list(available_bases))[:10],
                sorted(list(allowed_bases))[:10],
            )
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
                total_candidates,
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
                total_candidates,
                len(ligands_available),
            )
        allowed_size = len(allowed_bases) if allowed_bases is not None else 0
        missing_bases = sorted((allowed_bases or set()) - available_bases)[:10] if allowed_bases is not None else []
        logger.info(
            "%s action=select source=%s stage=%s allowed=%d candidates_bases=%d rescored=%d missing=%d missing_examples=%s",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            allowed_size,
            len(available_bases),
            len(ligands_available),
            len(missing_bases) if allowed_bases is not None else 0,
            missing_bases,
        )
        if allowed_bases is not None and control_bases is not None:
            noncontrol_allowed = allowed_bases - control_bases
            noncontrol_available = available_bases - control_bases
            if not noncontrol_allowed and noncontrol_available:
                logger.debug(
                    "%s action=select status=debug source=%s stage=%s reason=controls_only_selection noncontrol_available=%d examples=%s",
                    COMPONENT,
                    spec.source,
                    spec.stage_dir,
                    len(noncontrol_available),
                    sorted(list(noncontrol_available))[:10],
                )
        lig_path = _materialize_inputs(
            ph_root,
            combo_post_root,
            spec.stage_dir,
            ligands_available,
            overwrite,
            logger,
            run_mode=run_mode,
        )
        if lig_path is None:
            logger.error(
                "%s action=score status=failed source=%s stage=%s reason=materialize_failed",
                COMPONENT,
                spec.source,
                spec.stage_dir,
            )
            return False, None
    else:
        if score_csv is None:
            score_prefix = "dud_" if run_mode == "dud" else ""
            score_csv = run_root / pdb_id / variant / ph / f"{score_prefix}{spec.source}_docking_score_long.csv"
        if allowed_bases is not None and score_csv is not None:
            _log_score_csv_coverage(spec.source, combo, score_csv, allowed_bases, logger)
        ph_root = combo_post_root
        stage_dirs = stage_dirs_override or (spec.stage_dir,)
        ligands_available: List[Path] = []
        available_bases: Set[str] = set()
        total_candidates = 0
        for sd in stage_dirs:
            stage_pdbqts = _collect_stage_pdbqts(ph_root, sd)
            total_candidates += len(stage_pdbqts)
            ligands_available.extend(stage_pdbqts)
            available_bases.update(_pose_base_from_path(p) for p in stage_pdbqts)
        if allowed_bases is not None:
            ligands_available = [p for p in ligands_available if _pose_base_from_path(p) in allowed_bases]
        stage_counts = {}
        if allowed_bases is not None and not ligands_available and total_candidates > 0:
            logger.debug(
                "%s action=select status=debug source=%s stage=%s reason=filtered_empty candidate_bases=%s allowed_sample=%s",
                COMPONENT,
                spec.source,
                spec.stage_dir,
                sorted(list(available_bases))[:10],
                sorted(list(allowed_bases))[:10],
            )
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
        logger.debug(
            "%s action=score source=%s stage=%s candidates=%d ligands_unique=%d stage3=%d stage2=%d stage1=%d stage0=%d",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            total_candidates,
            len(ligands_available),
            stage_counts.get(3, 0),
            stage_counts.get(2, 0),
            stage_counts.get(1, 0),
            stage_counts.get(0, 0),
        )
        allowed_size = len(allowed_bases) if allowed_bases is not None else 0
        missing_bases = sorted((allowed_bases or set()) - available_bases)[:10] if allowed_bases is not None else []
        logger.info(
            "%s action=select source=%s stage=%s allowed=%d candidates_bases=%d rescored=%d missing=%d missing_examples=%s",
            COMPONENT,
            spec.source,
            spec.stage_dir,
            allowed_size,
            len(available_bases),
            len(ligands_available),
            len(missing_bases) if allowed_bases is not None else 0,
            missing_bases,
        )
        if allowed_bases is not None and control_bases is not None:
            noncontrol_allowed = allowed_bases - control_bases
            noncontrol_available = available_bases - control_bases
            if not noncontrol_allowed and noncontrol_available:
                logger.debug(
                    "%s action=select status=debug source=%s stage=%s reason=controls_only_selection noncontrol_available=%d examples=%s",
                    COMPONENT,
                    spec.source,
                    spec.stage_dir,
                    len(noncontrol_available),
                    sorted(list(noncontrol_available))[:10],
                )
        lig_path = _materialize_inputs(
            ph_root,
            combo_post_root,
            spec.stage_dir,
            ligands_available,
            overwrite,
            logger,
            run_mode=run_mode,
        )
        if lig_path is None:
            logger.error(
                "%s action=score status=failed source=%s stage=%s reason=materialize_failed",
                COMPONENT,
                spec.source,
                spec.stage_dir,
            )
            return False, None

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
        "run_mode": run_mode,
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
    post_root: Path,
    specs: List[StageSpec],
    combo: Tuple[str, str, str],
    logger: logging.Logger,
    *,
    run_mode: str = "fda",
    output_name: str = "scorch_scores_all.csv",
) -> Optional[Path]:
    pdb_id, variant, ph = combo
    combo_dir = post_root / pdb_id / variant / ph
    rows: List[Dict[str, str]] = []
    fields: List[str] = []

    for spec in specs:
        csv_name = spec.output_name if run_mode == "fda" else f"dud_{spec.output_name}"
        csv_path = combo_dir / csv_name
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

    out_path = combo_dir / output_name
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
        "%s action=aggregate status=ok rows=%d output=%s pdb_id=%s variant=%s ph=%s run_mode=%s",
        COMPONENT,
        len(rows),
        out_path,
        pdb_id,
        variant,
        ph,
        run_mode,
    )
    return out_path


def _compute_best_composites(rows: List[Dict[str, str]]) -> Dict[str, float]:
    def _as_float_local(val: Any) -> Optional[float]:
        try:
            if val is None:
                return None
            s = str(val).strip()
            if not s:
                return None
            out = float(s)
            return out
        except Exception:
            return None

    best: Dict[str, float] = {}
    for r in rows:
        lig_id = str(r.get("Ligand_ID", "") or r.get("ligand", "")).strip()
        if not lig_id:
            continue
        base = _pose_base_from_path(Path(lig_id))
        if not base:
            continue
        score_val = _as_float_local(r.get("SCORCH_score") or r.get("SCORCH_score_used"))
        cert_val = _as_float_local(r.get("SCORCH_certainty") or r.get("SCORCH_certainty_used"))
        if score_val is None or not math.isfinite(score_val):
            continue
        comp = score_val if cert_val is None or not math.isfinite(cert_val) else score_val * cert_val
        if not math.isfinite(comp):
            continue
        if base not in best or comp > best[base]:
            best[base] = comp
    return best


def _annotate_scorch_file(
    csv_path: Path,
    best_map: Dict[str, float],
    mu_decoy: Optional[float],
    sigma_decoy: Optional[float],
    n_decoys: int,
) -> None:
    rows: List[Dict[str, str]] = []
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]

    extra_fields = [
        "scorch_best_composite",
        "scorch_mu_decoy",
        "scorch_sigma_decoy",
        "scorch_n_decoys",
        "t_vs_decoys_scorch",
    ]
    for ef in extra_fields:
        if ef not in fields:
            fields.append(ef)

    for r in rows:
        lig_id = str(r.get("Ligand_ID", "") or r.get("ligand", "")).strip()
        base = _pose_base_from_path(Path(lig_id)) if lig_id else ""
        comp = best_map.get(base) if base else None
        r["scorch_best_composite"] = "" if comp is None else f"{comp:.6g}"
        r["scorch_mu_decoy"] = "" if mu_decoy is None else f"{mu_decoy:.6g}"
        r["scorch_sigma_decoy"] = "" if sigma_decoy is None else f"{sigma_decoy:.6g}"
        r["scorch_n_decoys"] = "" if n_decoys <= 0 else str(n_decoys)
        if comp is None or mu_decoy is None or sigma_decoy is None or sigma_decoy <= 0:
            r["t_vs_decoys_scorch"] = ""
        else:
            t_val = (comp - mu_decoy) / sigma_decoy
            r["t_vs_decoys_scorch"] = f"{t_val:.6g}"

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def annotate_scorch_t_scores(fda_csv: Path, dud_csv: Path, logger: logging.Logger) -> None:
    if not fda_csv.exists() or not dud_csv.exists():
        logger.info(
            "[t-score.schorch.skip] reason=missing_inputs fda_exists=%s dud_exists=%s",
            fda_csv.exists(),
            dud_csv.exists(),
        )
        return
    def _read(path: Path) -> List[Dict[str, str]]:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            return [dict(r) for r in reader]

    fda_rows = _read(fda_csv)
    dud_rows = _read(dud_csv)
    best_fda = _compute_best_composites(fda_rows)
    best_dud = _compute_best_composites(dud_rows)
    decoy_values = [v for v in best_dud.values() if v is not None and math.isfinite(v)]
    if not decoy_values:
        logger.info(
            "[t-score.schorch.skip] reason=no_decoys n_decoys=0 fda=%s dud=%s",
            fda_csv,
            dud_csv,
        )
        mu = sigma = None
        n_decoys = 0
    else:
        mu = sum(decoy_values) / len(decoy_values)
        variance = sum((v - mu) ** 2 for v in decoy_values) / len(decoy_values)
        sigma = math.sqrt(variance)
        n_decoys = len(decoy_values)
        if sigma <= 0 or not math.isfinite(mu) or not math.isfinite(sigma):
            logger.info(
                "[t-score.schorch.skip] reason=no_decoys_or_sigma0 n_decoys=%d mu=%s sigma=%s",
                n_decoys,
                mu,
                sigma,
            )
            mu = sigma = None
        else:
            logger.info(
                "[t-score.schorch] n_decoys=%d mu=%.6g sigma=%.6g fda=%s dud=%s",
                n_decoys,
                mu,
                sigma,
                fda_csv,
                dud_csv,
            )

    annotate_mu = mu if mu is not None else None
    annotate_sigma = sigma if sigma is not None else None
    annotate_n = n_decoys if mu is not None and sigma is not None else 0
    _annotate_scorch_file(fda_csv, best_fda, annotate_mu, annotate_sigma, annotate_n)
    _annotate_scorch_file(dud_csv, best_dud, annotate_mu, annotate_sigma, annotate_n)


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
        if _prep_missing(run_root, post_run_root, combos):
            _run_prep_for_scorch(args.run_id, repo_root, args.overwrite, logger)

    combos = discover_combos(run_root, post_run_root)
    if not combos:
        logger.warning("%s action=discover status=skip reason=no_combos run_id=%s", COMPONENT, args.run_id)
        return 0

    tasks: List[
        Tuple[
            StageSpec,
            Tuple[str, str, str],
            Path,
            Set[str],
            Set[str],
            str,
            Optional[Sequence[str]],
            Optional[Path],
        ]
    ] = []
    skipped_missing_receptor = 0
    skipped_missing_consensus = 0
    combos_with_tasks: Set[Tuple[str, str, str]] = set()
    combo_modes: Dict[Tuple[str, str, str], Set[str]] = {}
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
        if not consensus_csv.exists() or consensus_csv.stat().st_size == 0:
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
        mode_dirs = _discover_mode_dirs(combo, run_root, post_run_root, specs, logger)
        if mode_dirs["dud"]:
            logger.info(
                "[scorch.run] mode=dud stage_root=%s n_allowed=%d",
                mode_dirs["dud"][0],
                len(control_cache[pdb_id]),
            )
            for spec in specs:
                if spec.source in {"vina", "gnina"}:
                    stage_root = run_root / pdb_id / variant / ph
                else:
                    stage_root = post_run_root / pdb_id / variant / ph
                stage_dirs_override: Sequence[str] = stage_dir_candidates(spec.source, "dud", stage_root)
                score_csv, score_cols, higher_is_better = _score_csv_for_spec(run_root, combo, spec, "dud")
                allowed_bases, n_pool, k, controls_total = _select_top_bases_from_score_csv(
                    score_csv if score_csv is not None else Path(""),
                    top_fraction,
                    control_cache[pdb_id],
                    higher_is_better=higher_is_better,
                    logger=logger,
                    score_cols=score_cols,
                ) if score_csv is not None else (set(control_cache[pdb_id]), 0, 0, len(control_cache[pdb_id]))
                logger.info(
                    "[select.pool] mode=dud engine=%s n_controls=%d n_pool=%d k=%d",
                    spec.source,
                    controls_total,
                    n_pool,
                    k,
                )
                logger.info(
                    "[select.allowed] mode=dud engine=%s score_csv=%s n_candidates=%d n_allowed=%d",
                    spec.source,
                    score_csv if score_csv is not None else "None",
                    n_pool,
                    len(allowed_bases),
                )
                logger.info(
                    "%s action=select source=%s run_mode=dud pdb_id=%s variant=%s ph=%s score_csv=%s n_candidates=%d k=%d controls=%d allowed_total=%d",
                    COMPONENT,
                    spec.source,
                    pdb_id,
                    variant,
                    ph,
                    score_csv if score_csv is not None else "None",
                    n_pool,
                    k,
                    controls_total,
                    len(allowed_bases),
                )
                tasks.append(
                    (spec, combo, receptor, allowed_bases, control_cache[pdb_id], "dud", stage_dirs_override, score_csv)
                )
            combo_modes.setdefault(combo, set()).add("dud")
        if mode_dirs["fda"]:
            logger.info(
                "[scorch.run] mode=fda stage_root=%s n_allowed=%d",
                mode_dirs["fda"][0],
                len(control_cache[pdb_id]),
            )
            for spec in specs:
                score_csv, score_cols, higher_is_better = _score_csv_for_spec(run_root, combo, spec, "fda")
                allowed_bases, n_pool, k, controls_total = _select_top_bases_from_score_csv(
                    score_csv if score_csv is not None else Path(""),
                    top_fraction,
                    control_cache[pdb_id],
                    higher_is_better=higher_is_better,
                    logger=logger,
                    score_cols=score_cols,
                ) if score_csv is not None else (set(control_cache[pdb_id]), 0, 0, len(control_cache[pdb_id]))
                logger.info(
                    "[select.pool] mode=fda engine=%s n_controls=%d n_pool=%d k=%d",
                    spec.source,
                    controls_total,
                    n_pool,
                    k,
                )
                logger.info(
                    "[select.allowed] mode=fda engine=%s score_csv=%s n_candidates=%d n_allowed=%d",
                    spec.source,
                    score_csv if score_csv is not None else "None",
                    n_pool,
                    len(allowed_bases),
                )
                logger.info(
                    "%s action=select source=%s run_mode=fda pdb_id=%s variant=%s ph=%s score_csv=%s n_candidates=%d k=%d controls=%d allowed_total=%d",
                    COMPONENT,
                    spec.source,
                    pdb_id,
                    variant,
                    ph,
                    score_csv if score_csv is not None else "None",
                    n_pool,
                    k,
                    controls_total,
                    len(allowed_bases),
                )
                tasks.append(
                    (spec, combo, receptor, allowed_bases, control_cache[pdb_id], "fda", None, score_csv)
                )
            combo_modes.setdefault(combo, set()).add("fda")
        combos_with_tasks.add(combo)

    total_jobs = len(tasks)
    completed = 0
    failed_jobs = 0
    combos_attempted = len(combos_with_tasks)

    if args.jobs <= 1:
        for spec, combo, receptor, allowed_bases, control_bases, run_mode, stage_dirs_override, score_csv in tasks:
            ok, _ = _score_stage(
                spec,
                combo,
                run_root,
                post_run_root,
                receptor,
                args.threads,
                args.overwrite,
                logger,
                allowed_bases,
                control_bases=control_bases,
                run_mode=run_mode,
                stage_dirs_override=stage_dirs_override,
                score_csv=score_csv,
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
                    control_bases,
                    run_mode,
                    stage_dirs_override,
                    score_csv,
                ): (spec, combo, run_mode)
                for spec, combo, receptor, allowed_bases, control_bases, run_mode, stage_dirs_override, score_csv in tasks
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
                    spec, combo, run_mode = future_map[future]
                    logger.error(
                        "%s action=score status=failed reason=worker_exception source=%s stage=%s combo=%s run_mode=%s error=%s",
                        COMPONENT,
                        spec.source,
                        spec.stage_dir,
                        combo,
                        run_mode,
                        exc,
                    )

    for combo in sorted(combos_with_tasks):
        modes = combo_modes.get(combo, {"fda"})
        fda_all = None
        dud_all = None
        for mode in sorted(modes):
            out_name = "scorch_scores_all.csv" if mode == "fda" else "dud_scorch_scores_all.csv"
            agg_path = _aggregate_combo(post_run_root, specs, combo, logger, run_mode=mode, output_name=out_name)
            if mode == "fda":
                fda_all = agg_path
            elif mode == "dud":
                dud_all = agg_path
        chosen_all = fda_all
        if chosen_all is None:
            # fall back to any mode we aggregated
            for mode in sorted(modes):
                alt = "scorch_scores_all.csv" if mode == "fda" else "dud_scorch_scores_all.csv"
                candidate = (post_run_root / combo[0] / combo[1] / combo[2] / alt)
                if candidate.exists():
                    chosen_all = candidate
                    break
        if chosen_all and rerank_consensus_with_scorch and find_consensus_csv:
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
                    out_csv = chosen_all.parent / "consensus_reranked_scorch.csv"
                    scorch_inputs = [p for p in (fda_all, dud_all) if p is not None]
                    if not scorch_inputs:
                        scorch_inputs = [chosen_all]
                    rerank_consensus_with_scorch(
                        consensus_csv,
                        scorch_inputs,
                        out_csv,
                        logger,
                        overwrite=args.overwrite,
                    )
            except Exception as exc:
                logger.warning(
                    "%s action=rerank status=skip reason=exception combo=%s error=%s",
                    COMPONENT,
                    combo,
                    exc,
                )
        if fda_all and dud_all:
            annotate_scorch_t_scores(fda_all, dud_all, logger)
        elif chosen_all and not rerank_consensus_with_scorch:
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
