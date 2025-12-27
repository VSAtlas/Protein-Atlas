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
SCORCH_TOP_FRACTION_DEFAULT = 0.5
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


def _load_consensus_top_ligands(
    consensus_csv: Path, frac: float, logger: logging.Logger
) -> Tuple[List[Dict[str, str]], int]:
    rows: List[Dict[str, str]] = []
    try:
        with consensus_csv.open() as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            missing = {"ligand", "consensus_score", "best_engine"} - set(fieldnames)
            if missing:
                logger.error(
                    "%s action=select status=failed reason=missing_fields path=%s missing=%s",
                    COMPONENT,
                    consensus_csv,
                    ",".join(sorted(missing)),
                )
                return [], 0
            for row in reader:
                rows.append(row)
    except Exception as exc:
        logger.error(
            "%s action=select status=failed reason=read_error path=%s error=%s",
            COMPONENT,
            consensus_csv,
            exc,
        )
        return [], 0

    if not rows:
        return [], 0

    def _score(row: Dict[str, str]) -> float:
        try:
            return float(row.get("consensus_score", ""))
        except Exception:
            return float("-inf")

    rows_sorted = sorted(rows, key=_score, reverse=True)
    k = max(1, math.ceil(len(rows_sorted) * frac))
    return rows_sorted[:k], len(rows_sorted)


def _consensus_ligand_base(ligand: str) -> str:
    base = Path(ligand).name
    for ext in (".pdbqt", ".mol2", ".sdf"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    base = base.replace(".sanitized", "")
    return base


def _build_stage_stem_map(ph_root: Path, stage_dir: str) -> Dict[str, Path]:
    stage_map: Dict[str, Path] = {}
    for path in _collect_stage_pdbqts(ph_root, stage_dir):
        if path.stem not in stage_map:
            stage_map[path.stem] = path
    return stage_map


def _pick_vina_pose(base: str, vina_maps: Dict[str, Dict[str, Path]]) -> Optional[Path]:
    for stage_dir in ("stage3", "stage2", "stage1"):
        stage_map = vina_maps.get(stage_dir, {})
        stem = f"{base}_{stage_dir}"
        if stem in stage_map:
            return stage_map[stem]
        if base in stage_map:
            return stage_map[base]
    return None


def _pick_gnina_pose(base: str, gnina_maps: Dict[str, Dict[str, Path]]) -> Optional[Path]:
    for stage_dir in ("gnina_stage3", "gnina_stage2", "gnina_stage1"):
        stage_map = gnina_maps.get(stage_dir, {})
        stem = f"{base}_{stage_dir}"
        if stem in stage_map:
            return stage_map[stem]
        if base in stage_map:
            return stage_map[base]
    return None


def _extract_stage_num_from_stem(stem: str) -> int:
    match = re.search(r"stage(\d+)", stem)
    if match:
        try:
            return int(match.group(1))
        except ValueError:
            return -1
    return -1


def _normalize_prepped_stem_to_base(stem: str) -> str:
    base = stem.replace(".sanitized", "")
    base = re.sub(r"(__ledock_stage\d+|__dock6_stage\d+)$", "", base)
    if base.endswith(".mol2"):
        base = base[: -len(".mol2")]
    base = base.replace("__", "_")
    return base


def _build_prepped_pose_map(prepped_root: Path) -> Dict[str, Path]:
    pose_map: Dict[str, Tuple[int, Path]] = {}
    if not prepped_root.exists():
        return {}
    for path in prepped_root.rglob("*.pdbqt"):
        stem = path.stem
        base = _normalize_prepped_stem_to_base(stem)
        stage_num = _extract_stage_num_from_stem(stem)
        best = pose_map.get(base)
        if best is None or stage_num > best[0]:
            pose_map[base] = (stage_num, path)
    return {base: pair[1] for base, pair in pose_map.items()}


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
    ligands_override: Optional[List[Path]] = None,
    ph_root_for_hash: Optional[Path] = None,
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

    if ligands_override is not None:
        ligands_available = ligands_override
        if not ligands_available:
            logger.warning(
                "%s action=score status=skip source=%s stage=%s reason=no_selected_ligands",
                COMPONENT,
                spec.source,
                spec.stage_dir,
            )
            return True, None
        ph_root = ph_root_for_hash or ligands_available[0].parent
        lig_path = _materialize_inputs(ph_root, combo_post_root, spec.stage_dir, ligands_available, overwrite, logger)
        if lig_path is None:
            logger.error(
                "%s action=score status=failed source=%s stage=%s reason=materialize_failed",
                COMPONENT,
                spec.source,
                spec.stage_dir,
            )
            return False, None
    elif spec.source in {"vina", "gnina"}:
        ph_root = run_root / pdb_id / variant / ph
        ligands_available = _collect_stage_pdbqts(ph_root, spec.stage_dir)
        if not ligands_available:
            logger.warning(
                "%s action=score status=skip source=%s stage=%s reason=no_pdbqt path=%s",
                COMPONENT,
                spec.source,
                spec.stage_dir,
                ph_root / spec.stage_dir,
            )
            return True, None
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
        StageSpec("vina", "vina_selected", "scorch_scores_vina_selected.csv"),
        StageSpec("gnina", "gnina_selected", "scorch_scores_gnina_selected.csv"),
        StageSpec("ledock", "ledock_selected", "scorch_scores_ledock_selected.csv"),
        StageSpec("dock6", "dock6_selected", "scorch_scores_dock6_selected.csv"),
    ]
    spec_by_source = {spec.source: spec for spec in specs}

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

    tasks: List[Tuple[StageSpec, Tuple[str, str, str], Path, List[Path], Path]] = []
    skipped_missing_receptor = 0
    skipped_missing_consensus = 0
    skipped_no_selected = 0
    combos_with_tasks: Set[Tuple[str, str, str]] = set()
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
        if not consensus_csv.exists():
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

        selected_rows, total_rows = _load_consensus_top_ligands(consensus_csv, top_fraction, logger)
        if not selected_rows:
            logger.warning(
                "%s action=select status=skip reason=no_consensus_rows pdb_id=%s variant=%s ph=%s path=%s",
                COMPONENT,
                pdb_id,
                variant,
                ph,
                consensus_csv,
            )
            skipped_no_selected += 1
            continue

        ph_root = run_root / pdb_id / variant / ph
        vina_maps = {stage: _build_stage_stem_map(ph_root, stage) for stage in VINA_STAGE_DIRS}
        gnina_maps = {stage: _build_stage_stem_map(ph_root, stage) for stage in GNINA_STAGE_DIRS}
        ledock_root = post_run_root / pdb_id / variant / ph / "ledock_pdbqt"
        dock6_root = post_run_root / pdb_id / variant / ph / "dock6_pdbqt"
        ledock_map = _build_prepped_pose_map(ledock_root)
        dock6_map = _build_prepped_pose_map(dock6_root)

        selected_by_source: Dict[str, List[Path]] = {"vina": [], "gnina": [], "ledock": [], "dock6": []}
        missing_pose = 0
        for row in selected_rows:
            base = _consensus_ligand_base(row.get("ligand", ""))
            engine = str(row.get("best_engine", "")).strip().lower()
            pose_path: Optional[Path] = None
            if engine == "vina":
                pose_path = _pick_vina_pose(base, vina_maps)
            elif engine == "gnina":
                pose_path = _pick_gnina_pose(base, gnina_maps)
            elif engine == "ledock":
                pose_path = ledock_map.get(base)
            elif engine == "dock6":
                pose_path = dock6_map.get(base)
            if pose_path:
                selected_by_source[engine].append(pose_path)
            else:
                missing_pose += 1
                logger.debug(
                    "%s action=select status=skip reason=missing_pose engine=%s ligand=%s pdb_id=%s variant=%s ph=%s",
                    COMPONENT,
                    engine or "unknown",
                    base,
                    pdb_id,
                    variant,
                    ph,
                )

        found_counts = {k: len(v) for k, v in selected_by_source.items()}
        logger.info(
            "%s action=select status=ok pdb_id=%s variant=%s ph=%s consensus_rows=%d selected=%d frac=%.2f "
            "vina_found=%d gnina_found=%d ledock_found=%d dock6_found=%d missing_pose=%d",
            COMPONENT,
            pdb_id,
            variant,
            ph,
            total_rows,
            len(selected_rows),
            top_fraction,
            found_counts["vina"],
            found_counts["gnina"],
            found_counts["ledock"],
            found_counts["dock6"],
            missing_pose,
        )

        for engine, paths in selected_by_source.items():
            if not paths:
                continue
            spec = spec_by_source.get(engine)
            if not spec:
                continue
            ph_root_for_hash = ph_root if engine in {"vina", "gnina"} else (ledock_root if engine == "ledock" else dock6_root)
            tasks.append((spec, combo, receptor, paths, ph_root_for_hash))
            combos_with_tasks.add(combo)

    total_jobs = len(tasks)
    completed = 0
    failed_jobs = 0

    if args.jobs <= 1:
        for spec, combo, receptor, ligs, ph_root_for_hash in tasks:
            ok, _ = _score_stage(
                spec,
                combo,
                run_root,
                post_run_root,
                receptor,
                args.threads,
                args.overwrite,
                logger,
                ligs,
                ph_root_for_hash,
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
                    ligs,
                    ph_root_for_hash,
                ): (spec, combo)
                for spec, combo, receptor, ligs, ph_root_for_hash in tasks
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
        "%s action=summary status=%s combos=%d combos_scored=%d total_jobs=%d completed=%d failed=%d missing_receptor=%d missing_consensus=%d no_selected=%d",
        COMPONENT,
        "ok" if failed_jobs == 0 else "failed",
        len(combos),
        len(combos_with_tasks),
        total_jobs,
        completed,
        failed_jobs,
        skipped_missing_receptor,
        skipped_missing_consensus,
        skipped_no_selected,
    )
    return 0 if failed_jobs == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
