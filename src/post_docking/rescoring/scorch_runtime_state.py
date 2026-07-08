from __future__ import annotations

import argparse
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from config.output_paths import output_root, run_scoped_root

from post_docking.rescoring.scorch_types import ScorchTask

COMPONENT = "[scorch-rescore]"


def preserve_allowed_chunks(
    chunks: List[Set[str]],
    *,
    tail_mode: bool,
    min_split: int,
) -> List[Set[str]]:
    _ = tail_mode, min_split
    return chunks


def gpu_adaptive_chunk_size(
    *,
    allowed_count: int,
    base_chunk_size: int,
    free_cores: int,
    scheduler_queue_depth: int,
    local_queue_depth: int,
    jobs: int,
    want_threads: int = 1,
) -> int:
    _ = free_cores, scheduler_queue_depth, local_queue_depth, jobs, want_threads
    return max(1, min(max(1, int(base_chunk_size)), max(1, int(allowed_count))))


def preserve_task_for_rechunk(
    task: ScorchTask,
    suffix_seed: str,
) -> List[ScorchTask]:
    _ = suffix_seed
    return [task]


SCORCH_SCRIPT: Path | None = None
SCORCH_ENV: str = "scorch-env"
SCORCH_ENV_PREFIX: Path | None = None
SCORCH_ROOT: Path | None = None
SCORCH_USE_MICROMAMBA: bool = True
SCORCH_PYTHON: Path | None = None
SCORCH_TOP_FRACTION_DEFAULT = 0.10
SCORCH_TOP_FRACTION_KEY = "SCORCH_TOP_FRACTION"
SCORCH_DONE_DIR = "scorch"
SCORCH_DONE_SENTINEL = "_DONE"


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run SCORCH rescoring for a run-id")
    parser.add_argument("--run-id", required=True, help="Run identifier under docked/")
    parser.add_argument(
        "--repo-root",
        default=None,
        help="Repository root (default: directory containing rescoring_scorch.py)",
    )
    parser.add_argument(
        "--docked-root",
        default=None,
        help="Docked root (default: cfg DOCKED_DIR or <repo-root>/docked)",
    )
    parser.add_argument(
        "--post-docked-root",
        default=None,
        help="Post-docked root (default: cfg POST_DOCKED_DIR or <repo-root>/post_docked)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads passed to SCORCH (default: 1)",
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
    parser.add_argument(
        "--decoy-prefix",
        default=None,
        help="Override decoy prefix (disables TEST_MODE_ENABLE-driven multi-prefix)",
    )
    parser.add_argument(
        "--pdb-id",
        default=None,
        help="Restrict rescoring to a single PDB identifier",
    )
    parser.add_argument(
        "--variant",
        default=None,
        help="Restrict rescoring to a specific variant (e.g., APO/HOLO)",
    )
    parser.add_argument(
        "--ph",
        default=None,
        help="Restrict rescoring to a single pH label (e.g., ph_7_0)",
    )
    return parser.parse_args(argv)


def configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("rescoring_scorch")


def _resolve_roots(
    args: argparse.Namespace, cfg: Optional[Dict[str, object]] = None
) -> Tuple[Path, Path, Path, Path]:
    repo_root = (
        Path(args.repo_root).resolve()
        if args.repo_root
        else Path(__file__).resolve().parents[3]
    )
    cfg_map = cfg or {}
    docked_root = (
        Path(args.docked_root).expanduser().resolve()
        if args.docked_root
        else Path(str(cfg_map.get("DOCKED_DIR", output_root(repo_root, "docked"))))
        .expanduser()
        .resolve()
    )
    post_docked_root = (
        Path(args.post_docked_root).expanduser().resolve()
        if args.post_docked_root
        else Path(str(cfg_map.get("POST_DOCKED_DIR", output_root(repo_root, "post_docked"))))
        .expanduser()
        .resolve()
    )
    processed_root = run_scoped_root(
        Path(str(cfg_map.get("OUTPUT_DIR", output_root(repo_root, "processed_pdbs")))),
        getattr(args, "run_id", None),
    ).resolve()
    return repo_root, docked_root, post_docked_root, processed_root


def _default_scorch_script_for_repo(repo_root: Path) -> Optional[str]:
    repo_root = repo_root.expanduser().resolve()
    candidates = [
        repo_root.parent.parent / "tools" / "SCORCH" / "scorch.py",
        repo_root / "tools" / "SCORCH" / "scorch.py",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate.resolve())
    return shutil.which("scorch.py")


def _apply_scorch_config_defaults(
    cfg: Dict[str, object],
    *,
    repo_root: Path,
) -> None:
    if not (cfg.get("SCORCH") or cfg.get("SCORCH_SCRIPT")):
        script_default = _default_scorch_script_for_repo(repo_root)
        if script_default:
            cfg["SCORCH"] = script_default
    if not cfg.get("SCORCH_ENV"):
        cfg["SCORCH_ENV"] = os.environ.get("SCORCH_ENV") or "scorch-env"
    if not cfg.get("SCORCH_ENV_PREFIX"):
        env_prefix = os.environ.get("SCORCH_ENV_PREFIX")
        if env_prefix:
            cfg["SCORCH_ENV_PREFIX"] = env_prefix


def _current_python_matches_scorch_env() -> bool:
    current_python = Path(sys.executable).resolve()
    current_env = current_python.parent.parent if current_python.parent.name == "bin" else None
    if current_env is None:
        return False
    if SCORCH_ENV_PREFIX is not None and current_env == SCORCH_ENV_PREFIX.resolve():
        return True
    return current_env.name == str(SCORCH_ENV)


def _preflight(cfg: Dict[str, object], logger: logging.Logger) -> bool:
    global SCORCH_SCRIPT, SCORCH_ENV, SCORCH_ENV_PREFIX, SCORCH_ROOT, SCORCH_USE_MICROMAMBA, SCORCH_PYTHON

    script_cfg = cfg.get("SCORCH") or cfg.get("SCORCH_SCRIPT")
    env_prefix_cfg = cfg.get("SCORCH_ENV_PREFIX")
    env_cfg = cfg.get("SCORCH_ENV")

    if script_cfg:
        SCORCH_SCRIPT = Path(str(script_cfg))

    if env_cfg:
        SCORCH_ENV = str(env_cfg)
    if env_prefix_cfg and str(env_prefix_cfg).strip():
        SCORCH_ENV_PREFIX = Path(str(env_prefix_cfg).strip()).expanduser().resolve()
    else:
        SCORCH_ENV_PREFIX = None
    SCORCH_PYTHON = None

    SCORCH_USE_MICROMAMBA = shutil.which("micromamba") is not None
    if SCORCH_USE_MICROMAMBA and SCORCH_ENV_PREFIX is not None:
        prefix_python = SCORCH_ENV_PREFIX / "bin" / "python"
        if prefix_python.is_file() and os.access(prefix_python, os.X_OK):
            SCORCH_USE_MICROMAMBA = False
            SCORCH_PYTHON = prefix_python.resolve()
    if SCORCH_USE_MICROMAMBA and _current_python_matches_scorch_env():
        SCORCH_USE_MICROMAMBA = False
        SCORCH_PYTHON = Path(sys.executable).resolve()
    if not SCORCH_USE_MICROMAMBA:
        candidate_env_pythons: list[Path] = []
        if SCORCH_ENV_PREFIX is not None:
            candidate_env_pythons.append(SCORCH_ENV_PREFIX / "bin" / "python")
        mamba_root = os.environ.get("MAMBA_ROOT_PREFIX")
        if mamba_root:
            candidate_env_pythons.append(
                Path(mamba_root).expanduser().resolve()
                / "envs"
                / SCORCH_ENV
                / "bin"
                / "python"
            )
        current_python = Path(sys.executable).resolve()
        if current_python.parent.name == "bin":
            current_env_dir = current_python.parent.parent
            envs_dir = current_env_dir.parent
            if envs_dir.name == "envs":
                candidate_env_pythons.append(envs_dir / SCORCH_ENV / "bin" / "python")
        home = Path.home()
        candidate_env_pythons.extend(
            [
                home / "micromamba" / "envs" / SCORCH_ENV / "bin" / "python",
                home / ".mamba" / "envs" / SCORCH_ENV / "bin" / "python",
            ]
        )
        if SCORCH_PYTHON is None:
            for candidate in candidate_env_pythons:
                if candidate.is_file() and os.access(candidate, os.X_OK):
                    SCORCH_PYTHON = candidate.resolve()
                    break
        if SCORCH_PYTHON is None:
            SCORCH_PYTHON = Path(sys.executable).resolve()
            logger.warning(
                "%s action=preflight status=warn reason=missing_micromamba fallback=current_python",
                COMPONENT,
            )
        else:
            logger.warning(
                "%s action=preflight status=warn reason=missing_micromamba fallback=env_python path=%s",
                COMPONENT,
                SCORCH_PYTHON,
            )

    if not SCORCH_SCRIPT:
        logger.error(
            "%s action=preflight status=failed reason=missing_scorch_cfg",
            COMPONENT,
        )
        return False

    SCORCH_SCRIPT = SCORCH_SCRIPT.resolve()
    if not SCORCH_SCRIPT.exists():
        logger.error(
            "%s action=preflight status=failed reason=missing_scorch path=%s",
            COMPONENT,
            SCORCH_SCRIPT,
        )
        return False

    SCORCH_ROOT = SCORCH_SCRIPT.parent
    if (
        SCORCH_USE_MICROMAMBA
        and SCORCH_ENV_PREFIX is not None
        and not SCORCH_ENV_PREFIX.exists()
    ):
        logger.error(
            "%s action=preflight status=failed reason=missing_scorch_env_prefix path=%s",
            COMPONENT,
            SCORCH_ENV_PREFIX,
        )
        return False

    if not SCORCH_USE_MICROMAMBA:
        env_mode = (
            "current_python"
            if SCORCH_PYTHON is None
            or str(SCORCH_PYTHON) == str(Path(sys.executable).resolve())
            else "env_python"
        )
        env_value = str(SCORCH_PYTHON or Path(sys.executable).resolve())
    else:
        env_mode = "prefix" if SCORCH_ENV_PREFIX is not None else "name"
        env_value = (
            str(SCORCH_ENV_PREFIX) if SCORCH_ENV_PREFIX is not None else str(SCORCH_ENV)
        )
    logger.info(
        "%s action=preflight status=ok scorch_script=%s scorch_env_mode=%s scorch_env=%s scorch_root=%s",
        COMPONENT,
        SCORCH_SCRIPT,
        env_mode,
        env_value,
        SCORCH_ROOT,
    )
    return True
