from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Callable, Optional, Set, Tuple


def posebusters_missing(post_root: Path, combos: Set[Tuple[str, str, str]]) -> bool:
    if not post_root.exists():
        return True
    if combos:
        for pdb_id, variant, ph in combos:
            csv_path = post_root / pdb_id / variant / ph / "posebusters_all_stages.csv"
            if not csv_path.exists():
                return True
        return False
    any_pb = any(post_root.rglob("posebusters_all_stages.csv"))
    return not any_pb


def prep_missing(
    run_root: Path,
    post_root: Path,
    combos: Set[Tuple[str, str, str]],
    *,
    decoy_prefix: str,
    decoy_engine_stage_dirs_all: Callable[[str], Tuple[str, ...]],
) -> bool:
    def _dir_empty(path: Path) -> bool:
        return not path.exists() or not any(path.glob("*.pdbqt"))

    ledock_decoy_dirs = decoy_engine_stage_dirs_all("ledock")
    dock6_decoy_dirs = decoy_engine_stage_dirs_all("dock6")
    decoy_ledock_pdbqt = f"{decoy_prefix}_ledock_pdbqt"
    decoy_dock6_pdbqt = f"{decoy_prefix}_dock6_pdbqt"
    decoy_ledock_score = f"{decoy_prefix}_ledock_docking_score_long.csv"

    if combos:
        for pdb_id, variant, ph in combos:
            combo_run = run_root / pdb_id / variant / ph
            if _dir_empty(post_root / pdb_id / variant / ph / "ledock_pdbqt"):
                return True
            if _dir_empty(post_root / pdb_id / variant / ph / "dock6_pdbqt"):
                return True
            if any((combo_run / d).exists() for d in ledock_decoy_dirs):
                if _dir_empty(post_root / pdb_id / variant / ph / decoy_ledock_pdbqt):
                    return True
            if (combo_run / decoy_ledock_score).exists():
                if _dir_empty(post_root / pdb_id / variant / ph / decoy_ledock_pdbqt):
                    return True
            if any((combo_run / d).exists() for d in dock6_decoy_dirs):
                if _dir_empty(post_root / pdb_id / variant / ph / decoy_dock6_pdbqt):
                    return True
        return False
    any_pdbqt = any(post_root.rglob("*_pdbqt/*.pdbqt"))
    return not any_pdbqt


def run_pose_bust(
    run_id: str,
    repo_root: Path,
    docked_root: Path,
    post_docked_root: Path,
    overwrite: bool,
    max_workers: int,
    logger: Any,
    *,
    component: str,
    sys_executable: str,
    run_subprocess: Callable[..., subprocess.CompletedProcess[str]],
) -> bool:
    posebusters_chunk = 96
    cmd = [
        sys_executable,
        str(repo_root / "pose_bust.py"),
        "--run-id",
        run_id,
        "--docked-root",
        str(docked_root),
        "--post-docked-root",
        str(post_docked_root),
        "--posebusters-chunk",
        str(posebusters_chunk),
    ]
    if int(max_workers) > 0:
        cmd.extend(["--posebusters-max-workers", str(int(max_workers))])
    if overwrite:
        cmd.append("--overwrite")
    proc = run_subprocess(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.warning(
            "%s action=autofix status=skip step=pose_bust returncode=%s",
            component,
            proc.returncode,
        )
        return False
    logger.info(
        "%s action=autofix status=ok step=pose_bust run_id=%s", component, run_id
    )
    return True


def launch_pose_bust_async(
    run_id: str,
    repo_root: Path,
    docked_root: Path,
    post_docked_root: Path,
    overwrite: bool,
    max_workers: int,
    logger: Any,
    *,
    component: str,
    sys_executable: str,
    popen_subprocess: Callable[..., subprocess.Popen],
) -> int:
    workers = max(1, min(2, int(max_workers)))
    posebusters_chunk = 96
    cmd = [
        sys_executable,
        str(repo_root / "pose_bust.py"),
        "--run-id",
        run_id,
        "--docked-root",
        str(docked_root),
        "--post-docked-root",
        str(post_docked_root),
        "--posebusters-chunk",
        str(posebusters_chunk),
        "--posebusters-max-workers",
        str(workers),
    ]
    if overwrite:
        cmd.append("--overwrite")

    proc = popen_subprocess(  # noqa: S603
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        stdin=subprocess.DEVNULL,
        close_fds=True,
        start_new_session=True,
        cwd=str(repo_root),
    )
    logger.info(
        "%s action=autofix status=launched step=pose_bust mode=async run_id=%s pid=%d workers=%d",
        component,
        run_id,
        int(proc.pid),
        int(workers),
    )
    return int(proc.pid)


def run_prep_for_scorch(
    run_id: str,
    repo_root: Path,
    docked_root: Path,
    post_docked_root: Path,
    overwrite: bool,
    logger: Any,
    *,
    component: str,
    sys_executable: str,
    run_subprocess: Callable[..., subprocess.CompletedProcess[str]],
    decoy_prefix: Optional[str] = None,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
    max_workers: int = 1,
) -> bool:
    cmd = [
        sys_executable,
        str(repo_root / "src/post_docking/rescoring/prep_for_scorch.py"),
        "--run-id",
        run_id,
        "--docked-root",
        str(docked_root),
        "--post-docked-root",
        str(post_docked_root),
        "--workers",
        str(max(1, int(max_workers))),
    ]
    if decoy_prefix:
        cmd.extend(["--decoy-prefix", str(decoy_prefix)])
    if pdb_id:
        cmd.extend(["--pdb-id", str(pdb_id)])
    if variant:
        cmd.extend(["--variant", str(variant)])
    if ph:
        cmd.extend(["--ph", str(ph)])
    if overwrite:
        cmd.append("--overwrite")
    proc = run_subprocess(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.warning(
            "%s action=autofix status=skip step=prep_for_scorch returncode=%s",
            component,
            proc.returncode,
        )
        return False
    logger.info(
        "%s action=autofix status=ok step=prep_for_scorch run_id=%s", component, run_id
    )
    return True
