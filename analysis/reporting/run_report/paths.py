# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
def _run_data_dir(repo_root: Path, run_id: str) -> Path:
    return run_output_dir(repo_root, "data", run_id)

def _run_manifests_dir(repo_root: Path, run_id: str) -> Path:
    return run_output_dir(repo_root, "manifests", run_id)

def _run_post_docked_dir(repo_root: Path, run_id: str) -> Path:
    return run_output_dir(repo_root, "post_docked", run_id)

def _run_docked_dir(repo_root: Path, run_id: str) -> Path:
    return run_output_dir(repo_root, "docked", run_id)

