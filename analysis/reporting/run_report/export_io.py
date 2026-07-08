# ruff: noqa: F403, F405
from analysis.reporting.run_report.constants import *
from analysis.reporting.run_report.paths import _run_data_dir
def write_yaml(report: Dict[str, Any], out_path: Path) -> None:
    # Use a compact representation for lists of numbers (center/box) if possible?
    # PyYAML default dump is okay.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        yaml.dump(report, f, sort_keys=False, default_flow_style=False)

def _format_num(value: Any, fmt: str) -> str:
    try:
        num = float(value)
        if math.isfinite(num):
            return format(num, fmt)
    except (TypeError, ValueError):
        return ""
    return ""

def _stage_artifacts(
    output_dir: Path, artifacts: List[Tuple[str, Path]]
) -> List[Tuple[str, str]]:
    staged: List[Tuple[str, str]] = []
    if not artifacts:
        return staged
    artifacts_dir = output_dir / "artifacts"
    logger = logging.getLogger("run-report")
    for label, src in artifacts:
        if not src.exists():
            continue
        dest = artifacts_dir / src.name
        if dest.exists():
            try:
                if dest.resolve() == src.resolve():
                    rel_path = dest.relative_to(output_dir).as_posix()
                    staged.append((label, rel_path))
                    continue
            except ValueError:
                pass
            stem = src.stem
            suffix = src.suffix
            counter = 1
            while dest.exists():
                dest = artifacts_dir / f"{stem}-{counter}{suffix}"
                counter += 1
        try:
            artifacts_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
        except _REPORT_IO_ERRORS + (shutil.Error,) as exc:
            logger.warning(
                "%s action=stage_artifact_failed src=%s error=%s",
                COMPONENT,
                src,
                exc,
            )
            continue
        rel_path = dest.relative_to(output_dir).as_posix()
        staged.append((label, rel_path))
    return staged

def _resolve_heatmap_source(repo_root: Path, run_id: str) -> Optional[Path]:
    interactions_dataset = _run_data_dir(repo_root, run_id) / "dataset" / "interactions"
    heatmap_csv = _run_data_dir(repo_root, run_id) / "heatmap_input.csv"
    master_csv = _run_data_dir(repo_root, run_id) / "master_rows.csv"
    # Canonical preference: parquet dataset first, then CSV fallbacks.
    if interactions_dataset.exists():
        return interactions_dataset
    if heatmap_csv.exists():
        return heatmap_csv
    if master_csv.exists():
        return master_csv
    return None

