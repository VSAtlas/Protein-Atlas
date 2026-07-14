"""Runtime output path helpers with legacy-root fallback support."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping

_RUNTIME_NAMES = {
    "configs": "configs",
    "data": "data",
    "docked": "docked",
    "logs": "logs",
    "manifests": "manifests",
    "post_docked": "post_docked",
    "processed_pdbs": "processed_pdbs",
}


def repo_root_from_cfg(cfg: Mapping[str, Any]) -> Path:
    return Path(str(cfg.get("OVERALL_DIR") or ".")).expanduser()


def outputs_root(repo_root: Path) -> Path:
    return Path(repo_root).expanduser() / "outputs"


def output_root(repo_root: Path, name: str) -> Path:
    return outputs_root(repo_root) / _RUNTIME_NAMES.get(name, name)


def legacy_root(repo_root: Path, name: str) -> Path:
    return Path(repo_root).expanduser() / _RUNTIME_NAMES.get(name, name)


def preferred_root(repo_root: Path, name: str, *, prefer_existing: bool = False) -> Path:
    preferred = output_root(repo_root, name)
    legacy = legacy_root(repo_root, name)
    if prefer_existing and not preferred.exists() and legacy.exists():
        return legacy
    return preferred


def runtime_root(
    cfg: Mapping[str, Any],
    key: str,
    name: str,
    *,
    prefer_existing: bool = False,
) -> Path:
    explicit = cfg.get(key)
    if explicit:
        return Path(str(explicit)).expanduser()
    return preferred_root(repo_root_from_cfg(cfg), name, prefer_existing=prefer_existing)


def run_dir_candidates(repo_root: Path, name: str, run_id: str) -> list[Path]:
    return [
        output_root(repo_root, name) / str(run_id),
        outputs_root(repo_root)
        / str(run_id)
        / "outputs"
        / _RUNTIME_NAMES.get(name, name)
        / str(run_id),
        legacy_root(repo_root, name) / str(run_id),
    ]


def run_output_dir(repo_root: Path, name: str, run_id: str) -> Path:
    return first_existing(run_dir_candidates(repo_root, name, str(run_id))) or (
        output_root(repo_root, name) / str(run_id)
    )


def run_scoped_root(root: Path, run_id: str | None) -> Path:
    run_token = str(run_id or "").strip()
    base = Path(root).expanduser()
    if not run_token or base.name == run_token:
        return base
    return base / run_token


def first_existing(candidates: Iterable[Path]) -> Path | None:
    for path in candidates:
        if path.exists():
            return path
    return None
