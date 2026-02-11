from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any, Mapping, Sequence


def _cfg_get(cfg: object | None, key: str, default: object = None) -> object:
    if cfg is None:
        return default
    if isinstance(cfg, Mapping):
        return cfg.get(key, default)
    getter = getattr(cfg, "get", None)
    if callable(getter):
        try:
            return getter(key, default)
        except Exception:
            pass
    return getattr(cfg, key, default)


def resolve_tool(cfg: Mapping[str, Any] | object | None, key: str, fallback_cmd: str) -> dict[str, str]:
    """
    Resolve tool path with fixed prefer-config policy:
    configured value first, then PATH fallback.
    """
    configured = str(_cfg_get(cfg, key, "") or "").strip()
    if configured:
        p = Path(configured).expanduser()
        if p.is_absolute() or "/" in configured or "\\" in configured:
            # Configured executable/script paths must point to a file, not just an existing directory.
            if p.is_file():
                return {
                    "resolved_path": str(p),
                    "source": "config",
                    "configured_value": configured,
                }
            fallback = shutil.which(fallback_cmd)
            if fallback:
                return {
                    "resolved_path": str(fallback),
                    "source": "path_fallback",
                    "configured_value": configured,
                }
            return {
                "resolved_path": "",
                "source": "missing",
                "configured_value": configured,
            }
        found_cfg = shutil.which(configured)
        if found_cfg:
            return {
                "resolved_path": str(found_cfg),
                "source": "config_token",
                "configured_value": configured,
            }
        fallback = shutil.which(fallback_cmd)
        if fallback:
            return {
                "resolved_path": str(fallback),
                "source": "path_fallback",
                "configured_value": configured,
            }
        return {
            "resolved_path": "",
            "source": "missing",
            "configured_value": configured,
        }

    found = shutil.which(fallback_cmd)
    if found:
        return {
            "resolved_path": str(found),
            "source": "path",
            "configured_value": configured,
        }
    return {
        "resolved_path": "",
        "source": "missing",
        "configured_value": configured,
    }


def resolve_required_tools(
    cfg: Mapping[str, Any] | object | None,
    tool_specs: Sequence[tuple[str, str, bool]],
) -> list[tuple[str, str]]:
    missing_required: list[tuple[str, str]] = []
    for key, fallback_cmd, required in tool_specs:
        result = resolve_tool(cfg, key, fallback_cmd)
        if not str(result.get("resolved_path", "")).strip() and required:
            missing_required.append((key, fallback_cmd))
    return missing_required


def looks_like_pkgs_cache(prefix: Path) -> bool:
    return "/micromamba/pkgs" in prefix.as_posix()


def prefix_has_tool(prefix: Path, tool: str) -> bool:
    return (prefix / "bin" / tool).is_file()


def prefix_has_tools(prefix: Path, tools: Sequence[str]) -> bool:
    return all(prefix_has_tool(prefix, tool) for tool in tools)


def resolve_micromamba() -> Path | None:
    """
    Resolve the environment runner executable used for prefix-based launches.

    Preference order:
      1) MICROMAMBA_EXE env
      2) micromamba on PATH
      3) ~/micromamba/bin/micromamba
      4) CONDA_EXE env
      5) conda on PATH

    Note:
      Kept under legacy name for compatibility with existing call sites.
      Returned executable may be either micromamba or conda.
    """
    env_mm = os.environ.get("MICROMAMBA_EXE")
    if env_mm:
        env_path = Path(env_mm).expanduser()
        if env_path.is_file():
            return env_path
    found = shutil.which("micromamba")
    if found:
        return Path(found)
    home_candidate = Path.home() / "micromamba" / "bin" / "micromamba"
    if home_candidate.is_file():
        return home_candidate
    env_conda = os.environ.get("CONDA_EXE")
    if env_conda:
        conda_path = Path(env_conda).expanduser()
        if conda_path.is_file():
            return conda_path
    conda_found = shutil.which("conda")
    if conda_found:
        return Path(conda_found)
    return None


def resolve_conda() -> Path | None:
    env_conda = os.environ.get("CONDA_EXE")
    if env_conda:
        conda_path = Path(env_conda).expanduser()
        if conda_path.is_file():
            return conda_path
    conda_found = shutil.which("conda")
    if conda_found:
        return Path(conda_found)
    return None


def build_prefix_runner(prefix: Path) -> list[str]:
    runner = resolve_micromamba()
    if not runner:
        raise FileNotFoundError("neither micromamba nor conda found")
    return [str(runner), "run", "-p", str(prefix)]


def build_micromamba_runner(prefix: Path) -> list[str]:
    # Backward-compatible name used across the MMGBSA modules.
    return build_prefix_runner(prefix)


def resolve_ambertools_prefix(cfg: object | None) -> str | None:
    env_mmgbsa = os.environ.get("MMGBSA_AMBERTOOLS_PREFIX")
    if env_mmgbsa:
        return env_mmgbsa
    env_prefix = os.environ.get("AMBERTOOLS_PREFIX")
    if env_prefix:
        return env_prefix
    cfg_prefix = _cfg_get(cfg, "MMGBSA_AMBERTOOLS_PREFIX", None)
    if cfg_prefix:
        return str(cfg_prefix)
    fallback = _cfg_get(cfg, "AMBERTOOLS_PREFIX", None)
    if fallback:
        return str(fallback)
    return None


def common_ambertools_prefixes(cfg: object | None = None) -> list[Path]:
    candidates: list[Path] = []
    env_candidates = [
        os.environ.get("MMGBSA_AMBERTOOLS_PREFIX"),
        os.environ.get("AMBERTOOLS_PREFIX"),
        os.environ.get("CONDA_PREFIX"),
    ]
    for raw in env_candidates:
        if raw:
            candidates.append(Path(raw).expanduser())

    overall = _cfg_get(cfg, "OVERALL_DIR", None)
    if overall:
        base = Path(str(overall)).expanduser().resolve().parent
        candidates.append(base / "tools" / "envs" / "ambertools")

    mamba_root = os.environ.get("MAMBA_ROOT_PREFIX")
    if mamba_root:
        candidates.append(Path(mamba_root).expanduser() / "envs" / "AmberTools25")

    home = Path.home()
    candidates.extend(
        [
            home / "tools" / "envs" / "ambertools",
            home / "micromamba" / "envs" / "AmberTools25",
            home / "miniconda3" / "envs" / "AmberTools25",
            home / "anaconda3" / "envs" / "AmberTools25",
        ]
    )

    seen: set[str] = set()
    out: list[Path] = []
    for path in candidates:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        out.append(path)
    return out
