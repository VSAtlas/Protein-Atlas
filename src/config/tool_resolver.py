from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from config.value_access import cfg_get


@dataclass(frozen=True)
class ToolSelection:
    runner: list[str]
    exe: str
    source: str
    prefix: Path | None = None


def resolve_tool(cfg: Mapping[str, Any] | object | None, key: str, fallback_cmd: str) -> dict[str, str]:
    """
    Resolve tool path with fixed prefer-config policy:
    configured value first, then PATH fallback.
    """
    configured = str(cfg_get(cfg, key, "") or "").strip()
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
      2) MAMBA_EXE env
      3) micromamba on PATH
      4) ~/micromamba/bin/micromamba
      5) CONDA_EXE env
      6) conda on PATH

    Note:
      Kept under legacy name for compatibility with existing call sites.
      Returned executable may be either micromamba or conda.
    """
    env_mm = os.environ.get("MICROMAMBA_EXE")
    if env_mm:
        env_path = Path(env_mm).expanduser()
        if env_path.is_file():
            return env_path
    env_mamba = os.environ.get("MAMBA_EXE")
    if env_mamba:
        mamba_path = Path(env_mamba).expanduser()
        if mamba_path.is_file():
            return mamba_path
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


def repo_ambertools_roots(cfg: object | None) -> list[Path]:
    """
    Canonical repo-relative AmberTools prefixes: UNDER OVERALL_DIR parents and
    ``~/tools/envs/{ambertools-mpi,ambertools}`` (used by normalization and resolver).
    """
    roots: list[Path] = []
    overall_raw = cfg_get(cfg, "OVERALL_DIR", None)
    overall_str = str(overall_raw).strip() if overall_raw is not None else ""
    if overall_str:
        base = Path(overall_str).expanduser().resolve()
        for parent in (base, *base.parents):
            roots.append(parent / "tools" / "envs" / "ambertools-mpi")
            roots.append(parent / "tools" / "envs" / "ambertools")
    home = Path.home()
    roots.extend(
        [
            home / "tools" / "envs" / "ambertools-mpi",
            home / "tools" / "envs" / "ambertools",
        ]
    )
    return roots


def resolve_ambertools_prefix(cfg: object | None) -> str | None:
    env_mmgbsa = os.environ.get("MMGBSA_AMBERTOOLS_PREFIX")
    if env_mmgbsa:
        return env_mmgbsa
    env_prefix = os.environ.get("AMBERTOOLS_PREFIX")
    if env_prefix:
        return env_prefix
    cfg_prefix = cfg_get(cfg, "MMGBSA_AMBERTOOLS_PREFIX", None)
    if cfg_prefix:
        return str(cfg_prefix)
    fallback = cfg_get(cfg, "AMBERTOOLS_PREFIX", None)
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

    candidates.extend(repo_ambertools_roots(cfg))

    mamba_root = os.environ.get("MAMBA_ROOT_PREFIX")
    if mamba_root:
        candidates.append(Path(mamba_root).expanduser() / "envs" / "AmberTools25")

    home = Path.home()
    candidates.extend(
        [
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


def select_ambertools_tool(
    cfg: object | None,
    tool: str,
    *,
    explicit_prefix: str | Path | None = None,
    strict_explicit_prefix: bool = False,
    prefix_exe_style: str = "name",
) -> ToolSelection:
    prefix_value = explicit_prefix or resolve_ambertools_prefix(cfg)
    if prefix_value:
        prefix = Path(prefix_value).expanduser()
        if not looks_like_pkgs_cache(prefix) and prefix_has_tool(prefix, tool):
            exe = str(prefix / "bin" / tool) if prefix_exe_style == "absolute" else tool
            return ToolSelection(
                runner=build_micromamba_runner(prefix),
                exe=exe,
                source=f"prefix:{prefix}",
                prefix=prefix,
            )
        if strict_explicit_prefix:
            raise FileNotFoundError(f"could not locate {tool} under {prefix}")

    path_tool = shutil.which(tool)
    if path_tool:
        return ToolSelection(runner=[], exe=path_tool, source="PATH")

    for prefix in common_ambertools_prefixes(cfg):
        if not looks_like_pkgs_cache(prefix) and prefix_has_tool(prefix, tool):
            exe = str(prefix / "bin" / tool) if prefix_exe_style == "absolute" else tool
            return ToolSelection(
                runner=build_micromamba_runner(prefix),
                exe=exe,
                source=f"prefix:{prefix}",
                prefix=prefix,
            )

    raise FileNotFoundError(f"could not locate {tool}; set AMBERTOOLS_PREFIX or PATH")


def select_ambertools_tools(
    cfg: object | None,
    tools: Sequence[str],
    *,
    explicit_prefix: str | Path | None = None,
    strict_explicit_prefix: bool = False,
    prefix_exe_style: str = "name",
) -> tuple[list[str], dict[str, str], str]:
    prefix_value = explicit_prefix or resolve_ambertools_prefix(cfg)
    if prefix_value:
        prefix = Path(prefix_value).expanduser()
        if not looks_like_pkgs_cache(prefix) and prefix_has_tools(prefix, tools):
            executables = {
                tool: str(prefix / "bin" / tool) if prefix_exe_style == "absolute" else tool
                for tool in tools
            }
            return build_micromamba_runner(prefix), executables, f"prefix:{prefix}"
        if strict_explicit_prefix:
            missing = [tool for tool in tools if not prefix_has_tool(prefix, tool)]
            raise FileNotFoundError(
                f"could not locate required AmberTools executables under {prefix}: {','.join(missing)}"
            )

    path_tools = {tool: shutil.which(tool) for tool in tools}
    if all(path_tools.values()):
        return [], {tool: str(path_tools[tool]) for tool in tools}, "PATH"

    for prefix in common_ambertools_prefixes(cfg):
        if not looks_like_pkgs_cache(prefix) and prefix_has_tools(prefix, tools):
            executables = {
                tool: str(prefix / "bin" / tool) if prefix_exe_style == "absolute" else tool
                for tool in tools
            }
            return build_micromamba_runner(prefix), executables, f"prefix:{prefix}"

    raise FileNotFoundError(
        f"could not locate required AmberTools executables: {','.join(tools)}"
    )
