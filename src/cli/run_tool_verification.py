from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from cli.cli_utils import _cli_has
from cli.run_context import ConfigDict
from config.tool_resolver import resolve_micromamba, resolve_required_tools, resolve_tool

_TOOL_VERIFY_TARGETS: tuple[tuple[str, str, bool], ...] = (
    ("VINA_EXE", "vina", True),
    ("OPENBABEL_PATH", "obabel", False),
    ("P2RANK_PATH", "prank", True),
    ("SCORCH", "scorch.py", False),
)

_TOOL_FIX_HINTS: dict[str, str] = {
    "P2RANK_PATH": "Install P2Rank and set P2RANK_PATH to its install root or bin directory.",
    "MEEKO": "Install Meeko in the active docking environment.",
}


def _verify_scorch_env(cfg: ConfigDict) -> tuple[bool, str]:
    runner = resolve_micromamba()
    if not runner:
        return False, "missing_env_runner"
    runner_str = str(runner)

    env_prefix_raw = str(cfg.get("SCORCH_ENV_PREFIX", "") or "").strip()
    env_name = str(cfg.get("SCORCH_ENV", "scorch-env") or "scorch-env").strip()
    if not env_name:
        env_name = "scorch-env"

    if env_prefix_raw:
        env_prefix = str(Path(env_prefix_raw).expanduser().resolve())
        mode = f"prefix={env_prefix} runner={Path(runner_str).name}"
        cmd = [runner_str, "run", "-p", env_prefix, "python", "-V"]
    else:
        mode = f"name={env_name} runner={Path(runner_str).name}"
        cmd = [runner_str, "run", "-n", env_name, "python", "-V"]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except Exception as exc:
        return False, f"env_probe_exception mode={mode} error={exc}"

    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip().splitlines()
        summary = detail[0] if detail else ""
        if summary:
            return (
                False,
                f"env_probe_failed mode={mode} rc={proc.returncode} detail={summary}",
            )
        return False, f"env_probe_failed mode={mode} rc={proc.returncode}"
    return True, f"ok mode={mode}"


def _probe_meeko_runtime() -> tuple[bool, str]:
    try:
        import meeko
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        _ = MoleculePreparation, PDBQTWriterLegacy
        version = str(getattr(meeko, "__version__", "unknown"))
        return True, f"MEEKO ok version={version}"
    except Exception as exc:
        return False, f"MEEKO import probe failed err={exc}"


def _verify_tools_if_requested(cfg: ConfigDict, argv: list[str]) -> None:
    """
    Run tool verification only when explicitly requested by CLI flag.
    Config defaults remain false and do not auto-enable verification.
    """
    cfg.setdefault("TOOL_VERIFY_ON_START", False)
    cli_verify = _cli_has(argv, "--verify-tools")
    cfg["TOOL_VERIFY_ON_START"] = bool(cli_verify)
    if not cli_verify:
        return

    missing_required = resolve_required_tools(cfg, _TOOL_VERIFY_TARGETS)
    missing_optional: list[tuple[str, str]] = []
    scorch_tool_available = False
    for key, fallback_cmd, required in _TOOL_VERIFY_TARGETS:
        if required:
            continue
        result = resolve_tool(cfg, key, fallback_cmd)
        resolved_path = str(result.get("resolved_path", "")).strip()
        if key == "SCORCH":
            scorch_tool_available = bool(resolved_path)
        if not resolved_path:
            missing_optional.append((key, fallback_cmd))
    if scorch_tool_available:
        scorch_env_ok, scorch_env_reason = _verify_scorch_env(cfg)
        if not scorch_env_ok:
            missing_optional.append(("SCORCH_ENV", scorch_env_reason))

    if missing_required:
        print(f"{len(missing_required)} required tool checks failed", file=sys.stderr)
        for key, fallback_cmd in missing_required:
            result = resolve_tool(cfg, key, fallback_cmd)
            configured = str(result.get("configured_value", "") or "").strip() or "<unset>"
            source = str(result.get("source", "missing") or "missing")
            fix = _TOOL_FIX_HINTS.get(
                key,
                f"Set {key} in config.txt or install `{fallback_cmd}` on PATH.",
            )
            print(
                f"- {key}: failed to resolve `{fallback_cmd}` "
                f"(configured={configured}, source={source}). {fix} "
                "Re-run `atlas --verify-tools` and `atlas --doctor`.",
                file=sys.stderr,
            )
        sys.exit(2)
    runtime_failures: list[str] = []
    runtime_details: list[str] = []
    meeko_ok, meeko_reason = _probe_meeko_runtime()
    if meeko_ok:
        runtime_details.append(meeko_reason)
    else:
        runtime_failures.append(meeko_reason)
    if runtime_failures:
        print(f"{len(runtime_failures)} runtime probe failures", file=sys.stderr)
        for reason in runtime_failures:
            print(
                f"- {reason}. Fix the tool installation/path and re-run "
                "`atlas --verify-tools` (details) or `atlas --doctor` (overview).",
                file=sys.stderr,
            )
        sys.exit(2)
    print("Tools all successfully verified")
    for detail in runtime_details:
        print(f"[verify-tools] {detail}")
    if missing_optional:
        print(f"{len(missing_optional)} optional missing")
        for key, fallback_cmd in missing_optional:
            print(f"- {key} ({fallback_cmd})")
